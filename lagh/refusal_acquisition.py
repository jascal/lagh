"""Refusal-driven measurement acquisition on positive deterministic domains.

All discovery data are DESIGN data. Final observations are acquired once, after
freezing a candidate, and never inform a subsequent design. Experimental policy
is explicit; the existing acquisition baseline and hypothesis class are intact.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
import hashlib

import numpy as np
import sympy as sp

from .acquisition import _box_ladder, run_active_boxsearch
from .certify import Abstain, Certificate, check, epsilon
from .engine import Result
from .measurement_design import choose_measurement
from .passive import discover_passive


@dataclass(frozen=True)
class RefusalPolicy:
    admissible: tuple[tuple[float, ...], tuple[float, ...]]
    batch: int = 80
    final_points: int = 80
    max_rounds: int = 5
    max_queries: int = 1800
    probe_points: int = 257
    box_factor: float = 3.0
    floor_abs: float = 1e-12


@dataclass
class AcquisitionOutcome:
    result: Result
    queries: int
    initial_structural: bool
    initial_rivals: int
    history: list = field(default_factory=list)
    domain: dict | None = None
    final_passed: bool | None = None
    native_verdict: bool | None = None
    # Internal evidence for scoring/reproduction, never supplied to the chooser.
    X_design: np.ndarray | None = None
    y_design: np.ndarray | None = None
    X_final: np.ndarray | None = None
    y_final: np.ndarray | None = None

    @property
    def certified(self):
        return bool(self.result.certificate.certified and self.final_passed is True)


def fingerprint(*arrays):
    h = hashlib.sha256()
    for value in arrays:
        a = np.ascontiguousarray(value, dtype='<f8')
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def _sample(box, n, rng):
    return np.exp(rng.uniform(np.log(box[0]), np.log(box[1]), (n, box.shape[1])))


def _inside(box, admissible):
    return bool(np.all(box[0] >= admissible[0]) and np.all(box[1] <= admissible[1]))


def run_refusal_acquisition(oracle, X, y, *, initial_box, policy: RefusalPolicy,
                            seed: int, strategy: str = 'guided', initial_result=None,
                            on_plan=None) -> AcquisitionOutcome:
    """Count scalar oracle observations, including supplied initial rows.

    `initial_result`, if supplied, is a cached ordinary discover_passive result
    on exactly X,y. It is design evidence, never an accepted final certificate.
    `on_plan` persists each immutable-by-copy query plan BEFORE its oracle call.
    Native-ladder is an unchanged selection baseline plus the common final guard.
    Other fixed arms use the same batch/discovery machinery as guided.
    """
    if strategy not in ('guided', 'fixed-x10', 'fixed-ladder', 'native-ladder'):
        raise ValueError('unknown acquisition strategy')
    X, y = np.asarray(X, float).copy(), np.asarray(y, float).copy()
    admissible = np.asarray(policy.admissible, float)
    declared = np.asarray(initial_box, float)
    if (X.ndim != 2 or not len(X) or not X.shape[1] or y.shape != (len(X),)
            or not np.all(np.isfinite(X)) or not np.all(np.isfinite(y))
            or np.any(X <= 0) or admissible.shape != (2, X.shape[1])
            or declared.shape != admissible.shape
            or not np.all(np.isfinite(admissible)) or np.any(admissible <= 0)
            or np.any(admissible[0] >= admissible[1])
            or not np.all(np.isfinite(declared)) or np.any(declared <= 0)
            or np.any(declared[0] >= declared[1])
            or not _inside(declared, admissible)
            or np.any(X < declared[0]) or np.any(X > declared[1])):
        raise ValueError('finite initial observations and positive admissible boxes required')
    for n in (policy.batch, policy.final_points, policy.max_rounds,
              policy.max_queries, policy.probe_points):
        if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n <= 0:
            raise ValueError('positive integer sampling budgets required')
    if (not np.isfinite(policy.box_factor) or policy.box_factor <= 1
            or not np.isfinite(policy.floor_abs) or policy.floor_abs < 0):
        raise ValueError('finite box factor >1 and nonnegative floor required')
    if len(X) > policy.max_queries:
        raise ValueError('initial data already exceed declared query budget')
    design_seed, final_seed, probe_seed = np.random.SeedSequence(seed).spawn(3)
    rng = np.random.default_rng(design_seed)
    probe_rng = np.random.default_rng(probe_seed)
    if initial_result is None:
        initial_result = discover_passive(X, y, sigma=0., floor_abs=policy.floor_abs)
    r = deepcopy(initial_result.result)
    outcome = AcquisitionOutcome(r, len(X),
                                 r.certificate.abstain == Abstain.STRUCTURAL.value,
                                 len(r.rivals))
    syms = [sp.Symbol(f'x_{i}') for i in range(X.shape[1])]
    design_X, design_y = [X], [y]

    def demote(reason, note):
        c = deepcopy(outcome.result.certificate)
        c.certified, c.abstain = False, reason
        c.notes.append(note)
        outcome.result = Result(c, None, outcome.result.tier,
                                outcome.result.n_candidates, outcome.result.rivals)

    def query(points, role, **metadata):
        n = len(points)
        reserve = 0 if role == 'final' else policy.final_points
        if outcome.queries + n + reserve > policy.max_queries:
            raise ValueError('query would exceed budget including reserved final sample')
        plan = {'role': role, 'query_start': outcome.queries, 'count': n,
                'X': points.tolist(), 'X_hash': fingerprint(points), **metadata}
        # Persist before measurement. A consumer receives a copy, so it cannot
        # rewrite what the oracle is about to be asked.
        if on_plan is not None:
            on_plan(deepcopy(plan))
        outcome.history.append(plan)
        outcome.queries += n
        values = np.asarray(oracle(points.copy()), float)
        if values.shape != (n,) or not np.all(np.isfinite(values)):
            raise ValueError('oracle returned invalid observations; no rows may be dropped')
        outcome.history[-1]['y_hash'] = fingerprint(values)
        return values

    def finish_candidate():
        # Freeze before even GENERATING final inputs. No later design follows
        # either outcome. The old samples and the new samples must both pass.
        frozen = deepcopy(outcome.result)
        all_X, all_y = np.vstack(design_X), np.concatenate(design_y)
        old = check(frozen.expr, syms, all_X, all_y,
                    epsilon(all_y, floor_abs=policy.floor_abs))
        if not old['certified']:
            demote(Abstain.HELDOUT.value, 'candidate fails previously measured design rows')
            outcome.final_passed = False
            return
        box = np.array([all_X.min(axis=0), all_X.max(axis=0)])
        final_rng = np.random.default_rng(final_seed)
        fresh_X = _sample(box, policy.final_points, final_rng)
        outcome.X_final = fresh_X.copy()
        try:
            fresh_y = query(fresh_X, 'final', frozen_law=str(frozen.expr),
                            sampled_bounds=box.tolist())
        except (ValueError, FloatingPointError) as error:
            demote(Abstain.HELDOUT.value, str(error))
            outcome.final_passed = False
            return
        outcome.y_final = fresh_y.copy()
        final = check(frozen.expr, syms, fresh_X, fresh_y,
                      epsilon(fresh_y, floor_abs=policy.floor_abs))
        outcome.final_passed = bool(final['certified'])
        if not outcome.final_passed:
            demote(Abstain.HELDOUT.value, 'frozen candidate fails fresh final observations')
            return
        domain_X = np.vstack([all_X, fresh_X])
        domain = {'kind': 'finite observed rows, not the continuous bounding box',
                  'bounds': [domain_X.min(axis=0).tolist(), domain_X.max(axis=0).tolist()],
                  'rows': len(domain_X),
                  'design_hash': fingerprint(all_X, all_y),
                  'final_hash': fingerprint(fresh_X, fresh_y)}
        frozen.certificate.state_bounds = list(map(list, zip(*domain['bounds'])))
        frozen.certificate.domain_size = len(domain_X)
        frozen.certificate.notes.append(
            'acquisition: finite domain expanded to all prior and fresh rows; '
            'alpha is the original discovery statistic, not an adaptive family-wide bound')
        outcome.result, outcome.domain = frozen, domain

    if r.certificate.certified:
        finish_candidate()
    elif strategy == 'native-ladder':
        # Only an admissible PREFIX, preserving native ordering and boxes.
        prefix = 0
        for _, lo, hi in _box_ladder(*declared):
            if not _inside(np.array([lo, hi]), admissible):
                break
            prefix += 1
        def native_oracle(points):
            values = query(points, 'native-design')
            design_X.append(np.array(points, copy=True))
            design_y.append(values.copy())
            return values
        try:
            native = run_active_boxsearch(native_oracle, *declared, budget=200,
                                          max_boxes=prefix, floor_abs=policy.floor_abs,
                                          seed=seed, time_budget_s=None)
            outcome.result = native.active.result
            outcome.native_verdict = bool(outcome.result.certificate.certified)
            outcome.history.append({'role': 'native-summary', 'boxes': native.transforms,
                                    'native_holdout': native.heldout_box_ok})
            if outcome.native_verdict:
                finish_candidate()
        except (ValueError, FloatingPointError) as error:
            demote(Abstain.HELDOUT.value, str(error))
    else:
        ladder = list(_box_ladder(*declared))
        for step in range(policy.max_rounds):
            if outcome.queries + policy.batch + policy.final_points > policy.max_queries:
                break
            all_X = np.vstack(design_X)
            if strategy == 'guided':
                if X.shape[1] == 1:
                    probes = np.geomspace(*admissible[:, 0], policy.probe_points)[:, None]
                else:
                    probes = _sample(admissible, policy.probe_points, probe_rng)
                choice = choose_measurement(outcome.result.rivals, probes, bounds=admissible)
                if choice.point is None:
                    outcome.history.append({'role': 'stop', 'reason': choice.reason})
                    break
                q = np.asarray(choice.point)
                box = np.array([np.maximum(admissible[0], q/policy.box_factor),
                                np.minimum(admissible[1], q*policy.box_factor)])
                metadata = {'choice': asdict(choice),
                            'rivals': [str(e) for e in outcome.result.rivals]}
            elif strategy == 'fixed-x10':
                box = np.array([all_X.min(axis=0)/10, all_X.max(axis=0)*10])
                metadata = {'fixed': 'recover suggested_box'}
            else:
                if step >= len(ladder):
                    break
                name, lo, hi = ladder[step]
                box = np.array([lo, hi])
                metadata = {'fixed': name}
            if not _inside(box, admissible):
                outcome.history.append({'role': 'stop', 'reason': 'next fixed box is inadmissible',
                                        'box': box.tolist()})
                break
            new_X = _sample(box, policy.batch, rng)
            if strategy == 'guided':
                new_X[0] = q
            try:
                new_y = query(new_X, 'design', step=step, sampled_bounds=box.tolist(), **metadata)
            except (ValueError, FloatingPointError) as error:
                demote(Abstain.HELDOUT.value, str(error))
                break
            design_X.append(new_X.copy())
            design_y.append(new_y.copy())
            result = discover_passive(np.vstack(design_X), np.concatenate(design_y),
                                      sigma=0., floor_abs=policy.floor_abs)
            outcome.result = result.result
            outcome.history.append({'role': 'discovery', 'certified': result.certified,
                                    'abstain': result.result.certificate.abstain,
                                    'rivals': [str(e) for e in result.result.rivals],
                                    'hypotheses': result.result.n_candidates})
            if result.certified:
                finish_candidate()
                break
    outcome.X_design, outcome.y_design = np.vstack(design_X), np.concatenate(design_y)
    # The public outcome never leaves a failed/not-run final guard as certified.
    if outcome.result.certificate.certified and outcome.final_passed is not True:
        demote(Abstain.HELDOUT.value, 'no passed independent final guard')
    return outcome
