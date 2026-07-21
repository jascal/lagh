"""P2: active acquisition -- adaptive ranging, ledgers, and the two adopted
mechanisms (per-round micro-predictions; multi-objective utility with a
non-smuggled reference signal).

The policy is FROZEN at instrument registration; a run supplies only the target
declaration (oracle, box, noise declaration, budget). No per-target lever exists.

The two mechanisms, as adopted:

* MICRO-PREDICTIONS: before each active query batch, the current winner (and the
  best lower-tier law) predict the batch's values; predictions are recorded, then
  scored the moment observations arrive. Every round is a falsifiable
  mini-experiment -- the per-run predictions discipline moved into the loop.

* MULTI-OBJECTIVE UTILITY: query points scored by registered weights over
  (a) disagreement among surviving coherent classes, and (b) divergence of the
  winner from the BEST LOWER-TIER law -- an internal reference that injects no
  target knowledge. External "conjecture models" were rejected at design review:
  a supplied theory is per-target knowledge, and any recovery conditioned on it
  is not discovery.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp

from .base import eval_expr
from .certify import Abstain, Certificate, epsilon, sample_box
from .engine import Result, discover

MACHINE_FLOOR = 1e-12


@dataclass(frozen=True)
class Policy:
    """Frozen at instrument registration."""

    init_points: int = 40
    replicates: int = 2
    max_rounds: int = 8
    batch: int = 20
    ranging_iterations: int = 3
    ranging_probe: int = 16
    floor_fraction: float = 0.25      # >25% of outputs at floor -> contract
    target_signal: float = 0.75
    w_disagreement: float = 0.6
    w_reference: float = 0.4
    stall_rounds: int = 2


@dataclass
class Ledger:
    entries: list = field(default_factory=list)

    def record(self, kind: str, cost: int, **meta):
        self.entries.append({"kind": kind, "cost": cost, **meta})

    @property
    def spent(self) -> int:
        return sum(e["cost"] for e in self.entries)


@dataclass
class ActiveResult:
    result: Result
    box_initial: np.ndarray
    box_final: np.ndarray
    ranging_trajectory: list
    ledger: Ledger
    predictions: list            # per-round micro-prediction reports
    queries_used: int


def _sample(box: np.ndarray, n: int, rng) -> np.ndarray:
    lo, hi = box[0], box[1]
    return np.exp(rng.uniform(np.log(lo), np.log(hi), (n, len(lo))))


def _signal_profile(X: np.ndarray, y: np.ndarray, floor: float) -> dict:
    above = np.abs(y) > floor
    prof = {"fraction_above_floor": float(np.mean(above)), "octave_range": []}
    for j in range(X.shape[1]):
        # output dynamic range per input octave: which inputs kill the signal?
        lo, hi = X[:, j].min(), X[:, j].max()
        mid = np.sqrt(lo * hi)
        lo_half = above[X[:, j] <= mid]
        hi_half = above[X[:, j] > mid]
        prof["octave_range"].append(
            (float(np.mean(lo_half)) if lo_half.size else 0.0,
             float(np.mean(hi_half)) if hi_half.size else 0.0))
    return prof


def _contract(box: np.ndarray, prof: dict) -> np.ndarray:
    """Shrink along the inputs whose half kills the signal, toward the live half."""
    box = box.copy()
    for j, (lo_frac, hi_frac) in enumerate(prof["octave_range"]):
        mid = np.sqrt(box[0][j] * box[1][j])
        if hi_frac < lo_frac:          # upper half is dead
            box[1][j] = mid
        elif lo_frac < hi_frac:        # lower half is dead
            box[0][j] = mid
    return box


def run_active(oracle, box_lo, box_hi, *, budget: int = 200,
               policy: Policy = Policy(), sigma_declared: float | None = None,
               floor_abs: float = MACHINE_FLOOR, seed: int = 0) -> ActiveResult:
    """oracle: X (n,d) -> y (n,). Deterministic or noisy; replicates estimate sigma
    when it is not declared."""
    rng = np.random.default_rng(seed)
    box0 = np.array([np.asarray(box_lo, float), np.asarray(box_hi, float)])
    box = box0.copy()
    led = Ledger()
    trajectory = [box0.tolist()]

    # ---- adaptive ranging ----
    for it in range(policy.ranging_iterations):
        Xp = _sample(box, policy.ranging_probe, rng)
        yp = np.asarray(oracle(Xp), float)
        led.record("ranging", len(Xp), iteration=it)
        prof = _signal_profile(Xp, yp, floor_abs)
        if prof["fraction_above_floor"] >= 1 - policy.floor_fraction:
            break
        box = _contract(box, prof)
        trajectory.append(box.tolist())
    else:
        Xp = _sample(box, policy.ranging_probe, rng)
        yp = np.asarray(oracle(Xp), float)
        led.record("ranging", len(Xp), iteration="final-check")
        if _signal_profile(Xp, yp, floor_abs)["fraction_above_floor"] \
                < policy.target_signal:
            cert = Certificate(False, 0, 0, 0, box.tolist(), "",
                               abstain=Abstain.RANGE.value,
                               notes=["no sub-box with sufficient signal within "
                                      "the ranging budget", f"trajectory={trajectory}"])
            return ActiveResult(Result(cert, None, 0, 0), box0, box, trajectory,
                                led, [], led.spent)

    # ---- init + replicates ----
    X = _sample(box, policy.init_points, rng)
    y = np.asarray(oracle(X), float)
    led.record("init", len(X))
    rep = np.asarray(oracle(np.vstack([X[0]] * policy.replicates)), float)
    led.record("replicates", policy.replicates)
    sigma = (sigma_declared if sigma_declared is not None
             else float(np.std(np.concatenate([rep, [y[0]]]))
                        / max(abs(y[0]), 1e-30)))

    dim = X.shape[1]
    syms = list(sp.symbols([f"x_{i}" for i in range(dim)])) if dim > 1 \
        else [sp.Symbol("x_0")]
    predictions: list = []
    last_result: Result | None = None
    stall = 0

    for rnd in range(policy.max_rounds):
        idx = rng.permutation(len(X))
        a, b = int(0.6 * len(X)), int(0.8 * len(X))
        r = discover(X[idx[:a]], y[idx[:a]], X[idx[a:b]], y[idx[a:b]],
                     X[idx[b:]], y[idx[b:]], sigma=sigma, floor_abs=floor_abs)
        last_result = r
        if r.certificate.certified:
            return ActiveResult(r, box0, box, trajectory, led, predictions,
                                led.spent)
        if led.spent + policy.batch > budget:
            break

        # ---- multi-objective query selection ----
        cand_pts = _sample(box, 40 * policy.batch, rng)
        util = np.zeros(len(cand_pts))
        # (a) internal disagreement: winner-vs-best-failing spread is unavailable
        #     without certified rivals; use spread among top candidate laws when
        #     the abstention was structural (classes recorded in notes) -- else 0.
        # (b) reference divergence: winner at this tier vs best lower-tier law.
        exprs = []
        if r.expr is not None:
            exprs.append(r.expr)
        if r.certificate.law:
            try:
                exprs.append(sp.sympify(r.certificate.law))
            except Exception:                                 # noqa: BLE001
                pass
        lower = discover(X[idx[:a]], y[idx[:a]], X[idx[a:b]], y[idx[a:b]],
                         X[idx[b:]], y[idx[b:]], sigma=sigma,
                         floor_abs=floor_abs, max_tier=max(1, r.tier - 1))
        if lower.expr is not None:
            exprs.append(lower.expr)
        vals = [eval_expr(e, syms, cand_pts) for e in exprs]
        vals = [v for v in vals if v is not None]
        if len(vals) >= 2:
            V = np.vstack(vals)
            with np.errstate(all="ignore"):
                util = (policy.w_disagreement + policy.w_reference) * \
                       (np.nanmax(V, axis=0) - np.nanmin(V, axis=0))
            util[~np.isfinite(util)] = 0.0
        if util.max() <= 0:
            nxt = _sample(box, policy.batch, rng)      # stall fallback
            stall += 1
        else:
            nxt = cand_pts[np.argsort(-util)[:policy.batch]]
            stall = 0

        # ---- micro-predictions: registered BEFORE the queries are made ----
        pred_record = {"round": rnd, "predictions": []}
        for tag, e in [("winner", exprs[0] if exprs else None),
                       ("lower-tier", lower.expr)]:
            if e is None:
                continue
            pv = eval_expr(e, syms, nxt)
            if pv is not None:
                pred_record["predictions"].append(
                    {"law": tag, "expr": str(e)[:80], "values": pv.tolist()})
        ny = np.asarray(oracle(nxt), float)
        led.record("active", len(nxt), round=rnd)
        # score on arrival
        eps_b = epsilon(ny, sigma=sigma, floor_abs=floor_abs)
        for p in pred_record["predictions"]:
            hits = np.abs(np.asarray(p["values"]) - ny) <= eps_b
            p["hit_rate"] = float(np.mean(hits))
            del p["values"]
        predictions.append(pred_record)

        X = np.vstack([X, nxt])
        y = np.concatenate([y, ny])
        if stall > policy.stall_rounds:
            break

    return ActiveResult(last_result, box0, box, trajectory, led, predictions,
                        led.spent)
