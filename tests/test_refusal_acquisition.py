"""Gate counter-inputs for acquisition, using mocked discovery (no reach claim)."""
from copy import deepcopy

import numpy as np
import pytest
import sympy as sp

import lagh.refusal_acquisition as acq
from lagh.certify import Certificate
from lagh.engine import Result
from lagh.passive import PassiveResult

x = sp.Symbol('x_0')
X = np.linspace(.5, 3., 40)[:, None]
Y = X[:, 0]
POLICY = acq.RefusalPolicy(((.005,), (300.,)), batch=10, final_points=12,
                          max_rounds=3, max_queries=100)


def result(expr=None):
    if expr is None:
        c = Certificate(False, 0, 0, 8, [[.5, 3.]], '', abstain='structural')
        r = Result(c, None, 1, 2, (x, 2*x))
    else:
        c = Certificate(True, 0, 0, 8, [[.5, 3.]], str(expr), alpha_log10=-100)
        r = Result(c, expr, 1, 2)
    return PassiveResult(r, 1, None if expr is None else True)


def run(oracle, **kwargs):
    return acq.run_refusal_acquisition(oracle, X, Y, initial_box=[[.5], [3.]],
                                       policy=kwargs.pop('policy', POLICY), seed=7,
                                       initial_result=result(), **kwargs)


def test_expanded_domain_is_recomputed(monkeypatch):
    # Falsifies inheriting the original [.5,3] bounds after measuring at300.
    monkeypatch.setattr(acq, 'discover_passive', lambda *a, **kw: result(x))
    out = run(lambda X: X[:, 0])
    assert out.certified and out.queries == 62
    assert out.domain['bounds'][1] == [300.]
    assert out.result.certificate.state_bounds[0][1] == 300.
    assert out.result.certificate.domain_size == 62
    assert out.domain['rows'] == 62
    assert acq.fingerprint(out.X_final) != acq.fingerprint(out.X_design)


def test_final_disagreement_demotes_and_stops(monkeypatch):
    # Frozen x agrees with design, but fresh observations are x+1. No retry.
    monkeypatch.setattr(acq, 'discover_passive', lambda *a, **kw: result(x))
    calls = []
    def oracle(points):
        calls.append(len(points))
        return points[:, 0] + (len(calls) == 2)
    out = run(oracle)
    assert calls == [10, 12]
    assert not out.certified and not out.result.certificate.certified
    assert out.result.expr is None and out.result.certificate.abstain == 'heldout'
    assert out.final_passed is False


def test_future_responses_do_not_change_prequery_design(monkeypatch):
    # Two futures, identical initial/design data. Comparing plans INCLUDING
    # final inputs catches selection contaminated by future responses.
    monkeypatch.setattr(acq, 'discover_passive', lambda *a, **kw: result(x))
    plans = []
    for delta in (0., 100.):
        seen, records = [], []
        def oracle(points):
            seen.append(len(points))
            return points[:, 0] + (delta if len(seen) == 2 else 0.)
        run(oracle, on_plan=records.append)
        plans.append(deepcopy(records))
    assert plans[0] == plans[1]
    assert [r['role'] for r in plans[0]] == ['design', 'final']


@pytest.mark.parametrize('bad', [np.nan, np.inf, -np.inf])
def test_invalid_design_observation_cannot_be_filtered(monkeypatch, bad):
    # Constant candidate used to make NaN input witnesses too specific.
    # Here invalid output must be caught BEFORE any discovery/filtering call.
    monkeypatch.setattr(acq, 'discover_passive', lambda *a, **kw:
                        pytest.fail('invalid evidence reached discovery'))
    def oracle(points):
        y = points[:, 0].copy()
        y[3] = bad
        return y
    out = run(oracle)
    assert not out.certified and out.queries == 50
    assert len(out.history) == 1


@pytest.mark.parametrize('bad', [np.nan, np.inf])
def test_invalid_final_observation_demotes(monkeypatch, bad):
    monkeypatch.setattr(acq, 'discover_passive', lambda *a, **kw: result(x))
    calls = []
    def oracle(points):
        calls.append(len(points))
        y = points[:, 0].copy()
        if len(calls) == 2:
            y[0] = bad
        return y
    out = run(oracle)
    assert not out.certified and out.final_passed is False
    assert calls == [10, 12]


def test_old_observations_cannot_be_erased(monkeypatch):
    # A candidate x+1 matches new box data but contradicts all original rows.
    monkeypatch.setattr(acq, 'discover_passive', lambda *a, **kw: result(x+1))
    calls = []
    def oracle(points):
        calls.append(len(points))
        return points[:, 0] + 1
    out = run(oracle)
    assert calls == [10]     # stopped before final, not accepted on new box alone
    assert not out.certified and out.final_passed is False


def test_budget_reserves_final_measurements():
    # 61 leaves10 design +11 final slots;12 are required, so no query is legal.
    policy = acq.RefusalPolicy(((.005,), (300.,)), batch=10, final_points=12,
                               max_queries=61)
    out = run(lambda _: pytest.fail('budget exceeded'), policy=policy)
    assert out.queries == 40 and not out.certified


def test_no_rivals_cannot_trigger_fallback_queries():
    initial = result()
    initial.result.rivals = ()
    out = acq.run_refusal_acquisition(lambda _: pytest.fail('no informative query'),
                                      X, Y, initial_box=[[.5], [3.]],
                                      policy=POLICY, seed=7, initial_result=initial)
    assert out.queries == 40 and not out.certified
