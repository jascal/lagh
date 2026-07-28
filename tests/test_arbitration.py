"""Significance arbitration at the Müntz boundary (MUNTZ_ARBITRATION.md).

The rival shapes here are the MEASURED ones from the reach audit, not
hand-imagined whales: the rivals arbitration may dismiss are 60+ term
interpolations whose free parameters consume the whole certification sample
(h/n ~ 0.01), and the rivals it must NOT dismiss retain roughly half the
sample as held-out evidence (rational-d1: dof 34 vs 45 of n = 80).
"""
import numpy as np
import sympy as sp

from lagh.base import Candidate
from lagh.certify import arbitrate_significance, epsilon, free_dof


def _c(s, cx):
    return Candidate(expr=sp.sympify(s), complexity=cx, channel="linear")


def _classes(*cands):
    return [(c.expr, [c]) for c in cands]


def _dense(n_terms, seed):
    """A dense fractional-power fit: 2 free numbers per term (coefficient +
    rational exponent), the shape the linear channel proposes at the Müntz
    boundary."""
    x = sp.Symbol("x_0")
    rng = np.random.default_rng(seed)
    e = sp.Integer(0)
    for k, c in enumerate(rng.uniform(-1, 1, n_terms)):
        e += sp.Float(round(float(c), 6)) * x ** sp.Rational(k + 3, k + 2)
    return e


Y = np.random.default_rng(0).uniform(1, 10, 80)
EPS = epsilon(Y, sigma=0.0)


def test_muntz_contest_resolves_to_low_dof():
    """The measured winning contest: a 2-dof law against an interpolation with
    no held-out evidence left (mixed-4term-d2 / sparse6-d2)."""
    truth = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    whale = Candidate(expr=_dense(40, 1), complexity=113, channel="linear")
    assert free_dof(whale.expr) >= 0.9 * len(Y)
    r = arbitrate_significance(_classes(truth, whale), Y, EPS, 500)
    assert r is not None
    assert r[0][1][0].expr == truth.expr
    assert "domain claim" in r[1] and "evidence-starved" in r[1]


def test_constrained_rival_blocks_arbitration():
    """The rational-d1 counterexample that amended the rule: the margin is
    astronomical (the rival carries far more dof), but the rival still retains
    ~half the sample as held-out evidence -- a law the data constrains, so the
    contest is genuine and the abstain must stand."""
    truth = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    rival = Candidate(expr=_dense(16, 2), complexity=98, channel="linear")
    h = 1 - free_dof(rival.expr) / len(Y)
    assert 0.3 < h < 0.7                      # genuinely constrained, as measured
    assert arbitrate_significance(_classes(truth, rival), Y, EPS, 500) is None


def test_marginal_contest_stays_abstained():
    a = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    b = _c("3*x_0**2 - 1", 4)
    assert arbitrate_significance(_classes(a, b), Y, EPS, 500) is None


def test_single_class_is_noop():
    a = _c("x_0 + 1", 2)
    assert arbitrate_significance(_classes(a), Y, EPS, 500) is None


def test_two_whales_no_winner():
    w1 = Candidate(expr=_dense(40, 3), complexity=113, channel="linear")
    w2 = Candidate(expr=_dense(39, 4), complexity=110, channel="linear")
    assert arbitrate_significance(_classes(w1, w2), Y, EPS, 500) is None
