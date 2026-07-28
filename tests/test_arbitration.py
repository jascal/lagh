"""Significance arbitration at the Müntz boundary (MUNTZ_ARBITRATION.md)."""
import numpy as np
import sympy as sp

from lagh.base import Candidate
from lagh.certify import arbitrate_significance, epsilon


def _c(s, cx):
    return Candidate(expr=sp.sympify(s), complexity=cx, channel="linear")


def _classes(*cands):
    return [(c.expr, [c]) for c in cands]


Y = np.random.default_rng(0).uniform(1, 10, 200)
EPS = epsilon(Y, sigma=0.0)


def test_muntz_contest_resolves_to_low_dof():
    truth = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    whale = _c("0.54*x_0**(17/5) - 0.095*x_0**(4/3) + 0.086*x_0**(2/3)"
               " + 1.2*x_0**(1/2) - 0.3*x_0**(5/2) + 0.7*x_0 - 0.2", 113)
    r = arbitrate_significance(_classes(truth, whale), Y, EPS, 500)
    assert r is not None
    assert r[0][1][0].expr == truth.expr
    assert "domain claim" in r[1]


def test_marginal_contest_stays_abstained():
    a = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    b = _c("3*x_0**2 - 1", 4)
    assert arbitrate_significance(_classes(a, b), Y, EPS, 500) is None


def test_single_class_is_noop():
    a = _c("x_0 + 1", 2)
    assert arbitrate_significance(_classes(a), Y, EPS, 500) is None


def test_two_whales_no_winner():
    w1 = _c("0.54*x_0**(17/5) - 0.095*x_0**(4/3) + 0.086*x_0**(2/3)"
            " + 1.2*x_0**(1/2) - 0.3*x_0**(5/2) + 0.7*x_0 - 0.2", 113)
    w2 = _c("0.4*x_0**(7/3) - 0.1*x_0**(5/4) + 0.09*x_0**(3/5)"
            " + 1.1*x_0**(1/3) - 0.2*x_0**(9/4) + 0.6*x_0 - 0.1", 110)
    assert arbitrate_significance(_classes(w1, w2), Y, EPS, 500) is None
