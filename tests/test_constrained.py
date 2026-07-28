"""Constrained-input coherence (closure of the GAIA_P3 registered issue)."""
import numpy as np
import sympy as sp

from lagh.certify import input_constraints
from lagh.passive import discover_passive


def _circle_data(n=200, seed=0):
    th = np.random.default_rng(seed).uniform(0, 2 * np.pi, n)
    return np.column_stack([np.cos(th), np.sin(th)])


def test_constraint_detected_on_circle():
    X = _circle_data()
    syms = list(sp.symbols("x_0 x_1"))
    cons = input_constraints(X, syms)
    assert len(cons) == 1
    g = sp.expand(cons[0])
    x0, x1 = syms
    # proportional to x0^2 + x1^2 - 1
    ratio = sp.simplify(g / (x0**2 + x1**2 - 1))
    assert ratio.is_number


def test_no_false_positive_full_dim():
    X = np.random.default_rng(1).uniform(-2, 2, (200, 2))
    syms = list(sp.symbols("x_0 x_1"))
    assert input_constraints(X, syms) == []


def test_no_false_positive_statistical_correlation():
    rng = np.random.default_rng(2)
    x0 = rng.uniform(0, 1, 300)
    x1 = 0.6 * x0 + rng.normal(0, 0.2, 300)      # rho ~ 0.6, not exact
    X = np.column_stack([x0, x1])
    syms = list(sp.symbols("x_0 x_1"))
    assert input_constraints(X, syms) == []


def test_certifies_linear_law_on_circle():
    """On the unit circle the truth 2*x0 - 3*x1 has box-divergent rivals
    (any multiple of the constraint added); the closure must certify the
    canonical form with the domain-restriction note."""
    X = _circle_data(n=240, seed=3)
    y = 2.0 * X[:, 0] - 3.0 * X[:, 1]
    r = discover_passive(X, y, sigma=0.0)
    assert r.certified, r.result.certificate.abstain
    x0, x1 = sp.symbols("x_0 x_1")
    assert sp.simplify(r.result.expr - (2 * x0 - 3 * x1)) == 0
