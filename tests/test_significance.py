"""Significance accounting (DIRECTION_SIGNIFICANCE.md): alpha = |H| * q^h."""

import numpy as np
import sympy as sp

from lagh.certify import significance_log10
from lagh.engine import discover


def _splits(fn, n=200, seed=3, dim=1):
    rng = np.random.default_rng(seed)
    X = np.exp(rng.uniform(np.log(0.5), np.log(10), (n, dim)))
    y = fn(X)
    a, b = int(0.6 * n), int(0.8 * n)
    return X[:a], y[:a], X[a:b], y[a:b], X[b:], y[b:]


def test_alpha_hand_calculation():
    """40 held-out points, uniform eps giving q = 1e-3, |H| = 1000, dof = 2:
    alpha = 1e3 * (1e-3)^38 -> log10 = 3 - 114 = -111."""
    x0 = sp.Symbol("x_0")
    expr = sp.Rational(3, 2) * x0 ** 2          # dof = 2 numeric atoms (3/2, 2)
    y = np.linspace(0.0, 1.0, 40)
    eps = np.full(40, 5e-4)                      # q = 2*eps/range = 1e-3
    lg = significance_log10(expr, y, eps, 1000)
    assert abs(lg - (3 - 38 * 3)) < 1e-6


def test_certificate_carries_alpha():
    r = discover(*_splits(lambda X: 3 * X[:, 0] ** 2))
    assert r.certificate.certified
    assert r.certificate.n_hypotheses and r.certificate.n_hypotheses > 0
    assert r.certificate.alpha_log10 is not None
    # machine-precision certification over ~40 held-out points is hyper-significant
    assert r.certificate.alpha_log10 < -50


def test_alpha_degrades_with_wide_eps():
    """Same expr, looser eps -> larger (weaker) alpha."""
    x0 = sp.Symbol("x_0")
    expr = x0 ** 2
    y = np.linspace(0.0, 1.0, 40)
    tight = significance_log10(expr, y, np.full(40, 1e-6), 100)
    loose = significance_log10(expr, y, np.full(40, 1e-2), 100)
    assert tight < loose
