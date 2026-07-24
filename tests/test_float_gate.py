"""The exact-coefficient gate (CAP-E lesson): Float coefficients must be pinned by the
data at certification precision, and snap to exact rationals when one certifies."""

import numpy as np
import sympy as sp

from lagh.certify import check, epsilon, float_pinned
from lagh.engine import discover

X0 = sp.Symbol("x_0")


def _splits(fn, n=200, seed=3, lo=0.5, hi=10.0):
    rng = np.random.default_rng(seed)
    X = np.exp(rng.uniform(np.log(lo), np.log(hi), (n, 1)))
    y = fn(X)
    a, b = int(0.6 * n), int(0.8 * n)
    return X[:a], y[:a], X[a:b], y[a:b], X[b:], y[b:]


def test_unpinned_float_rejected():
    """Loose epsilon (the reverted-BE regime): a coefficient 3e-5 off still certifies,
    so perturbed neighbours certify too -> the gate must reject."""
    rng = np.random.default_rng(0)
    X = np.exp(rng.uniform(np.log(0.5), np.log(10), (60, 1)))
    y = X[:, 0] ** 2
    eps = 1e-3 * np.abs(y)                       # loose: 1e-3 relative
    expr = sp.Float(1.00003) * X0**2             # "certifies" but is not pinned
    assert check(expr, [X0], X, y, eps)["certified"]
    ok, _ = float_pinned(expr, [X0], X, y, eps)
    assert not ok


def test_pinned_float_accepted_and_snapped():
    """Tight epsilon: the coefficient is identified; the gate accepts and snaps it to
    the exact decimal rational."""
    rng = np.random.default_rng(1)
    X = np.exp(rng.uniform(np.log(0.5), np.log(10), (60, 1)))
    y = 6.674e-5 * X[:, 0] ** 2
    eps = epsilon(y)
    expr = sp.Float(6.674e-5) * X0**2
    ok, snapped = float_pinned(expr, [X0], X, y, eps)
    assert ok
    assert not snapped.atoms(sp.Float)           # now an exact rational
    assert check(snapped, [X0], X, y, eps)["certified"]


def test_no_floats_is_a_noop():
    expr = sp.Rational(3, 2) * X0**2
    ok, out = float_pinned(expr, [X0], np.ones((10, 1)), 1.5 * np.ones(10),
                           np.full(10, 1e-9))
    assert ok and out == expr


def test_discover_emits_exact_rationals_for_c3_winners():
    """End-to-end: the C3 power-law channel used to emit raw Float coefficients; the
    winner must now carry exact rationals only."""
    r = discover(*_splits(lambda X: 6.674e-5 * X[:, 0] ** 2 / 1.0))
    assert r.certificate.certified
    assert not r.expr.atoms(sp.Float)


def test_discover_still_recovers_previous_c3_form():
    """Regression: the gate must not break a previously-recovered power law."""
    r = discover(*_splits(lambda X: 2.5 * X[:, 0] ** sp.Rational(3, 2).__float__()))
    assert r.certificate.certified
