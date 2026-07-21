"""P1 acceptance: every curriculum tier recovers its class, refusals refuse.

The recovery targets are the shapes the predecessor measured (pilot + N1' smoke);
the refusal tests are the load-bearing half -- the product invariant is zero wrong
answers, and a discoverer that always answers has no product.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import sympy as sp
import pytest

from lagh import Abstain, discover

RNG = np.random.default_rng(0)


def splits(f, dim, lo, hi, n=60, sigma=0.0, seed=0):
    rng = np.random.default_rng(seed)
    X = np.exp(rng.uniform(np.log(lo), np.log(hi), (n, dim)))
    y = f(X)
    if sigma:
        y = y * (1 + sigma * rng.standard_normal(len(y)))
    i = rng.permutation(n)
    a, b = int(0.6 * n), int(0.8 * n)
    return (X[i[:a]], y[i[:a]], X[i[a:b]], y[i[a:b]], X[i[b:]], y[i[b:]])


def equivalent(expr, f, dim, lo, hi, tol=1e-6):
    P = np.exp(RNG.uniform(np.log(lo), np.log(hi), (200, dim)))
    syms = sp.symbols([f"x_{i}" for i in range(dim)])
    got = sp.lambdify(syms, expr, "numpy")(*P.T)
    want = f(P)
    return np.allclose(np.broadcast_to(got, want.shape), want,
                       rtol=tol, atol=tol * np.abs(want).max())


# -- recovery, one per tier ---------------------------------------------------

def test_c1_polynomial():
    f = lambda X: 2.1 * X[:, 0] - 0.5 * X[:, 0] ** 2        # noqa: E731
    r = discover(*splits(f, 1, 0.5, 4.0))
    assert r.certificate.certified and equivalent(r.expr, f, 1, 0.5, 4.0)
    assert r.tier == 1


def test_c2_rational_inverse_square():
    f = lambda X: 2 * X[:, 0] * X[:, 1] / X[:, 2] ** 2      # noqa: E731
    r = discover(*splits(f, 3, 0.5, 3.0))
    assert r.certificate.certified and equivalent(r.expr, f, 3, 0.5, 3.0)


def test_c3_powerlaw_fractional():
    f = lambda X: 3.0 * X[:, 0] / X[:, 1] ** 1.5            # noqa: E731
    r = discover(*splits(f, 2, 0.5, 3.0))
    assert r.certificate.certified and equivalent(r.expr, f, 2, 0.5, 3.0)


def test_c4_decay():
    f = lambda X: 3 * X[:, 0] * X[:, 1] * np.exp(-X[:, 1] * X[:, 2])  # noqa: E731
    r = discover(*splits(f, 3, 0.2, 3.0))
    assert r.certificate.certified and equivalent(r.expr, f, 3, 0.2, 3.0)


def test_c4_cos_squared():
    f = lambda X: 5 * X[:, 0] * np.cos(X[:, 1]) ** 2        # noqa: E731
    r = discover(*splits(f, 2, 0.1, 1.5))
    assert r.certificate.certified and equivalent(r.expr, f, 2, 0.1, 1.5)


def test_c5_bose_einstein_shape():
    f = lambda X: X[:, 0] / (np.exp(X[:, 0] / X[:, 1]) - 1)  # noqa: E731
    r = discover(*splits(f, 2, 0.3, 5.0))
    assert r.certificate.certified and equivalent(r.expr, f, 2, 0.3, 5.0)


def test_c5_sqrt_rational():
    def f(X):
        return np.sqrt(X[:, 0] / X[:, 1] - X[:, 2] ** 2 / 4)
    rng = np.random.default_rng(3)
    X = np.column_stack([np.exp(rng.uniform(np.log(2), np.log(4), 60)),
                         np.exp(rng.uniform(np.log(0.5), 0.0, 60)),
                         np.exp(rng.uniform(np.log(0.5), 0.0, 60))])
    y = f(X)
    i = rng.permutation(60)
    r = discover(X[i[:36]], y[i[:36]], X[i[36:48]], y[i[36:48]],
                 X[i[48:]], y[i[48:]])
    assert r.certificate.certified
    syms = sp.symbols(["x_0", "x_1", "x_2"])
    got = sp.lambdify(syms, r.expr, "numpy")(*X.T)
    assert np.allclose(got, y, rtol=1e-6)


# -- refusal: the load-bearing half ------------------------------------------

def test_pure_noise_never_answers():
    rng = np.random.default_rng(7)
    X = np.exp(rng.uniform(-1, 1, (60, 2)))
    y = rng.standard_normal(60) * 5.0
    i = rng.permutation(60)
    r = discover(X[i[:36]], y[i[:36]], X[i[36:48]], y[i[36:48]],
                 X[i[48:]], y[i[48:]])
    assert r.abstained
    assert r.certificate.abstain in {a.value for a in Abstain}


def test_vacuity_when_eps_swallows_signal():
    rng = np.random.default_rng(9)
    X = np.exp(rng.uniform(-1, 1, (60, 1)))
    y = 1e-15 * X[:, 0]                     # signal below the 1e-12 floor
    i = rng.permutation(60)
    r = discover(X[i[:36]], y[i[:36]], X[i[36:48]], y[i[36:48]],
                 X[i[48:]], y[i[48:]])
    assert r.abstained and r.certificate.abstain == Abstain.NOISE.value


def test_noise_never_confidently_wrong():
    """Under real noise: certified-and-close, or abstained. Never wrong."""
    f = lambda X: 2.1 * X[:, 0] - 0.5 * X[:, 0] ** 2        # noqa: E731
    r = discover(*splits(f, 1, 0.5, 4.0, sigma=0.01, seed=11), sigma=0.01)
    if r.certificate.certified:
        assert equivalent(r.expr, f, 1, 0.5, 4.0, tol=0.1)
    else:
        assert r.certificate.abstain in {a.value for a in Abstain}


def test_hardened_coefficients_are_exact_rationals():
    f = lambda X: 0.25 + 0.5 * X[:, 0]                      # noqa: E731
    r = discover(*splits(f, 1, 0.5, 4.0))
    assert r.certificate.certified
    assert all(a.is_Rational for a in r.expr.atoms(sp.Number))
