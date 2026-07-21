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


# -- P2: active acquisition ---------------------------------------------------

def test_active_recovers_and_ledgers():
    from lagh.acquisition import run_active
    f = lambda X: 2 * X[:, 0] * X[:, 1] / X[:, 2] ** 2      # noqa: E731
    r = run_active(f, [0.5] * 3, [3.0] * 3, seed=1)
    assert r.result.certificate.certified
    assert r.queries_used == r.ledger.spent > 0
    kinds = {e["kind"] for e in r.ledger.entries}
    assert "ranging" in kinds and "init" in kinds


def test_adaptive_ranging_contracts_dead_box():
    """Decay law over a huge box: most outputs underflow the floor (the N1'
    killer). Ranging must contract toward the live region and still recover."""
    from lagh.acquisition import run_active
    f = lambda X: 3 * X[:, 0] * np.exp(-X[:, 1] * X[:, 2])  # noqa: E731
    r = run_active(f, [0.1] * 3, [100.0] * 3, seed=2)
    contracted = np.array(r.box_final)
    assert (contracted[1] < 100.0 - 1e-9).any(), "box never contracted"
    assert r.result.certificate.certified or \
        r.result.certificate.abstain in ("structural", "range")


def test_range_abstention_when_no_signal_anywhere():
    from lagh.acquisition import run_active
    from lagh import Abstain
    f = lambda X: 1e-15 * X[:, 0]                            # noqa: E731
    r = run_active(f, [0.1] * 2, [100.0] * 2, seed=3)
    assert r.result.certificate.abstain in (Abstain.RANGE.value,
                                            Abstain.NOISE.value)


def test_micro_predictions_recorded_and_scored():
    from lagh.acquisition import run_active, Policy
    # a law needing escalation, so at least one active round fires
    f = lambda X: 5 * X[:, 0] * np.cos(X[:, 1]) ** 2         # noqa: E731
    r = run_active(f, [0.1] * 2, [1.5] * 2, seed=4,
                   policy=Policy(init_points=14))
    for rec in r.predictions:
        for p in rec["predictions"]:
            assert 0.0 <= p["hit_rate"] <= 1.0


# -- Ehrhart quasi-polynomial recovery (exact integer domain) ----------------

def test_ehrhart_counts_are_exact():
    from fractions import Fraction
    from lagh.adapters.ehrhart import lattice_count
    # unit interval [0,1]: L(t) = t+1
    assert [lattice_count([Fraction(1)], t) for t in (1, 2, 5)] == [2, 3, 6]
    # [0, 3/2]: L(t) = floor(3t/2)+1 -> t=1:2, t=2:4, t=3:5
    assert [lattice_count([Fraction(3, 2)], t) for t in (1, 2, 3)] == [2, 4, 5]


def test_ehrhart_recovers_and_validates_out_of_range():
    from fractions import Fraction
    from lagh.adapters.ehrhart import lattice_count
    from lagh.quasipoly import recover
    inter = [Fraction(5, 2), Fraction(7, 3)]           # period 6, degree 2
    ts = list(range(1, 49))
    Ls = [lattice_count(inter, t) for t in ts]
    r = recover(ts, Ls)
    assert r.certified
    # exact match far beyond the fit range -- interpolation cannot do this
    for t in range(60, 90):
        assert r.quasipoly(t) == lattice_count(inter, t)


def test_ehrhart_abstains_when_budget_too_small():
    from fractions import Fraction
    from lagh.adapters.ehrhart import lattice_count
    from lagh.quasipoly import recover
    inter = [Fraction(5, 2), Fraction(7, 3)]           # needs many t for period 6
    ts = list(range(1, 8))                             # far too few
    Ls = [lattice_count(inter, t) for t in ts]
    r = recover(ts, Ls)
    assert not r.certified and r.abstain == "range"


def test_ehrhart_never_confidently_wrong_on_random_polytopes():
    from lagh.adapters.ehrhart import lattice_count, random_simplex
    from lagh.quasipoly import recover
    rng = np.random.default_rng(99)
    for _ in range(12):
        dim = int(rng.integers(1, 4))
        inter = random_simplex(dim, rng)
        ts = list(range(1, 49))
        Ls = [lattice_count(inter, t) for t in ts]
        r = recover(ts, Ls)
        if r.certified:                                # if it answers, it is right
            for t in range(50, 80):
                assert r.quasipoly(t) == lattice_count(inter, t)


# -- C6 promotion: discover() escalates to the quasi-poly tier ----------------

def test_discover_escalates_to_c6_on_integer_lattice():
    from fractions import Fraction
    from lagh.adapters.ehrhart import lattice_count
    from lagh import discover
    inter = [Fraction(5, 2), Fraction(7, 3)]           # period 6, degree 2
    ts = np.arange(1, 49, dtype=float)
    Ls = np.array([lattice_count(inter, int(t)) for t in ts], float)
    i = np.random.default_rng(0).permutation(48)
    a, b = 30, 39
    r = discover(ts[i[:a], None], Ls[i[:a]], ts[i[a:b], None], Ls[i[a:b]],
                 ts[i[b:], None], Ls[i[b:]])
    assert r.certificate.certified and r.tier == 6
    for t in range(60, 80):
        assert r.expr(t) == lattice_count(inter, t)


def test_discover_still_prefers_float_tiers_for_smooth_laws():
    """A plain polynomial over integer x must NOT be over-promoted to C6."""
    from lagh import discover
    f = lambda X: 2 * X[:, 0] ** 2 - 3 * X[:, 0] + 1     # noqa: E731
    X = np.arange(1, 40, dtype=float)[:, None]
    y = f(X)
    i = np.random.default_rng(1).permutation(39)
    r = discover(X[i[:24]], y[i[:24]], X[i[24:31]], y[i[24:31]],
                 X[i[31:]], y[i[31:]])
    assert r.certificate.certified and r.tier <= 2      # C1/C2, not C6
