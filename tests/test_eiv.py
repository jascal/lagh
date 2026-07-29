"""Errors-in-variables epsilon: per-candidate bands, interval parameters.

The standard model bands y and treats X as exact. Weak-form rows (weakform.py)
break that: every column is a noisy functional of the SAME field, so the band
depends on the candidate's own coefficients and the column errors partly cancel.
These tests pin the three pieces that follow from it -- a callable epsilon, a
band that is tighter than the coefficient-bounded one AND still covers, and a
parameter INTERVAL where an exact rational is the wrong claim.
"""
import numpy as np
import sympy as sp

from lagh.certify import (band, check, epsilon, gated_atoms,
                          parameter_interval)
from lagh.weakform import PatchEpsilon, build, make_patches

X1 = sp.Symbol("x_0")


def heat(nu=0.1, nx=257, nt=81, sigma=0.0, seed=0):
    x = np.linspace(0.0, 2 * np.pi, nx)
    t = np.linspace(0.0, 1.0, nt)
    u = np.exp(-nu * t)[None, :] * np.sin(x)[:, None]
    if sigma:
        u = u + np.random.default_rng(seed).normal(0, sigma, u.shape)
    return u, x, t


TERMS = ["u_t", "u_xx", "u*u_x", "u_x", "u", "1"]
FEAT = [n for n in TERMS if n != "u_t"]


def patch_model(sigma, seed=0, coeff_max=2.0):
    u, x, t = heat(sigma=sigma, seed=seed)
    s = build(u, x, t, TERMS,
              make_patches(x, t, nx_half=24, nt_half=12, n_x=4, n_t=3),
              p=16, sigma=sigma)
    j = TERMS.index("u_t")
    cols = [TERMS.index(n) for n in FEAT]
    return PatchEpsilon(TERMS, "u_t", s.A[:, j], s.A[:, cols],
                        s.quad + s.roundoff, s.gram, sigma=sigma,
                        floor_abs=0.0, coeff_max=coeff_max), s


def test_band_resolves_arrays_and_callables():
    y = np.array([1.0, 2.0, 3.0])
    fixed = epsilon(y, floor_abs=1e-9)
    assert np.allclose(band(fixed), fixed)
    assert np.allclose(band(lambda expr: fixed, X1), fixed)
    # a callable band is a function OF THE CANDIDATE
    per = band(lambda expr: np.full(3, float(expr.subs(X1, 1))), 2 * X1)
    assert np.allclose(per, 2.0)


def test_check_uses_the_candidates_own_band():
    """The whole point: two candidates get two different bands."""
    y = np.array([1.0, 1.0])
    X = np.array([[1.0], [1.0]])

    def model(expr):
        # a band that is generous for the 1.5 law and tight for every other
        c = float(sp.expand(expr).coeff(X1)) if expr is not None else 0.0
        return np.full(2, 1.0 if abs(c - 1.5) < 1e-12 else 1e-9)

    assert check(sp.Float(1.5) * X1, [X1], X, y, model)["certified"]
    assert not check(sp.Float(1.4) * X1, [X1], X, y, model)["certified"]


def test_eiv_band_is_tighter_than_the_coefficient_bounded_one():
    m, _ = patch_model(1e-5)
    truth = sp.Rational(1, 10) * X1
    assert np.all(m(truth) < m(None))          # per-candidate beats bounding
    assert np.median(m(truth) / m(None)) < 0.5


def test_eiv_band_covers_the_noise_it_declares():
    """Over independent noise draws the true law's residual must stay inside the
    declared band -- an under-declared band would reject true laws, an
    over-declared one would admit impostors."""
    for seed in range(8):
        m, s = patch_model(1e-4, seed=seed)
        y = s.A[:, TERMS.index("u_t")]
        uxx = s.A[:, TERMS.index("u_xx")]
        assert np.all(np.abs(y - 0.1 * uxx) <= m(sp.Rational(1, 10) * X1))


def test_declared_noise_does_not_cost_patches():
    """The resolution gates must reject for UNDER-RESOLUTION, not for the noise
    the eps model already bands (measured: 16 of 24 patches lost before the
    ladder and aliasing tests were noise-corrected)."""
    _, clean = patch_model(0.0)
    _, noisy = patch_model(1e-4)
    assert noisy.rejected == 0
    assert len(noisy.A) == len(clean.A)


def test_parameter_interval_brackets_the_truth():
    y = np.array([1.0, 2.0, 3.0, 4.0])
    X = np.array([[1.0], [2.0], [3.0], [4.0]])
    eps = np.full(4, 0.05)
    iv = parameter_interval(sp.Float(1.0) * X1, [X1], X, y, eps,
                            sp.Float(1.0))
    assert iv is not None
    lo, hi = iv
    assert lo < 1.0 < hi
    # the widest coefficient that keeps every point within 0.05 is 1 + 0.05/4
    assert abs(hi - 1.0125) < 1e-3 and abs(lo - 0.9875) < 1e-3


def test_parameter_interval_reports_undetermined():
    """A parameter the data does not constrain returns None -- which is what
    keeps the interval certificate from claiming an unbounded family."""
    y = np.array([1.0, 1.0, 1.0])
    X = np.array([[0.0], [0.0], [0.0]])         # the column carries no signal
    eps = np.full(3, 0.5)
    assert parameter_interval(sp.Float(0.7) * X1 + sp.Float(1.0), [X1], X, y,
                              eps, sp.Float(0.7)) is None


def test_gated_atoms_catches_floats_in_costume():
    e = sp.Rational(90360785, 903607052) * X1 + sp.Rational(1, 10)
    got = {str(a) for a in gated_atoms(e)}
    assert "90360785/903607052" in got and "1/10" not in got
