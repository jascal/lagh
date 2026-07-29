"""The weak-form feature factory (lagh/weakform.py, docs/CASE_STUDY_PDE_DEV.md).

These tests check the thing the whole PDE arc rests on: that a patch integral of
the RAW field against analytic test-function derivatives reproduces the
integral of the derivative field, and that the DECLARED bound covers the actual
error. If the by-parts identity or the bound were wrong, every downstream
certificate would be a claim about numbers whose error was never declared.
"""
import numpy as np

from lagh.weakform import LIBRARY, build, make_patches


def heat_field(nu=0.1, k=1.0, nx=129, nt=65):
    x = np.linspace(0.0, 2 * np.pi, nx)
    t = np.linspace(0.0, 1.0, nt)
    u = np.exp(-nu * k ** 2 * t)[None, :] * np.sin(k * x)[:, None]
    return u, x, t


def advect_field(c=0.7, nx=129, nt=65):
    x = np.linspace(0.0, 2 * np.pi, nx)
    t = np.linspace(0.0, 1.0, nt)
    u = np.sin(x[:, None] - c * t[None, :]) + 0.3 * np.cos(
        2 * (x[:, None] - c * t[None, :]))
    return u, x, t


def _sys(u, x, t, terms, **kw):
    pa = make_patches(x, t, nx_half=16, nt_half=8, n_x=4, n_t=3)
    return build(u, x, t, terms, pa, **kw)


def test_by_parts_matches_the_derivative_field():
    """int phi*u_t computed from the RAW field must equal the integral of the
    analytic u_t -- the identity the whole approach is built on."""
    nu, k = 0.1, 1.0
    u, x, t = heat_field(nu, k)
    ut = -nu * k ** 2 * u                       # analytic u_t
    s_raw = _sys(u, x, t, ["u_t"])
    # the same patches, but integrating the analytic derivative directly:
    # int phi * u_t is the "u" term applied to the field u_t
    s_dir = _sys(ut, x, t, ["u"])
    assert s_raw.A.shape == s_dir.A.shape and len(s_raw.A) > 0
    err = np.abs(s_raw.A[:, 0] - s_dir.A[:, 0])
    scale = np.abs(s_dir.A[:, 0]).max()
    assert err.max() < 1e-5 * scale             # ~1e-9 absolute on ~1e-3 integrals


def test_declared_bound_covers_the_actual_error():
    """The declared quadrature bound must cover the measured error against the
    exact integral -- an under-declared bound is the failure mode that would
    make every weak-form certificate hollow."""
    nu, k = 0.1, 1.0
    u, x, t = heat_field(nu, k)
    ut = -nu * k ** 2 * u
    s_raw = _sys(u, x, t, ["u_t"])
    s_dir = _sys(ut, x, t, ["u"])
    actual = np.abs(s_raw.A[:, 0] - s_dir.A[:, 0])
    declared = s_raw.quad[:, 0] + s_raw.roundoff[:, 0] \
        + s_dir.quad[:, 0] + s_dir.roundoff[:, 0]
    assert np.all(actual <= declared)            # the bound COVERS the error
    assert np.all(declared < 1e3 * actual)       # ...and is not vacuously loose


def test_heat_law_holds_per_patch_within_the_declared_band():
    nu, k = 0.1, 1.0
    u, x, t = heat_field(nu, k)
    s = _sys(u, x, t, ["u_t", "u_xx", "u_x", "u"])
    ut = s.A[:, s.names.index("u_t")]
    uxx = s.A[:, s.names.index("u_xx")]
    res = np.abs(ut - nu * uxx)
    eps = s.declared_epsilon("u_t", coeff_max=1.0)
    assert np.all(res <= eps)
    assert np.abs(ut).max() > 1e3 * eps.max()   # not a vacuous patch family


def test_advection_law_holds_per_patch():
    c = 0.7
    u, x, t = advect_field(c)
    s = _sys(u, x, t, ["u_t", "u_x", "u_xx", "u"])
    ut = s.A[:, s.names.index("u_t")]
    ux = s.A[:, s.names.index("u_x")]
    res = np.abs(ut + c * ux)
    assert np.all(res <= s.declared_epsilon("u_t", coeff_max=1.0))


def test_burgers_law_holds_per_patch():
    """The exact traveling wave u = c - a*tanh(a(x-ct)/(2 nu)) solves
    u_t = nu u_xx - u u_x."""
    nu, a, c = 0.2, 1.0, 0.5
    x = np.linspace(-6.0, 6.0, 257)
    t = np.linspace(0.0, 2.0, 65)
    z = x[:, None] - c * t[None, :]
    u = c - a * np.tanh(a * z / (2 * nu))
    pa = make_patches(x, t, nx_half=24, nt_half=8, n_x=4, n_t=3)
    s = build(u, x, t, ["u_t", "u_xx", "u*u_x", "u_x", "u"], pa)
    ut = s.A[:, s.names.index("u_t")]
    uxx = s.A[:, s.names.index("u_xx")]
    uux = s.A[:, s.names.index("u*u_x")]
    res = np.abs(ut - (nu * uxx - uux))
    assert np.all(res <= s.declared_epsilon("u_t", coeff_max=1.0))


def test_a_wrong_coefficient_breaks_the_band():
    """The band must DISCRIMINATE: the true diffusivity fits, a perturbed one
    must not. Measured sensitivity for THIS patch family (p = 8, 33x17 points,
    no third-derivative term): nu(1 + 1e-4) is rejected, nu(1 + 1e-5) is not --
    so this family certifies ~4 digits of the coefficient. The campaign's p = 16
    family is far tighter (signal-to-band ~1e7)."""
    nu, k = 0.1, 1.0
    u, x, t = heat_field(nu, k)
    s = _sys(u, x, t, ["u_t", "u_xx", "u_x", "u"])
    ut = s.A[:, s.names.index("u_t")]
    uxx = s.A[:, s.names.index("u_xx")]
    eps = s.declared_epsilon("u_t", coeff_max=1.0)
    assert np.all(np.abs(ut - nu * uxx) <= eps)
    assert np.all(np.abs(ut - nu * 1.0001 * uxx) > eps)


def test_unresolved_patches_are_rejected_not_kept():
    """A field the grid cannot resolve must lose its patches, not acquire a
    quietly-wrong bound."""
    x = np.linspace(0.0, 2 * np.pi, 129)
    t = np.linspace(0.0, 1.0, 65)
    u = np.sin(31 * x[:, None]) * np.exp(-t[None, :])   # ~4 points per wavelength
    pa = make_patches(x, t, nx_half=16, nt_half=8, n_x=4, n_t=3)
    s = build(u, x, t, ["u_t", "u_xx"], pa)
    assert s.rejected > 0
    assert len(s.A) < len(pa)


def test_multiscale_family_breaks_the_constant_column():
    """A single-scale family makes the `1` column exactly constant, which the
    constrained-input detector reads (correctly) as a machine-exact input
    constraint. Pooling scales is what removes the degeneracy."""
    from lagh.weakform import multiscale_patches
    nu, k = 0.1, 1.0
    u, x, t = heat_field(nu, k, nx=257, nt=81)
    one = ["u_t", "u_xx", "1"]
    s1 = build(u, x, t, one, make_patches(x, t, nx_half=24, nt_half=12,
                                          n_x=4, n_t=3), p=16)
    sm = build(u, x, t, one, multiscale_patches(
        x, t, [(24, 12), (32, 16), (40, 20)], n_x=4, n_t=3), p=16)
    c1 = s1.A[:, s1.names.index("1")]
    cm = sm.A[:, sm.names.index("1")]
    assert (c1.max() - c1.min()) / abs(c1).max() < 1e-12      # degenerate
    assert (cm.max() - cm.min()) / abs(cm).max() > 0.1        # a real direction


def test_normalize_preserves_the_law_and_scales_the_bands():
    """Row normalization must not move the law: dividing a row by its own ∫φ
    cancels out of a linear relation, and the bands must ride along."""
    from lagh.weakform import multiscale_patches
    nu = 0.1
    u, x, t = heat_field(nu, 1.0, nx=257, nt=81)
    terms = ["u_t", "u_xx", "u_x", "1"]
    s = build(u, x, t, terms, multiscale_patches(
        x, t, [(24, 12), (32, 16)], n_x=4, n_t=3), p=16)
    n = s.normalize(by="1")
    assert "1" not in n.names and len(n.names) == len(s.names) - 1
    w = s.A[:, s.names.index("1")]
    # the law survives: u_t = nu*u_xx holds in both, within each one's own band
    for sys_ in (s, n):
        ut = sys_.A[:, sys_.names.index("u_t")]
        uxx = sys_.A[:, sys_.names.index("u_xx")]
        assert np.all(np.abs(ut - nu * uxx)
                      <= sys_.declared_epsilon("u_t", coeff_max=1.0))
    # and the rows really were divided by their own test-function integral
    assert np.allclose(n.A[:, n.names.index("u_t")],
                       s.A[:, s.names.index("u_t")] / w)
    assert np.allclose(n.gram, s.gram[:, :-1, :-1] / (w ** 2)[:, None, None])


def test_library_signs_are_the_by_parts_signs():
    for name, term in LIBRARY.items():
        assert term.sign == (-1) ** (term.ax + term.at), name
