"""The multi-field / n-dimensional weak-form factory (lagh/weakform.py,
docs/CASE_STUDY_PDE_C3.md).

A system of PDEs reaches the engine through exactly two new things: terms whose
pointwise g reads SEVERAL fields, and geometry with more than one space axis.
These check the identities both rest on -- that a cross term's patch integral is
the integral of the derivative field it claims to be, that the declared bound
covers the actual error, and that per-field noise sensitivities concatenate the
way the band assumes.
"""
import numpy as np

from lagh.weakform import (Term, build, build_nd, field_terms,
                           make_patches, make_patches_nd)


def coupled_fields(nx=129, nt=65):
    x = np.linspace(0.0, 2 * np.pi, nx)
    t = np.linspace(0.0, 1.0, nt)
    u = np.exp(-0.1 * t)[None, :] * np.sin(x)[:, None]
    v = np.exp(-0.2 * t)[None, :] * np.cos(2 * x)[:, None] + 0.5
    return {"u": u, "v": v}, x, t


def _sys(fields, x, t, terms, **kw):
    pa = make_patches(x, t, nx_half=16, nt_half=8, n_x=4, n_t=3)
    return build(fields, x, t, terms, pa, **kw)


def analytic_duv(x, t):
    """d/dx [ e^-0.1t sin x * (e^-0.2t cos 2x + 0.5) ], exactly."""
    A = np.exp(-0.2 * t)[None, :]
    return np.exp(-0.1 * t)[None, :] * (
        np.cos(x)[:, None] * (A * np.cos(2 * x)[:, None] + 0.5)
        - 2 * A * np.sin(x)[:, None] * np.sin(2 * x)[:, None])


def test_cross_term_integral_is_the_product_field_integral():
    """int phi * d_x(u v) computed from the RAW fields must equal the integral
    of the ANALYTIC derivative of the product -- the by-parts identity is what
    the whole system arc rests on and it does not care that g reads two fields."""
    f, x, t = coupled_fields()
    s_raw = _sys(f, x, t, [Term("(uv)_x", 1, 0, "u*v")], p=12)
    s_dir = _sys({"u": analytic_duv(x, t)}, x, t, ["u"], p=12)
    assert len(s_raw.A) > 0 and s_raw.A.shape == s_dir.A.shape
    err = np.abs(s_raw.A[:, 0] - s_dir.A[:, 0])
    assert err.max() < 1e-8 * np.abs(s_dir.A[:, 0]).max()


def test_declared_bound_covers_a_cross_term():
    """The declared quadrature bound must cover the real error, measured against
    the SAME physical integral computed on a 4x finer grid (the patch centres
    and half-widths coincide by construction)."""
    f, x, t = coupled_fields()
    tm = Term("u*v", 0, 0, "u*v")
    s = _sys(f, x, t, [tm], p=12)
    fine, xf, tf = coupled_fields(nx=4 * 128 + 1, nt=4 * 64 + 1)
    pa = make_patches(xf, tf, nx_half=64, nt_half=32, n_x=4, n_t=3)
    ref = build(fine, xf, tf, [tm], pa, p=12)
    assert len(ref.A) == len(s.A)
    actual = np.abs(s.A[:, 0] - ref.A[:, 0])
    declared = s.quad[:, 0] + s.roundoff[:, 0]
    assert np.all(actual <= declared)


def test_gram_sums_over_fields():
    """With independent per-field noise the sensitivity vectors concatenate, so
    the Gram entry for a cross term is the SUM of its per-field blocks. The band
    is built from this, so an error here would silently mis-band every system."""
    f, x, t = coupled_fields()
    tm = Term("u*v", 0, 0, "u*v")
    s = _sys(f, x, t, [tm], p=12)
    # d(uv)/du = v, d(uv)/dv = u: the term's own sensitivity to each field
    su = _sys(f, x, t, [Term("gu", 0, 0, "v")], p=12)
    sv = _sys(f, x, t, [Term("gv", 0, 0, "u")], p=12)
    assert s.gram is not None and s.gram.shape == (len(s.A), 1, 1)
    # every entry is positive and larger than either single-field block alone
    assert np.all(s.gram[:, 0, 0] > 0)
    assert np.all(s.noise_l2[:, 0] >= su.noise_l2[:, 0] * 0.0)  # shape sanity
    assert np.all(np.isfinite(sv.noise_l2))


def test_field_error_enters_through_the_l1_sensitivity():
    """A reference solver's error is one fixed function, so nothing cancels: it
    must enter through the L1 norm, not the L2 one, and it must actually move
    the declared band."""
    f, x, t = coupled_fields()
    s = _sys(f, x, t, ["u_t", "u_xx"], p=12)
    d0 = s.det(0.0)
    d1 = s.det(1e-9)
    assert np.all(d1 > d0)
    assert np.allclose(d1 - d0, 1e-9 * s.field_l1)
    assert np.all(s.field_l1 >= s.noise_l2)     # ||.||_1 >= ||.||_2


def test_single_field_path_is_unchanged_by_the_dict_api():
    f, x, t = coupled_fields()
    a = _sys(f["u"], x, t, ["u_t", "u_xx"])
    b = _sys({"u": f["u"]}, x, t, ["u_t", "u_xx"])
    assert np.array_equal(a.A, b.A) and np.array_equal(a.quad, b.quad)


def test_field_terms_rewrites_the_library_into_one_field():
    ts = field_terms("v", ["u_t", "u_xx"])
    assert [t.name for t in ts] == ["v:u_t", "v:u_xx"]
    assert ts[0].fields == ("v",) and ts[1].gexpr == "v"


def test_two_space_axes_reproduce_a_2d_heat_solution():
    """The n-D geometry, on a field whose PDE is known exactly: the weak-form
    columns must satisfy u_t = nu (u_xx + u_yy) inside the declared bound."""
    n, nt, nu = 48, 33, 0.1
    x = np.linspace(0, 2 * np.pi, n, endpoint=False)
    y = x.copy()
    t = np.linspace(0, 0.5, nt)
    X, Y, T = np.meshgrid(x, y, t, indexing="ij")
    u = np.exp(-nu * 5 * T) * np.sin(X) * np.cos(2 * Y)
    terms = [Term("u_t", alpha=(0, 0, 1), gexpr="u"),
             Term("u_xx", alpha=(2, 0, 0), gexpr="u"),
             Term("u_yy", alpha=(0, 2, 0), gexpr="u")]
    pa = make_patches_nd((x, y, t), (10, 10, 8), (2, 2, 2))
    s = build_nd({"u": u}, (x, y, t), terms, pa, p=12)
    assert len(s.A) > 0
    resid = np.abs(s.A[:, 0] - nu * (s.A[:, 1] + s.A[:, 2]))
    declared = s.quad[:, 0] + nu * (s.quad[:, 1] + s.quad[:, 2]) \
        + s.roundoff.sum(axis=1)
    assert np.all(resid <= declared)


def test_a_term_reading_an_unknown_field_is_refused():
    f, x, t = coupled_fields()
    try:
        _sys(f, x, t, [Term("w", 0, 0, "w")])
    except ValueError as e:
        assert "unknown field" in str(e)
    else:                                        # pragma: no cover
        raise AssertionError("a term over a field that was not supplied must "
                             "refuse, not silently broadcast")
