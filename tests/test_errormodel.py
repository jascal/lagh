"""Error provenance (lagh/errormodel.py, docs/DIRECTION_ERROR_PROVENANCE.md).

The verdict decides which band channel is correct -- L2 for a stochastic error
that cancels across a realization, L1 for a deterministic one that does not --
so the failure that matters is calling a deterministic error stochastic, which
UNDER-declares and admits impostors. These check that the discriminators fire,
that `undetermined` is returned rather than a default, and that the dangerous
direction does not happen.
"""
import numpy as np

from lagh.errormodel import EXPLAINS, characterize, modified_equation


def test_a_modified_equation_signature_is_found_and_named():
    """The residual of a stated law regressed on the next derivatives IS the
    discretization's signature (measured on PDEBench advection: u_xxx explains
    84% at c3 = 1.09e-7, u_xx 1.5%)."""
    rng = np.random.default_rng(0)
    uxxx = rng.normal(0, 1, 200)
    uxx = rng.normal(0, 1, 200)
    resid = 1.09e-7 * uxxx + 1e-9 * rng.normal(0, 1, 200)
    me = modified_equation(resid, {"u_xx": uxx, "u_xxx": uxxx})
    assert me["best"] == "u_xxx"
    assert me["best_variance_explained"] > 0.9
    assert np.isclose(me["terms"]["u_xxx"]["coefficient"], 1.09e-7, rtol=0.05)
    assert me["explains"]


def test_structured_error_selects_the_L1_channel():
    rng = np.random.default_rng(1)
    uxxx = rng.normal(0, 1, 200)
    r = characterize(2e-7 * uxxx, columns={"u_xxx": uxxx})
    assert r.verdict == "structured-deterministic"
    assert "L1" in r.recommended_channel
    assert "pipeline" in r.likely_reading
    assert not r.certified                       # never a certificate


def test_white_residual_selects_the_L2_channel():
    rng = np.random.default_rng(2)
    cols = {k: rng.normal(0, 1, 300) for k in ("u_xx", "u_xxx")}
    r = characterize(rng.normal(0, 1e-4, 300), columns=cols)
    assert r.verdict == "unstructured-stochastic"
    assert "L2" in r.recommended_channel
    assert r.modified_equation["best_variance_explained"] < EXPLAINS


def test_growth_in_time_is_a_structure_test():
    """Storage rounding is flat in t; an accumulated scheme error is not."""
    t = np.repeat(np.linspace(0.0, 2.0, 40), 5)
    rng = np.random.default_rng(3)
    r = characterize(rng.normal(0, 1, t.size) * (1e-6 + 1e-3 * t), times=t)
    assert r.verdict == "structured-deterministic"
    assert r.evidence["growth_late_over_early"] > 3.0


def test_nothing_measured_returns_undetermined_not_a_default():
    """A single trajectory with no replicates and no reference cannot separate
    the two; saying so is the correct answer."""
    r = characterize(np.random.default_rng(4).normal(0, 1, 50))
    assert r.verdict == "undetermined"
    assert "UNDETERMINED" in " ".join(r.notes)
    # ...and while the kind is unknown the CONSERVATIVE channel is recommended
    assert "L1" in r.recommended_channel


def test_the_stochastic_reading_is_reported_as_inconclusive_without_replicates():
    rng = np.random.default_rng(5)
    cols = {"u_xxx": rng.normal(0, 1, 200)}
    r = characterize(rng.normal(0, 1e-5, 200), columns=cols)
    assert r.verdict == "unstructured-stochastic"
    assert any("not a demonstration" in n for n in r.notes)


def test_a_free_fit_residual_is_a_lower_bound_and_refuses_to_be_a_declaration():
    """The dangerous direction, measured on PDEBench advection: with no stated
    law the fit ABSORBS the error being measured, so the residual is smaller and
    less structured than the truth's. u_xxx's share fell from 84% to 13.6%, the
    verdict flipped to stochastic, and the derived magnitude came out FOUR
    ORDERS below the independently measured solver error. Auto-setting a band
    that way manufactures confident-wrongs, so the code must refuse."""
    from types import SimpleNamespace

    from lagh.errormodel import characterize_rows

    rng = np.random.default_rng(7)
    n = 120
    # the columns must be CORRELATED the way real ones are: for Fourier mode k,
    # u_xxx = -k^2 u_x, so a fit on u_x alone can absorb part of a u_xxx term.
    # Independent columns cannot reproduce the effect at all -- measured, in a
    # first version of this test.
    ux = np.zeros(n)
    uxxx = np.zeros(n)
    for k in range(1, 9):
        c = rng.normal(0, 1, n)
        a = rng.uniform(0.2, 1.0) / k
        ux += a * c
        uxxx += -(k ** 2) * a * c
    # a field whose law is u_t = -0.7 u_x plus a dispersive pipeline term
    ut = -0.7 * ux + 1e-3 * uxxx
    names = ["u_t", "u_x", "u_xxx"]
    A = np.column_stack([ut, ux, uxxx])
    rows = SimpleNamespace(
        names=names, A=A,
        field_l1=np.ones((n, 3)) * 1e-3,
        gram=np.tile(np.eye(3), (n, 1, 1)),
        features=lambda: ["u_x", "u_xxx"])

    stated, req_s = characterize_rows(rows, "u_t", extensions=("u_xxx",),
                                      stated_law={"u_x": -0.7})
    assert stated.verdict == "structured-deterministic"
    assert req_s["usable_as_declaration"]

    free, req_f = characterize_rows(rows, "u_t", extensions=("u_xxx",))
    assert free.verdict == "undetermined"          # refuses to decide
    assert not req_f["usable_as_declaration"]
    assert "NOT A DECLARATION" in free.notes[0]
    # ...and the numbers it did compute are strictly smaller: that is the bias
    assert req_f.get("field_err", 0) <= req_s.get("field_err", 0)
