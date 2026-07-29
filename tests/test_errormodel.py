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
