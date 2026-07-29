"""The PDEBench adapter (experiments/pde/pdebench.py, docs/PDEBENCH_READINESS.md).

The adapter never fits anything; its whole job is to turn a downloaded array
into observations whose error is DECLARED. These check the two declarations that
would otherwise be silent -- float32 storage quantization, and a grid geometry
that rescales every derivative term without an error message.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from experiments.pde.pdebench import (FLOAT32_EPS, check_geometry,  # noqa: E402
                                      declared_noise, from_arrays,
                                      sigma_rep_for)


def field(nx=64, nt=17):
    x = np.linspace(0.0, 1.0, nx, endpoint=False)
    t = np.linspace(0.0, 0.5, nt)
    u = np.exp(-t)[None, :] * np.sin(2 * np.pi * x)[:, None]
    return u, x, t


def test_float32_storage_noise_is_computed_not_assumed():
    """A float32 array is the simulation's numbers ROUNDED; certifying against
    it without declaring that is the float32 lesson from the LLM-SRBench read."""
    u, x, t = field()
    u32 = u.astype(np.float32).astype(float)
    abs_, rel = sigma_rep_for(u)
    assert np.isclose(rel, 0.5 * FLOAT32_EPS)
    assert np.max(np.abs(u32 - u)) <= abs_        # the declaration COVERS it
    assert abs_ > 0


def test_declared_noise_names_what_is_measured_and_what_is_assumed():
    u, x, t = field()
    ds = from_arrays({"u": u}, x, t, field_err=1e-9)
    d = declared_noise(ds, extra_sigma=1e-5)
    assert d["sigma"] >= 1e-5 and d["sigma_rep_absolute"] > 0
    assert d["field_err_declared"] == 1e-9
    # the solver error CANNOT be measured after the fact for shipped data, and
    # the record has to say so rather than let a reader assume a ladder ran
    assert d["field_err_is_measured"] is False


def test_float32_coordinates_are_regularized_and_the_deviation_reported():
    """The coordinate VECTORS are float32 too, and their rounding is ~1e-6
    relative to the STEP -- which lands on the u_t column, because the by-parts
    weights carry 1/step^k with k per term."""
    u, x, t = field(nt=101)
    t32 = t.astype(np.float32).astype(float)
    raw_jitter = np.max(np.abs(np.diff(t32) - np.diff(t32)[0])) / np.diff(t32)[0]
    assert raw_jitter > 1e-7                      # the problem is real
    ds = from_arrays({"u": u}, x, t32)
    d = np.diff(ds.coords[-1])
    assert np.max(np.abs(d - d[0])) <= 1e-12 * d[0]          # ...and removed
    dev = ds.provenance["coord_deviation"][-1]
    assert 0 < dev["max_abs"] < 1e-6              # ...and reported, not hidden
    assert declared_noise(ds)["coord_regularized"] is True


def test_geometry_check_catches_a_nonuniform_axis():
    u, x, t = field()
    bad = x.copy()
    bad[10:] += 0.001                              # a subtly non-uniform grid
    ds = from_arrays({"u": u}, bad, t)
    rep = check_geometry(ds)
    assert not rep["ok"]
    assert any("uniform" in n for n in rep["notes"])


def test_geometry_check_flags_a_duplicated_periodic_endpoint():
    nx = 64
    x = np.linspace(0.0, 1.0, nx)                  # INCLUDES the right endpoint
    t = np.linspace(0.0, 0.5, 17)
    u = np.exp(-t)[None, :] * np.sin(2 * np.pi * x)[:, None]
    ds = from_arrays({"u": u}, x, t)
    rep = check_geometry(ds)
    assert any("endpoint" in n for n in rep["notes"])


def test_a_loaded_dataset_feeds_the_factory_unchanged():
    """The adapter's output must be exactly what build() already takes -- the
    point is a declared error model, not a new pipeline."""
    from lagh.weakform import build, make_patches
    u, x, t = field(nx=128, nt=65)
    ds = from_arrays({"u": u}, x, t, name="synthetic")
    pa = make_patches(ds.coords[0], ds.coords[-1], nx_half=16, nt_half=8,
                      n_x=3, n_t=2)
    s = build(ds.fields, ds.coords[0], ds.coords[-1], ["u_t", "u_xx", "u"], pa,
              p=12, sigma=ds.sigma_rep)
    assert len(s.A) > 0
    # u_t = -u for this field: the declared band must cover the residual
    resid = np.abs(s.A[:, 0] + s.A[:, 2])
    assert np.all(resid <= s.declared_epsilon("u_t", coeff_max=1.0,
                                              sigma=ds.sigma_rep))
    assert "synthetic" in ds.summary()
