"""PDEBench adapter: external field data -> a weak-form claim with a DECLARED
error model (docs/PDEBENCH_READINESS.md).

Nothing here fits or certifies. Its whole job is the step the program treats as
non-negotiable before any external data is touched: turn a downloaded array into
observations whose error is DECLARED rather than assumed. For PDEBench that has
three parts, and the first is the one the arc has already been bitten by:

1. **Storage quantization.** PDEBench ships float32 HDF5. A float32 value is not
   the number the simulation produced: it is that number rounded, and the
   rounding is a real, computable error of size eps32 * |u| / 2. Certifying
   against it without declaring it is exactly the float32-quantization lesson
   from the LLM-SRBench read (declare sigma_rep). `sigma_rep` below is that
   declaration, and because it is PROPORTIONAL to |u| rather than absolute, it
   is reported both ways: the weak-form band takes one scalar sigma, so the
   conservative choice is sigma_rep at the field's own peak.

2. **Solver error.** Every PDEBench field is the output of a numerical solver
   whose error is NOT distributed with the data, and no ladder can be run after
   the fact -- the coarse levels were never saved. This is declared, not
   measured, and that difference is stated in the certificate: a run against
   PDEBench data carries `field_err` as an ASSUMPTION with a stated value, or it
   carries zero and says so.

3. **Grid geometry.** PDEBench's 1-D files are (n_samples, n_t, n_x) with
   periodic x on [0, 1] (or [-1, 1] for some sets) -- so the axis order, the
   domain length and the endpoint convention all have to be converted, and a
   wrong L silently rescales every derivative term. The loader returns the
   coordinate vectors it derived, and the runner prints them.

The reader is deliberately dependency-light (h5py only, imported lazily) and
every function works on an in-memory array, so the whole path is testable
without the download.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

FLOAT32_EPS = float(np.finfo(np.float32).eps)      # 1.19e-7


@dataclass
class Dataset:
    """One PDEBench sample, in this program's coordinates."""
    fields: dict                  # name -> array shaped (nx, [ny,] nt)
    coords: tuple                 # (x, [y,] t), time LAST
    name: str = ""
    sigma_rep: float = 0.0        # declared representation noise (absolute)
    sigma_rep_rel: float = 0.0    # ...and the relative form it came from
    field_err: float = 0.0        # DECLARED (not measured) solver error
    provenance: dict = field(default_factory=dict)

    def summary(self) -> str:
        shp = "x".join(str(len(c)) for c in self.coords)
        return (f"{self.name}: fields {sorted(self.fields)} grid {shp} "
                f"sigma_rep {self.sigma_rep:.3e} "
                f"(rel {self.sigma_rep_rel:.3e}) field_err {self.field_err:.3e}")


def sigma_rep_for(u, dtype=np.float32) -> tuple:
    """(absolute, relative) representation noise for a stored field.

    A float32 array carries a rounding error of at most eps32/2 relative per
    value; the weak-form band takes ONE declared sigma for the field, so the
    honest scalar is the bound at the field's own peak -- an over-declaration
    everywhere else, which is the safe direction.
    """
    u = np.asarray(u)
    rel = 0.5 * (FLOAT32_EPS if np.dtype(dtype) == np.float32
                 else float(np.finfo(dtype).eps))
    return float(rel * np.max(np.abs(u))), float(rel)


def regularize_axis(c) -> tuple:
    """(uniform grid through the endpoints, max deviation, relative-to-step).

    PDEBench stores its COORDINATES in float32, so a perfectly uniform grid
    arrives with the storage rounding in it: on a 101-point t-axis over [0, 1]
    that is ~6e-8 absolute, which is 6e-6 RELATIVE TO THE STEP. That matters,
    because the by-parts weights carry 1/Δt^k with k differing per term, so the
    error does NOT cancel out of the relation -- it would land on the u_t column
    at 6e-6 of its own size, a hundred times the float32 field noise.

    Banding it would be honest but wasteful. Regularizing is honest AND tight:
    the grid IS uniform (the file's own coordinate vector asserts a linspace),
    so we reconstruct it in float64 and REPORT the deviation from what was
    stored. Nothing is invented -- the deviation is the storage rounding, and it
    is exactly what a reader needs to check that claim."""
    c = np.asarray(c, float)
    u = np.linspace(c[0], c[-1], len(c))
    dev = float(np.max(np.abs(c - u)))
    step = abs(float(u[1] - u[0])) if len(u) > 1 else 1.0
    return u, dev, dev / step if step else np.inf


def from_arrays(fields: dict, x, t, *, y=None, name="", dtype=np.float32,
                field_err=0.0, provenance=None, regularize=True) -> Dataset:
    """Build a Dataset from arrays already in memory -- the path the tests use,
    and the one a non-PDEBench source should also come through."""
    raw = (np.asarray(x, float),) + (
        () if y is None else (np.asarray(y, float),)) + (np.asarray(t, float),)
    coords, devs = [], []
    for c in raw:
        u, dev, rel_step = regularize_axis(c) if len(c) > 2 else (c, 0.0, 0.0)
        # ONLY absorb a deviation consistent with float32 storage rounding. A
        # genuinely non-uniform axis must survive to check_geometry and be
        # REFUSED there -- regularizing it would replace the grid the data lives
        # on with a grid we invented, which is the opposite of declaring it.
        tol = 8.0 * FLOAT32_EPS * max(1.0, float(np.max(np.abs(c))))
        ok = bool(regularize and len(c) > 2 and dev <= tol)
        coords.append(u if ok else c)
        devs.append({"max_abs": dev, "rel_to_step": rel_step,
                     "regularized": ok, "tolerance": tol})
    fields = {k: np.asarray(v, float) for k, v in fields.items()}
    peak = max(float(np.max(np.abs(v))) for v in fields.values())
    abs_, rel = sigma_rep_for(np.array([peak]), dtype)
    prov = dict(provenance or {})
    prov["coord_regularized"] = bool(regularize)
    prov["coord_deviation"] = devs
    return Dataset(fields=fields, coords=tuple(coords), name=name,
                   sigma_rep=abs_, sigma_rep_rel=rel, field_err=field_err,
                   provenance=prov)


# PDEBench file layouts, by the group/dataset names their HDF5 files use. Keyed
# by a short name the runner passes; `axes` says how to get from the stored
# shape to (space..., time).
LAYOUTS = {
    # 1D_Advection / 1D_Burgers / 1D_diffusion-reaction: tensor "tensor"
    # shaped (n_samples, n_t, n_x), plus "x-coordinate" and "t-coordinate"
    "1d_single": {"tensor": "tensor", "x": "x-coordinate", "t": "t-coordinate",
                  "fields": {"u": None}, "order": "sxt->xt"},
    # 1D_CFD: separate "Vx", "density", "pressure" of the same shape
    "1d_cfd": {"x": "x-coordinate", "t": "t-coordinate",
               "fields": {"rho": "density", "u": "Vx", "p": "pressure"},
               "order": "sxt->xt"},
    # 2D_diffusion-reaction / 2D shallow water: (n_samples, n_t, n_x, n_y, c)
    "2d_single": {"tensor": "data", "x": "x-coordinate", "y": "y-coordinate",
                  "t": "t-coordinate", "fields": {"u": 0}, "order": "sxyt->xyt"},
}


def load(path, *, layout="1d_single", sample=0, fields=None, dtype=np.float32,
         field_err=0.0, tmax_index=None, regularize=True) -> Dataset:
    """Read one sample out of a PDEBench HDF5 file.

    Deliberately explicit about the layout rather than sniffing it: a wrong axis
    order or a wrong domain length silently rescales every derivative term, and
    a silent rescale is the kind of error a certificate would happily certify.
    """
    import h5py                                             # lazy: optional dep
    spec = dict(LAYOUTS[layout])
    if fields:
        spec["fields"] = fields
    out, prov = {}, {"path": str(path), "layout": layout, "sample": int(sample)}
    with h5py.File(path, "r") as f:
        keys = list(f.keys())
        prov["keys"] = keys[:20]
        x = np.asarray(f[spec["x"]][()], float).ravel()
        t = np.asarray(f[spec["t"]][()], float).ravel()
        ykey = spec.get("y")
        y = (np.asarray(f[ykey][()], float).ravel()
             if ykey and ykey in f else None)
        for name, where in spec["fields"].items():
            if spec.get("tensor") and where is None:
                a = np.asarray(f[spec["tensor"]][sample], float)
            elif spec.get("tensor"):
                a = np.asarray(f[spec["tensor"]][sample], float)
                a = a[..., where] if a.ndim > (3 if y is not None else 2) else a
            else:
                a = np.asarray(f[where][sample], float)
            # stored as (n_t, n_x[, n_y]) -> (n_x[, n_y], n_t)
            out[name] = np.moveaxis(a, 0, -1)
        if tmax_index:
            t = t[:tmax_index]
            out = {k: v[..., :tmax_index] for k, v in out.items()}
    return from_arrays(out, x, t, y=y, name=f"{layout}:{sample}", dtype=dtype,
                       field_err=field_err, provenance=prov,
                       regularize=regularize)


def check_geometry(ds: Dataset, *, periodic=True, rel_tol=1e-6) -> dict:
    """Report what the loaded grid actually is, and refuse silently-wrong cases.

    Checks the things that would rescale derivatives without any error message:
    uniform spacing per axis, and whether the periodic endpoint is duplicated (a
    periodic grid stored WITH its right endpoint makes every spectral wavenumber
    slightly wrong -- measured on our own stage-1 fields, where it put 18244 of
    41634 forecast points outside the band with the law exact).

    Two things this deliberately does NOT do:

    * It does not demand exact uniformity. PDEBench stores its COORDINATES in
      float32 too, so a perfectly uniform grid arrives with ~1e-7 relative
      jitter in its spacing. That jitter is real -- the quadrature's h inherits
      it -- so it is MEASURED and reported (`nonuniformity_rel`) and folded into
      the declaration rather than treated as a broken grid.
    * It does not decide endpoint duplication from the coordinates, which cannot
      tell you the intended domain length. It asks the FIELD: on a duplicated
      grid u(x_0) = u(x_-1) identically.
    """
    coord_eps = 4.0 * FLOAT32_EPS         # float32-stored coordinates
    out = {"axes": [], "ok": True, "notes": []}
    for i, c in enumerate(ds.coords):
        d = np.diff(c)
        nonu = float(np.max(np.abs(d - d[0])) / abs(d[0])) if d[0] else np.inf
        uniform = bool(nonu <= max(rel_tol, coord_eps))
        out["axes"].append({"n": int(len(c)), "min": float(c[0]),
                            "max": float(c[-1]), "step": float(d[0]),
                            "span": float(c[-1] - c[0]), "uniform": uniform,
                            "nonuniformity_rel": nonu})
        if not uniform:
            out["ok"] = False
            out["notes"].append(
                f"axis {i} is not uniformly spaced (relative jitter {nonu:.2e}): "
                "the patch quadrature and the by-parts weights both assume a "
                "uniform grid")
        elif nonu > 0:
            out["notes"].append(
                f"axis {i} spacing carries {nonu:.2e} relative jitter from "
                "float32 coordinate storage; the quadrature step inherits it")
    if periodic and len(ds.coords) >= 2:
        for k, v in ds.fields.items():
            v = np.asarray(v, float)
            gap = float(np.max(np.abs(v[0] - v[-1])))
            scale = float(np.max(np.abs(v))) + 1e-300
            if gap <= max(8 * ds.sigma_rep, 1e-6 * scale):
                out["notes"].append(
                    f"field {k!r}: u(x_first) == u(x_last) to {gap:.2e}, so the "
                    "periodic right endpoint is DUPLICATED; the spectral verify "
                    "track needs it dropped (L = n*dx)")
                break
        # ...and the opposite failure, which this check MISSED entirely until a
        # transmissive-boundary file walked through it (2026-07-29): a field that
        # is not periodic at all. The raw gap cannot say -- on an
        # endpoint-excluded grid a periodic field's seam is dx*|u_x|, measured
        # larger on our own C2 fields than on PDEBench advection -- so the seam
        # is expressed in units of an ORDINARY interior step.
        seams = {}
        for k, v in ds.fields.items():
            v = np.asarray(v, float)
            step = np.abs(np.diff(v, axis=0))
            seams[k] = float(np.median(np.abs(v[0] - v[-1]))
                             / max(np.quantile(step, 0.99), 1e-300))
        out["periodicity_seam"] = seams
        worst = max(seams, key=seams.get) if seams else None
        if worst and seams[worst] > 2.0:
            out["notes"].append(
                f"field {worst!r} is NOT PERIODIC: its wrap seam is "
                f"{seams[worst]:.2f}x an ordinary interior step (periodic "
                "fields measure 0.14-0.91). Weak-form certification is "
                "unaffected -- test functions vanish inside the domain -- but "
                "the spectral verify track does not apply and will refuse")
    for k, v in ds.fields.items():
        if not np.all(np.isfinite(v)):
            out["ok"] = False
            out["notes"].append(f"field {k!r} contains non-finite values")
    return out


def declared_noise(ds: Dataset, extra_sigma: float = 0.0) -> dict:
    """The error model a run against this dataset declares, as one dict the
    runner can print verbatim into its results.

    Kept separate from the loader so that the declaration is a visible, quotable
    object rather than a set of keyword arguments scattered through a script."""
    return {"sigma": float(np.hypot(ds.sigma_rep, extra_sigma)),
            "sigma_rep_absolute": ds.sigma_rep,
            "sigma_rep_relative": ds.sigma_rep_rel,
            "sigma_extra_declared": float(extra_sigma),
            "field_err_declared": ds.field_err,
            "field_err_is_measured": False,
            "coord_regularized": ds.provenance.get("coord_regularized"),
            "coord_deviation": ds.provenance.get("coord_deviation"),
            "notes": [
                "the coordinate axes are float32 too, so they are rebuilt as "
                "exact linspaces and the deviation from what was stored is "
                "reported (coord_deviation): unregularized, that rounding lands "
                "on the time-derivative column at ~1e-6 of its own size, "
                "because the by-parts weights carry 1/step^k with k per term",
                "sigma_rep is the float32 storage quantization at the field's "
                "peak: computed, not assumed",
                "field_err is a DECLARED solver-error assumption -- PDEBench "
                "does not ship the coarse levels a tolerance ladder would need, "
                "so it cannot be measured after the fact",
                "any certificate from this data is conditional on that "
                "declaration and says so"]}
