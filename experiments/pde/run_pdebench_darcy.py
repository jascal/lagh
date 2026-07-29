"""PDEBench 2-D Darcy flow: the first STEADY-STATE claim in the arc (dev).

    .venv/bin/python experiments/pde/run_pdebench_darcy.py

Everything the weak-form arc has certified so far is an EVOLUTION equation: a
time-derivative target, time as the last axis, and a verify track that
integrates forward. Darcy has no time axis at all --

    -div(a(x) grad u(x)) = beta,

an elliptic boundary-value problem, shipped as (a, u) pairs on a 128x128 grid
with x- and y-coordinates and no t-coordinate. Steady states, equilibria and
constitutive laws are a large share of real science and this program has never
tested whether it can express such a claim.

TWO results are available here and they are different in kind.

**The reach boundary.** The general variable-coefficient equation is NOT in this
factory's reach, and the reason is precise rather than incidental. The library is
d^alpha(g(fields)) with g POINTWISE; `div(a grad u)` is divergence form at the
outer level but its integrand `a grad u` is not pointwise -- it pairs a field
with a DERIVATIVE of another field. By-parts once leaves int(grad phi . a grad u),
still a data derivative; moving it again gives int(div(a grad phi) u), which
needs grad a -- a derivative of measured data, which is exactly what the weak
form exists to avoid. Rearranging does not escape it:
(a u_x)_x = (a u)_xx - a_xx u - a_x u_x. So a variable coefficient that is ITSELF
DATA breaks the arc's central guarantee. That bears directly on the registered
next step "(c) variable coefficients" in DIRECTION_PDE.md: it needs either a
declared error model for grad a or a mixed formulation, and it is not free.

**What IS certifiable, on a stated domain.** PDEBench's coefficient is binary
(measured: a in {0.1, 1.0}), so wherever a is locally CONSTANT the equation
collapses to `a laplacian(u) = -beta`, i.e.

    u_xx + u_yy + beta * (1/a) = 0

and every term there is divergence form with a pointwise g. The patches where
that applies are selected FROM THE DATA -- a is constant on the window -- which
makes the restriction a domain statement the certificate carries, not a filter
that hides anything. Patches straddling a conductivity interface are reported as
excluded, with the count.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import pdebench as PB                        # noqa: E402
from lagh.pdesystem import (assemble, conjoin, discover_equation,  # noqa: E402
                            truth_check)
from lagh.weakform import Term, make_patches_nd                   # noqa: E402

OUT = Path("experiments/results/pdebench.json")
API = "https://darus.uni-stuttgart.de/api/access/datafile/{id}"
P_BUMP = 12

# The REGISTERED vocabulary for the constant-coefficient interior. Two spatial
# axes and NO time axis, so every alpha is stated in full.
TERMS = [
    Term("u_xx", alpha=(2, 0), gexpr="u"), Term("u_yy", alpha=(0, 2), gexpr="u"),
    Term("u_x", alpha=(1, 0), gexpr="u"), Term("u_y", alpha=(0, 1), gexpr="u"),
    Term("1/a", alpha=(0, 0), gexpr="1/a"), Term("a", alpha=(0, 0), gexpr="a"),
    Term("u", alpha=(0, 0), gexpr="u"), Term("a*u", alpha=(0, 0), gexpr="a*u"),
    Term("1", alpha=(0, 0), gexpr="1")]


def fetch(file_id, samples, out: Path):
    """Darcy ships (nu, tensor) with x/y coordinates and NO time axis, so it does
    not go through the general fetcher -- which assumes one."""
    import fsspec
    import h5py
    if out.exists():
        return out
    fs = fsspec.filesystem("http", client_kwargs={"trust_env": True})
    with h5py.File(fs.open(API.format(id=file_id), block_size=2 ** 20), "r") as h:
        a = np.asarray(h["nu"][:samples], np.float32)
        u = np.asarray(h["tensor"][:samples, 0], np.float32)
        x = np.asarray(h["x-coordinate"][()], np.float32)
        y = np.asarray(h["y-coordinate"][()], np.float32)
    out.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out, "w") as g:
        g.create_dataset("nu", data=a)
        g.create_dataset("tensor", data=u)
        g.create_dataset("x-coordinate", data=x)
        g.create_dataset("y-coordinate", data=y)
        g.attrs["source_file_id"] = int(file_id)
        g.attrs["source_doi"] = PB.__dict__.get("DOI", "doi:10.18419/darus-2986")
        g.attrs["notes"] = json.dumps(
            ["steady state: the file has x- and y-coordinates and NO "
             "t-coordinate; the tensor's singleton second axis is a placeholder"])
    return out


def constant_a_patches(a, x, y, half, counts, *, phase=None, min_dist=0.0):
    """Patches whose window lies INSIDE one conductivity region, optionally in a
    NAMED phase and at a stated distance from the nearest interface.

    Every part of that is a domain statement the certificate carries. Two of them
    were forced by measurement:

    * `phase`, because the stated law holds in one phase and not the other. With
      beta = 0.1 the high-conductivity interior gives an implied beta of 0.1000
      (IQR 7e-4) while the low-conductivity interior gives 0.135-0.163 at EVERY
      distance from an interface -- so it is not interface smearing, and the
      honest response is to name the phase rather than widen the band until both
      fit.
    * float32 comparison, because `a` arrives as float32 and float32(0.1) is not
      float64(0.1). An `a == 0.1` test silently keeps only the a = 1 patches --
      measured, in a first pass of this very diagnostic.
    """
    from scipy.ndimage import distance_transform_edt
    keep, drop = [], 0
    for pa in make_patches_nd((x, y), half, counts):
        w = a[pa.idx[0], pa.idx[1]]
        if float(w.max() - w.min()) > 0.0:
            drop += 1
            continue
        if phase is not None and not np.isclose(float(w.flat[0]), phase,
                                                rtol=1e-5):
            drop += 1
            continue
        if min_dist > 0:
            d = distance_transform_edt(a == w.flat[0])
            if float(d[pa.idx[0], pa.idx[1]].min()) < min_dist:
                drop += 1
                continue
        keep.append(pa)
    return keep, drop


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", type=int, default=133218)      # beta = 0.1
    ap.add_argument("--beta", type=float, default=0.1)
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--half", default="16,16")
    ap.add_argument("--counts", default="6,6")
    ap.add_argument("--field-err", type=float, default=0.0)
    ap.add_argument("--phase", type=float, default=1.0,
                    help="conductivity phase the claim is restricted to")
    ap.add_argument("--min-dist", type=float, default=8.0,
                    help="cells from the nearest interface")
    ap.add_argument("--tag", default="darcy_beta0.1")
    a_ = ap.parse_args(argv)

    t0 = time.time()
    path = fetch(a_.id, a_.samples + 1,
                 Path(f"data/pdebench/2D_Darcy_beta{a_.beta:g}_x{a_.samples+1}.h5"))
    import h5py
    with h5py.File(path, "r") as h:
        A = np.asarray(h["nu"][()], float)
        U = np.asarray(h["tensor"][()], float)
        x = np.asarray(h["x-coordinate"][()], float)
        y = np.asarray(h["y-coordinate"][()], float)
    # the coordinates are float32 like every other PDEBench axis
    x = np.linspace(x[0], x[-1], len(x))
    y = np.linspace(y[0], y[-1], len(y))
    vals = np.unique(A)
    print(f"== 2-D Darcy, steady state: no time axis. a takes {len(vals)} values "
          f"{vals[:4]}; grid {len(x)}x{len(y)}; {len(A)} samples")
    sigma = PB.sigma_rep_for(np.array([np.max(np.abs(U))]))[0]
    print(f"   declared sigma_rep (float32 storage, at the field peak): {sigma:.3e}")

    half = tuple(int(v) for v in a_.half.split(","))
    counts = tuple(int(v) for v in a_.counts.split(","))
    sols, dropped = [], 0
    for s in range(len(A)):
        fields = {"a": A[s], "u": U[s]}
        keep, drop = constant_a_patches(A[s], x, y, half, counts,
                                        phase=a_.phase, min_dist=a_.min_dist)
        dropped += drop
        if keep:
            sols.append((fields, (x, y), keep))
    print(f"   patches in phase a={a_.phase:g}, at least {a_.min_dist:g} cells "
          f"from an interface: {sum(len(k) for _, _, k in sols)}; excluded "
          f"(other phase, straddling, or too close): {dropped}")
    if len(sols) < 2:
        print("ABSTAIN: fewer than two samples retain any interior patch")
        return 0

    # the patch family differs PER SAMPLE, because which patches lie inside one
    # conductivity region depends on that sample's own a
    fams = iter([k for _, _, k in sols])
    rows = assemble([(f, c) for f, c, _ in sols], TERMS,
                    lambda c, _f=fams: next(_f), sigma=sigma,
                    field_err=a_.field_err, p=P_BUMP)
    if rows is None:
        print("ABSTAIN: every patch failed the resolution gate")
        return 0
    print(f"== rows {len(rows.A)} over {rows.n_solutions} samples, "
          f"{rows.rejected} rejected by the resolution gate, "
          f"vocabulary {rows.names}")

    truth = {"u_yy": -1.0, "1/a": -a_.beta}          # u_xx = -u_yy - beta/a
    tc = truth_check(rows, "u_xx", truth, sigma=sigma)
    print("== truth check (before any verdict):", tc)
    eq = discover_equation(rows, "u_xx", sigma=sigma, max_tier=3,
                           features=[n for n in rows.names if n != "u_xx"])
    print(f"== {'CERTIFIED' if eq['certified'] else 'ABSTAIN'} "
          f"{eq.get('abstain') or ''} {eq.get('law', '')}")
    if eq.get("certified"):
        print("   intervals:", eq["intervals"])
        print(f"   alpha <= 1e{eq['alpha_log10']:.0f}  "
              f"signal/band {eq['median_signal_to_band']:.3g}")
    res = {"path": str(path), "family": "2d_darcy", "beta": a_.beta,
           "phase": a_.phase, "min_dist_cells": a_.min_dist,
           "steady_state": True, "a_values": [float(v) for v in vals],
           "sigma_rep": sigma, "n_rows": int(len(rows.A)),
           "patches_excluded_straddling": int(dropped),
           "patches_rejected_resolution": int(rows.rejected),
           "vocabulary": rows.names, "stated_law": truth,
           "truth_check": tc, "equation": eq,
           "alpha_log10_total": conjoin([eq]).alpha_log10_total,
           "domain": ("patches lying strictly inside one conductivity region; "
                      "on a patch straddling an interface div(a grad u) is not "
                      "a laplacian and the stated law is false"),
           "reach_note": ("the GENERAL variable-coefficient equation is out of "
                          "reach: a grad u is not a pointwise g, and removing "
                          "grad u introduces grad a -- a derivative of measured "
                          "data"),
           "seconds": round(time.time() - t0, 1)}
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    prev[a_.tag] = res
    OUT.write_text(json.dumps(prev, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
