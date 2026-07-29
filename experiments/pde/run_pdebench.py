"""Point the weak-form instrument at a PDEBench file (docs/PDEBENCH_READINESS.md).

    .venv/bin/python experiments/pde/run_pdebench.py <file.h5> [--family burgers]

Does the pre-flight in order and REFUSES rather than guessing: load, check the
geometry, print the declared error model, assemble weak-form rows over several
SAMPLES (certifying on a held-out one), discover, then forecast-verify. Every
number that enters the band is either computed here or declared here; nothing is
assumed silently.

The registered term library per family lives in FAMILIES below. It is a
REGISTERED list, not a generated cross-product: |H| enters alpha directly.
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
from experiments.pde.verify import verify                         # noqa: E402
from lagh.pdesystem import (assemble, conjoin, discover_equation,  # noqa: E402
                            truth_check, weakest)
from lagh.weakform import LIBRARY, multiscale_patches             # noqa: E402

OUT = Path("experiments/results/pdebench.json")

# Registered libraries per PDE family. The target is always the time derivative.
FAMILIES = {
    "advection": ["u_t", "u_x", "u_xx", "u*u_x", "u", "1"],
    "burgers": ["u_t", "u_xx", "u*u_x", "u_x", "u", "1"],
    "diffusion-reaction": ["u_t", "u_xx", "u^3", "u^2", "u", "u_x", "1"],
    "generic-1d": ["u_t", "u_xx", "u_xxx", "u*u_x", "u_x", "u", "1"],
}
SCALES = [(24, 12), (32, 16), (40, 20)]
FAMILY = dict(n_x=4, n_t=3)
P_BUMP = 16


def load_samples(path, layout, n_samples, tmax_index, field_err):
    """The pre-flight: load each sample, CHECK its geometry, and refuse the run
    if any axis is not what the factory assumes."""
    out, reports = [], []
    for s in range(n_samples):
        ds = PB.load(path, layout=layout, sample=s, field_err=field_err,
                     tmax_index=tmax_index)
        rep = PB.check_geometry(ds)
        reports.append({"sample": s, "summary": ds.summary(), **rep})
        out.append(ds)
    return out, reports


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--family", default="generic-1d", choices=list(FAMILIES))
    ap.add_argument("--layout", default="1d_single")
    ap.add_argument("--samples", type=int, default=4)
    ap.add_argument("--tmax-index", type=int, default=None)
    ap.add_argument("--field-err", type=float, default=0.0,
                    help="DECLARED solver error; it cannot be measured from "
                         "shipped data and the record says so")
    ap.add_argument("--extra-sigma", type=float, default=0.0,
                    help="any declared measurement noise beyond float32 storage")
    ap.add_argument("--truth", default=None,
                    help="JSON dict of the generating law, when the filename "
                         "states it: runs the truth check BEFORE any verdict")
    ap.add_argument("--tag", default=None)
    a = ap.parse_args(argv)

    t0 = time.time()
    # one extra sample, reserved for the FORECAST: C2's discipline is that the
    # verify track runs on data no stage of the pipeline has seen, and the
    # certification holdout sample has already been seen by certification
    data, reports = load_samples(a.path, a.layout, a.samples + 1, a.tmax_index,
                                 a.field_err)
    data, forecast_sample = data[:a.samples], data[a.samples]
    noise = PB.declared_noise(data[0], extra_sigma=a.extra_sigma)
    print("== pre-flight")
    for r in reports:
        print("  ", r["summary"])
        for n in r["notes"]:
            print("     NOTE:", n)
    print("   declared error model:", json.dumps(noise, indent=None)[:400])
    if not all(r["ok"] for r in reports):
        print("REFUSED: the grid is not what the factory assumes (see notes)")
        return 2

    sigma = noise["sigma"]
    terms = [LIBRARY[n] for n in FAMILIES[a.family]]
    sols = [({k: v for k, v in ds.fields.items()}, ds.coords) for ds in data]
    rows = assemble(sols, terms, lambda c: multiscale_patches(
        c[0], c[-1], SCALES, **FAMILY), sigma=sigma, field_err=a.field_err,
        p=P_BUMP)
    if rows is None:
        print("ABSTAIN: every patch failed the resolution/aliasing gate -- the "
              "grid does not represent this field")
        return 0
    print(f"== rows {len(rows.A)} over {rows.n_solutions} samples, "
          f"{rows.rejected} patches rejected, vocabulary {rows.names}")

    res = {"path": str(a.path), "family": a.family, "declared_noise": noise,
           "geometry": reports, "n_rows": int(len(rows.A)),
           "patches_rejected": int(rows.rejected),
           "vocabulary": rows.names, "n_samples": rows.n_solutions}

    if a.truth:
        truth = json.loads(a.truth)
        res["truth_check"] = truth_check(rows, "u_t", truth, sigma=sigma)
        print("== truth check (before any verdict):", res["truth_check"])

    eq = discover_equation(rows, "u_t", sigma=sigma, max_tier=3)
    res["equation"] = eq
    cert = conjoin([eq])
    res["alpha_log10_total"] = cert.alpha_log10_total
    res["weakest"] = weakest([eq])
    print(f"== {'CERTIFIED' if eq['certified'] else 'ABSTAIN'} "
          f"{eq.get('abstain') or ''} {eq.get('law', '')}")
    if eq.get("certified"):
        print("   intervals:", eq["intervals"])
        print(f"   alpha <= 1e{eq['alpha_log10']:.0f}   "
              f"signal/band {eq['median_signal_to_band']:.3g}")
        # forecast-verify on a sample NO stage has seen
        ds = forecast_sample
        x, t = ds.coords[0], ds.coords[-1]
        ivs = {k: tuple(v) if v else (eq["coefficients"][k],) * 2
               for k, v in eq["intervals"].items()}
        try:
            # ETD, not explicit RK45: on a PDEBench-scale grid the explicit
            # scheme is diffusion-limited (512 points at nu = 0.2 needs ~250k
            # steps per trajectory and the run never finishes). The exponential
            # integrator solves the linear part exactly and declares its error
            # on a SUBSTEP ladder instead of a tolerance ladder.
            v = verify(ds.fields["u"], ds.fields["u"][:, 0], x, t, ivs,
                       sigma=sigma, scheme="etd", nsub=64)
        except Exception as e:                                 # noqa: BLE001
            v = {"verified": False, "refusal": f"verify raised: {e}"}
        res["verify"] = v
        print(f"   verify: {'OK' if v.get('verified') else 'FAIL'} "
              f"outside={v.get('n_outside')}/{v.get('n_points')} "
              f"{v.get('refusal') or ''}")
    res["seconds"] = round(time.time() - t0, 1)
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    prev[a.tag or f"{Path(a.path).stem}:{a.family}"] = res
    OUT.write_text(json.dumps(prev, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
