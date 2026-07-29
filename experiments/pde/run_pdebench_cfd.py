"""PDEBench 1-D CFD as a SYSTEM (dev; docs/CASE_STUDY_PDEBENCH.md).

    .venv/bin/python experiments/pde/run_pdebench_cfd.py <extract.h5>

Three shipped fields (density, Vx, pressure) and two conservation laws:

    rho_t + d_x(rho u)        = 0                      exact, no viscosity
    (rho u)_t + d_x(rho u^2 + p) = (4 eta/3 + zeta) u_xx

so this is the first PDEBench family that needs `lagh/pdesystem.py` rather than
the scalar path, and the first where the shipped fields contain shocks.

**The declared field error is REPORTED, not chosen.** CFD has no closed-form
solution and this program's integrator does not solve Euler, so the trick used
for advection (compare against the exact solution) and for Burgers (compare
against an independent solve) is unavailable. What is available is a statement
of the same fact from the other side: SCAN the declaration and report the
smallest field error under which the file's own stated law sits inside its own
band. That number is a lower bound on the shipped data's error, it is
falsifiable, and unlike a tuned band it is the RESULT rather than an input to
one.
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
from experiments.pde.run_pdebench import default_scales           # noqa: E402
from lagh.pdesystem import (assemble, conjoin, discover_equation,  # noqa: E402
                            truth_check, weakest)
from lagh.weakform import Term, multiscale_patches                # noqa: E402

OUT = Path("experiments/results/pdebench.json")
P_BUMP = 16

# The REGISTERED vocabulary for 1-D compressible flow. A list, not a generated
# cross-product over three fields: |H| enters alpha directly.
TERMS = [
    Term("rho_t", 0, 1, "rho"), Term("(rho u)_t", 0, 1, "rho*u"),
    Term("(rho u)_x", 1, 0, "rho*u"), Term("(rho u^2)_x", 1, 0, "rho*u**2"),
    Term("p_x", 1, 0, "p"), Term("rho_x", 1, 0, "rho"),
    Term("u_x", 1, 0, "u"), Term("u_xx", 2, 0, "u"),
    Term("rho_xx", 2, 0, "rho"), Term("p_xx", 2, 0, "p"),
    Term("rho", 0, 0, "rho"), Term("u", 0, 0, "u"), Term("p", 0, 0, "p"),
    Term("1", 0, 0, "1")]


def truth_for(eta, zeta):
    """The stated laws. Continuity is exact; the momentum equation carries the
    1-D viscous combination (4 eta / 3 + zeta) u_xx."""
    return {"rho_t": {"(rho u)_x": -1.0},
            "(rho u)_t": {"(rho u^2)_x": -1.0, "p_x": -1.0,
                          "u_xx": 4.0 * eta / 3.0 + zeta}}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--samples", type=int, default=4)
    ap.add_argument("--eta", type=float, default=0.01)
    ap.add_argument("--zeta", type=float, default=0.01)
    ap.add_argument("--field-err", type=float, default=None,
                    help="declare it directly; omitted means SCAN and report "
                         "the smallest declaration the stated law needs")
    ap.add_argument("--speed", type=float, default=2.0)
    ap.add_argument("--tag", default="cfd")
    a = ap.parse_args(argv)

    t0 = time.time()
    data, reports = [], []
    for s in range(a.samples):
        ds = PB.load(a.path, layout="1d_cfd", sample=s)
        rep = PB.check_geometry(ds)
        reports.append({"sample": s, "summary": ds.summary(), **rep})
        data.append(ds)
    noise = PB.declared_noise(data[0])
    print("== pre-flight")
    for r in reports:
        print("  ", r["summary"])
        for n in r["notes"]:
            print("     NOTE:", n)
    if not all(r["ok"] for r in reports):
        print("REFUSED: the grid is not what the factory assumes")
        return 2

    sigma = noise["sigma"]
    cx, ct = data[0].coords[0], data[0].coords[-1]
    scales = default_scales(len(cx), len(ct), speed=a.speed,
                            dx=float(cx[1] - cx[0]), dt=float(ct[1] - ct[0]))
    sols = [(d.fields, d.coords) for d in data]
    truth = truth_for(a.eta, a.zeta)
    res = {"path": str(a.path), "family": "1d_cfd", "declared_noise": noise,
           "patch_scales": [list(s) for s in scales], "eta": a.eta,
           "zeta": a.zeta, "stated_laws": truth, "geometry": reports}

    base = assemble(sols, TERMS, lambda c: multiscale_patches(
        c[0], c[-1], scales, n_x=4, n_t=3), sigma=sigma, field_err=0.0,
        p=P_BUMP)

    def build(fe):
        # only the deterministic band term depends on the declaration, so the
        # patch family is integrated ONCE and re-banded per candidate
        return None if base is None else base.rebanded(fe)

    # --- the scan: how much declared field error does the STATED law need?
    if a.field_err is None:
        need = {}
        for fe in [0.0] + [10.0 ** e for e in range(-6, 1)]:
            rows = build(fe)
            if rows is None:
                continue
            for tgt in truth:
                tc = truth_check(rows, tgt, truth[tgt], sigma=sigma)
                if tgt not in need and tc["truth_certifies"]:
                    need[tgt] = fe
                print(f"   field_err {fe:8.1e}  {tgt:12s} "
                      f"truth/band {tc['truth_max_ratio']:10.3g} "
                      f"signal/band {tc['signal_to_band']:9.3g}"
                      f"{'  <- certifies' if tc['truth_certifies'] else ''}"
                      f"{'  VACUOUS' if tc['vacuous'] else ''}", flush=True)
            if len(need) == len(truth):
                break
        res["field_err_required"] = need
        fe = max(need.values()) if len(need) == len(truth) else None
        print(f"== smallest declaration the stated laws need: {need}")
        if fe is None:
            print("== REFUSED: no declaration in the scan makes both stated "
                  "laws hold without going vacuous")
            res["verdict"] = "refused-no-declaration-works"
            _save(res, a.tag, t0)
            return 0
    else:
        fe = a.field_err

    rows = build(fe)
    print(f"== rows {len(rows.A)} over {rows.n_solutions} samples, "
          f"{rows.rejected} patches rejected, field_err {fe:g}")
    eqs = []
    for tgt in truth:
        tc = truth_check(rows, tgt, truth[tgt], sigma=sigma)
        eq = discover_equation(rows, tgt, sigma=sigma, max_tier=3)
        eq["truth"] = tc
        eqs.append(eq)
        print(f"   {tgt:12s} {'CERT   ' if eq['certified'] else 'ABSTAIN'} "
              f"{str(eq.get('abstain') or ''):12s} "
              f"truth/band {tc['truth_max_ratio']:.3g} "
              f"vacuous={tc['vacuous']} {str(eq.get('law', ''))[:60]}")
        if eq.get("certified"):
            print(f"        intervals {eq['intervals']}")
    cert = conjoin(eqs)
    res.update({"field_err_used": fe, "n_rows": int(len(rows.A)),
                "patches_rejected": int(rows.rejected), "equations": eqs,
                "alpha_log10_total": cert.alpha_log10_total,
                "weakest": weakest(eqs),
                "system_certified": bool(all(e.get("certified") for e in eqs))})
    _save(res, a.tag, t0)
    return 0


def _save(res, tag, t0):
    res["seconds"] = round(time.time() - t0, 1)
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    prev[tag] = res
    OUT.write_text(json.dumps(prev, indent=1))


if __name__ == "__main__":
    raise SystemExit(main())
