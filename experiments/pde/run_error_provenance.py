"""Score the error-provenance characterizer against known ground truth
(docs/DIRECTION_ERROR_PROVENANCE.md, predictions P1-P4).

Two cases where the answer is known independently:

  * PDEBench advection -- a shipped field whose deviation from the EXACT solution
    of its own stated law was measured today at 8.9e-4 to 2.0e-2 relative,
    dispersive, growing with t. Pipeline error, and P2 says the characterizer
    must say so and recover c3 within a factor of 2 of 1.09e-7.
  * this program's own C1 fields -- exact analytic solutions plus noise at a
    DECLARED sigma. Observational error, and P1 says no higher-derivative term
    may explain more than half the residual variance.

P3 (undetermined without evidence) and P4 (never call a deterministic error
stochastic) are scored on the same runs.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import fields as F                          # noqa: E402
from experiments.pde import pdebench as PB                       # noqa: E402
from experiments.pde.run_pdebench import default_scales          # noqa: E402
from lagh.errormodel import characterize                         # noqa: E402
from lagh.pdesystem import assemble                              # noqa: E402
from lagh.weakform import LIBRARY, multiscale_patches            # noqa: E402

OUT = Path("experiments/results/error_provenance.json")
HIGHER = ["u_xx", "u_xxx", "u_xxxx"]


def rows_for(sols, terms, sigma, scales, p=16):
    return assemble(sols, [LIBRARY[n] for n in terms],
                    lambda c: multiscale_patches(c[0], c[-1], scales,
                                                 n_x=4, n_t=3),
                    sigma=sigma, p=p)


def pdebench_case():
    data = [PB.load("data/pdebench/1D_Advection_beta0.7_x6.h5", sample=s)
            for s in range(4)]
    sigma = data[0].sigma_rep
    cx, ct = data[0].coords[0], data[0].coords[-1]
    sc = default_scales(len(cx), len(ct), speed=0.7,
                        dx=float(cx[1] - cx[0]), dt=float(ct[1] - ct[0]))
    rows = rows_for([(d.fields, d.coords) for d in data],
                    ["u_t", "u_x", "u_xx", "u_xxx", "u_xxxx", "u", "1"],
                    sigma, sc)
    n = rows.names
    resid = rows.A[:, n.index("u_t")] + 0.7 * rows.A[:, n.index("u_x")]
    cols = {k: rows.A[:, n.index(k)] for k in HIGHER}
    times = None            # assemble() pools solutions and does not keep the
                            # patch centres; the growth test is exercised in the
                            # unit tests instead
    scale = float(np.median(np.abs(rows.A[:, n.index("u_t")])))
    return characterize(resid, columns=cols, times=times, field_scale=scale)


def own_case(sigma=1e-4):
    """C1's own fields: exact analytic heat plus noise at a DECLARED sigma."""
    from experiments.pde.run_c0 import heat_solutions
    rng = np.random.default_rng(0)
    sols = []
    for u, x, t in heat_solutions(4):
        sols.append(({"u": u + rng.normal(0, sigma, u.shape)}, (x, t)))
    rows = rows_for(sols, ["u_t", "u_xx", "u_xxx", "u_xxxx", "u_x", "u", "1"],
                    sigma, [(24, 12), (32, 16), (40, 20)])
    n = rows.names
    resid = rows.A[:, n.index("u_t")] - 0.1 * rows.A[:, n.index("u_xx")]
    cols = {k: rows.A[:, n.index(k)] for k in ("u_xxx", "u_xxxx")}
    scale = float(np.median(np.abs(rows.A[:, n.index("u_t")])))
    return characterize(resid, columns=cols, field_scale=scale)


def main():
    res = {}
    print("== P2: PDEBench advection (known pipeline error, dispersive)")
    r = pdebench_case()
    res["pdebench_advection"] = r.__dict__
    print("  ", r.one_line())
    me = r.modified_equation
    if me.get("terms"):
        for k, v in sorted(me["terms"].items(),
                           key=lambda kv: -kv[1]["variance_explained"]):
            print(f"     {k:8s} explains {v['variance_explained']:6.1%}  "
                  f"coefficient {v['coefficient']: .4g}")
    print("     P2 met:", bool(me.get("best") == "u_xxx"
                               and 0.5e-7 < me["terms"]["u_xxx"]["coefficient"]
                               < 2.2e-7))

    print("== P1: our own C1 fields (declared sigma, observational)")
    r2 = own_case()
    res["own_c1_declared_sigma"] = r2.__dict__
    print("  ", r2.one_line())
    me2 = r2.modified_equation
    for k, v in (me2.get("terms") or {}).items():
        print(f"     {k:8s} explains {v['variance_explained']:6.1%}")
    print("     P1 met:", bool(not me2.get("explains")))

    print("== P3: nothing measured -> undetermined")
    r3 = characterize(np.random.default_rng(0).normal(0, 1, 50))
    res["no_evidence"] = r3.__dict__
    print("  ", r3.verdict, "|", r3.notes[0][:80])
    print("     P3 met:", r3.verdict == "undetermined")

    print("== P4: the dangerous direction -- deterministic never called stochastic")
    ok = (res["pdebench_advection"]["verdict"] == "structured-deterministic")
    print("     P4 met:", ok)
    OUT.write_text(json.dumps(res, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
