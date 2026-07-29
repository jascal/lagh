"""What declaration does a PDEBench file's stated law need in the WEAK FORM?

    .venv/bin/python experiments/pde/run_weakform_declaration.py <file.h5> \
        --family burgers --speed 2 --truth '{"u_xx": 0.0031831, "u*u_x": -1}'

The tool the STRATEGY rule "a declaration must say what quantity it bounds" asks
for. It assembles rows with EXACTLY the runner's configuration (same term
library, same patch family, same bump order, same sample count) and bisects the
declared field error to where the stated law exactly meets its band. That number
is what `run_pdebench.py --field-err` should be given; the POINTWISE deviation a
solver-error measurement reports belongs in `--forecast-err` instead, and on
advection the two differ by ~3900x.

It uses the stated law, so it is CIRCULAR as a discovery protocol and is a dev
instrument only -- `docs/CASE_STUDY_PDEBENCH.md` reason 2 for why this target is
dev. Where a second registered relation exists, prefer the declaration derived
from THAT one (`run_conservation_floor.py`): it is independent of the law being
certified, and the rule is to borrow the LOOSEST requirement among the relations
scanned, never the tightest.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import pdebench as PB                          # noqa: E402
from experiments.pde.declaration import required_declaration        # noqa: E402
from experiments.pde.run_pdebench import (FAMILIES, FAMILY, P_BUMP,  # noqa: E402
                                          default_scales)
from lagh.pdesystem import assemble                                 # noqa: E402
from lagh.weakform import LIBRARY, multiscale_patches               # noqa: E402

OUT = Path("experiments/results/weakform_declaration.json")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--family", default="generic-1d", choices=list(FAMILIES))
    ap.add_argument("--layout", default="1d_single")
    ap.add_argument("--samples", type=int, default=4)
    ap.add_argument("--speed", type=float, default=0.0)
    ap.add_argument("--truth", required=True,
                    help="JSON dict of the stated law, target u_t")
    ap.add_argument("--pointwise", type=float, default=None,
                    help="the pointwise solver-error measurement, for the "
                         "over-declaration ratio")
    ap.add_argument("--tag", default=None)
    a = ap.parse_args(argv)

    data = [PB.load(a.path, layout=a.layout, sample=s, field_err=0.0)
            for s in range(a.samples)]
    sigma = PB.declared_noise(data[0])["sigma"]
    cx, ct = data[0].coords[0], data[0].coords[-1]
    scales = default_scales(len(cx), len(ct), speed=a.speed,
                            dx=float(cx[1] - cx[0]), dt=float(ct[1] - ct[0]))
    terms = [LIBRARY[n] for n in FAMILIES[a.family]]
    base = assemble([(d.fields, d.coords) for d in data], terms,
                    lambda c: multiscale_patches(c[0], c[-1], scales, **FAMILY),
                    sigma=sigma, field_err=0.0, p=P_BUMP)
    law = json.loads(a.truth)
    print(f"rows {len(base.A)}  scales {scales}  sigma_rep {sigma:.3e}")
    exact, grid_point, scan = required_declaration(base, "u_t", law, sigma)
    for r in scan:
        if r["truth_over_band"] < 1e3:
            print(f"   field_err {r['field_err']:8.1e}  "
                  f"truth/band {r['truth_over_band']:10.3g}  "
                  f"signal/band {r['signal_over_band']:9.3g}"
                  f"{'   <- holds' if r['holds'] else ''}"
                  f"{'  VACUOUS' if r['vacuous'] else ''}")
    if exact is None:
        print("   -> the stated law NEVER meets its band on this grid")
        return 1
    print(f"\n   WEAK-FORM REQUIREMENT {exact:.3e}  ({exact/sigma:.0f}x sigma_rep)"
          f"   [a decadal scan would have said {grid_point:g}]")
    res = {"path": a.path, "family": a.family, "law": law, "n_rows":
           int(len(base.A)), "sigma_rep": sigma, "patch_scales":
           [list(s) for s in scales], "required_declaration": exact,
           "grid_point_that_holds": grid_point, "scan": scan}
    if a.pointwise:
        res["pointwise_measurement"] = a.pointwise
        res["over_declaration"] = a.pointwise / exact
        print(f"   the pointwise measurement {a.pointwise:g} would "
              f"OVER-DECLARE this band by {a.pointwise / exact:.0f}x")
    prev = json.loads(OUT.read_text()) if OUT.exists() else {}
    prev[a.tag or f"{Path(a.path).stem}:{a.family}"] = res
    OUT.write_text(json.dumps(prev, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
