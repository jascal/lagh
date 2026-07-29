"""Does the tightest-held relation measure a scheme's numerical FLOOR?
(docs/DIRECTION_ERROR_PROVENANCE.md item 5, route 2.)

REGISTERED BEFORE RUNNING.

The band's magnitude is the one thing this program cannot yet obtain without an
independent reference: the free-fit route is measured to under-declare by four
orders (it absorbs the error it is measuring), an exact solution exists for only
one family, and an independent solve is blocked wherever this program cannot
integrate the equations -- which is exactly the CFD case where it was needed.

Route 2 needs neither. PDEBench CFD showed twice that a scheme enforces its
relations UNEQUALLY: continuity required a 1e-4 declaration against momentum's
1e-2 (viscous) and 1e-5 (inviscid). The hypothesis under test generalizes that
observation into a measurement:

    H: among several registered relations of one system, the relation the
       generating scheme enforces MOST TIGHTLY requires a declaration at or near
       the data's numerical FLOOR (storage precision), while the others require
       declarations at the SCHEME-ERROR scale. The tightest therefore estimates
       the floor, and the spread across relations estimates the scheme error.

Scored on PDEBench advection, the one family whose scheme error is INDEPENDENTLY
known (2.75e-2 relative, measured against the exact shift solution). Two
relations that both hold exactly for the true solution:

    R1   u_t      = -beta u_x          the equation the solver integrates
    R2  (u^2)_t   = -beta (u^2)_x      a derived conservation law

PREDICTIONS, fixed before the run:

  * H HOLDS if the tighter of {R1, R2} requires a declaration well below the
    known 2.75e-2 -- ideally near sigma_rep = 5e-8 -- and the looser requires
    ~2.75e-2.
  * H FAILS if the tighter one also requires ~2.75e-2, i.e. every registered
    relation is violated at the scheme-error scale and none of them sees the
    floor. In that case route 2 does not measure what it was meant to and the
    honest output is to say so.

A third outcome is possible and is NOT a pass: both relations tight, which would
mean the declaration is not being driven by the scheme error at all and the
diagnostic is measuring something else.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import pdebench as PB                       # noqa: E402
from experiments.pde.run_pdebench import default_scales          # noqa: E402
from lagh.pdesystem import assemble, truth_check                 # noqa: E402
from lagh.weakform import Term, multiscale_patches               # noqa: E402

OUT = Path("experiments/results/conservation_floor.json")
KNOWN_SCHEME_ERROR = 2.75e-2          # measured independently, advection beta=0.7
BETA = 0.7

TERMS = [
    Term("u_t", 0, 1, "u"), Term("u_x", 1, 0, "u"),
    Term("(u^2)_t", 0, 1, "u**2"), Term("(u^2)_x", 1, 0, "u**2"),
    Term("u_xx", 2, 0, "u"), Term("u", 0, 0, "u"), Term("1", 0, 0, "1")]

RELATIONS = {
    "R1  u_t = -b u_x       (the integrated equation)": ("u_t", {"u_x": -BETA}),
    "R2  (u^2)_t = -b (u^2)_x (derived conservation)": ("(u^2)_t",
                                                        {"(u^2)_x": -BETA}),
}


def main():
    data = [PB.load("data/pdebench/1D_Advection_beta0.7_x6.h5", sample=s)
            for s in range(4)]
    sigma = data[0].sigma_rep
    cx, ct = data[0].coords[0], data[0].coords[-1]
    sc = default_scales(len(cx), len(ct), speed=BETA,
                        dx=float(cx[1] - cx[0]), dt=float(ct[1] - ct[0]))
    base = assemble([(d.fields, d.coords) for d in data], TERMS,
                    lambda c: multiscale_patches(c[0], c[-1], sc, n_x=4, n_t=3),
                    sigma=sigma, field_err=0.0, p=16)
    print(f"rows {len(base.A)}  sigma_rep {sigma:.3e}  "
          f"independently measured scheme error {KNOWN_SCHEME_ERROR:.3e}")
    res = {"sigma_rep": sigma, "known_scheme_error": KNOWN_SCHEME_ERROR,
           "n_rows": int(len(base.A)), "relations": {}}
    need = {}
    for label, (target, law) in RELATIONS.items():
        row = []
        smallest = None
        for fe in [0.0] + [10.0 ** e for e in range(-8, 0)]:
            tc = truth_check(base.rebanded(fe), target, law, sigma=sigma)
            row.append((fe, tc["truth_max_ratio"], tc["signal_to_band"],
                        tc["truth_certifies"], tc["vacuous"]))
            if smallest is None and tc["truth_certifies"] and not tc["vacuous"]:
                smallest = fe
        need[label] = smallest
        res["relations"][label] = {
            "required_declaration": smallest,
            "scan": [{"field_err": f, "truth_over_band": r,
                      "signal_over_band": s, "holds": h, "vacuous": v}
                     for f, r, s, h, v in row]}
        print(f"\n{label}")
        for f, r, s, h, v in row:
            if f == 0.0 or (smallest and f <= smallest * 10) or r < 100:
                print(f"   field_err {f:8.1e}  truth/band {r:10.3g}  "
                      f"signal/band {s:9.3g}"
                      f"{'   <- holds' if h and not v else ''}"
                      f"{'  VACUOUS' if v else ''}")
        print(f"   -> requires {smallest:g}" if smallest else "   -> never holds")

    vals = [v for v in need.values() if v]
    print("\n=== VERDICT")
    if len(vals) < 2:
        print("  INCONCLUSIVE: fewer than two relations ever hold")
        res["verdict"] = "inconclusive"
    else:
        tight, loose = min(vals), max(vals)
        floor_like = tight <= 0.1 * KNOWN_SCHEME_ERROR
        spread = loose / tight if tight > 0 else float("inf")
        print(f"  tightest requires {tight:g}, loosest {loose:g}, spread {spread:.0f}x")
        print(f"  independently measured scheme error {KNOWN_SCHEME_ERROR:g}")
        # H has TWO clauses and they must be scored separately: the tightest
        # must sit near the FLOOR (not merely below the scheme error), and the
        # SPREAD must estimate the scheme error. Scoring only the first is how a
        # half-failing hypothesis gets recorded as a pass.
        near_floor = tight <= 100 * res["sigma_rep"]
        spread_predicts = 0.1 <= (loose / KNOWN_SCHEME_ERROR) <= 10
        res["clause_tightest_near_floor"] = bool(near_floor)
        res["clause_spread_estimates_scheme_error"] = bool(spread_predicts)
        print(f"  clause 1 -- tightest near the floor ({res['sigma_rep']:.1e})? "
              f"{near_floor} (it is {tight/res['sigma_rep']:.0f}x sigma_rep)")
        print(f"  clause 2 -- spread estimates the scheme error? "
              f"{spread_predicts} (loosest {loose:g} vs known "
              f"{KNOWN_SCHEME_ERROR:g}, off by {KNOWN_SCHEME_ERROR/loose:.0f}x)")
        if near_floor and spread_predicts:
            res["verdict"] = "H holds"
            print("  H HOLDS")
        elif floor_like:
            res["verdict"] = "H fails as stated; the scan is still useful"
            print("  H FAILS AS STATED. The tightest relation is far below the "
                  "scheme error but is NOT the floor, and the spread does not "
                  "estimate the scheme error. What the run DOES establish is "
                  "below.")
        else:
            print("  H FAILS: the TIGHTEST relation already requires a "
                  "declaration at the scheme-error scale, so no registered "
                  "relation sees the floor -- route 2 does not measure what it "
                  "was meant to on this family")
            res["verdict"] = "H fails"
        res.update({"tightest": tight, "loosest": loose, "spread": spread})
    OUT.write_text(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
