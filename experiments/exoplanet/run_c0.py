"""Exoplanet Archive C0 (docs/CASE_STUDY_EXOPLANET_C0.md — predictions
frozen): the Archive's computed columns certify; the literature stratum and
real physics stay honest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.adapter import C0_ADQL, fetch  # noqa: E402
from experiments.gaia.run_p2 import half_ulp_rel, registered_floor  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402
from lagh.submit import submission  # noqa: E402

OUT = Path("experiments/results/exoplanet_c0.json")


def calc_mask(t, col):
    ref = [str(v) for v in t[col + "_reflink"]]
    return np.array(["Calculated" in v for v in ref])


def rec2(X0, X1, y, tag, out):
    ok = np.isfinite(X0) & np.isfinite(X1) & np.isfinite(y)
    X0, X1, y = X0[ok], X1[ok], y[ok]
    terms = [np.ones_like(y), X0, X1]
    floor = registered_floor(y, terms)
    r = recover(np.column_stack([X0, X1]).tolist(), y.tolist(), sigma=0.0,
                floor_abs=floor)
    out[tag] = {"n": int(ok.sum()), "floor_abs": floor,
                "certified": r.get("certified"),
                "law": (r.get("law") or "")[:140],
                "alpha_log10": r.get("alpha_log10"),
                "abstain": r.get("abstain")}
    return r


def main():
    out = {}
    t, snap = fetch(C0_ADQL, "c0")
    out["snapshot"] = str(snap)
    out["n_total"] = len(t)
    M = np.asarray(t["pl_bmasse"], float)
    R = np.asarray(t["pl_rade"], float)
    rho = np.asarray(t["pl_dens"], float)
    P = np.asarray(t["pl_orbper"], float)
    a = np.asarray(t["pl_orbsmax"], float)
    Ms = np.asarray(t["st_mass"], float)
    lum = np.asarray(t["st_lum"], float)
    S = np.asarray(t["pl_insol"], float)

    # P1: density identity on the Calculated subsample vs literature
    mc = calc_mask(t, "pl_dens") & (M > 0) & (R > 0) & (rho > 0)
    out["P1_n_calculated"] = int(mc.sum())
    rec2(np.log10(M[mc]), np.log10(R[mc]), np.log10(rho[mc]),
         "P1_density_calculated", out)
    ml = ~calc_mask(t, "pl_dens") & (M > 0) & (R > 0) & (rho > 0)
    rec2(np.log10(M[ml]), np.log10(R[ml]), np.log10(rho[ml]),
         "P1_density_literature", out)
    out["P1_const_expected"] = float(np.log10(5.51295))

    # P2: Kepler III on the Calculated pl_orbsmax subsample
    k = calc_mask(t, "pl_orbsmax") & (P > 0) & (a > 0) & (Ms > 0) \
        & np.isfinite(Ms)
    out["P2_n_calculated"] = int(k.sum())
    rec2(np.log10(P[k]), np.log10(Ms[k]), np.log10(a[k]),
         "P2_kepler_calculated", out)
    out["P2_const_expected"] = float(-(2.0 / 3.0) * np.log10(365.25))

    # P3: insolation identity on the Calculated pl_insol subsample
    i = calc_mask(t, "pl_insol") & (S > 0) & (a > 0) & np.isfinite(lum)
    out["P3_n_calculated"] = int(i.sum())
    rec2(lum[i], np.log10(a[i]), np.log10(S[i]), "P3_insol_calculated", out)

    # P4: population mass-radius must NOT certify
    ok = (M > 0) & (R > 0)
    sub = submission(np.log10(R[ok]).reshape(-1, 1), np.log10(M[ok]),
                     sigma=0.0)
    out["P4_mass_radius"] = {"track": sub["track"], "tag": sub["tag"],
                             "expr": (sub.get("expr") or "")[:90]}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
