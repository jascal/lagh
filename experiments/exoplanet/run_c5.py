"""Exoplanet Archive C5 (docs/CASE_STUDY_EXOPLANET_C5.md — predictions
frozen): the open sweep, the definitional triples, and the stress test.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.adapter import C0_ADQL, C5_ADQL, fetch  # noqa: E402
from experiments.exoplanet.run_c0 import half_step_rel  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402

OUT = Path("experiments/results/exoplanet_c5.json")

PHYS = ["pl_rade", "pl_bmasse", "pl_orbper", "pl_orbsmax", "pl_dens",
        "pl_insol", "pl_eqt", "pl_orbeccen", "pl_ratror", "pl_ratdor",
        "pl_trandep", "st_teff", "st_rad", "st_mass", "st_logg", "sy_dist"]


def scout_rel(x, y):
    best = np.inf
    ysc = float(np.max(np.abs(y))) + 1e-300
    bases = [np.column_stack([x, np.ones_like(x)]),
             np.column_stack([x, x**2, np.ones_like(x)])]
    with np.errstate(all="ignore"):
        if np.all(x > 0):
            bases.append(np.column_stack([np.log(x), np.ones_like(x)]))
            bases.append(np.column_stack([1.0 / x, np.ones_like(x)]))
        if np.all(x > 0) and np.all(y > 0):
            A = np.column_stack([np.log(x), np.ones_like(x)])
            c, *_ = np.linalg.lstsq(A, np.log(y), rcond=None)
            best = min(best, float(np.max(np.abs(np.exp(A @ c) - y))) / ysc)
    for A in bases:
        if not np.all(np.isfinite(A)):
            continue
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        best = min(best, float(np.max(np.abs(A @ c - y))) / ysc)
    return best


def main():
    out = {}
    t, snap = fetch(C5_ADQL, "c5")
    out["snapshot"] = str(snap)
    out["n_total"] = len(t)
    D = {c: np.asarray(t[c], float) for c in PHYS}

    census = {"certified": [], "conjectured": [], "abstained": []}
    details = []
    for xn in PHYS:
        for yn in PHYS:
            if xn == yn:
                continue
            x, y = D[xn], D[yn]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 50:
                census["abstained"].append(f"{xn}->{yn}")
                continue
            rel = scout_rel(x[ok], y[ok])
            if rel > 1e-6:
                census["conjectured" if np.isfinite(rel) else
                       "abstained"].append(f"{xn}->{yn}")
                continue
            sig = half_step_rel(np.abs(y[ok])) + half_step_rel(np.abs(x[ok]))
            r = recover(x[ok].reshape(-1, 1).tolist(), y[ok].tolist(),
                        sigma=float(sig))
            lab = "certified" if r.get("certified") else "abstained"
            census[lab].append(f"{xn}->{yn}")
            details.append({"pair": f"{xn}->{yn}",
                            "certified": r.get("certified"),
                            "law": (r.get("law") or "")[:120],
                            "alpha_log10": r.get("alpha_log10"),
                            "abstain": r.get("abstain")})

    # P1: registered triples on Calculated strata (flag where present)
    def calc(col):
        return np.array(["Calculated" in str(v) for v in t[col + "_reflink"]])

    for x0n, x1n, yn, tag in (
            ("pl_rade", "st_rad", "pl_ratror", "P1_ratror"),
            ("pl_orbsmax", "st_rad", "pl_ratdor", "P1_ratdor")):
        x0, x1, y = D[x0n], D[x1n], D[yn]
        m = calc(yn) & np.isfinite(x0) & np.isfinite(x1) & np.isfinite(y) \
            & (x0 > 0) & (x1 > 0) & (y > 0)
        sig = (half_step_rel(y[m]) + half_step_rel(x0[m])
               + half_step_rel(x1[m]))
        r = recover(np.column_stack([x0[m], x1[m]]).tolist(), y[m].tolist(),
                    sigma=float(sig))
        out[tag] = {"n": int(m.sum()), "sigma_rep": float(sig),
                    "certified": r.get("certified"),
                    "law": (r.get("law") or "")[:120],
                    "alpha_log10": r.get("alpha_log10"),
                    "abstain": r.get("abstain")}
    x, y = D["pl_ratror"], D["pl_trandep"]
    m = calc("pl_trandep") & np.isfinite(x) & np.isfinite(y) & (x > 0) \
        & (y > 0)
    sig = half_step_rel(y[m]) + 2 * half_step_rel(x[m])
    r = recover(x[m].reshape(-1, 1).tolist(), y[m].tolist(), sigma=float(sig))
    out["P1_trandep"] = {"n": int(m.sum()), "sigma_rep": float(sig),
                         "certified": r.get("certified"),
                         "law": (r.get("law") or "")[:120],
                         "alpha_log10": r.get("alpha_log10"),
                         "abstain": r.get("abstain")}

    out["P2_census"] = {k: len(v) for k, v in census.items()}
    out["P2_certified_list"] = census["certified"]
    out["P2_details"] = details

    # P3: stress on the C0 density identity
    t0, _ = fetch(C0_ADQL, "c0")
    M = np.asarray(t0["pl_bmasse"], float)
    R = np.asarray(t0["pl_rade"], float)
    rho = np.asarray(t0["pl_dens"], float)
    ok = np.isfinite(M) & np.isfinite(R) & np.isfinite(rho) & (M > 0) \
        & (R > 0) & (rho > 0)
    sig0 = half_step_rel(rho[ok]) + half_step_rel(M[ok]) \
        + 3 * half_step_rel(R[ok])
    r_mix = recover(np.column_stack([M[ok], R[ok]]).tolist(),
                    rho[ok].tolist(), sigma=float(sig0))
    out["P3a_full_mixture"] = {"n": int(ok.sum()),
                               "certified": r_mix.get("certified"),
                               "abstain": r_mix.get("abstain")}
    calc0 = np.array(["Calculated" in str(v) for v in t0["pl_dens_reflink"]])
    mc = calc0 & ok
    r_loose = recover(np.column_stack([M[mc], R[mc]]).tolist(),
                      rho[mc].tolist(), sigma=float(10 * sig0))
    law = r_loose.get("law") or ""
    out["P3b_loose_10x"] = {"certified": r_loose.get("certified"),
                            "law": law[:120],
                            "n_terms": law.count("+") + law.count("-")
                            if r_loose.get("certified") else None,
                            "abstain": r_loose.get("abstain")}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
