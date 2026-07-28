"""Exoplanet Archive C2 (docs/CASE_STUDY_EXOPLANET_C2.md — predictions
frozen): stellar-property conditioning, conjecture track only.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.adapter import C1_ADQL, fetch  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402

OUT = Path("experiments/results/exoplanet_c2.json")


def main():
    out = {}
    t, snap = fetch(C1_ADQL, "c1")
    out["snapshot"] = str(snap)
    R = np.asarray(t["pl_rade"], float)
    M = np.asarray(t["pl_bmasse"], float)
    Me1 = np.asarray(t["pl_bmasseerr1"], float)
    Me2 = np.asarray(t["pl_bmasseerr2"], float)
    Re1 = np.asarray(t["pl_radeerr1"], float)
    Re2 = np.asarray(t["pl_radeerr2"], float)
    prov = np.array([str(v) for v in t["pl_bmassprov"]])
    P = np.asarray(t["pl_orbper"], float)
    teff = np.asarray(t["st_teff"], float)
    met = np.asarray(t["st_met"], float)
    mstar = np.asarray(t["st_mass"], float)
    disc = np.array([str(v) for v in t["discoverymethod"]])
    merr = np.maximum(np.abs(Me1), np.abs(Me2)) / np.abs(M)
    rerr = np.maximum(np.abs(Re1), np.abs(Re2)) / np.abs(R)
    q = (prov == "Mass") & np.isfinite(M) & np.isfinite(R) & (M > 0) \
        & (R > 0) & (merr < 0.25) & (rerr < 0.08)

    # P1: giant vs small host metallicity
    gi = q & (M > 100) & np.isfinite(met)
    sm = q & (R < 4) & np.isfinite(met)
    d = float(np.mean(met[gi]) - np.mean(met[sm]))
    rng = np.random.default_rng(0)
    boots = [float(np.mean(rng.choice(met[gi], gi.sum()))
                   - np.mean(rng.choice(met[sm], sm.sum())))
             for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    out["P1_giant_metallicity"] = {
        "n_giant": int(gi.sum()), "n_small": int(sm.sum()),
        "delta_feh_dex": d, "ci95": [float(lo), float(hi)]}

    # P2: valley location vs stellar mass
    fv = (disc == "Transit") & np.isfinite(teff) & (teff >= 4700) \
        & (teff <= 6500) & (P < 100) & np.isfinite(R) & (R > 0.8) \
        & (R < 6.0) & (rerr < 0.08) & np.isfinite(mstar) & (mstar > 0)
    lms, lrf = np.log10(mstar[fv]), np.log10(R[fv])
    edges = np.quantile(lms, np.linspace(0, 1, 5))
    grid = np.linspace(np.log10(1.2), np.log10(2.6), 80)
    lc, vs = [], []
    for lo_, hi_ in zip(edges, edges[1:]):
        m = (lms >= lo_) & (lms < hi_)
        if m.sum() < 100:
            continue
        kde = np.array([np.sum(np.exp(-0.5 * ((lrf[m] - g0) / 0.04) ** 2))
                        for g0 in grid])
        w = (grid >= np.log10(1.4)) & (grid <= np.log10(2.4))
        lc.append((lo_ + hi_) / 2)
        vs.append(grid[np.flatnonzero(w)[int(np.argmin(kde[w]))]])
    lc, vs = np.array(lc), np.array(vs)
    slope = None
    if len(lc) >= 3:
        A = np.column_stack([lc, np.ones_like(lc)])
        (slope, _c), *_ = np.linalg.lstsq(A, vs, rcond=None)
    out["P2_valley_vs_mstar"] = {
        "n_bins": len(lc), "dlogRv_dlogMstar": float(slope),
        "valley_re_per_bin": [float(10 ** v) for v in vs],
        "edge_pinned": bool(np.any(np.isclose(vs, grid[0]))
                            or np.any(np.isclose(vs, grid[-1])))}

    # P3: volatile slope by metallicity half
    vw = q & (R >= 1.8) & (R < 4.0) & np.isfinite(met)
    zmed = float(np.median(met[vw]))
    slopes = {}
    for name, m in (("low_z", vw & (met < zmed)), ("high_z", vw & (met >= zmed))):
        lR, lM = np.log10(R[m]), np.log10(M[m])
        A = np.column_stack([lR, np.ones_like(lR)])
        (s, _c), *_ = np.linalg.lstsq(A, lM, rcond=None)
        r = recover(lR.reshape(-1, 1).tolist(), lM.tolist(), sigma=0.0)
        slopes[name] = {"n": int(m.sum()), "slope": float(s),
                        "recover_certified": r.get("certified"),
                        "recover_abstain": r.get("abstain")}
    out["P3_volatile_slope_by_met"] = slopes
    out["P3_abs_delta_slope"] = abs(slopes["high_z"]["slope"]
                                    - slopes["low_z"]["slope"])
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
