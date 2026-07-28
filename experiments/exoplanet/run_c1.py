"""Exoplanet Archive C1 (docs/CASE_STUDY_EXOPLANET_C1.md — predictions
frozen): mass-radius sequences and the radius valley, conjecture track only.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.adapter import C1_ADQL, fetch  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402

OUT = Path("experiments/results/exoplanet_c1.json")


def slope_window(R, M, lo, hi, tag, out, sigma_scout=0.02):
    w = (R >= lo) & (R < hi)
    lR, lM = np.log10(R[w]), np.log10(M[w])
    A = np.column_stack([lR, np.ones_like(lR)])
    (s, c), *_ = np.linalg.lstsq(A, lM, rcond=None)
    resid = lM - A @ np.array([s, c])
    r = recover(lR.reshape(-1, 1).tolist(), lM.tolist(), sigma=0.0)
    out[tag] = {"n": int(w.sum()), "slope": float(s),
                "scatter_dex": float(np.std(resid)),
                "recover_certified": r.get("certified"),
                "recover_abstain": r.get("abstain")}


def main():
    out = {}
    t, snap = fetch(C1_ADQL, "c1")
    out["snapshot"] = str(snap)
    out["n_total"] = len(t)
    R = np.asarray(t["pl_rade"], float)
    M = np.asarray(t["pl_bmasse"], float)
    Me1 = np.asarray(t["pl_bmasseerr1"], float)
    Me2 = np.asarray(t["pl_bmasseerr2"], float)
    Re1 = np.asarray(t["pl_radeerr1"], float)
    Re2 = np.asarray(t["pl_radeerr2"], float)
    prov = np.array([str(v) for v in t["pl_bmassprov"]])
    P = np.asarray(t["pl_orbper"], float)
    teff = np.asarray(t["st_teff"], float)
    disc = np.array([str(v) for v in t["discoverymethod"]])

    merr = np.maximum(np.abs(Me1), np.abs(Me2)) / np.abs(M)
    rerr = np.maximum(np.abs(Re1), np.abs(Re2)) / np.abs(R)
    mr = (prov == "Mass") & np.isfinite(M) & np.isfinite(R) & (M > 0) \
        & (R > 0) & (merr < 0.25) & (rerr < 0.08)
    out["n_mass_radius_window"] = int(mr.sum())

    # P1-P3: regime slopes
    slope_window(R[mr], M[mr], 0.0, 1.5, "P1_rocky", out)
    slope_window(R[mr], M[mr], 1.8, 4.0, "P2_volatile", out)
    g = mr & (M > 100)
    lM, lR = np.log10(M[g]), np.log10(R[g])
    A = np.column_stack([lM, np.ones_like(lM)])
    (s, c), *_ = np.linalg.lstsq(A, lR, rcond=None)
    rg = recover(lM.reshape(-1, 1).tolist(), lR.tolist(), sigma=0.0)
    out["P3_giant_plateau"] = {"n": int(g.sum()), "slope_dlogR_dlogM": float(s),
                               "recover_certified": rg.get("certified"),
                               "recover_abstain": rg.get("abstain")}

    # P4: radius valley on the Fulton-like sample
    fv = (disc == "Transit") & np.isfinite(teff) & (teff >= 4700) \
        & (teff <= 6500) & (P < 100) & np.isfinite(R) & (R > 0) \
        & (rerr < 0.08)
    out["P4_n_fulton_sample"] = int(fv.sum())
    lr = np.log10(R[fv & (R > 0.8) & (R < 6.0)])
    bins = np.arange(np.log10(0.8), np.log10(6.0), 0.05)
    h, edges = np.histogram(lr, bins=bins)
    centers = 10 ** ((edges[:-1] + edges[1:]) / 2)
    vwin = (centers >= 1.5) & (centers <= 2.2)
    vi = np.flatnonzero(vwin)[np.argmin(h[vwin])]
    vloc = float(centers[vi])
    lmax = int(h[:vi].max()) if vi > 0 else 0
    rmax = int(h[vi + 1:].max()) if vi + 1 < len(h) else 0
    contrast = min(lmax, rmax) / max(int(h[vi]), 1)
    out["P4_valley"] = {"location_re": vloc, "valley_count": int(h[vi]),
                        "flank_maxima": [lmax, rmax],
                        "contrast": float(contrast)}
    # (b) period dependence of the valley center
    pf, vloc_b = P[fv & (R > 0.8) & (R < 6.0)], []
    lp_edges = np.quantile(np.log10(pf), np.linspace(0, 1, 6))
    lp_c = []
    lrf = np.log10(R[fv & (R > 0.8) & (R < 6.0)])
    for lo, hi in zip(lp_edges, lp_edges[1:]):
        m = (np.log10(pf) >= lo) & (np.log10(pf) < hi)
        if m.sum() < 80:
            continue
        # smoothed distribution in the valley neighborhood
        grid = np.linspace(np.log10(1.4), np.log10(2.4), 60)
        kde = np.array([np.sum(np.exp(-0.5 * ((lrf[m] - g0) / 0.04) ** 2))
                        for g0 in grid])
        lp_c.append((lo + hi) / 2)
        vloc_b.append(grid[int(np.argmin(kde))])
    lp_c, vloc_b = np.array(lp_c), np.array(vloc_b)
    if len(lp_c) >= 3:
        A = np.column_stack([lp_c, np.ones_like(lp_c)])
        (sv, cv), *_ = np.linalg.lstsq(A, vloc_b, rcond=None)
        out["P4_valley_period_slope"] = {"n_bins": len(lp_c),
                                         "dlogRv_dlogP": float(sv),
                                         "valley_re_per_bin":
                                         [float(10 ** v) for v in vloc_b]}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
