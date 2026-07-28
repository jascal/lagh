"""Exoplanet Archive Phase 2, C3+C4 (docs/CASE_STUDY_EXOPLANET_PH2.md —
predictions frozen): density/irradiation + multi-planet architecture.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.adapter import PH2_ADQL, fetch  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402

OUT = Path("experiments/results/exoplanet_ph2.json")


def ols(x, y):
    A = np.column_stack([x, np.ones_like(x)])
    (s, c), *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(s)


def main():
    out = {}
    t, snap = fetch(PH2_ADQL, "ph2")
    out["snapshot"] = str(snap)
    out["n_total"] = len(t)
    R = np.asarray(t["pl_rade"], float)
    Re1 = np.asarray(t["pl_radeerr1"], float)
    Re2 = np.asarray(t["pl_radeerr2"], float)
    M = np.asarray(t["pl_bmasse"], float)
    rho = np.asarray(t["pl_dens"], float)
    S = np.asarray(t["pl_insol"], float)
    P = np.asarray(t["pl_orbper"], float)
    a = np.asarray(t["pl_orbsmax"], float)
    Ms = np.asarray(t["st_mass"], float)
    host = np.array([str(v) for v in t["hostname"]])
    npl = np.asarray(t["sy_pnum"], float)
    rerr = np.maximum(np.abs(Re1), np.abs(Re2)) / np.abs(R)

    # ---------- C3 ----------
    gi = np.isfinite(R) & (R > 8) & np.isfinite(S) & (S > 0) & (rerr < 0.08)
    hot, cool = gi & (S > 200), gi & (S < 200)
    d = float(np.mean(R[hot]) - np.mean(R[cool]))
    rng = np.random.default_rng(0)
    boots = [float(np.mean(rng.choice(R[hot], hot.sum()))
                   - np.mean(rng.choice(R[cool], cool.sum())))
             for _ in range(2000)]
    out["P1_inflation_offset"] = {
        "n_hot": int(hot.sum()), "n_cool": int(cool.sum()),
        "delta_re": d,
        "ci95": [float(v) for v in np.percentile(boots, [2.5, 97.5])]}
    coolf = gi & (S < 100)
    out["P2_threshold"] = {
        "cool_slope_re_per_dex": ols(np.log10(S[coolf]), R[coolf]),
        "n_cool": int(coolf.sum()),
        "hot_slope_re_per_dex": ols(np.log10(S[hot]), R[hot]),
        "n_hot": int(hot.sum())}
    sm = np.isfinite(R) & (R >= 1) & (R < 4) & np.isfinite(rho) & (rho > 0)
    s3 = ols(np.log10(R[sm]), np.log10(rho[sm]))
    r3 = recover(np.log10(R[sm]).reshape(-1, 1).tolist(),
                 np.log10(rho[sm]).tolist(), sigma=0.0)
    out["P3_density_radius"] = {"n": int(sm.sum()), "slope": s3,
                                "recover_certified": r3.get("certified"),
                                "recover_abstain": r3.get("abstain")}

    # ---------- C4 ----------
    multi = np.isfinite(P) & (npl >= 2)
    ratios, hills = [], []
    for h in np.unique(host[multi]):
        m = multi & (host == h)
        if m.sum() < 2:
            continue
        order = np.argsort(P[m])
        Ph = P[m][order]
        ah = a[m][order]
        Mh = M[m][order]
        Msh = Ms[m][order]
        for i in range(len(Ph) - 1):
            ratios.append(Ph[i + 1] / Ph[i])
            if (np.isfinite(Mh[i]) and np.isfinite(Mh[i + 1])
                    and np.isfinite(ah[i]) and np.isfinite(ah[i + 1])
                    and np.isfinite(Msh[i]) and Msh[i] > 0
                    and ah[i + 1] > ah[i] > 0):
                rh = (((Mh[i] + Mh[i + 1]) / (3 * Msh[i] * 332946.0))
                      ** (1.0 / 3.0)) * (ah[i] + ah[i + 1]) / 2.0
                hills.append((ah[i + 1] - ah[i]) / rh)
    ratios = np.array(ratios)
    hills = np.array(hills)
    out["P5_period_ratios"] = {
        "n_pairs": len(ratios), "median": float(np.median(ratios)),
        "frac_below_1.3": float(np.mean(ratios < 1.3))}
    out["P6_resonance_asymmetry"] = {
        "wide_2to1": int(np.sum((ratios >= 2.00) & (ratios <= 2.10))),
        "narrow_2to1": int(np.sum((ratios >= 1.90) & (ratios < 2.00))),
        "wide_3to2": int(np.sum((ratios >= 1.50) & (ratios <= 1.57))),
        "narrow_3to2": int(np.sum((ratios >= 1.43) & (ratios < 1.50)))}
    out["P7_hill_spacing"] = {
        "n_pairs": len(hills), "median_delta": float(np.median(hills)),
        "frac_below_8": float(np.mean(hills < 8))}
    idx = np.arange(min(len(ratios), 400))
    r8 = recover(idx.reshape(-1, 1).tolist(),
                 ratios[:len(idx)].tolist(), sigma=0.0)
    out["P8_ratio_sequence_recover"] = {"certified": r8.get("certified"),
                                        "abstain": r8.get("abstain")}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
