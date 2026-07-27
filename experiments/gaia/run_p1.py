"""Gaia Phase 1, C1+C2 (docs/CASE_STUDY_GAIA_P1.md — predictions frozen).

C1: Stefan-Boltzmann as a pipeline-definitional identity (P1), the loose-floor
    discipline check (P2), mass-luminosity as conjecture-only (P3).
C2: no deterministic kinematic law (P4), Stromberg asymmetric drift as a
    labeled conjecture (P5).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.gaia.adapter import C1_ADQL, C2_ADQL, fetch  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402
from lagh.submit import submission  # noqa: E402

OUT = Path("experiments/results/gaia_p1.json")


def quantization(col: np.ndarray) -> float:
    """Median relative decimal-rounding half-step of a positive column."""
    rels = []
    for v in col:
        s = f"{v!r}"
        if "e" in s or "E" in s:
            continue
        dec = len(s.split(".")[1]) if "." in s else 0
        rels.append(0.5 * 10.0 ** (-dec) / abs(v))
    return float(np.median(rels)) if rels else 1e-7


def c1(out):
    t, snap = fetch(C1_ADQL, "c1")
    out["c1_snapshot"] = str(snap)
    T = np.asarray(t["teff_gspphot"], float)
    L = np.asarray(t["lum_flame"], float)
    R = np.asarray(t["radius_flame"], float)
    M = np.asarray(t["mass_flame"], float)
    ok = np.isfinite(T) & np.isfinite(L) & np.isfinite(R) & (T > 0) & (L > 0) \
        & (R > 0)
    T, L, R, M = T[ok], L[ok], R[ok], M[ok]
    out["c1_n"] = int(ok.sum())
    # registered floor procedure: propagated log-rounding, 2x margin
    r_L, r_R, r_T = quantization(L), quantization(R), quantization(T)
    floor = 2.0 * (r_L + 2 * r_R + 4 * r_T) / np.log(10)
    out["c1_floor_measured"] = {"r_L": r_L, "r_R": r_R, "r_T": r_T,
                                "floor_abs": floor}
    X = np.column_stack([np.log10(R), np.log10(T)]).tolist()
    y = np.log10(L).tolist()
    r = recover(X, y, sigma=0.0, floor_abs=floor)
    out["P1_stefan_boltzmann"] = {
        "certified": r.get("certified"), "law": r.get("law"),
        "alpha_log10": r.get("alpha_log10"), "abstain": r.get("abstain")}
    # post-hoc constant decode (not part of the certificate claim)
    out["P1_const_expected_m4log10Tsun"] = -4 * float(np.log10(5772.0))
    # P2: deliberately loose floor -- must NOT certify a multi-term approximant
    r2 = recover(X, y, sigma=0.0, floor_abs=1e-2)
    law2 = r2.get("law") or ""
    out["P2_loose_floor"] = {
        "certified": r2.get("certified"), "law": law2,
        "abstain": r2.get("abstain"),
        "n_terms_if_certified": law2.count("+") + law2.count("-")
        if r2.get("certified") else None}
    # P3: mass-luminosity on the MS window -- conjecture only
    w = (M >= 0.5) & (M <= 2.0)
    sub = submission(np.log10(M[w]).reshape(-1, 1), np.log10(L[w]), sigma=0.0)
    slope = None
    if sub.get("expr"):
        import sympy as sp
        try:
            e = sp.sympify(sub["expr"])
            slope = float(e.diff(sp.Symbol("x_0")))
        except Exception:                                     # noqa: BLE001
            pass
    out["P3_mass_luminosity"] = {"track": sub["track"], "tag": sub["tag"],
                                 "expr": (sub.get("expr") or "")[:90],
                                 "slope": slope, "n": int(w.sum())}


def c2(out):
    t, snap = fetch(C2_ADQL, "c2")
    out["c2_snapshot"] = str(snap)
    from astropy import units as u
    from astropy.coordinates import ICRS, Galactic
    ic = ICRS(ra=np.asarray(t["ra"], float) * u.deg,
              dec=np.asarray(t["dec"], float) * u.deg,
              distance=(1000.0 / np.asarray(t["parallax"], float)) * u.pc,
              pm_ra_cosdec=np.asarray(t["pmra"], float) * u.mas / u.yr,
              pm_dec=np.asarray(t["pmdec"], float) * u.mas / u.yr,
              radial_velocity=np.asarray(t["radial_velocity"], float)
              * u.km / u.s)
    g = ic.transform_to(Galactic())
    g.representation_type = "cartesian"
    g.differential_type = "cartesian"
    U = g.velocity.d_x.to_value(u.km / u.s)
    V = g.velocity.d_y.to_value(u.km / u.s)
    W = g.velocity.d_z.to_value(u.km / u.s)
    bprp = np.asarray(t["bp_rp"], float)
    ok = np.isfinite(U) & np.isfinite(V) & np.isfinite(W) & np.isfinite(bprp)
    U, V, W, bprp = U[ok], V[ok], W[ok], bprp[ok]
    out["c2_n"] = int(ok.sum())
    # P4: no deterministic law among velocity components
    for name, x, y in (("U_to_V", U, V), ("U_to_W", U, W)):
        r = recover(x.reshape(-1, 1).tolist(), y.tolist(), sigma=0.0)
        out[f"P4_{name}"] = {"certified": r.get("certified"),
                             "abstain": r.get("abstain")}
    # P5: Stromberg -- 8 color bins, mean V vs sigma_U^2
    qs = np.quantile(bprp, np.linspace(0, 1, 9))
    mv, s2 = [], []
    for lo, hi in zip(qs, qs[1:]):
        m = (bprp >= lo) & (bprp < hi if hi < qs[-1] else bprp <= hi)
        if m.sum() < 30:
            continue
        mv.append(float(np.mean(V[m])))
        s2.append(float(np.var(U[m])))
    s2, mv = np.array(s2), np.array(mv)
    A = np.column_stack([s2, np.ones_like(s2)])
    (slope, icpt), *_ = np.linalg.lstsq(A, mv, rcond=None)
    corr = float(np.corrcoef(s2, mv)[0, 1])
    r5 = recover(s2.reshape(-1, 1).tolist(), mv.tolist(), sigma=0.0)
    out["P5_stromberg"] = {
        "n_bins": len(s2), "slope": float(slope), "corr": corr,
        "k_kms": float(-1.0 / slope) if slope < 0 else None,
        "binned_recover_certified": r5.get("certified"),
        "binned_recover_abstain": r5.get("abstain")}


def main():
    out = {}
    c1(out)
    c2(out)
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
