"""Gaia Phase 2, C3+C4 (docs/CASE_STUDY_GAIA_P2.md — predictions frozen).

C3: no cross-catalog Kepler law (P1), photocentric mass scale conjecture (P2).
C4: RRab pipeline metallicity as a certification target (P3), Leavitt law as
    conjecture (P4).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.gaia.adapter import (C3_ADQL, C4CEP_ADQL, C4RR_ADQL,  # noqa: E402
                                      fetch)
from lagh.mcp.core import recover  # noqa: E402
from lagh.submit import submission  # noqa: E402

OUT = Path("experiments/results/gaia_p2.json")


def half_ulp_rel(col: np.ndarray) -> float:
    """Storage-precision half-ulp (v2 procedure): float32 round-trip
    detection first (the P1 lesson), decimal-repr half-step fallback."""
    col = col[np.isfinite(col) & (col != 0)]
    if len(col) and np.all(col.astype(np.float32).astype(np.float64) == col):
        return 2.0 ** -25
    rels = []
    for v in col:
        s = repr(float(v))
        dec = len(s.split(".")[1]) if "." in s and "e" not in s else 17
        rels.append(0.5 * 10.0 ** (-dec) / abs(v))
    return float(np.median(rels)) if rels else 2.0 ** -25


def registered_floor(y, terms):
    """floor_abs = 2*(q_y + sum_j |beta_j| * q_term_j), beta from an OLS scout
    over the registered basis (logged calibration, before discovery)."""
    A = np.column_stack([t for t in terms])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    q_y = half_ulp_rel(y) * float(np.median(np.abs(y)))
    q_t = [half_ulp_rel(t) * float(np.median(np.abs(t))) for t in terms]
    return 2.0 * (q_y + float(np.sum(np.abs(beta) * np.array(q_t))))


def c3(out):
    t, snap = fetch(C3_ADQL, "c3")
    out["c3_snapshot"] = str(snap)
    A = np.asarray(t["a_thiele_innes"], float)
    B = np.asarray(t["b_thiele_innes"], float)
    F = np.asarray(t["f_thiele_innes"], float)
    G = np.asarray(t["g_thiele_innes"], float)
    plx = np.asarray(t["parallax"], float)
    P = np.asarray(t["period"], float) / 365.25
    u = (A**2 + B**2 + F**2 + G**2) / 2.0
    a0 = np.sqrt(u + np.sqrt(np.maximum(u**2 - (A * G - B * F) ** 2, 0.0)))
    a_au = a0 / plx
    ok = np.isfinite(a_au) & np.isfinite(P) & (a_au > 0) & (P > 0)
    a_au, P = a_au[ok], P[ok]
    out["c3_n"] = int(ok.sum())
    r = recover(np.log10(a_au).reshape(-1, 1).tolist(),
                np.log10(P).tolist(), sigma=0.0)
    out["P1_no_single_kepler"] = {"certified": r.get("certified"),
                                  "abstain": r.get("abstain")}
    m = a_au**3 / P**2
    out["P2_mass_proxy"] = {"median_msun": float(np.median(m)),
                            "q10": float(np.quantile(m, 0.1)),
                            "q90": float(np.quantile(m, 0.9))}


def c4(out):
    t, snap = fetch(C4RR_ADQL, "c4rr")
    out["c4rr_snapshot"] = str(snap)
    P = np.asarray(t["pf"], float)
    phi = np.asarray(t["phi31_g"], float)
    feh = np.asarray(t["metallicity"], float)
    ok = np.isfinite(P) & np.isfinite(phi) & np.isfinite(feh)
    P, phi, feh = P[ok], phi[ok], feh[ok]
    out["c4rr_n"] = int(ok.sum())
    # registered quadratic scout basis for the floor calibration
    terms = [np.ones_like(P), P, phi, P**2, P * phi, phi**2]
    floor = registered_floor(feh, terms)
    out["P3_floor_abs"] = floor
    X = np.column_stack([P, phi]).tolist()
    r = recover(X, feh.tolist(), sigma=0.0, floor_abs=floor)
    out["P3_rrab_metallicity"] = {
        "certified": r.get("certified"), "law": r.get("law"),
        "alpha_log10": r.get("alpha_log10"), "abstain": r.get("abstain")}
    if not r.get("certified"):
        # second track (logged): recover's abstain is a REACH gap (dim-2
        # cross-term quadratic support; registered instrument issue). The
        # calibration scout declares the form; verify certifies it or not.
        A = np.column_stack(terms)
        beta, *_ = np.linalg.lstsq(A, feh, rcond=None)
        keep = [i for i, b in enumerate(beta)
                if abs(b) * float(np.median(np.abs(terms[i]))) > floor]
        Ak = A[:, keep]
        bk, *_ = np.linalg.lstsq(Ak, feh, rcond=None)
        snapped = bk.copy()
        for i in range(len(bk)):
            for d in range(2, 11):
                trial = snapped.copy()
                trial[i] = round(bk[i], d)
                if np.max(np.abs(feh - Ak @ trial)) <= floor:
                    snapped[i] = trial[i]
                    break
        names = ["1", "x_0", "x_1", "x_0**2", "x_0*x_1", "x_1**2"]
        form = " + ".join(f"({float(snapped[j])!r})*{names[k]}"
                          for j, k in enumerate(keep))
        from lagh.mcp.core import verify
        v = verify(X, feh.tolist(), form, sigma=0.0, floor_abs=floor)
        out["P3_declared_form_verify"] = {
            "form": form, "certified": v.get("certified"),
            "strength": v.get("strength"), "law": v.get("law"),
            "note": v.get("note")}
        res = feh - Ak @ snapped
        out["P3_residual_decode"] = {
            "scout_max_abs_res": float(np.max(np.abs(res))),
            "scout_med_abs_res": float(np.median(np.abs(res)))}
    # P4: Leavitt law
    tc, snapc = fetch(C4CEP_ADQL, "c4cep")
    out["c4cep_snapshot"] = str(snapc)
    pf = np.asarray(tc["pf"], float)
    g = np.asarray(tc["int_average_g"], float)
    plx = np.asarray(tc["parallax"], float)
    ok = np.isfinite(pf) & np.isfinite(g) & np.isfinite(plx) & (plx > 0) \
        & (pf > 0)
    MG = g[ok] + 5 * np.log10(plx[ok] / 100.0)
    lp = np.log10(pf[ok])
    out["c4cep_n"] = int(ok.sum())
    sub = submission(lp.reshape(-1, 1), MG, sigma=0.0)
    slope = None
    if sub.get("expr"):
        import sympy as sp
        try:
            slope = float(sp.sympify(sub["expr"]).diff(sp.Symbol("x_0")))
        except Exception:                                     # noqa: BLE001
            pass
    out["P4_leavitt"] = {"track": sub["track"], "tag": sub["tag"],
                         "expr": (sub.get("expr") or "")[:90],
                         "slope": slope}


def main():
    out = {}
    c3(out)
    c4(out)
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
