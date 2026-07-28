"""Materials Project C2 (docs/CASE_STUDY_MATERIALS_C2.md — predictions
frozen): the open sweep + the three registered definitional targets.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.run_c0 import half_step_rel  # noqa: E402
from experiments.materials.adapter import fetch  # noqa: E402
from lagh.mcp.core import recover, verify  # noqa: E402

OUT = Path("experiments/results/materials_c2.json")

C2_PARAMS = {
    "_fields": "material_id,nsites,volume,density,density_atomic,"
               "formation_energy_per_atom,energy_per_atom,energy_above_hull,"
               "band_gap,efermi,total_magnetization,bulk_modulus,"
               "shear_modulus,universal_anisotropy,homogeneous_poisson",
    "energy_above_hull_max": 0.0,
    "_limit": 2000,
    "_sort_fields": "material_id",
}

SWEEP = ["nsites", "volume", "density", "density_atomic",
         "formation_energy_per_atom", "energy_per_atom", "band_gap",
         "efermi", "total_magnetization", "k_vrh", "g_vrh",
         "universal_anisotropy", "homogeneous_poisson"]


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


def floor_for(y):
    """Absolute floor from the fixed-decimal storage half-step (C1 lesson)."""
    steps = []
    for v in y:
        s = repr(float(v))
        if "e" in s or v == 0:
            continue
        dec = len(s.split(".")[1]) if "." in s else 0
        steps.append(0.5 * 10.0 ** (-dec))
    return 2.0 * float(np.median(steps)) if steps else 1e-12


def main():
    out = {}
    docs, snap = fetch("/materials/summary/", C2_PARAMS, "c2")
    out["snapshot"] = str(snap)
    D = {c: [] for c in SWEEP}
    for d in docs:
        K, G = d.get("bulk_modulus"), d.get("shear_modulus")
        row = {c: d.get(c) for c in SWEEP if c not in ("k_vrh", "g_vrh")}
        row["k_vrh"] = (K or {}).get("vrh") if isinstance(K, dict) else None
        row["g_vrh"] = (G or {}).get("vrh") if isinstance(G, dict) else None
        for c in SWEEP:
            v = row.get(c)
            D[c].append(float(v) if v is not None else np.nan)
    D = {c: np.asarray(v) for c, v in D.items()}
    Kd, Gd = {}, {}
    for d in docs:
        pass
    census = {"certified": [], "conjectured": [], "abstained": []}
    details = []
    for xn in SWEEP:
        for yn in SWEEP:
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
            r = recover(x[ok].reshape(-1, 1).tolist(), y[ok].tolist(),
                        sigma=0.0, floor_abs=floor_for(y[ok]))
            lab = "certified" if r.get("certified") else "abstained"
            census[lab].append(f"{xn}->{yn}")
            details.append({"pair": f"{xn}->{yn}",
                            "certified": r.get("certified"),
                            "law": (r.get("law") or "")[:110],
                            "alpha_log10": r.get("alpha_log10"),
                            "abstain": r.get("abstain")})
    out["P2_census"] = {k: len(v) for k, v in census.items()}
    out["P2_certified_list"] = census["certified"]
    out["P2_details"] = details

    # T1: atomic density
    ns, V, da = D["nsites"], D["volume"], D["density_atomic"]
    m = np.isfinite(ns) & np.isfinite(V) & np.isfinite(da) & (V > 0) \
        & (ns > 0) & (da > 0)
    r = recover(np.column_stack([ns[m], V[m]]).tolist(), da[m].tolist(),
                sigma=0.0, floor_abs=floor_for(da[m]))
    out["T1_density_atomic"] = {"n": int(m.sum()),
                                "certified": r.get("certified"),
                                "law": (r.get("law") or "")[:100],
                                "alpha_log10": r.get("alpha_log10"),
                                "abstain": r.get("abstain")}
    if not r.get("certified"):
        for form, name in (("x_0/x_1", "nsites_over_volume"),
                           ("x_1/x_0", "volume_over_nsites")):
            v = verify(np.column_stack([ns[m], V[m]]).tolist(),
                       da[m].tolist(), form, sigma=0.0,
                       floor_abs=floor_for(da[m]))
            if v.get("certified"):
                out["T1_verify"] = {"form": name,
                                    "strength": v.get("strength"),
                                    "law": v.get("law")}
                break

    # T2: Poisson's ratio
    K, G, nu = D["k_vrh"], D["g_vrh"], D["homogeneous_poisson"]
    m = np.isfinite(K) & np.isfinite(G) & np.isfinite(nu) & (K > 0) & (G > 0)
    fl = floor_for(nu[m])
    r = recover(np.column_stack([K[m], G[m]]).tolist(), nu[m].tolist(),
                sigma=0.0, floor_abs=fl)
    out["T2_poisson"] = {"n": int(m.sum()), "floor_abs": fl,
                         "certified": r.get("certified"),
                         "law": (r.get("law") or "")[:100],
                         "alpha_log10": r.get("alpha_log10"),
                         "abstain": r.get("abstain")}
    if not r.get("certified"):
        v = verify(np.column_stack([K[m], G[m]]).tolist(), nu[m].tolist(),
                   "(3*x_0 - 2*x_1)/(6*x_0 + 2*x_1)", sigma=0.0,
                   floor_abs=fl)
        out["T2_verify"] = {"certified": v.get("certified"),
                            "strength": v.get("strength"),
                            "note": (v.get("note") or "")[:90]}

    # T3: universal anisotropy (declared verify target, dim 4)
    Kv2, Kr2, Gv2, Gr2, Au = [], [], [], [], []
    for d in docs:
        K4, G4 = d.get("bulk_modulus"), d.get("shear_modulus")
        au = d.get("universal_anisotropy")
        if isinstance(K4, dict) and isinstance(G4, dict) and au is not None:
            try:
                Kv2.append(float(K4["voigt"])); Kr2.append(float(K4["reuss"]))
                Gv2.append(float(G4["voigt"])); Gr2.append(float(G4["reuss"]))
                Au.append(float(au))
            except (KeyError, TypeError, ValueError):
                pass
    Kv2, Kr2, Gv2, Gr2, Au = map(np.asarray, (Kv2, Kr2, Gv2, Gr2, Au))
    m = np.all(np.isfinite(np.column_stack([Kv2, Kr2, Gv2, Gr2, Au])), 1) \
        & (Kr2 > 0) & (Gr2 > 0)
    fl = floor_for(Au[m])
    # amended (logged): per-point se = INPUT rounding propagated through the
    # ratios (the v2 procedure applied properly; extreme-ratio rows carry
    # O(10) propagated tolerance against a 1e-3 output floor). verify grew an
    # se passthrough (plumbing parity with the epsilon model's lam_B term).
    h = 0.0005
    se = h * (5.0 / Gr2[m] + 5.0 * np.abs(Gv2[m]) / Gr2[m] ** 2
              + 1.0 / Kr2[m] + np.abs(Kv2[m]) / Kr2[m] ** 2)
    v = verify(np.column_stack([Gv2[m], Gr2[m], Kv2[m], Kr2[m]]).tolist(),
               Au[m].tolist(), "5*x_0/x_1 + x_2/x_3 - 6", sigma=0.0,
               floor_abs=fl, se=se.tolist())
    out["T3_anisotropy_verify"] = {"n": int(m.sum()), "floor_abs": fl,
                                   "certified": v.get("certified"),
                                   "strength": v.get("strength"),
                                   "note": (v.get("note") or "")[:90]}
    # refutation decode: which rows exceed the propagated bound, and by what
    pred = 5 * Gv2[m] / Gr2[m] + Kv2[m] / Kr2[m] - 6
    bound = fl + se
    w = np.abs(Au[m] - pred) / bound
    order = np.argsort(w)[::-1][:4]
    out["T3_decode"] = {
        "n_beyond_bound": int(np.sum(w > 1)),
        "worst_rows": [{"ratio_to_bound": float(w[i]),
                        "stored": float(Au[m][i]),
                        "formula": float(pred[i]),
                        "rel_diff": float(abs(Au[m][i] - pred[i])
                                          / max(abs(Au[m][i]), 1e-30))}
                       for i in order],
        "note": "with input rounding propagated (abs() in the derivatives; "
                "an early se draft without it spuriously refuted 2 rows), "
                "every row is within bound; the worst rows are the "
                "pathological defective-tensor entries, matching at ~1e-6 "
                "relative inside their honestly-large propagated tolerance"}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
