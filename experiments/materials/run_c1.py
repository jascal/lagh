"""Materials Project C1 (docs/CASE_STUDY_MATERIALS_C1.md — predictions
frozen): VRH identities certify; elastic scaling stays banded conjecture.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.run_c0 import half_step_rel  # noqa: E402
from experiments.materials.adapter import fetch  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402

OUT = Path("experiments/results/materials_c1.json")

C1_PARAMS = {
    "_fields": "material_id,composition,nsites,volume,density,"
               "bulk_modulus,shear_modulus",
    "energy_above_hull_max": 0.0,
    "_limit": 2000,
    "_sort_fields": "material_id",
}


def main():
    out = {}
    docs, snap = fetch("/materials/summary/", C1_PARAMS, "c1")
    out["snapshot"] = str(snap)
    out["n_fetched"] = len(docs)
    rows = []
    for d in docs:
        K, G = d.get("bulk_modulus"), d.get("shear_modulus")
        if not (isinstance(K, dict) and isinstance(G, dict)):
            continue
        try:
            r = {k: float(K[k]) for k in ("voigt", "reuss", "vrh")}
            r.update({"g_" + k: float(G[k]) for k in ("voigt", "reuss", "vrh")})
        except (KeyError, TypeError, ValueError):
            continue
        r["v_atom"] = float(d["volume"]) / float(d["nsites"])
        r["comp"] = d.get("composition") or {}
        rows.append(r)
    ok = [r for r in rows
          if all(np.isfinite(v) and v > 0 for k, v in r.items()
                 if k not in ("comp",))]
    out["n_with_elasticity"] = len(rows)
    out["n_clean_positive"] = len(ok)
    Kv = np.array([r["voigt"] for r in ok])
    Kr = np.array([r["reuss"] for r in ok])
    Kh = np.array([r["vrh"] for r in ok])
    Gv = np.array([r["g_voigt"] for r in ok])
    Gr = np.array([r["g_reuss"] for r in ok])
    Gh = np.array([r["g_vrh"] for r in ok])
    Va = np.array([r["v_atom"] for r in ok])

    # P1: VRH identities
    for tag, a, b, h in (("P1_K_vrh", Kv, Kr, Kh), ("P1_G_vrh", Gv, Gr, Gh)):
        sig = half_step_rel(h) + half_step_rel(a) + half_step_rel(b)
        r = recover(np.column_stack([a, b]).tolist(), h.tolist(),
                    sigma=float(max(sig, 1e-12)))
        out[tag] = {"n": len(h), "sigma_rep": float(sig),
                    "certified": r.get("certified"),
                    "law": (r.get("law") or "")[:100],
                    "alpha_log10": r.get("alpha_log10"),
                    "abstain": r.get("abstain")}

    # P2: bound ordering
    viol = int(np.sum(~((Kr <= Kh + 1e-9) & (Kh <= Kv + 1e-9)))) \
        + int(np.sum(~((Gr <= Gh + 1e-9) & (Gh <= Gv + 1e-9))))
    out["P2_ordering_violations"] = viol

    # P3: binary-oxide Birch/Anderson slope
    def is_binary_oxide(c):
        els = set(c)
        return len(els) == 2 and "O" in els
    m = np.array([is_binary_oxide(r["comp"]) for r in ok])
    lv, lk = np.log10(Va[m]), np.log10(Kh[m])
    A = np.column_stack([lv, np.ones_like(lv)])
    (s3, _), *_ = np.linalg.lstsq(A, lk, rcond=None)
    r3 = recover(lv.reshape(-1, 1).tolist(), lk.tolist(), sigma=0.0)
    out["P3_oxide_scaling"] = {"n": int(m.sum()), "slope": float(s3),
                               "recover_certified": r3.get("certified"),
                               "recover_abstain": r3.get("abstain")}

    # P4: Pugh ratio
    out["P4_pugh"] = {"median_G_over_K": float(np.median(Gh / Kh)),
                      "q10": float(np.quantile(Gh / Kh, 0.1)),
                      "q90": float(np.quantile(Gh / Kh, 0.9))}

    # P5: cross-family scaling must not certify
    lva, lka = np.log10(Va), np.log10(Kh)
    A = np.column_stack([lva, np.ones_like(lva)])
    (s5, _), *_ = np.linalg.lstsq(A, lka, rcond=None)
    r5 = recover(lva.reshape(-1, 1).tolist(), lka.tolist(), sigma=0.0)
    out["P5_cross_family"] = {"n": len(lva), "scout_slope": float(s5),
                              "recover_certified": r5.get("certified"),
                              "recover_abstain": r5.get("abstain")}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
