"""Materials Project C0 (docs/CASE_STUDY_MATERIALS_C0.md — predictions
frozen): the density identity certifies; band gaps must not; homogeneity
expectation tested.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.exoplanet.run_c0 import half_step_rel  # noqa: E402
from experiments.materials.adapter import C0_PARAMS, fetch  # noqa: E402
from experiments.materials.masses import MASSES  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402
from lagh.submit import submission  # noqa: E402

OUT = Path("experiments/results/materials_c0.json")
U_AMU = 1.66053906660  # g*A^3/cm^3 per amu/A^3 (expected constant decode)


def main():
    out = {}
    docs, snap = fetch("/materials/summary/", C0_PARAMS, "c0")
    out["snapshot"] = str(snap)
    out["n"] = len(docs)
    Mcell, V, rho, fe, gap, skipped = [], [], [], [], [], 0
    for d in docs:
        comp = d.get("composition") or {}
        try:
            m = sum(MASSES[el] * float(n) for el, n in comp.items())
        except KeyError:
            skipped += 1
            continue
        nsites = float(d.get("nsites") or 0)
        tot = sum(float(n) for n in comp.values())
        if tot > 0 and nsites > 0 and abs(tot - nsites) > 1e-9:
            m *= nsites / tot           # composition per formula unit -> cell
        Mcell.append(m)
        V.append(float(d["volume"]))
        rho.append(float(d["density"]))
        fe.append(float(d.get("formation_energy_per_atom") or np.nan))
        gap.append(float(d.get("band_gap") or 0.0))
    out["n_skipped_unknown_element"] = skipped
    Mcell, V, rho = map(np.asarray, (Mcell, V, rho))
    fe, gap = np.asarray(fe), np.asarray(gap)

    # P1: density identity via sigma_rep (v2 procedure, RAW columns)
    sig = half_step_rel(rho) + half_step_rel(Mcell) + half_step_rel(V)
    r = recover(np.column_stack([Mcell, V]).tolist(), rho.tolist(),
                sigma=float(max(sig, 1e-12)))
    out["P1_density_identity"] = {
        "n": len(rho), "sigma_rep": float(sig),
        "certified": r.get("certified"), "law": (r.get("law") or "")[:120],
        "alpha_log10": r.get("alpha_log10"), "abstain": r.get("abstain")}
    out["P1_expected_u"] = U_AMU
    # P3: homogeneity -- residual census of the analytic identity
    res = rho - U_AMU * Mcell / V
    rel = np.abs(res) / rho
    out["P3_homogeneity"] = {
        "median_rel": float(np.median(rel)),
        "q99_rel": float(np.quantile(rel, 0.99)),
        "max_rel": float(np.max(rel)),
        "frac_within_1e-4": float(np.mean(rel < 1e-4))}

    # P1 decode (registered clause): solve the PIPELINE's mass table from
    # rho*V/u = sum n_el * m_el (linear, overdetermined), snap to short
    # decimals, report disagreements with IUPAC 2021, re-certify with the
    # decoded table
    if not r.get("certified"):
        els = sorted({el for d in docs for el in (d.get("composition") or {})
                      if el in MASSES})
        idx = {el: i for i, el in enumerate(els)}
        A = np.zeros((len(Mcell), len(els)))
        row = 0
        for d in docs:
            comp = d.get("composition") or {}
            if any(el not in MASSES for el in comp):
                continue
            nsites = float(d.get("nsites") or 0)
            tot = sum(float(n) for n in comp.values())
            sc = nsites / tot if tot > 0 and nsites > 0 else 1.0
            for el, n in comp.items():
                A[row, idx[el]] = float(n) * sc
            row += 1
        bvec = rho * V / U_AMU
        m_fit, *_ = np.linalg.lstsq(A, bvec, rcond=None)
        m_snap = np.round(m_fit, 6)
        out["P1_mass_table_decode"] = {
            "post_fit_max_rel": float(np.max(np.abs(bvec - A @ m_fit) / bvec)),
            "elements_disagreeing_with_iupac2021": {
                el: {"pipeline": float(m_snap[idx[el]]),
                     "iupac2021": MASSES[el]}
                for el in els
                if abs(m_snap[idx[el]] - MASSES[el]) > 1e-3
                and (A[:, idx[el]] > 0).sum() >= 3},
            "note": "the decoded values are the ~2005-vintage IUPAC standard "
                    "atomic weights (the table pymatgen ships)"}
        M2 = A @ m_snap
        sig2 = half_step_rel(rho) + half_step_rel(M2) + half_step_rel(V)
        r2 = recover(np.column_stack([M2, V]).tolist(), rho.tolist(),
                     sigma=float(max(sig2, 1e-12)))
        out["P1_with_decoded_table"] = {
            "certified": r2.get("certified"),
            "law": (r2.get("law") or "")[:120],
            "alpha_log10": r2.get("alpha_log10"),
            "abstain": r2.get("abstain")}

    # P2: band gaps must NOT certify
    for xn, x in (("formation_energy", fe), ("density", rho)):
        ok = np.isfinite(x) & np.isfinite(gap)
        rr = recover(x[ok].reshape(-1, 1).tolist(), gap[ok].tolist(),
                     sigma=0.0)
        sub = submission(x[ok].reshape(-1, 1), gap[ok], sigma=0.0)
        out[f"P2_gap_vs_{xn}"] = {"certified": rr.get("certified"),
                                  "abstain": rr.get("abstain"),
                                  "track": sub["track"]}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
