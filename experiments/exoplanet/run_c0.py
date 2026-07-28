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
from lagh.mcp.core import recover  # noqa: E402
from lagh.submit import submission  # noqa: E402

OUT = Path("experiments/results/exoplanet_c0.json")


def calc_mask(t, col):
    ref = [str(v) for v in t[col + "_reflink"]]
    return np.array(["Calculated" in v for v in ref])


def half_step_rel(col):
    """Median relative decimal-rounding half-step of a RAW linear column.
    (Amended, logged: the first run fed LOG-transformed columns to the
    storage-precision detector -- the constant column round-trips float32 and
    yielded a nonsense 4e-8 floor. Archive columns are rounded to ~3
    significant digits: RELATIVE precision, so identities route through
    sigma_rep in linear space and the power-law channel, whose rational
    exponents pin where float perturbation cannot.)"""
    rels = []
    for v in col:
        s = repr(float(v))
        if "e" in s:
            continue
        dec = len(s.split(".")[1]) if "." in s else 0
        if abs(v) > 0:
            rels.append(0.5 * 10.0 ** (-dec) / abs(v))
    return float(np.median(rels))


def rec_monomial(X0, X1, y, sigma, tag, out):
    ok = np.isfinite(X0) & np.isfinite(X1) & np.isfinite(y) \
        & (X0 > 0) & (X1 > 0) & (y > 0)
    r = recover(np.column_stack([X0[ok], X1[ok]]).tolist(), y[ok].tolist(),
                sigma=float(sigma))
    out[tag] = {"n": int(ok.sum()), "sigma_rep": float(sigma),
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
    sig1 = half_step_rel(rho[mc]) + half_step_rel(M[mc]) \
        + 3 * half_step_rel(R[mc])
    rec_monomial(M[mc], R[mc], rho[mc], sig1, "P1_density_calculated", out)
    ml = ~calc_mask(t, "pl_dens") & (M > 0) & (R > 0) & (rho > 0)
    rec_monomial(M[ml], R[ml], rho[ml], sig1, "P1_density_literature", out)
    out["P1_const_expected_volumetric_Rearth"] = 5.51295

    # P2: Kepler III -- UNFULFILLABLE as registered (logged): zero rows in
    # PSCompPars carry a Calculated reflink for pl_orbsmax; the composite
    # table takes semi-major axes from literature only. Recorded, not forced.
    k = calc_mask(t, "pl_orbsmax")
    out["P2_n_calculated"] = int(k.sum())
    out["P2_kepler_calculated"] = {
        "n": int(k.sum()), "certified": False,
        "abstain": "unfulfillable: no Calculated stratum exists for "
                   "pl_orbsmax in PSCompPars"}

    # P3: insolation identity on the Calculated pl_insol subsample,
    # linear space: S = C * L * a^-2 with L = 10**st_lum
    i = calc_mask(t, "pl_insol") & (S > 0) & (a > 0) & np.isfinite(lum)
    out["P3_n_calculated"] = int(i.sum())
    L = 10.0 ** lum
    sig3 = half_step_rel(S[i]) + half_step_rel(a[i]) * 2 \
        + np.log(10.0) * 0.001   # st_lum rounded in dex: 1e-3 half-step
    rec_monomial(L[i], a[i], S[i], sig3, "P3_insol_calculated", out)
    if not out["P3_insol_calculated"]["certified"]:
        res = np.log10(S[i]) - (lum[i] - 2 * np.log10(a[i]))
        out["P3_residual_decode"] = {
            "median_dex": float(np.median(res)),
            "mad_dex": float(np.median(np.abs(res - np.median(res)))),
            "frac_within_1e-2_dex": float(np.mean(np.abs(res) < 1e-2))}

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
