"""Gaia C0 (docs/CASE_STUDY_GAIA_C0.md — predictions frozen): the definitional
photometric identity certified from live-archive data; real physics stays a
labeled conjecture."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.gaia.adapter import C0_ADQL, fetch  # noqa: E402
from lagh.mcp.core import recover  # noqa: E402
from lagh.submit import submission  # noqa: E402
from lagh.systems import discover_invariants  # noqa: E402

OUT = Path("experiments/results/gaia_c0.json")
SIGMA_MAG = 2e-4          # magnitudes published to ~1e-4; declared envelope


def main():
    t, snap = fetch(C0_ADQL, "c0")
    out = {"snapshot": str(snap), "n": len(t)}
    # P1: mag = zp - 2.5 log10(flux), per band
    for band, magn, fl in (("G", "phot_g_mean_mag", "phot_g_mean_flux"),
                           ("BP", "phot_bp_mean_mag", "phot_bp_mean_flux"),
                           ("RP", "phot_rp_mean_mag", "phot_rp_mean_flux")):
        lf = np.log10(np.asarray(t[fl], float))
        mag = np.asarray(t[magn], float)
        r = recover(lf.reshape(-1, 1).tolist(), mag.tolist(), sigma=SIGMA_MAG)
        out[f"P1_{band}"] = {"certified": r.get("certified"),
                             "law": r.get("law"),
                             "alpha_log10": r.get("alpha_log10"),
                             "abstain": r.get("abstain")}
    # P2: color identity bp_rp - BP + RP = 0
    cols = {"bp_rp": np.asarray(t["bp_rp"], float),
            "magBP": np.asarray(t["phot_bp_mean_mag"], float),
            "magRP": np.asarray(t["phot_rp_mean_mag"], float)}
    invs = discover_invariants(cols, sigma=SIGMA_MAG)
    lin = [iv for iv in invs if iv["n_terms"] == 3 and "**" not in iv["expr"]
           and "log" not in iv["expr"] and iv["expr"].count("*") <= 3]
    out["P2_color_identity"] = [
        {"expr": iv["expr"], "alpha_log10": iv["alpha_log10"]} for iv in lin[:2]]
    # P3: main-sequence CMD -- must NOT certify
    plx = np.asarray(t["parallax"], float)
    absG = np.asarray(t["phot_g_mean_mag"], float) + 5 * np.log10(plx / 100.0)
    bprp = np.asarray(t["bp_rp"], float)
    sub = submission(bprp.reshape(-1, 1), absG, sigma=0.0)
    out["P3_cmd"] = {"track": sub["track"], "tag": sub["tag"],
                     "expr": (sub.get("expr") or "")[:90]}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
