"""Catalysis C1 (docs/CASE_STUDY_CATALYSIS_C1.md — bands frozen in git
before this file was first run): the Abild-Pedersen rational-slope test.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.mcp.core import recover  # noqa: E402

OUT = Path("experiments/results/catalysis_c1.json")
DATA = Path(__file__).parent / "data" / "mamun_rebuilt_energies.csv"
SIGMA = 0.22
PAIRS = [("C1H1", "C1", 3, 4), ("C1H2", "C1", 1, 2), ("C1H3", "C1", 1, 4),
         ("H1N1", "N1", 2, 3), ("H1O1", "O1", 1, 2), ("H1S1", "S1", 1, 2)]


def main():
    import csv
    from collections import defaultdict
    surf = defaultdict(dict)
    for r in csv.DictReader(open(DATA)):
        surf[r["surface"]][r["species"]] = float(r["dE_star"])
    out = {"sigma_declared": SIGMA}
    intercepts = {}
    n_in_band = n_consistent = 0
    for hi, lo, p, q in PAIRS:
        gamma = p / q
        x = np.array([surf[k][lo] for k in surf if hi in surf[k] and lo in surf[k]])
        y = np.array([surf[k][hi] for k in surf if hi in surf[k] and lo in surf[k]])
        A = np.column_stack([x, np.ones_like(x)])
        (s, c), *_ = np.linalg.lstsq(A, y, rcond=None)
        res = y - A @ np.array([s, c])
        se = float(np.sqrt(np.sum(res**2) / (len(x) - 2)
                           / np.sum((x - x.mean())**2)))
        r = recover(x.reshape(-1, 1).tolist(), y.tolist(), sigma=SIGMA)
        in_band = abs(s - gamma) <= 0.10
        consistent = abs(s - gamma) <= 2 * se
        n_in_band += in_band
        n_consistent += consistent
        intercepts[hi] = float(c)
        out[f"{hi}_vs_{lo}"] = {
            "n": len(x), "gamma_exact": f"{p}/{q}", "slope": float(s),
            "se": se, "in_band_pm0.10": bool(in_band),
            "rational_consistent_2se": bool(consistent),
            "scatter_ev": float(np.std(res)),
            "intercept": float(c),
            "recover_certified": r.get("certified"),
            "recover_abstain": r.get("abstain")}
    out["P1_in_band"] = f"{n_in_band}/6"
    out["P2_rational_consistent"] = f"{n_consistent}/6"
    ivals = sorted(intercepts.values())
    out["P4_intercept_range_ev"] = float(ivals[-1] - ivals[0])
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
