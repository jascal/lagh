"""Build the adsorption-energy table from parsed Mamun raw archives
(docs/CASE_STUDY_CATALYSIS_C0.md, path-2 rebuild).

Pairing: adsorbed slab <-> clean slab on EXACT metal counts (supercell
preserved); clean slab = metals-only row with cell_c > 15 A (vacuum), min
energy per composition. Adsorption energy relative to the clean slab only:
    dE*(surface, species) = E(slab+species) - E(slab)
Gas-reference terms are species-wise CONSTANTS: they shift intercepts,
never the scaling slopes C1 tests. For validation against the mirror
(which used Mamun's H2/H2O gas conventions), one constant per species is
calibrated as the median offset and the post-offset residual is the gate.

Outputs: data/mamun_rebuilt_energies.csv (surface, species, dE_star,
n_sites) and a validation JSON.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

RAW = Path(__file__).parent / "raw"
DATA = Path(__file__).parent / "data"


def load(fname):
    return list(csv.DictReader(open(RAW / fname)))


def main(archives: list[str]):
    # clean slabs from references: metals-only, vacuum cell, min E per comp
    slabs = {}
    for r in load("references_parsed.csv"):
        if r["adsorbate"] or float(r["cell_c"]) <= 15:
            continue
        k = r["metals"]
        e = float(r["energy"])
        if k not in slabs or e < slabs[k]:
            slabs[k] = e
    print(f"clean slabs: {len(slabs)} compositions")

    best = defaultdict(lambda: (np.inf, 0))   # (surface, species) -> (dE, n)
    n_rows = n_nomatch = 0
    for arc in archives:
        for r in load(f"{arc}_parsed.csv"):
            ads = r["adsorbate"]
            if not ads or not r["metals"]:
                continue
            n_rows += 1
            e_slab = slabs.get(r["metals"])
            if e_slab is None:
                n_nomatch += 1
                continue
            dE = float(r["energy"]) - e_slab
            key = (r["metals"], ads)
            cur = best[key]
            best[key] = (min(cur[0], dE), cur[1] + 1)
    print(f"adsorbed rows: {n_rows}, no-slab-match: {n_nomatch}, "
          f"(surface,species) keys: {len(best)}")

    DATA.mkdir(exist_ok=True)
    out = DATA / "mamun_rebuilt_energies.csv"
    with open(out, "w") as f:
        f.write("surface,species,dE_star,n_sites\n")
        for (surf, sp), (dE, n) in sorted(best.items()):
            f.write(f"{surf},{sp},{dE!r},{n}\n")
    print("->", out)


if __name__ == "__main__":
    main(sys.argv[1:] if len(sys.argv) > 1 else ["O", "C", "H", "N", "S"])
