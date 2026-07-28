"""Extract the frozen Mamun energy table from the Ulissi-group CatHub mirror
(docs/CASE_STUDY_CATALYSIS_C0.md, amended acquisition path).

Source: ulissigroup/uncertainty_benchmarking (public research repo for Tran
et al. 2020), preprocessing/pull_data/catalysis-hub/cathub.pkl (git-LFS),
itself a pull of Catalysis-Hub's MamunHighT2019 (CC-BY, Mamun et al. 2019).
The 296 MB pickle is NOT committed; this script derives the small energy
table (surface key, adsorbate, energy) and freezes it as CSV with the
source SHA recorded. Surface key = metals-only composition of the slab
Atoms (alloys contain no adsorbate elements), read from raw arrays to
sidestep the old-ASE pickle incompatibility.
"""
from __future__ import annotations

import hashlib
import json
import pickle
import sys
from collections import Counter
from pathlib import Path

DATA = Path(__file__).parent / "data"
Z2EL = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
ADSORBATE_Z = set(Z2EL)


def slab_key(doc):
    """Reduced metals composition of the adsorbed-slab system.

    The systems list carries gas references, ELEMENTAL bulk references, and
    the slab+adsorbate — but no clean slab (measured: a metals-only filter
    returned pure-element keys). Take the LARGEST system and strip the
    adsorbate elements; Mamun alloys contain no H/C/N/O/S metals."""
    best = None
    for a in doc.get("systems") or []:
        try:
            nums = list(a.arrays["numbers"])
        except Exception:                                      # noqa: BLE001
            continue
        if best is None or len(nums) > len(best):
            best = nums
    if best is None:
        return None
    best = [z for z in best if int(z) not in ADSORBATE_Z]
    if not best:
        return None
    c = Counter(int(z) for z in best)
    # REDUCED composition: raw counts encode supercell size, and different
    # adsorbates were run in different supercells (measured: zero O/OH joins
    # on raw formulas). gcd-reduction collapses supercells to stoichiometry;
    # distinct orderings/structures at the same stoichiometry collapse too
    # (logged in the census as a known coarsening).
    from math import gcd
    from functools import reduce
    g = reduce(gcd, c.values())
    from ase.data import chemical_symbols
    return "".join(f"{chemical_symbols[z]}{n // g if n // g > 1 else ''}"
                   for z, n in sorted(c.items()))


def main(pkl_path: str):
    DATA.mkdir(parents=True, exist_ok=True)
    raw = open(pkl_path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    p = pickle.loads(raw)
    rows = []
    for ads, docs in p.items():
        for d in docs:
            if d.get("pubId") != "MamunHighT2019":
                continue
            k = slab_key(d)
            e = d.get("energy")
            if k is None or e is None:
                continue
            rows.append((k, ads, float(e), str(d.get("Equation", ""))))
    out = DATA / "mamun_mirror_energies.csv"
    with open(out, "w") as f:
        f.write("surface,adsorbate,energy,equation\n")
        for r in sorted(rows):
            f.write(f"{r[0]},{r[1]},{r[2]!r},\"{r[3]}\"\n")
    meta = {"source_repo": "ulissigroup/uncertainty_benchmarking",
            "source_file": "preprocessing/pull_data/catalysis-hub/cathub.pkl",
            "source_sha256": sha, "n_rows": len(rows),
            "underlying": "Catalysis-Hub MamunHighT2019 (CC-BY, "
                          "Mamun et al. Sci Data 6:76, 2019)"}
    (DATA / "mamun_mirror_energies.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
