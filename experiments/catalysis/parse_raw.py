"""Parse the Mamun raw QE archives into per-calculation tables
(docs/CASE_STUDY_CATALYSIS_C0.md, path-2 rebuild).

Streams each tar.gz member (no bulk extraction), reads the final structure +
energy with ASE's espresso-out parser, and records:
file id, author, metal composition (exact counts — supercell preserved for
clean-slab pairing), adsorbate composition, final energy (eV).

Usage: parse_raw.py <archive.tar.gz> <out.csv>
"""
from __future__ import annotations

import io
import sys
import tarfile
from collections import Counter

from ase.data import chemical_symbols
from ase.io import read

ADS = {"H", "C", "N", "O", "S"}


def comp_split(atoms):
    c = Counter(atoms.get_chemical_symbols())
    metals = {el: n for el, n in c.items() if el not in ADS}
    ads = {el: n for el, n in c.items() if el in ADS}
    fmt = lambda d: "".join(f"{el}{d[el]}" for el in sorted(d))  # noqa: E731
    return fmt(metals), fmt(ads)


def main(archive: str, out_csv: str):
    n_ok = n_fail = 0
    with open(out_csv, "w") as out:
        out.write("file,metals,adsorbate,natoms,cell_c,volume,energy\n")
        with tarfile.open(archive, "r:gz") as tf:
            for m in tf:
                if not m.isfile() or not m.name.endswith(".log"):
                    continue
                try:
                    raw = tf.extractfile(m).read()
                    atoms = read(io.StringIO(raw.decode("utf-8", "replace")),
                                 format="espresso-out")
                    e = float(atoms.get_potential_energy())
                    metals, ads = comp_split(atoms)
                    cc = float(atoms.cell[2][2])
                    vol = float(atoms.get_volume())
                    out.write(f"{m.name},{metals},{ads},{len(atoms)},"
                              f"{cc:.3f},{vol:.2f},{e!r}\n")
                    n_ok += 1
                except Exception:                              # noqa: BLE001
                    n_fail += 1
                if (n_ok + n_fail) % 2000 == 0:
                    print(f"  ...{n_ok} ok / {n_fail} fail", flush=True)
    print(f"DONE {archive}: {n_ok} parsed, {n_fail} failed -> {out_csv}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
