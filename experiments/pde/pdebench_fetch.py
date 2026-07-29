"""Fetch a small, reproducible EXTRACT of a PDEBench file (docs/PDEBENCH_READINESS.md).

    .venv/bin/python experiments/pde/pdebench_fetch.py --id 255666 --samples 6 \
        --out data/pdebench/1D_Advection_beta0.7_x6.h5

The published 1-D files are ~8 GB of 10000 samples each and the instrument needs
a handful. The DaRUS storage backend honours HTTP range requests and the
datasets are stored CONTIGUOUS (chunks=None), so h5py over an fsspec file reads
exactly the bytes of the samples asked for -- a few MB instead of a few GB.

The extract is written in the SAME layout as the source, so nothing downstream
knows the difference, and it carries its provenance as HDF5 attributes: the DOI,
the file id, the source filename, which sample indices were taken, and any
geometry quirk found on the way. A run whose data cannot be traced back to a
published file is not a result.

One quirk found on first contact and handled explicitly, never silently: the
1-D files ship a `t-coordinate` with one MORE entry than the tensor has time
slices (202 against 201). Truncating quietly would leave every patch integral
using a time grid the data does not have; the extract records the dropped value.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

DOI = "doi:10.18419/darus-2986"
API = "https://darus.uni-stuttgart.de/api/access/datafile/{id}"
CATALOG = ("https://darus.uni-stuttgart.de/api/datasets/:persistentId/"
           "versions/:latest/files?persistentId=" + DOI)


def catalog(cache: Path | None = None) -> list:
    """(size, id, path) for every published file, so a run can name its source."""
    import urllib.request
    if cache and cache.exists():
        d = json.loads(cache.read_text())
    else:
        with urllib.request.urlopen(CATALOG, timeout=120) as r:
            d = json.loads(r.read().decode())
        if cache:
            cache.write_text(json.dumps(d))
    out = []
    for f in d["data"]:
        df = f["dataFile"]
        out.append((df.get("filesize", 0), df["id"],
                    f"{f.get('directoryLabel', '')}/{df['filename']}"))
    return sorted(out)


def fetch(file_id: int, samples: int, out: Path, *, tmax_index=None,
          block=2 ** 20) -> dict:
    import fsspec
    import h5py
    url = API.format(id=file_id)
    fs = fsspec.filesystem("http", client_kwargs={"trust_env": True})
    notes = []
    with h5py.File(fs.open(url, block_size=block), "r") as h:
        keys = list(h.keys())
        if "tensor" not in h:
            raise SystemExit(f"unexpected layout {keys}: this fetcher handles "
                             "the 1-D single-field files")
        ten = h["tensor"]
        n_s, n_t, n_x = ten.shape
        x = np.asarray(h["x-coordinate"][()], np.float32)
        t = np.asarray(h["t-coordinate"][()], np.float32)
        dropped = []
        if len(t) != n_t:
            # DECLARED, not silently truncated: the published files ship one
            # more time coordinate than the tensor has slices
            dropped = [float(v) for v in t[n_t:]]
            notes.append(f"t-coordinate had {len(t)} entries against {n_t} "
                         f"time slices; kept the first {n_t}, dropped {dropped}")
            t = t[:n_t]
        keep = min(samples, n_s)
        data = np.asarray(ten[:keep], np.float32)
        if tmax_index:
            data = data[:, :tmax_index]
            t = t[:tmax_index]
            notes.append(f"time axis truncated to {tmax_index} slices")
    out.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out, "w") as g:
        g.create_dataset("tensor", data=data)
        g.create_dataset("x-coordinate", data=x)
        g.create_dataset("t-coordinate", data=t)
        g.attrs["source_doi"] = DOI
        g.attrs["source_file_id"] = int(file_id)
        g.attrs["source_shape"] = json.dumps([int(n_s), int(n_t), int(n_x)])
        g.attrs["samples_taken"] = json.dumps(list(range(keep)))
        g.attrs["notes"] = json.dumps(notes)
    return {"file_id": file_id, "out": str(out), "shape": list(data.shape),
            "source_shape": [n_s, n_t, n_x], "notes": notes,
            "x_span": [float(x[0]), float(x[-1])],
            "t_span": [float(t[0]), float(t[-1])]}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", type=int, help="DaRUS datafile id")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--samples", type=int, default=6)
    ap.add_argument("--tmax-index", type=int, default=None)
    ap.add_argument("--list", action="store_true",
                    help="print the published catalog instead of fetching")
    ap.add_argument("--grep", default="")
    a = ap.parse_args(argv)
    if a.list:
        for s, i, p in catalog():
            if a.grep.lower() in p.lower():
                print(f"{s/1e6:10.1f} MB  id={i:7d}  {p}")
        return 0
    if not (a.id and a.out):
        ap.error("--id and --out are required unless --list")
    info = fetch(a.id, a.samples, a.out, tmax_index=a.tmax_index)
    print(json.dumps(info, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
