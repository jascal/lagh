"""NASA Exoplanet Archive TAP adapter (Phase 0, docs/PROPOSAL_EXOPLANET.md):
query-driven acquisition with frozen-artifact snapshots, the Gaia adapter
pattern on the Archive's TAP service."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

DATA = Path(__file__).parent / "data"
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


def fetch(adql: str, tag: str):
    """Run an ADQL query against the Archive TAP; snapshot result + query to a
    frozen artifact. Idempotent: re-running with the same query reuses the
    snapshot (offline-reproducible)."""
    DATA.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256(adql.encode()).hexdigest()[:12]
    base = DATA / f"{tag}_{h}"
    meta = base.with_suffix(".json")
    csv = base.with_suffix(".csv")
    if csv.exists():
        from astropy.table import Table
        return Table.read(csv), csv
    import requests
    r = requests.get(TAP, params={"query": " ".join(adql.split()),
                                  "format": "csv"}, timeout=300)
    r.raise_for_status()
    csv.write_text(r.text)
    from astropy.table import Table
    t = Table.read(csv)
    meta.write_text(json.dumps({"adql": adql, "n_rows": len(t),
                                "columns": list(t.colnames)}, indent=1))
    return t, csv


C1_ADQL = """
SELECT pl_name, pl_rade, pl_radeerr1, pl_radeerr2,
       pl_bmasse, pl_bmasseerr1, pl_bmasseerr2, pl_bmassprov,
       pl_orbper, pl_dens, discoverymethod,
       st_teff, st_mass, st_met
FROM pscomppars
WHERE pl_rade IS NOT NULL AND pl_orbper IS NOT NULL
"""

PH2_ADQL = """
SELECT hostname, pl_name, sy_pnum, pl_orbper, pl_orbsmax,
       pl_rade, pl_radeerr1, pl_radeerr2,
       pl_bmasse, pl_bmasseerr1, pl_bmasseerr2, pl_bmassprov,
       pl_dens, pl_insol, st_mass, st_teff, discoverymethod
FROM pscomppars
WHERE pl_orbper IS NOT NULL
"""

C0_ADQL = """
SELECT pl_name, pl_bmasse, pl_rade, pl_orbper, pl_orbsmax, pl_dens,
       pl_insol, st_mass, st_rad, st_teff, st_lum, st_met, sy_pnum,
       discoverymethod,
       pl_dens_reflink, pl_orbsmax_reflink, pl_insol_reflink
FROM pscomppars
WHERE pl_bmasse IS NOT NULL AND pl_rade IS NOT NULL
  AND pl_dens IS NOT NULL AND pl_orbper IS NOT NULL
"""
