"""Gaia DR3 ADQL adapter (Phase 0, docs/PROPOSAL_GAIA.md): query-driven
acquisition with frozen-artifact snapshots. Each fetch is written to disk with
its exact ADQL before any discovery step touches it — the registration
discipline applied to a remote archive.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

DATA = Path(__file__).parent / "data"


def fetch(adql: str, tag: str):
    """Run an ADQL query; snapshot result + query to a frozen artifact.
    Returns (astropy table, artifact path). Re-running with the same query
    reuses the snapshot (idempotent, offline-reproducible)."""
    DATA.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256(adql.encode()).hexdigest()[:12]
    base = DATA / f"{tag}_{h}"
    meta = base.with_suffix(".json")
    csv = base.with_suffix(".csv")
    if csv.exists():
        from astropy.table import Table
        return Table.read(csv), csv
    from astroquery.gaia import Gaia
    job = Gaia.launch_job(adql)
    t = job.get_results()
    t.write(csv, format="csv", overwrite=True)
    meta.write_text(json.dumps({"adql": adql, "n_rows": len(t),
                                "columns": list(t.colnames)}, indent=1))
    return t, csv


C0_ADQL = """
SELECT TOP 400 source_id, random_index,
       phot_g_mean_mag, phot_g_mean_flux, phot_g_mean_flux_error,
       phot_bp_mean_mag, phot_bp_mean_flux,
       phot_rp_mean_mag, phot_rp_mean_flux,
       bp_rp, parallax, parallax_error
FROM gaiadr3.gaia_source
WHERE ruwe < 1.4 AND parallax_over_error > 20
  AND phot_g_mean_flux_over_error > 100
  AND phot_bp_mean_flux IS NOT NULL AND phot_rp_mean_flux IS NOT NULL
ORDER BY random_index
"""
