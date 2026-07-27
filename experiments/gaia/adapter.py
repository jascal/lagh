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
    try:
        job = Gaia.launch_job(adql)
    except Exception:                                          # noqa: BLE001
        # heavy queries (joins) exceed the sync endpoint's timeout
        job = Gaia.launch_job_async(adql)
    t = job.get_results()
    t.write(csv, format="csv", overwrite=True)
    meta.write_text(json.dumps({"adql": adql, "n_rows": len(t),
                                "columns": list(t.colnames)}, indent=1))
    return t, csv


C1_ADQL = """
SELECT TOP 500 g.source_id, g.random_index,
       ap.teff_gspphot, ap.lum_flame, ap.radius_flame, ap.mass_flame
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.astrophysical_parameters AS ap ON g.source_id = ap.source_id
WHERE g.random_index < 2000000
  AND g.ruwe < 1.4 AND g.parallax_over_error > 20
  AND ap.teff_gspphot IS NOT NULL AND ap.lum_flame IS NOT NULL
  AND ap.radius_flame IS NOT NULL AND ap.mass_flame IS NOT NULL
ORDER BY g.random_index
"""

C2_ADQL = """
SELECT TOP 3000 source_id, random_index, ra, dec, parallax,
       pmra, pmdec, radial_velocity, bp_rp
FROM gaiadr3.gaia_source
WHERE radial_velocity IS NOT NULL AND parallax > 5
  AND ruwe < 1.4 AND parallax_over_error > 20
  AND bp_rp IS NOT NULL
ORDER BY random_index
"""

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
