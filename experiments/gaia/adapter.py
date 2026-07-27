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

C3_ADQL = """
SELECT TOP 1000 source_id, nss_solution_type, period, eccentricity,
       a_thiele_innes, b_thiele_innes, f_thiele_innes, g_thiele_innes,
       parallax
FROM gaiadr3.nss_two_body_orbit
WHERE nss_solution_type = 'Orbital' AND period IS NOT NULL
  AND a_thiele_innes IS NOT NULL AND b_thiele_innes IS NOT NULL
  AND f_thiele_innes IS NOT NULL AND g_thiele_innes IS NOT NULL
  AND parallax IS NOT NULL AND parallax > 0
"""

C4RR_ADQL = """
SELECT TOP 800 source_id, pf, phi31_g, metallicity
FROM gaiadr3.vari_rrlyrae
WHERE best_classification = 'RRab' AND pf IS NOT NULL
  AND phi31_g IS NOT NULL AND metallicity IS NOT NULL
"""

C4CEP_ADQL = """
SELECT TOP 500 v.source_id, v.pf, v.int_average_g,
       g.parallax, g.parallax_over_error
FROM gaiadr3.vari_cepheid AS v
JOIN gaiadr3.gaia_source AS g ON v.source_id = g.source_id
WHERE v.pf IS NOT NULL AND v.int_average_g IS NOT NULL
  AND v.type_best_classification = 'DCEP'
  AND v.mode_best_classification = 'FUNDAMENTAL'
  AND g.parallax_over_error > 10 AND g.parallax > 0
"""

C5_ADQL = """
SELECT TOP 3000 source_id, random_index, ra, dec, l, b, parallax, pmra, pmdec
FROM gaiadr3.gaia_source
WHERE ABS(b) < 15 AND parallax BETWEEN 0.5 AND 2
  AND ruwe < 1.4 AND parallax_over_error > 10
  AND random_index < 2000000
ORDER BY random_index
"""

C6_ADQL = """
SELECT TOP 400 source_id, random_index,
       parallax, parallax_error, parallax_over_error, pmra, pmdec,
       phot_g_mean_flux, phot_g_mean_flux_error,
       phot_g_mean_flux_over_error, phot_g_mean_mag, ruwe,
       astrometric_sigma5d_max, bp_rp
FROM gaiadr3.gaia_source
WHERE ruwe < 1.4 AND parallax_over_error > 20
  AND phot_g_mean_flux_over_error > 100
  AND phot_bp_mean_flux IS NOT NULL AND phot_rp_mean_flux IS NOT NULL
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
