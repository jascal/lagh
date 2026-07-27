# Gaia DR3 — Phase 0 / C0 registration (proposal: docs/PROPOSAL_GAIA.md)

**Registered 2026-07-28, predictions before any discovery run.** Phase 0 per
the proposal: the ADQL adapter + one known relation end-to-end under full
discipline. Adjustments to the proposal recorded here: C0 anchors on an
EXACTLY-LAWFUL definitional identity (catalog magnitudes are computed from
fluxes with fixed zeropoints), not on empirical physics; the binaries stage
(C3) is reframed as mass-dependent-relation discovery since Kepler III is not
a single cross-catalog law.

## C0 sample (frozen query)

`gaiadr3.gaia_source`, TOP 400 by random_index, quality cuts `ruwe < 1.4`,
`parallax_over_error > 20`, `phot_g_mean_flux_over_error > 100`, columns:
G/BP/RP magnitudes + fluxes + errors, parallax, bp_rp. Snapshot committed
(frozen artifact).

## Registered predictions

- **P1 (the definitional identity):** `recover` on (log10 flux -> mag)
  certifies `G = zp_G − 2.5·log10(flux_G)` — structure exact, zp a fitted
  constant — at declared representation noise (mags published to ~1e-4).
  Same for BP and RP with their own zeropoints. α stated.
- **P2 (color definition):** the invariant scan over {bp_rp, mag_BP, mag_RP}
  certifies `bp_rp − mag_BP + mag_RP = 0` (exact by construction).
- **P3 (real physics stays honest):** the main-sequence color–magnitude
  relation (abs G vs bp_rp on the parallax-selected sample) does NOT certify —
  labeled conjecture only (astrophysical scatter, binaries, giants).
- **P4:** zero certified-wrong.

## Results

*(after the run; predictions frozen above)*
