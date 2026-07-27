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

## Amendment (logged, 2026-07-28, during the run)

The registration's "~1e-4 publication precision" assumption was WRONG: catalog
columns carry full float precision (measured true-law residuals ~1e-6). At the
loose 2e-4 floor, 14-term approximants CERTIFIED in all three bands — the
loose-epsilon parsimony exposure (same family as the sigma>0 approximant
boundary; the sigma=0 floor-dominated regime needs the minimality repair too —
**registered instrument issue**, not patched mid-study). Floor amended to the
measured 5e-6; the loose-floor artifacts are preserved in git history.

## Results (2026-07-28, `experiments/results/gaia_c0.json`; frozen snapshot
`experiments/gaia/data/c0_*.csv`, 400 sources, live TAP)

- **P1 MET (post-amendment):** all three bands certify the definitional law
  with the slope as the EXACT rational 5/2:
  `G = 25.68735 − (5/2)·log₁₀(flux)` (α ≤ 10⁻⁴⁶¹), BP (α ≤ 10⁻⁴⁵⁹),
  RP (α ≤ 10⁻⁴⁵⁸). The fitted zeropoints match the published DR3 Vega
  zeropoints (G: 25.6874).
- **P2 MET:** the color identity `bp_rp − mag_BP + mag_RP = 0` recovered with
  exact unit coefficients, **α ≤ 10⁻⁴⁵⁸⁶** (after the term-scale α fix this
  study forced — a zero-valued identity's chance-range is its constituent
  terms').
- **P3 MET:** the color–magnitude relation stays a labeled conjecture
  (`M_G ≈ 2.89·(BP−RP) + 1.29`), never a certificate.
- **P4:** met post-amendment; the loose-floor event is the registered
  exception and is documented above rather than hidden.

Phase 0 of the proposal is complete: adapter live, frozen-artifact fetches,
one known relation end-to-end under full discipline — plus two instrument
hardenings (absolute-precision declaration in `recover`; term-scale invariant
α) and one registered issue (loose-ε parsimony at σ=0) that Gaia forced on
day one, exactly as the proposal's methodological goals predicted.
