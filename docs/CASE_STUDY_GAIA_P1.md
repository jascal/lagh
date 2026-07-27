# Gaia DR3 Phase 1 (C1–C2) — registration and results

**Registered 2026-07-27, predictions before any fetch or run.** Governing
proposal: `PROPOSAL_GAIA.md`. Phase 0 (C0) closed in `CASE_STUDY_GAIA_C0.md`,
including the loose-ε parsimony issue it surfaced — now fixed at instrument
level (floor-dominated regime: refit-parsimony collapse + winner-level
coefficient gate + sticky re-split ambiguity; see `lagh/engine.py`,
`lagh/certify.py::refit_minimal`, `lagh/passive.py`). P2 below spends that fix
against the archive deliberately.

## C1 — astrophysical-parameter relations (new fetch: FLAME join)

**ADQL (frozen in `experiments/gaia/adapter.py::C1_ADQL`):** `gaia_source` ⋈
`astrophysical_parameters`, quality cuts `ruwe < 1.4`,
`parallax_over_error > 20`, all of `teff_gspphot`, `lum_flame`,
`radius_flame`, `mass_flame` non-null, ordered by `random_index`, TOP 500.

**Floor procedure (fixed here, C0 lesson):** before the discovery run, measure
each column's decimal quantization; `floor_abs` for a log-space target is the
propagated log-rounding `(r_L + 2·r_R + 4·r_T)/ln 10` with a 2× safety margin.
The measured value is recorded in results; no other amendment is anticipated.

- **P1 (Stefan–Boltzmann, pipeline-definitional).** FLAME derives its
  radius/luminosity chain through `L = R²·(T/T☉)⁴`. `recover` on
  (log₁₀ radius_flame, log₁₀ teff_gspphot) → log₁₀ lum_flame **certifies**
  `2·x₀ + 4·x₁ + c` — the exact exponents 2 and 4 — at the measured floor,
  with α stated. The constant should decode as `−4·log₁₀(5772)` (IAU nominal
  T☉; checked post-hoc, not part of the certificate claim). If a subpopulation
  breaks the identity (FLAME flag heterogeneity), the honest verdict is
  abstain; a flags-based cut is the only amendment allowed, and must be logged.
- **P2 (loose-floor discipline, in-field check of the C0 fix).** The same P1
  recovery run at a deliberately loose `floor_abs = 1e-2` must NOT certify a
  multi-term approximant: the only acceptable outcomes are a certificate of
  the same minimal 3-term form or an explicit abstain (parametric expected —
  coefficients are honestly unpinnable at a loose floor).
- **P3 (mass–luminosity is not a definitional identity).** On the
  main-sequence window `0.5 ≤ M/M☉ ≤ 2`, `recover` (log₁₀ mass_flame →
  log₁₀ lum_flame) must NOT certify; the fit scout reports a labeled
  conjecture with slope in the classical band **[3.0, 4.5]**.

## C2 — solar-neighborhood 6D kinematics (new fetch: RV sample)

**ADQL (frozen, `C2_ADQL`):** `gaia_source` with `radial_velocity` non-null,
`parallax > 5` (d < 200 pc), `ruwe < 1.4`, `parallax_over_error > 20`,
columns (ra, dec, parallax, pmra, pmdec, radial_velocity, bp_rp), ordered by
`random_index`, TOP 3000. Heliocentric UVW computed by the standard Galactic
transformation (deterministic feature construction, `astropy`).

- **P4 (no deterministic kinematic law).** `recover` on velocity components
  (U → V, U → W) must structurally abstain — phase-mixed scatter is not a law.
- **P5 (Strömberg asymmetric drift, labeled conjecture).** Binning by bp_rp
  color (8 bins), per-bin mean V vs per-bin σ_U²: a **negative linear trend**
  (redder/older → larger dispersion → larger lag), labeled empirical; the
  implied Strömberg constant `k = σ_U²/Δv_a` in the band **[50, 150] km/s**.
  `recover` on the 8 binned points (tiny-data mode) must NOT certify.
- **P6.** Zero confident-wrong across the phase.

**Execution amendment (logged, 2026-07-27, before any result):** the C1 join
timed out on the synchronous TAP endpoint (HTTP 408; no data returned). The
fetch now falls back to an async job, and the join is bounded by
`g.random_index < 2000000` (still random-index selection, same TOP 500).
Predictions unchanged.

## Results (2026-07-27, `experiments/results/gaia_p1.json`; frozen snapshots
`c1_85a589d51938.csv` (500 sources), `c2_5c940494e92b.csv` (3000 sources))

- **P1 MISSED — and the structural abstain is the correct verdict.** The FLAME
  trio does NOT satisfy Stefan–Boltzmann exactly as published: residuals of
  `log L − (2·log R + 4·log T − 4·log 5772)` have median |res| ≈ 2×10⁻³ dex
  with tails to 0.3 dex; only 0.8% of sources are within 10⁻⁵ dex. The
  registered flags-cut fallback was exercised (logged fetch of `flags_flame`,
  same window): both flag classes ('00', '10') sit at millidex residuals — no
  cut rescues certification. Decode: the implied Teff `5772·(L/R²)^¼` differs
  from `teff_gspphot` by ~6 K median (tails 700 K) — the published columns are
  per-parameter posterior point estimates, and point estimates of different
  marginals do not compose through a nonlinear identity even when each
  posterior sample satisfies it. The instrument refused to certify a relation
  the catalog genuinely does not satisfy at certification precision — the
  zero-confident-wrong invariant doing its job on real data. (The registered
  floor procedure also misread float32-quantized columns as full float64 —
  the measured floor 1e-9 was too tight, but the abstain stands at any floor
  below the millidex scatter, so no amendment changes the verdict.)
- **P2 MET — the loose-floor fix held in the field.** At the deliberately
  loose `floor_abs = 1e-2` the run returned a structural abstain: no
  multi-term approximant certified. This is the in-field validation of the
  C0 loose-ε closure.
- **P3 half-MET.** Conjecture-only as predicted (no certification). The slope
  prediction MISSED: 5.12 vs the registered band [3.0, 4.5]. Honest reading:
  `mass_flame` is itself derived from `lum_flame` through FLAME's isochrone
  chain, so this measures FLAME's internal mass–luminosity mapping, not the
  empirical (binary-calibrated) relation the classical band describes.
- **P4 MET:** U→V and U→W both structurally abstain — no deterministic law in
  phase-mixed kinematics.
- **P5 MET:** Strömberg asymmetric drift: per-color-bin mean V vs σ_U² has
  corr −0.70, negative slope, implied **k = 108.9 km/s** — inside the
  registered [50, 150] band; the 8-point binned `recover` honestly abstains
  (labeled conjecture only).
- **P6 MET:** zero confident-wrong across the phase.

**The sentence this phase earns:** on its first contact with derived
astrophysical columns, the instrument certified nothing — correctly: the
FLAME point estimates do not satisfy the law their own pipeline used, and the
refusal decodes into exactly why (marginal medians don't compose); meanwhile
the kinematic stratum yields the textbook Strömberg constant as a properly
labeled conjecture, and the loose-floor discipline check confirms the C0 fix
under live-archive conditions.
