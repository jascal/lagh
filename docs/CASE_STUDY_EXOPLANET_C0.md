# Exoplanet Archive Phase 0 / C0 — registration and results

**Registered 2026-07-27, predictions before any fetch or run.** Governing
proposal: `docs/PROPOSAL_EXOPLANET.md`. C0 = calibration & identity checks on
PSCompPars: the definitional stratum (columns the Archive computes from other
columns) must certify; the literature stratum must not; catalogue-consistency
findings are first-class results (the Gaia FLAME lesson).

**Provenance discipline (fixed here):** PSCompPars marks Archive-computed
values with a "Calculated Value" `_reflink`. Subsamples are selected by that
flag ONLY — never by agreement with the identity under test (that would be
circular). Floors by the v2 procedure (float round-trip quantization,
OLS-scout propagation, 2× margin); every certificate must clear the α ≤ 10⁻⁶
significance gate now built into certification.

## C0 sample (frozen query, `C0_ADQL`)

`pscomppars`, all rows with mass, radius, density, and period non-null
(~thousands), columns: planet mass/radius/period/semi-major axis/density/
insolation, stellar mass/radius/Teff/lum/met, multiplicity, discovery method,
and the three `_reflink` provenance columns for density, semi-major axis, and
insolation.

## Registered predictions

- **P1 (density is a computed column — the definitional anchor).** On the
  "Calculated Value" density subsample: `recover` on
  (log₁₀ M_⊕, log₁₀ R_⊕) → log₁₀ ρ certifies `x₀ − 3·x₁ + c` — exact
  exponents **1 and −3** — with c = log₁₀(5.513) ± (the Earth mean density
  constant, M_⊕/(4π/3 · R_⊕³), checked post-hoc against the IAU nominal
  values). On the literature-density complement: must NOT certify (mixture of
  independent measurements); abstain or conjecture only.
- **P2 (derived semi-major axis obeys Kepler III with the day/AU constant).**
  On the "Calculated Value" `pl_orbsmax` subsample: `recover` on
  (log₁₀ P_days, log₁₀ M_star) → log₁₀ a_AU certifies
  `(2/3)·x₀ + (1/3)·x₁ + c` — exact exponents **2/3 and 1/3** — with
  c = −(2/3)·log₁₀(365.25) (P in days → yr), checked post-hoc.
- **P3 (insolation identity).** On the "Calculated Value" `pl_insol`
  subsample: `recover` on (st_lum [already log₁₀ L/L☉], log₁₀ a_AU) →
  log₁₀ S certifies `x₀ − 2·x₁` — exact coefficients **1 and −2**, zero
  constant. (If the Archive's computation used a different chain, the honest
  abstain + residual decode replaces the certificate — the FLAME pattern.)
- **P4 (real physics stays honest).** The population mass–radius relation
  (all planets, log M vs log R) must NOT certify — labeled conjecture or
  abstention only; no registered slope band (the relation is known to be
  broken/multivalued across regimes — a certificate here would be a
  confident-wrong).
- **P5.** Zero confident-wrong; every certificate carries α ≤ 10⁻⁶.

**Execution amendments (logged, 2026-07-27, after the first run — floor/
error-model corrections; no prediction changed):** (1) the first run fed
LOG-transformed columns to the storage-precision detector, whose constant
column round-trips float32 and produced a nonsense 4×10⁻⁸ floor — the
detector must see RAW columns. Archive columns are rounded to ~3 significant
digits: RELATIVE precision, so the identities route through σ_rep in linear
space and the power-law channel, whose rational exponents pin where float
perturbation cannot (σ_rep composed from the per-column decimal half-steps).
(2) P3 runs in linear space (L = 10^st_lum) for the same reason.

## Results (2026-07-27, `experiments/results/exoplanet_c0.json`; frozen
snapshot `c0_0387e06ec923.csv`, 5843 planets)

- **P1 MET, both halves.** On the 4320-planet "Calculated Value" stratum:
  **`ρ = 5.4961·M/R³` certified, α ≤ 10⁻¹⁹¹⁸** — the strongest bound in the
  program's history — with the exact exponents 1 and −3. The constant
  decodes: 5.496 is Earth's mean density under the ~6378 km (equatorial/IAU
  nominal) Earth-radius convention, not the volumetric 6371 km (5.513) —
  the certificate reveals which constants the Archive's pipeline uses. The
  1523-planet literature complement **structurally abstains**, exactly as
  registered: independently measured densities are not a formula.
- **P2 UNFULFILLABLE — itself a catalogue finding.** Zero of 5843 rows carry
  a "Calculated Value" reflink for `pl_orbsmax`: PSCompPars takes semi-major
  axes from literature only, so the registered calculated-stratum Kepler
  test has no domain. Recorded, not forced.
- **P3 MISSED as a certificate; the abstain decodes (the FLAME pattern).**
  The calculated-insolation stratum satisfies `S = L/a²` with median offset
  10⁻⁸ dex and MAD 4×10⁻⁶ dex — but 1.3% of rows deviate beyond 10⁻² dex,
  and exhaustive certification correctly refuses: the composite table's
  per-column literature patchwork leaves a stratum whose inputs (st_lum, a)
  are not the ones the insolation was computed from. Point estimates that
  do not compose, third catalogue in a row.
- **P4 MET:** the population mass–radius relation stays a labeled conjecture
  (log-log slope 2.35); no certificate, as a multivalued relation demands.
- **P5 MET:** zero confident-wrong; the one certificate clears the
  significance gate by 1912 orders of magnitude.

**The sentence this phase earns:** on first contact with the composite
table, the instrument certified exactly the one column the Archive computes,
read the Archive's Earth-radius convention out of the constant, refused the
strata where the patchwork does not compose, and reported a registered
prediction as unfulfillable rather than bending it.
