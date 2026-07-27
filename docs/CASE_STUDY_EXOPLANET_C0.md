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

## Results

(after the runs)
