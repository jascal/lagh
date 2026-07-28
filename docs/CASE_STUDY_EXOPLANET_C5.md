# Exoplanet Archive C5 — open/stress exploration (campaign finale)

**Registered 2026-07-27, predictions before any fetch or run.** Governing
proposal: `docs/PROPOSAL_EXOPLANET.md`, final stage. The Gaia-C6 protocol
applied to the composite table: a blind pair sweep where the only
certificates permitted to survive are pipeline definitions, plus a stress
test of the instrument's own protections on this archive's data.

## Sample (frozen query, `C5_ADQL`)

`pscomppars`, rows with planet and stellar radius, 19 physical/system
columns plus the `_reflink` provenance for the three transit-geometry
columns the Archive computes.

## Protocol (fixed here)

The Gaia-C6 sweep: every ordered column pair (x → y) over the 16 physical
columns is scouted with the five-form OLS scout; pairs whose best form
reproduces y to ≤ 10⁻⁶ relative advance to `recover`. Monomial identities
route through σ_rep in linear space (the C0 lesson; composed per-column
decimal half-steps). Three dim-2 triples advance unconditionally, on their
"Calculated Value" strata where the flag exists:

- (pl_rade, st_rad) → pl_ratror — the radius-ratio definition; expected
  `ratror = R_p/(109.076·R★)` (R⊕ per R☉), exponents 1, −1.
- (pl_orbsmax, st_rad) → pl_ratdor — the scaled semi-major axis; expected
  `ratdor = 215.032·a/R★` (AU per R☉), exponents 1, −1.
- (pl_ratror) → pl_trandep — the transit depth; expected `depth ∝ ratror²`,
  exponent exactly 2 (constant fixes the units convention).

## Registered predictions

- **P1 (the registered definitional triples certify on their Calculated
  strata).** Exact exponents as above; constants decode the Archive's unit
  conventions (109.076, 215.032, and the depth unit). If a stratum is a
  literature mixture (the C0/P3 pattern), the honest abstain + residual
  decode replaces the certificate.
- **P2 (census under heavy abstention).** Certificates emerge ONLY for
  pipeline definitions; every astrophysical pair is a conjecture or an
  abstention; any unregistered certificate must decode into the pipeline
  documentation and is disclosed either way.
- **P3 (stress: the protections hold on this archive).** (a) The C0 density
  identity re-run on the FULL mixture (no provenance flag) must abstain.
  (b) The density identity re-run at 10× the measured σ_rep must NOT
  certify a multi-term approximant: acceptable outcomes are the true 2-term
  monomial or an explicit abstain — never a whale (the loose-ε closure,
  validated in the field on a second archive).
- **P4.** Zero confident-wrong; every certificate clears α ≤ 10⁻⁶.

## Results

(after the run)
