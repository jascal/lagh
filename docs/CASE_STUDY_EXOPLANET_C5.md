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

## Results (2026-07-27, `experiments/results/exoplanet_c5.json`; frozen
snapshot `c5_0a5089328824.csv`, 5957 planets)

- **P1 UNFULFILLABLE — the third catalogue-provenance finding.** Zero of
  5957 rows carry a "Calculated Value" reflink for ANY of the three
  transit-geometry columns: PSCompPars sources pl_ratror, pl_ratdor, and
  pl_trandep entirely from KOI tables and literature (top references:
  Morton 2016, Q1–Q17 DR24/DR25). The registered strata do not exist. The
  fallback mixture probe decodes why no certificate is possible there: the
  radius-ratio identity holds only to ~4 millidex MAD with 29% of rows off
  by more than 10⁻² dex — transit-fit ratios from one paper, radii adopted
  from another. The patchwork does not compose; the honest verdict is the
  abstention.
- **P2 MET:** the 240-pair open sweep certifies NOTHING — every pair is a
  scout-level conjecture, none reaches the exactness gate. No unregistered
  certificate appeared; heavy abstention exactly as designed.
- **P3 MET, both halves — the stress test passes on a second archive.**
  (a) The density identity on the full 5843-row mixture (no provenance
  flag): structural abstain. (b) At 10× the measured σ_rep, the certified
  law is `5.4961·M/R³` — the TRUE two-term monomial, not a whale: the
  loose-ε closure (refit-parsimony collapse + winner gate + sticky
  ambiguity) holds in the field on data it was not built against.
- **P4 MET:** zero confident-wrong; the only certificate C5 produced
  anywhere is the true density law under deliberate stress.

## Campaign closure

The exoplanet campaign (C0–C5, `PROPOSAL_EXOPLANET.md`) is complete. One
certificate — the Archive's own density formula at α ≤ 10⁻¹⁹¹⁸, constant
decoding its Earth-radius convention. Four provenance findings the
catalogue itself does not advertise: semi-major axes are never
Archive-computed; neither are any transit-geometry columns; the insolation
stratum has a 1.3% non-composing tail; the transit-geometry mixture
disagrees at the percent level. A population-science conjecture set with
every band hit or miss disclosed: the radius valley at 1.90 R⊕ sloping
−0.127 with period, the giant plateau flat to 1%, the metallicity
preference at 0.124 dex, the resonance asymmetry at 96:57 and 100:40, Hill
spacing at 17.6. And zero confident-wrong across six stages on a catalogue
built precisely of the heterogeneous literature patchwork the
zero-confident-wrong invariant exists to survive.
