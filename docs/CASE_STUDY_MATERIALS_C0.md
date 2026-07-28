# Materials Project C0 — registration and results

**Registered 2026-07-28, predictions before any fetch or run.** Governing
proposal: `docs/PROPOSAL_MATERIALS.md`. C0 = calibration & identity checks
on the `summary` endpoint.

## Sample (frozen query)

`summary` endpoint, stable materials (`energy_above_hull = 0`,
`theoretical = false` not required), fields: material_id, composition,
nsites, volume, density, formation_energy_per_atom, energy_per_atom,
band_gap. First 500 by material_id sort (deterministic).

## Registered predictions

- **P1 (the density identity — the definitional anchor).** The pipeline
  computes `density` from the cell composition and `volume`:
  ρ [g/cm³] = u·M_cell[amu] / V[Å³] with u = 1.66053906… (the atomic mass
  constant in g·Å³/cm³ normalization). Cell mass is deterministic feature
  construction from the returned composition and a standard atomic-mass
  table (logged; if the API's composition is per formula unit rather than
  per cell, rescaling by nsites/Σn is part of the construction and is
  noted). `recover` on (M_cell, V) → ρ certifies the exact monomial
  `u·x₀/x₁` — exponents 1, −1 — at the v2 floor (σ_rep if columns are
  rounded), with the constant decoding u to the pipeline's precision. If
  the atomic-mass table the pipeline uses differs from ours at the floor,
  the residual decode reports which elements disagree (a pipeline-decode
  finding, not a failure).
- **P2 (band gaps must NOT certify).** `recover` on
  (formation_energy_per_atom → band_gap) and (density → band_gap) must
  structurally abstain; the fit scout's conjectures stay labeled empirical.
  A certificate here would be a confident-wrong.
- **P3 (homogeneity expectation).** Because one documented pipeline
  computes every row, the P1 identity should hold across the ENTIRE
  sample — no patchwork tail (the registered contrast with PSCompPars). A
  tail, if found, is a finding about MP's ingestion history and is
  reported with its decode.
- **P4.** Zero confident-wrong; certificates clear α ≤ 10⁻⁶.

## Results

(after the run)
