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

## Results (2026-07-28, `experiments/results/materials_c0.json`; frozen
snapshot `c0_aec1e4a889cb.json`, 500 stable materials)

- **P1 MET — through the registered decode clause, which turned out to be
  the finding.** With the IUPAC-2021 mass table, `recover` correctly
  abstained: 95% of rows composed at float precision, but 5% deviated by up
  to 4×10⁻⁴. The registered residual decode then solved the PIPELINE's own
  mass table from the data (ρV/u = Σ n·m is linear in the masses;
  post-fit residual **3×10⁻¹⁴** — machine precision across all 500 rows):
  eleven elements disagree with IUPAC 2021 (Zn 65.409, S 32.065, Se 78.96,
  Mo 95.94, Ge 72.64, Cd 112.411, Yb 173.04, Li 6.941, B 10.811,
  Tl 204.3833, Hg 200.59) — **recognizably the ~2005-vintage IUPAC standard
  atomic weights, the table pymatgen ships.** With the decoded table:
  **certified, `ρ = 1.66053906888·M_cell/V`, α ≤ 10⁻⁷¹⁵** — exponents
  exactly 1 and −1, and the pinned constant is the atomic mass constant
  u = 1.66053906660×10⁻²⁴ g to nine significant figures.
- **P2 MET:** both band-gap cells structurally abstain; conjecture track
  only. The registered confident-wrong bait was not taken.
- **P3 MET, in the deep sense.** The homogeneity expectation held: there is
  NO patchwork tail. The apparent 5% tail was entirely the mass-table
  vintage — one documented pipeline, one internally consistent constant
  table, machine-precision composition once decoded. The registered
  contrast with PSCompPars is confirmed.
- **P4 MET:** zero confident-wrong.

**The sentence this stage earns:** the instrument refused to certify the
density identity against the modern atomic-mass table, and its refusal
residuals reconstructed, element by element, the twenty-year-old IUPAC
table the pipeline actually uses — then certified the identity at
α ≤ 10⁻⁷¹⁵ with the atomic mass constant read out to nine digits.
