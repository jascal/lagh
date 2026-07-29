# Materials Project C2 — open sweep under heavy abstention (campaign finale)

**Registered 2026-07-28, predictions before any fetch or run.** Governing
proposal: `docs/PROPOSAL_MATERIALS.md`, final MP stage. The Gaia-C6 /
exoplanet-C5 protocol, plus three registered definitional targets that the
C1 elastic work exposed.

## Sample (frozen query)

`summary` endpoint, stable materials, first 2000 by material_id, fields:
material_id, nsites, volume, density, density_atomic,
formation_energy_per_atom, energy_per_atom, energy_above_hull, band_gap,
efermi, total_magnetization, bulk_modulus, shear_modulus,
universal_anisotropy, homogeneous_poisson.

## Protocol (fixed here)

The established sweep: every ordered pair over the physical columns is
scouted with the five-form OLS scout; pairs reproducing y to ≤ 10⁻⁶
relative advance to `recover`, with the storage-precision floor measured
per column (absolute for fixed-decimal columns — the C1 lesson). Three
registered targets advance unconditionally:

- **T1 (atomic density).** (nsites, volume) → density_atomic is expected to
  be the monomial `x₀/x₁` (or its reciprocal — the orientation is recorded,
  not assumed) with unit coefficient.
- **T2 (Poisson's ratio).** (K_VRH, G_VRH) → homogeneous_poisson is the
  textbook rational function ν = (3K − 2G)/(2(3K + G)) — integer
  coefficients throughout. `recover` attempts it (the C2 rational channel's
  home game); on abstain, the declared form goes to `verify` (the
  established second track).
- **T3 (universal anisotropy).** A_U = 5·G_voigt/G_reuss + K_voigt/K_reuss
  − 6 (Ranganathan–Ostoja-Starzewski) — dim-4, integer coefficients 5, 1,
  −6. Registered as a `verify` target directly (dim-4 rational recovery is
  outside the current reach envelope, stated as such).

## Registered predictions

- **P1.** T1–T3 certify (recover or verify-pinned) on rows where the
  fields are jointly present, at measured storage floors; defective-tensor
  rows (the C1 nine) are NOT excluded — the identities are pipeline
  arithmetic and should hold on them too.
- **P2 (census).** Certificates emerge only for pipeline definitions;
  every astrophysical—er, physical—pair (band gap, efermi, magnetization
  against anything) is a conjecture or abstention; any unregistered
  certificate must decode into the pipeline documentation and is disclosed.
- **P3.** Zero confident-wrong; certificates clear α ≤ 10⁻⁶.

**Execution amendments (logged, 2026-07-28):** (1) `verify` grew a per-point
`se` passthrough — plumbing parity with the ε model's λ_B term, aligned with
the split permutation. (2) T3's tolerance is the v2 procedure applied
properly: per-point INPUT rounding propagated through the ratio derivatives
(with absolute values — an early draft without `abs()` spuriously refuted
two negative-moduli rows; caught before any result stood).

## Results (2026-07-28, `experiments/results/materials_c2.json`; frozen
snapshot `c2_86004c2b1f30.json`, 2000 stable materials)

- **T1 CERTIFIED OUTRIGHT — the program's strongest bound:**
  `density_atomic = volume/nsites` (orientation recorded: it is volume per
  atom, Å³), exact exponents ±1, unit coefficient, **α ≤ 10⁻⁵⁰⁷⁰** over all
  2000 rows.
- **T2 split by track:** `recover` structurally abstains on Poisson's
  ratio; the declared textbook form **ν = (3K − 2G)/(2(3K + G)) verifies
  pinned** — integer coefficients, no rival within the declared floors.
- **T3 verifies pinned** — the Ranganathan–Ostoja-Starzewski
  `A_U = 5·G_v/G_r + K_v/K_r − 6` holds on all 278 elastic cert rows once
  input rounding is properly propagated (extreme-ratio and
  defective-tensor rows carry honestly large tolerances; the worst
  agreement is 4×10⁻⁶ relative on a 10¹²-GPa pathological entry, within
  its bound). The initial refutation at an output-only floor was a floor
  mis-declaration, diagnosed and amended per the v2 procedure.
- **P2 MET:** the 156-pair open sweep certifies NOTHING — every pair a
  scout-level conjecture. No unregistered certificate; heavy abstention as
  designed.
- **P3 MET:** zero confident-wrong.

## Campaign closure

The Materials Project campaign (C0–C2, `PROPOSAL_MATERIALS.md`) is
complete. Five certified pipeline identities — density (with the decoded
2005-vintage mass table and the atomic mass constant to 9 digits), the two
VRH averages, volume-per-atom, and the verify-pinned Poisson and
anisotropy forms — at bounds up to **α ≤ 10⁻⁵⁰⁷⁰**; one upstream
data-quality report ([emmet#1499](https://github.com/materialsproject/emmet/issues/1499)
— eight materials with Reuss > Voigt plus one implausible magnitude; the
maintainer corrected our "unflagged" reading, see CASE_STUDY_MATERIALS_C1.md
P2); the banded elastic conjectures
all inside their registered bands; and zero confident-wrong across three
stages. The registered homogeneity contrast with PSCompPars held: a
computed database is internally consistent in a way a literature composite
is not — and where it isn't, the rows are enumerable by ID. Next campaign
per the proposal: Open Catalyst, with the Abild-Pedersen rational-slope
registration as its centerpiece.
