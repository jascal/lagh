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

## Results

(after the run)
