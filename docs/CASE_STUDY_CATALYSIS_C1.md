# Catalysis C1 — the rational-slope registration

**Registered 2026-07-28, bands frozen BEFORE any slope is computed** (the
C0 firewall has held: no energy–energy regression exists anywhere in the
program's artifacts for this data). Governing proposal:
`docs/PROPOSAL_CATALYSIS.md`. Data: the rebuilt table
`experiments/catalysis/data/mamun_rebuilt_energies.csv` (C0-validated).

## The hypothesis under test

Bond-counting (Abild-Pedersen 2007) predicts the scaling slope of a
hydrogenated adsorbate against its central atom is the exact rational
γ = (x_max − x)/x_max:

| pair | γ (exact) |
|---|---|
| CH* vs C* | **3/4** |
| CH₂* vs C* | **1/2** |
| CH₃* vs C* | **1/4** |
| NH* vs N* | **2/3** |
| OH* vs O* | **1/2** |
| SH* vs S* | **1/2** |

## Registered predictions

- **P1 (bands).** Each pair's OLS slope (most-stable-site ΔE*, all joinable
  surfaces) lands in **γ ± 0.10**. Six independent cells, scored
  separately.
- **P2 (exact-rational consistency).** |slope − γ| ≤ 2·SE(slope) for at
  least **4 of 6** pairs — the theory's exact rationals are statistically
  admissible, not merely near.
- **P3 (no certificates — registered confident-wrong bait).** `recover` on
  every pair (σ declared 0.22 eV from the C0 site-multiplicity spread)
  must abstain: DFT scatter is real, and an exact-law certificate on this
  data would be a confident-wrong. Zero certificates in C1.
- **P4 (intercepts are chemistry, not universal).** The six intercepts
  differ from each other by > 0.5 eV in at least one pair-of-pairs —
  scaling is per-species, not one law.
- **P5.** Zero confident-wrong.

## Results

(after the run)
