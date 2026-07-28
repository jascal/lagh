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

## Results (2026-07-28, `experiments/results/catalysis_c1.json`; rebuilt
table, 239–255 surfaces per pair; slopes computed for the first time AFTER
the bands froze)

| pair | γ | slope ± SE | in band ±0.10 | 2·SE-consistent |
|---|---|---|---|---|
| CH vs C | 3/4 | 0.832 ± 0.012 | ✓ | ✗ |
| CH₂ vs C | 1/2 | **0.531 ± 0.016** | ✓ | **✓** |
| CH₃ vs C | 1/4 | **0.225 ± 0.014** | ✓ | **✓** |
| NH vs N | 2/3 | 0.758 ± 0.015 | ✓ (barely) | ✗ |
| OH vs O | 1/2 | 0.715 ± 0.030 | ✗ | ✗ |
| SH vs S | 1/2 | 0.897 ± 0.026 | ✗ | ✗ |

- **P1 half-MET: 4/6 in band** (registered per-cell; OH and SH miss high).
- **P2 MISSED: 2/6** exact-rational-consistent vs the registered ≥ 4/6. The
  bond-counting rationals are statistically admissible only for the
  CH₂/CH₃ rungs; CH and NH run ~0.08 high of their rationals at small SE;
  OH (+0.21) and SH (+0.40) deviate decisively.
- **P3 MET: zero certificates** — all six `recover` calls structurally
  abstain at the declared 0.22 eV scatter. The confident-wrong bait was
  not taken on the campaign's own centerpiece.
- **P4 trivially met and honestly discounted:** intercepts span 129 eV, but
  in ΔE* space the intercepts carry the per-species gas-reference
  constants, so the registered criterion is confounded and claims nothing.
- **P5 MET:** zero confident-wrong.

**The sentence this stage earns — the honest version of the headline we
wanted:** on 25 GiB of rebuilt DFT with bands frozen in git before any
slope existed, bond counting's exact rationals survive cleanly on the
carbon ladder (CH₂ at 0.531 ± 0.016 against 1/2; CH₃ at 0.225 ± 0.014
against 1/4), bend for CH and NH, and break for OH and SH on bimetallic
alloys — deviations all toward stronger-than-predicted coupling, at
most-stable-site processing where adsorbate-induced reconstruction is
known to defeat simple valency counting. The instrument certified nothing,
banded everything, and reported the theory's misses with the same care as
its hits.
