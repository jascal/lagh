# R-gap closure: every curriculum gap bounded with a one-line cause

**Rewritten 2026-07-23** against the post-capability sweep (`experiments/results/
newtonbench_all.jsonl`: **87/108 recovered, 0 confident-wrong**). Supersedes the
2026-07-22 version (66/108 baseline) — notably its "UNITS/IRRATIONAL-IN-PRACTICE"
bucket was **wrong**: the 9 degree-unit snell cells were never effectively
out-of-class; symbolic `pi` carries the unit exactly and CAP-G recovers all 9.

## Buckets (21 abstains, zero unexplained)

| bucket | count | meaning |
|---|---:|---|
| **OUT-OF-CLASS (`^e` wedge)** | 11 | irrational exponent; no exact-rational form exists — correct-abstain, permanent |
| **ORACLE-PRECISION-LIMITED** | 3 | the benchmark's own float64 `exp(u)−1` at u≤1e-10 carries ≥1e-5 relative error; exactness is unmeasurable there — correct-abstain, permanent for THIS oracle |
| **OPEN-CAPABILITY (scoped)** | 6 | exact form exists; the needed grammar extension is named below, unbuilt |
| **DECLARED OUT-OF-SCOPE** | 1 | m2 hard_v0 `(I1+I2)^1.5` — fractional power of a sum; excluded by registered decision (NEWTONBENCH_GAP_PLAN.md) |

**Achievable-ceiling accounting:** 97 grammar-representable (88 algebraic + 9
angular) − 3 oracle-precision-limited = **94 measurable; 87/94 = 93% recovered.**
The addressable frontier is the 6 OPEN-CAPABILITY cells.

## Per-cell bounding (all 21)

### OUT-OF-CLASS — the `^e` wedge (11): correct-abstain by design
`m1_coulomb hard v0/v2` · `m5_decay hard v0/v1/v2` · `m7_malus medium v2, hard v1/v2`
· `m11_heat hard v0/v1/v2` — each ground truth carries an irrational exponent
(`^e`); an exact closed form does not exist; any "recovery" would be a
confident-wrong. These 11 are the zero-wrong line, not a gap.

### ORACLE-PRECISION-LIMITED (3): correct-abstain by measurement
| cell | cause |
|---|---|
| `m10_be easy v1` | u=C·√ω/T ∈ [1e-16, 4e-15]: float64 `exp(u)−1` in the oracle is ≥50% noise; y also saturates the 2^52 ceiling |
| `m10_be easy v2` | same mechanism, u ∈ 1e-16 range; most of the box returns NaN at the ceiling |
| `m10_be medium v2` | same mechanism (u=C·√ω·T^2.3 small over the box) |

### OPEN-CAPABILITY (6): exact form exists, extension named
| cell | cause (one line) |
|---|---|
| `m5_decay easy v2` | `N0·e^{−(λt)^{3/2}}` — exp of a CROSS-monomial inner; C4's inner pool has no x_i^p·x_j^q monomials (scoped extension: fractional cross-monomial inners) |
| `m5_decay medium v2` | `N0^{9/5}·e^{−(λt)^{3/2}}` — same inner, plus fractional prefactor power |
| `m6_underdamped medium v2` | `k·m^{−13/10} − b^2/4m^2` (verified against laws.py) — denom-10 fractional summand at dim-3; CAP-A2's registered bound is dim≤2 (scoped: lift to dim≤3 for negative denom-10 exponents) |
| `m6_underdamped hard v0` | `(k/m − b/2m^2)^{3/2}` (verified) — rational outer power 3/2; CAP-D's transform pair covers only {sqrt, square} (scoped: y^{p/q} outer-transform family for small p/q) |
| `m6_underdamped hard v2` | `k·m^{−13/10} − (b/2m)^{7/10}` (verified) — composite-BASE fractional power `(b/2m)^{7/10}`; exactly representable in sympy but the feature grammar has no (x_i/x_j)^{p/q} composite bases (scoped) |
| `m10_be hard v2` | `1/(−ln(...)−1)` — rational-of-LOG inner (a `log1p`-analogue transform with log inner; registered CAP-E stretch, unbuilt) |

### DECLARED OUT-OF-SCOPE (1)
| cell | cause |
|---|---|
| `m2_magnetic hard v0` | `(I1+I2)^{3/2}` — fractional power of a SUM; composite-power capability declined by registered decision (feature-inflation risk > single-cell gain) |

## R-gap verdict

Every non-recovered cell is bounded with a cause; 14 of 21 are *correct-abstain*
(wedge + precision), 6 are scoped-and-named extensions, 1 is a registered
exclusion. **R-gap is MET** on this dev set as "closed or explicitly bounded".
