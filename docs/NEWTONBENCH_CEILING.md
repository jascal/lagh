# NewtonBench: the achievable ceiling and the abstention wedge

**Established 2026-07-21 by full-matrix forensics** (108 cells = 12 modules × 3 difficulties × 3
versions), cross-referencing each cell's ground-truth law (`modules/m*/laws.py` docstrings) against
the frozen-instrument sweep (fractional-power grammar, `nb_frac` run). Zero confident-wrong across
all 108. This document supersedes the raw "55/108" headline.

## 1. NewtonBench is not 108 exact-recoverable laws

The difficulty/version axes do **not** vary noise on one physical law — they substitute a *different
synthetic law* per cell, many non-physical (the "malus" cells are `I₀·(sinθ+cosθ)²`, `I₀·tan²θ`,
`I₀·(cosθ/sinθ)^e` — none is Malus's actual `I₀cos²θ`). Classifying every cell by what an
**exact-rational** instrument can even represent:

| class | cells | can an exact-rational instrument recover it? |
|---|---:|---|
| **algebraic / elementary** | **88** | yes — polynomial / rational / power / trig-of-input, exact |
| **irrational-exponent** (`x^e`, `x^{np.exp(1)}`, `(…)^π`) | **11** | **no, by design** — the true form is irrational; any exact closed form is a lie |
| **inverse-trig** (`acos/asin/atan(...)`, all 9 of m4_snell) | **9** | only with an inverse-trig tier (exact, but a class we don't yet emit) |

**The achievable ceiling for the current grammar is 88, not 108.** With an inverse-trig tier it is 97.
The 11 irrational-exponent cells are **permanently** out of an exact instrument and *must* be abstained.

## 2. The scorecard against the real denominator

| class | recovered | correct outcome | note |
|---|---|---|---|
| algebraic (88) | **55** | recover | **33 addressable gaps** — trig products/ratios, rational-of-exp, transcendental products |
| irrational (11) | 0 | **abstain** | **11/11 abstained — the zero-wrong line.** A "pass" here is a confident-wrong exact claim |
| inverse-trig (9) | 0 | abstain (today) | recoverable only by adding an `arcsin/atan` tier (+9 ceiling) |
| **all 108** | **55** | — | **0 confident-wrong** |

**55/88 = 62.5 % of the achievable ceiling**, not 51 % of a mis-stated 108.

## 3. The wedge — why abstention is the benchmark claim, not a weakness

NewtonBench's own hard tier contains **11 cells whose ground truth is irrational** (`^e`). The
benchmark scores submissions by RMSLE + symbolic-equivalence. An LLM (or any regressor) that emits a
rational approximation to `(cosθ/sinθ)^e` can **pass on RMSLE while being exactly wrong** — it claims
a closed form that does not exist. lagh is the only instrument in the comparison that **abstains on
exactly those 11 cells and never on an algebraic one it can certify.** The composite guarantee
(`LLM + lagh ≥ LLM alone`) is sharpest here: on the irrational cells lagh subtracts a confident-wrong
answer and adds a machine-readable reason.

**The headline claim to make is therefore not "highest coverage" but:**
> *On NewtonBench, lagh recovers algebraic laws with an exact certificate, abstains with a reason on
> every irrational-exponent cell where an exact form does not exist, and produces zero confident-wrong
> answers over all 108 cells — a guarantee no RMSLE-scored baseline can state.*

## 4. Addressable algebraic gaps (33), by cluster — the coverage work that remains

| module | gaps | missing capability |
|---|---:|---|
| m6_underdamped_harmonic | 7 | transcendental **products** `e^{−γt}·cos(ωt)` (decay × oscillation) |
| m10_be_distribution | 6 | **rational-of-exponential** `1/(e^x − 1)` (C2 rational × C4 inner-exp) |
| m7_malus_law | 5 | **trig products / ratios** `tan²`, `sin²/cos³`, `(sinθ+cosθ)²` (sin·cos cross, sec) |
| m9_hooke_law | 4 | (inspect — likely product / higher-degree) |
| m8_sound_speed / m1_coulomb | 3 / 3 | hard-tier algebraic forms |
| m5_radioactive_decay / m0_gravity | 2 / 2 | — |
| m2_magnetic_force | 1 | — |

**Note the base laws already recover** (the current grammar gets `I₀cos²θ` via the half-angle
`x₀cos(2x₁)/2 + x₀/2`); the gaps are the *manglings* NewtonBench substitutes, which need
trig-product/ratio, rational-of-exp, and transcendental-product features. These are principled,
per-class grammar generalizations (like the fractional-power add), **not** per-target curve-fits —
and each must be registered with its predicted unlocked cells before the re-sweep.

## 5. Contamination note (unchanged)

NewtonBench stays **DEV-only** (`STRATEGY.md`): it is public and its laws are readable (as done here).
The zero-wrong invariant and the ceiling analysis are what transfer to a blind benchmark; the specific
33-gap grammar work is dev-set tuning and must be validated on held-out laws before any headline.
