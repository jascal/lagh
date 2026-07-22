# R-noise study — the zero-wrong invariant under measurement noise

**Run 2026-07-21.** Relative Gaussian noise injected at three SNR levels, correctly-declared sigma,
direct `discover` (max_tier=5) over all 108 NewtonBench-dev cells; each recovered law scored against
the **clean** oracle. Settles the R-noise readiness gate (`STRATEGY.md`).

## Result

| noise (σ_rel) | ≈SNR | recovered (exact) | abstained | certified ≠ exact truth | wrong **beyond** noise | wrong > 3× noise |
|---|---|---|---|---|---|---|
| 0.001 | 60 dB | 47 | 58 | 1 | 1 | **0** |
| 0.01 | 40 dB | 18 | 52 | 36 | 1 | **0** |
| 0.1 | 20 dB | 0 | 64 | 42 | 0 | **0** |

## Verdict — two properties, only one survives noise

1. **No gross fabrication (holds).** At every level, **zero** certified laws are wrong by more than 3×
   the injected noise, and at most one exceeds 1×. Every certificate is accurate to within the
   measurement noise. lagh does **not** produce a confidently-wrong-beyond-the-data answer under noise.

2. **Exact-symbolic recovery (clean-data only).** The "recover the *exact* closed form or abstain"
   guarantee does **not** survive noise: the exact-rational snap locks onto a *nearby* rational that
   fits the noisy data within the (correctly widened) ε=4σ, and lagh certifies it as exact instead of
   abstaining — 1 / 36 / 42 such cases at 0.1 / 1 / 10 %. Exact-recovery coverage collapses 47→18→0.

**Mechanism.** Coherence abstains on *functional under-determination* (rival materially-different laws).
Noise produces *parametric mis-determination* — one law with slightly-wrong coefficients/exponents and
no rival — which the abstention machinery does not catch. The widened C3 exponent grid (CAP-A) gives
the snap marginally more wrong rungs to land on under noise.

## Consequence for the headline claim

The standing **"zero confident-wrong across ~180 tasks" is a CLEAN-DATA result** (deterministic
oracles, machine-precision ε). It must be stated as such. Under noise the honest claim weakens from
*"exact law or abstain"* to *"a law accurate to within the stated noise, or abstain"* — which the study
supports (0 gross wrongs), but it is a different, weaker guarantee.

## R-noise gate — SETTLED by restriction

- For the **exact-symbolic** claim: register the blind read **clean / low-noise only** (σ_rel ≤ ~1e-3;
  even there 1 near-miss slips). Restrict the reserved-benchmark cells accordingly in the registration.
- The **numerical** claim ("certified accurate to within declared noise, else abstain") is defensible at
  all tested levels (no wrong-beyond-3×-noise), and is the honest claim if the blind benchmark is noisy.

## The fix to extend *exact-or-abstain* into the noisy regime (deferred, substantial)

A **parametric-uncertainty abstention**: propagate the fit's noise-induced covariance to each snapped
coefficient/exponent; if a *different* rational lies within that uncertainty band and also certifies,
**abstain** (the exact value is not pinned). Trades coverage for restored exactness under noise;
connects to the significance / `α = |H|·qʰ` counting direction (`DIRECTION_SIGNIFICANCE.md`). Until
built, noise is handled by restriction, not by the instrument.

---

## CORRECTED VERDICT (2026-07-22) — parametric gate active + STRUCTURAL scoring

The table above scored "confident-wrong" as clean rel-err > 1e-3. That **conflated two different
things**: a genuine STRUCTURAL error (wrong exponents/monomials) and an inherently noise-limited
COEFFICIENT (correct structure, a float constant off by ~the noise — which you *cannot* recover
exactly from noisy data). Re-scored by symbolic STRUCTURE (strip every additive term's multiplicative
coefficient, keep exponents & atoms), with `certify.pinned()` active:

| noise (σ_rel) | ≈SNR | recovered (correct structure) | abstained | STRUCTURAL confident-wrong |
|---|---|---|---|---|
| 0.001 | 60 dB | 46 | 19 | **0 / 65** |
| 0.01 | 40 dB | 41 | 24 | **0 / 65** |
| 0.1 | 20 dB | 30 | 32 | **3 / 65** |

(65 = cells where clean data recovers, i.e. there is a truth to match; 43 have no clean truth.)

**The "36 CW at 1%" was largely a SCORING artifact.** Every flagged case at ≤1% noise was verified
cell-by-cell to be coefficient-only — clean's exact rational vs noisy's float of the *same value*,
*identical* exponents/monomials. The instrument's real behaviour up to 1% noise is
**recover-correct-structure-or-abstain, with ZERO structural fabrication.** Coefficients are
noise-limited floats (inherent); the parametric gate pins the rational STRUCTURE.

**Residual at heavy (10%) noise: 3/65 gross FORM-substitutions** the parametric gate does *not* catch,
because they are CROSS-FORM, not nearby-rational:
- malus medium_v0: `(2sinθ+cosθ)²` → `4.6·x₀·√x₁`
- hooke hard_v0: polynomial → `exp(7.44·x^{1/6})`
- coulomb medium_v0: poly/`r²` → fractional-power product

The gate tests neighbour-*rationals of the winner*; at 10% noise a *different functional form* wins
outright (only one form certifies, so coherence sees no rival either). Catching these needs a distinct
mechanism — cross-tier rivalry under the noise-loosened ε, or a form-stability check — **not** the
parametric gate. Registered as the open heavy-noise failure mode.

**Corrected R-noise gate.** The **exact-STRUCTURE-or-abstain** guarantee holds to ~1% noise (0/65),
degrades gracefully (correct-structure 46→41→30, abstain 19→24→32), and breaks only for 3/65 cross-form
cases at 10%. The blind read can be registered **for SNR ≥ ~40 dB (≤1% noise) with the structural
guarantee intact**, tightening to clean-only only if the exact *coefficient* (not just structure) must
be certified. The parametric gate did its job (nearby-rational mis-snaps → abstain); the remaining
exposure is a heavier-noise, cross-form failure mode.
