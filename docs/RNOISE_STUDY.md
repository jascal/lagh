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
