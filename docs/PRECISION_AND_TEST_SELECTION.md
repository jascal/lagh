# Precision and pilot test selection

**Empirical result:** the registered covariance GLS estimator does not support
a universal 8% error budget. **Empirical result:** pilot-only test selection
is feasible but bought no certified reach on the registered Itô ladder.
The tweezers campaign remains closed; no raw record was reopened.

## Timescale result: outcome (b)

The legacy log-ACF window chooses which correlations to fit using their sign
and magnitude. The opt-in `covariance_timescale` instead fits A exp(-theta*t)+C
on fixed lags8–199, weighted by an empirical block covariance with a frozen
20% diagonal regularizer. C is a nuisance covariance offset, not an identified
contaminant amplitude. Historical `acf_timescale` defaults stay unchanged.

**Empirical.** All 12 regimes use the existing bead generator, with the detector
and anti-alias filter applied before sampling in filtered regimes. Unfiltered
controls use its same exact OU generator and physical truth. Contamination is
a 7Hz sine with twice the observed thermal variance. Each score includes full
and half estimates; an unreadable estimate counts as infinite error.

The initial pilot-size error is preserved in `timescale_initial.json`: eight
calibration records cannot give a finite rank-corrected 90% bound. The recorded
amendment increased calibration to12 (seeds20–31) and used fresh evaluation
seeds60–63, without changing estimator/window/regularization. No initial
validation seed was relabeled as fresh. Raw C3/C4 ratios are not recomputed.

**Empirical calibration rule.** For each calibration record define
E=max(|theta_hat/theta_true−1| for full and both halves,
max(theta_half)/min(theta_half)−1). With m=12, select sorted E at rank
ceil(.9*(m+1))=12. This derives a required tolerance from observed estimator
error rather than tuning a declared 8% threshold. Under exchangeability this
is the usual finite-rank calibration rule; applicability outside these
synthetic regimes is **open**. It is not a joint three-ratio uncertainty model.

| n | Filtered | Slow sine | GLS required tolerance | Evaluation max truth error | Evaluation max half disagreement | Evaluation max GLS/PSD disagreement |
|---:|:---:|:---:|---:|---:|---:|---:|
| 32768 | False | False | 245.0% | 35.3% | 74.7% | 8.3% |
| 32768 | False | True | 746.7% | 44.3% | 98.1% | 24.6% |
| 32768 | True | False | 287.1% | 42.2% | 57.0% | 13.4% |
| 32768 | True | True | 156.9% | 38.6% | 88.9% | 6.9% |
| 131072 | False | False | 46.8% | 13.1% | 15.4% | 5.0% |
| 131072 | False | True | 34.8% | 20.0% | 21.2% | 6.3% |
| 131072 | True | False | 35.7% | 13.1% | 20.1% | 11.8% |
| 131072 | True | True | 38.2% | 13.0% | 14.1% | 9.4% |
| 524288 | False | False | 9.3% | 5.9% | 3.7% | 1.8% |
| 524288 | False | True | 12.0% | 15.0% | 14.0% | 6.8% |
| 524288 | True | False | 13.5% | 3.7% | 4.3% | 3.3% |
| 524288 | True | True | 26.3% | 19.5% | 16.3% | 9.8% |

**Empirical comparison.** All three estimators, their per-seed full/half
readings, and errors are retained in `timescale.json`; poor methods/regimes
are not dropped. The GLS full-rate error is sometimes small while split scatter
is too large. No blanket claim that GLS improves the legacy estimator is made.
The registered broad 8% claim is withdrawn; this is outcome (b).

**Open.** ATTRIBUTION_TOL=.08 remains a historical compatibility threshold,
not a validated uncertainty guarantee. The derived table prices theta precision
only: using it as one scalar for theta, diffusion and variance would silently
weaken two independently uncalibrated residuals. No physical-attribution reach
is credited, and no old C3/C4 ratios are promoted. R-C3's unknown stage unit,
detector whiteness and cross-record rate causality require additional data.
The estimator exercise settles none of them.

## Restricted score-function experiment

**Open derivation / scope.** Write h=phi*w. Ignoring deterministic error,
Cauchy–Schwarz bounds (int h*a)^2 by int h²*b² times int a²/b².
Equality requires h proportional to a/b² where b²>0. This does not optimize
the deterministic term; a fixed state function w(x) generally cannot cancel a
time bump phi(t). At vanishing diffusion or outside the chosen function's
domain the formal score is not a usable test. OU with constant b already has
w proportional to x; GBM suggests w proportional to 1/x, hence f=log(x).
No unconstrained-optimality theorem is claimed for the implemented selector.

**Empirical protocol.** `choose_test` takes pilot trajectories only. It fits
drift on (1,x), nonnegative diffusion on (1,x²), then scores x, x²/2, log(x)
(the latter only for positive pilot paths) using a time-averaged martingale
proxy. The returned choice is immutable and carries a pilot hash. Certification
paths are generated after selection, from independent seeds; their hashes are
recorded. There is no validation-based fallback or retuning. The same vocabulary,
Itô rows, delta=.05 and existing coefficient-dependent bands certify both arms.
No pilot-fitted diffusion is inserted into the certification band.

| Case (4 seeds each) | Selected function | Baseline certificates | Selected certificates | Truth exclusions (either arm) |
|---|---|---:|---:|---:|
| OU | x**2/2 | 0/4 | 0/4 | 0 |
| GBM | log(x) | 0/4 | 0/4 | 0 |
| null | log(x) | 0/4 | 0/4 | 0 |

**Empirical.** All24 arms covered the true law's checked residuals and no
reported coefficient interval excluded truth. This is a small controlled
sample, not evidence of a new general soundness guarantee. Every interval and
band diagnostic is retained, including losses. All drift-x intervals include
zero: the selection buys neither certification nor resolved-mode reach.
For GBM, selected/baseline drift-x interval width ranges from 0.691 to 1.066 (median 0.963).
For null, selected/baseline drift-x interval width ranges from 0.027 to 0.925 (median 0.693).

**Open / failed prediction.** A measured certified reach gain is not delivered
by this restricted plug-in score. The broad gain claim is withdrawn in the
registration, not rescued by tuning windows or selecting on certification
residuals. Selection remains opt-in. Deterministic algebraic reach cannot gain
from a state-path test; its fresh unchanged-consumer run is reported separately
in `test_selection_reach.json`. Neither arm silently narrows the vocabulary to
the true drift or replaces the measured QV with pilot estimates.

## Residual evidence and reproduction

**Empirical.** `check()` now returns a verdict mapping with aligned `residual`,
`epsilon`, `predicted`, and `row_indices` attributes. Its bounded `measurement()`
is built from that evaluation. MCP no longer recomputes the full-data candidate;
engine terminal no-law refusals, passive full-data refusals, and acquisition
held-out failures attach that measurement to the certificate. It names the
candidate and domain, never a cause. Other kinds of refusal need not have
residual evidence. Invalid arithmetic is uncovered, not certified.

Use one process per test file. Reproduction commands (from repository root):

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_falsifiability
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_falsifiability_seeds
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_timescale
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_test_selection
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_selection_reach
```

C4's JSON contract is now explicit in its original registration: identical
structure/nonnumeric values, finite floats at rtol1e-12/atol1e-15. Compare with
`python -m experiments.compare_artifacts original.json rerun.json`. A 1e-5
perturbation fails and a one-ULP perturbation passes, as recorded in
`selection_artifact_witnesses.json`. This tolerates summation order, not model
error, and requires no renewed C4 analysis.

Run `experiments.run_timescale --initial` for the retained initial run.
`experiments.run_falsifiability --baseline` loads the pinned historical core
into a temporary directory and reproduces its eight failures without changing
the working tree. `experiments.run_selection_witnesses` and
`experiments.run_residual_consumers` reproduce the negative-input and consumer
artifacts; the latter explicitly reconstructs the rejected aliasing design for
the mutation baseline. Python/package versions and thread declarations are in
`experiments/results/research_environment.json`. Calibration seeds are reused
across regimes as common random numbers; no simultaneous90% coverage over all
regimes, or over theta/diffusion/variance jointly, is asserted.

## Final validation

**Empirical:** the complete unchanged deterministic reach audit certified35/36,
the same35/36 as the historical artifact, with no verdict changes. This is no
reach gain and no new out-of-domain soundness claim. Runtime telemetry is
excluded from the reproducible artifact.

205 tests passed, each file in its own process: evidence9, refusal7,
public surface55, instrument C4 25, instrument21, selection2, MCP17, passive5,
significance3, boxsearch3, Itô42, and the selected lagh core group16.
Five slow lagh tests were deselected; a full-suite result is not claimed.
Existing empty-array/overflow warnings remain. F/E9 lint and diff whitespace
checks passed. The initial/final timescale artifacts and selection artifact
were regenerated after the runner refactor and were byte-identical.
