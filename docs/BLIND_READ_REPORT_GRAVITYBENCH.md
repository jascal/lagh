# BLIND READ REPORT — Gravity-Bench-v1 (read executed 2026-07-27)

One read per variant, per `BLIND_READ_REGISTRATION_GRAVITYBENCH.md`. Logged
crash-fixes before any scoring information existed: HF-dataset loader (the repo
jsonl ships placeholder CSVs; all 412 first-pass solves crashed), `yrAUMsun`
unit marker, full-variant uniform 2400-row subsample (bounded runtime; uniform
cadence preserved). Scoring is mechanical per-instance thresholds from the
dataset itself.

## Scorecard vs the frozen baselines

| | lagh astronomer | frozen SOTA (o4-mini-high) |
|---|---|---|
| **Budgeted (100 observations)** | **94.66% (195/206)** | 49% |
| Full observations | 63.59% (131/206) | 74% |

**The headline: on the variant the benchmark was built around — plan your own
observations under a budget — the deterministic lagh astronomer scores 94.7%,
nearly DOUBLING the best published LLM-agent result, with zero LLM calls.**
The observation-planning thesis (active acquisition is lagh's home game) is
confirmed on a sealed benchmark: coarse scan → angle-based period → phase-dense
pass → strict-cadence acceleration triplets → digital-twin answers, all
deterministic, every quantity traceable to a fit. Modified-gravity exponents,
drag coefficients, Kepler booleans, sweep rates, and all COM-frame extremes
score at or near ceiling in the budgeted regime; one instance failed on an SVD
non-convergence (counted wrong, as registered).

## The full-observation shortfall (63.6% vs 74%), honestly

Below baseline, and the failure pattern is specific: min_* extremes
(acceleration/angular-velocity/momentum minima) and drag score ~0 in the full
variant while scoring near-perfectly in the budgeted one. Mechanism (post-read
diagnosis, no re-run): the full variant has NO planner, so the twin's epoch
state comes from the uniform-cadence fallback rather than the planner's strict
P/400 triplets — a subtly worse epoch propagates into the twin's extreme-value
minima while leaving broad maxima intact. The budgeted agent literally beats
the full-data agent BECAUSE its planner engineers better differentiation
points: the benchmark's own lesson, demonstrated by the instrument. A
triplet-refinement pass for dense tables is the obvious fix; per one-shot
discipline the 63.59% stands as read.

## Standing

- First blind-read WIN of the program (the LLM-SRBench loss is the other
  bookend; both reported as registered).
- The claim: *a fully deterministic, certificate-bearing instrument nearly
  doubles LLM-agent SOTA on budgeted scientific observation planning* — no
  tokens, reproducible, every answer traceable to a fitted, validated model.
- Gravity-Bench is now spent as a blind set (DEV from here), and the full-obs
  epoch-refinement fix is registered as dev work if the benchmark is revisited.
