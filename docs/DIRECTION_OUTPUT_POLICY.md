# Direction: two-track output policy for accuracy-scored benchmarks

Status: **pre-registered 2026-07-23**, before any blind read. Companion to
`STRATEGY.md` (which governs *which* benchmark and *when*); this governs *what is
submitted per problem*. Registered now so the policy cannot be tuned after seeing
blind results.

## The scoring reality this answers

Every reserved-blind candidate (LLM-SRBench, SRSD-Feynman, SRBench) scores
**symbolic/numeric accuracy only**: a wrong equation scores the same as no equation.
Abstention is unrewarded — the zero-confident-wrong invariant, as a *submission
policy*, leaves points on the table and cannot win. (REFUTE-style calibration
benchmarks are the exception; there, abstention IS the metric and track B below is
switched off.)

The resolution is to split what lagh *submits* from what lagh *claims*:

## The two tracks

Per problem, the submission is the first track that produces an expression:

- **Track A — certified.** `discover`/`discover_passive` certifies: submit the exact
  law. Tagged `proved`, certificate attached. This is the only track the zero-wrong
  claim covers.
- **Track B — labeled conjecture.** lagh abstains: submit the best *available*
  expression anyway, explicitly tagged `empirical` (never `proved`):
  1. the `fit` scout's conjecture (`lagh/mcp` `fit` — already built, diagnosis
     attached), else
  2. the engine's best non-certifying candidate (the `Certificate.law` field on a
     structural abstain — the best law that certified but was not coherent-unique), else
  3. the `characterize()` power-law snap when `class == power-law`, else
  4. nothing (only when no probe produced any expression — expected rare).
- On a calibration-scored benchmark, track B is disabled and the abstention +
  machine-readable reason is the submission itself.

## The claims, scoped

- **Headline:** "X% symbolic accuracy, of which Y points came with a machine-checked
  exact certificate; **zero wrong answers within the certified subset**."
- The zero-wrong invariant is a **track-A claim**. A track-B conjecture that scores
  wrong does not break it — it was labeled a conjecture at submission time. What
  WOULD break the program is a certified-track wrong answer; that remains the
  revert-the-capability tripwire.
- Track B never upgrades: no post-hoc "the conjecture was right, call it certified".
  Certification happens at discovery time or not at all.

## Why this is honest

The label travels with every answer. This is exactly the `proved`/`empirical`/`open`
tag discipline applied to benchmark submissions: the benchmark's accuracy metric sees
both tracks; anyone auditing the certified claim sees the partition. The alternative
(submit only certified answers) doesn't strengthen the claim — it just hides the
instrument's conjecture quality, which the dev record already measures separately.

## Pre-registered accounting for the blind read

At the blind read, report: total score; per-track counts (A submitted / A correct,
B submitted / B correct); the certified-subset wrong count (**must be 0**); and the
abstain-with-no-conjecture count. All four numbers are in the registration BEFORE
scoring, so the partition cannot be redrawn afterward.
