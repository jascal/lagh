# Direction: the passive-data regime (fixed datasets, no oracle)

Status: **registered 2026-07-23, BEFORE the sweep** (predictions below are scored, not
tuned toward). Companion to `STRATEGY.md` — this is the regime the reserved blind
benchmark actually is.

## Why this exists

Every symbolic-regression benchmark in the reserved-blind candidate set
(LLM-SRBench, SRSD-Feynman, SRBench) is **fixed-dataset**: the benchmark hands you
`(X, y)` sampled where *it* chose, and there is no oracle to query. Everything in
`lagh/acquisition.py` — adaptive ranging, replicates-based sigma, multi-objective
active queries, the box-search ladder — assumes a queryable oracle and therefore
**does not exist in the regime that will be scored**. All NewtonBench-dev numbers to
date (67/108, 0 CW) were measured in the active regime. The competitive question is
what survives without it.

## The passive entry point

`lagh/passive.py : discover_passive(X, y, sigma=0.0, ...)`:

1. Drop non-finite rows (a passive dataset may contain saturated/overflow cells).
2. **K deterministic re-splits** (default 3) of the same data into fit/select/certify
   (60/20/20, seeded permutations) — the passive substitute for the active loop's
   per-round re-splits. A split that happens to put all the high-signal points in
   `fit` is no longer fatal.
3. **Full-data exhaustive gate**: a law that certifies on the cert split must ALSO
   pass `check()` on **every point in the dataset** at the assembled epsilon. This
   is what keeps K re-splits sound: re-splitting multiplies hypothesis exposure ~K×,
   but the final gate is exhaustive on all n points — a wrong law would have to fit
   every observation to machine epsilon. Certification is never *granted* by a
   re-split, only *re-attempted*; the certifier itself is byte-identical.
4. No replicates exist, so sigma is **declared-only** in passive mode (0 for clean
   benchmarks; the benchmark's stated noise level otherwise). Estimating sigma from
   duplicated X rows, if a benchmark provides them, is future work — noted, not built.

What is **lost** relative to active, by construction: adaptive ranging (dead-signal
box contraction), the box-search ladder on abstain, active disambiguation queries,
replicate-estimated sigma. None has a passive equivalent; where the sweep shows one
was load-bearing, that is a finding about the regime, not a bug.

## The dev proxy

`experiments/run_newtonbench_passive.py`: for each of the 108 NewtonBench-dev cells,
generate ONE fixed dataset per sampling law — **n=250, seeded per cell** — over the
module's declared box, then `discover_passive` on it, scored by the same `dense_ok`
grid as the active sweep. Two sampling laws, both swept:

- `loguniform` — matches the instrument's own sampling assumption (isolates the
  *acquisition* variable: same distribution family as active, minus the feedback).
- `uniform` — what a benchmark most likely hands out (also stresses the low-decade
  coverage that log sampling gives for free).

## Registered predictions (to score, in `## Results` below, unchanged)

- **P1 (coverage):** passive/loguniform recovers **≥ 60 of the active 67** correct
  cells. The deficit candidates are cells where ranging/box-search was load-bearing
  (suspects: m5 decay far-tail, m10 BE scale, any cell the ladder rescued).
- **P2 (zero-wrong):** **0 confident-wrong in both variants.** The certifier is
  unchanged and the full-data gate only removes certifications.
- **P3 (sampling):** uniform recovers **≤** loguniform (worse low-decade coverage for
  power/log laws), with the gap concentrated in wide-box modules.
- **P4 (no false unlocks):** passive recovers **no cell the active sweep missed** —
  there is no mechanism for passive to see more.

## Amendment (2026-07-23, before any results were read)

The first sweep attempt crashed mid-run on `m10_be/easy/v1` (a dataset whose every
output exceeds the float ceiling -> zero finite rows -> empty splits; guard added).
Between the attempt and the re-run, three capabilities landed (the exact-coefficient
gate, CAP-C, CAP-G angular). Both regimes are therefore re-measured **on the same
instrument version**: the passive sweep AND a fresh active sweep, concurrently.
P1 is restated instrument-relative: *passive recovers ≥ 90% of what the same-day
active sweep recovers.* P2-P4 unchanged. The partial first-attempt log (66 cells,
all pre-gate) was discarded unread beyond the crash line.

## Results

*(to be filled by the sweep; predictions above stay frozen)*
