# Refusal-guided acquisition: scored result

**Empirical — clean-data false exactness is the central finding.** P2 certified
structurally wrong approximants on 8/16 noiseless initial 400-row draws
(`sigma=0`, machine-floor bands), before acquisition could act. P1 also retained
eight wrong arm outputs from two conditions sharing an initial draw, even after
fresh original-box checks. Local residual coverage does not establish the exact
generating form, including on clean data. The corrected acquisition entry guard
excludes initial certificates; it does not repair the ordinary engine.

**Open — the escalation rule.** The engine stops at the first tier with a
non-empty certifying set, whether it returns a winner or a structural refusal.
Later tiers therefore cannot contribute the truth to that comparison. The apex
diagnostic shows that the rational truth is already in the grammar. Whether and
when to continue escalation, and how to adjudicate later candidates without
weakening soundness, remain open. Escalating only on ambiguity would not address
the initial false-exact certificates. Neither the corrected scored zero-wrong
result nor significance against a random null settles this question. See the
[retained pilot evidence](ACQUISITION_PILOT_FINDINGS.md).

**Empirical result:** the implemented policy did not beat fixed x10 expansion.
It resolved the same 12/24 trials at the same 560 measurements. It did beat the
matched five-rung ladder's 640 measurements on those successes. The registered
stretch criterion required a strict gain against both; it **failed**.

**Empirical soundness:** zero structurally wrong terminal certificates across
all 96 corrected scored arms. This is conditional evidence from the manufactured
bank, not a universal checker guarantee. Earlier pilot failures remain recorded.

## Scored comparison

S1 was committed in 50db39d before any scored follow-up samples. Three frozen
initial rival sets (P2 seeds 30, 32, 35) each supplied a wide and a restricted
problem, with four fresh follow-up seeds 100–103. That is 24 paired trials per
arm and 96 arm runs; it is one rational truth family and three initial design
states, not 24 independent physical laws. Each initial state has at least two
representatives consistent with all 400 available observations.

| Strategy | Twin resolutions | Unresolved | Measurements per resolution | Total additional measurements, all trials |
|---|---:|---:|---:|---:|
| Refusal-guided |12/24|12/24|560|2320|
| Fixed x10 suggestion |12/24|12/24|560|1920|
| Matched five-rung ladder |12/24|12/24|640|3840|
| Native ladder, with common rival criterion |0/24|24/24|—|2832|

Successful totals include 400 initial observations, 80 new design observations
(or 160 for the matched ladder), and 80 independent final observations. Shared
initial data are charged to every trial as the cost of reproducing its starting
state; actual follow-up oracle calls are counted separately. The 6400 observations
used to screen the P2 bank are disclosed pilot setup, shared by every arm.

The native ladder needs a careful reading: it returned the correct rational
formula in all 24 raw runs. Those original-box observations did not eliminate
previously supported competitors, so it failed the common **twin-resolution**
criterion. Raw laws, verdicts and checks are preserved. This is not a claim that
native cannot recover the truth. It is a distinction between recovering the
right expression and measuring enough to reject the known alternatives.

Unresolved cases were retained. Under the registered 1801 penalty for unresolved
trials, mean censored costs are 1180.5 for both guided and x10, 1220.5 for the
matched ladder, and 1801 for native. Stopping unresolved is never scored as a
cheap resolution. Guided spent 400 more actual follow-up measurements than x10
across the unresolved controls, with no extra resolution.

All numerical forecasts in S1 matched:12 guided resolutions, 12 unresolved,
paired median differences 0 against x10 and -80 against the matched ladder,
zero wrong terminal certificates, and expanded domains for every resolved wide
case. The separate strict-gain criterion failed; the forecast tie is not a win.

## What was built and what it measures

**Empirical implementation:** `Result.rivals` retains structural refusal
representatives; passive re-splits preserve them and MCP exposes qualified
`design_evidence`. `choose_measurement` ranks a declared finite probe set by
finite rival disagreement per query cost. The opt-in `run_refusal_acquisition`
uses this choice, acquires a batch, and calls the unchanged discovery engine.

The corrected workflow requires a structural refusal and two full-design
survivors. It retains old competitors until measurements falsify them; their
absence from a later proposal list is insufficient. It declines to query when
predicted separation does not exceed the existing epsilon bands on any declared
probe. This is a finite-grid design heuristic, not a proof that the continuous
domain contains no distinguishing measurement.

All initial and intermediate discovery rows are DESIGN data, including the
engine's internal certification splits. A candidate and its design are frozen
before an independent final stream is generated. Any final failure terminates
that arm without feeding the observations back to selection. Accepted results
check every prior observation as well, and explicitly report the new finite
bounding ranges, row count and data hashes. The retained per-discovery alpha is
not presented as a new adaptive family-wide significance theorem.

The recorded winning guided choices all target x=300. Fixed x10 already reaches
an informative regime in one batch. Thus this bank supplies no evidence that
reading the particular refusal yields a better choice than simple expansion.
The ladder saving is avoiding its repeated original-box batch. In restricted
seed30 trials, edge disagreement at x=.5 (and once x=3) justified additional
probes, but none bought resolution. Restricted seeds 32 and 35 stopped without
new queries. No policy or tolerance was retuned after these scored outcomes.

**Open:** broader usefulness across other law families, heterogeneous query
costs and noisy data. This study rejects a gain claim for this implemented
heuristic on this bank; it does not establish that all refusal-guided design is
ineffective. The fixed x10 default remains unchanged.

## Apex: rational-d1

**Empirical:** the historical crc32-seeded draw still structurally refuses at 400
observations. Replaying the registered P1 refusal-model heuristic recovers the
exact `(2*x+1)/(x+3)` after measurement, using 560 observations and a new finite
domain extending approximately from 0.502 to 300. All observed rows and 80 fresh
final observations pass, and the returned expression is exactly equivalent to
the known truth. Fixed x10 also takes 560; matched ladder 640; native 598 including
its own holdout and the additional common final check.

There is an important limit to this result: both archived apex representatives
already miss original observations (2 and 3 misses on all 400 rows). They were
certifying split-local competitors, not two surviving full-data explanations.
The apex replay therefore demonstrates measurement-driven escape from the
engine's reach ordering, not the first falsification of those old models.
Its runner is pinned to 5f64d9a and clearly separated from the strict survivor
study. The original 35/36 passive reach artifact is not changed to 36/36; this is
a new acquisition regime and a new domain.

## Failures retained and guarantees scoped

P1's nonstructural seed10 inputs exposed false-exact initial certificates:
eight arm outputs accepted an approximant after original-box checks. P2 found
eight other initial draws with false-exact certificates and only 3/16 qualifying
twin states, missing its predicted minimum of 4. Every failed draw and the P1
policy are retained. P3 corrected the acquisition entry and rival-retention
rules; it did not repair these underlying engine behaviors or erase them from
the record. See [pilot findings](ACQUISITION_PILOT_FINDINGS.md).

The scoped headline paragraphs in README and INSTRUMENT_REPORT were completed
before the pilots. The new corrections have concrete before/after witnesses:
four new gate inputs fail on pinned P1 and all 15 gate cases pass afterward.
Scoring witnesses reject an x+100 replacement law, missing/duplicate arms,
a purported strict win at a tie, inherited old bounds, and a law different
from the one frozen before final observations. An unresolved 400-query trial
receives the registered 1801 censored cost. These checks are not credited as
new physical experiments.

## Reproduction and use

[Registration](REFUSAL_ACQUISITION_REGISTRATION.md) contains the sequence of
predictions, failed forecasts and amendments. The principal artifacts are:

- `experiments/acquisition_scored_protocol.json`: frozen code hash and S1 numbers.
- `experiments/results/acquisition_scored/`: all 96 outcomes and pre-query ledgers.
- `experiments/results/acquisition_scored_summary.json`: completeness, cost,
  soundness and domain audits for every scored trial.
- `experiments/results/acquisition_apex/`: the four qualified apex replays.
- `experiments/results/acquisition_reproduction.json`: actual fresh re-execution
  of one scored condition across all four arms. All four results and four plan
  ledgers reproduce byte-for-byte. It also records versions and byte-identical
  acquisition, certification and hypothesis-class sources against base e810b98.

The audit reconstructs initial and queried inputs, checks their hashes and all
oracle charges, verifies all 197 pre-query plans, reconstructs the declared
independent RNG stream for the 36 final samples, checks frozen-law identity,
rechecks accepted laws, and verifies that old compatible rivals were eliminated
by design measurements. No full-suite result is claimed:43 targeted test cases
passed, one test file per process. F/E9 lint and diff whitespace checks passed.

Run from the repository root, with the recorded package versions and thread
settings. Existing scored artifacts are verified and retained rather than
silently overwritten; use the reproduction runner for a genuinely fresh run.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_survivor_comparison --stage scored
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.score_refusal_acquisition
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.reproduce_acquisition
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_acquisition_gate_witnesses
.venv/bin/python -m pytest -q tests/test_refusal_acquisition.py
.venv/bin/python -m pytest -q tests/test_acquisition_scoring.py
```

For full historical P1/P2 reproduction, use their committed source revisions
(5f64d9a for P1, 577cefc for P2) in an isolated checkout with an empty stage output
directory. Each includes its registration and every retained original result;
the current runners intentionally reject a changed code hash. Tweezers data
were not reopened, no paid proposer was used, and no publication was pushed.

The library entry point is opt-in and currently restricted to finite,
deterministic observations on positive domains:

```python
from lagh.refusal_acquisition import RefusalPolicy, run_refusal_acquisition

policy = RefusalPolicy(admissible=((0.005,), (300.0,)))
outcome = run_refusal_acquisition(
    oracle, X_design, y_design,
    initial_box=[[0.5], [3.0]], policy=policy, seed=100,
)
# Credit a new law only when outcome.certified is True.
# outcome.domain names its measured finite domain; history records the evidence.
```
