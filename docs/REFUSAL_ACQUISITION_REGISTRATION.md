# Refusal-guided acquisition: preparation and pilot registration

Status: **open**. Written before new experiments. This is the pilot protocol,
not the scored-run registration. Numerical scored predictions must be appended
from the pilot and committed BEFORE evaluation. The incoming task is truncated;
the second precondition and missing constraints were initially pending clarification.
The scored protocol remains incomplete; no scored run has been performed.

## Complete instruction received

docs/STRETCH_GOAL_ACQUISITION.md supplies the missing text. Preconditions are
the scoped paragraph in README and INSTRUMENT_REPORT and a counter-input for
every new gate. Both documentation paragraphs are now present. An acquired
certificate must explicitly name its new finite domain. Apex success is optional;
the fixed-baseline comparison and measured soundness are mandatory.

## Pilot P0: apex design inspection (registered before execution)

**Open prediction P0:** the historical rational-d1 draw still returns a structural
refusal with at least two retained representatives. Its fixed input is 400 rows
uniform on [.5,3], RNG seed crc32('rational-d1'), ordinary passive defaults
(three re-splits, seed0), sigma0 and floor1e-12. Record every representative,
verdict and measured range. No target truth enters the selector.

On that refusal, score a fixed 257-point geometric grid in [.005,300] at unit
query cost using the new oracle-free selector. **Open prediction P0b:** its
selected query is outside the initial box. Grid evaluations are model
predictions, not measurements. Register all rival values at the selected point
before making any follow-up query. P0 does not claim certification or a gain.
If the initial result certifies or has fewer than two rivals, record failure
and retain the result; do not force a structural refusal. No final certification
sample is generated during this inspection. Full pilot acquisition settings and
scored predictions follow from this inspection, with amendments recorded before
the affected run. The task is no longer blocked on missing instructions.

## Scope

Use structural refusal representatives to select measurements on a declared
admissible domain, maximizing finite rival disagreement per declared query cost.
Keep the hypothesis class, certification gates and fixed baseline policies
unchanged. Rational-d1 is the apex, with truth (2*x+1)/(x+3), and must go through
ordinary discovery after acquisition: neither historical rival is the truth.
Rejecting one or even every old rival is not itself a certificate.

**Empirical historical baseline:** the Müntz and reach artifacts report 35/36,
with rational-d1 structurally unresolved. The current base is e810b98 (merged
PR14). The fixed comparators are recover's x10 box suggestion and acquisition's
five-rung _box_ladder, in their existing order. Any matched-budget experiment
must be labelled separately from a literal unchanged-baseline run.

## Information boundary

1. Register admissible input bounds, query costs, budgets, candidate probe grid,
   splits and seeds before oracle queries. No target-specific truth is passed
   to the selector. Truth is available only to the generator and scoring code.
2. Initial observations and every observation used to produce the refusal are
   DESIGN evidence. This includes the engine's internal certification rows:
   their resulting rivals influence measurement selection.
3. The selector accepts only rival expressions, admissible probe inputs and
   declared costs. It cannot inspect future oracle responses or final check
   samples. Save the selected design and predicted rival values before query.
4. Acquire fitting/selection evidence on that design and use the unchanged
   engine. A final acceptance also requires a fresh sample generated only
   AFTER the candidate and design are frozen. Keep final evidence sealed from
   later design updates. Stop on a failed final guard; do not tune against it.
5. Count every oracle output: initial data, ranging, replicates, attempted
   acquisitions, failed attempts and final checks. Report actual query counts,
   separately from the existing baseline ledger, which omits its 60-row guard.
   Independent truth-scoring probes are reported separately, never fed back.
6. A candidate must also remain consistent with prior measured evidence. A
   narrower new box cannot silently erase old contradictions. Multiple testing
   exposure and the scope of any significance claim must be explicit.

## Pilot questions and kill criteria

**Open predictions:** actual rival expressions can be retained without changing
existing verdicts; finite disagreement can select a design without an oracle;
new observations can sometimes eliminate a structural ambiguity at lower cost.

Before crediting each test, exhibit its failing input. Seed witnesses: absent
rivals; identical rivals; a rival undefined at the proposed query; nonpositive
or nonfinite cost; a design outside admissible bounds; a final sample with
changed target values; leakage of final responses into design; a selector that
ignores rival content; both old twins wrong. Undefined predictions are not
infinite information. No finite informative query means explicit refusal.

Pilot families and sample sizes, then scored N, resolve-count and query-count
predictions remain to be fixed. Include known-truth distinguishable cases and
controls whose rivals remain indistinguishable on their admissible domain.
Do not remove failed pilot regimes from the report. A test of twin elimination
is distinct from end-to-end discovery and certification.

**Kill:** no independent design/certification split -> stop the dependent work.
Any confident-wrong -> fail the soundness prediction and retain the artifact.
No lower query cost -> report a failed gain prediction; do not redefine cost or
remove unresolved cases. An untestable prediction is withdrawn. No grammar,
tolerance or baseline change may rescue a failed comparison. Tweezers stays
closed; no paid proposer or publication push.

## Preparation evidence

**Empirical:** 13 oracle-free unit cases in tests/test_measurement_design.py
passed, each with its concrete counter-input stated in the test. Opposite
rival shapes choose opposite endpoints of the same probe set; a high query
cost changes the choice; identical/absent/undefined rivals cannot manufacture
information. Invalid costs and inadmissible probes raise ValueError. Mutation
of a caller's probe array cannot change a stored choice.

This is unit coverage only. No end-to-end acquisition, rational-d1 replay,
pilot, scored run, or sample-efficiency improvement has been measured.

Reproduce the unit cases in one process:
`.venv/bin/pytest -q tests/test_measurement_design.py`.

### Refusal transport counter-inputs

**Empirical:** five mocked-split witnesses fail against preparation commit
cef9feef213c860e67c6eb8ec9bf4f9300032ea2 and pass after transport repair.
The artifact is experiments/results/refusal_transport_witnesses.json.

- A structural first split followed by a rival-free last refusal previously
  erased the first pair. It now retains it.
- The same first split followed by a successful split remains an abstention
  under the unchanged sticky guard, and now retains the refusal's pair.
- Two overlapping pairs from separate splits are retained in first-seen order
  without syntactic duplicates. They are not claimed to be a single exhaustive
  coherence partition or laws that pass on every observed row.
- Passive and active MCP recover now expose these expressions in
  design_evidence, tagged empirical and explicitly requiring fresh
  certification. They remain absent from the returned law field. Passive
  recover's fixed suggested_box remains [lo/10, hi*10].

Reproduce before/after witnesses (one test file per process):
`.venv/bin/python -m experiments.run_refusal_transport_witnesses`.
Their failure inputs and expected behavior are in tests/test_refusal_transport.py;
the pinned run actually exhibits the loss, rather than merely naming a possible
mutation. No acquisition or discovery-success claim follows from mocked splits.
Independent final validation and the acquisition comparison remain unimplemented.

**Empirical regression check:** the 17 existing MCP tests pass (one process),
including wrong-form, below-floor, weak-significance and unpinned-coefficient
counter-inputs. Three existing overflow warnings occur in scout/abstain cases.
This is regression coverage of those cases, not a global soundness guarantee.
The five existing passive tests also pass in a separate process, including
the wrong-law full-data demotion witness and the irrational-power abstention.
The latter emits an existing overflow warning. Together with the five repaired
transport cases, this checkpoint has 27 passing cases; no full-suite run is
claimed.

## P0 result and acquisition pilot P1 (before execution)

**Empirical P0:** the historical draw still refuses with two representatives,
dof34 and45. Unit-cost maximum disagreement selects x=300, outside [.5,3];
its rival predictions are about 6.7e6 and -1.6e127. Both P0 predictions hold.
Large extrapolated disagreement is a selection heuristic, not reliable physics.
Artifact: acquisition_pilot_p0.json; no follow-up measurement was made in P0.

**P1 protocol, open predictions:** four manufactured cases from the existing
rational reach family: (2x+1)/(x+3), (3x+2)/(x+4), (x+1)/(x+5), and a restricted
copy of the first. Each starts with400 uniform design rows on [.5,3], seed10.
The first three allow [.005,300]; the last permits only [.5,3]. All use sigma0,
floor1e-12, unchanged passive discovery with three re-splits. Initial non-twin
cases remain in the report, explicitly ineligible for a twin-resolution gain.

Four arms share the same initial observations and charge them as400 measurements:

1. Guided:257 geometric probe points, unit cost, maximum finite rival spread.
   Acquire80 points in [q/3,3q] intersected with admissible bounds, including q
   itself; other79 are log-uniform. Append to design observations and rediscover.
2. Fixed x10: use [observed_min/10,observed_max*10] as in recover,80 log-uniform
   measurements per attempt, append and rediscover. Stop when the whole next box
   would exceed admissible bounds; do not silently clip the baseline.
3. Matched fixed ladder: exactly _box_ladder's five boxes and order,80
   log-uniform observations per box, append and rediscover. This matched wrapper
   isolates WHERE to measure; it is not a claim to run native run_active.
4. Native fixed ladder: run_active_boxsearch unchanged, budget200 per box,
   max_boxes5, default Policy, no wall-clock scoring cutoff. Only admissible
   prefix rungs are allowed. Count every actual oracle call, including existing
   holdouts. Report its native verdict separately from the shared final guard.

Guided and matched arms have at most5 acquisition rounds. Global ceiling1800
observations per arm includes initial data and a reserved80-row final guard.
The native ladder's predeclared ceiling (initial400 + five200 budgets + five60
holdouts + final80) is1780; actual run_active round ledgers are also retained.
Only deterministic finite positive-domain oracles are in this first experiment.

Selection sees design data only. All internal discovery splits are treated as
design. Freeze the first returned candidate, then generate80 fresh log-uniform
points across the measured design bounds using a separate SeedSequence child.
Check that frozen candidate on ALL acquired/design observations and these fresh
rows. Stop after any final guard failure: never return final responses to the
selector or retry a law. Recompute finite-domain bounds and row count explicitly.
No new noisy-data or adaptive significance theorem is claimed; report search
exposure rather than treating a per-discovery alpha as a family-wide guarantee.

**Pilot predictions:** at least one initial case supplies two rivals; at least
one unrestricted arm resolves; the restricted case remains unresolved. These
can fail. P1 is for measuring query counts and choosing scored N/thresholds;
no efficiency prediction is promoted from this pilot. Scored seeds and numbers
will be separately frozen and committed. Full truth is used only for oracle
responses and post-run scoring (exact rational equivalence plus an independent
257-point probe); it is never given to the selector. Count every structurally
wrong certified law, even if it agrees locally. No failed arm is removed.

New gate witnesses registered before tests: invalid/nonfinite oracle rows cannot
vanish in discovery's finite-row filter; a law agreeing only in the selected
box must fail prior-data checks; a changed fresh target must demote and stop;
a fresh sample outside the new claimed bounds exposes inherited-domain bugs;
budget below batch+reserved-final must cause zero new queries; same design and
different future final responses must produce byte-identical pre-query plans.
