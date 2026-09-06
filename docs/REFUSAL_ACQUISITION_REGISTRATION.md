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

## P1 result and P2 survivor screen (registered before execution)

P1 results and failures are in ACQUISITION_PILOT_FINDINGS.md. The positive
soundness reading is withdrawn: eight terminal certificates from nonstructural
initial inputs were approximants, and the intended restricted twin did not
exist. The two structural inputs tied fixed x10. Do not present this as a gain
over both fixed strategies. All16 outcomes remain recorded.

**P2 screen, open prediction:** of16 datasets from truth (3x+2)/(x+4), each400
uniform rows on [.5,3] at seeds20 through35, at least4 produce a structural
refusal with at least2 representatives passing check on ALL400 rows at the
unchanged epsilon (sigma0, floor1e-12). Ordinary passive defaults remain frozen.
Record nonstructural outcomes, all per-rival miss counts and exact truth of any
initial certificate. This is a DESIGN-only pilot screen,6400 oracle observations,
no follow-up or final observations. No seed is discarded from the artifact.

If fewer than4 qualify, score that failure before any amendment. If none
qualify, stop the strict survivor-policy experiment and report the absence of
the required rival set; do not relabel already-refuted classes as survivors.
Qualified frozen design states can be used as manufactured twin problems with
fresh, separately registered acquisition/final seeds. Conditioning on their
initial refusal is explicit; these will not be a random population sample.

Proposed corrections, not yet executed: require structural entry; retain only
full-design-compatible representatives for strict disagreement; require new
evidence to exclude old incompatible explanations before crediting resolution.
Every correction needs its own counter-input and pre-run amendment. No old P1
artifact is overwritten, and no core hypothesis class or tolerance changes.

## P3 corrected survivor-policy pilot (registered before execution)

P2 delivered3 qualifying draws, seeds30,32,35, rather than the predicted4.
The numerical prediction FAILED. Use all three as the frozen twin bank: same
truth (3x+2)/(x+4), different initial400-row draws and fitted competitors. For
each, run both unrestricted [.005,300] and restricted [.5,3] acquisition domains.
That is6 problem conditions, not6 independent physical laws. P3 follow-up seed99;
the future scored follow-up seeds100,101,102,103 remain untouched.

All P1 policy numbers and fixed selection policies remain unchanged. The
following common acceptance/entry corrections are explicit amendments:

- Require an actual structural initial refusal with at least2 representatives
  compatible with every initial row. Nonqualifying inputs return an uncovered
  acquisition claim, not the incoming certificate. Preserve the incoming status
  in the experiment artifact. This fixes P1's out-of-scope acceptance branch;
  it does not repair or hide the underlying discovery behavior.
- Use only currently full-design-compatible rivals for guided disagreement.
  At later rounds, retain any old representatives still compatible with all
  observed design rows. New measurements, not disappearance from a later
  proposal list, must eliminate an old competitor.
- Query only if finite predicted separation exceeds the sum of the two
  existing epsilon bands (evaluated at the extreme rival predictions). This
  is a deterministic design heuristic using the unchanged machine/floor model,
  not a new statistical significance threshold. No informative finite probe
  means refusal; it is not proof about an unprobed continuum.
- Before freezing an engine candidate, require every surviving old rival to
  be that same symbolic expression or to have failed on measured design rows.
  Otherwise continue the matched strategy using DESIGN evidence only; no final
  sample has yet been drawn. Native ladder is left byte-identical and its raw
  verdict recorded: if its returned law does not eliminate the prior rivals,
  the wrapper reports no twin resolution. It does not rerun native on the final
  data or silently replace its stopping policy.
- The common independent80-row final guard and explicit expanded finite domain
  remain as P1. Any final failure terminates the arm without reuse.

**New counter-inputs before credit:** an initially certified x on exact x data
must not count as a refusal resolution; a pair(x,2x) on exact x data is not two
survivors; a compatible pair(x, x+1e-16*x**8) on [.5,3] supplies a valid unit
fixture; a later proposal list containing only x does not eliminate its old
compatible twin on that restricted domain; an all-data-rejected rival cannot
determine the next query. The same pair has detectable separation at x=300;
the admissible-domain change must change queryability. A future changed final
response must still not change any pre-query plan. These are mocked unit
witnesses, not manufactured scientific successes.

**P3 numerical predictions (open):** guided resolves3/6 conditions, with at
least2/6 remaining unresolved; successful guided and x10 cases use560 total
queries, so their paired median difference is0 (NO strict gain predicted).
Matched ladder requires640 on successful cases. Native raw successes are not
automatically twin resolutions; record both. Zero wrong terminal certificates
is required in all corrected arms. Scored numbers will be finalized from P3
before running any scored follow-up seed.

**Separate apex diagnostic:** replay the historical crc32 rational-d1 draw
using the original P1 refusal-model heuristic, loaded from pinned5f64d9a,
with follow-up seed99, all four arms and the existing80-row final guard. Predict
guided obtains the exact rational at560 queries and ties x10; this can fail.
This is explicitly NOT a surviving-twin experiment (both initial models already
miss old rows). It tests measurement-driven escape from the engine's reach
ordering without changing its grammar. Record its new finite domain, exact
truth check, cost and that qualification. Do not relabel it a36/36 passive audit.

## Scored S1 registration (after P3, before any scored follow-up samples)

P3 results are in ACQUISITION_PILOT_FINDINGS.md. Freeze the current study code
hash in experiments/acquisition_scored_protocol.json. Run the same6 manufactured
problem conditions (three frozen P2 initial states, each with wide/restricted
admissibility) at follow-up seeds100,101,102,103. N=24 paired acquisition trials
per arm,96 total arm runs. There are only3 distinct initial rival sets and one
rational truth family; the four follow-up replicates do not create24 independent
physical laws. Selection/discovery/final rules and every numerical policy
parameter remain exactly P3. No scored seed has been used in P0-P3.

**S1a, empirical prediction:** guided, fixed x10 and matched ladder each resolve
12/24 trials;12/24 remain unresolved. In particular, at least8/24 guided trials
remain unresolved. Successful totals are560 for guided/x10 and640 for matched
ladder. These numbers come directly from P3, not from the scored seeds.

**S1b, empirical prediction:** paired median(guided queries - fixed x10 queries)
is0 on jointly resolved trials; against matched ladder it is-80. Thus no strict
gain over x10 is predicted. The stretch criterion is nevertheless explicit:
guided must resolve at least as many trials as each fixed comparator AND have
strictly smaller paired median measurement cost than BOTH, with zero wrong
certificates. A tie fails that criterion. Do not redefine a tie as a win.

**S1c:** zero structurally wrong terminal certificates across ALL96 corrected
arms. Raw native laws and raw initial statuses are separately retained; their
correct formula outputs do not imply known rivals were eliminated. Native is
predicted to have0/24 twin resolutions while retaining its raw successful laws.
Any wrong terminal certificate fails soundness; preserve the result and withdraw
the gain claim. No core correction or threshold retuning occurs on scored data.

**S1d:** every resolved wide case names a finite domain extending beyond the
initial [.5,3] box; every result counts all actual follow-up oracle observations,
including final and native holdouts. The400 shared initial rows are charged to
each trial as the cost of reproducing its starting state; new actual oracle
calls are recorded separately. Manufacturing the P2 bank cost6400 design rows
and is reported as pilot setup, not silently charged only to one arm.

Report all unresolved trials, not only successful pairs. For a complementary
aggregate, charge unresolved trials the common ceiling+1=1801 (censored cost);
never count stopping early unresolved as cheap resolution. Report both actual
spent queries and this declared penalty. This finite controlled comparison
cannot establish a universal soundness guarantee or a general optimality result.

**Witnesses for scoring:** replacing a true terminal law with x+100 must mark
it wrong; a missing arm or duplicated trial must invalidate completeness; a
fixed-x10 tie must fail the strict-gain criterion; an unresolved400-query run
must receive1801 rather than outperforming a resolved560-query run. These
scoring counter-inputs must be exhibited before crediting the scored summary.

## S1 scored result

**Empirical:** all96 arms completed at the frozen code hash. Guided/x10/matched
ladder each resolved12/24; native had24 raw correct formulas but0 twin
resolutions. Guided and x10 used560 per resolution, matched ladder640. Zero
corrected terminal certificates were wrong. All numerical S1 forecasts held;
the separate strict-gain criterion failed because guided tied x10. Guided also
spent400 extra actual follow-up observations on unresolved controls relative
to x10. No scored-data tuning or replacement experiment followed.

See REFUSAL_ACQUISITION_RESULTS.md and acquisition_scored_summary.json for the
complete account, including native's raw-versus-resolution distinction, the
qualified apex recovery and the retained initial-checker failures. The summary
audits all96 trials,197 pre-query records and36 independent final samples.


## Post-result interpretation correction (2026-09-06)

**Empirical:** P2's 8/16 false-exact initial certificates occurred on clean,
noiseless data. The approximant-impostor boundary limits exact-form claims in
both clean and declared-noise regimes. The first non-empty certifying tier
terminates escalation, preventing later-tier truth from entering that comparison.
P1's eight wrong arm outputs remain retained; S1's zero wrong terminal
certificates is conditional on its corrected refusal-only entry and bank.

**Open:** when to continue escalation past a certifying tier and how to
adjudicate later candidates at unchanged soundness. No escalation repair or
validation is claimed. This is a post-result scope correction, not a change
to frozen predictions or scoring. The full-data apex replay also withdraws
permanent-indistinguishability language for its archived split-local rivals.
See [results](REFUSAL_ACQUISITION_RESULTS.md).
