# Refusal-guided acquisition: preparation and pilot registration

Status: **open**. Written before new experiments. This is the pilot protocol,
not the scored-run registration. Numerical scored predictions must be appended
from the pilot and committed BEFORE evaluation. The incoming task is truncated;
the second precondition and missing constraints remain pending clarification.
The scored protocol remains incomplete; no scored run has been performed.

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
