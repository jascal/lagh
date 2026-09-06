# Refusal measurement audit — registration, 2026-09-05

Status: **open**, registered before executing the new experiment.

First bounded deliverable of the research handoff: audit every `certify.Abstain`
member, distinguish live producers from reserved vocabulary, and retain the
exhaustive residual measurement at the public `verify` refusal. The optimal
weak test, approximant boundary and LawSystemBench v1.3 remain separate open
research deliverables; this experiment makes no reach or theory claim.

Predictions:

1. On deterministic `verify('x_0', X, y)` with y = 2x + 0.1x², 100
   equally spaced x in [1, 2], floor 1e-6, the structural refusal will retain
   signed residuals and their declared bands for exactly the checked split.
2. A second case with the mismatch confined to a fitting row will reach the
   full-data refusal and identify its domain as all finite supplied rows.
3. A clean 2x control keeps its original verdict and law. No refusal upgrades,
   and the existing public-surface regression file passes.
4. The residual alone cannot distinguish an additive instrument offset from
   an omitted physical contribution. Both explanations produce identical
   observations; attribution remains explicitly unresolved. No confidence
   level is invented for an epsilon envelope.

Kill criteria: discard the diagnostic if it changes any decision, evaluates a
new candidate, or attributes cause from residual alone. Nonfinite arithmetic
must omit the diagnostic, never replace the original refusal with an exception.

Artifact: `experiments/results/refusal_measurements.json`. Tag: **empirical**.

## Executed result

P1/P2 met on the registered fixtures (20 and 100 rows respectively); P3 clean
control certified `2*x_0`. Targeted tests: `test_refusal.py` 3 passed;
`test_tool_surface.py` 55 passed (two existing empty-array warnings). P4
remains explicitly unattributed; no physical attribution experiment was run.
No kill criterion fired. Full research program and full-suite validation are
not claimed by this bounded result. See `REFUSAL_MEASUREMENT_AUDIT.md`.

## Review amendment — registered before revised execution

Claude's PR review found that the partial union hull can hide a detected
nonzero discrepancy and falsely compose unrelated datasets. Remove `partial`
entirely; ship only empirical measurements. Preserve the original top-level
`law` behavior. Predictions: finite rows survive an undefined neighbor; malformed
index alignment reports an omission; scalar epsilon broadcasts; at most 64 rows
are returned, prioritizing largest band exceedances, with complete counts and
maximum excess over every measurable row. The measurement names its candidate.
The fixed full-domain fixture must assert its intended gate/domain, so a changed
split fails the experiment instead of silently changing its meaning.

Revised execution: 7 focused measurement tests, 55 public-surface tests, and
10 MCP verify/server tests passed (7 unrelated MCP tests deselected). Public
surface retains two empty-array warnings. Ruff F,E9 and diff whitespace checks
passed. The 20,000-row synthetic diagnostic returned 64 rows, counted all 100
misses and 19,936 elisions, and serialized below 6,000 characters. Fixture
artifact regenerated; clean control unchanged. No partial determination ships.

## Falsifiability correction (2026-09-05)

P4 is **WITHDRAWN as an empirical prediction** and restated as an identifiability
limitation: identical observations cannot discriminate their identical-input
explanations. It supplies no independent physical validation. P3's diagnostic
noninterference requirement remains; rejecting invalid/nonfinite check inputs
is a separately registered soundness correction, not a diagnostic upgrade.
