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
