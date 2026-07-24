# lagh

*(Scottish Gaelic: law)*

A general-purpose **certified law discoverer**: given data from a system with an exact
underlying law, return either a symbolic law with a machine-checked certificate over a
stated domain, or a machine-readable refusal saying why no law is certifiable — and
never a confident wrong answer.

**Design:** [`wyly/docs/DISCOVERER.md`](../wyly/docs/DISCOVERER.md). Every mechanism
here is either something the predecessor project measured working, or the registered
fix for a failure it measured. The predecessor's record across 114 scored benchmark
tasks: **zero wrong submissions** — that invariant is the product definition.

## Status: P1 (core extraction + curriculum C1–C5)

```python
from lagh import discover
r = discover(X_fit, y_fit, X_sel, y_sel, X_cert, y_cert, sigma=0.0)
r.certificate.certified   # True -> r.expr is the law, exact rationals
r.certificate.abstain     # else: domain|structural|noise|surrogate|numerical|range
```

- `lagh/certify.py` — the fixed honesty core: exhaustive check, four-term ε model,
  null-law vacuity, set-coherence, the exact-coefficient gate (unpinned Float /
  huge-denominator coefficients never certify). Never per-target, never per-class.
- `lagh/classes/` — the math-class curriculum, complete at registration:
  C1 polynomial (+ fractional/trig-product/high-degree features) · C2 rational ·
  C3 power-law · C4 inner-scaled transcendentals · C5 target transforms (incl.
  bose/fermi) · C6 quasi-polynomial · C7 Lévy · C8 angular/inverse-trig ·
  C9 generalized monomial (e^x / ln x factors).
- `lagh/engine.py` — MDL-ordered tier escalation; escalates **only** when the
  certifying set is empty (ambiguity is a verdict, not a reason to escalate).
- `lagh/acquisition.py` — P2 active loop: adaptive ranging, ledgers, box-search.
- `lagh/passive.py` — the fixed-dataset regime (no oracle): K re-splits + a
  full-data exhaustive gate (`docs/DIRECTION_PASSIVE.md`).
- `lagh/submit.py` — two-track benchmark submission: certified law, else labeled
  conjecture (`docs/DIRECTION_OUTPUT_POLICY.md`).
- `machine/` — the research loop as a verified Orca state machine; the LLM is a
  bounded k-form proposer at exactly one state (`docs/DIRECTION_MACHINE.md`).

```bash
.venv/bin/python -m pytest tests/ -q
```
