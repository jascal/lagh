# BLIND READ REPORT — LLM-SRBench (read spent 2026-07-24/25)

One read, per `BLIND_READ_REGISTRATION.md` (frozen pre-download; two logged
amendments: 600s/problem cap, 240-problem count). Scores computed once by the
registered conservative sympy judge. **This is a loss on the headline metric,
reported as registered — win, lose, or mixed.**

## Scorecard vs the frozen SOTA

| | lagh (this read) | frozen SOTA (LLM-SR, GPT-4o-mini) |
|---|---|---|
| LSR-Transform SA | **0.9%** (1/111) | 31.53% |
| LSR-Synth SA | **0.0%** (0/129) | 20.24% (LaSR 28.12%) |
| ALL SA | **0.42%** (1/240) | — |
| ALL Acc₀.₁ (ID test) | 11.25% (Transform 19.8%) | not the headline metric |

Pre-registered accounting: **track A submitted 0 / wrong 0** (the invariant is
intact but vacuous — no certificate was ever issued); track B submitted 95, SA-
correct 1 (`2·x_0/(x_1·x_2)`, a genuine blind recovery by the log-log probe);
abstain-with-nothing **145/240 (60%)** — the policy predicted this would be
rare; it was the modal outcome. 80/240 problems hit the 600s cap. 14 gt
expressions failed the conservative judge's parser (functional notation like
`A(t)`); by construction this can only under-count our SA.

## Root cause 1 (dominant, structural): the benchmark's data is float32-quantized

Post-read diagnosis: the ground-truth equations do NOT reproduce their own train
data at float64 precision — median relative residual ≈ 4-9 × 10⁻⁸ (= float32
epsilon; the HDF5 arrays are float32-born, cast to float64 by the loader). lagh's
`sigma=0` epsilon is ~2.2 × 10⁻¹³ relative. **Exact certification was impossible
on every one of the 240 problems by data representation alone.** Zero
certificates is therefore not a grammar verdict — it is the honesty core
refusing to claim machine-precision exactness on quantized data, on all 240,
correctly.

Lesson (the big one): **representation precision is DECLARED noise, not
statistical noise.** The array dtype is metadata, readable before discovery
without contamination. Running with σ_rep ≈ 1.2e-7 (declared, relative) would
have engaged the structural-certification regime (`RNOISE_STUDY.md`: exact
structure to ~1% noise, zero fabrication) on every problem. The frozen protocol
said `sigma=0`; the read is spent under it; this is the first change any future
benchmark run makes.

## Root cause 2 (Transform half): inversion orthogonality

LSR-Transform is built by solving Feynman equations for a pivot input —
manufacturing quadratic-root shapes, nested radicals, inverse-trig of compound
arguments. The forward-form curriculum (C1–C9) does not span inversion images;
CAP-G covers one such family. 0 certificates + 71 timeouts on this half was
predicted from the public construction before any score existed. The scoped
capability family is **inversion closure**: for each forward family in the
grammar, add its solve-for-pivot forms.

## Root cause 3: the conjecture floor is a floor

Track B was the `fit` scout's log-log probe (or a timeout fallback of the same
shape). It still put 27/95 conjectures within 10% max-relative-error on ID test
data, and landed one exact SA hit — but a monomial probe cannot express sums,
which is most of the benchmark. The k-form LLM proposer (already built for the
machine) was deliberately not in the frozen submission path; wiring it in is the
obvious upgrade.

## What survives this read

- **Zero confident-wrong remains unbroken** across the program — trivially here
  (no certificates), substantively on dev (0/108 clean, 0/103 under noise).
- The two-track policy worked as designed: every submission labeled, the
  certified partition empty rather than polluted, the accounting exactly as
  pre-registered.
- The instrument diagnosis is crisp and actionable: declared-representation
  epsilon (small, high-leverage), inversion closure (medium), LLM-proposer
  track B (medium). None of it is benchmark-specific tuning.

## Standing of the claims

Per `STRATEGY.md`: this was the reserved blind read and its numbers are the
headline. **lagh does not beat LLM-SRBench SOTA under the frozen protocol.** The
honest claim after this read: *a certified law discoverer with a proven
zero-wrong record on clean data, whose exactness bar correctly refused a
float32-quantized benchmark, and whose competitive submission mode does not yet
exist.* LLM-SRBench is now spent as a blind set (opened, diagnosed); if a future
run targets it, it is a DEV benchmark and must be labeled as such.
