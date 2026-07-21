# Direction (registered 2026-07-21): certified significance on random / pseudorandom domains

**Status:** design direction, scheduled *after* the NewtonBench-dev → blind-benchmark arc
(`STRATEGY.md`). Recorded now so it is not reverse-justified later.

## The gap it closes

Today's certificate says *"the law fits every point of the stated domain."* That is a statement
about the **sample**. On a thin domain a spurious law can fit by chance — this was the `|D|=30`
MDBench failure and the agreeing-impostor risk. The zero-wrong record (0/174) is **empirical**;
it has no *quantified* guarantee attached. This direction attaches one:

> **Significance certificate:** *"the law fits, and under the null (a random target of the same
> range and size) a spurious certification of this hypothesis-class complexity occurs with
> probability < α,"* with α **computed**, not hoped.

## Why it is analytically tractable (the good news)

For exact / tolerance certification the significance is a counting bound, not a fitted statistic
(so it is *not* the significance-gated abstention the Tübingen post-mortem warns about — there is
no test-set selection, the bound is a-priori):

- a candidate with `dof` free parameters is pinned by `dof` points; the remaining
  `h = |D| − dof` held-out points must each match by chance;
- per-point chance-match probability under a null of value-range `R` (or tolerance `ε` over
  `[y_min, y_max]`): `q ≈ 1/R` (discrete) or `q ≈ 2ε/(y_max−y_min)` (continuous);
- union bound over the `|H|` candidate forms actually searched:
  **`α ≤ |H| · q^h`.**

For a degree-2 law certified over 30 held-out integer points in range 100, with `|H| ~ 10⁴`
searched forms: `α ≤ 10⁴ · 100⁻³⁰ ≈ 10⁻⁵⁶`. **Exact certification is already hyper-significant,
and lagh can print the number** — turning "fits the sample" into "fits, p < 10⁻⁵⁶."

## Random / pseudorandom domains do three jobs

1. **Null calibration (validate the bound).** Run recovery on many *true-random* targets
   (`os.urandom`, a crypto RNG). The observed false-certification rate must be `≤ α`. If it
   **exceeds** α, the effective `|H|` was undercounted (the transform/inner-scale tiers multiply
   the search) — a genuine, quantified diagnostic that bounds the agreeing-impostor risk instead
   of only observing it hasn't bitten yet.

2. **Graded recovery targets (a real exact-but-undocumented benchmark).** Pseudorandom generators
   form a difficulty gradient with a crisp success criterion:
   - **weak LCG** `x_{n+1} = (a·x_n + c) mod m` — the recurrence is an exact integer law (C6 /
     modular family); the specific `(a,c,m)` of a random instance is written nowhere → recover and
     certify. Exact-but-undocumented, textbook example.
   - **truncated LCG** — harder; partial state hidden.
   - **crypto PRNG** (AES-CTR, ChaCha) — the recurrence is infeasible to recover → **abstain**,
     and the abstention is now *meaningful and certified*: "no law in class C certifies at
     significance α."

3. **A dual certificate that is itself a product.** Every run ends in one of:
   - *recovered the generator, significant at α*, or
   - *no generator recoverable — indistinguishable from random-in-class-C at level α.*

   The second is a **certified randomness statement relative to a hypothesis class** — a real,
   citable object (crypto/randomness-testing adjacent) that no SR tool emits. It is the abstention
   thesis at its strongest: refusal with a *quantified* guarantee.

## Honest subtleties (recorded up front)

- **`|H|` must be upper-bounded, not guessed.** Count the candidates actually generated per run
  (the engine already enumerates them); the inner-scale grids and transforms inflate it and must be
  included, or empirical null calibration will (correctly) reject the bound.
- **The null must be stated.** `α` is *relative to a null distribution* (uniform-in-range, or the
  target's own marginal). A structured pseudorandom sequence is not null — that is the point: a
  recoverable PRNG *is* structured, a good one is not.
- **Independence of held-out points** is assumed by `q^h`; for autocorrelated domains (trajectory
  tails — the MDBench sin) the effective `h` is smaller and must be discounted. This is the same
  "thin tube" caveat, now made quantitative.

## Where it slots

After the blind SR benchmark. It is not a benchmark win itself; it is the **theoretical upgrade to
the win's foundation** — converting lagh's headline from "never wrong on our tests" to "certifies
at a stated, validated significance level, and its refusals are certified randomness statements."
That is the sentence that makes the instrument publishable *and* a benchmark contender at once.
