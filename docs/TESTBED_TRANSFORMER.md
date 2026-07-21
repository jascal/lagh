# Test-bed registration: transformer-idiom recovery

**Registered 2026-07-21, before any recovery is run.** The frozen instrument (`lagh`, C1–C6
+ policy v2) is unchanged; this supplies the target declaration. Amendments dated and appended.

**One-line:** train a tiny transformer on modular addition mod `p` until it computes `(a+b) mod p`
**exactly**, then recover the idiom it computes — the modular structure — from black-box argmax
queries, **without reading its weights**; and use an **undertrained checkpoint as a negative
control** where no clean idiom exists.

## 1. Why this target — it closes the original THESIS.md arc

The founding claim of the whole program (`wyly/THESIS.md`): a trained transformer is *exact*
(deterministic given weights), *real* (we did not author its algorithm), and *undocumented*
(nobody wrote the closed form). `rosetta`/`threx` recovered one computed idiom (`THINGS[i+j]`) and
certified it. This is that class, with the new instrument: recover what a trained model computes,
certified, from queries alone.

| criterion | transformer idiom |
|---|---|
| exact law | ✓ a grokked model computes `(a+b) mod p` exactly (verified: 100% on the full table) |
| closed-form-expressible | ✓ the idiom is a **period-`p` quasi-polynomial** — exactly C6's class |
| undocumented | ✓✓ the closed form lives in distributed weights, never written; the learner reads only argmax outputs |
| queryable | ✓ forward pass, deterministic |
| noise / floor | ✓ integer (argmax) — the exact-integer C6 path, no ε |
| sealed GT | ✓ verify on held-out `a` and on `a` far beyond the queried range |
| abstention meaningful | ✓✓ the **undertrained control** has an error-corrupted map with no clean idiom — the honest output there is abstention, and it is the point |

## 2. The idiom and how it maps to C6

Fix `b = b₀`; vary the integer input `a` over `1..T_MAX`. The grokked model outputs
`(a + b₀) mod p`. Within residue class `a ≡ r (mod p)` this is `(r + b₀) mod p` = **a constant**.
So the idiom is a **period-`p`, degree-0 quasi-polynomial** — recovered by C6, which also recovers
**`p` itself as the period.** Recovering the period *is* discovering the modulus the model learned.

The float tiers C1–C5 run first and fail (a periodic sawtooth is no polynomial), so the engine
escalates to C6 — the natural, unforced path.

## 3. Targets

- **p = 11**, verified checkpoints: **grokked** (4000 steps, 100% exact) and **undertrained**
  (50 steps, 87.6% exact — the negative control).
- Recover on slices `b₀ ∈ {0, 3, 7}`, `a ∈ 1..48` (`T_MAX = 48`), for each checkpoint.
- Period search covers `p` (`P_MAX = 12 ≥ 11`).

## 4. Scoring

- **Certified** ⇔ C6 recovers a quasi-polynomial matching the model's argmax exactly on held-out
  `a`.
- **Correct (idiom)** ⇔ recovered period `= p = 11` **and** the quasi-polynomial matches the
  model's output on `a ∈ 49..88`, far beyond the queried range.
- **Zero-wrong invariant** carries over: a certified recovery whose extended-range check fails is
  confident-wrong and caps the bed. Standing record: **0 / 168**.

## 5. Predictions, registered before the run

- **PT-1:** the grokked model's idiom is recovered on all three slices, period `= 11`, validated
  out of range.
- **PT-2:** the **undertrained model abstains** on most/all slices — its map has ~12% errors, so
  no period-`p` quasi-polynomial reproduces it exactly, and exact-integer certification cannot
  admit an approximate one. This is the grokked-vs-ungrokked contrast, mirroring
  Ehrhart-vs-econ-sae.
- **PT-3:** **zero confident-wrong** — a certified idiom is always the true modular structure, and
  the undertrained control is refused, not fabricated.

## 6. What this can and cannot claim

A certified grokked recovery is a real result: *the idiom a trained transformer computes, recovered
and certified from black-box queries without reading its weights, including the modulus it learned.*
It is a **1-D-slice** recovery (C6 is 1-D in the dilation parameter); recovering the full 2-input
map is future work (a multivariate quasi-polynomial tier). The undertrained abstention is equally
the product. It is **not** a claim about large models or non-arithmetic idioms — it is the
minimal, clean instance of the founding thesis, with the new certified-abstaining instrument.


---

## VERDICT (2026-07-21). FLAWED DESIGN — reported as such, not as a success.

Surface result: grokked 3/3 "recovered period 11", undertrained 3/3 also "recovered period 11",
zero confident-wrong. But **the negative control (PT-2) did not abstain, and the reason exposes a
design flaw that voids the demonstration:**

**The oracle reduces `a mod p` before every query, so the model only ever sees `p` distinct inputs
per `b₀`.** Every 1-D slice is therefore a period-`p` constant-per-class lookup **by construction** —
whether the model is right or wrong. A corrupted (undertrained) model is *still* a deterministic
period-`p` map, just with some classes holding the wrong constant, so C6 certifies it too (it
faithfully matches the model, including its errors). **grokked and undertrained are
indistinguishable this way** — both are period-`p` lookups, and the intended contrast is impossible.

Verified: the undertrained model errs on 9 of 11 slices (e.g. `b0=1`: 3/11 wrong), yet those
corrupted slices still recover as clean period-11 laws that match the model exactly — because within
each residue class every `a` maps to the *same* model input, so the class is constant even when
wrong.

**The deeper lesson (this is the value):** a finite-vocabulary classifier **cannot be queried
outside its training grid** (no embedding exists for `a ≥ p`), so the real idiom-vs-memorization
question — *does the model compute the algorithm or store a table?* — **cannot be posed via
extrapolation**, and within the grid any deterministic map is a trivially-recoverable lookup. The
mod-`p` transformer is **fundamentally unsuited** to demonstrating certified idiom recovery, and
`oracle_fn`'s `a mod p` reduction made this concrete.

**PT-1 technically holds but is vacuous; PT-2 FALSIFIED; PT-3 (zero confident-wrong) holds** — and
that last one is the only real positive: lagh recovered every model's actual computed map exactly,
never fabricating, even when the model itself was wrong. The instrument did its job; the *target*
does not test what was claimed.

**Disposition: the mod-`p` transformer is RETIRED as an idiom target.** The founding thesis
(`wyly/THESIS.md`, `rosetta`/`threx`) is genuinely served by a **continuous-input** idiom that can
be queried off-distribution — index arithmetic over model activations (`THINGS[i+j]`), where a
grokked idiom extrapolates and a memorized table does not. That requires either a multivariate
quasi-polynomial tier or the activation-probing setup, and is left as **future work, honestly
scoped** — not claimed here on a flawed proxy.

**Record:** zero-wrong invariant stands (lagh never fabricated), now **0 / 174** scored tasks, but
this bed contributes **no positive idiom-recovery claim** — a null result from a flawed target
design, reported as one.
