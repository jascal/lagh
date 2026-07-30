# Level 0 of the stochastic suite: what the coverage factor actually buys

**Run 2026-07-29.** `experiments/stochastic/run_level0.py` →
`experiments/results/stochastic_level0.json`. Machinery: `lagh/ito.py` (the Itô
weak form), `lagh/certify.py` (`coverage_factor`, `admissible_interval`), scored
through the frozen checker `lagh/stochcheck.py`. Predictions S1–S7 were registered
in [`DIRECTION_STOCHASTIC.md`](./DIRECTION_STOCHASTIC.md) **before** any of this
existed; this document reports what happened to them, including where they failed.

This is a **DEV** result by construction (`STRATEGY.md`): the suite is
self-authored, the systems are simulated from exact transition laws, and its
numbers drive development. Nothing here is a benchmark claim.

## The deliverable: κ holds, and with a large margin

`certify.coverage_factor(n, δ) = √(2 ln(2n/δ))` was registered as the answer to
"what does κ mean when the residual is intrinsic". Level 0 asks the direct
question: over 40 independent replicate runs, does the exhaustive band hold on
**every row** at the TRUE law?

| declared δ | κ | replicates failing | empirical rate | worst max\|resid\|/band |
|---|---|---|---|---|
| 0.50 | 3.15 | 0 / 40 | 0.000 | 0.973 |
| 0.20 | 3.43 | 0 / 40 | 0.000 | 0.894 |
| 0.05 | 3.81 | 0 / 40 | 0.000 | 0.805 |
| 0.01 | 4.21 | 0 / 40 | 0.000 | 0.728 |

**The band holds at every budget, and the margin is real** — even at δ = 0.5 the
worst row across 40 runs reached only 0.97 of its band. So the coverage factor is
sound and *conservative*, which is the direction that costs reach rather than
correctness. The conservatism is expected and its size is accounted for: the
exponential martingale inequality is ~10× looser than an exact Gaussian tail at
κ ≈ 4 (`certify.coverage_budget`), and the band carries three further terms
(the quadrature bound, the quadratic-variation estimator's own error, and the
machine term) on top of the martingale.

## The five tasks, scored by the frozen checker

Every task's verdict was submitted as an **abstention token carrying its partial
determination** — the shape this arc produces — and scored by
`stochcheck.score_task`. S5 and S4 are read off this table.

| task | verdict | CW | covered | informative | abstained-ok | missed | resolved |
|---|---|---|---|---|---|---|---|
| `L0-ou` | `ABSTAIN[structural]` | 0 | 4 | 0 | 0 | 0 | 0 |
| `L0-gbm` | `ABSTAIN[structural]` | 0 | 4 | 1 | 0 | 0 | 1 |
| `L0-ode-obs` | `ABSTAIN[resolution]` | 0 | 2 | 0 | 0 | 2 | 0 |
| `L0-null-noise` | `ABSTAIN[noise]` | 0 | 2 | 0 | 2 | 0 | 0 |
| `L0-null-coarse` | `ABSTAIN[resolution]` | 0 | 4 | 0 | 0 | 0 | 0 |
| **totals** | | **0** | 16 | 1 | 2 | 2 | 1 |

* **S5 holds: zero confident-wrong** — after one was found and fixed (below).
* **S4 holds on both nulls, with the registered reason.** Pure noise abstains on
  `noise` (the zero law certifies — vacuity), and the Δt-unidentifiable null
  abstains on `resolution`, not with a wide interval. The reason match was scored
  against the task's registered `null_reason`: `reason_correct = True` for both.
* Nothing certified, and 16 covered components are honest partial determination
  rather than reach: see "What Level 0 does NOT show".

`exceeded_expectation` fired 6 times — components the registered expectation said
would abstain and which came back with covering bounds instead. Recorded as a
finding about the expectations, never as credit.

## S7 confirmed, where the drift accumulates

| n disjoint windows (GBM) | κ | median q | log₁₀ α |
|---|---|---|---|
| 10 | 3.66 | 0.0019 | −24.6 |
| 20 | 3.84 | 0.0019 | −51.5 |
| 40 | 4.02 | 0.0020 | −105.0 |
| 80 | 4.19 | 0.0021 | −212.3 |

8× the windows buys 8.6× the significance exponent while the band widens 14% and
q barely moves. **Patch count is a strictly winning resource** — the opposite of
the usual multiple-comparisons intuition, and the reason the regime is workable.
On stationary OU the same measurement returns q = 1 and α = 10⁰ at every window
count: vacuous, because the band exceeds the target's range.

## The finding that cost the most: the f-family is not optional

With the plain `dX` form (test function `f = x`) **an OU drift cannot be certified
at any window size**, and this is structural, not tuning. For any 1-D diffusion
with a stationary law, `E[a(X)] = 0` — so the drift's φ-weighted integral is pure
*fluctuation*, of the same scale as the martingale it is competing with. Measured:
|y| / band falls from 0.14 to 0.025 as the window grows from 0.1 to 40. The zero
law certifies everywhere; the honest verdict is vacuity, at every scale.

The generator's `f`-family fixes it. With `f = x²/2` the drift column becomes
`∫φX²dt` and `E[−θX²] = −θ·Var(X) ≠ 0`: the drift accumulates while the martingale
grows as the root of the window. Measured |y| / band then *grows* as √L and crosses
1 near L = 50, matching √(θL/2)/κ. The identifiability condition is

> **θ·L > 2·κ(n, δ)²** — a mean-reverting drift needs a window long against its own
> relaxation time, and no amount of sampling rate substitutes for it.

Confirmed across the window sweep: at θL = 4 and 16 the verdict is
`ABSTAIN[noise]` (vacuous, correctly); at θL = 64 and 128 it is
`ABSTAIN[structural]` with a partial determination.

## S1 is half false, and the reason is the coverage factor itself

**S1 predicted the interval half-width scales as the CLT rate σ/√(N·T).** Measured
separately in the two variables, because they are not the same question:

| rows (window fixed at L = 64) | κ | joint width on the `x` coefficient |
|---|---|---|
| 36 | 3.81 | 5.44 |
| 72 | 3.99 | 3.75 |
| 144 | 4.16 | 3.39 |
| 288 | 4.32 | 3.40 |
| 576 | 4.48 | **3.39** |

| window L (rows fixed at 144) | width | width·√L |
|---|---|---|
| 16 | 6.36 | 25.4 |
| 32 | 5.45 | 30.8 |
| 64 | 3.39 | 27.1 |
| 128 | 3.82 | 43.2 |

**Adding rows stops helping entirely.** From 144 to 576 rows — a 4× increase — the
width does not move (3.39 → 3.39). That half is unambiguous: the row count is held
against one fixed path set, so the five points differ only in how many rows are
added.

Adding *window length* does help, and 1/√L is consistent with the first three
points (width·√L flat at 25–31). The L = 128 row sits well above that trend, and
the sweep varies T alongside L so its four points are different realizations rather
than nested subsets of one — so the **exponent here is indicative, not measured**.
What is solid is the sign: window length buys precision and row count does not.

The mechanism is exactly the coverage factor. An exhaustive per-row band admits a
coefficient only if **every** row accepts it, so the admissible interval is an
intersection of n intervals whose centres scatter by roughly the per-row scale.
The intersection's half-width goes like (κ(n,δ) − √(2 ln n)) ≈ ln(2/δ)/κ times that
scale — and κ(n,δ)² = 2 ln n + 2 ln(2/δ) is built so that the √(2 ln n) term
tracks the extreme of n samples. **The union bound's cost is the averaging gain,
to leading order.** Precision comes only from making each row more informative.

So: S1's T half holds, via window length. **S1's N half is false**, and the same
formula that makes the certificate sound is what makes it false. This is a reach
limit of the exhaustive doctrine under intrinsic noise, and it is not a soundness
problem — the intervals cover throughout.

For scale, the OU MLE for θ has asymptotic sd √(2θ/T) — ±5.6% of θ = 1 from one
trajectory of length 640 — where the certificate's joint bound is ±169% (half of
the 3.39 width). The ~30× gap is attributable: κ ≈ 4 instead of
~1 for a guaranteed rather than a confidence band (4×), windows of 64 rather than
the whole 640 (√10 ≈ 3.2×), and the joint 4-term library rather than a known form
(~5× from collinearity). A guaranteed exhaustive band is expensive, and now the
price is measured rather than guessed.

## S6: realized quadratic variation cannot tell the two noises apart

The `ode_obs` system has **no process noise** — deterministic decay observed with
σ_obs = 10⁻³. Realized quadratic variation nonetheless reports a diffusion, and
reports exactly the predicted spurious value 2σ_obs²/Δt:

| Δt | implied b² | 2σ_obs²/Δt | ratio | true b² |
|---|---|---|---|---|
| 4×10⁻³ | 5.02×10⁻⁴ | 5.00×10⁻⁴ | 1.004 | 0 |
| 2×10⁻³ | 1.003×10⁻³ | 1.00×10⁻³ | 1.003 | 0 |
| 1×10⁻³ | 1.991×10⁻³ | 2.00×10⁻³ | 0.995 | 0 |
| 5×10⁻⁴ | 3.992×10⁻³ | 4.00×10⁻³ | 0.998 | 0 |

Confirmed to under 0.5% at every rate, and it **diverges as Δt → 0** — the
estimator gets worse the better you sample. This is the sharpest possible statement
of why Level 2 exists: separating process noise from measurement noise cannot be
done by the estimator that makes the martingale band honest.

## The confident-wrong: one measured, one mechanism, and the fix

**The first Level 0 run scored one confident-wrong**, on `L0-ode-obs`: the joint
bound on the drift coefficient came back [−1087.5, −902.8] where the truth is
−1.0. Found by the frozen checker, on its first contact with a real producer.

The mechanism is a two-consumer problem, which is the rule the checker's
declarations encode:

> Realized quadratic variation of a noisily-observed path estimates
> `[X] + 2nσ_obs²`. In the **band's scale** that contamination is conservative — a
> wider band loses laws and never admits them. In the **Itô correction on the
> target** it is a systematic offset, and nothing bands it. One measured quantity,
> two consumers, safe for one and unsafe for the other.

Three fixes, all in `lagh/ito.py`:

1. Both quadratic-variation functionals are **debiased** with the declared σ_obs,
   and the debias carries its own error into the band.
2. A row whose raw quadratic variation is more than 50% explained by that
   declaration is **refused** — past that, the process diffusion is a small
   difference of two large numbers.
3. The observational channel's **Gram is built**, so σ_obs enters the band as a
   first-order term. It previously did not: `ItoBand` accepted `sigma_obs` and
   nothing constructed a Gram, so the declaration was a silent no-op. And
   declaring σ_obs to `certify_drift` while `build_rows` never saw it now **raises**
   — a declared error that does not reach the assembler is an error, not a default.

After the fix, `L0-ode-obs` refuses 105 of 144 rows (the declaration explains a
median 99% of their raw quadratic variation), reports `ABSTAIN[resolution]`, and
every bound covers the truth.

## Two defects found in code that predates this run

**`invariant_content` overstated its claim.** Its docstring said the reported range
covered "every law in the declared vocabulary consistent with these observations",
and argued it "cannot weaken the zero-confident-wrong record". Both were wrong: it
ranges over the certifying laws the SEARCH FOUND, and that can exclude a law that
certifies. Measured here — the true drift −5x certifies every row inside its band
while the reported range for that coefficient was [−0.66, 0], six candidates read
and none near the truth. The claim text is corrected, `over` now says what it
ranges across, and `certify.admissible_interval` answers the bounding question
directly by LP over the consistent polytope.

Per the `STRATEGY.md` gating rule, the LP is **scoped to the Itô path** in this
session. The PDE/PDEBench campaigns that predate it keep the machinery they were
run with; whether their reported partial-determination ranges should be re-derived
by LP is a separate question with its own validation, and the honest statement
meanwhile is that those ranges are reports on a certifying set, which is now what
they say.

**A martingale band is not clean data.** Passing `sigma=0` to `discover` alongside
an `eps_model` told the engine the data was clean, with two consequences: the drift
came back as the exact rational `495334·x/606799` — false exactness for a
diffusion-scale coefficient, exactly what `RNOISE_STUDY.md`'s parametric gate
exists to refuse — and the rival search was narrowed, so a law certified that the
wider search does not support. The band is now translated into the relative scale
it corresponds to and declared, putting the exactness and parametric gates in
their noise regime.

## What Level 0 does NOT show

* **The drift's symbolic FORM did not certify** on any Level 0 system with the
  declared 4-term polynomial library. The registered headline — "at σ > 0, certify
  the drift's symbolic form with calibrated intervals" — is **not met at Level 0**.
  The abstains are honest: the LP proves the admissible set genuinely contains
  materially different laws, so this is under-determination and not a search
  failure.
* **What IS produced is partial determination**: sound joint bounds that cover the
  truth on every component, with the `x` term RESOLVED (established as present) on
  GBM at b ≤ 0.02. Given the direction's own title — certified discovery on
  *partially determined* systems — that is the expected output shape, but it should
  be read as reach short of the registered claim, not as the claim met.
* **S2** (the drift/diffusion sampling-rate asymmetry) needs a diffusion claim,
  which is Level 1. The S6 table is consistent with its diffusion half.
* **S7 holds on GBM and is unmeasurable on OU** (above). The prediction is
  confirmed where α means anything and the amendment is registered in the direction
  doc.
* **One more interface gap, found and closed.** The checker silently dropped the
  record on an abstention submission, so it could not score
  "ABSTAIN[structural] + here is what every consistent law agrees on" — the primary
  output of the partial-determination arc. Now the reason is scored, the record's
  components are scored, and the token speaks only for what the record does not
  mention. The field shapes did not change, so the freeze holds.

## Reproduce

```bash
.venv/bin/python experiments/stochastic/run_level0.py     # ~2 min
.venv/bin/pytest tests/test_ito.py tests/test_stochcheck.py
```
