# Level 1: the diffusion wants the opposite noise, and quadratic variation is its estimator

**Run 2026-07-29/30,** continuing from [`CASE_STUDY_STOCHASTIC_L0.md`](./CASE_STUDY_STOCHASTIC_L0.md).
Machinery: `lagh/ito.py` (`build_rows(diff_names=...)`, `build_qv_rows`,
`certify_diffusion`), `lagh/certify.py` (`admissible_functional`). DEV, by
construction, as the whole suite is.

Level 1's registered targets are *symbolic drift; symbolic diffusion or a reasoned
abstention on it; any almost-sure invariants* on four systems. This document covers
the first two increments: making the diffusion a **claim** rather than a
measurement, what measuring it revealed, and then certifying it the way the
measurement said it had to be done. The 2-D systems, the invariants target and the
remaining nulls are not in it, and §*What is not here* says so precisely.

## The mechanism: b² in the design matrix

At Level 0 the Itô correction was a MEASURED quantity moved to the left-hand side.
Given a library {h_j} for b², it becomes ordinary `dt` columns instead:

    -∫φ' f(X) dt  =  Σ_k c_k ∫φ f'(X) g_k(X) dt
                   + Σ_j d_j · ½∫φ f''(X) h_j(X) dt  +  M

so drift and diffusion are identified **jointly, from the same rows**, and the
columns are ordinary weak-form terms — no new `weakform` capability was needed
beyond what step 3 already landed. Columns are named `drift:<g>` / `diffusion:<h>`,
which is the frozen checker's component vocabulary.

**This also removes the one unsafe consumer of realized quadratic variation.** The
measured correction was where observation-noise contamination became a systematic
offset on the *target* — the Level 0 confident-wrong. With b² in the design matrix,
realized QV only sets the BAND, where contamination is conservative. That is a
soundness improvement obtained by changing what is claimed, not by adding a guard.

## The finding: drift and diffusion want opposite noise intensities

GBM, `dS = μS dt + b S dW`, μ = 0.8, T = 6, L = 1, 80 rows, δ = 0.05. Relative
width = the joint bound's width divided by the true coefficient.

| b | `drift:x` rel. width | `diffusion:x²` rel. width | resolved |
|---|---|---|---|
| 0.02 | **1.5** | 6388 | drift |
| 0.05 | 3.7 | 2417 | — |
| 0.10 | 7.0 | 1181 | — |
| 0.20 | 13.4 | 532 | — |
| 0.50 | 29.4 | 174 | — |
| 1.00 | 20.8 | 47.7 | — |
| 2.00 | 22.1 | **19.7** | — |

Every bound covers the truth at every b. The two columns are monotone in **opposite
directions**, and the reason is elementary once seen:

* the **drift** contributes `μ·∫φf'g dt` against a band `κ√⟨M⟩ ∝ b`, so its
  signal-to-band goes as **1/b** — less noise is better;
* the **diffusion** contributes `½·b²·∫φf''h dt` against the same band, so its
  signal-to-band goes as **b** — *more* noise is better, because the thing being
  measured **is** the noise.

**So there is no single noise intensity at which both are determined**, and on this
system neither is determined at the crossover (both ≈ 20× relative width near
b ≈ 1–2). That is a structural tension between the two Level 1 targets, not a
tuning problem, and it is the sharpest statement this increment produced.

### What it means for S2, which was registered in the other picture

**S2** predicted: shortening Δt at fixed T tightens the *diffusion* interval and not
the drift's; extending T at fixed Δt tightens the *drift* and not the diffusion. In
the weak form the diffusion barely depends on Δt at all — nothing here estimates b
from quadratic variation, so the classic Δt → 0 demand does not apply to it. The
asymmetry S2 names is real but it lives in **which estimator you use**, not in two
intervals from one design:

| | efficient estimator | its precision |
|---|---|---|
| drift | the weak form (this) | improves with horizon T |
| diffusion | realized quadratic variation | relative error √(2/m), m = T/Δt |

At m ≈ 10⁴ increments quadratic variation pins b² to about 1.4%, where the weak
form's joint bound is 2000× the truth. **For the diffusion the weak form is a very
poor estimator and QV is the efficient one; for the drift it is the reverse, since
QV says nothing about drift at all.** S2 should be re-registered as a statement
about estimators before it is scored.

**The design conclusion:** the diffusion's symbolic form should be certified from
quadratic variation, not from the weak-form design matrix. **Built and measured
below** — the increment the measurement chose over the curriculum's ordering.

## The diffusion, certified from quadratic variation

Same weak form, different target. With a state-weight family {w} and a library
{h_j} for b²,

    ∫φ w(X) d[X]  =  Σ_j d_j ∫φ w(X) h_j(X) dt

where the left side is observable as `Σ_i (φw)_i (ΔX_i)²` and the right side is
ordinary `dt` columns. `w` plays exactly the role `f` plays for the drift.
`lagh.ito.build_qv_rows` / `certify_diffusion`.

**Same discipline, and the band's scale is again measured, not declared.** The
estimator's variance is `Var(Σ w(ΔX)²) = 2Σw²(b²Δt)²`, and `E[(ΔX)⁴] = 3(b²Δt)²`,
so `Var = (2/3)Σw²(ΔX)⁴` — a fourth moment of the increments, needing no knowledge
of b. The same property that made the martingale scale measurable.

**The coupling runs both ways, and the second direction is new.** The drift's band
needed the diffusion (via ⟨M⟩); this estimator needs the drift, because
`E[(ΔX)²] = b²Δt + a²Δt²`. That O(Δt) **drift leakage** enters the band as a
deterministic term with coefficient 1, against a declared bound on |a| that must be
verified. Measured share of the band: 0.3% on GBM at b=0.2, 22% on OU at b=1.4 with
`drift_max = 5`.

One caveat stated rather than assumed: the increments `w((ΔX)² − E)` are martingale
differences with **chi-square** tails, not a continuous martingale, so the applicable
inequality is Bernstein's rather than the exponential martingale bound. Its
sub-exponential correction is `O(√(Δt/L))` — negligible when a window holds many
increments — and `QvBand.bernstein_correction` returns its measured size.

### The result: 2600× tighter, and OU's diffusion certifies

Same truth, same data, two estimators:

| | `diffusion:x²` bound (GBM, b = 0.2, truth 0.04) | relative width |
|---|---|---|
| inside the drift's design matrix | [−10.71, +10.78] | 532 |
| from quadratic variation | [**+0.0360, +0.0443**] | **0.207** |

Both cover. Only the second **resolves** (excludes zero), and it is 2600× tighter.
The other two components come back bounded near zero — `diffusion:1` in
[−0.098, +0.071], `diffusion:x` in [−0.045, +0.059] — so the diffusion's form is
effectively determined even where the engine still reports a structural abstain.

**On OU the diffusion CERTIFIES**: b² = 1.9603 against a truth of 1.96, joint bound
[1.855, 2.097], resolved, with `diffusion:x` and `diffusion:x²` bounded inside
±0.084 of zero. That is the first certified diffusion in the arc.

The relative width is now roughly **independent of b** (0.21 at b = 0.2 through 0.33
at b = 0.02, the residual growth being drift leakage), which is the signature of the
QV estimator: relative error √(2/m) with m increments, independent of the noise
intensity. At m = 10⁴ that floor is 1.4%; the measured 21% is that floor times κ ≈ 4,
times two for a two-sided bound, times the {1, x, x²} collinearity.

One cosmetic honesty note: the certified OU law prints as an exact rational
(`459873/234589`). The parametric gate passes it legitimately — no materially
different neighbour certifies — but the coefficient is determined to ±6%, and the
**determination record says so** while the law string does not. Read the record.

### The diffusion's own identifiability condition

OU's diffusion is a constant and certified; GBM's `x²` resolved over a 170× state
range. **CIR** (`dX = θ(m−X)dt + b√X dW`, so `b² = b²x`) is the case where the form
is genuinely state-dependent and the state range is modest, and it gives the
diffusion's counterpart to the drift's `θ·L > 2κ²`:

| b | Feller `2θm/b²` | sd/mean | `diffusion:x` bound (truth b²) | rel. width | resolved |
|---|---|---|---|---|---|
| 0.80 | 3.12 | 0.548 | [−0.044, +1.203] (0.640) | 1.95 | — |
| 1.00 | 2.00 | 0.710 | [**+0.499, +1.362**] (1.000) | 0.86 | **yes** |
| 1.20 | 1.39 | 0.837 | [+0.610, +2.232] (1.440) | 1.13 | yes |
| 1.35 | 1.10 | 0.884 | [+0.931, +2.763] (1.823) | 1.01 | yes |

Every bound covers at every spread. **The form is identifiable in proportion to the
state's RELATIVE SPREAD** — what separates `{1, x, x²}` from one another is how much
of the state axis the process actually visits, and the crossing here sits between
sd/mean 0.55 and 0.71. The same "visited range" theme as the drift's, one level
down: GBM's 170× range makes `x²` unmistakable, CIR's factor-of-two spread barely
separates `x` from a constant.

**A generator bug worth recording, because the instrument would have been blamed.**
The exact CIR transition scales a noncentral chi-square by
`c = 4θ/(b²(1−e^{−θΔ}))`; written with a 2 instead of a 4, the stationary mean comes
out at 2m rather than m. It was caught by printing the mean before trusting any
verdict, and `tests/test_ito.py` now checks both stationary moments. An exact
sampler with a wrong constant is a silent generator defect that looks exactly like
an instrument defect.

### A claim I asserted and the measurement changed

I wrote that the w-family is "what makes a state-dependent b² identifiable at all".
That is only true for a **stationary** process. Measured:

| process | `("1",)` | `("1","x","x²")` |
|---|---|---|
| OU (stationary, every window sees the same state distribution) | 0.700 | **0.242** |
| GBM (state grows exponentially, windows sample different regions) | 0.159 | 0.168 |

So the w-family supplies state variation when the windows cannot; when the process
is non-stationary the windows already provide it and the extra rows only cost a
little κ. The docstring now says the measured thing.

## Two identifiability findings on the 1-D systems

**A bistable system hides its own drift.** The double well `a(x) = x − x³` was the
harder case and not for the reason first suspected. Its coefficients are wide
(width 18 on a truth of 1.0 at b = 1.4, L = 64) and shrinking the library to the
true support `{x, x³}` barely helps (width 16), so it is not library size. Nor is it
coefficient trade-off: the new `admissible_functional` bounds the drift **as a
function**, and `a(x)` is undetermined too (±4.9 at x = ±1, where its true value is
0, against a function whose range is ±6). The reason is structural — **the wells sit
where the drift vanishes, so a bistable process spends its time exactly where its
drift carries least information.** Density and information are anti-correlated by
construction.

**Multiplicative noise can disconnect the state space.** With `b(x) = b·x` the
volatility vanishes at the origin, so the origin is unreachable and every
trajectory is trapped in one well: measured zero well-crossings over T = 400 across
8 trajectories, against ~21 000 crossings for the additive case at comparable
noise. The drift is then determinable only on the visited well, which is the
registered third null ("trajectories do not visit enough of state space") arriving
as a property of a *system* rather than as a null. A domain-qualified claim
(`certify.domain_qualifier`) is the right output shape and is not yet wired here.

## A capability this needed: bounding the law, not the coefficients

`certify.admissible_functional` gives the joint range of any **linear functional**
of the coefficients over the consistent polytope — one LP per query. It exists
because per-component bounds can badly understate what the data determined, and it
is what turns "these four coefficients are each loose" into "the drift function is
in this envelope". Here it happened to confirm the pessimistic reading, which is
worth stating: it was built expecting to find hidden determination and it found
none.

It also matters for the frozen checker's semantics. Interval coverage is scored per
component (`STOCHASTIC_CHECKER.md` §5), so a system whose *function* is determined
while its coefficients are not would score as unresolved. That is the registered
scoring and this increment does not change it — but the case is now known to be
possible, and a functional-claim component kind is the obvious extension if a real
system exhibits it.

## Level 1 scored through the frozen checker

Four systems, drift and diffusion per component, one submission each carrying its
partial determination. `experiments/stochastic/run_level1.py` →
`experiments/results/stochastic_level1.json`.

| task | drift | diffusion | CW | resolved |
|---|---|---|---|---|
| `L1-dw-additive` | `ABSTAIN[noise]` | **CERTIFIED** b² = 1.9536 (truth 1.96) | 0 | `diffusion:1` |
| `L1-dw-multiplicative` | `ABSTAIN[noise]` | **CERTIFIED** b² = 0.4890·x² (truth 0.49) | 0 | — |
| `L1-cir` | `ABSTAIN[noise]` | `ABSTAIN[structural]` | 0 | `diffusion:x` |
| `L1-gbm-mult` | `ABSTAIN[noise]` | `ABSTAIN[parametric]` | 0 | `diffusion:x²` |
| **totals** | | | **0** | 3 |

**S5 holds on Level 1**: zero confident-wrong across nonlinear drift, state-dependent
diffusion, multiplicative noise, and a system whose state space is disconnected.
Every one of the 28 component bounds covers its truth. Two diffusions certify —
**the first certified laws in the stochastic arc.**

Three things this table says that its headline does not:

* **28 covered, only 2 informative.** Twenty-six of the covering bounds are VACUOUS
  at the task's own coefficient scale. The vacuity axis was added during the freeze
  because a hypothetical entrant could game coverage with [−1e9, 1e9]; it turns out
  to be what makes *our own* run readable. Without it the table reads "28/28 covered,
  zero wrong", which sounds like a triumph and means almost nothing.
* **No drift resolved anywhere.** Consistent with everything measured: the wells hide
  the double well's drift, and CIR and GBM need a longer horizon than these
  configurations give.
* **`abstained_correctly` is 0, and that is an interface observation.** Once a
  producer ALWAYS emits partial determination, the record speaks for every component
  and the abstention token speaks for none — so the abstention-correctness axis stops
  being exercised. At Level 0 it fired because the nulls' records were empty. Not
  wrong, but worth knowing before that axis is used to rank anything.

One subtlety worth stating: `L1-dw-multiplicative` **certifies** b² = 0.4890x² —
correct against a truth of 0.49 — while its joint bound on that coefficient,
[−0.324, +1.117], straddles zero and is therefore unresolved. Not a contradiction:
the certificate says *this* law fits every row, the bound says other laws do too,
including some with a zero coefficient. The disconnected state space is why. The
record reports the wider truth, which is the point of having both.

## The declaration that was under-declared, and the fix that made things better

The first version of this run hand-declared the bound on |a| that the diffusion's
band needs. Checked against the truth afterwards, **all four were under-declared** —
`L1-gbm-mult` by 96× (2.0 declared against 191.9 needed). Under-declaring is the
impostor-admitting direction, so that table was not sound as declared, and the one
diffusion certificate in it rested on a violated declaration.

**This was the third instance of one pattern in this arc**: a band input that no code
verifies. The others were `sigma_obs` accepted while nothing built its Gram, and
realized quadratic variation serving two consumers with opposite safety.

The fix is not better numbers — it is removing the declaration. The bound is now
DERIVED from the drift's own admissible envelope via
`certify.admissible_functional`, one LP per query state:

> measured QV → drift band → drift envelope → diffusion band

one direction, no circularity, and with the right failure mode: an undetermined drift
gives a wide envelope and a wider diffusion band, and an unbounded envelope makes the
run decline rather than pick a number. `drift_max` is gone from the config; there is
nothing left to under-declare.

**And the sound version is strictly better than the unsound one.** Naively the
correct envelope is far wider than the numbers it replaced (67, 31, 137, 50312
against 8, 8, 3, 2) and the leakage term is quadratic in it — the first attempt at the
derived bound made every task abstain on vacuity. What rescues it is that the bound
must be **pointwise**: the leakage is `Σ|w_i| a(X_i)² Δt²`, and a process spends its
time where its own drift is SMALL — most of all a bistable one, whose wells sit at
a(x) = 0. Bounding a² by its worst value anywhere over-declares by the square of a
ratio. With the envelope evaluated at the states actually visited, the leakage share
falls to 9–29% of the band and **two diffusions certify where one did before**. Same
distinction `weakform` draws with `field_l1`, and the same lesson: the honest
propagation is also the powerful one.

## Multi-field Itô: Van der Pol, and a component that carries no noise

The migration Level 1 needed. `lagh.ito.ito_terms_nd` emits the multi-dimensional
identity in `weakform.Term`, and `build_nd` assembles it — one row set spanning both
components, which is the same shape `DIRECTION_PDE_SYSTEMS.md` already uses for a
coupled PDE. For dX_i = a_i dt + Σ_j b_ij dW_j and a smooth f of the whole state:

    -∫φ' f(X) dt = Σ_i ∫φ (∂f/∂x_i) a_i dt
                 + ½ Σ_ik ∫φ (∂²f/∂x_i∂x_k) d[X_i, X_k]  +  M

Two things the scalar case hid, both of which needed new machinery:

* **Cross-variation.** `Term.measure` gained `d[u,v]`, integrating against
  `Σ_i w_i Δu_i Δv_i`. The Hessian's off-diagonal entries appear **twice** in Σ_ik,
  so the emitted term carries that factor 2 explicitly rather than leaving a reader
  to notice. Its variance estimator differs from the diagonal's and is stated:
  `Var(du_i du_k) = v_i v_k + c²` while `E[(du_i du_k)²] = v_i v_k + 2c²`, so the
  fourth-moment sum is an upper bound in the cross case (used as is) and the
  diagonal's 2/3 factor makes it exact.
* **A per-field martingale, which collapses.** ⟨M⟩ = ∫φ² Σ_ik ν_i ν_k d[X_i,X_k]
  is `Σ_n (φ_n Σ_i ν_i,n Δu_i,n)²` — the square of the *summed* increment, one
  scalar per sample. No d² loop, no cross terms to enumerate, and the scalar case is
  the same formula with one field.

Validated on **Van der Pol with additive noise**, `dx = y dt`,
`dy = (μ(1−x²)y − x)dt + b dW`: the true law is inside the band on every row, for
both the deterministic and the stochastic test function.

### A component's quadratic variation says whether it carries noise

Van der Pol's `x` component has **no noise**, so its weak form is a deterministic
identity: `−∫φ'x dt = ∫φ y dt` up to quadrature. Measured residual at the truth:
1.8×10⁻⁴. Its martingale band: 0.35 — **2000× too wide**.

The reason is worth having. A differentiable path has zero quadratic variation, and
the realized estimator sees the O(Δt) residue instead: `Σ(Δx)² ≈ ∫y²dt·Δt`. So the
band treats a discretization residue as if it were a martingale variance. Which
makes the diagnostic measurable:

| Δt | qv, `x` (noiseless) | qv, `y` (driven) |
|---|---|---|
| 4×10⁻³ | 2.81×10⁻² | 1.94 |
| 2×10⁻³ | 1.59×10⁻² | 2.04 |
| 1×10⁻³ | 6.18×10⁻³ | 1.60 |

**A component's realized quadratic variation scales as O(Δt) when it carries no
noise and O(1) when it is driven, so halving Δt separates them.** That is the
multi-dimensional form of the error-provenance question — *which components are
noise-driven* — and it is measured rather than declared.

The band is conservative here, never unsound: reach lost, not correctness.

### And it needs no declaration at all — the two parts separate

The Δt-scaling test above needs two grids. There is an equivalent test on ONE, and it
does better than diagnose: it **decomposes**. Summing squared increments over all
offsets at lag s,

    Σ_i (u[i+s] − u[i])²  ≈  α·s  +  β·s²

because a martingale's increment variance grows like s while a differentiable path's
*squared* increment grows like s². So **α is the martingale part and β the smooth
residue**, and two strides suffice. Measured against known truth:

| path | recovered α | true martingale QV |
|---|---|---|
| pure smooth | 2.75×10⁻⁴ | 0 (total QV 0.623) |
| smooth + 0.01·BM | 0.0203 | 0.0200 |
| smooth + 0.05·BM | 0.5013 | 0.500 |
| smooth + 0.5·BM | 49.83 | 50.0 |

Within 0.3% on the mixed paths, and ~2000× below the total for a purely smooth one —
which is exactly the over-declaration factor measured on Van der Pol's x row. **So
the fifth instance of the recurring pattern did not need a declaration after all: the
quantity is measurable.**

Because this *shrinks* a band, two safety rules are built in rather than added later:
the returned scale is `α + κ·se(α)`, an upper estimate and never the point estimate;
and a poor two-term fit **falls back to the total**, which is always sound. The
two-term model is an approximation — a curved path contributes O(s³), a
state-dependent b contributes O(s^1.5) — so the fallback is what stops the tightening
from ever being a claim the data cannot support.

On Van der Pol, `build_nd(martingale_decompose=True)` does what it should and only
that:

| row | band, raw total | band, decomposed | truth covers |
|---|---|---|---|
| `f = x` (component noise-FREE) | 3.14×10⁻¹ | **2.61×10⁻²** (12× tighter) | yes, 43× margin |
| `f = y²/2` (component DRIVEN) | 5.056 | 5.048 (unchanged) | yes |

It tightens where there is a smooth residue to remove and leaves a genuinely driven
row alone. The gain is 12× rather than the full 2000× because the fallback fires on
some patches and the inflation is deliberately conservative — the safety machinery
costing reach, which is the trade it exists to make.

## Van der Pol certified — the arc's first certified DRIFT

`ito.build_rows_nd` completes the step-3 migration: it packs multi-field Itô rows into
the same `ItoRows` the scalar path produces, so `certify_drift` — with its coherence,
significance, parsimony, holdout and admissible-bound machinery — applies to a 2-D
state with no further work. That was the point of putting the vocabulary in
`weakform` rather than growing a second assembler.

**Identifiability first, per equation**, because Van der Pol on its limit cycle is
stationary in the long run and that is what made stationary drifts vacuous:

| test function | L = 4 | L = 16 | L = 64 |
|---|---|---|---|
| `f = x` (component **noise-free**) | 9.0 | **32.2** | 3.9 |
| `f = y²/2` (component driven) | 0.49 | 0.12 | 0.24 |
| `f = x·y` (both) | 0.56 | 0.015 | 0.001 |

(median \|target\|/band). The split is stark and it is the one the QV decomposition
predicts: with the O(Δt) residue removed from the noise-free row's band, that
equation carries 30× more signal than its band, while the driven equation's own
cancellation leaves it vacuous.

So:

> **`dx = y dt` CERTIFIES** — coefficient recovered as 1.0000, joint bound
> [0.9700, 1.0303] (±3%), resolved, α ≤ 10⁻¹⁵, over a declared 6-term bivariate
> library. The arc's first certified **drift**.

The driven equation abstains on vacuity at the same configuration, and every one of
its bounds covers the truth. **In a partially-driven system the noise-free components'
equations are far more determinable than the driven ones** — which is obvious in
hindsight and was not obvious enough to state before the machinery could measure it.

### What the decomposition actually buys here — a claim I had to correct

I expected to find the certificate *depends* on the stride decomposition. It does
not. Measured:

| | certified | signal/band | `x:y` bound | width |
|---|---|---|---|---|
| raw total band | yes | 3.6 | [0.797, 1.202] | 0.405 |
| decomposed | yes | 29.1 | [0.970, 1.030] | **0.060** |

Even the un-decomposed band leaves signal-to-band at 3.6, so the law certifies either
way. What the decomposition buys is **precision**: the same certified law, with its
coefficient pinned to ±3% instead of ±20%. A 6.7× tighter interval, not a certificate
that would otherwise not exist. The test says so in its name.

## What is not here

* ~~Van der Pol with additive noise.~~ **DONE** (above): the multi-field machinery
  works, the identifiability arithmetic is measured per equation, and `dx = y dt`
  certifies. What is still open is Van der Pol inside the SCORED task set — the
  checker's `drift[i]:term` component convention exists for per-component drifts and
  is not yet wired — and the driven equation, which needs a configuration where its
  own cancellation does not dominate. `ito.build_rows` also still has its scalar
  assembler beside `build_nd`; the equality test binds them and collapsing them is
  now tidying, not a capability gap.
* **The invariants target.** Untouched. The checker scores invariants up to affine
  reparametrization already, so this is generator work rather than instrument work.
* ~~A certified diffusion.~~ **DONE** (above): OU's constant b² certifies, GBM's x²
  resolves at 0.21 relative width, and CIR's state-dependent x resolves once the
  state spreads. What is still open is a state-dependent b² whose form **certifies**
  rather than merely resolving — CIR's best relative width is 0.86, so the engine
  correctly declines to crown a form.
* ~~Level 1 scored through the checker.~~ **DONE** (above): S5 holds, zero
  confident-wrong across four systems, two diffusions certified.

## Reproduce

```bash
.venv/bin/python experiments/stochastic/run_level1.py   # ~40 s
.venv/bin/pytest tests/test_ito.py                     # 28 green
```
