# Direction: error provenance — is this a measurement or a simulation?

**Scoped 2026-07-29** (user direction), out of the PDEBench dev pass, which
supplied both the motivating failure and most of the diagnostics.

Handed a field, this program should produce the best certified symbolic law it
can whether the field came from an instrument or from a solver — and it should
say WHICH it is looking at, because the two carry different errors and the
difference is not cosmetic. It decides which band channel is correct, and
getting that wrong fails in both directions.

## The crux: the error kind chooses the channel

| | observation error | simulation error |
|---|---|---|
| nature | stochastic, one realization | deterministic, reproducible |
| in time | flat | **grows** (PDEBench advection: 0 → 2.75e-2 by t = 2) |
| in wavenumber | white, or an instrument spectrum | **structured** (β_eff 0.70001 → 0.70035 with k) |
| vs the field's own content | independent | **scales with it** (8.9e-4 at k_max = 4 → 2.0e-2 at k_max = 34) |
| across conserved quantities | uniform | **asymmetric** — the scheme's invariants (PDEBench CFD: mass 1e-4, momentum 1e-2) |
| correct band channel | `KAPPA·σ·sqrt(a'Ga)` — **L2, cancels across the realization** | `field_err·‖w·g'‖₁` — **L1, nothing cancels** |

Using L2 for a deterministic error **under-declares** and admits impostors;
using L1 for a stochastic error **over-declares** by roughly √n and loses laws
that are really there. Both the channel and its magnitude are set **by hand**
today, and that hand is what this direction removes.

**Corrected 2026-07-29, and the correction is load-bearing:** this section first
read that PDEBench's ±17% certificates were wide *because the L1 channel was
chosen*. They are not. L1 was the right channel — the error there is
deterministic, as the modified-equation fit below then confirmed — and the width
came from the MAGNITUDE poured into it, over-declared ~3900× (the Route 2 run
below). On the same L1 channel, a declaration scan already reached β ± 0.02% at
3e-5. So PDEBench is evidence for *automating the magnitude*, not
for automating the channel, and it is a mis-example of this section's own thesis.
A live case of the channel error is still wanted; none is in hand.

## The sharp instrument: recover the modified equation

A scheme integrating `u_t + βu_x = 0` does not solve that equation. It solves its
MODIFIED equation,

    u_t + β u_x = c₂ u_xx + c₃ u_xxx + ...

whose leading coefficients are the discretization's signature: c₂ is numerical
dissipation, c₃ numerical dispersion. So the residual of a stated law, regressed
on the next derivatives in the hierarchy, IS that signature — and it is measured
with machinery this program already has, since those derivatives are just more
weak-form columns.

**Measured 2026-07-29 on PDEBench advection (β = 0.7), and it works:** the
residual of `u_t + 0.7u_x = 0` is explained at **84% of its variance by `u_xxx`
alone**, coefficient **1.09e-7**, against 1.5% for `u_xx` and 0.9% for `u_xxxx`.
Purely dispersive, no dissipation — which independently reproduces the spectral
measurement (amplitudes preserved to 1.0000, phase speed rising with k). The two
agree quantitatively: a `c₃u_xxx` term gives phase speed `β + c₃k²`, predicting
`1.088e-7 × (2π·14)² = 8.4e-4` at mode 14 against a spectrally measured 3.4e-4.
Two unrelated diagnostics, same scheme property, same order.

So the honest description of that file is `u_t = −0.7u_x + 1.09e-7 u_xxx`: the
solver's equation, not advection.

### ...and its measured limit, which is the interesting part

Certifying the modified equation does NOT separate it from the physical law on
this data. Scanning the declaration:

| declared field error | plain `u_t = −0.7u_x` | modified (+ c₃u_xxx) | verdict |
|---|---|---|---|
| 1e-4 | 0.071 | 0.038 | both certify |
| 1e-5 | 0.71 | 0.38 | both certify |
| 3e-6 | 2.35 | 1.25 | neither |

(The plain column is the same measurement Route 2 later bisected: it crosses 1 at
**7.06e-6**, and the two runs agree to three digits. What tightens the band is
lowering the declaration to what the weak form actually requires — the modified
vocabulary buys the *diagnosis*, not the tightening.)

The two residual ratios differ by a constant ~1.85×, so they pass and fail
together and no band separates them. The mechanism is structural: the weak form
integrates `u_xxx` against `φ_xxx`, which a smooth bump window damps hard, so the
term's contribution to a patch integral is small — while the 16% of the residual
it does NOT explain is comparable. Parsimony then correctly drops the term
(measured: the engine certifies `u_t = −0.70003u_x` alone at every band where
anything certifies).

What that buys anyway is large: moving from an opaque L1 declaration (0.0275) to
the modified-equation vocabulary at 5e-4 tightened β from **±17% to ±0.3%**, and
at 3e-5 to **±0.02%** with α from 1e-16 to 1e-87. The diagnosis pays even when
the diagnosed term is not itself certifiable.

## What the verdict is a claim ABOUT

The labels `observation` and `simulation` are operational, not ontological, and
the characterizer must say so in its own output. What is measured is STRUCTURE —
deterministic reproducibility, growth in time, monotone drift in wavenumber,
asymmetry across conserved quantities — and "simulation error" is shorthand for
that structure, never a claim about where the numbers came from. Nothing in this
program can distinguish a physical process from a sufficiently good simulation of
one, and it does not need to: a certificate claims that over a stated domain,
within a declared band, a relation holds. That claim is untouched by what the
substrate turns out to be.

The distinction is not idle, though, because the same diagnostic is used in
physics for exactly the substrate question. The measurement made here —
fitting `β_eff(k) = β + c₃k²` and attributing the drift to a `c₃u_xxx` term — is
formally the modified-equation fit that Lorentz-invariance-violation searches run
against gamma-ray burst arrival times (`v(E) = c(1 ± E/E_QG)`), where a
wavenumber-dependent propagation speed would be evidence of a discrete
substrate rather than of a discrete SOLVER. Same logic, same fit, different
subject; the experimental bounds there are severe (E_QG ≳ 10¹⁹ GeV), which is
what makes the question empirical rather than metaphysical.

So the discipline is the one this arc has already been bitten by twice: the
mundane explanation is checked against its own band FIRST. The c₃ = 1.09e-7
measured here is a finite-volume scheme's truncation error, fully accounted for
by the scheme, and reading anything else into it would be the interpretive form
of reading a construction bug as a finding.

## What to build

1. **`characterize_error(fields, coords, law=None, replicates=None)`** — measures
   the diagnostics in the table above and returns an `empirical` verdict, never a
   certificate: `structured-deterministic` / `unstructured-stochastic` / `both` /
   **`undetermined`**, with the evidence and the recommended channel and
   magnitude. The verdict names the STRUCTURE it measured; `simulation` and
   `observation` are the usual causes of each and are reported as the likely
   reading, not as the finding.
2. **UNDETERMINED must be first-class.** A single trajectory with no replicates
   and no reference cannot distinguish a deterministic error that happens to look
   irregular from a stochastic one. Saying so is the correct answer and is the
   same discipline as the state certificate's undetermined modes.
3. **The modified-equation regression** as a reportable diagnostic, with its
   variance-explained and the coefficient — plus the honest note that explaining
   the residual is not the same as certifying the term.
4. ~~**Wire the verdict to the band**~~ — **ATTEMPTED AND REFUSED**, see "The
   automation does NOT happen" below. The free-fit route under-declares by four
   orders and mislabels the kind, because the fit absorbs the error being
   measured. `characterize_rows` now returns `undetermined` and
   `usable_as_declaration: False` without a stated law.

5. **THE ACTUAL NEXT STEP: an INDEPENDENT reference for the magnitude.** The
   channel can be chosen automatically wherever one exists; the whole problem is
   that PDEBench CFD had none. Four candidate routes, cheapest first, none yet
   built:

   * **Multi-resolution differencing.** Where a source ships the SAME simulation
     at two resolutions, their difference on the coarse grid IS the
     discretization error — no model, no closed form, no assumption. Check what
     PDEBench actually ships per family before assuming it is available.
   * ~~**The tightest conservation law as an internal floor.**~~ **RUN AND
     FALSIFIED** — see the Route 2 section below. The observation that motivated
     it stands (PDEBench CFD needed 1e-4 for continuity against 1e-2 for
     momentum, and 1e-5 vs 1e-4 at near-zero viscosity), but on advection, where
     the scheme error is independently known, the tightest relation is not the
     floor (142× σ_rep) and the spread across relations is 1.47× against a scheme
     error 2643× away. It survives only as a source of an independent
     declaration, under the soundness rule recorded there.
   * **An independent high-accuracy solve**, which worked for Burgers (2.44e-3,
     with the reference's own ladder error reported at 9.5e-5) and is blocked
     wherever this program cannot integrate the family — Euler, and any 2-D
     system.
   * **A learned operator**, per the section below: strongest for the
     deterministic/stochastic SPLIT (σ without replicates) and useless for error
     shared across the training distribution.

   Route 2 was the one tried first — measurable today, on data already fetched,
   and the only one needing neither a second resolution nor an integrator. It
   failed as stated. **Route 1 (multi-resolution differencing) is now the one to
   try next**, and the first step is cheap: check what PDEBench actually ships per
   family before assuming a second resolution exists.

### Route 2 RUN 2026-07-29 — H fails as stated, and corrects an earlier error

`experiments/pde/run_conservation_floor.py`, registered before running, scored on
PDEBench advection — the one family whose scheme error is INDEPENDENTLY known
(2.75e-2). Two relations that both hold exactly for the true solution:
`u_t = −βu_x` (the equation the solver integrates) and `(u²)_t = −β(u²)_x`
(derived conservation).

| relation | smallest declaration it needs | in σ_rep |
|---|---|---|
| R1 `u_t = −βu_x` | **7.06e-6** | 142× |
| R2 `(u²)_t = −β(u²)_x` | **1.04e-5** | 210× |

**Both clauses of H are false.** The tightest is 142× σ_rep, so it is NOT the
floor. The spread is **1.47×**, and the loosest sits **2643× below** the known
scheme error, so the spread does NOT estimate it. Route 2 does not measure what
it was designed to measure.

Two qualifications on that verdict, because it was first reported harder than the
data supports:

* **The first run's numbers were grid artifacts.** The decadal scan reported the
  smallest grid POINT that holds, so every requirement was quantized up to the
  next decade and the spread — a ratio of two such numbers — could only ever come
  back a power of ten. "The spread is 10×" was one grid step. The requirement is
  now bisected to where the truth exactly meets its band, which is assumption-free
  (the band is affine in `field_err`, so the bisection is monotone and cheap).
  The corrected spread of 1.47× fails clause 2 *harder*: two relations needing
  the same declaration to within a factor of 1.5 carry even less information
  about the scheme error than a decade of apparent spread did.
* **Clause 1 fails by 1.42×, not decisively.** The threshold fixed before the run
  was 100× σ_rep and the measurement is 142×. That is a fail, and it is also
  inside the original grid's own resolution of its threshold — worth saying,
  because "142× σ_rep, not the floor" reads as a rout and this is a near miss.
  Whether a relation exists on this family that DOES see the floor is open; both
  tried here land in the same narrow band.

**What the run did establish, and it matters more:** a POINTWISE solver-error
measurement is the wrong magnitude for a WEAK-FORM band. The 2.75e-2 is an
*accumulated* pointwise deviation over a whole trajectory; the weak-form residual
is a *local* violation over a patch. They are different quantities and this
program conflated them — every PDEBench certificate was declared at 2.75e-2 when
the weak form only required **7.06e-6**, an over-declaration of **~3900×**. That
is why β came back at ±17% instead of far tighter. It corrects
`CASE_STUDY_PDEBENCH.md`'s numbers, and one of its conclusions with them: the
certificates were sound and needlessly weak, and they were weak for a property of
THIS METHOD, not — as that case study said before it was corrected — a property
of the benchmark. Both corrections are now carried there.

**Route 2 survives in weakened form.** R2 gives an INDEPENDENT, non-circular
declaration (1.04e-5) for certifying R1 — no exact solution, no reference solve,
no second resolution — and it is 2643× tighter than the pointwise number. That is
exactly the situation CFD is in. What it cannot claim is to measure the floor or
the scheme error.

### The re-run, and the second half of the same error (2026-07-29)

Every PDEBench certificate was re-run at the corrected declarations. β ± 17% →
**± 0.0066%** (α 1e-16 → 1e-115), ν ± 25% → **± 0.33%** (α 1e-43 → 1e-106),
truth inside every interval, every forecast still 0 / 205824. Advection took its
band from R2's requirement (1.040e-5, a relation independent of the one being
certified); Burgers, which has no second relation in the registered vocabulary,
took 2× its own scanned requirement and says so.

**And the conflation had a mirror image that the first write-up read as evidence
AGAINST tightening.** `advection_modified_tight`, run by hand at 3e-5, certified
β to ± 0.019% and then FAILED its forecast at 26990 / 205824 points. That looked
like the tight band overreaching. It was the same confusion pointing the other
way: `run_pdebench.py` passed ONE number to both consumers, and the forecast
compares u pointwise against a trajectory integrated from t = 0, so it needs the
*accumulated pointwise* deviation (2.75e-2) exactly as much as the band needs the
*local weak-form* one (7.06e-6). Each was wrong in the other's slot. With
`--field-err` and `--forecast-err` separated, the same law certifies at
± 0.0066% **and** verifies at 0 / 205824.

So the rule generalizes past the band: **a declared error belongs to a consumer,
not to a dataset.** One field can owe two different numbers to two different
checks, and the question to ask of each is what it integrates against.

**The weakened form needs a soundness rule, and this family is the easy case.**
Borrowing R2's requirement to certify R1 is sound here only because R2 happens to
be the LOOSER of the two, by 1.47×. Nothing guarantees that ordering, and the
counterexample is in this same document: PDEBench CFD needed 1e-4 for continuity
against 1e-2 for momentum, so borrowing continuity's requirement to declare
momentum would under-declare by 100× — the impostor-admitting direction. The rule
is therefore: **scan every registered relation, borrow the LOOSEST requirement,
and treat it as a declaration only for relations scanned alongside it.** Even
then it is a lower bound on the error rather than a measurement of it, since all
of them can be blind to a term orthogonal to the whole registered vocabulary.

## Naming

**Observational error** and **pipeline error** are the pair to use. "Pipeline"
because a solver is one, and so is a data-assimilation system, a photometric
reduction, or a counting-and-binning step — the program already says "pipeline
decode" for exactly this in `DIRECTION_PDE.md`. Where the pipeline is
specifically a discretization, **scheme error** is the narrower term. The
measured verdict stays `structured-deterministic` / `unstructured-stochastic`,
because that is what is measured; the two names above are the usual causes.

## A learned operator (FNO) as an instrument for this

Registered here rather than in `PDEBENCH_READINESS.md`, because error provenance
turns out to be the better use for one than equation proposal is. Four roles, in
order of value, and none of them ever touches the band:

1. **A reference where none exists.** The measurement needs something to compare
   against: advection had a closed form, Burgers an independent high-accuracy
   solve, and CFD had NEITHER -- which is why that run could only report the
   declaration required rather than measure the error. A learned operator
   supplies a reference for families this program cannot integrate.
2. **Resolution extrapolation → the discretization error itself.** An FNO is
   parameterized in mode space and is resolution-transferable, so training across
   the resolutions a benchmark ships and extrapolating the learned operator to
   h → 0 estimates the scheme's error directly. Richardson extrapolation in
   operator space, no closed form needed. This is the role that would unblock
   the CFD case.
3. **Symbol probing → the modified equation, better conditioned.** Perturbing a
   trained operator with single Fourier modes at small amplitude and reading the
   linearized response recovers lambda-hat(k) directly, giving beta_eff(k) or
   nu_eff(k). The same quantity the weak-form residual regression measured
   (c3 = 1.09e-7), but fitted over all 10^4 samples rather than 4.
4. **The deterministic/stochastic split, i.e. sigma WITHOUT replicates.** This is
   the case the characterizer would otherwise have to call `undetermined`. A
   learned operator is a deterministic function of its input, so what it predicts
   ON HELD-OUT SAMPLES is reproducible-given-input (pipeline), and the
   irreducible residual is what no deterministic operator can recover
   (observational). Held-out is load-bearing -- enough capacity memorizes noise.

**Limits, which matter more than the promises:**

* It cannot detect error shared by the whole training distribution. If every
  sample carries the same dispersion, the operator learns it as SIGNAL. It
  separates reproducible-given-input from not, never true from false -- and the
  advection dispersion measured here would be invisible to it for exactly that
  reason.
* It carries no declarable bound of its own, so it may inform which CHANNEL to
  use and never the magnitude that goes into it.
* On shocks it rings, so its residual there is its own artifact rather than the
  data's -- which is the regime it would otherwise be most wanted for.

## Results (2026-07-29, `experiments/pde/run_error_provenance.py`)

Built as `lagh/errormodel.py` and scored against ground truth this session
produced independently. P1-P4 met; P5 and P6 are not yet run.

| case | known to be | verdict | evidence |
|---|---|---|---|
| PDEBench advection | pipeline error (dispersive, measured) | **structured-deterministic → L1** | `u_xxx` explains **84.0%** at c₃ = **1.088e-7**; `u_xx` 1.5%, `u_xxxx` 0.9% |
| this program's C1 fields | observational (declared σ) | **unstructured-stochastic → L2** | best term explains **1.1%** |
| no times, no columns, no replicates | — | **undetermined** | nothing measured that could separate them |

The two real cases are separated by a factor of ~76 in variance-explained, so
the discriminator is not marginal. **P4 (the only forbidden direction) holds**:
the deterministic case was never called stochastic. The reverse error is
permitted and remains possible — it merely widens the band.

One honesty property worth keeping: when the tests run and find no structure, the
report says `unstructured-stochastic` but NOTES that this is evidence rather than
a demonstration — a deterministic error orthogonal to every column tested would
look identical. Only replicates settle it.

## The automation does NOT happen, and the measurement says why

Attempted 2026-07-29: wire the characterizer into the runners so the band's
channel and magnitude stop being hand-set. **It fails, in the dangerous
direction, and the failure is the useful result.**

Characterizing a residual needs a law, and the law is what discovery is for. The
apparently law-free route is to take the residual of the BEST FIT over the
registered library. Measured on PDEBench advection, where the solver error is
independently known to be 2.75e-2:

| residual taken from | u_xxx explains | verdict | derived field_err |
|---|---|---|---|
| the STATED law (`u_t = −0.7u_x`) | **84.0%** | structured-deterministic | — |
| the best FIT over the library | **13.6%** | unstructured-stochastic ✗ | **2.6e-6** ✗ |

The fit moves β to −0.70003 to soak up the dispersion's mean effect, leaving a
residual that is both smaller and less structured. So the free-fit route
under-declares by **four orders of magnitude** and mislabels the kind — and
under-declaring is the direction that admits impostors. Auto-setting a band this
way would have made every PDEBench certificate four orders too tight.

`characterize_rows` therefore REFUSES: with no stated law it returns
`undetermined`, marks `usable_as_declaration: False`, and keeps the computed
magnitudes only as labelled lower bounds. The channel choice can be automated
only where an INDEPENDENT reference exists — an exact solution, a separate solve,
or replicates — which is precisely the CFD case where none was available. The
human declaration is not removable by this route, and saying so is better than an
automation that quietly manufactures confident-wrongs.

## Registered predictions

Left in the words they were registered in. The shipped characterizer returns
`structured-deterministic` / `unstructured-stochastic` where these say
`simulation` / `observation` — a rename made after registration (see "Naming"),
and the mapping is one-to-one, so the scoring in the results table above is
unambiguous. Recorded here rather than silently rewritten: editing a registered
prediction to match what was built is how a prediction stops being one.

- **P1.** On our own C1/C2 fields with declared σ (exact analytic solutions plus
  known noise), the characterizer returns `observation`: no higher-derivative
  term explains more than 10% of the residual variance, and the residual is flat
  in t.
- **P2.** On PDEBench advection it returns `simulation`, names dispersion, and
  recovers c₃ within a factor of 2 of 1.09e-7 across independent sample sets.
- **P3.** On a single trajectory with no replicates and no stated law it returns
  `undetermined` rather than guessing.
- **P4.** Zero confident-wrong in the DANGEROUS direction: it never returns
  `observation` for a field whose error is deterministic. The reverse error
  (calling observation error `simulation`) is merely conservative — it widens the
  band — and is permitted.
- **P5.** On PDEBench CFD the modified-equation terms differ between the
  continuity and momentum equations, consistent with the measured 100× asymmetry
  in the declaration each needs.
- **P6 (learned operator).** Trained on a family with declared σ and held out by
  sample, an operator's irreducible residual estimates that σ within a factor of
  2 — recovering the noise scale with no replicates. On PDEBench advection the
  same procedure returns a residual floor far BELOW the measured 2.75e-2
  deviation, because the dispersion is common to the training distribution and is
  learned as signal. Both halves are the prediction: it finds observational error
  and is blind to distribution-wide pipeline error, and that asymmetry is exactly
  why it informs the channel and never the magnitude.

## ANSWERED 2026-07-30: the central question now has ground truth

This document's central question — separating process noise from measurement noise —
is recorded above as scoreable on exactly two cases, with ground truth for the
separation existing nowhere in the repo. The stochastic suite's Level 2 constructed
it, and the answer turned out to be a measurement rather than a declaration.

**The three sources have three distinct STRIDE exponents.** Summing squared increments
over all offsets at lag s,

    Σ_i (u[i+s] − u[i])²  ≈  c + α·s + β·s²

because iid observation error gives `E[(e_{i+s}−e_i)²] = 2σ²` at every lag (constant
in s), a martingale's increment variance grows linearly in the lag, and a
differentiable path's squared increment grows as s². One polynomial fit separates all
three, so σ_obs is measured: median relative error **0.32%** over 45 constructed cases
(`weakform.qv_three_way`, `experiments/stochastic/run_level2_separation.py`).

Three consequences for this direction specifically:

* **The L1-vs-L2 channel choice can be measured rather than reasoned.** The channel
  question this document opens with — is the error unstructured/stochastic (L2) or
  structured/deterministic (L1) — is answered by which stride term dominates, with
  `dominant` reported per fit.
* **A declared σ becomes redundant rather than merely unverified.** The confident-wrong
  the stochastic arc measured came from a σ_obs nobody checked; the measurement
  retires that declaration at its root.
* **And it has a stated regime of validity, with a refusal outside it.** σ_obs is
  recoverable only while its term is at least ~1/10 of the process term. Outside that,
  the fit OVER-estimates σ_obs by up to 300%, which would tighten a band — so the
  boundary is measured (`SEP_MIN_FRAC`) and the fit refuses beyond it. That direction
  matters: over-declaring loses laws, under-declaring admits impostors.

What this does NOT do: it says nothing about a structured error that happens to scale
linearly in the lag, and nothing about autocorrelated observation error — for which
the constant-in-s signature is exactly wrong. The second is registered as the expected
failure on the intended highD read (computer-vision tracking error is correlated
across frames) with a testable remedy: restrict the fit to strides beyond the
correlation length.

## Why this matters beyond PDEBench

Every external data source this program has touched or plans to touch falls on
one side or the other, and several are mixed: Gaia and the exoplanet archive are
instruments with pipelines, ERA5 and PDEBench are model output, traffic data is
an instrument with a discretization (counting and binning). The two-strata
treatment the traffic case study already plans — certifiable conservation
against a conjectured closure — is the same idea one level down. An instrument
that can say which kind of error it is looking at can pick its own channel, and
an instrument that cannot has to be told, once per dataset, by someone who might
be wrong.
