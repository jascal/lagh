# PDEBench — DEV — registration and results

> **PDEBench is a DEV target, not a blind read** (decided 2026-07-29, user call,
> during this pass). Everything below is a dev result and is tagged `empirical`;
> no score, ranking or win/lose claim is made or is available to be made. The
> reasons are in "Why this is dev and not a read" below, and they are reasons of
> KIND — the same session that produced these numbers also produced the evidence
> that a read is not the right instrument here.

**Run 2026-07-29**, against the published PDEBench data (DaRUS
`doi:10.18419/darus-2986`), following the pre-flight registered the night before
in `PDEBENCH_READINESS.md`. Runner `experiments/pde/run_pdebench.py`, loader
`experiments/pde/pdebench.py`, extractor `experiments/pde/pdebench_fetch.py`,
results `experiments/results/pdebench.json`.

Four families in this pass, and they produced four different verdicts — which is
the most useful thing about the pass:

| family | verdict |
|---|---|
| 1-D advection (β = 0.7) | **CERTIFIED**, β = 0.70002, forecast verified |
| 1-D Burgers (ν = 0.01) | **CERTIFIED**, both coefficients, forecast verified |
| 1-D reaction–diffusion | **ABSTAIN[noise]** — vacuous: the field is frozen over ~90% of the record |
| 1-D CFD (η = ζ = 0.01) | **ABSTAIN[structural]** — the stated laws hold but 15 / 662 classes fit |
| 2-D Darcy (β = 0.1, steady state) | **ABSTAIN** — no declaration is both true and informative; β = 0.1000 measurable in one phase |
| 1-D CFD transmissive (non-periodic) | **ABSTAIN[structural]** + the verify track now REFUSES instead of silently forecasting |

Zero confident-wrong across all five. Each refusal names a DIFFERENT mechanism —
a swallowed signal, an under-determined one, and a required declaration that
exceeds the signal — and each was reached through the truth check rather than
around it. Darcy also exercised the first non-evolution geometry in the arc.

## What was actually run

| | advection β=0.7 | Burgers ν=0.01 |
|---|---|---|
| source | DaRUS id 255666 | id 281363 |
| grid | 1024 × 201, periodic, L = 1 | same |
| samples | 4 fitted/certified + 1 held out for the forecast | same |
| patch rows | 120 (24 patches rejected) | 143 (1 rejected) |
| declared σ (float32 storage) | 5.0e-8 | 5.0e-8 |
| **declared field error** | **2.75e-2 (measured)** | **2.44e-3 (measured)** |
| law | `u_t = −(157943/225626)·[u_x]` | `u_t = 0.0032124·[u_xx] − 1.00096·[u*u_x]` |
| truth | β = 0.7 | ν/π = 0.0031831, −1 |
| intervals | β ∈ [0.587, 0.830] ✓ | [0.002422, 0.004029] ✓, [−1.0565, −0.9494] ✓ |
| α | ≤ 1e-16 | ≤ 1e-43 |
| signal / band | 2.22 | 2.27 |
| forecast verify | OK, 0 / 205824 | OK, 0 / 205824 |

The recovered β is **0.70002** against a stated 0.7 — agreement to 3e-5, from
patch integrals of the raw field with no differentiated data anywhere in the
pipeline.

## The finding: PDEBench's own numerical error is four orders above its storage noise

The first run **abstained**, and the truth check said why before any verdict was
read: with a band assembled from float32 storage noise alone, the file's own
stated law misses its own band by **1185×**. The pre-flight's second declaration
(“solver error — DECLARED, and it cannot be measured after the fact”) is not a
formality; on this data it is the dominant term.

For advection it *can* be measured, because the stated law has an exact
solution: `u(x, t) = u0(x − βt)`, a circular shift computable spectrally from the
file's own first time slice. Measured across 6 samples
(`experiments/pde/pdebench_solver_error.py`):

| sample | spectral content of IC (k with ≥1e-3 of peak) | deviation at t = 2 | relative |
|---|---|---|---|
| 0 | k ≤ 4 | 7.4e-4 | 8.9e-4 |
| 2 | k ≤ 2 | 8.6e-4 | 8.2e-4 |
| 3 | k ≤ 3 | 3.1e-3 | 2.2e-3 |
| 4 | k ≤ 4 | 6.4e-3 | 3.5e-3 |
| 5 | k ≤ 6 | 5.6e-3 | 7.2e-3 |
| **1** | **k ≤ 34** | **2.75e-2** | **2.0e-2** |

It is **dispersion, not dissipation**: mode amplitudes are preserved to 1.0000
and mass to 1e-7, while the per-mode phase speed runs **0.70001 at k=2 to
0.70035 at k=14** — higher modes travel slightly faster, and the pointwise
deviation accumulates with time and with spectral richness. Against a float32
storage noise of 6e-8, that is four orders of magnitude.

For Burgers there is no closed form, so the same measurement is made against an
**independent high-accuracy solve** of the file's own stated law (spectral
exponential integrator, started from the file's first time slice, its own error
bounded by a substep ladder and reported): deviation **2.44e-3** (1.3e-3
relative) with the reference's ladder error at 9.5e-5 — comfortably below the
deviation being measured, so the measurement stands rather than being refused.

**Consequence for anyone certifying against these files:** the band is set by
that undeclared term, not by the storage precision. It is what holds these
certificates to β ± 17% and ν ± 25% instead of the 1e-4 the same instrument
reaches on its own exact fields (C1). The certificate is honest and weak, and it
is weak for a reason that is a property of the benchmark rather than of the
method.

## What this target actually is: model output, so the exercise is pipeline decode

Framing corrected 2026-07-29 (user observation), and it improves the results
rather than excusing them.

**For PDEBench's own task, the solver error measured below is not an error.**
Surrogates are scored by RMSE against the same solver's output, so systematic
solver error cancels exactly out of their metric — it is part of the target
function they intend to be learned. It is invisible from that side and dominant
from this one because the two sides are asking different questions, not because
the data is careless. Generating 10⁴ trajectories per family at a resolution
where the error would be negligible is not a trade anyone would make.

It is also not *realistic* error, and the distinction is worth keeping: real
measurement error arrives characterized and shipped — an error budget, a
calibration, an uncertainty column — while this arrives undocumented, and its
structure (a phase speed rising monotonically with wavenumber) is a scheme
signature rather than an instrument's. A benchmark of that shape is a poor
substitute for real systematics; if this program wants those, it should inject
them itself, declared and controlled.

The productive reading is the one `DIRECTION_PDE.md` already reached for ERA5,
which it rejected as a discovery target because reanalysis is model output (the
FLAME circularity) while noting it "becomes interesting later, as a
pipeline-decode exercise (certifying which balances the assimilation system
enforces)". PDEBench is model output too, so **that is what these runs are**: not
recovering physics from a benchmark, but decoding what the generating solver
enforces and how tightly. Under that reading the results below are positive
findings rather than complaints:

- continuity needs a 1e-4 declaration and momentum 1e-2 — the scheme enforces
  **mass conservation about 100× more tightly** than the momentum balance;
- advection's effective speed rises 0.70001 → 0.70035 with wavenumber — the
  scheme's numerical dispersion, quantified per mode;
- reaction–diffusion ships pure diffusion at the stated ν, with the reaction
  absent.

Each is a statement about a generating pipeline, and each needed a declared
error model to be sayable at all.

## 1-D CFD: the stated laws hold, and the data does not determine them

`1D_CFD_Rand_Eta0.01_Zeta0.01_periodic` (id 164672), three shipped fields
(density, Vx, pressure) run through `lagh/pdesystem.py` as a two-equation system
— the first PDEBench family needing the system path, and the first with shocks.

**The declared field error is reported, not chosen.** CFD has no closed form and
this program's integrator does not solve Euler, so neither the advection trick
(exact solution) nor the Burgers one (independent solve) is available. Instead
the run scans the declaration and reports the smallest one under which each
stated law sits inside its own band:

| equation | smallest declaration needed |
|---|---|
| continuity `ρ_t = −∂ₓ(ρu)` | **1e-4** |
| momentum `(ρu)_t = −∂ₓ(ρu²+p) + (4η/3+ζ)u_xx` | **1e-2** |

Two equations of the same system, on the same rows, differing by **100×** in how
much declared error their own stated form requires. Read as pipeline decode: the
generating scheme enforces mass conservation about two orders more tightly than
it enforces the momentum balance, which is what a conservative finite-volume
method with an approximated viscous closure would do. A single tuned band would
have hidden that entirely.

**Verdict: ABSTAIN[structural] on both equations**, and it is the correct one.
The stated laws hold comfortably inside the band they need (truth/band 0.003 and
0.507, neither vacuous) — but 15 materially different classes certify for
continuity and **662** for momentum. The data admits the stated law and does not
single it out.

The proximate cause is thinness rather than looseness alone: the resolution gate
rejected **102 of 144 patches** (71%) as unresolved, which is the honest response
to shocks, leaving 42 rows over 3 samples and only **5 certification rows**
against 12 declared features. Five points cannot discriminate 12 columns, so
supports interpolate and coherence correctly finds a crowd. **Large parts of
these files are outside the instrument's reach, and the report says so rather
than certifying on the thin remainder.**

### What the abstain was throwing away

Added 2026-07-29, after the pass. `ABSTAIN[structural]` reported nothing about
the coefficients — and the truth check knew the stated law sat at 0.003 of its
band, so something was determined. `certify.invariant_content` now computes what
EVERY law consistent with the data at the declared band agrees on, and reports it
alongside the abstain. On this same CFD run:

| equation | consistent laws | recovered range | truth |
|---|---|---|---|
| `ρ_t` | 485 | `(ρu)_x ∈ [−1, 0]` | −1 ✓ |
| | | `p_xx`, `ρ_xx`, `p_x` pinned to spans 1e-6 – 6e-5 | ≈ 0 ✓ |
| `(ρu)_t` | 219 | `(ρu²)_x ∈ [−1.06, 0]` | −1 ✓ |
| | | `p_x ∈ [−1.19, 0]` | −1 ✓ |
| | | `u_xx ∈ [−1.27, 0.35]` | 0.0233 ✓ |
| | | `u ∈ [−42, 1.4]` — genuinely unconstrained | — |

**Every true coefficient lies inside its recovered range**, and the terms that
are NOT determined are named rather than implied. The claim is about the
vocabulary, the data and the band — not about nature — so it needs no assumption
that the truth is in the certifying set and cannot weaken zero-confident-wrong:
it is strictly weaker than any member that already certifies, and it is reported
ALONGSIDE the abstain, never as a certificate.

### The defect this found: cost scales with how loose the band is

The run took 1553 s, almost all of it in `certify.coherent`, which clusters the
certifying set by PAIRWISE divergence on the probe box — quadratic in a set that
a loose declaration inflates (662 classes implies a far larger certifying set
behind them). The dev campaigns never hit this because their bands were tight and
a handful of candidates ever certified. **Cost rising with the messiness of the
data is backwards**: the harder the case, the more expensive the refusal, and the
answer at the end is an abstain that the certifying-set size alone would have
predicted. Two fixes worth registering, neither taken yet: short-circuit the
clustering once the certifying set exceeds a reported bound (thousands of
materially-different certifying laws IS a structural abstain), and make the
proposal budget aware of the certification split (5 rows cannot support 8191
enumerated supports, and enumerating them anyway is what generated the crowd).

Also fixed on the way: `PatchEpsilon` re-differentiated every candidate
symbolically — one `sympy.diff` plus one `lambdify` per feature per candidate,
~8000 × 13 for this vocabulary — where every weak-form law is linear in its
columns and the gradient IS the coefficient.

## 2-D Darcy: the first steady-state run, and a third way to fail

`2D_DarcyFlow_beta0.1` (id 133218). Every claim the weak-form arc had made until
now was an EVOLUTION equation — a time-derivative target, time as the last axis,
a verify track that integrates forward. Darcy has no time axis at all: it ships
`nu (10000,128,128)` and `tensor (10000,1,128,128)` with x- and y-coordinates
only, for `−∇·(a∇u) = β`. Steady states, equilibria and constitutive laws are a
large share of real science and this program had never tested whether it could
express such a claim.

**The machinery works unchanged.** `build_nd` over two SPATIAL axes, patches,
weak-form columns, bands and the discovery path all run with no modification —
the geometry was already n-dimensional and nothing in it required the last axis
to be time. That is a capability confirmation worth having.

**The reach boundary, stated precisely.** The general variable-coefficient
equation is NOT in this factory's reach. The library is `∂^α(g(fields))` with g
POINTWISE; `∇·(a∇u)` is divergence form at the outer level but its integrand
`a∇u` pairs a field with a DERIVATIVE of another field. By-parts once leaves
`∫∇φ·a∇u`, still a data derivative; moving it again gives `∫∇·(a∇φ)u`, which
needs `∇a` — a derivative of measured data, exactly what the weak form exists to
avoid. Rearranging does not escape it: `(a u_x)_x = (au)_xx − a_xx u − a_x u_x`.
**A variable coefficient that is itself DATA breaks the arc's central guarantee**,
which bears directly on the registered next step "(c) variable coefficients": it
needs a declared error model for `∇a` or a mixed formulation, and it is not free.

**What is measurable, and it is exact where it should be.** PDEBench's
coefficient is binary (a ∈ {0.1, 1.0}), so wherever a is locally constant the
equation collapses to `∇²u = −β/a`, whose every term is divergence form with a
pointwise g. Measured on patches lying inside one phase, binned by distance to
the nearest interface:

| distance to interface | implied β, a = 1.0 | implied β, a = 0.1 |
|---|---|---|
| 0–2 cells | 0.1250 | 0.1430 |
| 2–4 | 0.1008 | 0.1634 |
| 4–8 | 0.1008 | 0.1350 |
| **≥ 8 (deep interior)** | **0.1000** | 0.1492 |

The high-conductivity interior recovers **β = 0.1000**, matching the filename to
four digits with an IQR of 7e-4. The low-conductivity phase is ~49% off at EVERY
distance, so it is not interface smearing; the natural reading is that `u` is
stored downsampled from a finer solve, which is harmless where u is smooth
(a = 1) and biases second derivatives where it is steep (a = 0.1, curvature 16×
larger).

**And a third distinct way to fail.** Restricting to the a = 1 interior, no
declaration works:

| declared field error | truth/band (max) | truth/band (median) | signal/band |
|---|---|---|---|
| 0 | 172 | 0.701 | 609 |
| 1e-6 | 5.56 | 0.003 | 4.44 |
| 1e-5 | 0.573 | 0.0003 | **0.446** |

At 1e-6 the claim is informative but the stated law misses; at 1e-5 the law holds
but the band has swallowed the signal. **The window between "the law holds" and
"the claim is not vacuous" is empty** — a few outlier patches (median 0.70, max
172) force a declaration larger than the signal itself. That is neither the
vacuity of reaction–diffusion (where the signal was genuinely zero) nor the
under-determination of CFD (where many laws fit); it is a fourth verdict, and the
scan is what makes it sayable.

**A lesson from a bug in the diagnostic, not the instrument.** `a` arrives as
float32 and `float32(0.1) ≠ float64(0.1)`, so a first pass comparing `a == 0.1`
silently kept only the a = 1 patches and reported a clean result for a domain
half the size of the one it claimed. The phase selector now compares with a
tolerance. A silent filter that agrees with your hypothesis is the most dangerous
kind.

## Non-periodic CFD: a capability that applied itself outside its domain

`1D_CFD_Rand_Eta1.e-8_Zeta1.e-8_trans_Train` (id 133155). Run against a
registered expectation for once — `PDEBENCH_READINESS.md` already said
non-periodic runs "get weak-form certification without the verify track, and
must say so". **Scored: it failed.** `check_geometry` returned `ok=True` with no
notes, and the spectral forecast would have run happily on a field whose wrap
seam is 4.6× an ordinary interior step. That is defect #9 and a more serious
class than the previous eight: not a wrong number, but a capability applying
itself where it does not hold and saying nothing.

The discriminator needed care. The RAW wrap gap does not work: on an
endpoint-excluded grid a periodic field's seam is `dx·|u_x|`, and this program's
own C2 fields measure 2.6e-2 – 3.0e-2 there — LARGER than PDEBench advection's
9.7e-3. Expressed in units of an ordinary interior step it separates cleanly:

| field | seam / interior step |
|---|---|
| our C2 heat / Burgers | 0.91 / 0.65 |
| PDEBench advection / Burgers / periodic CFD | 0.74 / 0.14 / 0.16 |
| **PDEBench CFD transmissive** | **4.63** |

`verify()` now refuses with a named refusal that states what remains valid, and
the pre-flight names the condition. Confirmed selective: the transmissive file
refuses, advection still verifies 0/205824, and the C2 campaign is identical
across all 36 entries.

**An unplanned cross-check came with it.** The declaration each equation needs
reverses between the two CFD files:

| | continuity | momentum |
|---|---|---|
| periodic, η = ζ = 0.01 | 1e-4 | **1e-2** |
| transmissive, η = ζ = 1e-8 | 1e-4 | **1e-5** |

At near-zero viscosity momentum needs 1000× LESS. That independently supports
the earlier reading that the periodic file's 100× asymmetry was the **viscous
closure** rather than the momentum balance — a prediction nobody made in advance,
confirmed by a file chosen for an unrelated reason.

Certification itself: 34 rows, 110 of 144 patches rejected (77%, against 71% for
the viscous file — consistent with inviscid shocks), both equations structural.

## Why this is dev and not a read

A blind read freezes SOTA and protocol before download, runs once, and reports
win, lose or mixed. Four things make that impossible here, and three of them
were demonstrated by this pass rather than argued in advance:

1. **The band's dominant term is chosen by the analyst, per file.** The declared
   field error is what holds these certificates to β ± 17%, and it was set after
   the truth check showed σ_rep alone could not work. A protocol in which the
   binding parameter is picked once the answer is in view is a dev loop, however
   carefully each individual step is declared.
2. **Recovering the truth requires the truth.** The measurement that sets the
   field error uses the file's own stated law (exactly for advection, via an
   independent solve otherwise). Using the answer to build the band and then
   reporting recovery of the answer is not a test of discovery.
3. **Every family costs a convention hunt.** β is signed as `u_t + β u_x = 0`;
   Burgers' `Nu0.01` is ν/π, not ν; the reaction–diffusion family did not match
   `u_t = ν u_xx + ρu(1−u)` at its filename's parameters at all (the field
   homogenizes by t ≈ 0.25 and the mean then stops moving, which no logistic
   growth at ρ = 1 does). Those are not discoveries about physics, they are
   discoveries about a file format. Scoring them would be scoring the analyst.
4. **There is no SOTA to be blind against.** PDEBench's leaderboard scores
   forecasting surrogates. Nothing on it answers "does a law certify, and over
   what domain", so there is no frozen number a read could win or lose against.
   The honest comparators (PDE-FIND, WSINDy) are not on that board and would
   have to be run here, which is itself a dev exercise.

The pass is still worth its cost, and what it bought is worth stating plainly:
three real defects in this program's own code, one hard measurement about the
benchmark's data (below), and a working end-to-end path from a published file to
a certificate. Those are dev deliverables and they are exactly what a dev target
is for.

**What a genuine evaluation would need**, if one is ever wanted: a band rule
fixed in advance and applied unchanged to every file (for example, "declare the
field error from an independent solve, refuse where the reference is not clearly
better"), conventions taken from the benchmark's own metadata rather than
reverse-engineered per family, abstention counted as a permitted answer rather
than a failure, and a comparator run under the same rule. That is a different
project from this one, and it should be registered as such rather than grown out
of this pass.

## What this pass does NOT claim

- Not a comparison against PDEBench's leaderboard. That board scores forecasting
  surrogates; this instrument answers a different question (does a law certify,
  and over what domain), so the numbers are not commensurable and no ranking is
  implied.
- Not a claim about the whole benchmark. Two 1-D families, 4+1 samples each. The
  2-D families, the CFD sets and the non-periodic cases are untouched.
- Not a claim that PDEBench is wrong. A 1e-3–1e-2 numerical error is ordinary for
  the schemes and resolutions involved; what is worth stating is that the files
  do not ship it, so anyone who certifies against them at storage precision is
  claiming more than the data supports.

## Three bugs the real data found

All three were in this program, not in the benchmark, and all are fixed:

1. **The forecast band omitted the declared field error** — 58736 of 205824
   points outside while the certified law was right to 3e-5. A declared field
   error belongs in the forecast band as well as the certification band, with
   coefficient 1 (the `hard`-channel convention: a computed bound is not
   multiplied by a coverage factor).
2. **Three corners do not bound a translation family.** The forecast envelope
   sampled `{lo, mid, hi}` of each parameter interval; for advection the
   forecast is not monotone in β, so the pointwise min/max of three shifted
   profiles left 45309 points outside an "envelope" every member of the family
   respected. Intervals are now SAMPLED, and the sample count is reported.
3. **The published files ship a `t-coordinate` one entry longer than the
   tensor** (202 against 201). Truncating quietly would leave every patch
   integral using a time grid the data does not have; the extractor records the
   dropped value.

And one that was not a bug but a real constraint: the verify track's explicit
RK45 is diffusion-limited, so a 512-point forecast at ν = 0.2 never finishes.
The exponential integrator added for this run does it in 0.09 s
(`CASE_STUDY_PDE_C2.md`, amended).

## Method notes worth carrying forward

- **No bulk download.** The published 1-D files are ~8 GB of 10000 samples; the
  instrument needs a handful. The DaRUS backend honours HTTP range requests and
  the datasets are stored contiguous, so `pdebench_fetch.py` reads exactly the
  samples asked for (4.9 MB) and writes an extract carrying its provenance —
  DOI, file id, source shape, sample indices, and any geometry quirk found.
- **The patch family has to be sized to the grid and to the physics.** A family
  sized for the dev campaign's 257×81 fields is wrong for 1024×201, and at
  β = 0.7 a signal crosses a narrow patch inside one time window, which the
  resolution gate correctly throws away. `default_scales` sets half-widths as
  fractions of the grid and caps the time window by the declared propagation
  speed.
- **The truth check earned its place again.** Both refusals in this pass were
  diagnosed by asking whether the TRUE law sits inside its own band before
  reading anything into the abstain.

## 1-D reaction–diffusion: the record is 90% dead, and the reaction is absent

`ReacDiff_Nu0.5_Rho1.0` (id 133177), resolved after the reclassification.
Verdict: **ABSTAIN[noise]** — vacuous — and that is the correct answer, for
reasons that are a property of the file:

- The initial condition is **binary**: min 0.000, max 1.000 per slice.
- ν = 0.5 is **confirmed exactly**, measured per Fourier mode at early times
  (ν_eff = 0.4934 – 0.4996 for k = 1, 2 at t = 0.02–0.05), the same
  mode-resolved method that measured advection's dispersion.
- Diffusion at that strength erases the field almost immediately: by t = 0.10
  the slice is **spatially constant** (min = max = 0.5169), and it then stays at
  that value, unchanged to four decimals, through t = 1.00. **About 90% of the
  time record carries no information at all.**
- **No logistic reaction is present.** At ρ = 1 a uniform field at u = 0.517
  grows to ≈ 0.73 by t = 1; the shipped field does not move. Nor is the deficit
  a units mismatch that a rescaled ρ would fix, because the field is *static*,
  not slow. Fitting the homogenized regime returns ρ = 0.000.
- Consistent with pure diffusion: the reference deviation is 3.4% for
  `u_t = 0.5 u_xx` against 21% for the stated law with its reaction.

The engine's verdict on the assembled rows is vacuity: median |target| = 1.9e-14
against a median band of 37. Nothing can be certified from a frozen field, and
the instrument said so.

### The defect this found in our own discipline

`pdesystem.truth_check` — the gate that exists so a null is never read as a
finding — reported **`truth_certifies: True`** on this data. It is technically
right and practically a false green light: when the band swallows the target,
*every* law sits inside it, so the truth doing so is not evidence about the
truth. The check now tests vacuity first and reports `vacuous`, `signal_to_band`
and a note; `truth_certifies` is False when the target is swallowed. Fixed with
a regression test.

That is the third defect this dev target has found in this program's code, and
the first one in the honesty machinery itself.

## Next, as a dev target

Not a queue of families to score. The useful work here is whatever exercises the
instrument on data it did not generate: the 2-D families (`2D_diff-react` is a
1000-group layout and is a reaction–diffusion SYSTEM in 2-D, which would exercise
C3's stage 2 and stage 4 together), and the non-periodic 1-D CFD sets, where
weak-form certification runs without the spectral verify track and has to say so.
The FNO-as-proposer extension in `PDEBENCH_READINESS.md` remains the more
interesting direction, and it does not need a benchmark score to be worth doing.
