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

Two families in this pass: **1-D advection** (β = 0.7) and **1-D Burgers**
(ν = 0.01). Both certify their true support, both report intervals containing
the truth, both forecast-verify on a sample no stage of the pipeline saw, and
zero confident-wrong.

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
