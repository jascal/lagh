# PDEBench — first pass — registration and results

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

## Next

The registered next steps are the remaining 1-D families
(diffusion–reaction, CFD — the last is non-periodic, so weak-form
certification runs without the spectral verify track and must say so), then the
2-D families using the geometry C3 stage 4 exercised. The FNO-as-proposer
extension registered in `PDEBENCH_READINESS.md` can now be scored against this
pass's ground truth.
