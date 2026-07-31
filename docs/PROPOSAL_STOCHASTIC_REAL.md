# Proposal: the stochastic arc's real-data read — a trapped bead

**Drafted 2026-07-30**, replacing highD as the intended real-data read for
`DIRECTION_STOCHASTIC.md`. Registered before any analysis of the data; the
acquisition of one bead file preceded this document, no measurement did.

> **C0 RAN THE SAME DAY AND REFUSED THE READ.** The dataset cannot support the
> diffusion target: its detection noise exceeds the bead's entire thermal motion
> by 338×, and no averaging closes the gap. The scout cost one afternoon and
> killed the read before it was registered, which is what a scout is for. It
> also returned two findings worth more than the read would have been — a
> one-directional guard in the frozen checker, and an amendment to a registered
> prediction. See §9. Sections 2–8 are left as written, pre-measurement.

## 1. Why highD was dropped

`DIRECTION_PDE.md:88` names German highway drone data (highD) as the arc's
real-data case study, and `DIRECTION_STOCHASTIC.md:457` inherits it as the
target Level 2 was scoped for. Two reasons it is the wrong instrument for
*this* direction, one practical and one technical:

* **Access.** highD is granted by manual review against an *official
  university/company address*. That gate is a real risk for an independent
  project, and nothing about the arc should be blocked behind someone else's
  approval queue.
* **25 fps is the wrong Δt for a quadratic variation.** The diffusion target
  is a Δt → 0 limit. At 40 ms per frame, with computer-vision tracking error
  that `DIRECTION_STOCHASTIC.md:486` already registers as autocorrelated, the
  measurement channel dominates the realized QV in precisely the band where
  the QV would be measured. The dataset's own error would be the answer.

highD stays the right data for the **traffic/PDE** case study, where the
target is a weak-form conservation over road–time patches and the counting
error is a declared binning term rather than a QV contaminant. It is not
retired; it is returned to the arc that wants it.

## 2. The data

**Zenodo 3333905** — `doi:10.5281/zenodo.3333905`, Skidmore, Rajasekharan &
Farrell — *calibration of an optical tweezers within the linear Hookian
region*. **CC-0**, direct anonymous HTTP, no form, no account, no
institutional address.

| property | value |
|---|---|
| detector | quadrant photodiode (QPD), 3 channels: ΔX, ΔY, Σ |
| rate | **200 kHz** (Δt = 5 µs), ~60 s per record |
| replicates | **153 beads over 69 days**, 135 files (~1.1 GB HDF5 each) |
| drive | AOD square wave, commanded ±800 and ±500 nm, 80 ms period, ×4 each |
| bead | sulphate-coated fluorescent, nominal radius 2 µm, in buffered saline |
| container | HDF5 with per-dataset unit and provenance annotations |

The authors' own analysis script (`calibrateopticaltweezers.py`, in the record)
documents the layout, so the format is decoded from the producer's code rather
than inferred:

```
/calibration/Stokes_Faxen_coeffiecent        β  [pN·s/nm]   drag, from R, h, η
/calibration/gain                            G  [V/V]
/measurement_datum/
   planned_displacement_of_trapped_bead/<±800|±500>/
       displacement_x, displacement_y, sum_signal_in_light   [V] @ 200 kHz
   mean_background_signal_in_darkness/
       x_signal_darkness, y_signal_darkness, sum_signal_darkness  [V]
/mean_background_corrected_datum/
   optical_tweezer_calibration_parameters_for_experiment/
       QPD_slope_of_line                     [V/V/nm]  volts→nm
   <±800|±500>/ reciprocal_time_constant     1/τ [s⁻¹]
                spring_constant              k  [pN/nm]
```

## 3. The property that decides it: the circularity is broken

> **AMENDED 2026-07-30 by C1** (`CASE_STUDY_TWEEZERS_C1.md` §1). This section is
> right about the *fluorescence* record and right in principle about the C-Trap, but
> the break lives in the **active** record only. On the C-Trap's PASSIVE record the
> displacement scale `Rd` is *defined* as `sqrt((k_B T/gamma_0)/D_volts)` — verified
> identically on 8 calibration items — so b² there is circular by construction. C1
> ran on the passive record and found it. Read §3 as a property of a *driven*
> calibration, never of a passive one.

Optical-trap calibration is normally circular for our purposes — the volts→nm
scale and the trap stiffness are both derived from the *same Brownian
statistics* we would be certifying, by equipartition or by a Lorentzian fit to
the power spectrum. Certifying a diffusion against a scale calibrated from that
diffusion is the FLAME circularity that got reanalysis rejected in
`DIRECTION_PDE.md`.

This dataset does not do that. Three quantities are fixed by routes that never
touch the passive fluctuations:

1. **volts→nm** (`QPD_slope_of_line`) from the **commanded** AOD displacements —
   ±800 and ±500 nm are instrument inputs, an external length standard.
2. **the drag γ** (`Stokes_Faxen_coeffiecent`) from bead radius, measured height
   above the dish and solution viscosity, through the Faxén wall correction.
3. **the stiffness k** (`spring_constant`) from the **driven step response** —
   the reciprocal time constant of the exponential relaxation after a commanded
   square-wave edge, times β.

So the read has an answer key that was produced by a different measurement.

## 4. The physics, in the arc's own vocabulary

A bead in a harmonic trap is an Ornstein–Uhlenbeck process — the Level 0
calibration system, in the world:

    dx = −(k/γ) x dt + √(2k_B T / γ) dW

with, per §3, **k and γ both known externally**. Three consequences:

* the drift coefficient θ = k/γ is exactly the authors' `reciprocal_time_constant`,
  measured by them from the driven response and by us from passive fluctuation —
  a head-to-head against an independent instrument;
* the diffusion b² = 2k_BT/γ, so a certified b² **decodes Boltzmann's constant**,
  `k_B = b²γ / 2T`. This is the arc's analogue of the materials campaign's
  `u` to 9 digits, and it is Perrin's measurement;
* equipartition `Var(x) = k_B T / k` is an invariant relating the two, and the
  invariant-discovery target of Level 1 has a real instance.

## 5. Strata

**Certifiable stratum.** θ and b² from the passive segments; the equipartition
invariant; k_B as the decoded constant. Scored against §3's external values.

**Conjecture stratum (banded, must not certify).**
* The trap is harmonic only in the linear Hookian region the title names — the
  restoring force is anharmonic at large excursion, so a single global θ is a
  conjecture outside it, and the domain qualifier (`certify.domain_qualifier`)
  is the honest form.
* Cross-session stiffness: k varies with alignment and laser power between
  sessions. Any law claimed across beads is banded, never certified.

**Must-not-certify / nulls.**
* `mean_background_signal_in_darkness` — the QPD recorded with **no bead**.
  Real instrument noise, no process. It must certify no drift and no diffusion.
  This is a *measured* null, not a synthetic one, and the arc has not had one.
* A single bead must refuse a cross-bead law, per the multi-solution holdout.
* Sessions the authors annotate as **TLED, signal-to-noise significantly
  decreased** — an external label on the measurement channel.

## 6. Registered predictions

The arc's S-series is registered in `DIRECTION_STOCHASTIC.md`; these are the
real-data R-series, frozen here before measurement.

* **R1 — the three-way separation finds all three terms on a real QPD channel.**
  `qv_three_way`'s `c + α·s + β·s²` should show QPD noise in `c`, Brownian
  motion in `α`, and the driven relaxation in `β`. At 200 kHz the per-sample
  Brownian step is ~1 nm (γ ≈ 3.8e-8 kg/s at R = 2 µm), comparable to QPD
  localization noise, so the prediction is that the run lands in the
  **separable** regime (`c ≥ 0.1·α`) rather than being refused as buried.
  A refusal is an acceptable outcome and must be reported as one.

* **R2 — the registered `s > τ` prediction gets its first real-data test.**
  `DIRECTION_STOCHASTIC.md:486` registers that autocorrelated measurement error
  breaks the constant-in-`s` signature, and predicts that restricting the fit to
  strides beyond the correlation length τ recovers it. A QPD's anti-alias filter
  correlates the detector noise over a few samples, so this dataset supplies the
  violating case. Prediction: `c` fitted from `s ≥ 4` differs systematically from
  `c` fitted from `s ≥ 1`, and the `s ≥ 1` fit is the biased one.

* **R3 — the square-wave transitions masquerade as process noise.** Each
  commanded edge is straddled by `s` increments, so ~750 jumps contribute a term
  **linear in `s`** — the same exponent as a martingale. Prediction: `α`
  measured over the full record exceeds `α` measured over the flat inter-edge
  segments, and the excess scales with jump amplitude (±800 vs ±500). If this
  holds, edge-excision is a required preprocessing step and must be declared.

* **R4 — θ recovered from passive fluctuation covers the authors'
  `reciprocal_time_constant`.** Interval claim; must cover.

* **R5 — b² certifies and k_B is covered.** The certified interval for
  `k_B = b²γ/2T` must contain 1.380649e-23 J/K. A confident-wrong here is the
  most expensive outcome in the campaign and dominates the ranking.

* **R6 — α (significance) is vacuous on the passive segments without the
  f = x²/2 family.** Level 0 measured that a stationary drift makes the per-row
  chance-match saturate. A trapped bead is stationary by construction, so the
  structural limit should reproduce on real data. The driven segments are the
  accumulating regime where S7 is testable.

* **R7 — the dark record certifies nothing.** No drift, no diffusion, and the
  measured `σ_obs` from the dark channel should agree with the `c` term of the
  light channel to within its own spread — the same quantity, two routes.

* **R8 — the TLED-annotated sessions measure a larger `σ_obs`.** An ordering
  prediction against an external label, scoreable without any physics.

## 7. Curriculum

* **C0 — format decode, provenance and the cold three-way probe.** Decode the
  HDF5 against the producer's script; record the external answer key (β, k, 1/τ,
  QPD slope, gain) without using it; measure the **storage floor** (σ_rep — the
  quantization of the stored QPD values, which is itself an iid observation error
  and a confound for `c`); run `qv_three_way` cold on a raw channel and on the
  dark channel. Answers R1/R2/R3 in their diagnostic form and decides whether the
  read is registered at all.
* **C1 — θ and b² certified on the passive segments** of one bead, against §3's
  external values. R4, R5, R6.
* **C2 — the nulls and the holdout.** Dark record, single-bead refusal,
  cross-session banding. R7, R8.
* **C3 — the Δt axis.** Decimate 200 kHz to build the sampling-rate sweep on
  *real* data — `S2`'s drift/diffusion asymmetry, which
  `CASE_STUDY_STOCHASTIC_L1.md` flagged as registered in the wrong picture and
  needing re-registration before it is scored.

## 8. Discipline

As established: predictions frozen in git before measurement; frozen-artifact
acquisition with the exact URL, DOI and SHA-256 recorded; the floor procedure on
raw stored columns; every certificate behind its significance gate or reported as
an abstention; the declaration audit's per-consumer accounting (`drift-band`,
`diffusion-qv`, `observation`) applied to a real instrument for the first time;
misses reported as misses.

**Tag on arrival: `open`.** Nothing here is measured yet. C0 is what moves R1–R3
to `empirical`.

## 9. C0 outcome — MEASURED 2026-07-30

`experiments/tweezers/run_c0.py`, one bead (07-02-2013 #1, sha256 `8ede4bef…`,
4.1 s, 4 × 10.4 M samples). `experiments/results/tweezers_c0.json`.

### The refusal

| quantity | OU says | measured |
|---|---|---|
| thermal rms | 6.39 nm | per-sample detection noise **2159 nm** |
| per-sample Brownian step | 1.03 nm | — (**2104×** below the noise) |
| b² | 2.106e5 nm²/s | `qv_three_way` returned 1.15e10, **5.4e4× high** |
| lag-1 autocorrelation | > 0, decaying over 77 samples | **−0.063** |

Averaging does not rescue it, and this is measured rather than argued. For white
noise `m·Var[block mean]` is flat in `m`; here it **rises 118×** between m = 1 and
m = 5000, so the error is not white, and the block variance never comes within
**3396×** of the OU's own predicted share at any averaging scale. The process
would need ~1e5 samples of averaging against a correlation time of 77 samples.

**Root cause, which §2 should have caught and did not:** detection is
*fluorescence* from a photon-starved bead — the sum channel is 1.91 mV — not
back-focal-plane interferometry. The authors' method recovers signal by averaging
~750 **identical** driven responses. A Brownian path is not repeatable, so the
same instrument that calibrates a trap beautifully cannot deliver a diffusion.

### Two findings that outlive the read

**F1 — the frozen checker's separability guard is one-directional.**
`qv_three_way` reported `separable=True`, `dominant="observation"`, and handed
back a b² wrong by 5.4e4×. `SEP_MIN_FRAC` gates the case where the *observation*
term is buried under the *process* term; when the reverse holds it passes
trivially. Worse, the measured graceful degradation recorded in
`DIRECTION_STOCHASTIC.md` (720× observation noise → α 40% high) was established
under **iid** observation error and does not transfer: here, at only 162×, α is
wrong by four to five orders of magnitude, because the α term is fitting the
noise's own stride structure rather than the process. **A symmetric bar gating
the b² consumer is the missing guard**, and this is the recurring class again —
a declared safeguard that silently does nothing in the direction it was not
calibrated for.

**F2 — amendment to `DIRECTION_STOCHASTIC.md:486` and to R2.** The registered
prediction assumes autocorrelated observation error has a finite correlation
length τ, so restricting the fit to strides `s > τ` recovers `c`. This
instrument's error is **periodic** — ACF +0.249 at lag 5 and +0.243 at lag 10, a
40 kHz line — which has no finite correlation length for `s > τ` to step past.
The `s ≥ 4` fit moved `c` by 1.7%, confirming the mitigation does nothing here.
The prediction must be re-registered as conditional on *decaying* correlation,
with periodic error a separate case needing a different treatment.

### Prediction verdicts

* **R1 — FALSE in the sense that matters.** All three terms are present and the
  fit is clean, but the α term is not the bead. Being in the "separable" regime
  was necessary and nowhere near sufficient.
* **R2 — AMENDED, see F2.**
* **R3 — TRUE.** α over the full record exceeds α over the edge-excised flat
  segments by **1.277×**: the square-wave transitions do contribute a term linear
  in `s` and masquerade as process noise. Edge excision is required preprocessing
  and is now implemented in `stride_sums`.
* **R4, R5 — REFUSED, not missed.** The AR(1) drift estimator returns a negative
  lag-1 autocorrelation, so it declines rather than reporting a number; k_B is
  not attempted. This is the correct behaviour and cost zero confident-wrongs.
* **R6 — not reached.**
* **R7 — UNAVAILABLE.** `mean_background_signal_in_darkness` is a **scalar** (a
  mean, subtracted in real time), not a dark time series. The measured null the
  proposal wanted does not exist in this container.
* **R8 — not reached.**

### Provenance findings

* **`spring_constant` is derived as `(1/τ)·β` to machine precision** (max
  deviation 0.0e+00 across all four amplitudes). The answer key therefore carries
  **one** independent drift measurement, not two — §3's "three externally-fixed
  quantities" is really two plus a restatement.
* Storage floor: LSB 3.216e-6 V → σ_rep 1.85 nm, well below the 2159 nm
  detection noise, so quantization is not the limit here. Recorded for the next
  dataset, where it may be.

### What the next candidate must show

The scout's value is a sharpened filter. A replacement must have:

1. **Back-focal-plane interferometry or equivalent** — sub-nm detection at
   ≥10 kHz. Fluorescence detection is disqualified on this evidence.
2. **A published per-sample noise floor**, or a passive record long enough to
   measure one, *below* the thermal rms — checkable before download from the
   paper's Allan or PSD plot.
3. **A passive (undriven) segment.** This record is driven throughout; the
   thermal-only stretch had to be manufactured by excision.
4. **A genuine dark/blank time series**, not a scalar.

Criteria 1 and 2 are cheap to check from a paper figure, which makes the next
screening round much faster than this one.

### F1 discharged — the symmetric guard, 2026-07-30

`lagh/weakform.py`, `qv_three_way` now returns **`process_resolved`**,
**`process_significance`** and a `process_note`, on the same contract as
`separable`/`separable_note`: the number is still returned, the flag says whether a
caller may read it as a point estimate.

The test introduces **no new constant**. A process term smaller than its own
`KAPPA`-sigma is not resolved:

    process_resolved  ==  alpha > KAPPA * se(alpha)

which self-calibrates to how well the stride design constrains α — a fixed second
ratio bar cannot. `process_scale` is unchanged and remains a sound upper bound for a
band; only the point-estimate reading is gated. The note also states the **domain**
the graceful-degradation figure was measured in (iid), which was the second half of
F1.

Validated in both directions:

* **Costs nothing already working.** `run_level2_separation.py` re-run over its 45
  cases: **24/24** cases that recovered b² within 10% are flagged `resolved`, and
  **0** recovered cases are flagged unresolved. The run's own summary is unchanged
  (median 1.19%, worst 7.75%, 24/24 within 10%).
* **Fires on the record that motivated it.** The tweezers passive channel reports
  `process_significance` **0.145** — α is 7× *smaller* than its own 4σ error — and
  `process_resolved: False`, while the old `separable` flag still reads True.

`tests/test_ito.py::test_the_separation_refuses_when_the_PROCESS_term_is_buried`
locks it, including the assertion that the guard does not fire on the documented
720× graceful-degradation case.

### The replacement candidate — SCREENED AND PASSED, 2026-07-30

**Zenodo 14726586** (`CC-BY-4.0`) — LUMICKS C-Trap, Pylake tutorial force
calibration. Back-focal-plane interferometry at **78125 Hz** (100 kHz on one file).
110 MB total. `experiments/tweezers/screen_bfp.py`,
`experiments/results/tweezers_bfp_screen.json`.

The record set answers all four §9 criteria as separate files —
`passive_calibration.h5`, `near_surface_active_calibration.h5`, `noise_floor.h5`,
`fast_measurement_25.h5` — and the screen is the new guard itself.

| | fluorescence record (refused) | C-Trap passive (passed) |
|---|---|---|
| thermal rms, expected vs observed | 6.4 nm vs **2159 nm** | 5.52 nm vs **5.30 nm** (0.96×) |
| `process_resolved` | **False** (sig 0.145) | **True** (sig 4.7–11.9) |
| b²/truth | **5.4e4** | **0.65–0.83** |
| measured detector floor | — | **0.249 nm** vs 4.09 nm thermal |

**The two records exercise opposite directions of the same guard pair**, which is the
cleanest validation either guard could get: the fluorescence record is
observation-dominated (`separable=True`, `process_resolved=False`); the C-Trap record
is process-dominated, so the *original* `SEP_MIN_FRAC` bar is the binding one
(`separable=False` — σ_obs correctly refused as buried) while the new one passes. Two
real instruments, one on each side.

**The answer key is richer than the last one**: κ, γ₀, f_c, Rd, temperature,
viscosity and bead diameter, each with fitted errors — and, decisively, the
photodiode's own filter parameters, **`alpha` ≈ 0.42 and `f_diode` ≈ 8.1–14.3 kHz**.
The correlated measurement error is not merely present here, it is *published and
parameterised*, so R2's `s > τ` question has a ground-truth correction to be scored
against rather than only a violating case.

The remaining b² shortfall is already diagnosed rather than mysterious: the observed
per-sample increment (0.87 nm) sits *below* the free-Brownian prediction (1.60 nm),
and b²/truth **rises with stride** (0.647 on s ≤ 8, 0.799 on s ≤ 32). That is diode
filtering, going the direction physics requires. Recovering b² through a known
first-order detector filter is the natural C1 and the reason to prefer this dataset:
the failure mode is understood, quantified and independently calibrated.

Companion record **11105579** (60 files, active calibration, dual-trap, bead1–4 with
technical replicates) supplies the multi-trajectory holdout when C2 needs it.

Two cautions carried forward: the trap correlation time is only **25 samples**
(6 on the near-surface active record), so `QV_STRIDES_3`'s max stride of 32 exceeds it
and the ladder must be shortened per record; and the exported force channel's `Rf` may
come from a different calibration item than the κ used to convert it, a ≤5% scale
ambiguity that is irrelevant to the screen and must be resolved before C1.

### What survives here

The **driven** read: a deterministic relaxation observed through a stochastic
channel is Level 0's third rung, and θ has an external answer (1/τ = 2581 ± 128
s⁻¹) this record can be scored against. That is a real-data instance of a rung
the arc has only ever run in simulation. It is a smaller claim than §4's and is
**not** the stochastic read — recorded so the acquisition is not wasted, not
promoted to fill the gap.

