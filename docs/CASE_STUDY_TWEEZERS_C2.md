# Case study: C2 on a trapped bead — the instrument, declared

C1 (`CASE_STUDY_TWEEZERS_C1.md`) measured what a real instrument does to the
stochastic arc's estimators on the LUMICKS C-Trap passive record: the diffusion read
32% of truth, the drift 46% low, and one axis carried a contaminant the calibration
never saw. It listed four things C2 must do. This run does them, in C1's order, on
both records (`experiments/tweezers/run_c2.py`, artifact
`experiments/results/tweezers_c2.json`, 1.8 s). Library: `lagh/instrument.py`,
`ito.build_qv_rows(band_loss=)`, tests `tests/test_instrument.py`.

## 1. Gate axes on the ACF timescale, before anything is claimed

`instrument.axis_gate`: the record's autocorrelation decay θ against the
calibration's own 2πf_c, within a factor 1.25. A timescale is scale-free and is read
at lags beyond the detector's filter, so neither of the instrument's distortions
touches it.

| record | axis | θ_acf / 2πf_c | gate |
|---|---|---|---|
| passive | Force 1x | **0.040** | refused (C1's contaminant, 1.9× thermal variance) |
| passive | Force 1y | 1.008 | passed |
| active | Force 1x (drive removed) | 0.692 | refused |
| active | Force 1y | 0.697 | refused |

Nothing below is claimed on a refused axis. The active record's two axes agree with
each other and disagree with their calibration by the same 30%: that is not a
contaminant on one axis, it is the near-surface record not being the bulk
Ornstein–Uhlenbeck process the gate is calibrated for (frequency-dependent drag; the
calibration fits with "hydrodynamic correction enabled", the ACF reads a plain
exponential). The gate refuses it honestly and the question is left **open**.

## 2. Declare the band loss as a measured input

`instrument.band_loss`: the fraction of the process's quadratic variation the RECORD
can contain, a product of three factors, none fitted alongside the claim:

| factor | what it is | Force 1y |
|---|---|---|
| r_nyquist | geometric — the sampled increments see the process through w(f) = 4 sin²(πf/f_s), whose full-band integral is the true QV; the record holds only \|f\| < f_N | 0.768 |
| r_diode | the instrument's DECLARED detector model α² + (1−α²)/(1+(f/f_diode)²) | 0.444 |
| r_antialias | MEASURED — recorded PSD over the declared model above the calibration's 23 kHz fit ceiling, clipped to [0, 1] | 0.937 |
| **r_total** | | **0.320 ± 1e-4** (two halves of the record) |

C1 attributed the miss to the anti-alias filter. Measured, the anti-alias filter is
6% of it; the detector's own declared plateau (α² = 0.19 at high frequency) is the
bulk, and a quarter of the process's variation lies above Nyquist where no record can
hold it. The truncated denominator was itself a 23% error until it was noticed.

The raw QV reads 0.323 of 2k_BT/γ₀; divided by the declared loss it reads **1.009**
(Force 1x, refused above but computable: 0.325 → 1.022).

## 3. The diffusion, with the loss declared

`build_qv_rows(band_loss=(r, r_se))` rescales every row (y/r, Var/r²) and puts the
loss's own uncertainty in the band at coefficient 1 — a declared input's error, like
the discretization term, never multiplied by κ. Force 1y, 312 rows, 4 holdout
segments:

| | signal / band | b² mid / bulk | interval covers 2k_BT/γ₀ |
|---|---|---|---|
| raw (C1) | 1.42 | 0.328 | **no** |
| loss declared | 3.05 | 1.026 | **yes** |

The C1 coverage miss is closed by a declaration measured from the record and the
instrument's own model. Tag: **empirical**, and still **circular** with respect to
physics on this record — Rd is defined by the thermal motion (C1 §1), so this is the
QV estimator agreeing with the PSD fit once both see the same band, not b² measured
against nature.

## 4. The b²-free drift, scored

`ito.build_lag_rows` (the lagged form of commit 14511a5), one lag at a time,
inverted analytically. Force 1y:

| lag | inverted θ / 2πf_c |
|---|---|
| 16 | 0.968 |
| 32 | 1.057 |

against 0.54 for the single-time form in C1. A scored estimator, not a certificate —
the ridge recorded in 14511a5 stands.

## 5. The active record: a scale that does not pass through the thermal motion

> **SUPERSEDED by `CASE_STUDY_TWEEZERS_C3.md` (2026-09-05).** The drag figure below
> (1.48× bulk) was computed from a drive scale that used the CALIBRATION's corner
> frequency. At f_drive ≪ f_c the driven response determines only the product Rd·f_c,
> so supplying that corner assumes the quantity under test and returns it: measured in
> simulation, the error is exactly the drag ratio. Two scale-free routes put the drag
> at **1.353 ± 0.079 × bulk**, and the gate failure below is a measurement of κ/γ
> rather than a reason to refuse the record. The section is kept as written.

`instrument.drive_scale`: lock-in of the bead's volts and the recorded nanostage
position at the drive line (38.15 Hz, refined below the bin width), with the
fluid-drag transfer |x_bead/x_stage| = f_d/√(f_d²+f_c²) = 0.0197. Force 1x carries
72% of its power at the drive; Force 1y is not driven.

| | value |
|---|---|
| Rd from the drive | **3.755** stage-units / V |
| Rd stored (defined from D, verified identical to 1e-16 on this file too) | 4.568 µm / V |
| ratio | 0.822 |

If the stage reads in µm (the file carries no unit attribute; Bluelake's convention),
this is the first non-circular displacement scale on this instrument, and the ratio
is what a near-surface drag of (4.568/3.755)² = **1.48 × bulk Stokes** would produce
through the passive calibration's definition of Rd. Tag: **empirical**, and **open**:
the axis fails the gate, so no stochastic claim is built on it; the declared-loss QV
with this scale reads 0.51 of bulk (γ_eff/γ₀ ≈ 2.0), which disagrees with 1.48 by
30%, and the file does not record the height above the surface that Faxén's law would
need to adjudicate. C1 §7's item 1 is therefore reached but not closed: the scale is
in hand; the record it came from is not a clean OU at the calibration's timescale.

## 6. Verdicts

| target | verdict |
|---|---|
| Force 1x (passive) | refused by the gate, as in C1 |
| b² (passive, 1y), loss declared | interval covers; empirical, circular w.r.t. physics |
| θ (passive, 1y), b²-free form | 0.97 / 1.06 of 2πf_c, scored |
| Rd (active, 1x) | 3.755 stage-units/V, non-circular; empirical |
| any stochastic claim on the active record | refused by the gate on both axes; open — **reopened and answered in C3**: the common-mode refusal was a measurement of κ/γ, attributed to the drag at 1.353 ± 0.079 × bulk |

Zero confident-wrong. Nothing here is `proved`.

## 7. What was learned that outlives the read

* **A band loss decomposes, and the parts are not where C1 put them.** The
  detector's declared plateau and the sampling band dominate; the anti-alias filter is
  a 6% effect. The instrument had declared the largest term all along.
* **The denominator of a loss is the full-band process, not the record's band.** A
  sampled increment sees every frequency; truncating at Nyquist was a 23% error that
  the simulation could not have shown, because a simulation that filters an
  already-sampled series leaves the aliased part unfiltered. The test now oversamples,
  filters, then decimates — the instrument's order.
* **The gate is doing physics.** Two axes of one record disagreeing with their
  calibration by the same factor is a statement about the record, and the honest
  output is a refusal with the factor attached.
