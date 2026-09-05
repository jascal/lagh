# Case study: C3 on a trapped bead — the gate's ratio was the measurement

C2 (`CASE_STUDY_TWEEZERS_C2.md`) gated both axes of the near-surface active record
OUT: their autocorrelation timescale was 0.70 of the calibration's, so nothing was
claimed on them, and the drag was left open with two readings 30% apart. C3 asks what
that 0.70 is. It is a measurement, the apparatus that was refusing it can be made to
report it, and **C2's drag figure was its own assumption coming back**.

`experiments/tweezers/run_c3.py`, artifact `experiments/results/tweezers_c3.json`
(1.5 s). Library: `lagh/instrument.py`, tests `tests/test_instrument.py` (21).

## 1. What the calibration's numbers actually are

Verified on both files rather than read from documentation: `gamma_0` is **exactly**
6πηa (ratio 1.000000 at six digits, both records), and `kappa` is **defined** as
2π·gamma_0·f_c (ratio 1.000000). So the calibration's timescale 2πf_c is a statement
about a *reference* drag, and a record whose real drag differs cannot agree with it.
The passive record has `Hydrodynamic correction enabled = 0`; the active one has it
`= 1` and a fit range starting at 10 Hz rather than 100 Hz.

## 2. A common-mode gate failure is a measurement, and equipartition is what proves it

Equipartition fixes the position variance from κ and temperature **alone** — the drag
does not appear — while the timescale is κ/γ. So a record whose drag is not the
calibration's reference has healthy variances on every axis and one common timescale
factor. Contamination moves both. That is the discriminator, and one axis cannot
supply it (`instrument.gate_record`, which refuses to call a single off axis anything
but `unresolved`).

| record | axis | θ/2πf_c | equipartition | verdict |
|---|---|---|---|---|
| passive | Force 1x | 0.019 | **3.80** | contaminated — C1's slow contaminant |
| passive | Force 1y | 1.011 | 1.017 | clean (the control) |
| active | Force 1x | 0.701 | 1.012 | — |
| active | Force 1y | 0.706 | 1.001 | — |
| active | *record* | | | **common-mode**, ratio 0.704, axis spread **0.8%** |

The equipartition ratio divides out the record's own **variance retention**, which is
not its quadratic-variation retention: 0.90 against 0.38 on this record. One measured
instrument response, two consumers, and an estimator that used one where it needed the
other would be wrong by their ratio (`instrument.Retention`).

## 3. The gate stops there, on purpose: κ/γ is one number and they are two

C2 read the 0.70 as a drag and quoted 1.48× bulk. That attribution needs a second
measurement, and the second measurement C2 used was **degenerate**.

The bead's response to a stage driven at f_d is |x_bead/x_stage| = f_d/√(f_d²+f_c²),
so at f_d ≪ f_c the lock-in determines only the **product** Rd·f_c. C2 supplied the
corner from the calibration — the quantity under test. Measured in simulation on a
bead whose drag is known:

| true drag | Rd from the declared corner | Rd from the record's corner | truth |
|---|---|---|---|
| 1.00× | 4544 | 4423 | 4568 |
| 1.44× | 3159 | 4421 | 4568 |
| 2.30× | 1981 | 4503 | 4568 |

The declared-corner scale is wrong by **exactly the drag ratio**: it returns the
assumption. `drive_scale` now takes the corner from the record by default and reports
the `degeneracy_factor` a declared corner would have introduced (0.701 on this record).

## 4. Three scale-free ratios name what moved

Read in the **detector's own volts**, so no length is assumed, three ratios against the
calibration's own reported numbers respond differently to the three things that could
be wrong:

| what moved | θ ratio | b² ratio | Var ratio |
|---|---|---|---|
| drag γ → F·γ | 1/F | 1/F | 1 |
| stiffness κ → c·κ | c | 1 | 1/c |
| scale Rd → s·Rd | 1 | 1/s² | 1/s² |

b² comes from the increments (`instrument.realized_diffusion`), each stride corrected
by the retention **at its own weight** 4sin²(πfs/f_s) — so "step past the filter"
becomes a measured correction whose flatness in the stride is the self-check — and is
compared against the calibration's own fitted `D`. Each hypothesis has one free
parameter and one residual left over.

| record / axis | θ | b² | Var | verdict |
|---|---|---|---|---|
| passive 1y | 1.011 | 1.016 | 1.017 | **consistent** — nothing moved (the control) |
| passive 1x | 0.019 | 0.885 | 3.797 | **unattributed** — no single change explains it |
| active 1x | 0.701 | 0.785 | 1.012 | **drag**, F = 1.348, residual 5.8% |
| active 1y | 0.706 | 0.767 | 1.001 | **drag**, F = 1.359, residual 4.2% |

**γ = 1.353 ± 0.079 × bulk Stokes**, from two axes that never share a displacement
scale, with the stiffness and scale hypotheses leaving 5–10× the residual. Faxén's
parallel correction then puts the bead centre at **1.11 µm**, a gap of 0.59 µm above
the glass — physical, and the record is named `near_surface`.

## 5. The one route that needs a length is the one that disagrees

With the corner taken from the record and the nanostage read in µm, the drive route
gives γ = 1.03× bulk and κ = 0.72 of the calibration's — against 1.35 from the two
scale-free routes. Inverting the disagreement rather than splitting it:

| | drive route | implied by the scale-free drag |
|---|---|---|
| Rd | 5.356 µm/V | **4.661 µm/V** — 1.020 of the calibration's 4.568 |
| stage amplitude | 0.527 (recorded) | **0.459** |

So three numbers agree — the two thermal ratios and the calibration's own Rd — and the
odd one out is the **recorded nanostage amplitude**, 13% high, or equivalently the
transfer at the drive frequency. That is a provenance finding of the same family as
C1's applied-vs-derived calibration mix-up, and it is where the remaining
disagreement lives. Tag: **empirical**.

## 6. Verdicts

| target | verdict |
|---|---|
| the active record's κ/γ against the calibration's reference | **0.704 ± 0.008**, scale-free, both axes |
| what moved | **the drag**, on both axes independently, alternatives 5–10× worse |
| γ / bulk Stokes | **1.353 ± 0.079** (empirical) |
| bead centre height (Faxén, parallel) | **1.11 µm**, gap 0.59 µm (empirical) |
| C2's 1.48× | **superseded** — computed from the calibration's own corner, which the degeneracy makes circular |
| the nanostage's recorded amplitude | 13% above what the drag implies — **open** |
| passive Force 1x | contaminated, as in C1 and C2 |

Zero confident-wrong. Nothing here is `proved`.

## 7. What outlives the read

* **A gate that only refuses is throwing away a measurement.** The ratio a gate
  computes to decide pass/fail is a number about the record; what was missing was the
  second observable that says whether it is a fault or a property. Equipartition is
  that observable here because it is the one thing the drag does not enter.
* **A one-parameter estimator is degenerate until something else fixes the other
  parameter, and the tempting thing to fix it with is the quantity under test.** The
  measured cost was a factor of two, and the error returns the assumption — the worst
  shape a defect can have, because it looks like agreement.
* **Attribute by signature, not by assumption.** Three ratios and three hypotheses
  each with one free parameter leaves a residual that can refuse; `unattributed` on the
  contaminated axis and `consistent` on the control are the same machinery working.
* **The route that needs an undeclared input is the route to distrust.** Both thermal
  ratios are scale-free and agree; the drive route needs the stage's unit and is the
  one that disagrees. Localising a disagreement to its least-supported input is more
  useful than averaging it away.

## 8. What would close the stage question

A **second drive frequency** determines Rd and f_c together with no appeal to the ACF,
and their consistency tests both. A **declared stage unit** would settle it outright.
Neither is in this file, and no claim here rests on the stage channel.


## C4 follow-up

`CASE_STUDY_TWEEZERS_C4.md` attributes passive 1x's second-half deviation to
slow additive contamination under an independently fitted thermal model,
while rejecting a stationary whole-record interpretation. Its measured
noise-floor residual does not validate a scalar iid sigma_obs. R-C3 and the
13% stage-amplitude discrepancy remain open; no stage unit was introduced.
