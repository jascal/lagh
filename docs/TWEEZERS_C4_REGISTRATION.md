# Tweezers C4 — registered 2026-09-05 before new numerical reads

Scope: empirical attribution and instrument measurements. No new certificates,
no stage length unit, no promotion to proved. C3 and instrument.py read first.
The supplied C3 triple (0.019, 0.885, 3.797) is already known. Producer docs
identify noise_floor.h5 as a trapped-bead record with a high-frequency floor,
not a dark detector record. screen_bfp.py previously included its Force 1x;
this is a registered follow-up, not a blind first read.

## 1. Slow additive contamination

Add hypothesis (free apparent theta < 1, short-lag diffusion = 1, variance =
1+v with v>0), retaining tolerance 0.08 and unique-passing-hypothesis rule.
The theta coordinate is unpredicted, NOT a second fitted physical parameter;
this signature alone cannot yield a contaminant timescale. P1: the old C3 triple
still refuses (11.5% diffusion deficit); no tolerance amendment to make it pass.

P2: estimate the THERMAL rate from the existing in-band Retention fit for the
short-lag inversion on passive 1x, rather than using the contaminant-dominated
single-exponential ACF. Freeze strides 2,3,4,6,8 with theta_thermal*tau <= 0.3.
The thermal rate comes from a disjoint first half; test increments on the second.
Predict b²/calibration within 8% and slow-contaminant attribution. Report both
old and revised routes; do not apply this thermal assumption to drag controls.

P3: conditional on attribution, estimate recorded contaminant variance as total
variance minus the thermal model's retained variance (do not divide slow variance
by thermal retention). Measure its covariance by subtracting thermal covariance,
then report first 1/e crossing within 0.5 s; do not label it an OU decay constant.
Report two-half spread as sensitivity, not confidence coverage. Predict positive
excess and timescale > 10 / theta_thermal. Check eight contiguous blocks for
nonstationarity and use their range as a stability diagnostic. If no crossing,
report a lower bound only; if subtraction is negative or block behavior unstable,
refuse a stationary single-amplitude/timescale interpretation. Synthetic independent
slow-OU and slow-line contamination must distinguish amplitude from decay shape;
clean, drag, stiffness, scale and mixed-change controls must not be misattributed.

## 2. Detector noise and sampling rate

Read four noise_floor axes with their own applied calibration and fs metadata.
P4: a positive high-frequency excess above thermal*diode PSD is reproducible in
held-out halves. Fit thermal parameters only in [100,2300] Hz (producer's stated
floor-avoiding window). Estimate a constant one-sided additive PSD from 30–40 kHz
on the first half; test 20–30 and 40–48 kHz on the second half (10% agreement bar
for band-average excess). A passing white floor gives sigma_obs=sqrt(N*fs/2),
conditional on whiteness across the band. Otherwise report band-limited noise
power and refuse a unique iid sigma_obs. No clipping additive excess into a
multiplicative response. No transferring sigma between recordings without matching
acquisition/noise provenance. Noise on a trapped-bead record is model-conditional.

P5: test response without learning the target band from itself. Freeze in-band
thermal fit and diode parameters on half 1, predict half-2 increment moments at
strides 1,2,4,8,16 with the declared diode/sampling response, with and without the
independently estimated floor. Compare held-out measured/model ratios and report
10% agreement/failure; do not tune cutoffs after reading. Quantify r_nyquist and
r_diode at 100 kHz and the existing 78.125 kHz records separately. Different trap,
diode, power or filter settings preclude attribution of their difference to fs
alone. Also decimate 100 kHz by 2 with NO new filtering: covariance at equal
physical lags must agree within 2%, while an unaliased Nyquist-truncation model
may fail. This deliberately tests the Retention sampling assumption; it is not
an independent second acquisition and cannot validate unknown analog response.

## 3. Stage disagreement

R-C3 stays open regardless of P1–P5. Neither an additive contamination measurement
nor a noise PSD declares the nanostage unit or supplies a second drive frequency.
Do not run the old C3 main path that assumes micrometres for the stage.

Artifacts: experiments/results/tweezers_c4.json, reproducible run_c4.py; synthetic
checks in test_instrument.py or a separate per-process instrument C4 test file.
Record input hashes and calibration provenance. Kill criteria: any attribution
requiring threshold relaxation, returning the calibration's assumption as an
independent answer, or reusing a validation band to fit its own correction is
rejected. All failed predictions remain in the artifact and report.

Pre-execution identifiability clarification: the slow signature overlaps an exact
stiffness decrease (c,1,1/c). Both must pass there, so the unique-verdict rule now
refuses that formerly attributed case. Preserve drag/scale/clean controls;
accept lost stiffness attribution as the explicit cost of the wider hypothesis
class. Do not invent a theta constraint to hide this overlap. The added hypothesis
is conditional on positive excess variance and slower apparent theta.

## Measured amendment after first execution (initial artifact retained)

`tweezers_c4_initial.json`: 1x half 2 attributes at b² ratio 1.059, but half 1
has no excess. Its four subblocks range from negative excess to 3.95e-6 V².
A stationary single-amplitude/decay interpretation is killed; report scoped
amplitudes and crossing bounds. Half 2's covariance does not cross 1/e by 0.5 s.

The initial signature also called control 1y's 3.8% excess a contaminant,
although that is below the EXISTING 8% material tolerance. Require Var > 1+tol
before offering the contaminant hypothesis. This is a measured false-positive
fix, not a relaxed diffusion bar. Prediction: 1x half 2 stays attributed; the
control loses the spurious contamination verdict. Both halves remain reported,
including their unstable apparent-ACF fits. No inference of a physical OU rate
from those fits. Eight-block summaries are sensitivity ranges, not confidence
intervals. No change to the noise predictions, whose samples are still unread.

## Noise provenance amendment (stored-diode read retained)

`tweezers_c4_stored_diode.json`: all four white-floor predictions fail; a floor
learned at 30–40 kHz is reproducible there but is not flat over 20–48 kHz. A
single iid sigma is refused. Half-2 moment predictions fail too. The stored
calibration has backing=0, chi²/dof 73–101 and precedes the record by minutes;
its power 0.118/0.256 V differs from diagnostics 0.362/0.782 V. It is not this
record's thermal answer key. The producer specifies the diode model as
max - delta*exp(-rate*mean(power)), applied to the record's Diagnostics trap
power. Source inspected: lumicks/pylake force_calibration/calibration_models.py,
DiodeCalibrationModel / diode_params_from_voltage (GitHub main, 2026-09-05).

Re-evaluate ONLY the declared diode parameters at measured diagnostic power,
then repeat the frozen fit windows, splits and bars. Keep stored-diode results
as failed initial predictions. Predict white-noise refusal persists; no fit of
the diode to validation residuals. Also report each band excess on both halves,
with half-difference as repeatability (not coverage), and aggregate 20–48 kHz
excess power. Do not integrate a rejected white-floor model across the entire
Nyquist band as if it were a measured sigma. For the contaminant amplitude,
report training-block sensitivity of the independent thermal subtraction alongside
the much larger variation across evaluation blocks. These are empirical ranges.

Further fixed diagnostic before execution: since the white-floor moment model
has failed, isolate decimation geometry using the FIRST-half measured PSD itself,
then predict SECOND-half increments with either its full native frequency band
(alias-retaining) or truncation at 25 kHz (alias-discarding). Predict <10% error
for the former and a positive lost-tail contribution for the latter. This tests
held-out response repeatability and alias accounting, not separation of latent
thermal and electronic spectra. Preserve the failed thermal+white prediction.


## Final score

P1 met; P2 met on second-half 1x only; P3 stationary interpretation killed,
windowed amplitude/crossing bound retained. P4 white-floor test failed on four
axes despite reproducible 30–40 kHz residual. P5 thermal+white prediction failed;
equal-lag decimation passed and the subsequently registered alias-retaining
held-out PSD prediction agreed within 0.87%. R-C3 stays open.

11 C4 tests and 21 existing instrument tests pass in separate processes.
Ruff F,E9 and whitespace checks pass. Figure and final numeric artifact generated
from the script. See CASE_STUDY_TWEEZERS_C4.md for domains, measured amendments,
uncertainty limits and retained failed initial artifacts.

## PR #13 review correction — registered before revised execution

The pre-review C4 artifact is preserved in commit ca424ca. Withdraw the claim
that unfiltered decimation independently validates aliasing or Retention._full:
equal-lag increments are a subsample consistency check, and the trained native-
band PSD comparison tests half-to-half spectral repeatability only. The discarded
spectral integral is a counterfactual truncation, not observed above-Nyquist power.

The original 40–48 kHz test overlaps the acquisition rolloff (~43 kHz, identified
in review), while the subtracted model omits that transfer. The previous rejection
cannot determine detector-floor whiteness. New exploratory passband tests are
20–30 and 40–43 kHz, unchanged 10% bar; 30–40 kHz is training-band repeatability
ONLY, never ANDed into validation. Report 43–45, 45–48, 48–50 kHz separately as
rolloff diagnostics. This is a newly registered protocol, not a rescue of the
original prediction. Predict unresolved detector whiteness and null sigma_obs,
regardless of the composite thermal+floor passband test. The transfer is not
independently known, so do not treat Retention.measured (clipped from this noisy
PSD) as an independent correction. The old test-set inclusion is a protocol
violation; disjoint temporal halves are not reuse of the same observations, but
the former 'none of these test bands fits its correction' wording is withdrawn.

Remove theta<1 from admission of the free-theta signature; separate and report
the existing 0.08 variance-materiality floor from the 0.08 residual tolerance.
Test behavior on both sides of theta=1 and the materiality boundary; no hysteresis
is asserted without an estimator uncertainty model. Restore axis_gate on each
half before interpretation. Report apparent/thermal/calibration theta ratios
and half-to-half spread. If apparent theta is unstable beyond the tolerance,
label the slow signature as compatibility, not identified physical cause; keep
window-specific covariance diagnostics, never a promoted law certificate.

Guard missing ACF, empty strides, failed thermal fits and negative sensitivity
subtractions. Validate diode provenance with explicit exceptions under python -O.
Name the excess as a within-window, sample-mean-removed statistic, not population
variance; explicitly report the unbounded/unknown slow-process mean-removal bias.
Show thermal-clipping sensitivity by removing that attenuation as an alternative
subtraction; do not call the narrow thermal-fit range total uncertainty. Reuse
observed covariance across sensitivity models. Predict that these domain limits
survive even if the numeric 2.402 mV summary changes little.

C3 remains a historical frozen artifact from its recorded commit; add an explicit
current-classifier replay artifact rather than overwrite its stage-unit-assuming
analysis or promise byte-identical reexecution with a changed classifier. R-C3
stays open; no new stage unit or drive frequency.


Revised validation: 25 C4 tests and 21 instrument tests passed, each file in
its own process. New checks cover training-band exclusion from validation,
rolloff exclusion, absent ACF/empty strides/failed fits, unresolved sensitivity,
stiffness increases, free-theta admission, mean-removal matching, explicit diode
provenance under python -O, direct voltage conversion, and unresolved plotting.
Ruff F,E9 and diff whitespace checks passed. The revised real-data artifact and
figure were regenerated; no full-suite or new certified-physics claim is made.

## Closed-campaign falsifiability disposition (2026-09-05)

The PR13 withdrawals remain binding. P5 equal-lag subsampling is **WITHDRAWN
as an independent response test** (paired increments are identical). P4's
original stopband statistic is **WITHDRAWN as a whiteness test** (unknown
transfer confounds it). Training-band agreement is repeatability only. No
threshold change or new data analysis rescues these claims. R-C3, detector
whiteness and cross-record sampling causality remain **open**. The complete
prediction inventory and uncovered branches are in FALSIFIABILITY_AUDIT.md.

Artifact reproducibility is numerical: identical structure and nonnumeric
values, finite floats compared with rtol=1e-12, atol=1e-15. Last-bit summation
differences are permitted; this is not a scientific agreement tolerance.
