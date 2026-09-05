# Tweezers C4 — the contaminated axis is a late excursion; the noise record rejects an iid floor

**Empirical.** Passive Force 1x's second half now matches additive slow
contamination at the unchanged 8% signature tolerance. Its recorded excess RMS
is **2.402 mV**, conditional on an independent thermal fit, over **8.017–16.034 s**.
The excess covariance stays above 1/e through **0.500 s**. This is a crossing-time
lower bound for that window, not a stationary OU timescale. Block observations
show why: the signal changes level late in the record.

The 100 kHz noise-floor recording supplies a reproducible high-frequency spectral
residual, but **no validated iid sigma_obs**. Its colored residual and acquisition
response defeat the registered white-floor model. A separate held-out decimation
experiment isolates alias accounting. **R-C3 stays open.**

Registration: [TWEEZERS_C4_REGISTRATION.md](TWEEZERS_C4_REGISTRATION.md).
Code: `experiments/tweezers/run_c4.py`; final artifact:
`experiments/results/tweezers_c4.json`. Failed initial reads are retained as
`tweezers_c4_initial.json` and `tweezers_c4_stored_diode.json`. Input byte counts,
SHA-256 hashes, applied/reference calibration items, timestamps, sample rates,
fit/validation halves and numeric diagnostics travel in the artifacts.

![C4 measurements](../experiments/results/tweezers_c4.png)

## 1. A signature with a prediction left over

`instrument.attribute_deviation` adds `(free apparent theta < 1, 1, 1+v)`:
variance estimates v, while short-lag diffusion must remain unchanged. Apparent
theta is deliberately unpredicted; it cannot supply the contaminant's rate.
The old C3 triple **(0.0185, 0.8851, 3.7967) still refuses**, as registered: the
11.5% diffusion deficit does not pass 8%. No bar was widened.

A contaminated record's apparent ACF rate is not its thermal rate. Using it in
the finite-lag OU inversion undercorrects the thermal increments. C4 instead fits
the thermal PSD in the calibration's original band on one half and evaluates
short-lag diffusion on the other. Strides are frozen at 2,3,4,6,8, additionally
restricted to theta_thermal*tau <= 0.3. This is conditional on a clean thermal
fit band, not a new independent validation of the stored calibration physics.

| passive axis / evaluation half | apparent theta ratio | short-lag b² ratio | Var ratio | verdict |
|---|---:|---:|---:|---|
| 1x, first | 1.166 | 1.010 | 0.986 | unattributed |
| 1x, second | 0.0139 | **1.059** | **4.938** | **slow-contaminant**, residual 5.88% |
| 1y, first | 1.215 | 1.018 | 0.989 | unattributed |
| 1y, second | 0.558 | 1.027 | 1.038 | unattributed |

The ACF slopes are themselves unstable across halves, including the control;
we do not relabel them as physical drift changes. Replaying the original C3
triples preserves its full-record control (`consistent`) and both active `drag`
verdicts; the active data and stage are not reanalyzed here.

Two costs are explicit. First, the initial new signature falsely called the
control's 3.8% variance excess contamination. The retained failed artifact
motivated requiring excess **greater than the existing tolerance**, not merely
positive. Second, a stiffness decrease `(c,1,1/c)` also fits the slow signature:
that formerly unique stiffness attribution is now ambiguous. The classifier
refuses rather than inventing a theta constraint to separate identical evidence.
The library's stiffness test now checks that honest loss of identifiability.

## 2. Amplitude is measurable; one stationary timescale is not

Subtract the independent first-half thermal model's **retained covariance** from
the second half's observed covariance. Do not divide the entire inflated variance
by thermal retention: the slow excess has a different frequency weighting.

| second-half quantity, detector volts | measurement |
|---|---:|
| observed variance | 7.2483e-6 V² |
| thermal retained variance | 1.4795e-6 V² |
| recorded excess variance | **5.7687e-6 V²** |
| recorded excess RMS | **2.4018 mV** |
| RMS range using four training-block thermal fits | **2.3817–2.4128 mV** |
| first 1/e crossing of excess covariance | **not observed through 0.499994 s** |
| minimum normalized excess covariance in that lag window | **0.864** |

The RMS range measures sensitivity to the thermal subtraction, **not confidence
coverage or stationary sampling uncertainty**. Variation across evaluation blocks
is much larger: the four excess variances range from negative subtraction noise
to **3.9503e-6 V²**. The eight-block observations show a near-constant mean through
12 s, a transition in 12.025–14.030 s, and a last-block mean **4.921 mV** above
the first block. This supports a late additive excursion/level shift; it does not
identify an electronic, mechanical or biological cause.

The registered stationary interpretation is killed. No single contaminant OU
rate, displacement amplitude, or whole-record stationary variance is reported.
The 0.5 s lower bound concerns the empirical covariance of the stated window;
it is not a lower confidence bound for a latent process parameter. Synthetic
independent slow-OU and slow-sinusoid tests recover the known variance and first
crossing, while retaining the line's negative covariance lobe instead of calling
both signals exponential decay.

## 3. What noise_floor.h5 actually declares

The file contains four trapped-bead force channels, each **1,000,139 samples at
100 kHz**. It is not a dark/no-bead record. The prior `screen_bfp.py` and committed
screen artifact already included Force 1x, so this is not a blind first read.
The [producer tutorial](https://lumicks-pylake.readthedocs.io/en/v1.8.0/tutorial/force_calibration/diode_model.html)
identifies the visible noise floor and describes using the stored diode model
at the record's diagnostic trap power.

The sole calibration item precedes the recording; its own backing is zero,
chi²/dof is 73–101, and its reported corner frequencies are not a truth key for
these samples. Stored trap powers (0.118/0.256 V) differ from the recorded
Diagnostics means (0.362/0.782 V). The independently stored power-dependent model,
`max - delta exp(-rate*mean(power))`, changes the declared diode parameters:

| trap | stored alpha / diode Hz | at recorded power alpha / diode Hz |
|---|---|---|
| 1 | 0.40746 / 14299.6 | **0.44893 / 14829.5** |
| 2 | 0.48532 / 13555.2 | **0.54511 / 14621.9** |

This uses the [producer implementation](https://github.com/lumicks/pylake/blob/main/lumicks/pylake/force_calibration/calibration_models.py),
not a diode fit to held-out residuals. Evaluating at the calibration's recorded
power reproduces its stored diode parameters. The first, stored-power read
remains in its own artifact; updating this provenance does not rescue the
white-floor prediction.

## 4. A measured spectral excess, with its domain and failure to be white

Fit thermal Lorentzian*declared-diode PSD on half 1 in the frozen **100–2300 Hz**
window. Estimate constant additive one-sided PSD N from **30–40 kHz** on half 1;
test mean excess in 20–30, 30–40 and 40–48 kHz on half 2. None of these test
bands fits its own correction. Values below are conditional on that thermal
model and recorded-power diode declaration, in **V²/Hz**.

| axis | training N at 30–40 kHz | half-difference there | held-out excess / N at 20–30, 30–40, 40–48 kHz |
|---|---:|---:|---|
| 1x | 2.0167e-13 | 8.03e-17 | 0.837, 0.999, 0.365 |
| 1y | 1.9343e-13 | 7.03e-16 | 0.805, 1.007, 0.365 |
| 2x | 1.0002e-13 | 4.65e-16 | −0.439, 1.009, −0.034 |
| 2y | 8.6227e-14 | 6.78e-16 | −0.796, 1.016, −0.202 |

The half-differences describe repeatability, not model uncertainty. Reproducible
30–40 kHz excess is useful, but a white additive detector floor alone cannot
explain these spectra. Negative residuals are **kept**, never clipped into
positive noise power: on trap 2, the frozen thermal-response model itself
overshoots parts of the observed spectrum. Additive noise and response loss are
not separately identified over the full band by this record.

Therefore all four `sigma_obs_V` fields are **null**. White-equivalent amplitudes
are retained only as explicitly conditional counterfactuals. The signed integrated
band residual is not called a detector variance. Neither those counterfactuals
nor a new sigma are injected into Ito rows or transferred to the other files.
This is a measured failure of the iid-noise recipe, not a claim that detector
noise can never be identified. A detector-only record or independent transfer
measurement would distinguish the remaining explanations.

## 5. Retention tested on withheld observations, and on changed sampling

The declared thermal*diode response plus that constant floor predicts held-out
increment variances too high: observed/predicted is **0.728–0.901** across four
axes and strides 1,2,4,8,16. Almost all miss the registered 10% bar. The existing
Retention factor learned above the fit ceiling is recorded for diagnosis, but
is **not** used as evidence for these failed independent predictions: clipping a
noisy PSD/model ratio to [0,1] cannot identify a multiplicative response in the
presence of additive noise.

At 100 kHz the geometric Nyquist fraction is about **0.770**, versus **0.768**
on the passive 78.125 kHz record. Diode fractions are **0.544–0.597**, versus
**0.439–0.444** there. These are different beads, trap powers, diode settings and
acquisitions. The difference is not attributable to sampling rate alone, and the
new file does not independently validate the old record's entire Retention model.

A controlled digital intervention does isolate sampling geometry: decimate the
100 kHz held-out half to 50 kHz **without new filtering**, then compare increments
at identical physical lags. They agree within **0.15%**. Predict those increments
from the other half's measured PSD, either retaining its entire native band
(which will alias) or discarding everything above the new Nyquist:

| axis | observed / alias-retaining prediction at 20 µs | observed / alias-discarding prediction | discarded spectral contribution |
|---|---:|---:|---:|
| 1x | 0.9975 | 1.3359 | 25.33% |
| 1y | 1.0002 | 1.3224 | 24.36% |
| 2x | 0.9950 | 1.1551 | 13.86% |
| 2y | 1.0016 | 1.1485 | 12.79% |

Across all four lags, alias-retaining predictions agree within **0.87%**. This
is held-out response repeatability and a test of the sampling operation, not a
second independent analog acquisition. It demonstrates why an alias-free
Nyquist-truncation correction cannot be applied unchanged to unfiltered
decimation. It does not identify the original analog anti-alias transfer or
prove that every native acquisition aliases in this way.

## 6. R-C3 is unchanged

C4 identifies a **passive-record**, window-specific slow excursion and quantifies
where the noise-floor record defeats a white-noise/response model. Neither
measurement supplies the active stage's length unit or a second drive frequency.
The 13% stage-amplitude discrepancy remains **open**, as do its unit and
transfer-function explanations. C4 neither settles that discrepancy nor revises
the C3 active drag estimate. No new stage-based length or drag is computed.

## 7. Reproduction and prediction score

Run `.venv/bin/python experiments/tweezers/run_c4.py` and then
`.venv/bin/python experiments/tweezers/plot_c4.py` (optional figure).
`OPENBLAS_NUM_THREADS=1` avoids excessive thread overhead. No paid proposer,
new dataset download, or certificate is involved.

P1 met: old triple refuses. P2 met only on the contaminated second half; whole-axis
stationarity was not supported. P3 stationary interpretation killed; scoped
amplitude and crossing bound survive. P4 reproducible band excess met, white
floor failed on all four axes. P5 thermal+white response failed; equal-lag
sampling check and subsequent independently trained alias prediction passed.
R-C3 remains open. All results here are **empirical** or explicitly **open**;
none is `proved`.

Validation: `test_instrument_c4.py` 11 passed and `test_instrument.py` 21
passed, each file in its own process. Ruff F,E9 and diff whitespace checks
passed. No full-suite or new zero-wrong campaign count is claimed.
