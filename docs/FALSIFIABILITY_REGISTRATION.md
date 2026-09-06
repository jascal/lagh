# Falsifiability, estimator and test selection — prospective registration

Status: **open**. Registered before new witness execution, 2026-09-05.
The C1–C4 campaign is closed. Only synthetic inputs and frozen artifacts are
used below. No new physical attribution or stage-unit inference is permitted.

## A. Audit first

A1 (**open**): each empirical prediction in the four pre-existing registration
files and each decision gate in certify.py/instrument.py has a specified
counter-input. Run witnesses before crediting regression tests. Historical
protocol commitments are audited as provenance, not re-run blind experiments.
A missing executed witness is **no coverage**, not a pass. Algebraic identities
are **tautologies**, not independent validation. Kill: withdraw or explicitly
restate any unfalsifiable empirical prediction in its original registration.

Seed interventions: change the true drag while supplying the assumed corner;
add a high-frequency tone before unfiltered subsampling; apply a stopband to
white noise while omitting it from the subtraction model; add a spectral bump
outside the noise-fitting band. These must defeat the respective scientific
claims, even when the old statistic still passes.

A2 (**open**): exhaustive checks refuse nonfinite observations, nonfinite or
negative bands, and empty domains. Witnesses are y=NaN, epsilon=NaN/Inf/-1,
and zero rows. Missing axis variance or an unreadable axis must not support a
whole-record common-mode claim. Faxen at exactly bulk must return no finite
height. Fixes are reject-only; residual reporting must not change finite-input
verdicts and must reuse the evaluated candidate. Kill any diagnostic which
adds a candidate evaluation or a cause claim.

## B. Timescale

B1 (**open**): fit covariance directly on a fixed lag window, including a
constant slow component and free amplitude, rather than selecting positive
log-ACF samples. Compare to the legacy slope and the PSD fit on synthetic
bead() records. Fixed lags 8–199; covariance model A exp(-theta*t)+C;
no reference theta enters fitting. Score n=32768,131072,524288, clean and
slow-sine contamination (variance twice thermal), with/without detector
filtering. Calibration seeds 20–27; untouched evaluation seeds 40–43.
Report truth error, half disagreement and ACF/PSD disagreement separately.

Kill the universal 8% precision claim if any registered regime's evaluation
90th-percentile relative error or half disagreement exceeds 8%. Then derive
an empirical tolerance by the finite-sample 90% upper calibration quantile of
maximum truth/half errors, separately by regime. Never silently widen the
physical cause classifier: a theta error bound is not a joint three-ratio
error model. Report frozen C3/C4 sensitivity only, no reconstructed estimates.
An estimator which cannot read a record is a failure, not dropped from scoring.

## C. Pilot-selected weak test

C1 (**open**): compare f=x²/2 with pilot selection among x, x²/2 and log(x)
(the latter only on positive pilot paths). Fit drift/diffusion on independent
pilot trajectories and score estimated martingale signal/band there. Freeze
selection before opening certification paths. Use separate seeds 100–103
(pilot) and 200–203 (certification), OU and GBM and a zero-drift control.
No certification residual may select a test. Keep existing bands, coverage
budgets and vocabulary. Report both interval reach and wrong exclusions.
Compare the deterministic reach audit unchanged: state tests cannot improve
ordinary algebraic rows. Kill a general reach claim without an observed gain;
kill selection if certification data enter the choice. A loss is a result.

## D. Serialization

D1 (**open**): C4 artifact numerical reproducibility means recursive identical
keys/types/non-numerical values and finite floats within rtol=1e-12, atol=1e-15.
One-ULP perturbation should pass; a 1e-5 relative perturbation of a unit-sized
value must fail. This defines artifact comparison, not scientific tolerance.

## A1 disposition, before estimator execution

**WITHDRAWN:** universal executed counter-input coverage. The delivered inventory
explicitly names unexecuted protocol/branch witnesses as **no coverage**. A1 is
restated as an inventory plus the bounded executed witness set in
FALSIFIABILITY_AUDIT.md. Historical identities and unavailable provenance are
not retroactively credited as tests. Eight A2 baseline failures have concrete
inputs and are repaired before proceeding to B.

B implementation detail, frozen before synthetic scoring: covariance weights
come from nonoverlapping 2048-sample blocks, centered by the record mean; use
80% empirical covariance plus 20% diagonal regularization (singular otherwise
on short records). Fit amplitude, offset and theta with GLS on lags8–199.
Bounds theta*dt in [1e-6,2], normalized amplitude [0,10], offset [-2,2].
No choice among regularizers or windows is made on evaluation seeds.
Unfiltered controls use bead's exact OU generator and physical truth; filtered
controls call bead(driven=False) with its default diode and anti-alias filter.

## B measured amendment: inadequate calibration count

The initial synthetic run is retained as `timescale_initial.json`. The 8%
universal claim failed (e.g. filtered clean n=524288: half disagreement8.26%;
shorter/contaminated regimes worse). Eight calibration records also cannot
supply a finite split-conformal 90% bound: ceil(.9*(8+1))=9, beyond the sample.
That protocol error is not solved by declaring its maximum a 90% bound.
Before new draws, expand calibration to seeds20–31 (12 records/regime) and
move final evaluation to unused seeds60–63. Keep every estimator setting and
all regimes unchanged. The 90% rank is12; the measured maximum is the resulting
bound under exchangeability within this synthetic regime, not on real records.
Predict that at least one regime still needs more than8%; kill any universal
8% reach claim accordingly. Retain the initial evaluation, but never relabel it
as newly held out. No calibration-derived theta bound will become a joint
attribution tolerance without three-observable calibration.

C implementation detail before execution: pilot least-squares drift on (1,x),
nonnegative least-squares squared increments/dt on (1,x²); candidate score is
|mean(w*a_hat)|/sqrt(mean(w²*b²_hat)). This is a time-averaged proxy restricted
to the three state tests, not the unconstrained optimum with a variable bump.
Use OU(theta=1,b=.5,T=256,dt=.002), GBM(mu=.015,b=.2,T=128,dt=.002),
and null GBM(mu=0, otherwise identical), eight trajectories each. Pilot and
certification sets generated separately. Half-window4000 for OU,8000 for GBM;
vocabulary(1,x,x²), delta=.05. Score certificate, resolved drift:x, interval
width, true-law band misses, and truth exclusions on exactly the same validation
paths for baseline and selected test. Baseline x²/2 only. All thresholds and
bands stay unchanged. Null/nonfinite/constant pilot is a failing selector input;
a wrong drift on synthetic noiseless rows must fail the residual check.
The existing deterministic reach cells use unchanged discovery (selector cannot
consume those rows); report a fresh run separately from historical timings.

## Final B/C dispositions

B1's universal8% precision prediction is **WITHDRAWN**. Twelve independently
calibrated regime tolerances are reported in PRECISION_AND_TEST_SELECTION.md;
some exceed8% substantially. This is outcome(b), with no silent change to the
three-ratio classifier and no real-data reach promotion. Initial results and
the sample-size amendment remain visible.

C1's certified-reach-gain prediction is **WITHDRAWN** for the registered ladder:
0/12 baseline and0/12 selected certifications, no resolved drift-x intervals,
and no truth exclusions in either arm. Selection was possible on pilot data
alone; the scientific gain failed. No window, band or vocabulary is tuned to
rescue it. The restricted proxy is not advertised as the full optimum.

Residual verification amendment before mutation test: retained check evidence
must snapshot its input values. Counter-input: check law1 against y=[1,3], then
mutate caller y to [1,1] and epsilon to100 before requesting the measurement.
The stored verdict still has one miss; the measurement must retain residual2,
epsilon.1 and the original row labels. A changed measurement kills the retained
record. This is an evidence-integrity fix, with no new candidate or verdict.
