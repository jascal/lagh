# PR14 review revision

**Empirical:** the review found real defects and gaps in the first pass. The
counter-input families were registered before execution; the initial revision
is available at `150c617`. `pr14_review_baseline.json` is reproducible against
that pinned revision, and `pr14_review.json` records the corrected behavior.
No raw tweezers record was opened and no scientific prediction was retuned.

| Finding | Counter-input / action | Disposition |
|---|---|---|
| B1 malformed check domain | 1-D, scalar, 3-D, wrong column count, nonnumeric matrix | Refuse without evaluating the law; no AxisError |
| B2 too-specific Faxen witness | 1, nextafter(1,2), 1+1e-12, 1+1e-9, 1.5 | Stable scaled-series inversion; missing precision no longer means measured height |
| B3 invalid uncertainty | NaN, infinity, negative uncertainty | Explicit refusal before inversion |
| Clean without equipartition | both ratios1, no var_ok | Unresolved; a healthy timescale alone cannot grant clean |
| Write-only tool evidence | failed x against2x through passive and active recover | Both tool refusal payloads expose retained evidence |
| Invalid row labels | short, duplicate, negative, floating row indices | ValueError at check boundary; no deferred silent omission |
| Misleading verify note | exact2x with NaN band | Invalid band note, not form divergence |
| Incomplete GLS refusal schema | zero/short record | Normal metadata keys with null unavailable estimates |
| Rounded artifact comparison | one-ULP neighbors at1.234567890125 round to different12-digit values | Explicit quantized profile; a1e-5 perturbation still fails |
| Other engine refusals | fixed-proposer parametric, structural and wrong integer candidate | Existing checked evidence attached; candidate and domain named |
| Asserted versus measured seeds | declared white spectrum/filter and constructed band ratios | Inputs labeled declarations; verdict flags computed from those inputs; identities identified |
| Repeated residual / copies | mutate caller arrays after check | One subtraction; retain residual, band and row-label snapshots, not redundant observation/prediction copies |
| Dead empty-domain disjunct | early empty-y return | Redundant later branch removed |

## Faxen: a neighborhood and an uncertainty contract

**Empirical limitation of the original audit:** the witness at exactly1 did
not cover the open interval between1 and the old root bracket. A too-specific
witness can make an incomplete fix look complete. The revised witness sweeps
both that neighborhood and uncertainty states, rather than adding a special
case for the review's1e-12 input.

**Model-conditional calculation:** write delta=(F-1)/F and lambda=delta*z.
The root equation becomes

```
z*(9/16 - lambda²/8 + 45*lambda³/256 + lambda⁴/16) = 1.
```

It remains well scaled near bulk. Over lambda<=.6, delta/lambda is greater
than1/2, so z can be bracketed in[0,2]; the old hard lower limit lambda=1e-9
is unnecessary. No guessed physical-height cutoff or instrument uncertainty
is introduced. The numerical tests also check the independent leading-order
relation h*(F-1)/radius ->9/16 near bulk.

**Open without uncertainty:** a ratio arbitrarily close to1 can legitimately
have an enormous *conditional mathematical inverse*. It cannot establish a
measured height without precision information. For *every* ratio, omitted
`drag_ratio_se` now gives `height_um=None`, `height_determined=False`, and the
conditional calculation in `conditional_height_um` when it exists. This is an
intentional tightening of the API, not a claim that all previous outputs stay
unchanged. Exact synthetic inversion tests explicitly pass `drag_ratio_se=0`.
Finite supplied uncertainty reaching bulk gives a lower bound only; intervals
outside the series domain do not yield a determined height. No added confidence
level is attached to the caller's uncertainty envelope.

## Two corrections to the review interpretations

**Empirical artifact already committed.** The full rerun is
[`experiments/results/test_selection_reach.json`](../experiments/results/test_selection_reach.json),
introduced by commit150c617, with runner `experiments/run_selection_reach.py`.
It contains all36 cells and35 certifications, with no verdict changes relative
to July's `reach_audit.json`. The old artifact was deliberately left untouched.
This is evidence from the initial PR run, not a claim that the entire reach
audit was rerun after these review fixes.

**Open statistical power, not impossibility.** A zero-success baseline does
not mathematically prevent the selected arm from certifying. A result changing
0 to1 would have been observable. However, all four OU selections are the
baseline x²/2, so those four pairs do not test a weight intervention. Only the
eight GBM/null pairs change functions. Both arms failed to resolve drift-x or
certify in the registered short ladder. There was no prospective power study
or positive-control demonstration of sensitivity to a weight change, so these
results cannot establish that the broader score-function approach is ineffective.
They establish no observed certification gain on these inputs. The paired
interval-width changes (including null median ratio0.693) are informative
continuous measurements; they are not rescued certification gains or formal
power estimates. The failed gain prediction remains withdrawn. No new tuning
or scientific run was performed to rehabilitate it.

## Evidence and serialization scope

**Empirical:** `CheckResult` owns the residual computed once by `check()`, a
band snapshot and original row labels. `residual_evidence()` builds the bounded
payload directly from that residual. The standalone `residual_measurement()`
API still accepts observed/predicted arrays for existing callers. Retained
observation/prediction attributes from the initial PR are removed; they had no
consumer and duplicated data. Malformed labels now raise before evaluation.
An unavailable measurement carries an explicit omission through certificates
and recover. On structural or parametric refusals a checked candidate can fit
all rows: its residual is diagnostic evidence, not an explanation of the
ambiguity. The payload explicitly identifies that distinction. No candidate is
evaluated just to populate a diagnostic.

**Numerical comparison contract:** raw C4 floats retain rtol1e-12, atol1e-15.
For artifacts deliberately rounded to p significant digits, each storage error
is at most .5*10^(1-p) relative; the optional profile permits the two rounding
errors in addition to the raw tolerance. Use

```
.venv/bin/python -m experiments.compare_artifacts first.json second.json --significant-digits 12
```

The profile is explicit; it does not silently loosen C4's raw-float contract or
any scientific gate. Nonfinite numbers, changed keys/types and material changes
still fail. The seed artifact now labels constructed inputs separately from
computed observations, and calls construction identities identities.

Reproduce the review witnesses without modifying the working tree's core:

```
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m experiments.run_pr14_review --baseline
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m experiments.run_pr14_review
OPENBLAS_NUM_THREADS=1 .venv/bin/pytest -q tests/test_pr14_review.py
```

Other touched artifacts have their existing runners. Tests run one file per
process. The original research results and registrations remain visible.

## Validation of this revision

**Empirical:** 233 tests passed, each file in its own process (the core file
was split into two processes). Counts: review24, retained evidence9, refusal7,
instrument21, instrument C4 25, public surface55, MCP17, passive5, Itô42,
boxsearch3, submission3, significance3, selection2, and core17 including both
noise guarantees. Four other slow core tests were not run; no full-suite claim
is made. Existing numerical/empty-array warnings remain. The instrument
refusal-message regression was corrected and its file rerun successfully.
F/E9 lint and whitespace checks passed.

The review artifact contains46 counter-input observations. Its pinned baseline
has13 exceptions; the corrected artifact's only exceptions are the three
intentional row-label ValueErrors. The original44 refusal witnesses still meet
their expected dispositions. Array/copy reduction is an implementation change;
no new end-to-end speedup measurement is claimed.
