# Reproducing the lagh Gravity-Bench result

Result being reproduced (one-shot blind read, 2026-07-27 — registration and
report in [`docs/BLIND_READ_REGISTRATION_GRAVITYBENCH.md`](../../docs/BLIND_READ_REGISTRATION_GRAVITYBENCH.md)
and [`docs/BLIND_READ_REPORT_GRAVITYBENCH.md`](../../docs/BLIND_READ_REPORT_GRAVITYBENCH.md)):

| variant | score |
|---|---|
| budgeted (100 observations, planned) | **94.66%** (195/206) |
| full observations | 63.59% (131/206) — below the 74% baseline; mechanism in the report |
| full observations + dev epoch fix (labeled DEV, post-read) | 91.75% |

> **The code in this repo no longer reproduces those numbers exactly, and that
> is deliberate.** A mass-estimate bug was found and fixed on 2026-07-29, after
> the read (details below). **To reproduce the read as it was run, check out
> `cf54706`** — the commit that recorded it. Current `master` scores are DEV
> numbers, reported separately and never as the read.

## What this pipeline is (and is not)

A **fully deterministic, non-LLM baseline**: a fixed observation-planning
policy plus a fitted "digital twin" that answers all 50 task types. It does
NOT run inside the GravityBench agent harness — it drives the same
observation data (the HF dataset's simulation tables, interpolated exactly as
the harness's `Observe` tool does) and is scored mechanically against the
dataset's own per-instance thresholds (`budget_obs_threshold_percent` /
`full_obs_threshold_percent`) versus `true_answer`. It is therefore a
different baseline *class* from the paper's LLM agents — closer in kind to
the expert-solution reference. If our threshold interpretation differs from
the official scoring, we want to know.

## Steps

```bash
git clone https://github.com/jascal/lagh && cd lagh
python3.12 -m venv .venv
.venv/bin/pip install numpy scipy sympy pandas datasets

# the one-shot read (both variants, ~1-2 h; incremental journal, resumable)
.venv/bin/python experiments/gravitybench/run_blind_read.py

# expected: experiments/results/gravitybench_read_scores.json ->
#   {"budget": {"n": 206, "correct": 195, "pct": 94.66},
#    "full":   {"n": 206, "correct": 131, "pct": 63.59}}

# the post-read DEV fix for the full variant (native-cadence epoch triplets)
.venv/bin/python experiments/gravitybench/dev_full_retest.py
# expected: 189/206 = 91.75%
```

Determinism: all sampling is seeded; observation times are chosen by a fixed
policy; no network calls except the HF dataset download. The committed
journals (`experiments/results/gravitybench_read.jsonl`, `..._scores.json`)
are the exact artifacts behind the report.

Development provenance: the agent was built and validated ONLY on self-built
synthetic orbits (`integrator.py`, `battery.py` — 8 hardening rounds,
committed) before the benchmark data was first opened; the registration doc
predates the download in git history.

## Post-read fix, 2026-07-29: the snapped exponent kept a stale intercept

Found by a test that had been committed RED — `test_twin_end_to_end` asserts the
twin reproduces the observations it was fitted from to 5%, and it had been
returning 11.2% since `4f096aa` ("astronomer hardened"), which is the commit
that introduced explicit triplet metadata. It passed at `b9a211a` and failed at
`4f096aa`; nobody re-ran it in between.

`system_id` fits `log|a_rel| = log(GM) + p log r`, then SNAPS `p` to the
Newtonian −2 when it is within tolerance. It was not re-fitting the intercept
after moving the slope, so `GM` continued to describe the force law at the free
`p` while everything downstream used `p_used`. The bias is `r_mid^(p − p_used)`:
on the dev orbit `p_raw = −1.999828`, `r_mid = 1.5e11`, bias **0.444%**, and the
total mass came back **0.445% low**.

That size of error is invisible to every direct answer — the tolerances here are
3–10% — and fatal to the twin, which is *integrated forward*. A 0.44% mass error
is a 0.22% period error, which over a four-period window is a 5.6% phase drift
and 11% at periastron on an `e = 0.3` orbit. Re-fitting the intercept at the
snapped exponent takes the mass error to **0.0029%** and the twin validation
from **0.112 to 0.00006**.

Effect on the synthetic battery — a Pareto improvement, which is what the
gating rule asks of a change like this:

| | before | after |
|---|---|---|
| cases improved / unchanged / regressed | — | **40 / 3 / 0** |
| twin-validation median | 0.1120 | **0.0023** |
| battery verdict | 42 / 43 (drag 17.3%) | **42 / 43 (drag 17.3%, unchanged)** |

`run_blind_read.py` grew `--out` / `--scores` at the same time, so a post-fix
re-run cannot append to the recorded read: the runner skips rows already in its
journal, so pointing it at the committed file would have silently restated a
one-shot result instead of producing a new one.

### What the fixed code scores — DEV, not a read

Both variants re-run against the same dataset
(`gravitybench_read_postfix{,_scores}.json`):

| variant | read as executed (`cf54706`) | post-fix DEV |
|---|---|---|
| budgeted (100 obs) | 94.66% (195/206) | **94.66% (195/206)** — bit-identical, no instance flipped either way |
| full observations | 63.59% (131/206) | **86.89% (179/206)** — 50 newly correct, 2 newly wrong |

The budgeted variant is untouched, which fits the mechanism: its answers are
read off a twin fitted to a short, dense, well-conditioned observation plan,
and none of them sat within 0.44% of a threshold. The full variant is where the
bias was doing its damage — median twin validation across the read drops from
**0.538 to 0.0074**.

**The two regressions are worth more than the 48 net gains**, because they are
a case of being right for the wrong reason. Both are extremes
(`max_velocity_star1`, `max_angular_velocity_star1`) on one scenario,
`9p6_M_3p1_M_Proper_Motion2`, whose twin validation improved from **4.09 to
0.62** — that is, from 409% max reconstruction error to 62%. It was never a
usable twin and still isn't; the old answers happened to land inside the 5%
threshold (2.2% and 0.4% error) and the new ones land outside (6.3% and 8.6%).
Nothing about that scenario got worse except the luck. **The instrument already
computes the number that would have flagged it** — `twin_validation` is
recorded per instance and never gates the answer. Gating on it would convert
both of these from a wrong answer to an abstain, and is the registered
follow-up. (Wrong ANSWERS: this pipeline is Track B of
`docs/DIRECTION_OUTPUT_POLICY.md` and emits no certificates, so nothing here is
a confident-wrong in this program's sense of the term — a false claim carrying
an alpha. It never calls the certified engine at all.)

### The gate, built 2026-07-29

`Twin.gated_answer` returns `(answer, validation, refusal)` with `answer = None`
exactly when `refusal` is set, so a caller cannot read a number without being
told the model behind it does not fit. `TWIN_VALIDATION_MAX = 0.05`.

**The bar comes from the synthetic battery alone and was fixed before the
benchmark effect was measured.** Across the 43 dev cases, the 42 whose worst
answer error is ≤ 15% all validate at ≤ 0.0231; the one that fails (drag, 17.3%)
validates at 1.7562. A 76× gap with nothing in it. 0.05 sits 2.2× above the
worst passing dev case, 35× below the failing one, and is the bar
`test_twin_end_to_end` already used for the same quantity. Under it the dev
battery abstains on 1 of 43 cases and **no wrong answer survives**.

It gates **every** task. The narrower reading — gate only the answers read off
the integrated trajectory, since masses and the exponent come from fits — is
defensible and was rejected: this pipeline's claim is "every answer traceable to
a fitted, *prediction-validated* model", and a twin that cannot reproduce its own
fitting data has validated nothing it emits.

**It costs score, and that is the trade being made.** Abstain scores the same as
wrong here, so the gate can only lose points:

| variant | DEV score | | wrong answers submitted | |
|---|---|---|---|---|
| | ungated | gated | ungated | gated |
| budgeted | 94.66% | **91.26%** | 9 | **5** |
| full | 86.89% | **81.07%** | 26 | **9** |

13 of 206 budgeted instances abstain (7 of them would have been right) and 30 of
206 full instances (12 would have been right). Confident-wrongs fall 44% and
65%. The gate does **not** eliminate them: 5 and 9 survive with good validation,
which is the honest limit — those errors are invisible to this diagnostic, and
a twin reproducing its observations is necessary for a right answer, not
sufficient.

The other post-read DEV path, `dev_full_retest.py` (native-cadence epoch
triplets), is **byte-identical before and after** — still 189/206 = 91.75%, so
the expected output quoted above is unchanged. It builds its epoch state a
different way and never consumed the biased mass. Note it still beats the
mass-fix-only full number (91.75% vs 86.89%): the two fixes address different
things and the epoch one is doing more work on this variant.
