# Reproducing the lagh Gravity-Bench result

Result being reproduced (one-shot blind read, 2026-07-27 — registration and
report in [`docs/BLIND_READ_REGISTRATION_GRAVITYBENCH.md`](../../docs/BLIND_READ_REGISTRATION_GRAVITYBENCH.md)
and [`docs/BLIND_READ_REPORT_GRAVITYBENCH.md`](../../docs/BLIND_READ_REPORT_GRAVITYBENCH.md)):

| variant | score |
|---|---|
| budgeted (100 observations, planned) | **94.66%** (195/206) |
| full observations | 63.59% (131/206) — below the 74% baseline; mechanism in the report |
| full observations + dev epoch fix (labeled DEV, post-read) | 91.75% |

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
