# BLIND READ REGISTRATION — LLM-SRBench

**Registered 2026-07-24, BEFORE the dataset is downloaded or any problem is
inspected.** User authorization for the one-shot read: given 2026-07-24
("Proceed with 1 then 2"). Everything below is frozen; nothing may be changed
after the data is opened. Per `STRATEGY.md` this is the single reserved blind
test; its numbers are the headline, read exactly once.

## 1. Frozen public SOTA (source: arXiv 2504.10415 / ICML 2025 oral, fetched
2026-07-24 — metadata only, dataset untouched)

- Benchmark: **LLM-SRBench**, 239 problems = **LSR-Transform 111** (Feynman
  equations re-represented) + **LSR-Synth 128** (chem 36, bio 24, phys 43,
  matsci 25). Data: HuggingFace `nnheui/llm-srbench`.
- Official metric: **symbolic accuracy** judged by GPT-4o (mathematical
  equivalence after removing parameters/constants; 94.6% human agreement),
  plus numeric **Acc_τ (τ=0.1)** and **NMSE** on ID test sets (OOD test sets
  exist for LSR-Synth only).
- **Frozen SOTA:** best published system = LLM-SR (GPT-4o-mini backbone):
  **31.53% SA on LSR-Transform, 20.24% SA on LSR-Synth** (LaSR: 6.31 / 28.12;
  best direct prompting: 7.21 / 0.0).

## 2. Frozen run protocol

Per problem (all 239; no problem skipped or added):

1. Load `train` samples only for discovery. **`gt_equation` is never read by the
   discovery path** — it is touched exactly once, by the scoring function, after
   all submissions are produced.
2. **Submission = `lagh.submit.submission(X_train, y_train, sigma=0)`** (the
   pre-registered two-track policy, `DIRECTION_OUTPUT_POLICY.md`):
   - Track A (certified): `discover_passive` certificate — K=3 re-splits +
     full-data exhaustive gate, `sigma=0` (the benchmark declares no noise
     model; if the data is in fact noisy, certification honestly fails and the
     problem falls to track B — no certificate is ever claimed on data that
     does not support exactness).
   - Track B (labeled conjecture): the `fit` scout / best non-certifying
     candidate, tag `empirical`.
   - Datasets with more than 600 train rows are deterministically subsampled
     (seed 0) to 600 for the discovery call; a track-A certificate must then
     ALSO pass the exhaustive check on the full train set or it is demoted to
     track B. Bounded-runtime rule, registered here.
3. Scoring, computed once, after all 239 submissions exist:
   - **Primary: symbolic accuracy** under a CONSERVATIVE local judge: sympy
     structural equivalence after stripping numeric constants (the official
     criterion, implemented symbolically: `simplify(pred_struct - gt_struct) == 0`
     or equivalent-set-of-terms match). This can only UNDER-count relative to
     the official GPT-4o judge (it recognizes fewer equivalences), so the
     reported SA is a **lower bound** on the official-protocol SA. No LLM judge
     is used (deterministic, reproducible, no leakage).
   - Secondary: **Acc_0.1** and **NMSE** on the provided ID test split, and OOD
     where present.
4. **Pre-registered accounting** (`DIRECTION_OUTPUT_POLICY.md`): SA overall and
   per category; track-A submitted / track-A SA-correct; track-B submitted /
   track-B SA-correct; abstain-with-nothing count; **track-A structural-wrong
   count (the invariant: must be 0)**.

## 3. Noise clause (R-noise gate)

`sigma=0` is declared for every problem (the benchmark publishes no noise
model). Consequence: the certified track claims exactness only where the train
data supports it at machine precision; noisy problems produce track-B
conjectures, never certificates. The passive R-noise re-confirmation
(`experiments/run_rnoise_passive.py`, running at registration time) informs the
*interpretation* of abstains but does not change this protocol.

## 4. Honesty commitments

- One read. No re-runs after seeing scores; no grammar changes, no threshold
  changes, no re-splits. Bug fixes required to *complete* the run (crashes)
  are allowed but logged in the report; scoring-affecting changes are not.
- The report states both tracks, the certified partition, and the comparison
  against the frozen SOTA above — win, lose, or mixed.
- NewtonBench remains dev; this read's numbers are the only headline claims.
