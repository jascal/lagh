# LLM-SRBench as a DEV benchmark — registered protocol + predictions

**Status: DEV, permanently** (the blind read was spent 2026-07-25,
`BLIND_READ_REPORT.md`). Everything here is dev iteration under `STRATEGY.md`
I7-at-benchmark-level: numbers are dev metrics, never a headline. Blind-read
artifacts (`blind_llmsrbench_*.jsonl/json`) stay frozen; dev outputs are
separate files.

## Measured data-quality facts (drive the protocol)

- Arrays are float32-born: gt self-residual median ≈ 5e-8.
- **Input quantization amplifies through gradients**: per-problem max gt
  self-residual p90 ≈ 1e-4, worst ≈ 3e-3. A single tight sigma cannot cover the
  tail without being wrong for the head.
- The conservative sympy judge failed to parse ~14 gt expressions (functional
  notation `A(t)`); fixed for dev (AppliedUndef → Symbol before mapping). Blind
  scores stay as scored.

## Dev protocol v1

- **σ_rep = 1e-4 declared globally** (the amplification envelope; a property of
  the benchmark's data format, not per-target tuning). RNOISE says the
  structural guarantee holds far past this (to ~1%), and the SA judge strips
  coefficients anyway, so structure is the claim.
- Two-track submission as before, plus the **k-form LLM proposer on abstain**
  (machine/llm.py, k=3, characterization context): each proposed form is
  `verify`d against the data at σ_rep — a verified proposal is a *certificate*
  (the checker is the same sound core), else the best proposal becomes the
  track-B conjecture in place of the log-log probe.
- 600 s/problem cap, 8-way process parallelism, 600-row discovery subsample
  with full-train gate at σ_rep.

## Registered predictions (before the run)

- **P1:** declaring σ_rep unlocks certificates on a nontrivial subset of
  LSR-Synth (forward skeletons intersect C1–C9); prediction: **≥ 15 of 129**
  synth problems certify, vs 0 at sigma=0.
- **P2:** LSR-Transform stays certificate-poor (**< 10 of 111**) — inversion
  orthogonality is a grammar fact, not an epsilon fact.
- **P3:** zero *structural* confident-wrong among certificates (the RNOISE
  structural regime; coefficient exactness is not claimed at σ_rep).
- **P4:** the LLM proposer adds SA points on Transform specifically (it can
  guess inverse forms the grammar lacks; the sound checker filters). Prediction:
  **≥ 5 problems** where a proposer form verifies or scores SA-correct.
- **P5:** overall SA beats the blind read's 0.42% by an order of magnitude;
  whether it approaches the 31.5%/20.2% SOTA numbers is the open question the
  sweep answers (no prediction).

## Amendment v1.1 (after 16 problems of v1, before any scoring)

v1's single 600s cap starved the proposer: sigma-widened discovery is slower
(more certifying candidates per tier), so 14/16 early problems timed out
BEFORE the LLM stage and emitted bare abstains. v1.1 = per-stage caps
(discovery 350s; proposer+verify 200s; log-log fallback always reachable).
The 16 v1 rows were discarded and the sweep restarted clean. Positive early
signal kept for the record: 2 grammar certificates in those 16 (P1 direction).

## Results (v1.1 sweep 2026-07-25; judge v3 = as_independent constant strip +
locals-dict gt parsing; scores in `dev_llmsrbench_v1_scores_v2.json`)

| | SA | Acc₀.₁(ID) | certified | cert struct-wrong | notes |
|---|---|---|---|---|---|
| **ALL** | **12.5%** | 16.25% | 36 | 6 | blind read was 0.42% |
| **LSR-Transform** | **27.0%** | 27.9% | **30 (all SA-correct)** | **0** | frozen SOTA 31.53% |
| LSR-Synth (4 domains) | 0% | ~6% | 6 | 6 | novel-term wall |

Certified channels: grammar 16 (incl. the 1 judge-error), **llm-verified 19** —
grok proposes, the sound checker certifies at σ_rep. Conjecture track: 108
submissions, 0 SA-correct (the log-log floor and raw LLM guesses earn numeric
Acc points but no structural hits).

**Prediction scoring:**
- **P1 FAILED** (≥15 synth certificates; got 6, all structurally wrong): synth's
  novel terms are the real wall — not representation, not epsilon.
- **P2 FAILED in the good direction** (<10 Transform certs; got 30): the
  llm-verified channel was not in the prediction's model of the world.
- **P3 FAILED as stated, with structure:** Transform certificates are **30/30
  structurally correct — zero fabrication where a full structure was reachable.**
  All 6 synth certificates are structurally wrong the same way: a certified
  skeleton that omits the novel term, which stays inside the σ_rep envelope even
  on OOD (all 6 pass Acc₀.₁ OOD). This is the RNOISE **asymptotic-degeneracy
  class at benchmark scale** — the certificate is numerically honest and
  symbolically incomplete. Naming + gating this (e.g. residual-structure probe
  before certifying at σ>0) is now a registered instrument direction.
- **P4 MET, exceeded** (≥5 proposer wins on Transform; got 19 verified, all
  SA-correct).
- **P5 MET** (12.5% = 30× the blind 0.42%).

**Read on the SOTA gap:** Transform 27.0% vs LLM-SR's 31.53% — within 4.5 points
of the best published system, with a claim it cannot make (30 machine-checked
certificates, zero structurally wrong). The gap is entirely conjecture-side;
synth is a grammar/term problem, not an epsilon problem. Next (no-LLM per user
directive): mine the 19 verified forms into grammar capabilities; census the
synth novel terms; build the residual-structure gate.
