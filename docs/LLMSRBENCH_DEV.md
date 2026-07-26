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
| **ALL** | **13.75%** | 16.25% | 36 | **3** | blind read was 0.42% |
| **LSR-Transform** | **27.0%** | 27.9% | **30 (all SA-correct)** | **0** | frozen SOTA 31.53% |
| LSR-Synth (4 domains) | ~2.3% | ~6% | 6 | 3 | novel-term wall, softened |

*(Judge v4: three successive judge defects fixed against the benchmark's own
data quirks — float exponents `x**0.333333333333333` rationalized; `A(t)`
functional notation collapsed; a literal `_f`-suffix ARTIFACT IN THE BENCHMARK'S
OWN gt STRINGS stripped. Each fix only ever upgraded our SA — the conservative
lower-bound property is preserved. 3 certificates remain structurally wrong,
all grammar-channel synth: CRK13, BPG3, PO18 — the true asymptotic-degeneracy
residue, smaller than first reported.)*

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

## Registered capabilities from proposal mining (2026-07-25, before code)

The 19 grok-verified certificates, clustered — each cluster becomes a registered
change with predicted cells, scored by re-probing those problems GRAMMAR-ONLY:

- **CAP-S — cost-aware cheap pre-pass (instrument scheduling, no grammar change).**
  12/19 verified forms are plain C3 power-law monomials the tier loop never
  reached: at dim≥4 the time budget dies inside C2's implicit enumeration before
  tier 3 runs. Pre-pass = {C3, C9, C8, C3-under-C5-transforms} candidates
  (≈50 cheap log-fits), full certification + gate + coherence + pinned; unique
  certifying class → return; anything else falls through to the normal loop.
  **Predicts ≥10 of the 12 monomial cells recover lagh-alone; zero regressions
  on the NewtonBench-dev 108 (re-run to verify).**
- **CAP-R — ratio features `x_i/x_j`** (dim≤5, complexity 3): unlocks
  affine-denominator rationals via the C2 implicit pass. **Predicts 4 cells**
  (`x_0/(x_1(x_2+1))` ×2, `x_0/(x_1x_2)−1`, `1+x_1x_3/(x_0x_2)`).
- **CAP-B bound lift 2→3 with an angle-plausibility guard** (trig features only
  for columns with range ⊂ (0, 2π)): **predicts the 3 sec-form cells.**

## Synth census (2026-07-25, judge-v4 parser) + registered capabilities round 2

Ground-truth term families (129 synth problems; 32 chem/bio parse-fails remain
— deeper `A(t)`-notation variants, census-only issue):

- **phys_osc (44):** sparse LINEAR combinations over a small nonlinear library —
  `sin(x_j)`, `x^3`, `x`, `exp(-|x|)`, `log(|x|+1)`, `|x|^{1/3}`, `x_i(1-x_j^2)`
  products, with per-problem float constants. This is lagh's linear channel with
  missing library terms.
- **matsci (25):** Arrhenius-family `x_0^p * exp(-E/x_1)` (+ shifted `(x_1-c)`
  factors) with CONTINUOUS fitted exponents — outside the exact-rational class
  by construction; winnable only as labeled conjectures (the official judge
  strips parameters).
- **bio/chem:** additive mixes of logistic rationals `x^2/(cx+1)`, `x^{1/3}`,
  `sqrt(x)`, exp terms — partially in-grammar today (2 chem certificates).

Registered, predictions before code:

- **CAP-P — damped/saturating library features** (dim≤3): `exp(-|x_j|)`,
  `log(|x_j|+1)`, `|x_j|^{1/3}`, and `x_i ×` each. **Predicts ≥15 of 44
  phys_osc problems certify** (gt = sparse linear combo once the features
  exist) **with zero structural wrongs.**
- **CAP-T — continuous-parameter conjecture mode (track B, no LLM):** when
  certification abstains, emit the UNSNAPPED c9/c3 fit (float exponents kept)
  as the labeled conjecture; add `1/x_j` to c9's factor pool (Arrhenius).
  **Predicts first nonzero matsci SA** (judge treats float exponents as
  parameters, official-protocol style).

## Grammar-only sweep v2 (2026-07-25, ZERO LLM calls, judge v4)

| | SA | cert (SA-ok / wrong) | vs v1.1 (with grok) |
|---|---|---|---|
| **ALL** | **12.08%** | 26 (24 / **2**) | 13.75%, 36 certs |
| **LSR-Transform** | **23.42%** | **21 (21 / 0)** | 27.0%, 30 certs |
| chem / bio | 5.6% / 4.2% | 5 (3 / 2) | first synth SA points |
| matsci / phys_osc | 0 / 0 | 0 | unchanged |

**The round's verdict: the instrument now reaches 88% of the LLM-composite's SA
autonomously** (12.08 vs 13.75), with zero fabrication on Transform (21/21).
The grok proposer's residual value is +1.7 SA points. Mining scorecard:
- **CAP-S/R/B/N/G2 + gate fix: PASS** — grammar-channel Transform certificates
  11 → 21, all structurally correct, most in <1 s via the pre-pass.
- **CAP-T: PARTIAL** — 5 Transform conjecture SA hits (first-ever conjecture-track
  SA); the matsci prediction FAILED (0 SA — continuous exponents fit numerically,
  acc01 8%, but the judge's structure signature keeps distinct float exponents).
- **CAP-P: FAILED as registered** (0/44 phys_osc certificates vs ≥15 predicted).
  Post-hoc diagnosis identifies the likely blocker, same family as the _FIT_TOL
  bug: the linear channel's `PREFILTER_REL = 1e-6` validation gate rejects
  every candidate whose residual sits at the declared σ_rep = 1e-4 — correct
  multi-term laws are filtered before certification on quantized data. A
  **sigma-aware prefilter** is the registered next lever (engine + c2), with
  CAP-P re-scored after it.

**Read on the SOTA gap:** Transform 27.0% vs LLM-SR's 31.53% — within 4.5 points
of the best published system, with a claim it cannot make (30 machine-checked
certificates, zero structurally wrong). The gap is entirely conjecture-side;
synth is a grammar/term problem, not an epsilon problem. Next (no-LLM per user
directive): mine the 19 verified forms into grammar capabilities; census the
synth novel terms; build the residual-structure gate.
