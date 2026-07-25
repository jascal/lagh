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

## Results

*(filled after the sweep; predictions frozen above)*
