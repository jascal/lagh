# Benchmark strategy — dev/test separation at the benchmark level

**Decided 2026-07-21.** The program's `I7` ("test read once") is now applied one level up, to
whole benchmarks, because the only viable deliverable is beating a high-profile benchmark and a
win is worthless if the benchmark was used during development.

## The decision

1. **NewtonBench is a DEVELOPMENT set, permanently. No win claim will ever be made from it.**
   Justification: wyly's N1/N1′ already read its easy-vanilla cells, and `lagh`'s C4/C5/C6/v2
   capabilities were built partly informed by which families failed there. It is contaminated as a
   test by construction. It is therefore used *freely* from here on — all difficulties, all cells,
   iterate as much as useful — to drive and measure capability. Its numbers are **dev metrics**,
   never a headline.

2. **A brand-new, never-touched benchmark is RESERVED as the single blind test, read exactly once,
   only when development is declared complete** (readiness bar below). Until then it is not
   cloned, not inspected, not sampled — inspection is contamination.

3. **Candidate blind benchmarks are named by reputation/metadata only, and not opened.** Selection
   among them happens at the moment of the blind read, from public leaderboard/SOTA facts frozen
   into the run registration first (as `GOAL.md` R4 required for NewtonBench).

   - **LLM-SRBench** (ICML 2025) — symbolic regression designed explicitly to defeat memorization,
     multi-domain, published LLM-SR baselines. *Primary candidate:* recent, high-profile,
     anti-contamination by design, so a genuine-discovery win is meaningful. Untouched.
   - **SRSD-Feynman** (Symbolic Regression for Scientific Discovery) — rigorous exact-recovery
     protocol over the Feynman equations. Untouched.
   - **SRBench** black-box track — older, more saturated; fallback only.

   None of the above has been read by this program. This document is the commitment that they stay
   sealed.

## Why lagh can win a blind SR benchmark (the thesis of the attempt)

lagh's demonstrated edge is **exact certified recovery when the law is in C1–C6, with zero
confident-wrong** (0/174 across five beds). On an SR benchmark scored by *symbolic accuracy* /
*exact recovery rate*, an incumbent that fits within tolerance but not exactly loses a point that
lagh's exact recovery wins — and lagh never emits a confident-wrong to lose points the other way.
The edge is strongest on clean / low-noise cells whose true law is elementary; it is honestly
absent where the law is outside the curriculum (lagh abstains there, scoring 0 like a wrong answer
but never *below* by fabricating).

## Readiness bar (pre-registered, so development terminates)

The blind read is spent only when **all** hold on NewtonBench-as-dev (vanilla, all difficulties):

- **R-cap:** ≥ 10/12 modules recovered on ≥ 2 of 3 difficulties (broad curriculum coverage).
- **R-gap:** every known curriculum gap is either closed or explicitly bounded with a one-line
  cause (today's open gap: joint-two-parameter products `e^{−bt}·cos(ωt)`, underdamped harmonic).
- **R-zero:** the zero-confident-wrong invariant still holds across the full dev sweep.
- **R-noise:** a decision on the noisy regime — either lagh clears a stated noisy-dev bar, or the
  blind read is declared *clean/low-noise only* and the benchmark cells are restricted accordingly
  in the registration (MDBench proved noise is where this method is weakest; the blind claim must
  not be silently exposed to it).

When R-cap ∧ R-gap ∧ R-zero ∧ R-noise, freeze the blind benchmark's published SOTA into a run
registration, then read it **once**. Not before.

## What NewtonBench-dev is telling us so far (2026-07-21, partial)

m0–m3 (gravity, coulomb, magnetic, fourier) all versions recovered, zero wrong — the
power-law/rational/linear families are solid. Full sweep and per-difficulty coverage pending;
the underdamped-harmonic gap is the known risk to R-cap.

## Amendment 2026-07-21 — no August competition target; deadline framing dropped

We did **not** register for AISB/NLPCC Task 9 (registration closed 2026-05-25), so its 2026-08-01
submission is unavailable. A fresh web survey confirms the earlier finding: **the LLM-science /
symbolic-discovery space is paper-benchmarks, not live competitions** — no open submission server
with a near-term deadline exists that fits. The "August deadline" therefore has **no valid
target**, and the deadline framing is dropped. This *confirms* the disciplined path: no artificial
deadline; readiness bar → reserved blind benchmark → self-scored blind read, with `I7` standing in
for the extinct held-out server.

Sealed candidate leads from the survey (named by metadata only, UNOPENED):
- **LLM-SRBench** (ICML 2025, 239 problems, 4 domains, anti-memorization) — primary reserved blind.
- **Gravity-Bench-v1** — gravitational-physics discovery *for agents*; queryable, agent-shaped →
  a direct fit for the tool-shape (LLM + lagh active recovery). Strong secondary.
- **SURFACEBENCH** — symbolic surface discovery, geometry-aware. Symbolic, in-class-adjacent.
- **REFUTE** (2026.06) — scores missing-evidence **refusal**, **calibration**, planted-flaw
  detection. The one family where *abstention is the scored metric*: a native fit for the
  zero-wrong + significance-α story. Flags a whole benchmark class (epistemic-calibration) where
  the federation wins by design rather than by bolt-on.

Decision: keep developing lagh against NewtonBench-dev to the readiness bar; build the tool-shape
(MCP) so lagh/i-orca become LLM-orchestrated; then a single blind read of a reserved benchmark on
our schedule. The externally-adjudicated live win remains desirable but is not currently available;
watch for a new competition (esp. calibration-scored) rather than force a dead one.
