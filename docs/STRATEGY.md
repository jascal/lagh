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

## FULL SWEEP + READINESS ASSESSMENT (2026-07-21, `experiments/results/newtonbench_all.jsonl`)

After CAP-B (trig), the C2 power-denominator bug fix, and CAP-D (`sqrt`-transform):
**64/108 recovered (was 55), 64/88 algebraic = 73% of the exact ceiling, confident-wrong 0/108.**
Per-difficulty **module** coverage (any version): **easy 11/12, medium 11/12 (was 8/12), hard 6/12.**

| gate | status | evidence |
|---|---|---|
| **R-cap** (≥10/12 on ≥2 difficulties) | **✅ MET** | easy 11/12 **and** medium 11/12 |
| **R-zero** (zero confident-wrong, full sweep) | **✅ MET** | 0/108 |
| **R-gap** (each gap closed or bounded w/ cause) | **≈ needs a bounding pass** | remaining abstains classify cleanly (see below); doc them one-line each |
| **R-noise** (explicit noisy-regime decision) | **◐ SETTLED** (`RNOISE_STUDY.md`, corrected 2026-07-22) | with `certify.pinned()` + STRUCTURAL scoring: exact-structure-or-abstain holds to ~1% noise (**0/65 structural CW** at 60 & 40 dB), 3/65 cross-form failures only at 10 dB. The original '36 CW' was a coefficient-representation scoring artifact. Blind read registrable for **SNR ≥ ~40 dB with the structural guarantee**; clean-only only if the exact coefficient must be certified |

**Remaining abstains, all accounted for (R-gap material):** 11 irrational-`^e` cells (out of class BY
DESIGN — the wedge, correct-abstain); 9 inverse-trig (all m4_snell — needs an `arcsin/atan` tier,
+9 ceiling); the rest are CAP-A fractional-grid (hooke, frac-underdamped), CAP-C (gravity/coulomb
high-degree), CAP-E (BE rational-of-exp), CAP-F (sound transcendental×power) — each a registered,
scoped capability in `NEWTONBENCH_GAP_PLAN.md`. The old "underdamped is the R-cap risk" note is
**superseded**: underdamped was a C2 reach bug (now 4/9), not the blocker.

**Bottom line:** two of four gates (R-cap, R-zero) are now MET on dev. The blind read remains sealed
— it is gated on **R-noise (an explicit decision)** and an R-gap bounding pass, and the blind read is
a deliberate one-shot the user authorizes, never an automatic consequence of hitting coverage.

## RE-SWEEP 2026-07-23 — post-capability baseline (supersedes the table above)

After the exact-coefficient gate, the escalating-snap reach fix, CAP-C, CAP-G
(angular/inverse-trig), CAP-E (bose/fermi), CAP-F (generalized monomial), CAP-A2:

**87/108 recovered, confident-wrong 0/108. Module coverage: easy 12/12,
medium 12/12, hard 10/12.** Ceiling accounting: 97 grammar-representable − 3
oracle-precision-limited = 94 measurable → **87/94 = 93% of the achievable ceiling.**

| gate | status | evidence |
|---|---|---|
| **R-cap** | **✅ MET (exceeded)** | easy 12/12, medium 12/12, hard 10/12 |
| **R-zero** | **✅ MET** | 0/108 across the full re-sweep |
| **R-gap** | **✅ MET** | `RGAP_CLOSURE.md` rewritten 2026-07-23: all 21 abstains bounded (11 wedge, 3 oracle-precision, 6 scoped-open, 1 registered exclusion), zero unexplained |
| **R-noise** | ◐ unchanged | `RNOISE_STUDY.md` decision stands (structural guarantee to ~40 dB); must be re-confirmed on the post-capability instrument before registration |

**New since the last table:** the PASSIVE regime (fixed datasets, no oracle — what
the blind candidates actually are) is built (`DIRECTION_PASSIVE.md`, predictions
pre-registered, sweep in flight), and the submission policy for accuracy-scored
benchmarks is pre-registered (`DIRECTION_OUTPUT_POLICY.md`: certified track +
labeled-conjecture track). The blind read remains sealed and user-authorized; the
remaining pre-read work is the passive-sweep verdict + an R-noise re-confirmation.

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

---

## Reconnaissance: a third category (added 2026-07-29, after PDEBench)

The taxonomy above has two members: a benchmark is **DEV** (used freely, its
numbers are dev metrics, never a headline) or **BLIND** (sealed, read once,
scored). PDEBench fitted neither, and the cost of not noticing that in advance
was paid in a single long session.

A **RECONNAISSANCE** pass is first contact with a data source whose character is
not yet known. Its purpose is to characterize the SOURCE, not to score the
instrument. What makes it a category rather than an excuse is that it carries its
own registration, and the registration is about the data:

> Before opening the files, state what the source is expected to BE — its error
> scale and kind, whether its stated laws are expected to hold, what would count
> as the pass failing — and state what follows from each answer.

`PDEBENCH_READINESS.md` registered three declarations (float32 σ_rep, an
undeclarable solver error, geometry checks) and all three earned their place. But
they were declarations about **mechanics** — how to load, what to declare — not
predictions about the target. Nothing was registered that the first run could
FALSIFY. So when the truth check fired at 1185× on the opening run, it arrived as
a surprise to be interpreted rather than as a prediction being cleanly refuted,
and everything downstream was mid-course correction: the DEV reclassification,
the pipeline-decode reframing, inventing "report the required declaration" for
CFD. Each was right; each cost more than it should have.

The single line that was missing:

> **The shipped fields' deviation from their own stated laws will be at or below
> storage precision.** If it is not, the target is model output whose generating
> error dominates, the exercise is pipeline decode rather than discovery, and no
> score is available from it.

That is falsifiable in the first hour, and it was false by four orders of
magnitude. The whole framing of the pass would have been its *starting* position.

### The gating rule this forces

A reconnaissance pass finds defects — PDEBench found eight — and the temptation
is to fix them in place and keep going. But **a pass with nothing registered
cannot tell you whether a change made during it is an improvement**: there is no
baseline anyone committed to in advance, so "the new number looks better" is not
evidence. Therefore:

> **An engine change discovered while mining an uncharacterized target is SCOPED
> to that target's path until it has been validated against the campaigns that
> predate it.** Prove identity by construction where possible; measure it where
> not.

Measured on 2026-07-29, and the rule is not decorative. The coherence early exit
(a 119× speedup on the case that motivated it) was first wired into the main
`discover` path, where an audit found it could turn an arbitrated winner into a
structural abstain on an open-ended library — `arbitrate_significance` scores a
class by its min-complexity representative while the exit tested the member that
opened the class, and complexity is not dof. Conservative, never a false
certificate, but a **reach regression**, which is not a Pareto improvement. Gated
behind `linear_basis`, the change keeps its speedup and every campaign that
predates it is bit-identical by construction rather than by testing. The one
change that could NOT be gated (reading a linear candidate's coefficients instead
of differentiating it symbolically) was validated the other way: the C1 ladder
re-run bit-identical across all 15 rungs.
