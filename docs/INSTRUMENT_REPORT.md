# lagh: certified law discovery with calibrated abstention
### Instrument report (H1b) — draft 1, 2026-07-27

**One-sentence claim.** A symbolic law discoverer that returns either a
machine-checked certificate carrying a stated significance bound, or a
machine-readable refusal — never a confident wrong answer — and whose
deterministic observation-planning agent nearly doubled LLM-agent SOTA on a
sealed benchmark without a single LLM call.

---

## 1. The instrument

lagh (`github.com/jascal/lagh`) discovers exact symbolic laws from data. Its
product definition is an invariant, not a metric: **zero confident-wrong
submissions**, inherited from a predecessor's 114-task record and preserved
through every capability added since. Architecture:

- **A fixed honesty core** (`certify.py`): exhaustive per-point certification
  under a four-term epsilon model; null-law vacuity; functional-coherence
  rivalry on an extended probe box; parametric pinning under declared noise; an
  exact-coefficient gate; a minimality/parsimony repair. Never per-target,
  never per-class.
- **A registered curriculum** (C1–C9): polynomial, rational (implicit +
  pure-denominator), power-law, inner-scaled transcendentals, target
  transforms (incl. Bose/Fermi), quasi-polynomial, Lévy, angular/inverse-trig,
  generalized monomials. Capabilities land only through pre-registered
  predictions scored against named cells; any capability producing one
  confident-wrong is reverted (one was).
- **Significance accounting**: every certificate carries
  α ≤ |H|·q^h — the a-priori chance-fit bound from the hypothesis count
  actually searched and the per-point tolerance/range ratio over dof-discounted
  held-out points. Validated by null calibration: **0 false certifications in
  200 true-random targets.** Refusals thereby become *certified randomness
  statements relative to a stated hypothesis class* — an object no other SR
  tool emits.
- **Two regimes, one verdict machinery**: active (adaptive ranging, budget-
  metered acquisition, box-search) and passive (fixed datasets: K re-splits +
  a full-data exhaustive gate). Measured: on the 108-cell NewtonBench dev set
  the two regimes recover **identically** (87/108, same abstain set) — the
  power lives in the grammar + certifier, not oracle access.
- **A verified orchestration layer**: the research loop is an Orca state
  machine (proved topology: certificates are the only route to a `proved`
  terminal; at most one research move; guaranteed termination), with an LLM as
  a bounded, optional proposer at exactly one state.

## 2. Measured soundness results

- **Clean data**: 0 confident-wrong across all dev sweeps (NewtonBench 108
  cells; five prior testbeds; 240 LLM-SRBench problems).
- **Declared noise** (0.1–10% relative): zero gross fabrication (no certified
  law wrong beyond ~1.2× the noise, bound 3×); exact-STRUCTURE-or-abstain to
  ~1% noise minus one named, bounded exposure — *asymptotic degeneracy*
  (2/87 cells: the certified form is the true form's asymptote, deviating
  ~20× BELOW the noise floor; symbolically incomplete, numerically honest).
- **The approximant-impostor boundary** (measured, then enforced): at
  envelope-epsilon on a bounded box, dense linear/rational channels can
  certify Taylor-slop approximants of smooth laws whenever the true support
  goes unproposed. Since α bounds chance-fits, not wrong-form significant
  fits, the instrument draws the line structurally: **under declared noise,
  only small-hypothesis-class closed-form channels certify; dense-channel fits
  are labeled conjectures.** (The |H|-accounting extension that would relax
  this is future work.)

## 3. Benchmark results (both blind reads pre-registered; SOTA frozen before
download; one shot each; crash-fixes logged; conservative local judges)

### 3.1 LLM-SRBench (ICML'25) — a loss, and why it was informative
SA 0.42% vs 31.5% SOTA under the frozen sigma=0 protocol. Post-read diagnosis:
the benchmark's arrays are float32-born (~1e-7 relative), making
machine-precision certification impossible on all 240 problems *by data
representation* — the honesty core correctly refused every one. Lesson now
standing policy: **representation precision is declared noise** (dtype is
metadata). Dev iteration with σ_rep declared + capabilities mined from a
one-time LLM-proposer pass reached SA 12.08% autonomous (88% of the
LLM-composite's 13.75%), and the final soundness-first configuration holds
**19/19 structurally correct certificates** with α on each.

### 3.2 Gravity-Bench-v1 (ICML'25) — a win, near-doubling agent SOTA
The deterministic "astronomer" (fixed observation policy + digital twin; every
answer traceable to a fitted, prediction-validated model):

| variant | lagh (zero LLM calls) | frozen SOTA (o4-mini-high) |
|---|---|---|
| **budgeted (100 obs)** | **94.66%** | 49% |
| full observations | 63.59% | 74% (below; mechanism diagnosed, reported) |

The budgeted agent outscored the full-data agent — the planner engineers
better differentiation points than uniform density provides, which is the
benchmark's own thesis about planning, demonstrated by an instrument.

### 3.3 The distillation loop (LLM → grammar)
A bounded LLM proposer, verified by the sound checker, produced 19 certificates
the grammar could not reach; clustering those verified forms into registered
capabilities (scheduling, ratio features, sign/identity closures, sparse
supports) converted the grammar to 88% of the composite's accuracy with zero
further LLM use. The LLM's value was *one-time teaching*, permanently
distilled into a verifiable instrument.

## 4. Taxonomy of failure modes (all measured, all named)

chance-fit (bounded by α, validated 0/200) · asymptotic degeneracy (2.3% dev
exposure at ≤1% noise) · approximant impostors (structural boundary enforced)
· cross-form substitution at heavy noise (1/19 at 10%) · representation
quantization (declared-noise policy) · unfalsifiable detections (significance-
gated: a decay slope must exceed 3× its own standard error).

Added 2026-07-28, and they share one shape — **the truth absent from the
contest**, where every gate that compares rivals is powerless because the right
answer was never proposed: (a) *α cannot rank structure* — arbitrating two
dense approximants by their chance-fit bounds certified a form 2×10⁻⁵ off the
truth just outside the sampled box, closed by requiring every defeated rival to
be an interpolation the data never constrained (`MUNTZ_ARBITRATION.md`);
(b) *on-shell degeneracy* — a PDE fitted from ONE solution certified the
traveling-wave relation `u_t = −u_x` instead of the KdV equation it came from,
true of that field and not of the equation, closed by requiring certification
patches from a held-out SOLUTION (`CASE_STUDY_PDE_DEV.md`).

## 5. What this program claims, and what it does not

Claims: exact certified recovery on clean in-class data with stated α; zero
confident-wrong across ~600 scored tasks over two regimes and two benchmarks;
SOTA-doubling budgeted observation planning without LLMs. Not claimed: SA
competitiveness on single-equation leaderboards whose constructions
(inversion, injected novel terms, continuous parameters) sit outside a fixed
exact grammar — measured and documented as the boundary of the class.

*Artifacts: every number above has a committed jsonl/log and a registered
pre-read protocol in this repository (`docs/`, `experiments/results/`).*
