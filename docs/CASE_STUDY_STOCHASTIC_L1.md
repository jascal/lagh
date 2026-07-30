# Level 1, first increment: the diffusion is a claim, and it wants the opposite noise

**Run 2026-07-29,** continuing from [`CASE_STUDY_STOCHASTIC_L0.md`](./CASE_STUDY_STOCHASTIC_L0.md).
Machinery: `lagh/ito.py` (`build_rows(diff_names=...)`), `lagh/certify.py`
(`admissible_functional`). DEV, by construction, as the whole suite is.

Level 1's registered targets are *symbolic drift; symbolic diffusion or a reasoned
abstention on it; any almost-sure invariants* on four systems. This document covers
the first increment — making the diffusion a **claim** rather than a measurement,
and what measuring it revealed. The two 2-D systems, the invariants target and the
remaining nulls are not in it, and §*What is not here* says so precisely.

## The mechanism: b² in the design matrix

At Level 0 the Itô correction was a MEASURED quantity moved to the left-hand side.
Given a library {h_j} for b², it becomes ordinary `dt` columns instead:

    -∫φ' f(X) dt  =  Σ_k c_k ∫φ f'(X) g_k(X) dt
                   + Σ_j d_j · ½∫φ f''(X) h_j(X) dt  +  M

so drift and diffusion are identified **jointly, from the same rows**, and the
columns are ordinary weak-form terms — no new `weakform` capability was needed
beyond what step 3 already landed. Columns are named `drift:<g>` / `diffusion:<h>`,
which is the frozen checker's component vocabulary.

**This also removes the one unsafe consumer of realized quadratic variation.** The
measured correction was where observation-noise contamination became a systematic
offset on the *target* — the Level 0 confident-wrong. With b² in the design matrix,
realized QV only sets the BAND, where contamination is conservative. That is a
soundness improvement obtained by changing what is claimed, not by adding a guard.

## The finding: drift and diffusion want opposite noise intensities

GBM, `dS = μS dt + b S dW`, μ = 0.8, T = 6, L = 1, 80 rows, δ = 0.05. Relative
width = the joint bound's width divided by the true coefficient.

| b | `drift:x` rel. width | `diffusion:x²` rel. width | resolved |
|---|---|---|---|
| 0.02 | **1.5** | 6388 | drift |
| 0.05 | 3.7 | 2417 | — |
| 0.10 | 7.0 | 1181 | — |
| 0.20 | 13.4 | 532 | — |
| 0.50 | 29.4 | 174 | — |
| 1.00 | 20.8 | 47.7 | — |
| 2.00 | 22.1 | **19.7** | — |

Every bound covers the truth at every b. The two columns are monotone in **opposite
directions**, and the reason is elementary once seen:

* the **drift** contributes `μ·∫φf'g dt` against a band `κ√⟨M⟩ ∝ b`, so its
  signal-to-band goes as **1/b** — less noise is better;
* the **diffusion** contributes `½·b²·∫φf''h dt` against the same band, so its
  signal-to-band goes as **b** — *more* noise is better, because the thing being
  measured **is** the noise.

**So there is no single noise intensity at which both are determined**, and on this
system neither is determined at the crossover (both ≈ 20× relative width near
b ≈ 1–2). That is a structural tension between the two Level 1 targets, not a
tuning problem, and it is the sharpest statement this increment produced.

### What it means for S2, which was registered in the other picture

**S2** predicted: shortening Δt at fixed T tightens the *diffusion* interval and not
the drift's; extending T at fixed Δt tightens the *drift* and not the diffusion. In
the weak form the diffusion barely depends on Δt at all — nothing here estimates b
from quadratic variation, so the classic Δt → 0 demand does not apply to it. The
asymmetry S2 names is real but it lives in **which estimator you use**, not in two
intervals from one design:

| | efficient estimator | its precision |
|---|---|---|
| drift | the weak form (this) | improves with horizon T |
| diffusion | realized quadratic variation | relative error √(2/m), m = T/Δt |

At m ≈ 10⁴ increments quadratic variation pins b² to about 1.4%, where the weak
form's joint bound is 2000× the truth. **For the diffusion the weak form is a very
poor estimator and QV is the efficient one; for the drift it is the reverse, since
QV says nothing about drift at all.** S2 should be re-registered as a statement
about estimators before it is scored.

**The design conclusion:** the diffusion's symbolic form should be certified from
quadratic variation — regressing b²(x) = Σ d_j h_j(x) on the *local* realized QV as
a function of state, with its own band — not from the weak-form design matrix. The
`measure="d[u]"` terms step 3 landed are what express that, so the capability is
already in the vocabulary; what is missing is the regression and its band. That is
the next Level 1 increment, and it is now driven by a measurement rather than by
the curriculum's ordering.

## Two identifiability findings on the 1-D systems

**A bistable system hides its own drift.** The double well `a(x) = x − x³` was the
harder case and not for the reason first suspected. Its coefficients are wide
(width 18 on a truth of 1.0 at b = 1.4, L = 64) and shrinking the library to the
true support `{x, x³}` barely helps (width 16), so it is not library size. Nor is it
coefficient trade-off: the new `admissible_functional` bounds the drift **as a
function**, and `a(x)` is undetermined too (±4.9 at x = ±1, where its true value is
0, against a function whose range is ±6). The reason is structural — **the wells sit
where the drift vanishes, so a bistable process spends its time exactly where its
drift carries least information.** Density and information are anti-correlated by
construction.

**Multiplicative noise can disconnect the state space.** With `b(x) = b·x` the
volatility vanishes at the origin, so the origin is unreachable and every
trajectory is trapped in one well: measured zero well-crossings over T = 400 across
8 trajectories, against ~21 000 crossings for the additive case at comparable
noise. The drift is then determinable only on the visited well, which is the
registered third null ("trajectories do not visit enough of state space") arriving
as a property of a *system* rather than as a null. A domain-qualified claim
(`certify.domain_qualifier`) is the right output shape and is not yet wired here.

## A capability this needed: bounding the law, not the coefficients

`certify.admissible_functional` gives the joint range of any **linear functional**
of the coefficients over the consistent polytope — one LP per query. It exists
because per-component bounds can badly understate what the data determined, and it
is what turns "these four coefficients are each loose" into "the drift function is
in this envelope". Here it happened to confirm the pessimistic reading, which is
worth stating: it was built expecting to find hidden determination and it found
none.

It also matters for the frozen checker's semantics. Interval coverage is scored per
component (`STOCHASTIC_CHECKER.md` §5), so a system whose *function* is determined
while its coefficients are not would score as unresolved. That is the registered
scoring and this increment does not change it — but the case is now known to be
possible, and a functional-claim component kind is the obvious extension if a real
system exhibits it.

## What is not here

* **Van der Pol with additive noise** — 2-D, and the assembler is still scalar. The
  migration of `ito.build_rows` onto `build_nd` is the prerequisite and remains
  Level 1's other half. Multi-dimensional Itô also needs cross-variation
  `d[u,v]` and a per-field martingale sensitivity list, neither of which exists yet;
  `Term.measure` covers only the diagonal `d[u]` today.
* **The invariants target.** Untouched. The checker scores invariants up to affine
  reparametrization already, so this is generator work rather than instrument work.
* **A certified diffusion.** Nothing here certified b², and the honest reading is
  that the weak form is the wrong tool for it (above), not that b² is
  undetermined — quadratic variation determines it to ~1.4% on the same data.
* **Level 1 scored through the checker.** The Level 0 run does this; these
  measurements were taken directly, so there is no scored table yet and no S5
  statement for Level 1.

## Reproduce

```bash
.venv/bin/pytest tests/test_ito.py         # 22 green, incl. the joint-claim path
```
