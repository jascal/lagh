# Test-bed registration: econ-sae emergent-law discovery

**Registered 2026-07-21, before any target measurement.** The frozen instrument
(`lagh`: curriculum C1–C5, honesty core, acquisition policy) is unchanged; this document
supplies only the *target declaration* per the registration discipline. Amendments follow
the correction discipline: dated, appended, never edited in place.

**One-line:** point `lagh` at the `econ-sae` stock-flow-consistent macroeconomy as a black-box
oracle and ask it to recover the **emergent steady-state response laws** — functions of shock
inputs that no closed form in the source computes — certifying each or abstaining.

## 1. Why this target (against the six criteria)

| criterion | econ-sae |
|---|---|
| **exact law** | ✓ deterministic simulator; six accounting identities hold to machine precision (README) |
| **undocumented** | ✓ **for emergent aggregates only** — see the G7 firewall in §2 |
| **queryable** | ✓ `Simulator.step(shock) → (state, txns, dict)`; a run is init → N steps → read a steady-state aggregate |
| **noise / floor** | ✓ machine-precision determinism; replicates confirm σ≈0; floor = 1e-12 |
| **sealed GT** | ✓ the emergent law is computed by *running the sim*, never read from a formula; verification is a held-out shock grid |
| **abstention meaningful** | ✓ many emergent responses are genuinely non-closed-form or only locally analytic; `structural` / `range` refusal is the honest output and the expected majority outcome |

econ-sae is the first **zero-exposure, principled "exact-but-undocumented"** target — the original
`THESIS.md` class — living in the workspace, queryable today.

## 2. The G7 firewall (the discipline MDBench lacked)

The simulator contains **embedded closed-form parametric rules**: the Taylor monetary rule
(`i = base + π_w·(inflation−target) + y_w·gap`, source lines ~206, ~292), the balanced fiscal
rule, the MPC consumption rule, the sentiment-threshold multiplier. **Recovering any of these is
transcription, not induction — they are disqualified as targets**, exactly as CLRS was.

- **DISQUALIFIED (embedded, written-down):** the interest-rate response to inflation (Taylor), the
  transfer rule, per-household consumption as `mpc·disposable`, the sentiment step function. If a
  target's closed form appears in `simulator/core.py`, it is out.
- **VALID (emergent, no closed form in source):** a steady-state *aggregate* that is the fixed
  point of the whole multi-agent interaction and is computed by no single formula — the micro-rules
  are known; their *aggregate consequence* is not solved anywhere.

A target passes G7 iff a competent reader of `core.py` **cannot write its closed form** from the
source. That test is applied per candidate and recorded.

## 3. Registered target set (named before measurement)

Each target is a scalar `y = g(shock_params)`: run the sim to a steady state under a shock vector,
read one aggregate. Candidates, all `[VALID]` pending the §2 test at scoring time:

| id | shock inputs (the box) | emergent aggregate `y` | why not closed-form |
|---|---|---|---|
| **E1** | government spending multiplier `g ∈ [0.5, 2.0]` | steady-state GDP ratio vs baseline (the fiscal multiplier) | the multiplier is the fixed point of the full circular flow; no line computes it |
| **E2** | productivity shock `a ∈ [0.8, 1.25]` × sector | steady-state aggregate price level | emerges from firm pricing × household demand interaction |
| **E3** | credit-spread shock `s ∈ [0.0, 0.1]` | steady-state aggregate leverage | bank↔firm↔household credit loop fixed point |
| **E4** | two-input: `(g, a)` jointly | steady-state Gini across HH cohorts | distributional outcome of the coupled shocks |

- **Box:** log-uniform over each stated range (the registered `lagh` sampler). Multi-input targets
  sample the product box.
- **Steady state:** run `N_burn = 200` steps, then average the aggregate over `N_avg = 50` steps
  (frozen; a target that has not converged by then is reported, not silently averaged).
- **Determinism:** the sim is seeded per query; replicates at a fixed shock confirm σ. If σ>0
  (stochastic shocks), it feeds ε exactly as declared.

## 4. Scoring

- **Certified** ⇔ `lagh` returns a law and its exhaustive check passes on a held-out shock grid
  (20 points not used in fit/select), `nmiss = 0`.
- **dev-SA analogue:** a recovered law is *correct* iff it matches a dense-grid reference (200
  shock points, run fresh) to relative 1e-3 — the reference is the simulator itself, never a
  formula, so this is honest even when no closed form exists.
- **The zero-wrong invariant is the product** and carries over: a certified law that fails the
  dense-grid reference is a confident-wrong and caps the test bed. The record to defend is 0/all.
- **Abstention is a first-class, expected outcome** and is reported by reason, not scored as a
  loss — this is not a benchmark with a no-abstention protocol.

## 5. Predictions, registered before the run

- **PE-1:** E1 (fiscal multiplier) is recovered by C1–C2 (it is smooth and low-order in `g` over
  this range) — the tractable case.
- **PE-2:** E4 (Gini vs two coupled shocks) **abstains `structural`** — a distributional fixed
  point over a 2-D shock box is unlikely to be a low-complexity closed form, and honest refusal is
  the correct output.
- **PE-3:** **zero confident-wrong across all four** — every miss is an abstention. This is the
  invariant, now tested on a genuinely emergent (not manufactured) system.

## 6. What this can and cannot claim

A certified E1–E3 is a real result: *a closed-form emergent law of a multi-agent economy that the
simulator's own source does not contain, recovered from queries and certified over a stated shock
domain.* It is **not** a claim about real economies (econ-sae is a model), and it is **not** the
Standard-Model or bio ground-truth vocabulary (those are documented and fail G7 — recorded in the
menu discussion, not targeted here). Serving, if reached, is a `lagh` package per the design.

## 7. Build note

econ-sae is a published-package dependency, queried as an oracle — **not modified**. The adapter is
a thin `oracle(shock_matrix) → aggregate_vector` wrapper over `Simulator.step`, living in
`lagh/adapters/econsae.py`, that runs burn+avg per row. No econ-sae source is read beyond the
`step`/observable interface already inspected for this registration; the target law bodies (the
emergent aggregates) are unread by construction — they don't exist as source.
