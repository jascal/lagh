# Direction: characterization-on-abstain + the research playbook

Status: **proposal + first implementation** (2026-07-22). Companion to
`DIRECTION_TOOLSHAPE.md` (the `recover`/`verify`/`fit` verbs) and `RNOISE_STUDY.md`
(the parametric gate). Motivated by the full-108 composite run (Grok-4.5 + lagh:
67/108, 0 confident-wrong, but 41 abstain cells timed out because a *bare* abstain
gave the model nothing to act on).

## The problem a bare abstain creates

lagh's abstention is sound but **mute**: it returns a reason enum (`structural`,
`range`, `parametric`, …) and stops. In the composite loop that muteness costs twice:

1. **It wastes the LLM.** With no diagnosis, the model either gives up (a MISS — value
   left on the table) or thrashes — re-samples and re-calls `discover` hoping to break
   through, blowing the wall-clock (the 41 snell/BE/decay timeouts). A bare abstain is
   an *invitation to thrash*.
2. **It hides what lagh actually knows.** After discover exhausts, lagh has measured a
   great deal about the function — its trend, its power-law slope, whether that slope is
   a clean rational or an irrational reaching for `e`. Throwing that away is the
   difference between "I can't" and "I can't, *and here's precisely why, and here's what
   would help*."

## The fix: a graceful-degradation ladder, every rung zero-wrong

Replace certify-or-silence with a three-rung ladder, tagged by the existing tag
discipline (`proved` / `empirical` / `open`):

| rung | tool result | tag | guarantee |
|------|-------------|-----|-----------|
| **certify** | exact law + `strength` (pinned/consistent) | `proved` | exhaustive `|f−y|≤ε` |
| **characterize** | sound, UNCERTIFIED structural diagnosis | `empirical` | grounded observations, labeled as a guess |
| **abstain** | reason enum only | `open` | nothing claimed |

The middle rung is new. It is **strictly more information than a bare abstain at exactly
the same safety** — because it never carries a `certified` field and never a `law`, only
a `characterization`. It cannot become a confident-wrong: there is no certificate to be
wrong *about*.

## The characterization contract (tool-side, `lagh/characterize.py`)

`characterize(X, y, sigma, abstain_reason)` returns a dict that **only reports what the
samples measure**, each observation heuristic and labeled:

- **`shape`** — per-input trend (rank-correlation sign: increasing / decreasing / none),
  output sign, and boundedness (does `|y|` stay bounded as inputs grow).
- **`power_law`** — the centerpiece. A log-log least-squares probe (valid when `X>0,
  y>0`) yields an exponent per input, each snapped to a small rational. Then the
  **exponent-rationality flag**:
  - all exponents pin to small rationals, low residual → clean power law.
  - an exponent does **not** pin, but sits within `5e-3` of a named constant
    (`e`, `π`, `√2`, `1/e`, `ln2`, `golden`) → **`irrational_hint`**. This is the
    NewtonBench **wedge** turned into a *positive diagnosis*: "exponent ≈ e; no exact
    rational closed form exists — any exact algebraic law here would be a lie."
  - high log-log residual → not a pure power law (additive/polynomial/periodic).
- **`class`** — synthesis: `power-law` | `irrational-power` | `exponential`
  (log-linear fit) | `non-algebraic` (bounded output + poor power-law fit → possible
  trig/inverse-trig/saturating) | `unresolved`.
- **`why`** — one human-readable line.
- **`research`** — the **playbook pointer**: the single next move (see below).

**Hard invariant:** the return has `certified: False` and `tag: "empirical"`, no `law`
key. It lives entirely in the `fit` lane. The skill must relay it as a hedge, never as
a result. Soundness rule: only claim measurables; never snap to a specific wrong form.

## The research playbook (skill-side) + the stop rule

A bare abstain makes the model thrash; a *diagnosed* abstain lets it act — or stop. The
skill gets a decision tree keyed on `research.move`:

| `research.move` | when | what the LLM does |
|-----------------|------|-------------------|
| `acquire_divergent` | under-determined (`structural`; rivals fit the thin range) | sample a **wider / asymptotic** regime where the rival forms separate, then `recover`. *(This is the underdamped-v0 gain path.)* |
| `acquire_more_data` | `range` / thin / needs-data | sample more high-leverage points, widen 2–10×, `recover`. |
| `declare_and_verify` | `non-algebraic` / bounded-trig | declare a form from priors (physics hint, trig) and `verify` — accept a `consistent` hedge. |
| `report_and_stop` | `irrational-power`, or nothing pins after the moves | **report the characterization as a hedged answer and abstain. Do NOT keep sampling.** |

`report_and_stop` **is the timeout fix**: on a genuinely hopeless cell (irrational
exponent, out-of-grammar), the diagnosis tells the model to quit *directed* instead of
thrashing to the wall. And a bounded budget — **at most K research moves, then
`report_and_stop`** — caps the cost even when the model wants to keep trying.

## The spine: bold proposer, sound checker

The reason we can push the model to *explore aggressively* is that lagh is a sound
checker and the model is an untrusted proposer. Every form it invents is run through
`verify`/`recover`; a wrong guess is rejected, never certified. So the honesty
constraint lives at the **output** (report only what a tool certified, plus the labeled
hedge), not at the exploration. The instruction to the model is: **explore boldly;
report only what's certified.**

## Guardrails

1. **Bounded** — the K-move stop rule; `report_and_stop` on hopeless diagnoses. (Fixes
   the timeout.)
2. **Honest** — the characterization is `empirical`, never a certificate; a `fit`
   conjecture is never an answer.
3. **General, not benchmark-specific** — the playbook is a *research method* ("on
   under-determination, sample where rivals diverge"), never encoded answers ("snell is
   arcsin"). This is the line guarded hardest: benchmark-specific hints would collapse
   the honesty story into overfitting.

## Honest scope

On NewtonBench this does **not** move the exact-CORRECT count much: ~20 of the abstains
are irrational-`^e` (impossible by design) or inverse-trig (out of grammar) — those
become *richer hedges* (`verify`→`consistent`, or a characterization), not new
CORRECTs. The extra exact recoveries come from the under-determined / reach-bounded
cells (fewer). The real payoff is (a) killing the timeout so abstains are fast and
clean, and (b) open/real problems, where "under-determined, fixable by smarter
sampling" is far more common than NewtonBench's adversarial out-of-class walls — there
the playbook makes the LLM a genuine force multiplier.

## Validation (2026-07-22, `lagh/characterize.py`, `tests/test_characterize.py`)

Ran the characterizer over all 42 NewtonBench-dev abstain cells + the stress cells:

- **Snell → `non-algebraic` → `declare_and_verify`** (bounded, log-log residual ~0.49 —
  correctly reads as trig/inverse-trig, not power law).
- **underdamped-v0 → `acquire_divergent`** — *exactly the gain path* Grok used to recover
  it via a wider hand-sample.
- **Class distribution over the 42:** unresolved 19, additive-or-mixed 13, non-algebraic
  8, power-law 2 — every abstain now carries a sound structural read + a next move.
- **Invariant held on every real cell:** `certified: False`, no `law` key.

**A soundness bug caught and fixed here:** the first pass flagged 3 cells as
irrational-`^e` "≈ ln2 / √2" — all **false positives**: the log-log "exponent" of a
compound/multi-term law is meaningless and coincidentally landed near a named constant
(the flagged set even *changed* with the random sample). Fix: the irrational-hint only
counts when the monomial fit residual is low (a genuine clean irrational monomial like
`x**e`, residual ≈ 0, is still caught — unit-tested). This is the "under-claim rather
than mislead" rule made concrete.

**Honest scope confirmed:** the wedge detector fires only on *clean* irrational monomials.
NewtonBench's 11 `^e` cells are **compound** (not pure `x**e`), so they read as
additive/non-algebraic, not the specific "irrational-power" — a sound hedge, just not the
sharp label. Zero false wedges is the win; catching compound irrationals is future work
(needs residual-structure analysis, not a single log-log slope).

## Test claim still to close (live)

A diagnosed, playbook-guided abstain is **both** faster (kills the 150s thrash) **and**
yields **more gains**, with **zero** new confident-wrong — to be measured end-to-end on a
composite re-run once the guidance is deployed. The tool-side pieces are validated above;
the end-to-end composite measurement is the remaining step.
