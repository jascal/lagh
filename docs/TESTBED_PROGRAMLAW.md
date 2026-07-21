# Test-bed registration: program-law recovery in the wild

**Registered 2026-07-21, before any target function is scored.** The frozen instrument
(`lagh`, C1–C6) is unchanged; this supplies the target declaration and — critically — a
**mechanical function-selection protocol** so functions cannot be cherry-picked for
recoverability. Amendments dated and appended.

**One-line:** treat a real numeric function as a black-box oracle — call it, never read it —
and recover a certified closed form or abstain; read the source *only after* submission, to
score. This is the product demo: *"point lagh at a function you can call but haven't read."*

## 1. Why this target and what is genuinely new

Ehrhart proved the honesty core in an exact-integer clean room where the law's *form* was
theorem-given. Program-law recovery is the first target where the **functional form is unknown**
(any of C1–C6, or none) *and* the ground truth is a real artifact rather than a generated one. It
directly tests the whole curriculum's reach on functions written for other purposes.

| criterion | program laws |
|---|---|
| exact law | ✓ deterministic code (replicates confirm σ=0) |
| closed-form-expressible | **unknown per function — that is the test.** Many functions loop / branch / table-look-up and have NO closed form in C1–C6; the honest output there is abstention (the econ-sae lesson operationalized) |
| undocumented | ✓ the closed form of an arbitrary numeric function is written nowhere but its source, which the learner never reads |
| queryable | ✓ free, unlimited, deterministic |
| noise / floor | ✓ float semantics already modelled (machine-precision ε term) |
| sealed GT | ✓✓ **cleanest possible** — verification is reading the source *after* submission |
| abstention meaningful | ✓ non-closed-form functions are the majority; refusing on them is the product |

## 2. The function-selection protocol (frozen — the anti-cherry-pick firewall)

The danger is obvious: pick functions I know lagh recovers. So selection is **mechanical and
declared before any function is read for recoverability**:

- **Source pool:** the numeric standard library + numpy, restricted to functions with signature
  `float^k → float`, `k ∈ {1, 2, 3}`, `k` known from the signature only.
- **Candidate list, frozen below** — chosen by *name and arity alone*, before probing any for
  recoverability, spanning known-easy, known-hard, and unknown:

  | id | function | arity | a-priori expectation (not a selection criterion) |
  |---|---|---|---|
  | F1 | `math.hypot(x, y)` | 2 | recoverable (√(x²+y²), C5) |
  | F2 | `math.log2(x)` | 1 | recoverable (log/ln ratio) |
  | F3 | `math.gamma(x)` | 1 | **abstain expected** — no elementary closed form |
  | F4 | `math.erf(x)` | 1 | **abstain expected** — non-elementary |
  | F5 | `math.atan2(y, x)` | 2 | branch/quadrant — abstain expected off a single sheet |
  | F6 | `numpy.logaddexp(x, y)` | 2 | recoverable (`log(e^x+e^y)`, C4/C5) |
  | F7 | `math.expm1(x)` | 1 | recoverable (`e^x − 1`) |
  | F8 | `math.copysign(x, y)` | 2 | discontinuous — abstain expected |
  | F9 | `math.remainder(x, y)` | 2 | periodic/branch — abstain expected |
  | F10 | `numpy.hypot`-of-3 `√(x²+y²+z²)` | 3 | recoverable (C5) |

  The a-priori column is **recorded, not used to select** — the point is a mix, and the
  predictions (§4) are scored against it.
- **Domain box:** each function's declared safe positive-real domain, log-uniform, stated per
  function in the runner (e.g. gamma on `[0.5, 6]`, log2 on `[0.1, 100]`). Adaptive ranging on.
- **Budget:** the registered acquisition policy, ≤ 200 queries.

**No function's source is read until after its result is committed.** For stdlib/numpy the "source"
is the mathematical definition, known — so scoring uses a **held-out dense reference** (a fresh
grid, 200 points, matched to relative 1e-9) rather than source inspection, which is stricter and
oracle-only.

## 3. Scoring

- **Certified** ⇔ `lagh` returns a law certified over its held-out query domain.
- **Correct** ⇔ it also matches the oracle on a fresh dense grid to relative 1e-9.
- **Zero-wrong invariant** is the product: a certified law failing the dense grid is
  confident-wrong and caps the bed. Standing record to defend: **0 / 148**.
- Abstention reported by reason; expected for F3–F5, F8, F9 and *not* scored as failure.

## 4. Predictions, registered before any run

- **PP-1:** F1, F2, F6, F7, F10 recovered (elementary closed forms in C1–C6). F2 (`log2`) tests
  whether the constant `1/ln2` is recovered as an exact rational or forces abstention — recorded
  either way.
- **PP-2:** F3, F4 (gamma, erf) **abstain** — no elementary closed form; refusal is correct, and a
  *certified* answer here would be a confident-wrong (they are the invariant's real test).
- **PP-3:** F5, F8, F9 (atan2, copysign, remainder) **abstain** — branch/discontinuity/periodicity
  outside a single analytic sheet.
- **PP-4:** **zero confident-wrong across all ten.** This is the first test of the invariant where
  the target's recoverability is genuinely unknown per function and non-elementary functions are
  deliberately included as traps.

## 5. What this can and cannot claim

A certified F1/F6/F10 is a real result: *the closed form of a function recovered from black-box
calls, certified, without reading its source.* An abstention on F3/F4 is equally the product: the
instrument declining to fabricate an elementary form for a non-elementary function. What it is
**not**: a claim that lagh recovers arbitrary programs (it recovers what its curriculum reaches);
the traps are included precisely to bound that claim honestly.


---

## VERDICT (2026-07-21). 3 recovered, 7 abstain, zero confident-wrong. Two prediction
## misses, both diagnosed to specific frozen-instrument limits.

| func | outcome | law / reason |
|---|---|---|
| F1 hypot | **CERT ✓** | `sqrt(x_0**2 + x_1**2)` (C5) |
| F2 log2 | **CERT ✓** | `1206321·log(x_0)/836158` = `log(x)/ln2` — `1/ln2` recovered as exact rational |
| F3 gamma | abstain `structural` | correct — no elementary closed form |
| F4 erf | abstain `structural` | correct — non-elementary |
| F5 atan2 | abstain `structural` | correct — branch/quadrant |
| F6 logaddexp | abstain `structural` | **predicted-recoverable MISS — curriculum gap** |
| F7 expm1 | **CERT ✓** | `exp(x_0) - 1` (C1) |
| F8 copysign | abstain `structural` | correct — discontinuous |
| F9 remainder | abstain `structural` | correct — periodic/branch |
| F10 hypot3 | abstain `structural` | **predicted-recoverable MISS — acquisition budget** |

**Zero confident-wrong (PP-4 confirmed).** PP-2 (gamma, erf abstain) and PP-3 (atan2, copysign,
remainder abstain) **confirmed** — the non-elementary traps were declined, not fabricated. This is
the invariant's hardest test yet: five deliberate traps whose recoverability was unknown per
function, and the instrument refused all five without a single confident-wrong.

**PP-1 partially falsified — informatively, in two different ways:**

1. **F10 hypot3 — acquisition-policy limit, not curriculum.** Bare `discover` (90 points) recovers
   `sqrt(x_0**2 + x_1**2 + x_2**2)` exactly. The full run abstained because the registered
   acquisition `init_points = 40` is too few for a 3-arity square-transform target (24 fit points
   vs. a 10-term quadratic library). **The curriculum reaches it; the frozen policy's init budget
   does not.**
2. **F6 logaddexp — genuine curriculum gap.** `log(e^x + e^y)` needs an **exp-of-target transform**
   (`exp(y) = e^x + e^y`, a clean C4 sum), and C5's transform set is `{square, 1/y, log, sin}` —
   **exp is absent.** Bare discover fails too. This is the honest curriculum boundary the traps
   were meant to probe, hit by a function that *does* have a closed form we simply don't reach yet.

**Both misses are honest under the frozen instrument** and neither is retuned to force the number.
Two disclosed, exposure-recorded improvements for a *future* registered instrument version (not
applied here):
- add `exp` to the C5 transform set (generic — closes F6 and the whole log-sum-exp family);
- raise or adapt `init_points` with arity (generic — the active loop should escalate data when
  nothing certifies, rather than abstain on thin init).

**Record:** zero-wrong invariant now **0 / 158** scored tasks. Program-law recovery is the first
target with genuinely-unknown per-function recoverability, and the instrument sorted ten functions
into recovered / honestly-refused / diagnosed-limit with no fabrication.
