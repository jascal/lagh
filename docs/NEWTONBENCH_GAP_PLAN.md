# NewtonBench algebraic-gap plan — REGISTERED 2026-07-21, before any grammar change

The 33 addressable algebraic gaps (`NEWTONBENCH_CEILING.md`) decompose into **6 capabilities**, not
3. Each is a *structural* grammar class (never per-target), registered here with the exact cells it
predicts it will unlock **before** the code exists — a prediction I score against, not tune toward.
Zero-wrong is the gate: any capability that produces one confident-wrong cell is reverted, not kept.

Ground rule: a capability is scored PASS only if (a) it unlocks its predicted cells, (b) it breaks no
previously-recovered cell, (c) it adds zero confident-wrong across the whole 108. Register → implement
→ isolated verify (predicted module + a regression module) → targeted re-sweep → score.

## CAP-B — input × trig-monomial products / trig ratios  ✅ SCORED PASS (2026-07-21)
**Result under the real NewtonBench harness (noise + active acquisition):** malus **1/9 → 6/9** —
recovered easy_v0/v1/v2, medium_v0/v1, hard_v0 (the 5 predicted + easy_v0); the 3 `^e` cells abstain
structural; **0 confident-wrong.** Prediction met exactly. Registered bound: emitted only for dim≤2
(the features are admissible on any input, so dim≥3 blows up the C2 search; all trig gaps are dim-2).
NewtonBench algebraic 55→60 / 88. Below is the pre-registration, unchanged.

Feature: `x_i · sin^a(x_j) cos^b(x_j)`, small `a` and signed `b` (covers `(α sinθ+β cosθ)²` →
{sin², cos², sincos}, `tan²=sin²cos^{-2}`, `cot²=cos²sin^{-2}`, `sin²/cos³`). Guarded finite.
**Predicts (5):** m7_malus easy_v1 (`I₀tan²θ`), easy_v2 (`I₀cot²θ`), medium_v0 (`I₀(2sinθ+cosθ)²`),
medium_v1 (`I₀sin²/cos³`), hard_v0 (`I₀(2sinθ+1.5cosθ)²`). Must still abstain: malus medium_v2/hard_v1/
hard_v2 (`^e`, irrational).

## CAP-A — wider fractional-power grid + input×fractional products  ◐ PART 1 SCORED (+2)
> **Part 1 (2026-07-21): extend C3 exponent denom caps 4→{3,5,10}.** Data-driven exponent snapping
> from the log-log slope, each cap one checked candidate (no feature inflation). **Scored under noise:
> hooke 5→6/9 (easy_v2 `x^{17/5}`), sound 6→7/9 (hard_v2 `(RTM^1.5)^{-2.8}`), zero confident-wrong.**
> Total 64→66/108. **Part 2 DEFERRED** — the multi-term hooke *sums* (medium_v2/hard_v1/hard_v2) need
> additive C1 fractional features (`x^{17/5}`, `x^{-3/10}`, `x^{-10/3}` as summands), which inflates
> the feature set combinatorially; not added until an exact-coefficient gate is in place (CAP-E lesson).

Add denominators 5, 10, 3 to `FRAC_EXP` (`x^{17/5}=x^{3.4}`, `x^{-3/10}=x^{-0.3}`, `x^{-10/3}`), and the
`x_i·x_j^{p/q}` products over the same grid. **Risk: inflates the feature set globally → impostor/zero-
wrong exposure; the extended-coherence probe must hold.** Verify zero-wrong on ≥3 unrelated modules.
**Predicts (5):** m9_hooke easy_v2 (`2k x^{3.4}`), medium_v2 (`+K₂x^{0.5}`), hard_v1 (`x^{0.5}+x³+x^{-0.3}`),
hard_v2 (`x^{3.4}+x^{0.5}+x^{-10/3}`); m8_sound hard_v2 (`(RTM^{1.5})^{-2.8}` = product of denom-5 powers).

## CAP-C — richer monomial products (poly × negative-power, higher cross-degree)
`x_i²x_k²`, `x_i²x_j·x_k^{-2}`, and degree-≤3 sum-expansions × `r^{-2}`. **Predicts (5):** m0_gravity
hard_v0 (`(m1+m2)²/r^{1.5}`), hard_v2 (`(m1²+m2²)r²`); m1_coulomb medium_v0 (`q1q2(q1+q2)/r²`), medium_v1
(`(q1+q2)³/r²`), hard_v1 (`q2²(q1+q2)³/r²`).

## CAP-D — rational base raised to a power, via a `sqrt` target-transform
> **RESOLVED 2026-07-21 — it was REACH, not box. ✅ implemented.** The first probe misread the
> "structural" abstain as multiplicity; it was the *terminal no-candidate* case. Root cause: **C2's
> pure-term-denominator pass returned 0 candidates** for `k/m − b²/4m²` even though `x_1²`, `x_0x_1`,
> `x_2²` were all present. Two bugs, both fixed:
> 1. **C2 `denom_idx` string bug** — excluded any term whose name contains `*`, but `x_j**2` contains
>    `**`, so integer-power denominators were never tried. Fixed to `free_symbols == 1` (single
>    variable). *General reach fix* for every `P/x_j**p` rational — not underdamped-specific.
> 2. **CAP-D `sqrt`-forward transform** added (fit `√y`, invert by squaring) so `(rational)²` reduces
>    to the C2 base under the root.
> Result: `easy_v2` recovers `(x_0x_1 − x_2²/4)/x_1²` exactly; `(rational)²` variants recover via
> `sqrt`→C2. Box-widening was NOT the fix — the abstains were never box-shaped. Scoring the 7 cells
> under noise: targeted re-sweep pending.

Add the missing `sqrt` transform (forward `√y`, invert `·²`) so `(rational)²` reduces to a C2 rational
base. **Predicts (up to 5):** m6_underdamped easy_v1/v2, medium_v1, hard_v1 (the squared/plain rationals
`(k/m − (b/2m)²)^{1 or 2}`). Fractional-power variants (hard_v0 `^1.5`, hard_v2/medium_v2 `m^{1.3}`) are
weaker — may need CAP-A too; registered as *stretch*, not core.

## CAP-E — rational-of-exponential `1/(e^{u}−1)` (Bose–Einstein)  ❌ REVERTED — UNSOUND AT SCALE
> **Attempted + reverted 2026-07-21.** A `bose` transform (forward `log(1/y+1)`, invert `1/(e^e−1)`)
> recovered all four BE forms *in isolation* (nice `C=1`) but **produced a confident-wrong on the real
> oracle** (medium_v1). Two coupled causes: (1) `log(1/y+1)` is numerically ill-conditioned for large
> `y` — the real BE outputs are ~1e10, so `1/y+1 ≈ 1` and catastrophic cancellation gives a ~5-digit
> `τ`; (2) the power-law tier emitted the resulting coefficient (~1.05e-14) as a **raw float, unsnapped**,
> and at output scale ~1e10 the relative `eps` was loose enough for that float form to certify while
> being 3e-5 off. The real BE coefficients are ~1e-14 — **un-representable as exact rationals** — so the
> transform recovered *zero* real BE cells and only manufactured a false certificate. **Reverted; the
> invariant is the product.** Re-introducing BE safely needs BOTH `np.log1p(1/y)` stability AND an
> exact-coefficient gate (reject certified candidates carrying un-snappable `Float` coefficients) — a
> separate, carefully-verified task. **LESSON: target-transforms that compress the signal (`log`, `inv`,
> `bose`) can spuriously certify at extreme output scale; prefer additive polynomial/power features
> whose coefficients are exact rationals by construction (lstsq→snap), and add an exact-coefficient
> gate before any more transform capabilities.**

### (original CAP-E pre-registration, superseded by the revert above)
The `inv` transform already yields `e^u − 1`; needs C4 inner `e^{C·ω^p/T^q}` with fractional inner
powers. **Predicts (up to 5):** m10_be easy_v1/v2, medium_v1/v2, hard_v1 (`1/(e^{…}−1)` with rational
inner exponents). hard_v2 (`1/(−ln(...)−1)`) is a different inner (log) — stretch.

## CAP-F — transcendental(input) × power products `e^{γ}`, `ln(γ)` × monomials
**Predicts (2):** m8_sound hard_v0 (`√(e^γ RT²/M^{1.5})` = `e^{γ/2}R^{0.5}TM^{-0.75}`), hard_v1
(`ln(γ)TR/M^{1/3}`). m5_decay easy_v2 (`e^{−(λt)^{0.5}}`) is exp-of-fractional-inner — stretch.

## Out of scope (correctly abstained, NOT gaps)
- 11 irrational-exponent cells (`^e`) — the boundary; must stay abstained (the wedge).
- m2_magnetic hard_v0 `(I1+I2)^{1.5}` — fractional power of a **sum**, not a monomial; genuinely hard,
  left abstained unless a composite-power capability is later justified.

## Predicted totals if all core capabilities land
CAP-B 5 + CAP-A 5 + CAP-C 5 + CAP-D 4 + CAP-E 5 + CAP-F 2 = **26 of 33**, lifting algebraic 55→~81 / 88.
Ceiling stays 88 (inverse-trig tier is separate, +9 → 97). Every number above is a prediction to be
scored, not a claim.
