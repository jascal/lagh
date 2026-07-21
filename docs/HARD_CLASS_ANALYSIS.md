# NewtonBench-dev hard difficulty: the genuine new function class (2026-07-21)

Read from the hard law bodies (dev, read freely). The hard mutations are NOT box-degenerate
versions of in-class laws — they introduce a genuine new class that splits in two.

## Half 1 — fractional-RATIONAL exponents (tractable, exact-recoverable)

Rational exponents C3 currently misses because it snaps only to denominators {1, 2, 4}:

| module | hard law fragment | exponent |
|---|---|---|
| sound | `T**2 / M**1.5`, `(...)**-2.8` | 3/2, -14/5 |
| harmonic | `k/m**1.3 - (b/2m)**0.7` | 13/10, 7/10 |
| snell | `n2*sin(θ)/n1**2.5` | 5/2 |
| coulomb | `q2**2 * (q1+q2)**3` | integer 3 (in-class already) |

**Fix:** extend C3's exponent search to finer rationals (denominator ≤ ~12). Preserves EXACT
recovery; recovers this half outright. Small, principled (`snap_small` already does this for C4
inner scales — reuse the pattern).

## Half 2 — IRRATIONAL exponents + new operators (a fundamental boundary)

| module | hard law fragment | why unreachable |
|---|---|---|
| coulomb | `/ r**e` (Euler's e) | **irrational exponent** — no exact rational-closed-form |
| decay | `N0**1.2 * e**(-λ**e * t**1.5)`, `(λt)**(e+1.5)` | nested irrational exponents |
| sound | `(e**gamma) * ...`, `ln(gamma)*...` | irrational base/exponent |
| snell | `atan((n1/n2)**2 * tan(θ))` | **tan / arc-trig** — new operators |

**`x^e` is not an exact rational-closed-form (C2′).** No finite exact certificate exists for it
over a wide box. lagh MUST abstain — and that abstention is correct, not a defect. This caps hard
R-cap by construction: a fraction of hard vanilla cells are outside the certifiable class.

- `tan` / `atan` are addable as C1/C5 operators (a real, bounded capability add).
- **Irrational exponents are the boundary.** lagh does not approximate; it certifies exact forms
  or abstains. This is the §3 distinction (`wyly/GOAL.md`): NewtonBench's SA is LLM-judged
  structural equivalence *with constants disregarded* — it would accept `x^e` as "a power law".
  lagh's certificate demands the exponent exactly. Different objects.

## Strategic consequence (the tool-shape argument, sharpened)

On irrational-exponent hard cells, **lagh abstains (exact-honest) while an LLM recovers the
structure** (constants-disregarded). Neither alone covers hard; the **LLM + lagh composite** does —
lagh certifies the exact-rational cells (including Half 1 after the C3 extension), the LLM handles
the irrational-structure cells, and the composite guarantee (`≥ LLM_alone`) holds because lagh
never fabricates an exact `x^e`. This is the strongest concrete case for `DIRECTION_TOOLSHAPE.md`:
hard difficulty is exactly where a text-blind exact solver and a text-fluent approximate one are
complementary rather than competing.

## Plan

1. **C3 finer-rational exponents** (denominator ≤ 12) — recovers Half 1, preserves exact recovery.
2. **tan/atan operators** — bounded add, recovers the snell arc-trig subclass where the argument
   is in-class.
3. **Irrational exponents: accept as the boundary.** Report the abstention rate on hard as a
   measured property, not a gap to close — and route those cells to the LLM in the composite.
4. R-cap is re-evaluated with (1)+(2); if hard still < 10/12 due to irrational exponents, that is
   the honest ceiling for lagh-alone, and the composite is the path.

## Design principle: DO NOT add an irrational-exponent class (load-bearing)

The rational/irrational exponent line is the exact-certifiable boundary, not a limitation to
overcome. Adding irrational support (`x^e`) is refused, for five compounding reasons:

1. **Breaks "certified."** `e` is non-terminating; supporting it means approximating it, which turns
   the certificate from "holds exactly" into "fits within tolerance for this approximation" — i.e.
   approximate SR, the thing lagh exists not to be. It would rebuild PySR and discard the only
   distinguishing claim.
2. **Destroys the significance guarantee.** `α ≤ |H|·q^h` needs `|H|` finite/countable (bounded-
   denominator rationals). A continuous real exponent is uncountably many hypotheses; `|H|→∞`,
   significance collapses — the multiple-testing failure the min-domain guard just fought.
3. **Forces the nonlinear search we refused.** A real-exponent fit is a search, not a solve, and a
   continuous-parameter residual-minimizer ALWAYS fits within tolerance → it never abstains,
   deleting the product.
4. **The abstention is a TRUE statement.** On a genuine `x^e`, "not an exact rational-closed-form"
   is correct. Supporting it replaces true refusal with fabricated certification of an
   approximation — the confident-wrong direction.
5. **No principled stop past the line.** Irrational exponents → irrational coefficients →
   transcendental compositions → algebraic numbers → …; rational/irrational is the only principled
   boundary where "certify" has content.

**The coincidence is the design cohering:** fractional-RATIONAL exponents stay in-class (add them,
exact); irrational ones are out (abstain). The composite covers `x^e` via the LLM's
constants-disregarded structure judgment. lagh's abstention there is the value working, not failing.
