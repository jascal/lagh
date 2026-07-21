# Test-bed registration: Ehrhart quasi-polynomial discovery

**Registered 2026-07-21, before any target measurement.** New capability (a modular /
quasi-polynomial tier) is registered here as part of the instrument, complete before the
run. Amendments dated and appended.

**One-line:** given a random rational polytope `P`, query the exact lattice-point count
`L_P(t) = #(t·P ∩ ℤ^d)` as an integer oracle, and recover its **Ehrhart quasi-polynomial**
— or abstain when the query budget cannot determine it.

## 1. Why this target — it fixes what every prior target strained

| criterion | Ehrhart |
|---|---|
| **exact law** | ✓✓ **integers.** No floats anywhere. The ε/noise/floor machinery (Corrections 5–7, 11, the N1′ machine-precision term) simply does not apply — certification is exact integer equality. |
| **closed form provably exists** | ✓✓ **Ehrhart's theorem** guarantees `L_P(t)` is a quasi-polynomial of degree `dim P`, period dividing the lcm of vertex denominators. This is the criterion econ-sae failed: the law is guaranteed to be in the class. |
| **undocumented** | ✓ nobody has written the quasi-polynomial of a *random* rational polytope; recovering it is induction, and computing it is a genuine (#P-hard in general) problem |
| **queryable** | ✓ `L_P(t)` for any positive integer `t` by exact enumeration |
| **sealed GT** | ✓ verify on held-out `t`, and on `t` far beyond the fit range; the "reference" is more counting, never a formula |
| **abstention meaningful** | ✓ a high-period quasi-polynomial needs many `t` per residue class; under budget the honest output is "underdetermined → abstain" |

## 2. The new capability (registered, generic — not per-target)

**Quasi-polynomial tier.** A quasi-polynomial is `L(t) = Σ_i c_i(t)·t^i` where each `c_i` is
periodic in `t` with some period `p`. Recovery, in **exact `Fraction` arithmetic** (the full C2′
ideal — no float ever touches it):

- search `(period p, degree d)` over `p ∈ {1..P_MAX}`, `d ∈ {0..D_MAX}`;
- split the fit `t` by residue class `t mod p`; per class, exact Lagrange-interpolate a
  degree-`d` rational polynomial through `d+1` points;
- **certify** iff every residue class's interpolant reproduces **every held-out `t` in that class
  exactly** (integer equality, `nmiss = 0`);
- **parsimony**: smallest `(p, d)` lexicographically — a true period-`p` law also "fits" period
  `2p` with duplicated classes, so the smallest certifying period is the true one;
- **coherence**: a degree-`d` polynomial is fixed by `d+1` points, so two quasi-polynomials that
  certify on ≥`d+1` held-out points per class are equal — coherence is automatic when the held-out
  set is adequate, and its inadequacy is exactly the abstention trigger;
- **abstain `range`** (budget analog): if the smallest true period needs more `t` per class than
  the budget provides, no `(p,d)` certifies → honest refusal.

## 3. Targets (generated before scoring; law bodies never written by hand)

Random rational simplices, dim 1–3, vertex denominators in {2, 3} (so periods ∈ {1,2,3,6}). The
oracle enumerates the bounding box in exact arithmetic. Ten polytopes per dimension, seeded.

- **H1** dim-1 intervals `[0, a/b]` — degree-1 quasi-polynomials, period `b`.
- **H2** dim-2 rational triangles — degree-2, period `lcm` of leg denominators.
- **H3** dim-3 rational simplices — degree-3.

Budget: up to `T_MAX = 48` distinct `t` values per polytope (the "queries"); fit/certify split
disjoint. `P_MAX = 12`, `D_MAX = 4`.

## 4. Scoring

- **Certified** ⇔ the recovered quasi-polynomial matches the oracle exactly on a held-out `t` set,
  `nmiss = 0`, `|D|` reported.
- **Correct (dense reference)** ⇔ it also matches the oracle exactly on `t ∈ {T_MAX+1 .. T_MAX+40}`
  — far beyond the fit range, so a merely-interpolating fit is caught.
- **Zero-wrong invariant** carries over and is the product: a certified quasi-polynomial that
  fails the extended-range reference is confident-wrong and caps the bed. Record to defend:
  the program's standing **0 / 118**.

## 5. Predictions, registered before the run

- **PH-1:** H1 and H2 recovered at essentially 100% within budget — periods ≤ 6 need ≤ `6·(d+1)`
  points, well under `T_MAX = 48`.
- **PH-2:** H3 (dim-3) mostly recovered; any misses are **`range` abstentions**, not wrong
  submissions — period 6 × degree 4 needs 30 points/class × 6 = beyond budget for the worst case.
- **PH-3:** **zero confident-wrong.** The extended-range reference cannot be gamed by
  interpolation, and exact-integer certification cannot admit an approximate law.
- **PH-4:** this is the **first target where no abstention is `noise`/`numerical`** — exact
  integers remove those failure modes entirely; only `range` (budget) and `structural`
  (underdetermined) can fire.

## 6. What this can and cannot claim

A certified H1–H3 is a real result: *the Ehrhart quasi-polynomial of a random rational polytope,
recovered from exact lattice-point queries and certified over a stated `t`-domain, in exact
arithmetic with no tolerance.* It is a clean-room demonstration that the honesty core (exhaustive
check, parsimony, coherence, first-class abstention) transfers to a **noise-free integer** domain
where it is at its strongest. It is **not** a claim about #P-hardness or general polytopes (small
dim, small denominators, bounded budget).


---

## VERDICT (2026-07-21). 30/30 recovered, zero abstain, zero confident-wrong.

| target | recovered (ext-range validated) | abstain | confident-wrong |
|---|---|---|---|
| H1 dim-1 | **10/10** | 0 | 0 |
| H2 dim-2 | **10/10** | 0 | 0 |
| H3 dim-3 | **10/10** | 0 | 0 |

Every recovery matches the exact oracle on `t = 49..88` — far beyond the fit range `1..48` —
so no result is interpolation. **PH-3 confirmed** (zero confident-wrong). **PH-4 confirmed**
(no `noise`/`numerical` abstention is even possible — pure integers). **PH-1 confirmed.**
**PH-2 over-delivered:** dim-3 recovered 10/10, not "mostly with `range` misses" — the
per-class self-split (below) uses the budget efficiently enough that period-6 degree-3 fits.

**First outright success in the program.** Set against the econ-sae verdict (0/4, all honest
abstention), the same instrument and same honesty core gives:

- **jagged, no closed form (econ-sae)** → refuses every time;
- **guaranteed closed form (Ehrhart)** → recovers every time.

That contrast *is* the product: abstention is meaningful because the instrument recovers exactly
when a law exists in its class and refuses exactly when one does not. Neither result is the
instrument being trivially cautious or trivially confident.

**Two defects found and fixed at cause (both over-abstained, never mis-answered — the safety
property):**
1. dim-3 Fraction enumeration was 12.5 s/query; cleared denominators to a pure-integer bounded
   count with the last axis closed-formed → 2 ms/query (6000×), still exact.
2. a *strided* fit/certify split aligned with period 2 and starved both residue classes; replaced
   by a **per-class self-split** (first `d+1` points of each class interpolate, the rest certify),
   which cannot starve a class and uses the budget maximally. Lifted H2/H3 from 7–8/10 to 10/10.

**The honest caveat — what this is NOT.** The recovery *uses* Ehrhart's theorem: it searches
exactly the class the theorem guarantees (quasi-polynomial, degree ≤ dim, small period). So the
induction is *"find the period, degree, and exact coefficients,"* not *"discover that the law is a
quasi-polynomial."* The functional *form* is theorem-given, not discovered — unlike the physics
targets where the form itself was unknown among C1–C5. This is a **certified-recovery** result in a
noise-free domain, deliberately, and it is claimed as exactly that: a clean-room demonstration that
the honesty core (exhaustive exact check, parsimony over period, first-class abstention) is at its
strongest when the arithmetic is exact. It is not a claim of open-form-discovery.

**Record:** the program's zero-wrong invariant now stands at **0 / 148** scored tasks, and Ehrhart
is its first domain of outright recovery.
