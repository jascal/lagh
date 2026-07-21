# Test-bed registration: infinitely divisible laws via the Lévy exponent

**Registered 2026-07-21, before any run.** The frozen instrument (`lagh`, C1–C6) is reused; the new
piece is a **statistical certificate** — the significance layer of `DIRECTION_SIGNIFICANCE.md`,
pulled forward because a stochastic target requires it. Amendments dated, appended.

**One-line:** given samples from an infinitely divisible distribution, recover its **Lévy exponent**
`ψ(u)` — an *exact closed form that exists under genuine randomness* — from the empirical
characteristic function, certified at a stated significance, or abstain.

## 1. The target and why it fits

A distribution is infinitely divisible iff it is the marginal of a Lévy process; Lévy–Khintchine
gives every one an exact `φ(u)=exp(ψ(u))`. The exponent `ψ` is the law. **The exactness survives
the randomness** — this is the "exact but undocumented" class in a stochastic setting (a
mixed/compound Lévy exponent is written nowhere yet is an exact closed form).

**Symmetric case (v1):** `ψ(u)` is real and even, so the real observable `L(u)=log|φ̂(u)|=ψ(u)`
maps onto the existing real-function engine:

| law | `L(u)` | tier |
|---|---|---|
| Gaussian(σ) | `−σ²u²/2` | C1 |
| symmetric α-stable | `−c·\|u\|^α` (α rational) | C3 |
| symmetric compound Poisson (jump ±a, rate λ) | `λ(cos(a u) − 1)` | C4 |
| variance-gamma / Gamma-symmetrized | `−c·log(1 + u²/b²)` | C5 |
| **independent sum** of the above | **additive** | native |

## 2. The statistical certificate (the new, load-bearing piece)

You cannot exact-certify random data pointwise: `φ̂(u)=(1/n)Σe^{iuX_k}` carries `O(1/√n)` error. So:

- observation: `L̂(u)=log|φ̂(u)|` on a grid of `u`;
- **statistical `ε(u)` by bootstrap**: resample the data `B` times, `se(u)=std_b(L̂_b(u))`; feed
  `se` to the engine's existing `se` argument. Certify `|L_rec(u) − L̂(u)| ≤ λ_stat·se(u)` at every
  grid `u`. `λ_stat = 4` (z-score, matching κ) → `α`-level significance.
- This is a **proved-over-the-observation-set** claim in the `L2` sense: certified at significance
  `α` for this sample, **not** "proved the true law" — the honesty distinction, restated for the
  stochastic domain.

## 3. Targets (generated before scoring; exponents never hand-written into the recovery)

Symmetric ID laws, `n = 20000` samples, CF grid `u ∈ (0, U]`:

- **L1** Gaussian — `−σ²u²/2` (C1). Sanity.
- **L2** symmetric compound Poisson — `λ(cos(au)−1)` (C4).
- **L3** variance-gamma-like — `−c·log(1+u²/b²)` (C5).
- **L4** symmetric stable, **rational α** (e.g. 3/2) — `−c|u|^{3/2}` (C3). In-class.
- **L5** symmetric stable, **irrational α** (e.g. √2) — `−c|u|^{√2}`. **Out of class → abstain.**
- **L6** **mixed** = Gaussian + compound Poisson (undocumented sum) — additive `L(u)`.

## 4. Scoring

- **Certified** ⇔ `L_rec` certifies against `L̂` at `λ_stat·se` over the grid, min-domain guard met.
- **Correct** ⇔ `L_rec` matches the *true* `ψ` (known for the generator) to a stated tolerance on a
  fresh CF grid — the reference is the true exponent, computed, never fed to the recovery.
- **Zero-wrong invariant** carries over and is the product. A certified `L_rec` that fails the true
  exponent is confident-wrong. Standing record: **0 / ~180**.

## 5. Predictions, registered before the run

- **PL-1:** L1–L4, L6 recovered (Gaussian, compound Poisson, VG, rational-α stable, and the mixed
  sum) at significance; each `L_rec` matches the true `ψ`.
- **PL-2:** **L5 (irrational α = √2) abstains** — `|u|^{√2}` is not an exact rational-closed-form;
  the rational/irrational boundary transfers intact from the deterministic domain.
- **PL-3:** **zero confident-wrong.** The statistical `ε` + min-domain guard must not admit a
  spurious exponent; this is the first test of certification under genuine sampling noise.
- **PL-4:** the recovered mixed exponent (L6) is the additive sum, demonstrating superposition
  recovery — a genuinely *undocumented* Lévy exponent recovered and certified.

## 6. What this can and cannot claim

A certified recovery is: *the Lévy exponent of an infinitely divisible law, recovered from its
empirical characteristic function and certified at significance α over a stated CF grid, in exact
symbolic form.* It is **not** a claim about the true law with probability 1 (finite sample, stated
α), nor about non-symmetric ID laws (v1 is symmetric; the imaginary part / drift is v2), nor about
processes that are not infinitely divisible (those abstain — correctly).

---

## PARTIAL VERDICT (2026-07-21). Pipeline works; zero-wrong holds; needs a Lévy grammar.

- **Mechanics work:** CF-domain ranging (noise-floor `~1/√n`) + bootstrap statistical `ε` correctly
  turn a stochastic process into a certifiable real function `L(u)`.
- **PL-3 CONFIRMED (the invariant):** zero confident-wrong across all six targets — after the
  extended-coherence fix caught a masked spurious `x^(4/3)` recovery for Gaussian.
- **PL-1 NOT met yet:** all six currently *abstain*. Cause established, not guessed: the general
  fractional-power library (C1–C6) is too rich for the narrow achievable CF range — `x^(4/3)`,
  `x^(2/3)` etc. certify as impostors and the exponent is under-determined. More data does NOT fix
  it: `n` 20k→200k widens the reliable range only 1.83→2.16 (sqrt-log). The over-abstention is
  honest, not a bug.
- **The fix is a Lévy-specific grammar** (next build): restrict this target's curriculum to the
  Lévy-Khintchine exponent forms — `−a u²` (Gaussian), `λ(cos(bu)−1)` (compound Poisson),
  `−c log(1+u²/b²)` (VG/Gamma), `−c|u|^α` with **rational** α (stable) — plus their sums. This
  removes the impostors and makes the true form unique, exactly as C6 is a specialized tier for
  integer lattices. Irrational-α stable (L5) then abstains for the *right* reason (irrational
  boundary), not under-determination.
- **Net:** the honesty core generalized to the stochastic domain (statistical `ε`, ranging, the
  extended-coherence improvement that helps everywhere), and the remaining work is a focused
  domain grammar, not a soundness question. The zero-wrong record stands.
