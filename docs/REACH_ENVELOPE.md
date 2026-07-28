# The measured reach envelope

**Audit registered and run 2026-07-28** (`experiments/reach/audit.py`,
results `experiments/results/reach_audit.json`): 36 cells of
(dimension × form family), each a clean synthetic recovery at machine floor,
n = 400, inputs uniform on [0.5, 3], per-cell stable seeds. The point: every
reach gap should be a **stated cap**, not a benchmark-day surprise (the RRab
lesson — proposal reach, not certification, is the binding constraint).

## Verdict: 33/36 certified — 35/36 after significance arbitration (2026-07-28)

| family | cells | result |
|---|---|---|
| linear, 2–7 terms, dims 1–6 | 6 | all certify (≤ 74 s) |
| complete quadrics, dims 2–5 | 4 | all certify — CAP-T's d ≤ 4 cover plus support search carries d = 5 |
| cubics incl. cross terms, dims 1–3 | 3 | all certify |
| monomials incl. fractional powers, dims 1–4 | 5 | all certify |
| rationals (affine denominators), dims 1–3 | 3 | d2/d3 certify; d1 abstains (boundary-class, below) |
| transcendental sums/products | 7 | 6 certify; mixed-4term abstains (boundary-class) |
| sparse sums, sizes 3–6 | 4 | sizes 3–5 certify; size 6 abstains (boundary-class) |
| compositions (gaussian, inner-linear trig, √Σx², 1/Σx²) | 4 | all certify (the Σx² forms are slow: 13–21 min) |

## Two gaps found and CLOSED (commit-level, suite 76 green)

1. **Negative integer powers were absent from the library** — x^(−1/2) and
   x^(−3/2) existed, plain 1/x and 1/x² did not; `3/x` in a sparse sum was
   unreachable at any dimension. Closed: p ∈ (−1, −2) features, guarded.
2. **The rational channel's enumeration died before its second
   denominator** — exhaustive C(n,2) numerator pairs (~25k per denominator
   at dim 3) against a 40k attempt budget meant `(x₀+x₁)/(x₂+2)` was
   unreachable because x₂ was never *tried*. Closed: complexity-ordered
   denominators; pair pool = top-14 singles ∪ the low-complexity core (the
   top-singles ranking measurably favors composite features — {x₀, x₁}
   never entered the old pool).

## The boundary finding: Müntz twins at machine floor

The three remaining abstains are not reach failures. Diagnosis (worked
example in the audit history): on a smooth target over a narrow positive
range, an ~8-term combination of non-integer powers (x^17/5, x^4/3, …) can
fit within MACHINE epsilon on a data split — Müntz-type approximation of
analytic functions is exponentially good — while diverging percent-level on
the extended probe box. That is a genuine rival class the data cannot
refute, the sticky-ambiguity rule vetoes certification, and the abstain is
the correct verdict: the approximant-impostor boundary, previously mapped
under declared noise, exists at machine floor too.

**Amended 2026-07-28 after the arbitration validation** (`MUNTZ_ARBITRATION.md`,
Results): these three cells are *two* populations, not one, and the difference
is measurable at the contest — the defeated rival's held-out fraction
h/n = (n − dof)/n.

- `mixed-4term-d2`, `sparse6-d2`: the rival is a 66/68-term **interpolation**
  of the 80 certification points (h/n = 0.01, −0.01) — no held-out evidence at
  all. Arbitration dismisses it; both cells now certify their true laws
  (agreement ~1e-13 out to a wide box). **Reach envelope: 35/36.**
- `rational-d1`: BOTH rivals are dense fractional-power twins that each retain
  ~half the sample as evidence (h/n 0.57 and 0.44), and *neither is the truth*.
  Genuinely indistinguishable on the stated domain; the abstain is correct and
  permanent under this instrument's evidence. Arbitrating it by α margin was
  measured to certify an approximant deviating 2e-5 just outside the box —
  which is why the arbitration rule now also requires every defeated rival to
  be evidence-starved.

Its signature: **draw-conditional incidence** (~10% of smooth-cell draws).
Between the two full audit passes, `cubic-full-d1` flipped abstain→certify
and `mixed-4term-d2` flipped certify→abstain on fresh draws of the same
laws. The instrument's behavior is correct on both sides: when the twin
materializes, refusing to choose is the sound verdict; when it does not,
the certificate stands; zero confident-wrong either way. The road through
this boundary is significance-based arbitration (|H|·q^h accounting across
rival classes) — the program's stated eventual direction, not patched here.

## Standing caps (stated, not silent)

- Sparse sums of size ≥ 6 over the wide basis sit at the proposal-reach
  frontier (top-6-singles quintuples and CAP-Q's size-4 exhaustive are the
  current mechanisms); size-6 recovery is draw-dependent.
- Σx²-composition cells certify but slowly (13–21 min) — the escalation
  cost, not a reach failure.
- **Reach ORDERING** (diagnosed 2026-07-28): certification ends escalation,
  even ambiguous certification. When the dense channel certifies rival
  approximants at tier 1, the engine never reaches the later-tier channel that
  holds the true form — `rational-d1`'s `(2x+1)/(x+3)` is well inside the
  rational channel's reach, and the cell still abstains. A cap on the ORDER of
  the escalation, not on its reach; closing it means letting ambiguity (as
  opposed to certification) escalate, which is a separate registered change.
- Presentation polish: certified laws are sometimes trig-disguised
  equivalents (sin²+cos² foldings of exact laws); `reduce_to_minimal`'s
  trigsimp runs only under declared noise today.
