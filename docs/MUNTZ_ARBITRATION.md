# Significance arbitration at the Müntz boundary — registration

**Registered 2026-07-28, before implementation.** The road through the
approximant-impostor boundary named in `REACH_ENVELOPE.md` and
`DIRECTION_SIGNIFICANCE.md`: when materially-different classes certify, the
program's own α currency can sometimes arbitrate — soundly, as a domain
claim.

## The principle

Rival classes A and B both certify (every point within ε); they diverge only
on the extended probe — i.e., OUTSIDE the certificate's stated domain. The
per-class chance-fit bounds α ≤ |H|·q^h share |H| and q but differ in
h = n − dof: a high-dof approximant burns its evidence on its own
parameters. The α ratio between class representatives is q^(Δdof) — for the
measured Müntz case (q ≈ 6×10⁻¹³, Δdof ≈ 12) the whale's bound is ~146
orders of magnitude weaker; for the historical Gaia loose-floor case, ~50.

**Arbitration rule (registered):** among certifying classes, if exactly one
class's representative α_log10 is smaller than EVERY rival's by at least
**MARGIN = 30 orders of magnitude**, that class wins; the winner then passes
every existing gate (pinning, floor-dominated winner gate, significance
gate) unchanged, and the certificate carries an arbitration note naming the
margin and the defeated rivals' bounds. If no class clears the margin over
all rivals, the structural abstain stands exactly as today.

## Why this is sound under zero-confident-wrong

The certificate is a DOMAIN claim; rivals agree on the domain by
construction (they all certify). The risk arbitration must not take is
form-overclaim beyond the domain — hence the astronomical margin: a rival
whose chance-fit bound is ≥ 10³⁰ times weaker is, within the program's own
accounting, overwhelmingly an artifact of its richer family. Conservative
bias notes: dof counting excludes ±1/0 coefficients, which favors simple
integer laws — the direction that protects the truth in every measured
case; marginal contests (Δdof small, q moderate) fail the margin and
abstain, exactly as today.

## Registered predictions

- **P1.** The three boundary cells in the current reach matrix
  (`rational-d1`, `mixed-4term-d2`, `sparse6-d2`) certify their true laws
  after arbitration; the full 36-cell audit reaches **36/36** with no cell
  regressing.
- **P2.** The full suite stays green; NO existing certification changes its
  law (arbitration only converts abstains, never alters winners).
- **P3.** Null validation stays **0 false certifications / 200
  true-random targets** (arbitration operates only on certifying classes;
  nulls produce none).
- **P4.** The historical Gaia loose-floor scenario (floor 2e-4, whale vs
  truth classes) resolves to the truth under arbitration — consistent with
  the eventually-amended answer, now reachable without the amendment.
- **P5.** Zero confident-wrong throughout.

## Results

Implemented 2026-07-28 (`certify.arbitrate_significance`, applied in
`engine.discover` after the constraint branch). Validation **falsified the
registered rule as stated** and the rule was amended the same day, before any
campaign consumed it. Everything below is measured.

### The counterexample that amended the rule

Every contest arbitration actually meets, instrumented at the call site
(σ = 0, tier 1, n = 80 certification points, per-cell audit seeds):

| cell | winner dof, h/n | defeated rival dof, h/n | α margin | registered rule | amended rule |
|---|---|---|---|---|---|
| `mixed-4term-d2` | 2, **0.97** | 79, **0.01** | 960 | certify | certify — and the winner IS the truth |
| `sparse6-d2` | 1, **0.99** | 81, **−0.01** | 971 | certify | certify — and the winner IS the truth |
| `rational-d1` | 34, **0.57** | 45, **0.44** | 125 | certify an **approximant** | abstain |

(h = n − dof: the points left over after the form's own free parameters are
pinned — the evidence the certificate actually rests on.)

`rational-d1` is the failure. At tier 1 the dense channel proposes two
fractional-power twins — *neither of them the truth* `(2x+1)/(x+3)`, which
lives in the rational channel the escalation never reaches, because something
certified first and the engine stops escalating on ambiguity, not on reach.
The margin rule then crowns the cheaper twin at 125 orders. On the audit's own
draw the passive full-data gate happened to catch it; on a seed-0 draw of the
same law it certified through: a law agreeing with the truth to 1×10⁻¹² inside
the sampled box, **1.9×10⁻⁵ at 1.5× the box and 7×10⁻² on a wide box**. That is
a form-overclaim — the one risk the registration named, arriving by a route the
registration did not anticipate: *all* rivals dense, the truth absent from the
contest entirely, so nothing in the α ranking has any purchase on structure.

### The amendment

Arbitration now requires **both** conditions:

1. **Margin** — one class's α beats every rival's by ≥ 30 orders (as
   registered), and
2. **Rival evidence** — every defeated rival is *evidence-starved*:
   h/n < `ARBITRATION_RIVAL_EVIDENCE_MAX` = **0.10**. It must be an
   interpolation of the certification sample, not a law the sample constrains.

The measurement picks the bar: dismissible rivals sit at h/n = 0.01 and −0.01
(66- and 68-term fits of 80 points — zero held-out evidence, the H1a mechanism
`_significance_gate` was built for), genuine rivals at 0.44. Any bar in
(0.02, 0.44) separates them; 0.10 keeps a 4× margin on the side that matters.

Why this is the sound reading, and the margin alone was not: α bounds **chance**
agreement, never structure — `significance_log10`'s own docstring says so
("structure uniqueness is coherence's job, and the approximant-impostor
boundary is unaffected"). What α licenses structurally is *dismissal of a rival
the data never constrained*: a 66-term fit of 80 points is a restatement of the
sample, not a competing law. It does not license ranking two rivals that each
hold out half the sample — that is structure inference, and α cannot do it.
Note that the registration's own soundness argument ("a rival whose bound is
≥ 10³⁰ weaker is an artifact of its richer family") silently assumed the truth
was *in* the contest; `rational-d1` is the case where it is not.

### Predictions scored

- **P1 — MISSED as stated.** Not 36/36: the audit re-run under the amended
  rule gives **35/36**. `mixed-4term-d2` and `sparse6-d2` convert
  abstain→certify and their certified laws agree with the truths to ~1×10⁻¹³
  out to a wide box (they are trig-folded restatements). `rational-d1` still
  abstains, and now with a diagnosis: its rivals are genuine, so refusing to
  choose is the correct verdict — an abstain no arbitration rule should ever
  convert. No cell regressed.
- **P2 — met.** Full suite green (**88 tests**, including the amended
  arbitration set and the new weak-form set); no existing certification changed
  its law; the audit shows **no cell regressed** against the committed
  baseline, only the two gains. The Gaia C0 photometric certificate is
  byte-identical at floor 5e-6 (`12660019/492850 − 5·x₀/2`, α ≤ 10⁻⁴⁶¹·¹⁴⁵).
- **P3 — MET, measured 2026-07-29.** `experiments/run_null_calibration.py --out
  experiments/results/null_calibration_amended.jsonl`, 200 fresh OS-seeded
  trials under the amended rule: **0 / 200 false certifications**, and no trial
  emitted a law at all. The argument this replaces — that the amendment can only
  make certification stricter, so a pre-amendment 0/200 implies an amended one —
  was sound, and is now a measurement rather than an argument. The
  pre-amendment run is kept beside it (`null_calibration.jsonl`) rather than
  overwritten, which is what the new `--out` flag is for.

  One observation that is NOT a controlled comparison and should not be read as
  one: the amended run took 117 min against the original's 190 (mean 35.1 s vs
  57.1 s per trial). Several engine changes landed between the two — the
  coherence early exit chief among them — and the machine was loaded
  differently. It is consistent with the engine having got faster on
  no-certifying-class inputs; it does not measure that.
- **P4 — NOT APPLICABLE as registered.** Measured on the frozen C0 snapshot at
  floor 2e-4: **parametric abstain**, and *no arbitration contest ever forms* —
  the whale/truth rivalry the prediction assumed was already dissolved upstream
  by refit-parsimony collapse (the loose-ε closure, C0 issue). Arbitration is a
  no-op on that scenario; the honest verdict there remains "the slope is not
  pinned at a loose floor". The prediction described an engine state that no
  longer exists.
- **P5 — MISSED for the registered implementation, restored by the
  amendment.** The registered rule *did* produce a confident-wrong (the
  `rational-d1` seed-0 draw above). It was caught by this validation, before
  any campaign or benchmark consumed it, and the seed-0 draw now abstains. The
  program's zero-confident-wrong record holds for shipped results; the
  registration-then-validation discipline is what kept it.

### What this changes about the boundary map

The three abstaining cells were described in `REACH_ENVELOPE.md` as one
population ("Müntz twins"). They are **two**, and only one is arbitrable:

- **Interpolation rivals** (dof ≈ n, h/n ≈ 0): the sample cannot see them as
  laws at all. Arbitration dismisses them; `mixed-4term-d2` and `sparse6-d2`.
- **Constrained twins** (dof ≈ n/2, both significant, both with real held-out
  evidence): genuinely indistinguishable on the stated domain. The abstain is
  correct and is **permanent under this instrument's evidence** — only more
  data, a wider box, or an oracle query outside the box can break the tie.
  `rational-d1`.

A second, separate finding falls out of the same diagnosis: **dense-channel
certification at tier 1 pre-empts escalation to the channel that holds the
truth.** `rational-d1`'s law is trivially within reach of the rational channel,
but the engine never gets there — certification (even ambiguous certification)
ends escalation. That is a reach-ordering cap, not a reach failure, and it is
now stated in `REACH_ENVELOPE.md` rather than lurking.
