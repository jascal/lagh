# Escalation rule — registration

**Registered 2026-10-04, before any scored run.** This addresses the open boundary
left by [refusal-guided acquisition](REFUSAL_ACQUISITION_RESULTS.md): on clean
data, the engine stops at the first tier with a non-empty certifying set, so a
later tier's truth never enters the comparison.

## The measured failure (diagnosis, pilot)

**Empirical (pilot, disclosed):** on P2 seed 20, `(3x+2)/(x+4)` sampled at 400 points on
[.5, 3] with `sigma=0`, the first split's **tier 1** certifies a 22-term dense
linear-channel fit. It mixes fractional powers, logs, `exp(-|x|)` and
`sin²`/`cos²`, with 11-digit coefficients, and its maximum error is
4.9×10⁵ on the [.005, 300] probe. The truth `(3/4·x+1/2)/(1+x/4)` certifies on
the same certification rows (0 misses). It is a tier-2 implicit rational, so
the engine never proposes it: escalation stops at tier 1. The other seven
false-exact P2 draws carry the same signature: a dense fractional-power sum
certified at the first tier.

## The rule under test

`discover(..., escalation="pool")`: judge the certifying set of the **highest
tier** (`_tier_candidates` already includes every lower tier's candidates),
instead of returning at the first non-empty tier. Everything downstream is
unchanged: coherence, constrained-input coherence, significance arbitration,
pinning gates, the significance gate and passive re-split stickiness. A
lower-tier truth and a higher-tier approximant now meet in one coherence pass.
They either share a class, and the minimum-complexity member wins, or they form
rival classes and the verdict is a structural refusal. The rule is opt-in;
`escalation="first"` stays the default, and the default path is unchanged.

**Pilot evidence (10 of the 16 P2 seeds, run during diagnosis; disclosed):**

| P2 seeds | `first` | `pool` |
|---|---|---|
| 20, 25, 29, 33, 34 | wrong approximant certified | exact truth certified |
| 22, 23, 24 | wrong approximant certified | structural refusal, truth among rivals |
| 21, 30 | structural refusal | structural refusal, truth among rivals |

`pool` took 142–200 s per draw. These ten seeds are not fresh evidence for the
P2 predictions below. Seeds 26, 27, 28, 31, 32 and 35 are.

## Registered measurements

Runner: `experiments/run_escalation_study.py`. Each of 55 cases runs under both
rules: p2 (16), p1 (3 distinct seed-10 datasets) and reach (the 36 unchanged
`reach/audit.py` cells). Each certified law is scored on 400 fresh in-box
points (**domain claim**) and on 400 points of the extended box [.25, 6]^d
(**form claim**). Exactness is checked symbolically when the truth is a sympy
expression. A certificate is **wrong** if its in-box relative error exceeds
1e-6 (wrong domain), or if its form is wrong: inexact against a symbolic truth,
or extended-box relative error above 1e-6 otherwise. Null calibration:
`run_null_calibration.py --escalation pool`, 200 OS-seeded trials.

## Predictions

- **E1 (P2, 16 seeds).** `pool`: 0 wrong certificates. At least 5 exact-truth
  certificates (the pilot shows 5; on the 6 unseen seeds, no wrong certificate).
  `first` reproduces the recorded 8 wrong certificates (sanity check).
- **E2 (P1, 3 datasets).** `pool`: 0 wrong certificates.
- **E3 (reach, 36 cells).** `pool`: 0 wrong certificates; at least 30/36
  certified. Regressions are expected where a higher-tier approximant becomes a
  rival of a lower-tier truth. Count, not cause, is registered.
- **E3′ (measurement, no prediction).** The `first` rule's wrong-form count on
  reach. The reach audit has only ever recorded whether a cell certified, never
  whether the law was right, so this is new.
- **E4 (null).** `pool`: 0 false certifications in 200 trials. Pooling enlarges
  |H|, and α counts every candidate, so the bound should still hold.
- **E5 (default unchanged).** The full suite stays green. `first` on reach
  reproduces `reach_audit.json`'s 35/36 certified set.
- **E6 (cost, measurement).** Median seconds per case under each rule, per bank.

## Decision rule (registered)

`pool` becomes **eligible** to replace `first` as the default only if all of
the following hold:
- 0 wrong certificates on every bank;
- null 0/200;
- at most 2 reach cells lose certification relative to `first`;
- median reach runtime under `pool` at most 5× `first`.

Even when eligible, switching the default is a separate change, and every
campaign result that depends on the engine is re-scored first. If any
condition fails, `pool` stays opt-in and the failure is reported. The P2
entry guard in `refusal_acquisition` stays in place either way.

**Scope of the rule.** `pool` changes only the main tier loop (tiers 1–5).
Three paths outside that loop are unchanged and still return early:
- the C7 Lévy grammar, on its own domain, before the loop;
- the dim ≥ 3 CAP-S pre-pass, which returns a unique certifying class of
  closed-form monomials/angulars before the loop;
- the C6 exact-integer tier, after the loop.

A first-certifying shortcut therefore survives in the pre-pass. No failure has
been measured there, so this registration does not touch it. It also gives no
guarantee for a truth outside the grammar; that boundary stays **open**.

## Results — `pool` (scored 2026-10-04)

**Empirical: `pool` fails the registered decision rule. It produced a new wrong
certificate.** Artifacts: `experiments/results/escalation/{first,pool}/`,
`summary.json`, `null/`.

| bank | `first` certified | `first` wrong | `pool` certified | `pool` wrong | median s (`first` → `pool`) |
|---|---:|---:|---:|---:|---:|
| P2 (16) | 8 | **8** | 6 (all exact) | 0 | 18.9 → 263.4 |
| P1 (3) | 1 | **1** (rational-a) | 1 (rational-b, exact) | 0 | 19.5 → 260.3 |
| reach (36) | 35 | 0 | 31 | **1** (sparse5-d2) | 12.7 → 269.2 |
| null (200) | — | — | 0 false certs | — | — |

- **E1 holds.** `pool` has 0 wrong and 6 exact certificates on P2: seeds 20, 25,
  27, 29, 33 and 34. Seed 27 was unseen and refused under `first`. `first`
  reproduced its 8 wrong certificates.
- **E2 holds.** `pool` has 0 wrong on P1. rational-a goes from wrong to refused,
  rational-b from refused to exact, and rational-c refuses under both rules.
- **E3 fails.** `pool` certifies 31/36 (≥ 30 holds) but issues one wrong
  certificate. On `sparse5-d2`, `first` certifies the exact truth
  `x0² − x1 + log x1 + sin x0 + 2/x0` at tier 1. `pool` certifies a 29-term
  approximant: in-box relative error 9×10⁻¹⁴, extended-box 6×10⁻⁶.
- **E3′ (measured).** `first` has 0 wrong-form certificates on its 35 reach
  certificates.
- **E4 holds.** 0/200.
- **E5 holds for the reach bank.** `first` reproduces `reach_audit.json`
  exactly: same 35 certified cells, same laws. The full suite has not yet been
  run against this commit.
- **E6.** `pool` is about 21× slower on the median reach case.
- **Decision.** `pool` is not eligible. Three conditions fail: one wrong
  certificate; 4 reach cells lost (sparse3-d1, sparse4-d1, trig-prod-d2,
  trig-sum-d1, all now structural refusals); and cost 21× against a 5× limit.

**Cause, verified.** The rule's premise was false. `_tier_candidates(t)`
includes the lower tiers' *terms*, not their *proposals*. On the
`sparse5-d2` split, tiers 1–3 propose the exact truth among 20, 20 and 26
candidates. Tiers 4 and 5 propose 36 and 73 candidates and none of them is
the truth: the larger library changes which sparse supports the greedy
channels emit. Judging only the top tier therefore discards lower-tier truths,
which is the mirror image of the failure it was meant to fix.

## Amendment A1 — `accumulate` (registered 2026-10-04, before its scored run)

`discover(..., escalation="accumulate")` runs every tier in order, keeps each
tier's own certifying candidates, and judges their **union** once, after the
last tier, using the same downstream machinery. |H| (`total`) already
accumulates across tiers; duplicates are counted again, which overstates α.
That direction is conservative. The `pool` rule and its results are retained
unchanged. Same runner, banks, scorer and null protocol.

Predictions:

- **A1 (P2).** 0 wrong. Exact certificates ≤ 6, `pool`'s count. A larger
  certifying set can add rival classes, for example the tier-1 approximant,
  which `pool` never saw.
- **A2 (P1).** 0 wrong.
- **A3 (reach).** 0 wrong. `sparse5-d2` either certifies the exact truth or
  refuses. Certified count ≥ 30/36.
- **A4 (null).** 0 false certifications in 200 trials.
- **A5 (cost).** Median per bank within 1.3× of `pool`, because the lower tiers
  are cheap next to tier 5. So the 5× default-eligibility condition is
  **expected to fail**. If A1–A4 hold, `accumulate` is the sound opt-in rule,
  and making it fast enough to be the default is a separate, open problem.
