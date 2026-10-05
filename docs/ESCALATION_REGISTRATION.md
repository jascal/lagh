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
