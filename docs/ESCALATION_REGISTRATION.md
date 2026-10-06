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
  exactly: same 35 certified cells, same laws. With `accumulate` in place, the
  full suite passes: 36 files, 423 tests, one process per file. The 3 new
  tests in `tests/test_escalation.py` also pass.
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

## Results — `accumulate` (scored 2026-10-04)

**Empirical: `accumulate` issued zero wrong certificates on every bank and kept
null 0/200. It is not eligible as the default: it loses 5 reach cells and is
about 31× slower.**

| bank | `first` cert / wrong | `pool` cert / wrong | `accumulate` cert / wrong | median s (`accumulate`) |
|---|---:|---:|---:|---:|
| P2 (16) | 8 / **8** | 6 / 0 | 5 / 0 (all exact) | 289.8 |
| P1 (3) | 1 / **1** | 1 / 0 | 0 / 0 | 239.1 |
| reach (36) | 35 / 0 | 31 / **1** | 30 / 0 | 393.9 |
| null (200) | — | 0 | 0 | — |

- **A1 holds.** 0 wrong. 5 exact (seeds 20, 25, 29, 33, 34), which is ≤ 6.
  Seed 27 drops from `pool`'s exact certificate to a refusal: the larger union
  adds a rival class.
- **A2 holds.** 0 wrong. All three P1 datasets refuse; rational-b loses
  `pool`'s exact certificate the same way.
- **A3 holds.** 0 wrong. `sparse5-d2` certifies the exact truth
  `x0² − x1 + log x1 + sin x0 + 2/x0`, the same law `first` returns. 30/36
  certified, exactly the registered floor.
- **A4 holds.** 0/200.
- **A5 fails for reach.** Median cost relative to `pool`: P2 1.10× and P1 0.92×
  (both within 1.3×), but reach 1.46×. Against `first` on reach, `accumulate`
  is 31× slower.
- **Decision.** Not eligible. Lost reach cells: linear-7term-d6, sparse3-d1,
  sparse4-d1, trig-prod-d2, trig-sum-d1 (limit 2). Cost is 31× against a 5×
  limit.

**What the losses are.** Every lost cell now refuses with the true law among
its rivals. In four cells the other rival is a dense fractional-power
approximant with 11-digit coefficients, the same signature as P2's false
certificates. linear-7term-d6 is instead paired with a `sqrt((…)²)` twin of the
linear truth. Each is a structural refusal, which is sound but costs reach.
The approximants keep held-out evidence (h/n ≥ 0.10), so registered
significance arbitration (`MUNTZ_ARBITRATION.md`) correctly declines to
dismiss them.

## Status after this registration

**Empirical:** clean-data false exactness is caused by the first-certifying-tier
stop. Judging the union of every tier's certifying candidates removes it on
all measured banks (9 → 0 wrong on P2/P1, 0 on reach, 0/200 null). It is
available opt-in as `escalation="accumulate"`; the default stays `"first"`.

**Open:**
1. A sound way to dismiss the dense approximants that `accumulate` meets.
   They are evidence-bearing fits of the wrong form, and |H|·q^h does not
   separate them from the truth, so this needs a structural criterion, not a
   significance margin.
2. Cost. Every tier always runs, and tier 5 dominates.
3. The dim ≥ 3 pre-pass shortcut remains untested.

The P2 entry guard in `refusal_acquisition` stays in place.

## Amendment A2 — the joint coefficient gate (registered 2026-10-05, before its scored run)

**Diagnosis.** Every approximant that blocks or fools escalation is a dense,
nearly collinear sum whose coefficients are huge-denominator rationals. These
are gated atoms, but each passes `float_pinned`, which perturbs one atom at a
time by ±10⁻⁵ and ±10⁻⁴ relative. At machine floor, any single-coordinate move
breaks the fit. That does not show the coefficient **vector** is determined. A
near-collinear support has a direction along which all coefficients move
together while every prediction stays in the band. This is the same
marginal-vs-joint distinction the state certificates already had to learn.

**Rule.** `certify.joint_pinned` keeps the marginal gate's perturbation sizes
(10⁻⁵ and 10⁻⁴ relative, the maximum over atoms). It takes them along the
least-determined direction instead: the smallest right singular vector of the
Jacobian of the predictions with respect to each gated atom's relative change,
scaled by the band. The moved law goes through the real `check`. If it still
certifies, the candidate is rejected exactly as a `float_pinned` failure is.
There is no new tolerance. The gate only removes candidates; candidates with
fewer than two gated atoms are left to the marginal gate. It applies on the
same path as the marginal gate (clean data, not floor-dominated, per
candidate). Opt-in: `discover(..., coefficient_gate="joint")`, default
`"marginal"`.

**Pilot (disclosed; not evidence for the predictions below):**
- P2 seed 20, first split, `first` + joint gate, `max_tier=2`: the tier-1
  approximant is rejected and escalation certifies the exact truth at tier 2
  in 12 s.
- reach `trig-sum-d1`, `accumulate` + joint gate: certifies the exact
  `2 sin x0 + cos x0` with no rivals.

**Arms.** `first+joint` and `accumulate+joint` on the same 55 cases, with the
same scorer, plus 200 null trials for each arm.

- **J1 (P2).** `first+joint`: 0 wrong; at least 5 exact certificates.
- **J2 (P1).** `first+joint`: 0 wrong.
- **J3 (reach).** `first+joint`: 0 wrong; at least 34/36 certified.
- **J4.** `accumulate+joint`: 0 wrong on every bank; reach at least 33/36. The
  four approximant-blocked cells may recover. linear-7term-d6's `sqrt((…)²)`
  twin has no gated atoms and is expected to keep refusing.
- **J5 (null).** 0/200 for both arms.
- **J6 (cost).** `first+joint` median at most 2× `first` on every bank.
- **J7 (measurement).** The full suite with the default temporarily set to
  `"joint"`, to find which existing tests encode certificates that the gate
  would remove. The default itself is not changed here.

**Decision rule.** The joint gate is eligible to become the default (with
`first`) only if J1–J3, J5 and J6 hold and J7 shows no failures. Even then the
default changes in a separate change, after the campaign certificates in
`CERTIFICATES.md` with two or more gated atoms are re-scored. Under that
outcome, `accumulate` would remain the slower, more conservative opt-in.

## Results — joint coefficient gate (scored 2026-10-05)

**Empirical: the joint gate issued zero wrong certificates under both rules and
kept null 0/200. Under the default `first` rule it certifies the exact truth on
every P2 and P1 witness and 36/36 reach cells. As registered, it is not
default-eligible: J6 (cost) fails, and J7 failed on a crash in the gate, now
fixed.**

| bank | `first` | `first+joint` | `accumulate` | `accumulate+joint` |
|---|---|---|---|---|
| P2 (16): cert / wrong / exact | 8 / **8** / 0 | **16 / 0 / 16** | 5 / 0 / 5 | 16 / 0 / 16 |
| P1 (3) | 1 / **1** / 0 | **3 / 0 / 3** | 0 / 0 / 0 | 3 / 0 / 3 |
| reach (36): cert / wrong | 35 / 0 | **36 / 0** | 30 / 0 | 34 / 0 |
| null (200) false certs | — | 0 | 0 | 0 |
| median s, P2 / P1 / reach | 18.9 / 19.5 / 12.7 | 68.2 / 78.1 / 21.2 | 289.8 / 239.1 / 393.9 | 373.0 / 340.2 / 287.4 |

- **J1 holds.** All 16 P2 draws certify the exact truth at tier 2. That includes
  the 8 former false certificates and the 8 former refusals.
- **J2 holds.** All three P1 datasets certify their exact truths at tier 2.
- **J3 holds.** 36/36, 0 wrong. `rational-d1`, the cell that kept the passive
  reach audit at 35/36, now certifies `(2x/3+1/3)/(x/3+1)` = `(2x+1)/(x+3)`
  exactly, with extended-box error 3×10⁻¹⁶. This is a new rule, so the
  historical 35/36 artifact is not amended.
- **J4 holds.** `accumulate+joint` has 0 wrong and reaches 34/36. sparse3-d1,
  sparse4-d1, trig-sum-d1 and rational-d1 recover. trig-prod-d2 and
  linear-7term-d6 keep refusing. Each pairs the truth with a `sqrt(f²)` twin,
  which has no gated atoms. The twin equals the truth wherever f > 0, and the
  probe leaves that region.
- **J5 holds.** 0/200 for both arms. Null median per trial: 46.2 s
  (`first+joint`), 29.2 s (`accumulate+joint`), 31.0 s (historical baseline).
- **J6 fails.** `first+joint` over `first`: P2 3.6×, P1 4.0×, reach 1.67×,
  against 2× on every bank. On P2/P1 most of the increase is the tier-2
  escalation the correct answer needs: `first` was fast because it stopped at a
  wrong tier-1 certificate. The registered criterion compared against that
  cost, but it stands as written.
- **J7 failed as registered.** The suite run with `"joint"` as the default had
  2 failures out of 426:
  - `test_escalation::test_first_rule_stops_at_the_tier1_approximant`, which
    documents the default rule's recorded failure, so it is expected to fail
    under a default that removes it;
  - `test_pdesystem::test_both_equations_certify_the_true_support_over_shared_rows`,
    which was a **crash in `joint_pinned`**. It assumed an array band; weak-form
    rows pass a per-candidate callable `PatchEpsilon`. Fixed by resolving the
    band through `certify.band(eps, expr)`, as the other gates do. With the fix
    and the joint default, `test_pdesystem.py` passes 7/7. For array bands the
    fix is a no-op, so the scored arms above are unaffected.
- **Decision.** As registered, not eligible: J6 and J7 failed. The gate stays
  opt-in and the default is unchanged.

## Status after A2

**Empirical:** the joint gate closes clean-data false exactness on every
measured bank, under the default escalation rule. It turns each wrong
certificate into the exact truth instead of a refusal, and adds one reach cell.
Its cost is the escalation that the correct answer requires.

**Open:**
1. Whether to adopt it as the default. That needs a new registration with a
   cost criterion that doesn't credit the old rule for stopping at wrong
   answers, a clean J7 rerun, and re-scoring of the campaign certificates
   with two or more gated atoms (`CERTIFICATES.md`).
2. `sqrt(f²)` twins (trig-prod-d2, linear-7term-d6) under `accumulate`.
3. Noise. The gate is registered only for the clean-data per-candidate path;
   the declared-noise and floor-dominated winner gates are unchanged.
4. The dim ≥ 3 pre-pass.

## Amendment A3 — switching the default to the joint gate (registered 2026-10-05, before any scored run)

**Change under test.** Make `coefficient_gate="joint"` the default of
`discover` and `discover_passive`. Escalation stays `"first"`. Nothing else
changes.

**Why a new registration.** A2 left the gate opt-in on two registered failures:
- **J6's cost criterion** compared against the old rule's time, which was
  short because it stopped at wrong tier-1 certificates. Rewriting that
  criterion after seeing the J6 data would be moving the goalposts, so the new
  criterion below is evaluated **only on fresh banks** that no earlier run has
  seen.
- **J7** failed on a crash in the gate that has since been fixed. It is rerun
  here.

### Banks

1. **Fresh rational (`fr`, 24 draws).** `(a·x+b)/(x+c)` on [.5, 3], 400 rows,
   seeds 100–123. a, b, c are integers 1..5 from `rng(seed)`, redrawn
   deterministically when a·c = b. Runner: `experiments/run_default_switch.py`.
2. **Fresh reach (`frch`, 36 cells).** Every `reach/audit.py` cell re-drawn with
   seed `crc32(name)+1`. Same runner.
   Both fresh banks run both gates and use the A2 scorer: fresh in-box points
   for the domain claim, the extended box [.25, 6]^d for the form claim.
3. **Campaigns (14 scripts).** Every real-data campaign under `experiments/`
   (gaia c0, p1, p2, p3; solar; macro; exoplanet c0, c1, c2, c5, ph2;
   materials c0, c1, c2). This covers every certificate in `CERTIFICATES.md`.
   Driver: `experiments/run_default_switch_campaigns.py`. Each arm runs in its
   own clean worktree of the registered commit; the joint arm differs only by
   the one-line default flip in `engine.py` and `passive.py`. Both arms run with
   4 jobs, and the campaigns run alone, not alongside the fresh banks, because
   `recover` has a 45 s wall-clock budget. Compared fields: every
   `certified`, `law`, `alpha_log10` and `abstain` in every result file.
   **Fresh marginal vs fresh joint, not the committed artifacts**: the engine
   has changed since July, so marginal-vs-committed drift is reported
   separately and not charged to the gate.
4. **Null.** 200 new OS-seeded trials with `--coefficient-gate joint`.
5. **Suite.** The full suite in a worktree with the default flipped (J7 rerun).

### Predictions

- **D1 (fresh rational).** Joint: 0 wrong; at least 20/24 exact.
  Marginal: measured, no prediction.
- **D2 (fresh reach).** Joint: 0 wrong; no cell that marginal certifies
  correctly is lost.
- **D3 (campaigns).** Every certificate the fresh marginal arm issues is
  issued by the joint arm with the same law. Any new joint certificate is
  reported and inspected; a new certificate that a campaign document records
  as an abstain, conjecture or wrong counts as a failure.
- **D4 (null).** 0/200.
- **D5 (suite).** With the joint default, exactly one test fails:
  `test_escalation::test_first_rule_stops_at_the_tier1_approximant`, which
  documents the marginal default's recorded failure.
- **D6 (cost, fair).** On each fresh bank, over the cases where marginal did
  **not** issue a wrong certificate, the median per-case ratio of joint time
  to marginal time is ≤ 2×. On campaigns, every script's wall time under joint
  is ≤ 3× marginal, and no D3 difference is caused by a time-budget
  truncation.

### Decision rule

If D1–D6 all hold, open the switch PR:
- flip both defaults;
- re-point the documenting test at `coefficient_gate="marginal"`;
- add a ledger row to README;
- add a note to `CERTIFICATES.md`.

If any prediction fails, the default stays `marginal` and the failure is
recorded. A D3 failure is decisive on its own: the gate may not cost a
certificate the campaigns stand on.

**Protocol deviation (2026-10-05, before any result was compared).** The first
campaign run was discarded unscored: the machine suspended three times during
it, at 09:43–09:51, 10:24–10:37 and 10:58–12:15. Each was a logind suspend on
AC power with no lid event. A sleep stretches wall times and can exhaust
`recover`'s wall-clock budget, so both arms were contaminated. Only exit codes
and times had been seen; no law or certificate comparison had been made.

The rerun:
- runs under `systemd-inhibit --what=sleep:idle`;
- records each run's `suspended_seconds`, the drift between `CLOCK_BOOTTIME`
  and `CLOCK_MONOTONIC`;
- treats any run with non-zero `suspended_seconds` as void, and reruns it.

An unrelated job from another session (`pic`) was loading the CPU during the
run. It is outside this study's control and is noted, not corrected for.

## Results — A3 default switch (scored 2026-10-05)

**Empirical: the switch is rejected under the registered rule. D3 fails: the
joint gate costs the Gaia P3 frame-rotation certificate (α ≤ 10⁻⁷⁶⁸⁰). The
default stays `marginal`.** Every run used below recorded 0 s suspended.
Void runs are kept under `experiments/results/default_switch/void/`, and both
arms of each void run were rerun side by side.

| prediction | result |
|---|---|
| **D1** fresh rational (24) | **holds.** joint: 24 certified, all 24 exact, 0 wrong. marginal: 14 certified, 5 exact, **9 wrong** |
| **D2** fresh reach (36) | **holds.** joint: 36 certified, 0 wrong, no correct certificate lost. marginal: 33 certified, **1 wrong** (`rational-d1` on fresh data) |
| **D3** campaigns (14 scripts, 13 result files) | **fails.** 19 certificates under marginal, 18 under joint. The one difference: `gaia_p3/P1_frame_rotation` is certified under marginal and structurally refused under joint. Every other certificate is identical, law and α included |
| **D4** null | **holds.** 0/200 |
| **D5** suite with joint default | **holds.** 428 passed; the only failure is the predicted `test_first_rule_stops_at_the_tier1_approximant` |
| **D6** cost | **fails.** Fair median ratio: fresh rational 2.26× (limit 2×), fresh reach 1.00×. Campaigns: gaia p3 3.35× (limit 3×), and it is the D3 case, searching every tier before refusing. The other 13 scripts ran at 0.49–1.03× |

macro's result file was unchanged from the committed artifact in both arms,
so it is identical between them. Campaign wall times are confounded by load:
another session's CPU job ran during the marginal arm and had stopped before
most of the joint arm, so ratios below 1 are load, not speed. The fresh banks
had all finished before the machine's power profile was changed at 14:30:43.

**Cause of the D3 failure, verified.** The Gaia inputs are direction cosines,
which satisfy x0² + x1² + x2² = 1 exactly. At tier 1 the engine never
proposes the bare three-term linear law. It proposes three certifying
candidates that are the truth plus multiples of the constraint, for example
`−0.2169·x0²·x1 − 0.2169·x1³ − 0.2169·x1·x2² + …`, which equals
`−0.2169·x1` on the sphere. Under marginal, the constrained-input path
certifies one of them and reduces it modulo the constraint. The constraint
ideal is an **exact** joint flat direction, so the joint gate correctly finds
each candidate's coefficients undetermined and rejects all three. Nothing
certifies through tier 7. The full-data three-term law itself passes the joint
gate by a factor of about 3×10⁷.

**What the fresh banks add, empirical.** The shipped default issues wrong
certificates on 9 of 24 new clean rational draws and on fresh-seed
`rational-d1`. Clean-data false exactness is not confined to the P2/P1
witnesses; it is common on this family. Under the joint gate those cases
certify the exact truth.

**Open, and the next registration (A4).** Apply the joint test modulo
machine-exact input constraints: reduce each candidate through
`reduce_mod_constraints` when `input_constraints` finds any, and test joint
pinning on the reduced law, so the constraint ideal does not count as
non-identification. It must be registered with the frame rotation as a named
cell, and the cost criterion will have to deal with the 2.26× fresh-rational
ratio. That ratio is the price of escalating to the tier where the right
answer is.

## Amendment A4 — the joint gate modulo input constraints (registered 2026-10-05, before any A4 run)

**Rule.** `coefficient_gate="joint_modulo"`. When `input_constraints` finds
machine-exact polynomial constraints on the inputs (computed once per
`discover`, on all rows), each candidate is reduced with
`reduce_mod_constraints`, float dust is swept with `reduce_to_minimal`, the
same canonicalization the constrained-input path already applies to winners,
and `joint_pinned` tests that reduced law. The candidate itself enters the
certifying set unchanged, so everything downstream is as before. Without
constraints the rule is exactly A2's `joint`.

The certificate on constrained inputs is already a domain-restricted claim,
made modulo the constraint ideal. Identification is asked modulo the same
ideal, so the ideal is not counted as an undetermined direction.

**Pilot (disclosed):**
- Gaia P3 frame rotation under `joint_modulo`: certifies the same law at
  tier 1 in 10 s (α 10⁻⁷⁶⁴⁸ with the pilot's floor, which omits the campaign's
  ulp term).
- P2 seed 20 (`max_tier=2`): exact truth at tier 2.

### Banks (all new or re-run; nothing from A3 is re-scored)

1. **`fr2`**: 24 new rational draws, seeds 200–223, same generator as `fr`.
2. **`frch2`**: the 36 reach cells on seed `crc32(name)+2`.
   Both run marginal and `joint_modulo`, interleaved case by case so the two
   gates share load.
3. **Campaigns**: all 14 scripts. The marginal and `joint_modulo` worktrees
   run **concurrently**, each script beside its twin, 4 jobs per arm. Output:
   `experiments/results/default_switch_a4/campaigns/`.
4. **Null**: 200 new trials, `--coefficient-gate joint_modulo`.
5. **Suite**: the full suite with `joint_modulo` as the default.

Every run sits under `systemd-inhibit` and records `suspended_seconds`. A run
that slept is void, and both arms of it are rerun together.

### Predictions

- **E1 (`fr2`).** `joint_modulo`: 0 wrong, at least 20/24 exact.
- **E2 (`frch2`).** `joint_modulo`: 0 wrong; no correct marginal certificate
  lost.
- **E3 (campaigns).** Every marginal certificate is issued under
  `joint_modulo` with the same law and the same α. This names the **Gaia P3
  frame rotation** explicitly. Any new certificate is reported and inspected.
- **E4 (null).** 0/200.
- **E5 (suite).** Only
  `test_escalation::test_first_rule_stops_at_the_tier1_approximant` fails.
- **E6 (cost).** A3's single fair ratio mixed two different things, so it is
  split on fresh data:
  - **same verdict** (both gates certify the same law, or both refuse): median
    `joint_modulo`/marginal ratio ≤ 2×. This is the gate's own overhead;
  - **changed verdict** (marginal refuses or is wrong, `joint_modulo`
    certifies): no ratio limit, since this is the cost of reaching the tier
    where the right answer is. Absolute caps instead: median ≤ 120 s, and no
    fresh case over 600 s under `joint_modulo`;
  - **campaigns:** every script ≤ 3× its concurrent marginal twin.

### Decision rule

If E1–E6 all hold, open the switch PR:
- default `coefficient_gate="joint_modulo"`;
- point the documenting test at `"marginal"`;
- add a README ledger row;
- add a `CERTIFICATES.md` note.

Otherwise the default stays `marginal` and the failure is recorded. An E3
failure is decisive on its own.

## Results — A4 (scored 2026-10-05)

**Empirical: every soundness and campaign prediction holds; the switch is
rejected under the registered rule on cost alone. E6's same-verdict ratio on
`fr2` is 2.10×, against a 2× limit.** No A4 run recorded any suspended time.

| prediction | result |
|---|---|
| **E1** `fr2` (24 new rationals) | **holds.** `joint_modulo`: 24 certified, all 24 exact, 0 wrong. marginal: 13 certified, 5 exact, **8 wrong** |
| **E2** `frch2` (36) | **holds.** `joint_modulo`: 36/36, 0 wrong, none lost. marginal: 34, **1 wrong** (`rational-d1` again) |
| **E3** campaigns | **holds.** 19 certificates in each arm and 0 differences across all 13 result files: every law and α identical, Gaia P3 frame rotation included |
| **E4** null | **holds.** 0/200 (median 32.7 s per trial) |
| **E5** suite with `joint_modulo` default | **holds.** 428 passed; the only failure is the predicted documenting test |
| **E6** cost | **fails.** Same verdict: `fr2` **2.10×** over 5 cases, `frch2` 1.00× over 33. Changed verdict: median 30.0 s and max 65.1 s (`fr2`), max 402.4 s (`frch2`), within the caps. Campaigns: worst script 1.03×, with each arm run concurrently beside its twin |

**Where the 2.10× comes from.** The five same-verdict `fr2` cases all certify
the same exact rational at tier 2 under both gates. Two show no overhead (1.01×,
0.99×). Three roughly double: 7.2 → 15.6 s, 11.7 → 24.6 s and 8.9 → 19.1 s.
The truth has no gated atoms, so it returns immediately; the added seconds
are `joint_pinned` on rival tier-2 candidates that carry gated atoms, across
three passive re-splits. That is a re-lambdify and up to four `check` calls
per such candidate. It is the gate's own overhead, and it never changes the
verdict on these cases.

**Status.** The default stays `marginal`. `joint_modulo` is available opt-in
and has the strongest record of any rule measured here:
- 0 wrong certificates on every bank;
- 24/24 and 24/24 exact on two independent fresh rational banks (marginal is
  wrong on 9/24 and 8/24);
- 36/36 reach on both fresh draws;
- every campaign certificate unchanged.

**Open:** the overhead of `joint_pinned` on cheap cases, which is
engineering, not soundness.

## Diagnosis after A4 (2026-10-05): the E6 overhead is not the gate

A performance-only A5 was planned, on the reading that the 2.10× was
`joint_pinned` overhead. Profiling **falsified that reading before any A5
change was made.**

- On `fr2-seed205`, `joint_pinned` ran 9 times for 0.2 s in total. The time
  is in roughly 2000 ordinary `check`/`lambdify` calls.
- Per split, without passive aggregation:
  - seed 205, marginal: splits 0 and 2 certify at **tier 1** (2.5 s and
    2.2 s); split 1 at tier 2.
  - seed 205, `joint_modulo`: all three splits certify at tier 2 (7.2 s,
    3.0 s, 6.6 s).
  - seed 207: the same pattern on two splits (marginal 3.8 s and 5.3 s at
    tier 1; `joint_modulo` 9.3 s and 11.8 s at tier 2).
  - seed 201: tier 2 on every split under both rules, and no overhead
    (2.4–2.6 s against 2.5–3.7 s).
- The marginal tier-1 split certificates are approximants. The passive
  full-data gate rejects them, so the final verdict is the same tier-2 truth.

**Empirical:** E6's "same verdict" bucket classified cases by **final**
verdict and so included cases where marginal certified a wrong law on
individual splits. The measured 2.10× is the cost of escalating those splits
to the correct tier, the same cost E6 meant to exempt under "changed
verdict", not gate overhead. A performance-only change to `joint_pinned`
cannot address it; the planned A5 is withdrawn. A criterion that classifies
by per-split verdict would be a redefinition made after seeing data, so it
could only be judged on a further fresh bank.

## Owner decision (2026-10-05)

The repository owner accepts E6's cost on the A4 evidence and authorizes
switching the default to `coefficient_gate="joint_modulo"`. This is
**recorded as the owner's judgment, not as a registered pass**: under the
registered rule A4 fails E6. The judgment rests on:
- A4's soundness and campaign results: 0 wrong on every bank; every campaign
  certificate identical; null 0/200; suite clean;
- the diagnosis above, that the extra time is escalation to the correct tier
  on splits where marginal certified a wrong law.

A5 below is run as independent measurement of that diagnosis, and its result
is reported in the switch PR whichever way it comes out.

## Amendment A5 — per-split cost classification (registered 2026-10-05, before any A5 run)

**Change.** The cost criterion only; the gate code is A4's `joint_modulo`,
unchanged. A fresh case counts as **changed** if the final verdicts differ
**or** marginal certified a wrong law on any of the three passive splits.
Per-split verdicts are recomputed outside the timed call, using passive's own
split procedure, and each split's certified law is scored with the A2 scorer.
Every other case counts as **same**. The A4 campaign (E3, campaign cost),
null (E4) and suite (E5) results stand for A5 because the code is identical.

**Banks.** `fr3`: 24 rationals, seeds 300–323, same generator. `frch3`: the 36
reach cells on seed `crc32(name)+3`. Both gates, interleaved case by case,
under `systemd-inhibit`, with suspend-voided runs rerun in both arms.

**Predictions.**
- **F1.** `joint_modulo`: 0 wrong on both banks; at least 20/24 exact on `fr3`;
  no correct marginal certificate lost on `frch3`.
- **F2 (cost).** Same cases: median `joint_modulo`/marginal ≤ 2× on each bank.
  Changed cases: median ≤ 120 s, and no case over 600 s.
- **F3 (diagnosis).** On `fr3`, every case whose ratio exceeds 1.5× has a
  marginal split that certified a wrong law.

## Results — A5 (scored 2026-10-05)

**Empirical: F1 and F2 hold; F3 fails on 1 of 21 cases.** None of the 120
runs recorded suspended time.

| prediction | result |
|---|---|
| **F1** | **holds.** `fr3`: `joint_modulo` 24 certified, all 24 exact, 0 wrong. Marginal: 15 certified, 8 exact, **7 wrong**, and a wrong law certified on at least one split in 21/24 cases. `frch3`: `joint_modulo` 36/36, 0 wrong. Marginal: 34, **1 wrong** (`rational-d1`, for the fourth fresh draw running). No correct marginal certificate lost |
| **F2** | **holds.** Same cases: `fr3` 0.99× (2 cases), `frch3` 1.00× (33). Changed cases: median 35.0 s (`fr3`) and 12.8 s (`frch3`); maximum 93.5 s and 423.4 s, against caps of 120 s and 600 s |
| **F3** | **fails, 20 of 21.** 21 `fr3` cases exceed 1.5×. 20 have a marginal split that certified a wrong law. The exception, `fr3-seed322` (2.18×), never certified a wrong law: marginal splits 0 and 2 stopped at **tier 1 with a structural refusal** among dense fractional-power rivals, and `joint_modulo` removes those rivals and certifies the exact truth at tier 2 on every split |

**Reading.** The registered diagnosis named the mechanism too narrowly. The
measured cost is escalation past tier 1 whenever tier 1's certifying set was
made of approximants: they either certify wrongly or force a refusal. On all
60 A5 cases with an unchanged outcome, the gate costs 0.99–1.00×.

**Status.** The default switch rests on the owner decision recorded above,
which A5 is consistent with. The record shows E6 failing under A4's
classification and F3 failing as worded; neither is rewritten.

## Review corrections (2026-10-05, PR #17)

1. **A4's rule text overstated what `reduce_to_minimal` does.** It is not dust
   sweeping. It drops any term whose removal still certifies on `X_all_m`,
   the certify split included, so the law the joint gate tests was chosen
   using the certify rows (a second use of that split), and it can remove more
   than the constraint ideal. The path runs only when
   `reduce_mod_constraints` changes the candidate. In principle a dense
   candidate on constrained inputs could pass because the minimal law passes,
   while the candidate itself enters the certifying set. Campaign identity (E3)
   is the empirical guard. Nothing proves the gate target is the constraint
   quotient, and the unit test exercises `reduce_mod_constraints` only, not
   the engine path. **Open:** restrict the step to machine-scale coefficient
   dust, or show the drop set lies in the ideal; either change requires
   re-scoring.
2. **The default flip itself was not re-scored.** Every A4/A5 artifact was
   produced through the opt-in path (`coefficient_gate="joint_modulo"`,
   introduced in 8574a7e and unchanged through 95da103). eea8b2d changes only
   the default strings in `engine.py` and `passive.py`, plus docs and tests.
   The artifacts carry over **only** for that branch as it stands. Any later
   edit to the `joint_modulo` branch of the gate must be re-scored, not
   credited with these results.
3. **F3 stays failed as worded.** `fr3-seed322`'s tier-1 structural refusal
   among approximants is the same mechanism that made those splits cheap
   under marginal. It is not a near-pass.
4. **Scope of the closure.** It covers the clean, non-floor-dominated,
   per-candidate path only. The declared-noise and floor-dominated winner
   gates, `sqrt(f²)` twins under `accumulate`, and the dim ≥ 3 pre-pass
   (where this gate does not run) remain open. "Removed the dense-approximant
   class on every measured bank" is a statement about those banks, not a
   closure of false exactness.

## Amendment A6 — the joint test in the exact quotient space (registered 2026-10-05, before any A6 run)

**Why.** The PR #17 review found that `joint_modulo` does not test the
constraint quotient. Its `reduce_to_minimal` step drops any term that still
certifies on all rows, certify split included. A "dust only" sweep was tried
first and abandoned before any scored run. On the Gaia frame-rotation
candidates, the post-reduction residue reaches 2.6–3.9×10⁻¹³ of the law's
scale, above `MACHINE_REL` (2.2×10⁻¹³) and above the band in absolute terms.
It is the candidate's own fit noise along the constraint direction, not
rounding, so any dust threshold would be an arbitrary cut in a 12-order gap.

**Rule.** `coefficient_gate="joint_quotient"` (opt-in; `joint_modulo` stays
the default and stays reproducible). `joint_pinned(..., constraints=…)`
reduces each gated atom's derivative term modulo the constraint ideal in
exact arithmetic. The exact null space of those remainders is the set of
ideal directions, and the least-determined direction is sought only in their
orthogonal complement. The candidate is never reduced; no term is dropped;
`y` is not used beyond `check`. The rule applies when the law is linear in its
gated atoms. Otherwise, and whenever there are no constraints, it is exactly
A2's joint test.

**Pilot (disclosed):**
- A padded sphere law, truth + 0.2169·x1·(|x|² − 1): plain joint rejects it,
  `joint_quotient` passes it; the bare truth passes both.
- Gaia frame rotation, end to end: same law and α as `joint_modulo`, 10 s.
- Constrained `1/(2+x0)` on the sphere: marginal, `joint_modulo` and
  `joint_quotient` all give the exact truth.

**Observation, not part of the rule.** On the plane x0 + x1 + x2 = 1,
`input_constraints` returns two arbitrary quadratic combinations with snapped
rational coefficients instead of the linear constraint, because the quadratic
feature matrix contains several multiples of it. This predates A6 and also
feeds the constrained-input coherence path. Recorded as **open**.

### Banks

1. **`cs` (new):** 4 varieties (sphere, circle arc, hyperbola, plane) × 6 laws
   (float-coefficient linear, integer linear, rational, product, exp bait,
   sqrt bait), 400 rows each. Runner: `experiments/run_quotient_study.py`.
   Gates: marginal, `joint_modulo`, `joint_quotient`. Scored **on the
   variety** only: fresh in-region points for the domain claim, fresh points
   from a wider region of the same variety for the form claim.
2. **Campaigns:** all 14 scripts, `joint_modulo` and `joint_quotient`
   worktrees concurrent. Output `experiments/results/joint_quotient/campaigns/`.
3. **Unconstrained banks and null:** the code path is identical whenever
   `input_constraints` returns nothing. Verified directly, with no
   rediscovery: on every `fr*`/`frch*` case and on 200 null-style input
   draws, it must return `[]`.
4. **Suite** with `joint_quotient` as the default.

### Predictions

- **G1.** `joint_quotient`: 0 wrong on `cs`. marginal and `joint_modulo`:
  measured; any wrong certificate there is reported as evidence about the
  hole.
- **G2.** Campaigns: `joint_quotient` matches `joint_modulo` exactly (19
  certificates, same laws and α).
- **G3.** `input_constraints` returns `[]` on all 180 fresh unconstrained
  cases and all 200 null-style draws.
- **G4.** Suite: only the documenting test fails.
- **G5 (cost).** `cs` median `joint_quotient`/`joint_modulo` ≤ 1.5×; every
  campaign script ≤ 1.5× its concurrent twin.

**Decision rule.** If G1–G5 hold, open a PR making `joint_quotient` the
default; `joint_modulo` stays as the recorded former default. Otherwise the
default stays `joint_modulo` and the failure is recorded.

## Results — A6 (scored 2026-10-05)

**Empirical: G1, G2, G3 and G5 hold. G4 fails as worded, in the safe
direction. `joint_quotient` also loses one correct certificate that
`joint_modulo` issues, on the variety where constraint detection is
malformed.** No A6 run recorded suspended time.

| prediction | result |
|---|---|
| **G1** `cs` (24) | **holds.** `joint_quotient`: 22 certified, 0 wrong. `joint_modulo`: 23 certified, 0 wrong. marginal: 23 certified, **4 wrong**: circle exp-bait (form error 1.7×10⁻¹), plane rational (5.4×10⁻²), plane exp-bait and plane sqrt-bait (domain wrong) |
| **G2** campaigns | **holds.** 19 certificates in each arm, 0 differences across 13 result files |
| **G3** | **holds.** `input_constraints` returns `[]` on all 180 fresh unconstrained cases and all 200 null-style draws, so `joint_quotient` takes exactly the A2/A4 code path there |
| **G4** suite with `joint_quotient` default | **fails as worded.** It predicted the documenting test would fail; **no test failed** (431 passed), because PR #17 had already pinned that test to `"marginal"`. A stale prediction, not a regression |
| **G5** cost | **holds.** `cs` median per-case ratio `joint_quotient`/`joint_modulo` 1.00×; worst campaign script 1.02× |

**The difference between the two gates.** On `cs-plane-float-linear`,
marginal and `joint_modulo` certify the same 10-term quadratic
representative with 11-digit coefficients. It agrees with the linear truth on
the plane to 1.6×10⁻¹² (in-region and extended): a correct domain-restricted
certificate with an unreduced representative. `joint_quotient` refuses
structurally: a split certifies but fails the full-data gate. This is the
variety where `input_constraints` returns two snapped quadratic combinations
instead of x0 + x1 + x2 − 1, so the ideal that `joint_quotient` excludes is
not the variety's ideal.

**What this says about the review's hole.** On this bank it did not occur
empirically: `joint_modulo` issued no wrong certificate on any constrained
case. `joint_quotient` closes it by construction and costs one cell, and that
cost traces to the constraint-detection defect, not to the quotient rule.

**Decision.** Under the registered rule (G1–G5 must all hold), G4's wording
fails and the default stays `joint_modulo`. G4's miss is in the safe
direction, so the substantive trade is the plane cell. **Open, and the
natural next step:** make `input_constraints` return minimal-degree
generators (the linear constraint on the plane), then re-register the switch
with the plane cell named.

## B1 — graded constraint detection (registered 2026-10-05, before any B1 run)

**Defect (measured in A6).** `input_constraints` SVDs the full quadratic
feature matrix at once. A linear constraint l = 0 produces a (d+1)-dimensional
null space (l and each x_i·l), and the function returned two arbitrary,
snapped mixtures of it, never l itself. On the plane x0 + x1 + x2 = 1 that
cost `joint_quotient` a correct certificate. The same detector feeds the
default engine's constrained-input coherence.

**Rule.** `constraint_detection="graded"` (opt-in; default stays `"flat"`).
- Find linear constraints first, from the null space of (1, x_i).
- Then keep only the quadratic null directions outside the span of
  {1, x_i}·l for each linear l, using a rank-revealing basis.
- Report each degree's constraints in reduced row echelon form, so rational
  constraints come out with rational coefficients.
- Same tolerance, rationalization and cap as before.

Pilot (disclosed):
- plane → `−x0 − x1 − x2 + 1`;
- sphere, circle and hyperbola unchanged (up to sign);
- a line in 3-D → two linear constraints;
- x2 = x0² + x1 → `x0² + x1 − x2`;
- log-uniform random inputs → `[]`.

### Banks

1. **Detection directly:**
   - each `cs` variety, plus the line and quadric pilots, as unit tests;
   - all 180 fresh unconstrained cases and 200 null-style draws, which must
     return `[]`.
2. **`cs`:** four arms, all run in this study under the same load:
   `joint_modulo` and `joint_quotient`, each with `flat` (rerun, tagged
   `+flatB1`) and `graded`.
3. **Campaigns:** all 14 scripts, default gate (`joint_modulo`), `flat` and
   `graded` worktrees concurrent. Named cells: the Gaia P3 frame rotation
   and every certificate whose inputs carry an exact linear relation.
4. **Suite** with `graded` as the default.

### Predictions

- **H1.** Graded detection returns the minimal generators on every `cs`
  variety and on the pilot shapes, and `[]` on all 180 + 200 unconstrained
  inputs.
- **H2.** Campaigns: every certificate under `flat` is issued under `graded`
  with the same law and α.
- **H3.** `cs`, both gates: 0 wrong under `graded`, and no cell certified
  under `flat` is lost under `graded` for the same gate.
  `cs-plane-float-linear` certifies under `joint_quotient+graded`.
- **H4.** Suite with `graded` default: no failures.
- **H5 (cost).** `cs` median per-case ratio `graded`/`flat` ≤ 1.5× per gate;
  every campaign script ≤ 1.5× its concurrent twin.

**Decision rule.** If H1–H5 hold, open a PR making `graded` the default
detection. Re-registering the `joint_quotient` default switch, with the plane
cell named, follows as a separate step.

## Results — B1 (scored 2026-10-05)

**Empirical: H1–H5 all hold. Graded detection becomes the default.** No B1
run recorded suspended time.

| prediction | result |
|---|---|
| **H1** | **holds.** The minimal generators on every `cs` variety and pilot shape (unit tests, 4 pass); `[]` on all 180 fresh unconstrained cases and all 200 null-style draws |
| **H2** campaigns | **holds.** 19 certificates under each detection, 0 differences across 13 result files |
| **H3** `cs` | **holds.** 0 wrong in all four arms. No cell certified under `flat` is lost under `graded` for either gate. **`cs-plane-float-linear` certifies under `joint_quotient+graded`** (extended-region error 1.6×10⁻¹²): the cell `joint_quotient` lost in A6 is restored |
| **H4** suite with `graded` default | **holds.** 435 passed, 0 failed |
| **H5** cost | **holds.** `cs` median per-case ratio `graded`/`flat` 1.00× for both gates; worst campaign script 1.01× |

**Changed representatives, not changed laws.** On `cs-plane-int-linear` and
`cs-plane-product`, both gates under `graded` return a different
representative of the same law on the plane:
- `2·x0 − x1 + x2` becomes `−3·x1 − x2 + 2`;
- `x0·x1 + 3` becomes `−x1² − x1·x2 + x1 + 3`.

Each is the old law with x0 = 1 − x1 − x2 substituted (extended-region error
around 10⁻¹⁶). The engine now canonicalizes modulo the true constraint and
names `−x0 − x1 − x2 + 1 = 0` in the certificate, where under `flat` it named
a snapped quadratic that is not the variety's ideal. The cost is readability:
the elimination form is less natural than the input's own form.

**For the next registration (joint_quotient as default gate):** under
`graded`, `joint_modulo` and `joint_quotient` produce identical results on all
24 `cs` cells (median ratio 1.02×). That is evidence, not a registered pass.

## Review of PR #18 (2026-10-06): one fix landed, four preconditions named

Verdict: merge for B1, keep `joint_quotient` opt-in.

**Fixed in this PR: the certified payload did not carry the constraint.** The
review asked to confirm that the named constraint reaches the serialized
certificate. **It did not.** `mcp.core.recover`'s certified payload dropped
the certificate's notes entirely, so on the plane it returned
`law: "-3*x_1 - x_2 + 2"` with only bounds. That is an affine chart presented
as an ambient law. It predates B1 (the Gaia frame rotation's constraint never
reached the payload either), but B1 puts chart representatives on more cells.

`Certificate` now has a first-class `constraints` field, set by the engine
whenever it issues a domain-restricted certificate. The certified `recover`
payload serializes `constraints`, a `domain_restriction` statement and every
note. Unconstrained laws have no `constraints` key. Tested in
`tests/test_constraint_detection.py`. Verdicts are unchanged: only the payload
gains fields.

**Also added:** the A6 pilot is now a unit test
(`test_joint_quotient_excludes_exactly_the_ideal_directions`). The padded
sphere law fails plain joint and passes `joint_quotient`; the bare truth passes
both.

**Preconditions for re-registering `joint_quotient` as the default gate:**

1. **Gröbner basis.** `sp.reduced` is given the raw generator list. A single
   linear, a single quadratic, or an RREF set of linears is already a Gröbner
   basis. A linear plus a kept quadratic need not be: LM(l) can divide a term
   of the quadratic, and the remainder is then not unique. `cs` never mixes
   degrees. Reduce modulo `groebner(constraints)` instead.
2. **The float step after the exact null space.** Exact null vectors are cast
   to float, divided by v, then rank-cut at 1e-12. Over-ranking would exclude
   a real direction, and a non-identified candidate could pass. Make it a
   named prediction, with fixtures whose ideal directions are nearly
   degenerate.
3. **`max_constraints=2` truncates a graded basis.** Linears are appended
   first, so a quadratic beside two linears, or a third linear, is dropped.
   The cap was inherited from the flat detector, but graded detection is what
   makes a multi-generator return meaningful. This goes with the dim ≥ 3
   work.
4. **Three unshared cutoffs** in graded detection: null space at
   `tol_rel·s₀` (1e-10), implied span at 1e-9, new quadratic directions at
   1e-6 absolute when small. That is conservative against false constraints,
   but a weak genuine quadratic can vanish, and a slightly-off float RREF row
   could leave a residual of l above 1e-6 and reintroduce a quadratic
   mixture, the bug B1 fixes. Needs an adversarial fixture with a poorly
   scaled linear constraint.

## C1 — the four preconditions, and joint_quotient as default (registered 2026-10-06, before any C1 run)

**The fixes**, in code on this branch:
1. **Gröbner basis.** `_ideal_directions` and `reduce_mod_constraints` reduce
   modulo `groebner(constraints, order='lex')`, not the raw list.
2. **No float rank cut.** The ideal directions come from an exact null
   space, so their number is exact, and diagonal scaling by v cannot change
   rank. The complement comes from a complete QR with that rank. Over-ranking
   cannot occur.
3. **No truncation of a graded basis.** `max_constraints` now applies only to
   the legacy flat detector.
4. **One tolerance.** Nullity is decided only by `tol_rel` (1e-10). The
   implied span is built from the exact rationalized linear constraints and
   its rank is exact. The number of new quadratic constraints is a count,
   dim(quadratic null) − rank(implied). If any linear row has no exact
   rational form, the quadratic stage is skipped. The detector returns the
   **reduced lex Gröbner basis** (unique and readable); if the snapped
   generators are inconsistent (basis {1}) it returns no constraints.

**Fixtures** (`tests/test_constraint_detection.py`, `tests/test_escalation.py`;
19 pass on this branch). Run against master's code:
- **`test_two_linears_and_a_quadratic_all_returned` fails on master:** 2
  generators, the quadratic truncated.
- **`test_joint_quotient_on_mixed_degree_constraints` fails on master.** With
  the new detector but the raw-list reduction, it still fails; with the
  Gröbner reduction it passes. So finding 1 matters on its own.
- The poorly-scaled linear, non-ideal flat direction and ill-scaled
  coefficient fixtures pass on master too. They are regression guards, not
  evidence for a fix.

**Change under test.** The default `coefficient_gate` becomes
`"joint_quotient"`, with fixes 1–4 in place. **Baseline:** master at 34bfddc,
as committed (`joint_modulo`, the pre-fix graded detector). **Candidate:**
this branch with the gate default set to `joint_quotient`.

### Banks

1. **Detection:** all 180 fresh unconstrained cases and 200 null-style draws
   must return `[]` under the new detector.
2. **`cs`** (24 cells): baseline (from a master worktree) and candidate, run
   concurrently, case by case.
3. **Campaigns:** all 14 scripts, baseline worktree at 34bfddc and candidate
   worktree, concurrent.
4. **Suite** with the candidate default.

### Predictions

- **C-1.** The fixtures above pass on the candidate.
- **C-2.** `[]` on all 380 unconstrained inputs.
- **C-3 (`cs`).** Candidate: 0 wrong. No cell the baseline certifies is lost.
  Any law difference must be a representative of the same law on the variety
  (extended-region error ≤ 1e-6); each one is reported.
- **C-4 (campaigns).** The same 19 certificates with the same laws and α.
  Constraint generators may change form (reduced Gröbner basis); each change
  is reported.
- **C-5.** The suite with the candidate default: no failures.
- **C-6 (cost).** `cs` median per-case ratio ≤ 1.5×; every campaign script
  ≤ 1.5× its concurrent twin.

**Decision rule.** If C-1 to C-6 hold, open a PR that makes `joint_quotient`
the default gate, carrying fixes 1–4. Otherwise record the failure; the
default stays `joint_modulo`, and any fix that held still lands on its own
evidence.

## Results — C1 (scored 2026-10-06)

**Empirical: C-1 to C-6 all hold. `joint_quotient` becomes the default gate,
with fixes 1–4.** No C1 run recorded suspended time.

| prediction | result |
|---|---|
| **C-1** | **holds.** All fixtures pass on the candidate (in the suite below) |
| **C-2** | **holds.** `[]` on all 180 fresh unconstrained cases and all 200 null-style draws |
| **C-3** `cs` (24) | **holds.** Baseline (master 34bfddc, `joint_modulo`) and candidate (`joint_quotient`, fixes 1–4) both certify 23, **0 wrong**, nothing lost or gained, and **no law differences** |
| **C-4** campaigns | **holds.** 19 certificates in each arm, 0 differences in `certified`, `law`, `alpha_log10` or `abstain` across 13 result files. Not compared: the text of constraint generators in notes, which may now read in reduced Gröbner form |
| **C-5** suite with candidate default | **holds.** 442 passed, 0 failed |
| **C-6** cost | **holds.** `cs` median per-case ratio 1.00×; worst campaign script 1.03× (macro) |

**What changes for users.** The default exact-coefficient gate now tests joint
pinning with the constraint ideal excluded exactly, so the #17 concern (a
`reduce_to_minimal` target chosen with the certify rows) no longer applies to
the default path. `joint_modulo` stays available and reproduces A4/A5.
Detection returns a reduced Gröbner basis with no truncation, and
constrained-coherence canonicalization reduces modulo a Gröbner basis.

**Still open:**
- gates whose law is non-linear in its gated atoms get no ideal exclusion
  (they fall back to the plain joint test, so they refuse more, never accept
  more);
- the declared-noise and floor-dominated paths;
- the dim ≥ 3 pre-pass.
