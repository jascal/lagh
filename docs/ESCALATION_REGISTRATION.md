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
