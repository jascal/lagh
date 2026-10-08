# Certification boundaries — BND1 outcome

2026-10-08. **Empirical, bounded study.** Registration and baseline runner were
committed at `4701e5b`; baseline evidence and repair predictions at `6e3000d`;
production repair at `20c9bd0`. The repair and fresh confirmation share identical
production sources and environment. No hypothesis vocabulary, tolerance or
significance threshold was relaxed.

The study found and repaired a bypass of the selected coefficient gate and
interval/reporting defects. It also measured a limit that a gate change
cannot erase: noisy/floor-limited observations support fits that differ from
the generating law. Those outputs now explicitly claim finite-row band
consistency. **The change is not an exact-recovery gain.**

## Results and stopping rule

All three registered arms completed, 120 discovery attempts each. Each row had
240 observations, a fixed 60/20/20 split, a 20-second discovery limit, and a
separate diagnostic budget. Four workers, one BLAS thread per worker. Baseline
and replay use seeds 710/711; confirmation uses 810/811. The full bank and its
fixed stopping rule are in [the registration](CERTIFICATION_BOUNDARIES_REGISTRATION.md).
These are twelve designed families (including the negative control) repeated
across regimes and seeds, not 120 independent physical laws per arm.

| Regime (24 attempts per arm) | Baseline/replay: certified / refused / timed out | Fresh: certified / refused / timed out | Generator mismatches: baseline/replay → fresh | Median discovery seconds: baseline → replay → fresh |
|---|---:|---:|---:|---:|
| Clean | 18 / 2 / 4 | 20 / 2 / 2 | 0 → 0 | 3.14 → 3.16 → 4.74 |
| Relative noise 1e-6 | 5 / 11 / 8 | 7 / 9 / 8 | 5 → 7 | 6.56 → 6.71 → 6.32 |
| Relative noise 1e-3 | 4 / 14 / 6 | 4 / 14 / 6 | 4 → 4 | 3.84 → 3.90 → 3.98 |
| Absolute floor 1e-6 | 13 / 9 / 2 | 13 / 9 / 2 | 1 → 1 | 1.37 → 1.38 → 1.56 |
| Absolute floor 1e-3 | 0 / 18 / 6 | 0 / 18 / 6 | 0 → 0 | 3.86 → 3.88 → 3.89 |
| Total | **40 / 54 / 26** | **44 / 52 / 24** | **10 → 12** | |

There were no exceptions, missing rows, or missing diagnostic scores. No null
target certified (10 per arm; repeated baseline/replay inputs are not independent
null trials). Every generating law fit its declared observed-row band. Every
returned law also passed an independent recheck on all supplied observations.

The replay has **zero law or verdict differences** from baseline. Its aggregate
discovery time was 809.1 seconds versus 807.5 seconds baseline (summed per-job
elapsed time, not wall time or a controlled microbenchmark). The audit reports
refusal and timeout rates separately over all 24 attempts per regime, as well as
scoring runtime. Timeouts are censored outcomes, never counted as refusals or
successful exact recovery. Clean rational cases remain timed out at this budget;
several multidimensional/noisy cases do too. Their reach remains unmeasured here.

## False exactness versus finite-data consistency

The mismatch score preserves polynomial-support differences and coefficient
errors above 1e-8, and otherwise a normalized extended-region error above 1e-8.
Symbolic equality is recorded separately: passing the numerical threshold does
not prove exact equivalence. Original-region and extended-region discrepancies
are separate; extension never enlarges the certificate's finite domain.

The ten baseline mismatches include measured coefficient shifts, a small
quadratic term disappearing under noise, and nearly collinear declared columns
collapsing to a simpler fit. Three already had a conditional parameter record;
seven lacked such qualification. Thus “ten wrong certificates” would conflate
different claims. All ten were compatible with the supplied observation bands.
The public `pinned` wording and absence of machine-readable scope could nevertheless
present these fits as exact identification, which the indistinguishable-coefficient
witness directly refutes.

After repair, all 22 noisy/floor certificates in replay and all 24 in fresh
confirmation carry `claim.kind = finite-data-consistency` and
`exact_form_identified = false`. Public recover/verify use `strength=consistent`
for those regimes. Clean `gate-qualified-fit` records operational pinning and
also explicitly disclaims proof of the generator. The mismatch counts remain
10 and 12; correcting the claim does not turn them into exact discoveries.

No confidence probability for exact structure is claimed. The existing alpha
still concerns chance agreement under its declared null, not exact-form error
or adaptive family-wide calibration.

## Supported coefficient ranges and interval defects

The natural-bank diagnostics measured numeric-atom slices (which can include
an exponent). They are separate from the independent coefficient-position
slices returned by the PDE adapter.

1. **Unchecked endpoint:** baseline bisection set its near endpoint to half the
   first rejected step, even when that half-step had never passed. Independent
   replay found 82 infeasible endpoints across 16 certified cases. Repair starts
   from the checked central value. All **170 replay** and **194 fresh** diagnostic
   endpoints passed independent re-evaluation. Median clean slice width changed
   from 7.0e-10 to 9.3e-13; the wider historical numbers were not supported.
2. **Conditional is not marginal or simultaneous:** the controlled two-coefficient
   witness has valid individual endpoints but an invalid joint corner, and a
   feasible alternative vector outside both slices. Natural replay still has
   28 cases with an infeasible simultaneous corner; confirmation has 32. Those
   are expected limitations, not newly claimed boxes. Metadata explicitly sets
   `joint_box=false` and `marginal_coverage=false`.
3. **Unit/repeated coefficient positions:** the PDE adapter previously reported
   unit coefficients as exactly [1,1], and shared numeric atoms could alias
   distinct positions. Each linear coefficient position is now perturbed
   independently, including ±1, with callable bands evaluated on the actual
   perturbed expression.
4. **Independent midpoint substitution:** the engine no longer replaces several
   coefficients with separately computed slice midpoints. It retains the checked
   law and reports slices around it, avoiding an unchecked joint replacement.

The diagnostic can match a generator coordinate unambiguously only when a single
replacement establishes symbolic equality. All 54 such replay matches and all
62 fresh matches were inside their slices. This selected, conditional statistic
is **not marginal coverage** for the general fitted vector. A missing bound is
a failed/limited local search, not proof of global unboundedness. For nonlinear
parameters or candidate-dependent bands, bisection certifies its sampled endpoints;
it is not a proof that every interior value is feasible.

The existing PDE regression supplied a second, natural demonstration: the full
true vector (0.1, 0.5) fits, while replacing one fitted coefficient alone by its
true value fails with the other fitted coefficient held fixed. Its corrected
slices are approximately [0.0999999999857, 0.0999999999945] and
[0.4999999999272, 0.4999999999717]. The old test's marginal-coverage assertion
therefore depended on the invalid endpoint inflation. The retained
`pde_conditional_witness.json` records both facts. The regression now checks
true-vector feasibility, support recovery, fitted-center inclusion, and actual
endpoint feasibility; it does not widen the slices to make the old assertion pass.

## Execution-path contract

| Path | Supported claim | Evidence / limitation |
|---|---|---|
| Clean ordinary tiers | Gate-qualified candidate on checked rows; constraint-qualified where named | Existing quotient regressions plus 18 replay / 20 fresh clean certificates. No universal exactness theorem. |
| Clean multidimensional pre-pass | Same selected coefficient gate as ordinary tiers | Controlled proposal passes marginal and fails joint; legacy arm remains reproducible. Quotient constraints travel with the result. Natural clean pre-pass monomials retain reach. |
| Noisy ordinary tiers | Finite-row band consistency; dense-channel restrictions preserved | Nine replay / eleven fresh noisy certificates; every one differs from the generator under the registered score. No exact-structure claim. |
| Noisy multidimensional pre-pass | Finite-row band consistency | Noisy monomial coefficients differ from the generator. No blind application of a clean gate to uncertain coefficients. |
| Loose-floor tier loop | Finite-row band consistency after existing parsimony/winner gates | Thirteen certificates per repaired arm at 1e-6; one generator mismatch each. At 1e-3: 18 refusals and 6 timeouts. Pre-pass remains disabled for loose floors. |
| Declared linear/basis and PDE error bands | Consistent law plus separately scoped conditional ranges | Independent/near/exact-collinear banks, endpoint witnesses, unit/repeated positions, callable-band regression. No joint coverage from conditional slices. |
| Passive, tiny-data and active recover; verify | Preserve the applicable scope in the public result | Transport regressions for all recover modes; real noisy/floor verify regressions. Active transport is tested, not a new acquisition campaign. |
| PDE and benchmark submission | Preserve consistency scope; no “exact certificate” wording | PDE result includes claim and interval scope; submission includes claim. No benchmark score upgrade. |

Specialized C6 integer and C7 Lévy grammars are not new discovery banks in BND1
(max tier was 3); their successful engine results also receive the conservative
claim record, and existing regressions remain part of validation. No new reach or
soundness claim for those grammars is inferred from this study.

## Artifacts and reproduction

`experiments/results/certification_boundaries/` contains the three append-only
arms, source/environment manifests, original and repaired mechanism witnesses,
the original failing regressions, the PDE conditional witness, and
`audit_verified.json`. The audit regenerates input hashes, independently checks
returned laws and slice endpoints, verifies all 120 expected identities per arm,
and checks that repair/confirmation used identical sources. Tampering tests reject
a missing case, a substituted wrong law despite saved success flags, and an
overstated exactness claim.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.run_certification_boundaries --output /tmp/bnd1-fresh-replay
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m experiments.score_certification_boundaries --output /tmp/bnd1-audit.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/pytest -q tests/test_certification_boundaries.py tests/test_boundary_audit.py tests/test_pdesystem.py
```

The first command is a reproduction of seeds 710/711, not a new independent
study. Use isolated checkouts at the recorded revisions for historical baseline
reproduction. Existing artifacts are never overwritten. The bounded campaign
ends at the three registered arms; broader noisy identifiability and timeout
reach remain open, explicitly outside the supported claims.

## Final validation and completion audit

**455 tests passed, zero outstanding failures or skips**, across all 40 current
test files. The initial 39-file run had 451 passes and the single obsolete PDE
coverage assertion described above; the corrected PDE file passed 7/7 and the
new artifact-audit file passed 3/3. This is an aggregate of the full run and
targeted reruns, not a claim that the original full invocation exited green.
`validation/` retains that failure, corrected output, per-file source hashes,
and the effective results. Bug-class Ruff (`E9,F63,F7,F82`) and whitespace checks
passed. Current production/registration hashes match both scored candidate arms.

| Goal requirement | Authoritative evidence |
|---|---|
| Register before measurement; bounded stopping rule | Registration commit `4701e5b`, 120-case manifests, exactly three completed discovery arms |
| Cover noise, floors and multidimensional pre-pass | Five-regime bank plus controlled pre-pass witness; route contract above |
| Separate false exactness, coefficient support, refusal and runtime | Preserved generator mismatches and symbolic equality, conditional slices and joint-corner diagnostics, all-case status rates and separate discovery/scoring times in `audit_verified.json` |
| Change gates only after measured evidence | Baseline pre-pass witness and failing regression precede repair `20c9bd0`; noisy gates unchanged |
| Limit unsupported claims and preserve public scope | Claim fields and transport tests; consistency for uncertainty, explicit slice limits; no global exactness/coverage assertion |
| Verify repairs and retain negative results | Same-input replay unchanged, fresh-seed confirmation, independent checks of all hashes/laws/endpoints, original failures retained |
| End the finite campaign | 360 attempts completed; no additional seed search or capability iteration; timeout reach and general noisy identifiability remain documented limits |
