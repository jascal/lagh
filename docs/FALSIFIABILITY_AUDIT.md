# Falsifiability audit

Status: **empirical** for executed witnesses; **open** for coverage explicitly
missing below. No theorem is promoted to `proved` by a passing numerical check.
Scope is the four registrations existing at the start of this pass, plus the
prospective registration for this work. A counter-input tests the discriminator;
it does not retroactively falsify an observation of a different, frozen input.

Reproduce: `OPENBLAS_NUM_THREADS=1 .venv/bin/python -m experiments.run_falsifiability`.
The baseline artifact records the pre-fix behavior of core revision
`582227fc59f22d6d1a5d09fbec7a575a8a2224df`. Of 44 specified refusal
witnesses, 36 initially refused as expected. Eight exposed acceptance or crashes:
NaN target, NaN/infinite band, empty check, empty significance, missing variance,
missing axis, and exactly-bulk Faxen. All 44 now give the specified refusal.
These counts establish witness coverage, not universal correctness.

## Executed gate counter-inputs

| Gate | Input making the positive claim fail | Required result | Evidence |
|---|---|---|---|
| `check/wrong value` | x=linspace(1,2,40), y=2x, law=x, eps=1e-6 | 40 misses | empirical; witnessed |
| `check/undefined` | law=log(-x), x in [1,2] | 40 uncovered | empirical; witnessed |
| `check/NaN target` | NaN target on 40 exact x rows | refuse | empirical; witnessed |
| `check/NaN band` | NaN band on 40 exact x rows | refuse | empirical; witnessed |
| `check/infinite band` | infinite band on 40 exact x rows | refuse | empirical; witnessed |
| `check/negative band` | negative band on 40 exact x rows | refuse | empirical; witnessed |
| `check/empty` | X=(0,1), y=eps=[] | refuse | empirical; witnessed |
| `vacuous` | y=x, eps=3 | zero law covers: vacuous | empirical; witnessed |
| `significance/empty` | y=[], eps=[] | no evidence: log bound 0 | empirical; witnessed |
| `significance/evidence` | y=x, eps=10 | no evidence: log bound >=0 | empirical; witnessed |
| `float_pinned` | law=2.1*x, y=2.1*x, eps=.1 | unidentified coefficient | empirical; witnessed |
| `pinned` | x in [1,1.001], law=x**(3/2), sigma=.2, eps=.1; probe [1,2] | rival exponent fits narrow domain but diverges on probe | empirical; witnessed |
| `minimal` | law=x+1/1000, y=x, eps=.01 | junk term removable | empirical; witnessed |
| `parameter_interval` | law=2*x, y=2*x, eps=100 | unbounded at max_rel=.5 | empirical; witnessed |
| `admissible_interval/infeasible` | A=[[1],[1]], y=[0,1], eps=0 | infeasible | empirical; witnessed |
| `admissible_interval/unbounded` | A=[[0]], y=[0], eps=0 | no coefficient bound | empirical; witnessed |
| `admissible_interval/coefficient budget` | A=[[1]], y=[2], eps=lambda cap:cap, cap=1, iters=1 | coefficient declaration unverified | empirical; witnessed |
| `admissible_functional` | A=[[1],[1]], y=[0,1], eps=0, W=[[1]] | infeasible | empirical; witnessed |
| `max_divergence/overlap` | identical x laws, only 7 probes | insufficient overlap: infinity | empirical; witnessed |
| `arbitrate_significance/margin` | rivals x and 2x on 40 rows | no winner | empirical; witnessed |
| `arbitrate_significance/rival evidence` | x versus polynomial with 14 coefficients, 40 rows, eps=1e-6 | margin alone cannot defeat evidence-bearing rival | empirical; witnessed |
| `coherent` | laws x and 2x, probe [1,2], yscale=1 | two classes | empirical; witnessed |
| `input_constraints` | seed 0 independent uniform 40x2 design | no exact constraint | empirical; witnessed |
| `conjoin_determination/domain` | qualified domains x>0 and x<0 | refused | empirical; witnessed |
| `conjoin_determination/contradiction` | a in [1,2] AND [3,4] | contradiction | empirical; witnessed |
| `determination/resolution` | a in [-1,1] | not resolved | empirical; witnessed |
| `acf_timescale` | zeros(1000), dt=.001 | unreadable | empirical; witnessed |
| `axis_gate` | zeros(1000), dt=.001, theta_ref=1 | refuse | empirical; witnessed |
| `gate_record/unreadable` | {'a': {'ratio': None}} | unreadable | empirical; witnessed |
| `gate_record/contaminated` | {'a': {'ratio': 0.1, 'var_ok': False}} | contaminated | empirical; witnessed |
| `gate_record/single` | {'a': {'ratio': 0.7, 'var_ok': True}} | unresolved | empirical; witnessed |
| `gate_record/spread` | {'a': {'ratio': 0.5, 'var_ok': True}, 'b': {'ratio': 0.7, 'var_ok': True}} | unresolved | empirical; witnessed |
| `gate_record/missing variance` | {'a': {'ratio': 0.7}, 'b': {'ratio': 0.7}} | unresolved | empirical; witnessed |
| `gate_record/missing axis` | {'a': {'ratio': 0.7, 'var_ok': True}, 'b': {'ratio': 0.7, 'var_ok': True}, 'c': {'ratio': None}} | unresolved | empirical; witnessed |
| `retention/fit bins` | zeros(20), fs=100, fit_range=(10,11) | ValueError insufficient bins | empirical; witnessed |
| `faxen_height/0.9` | radius=1, ratio=0.9 | no finite height | empirical; witnessed |
| `faxen_height/3.0` | radius=1, ratio=3.0 | no finite height | empirical; witnessed |
| `faxen_height/1.0` | radius=1, ratio=1.0 | no finite height | empirical; witnessed |
| `faxen_height/1.02` | radius=1, ratio=1.02, se=.1 | no finite height | empirical; witnessed |
| `attribute_deviation/(0, 1, 1)` | (0, 1, 1) | unreadable | empirical; witnessed |
| `attribute_deviation/(0.019, 0.885, 3.797)` | (0.019, 0.885, 3.797) | unattributed | empirical; witnessed |
| `attribute_deviation/(0.5, 1, 2)` | (0.5, 1, 2) | unattributed | empirical; witnessed |
| `attribute_deviation/(0.1, 1, 1.04)` | (0.1, 1, 1.04) | unattributed | empirical; witnessed |
| `attribute_deviation/(2, 3, 4)` | (2, 3, 4) | unattributed | empirical; witnessed |

## Remaining branches and non-gates

| Surface | Counter-input / limit | Status |
|---|---|---|
| coverage_factor | delta=0 or 1 raises; the inverse budget is algebra, not coverage validation | open: no new tail-distribution coverage experiment |
| coverage_budget | unknown tail raises; kappa=4, n=2000 yields budget >1 | arithmetic diagnostic; no distributional validation |
| epsilon, band | negative/nonfinite declarations reach check's invalid-band refusal | assembly helpers, not independent gates |
| free_dof, free_atoms, gated_atoms | integer x has no fitted atoms; 2.1*x does | bookkeeping, not tests of identifiability |
| sample_box | a sampled probe cannot exclude an unsampled localized discrepancy | no coverage of the entire continuous box |
| reduce_to_minimal | x+x² with y=x+x², eps=1e-6 must retain both terms | open: witness specified, not newly executed |
| refit_minimal | x+x²+x³ on [1,2] at eps=1e-12 cannot discard support | open: witness specified, not newly executed |
| reduce_mod_constraints | a nonpolynomial input remains unchanged | algebraic rewrite, not new evidence |
| invariant_content | search set {x} omits another certifying law 2x at eps=2 | no coverage over unsearched laws; existing qualifier essential |
| domain_qualifier | different predicate strings fail conjoin (above) | labels alone establish no physical domain coverage |
| Retention.fraction, _full, diode_response | forward response parameters fixed; wrong diode can fit residual compensation | model evaluation/identity, no independent response validation |
| band_loss half spread | two identical halves have zero spread even if both use a wrong diode | repeatability, not accuracy; short-record 5% fallback is declared |
| drive_frequency, lockin | zero reference/zero response has no identifiable drive | open: missing signal-quality coverage; not an independent calibration gate |
| drive_scale corner gate | zero residual and no declared corner gives scale=None | open: specified counter-input; existing synthetic drag witness below |
| equipartition_ratio, local_drag | multiply declared length by 2: stiffness/drag change by 1/4 | identities conditional on units, no validation of those units |
| realized_diffusion | strides=[] gives no plateau; flatness with ret fitted on the same spectrum is not independent response evidence | open: no new flatness witness executed |
| axis_gate variance annotation | normal theta with var_ratio=4 can pass the timescale gate; gate_record must refuse contamination | timescale-only pass, not a clean-axis claim |
| attribution uniqueness/materiality/tolerance | frozen refusal, stiffness overlap, sub-floor excess, mixed-change witnesses above | empirical; no joint uncertainty coverage for 8% |

## Registered predictions and four seed defects

Every prediction and amendment is accounted for below. Protocol promises cannot
be tested by substituting a physical input: their falsifier is a provenance
violation. That distinction prevents counting a metadata audit as scientific
validation. No sealed benchmark is re-read.

| Registration / prediction | Input or intervention defeating the claim | Disposition |
|---|---|---|
| LLM §1 frozen baselines | a different authoritative historical benchmark count/score | historical declarations, not predictions; 239 already amended to 240 |
| LLM §2 train-only, GT-after-submissions | a submission timestamp later than GT access | open: no new provenance replay |
| LLM §2 A full-data/K3, B empirical, >600 cap | one omitted full-data row with residual >eps; a conjecture tagged proved | core wrong-row witness; protocol coverage open |
| LLM §2 local judge lower-bounds GPT judge | GPT rejects a pair which stripping constants accepts (e.g. x and x² if the exponent is stripped) | WITHDRAWN: no coupling to the stochastic judge, no lower-bound theorem |
| LLM §2 accuracy/NMSE/accounting, A-wrong=0 | replace submitted 2x by x² on x in [1,2], score wrong; any A wrong kills invariant | metric falsifiable; historical score not rerun |
| LLM §3 noisy data never certify | bounded perturbation < machine floor is noisy yet passes | WITHDRAWN universal exactness claim; finite-band agreement only |
| LLM §3b timeout demotes/count240 | timeout path emitting A without a certificate; missing row | open: protocol witness not executed |
| LLM §4 one read/no tuning/report both | changed grammar after scores; omitted B failures | open: provenance commitments, not physical tests |
| Gravity §1 baselines | mismatched historical source | declarations, not a claimed future score |
| Gravity §2 fixed planner/budget | request 101 observations under cap100 | open: no new protocol execution |
| Gravity §2 playbook/fit provenance | mass value without source fit, or alpha on noncertifying fit | open: no new protocol execution |
| Gravity §2 synthetic-only development | scenario parameters in development inputs | open: provenance coverage |
| Gravity §3 dtype noise first | float32 observations declared float64 epsilon after seeing errors | open: protocol coverage |
| Gravity §4 metrics/types/unexpressible=0 | a fabricated correct answer for an unexpressible type or omitted miss | falsifiable accounting; no benchmark reread |
| Gravity §5 walk-away/one-read | unsupported free-text task processed after opening scores | open: historical protocol, no new coverage |
| Refusal P1/P2 signed residual and domain | perturb fitting-only row; scramble the reported original row indices | existing fixtures; new check keeps indices at evaluation boundary |
| Refusal P3 no decision/law change | clean 2x returns changed law or upgraded refusal | finite-input regression required; invalid-band soundness fixes explicitly exempted |
| Refusal P4 identical observational causes | no input separates two literally identical observation arrays | tautology/identifiability limitation, WITHDRAWN as empirical prediction |
| Refusal amendment finite peers/invalid indices/scalar eps | [finite,NaN] predictions; duplicate indices; scalar .1 | falsifiable implementation contract; focused regression file |
| Refusal amendment64/counts/priority/name/no partial | 20000 rows, late largest miss, missing candidate or partial key | falsifiable implementation contract; existing fixtures |
| C4 P1 | (.019,1,3.797) permits slow compatibility, whereas frozen .885 must refuse | executed frozen-triple witness; no historical reinterpretation |
| C4 P2 | short-lag diffusion .885 while Var=3.797 (above) | falsifiable; original second-half-only result preserved |
| C4 P3 positive/stationary excess and crossing | zero excess; independent low-rate sinusoid never crossing within window | stationary interpretation already killed; window statistic remains conditional |
| C4 P3 clean/drag/stiffness/scale/mixed controls | (.5,1,2) stiffness/slow overlap; (.1,1,1.04) sub-floor excess | executed; formerly unique stiffness reading withdrawn |
| C4 P4 white excess | flat detector noise passed through a low-pass filter while model omits filter | old statistic measures transfer: WITHDRAWN as whiteness test, no coverage |
| C4 P5 thermal+floor moment prediction | add independent tone outside thermal/noise fitting windows | falsifiable composite prediction, historically failed |
| C4 P5 equal-lag decimation | alternating/high-frequency component before subsampling | equality survives wrong response: tautology at paired samples; independent validation WITHDRAWN |
| C4 initial amendment materiality | Var=1.038 vs 1.09 with D=1 and free theta=.1 | sub-floor witness executed; threshold declaration is not an error estimate |
| C4 stored-diode/power amendment | change independent diagnostic power rather than fit diode on residual | falsifiable provenance/model choice; raw records closed |
| C4 alias-retaining half-PSD prediction | second half adds a large tone absent from first half | falsifiable repeatability; cannot identify sampling causality |
| C4 PR13 train-band exclusion | bump only 30–40kHz validation repeatability while held-out bands agree | same-band vote is protocol defect; exclusion required, not whiteness proof |
| C4 PR13 40–43passband/rolloff | response attenuates even inside presumed passband | no transfer coverage: whiteness and sigma remain open regardless of composite fit |
| C4 PR13 free theta / materiality | theta=1±small at fixed Var>1.08,D=1 | no theta admission; specified boundary fixture remains in C4 tests |
| C4 PR13 axis gate / instability | unreadable theta; half rates separated >8% | whole-record missing-evidence witnesses executed; no physical-cause upgrade |
| C4 PR13 missing fits/strides/negative subtraction | fit exception, [], or thermal variance above total | existing guarded synthetic fixtures; no new raw analysis |
| C4 PR13 provenance under -O | invalid declared diode provenance with assertions disabled | existing subprocess fixture, not a numerical theorem |
| C4 PR13 centering/clipping/unknown bias | near-constant slow offset removed by sample centering | population amplitude unidentifiable; no coverage, claim remains window-scoped |
| C4 R-C3 and frozen replay | two possible stage units giving same stored numbers | nonidentifiability, not empirical prediction; remains open |
| Seed C2 circular corner | true drag2.3, supplied bulk-referenced corner | existing synthetic bead witness; assumed-corner result is not independent drag evidence |
| Current A1 | any row in this table lacking an executable witness | universal executed-coverage prediction FAILED; no coverage is named, not counted as pass |
| Current A2 | NaN y/eps, missing axes, bulk Faxen | eight baseline defects; repaired witnesses above |
| Current B1 | short or contaminated record with truth/half error >8% | prospective experiment; remains open until scored |
| Current C1 | pilot chooses worse function on independent certification paths; any path leakage | prospective experiment; remains open until scored |
| Current D1 | perturb unit float by1e-5 | prospective artifact-comparison witness |

A1 is restated in its registration: this is an inventory plus an executed
bounded witness set, not complete branch or scientific coverage. Four seed
failures have different counterfactuals but the same lesson: a repeatability
statistic or identity cannot validate the mechanism it assumes.

The four seed interventions were additionally executed during final verification
in `falsifiability_seeds.json` (synthetic only): supplied bulk corner gives
scale1.986 versus truth4.568; equal-lag paired error is exactly0 despite the
high-frequency tone; a white detector plus omitted stopband leaves61.75%
apparent nonwhiteness; changing training-repeatability alone flips the old
vote but leaves disjoint validation unchanged. Changing a disjoint band to
twice its prediction defeats the revised vote. These are empirical witnesses,
not new C2/C4 analyses. B and C dispositions are now recorded in the prospective
registration and PRECISION_AND_TEST_SELECTION.md.

Retained-row evidence added one intermediate-design failure: mutating caller
arrays changed the later diagnostic although the verdict still reported a miss.
The mutation input was exhibited before its fix; CheckResult now snapshots those
values. `residual_mutation_baseline.json` and `residual_consumers.json` preserve
the rejected aliasing design and corrected behavior. The reproducible consumer
script exercises actual refusal paths using an explicitly fixed wrong proposer.
Original44 witnesses and this dataflow test are separate counts.

## Review: witness specificity is itself a failure mode

**Empirical:** the equality-only Faxen witness hid a surviving crash at1+1e-12.
The PR14 revision replaces that limited interpretation with a neighborhood and
uncertainty-state sweep, and separates conditional inversion from a measured
height. Other newly exhibited inputs include malformed check domains, NaN
uncertainty, clean-looking axes with missing variance evidence, and rounded
storage boundaries. See `PR14_REVIEW_REVISION.md` and the pinned baseline/current
review artifacts. The44 original witnesses are not claimed to exhaust these
input families. No absent coverage is retroactively credited.
