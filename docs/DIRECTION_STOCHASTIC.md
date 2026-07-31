# Direction: certified discovery on partially determined systems

**Scoped 2026-07-29 from a user design.** The catalogue, the three-regime scoring
and the zero-confident-wrong ranking are the user's; this document adds the
enabling mechanism, the one decision that has to be registered before any code,
the strategic classification, and three things the minimal suite was missing.

Every benchmark this program has met assumes a deterministic law plus
observational error. A stochastic system is different in kind: the law is
*partially* determined, and the honest output is a skeleton plus a noise
structure plus an explicit statement of what the data cannot pin. That is the
shape lagh already produces — exact / interval / abstain — so a suite scored on
those regimes tests the discipline directly rather than through a proxy.

## The enabling insight: this is the weak form again

An SDE has **no pathwise law to certify**. A realization of
`dX = a(X)dt + b(X)dW` satisfies no deterministic pointwise relation; the
increment is random. But the Itô/generator form converts it into precisely the
structure the PDE arc finished building. Testing against a test function φ and a
smooth f, and integrating over a window,

    ∫ φ df(X)  =  ∫ φ·[a f' + ½ b² f''] dt  +  ∫ φ f' b dW

and the last term is a **martingale with mean zero and variance
∫ φ² f'² b² dt** — a COMPUTABLE error scale, not an assumed one. So drift and
diffusion sit in a linear system whose error is declared, which is exactly the
errors-in-variables shape `weakform.PatchEpsilon` assembles. Even the
circularity carries over unchanged: the error scale depends on `b`, which is
also unknown, and the band is already assembled PER CANDIDATE as the quadratic
form `a'Ga` for that reason.

Consequences worth stating plainly:

* Level 0–1 is mostly new FIELDS, not new engine. The generator columns are
  `∂^α(g(fields))` terms with the same by-parts discipline.
* No differentiated data enters, for the same reason it does not in the PDE arc:
  the derivatives live on φ.
* The multi-solution holdout transfers directly — independent trajectories are
  exactly the "several solutions" the PDE curriculum already requires, and a
  single trajectory should refuse for the same reason it does there.

## The decision that must be registered first

lagh certifies `|pred − y| ≤ ε` EXHAUSTIVELY, over every point. A martingale
increment is unbounded, so no finite ε holds for every point almost surely. The
certificate's kind therefore changes, and the choice is a registration decision
rather than an implementation detail. Three options:

1. **Coverage-factor calibration.** lagh ALREADY makes a probabilistic band
   claim — `KAPPA = 4` in `certify.epsilon` is 4σ coverage, and `alpha` bounds
   chance agreement. So this is not a new epistemology; it is the question of
   what κ means when the dominant term is intrinsic rather than observational,
   and how a union bound over n patches enters α. Continuous with everything
   that exists.
2. **Certify only almost-sure objects** — quadratic variation, martingale
   residuals, generator identities. Genuinely pathwise, materially narrower.
3. **A distinct statistical-certificate class**, labelled and never compared
   against a deterministic α — the discipline `statecert.alpha_kind` already
   uses for state certificates.

**Registered choice: (1), with (2) as the fallback where (1) cannot be
calibrated, and (3) deferred out of the minimal suite.** A run that cannot state
its coverage refuses.

### ANSWERED 2026-07-29: what κ means for an intrinsic term

The question option (1) left open — what κ means when the dominant term is
intrinsic, and how a union bound over n patches enters α — is settled in
`certify.coverage_factor`, and in full in `STOCHASTIC_CHECKER.md` §1. The short
form: the tool is the **exponential martingale inequality**, not a Gaussian
quantile. For any V > 0,

    P( |M_p| ≥ κ√V  and  ⟨M_p⟩ ≤ V )  ≤  2 exp(−κ²/2)

which is pathwise and assumes no Gaussianity, no CLT and no independence.
Exhaustiveness over n rows costs a union bound, and inverting for a declared
false-abstain budget δ gives **κ(n, δ) = √(2 ln(2n/δ))**. Four consequences:

* **The continuity is numerical, not rhetorical.** κ(300, 0.05) = **4.334**
  against the existing constant `KAPPA = 4`. What the constant never had is the
  budget: κ = 4 over 300 rows buys δ = 0.20 under the martingale bound, 0.019
  under an exact Gaussian tail — and past ~1491 rows it states **no coverage at
  all**, which is worth knowing since weak-form runs carry thousands of patches.
* **κ becomes a function of the run** and is reported with (n, δ). The checker
  refuses an answer whose κ does not meet its own stated budget at its own row
  count, which is what makes this answer load-bearing rather than decorative.
* **The band's scale is measured and never by the candidate.** V bounds
  ⟨M_p⟩ = ∫φ²f'²b²dt, which depends on `b` — one of the objects being
  discovered. A candidate that set its own band could widen it at will: the
  loose-ε failure mode, self-serving. So V comes from realized quadratic
  variation (a functional of the DATA), with only the φ²f'² weights coming from
  the test function, and a residual-derived scale is refused as circular.
* **The union bound enters α through the DISJOINT patch count, and the trade is
  favourable.** q^h assumes independent held-out rows and overlapping patches
  share the driving path, so the honest count is the number of patches with
  disjoint support and it must be stated. Then κ grows like √(2 ln n) — 4.33 at
  300 rows, 5.08 at 10 000 — while log α gains linearly in h. See S7.

## Strategic classification: a DEV instrument, not a win vehicle

`STRATEGY.md` separates dev targets from the sealed blind read, and this suite
is unambiguously the former — by construction, not by concession. A benchmark
whose primary ranking penalises confident-wrongs and rewards correct abstention
is one **only this program can play**: SINDy-for-SDEs, Kramers–Moyal estimators
and sparse-identification-for-SDEs all emit point estimates, so "abstention
quality" is unscoreable for them and lagh wins by definition. That is a test
suite, not a competition.

The precedent is already in the repo: **LawSystemBench** — self-authored,
registered spec-before-generator, used to drive development, never a headline.
This follows it exactly, including the rule that its numbers are dev metrics.

## Curriculum

Levels 0–2 as proposed, with three changes, each for a measured reason.

**Level 0 — calibration.** OU (linear drift, constant diffusion), geometric
Brownian motion, and a deterministic ODE with observational noise only.
*Renamed from "exact recovery at σ = 0":* at σ = 0 an SDE is an ODE, so that rung
tests existing machinery rather than stochastic discovery. The headline claim is
the harder one — **at σ > 0, certify the drift's symbolic form with calibrated
intervals even though no pathwise law holds.**

**Level 1 — core SDE discovery.** Double-well with constant diffusion,
double-well with multiplicative diffusion, Van der Pol with additive noise, a
simple multiplicative-noise SDE. Targets: symbolic drift; symbolic diffusion or
a reasoned abstention on it; any almost-sure invariants.

**Level 2 — partial determinism.** Piecewise-deterministic Markov process
(deterministic flow + Poisson jumps), regime switching (two skeletons, stochastic
switching), and a deterministic system observed through a stochastic channel.
Targets: the continuous skeleton(s); jump/switching structure or an abstention
with a reason; **and the separation of process noise from measurement noise**.

**Level 3 — generator/Fokker–Planck-only identifiability.** Deferred.

### Added: nulls

The proposal has abstention *targets* but no task where **no law exists**. Every
lagh campaign carries nulls and they are what make an abstention rate meaningful:

* pure noise (no drift, no structure) — must certify nothing;
* a drift that is unidentifiable AT THE GIVEN Δt (below);
* a system whose trajectories do not visit enough of state space to determine
  the nonlinearity — the stochastic analogue of C0's single-solution refusal.

### Added: sampling rate as a first-class axis

Drift and diffusion have OPPOSITE data requirements, and this is the sharpest
identifiability structure in the whole suite: the drift needs a long horizon T,
the diffusion needs fine sampling Δt (quadratic variation is a Δt → 0 limit, and
the drift estimate carries an O(Δt) discretization bias). So there are regimes
where one is determined and the other is not, and an instrument with calibrated
intervals must show it. This is the direct analogue of the patch-resolution gate
in the PDE arc, and it makes abstention physically meaningful rather than
synthetic.

## Task format and scoring

Per system: several independent trajectories at several noise intensities, clean
and observationally-noisy versions, the ground-truth decomposition (skeleton,
noise law, invariants), fixed seeds, and a frozen checker accepting either
symbolic expressions with certificates or explicit abstention tokens with
structured reasons. State dimension ≤ 3.

Scored on: exact recovery, interval COVERAGE of the true coefficients,
abstention correctness, certificate validity under an independent checker, and
invariant discovery as a bonus. **Any confident-wrong dominates the ranking.**

## Registered predictions

- **S1.** OU: the drift certifies with intervals containing the truth wherever
  the martingale band is not vacuous, and the interval half-width scales as the
  CLT rate `σ/√(N·T)` — measurable across trajectory count and horizon.
  **MEASURED 2026-07-29 (`CASE_STUDY_STOCHASTIC_L0.md`): HALF FALSE.** The
  intervals contain the truth throughout. The width improves with the WINDOW
  LENGTH — the T half, consistent with `1/√L` though the sweep is single-seed so the
  exponent is indicative — and **does not improve at all** in the row count: 3.39 at
  144 rows, 3.39 at 576, on one fixed path set. The reason is the coverage factor
  itself, and it
  is structural rather than incidental: an exhaustive per-row band admits a
  coefficient only if every row accepts it, and `κ(n,δ)² = 2 ln n + 2 ln(2/δ)`
  grows exactly fast enough to cancel the √n averaging gain. The union bound's
  cost IS the averaging gain. A reach limit of the exhaustive doctrine under
  intrinsic noise, not a soundness problem.
- **S2.** The drift/diffusion asymmetry appears as stated: at fixed sample count,
  shortening Δt at fixed T tightens the DIFFUSION interval and not the drift's;
  extending T at fixed Δt tightens the DRIFT interval and not the diffusion's.
  This is the prediction that tests whether the intervals are calibrated for the
  right reason rather than merely being the right size.
- **S3.** Additive and multiplicative diffusion yield materially different
  certified diffusion forms on the two double-well tasks; where trajectories do
  not visit enough of state space, the diffusion abstains while the drift may
  still certify.
- **S4.** Every null certifies nothing, including the Δt-unidentifiable drift —
  which must abstain with a resolution reason, not a wide interval.
- **S5.** Zero confident-wrong across the suite.
- **S6.** Level 2 supplies what PDEBench structurally could not: replicates,
  known σ, and a known process-vs-measurement decomposition. The error-provenance
  characterizer (`lagh/errormodel.py`) separates them, scoring P1–P6 of
  `DIRECTION_ERROR_PROVENANCE.md` on more than two cases.
- **S7 (added 2026-07-29 with the κ answer).** Patch count is a strictly winning
  resource: the coverage factor costs √(2 ln n) while α's exponent gains linearly
  in the disjoint held-out count, so certifying over MORE patches tightens α
  faster than it loosens the band. **AMENDED after Level 0: measurable only where
  the drift ACCUMULATES.** On a stationary drift the martingale band exceeds the
  target's own range, so the per-row chance-match `q = min(1, 2ε/range(y))`
  saturates at 1 and α is vacuous — the significance bound says nothing in the
  regime where the drift is pure fluctuation. That is itself a result, and it puts
  S7's test on GBM rather than OU. Measurable by holding Δt and T fixed and
  varying the patch family alone — the band widens 17% from 300 to 10 000 rows
  while log α improves ~32× (the per-row chance-fit contribution loses ~3% to the
  wider band; the row count multiplies by 33). This is the opposite of the usual
  multiple-comparisons intuition and is the reason the regime is workable; if it
  fails, the independence discount (`n_disjoint`) is the thing that is wrong.

## Why this is worth doing now

`DIRECTION_ERROR_PROVENANCE.md` is currently scoreable on exactly two cases —
PDEBench advection (pipeline error) and this program's own C1 fields
(observational). Its central question, separating process noise from measurement
noise, has no ground truth anywhere in the repo. Level 2 constructs that ground
truth by design. That alone justifies Levels 0–2 independently of the benchmark
framing.

## Sequencing

0. **Partial determination as a first-class verdict — STARTED 2026-07-29.**
   `Certificate.partial` and `certify.invariant_content` are in: on a structural
   abstain the verdict now carries what EVERY law consistent with the data at
   the declared band agrees on, per coefficient, with the terms that are required
   in all of them, excluded from all, and the tightest/loosest spans. Measured on
   PDEBench CFD, which previously reported nothing: every true coefficient lies
   inside its recovered range, several are pinned to spans of 1e-6, and the
   genuinely unconstrained ones are named. This is a prerequisite rather than
   cleanup — the per-component scoring registered above cannot be expressed by
   `certified: bool`, so the frozen checker's interface depends on it.

   **THREE of five dimensions are now on the shared vocabulary** (2026-07-29):
   STRUCTURE (`invariant_content` over the certifying set), INTERVAL (the
   certified law's own ranges), and — added after the PDEBench re-run — MODE:
   `StateCertificate.partial` is a `determination(status="state")` record keyed
   by mode label, and a test asserts it cannot drift from the `modes` dict it
   duplicates (kind, bounds and `resolved` component-wise). The mode dimension
   was cheap precisely because `determination`'s vocabulary was lifted from
   C4's state certificates in the first place; this closes the loop.

   Remaining, and neither is code-shaped:

   * **DOMAIN — DECIDED AND BUILT 2026-07-29.** A domain restriction ("this
     holds where the conductivity is in its high phase") qualifies EVERY
     component at once rather than being one among them, so it does not fit the
     `(name, lo, hi)` entry shape. It is now a QUALIFIER on the record
     (`certify.domain_qualifier`), and the rule that gives it teeth is in
     `conjoin_determination`: records established on DIFFERENT domains REFUSE to
     conjoin, because the conjunction is defined only where both were
     established and this layer cannot intersect two predicates — guessing that
     two differently-worded regions are the same set is the unsound direction.
     Same-domain records intersect component-wise, and an EMPTY intersection is
     reported as a contradiction rather than dropped. The decision arrived via
     `DIRECTION_PDE.md` (c): the honest domain restriction was chosen as the
     variable-coefficient route, and this qualifier is what it costs. Darcy is
     the live case (β = 0.1000 in one phase, ~49% off in the other) and is now
     also the first PRODUCER: both phases were re-run with
     `discover_equation(..., qualifier=)`, and `run_darcy_domains.py` checks the
     refusal on that output rather than on constructed records. The interface
     step 1 freezes has therefore been exercised end-to-end by something real,
     which was the point of doing it before the freeze rather than after.
   * **STRATUM is a campaign convention** (certifiable conservation vs
     conjectured closure, per the traffic plan) rather than a per-verdict
     record, and is the lowest priority of the five.
1. Register the certificate-kind decision (above) and the checker's interface.
   **DONE 2026-07-29 — `docs/STOCHASTIC_CHECKER.md`, `lagh/stochcheck.py`,
   `certify.coverage_factor`/`coverage_budget`, 32 tests.** Both halves:
   * The certificate-kind decision is registered above, and the κ question it
     left open is answered above (§ *ANSWERED*) and in full in
     `STOCHASTIC_CHECKER.md` §1. The martingale inequality is itself measured in
     simulation rather than trusted.
   * The checker interface is **frozen** against `certify.determination`, with
     `test_the_interface_is_frozen` as the lock — every vocabulary tuple and
     dataclass field list asserted verbatim. Component claims travel only in a
     determination record, so lagh's own verdicts are submissions by
     construction; five outcomes with `confident-wrong` first in the ranking key;
     silence is not abstention; a procedural refusal voids credit but not
     exposure; domain-qualified answers are scored on their own domain.
   * One hole found and closed while freezing: a key rewarding COVERAGE alone can
     be topped by submitting [-1e9, 1e9] everywhere — covered on every component
     with zero confident-wrongs. `certify.vacuous`'s doctrine now applies per
     component (an interval must be tighter than 2·TAU of the task's coefficient
     scale to count as a determination), and determinations rank ahead of
     coverage. A vacuous answer is not wrong, just worth nothing.
   The constraint carried in from the PDEBench correction — **a declared error
   belongs to a CONSUMER, not to a dataset** (`STRATEGY.md`,
   `DIRECTION_ERROR_PROVENANCE.md`) — is discharged as validation: three
   consumers (`drift-band`, `diffusion-qv`, `observation`), each stating the
   quantity and units it bounds, and a submission that gives one magnitude two
   jobs is refused. The declaration audit reports declared/true per consumer and
   flags a decade either way, so the 3900× over-declaration that took a session
   to find is now a column in the score sheet.
2. Level 0 — three systems, mostly existing machinery, calibrate κ.
   **DONE 2026-07-29 — `docs/CASE_STUDY_STOCHASTIC_L0.md`, `lagh/ito.py`,
   `experiments/stochastic/`.** The deliverable landed: the band holds on every row
   in 0/40 replicate failures at every declared δ from 0.5 to 0.01, conservatively
   (the worst row across 40 runs reached 0.97 of its band at δ = 0.5). Five headline
   results beyond it, three of them corrections:
   * **The generator's f-family is not optional.** With `f = x` alone a stationary
     drift is vacuous at EVERY window size, because `E[a(X)] = 0` under any
     stationary law. `f = x²/2` makes it accumulate, and the identifiability
     condition is **θ·L > 2κ(n,δ)²**.

     **AMENDED 2026-07-30 — the f-family has a price, and it is only visible on a
     real instrument** (`CASE_STUDY_TWEEZERS_C1.md` §3). A nonlinear f imports the
     quadratic variation: `d(x²/2) = x dx + ½d[x]`, so with `f = x²/2` the identity
     reads `θ ≈ b²/(2⟨x²⟩)` for a stationary path — **the drift is determined
     through b²**. Every real instrument anti-alias-filters before storage, which
     destroys quadratic variation, and that bias lands on the drift. Measured on a
     C-Trap bead: b² 0.554 of truth, weak-form θ 0.539 — the same factor.

     **The trade is structural, not a tuning failure.** For a stationary scalar
     diffusion a single-time Itô weak form either fails to accumulate (f linear ⇒
     `E[f'a] = E[a] = 0`) or imports `d[X]` (f nonlinear ⇒ `f'' ≠ 0`). There is no
     b²-free accumulating f. So the fix leaves the single-time form:
     **`ito.build_lag_rows`** — state-weighted Itô increments at a lag,
     `y = Σφ w(X_i)(X_{i+h} − X_i)`, whose left-endpoint sum IS the Itô integral and
     forms no quadratic-variation term at all, while `w = x` supplies the
     accumulation `f = x²/2` was adopted for.

     Measured (`run_lagform.py`, an OU matched to the bead, 23 kHz anti-alias
     stand-in reproducing the instrument's b² attenuation at 0.539 against the real
     0.554):

     | | θ / truth |
     |---|---|
     | single-time weak form, filtered | **0.533** |
     | lagged form, lag 32, inverted | **0.981** |
     | the real bead, Force 1y, lag 16, inverted | **0.980** (single-time: 0.341) |

     The lag turns an UNKNOWN instrument bias into a KNOWN estimator bias:
     `θ_apparent(h) = (1 − e^{−θh})/h` exactly for a linear drift, which inverts to
     within 1% across the ladder, while the filter's own effect
     (filtered/clean) decays to under 2% from lag 16 on. Force 1x stays broken at
     0.448, correctly — C1's ACF gate had already flagged it as carrying a
     low-frequency contaminant, and a better estimator does not rescue a corrupted
     axis.

     **The h-columns were built, and they fix the ESTIMATE but not the
     CERTIFICATE.** `build_lag_rows(bias_names=, bias_orders=, generator_max=)` puts
     the conditional-mean bias into the design instead of bounding it in `quad`,
     with order k carrying `h^k/(k+1)!` — the exact expansion of `(1−e^{−z})/z` —
     and rescaled by `rate^k` so a search does not have to span θ and θ² at once.
     Measured (`run_lagform.py`):

     | bias orders | clean | filtered |
     |---|---|---|
     | none | 0.743 | 0.731 |
     | 1 | 0.946 | 0.913 |
     | 2 | **0.989** | 0.933 |
     | 3 | 0.990 | 0.928 |

     Order 3 = order 2, so the expansion has CONVERGED and what remains on the
     filtered arm is the measurement channel rather than the estimator.

     **It still does not certify, and the reason is structural rather than
     statistical.** Signal-to-band is not the obstacle: it scales as √L exactly as
     Level 0's window-length finding predicts (0.96 → 2.27 → 4.07 for 5× and 15×
     windows) and reaches **6.97** at 14 trajectories × 9 s — where the joint
     admissible interval is still unbounded. The cause is DEGENERACY. With `w = x`,
     the drift column `∫φx²dt` and every bias column are the SAME integral times a
     scalar function of h, so the columns are near-parallel and the admissible set is
     a long thin ridge, unbounded in the compensating direction however good the band
     is. More data does not help; it is a rank problem, not a noise problem.

     **And the ridge exists because the fit throws away a constraint it has.** The
     bias coefficients are not free: `La = a a' + ½b²a''` is DETERMINED by the drift
     (and the diffusion). Fitting them as independent parameters discards that and
     manufactures the degeneracy. Two honest routes out, neither of them more columns:

     * **Constrain the bias coefficients to the drift's** — the correct model, but
       nonlinear in the coefficients, so it leaves `certify_drift`'s linear contract
       and needs its own band treatment.
     * **Invert the h-dependence outside the fit**, which is what already works:
       `θ = −ln(1 − θ_app h)/h` recovers **0.980 on the real bead** (single-time
       form: 0.341) and 0.981 in simulation. It is a measured point estimate with a
       measured operating window, not a certificate.

     So the b²-free drift is, at this stage, a **scored estimator and not a certified
     one**, and the next lever is neither columns nor data: signal-to-band is
     `∫φ w a dt / κ√(∫φ²w²b²dt)`, a ratio maximisable over `(φ, w)` that nothing has
     yet been chosen to optimise — `w = x` came from wanting accumulation and the
     bump φ came from the PDE arc. Failing that, the registered fallback is option 3
     of the certificate-kind decision (a distinct statistical-certificate class),
     which two independent estimators have now hit the same wall without.
   * **S1 is half false** (above), and the coverage factor is why.
   * **One confident-wrong, found by the frozen checker on its first real
     producer, and fixed.** Realized quadratic variation of a noisily-observed path
     estimates `[X] + 2nσ_obs²`: conservative in the band's scale, a systematic
     offset in the Itô correction on the target. One measured quantity, two
     consumers, safe for one and unsafe for the other — the declaration rule
     arriving in a new form. Debias + refuse the rows the declaration explains +
     build the observational Gram (which had been a silent no-op).
   * **Two defects in code that predates the run**: `invariant_content`'s claim
     overstated what it computes (a range over the certifying set the SEARCH found
     is not a bound over the laws that exist — measured excluding a law that
     certifies), now corrected and joined by `certify.admissible_interval`, an LP
     bound that does hold; and passing `sigma=0` with a martingale band told the
     engine the data was clean, producing a false exact rational for a
     diffusion-scale coefficient and a narrowed rival search.
   * **The registered headline is NOT met at Level 0.** The drift's symbolic FORM
     did not certify on any Level 0 system with the declared 4-term library. The
     abstains are honest — the LP proves the admissible set contains materially
     different laws — and the output is partial determination with sound covering
     bounds. Read as reach short of the claim, not the claim met.
3. The Itô/generator weak form in `lagh/weakform.py` (new terms, same discipline).
   **DONE 2026-07-29 for the vocabulary; the assembler migration is deliberately
   deferred.** What landed, all of it opt-in so no campaign predating it changes:
   * **`Term.measure`.** `"dt"` is every deterministic term; `"d[<field>]"`
     integrates against that field's REALIZED QUADRATIC VARIATION. That single
     addition is what lets the Itô correction be a library term — because
     `b² dt IS d[u]`, so the term can be written down without modelling the
     diffusion at all. `lagh.ito.ito_terms(f, library)` emits the whole identity
     (target, correction, drift columns) in that vocabulary.
   * **The martingale scale is NOT a term, and cannot be.** `∫φ²f'²d[u]` is
     QUADRATIC in the test function while the weak form is linear in it. It belongs
     where `noise_l2` belongs — a per-patch scale a coverage factor multiplies —
     and `build_nd(martingale=(field, gexpr))` measures it into `WeakSystem.qv`.
     The only difference from `noise_l2` is that one multiplies a DECLARED sigma
     and this one is measured, which is the whole of §1c.
   * **`build_nd(rough=True)`** relaxes the three deterministic resolution gates,
     each for a stated reason: a Brownian path is aliased by construction, its
     quadrature converges at O(h) so the observed ladder order sits near 1, and a
     bound relative to the term's own scale is the wrong comparison when the band
     is made of the martingale. The bound is declared and the CALLER gates on it,
     because only the caller knows the coverage factor.
   * **A `d[]` term does not refine.** Subsampling a realized-QV estimator computes
     a different estimator rather than resolving the same integral, so its declared
     bound is its own statistical spread and its ladder order is reported as NaN.
   Validated to machine precision against the hand-rolled assembler
   (`test_the_term_vocabulary_reproduces_the_hand_rolled_rows`, rel ≤ 1e-12 on the
   target, every column, the martingale scale and the correction's bound), and all
   55 weakform-dependent tests are unchanged.

   **What is deferred and why.** `ito.build_rows` keeps its own scalar assembler;
   the equality test is the bridge that stops the two drifting. Delegating it onto
   `build_nd` buys nothing until MULTI-FIELD state is needed — which is exactly what
   Level 1 (Van der Pol, 2-D) and Level 2 (up to 3-D) require, and `build_nd`
   already has that geometry. So the migration is sequenced as the first act of
   Level 1, driven by a requirement rather than by tidiness.
4. Level 1. **FIRST INCREMENT DONE 2026-07-29 —
   `docs/CASE_STUDY_STOCHASTIC_L1.md`.** The diffusion became a CLAIM:
   `build_rows(diff_names=...)` puts b² in the design matrix as ordinary dt columns
   (`½∫φf''h_j dt`), so drift and diffusion are identified jointly and the measured
   Itô correction — the one unsafe consumer of realized quadratic variation —
   disappears. `certify.admissible_functional` bounds any linear functional of the
   coefficients, i.e. the drift AS A FUNCTION. Three findings:
   * **The drift and the diffusion want OPPOSITE noise intensities.** The drift's
     signal-to-band goes as 1/b and the diffusion's as b, because the thing the
     diffusion measures IS the noise. Measured on GBM across b ∈ [0.02, 2]: the two
     relative widths move monotonically in opposite directions and neither is
     determined at the crossover. A structural tension between two Level 1 targets.
   * **S2 was registered in the wrong picture and should be re-registered before it
     is scored.** The Δt/T asymmetry is real but lives in WHICH ESTIMATOR is used:
     the weak form for the drift, realized quadratic variation for the diffusion
     (relative error √(2/m)). In the weak form the diffusion barely depends on Δt.
     So the next increment certifies b² from quadratic variation — regressing
     b²(x) on the local realized QV, which the `measure="d[u]"` terms already
     express — rather than from the design matrix.
   * **A bistable system hides its own drift**: the wells sit where the drift
     vanishes, so the process spends its time exactly where the drift carries least
     information. Confirmed as a statement about the drift FUNCTION, not just its
     coefficients. And multiplicative noise vanishing at the origin DISCONNECTS the
     state space (zero well-crossings measured), making the registered third null a
     property of a system rather than a null.
   **SECOND INCREMENT DONE — the diffusion certified from quadratic variation.**
   `build_qv_rows` / `certify_diffusion`: the same weak form with the target
   `∫φ w(X) d[X] = Σ_j d_j ∫φ w h_j dt`, the band's scale again MEASURED (from the
   fourth moment of the increments, since `Var((ΔX)²) = (2/3)E[(ΔX)⁴]`).
   * **2600× tighter, and OU's diffusion CERTIFIES** — b² = 1.9603 against a truth
     of 1.96, joint bound [1.855, 2.097], resolved. GBM's x² coefficient resolves at
     0.21 relative width against 532 inside the drift's design matrix. Both routes
     cover; only this one is worth having. **The first certified diffusion in the
     arc.**
   * **The coupling runs both ways.** The drift's band needed the diffusion via
     ⟨M⟩; this estimator needs the DRIFT, because `E[(ΔX)²] = b²Δt + a²Δt²`. That
     O(Δt) drift leakage is a deterministic band term against a declared bound on
     |a| (0.3% of the band on GBM, 22% on OU at `drift_max = 5`).
   * **A caveat stated rather than assumed:** these increments have CHI-SQUARE
     tails, so the applicable inequality is Bernstein's, not the continuous-
     martingale bound. The sub-exponential correction is O(√(Δt/L)) and
     `QvBand.bernstein_correction` reports its measured size.
   * **A claim I asserted and the measurement changed:** the w-family is what makes
     b²'s state dependence identifiable only for a STATIONARY process (OU: bounds
     tighten 2.9×). Where the process grows, the windows already sample different
     state regions and w adds nothing (GBM: marginally wider).
   **THIRD INCREMENT — MULTI-FIELD ITÔ, and a scored Level 1 table.** S5 holds:
   zero confident-wrong over four systems, all 28 component bounds covering, two
   diffusions CERTIFIED. Then the 2-D machinery: `Term.measure` gained
   cross-variation `d[u,v]`, `build_nd(martingale=)` takes a per-field sensitivity
   list (⟨M⟩ collapses to the square of the SUMMED increment — no d² loop), and
   `ito.ito_terms_nd` emits the multi-dimensional identity. Validated on Van der Pol:
   the truth is inside the band on the deterministic and the stochastic row alike.
   * **A component's quadratic variation says whether it carries noise.** Van der
     Pol's `x` is noise-free, so its weak form is a deterministic identity — residual
     1.8e-4 against a martingale band of 0.35, over-declared 2000×, because a
     differentiable path's realized QV is the O(Δt) residue rather than zero. Halving
     Δt halves it while a driven component's stays flat, so **the diagnostic is
     measurable**: this is the multi-dimensional form of the error-provenance
     question. Conservative, never unsound.
   * **The noiseless-component band needed NO declaration** — the fifth instance of
     the recurring pattern turned out to be measurable. The martingale part of a
     realized QV separates from the smooth residue by stride scaling
     (`Σ(u[i+s]−u[i])² ≈ α·s + β·s²`, α accurate to 0.3%), so
     `build_nd(martingale_decompose=True)` tightens a noise-free row 12.5× with
     coverage intact. Safety is intrinsic: it returns α + κ·se(α) and falls back to
     the total on a poor fit.
   **FOURTH/FIFTH INCREMENTS also gave the arc's first certified DRIFT.**
   `ito.build_rows_nd` packs multi-field rows into the `ItoRows` `certify_drift`
   already consumes, so the 2-D case inherits every gate. On Van der Pol
   **`dx = y dt` certifies** at coefficient 1.0000, bound [0.970, 1.030], α ≤ 1e-15,
   over a 6-term bivariate library — while the DRIVEN equation is vacuous at the same
   configuration (median signal/band 32 against 0.12). **In a partially-driven system
   the noise-free components' equations are far more determinable.**
   **SIXTH INCREMENT — Van der Pol SCORED**, and the first producer for the frozen
   interface's component INDEX: a 2-D system has one drift equation per component, so
   `part[index]:term` plus several submissions per task (both in the interface since
   the freeze) is what expresses it. Zero confident-wrong across 18 components, all
   covering, 6 informative, 2 resolved — `drift[0]:y` at [0.9827, 1.0177] with the
   other five terms of that equation excluded to ±0.05, so its FORM is determined, and
   `diffusion[1]:1` resolved at [0.123, 0.369] against 0.25. `diffusion[0]` abstains on
   vacuity, correctly: the x component carries no noise, so zero and small-nonzero are
   genuinely indistinguishable at that band and vacuity is the honest form of "there is
   none".
   Still open in Level 1: the DRIVEN Van der Pol equation, vacuous at every window
   length tried (0.49 at L = 4 the least bad), which suggests it needs a different f
   family rather than a different window; the invariants target; and a state-dependent
   b² whose FORM certifies.
5. Level 2, as the error-provenance testbed. **SEPARATION INCREMENT DONE
   2026-07-30** — scoped to the process-vs-measurement question rather than the full
   PDMP/regime-switching curriculum, on the reasoning that this is the capability the
   intended real-data read (highD) will demand first.

   **The three noise sources have three distinct STRIDE exponents, so one polynomial
   fit separates all three:**

       Σ_i (u[i+s] − u[i])²  ≈  c + α·s + β·s²

   observation noise is iid so `E[(e_{i+s}−e_i)²] = 2σ²` at EVERY lag (constant in
   s); process noise is a martingale, whose increment variance grows linearly; a
   smooth drift's squared increment grows as s². So **σ_obs is MEASURED, not
   declared** — retiring at its root the declaration that produced this arc's only
   confident-wrong. `weakform.qv_three_way`, `run_level2_separation.py`, 45 cases:

   | | result |
   |---|---|
   | σ_obs, separable regime (c ≥ 0.1·α) | median rel err **0.32%**, 14/15 within 10% |
   | σ_obs, buried | **6 cases REFUSED** — worst error had they been used: **301%** |
   | process b² | median **1.2%**, **24/24** within 10% |
   | null: no observation noise | worst spurious σ_obs 2.1e-3 |
   | null: no process noise | worst spurious b² 1.8e-3 |

   **The refusal bar is the necessary half, not a caveat.** The buried cases fail by
   OVER-estimating σ_obs, which subtracts too much from the quadratic variation and
   makes a debiased band too tight — the impostor-admitting direction. The boundary
   (`SEP_MIN_FRAC = 0.1`) is measured, not chosen.

   This answers `DIRECTION_ERROR_PROVENANCE.md`'s central question, which that
   document records as having ground truth NOWHERE in the repo. It now has 45 cases.

   **Registered before the highD read**: the mechanism assumes iid observation error,
   which is why the term is constant in s. Computer-vision tracking error is
   AUTOCORRELATED across frames, which breaks that signature and will bias `c`. It is
   testable directly — if the error is correlated to lag τ, the constant term appears
   only for strides beyond τ, so restricting the fit to s > τ should recover it. That
   prediction is registered here rather than discovered during the read.

   Still open in Level 2: the PDMP and regime-switching systems, and the jump/switching
   structure targets.

Level 3, SPDEs, higher dimension, partial observation and non-Markovian noise are
explicitly out of the minimal suite.
