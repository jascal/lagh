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
   `certified: bool`, so the frozen checker's interface depends on it. Remaining:
   the same treatment for the other four dimensions (interval / mode / domain /
   stratum) under ONE vocabulary, so verdicts compose.
1. Register the certificate-kind decision (above) and the checker's interface.
2. Level 0 — three systems, mostly existing machinery, calibrate κ.
3. The Itô/generator weak form in `lagh/weakform.py` (new terms, same discipline).
4. Level 1.
5. Level 2, as the error-provenance testbed.

Level 3, SPDEs, higher dimension, partial observation and non-Markovian noise are
explicitly out of the minimal suite.
