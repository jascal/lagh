# Certification boundaries — BND1 registration

Registered 2026-10-07 before executing the bank. Baseline: `049238d`.
This is a finite engineering/research study, not a new universal soundness proof.
The objective is to give the noisy, loose-floor and multidimensional pre-pass
paths explicit, tested claims, and to change gates only after a measured witness.

## Bank and fixed budget

`experiments/run_certification_boundaries.py` defines 120 discovery cases:

- Four one-dimensional laws: affine, rational, fractional power, and a linear
  law with a small quadratic contribution.
- Four three-dimensional laws: monomial, monomial with a small extra term,
  affine, and affine on the unit sphere.
- Three declared linear libraries: independent, nearly collinear, and exactly
  collinear columns. These exercise the existing interval-parameter path.
- One independent-random-target negative control.

Each has seeds 710 and 711 and five regimes: clean, relative Gaussian noise
at 1e-6 and 1e-3, deterministic bounded error under absolute floors 1e-6 and
1e-3. Exactly 240 rows, split 60/20/20; max tier 3 (declared linear: tier 1).
Inputs are in [.5,3], except sphere inputs normalized from Gaussian vectors.
The floor error is a fixed sinusoid with amplitude floor/4. Noise realizations
are not clipped: truth-band failures are recorded and coverage is conditional
on the truth satisfying the supplied finite-row band.

Each discovery has a 20-second elapsed-time limit, single-threaded BLAS.
Timeouts and errors remain in the denominator, distinct from refusals. Four
independent worker processes maximum. No retry at a larger budget, no seed
replacement. Scoring uses 512 new original-region inputs and 512 extended-region
inputs [.25,6], on the same variety when constrained. These are diagnostics,
never new certification or selection data.

## Measurements (reported separately)

1. Finite-row band validity of the truth and returned law; fresh in-region
   and extended-region error (the latter does not change the certified domain).
2. Structural and coefficient agreement with the known generating expression:
   equality of monomial supports after cancellation/expansion where polynomial;
   otherwise symbolic equivalence or normalized numerical discrepancy >1e-8.
   Numerical agreement alone is explicitly weaker than exact-form equality.
   A mismatch is scored as a false *exactness claim* only if output claims
   exact identification; a scoped consistency result still retains the mismatch.
3. Returned parameter intervals, plus diagnostic coordinate slices computed
   with all other parameters fixed. Endpoint feasibility, truth inclusion where
   coordinates can be matched, and simultaneous-corner feasibility are separate.
   Coordinate slices are neither joint boxes nor marginal projections of a
   joint feasible set. No confidence/coverage guarantee is inferred from them.
4. Certificates, refusals by reason, timeouts, errors, proposal count, route,
   discovery runtime and scoring runtime. Report rates over ALL registered rows.

## Mechanism witnesses and predictions

Run a separate, labeled diagnostic battery (not credited as natural discovery):

- Correlated linear coefficients: coordinate-slice endpoints each fit while a
  simultaneous corner fails. Prediction: slices cannot license a joint box.
- Offset fitted coefficients on correlated columns: a feasible true coefficient
  vector can lie outside a coordinate slice with other coefficients fixed.
- A controlled near-collinear candidate in the 3-D pre-pass, with candidate
  generation replaced only for this witness: marginal gate passes, joint gate
  fails. Prediction: baseline pre-pass does not honor the configured joint gate.
- Noise/loose-floor observations compatible with two generating coefficients:
  prediction: an exact generating value cannot be established by these data.

Natural-bank predictions: some noisy/floor certificates differ from the exact
generator; null controls never certify; clean ordinary affine/power cases retain
reach; interval records require conditional wording. These are scored predictions,
not permission to erase failures or silently broaden error declarations.

## Decision and stopping rules

Freeze this registration and runner in git before any scored run. Preserve all
baseline rows, software versions, source hashes and input hashes. Review the
entire bank, including failed predictions, before changing production behavior.

One repair round is allowed: (a) route an existing gate through an uncovered
path only after a rejecting witness; (b) correct overstated exactness/interval
claims and unsafe interval handling without suppressing supported consistency
results. Do not apply the clean joint test blindly to noisy data. Do not widen
grammar, relax epsilon or revise the bank to purchase successes.

Before repair scoring, register the concrete changes and expected effects.
Replay the same bank and mechanism witnesses on the candidate. For fresh
confirmation run the same bank on seeds 810 and 811, once. Thus at most 360
natural discovery runs, plus targeted regressions and the existing suite.
If a path still lacks justified exactness, expose finite-data consistency or
refuse the stronger claim, with the limitation carried through public payloads.
Unresolved scientific identifiability is a valid measured boundary; silently
claiming it is solved is not. Completion requires an outcome report, a route
claim matrix, retained evidence, meaningful regression tests and required checks.

## BND1 baseline and repair registration (before production edits)

Baseline at registration commit `4701e5b`: all 120 rows completed, 40 certified,
54 refused, 26 timed out, zero exceptions or missing diagnostic scores. Ten
certified expressions mismatch the generator (9 noisy, 1 loose floor), while
every returned law fits all supplied observations. All true laws fit their
declared bands. No null target certified. Clean rational and perturbed 3-D
monomial cases time out: this bounded study does not establish their reach.

The interval diagnostic returns infeasible endpoints on 16 certified cases.
`parameter_interval` initializes its near endpoint at half the first rejected
step without checking it. The controlled pre-pass passes both gate selections
although the same snapped candidate fails `joint_pinned`. Coordinate slices
have feasible individual endpoints but an infeasible simultaneous corner, and
can exclude a feasible true vector. Additional source-directed witness:
`pdesystem.intervals_for` reports [1,1] for both coefficients of x0+x1 under a
.01 band although 1.001*x0+.999*x1 fits. Repeated numeric atoms also alias
different coefficient positions. The new regression file fails three tests
and passes the slice-semantics witness on the baseline.

One repair round is now registered:

1. Apply the configured clean coefficient gate to the 3-D pre-pass, including
   quotient constraints and their serialized domain restriction. Preserve the
   noisy pre-pass policy; do not substitute a clean pinning test for uncertainty.
2. Carry an explicit claim kind from engine certificates through public
   recover (active/passive/tiny), verify, PDE and submission payloads. Noisy,
   loose-floor and declared-error results claim finite-row band consistency,
   not exact structure or coefficients. Clean pinning is an operational gate
   result, not proof of the generating law. Preserve underlying certification
   metrics: a generator mismatch remains a mismatch after a claim correction.
3. Start interval bisection from a verified feasible point. Keep conditional
   slice provenance explicit, never call slices a simultaneous box or marginal
   coverage. Do not replace all coefficients by independently chosen slice
   midpoints without rechecking the resulting law. Give unit/repeated linear
   coefficients independent positions in the PDE adapter's interval diagnostic.

Predictions: the pre-pass rejecting witness is repaired while the explicit
legacy marginal arm remains reproducible; every reported slice endpoint fits;
PDE unit/repeated coefficients have finite nonzero-width slices; the natural
bank has no null certificates or new observed-row failures; clean baseline
certificates retain their laws on cases completing in both arms. Runtime,
timeouts, natural mismatches and refusal differences are reported without a
post-hoc speed or reach win criterion. Noise/floor mismatch certificates must
carry the consistency claim. Unresolved paths retain explicit limitations.

The fixed bank replay and fresh 810/811 confirmation remain the last two
natural discovery arms. Targeted wrapper, callable-band and interval fixtures
are regression checks, not additional discovery success counts. Any failure of
these predictions is retained and explained; no second capability campaign is
authorized by this registration.
