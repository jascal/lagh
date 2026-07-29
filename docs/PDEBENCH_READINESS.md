# PDEBench readiness

> **Superseded in one respect (2026-07-29): PDEBench is a DEV target, not a
> blind read.** The pre-flight below is what the first pass ran and it stands as
> written; what did not survive contact is the assumption that a run against it
> could be scored. See `CASE_STUDY_PDEBENCH.md` § "Why this is dev and not a
> read" — the short version is that the band's dominant term has to be chosen
> per file, choosing it requires the answer, every family costs a convention
> hunt, and the benchmark's own leaderboard scores forecasting surrogates, so
> there is no frozen number to be blind against.

**Written 2026-07-28, after C3 (systems) and C4 (state certificates).** This is
the pre-flight list for pointing the weak-form instrument at PDEBench, and the
statement of what it can and cannot claim there. It is deliberately written
BEFORE any file is downloaded, so that the error model is a registration rather
than a reaction to results.

## Retrospective: this was a RECONNAISSANCE pass and did not know it

Added 2026-07-29, after the pass. The three declarations below are about
MECHANICS — how to load, what to declare, what to refuse — and all three earned
their place. What this document never registered was a prediction the first run
could FALSIFY, and that omission shaped the whole session: the truth check firing
at 1185× arrived as a surprise to be interpreted rather than as a clean
refutation, so the DEV reclassification, the pipeline-decode reframing and the
"report the required declaration" device for CFD were all discovered mid-run.

The line that belonged here, in the pre-flight, before any file was opened:

> **The shipped fields' deviation from their own stated laws will be at or below
> storage precision (σ_rep ≈ 6e-8 relative).** If it is not, this target is model
> output whose generating error dominates the band, the exercise is pipeline
> decode rather than discovery, and no score is available from it — in which case
> the pass proceeds as dev, and its deliverables are defects found and
> capabilities exercised, not results.

Measured: **8.9e-4 to 2.0e-2 relative**, four orders above storage precision.
False in the first hour, and every later conclusion follows from it.

`STRATEGY.md` now carries reconnaissance as a category, with the gating rule for
engine changes found during one.

## What is ready

| capability | state | where |
|---|---|---|
| scalar weak-form laws, clean | certified C0 | `experiments/pde/run_c0.py` |
| declared-noise ladder, interval parameters | certified C1/C1b | `run_c1.py` |
| forward-integration verify | certified C2 | `verify.py`, `run_c2.py` |
| **systems** (multi-field terms, per-equation targets, conjoined α) | certified C3 | `lagh/pdesystem.py`, `run_c3_systems.py` |
| **2-D space** (x, y, t geometry) | certified C3 stage 4 | `lagh/weakform.py::build_nd` |
| **state certificates** (initial-condition inversion) | certified C4 | `lagh/statecert.py`, `run_c4_state.py` |
| declared FIELD error from a reference solver | measured C3 | `systems_fields.py::solve_declared` |
| vector-RHS verify for systems | built, unexercised on PDEBench | `verify.py::verify_system` |
| loader + declared error model for external files | built, tested on synthetic | `experiments/pde/pdebench.py` |

The PDEBench 1-D sets (advection, Burgers, diffusion–reaction, CFD) map onto
stages the arc has already certified; the 2-D sets (diffusion–reaction, shallow
water, incompressible NS) map onto stage 4's geometry.

## The three declarations a PDEBench run must make

**1. Storage quantization — computed, not assumed.** PDEBench ships float32
HDF5. A float32 value is not the number the simulation produced; it is that
number rounded, and the rounding is a real error of size `eps32·|u|/2 ≈
6e-8·|u|`. This is the float32 lesson from the LLM-SRBench blind read (declare
σ_rep) and it applies verbatim. `pdebench.sigma_rep_for` computes it; because
the weak-form band takes one scalar σ per field, the honest scalar is the bound
at the field's own peak — an over-declaration everywhere else, which is the safe
direction.

**2. Solver error — DECLARED, and it cannot be measured.** Every PDEBench field
is a solver's output, and the coarse levels a tolerance ladder would need were
never saved. C3 measures this term for its own fields (tolerance ladder plus
spatial-resolution ladder) precisely because PDEBench will not let it be
measured. So a PDEBench certificate is conditional on a stated `field_err`, the
statement travels with the result (`pdebench.declared_noise` returns
`field_err_is_measured: False`), and the sensitivity of the verdict to that
number is part of the report, not a footnote. **A run that sets it to zero is
claiming the data is exact and must say so in those words.**

> **Amended 2026-07-29, from the run this declaration governed.** Two things were
> missing here and both bit. (a) Solver error CAN be measured where an exact
> solution exists — advection has one, and `pdebench_solver_error.py` measured
> 2.75e-2. (b) Having measured it, the run declared it, and **that was the wrong
> quantity**: a pointwise deviation accumulated over a trajectory is not a local
> weak-form residual over a patch, and it over-declared by ~3900× (the weak form
> requires 7.06e-6 on that file). The declaration this section should have
> demanded is therefore sharper: **state the field error, AND state which
> quantity it is a bound on — pointwise or weak-form, local or accumulated —
> because only the second is what the band consumes.** The direction of the
> error was safe; four orders of needless width was not.

**3. Coordinate storage — regularized, and the deviation reported.** Found in
the dry run, and it would have been a real error: PDEBench stores its
COORDINATE vectors in float32 as well as its fields. A perfectly uniform
101-point t-axis over [0, 1] therefore arrives with ~6e-8 absolute rounding in
it, which is **6e-6 relative to the step** — and that does not cancel out of the
relation, because the by-parts weights carry `1/step^k` with k differing per
term, so it lands on the u_t column at 6e-6 of its own size, a hundred times the
float32 FIELD noise. The loader rebuilds each axis as an exact linspace through
its endpoints and reports the deviation from what was stored
(`coord_deviation`); nothing is invented, and the deviation is exactly what a
reader needs to check the claim.

**4. Grid geometry — checked, because a wrong L is silent.** PDEBench stores
`(n_samples, n_t, n_x[, n_y])` with periodic space; the axis order, the domain
length and the endpoint convention all have to be converted, and a wrong L
rescales every derivative term without any error message — a certificate would
happily certify the rescaled law. `pdebench.check_geometry` reports uniform
spacing per axis, the span, and whether the periodic endpoint is duplicated —
the last from the FIELD (`u(x_first) == u(x_last)`), because the coordinates
cannot tell you the intended domain length — and refuses rather than guessing.
Our own stage-1 fields tripped exactly this trap during C3's verify track: an
endpoint-inclusive grid put 18244 of 41634 forecast points outside the band with
the law exact to machine precision.

## What a PDEBench certificate may and may not claim

- It is a claim over the STATED PATCH FAMILY on the sampled solutions, at the
  declared σ, conditional on the declared solver error.
- Multi-solution holdout is non-negotiable: PDEBench ships many samples per PDE,
  so certifying on patches from a HELD-OUT sample is available and required. A
  single-sample claim is an on-shell statement about that sample (C0 N5, C3 Y3).
- Where the recovered coefficients are the benchmark's own generating parameters
  (β, ν, and the CFL/η settings recorded in the filenames), agreement is a
  measurable check and disagreement beyond the reported interval is a finding
  about the data or the declaration, not a tuning target.
- **Nothing about "beating" a baseline.** PDEBench's own leaderboard is a
  forecasting-error benchmark for surrogates; this instrument answers a
  different question (does a law certify, and over what domain), so the honest
  comparison is against PDE-FIND/WSINDy-style recovery on the same files, with
  abstentions reported as abstentions.

## Pre-flight checklist

1. `h5py` present in the venv (the loader imports it lazily; nothing else needs
   it).
2. For each file: `load(...)` → `check_geometry(...)` → print
   `declared_noise(...)` into the results record BEFORE any fit.
3. Confirm the patch family is resolved on the actual grid: the factory's
   aliasing and ladder gates reject patches the grid cannot represent, and the
   REJECTION COUNT is the first number to read. PDEBench's 1-D grids are
   typically 1024 in x and 41–201 in t; the t direction is the tight one, and
   the C4 floor study measured what a coarse time grid does to the band.
4. Declare the term library per PDE family as a REGISTERED list before running
   (C3's discipline: |H| enters α directly).
5. Run the truth check where the generating law is known from the filename —
   `pdesystem.truth_check` — before reading any abstain as a finding.
6. Verify track: `verify.py` for scalars, `verify_system` for systems, from a
   held-out sample.

## After the run: an FNO as a PROPOSER (user suggestion, 2026-07-28)

Fourier Neural Operators do well on PDEBench's own forecasting leaderboard, and
a surrogate that forecasts well is a candidate PROPOSER — the same slot the
bounded LLM proposer occupies in `machine/`. It is admissible with no risk to
soundness for the structural reason the split already enforces: proposers see
`X_fit`/`X_sel`, never the certification split, so a proposer cannot manufacture
a certificate no matter how it was trained.

The concrete mechanism is better than "suggest a form": probe a trained operator
with small single-mode perturbations and read the LINEARIZED response. For a
linear PDE that is `exp(λ(k)Δt)`, and λ(k) is linear in the library's
coefficients (`−c_xx k² + i c_x k − i c_xxx k³ + …`), so a least-squares fit of
the probed symbol is an estimate of the operator rather than a guess at a form.

Where it would genuinely earn its keep is the case this instrument cannot reach:
the library is divergence-form-only by construction, so on files where an FNO
forecasts well and our vocabulary abstains, the probed symbol is evidence about
**what the vocabulary is missing** — a coverage instrument for negative results,
which is worth more than proposing laws we already recover.

Sequencing, and it is not optional: the PDEBench pass is REGISTERED, so a
proposer introduced mid-campaign would select the hypothesis using the data the
certificate then rests on. Run the registered pass first; then register the FNO
proposer as its own extension with its own predictions, scored against the
ground truth this pass establishes (e.g. "the proposed support contains the true
support at rate X on the families where we certify", plus the null: on
smooth-random fields it must still produce nothing certifiable). Two costs
stated plainly: a torch dependency, and an artifact whose error is not
declarable — fine for a proposer, and it must never touch the band.

## Known gaps, stated rather than discovered

- **Non-periodic boundaries.** PDEBench's 1-D CFD sets include non-periodic
  cases; the weak form itself is fine there (the test function vanishes inside
  the domain), but the spectral VERIFY track assumes periodicity, and the state
  certificate needs a basis diagonalizing the forward operator. Those runs get
  weak-form certification without the verify track, and must say so.
- **Non-divergence-form terms** (`u·u_xx` and friends) remain out of reach by
  construction — stated in the C0 registration, unchanged.
- **Per-field σ.** The band takes one declared σ; a system whose fields have
  genuinely different scales must declare the largest (conservative) until
  per-field Gram blocks exist.
- **Shocks.** PDEBench's Burgers at low viscosity and its CFD sets contain
  shocks. C4 measured what happens: the patch gate rejects the unresolved
  patches and the certificate refuses on the resolution bound. That is the
  correct behaviour and it means large parts of those files are outside the
  instrument's reach — expected, and to be reported as coverage, not hidden.
