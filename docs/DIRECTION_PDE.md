# Direction: PDE support (the weak-form arc)

**Scoped 2026-07-28; C0 BUILT AND RUN the same day.** The design below stands
as written — see `CASE_STUDY_PDE_DEV.md` for the registration and the measured
results, and `lagh/weakform.py` for the factory. Status: pieces 1–3 built
(weak-form factory with a declared patch band, patch-ε via the new
deterministic `hard` channel, multi-solution holdout), heat/advection/Burgers
certified at exact rational coefficients from held-out solutions, single-
solution data refused structurally. Pieces 4–5 (forward-integration verify
track, the rest of the curriculum) and everything under declared noise remain
`open`.

## What a PDE claim is, in this program's terms

A PDE law is a relation among derivative fields — `u_t = ν·u_xx − u·u_x`.
The engine's existing machinery transfers almost unchanged for the
CANDIDATE side: a term library over derivative combinations (SINDy's
candidate space) driven by the linear channel, OMP/CAP-Q/CAP-T support
proposal, `float_pinned`, `refit_minimal`, the significance gate, and
`systems.py` for coupled fields. What does NOT transfer is the data model.
Two hard problems, and both already have program precedents:

### Problem 1 — derivatives are not data

lagh certifies |pred − y| ≤ ε per point under a DECLARED error model; u_t
and u_xx estimated by finite differences amplify noise like σ/hᵏ plus
truncation — neither measured nor declared. **The honest route is the weak
form**: integrate the candidate law against smooth compactly-supported test
functions, move all derivatives onto the test functions by parts. Each
"data point" becomes a patch integral of the RAW field against known
analytic test-function derivatives, and its error bound is computable:
quadrature error + measurement noise through a known linear functional.
Certification becomes **per-patch**: "every test-function residual within
ε_patch over the stated patch family." The four-term ε model, vacuity, and
α accounting (patches as n, dof from the form) carry over directly. This is
the σ_rep/float32 lesson generalized: never certify against numbers whose
error was not declared. (Precedent: WSINDy uses weak forms for robustness;
we use them for *declarability*.)

### Problem 2 — single-solution under-determination

On one solution, the field satisfies on-shell constraints (a traveling wave
has u_t ≡ −c·u_x identically), so rival PDEs differing by multiples of an
on-shell-zero quantity are indistinguishable — this is EXACTLY the
constrained-input coherence problem (closed 2026-07-28 for algebraic input
constraints) fused with the LawSystemBench single-trajectory
identifiability finding (H3). The fix is the same one H3 found: **multiple
solutions** (diverse initial conditions and forcings), with coherence
probing across solutions rather than on an off-manifold box. A structural
abstain on single-IC data is the CORRECT verdict and should be a registered
prediction of the dev campaign, not a surprise. The closed
constrained-input machinery (`input_constraints` / manifold-probe
coherence / canonicalization modulo the constraint ideal) is the template —
here the "constraint" is the on-shell relation itself, detectable the same
way over the derivative-feature matrix.

## The new pieces to build

1. **Weak-form feature factory** (the big piece): given field snapshots
   u(x, t) on a grid, a family of test functions φ (tensored bump/Legendre
   windows at multiple scales/locations), produce per-patch integrals
   ⟨∂ᵅφ, u⟩ and products (⟨φ, u·u_x⟩ via integration-by-parts-compatible
   forms), each with a DECLARED bound: quadrature truncation (h-order,
   computable from φ smoothness) + noise propagation (‖φ‖·σ). Output shape:
   a design matrix over patches — then the existing linear channel does the
   discovery.
2. **Patch-ε model**: `epsilon()` already supports per-point `se` (the
   anisotropy plumbing) — patch bounds enter there; nothing new in the
   checker.
3. **Multi-solution coherence**: probe = patches from held-out solutions
   (other ICs), the manifold-probe pattern.
4. **PDE verify track** (strongest falsification): integrate the candidate
   law forward with a high-order solver from held-out ICs; certify forecast
   agreement within a declared tolerance — the verify-the-declared-form
   pattern (RRab/Poisson precedent) for dynamics.
5. **C-PDE curriculum**: constant-coefficient linear (heat, advection,
   wave) → quadratic nonlinearity (Burgers, KdV) → variable coefficients →
   conservation-law form detection.

## The campaign plan (established discipline throughout)

- **Dev campaign: PDEBench / the classic PDE-FIND suite** (Burgers, KdV,
  heat, NS vorticity — known ground truth). This is the NewtonBench role:
  registration, null runs (weak-form nulls: random fields must certify
  nothing), escalating noise (the R-noise treatment for the patch-ε model),
  single-IC identifiability predictions (must abstain), multi-IC recovery
  predictions.
- **Real-data case study: highway traffic (NGSIM, or the cleaner German
  HighD drone data)** — chosen over ERA5 (reanalysis is model output — the
  FLAME circularity at continental scale) and lab data (none public at
  quality). Traffic has the program's signature two-strata structure:
  - *Certifiable stratum*: car conservation ρ_t + q_x = 0 is EXACT —
    vehicles are countable; in weak form (cars-in − cars-out = accumulation
    over a road-time patch) it certifies against trajectory data with a
    declared counting/binning error, no closure assumptions. A certified
    conservation law from real traffic, with α, would be a first anywhere.
  - *Conjecture stratum*: the flux closure q(ρ) — the "fundamental
    diagram" — is empirical, scattered, hysteretic, multivalued near
    capacity. It must NOT certify; banded conjectures (free-flow speed,
    capacity drop, backward wave speed −20…−15 km/h) mirror the
    exoplanet/catalysis treatment.
- ERA5 becomes interesting only later, as a pipeline-decode exercise
  (certifying which balances the assimilation system enforces).

## Effort estimate

Weak-form factory + patch-ε declaration: the dominant piece (roughly the
size of the catalysis rebuild). Curriculum + dev campaign: the standard
registration mechanics. Engine core: nearly untouched. Multi-solution
coherence: small, template exists. Verify track: moderate (needs a
reference solver; `scipy.solve_ivp` + method-of-lines is enough for the
curriculum PDEs).

## Pickup checklist

1. ~~Read this doc, `REACH_ENVELOPE.md`, `MUNTZ_ARBITRATION.md`, the memory
   index.~~ Done 2026-07-28.
2. ~~Weak-form factory `lagh/weakform.py` + registration
   `docs/CASE_STUDY_PDE_DEV.md`.~~ Done: factory, tests, C0 registered and
   run. One registered prediction was falsified (single-solution KdV
   certified the on-shell relation under a row split) and closed by the
   solution-holdout gate — read that section before extending anything.
3. ~~Locally-generated solutions before any download.~~ Done: C0 is entirely
   analytic (`experiments/pde/fields.py` — Fourier heat, translation
   advection, Cole–Hopf Burgers, tanh wave, KdV soliton).

### (a) the errors-in-variables ε model — **BUILT 2026-07-28**, see
`CASE_STUDY_PDE_C1.md` for the registration, the three engine changes it
required (per-candidate ε, noise-corrected resolution gates,
interval-parameter certificates) and the measured noise ladder. The design as
sketched below is what was built; the one thing it did not anticipate is that
the *claim shape* had to change too — under noise a physical coefficient is
determined to an INTERVAL, and demanding an exact rational there was the wrong
question, not a stricter one.

### Original design sketch (kept for provenance)

Under declared σ the weak-form design matrix is noisy too, and C0's band
handles that by bounding the coefficients (`coeff_max`) and summing per-column
errors. That is honest but loose, and for the stochastic part it is also
*wrong in shape*: every column is a linear functional of the SAME noise
realization, so the errors do not add in absolute value, they cancel. For a
candidate law y = Σ c_k X_k the residual carries

    δ_y − Σ c_k δ_k  =  Σ_i (w_y,i − Σ_k c_k w_k,i g'(u_i)) e_i

whose band is `4σ · ‖w_y − Σ_k c_k w_k g'(u)‖₂` — one norm over the combined
weight vector, not a sum of norms. This is both tighter and coefficient-aware,
and it needs one engine extension: ε assembled **per candidate** rather than
once, i.e. `check()` taking a band *function* of the expression under test
alongside today's fixed vector. Deterministic parts (quadrature, roundoff) do
not share a realization and keep the Σ|c_k|·q_k combination they have now.

Registering the shape here because it is the piece that decides whether noisy
PDE certification is possible at all: with the loose sum-of-columns band, a
realistic σ swamps the signal well before the noise ladder gets interesting.

**Next (in order):** ~~(a) the σ > 0 ladder~~ — done, C1: heat/advection/Burgers
certify with interval parameters through σ = 1e-4 (advection to 1e-3), every
interval containing the truth, zero confident-wrong. ~~(a′) multi-scale patch families~~ — done, C1b: pooled scales plus ROW
NORMALIZATION by each patch's own ∫φ (a pooled family is only as sound as its
worst-conditioned member, and an un-normalized varying `1` column turns the
patch scale into a state variable). Parameters 1.4–3.4× more tightly
determined. ~~(b) the verify track~~ — done, C2: every one of 18 certificates forecasts
inside its declared band from a fifth, unseen initial condition, and a
10x-interval-wrong law fails at ~95% of points. It forced the
initial-condition-noise term into the band and the law-verify / data-verify
split. ~~systems of PDEs~~ — done, C3 (`CASE_STUDY_PDE_C3.md`): the whole
curriculum certifies, linear coupled through **2-D Navier–Stokes vorticity**,
with a conjoined union-bound α; it added multi-field terms, n-dimensional patch
geometry, a MEASURED field-error channel for solver-produced fields, and the
`linear_basis` engine declaration. ~~state certificates~~ — done, C4
(`CASE_STUDY_PDE_C4.md`): initial conditions inverted mode by mode with
UNDETERMINED as a first-class outcome, per-mode intervals as exact LP
projections of the feasible set, and the exponential ill-posedness located where
it actually lives (back-propagation, not the reading). Conservation-law form
came free with C3 stage 3 (shallow water certifies `h_t = −(hu)_x` with
coefficient exactly −1).

~~(d) PDEBench~~ — run 2026-07-29 as a RECONNAISSANCE/dev pass, not a scored
read: `CASE_STUDY_PDEBENCH.md`. Six families, six verdicts, zero confident-wrong,
nine defects found in this program's own code. Its durable output is the
`error provenance` direction and two process rules now in `STRATEGY.md`.

**Next (in order):**

**(c) VARIABLE COEFFICIENTS — and Darcy priced it.** The registered slot said
"variable coefficients" as though it were a library extension. It is not. For a
coefficient field that is ITSELF DATA the weak form's central guarantee breaks:
the library is `∂^α(g(fields))` with g POINTWISE, and `∇·(a∇u)` has an integrand
`a∇u` pairing a field with a DERIVATIVE of another field. By-parts once leaves
`∫∇φ·a∇u`, still a data derivative; moving it again gives `∫∇·(a∇φ)u`, which
needs `∇a`. Rearranging does not escape it — `(a u_x)_x = (au)_xx − a_xx u −
a_x u_x`. Measured on PDEBench 2-D Darcy, where the piecewise-constant case IS
reachable (`∇²u = −β/a` on patches interior to one phase, β = 0.1000 recovered)
and the general case is not. So the work is one of:
  * a DECLARED error model for `∇a` — spectral differentiation of a band-limited
    coefficient with a stated aliasing bound, which the factory's own resolution
    gates could check; or
  * a MIXED formulation carrying the flux `q = a∇u` as a field, which certifies
    `∇·q + β = 0` exactly but needs q measured, not derived; or
  * the honest restriction: certify only where `a` is locally constant, and
    report the domain — which is what the Darcy run did.
The choice is a registration decision, not an implementation detail.

**DECIDED 2026-07-29 (user): the honest restriction.** Certify where `a` is
locally constant and REPORT THE DOMAIN. The two rejected routes each buy reach
with something this program does not want to spend: a declared `∇a` error model
adds another hand-set number in a slot where a wrong one is invisible (the
session that made this decision had just spent a correction on exactly that), and
the mixed formulation needs `q = a∇u` MEASURED, which the only live dataset does
not ship. What the chosen route costs is a vocabulary, not an approximation:
`certify.domain_qualifier` + `conjoin_determination`, now built — a restriction
that scopes an entire determination record, with the rule that records made on
different domains REFUSE to conjoin (their conjunction is defined only where both
were established, and nothing here can intersect two predicates). That is the
same qualifier the DOMAIN dimension of partial determination needed, so one
design closes both open items; see `DIRECTION_STOCHASTIC.md` step 0. Reach stays
piecewise-constant `a`, stated rather than discovered.

**Applied the same day.** `discover_equation(..., qualifier=)` carries the domain
on BOTH paths — a domain-restricted abstain is still a domain-restricted
statement — and where there is no partial content at all it emits a record with
the domain and no components, because "nothing determined, over this region" is
the honest reading and `None` was being read as a claim about the whole field.
Darcy re-run in both phases (`darcy_beta0.1_qualified`,
`darcy_beta0.1_lowphase_qualified`): verdicts unchanged (both
ABSTAIN[structural]), and `run_darcy_domains.py` confirms on that real output
that the two phases refuse to conjoin. Every campaign predating the kwarg passes
`None` and is unaffected by construction, which is what the gating rule asks for.

**(e) THE TRAFFIC CASE STUDY** — unchanged, and now better prepared: its
two-strata structure (exact conservation, conjectured closure) is the same shape
as the CFD result, where continuity needed a 1e-4 declaration and momentum
1e-2.

One general question the ladder raised and did not settle: **should
interval-parameter certificates be the default under any declared noise**, not
just for a declared basis? An exact rational is the right claim for a
definitional identity and the wrong one for a measured coefficient, and that
distinction is not specific to PDEs.
