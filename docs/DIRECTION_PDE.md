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
split. Then, both scoped 2026-07-28 after POSITIVE scoping probes and each with its
own direction doc: **systems of PDEs** (`DIRECTION_PDE_SYSTEMS.md` — the linear
coupled case already certifies with no factory change) and **state
certificates**, i.e. inverting for initial conditions (`DIRECTION_PDE_STATE.md`
— a different KIND of claim from a law, deliberately kept separate; the
resolution curve measures as sigma*exp(+nu k^2 T), the exact inverse of the
forward decay). Then (c) variable coefficients and conservation-law form; then
(d) PDEBench, and only then the traffic case study.

One general question the ladder raised and did not settle: **should
interval-parameter certificates be the default under any declared noise**, not
just for a declared basis? An exact rational is the right claim for a
definitional identity and the wrong one for a measured coefficient, and that
distinction is not specific to PDEs.
