# Direction: PDE support (the weak-form arc)

**Scoped 2026-07-28, not yet started.** This document is the pickup point:
it records the full design so the arc can begin in a fresh context with no
re-derivation. Status of everything here: `open` (design), except where
marked as established program results it builds on.

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

## Pickup checklist (fresh context)

1. Read this doc, `REACH_ENVELOPE.md`, `MUNTZ_ARBITRATION.md` (results
   section — confirm the arbitration arc landed), and the memory index.
2. Start with the weak-form factory as `lagh/weakform.py` + a registration
   doc `docs/CASE_STUDY_PDE_DEV.md` (predictions before any PDEBench run:
   null-certification zero; heat/advection certify multi-IC; single-IC
   Burgers abstains structurally; noise ladder bands).
3. PDEBench data: start with locally-generated solutions (method of lines,
   exact ICs) before downloading anything — the dev campaign's C0 needs no
   external data at all.
