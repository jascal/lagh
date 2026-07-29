# Direction: state certificates — inverting for initial conditions

**BUILT AND RUN 2026-07-28 — see `CASE_STUDY_PDE_C4.md` for the registration and
the measured results.** All four pieces below exist (one-sided-in-time test
functions with Romberg-in-time quadrature and analytic IC columns in
`lagh/weakform.py`; assembly, semantics and per-mode LP intervals in
`lagh/statecert.py`; the curriculum in `experiments/pde/run_c4_state.py`; the
forecast track in `verify.py::verify_state`). W1, W2, W4 and W5 met; W3 met in
direction. Two of the predictions needed their own geometry corrected before
they meant anything — the `exp(+νk²T)` curve lives in the BACK-PROPAGATION, not
in a window that starts at t = 0, and the Cole-Hopf family cannot form a shock —
and both corrections are in the case study.

**Scoped 2026-07-28 after a positive scoping probe.** Everything this program
has certified is a LAW over a stated domain. Recovering an initial condition is
a claim about one system's particular history, which is a different kind of
claim, and the decision to treat it as a separate category rather than sliding
into it is deliberate (user call, 2026-07-28).

## The mechanism: the boundary term the weak form has been throwing away

Every test function used so far vanishes at the patch boundary, which is what
kills the boundary terms in the by-parts identity. Drop that **in time only** —
a φ that is nonzero at t = 0 — and

    ∫∫ φ u_t dx dt  =  [∫ φ u dx]₀ᵀ − ∫∫ φ_t u dx dt

so the initial condition enters as `−∫ φ(x,0) u(x,0) dx`: a linear functional of
the unknown, with an ANALYTIC weight, inside the same linear system, under the
same declared band. No differentiated data, no new error model. Expand the
initial condition in a basis and each mode coefficient is a parameter — which is
exactly what `certify.parameter_interval` already measures.

The pleasant asymmetry: the initial-condition columns are `∫ φ(x,0)·b_j(x) dx`,
known analytic constants carrying **no data error at all**. The band is simpler
here than in the law case, not harder.

## Probe result (2026-07-28, `experiments/pde/probe_ic.py`)

The probe skips the factory: the law is known and the domain periodic, so the
forward operator is diagonal in Fourier space and mode k arrives multiplied by
d_k. It then asks the production interval search what the declared band permits
per mode. Prediction under test: half-width ≈ ε / d_k.

| system | d_k range | measured half-width / σ | half-width · d_k / σ |
|---|---|---|---|
| advection (unitary) | 1.0 | 1.12 – 1.62, **flat in k** | 1.12 – 1.62 |
| heat, T = 0.3, ν = 0.1 | 0.970 → 0.0133 (k = 1 → 12) | 1.62 → **148** | **1.33 – 1.97** |

So the resolution curve is `σ·exp(+νk²T)` up to an O(1) constant that stays
inside a factor of 1.5 across four decades of decay — the exact inverse of the
forward decay, measured, not asserted. Widths are linear in σ (identical
half-width/σ at σ = 1e-6 and 1e-4), consistent with C1b.

**A gotcha the probe found, worth not rediscovering:** `certify.epsilon` treats
`sigma` as RELATIVE to |y| by default (`prop=None` → `KAPPA·σ·|y|`). Absolute
field noise must be declared with `prop=ones`, or the band collapses wherever
the field crosses zero and nothing certifies — which is exactly what the first
probe run did, uniformly, and it looks identical to a real negative result.

## What a state certificate says

Draft semantics, to be fixed before anything runs:

- **Claim.** "Over the stated observation window and patch family, every initial
  condition whose mode coefficients lie in the reported intervals reproduces the
  observations within the declared band; modes outside the reported set are NOT
  determined." Undetermined modes are named, not silently dropped — an
  ill-posed inversion's honest output is a resolution statement.
- **Domain.** The observation window and the basis. A state certificate says
  nothing about other times, other solutions, or modes above the reported cut.
- **α.** With a fixed basis and a known law there is no search: |H| ≈ 1, so α
  reduces to a pure chance-agreement bound rather than a
  multiple-comparisons-corrected one. That is defensible but it is a DIFFERENT
  quantity from the α on a law certificate and must be labelled as such.
- **Resolution bound, stated up front.** dof = number of basis modes, so
  h = n − dof must stay positive with margin: at most as many modes can be
  certified as there are independent patch equations. This is the same mechanism
  as the 35-term-interpolation lesson, and it belongs in the registration, not
  in the results.

## Build plan

| piece | work | notes |
|---|---|---|
| one-sided-in-time test functions + boundary column in `weakform` | ~2–3 h | new ψ shape; IC columns are analytic constants |
| state-certificate assembly + semantics doc | ~1–2 h | reuses `check`, `parameter_interval`; no engine change |
| curriculum campaign (advection / heat / Burgers) | ~1–2 h | fields exist; runner mirrors `run_c1` |
| verify (recovered IC forward through `verify.py`) | ~0.5 h | nearly free |

**Total ≈ one session the size of the C0–C2 work.** Highest-risk piece: the
ill-conditioning. The per-mode interval search must happen in a basis where the
operator is diagonal, which is Fourier on a periodic domain — so anything
non-periodic is a later problem, not a first one.

## Registered predictions (for when it runs)

- **W1.** Advection: every mode in the reported set certifies, half-width flat
  in k within a factor of 2, and the true amplitude inside every interval.
- **W2.** Heat: half-width(k)/σ tracks `exp(+νk²T)` within a factor of 3 up to
  the cut; above the cut `parameter_interval` returns unbounded and the mode is
  reported UNDETERMINED rather than estimated.
- **W3.** The cut moves the way the physics says: later T or larger ν pushes it
  to lower k, and smaller σ pushes it up, roughly as `k_cut ≈ √(ln(1/σ)/(νT))`.
- **W4.** Burgers after shock formation: refusal, not a wide interval — the
  information is destroyed, not merely diluted.
- **W5.** Zero confident-wrong: no reported interval excludes the true
  amplitude.

## Explicitly out of scope for the first pass

Joint law-and-initial-condition inversion. Both are linear in the same system,
so it is tempting, but the identifiability is materially worse and it deserves
its own registration rather than being folded in.
