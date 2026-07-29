# PDE dev campaign C3 — systems of PDEs — registration

**Registered 2026-07-28, before implementation**, from the scoping probe and the
direction doc `DIRECTION_PDE_SYSTEMS.md` (predictions Y1–Y5 were registered
there and are restated below unchanged). C0–C2 certified SCALAR weak-form laws:
exact coefficients on clean fields, interval parameters under declared noise,
and forward-integration verification from unseen initial conditions. Coupled
fields are where real PDE discovery lives — reaction–diffusion, shallow water,
Navier–Stokes — and the arc's intended real-data case study (traffic) is a
system in disguise.

## What C3 claims

One row set, one target per equation, features spanning every field, and ONE
conjoined certificate whose α is a union bound over the equations.

Nothing in the weak-form machinery was single-field by nature: the by-parts
identity does not care whether g is `u²/2` or `u·v`; with independent per-field
noise the sensitivity vectors concatenate, so the Gram is block-diagonal over
fields and the band is still the single quadratic form `a'Ga`; the solution
holdout, multi-scale pooling, row normalization, interval parameters and the
verify track are all field-agnostic.

Four things are genuinely new and all four are declared:

1. **Multi-field terms** — `g` reads several fields (`u*v`, `u²v`, `h*u`), so
   the Gram sums per-field sensitivity blocks (`lagh/weakform.py`).
2. **n-dimensional geometry** — patches over (x, …, t), which stage 4 (2-D
   Navier–Stokes) needs and which the 2-D half of PDEBench will need.
3. **A declared FIELD error** — stages 2–4 need a reference solver, and a
   solver's error is a new undeclared quantity entering the claim. Each solved
   field arrives with a MEASURED bound (tolerance ladder + spatial-resolution
   ladder), propagated into the band through the L1 sensitivity, because a
   solver error is one fixed function and nothing cancels across it.
4. **The conjoined certificate** — cross-equation constant agreement plus a
   union-bound α (`lagh/systems.py`, `lagh/pdesystem.py`).

**The vocabulary of each stage is a REGISTERED list, never a generated
cross-product.** Cross terms multiply fast and |H| enters α directly; that is a
registration discipline, not a code problem.

## Curriculum

1. **Linear coupled** — `u_t = a u_xx + b v`, `v_t = c v_xx + d u`, exactly
   solvable per Fourier mode (2×2 propagator), so a stage-1 failure is the
   instrument's and never an integrator's.
2. **Reaction–diffusion** — FitzHugh–Nagumo (nonlinear in ONE field: the cubic)
   and Brusselator (the first genuine CROSS term, `u²v`). Split deliberately, so
   a failure can be attributed to "nonlinear" or to "cross-coupled".
3. **Shallow water** — a conservation-law system with a real flux; the first
   place conservation-form detection means something.
4. **Navier–Stokes vorticity–streamfunction (2-D space)** — where the closed
   constrained-input machinery should earn its keep: `∇²ψ = −ω` is a
   machine-exact input constraint of exactly the kind that closure handles.

## Registered predictions

- **Y1.** Every equation of a certified system has the true support, and every
  reported parameter interval contains the truth — the C1b standard, per
  equation.
- **Y2.** The conjoined certificate's union-bound α is dominated by the weakest
  equation, and is reported as such rather than as a product of the strongest.
- **Y3.** Single-solution data refuses for systems too, and MORE often than for
  scalars.
- **Y4.** Nonlinear coupling needs spectrally richer initial data than the
  linear case for the same identifiability — measurable as the number of modes
  at which the abstain turns into a certificate.
- **Y5.** Zero confident-wrong across the curriculum.

## Results

Run 2026-07-28 (`experiments/pde/run_c3_systems.py`, fields
`experiments/pde/systems_fields.py`, results
`experiments/results/pde_c3.json`). Four solutions per stage, three pooled patch
scales, certification on patches from a HELD-OUT solution, σ ∈ {0, 1e-6, 1e-4}.

| stage | equations | law recovered | σ range |
|---|---|---|---|
| 1 linear coupled | 2/2 | `u_t = [u:u_xx]/10 + [v:u]/2`, `v_t = −3[u:u]/10 + [v:u_xx]/20` | 0 – 1e-4 |
| 2a FitzHugh–Nagumo | 2/2 | `u_t = [u:u_xx]/20 + [u:u] − [u:u^3]/3 − [v:u] + 1/2` | 0 – 1e-4 |
| 2b Brusselator | 2/2 | `u_t = [u:u_xx]/50 − 27[u:u]/10 + [u^2*v] + 1` | 0 – 1e-4 |
| 3 shallow water | 2/2 | `h_t = −[(hu)_x]`, `(hu)_t = −[(hu^2)_x] − 981[(h^2)_x]/200` | 0 – 1e-4 |
| 4 Navier–Stokes (2-D) | 1/1 | `w_t = [w_xx]/50 + [w_yy]/50 − [(uw)_x] − [(vw)_y]` | 0 – 1e-6 |

- **Y1 — met.** Every certified equation carries the true support and every
  reported interval contains the truth, at every σ, across all four stages.
  Shallow water recovers the conservation form with coefficient exactly −1 and
  `g/2` as the exact rational 981/200; the 2-D vorticity equation recovers
  ν = 1/50 in both Laplacian components independently.
- **Y2 — met.** The conjoined α is the union bound (log-space `max + log₁₀ n`),
  and the runner reports WHICH equation dominates it rather than the best one.
- **Y3 — met, and it is the holdout that does the work.** One solution abstains
  `single-solution` for both equations; the same data under a ROW split (what a
  naive pipeline does) certifies both equations happily. The gap between those
  two lines is the whole argument for the solution holdout.
- **Y4 — met, and it produced the most interesting result of the stage.**
  See below.
- **Y5 — met.** 40 certified equations across the curriculum, zero
  confident-wrong (no certified law disagrees with the truth by more than its own
  declared band anywhere on its certified domain) and zero SILENT degeneracies
  (every support that differs from the truth carries the certificate's own
  domain-restriction note).

### Y4, and what "spectrally poor data" actually does

With single-mode initial data every solution in the family satisfies
`u_xx ≡ −u` **exactly**, so the degeneracy holds across the whole family and the
solution holdout cannot break it — holding out a solution only helps when the
held-out solution breaks the relation. The linear system then certifies

    u_t = −[u:u]/10 + [v:u]/2

which is not the general PDE. It is also not wrong: the engine's
constrained-input detector fired, the certificate is issued **domain-restricted**
("inputs satisfy … machine-exact"), and the certified law agrees with the truth
to 1.3e-3 of its own band everywhere the domain claim applies. At two or more
modes the degeneracy disappears and the true support comes back.

This forced a scoring change that is worth keeping: **support equality is the
wrong test on its own.** The campaign now scores three separate things —
support equality, agreement with the truth over the certified domain, and
whether the certificate flagged its own domain restriction — and reserves
"confident-wrong" for a law that DISAGREES beyond its band. A support that
differs while agreeing everywhere the claim applies is an under-determination,
and a *silent* one only if the certificate failed to say so.

Brusselator, by contrast, certifies its true support even at one mode: the
nonlinear cross term breaks the degeneracy that sinks the linear system.

### The verify track, at system scale

`experiments/pde/run_c3_verify.py` (results
`experiments/results/pde_c3_verify.json`) asks C2's question of systems:
integrate the certified system forward, in the weak form's own vocabulary, from
an initial condition no stage of the pipeline has seen, and ask whether every
evolved quantity stays inside the declared band. **All 12 rungs (4 stages × 3 σ)
verify, law and data, with zero points outside**, and the control — every
coefficient pushed 10 interval half-widths off centre — fails at ~98% of points
on every rung. The track carries information.

Two things it forced, and both are the C2 lessons re-earned at system scale:

- **Law-verify must start from the CLEAN initial state.** The first version
  forecast both claims from the measured start, so the law claim inherited the
  initial condition's noise without a band for it and failed at every σ > 0
  while the law was exact. Law-verify tests the law; data-verify tests the law
  as a predictor of measurements, and carries 4σ plus a MEASURED
  initial-condition-noise term (re-integration from perturbed starts).
- **The forecast grid must be genuinely periodic.** The stage-1 fields are
  generated on an endpoint-INCLUSIVE grid, which the weak-form patches do not
  care about and the spectral integrator very much does: 18244 of 41634 points
  fell outside the band at σ = 0 with the law exact. This is the same trap
  `pdebench.check_geometry` exists to catch on external data, met here on our
  own.

The envelope is a per-parameter triangle bound rather than the scalar track's
full corner product (a 2-equation system with 3 parameters each has 3⁶ = 729
corners); it upper-bounds the corner envelope exactly for a linear-in-parameters
response and conservatively otherwise, and it is reported as `envelope_method`
rather than silently substituted.

### Nulls

Two smooth fields solving no system, and an i.i.d. pair, both at σ = 0 and
σ = 1e-4: structural abstain everywhere, and the i.i.d. pair loses every patch
to the resolution gate ("the grid does not represent this field"). The truth
check reports the smooth null's truth/band at 4e12 — the registered system is
nowhere near holding, which is the point of running it.

### The engine change this stage forced

`discover(..., linear_basis=True)`. A weak-form PDE law IS a linear combination
of the declared library's columns, and with 8 declared columns and a 5-term true
support the greedy proposal channels (STLSQ, OMP, size-5-from-top-6-singles)
proposed **nothing at all** — the engine returned zero candidates while the
truth sat four orders inside its own band. Under `linear_basis` the supports are
enumerated exhaustively over the library (budgeted, and the budget is reported),
and the engine builds no products or powers of patch integrals — which is also
the direct fix for the C1b `u_xx*[1]^(3/2)` failure mode. It is a declaration
about the claim, not a search-budget knob, and it is off by default.

### The discipline that kept a null honest

`pdesystem.truth_check` runs BEFORE any abstain is read as a finding: it asks
whether the TRUE law sits inside its own declared band on these rows. The system
scoping probe bought this the hard way — its first run abstained structurally at
every σ including 0, which looks exactly like an identifiability finding and was
a construction bug (independent random phases made the coupling term spatially
orthogonal, so the system did not hold at all; the truth missed its own band by
0.85 against a band of 4e-10). Every stage now reports truth/band alongside its
verdict.

### Open

Systems under a MEASURED rather than declared field error (a PDEBench field
cannot have its solver ladder re-run after the fact); per-field σ that differ
(the band takes one declared σ and the conservative reading is the largest);
shocks in the shallow-water stage; and 3-D geometry, which the factory supports
but nothing has exercised.
