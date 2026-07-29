# Direction: systems of PDEs

**BUILT AND RUN 2026-07-28 — see `CASE_STUDY_PDE_C3.md` for the registration and
the measured results.** All five pieces below exist (`lagh/weakform.py`
multi-field + n-D geometry, `lagh/pdesystem.py` driver + conjoined certificate,
`experiments/pde/verify.py::verify_system`, `experiments/pde/run_c3_systems.py`
curriculum). The whole curriculum certifies: linear coupled, FitzHugh-Nagumo,
Brusselator, shallow water, and 2-D Navier-Stokes vorticity. Y1-Y5 all met, with
one correction the run forced (Y4: spectrally poor data creates a degeneracy the
solution holdout cannot break, which the constrained-input path correctly
reports as a domain-restricted certificate rather than a wrong one) and one
engine change (`discover(linear_basis=True)`, without which the greedy proposal
channels found nothing for a 5-term support over 8 declared columns).

The design below stands as written.

**Scoped 2026-07-28 after a positive scoping probe.** Coupled fields are where
real PDE discovery lives — reaction–diffusion, shallow water, Navier–Stokes —
and traffic, the arc's intended real-data case study, is already a system in
disguise (`ρ_t + q_x = 0` plus a closure).

## What transfers unchanged

Nothing in the weak-form machinery is single-field by nature:

- **The by-parts identity.** `∫φ ∂^α(g) = (−1)^|α| ∫(∂^αφ) g` does not care
  whether g is `u²/2` or `u·v`.
- **The band.** With independent per-field noise the sensitivity vectors simply
  concatenate, so the Gram becomes block-diagonal over fields and the band is
  still the one quadratic form `a'Ga`. Correlated field noise would populate the
  off-diagonal blocks; the formula is unchanged.
- **Solution holdout, multi-scale pooling, row normalization, interval
  parameters, the verify track** — all field-agnostic.
- **`lagh/systems.py`** already provides the system-level parts a per-equation
  loop cannot: cross-equation constant agreement, invariant discovery, and a
  union-bound α over the conjoined certificate.

## Probe result (2026-07-28, `experiments/pde/probe_system.py`)

The probe needs NO factory extension, because for a LINEAR system every feature
is a single-field term of one field or the other: run the existing `build()` on
each field over the SAME patches and concatenate. It tests the genuinely new
part — two targets sharing one row set, each equation's features spanning both
fields.

Exactly-solvable system (2×2 propagator per Fourier mode, so no integrator error
enters):

    u_t = 0.1 u_xx + 0.5 v
    v_t = 0.05 v_xx − 0.3 u

| σ | u-equation | v-equation |
|---|---|---|
| 0 | `[u:u_xx]/10 + [v:u]/2` — exact rationals | `−3[u:u]/10 + [v:u_xx]/20` — exact |
| 1e-6 | correct support, coefficients to 7 digits | correct support, 7 digits |
| 1e-4 | correct support, 4–5 digits | correct support, 4–5 digits |

144 rows over 4 initial conditions, 8 columns, certifying on a held-out
solution. **Both equations of a coupled system certify their true support, with
zero changes to the factory or the engine.**

**A probe bug worth recording, and the discipline that caught it.** The first
run abstained structurally at every σ including 0. Before reading that as an
identifiability finding, the check was: does the TRUTH certify? It did not —
residual 0.85 against a band of 4e-10 — so the fields were wrong, not the
instrument. The cause: the 2×2 propagator acts on the amplitude pair of the
SAME basis function, and independent random phases for u and v make the coupling
term `B·v` spatially orthogonal to u, so the system does not hold at all. One
shared phase per mode fixed it. A null result that is really a construction bug
looks exactly like a finding; checking the truth against its own band first is
what separates them.

## What is genuinely new, and what it costs

| piece | work | notes |
|---|---|---|
| multi-field library (`u·v`, `u²v`, `∂_x(uv)`) | ~2–3 h | `Term.gexpr` over several symbols; `build()` takes a dict of fields; the Gram spans concatenated sensitivities |
| system driver: per-equation targets over shared rows | ~1–2 h | one target per equation, features spanning all fields, solution holdout unchanged |
| conjoined certificate + union-bound α | ~1 h | `lagh/systems.py` has the mechanism; wire the PDE case to it |
| verify with a vector RHS | ~1 h | `verify.py` generalization; the IC-noise term applies per field |
| curriculum + registration | ~1–2 h | below |

**Total ≈ one session**, and the linear case already works today.

The library growth is the thing to watch: cross terms multiply fast and |H|
enters α directly, so the vocabulary must be a REGISTERED list per curriculum
stage, not a generated cross-product. That is a registration discipline, not a
code problem.

## Curriculum

1. **Linear coupled** (done in the probe): certifies exactly, no new code.
2. **Reaction–diffusion** — FitzHugh–Nagumo or Brusselator. First nonlinear
   coupling, and the first stage needing a reference solver, which means a
   DECLARED solver error (the C2 tolerance-ladder treatment) rather than an
   assumed-exact field.
3. **Shallow water** — a conservation-law system with a real flux; the first
   place conservation-form detection means something.
4. **Navier–Stokes vorticity–streamfunction** — where the closed
   constrained-input machinery should earn its keep: `∇²ψ = −ω` is a
   machine-exact input constraint of exactly the kind that closure handles.

## Registered predictions (for when it runs)

- **Y1.** Every equation of a certified system has the true support, and every
  reported parameter interval contains the truth — the C1b standard, per
  equation.
- **Y2.** The conjoined certificate's union-bound α is dominated by the weakest
  equation, and is reported as such rather than as a product of the strongest.
- **Y3.** Single-solution data refuses for systems too, and MORE often than for
  scalars: a coupled system on one trajectory can satisfy relations that hold in
  both equations at once.
- **Y4.** Nonlinear coupling (stage 2) needs spectrally richer initial data than
  the linear case for the same identifiability — measurable as the number of
  modes at which the abstain turns into a certificate.
- **Y5.** Zero confident-wrong across the curriculum.
