# PDE dev campaign C2 — the verify track — registration

**Registered 2026-07-28, before implementation.** C0 certified weak-form PDEs on
exact solutions; C1/C1b certified them under declared noise with interval
parameters. Both certify a RELATION AMONG PATCH INTEGRALS. This stage asks the
harder question a physicist would ask: **does the certified law, integrated
forward from an initial condition it never saw, reproduce the field?**

`DIRECTION_PDE.md` piece 4 called this the strongest falsification available,
and it is the same pattern as the RRab/Poisson verify-the-declared-form track:
the instrument states a form, then something outside the fitting machinery
checks it.

## What C2 claims, and what it must declare

For a law certified at noise σ, with parameter intervals from the certificate:

1. **Integrate forward** from a held-out initial condition — method of lines
   with spectral (FFT) spatial derivatives on the periodic domain, `solve_ivp`
   with tight tolerances. This is a NEW numerical error entering the claim, so
   it is declared like every other: a tolerance ladder (rtol, rtol/10) with the
   difference as the declared solver bound. Anything that does not converge on
   the ladder is refused, not banded.
2. **Propagate the parameter interval.** The certificate does not claim a value,
   it claims an interval, so the forecast is an INTERVAL too: integrate at the
   interval endpoints (and the centre) and take the envelope. A verify claim
   that ignored this would be claiming more than the certificate it tests.
3. **Compare against the true field** on the held-out solution, over the whole
   time window, with the total declared band = solver bound + interval envelope
   + the field's own declared noise σ.

Verdict per system: FORECAST-VERIFIED (every point inside the band), or a
refusal naming which of the three parts was exceeded.

## Registered predictions

- **V1 — the certified laws forecast.** For heat, advection and Burgers, at
  every σ where C1b certified, the forward integration from a held-out IC stays
  inside the declared band over the full window. A failure here would mean a
  law can certify on patch integrals while being wrong as dynamics, which would
  be the most important negative result of the whole arc.
- **V2 — the track discriminates.** A law perturbed by 10× its own certified
  interval half-width FAILS the forecast check for every system. If a 10×-wrong
  coefficient forecasts just as well, the verify track carries no information
  and should be reported as such rather than kept as decoration.
- **V3 — the forecast band grows with σ**, roughly linearly, because it is
  dominated by the parameter interval, which C1b measured to be linear in σ.
- **V4 — the solver error is not the binding term** at any σ > 0: the declared
  solver bound stays below the interval envelope. If it dominates, the ladder
  is not tight enough and the comparison is measuring the integrator, not the
  law.
- **V5 — zero confident-wrong**: no system reported FORECAST-VERIFIED whose
  forecast in fact leaves the band, and no refusal for a law that is right.

## What would end the arc

V1 failing for a law that certified: it would mean weak-form certification is
not evidence about dynamics, and the honest response would be to demote every
PDE certificate in C0/C1 to a statement about patch integrals only.

## Results

Run 2026-07-28 (`experiments/pde/verify.py`, `experiments/pde/run_c2.py`,
results `experiments/results/pde_c2.json`). Every certificate C1b produced — 18
rungs across heat, advection and Burgers — integrated forward from a **fifth
initial condition** that no stage saw (not fitted, not certified against, not
the held-out solution), on a periodic grid the patch family never used.

| system | law-verify | data-verify | law error at σ=1e-4 | 10×-wrong control |
|---|---|---|---|---|
| heat | 6/6 rungs | 6/6 | 4.9e-6 | fails, 5120/5376 points |
| advection | 6/6 | 6/6 | 5.7e-6 | fails, 5120/5376 |
| Burgers | 6/6 | 6/6 | 1.2e-5 | fails, ~5120/5376 |

- **V1 — met.** Every certified law forecasts inside its declared band, at every
  σ, on a solution it never saw. The weak-form certificates are evidence about
  dynamics, not only about patch integrals.
- **V2 — met, decisively.** Pushing every coefficient 10 interval half-widths
  off centre fails at ~95% of points on **every** rung. The track carries
  information; it is not decoration.
- **V3 — met.** The forecast band is linear in σ: the interval envelope runs
  4.5e-9 → 4.5e-5 and the initial-condition term 3.9e-8 → 3.9e-4 as σ runs
  1e-8 → 1e-4, matching C1b's linear interval widths.
- **V4 — met.** The declared solver bound (~1e-10, uniform) is the smallest term
  at every σ > 0, three to six orders below the interval envelope. At σ = 0 it
  is the only term, which is the correct degenerate case.
- **V5 — met.** No rung reported verified whose forecast leaves its band, and
  no correct law was refused.

### The declaration this stage forced

The first C2 run FAILED advection at 1–5 points of 5376 on every noisy rung
while heat and Burgers passed. Diagnosis, and it was a band error rather than a
law error: the forecast from the clean initial condition matches the clean field
to **3e-11**, so the law was exact; what the band omitted is that **a forecast
started from measured data carries that data's noise for the whole window**.
Advection is non-dissipative, so the initial condition's noise never damps —
measured at 2.8σ, against a band carrying only the target's 4σ. Heat passing was
luck, not correctness: its own term measures 3.9σ over the same window.

Two consequences, both now in the code:

- The band carries an `ic_noise_bound`, MEASURED by re-integrating from
  perturbed initial conditions. It is computable from measured data alone: σ is
  declared, so the perturbation can be simulated without knowing the clean
  field.
- **Law-verify and data-verify are kept apart**, because they are different
  claims: "the law reproduces the field" (clean IC, clean field, band = solver +
  interval) and "the law predicts measurements" (measured IC, measured field,
  band additionally carrying the target's noise and the propagated IC noise).
  The second is strictly harder, and a dev campaign that reports only the first
  is quietly claiming the easier one.

Also measured, and fixed the same way as the patch ladder: the solver's
tolerance-ladder difference **under-covers pointwise** (~6% of points on exact
heat, by up to 1e-12) because adaptive stepping makes the per-point error
fluctuate while the ladder difference does not, and part of the residual is
float roundoff the ladder cannot see. The bound is declared uniformly over the
domain, floored at the program's machine term.

### Amended 2026-07-28 (C3/PDEBench readiness): the integrator

The explicit RK45 method-of-lines used above is DIFFUSION-LIMITED — its stable
step goes like `1/(ν k_max²)`, so on a 512-point grid at ν = 0.2 it needs ~250k
steps per trajectory and a verify run does not finish. Measured while preparing
the PDEBench pass: the same forecast that hangs there completes in 0.09 s under
an exponential integrator (ETD-RK2), which solves every LINEAR term exactly per
substep and steps only the nonlinearity.

`verify.integrate(scheme=...)` now carries both. `"direct"` is the default, so
the C2 results above stand exactly as run; `"etd"` is what PDEBench-scale grids
use, and it declares its error the same way everything else here does — by a
ladder, in substeps rather than in tolerance. On the C2 fields the two agree to
1e-9 or better, and on the linear parts the exponential scheme is the more
accurate of the two (heat: 3.8e-15 against 1.3e-11).

An integrating factor unwrapped from t₀ rather than per substep was tried first
and is simply wrong for diffusion: it carries `exp(+ν k² t)`, which overflows
within a fraction of a time unit at any real resolution.

### Open

Longer windows (this ran to t = 0.4, where Burgers has not yet steepened);
shock formation, where the spectral integrator itself would need a declared
treatment; and the verify track against a law certified from REAL data, which
is the point at which it stops being a self-check.
