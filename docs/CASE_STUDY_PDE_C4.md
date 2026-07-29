# PDE dev campaign C4 — state certificates — registration

**Registered 2026-07-28, before implementation**, from `DIRECTION_PDE_STATE.md`
(predictions W1–W5 were registered there and are restated below unchanged).

Everything this program has certified so far is a LAW over a stated domain.
Recovering an initial condition is a claim about one system's particular
history — a different KIND of claim, and the decision to treat it as a separate
category rather than sliding into it is deliberate (user call, 2026-07-28).

## The mechanism

Every test function used through C3 vanishes on the whole patch boundary, which
is what kills the boundary terms in the by-parts identity. Drop that **in time
only** — a φ that is nonzero at the window's initial time — and

    ∫∫ φ u_t = [∫ φ u dx]_{t0}^{T} − ∫∫ φ_t u

so with φ(·, T) = 0 the initial condition enters as a linear functional with an
ANALYTIC weight, inside the same linear system, under the same declared band.
With the law KNOWN, the system rearranges to

    Σ_j a_j B_ij = y_i,   B_ij = ∫ φ_i(x, t0) b_j(x) dx

where the design matrix B is exact — known analytic constants carrying no data
error at all — and only the target y carries the weak-form band. The band is
SIMPLER here than in the law case, not harder.

## What a state certificate says

- **Claim.** Over the stated observation window and patch family, every initial
  condition whose mode coefficients lie in the reported intervals reproduces the
  observations within the declared band. Modes outside the reported set are NOT
  determined — they are named, not silently dropped.
- **Domain.** The observation window and the basis. Nothing about other times,
  other solutions, or modes above the reported cut.
- **α.** With a fixed basis and a known law there is no search: |H| = 1, so α is
  a pure chance-agreement bound rather than a multiple-comparisons-corrected
  one. Defensible, but a DIFFERENT quantity from a law certificate's α, and
  labelled `alpha_kind` so the two are never compared.
- **Resolution bound, stated up front.** dof = number of basis modes, so
  h = n_rows − dof must stay positive with margin (`MIN_HELDOUT = 8`). At most
  as many modes can be certified as there are independent patch equations —
  the 35-term-interpolation lesson, in the registration rather than the results.

Two refinements the implementation added, both reported per mode:

- **DETERMINED vs RESOLVED.** Determined means the interval is bounded;
  resolved means it EXCLUDES ZERO, i.e. the certificate can tell the mode is
  there at all. A bounded-but-zero-straddling interval is a real result and must
  not read as a recovered mode.
- **Joint vs conditional intervals.** The reported interval is the exact
  projection of the whole feasible set `{a : |Ba − y| ≤ ε}` onto that
  coordinate, by two linear programs. `certify.parameter_interval`'s bisection
  holds every other mode FIXED and answers a conditional question; both are
  reported, and where they differ the joint one is the claim.

## Registered predictions

- **W1.** Advection: every mode in the reported set certifies, half-width flat
  in k within a factor of 2, and the true amplitude inside every interval.
- **W2.** Heat: half-width(k)/σ tracks `exp(+νk²T)` within a factor of 3 up to
  the cut; above the cut the mode is reported UNDETERMINED rather than estimated.
- **W3.** The cut moves the way the physics says: later T or larger ν pushes it
  to lower k, smaller σ pushes it up.
- **W4.** Burgers after shock formation: refusal, not a wide interval — the
  information is destroyed, not merely diluted.
- **W5.** Zero confident-wrong: no reported interval excludes the true amplitude.

## Results

Run 2026-07-28 (`experiments/pde/run_c4_state.py`, module `lagh/statecert.py`,
results `experiments/results/pde_c4.json`). Fourier basis {1, cos kx, sin kx},
one-sided-in-time patches at three pooled scales, NX = 256, NT = 161.

- **W1 — met.** Advection certifies every mode at σ ∈ {1e-6, 1e-5, 1e-4}, with
  half-widths flat in k (0.62σ – 0.97σ at σ = 1e-4, a spread of 1.55× across
  k = 1…8, against the registered factor of 2) and every true amplitude inside
  its interval. The eight modes actually present come back
  RESOLVED; the eight sine modes, which are absent, come back determined but NOT
  resolved — their intervals straddle zero, which is the correct statement.
- **W2 — met, but only after the prediction's own geometry was corrected.**
  A window that starts at t = 0 sees every mode at full amplitude, so its
  resolution is FLAT in k — the registered `exp(+νk²T)` curve does not appear
  there, and reporting it as absent would have been a misreading. The
  exponential ill-posedness is entirely in the PROPAGATION: certify the state at
  a later t₀ and carry it back to t = 0 through the known (diagonal) law, and
  the interval is divided by `exp(−νk²t₀)` exactly. Measured at ν = 0.3,
  t₀ = 0.3, σ = 1e-4: half-widths run 4.0e-4 (k=1) → 1.1e-1 (k=8), a ratio of
  **278** against the predicted `exp(ν t₀ (k²−1))` = **290** — inside 1.05×,
  well within the registered factor of 3.
- **W3 — met in direction.** `k_cut` (the largest k whose back-propagated
  interval excludes zero) falls as ν·t₀ grows and as σ grows, and rises as σ
  falls; with the basis capped at k = 12 the cut sits at the basis edge for the
  gentler rungs, so the measurement bounds the trend rather than tracing the
  full `√(ln(1/σ)/(νT))` curve.
- **W4 — met, with the mechanism named, and only after the field family was
  fixed.** The Cole–Hopf family used elsewhere in the arc CANNOT form a shock:
  with a fixed positive potential its initial profile is ν times a fixed shape,
  so lowering ν scales the field down instead of steepening it (measured: the
  basis-truncation error is identical at ν = 0.02, 0.01 and 0.005). Testing a
  shock prediction on that family would have tested nothing. On a field that
  really steepens (u₀ = −sin x, ν = 0.02 and 0.005, spectral with dealiasing and
  a declared solver error): a window starting before the front certifies; a
  window starting AFTER it loses 36 of 90 patches to the resolution gate — the
  grid no longer represents the field — leaving fewer rows than unknowns, and
  the certificate refuses with `resolution`. A refusal, not a wide interval.
- **W5 — met.** No reported interval excludes a true amplitude, in any rung, in
  either the forward or the back-propagated form.

### The refusals, and what each one names

Across 41 rungs there are 7 abstains and every one of them names a mechanism.
Four are Cole–Hopf Burgers rungs whose DECLARED BASIS cannot represent the
state, and the truth check proves it is the basis and not the instrument:

| rung | verdict | truth/band | basis truncation |
|---|---|---|---|
| ν = 0.2, k ≤ 8 | abstain | 15.9 | 1.1e-4 |
| ν = 0.2, k ≤ 16 | **certifies** | 0.14 | 6.7e-9 |
| ν = 0.05, k ≤ 8 | abstain | 633 | 3.0e-3 |
| ν = 0.05, k ≤ 16 | **certifies** | 0.78 | 6.8e-6 |
| ν = 0.02, k ≤ 16 | abstain | 354 | 2.9e-3 |
| ν = 0.02, k ≤ 32 | **certifies** | 0.41 | 5.7e-6 |

The Cole–Hopf profile is a RATIO, so its Fourier spectrum is infinite and a
truncated basis is a modelling choice with a measurable cost. Declaring a bigger
basis fixes it — until the resolution bound says it may not be declared, which
is the other two abstains (the shock rungs) and the over-parameterized k ≤ 40
rung (81 unknowns against 30 rows).

### The verify track

Every certified state was integrated forward from its recovered amplitudes under
the known law, with the per-mode intervals propagated as an envelope (per-mode
triangle bound) and the solver's own tolerance-ladder error declared. Zero
points outside the band across all rungs (41216 points per rung).

### Two measured facts about the instrument, both worth not rediscovering

**The band has a quadrature floor, and the floor is visible.** A one-sided
window does not vanish at t₀, so it loses the Euler–Maclaurin cancellation that
makes the two-sided rule spectrally accurate. The rule there is composite
trapezoid, Richardson-extrapolated in TIME ONLY (extrapolating across a joint
space–time refinement corrects with junk — a 21-point space window subsampled to
6 points is not a quadrature — and that mistake moved the u_t integral by 1.3%
and rejected every state patch), with the space direction carrying its own
separate, non-extrapolated bound. What remains is a declared O(h⁴) term, and at
a coarse time grid it DOMINATES the band:

| NT | half-width/σ at σ = 1e-6 | at σ = 1e-4 |
|---|---|---|
| 81 | 17 – 32 | 1.8 – 2.9 |
| 161 | 2.4 – 4.2 | 1.2 – 1.9 |
| 321 | 1.0 – 1.9 | 0.8 – 1.4 |

So the small-σ rungs are quadrature-limited until the time grid resolves the
window, and only then does resolution become linear in σ. Reporting the σ = 1e-6
rung at NT = 81 as "the instrument's resolution" would have been an artifact of
the grid.

**Least squares is the wrong solver for a per-row band.** The plain lstsq fit
violated the band at 9.8× on some rows while the TRUE amplitudes sat at 0.97×:
the band spans a decade across patch scales, and an L2 fit happily trades a
loose row's slack for a tight row's violation. The certificate's question is
feasibility — is there an initial condition in the declared basis consistent
with the observations? — so it is answered by a linear program, and the
per-mode intervals are LP projections. The program must also be scaled (rows by
the band, columns by their norm): unscaled, the solver failed on exactly the
hard rungs, the code fell back to least squares, and a state whose truth sat at
0.57 of its band came back as "no state explains the observations". A spurious
refusal is still a wrong answer.

### Open

Non-periodic domains (the interval search needs a basis where the forward
operator is diagonal, which is Fourier on a periodic domain — anything else is a
later problem, and it is stated rather than discovered); joint law-and-initial-
condition inversion, explicitly out of scope for this pass and deserving its own
registration; and state certificates for SYSTEMS, where the boundary term
carries one initial condition per field.
