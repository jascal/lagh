# PDE dev campaign C1 — the noise ladder — registration

**Registered 2026-07-28, before the ladder is run.** C0
(`CASE_STUDY_PDE_DEV.md`) certified weak-form PDEs on exact analytic solutions
and deliberately deferred σ > 0, because the design matrix columns are noisy
functionals of the same field as the target — an errors-in-variables problem the
ε model could not express. C1 is that case, and it required three engine
changes, all built before this registration and each traceable to a measurement
made while building them.

## The three engine changes C1 rests on

1. **Per-candidate ε** (`certify.band`, `engine.discover(eps_model=)`). ε may
   now be a CALLABLE assembled for the law under test. The weak-form residual is
   `δ_y − Σ c_k δ_k` where every δ is a functional of ONE noise realization, so
   its band is a single norm `κσ‖ν_y − Σ c_k ν_k‖₂` over the combined
   sensitivity vector — computed as the quadratic form `a'Ga` from a per-patch
   Gram matrix (`weakform.PatchEpsilon`). The C0 band had to BOUND the
   coefficients instead (`coeff_max`); measured, the per-candidate band is
   **6.9× tighter** on heat rows and covered the true law's residual in 960/960
   patch-draws. Every gate that re-checks a modified law — `float_pinned`'s
   perturbations, `reduce_to_minimal`'s drops — now gets the band of the law it
   is actually checking.
2. **Noise-corrected resolution gates** (`weakform`). The patch ladder and the
   aliasing test are computable functionals of the declared noise, so both now
   subtract what σ explains before judging resolution. Without it the gates
   rejected patches for carrying noise ε already bands: **16 of 24 heat patches
   lost at σ = 1e-4, all 24 at σ = 1e-2**; with it, 24/24 and 23/24.
3. **Interval-parameter certificates** (`certify.parameter_interval`,
   `engine.discover(declared_basis=)`). Two findings forced this:
   - the under-noise dense-channel skip (empirical-only, PO12/PO40) makes a
     weak-form library uncertifiable at any σ > 0, because the PDE terms come
     through the linear channel. That skip's justification is about UNREGISTERED
     bases; `declared_basis=True` is the caller declaring a registered term
     library, and it is opt-in, default off.
   - with the channel open, the fitted coefficient arrives snapped to
     `90360785/903607052` for a true `1/10` — a float in a costume — and the
     exact-value gate then abstains. But an exact rational is the right claim for
     a definitional identity and the WRONG one for a diffusivity. The engine now
     measures what the data determines: the interval over which the law still
     certifies. Measured on heat rows: ν = 0.1 ± 1.5e-8 / 2.1e-6 / 2.6e-4 at
     σ = 1e-8 / 1e-6 / 1e-4 — linear in σ, always containing the truth. An
     unbounded interval still abstains.

**Already measured before this registration** (disclosed, so nothing below is
scored as a prediction it isn't): the heat interval widths above, the 6.9×
tightening, the patch-survival counts, and one end-to-end heat run at σ = 1e-6
which certified `u_t = 0.09999996·u_xx` with the interval [0.09999690,
0.10000302]. Advection, Burgers, the nulls, and every other rung are unrun.

## Registered predictions

Ladder: σ ∈ {0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3} on the same analytic fields
as C0 (heat, advection, Burgers), pooled over 4 initial conditions, certifying
on a held-out solution, patch family A (p = 16).

- **M1 — zero confident-wrong across the whole ladder.** Wherever a law
  certifies, its SUPPORT is the true one (no extra terms, none missing) and
  every certified parameter interval CONTAINS the true value. This is the claim
  the arc lives or dies by.
- **M2 — graceful degradation, not a cliff.** As σ rises each system passes
  through: exact-rational certificate → interval certificate → abstain. No
  system jumps from certifying to wrong.
- **M3 — interval width scales linearly in σ**, within a factor of 3 across the
  ladder (it is a linear functional of the noise; anything else means the band
  is not doing what it claims).
- **M4 — Burgers degrades first.** Its nonlinear term contributes least to the
  target, so it loses certification at a lower σ than heat or advection.
- **M5 — nulls stay silent at every rung.** The i.i.d. and smooth-random nulls
  certify nothing at every σ, and their patches are rejected or abstained, never
  certified.
- **M6 — the single-solution refusal survives noise.** A single traveling wave
  still returns `single-solution` at every σ (structural, not luck).

## What would end the arc

Any certified law at any rung whose interval excludes the truth, or whose
support is wrong (M1). Unlike C0 this is not guaranteed by construction: the
band is wider, more candidates certify (measured: 1816 of 2118 candidate checks
passed at σ = 1e-6), and coherence is carrying more weight than before.

## Results

Run 2026-07-28 (`experiments/pde/run_c1.py`, results
`experiments/results/pde_c1.json`). Three systems × seven σ, pooled over 4
initial conditions, certifying on a held-out solution.

| σ | heat (ν = 1/10) | advection (c = −7/10) | Burgers (ν = 1/5, −1) |
|---|---|---|---|
| 0 | exact `1/10` | exact `−7/10` | exact `1/5`, `−1` |
| 1e-8 | interval ±2.5σ | interval ±5.0σ | intervals ±6.5σ, ±93σ |
| 1e-7 | ±3.06σ | ±5.28σ | ±6.45σ, ±92.7σ |
| 1e-6 | ±3.06σ | ±5.28σ | ±6.45σ, ±92.7σ |
| 1e-5 | ±3.06σ | ±5.28σ | ±6.45σ, ±92.7σ |
| 1e-4 | ±3.06σ | ±5.28σ | ±6.45σ, ±92.7σ |
| 1e-3 | **abstain** (structural) | ±5.28σ | **abstain** (structural) |

- **M1 — met, and it is the result that matters.** Every certifying rung has
  the TRUE support (no term added, none missing) and **every reported parameter
  interval contains the true value** — 0 violations across 18 certificates and
  25 intervals. Zero confident-wrong on the whole ladder.
- **M2 — met.** Every system walks exact → interval → abstain. No system jumps
  from certifying to wrong; the two that stop certifying do so structurally
  (rival classes at that band), which is the honest verdict, not a cliff.
- **M3 — met, more sharply than registered.** Interval half-width / σ is
  CONSTANT per system to three significant figures across four decades: heat
  3.06, advection 5.28, Burgers 6.45 (ν) and 92.7 (the u·u_x coefficient). The
  band is exactly the linear functional it claims to be. The ratio also reads
  as an identifiability constant per coefficient: Burgers' nonlinear
  coefficient is ~14× less determined than its diffusivity by the same data.
- **M4 — MISSED.** Burgers was predicted to degrade first; it and heat both
  abstain at σ = 1e-3 while ADVECTION survives the whole ladder. The reason is
  visible in M3: what decides survival is the per-coefficient identifiability
  constant, not whether the system is nonlinear. Advection's single coefficient
  is the best-determined of the three.
- **M5 — met.** Nulls certify nothing at every σ: the i.i.d. field loses every
  patch at the resolution gate, the smooth random field abstains structurally.
- **M6 — met.** The single traveling wave returns `single-solution` at σ = 0
  and σ = 1e-4 — structural, as designed, not luck.

### What the ladder exposed in the engine (both fixed here)

- **A modest-denominator rational is an exactness claim too.** At σ = 1e-3 the
  advection winner was `−72071/102974`-shaped: the float gate only inspects
  Floats and denominators > 10⁶, so it passed with a value 1.4e-4 off the truth
  and the certificate said nothing about its own precision. Under a declared
  basis the interval is now reported for EVERY free parameter, whatever costume
  it arrives in — the exact/interval boundary is now set by the data rather
  than by `float_pinned`'s σ-scaled perturbation.
- **A single-scale patch family makes a degenerate column.** With every patch
  the same size, the `1` term's column is exactly constant, and the
  constrained-input machinery correctly detects it as a machine-exact input
  constraint and switches to the domain-restricted path. Nothing was
  mis-certified, but the mechanism is an artifact of the patch family, not
  physics: a multi-scale family is the fix, and it is the next change (also the
  right move for identifiability, since the C0 registration already called for
  "multiple scales").

## C1b — multi-scale patch families (registered 2026-07-28, before the re-run)

The C1 ladder ran on a single-scale family, which makes the `1` term's column
exactly constant (measured spread 2.2e-16, i.e. 1.3e-14 relative) — a genuine
machine-exact input constraint that the constrained-input machinery correctly
detected, switching the engine to its domain-restricted path at two rungs. The
degeneracy is an artifact of the patch family, and pooling scales removes it
(`weakform.multiscale_patches`): measured, the same column then spans 75% of its
own magnitude with signal-to-band essentially unchanged (1.05e7 vs 1.09e7).

Scales pooled: (16, 8), (24, 12), (32, 16) half-widths in grid cells, 12 patch
centres each, 36 rows per solution against the previous 24.

- **S1 — the domain-restricted path stops firing.** No input constraint is
  detected at any rung (it fired at heat σ = 1e-5 and Burgers σ = 1e-4).
- **S2 — no regression.** Every rung that certified the true support under the
  single-scale family still certifies the true support.
- **S3 — the interval story survives.** Half-width/σ stays constant in σ per
  system (the linearity result), and the constants move by less than 2× in
  either direction. Adding rows adds constraints, so I expect tightening rather
  than widening, but the pooled small patches carry less signal per row, which
  is why the prediction is a bound and not a direction.
- **S4 — zero confident-wrong.** Every reported interval contains the truth.

### C1b results: three attempts, and the two failures are the content

**Attempt 1 — scales (16,8), (24,12), (32,16): every heat and advection rung
went to a structural abstain, including σ = 0.** Diagnosis: the (16,8) patches
carry a band ~2000× looser than the others (signal-to-band 5.3e3 vs 1.05e7,
because a smaller support means fewer points across the bump and a worse
quadrature). Rivals certify at the loosest row's band, so **a pooled patch
family is only as sound as its worst-conditioned member** — the loose-ε lesson,
one level up.

**Attempt 2 — resolved scales only, (24,12), (32,16), (40,20): heat still
abstained, now with 21 materially different certifying classes** on clean data.
Diagnosis, and it was a modelling error of ours rather than an engine one: with
multiple scales the `1` column VARIES, so it is a genuine input, and the general
library duly built laws like `u_xx·[1]^(3/2)`, `u_xx·√[1]` and `u_xx/[1]`.
Because `[1]` takes exactly one value per scale, any function hitting the right
value at three points fits — the patch scale had become a state variable, which
a PDE law must never depend on.

**Attempt 3 — normalize each row by its own ∫φ and drop that column
(`WeakSystem.normalize`).** Every column becomes a patch AVERAGE, commensurable
across scales; the linear relation is unchanged (row scaling cancels); a source
term becomes the intercept the engine already proposes; bands scale with their
rows (the Gram by the square). 144 rows, 4 features.

| σ | heat | advection | Burgers |
|---|---|---|---|
| 0 | exact `1/10` | exact `−7/10` | exact `1/5`, `−1` |
| 1e-8 … 1e-4 | ±0.89σ | ±1.77σ | ±2.67σ, ±67.3σ |
| 1e-3 | abstain | abstain | abstain |

- **S1 — met.** No input constraint is detected at any rung; the
  domain-restricted path no longer fires.
- **S2 — MISSED at one rung.** Advection at σ = 1e-3 certified under the
  single-scale family and now abstains structurally. Everywhere else the true
  support still certifies, and heat/Burgers/advection now all reach σ = 1e-4.
- **S3 — MISSED on the bound, met on the substance.** Half-width/σ stays
  constant in σ (the linearity result survives), and the constants TIGHTEN, as
  expected — but by more than the 2× I registered: heat 3.06 → 0.89 (3.4×),
  advection 5.28 → 1.77 (3.0×), Burgers 6.45 → 2.67 (2.4×) and 92.7 → 67.3
  (1.4×). More rows, better conditioned, and the scale nuisance removed.
- **S4 — met.** Every interval at every rung contains the truth; zero
  confident-wrong across the re-run.

Net: same verdict shape as C1, one rung lost, parameters determined 1.4–3.4×
more tightly, and the constant-column artifact gone.

### Open

Multi-scale patch families; a σ-ladder on the harder curriculum (KdV,
variable coefficients); the verify track; and the standing question of whether
interval-parameter certificates should be the general default under declared
noise rather than a declared-basis-only path.
