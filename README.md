# lagh

*(Scottish Gaelic: law)*

A **certified symbolic law discoverer**: given data from a system with an exact
underlying law, lagh returns either a symbolic law with a machine-checked
certificate — carrying a stated significance bound α — or a machine-readable
refusal saying why no law is certifiable. **Never a confident wrong answer.**

## Headline results

- **Gravity-Bench-v1 (ICML'25), budgeted variant: 94.7% vs 49% LLM-agent SOTA**
  — a fully deterministic observation-planning agent (no LLM calls), one-shot
  blind read, pre-registered protocol
  ([report](docs/BLIND_READ_REPORT_GRAVITYBENCH.md)).
- **Zero confident-wrong across ~600 scored tasks** spanning two data regimes
  (active/oracle and passive/fixed-dataset) and two public benchmarks.
- **Significance-bearing certificates**: every certification carries
  α ≤ |H|·q^h, null-validated at 0 false certifications / 200 true-random
  targets. Refusals are certified randomness statements relative to a declared
  hypothesis class.
- An honest loss, fully reported: LLM-SRBench under the frozen blind protocol
  ([report](docs/BLIND_READ_REPORT.md)) — whose diagnosis (float32-quantized
  benchmark data; representation precision is *declared noise*) became
  standing policy.

- **Twenty certificates from four real-archive campaigns**
  ([docs/CERTIFICATES.md](docs/CERTIFICATES.md)): Gaia DR3 (photometric
  zeropoints with the exact 5/2 slope; the DR3 RRab metallicity pipeline
  formula, rediscovered with Nemec's exact coefficients, α ≤ 10⁻¹⁰⁰⁹), the
  NASA Exoplanet Archive (the composite density formula at α ≤ 10⁻¹⁹¹⁸,
  its constant decoding the Archive's Earth-radius convention), the
  Materials Project (`volume/nsites` at **α ≤ 10⁻⁵⁰⁷⁰**; the density
  identity whose refusal residuals reconstructed the pipeline's
  2005-vintage atomic-mass table element by element, then certified with
  the atomic mass constant read out to nine digits), Kepler III three
  times over from NSSDC tables, and the GDP accounting identity from raw
  FRED columns.
- **A registered theory test on 25 GiB of rebuilt DFT** (catalysis
  campaign, [docs/CASE_STUDY_CATALYSIS_C1.md](docs/CASE_STUDY_CATALYSIS_C1.md)):
  bond-counting's exact rational scaling slopes, bands frozen in git
  before any slope was computed — confirmed on the carbon ladder
  (CH₂: 0.531 ± 0.016 vs 1/2; CH₃: 0.225 ± 0.014 vs 1/4), broken for
  OH (+0.21) and SH (+0.40) on bimetallic alloys; zero certificates
  granted on scattered chemistry, exactly as registered.
- **Refusals that decode into findings**: Saturn's Kepler abstain is J₂
  oblateness plus Hyperion's resonance, quantitatively; Gaia's
  Stefan–Boltzmann abstain is FLAME's marginal medians failing to compose;
  the frame-rotation abstain is the unit-sphere constraint on direction
  cosines; the exoplanet abstains map PSCompPars's literature patchwork;
  the Materials Project census pulled eight materials with impossible elastic
  bound ordering (Reuss > Voigt) out of the database by ID, plus one with an
  implausible magnitude (reported upstream, emmet#1499 — where the maintainer
  corrected our reading: those entries ARE flagged, on the `elasticity`
  collection whose per-property `deprecated` we had not read, not on the
  aggregated `summary` endpoint we swept).

The full program narrative, measured soundness results, and the taxonomy of
failure modes: [docs/INSTRUMENT_REPORT.md](docs/INSTRUMENT_REPORT.md).

## The failure → hardening ledger

Every protection in the engine is traceable to a measured failure — none is
speculative. The ledger, in the order the field taught it:

| measured failure | hardening |
|---|---|
| a coefficient 3×10⁻⁵ off the truth certified at loose ε (reverted-transform output scale) | the exact-coefficient gate: perturbing every coefficient must break certification (`float_pinned`) |
| dyadic-garbage overfits certified at floor-dominated ε and poisoned coherence | per-candidate gating on clean data; gate moves to the winner under noise |
| quantized (float32) benchmark data rejected every true multi-term law | representation precision is *declared noise* (σ_rep); prefilters may not out-tighten ε |
| a 14-term approximant certified on Gaia photometry at a loose declared floor — the gate had *inverted* (it rejects the honestly-unpinned truth, passes delicately-canceling whales) | the floor-dominated regime: refit-parsimony collapse (`refit_minimal`), winner-level gate, sticky re-split ambiguity |
| the DR3 metallicity quadratic was unreachable — size-5 supports existed only via fragile top-6 singles | CAP-T: complete quadric/cubic polynomial supports proposed unconditionally |
| a 35-term interpolation of 12 binned points "certified" while its own α bound said 1 | **significance is part of certification**: α ≤ 10⁻⁶ or the certificate demotes to an abstain |
| a certified-looking result at an amended floor was split luck — one re-split dissolved a genuine 3-class ambiguity | structural ambiguity is sticky across re-splits |
| a constant over 4 overflow-artifact points was the program's only confident-wrong | the minimum-domain guard (evidence floors under \|H\|·q^h) |
| the IAU frame rotation abstained: rival classes equal the true law modulo the unit-sphere constraint the inputs satisfy — the probe box asked an off-manifold question the domain claim never made | constrained-input coherence: machine-exact constraint detection (SVD null-space, 10⁻¹⁰ tolerance), manifold-probe coherence, and winner canonicalization modulo the constraint ideal — the frame rotation now certifies at α ≤ 10⁻⁷⁶⁸⁰ with the constraint named in the certificate |
| a forward forecast from a certified PDE reproduced the field everywhere except a handful of tail points, and only on the non-dissipative system — while the law itself was exact to 3e-11 | a forecast started from MEASURED data carries that data's noise for the whole window, and without dissipation it never damps (measured 2.8σ); the verify band now measures that term by re-integrating from perturbed initial conditions, and "the law reproduces the field" is reported separately from "the law predicts measurements" |
| pooling patch scales to remove a degenerate column made every heat and advection rung abstain — the smallest patches carry a ~2000× looser band and rivals certify at the loosest row's — and once the scales were comparable, a varying `∫φ` column turned the TEST-FUNCTION SCALE into a state variable, so 21 rival classes fit laws like `u_xx·[1]^(3/2)` on clean data | a pooled family is only as sound as its worst-conditioned member (scales must be comparably resolved), and rows are normalized by their own `∫φ` so every column is a patch average and a source term is the intercept, not a scale-valued input |
| under declared noise a weak-form coefficient certified as `−72071/102974` — an exactness claim the data could not support (1.4e-4 off the truth), because the float gate only inspects Floats and denominators > 10⁶ | **interval-parameter certificates**: where an exact value is not what the data determines, measure what IS — the interval over which the law still certifies — and claim that; reported for every free parameter, and an unbounded interval still abstains |
| a weak-form PDE system built from ONE solution certified `u_t = −u_x` on a KdV soliton — the on-shell traveling-wave relation, true of that field and not the equation it came from | a PDE claim must certify on patches from a **held-out solution**; single-solution data refuses structurally (`abstain: single-solution`) instead of certifying whatever the library's simplest fit happens to be |
| Müntz twins at machine floor abstained three reach cells; arbitrating them by α margin alone then **certified an approximant** on `rational-d1` (2×10⁻⁵ off the truth just outside the box) — α ranks dof, and the truth was not even in that contest | significance arbitration with an **evidence bar**: a rival may be dismissed only when it is an interpolation of the sample (held-out fraction h/n < 0.10), never when it retains evidence — reach 33/36 → 35/36, and the genuinely-twinned cell keeps its abstain |
| a coupled PDE system with a 5-term true support over 8 declared columns proposed **nothing at all** — STLSQ, OMP and size-5-from-top-6-singles all missed it — while the truth sat four orders inside its own band | `discover(linear_basis=True)`: when the caller declares that the claim is a linear combination of the declared columns (which a weak-form PDE always is), supports are enumerated exhaustively over the library and the engine builds no products or powers of patch integrals |
| single-mode initial data made a coupled system certify `u_t = −u/10 + v/2` — because `u_xx ≡ −u` holds across the WHOLE solution family, so the multi-solution holdout cannot break it | support equality is not the test: score agreement with the truth **over the certified domain** separately from support, and reserve confident-wrong for a law that disagrees beyond its band. The constrained-input path already reported it as a domain-restricted certificate; what was missing was a scorer that could tell the two apart |
| a state certificate's least-squares fit violated the band at 9.8× while the TRUE amplitudes sat at 0.97× — the per-row band spans a decade across patch scales, and L2 trades a loose row's slack for a tight row's violation | a certificate's question is **feasibility**, so it is answered by a linear program (and the per-mode intervals are exact LP projections of the feasible set); the program is scaled by the band and the column norms, because unscaled it failed on exactly the hard rungs and a state whose truth sat at 0.57 of its band came back as "no state explains the observations" |

One registered open boundary, honestly held: the approximant-impostor class
under declared noise (dense channels are empirical-only by design —
|H|·q^h accounting is the eventual road through). A second, newly stated:
**constrained twins** — rivals that each retain real held-out evidence and
diverge only outside the sampled box — are indistinguishable under this
instrument's evidence, and the abstain there is permanent, not a gap
(`MUNTZ_ARBITRATION.md`).

## Error-model repertoire (hard-won, both directions)

Real archives taught the ε model its full shape: relative representation
noise (σ_rep) for float32-rounded columns (Gaia RUWE-class columns, the
`_over_error` ratios, archive 3-significant-digit rounding); absolute
floors for fixed-decimal storage (Gaia magnitudes, MP's 3-decimal GPa
moduli — the same lesson, inverted); per-point propagated `se` for derived
columns whose input rounding amplifies through ratios (the anisotropy
formula's extreme-ratio rows); and reference-offset-free ΔE\* spaces where
unknown per-species constants provably cannot touch the slopes under test
(the catalysis rebuild, whose deposit ships no gas references).

## Use

```python
from lagh import discover                      # core engine (split-level API)
from lagh.passive import discover_passive      # fixed datasets (benchmarks)
from lagh.submit import submission             # two-track benchmark submission

r = discover_passive(X, y, sigma=0.0)
r.certified                    # True -> exact law with certificate
r.result.expr                  # sympy expression, exact rationals
r.result.certificate.alpha_log10   # log10 of the chance-fit bound
```

`lagh/mcp/` exposes the three verbs over MCP: `recover` (bounded — certificate
or reasoned abstention), `verify` (check a declared form), `fit` (an unbounded
scout whose output is labeled conjecture, never certificate).

**SRBench**: `lagh.sklearn.LaghRegressor` implements the SRBench estimator
contract (`fit`/`predict`/`model()`, `max_time`, `random_state`); the
ready-to-PR method directory is `srbench-submission/algorithms/lagh/`. The
conjecture track always answers; certificate status rides along as estimator
metadata (`track_`, `tag_`, `alpha_log10_`).

## Layout

| path | what |
|---|---|
| `lagh/certify.py` | the fixed honesty core: exhaustive check, ε model, vacuity, coherence, pinning, coefficient gate, minimality repair, significance |
| `lagh/classes/` | the registered curriculum C1–C9 |
| `lagh/engine.py` | MDL-ordered escalation; escalates only on an empty certifying set |
| `lagh/acquisition.py` / `lagh/passive.py` | active and passive regimes |
| `lagh/weakform.py` | PDE claims in weak form: patch integrals of the RAW field against analytic test-function derivatives, each with a declared bound (no differentiated data ever enters a certificate). Multi-field terms and n-dimensional geometry, so systems and 2-D fields use the same factory |
| `lagh/pdesystem.py` | systems of PDEs: one row set, one target per equation, features spanning every field, and a conjoined certificate whose α is a union bound dominated by its weakest equation |
| `lagh/statecert.py` | state certificates — inverting for an initial condition, a claim about one system's history rather than a law, with UNDETERMINED as a first-class per-mode outcome |
| `machine/` | the research loop as a verified Orca state machine; optional bounded LLM proposer |
| `experiments/gravitybench/` | the deterministic astronomer (planner + digital twin) |
| `docs/` | registrations, blind-read reports, direction docs — every claim's provenance |

## Method discipline

Every claim is tagged `proved` / `empirical` / `open`. Capabilities land only
through pre-registered predictions scored against named cells; blind reads
freeze SOTA and protocol before download, run once, and report win, lose, or
mixed. Both practices are documented in `docs/` with the artifacts inline.

## Development

```bash
.venv/bin/pytest tests/ -q
```

## License

Apache 2.0 — see [LICENSE](LICENSE).
