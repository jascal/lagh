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

- **Fourteen certificates from real public archives**
  ([docs/CERTIFICATES.md](docs/CERTIFICATES.md)): Gaia DR3 (photometric
  zeropoints with the exact 5/2 slope; the DR3 RRab metallicity pipeline
  formula, rediscovered with Nemec's exact coefficients, α ≤ 10⁻¹⁰⁰⁹), the
  NASA Exoplanet Archive (the composite density formula at **α ≤ 10⁻¹⁹¹⁸**,
  its constant decoding the Archive's Earth-radius convention), Kepler III
  three times over from NSSDC tables, and the GDP accounting identity from
  raw FRED columns.
- **Refusals that decode into findings**: Saturn's Kepler abstain is J₂
  oblateness plus Hyperion's resonance, quantitatively; Gaia's
  Stefan–Boltzmann abstain is FLAME's marginal medians failing to compose;
  the frame-rotation abstain is the unit-sphere constraint on direction
  cosines; the exoplanet abstains map PSCompPars's literature patchwork.

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

Two registered open issues, honestly held: constrained-input coherence (the
IAU frame rotation cannot certify while the probe box leaves the data
manifold — every rival class equals the true law modulo the constraint) and
the approximant-impostor boundary under declared noise (dense channels are
empirical-only by design).

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
