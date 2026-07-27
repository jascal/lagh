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

The full program narrative, measured soundness results, and the taxonomy of
failure modes: [docs/INSTRUMENT_REPORT.md](docs/INSTRUMENT_REPORT.md).

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
