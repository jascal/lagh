# lagh MCP server — `recover` / `verify` / `fit`

The tool-shape from [`docs/DIRECTION_TOOLSHAPE.md`](../../docs/DIRECTION_TOOLSHAPE.md),
made runnable. An LLM parses a problem, samples the black-box oracle, and calls these
three tools; lagh returns a **law checked over a stated finite domain or a
reasoned abstention**. Certification is conditional on the declared error model
and checker assumptions. The `proved` tag records passage through those checks,
not universal truth or a guarantee against checker defects. See the
[falsifiability audit](../../docs/FALSIFIABILITY_AUDIT.md) for measured failures.

## Install & run

```bash
pip install -e ".[mcp]"     # the [mcp] extra pulls the MCP SDK; core logic needs neither
lagh-mcp                    # or: python -m lagh.mcp.server
```

`import lagh.mcp` (the tool *logic*) works without the SDK — only the transport needs it.
Agent discovery config is in [`.mcp.json`](../../.mcp.json).

## The three acts

Named for what they **do**, never for the claim they make (Orca convention:
`verify_machine`, never `certify_machine`). Certification is the typed *return value*,
with a `proved` / `open` tag — not a tool name.

| tool | bound | act | returns |
|---|---|---|---|
| **`recover(X, y, sigma)`** | bounded | discover an exact law | `Certificate` \| `Abstention` |
| **`verify(X, y, form, sigma)`** | bounded | check a *declared* form (restricted math grammar, never evaluated as Python) | `Certificate` \| `Abstention` |
| **`fit(X, y, sigma)`** | **UNBOUNDED** | best-guess + identifiability diagnosis | `Conjecture[]` + `Diagnosis` — **no `certified` field** |

- **Certificate strength** (`recover`/`verify`): `pinned` records passage through
  the clean-data operational gates; it does not prove the generating form.
  `consistent` covers declared noise, loose floors, other declared error, and
  declared irrational constants. The machine-readable `claim` states the scope.
  Neither strength implies exact structure or exact coefficients from finite data.
- Parameter slices, when supplied in `partial`, hold the other parameters fixed.
  They are not marginal uncertainty intervals or a simultaneous parameter box;
  read their `scope` before using them. A missing bound can be a search limit.
- **The wall is structural.** `fit`'s result has no `certified` key — a guarantee
  cannot be read off it by accident. `fit` is a scout; to make a conjecture
  trustworthy, feed it to `recover` or `verify`.

## Active acquisition (two modes)

lagh's strength is *adaptive* querying — ranging to find the informative box, then
budget-metered multi-objective queries. That needs a **live oracle**, which JSON can't
carry, so `recover` has two modes:

- **In-process (library):** hand lagh the oracle and box and it drives the whole loop:
  ```python
  from lagh.mcp import recover
  recover(oracle=lambda X: 3*X[:,0]**2, box=[[0.5],[4.0]])          # → certificate + `acquisition` provenance
  recover(oracle=my_oracle, box=[lo, hi], box_search=True)          # broaden-the-box ladder on abstain, held-out-box guarded
  ```
  The certificate comes with `acquisition` = {queries_used, budget_spent, box_final, …}.
- **Over the wire (caller-driven):** the data-only tool can't call back to your oracle,
  so on a thin/under-determined abstain it returns `next_action:"acquire"` +
  `suggested_box` (10× wider). Re-sample that box, call `recover` again — the acquisition
  loop, driven by you.

Structural refusals can also carry `design_evidence`: an empirical list of rival
expressions retained across discovery splits. These are candidates for planning
a distinguishing measurement, not certified output laws or an exhaustive list
of explanations. The truth need not be among them. Once used to choose a query,
all observations behind these rivals are design evidence; a subsequent claim
requires fresh certification observations. This metadata does not change the
fixed `suggested_box` or automatically query a new measurement.

## The loop `fit` is built for

`fit`'s real product is the **diagnosis**, not the guess. It names the next move:

```
fit(X, y)  →  diagnosis: "pinned to 5/2"          →  recover(X, y)          # get the certificate
           →  diagnosis: "continuum ~2.718 (≈e)"  →  verify(X, y, "x_0**E") # → consistent
           →  diagnosis: "under_determined"        →  acquire data, retry
```

so the scout → declare → check loop is driven by the caller, guided by `next_action`.

## Honest boundaries (carried in the tool docstrings)

- The certificate is over the **domain you sampled**, not a claim about the world —
  check the returned `bounds` match your question.
- Under noise, `recover` checks **finite-row band consistency** or abstains.
  Significant approximants and unresolved small terms can pass; significance
  against a random null does not establish exact structure. See the
  [boundary study](../../docs/CERTIFICATION_BOUNDARIES_REGISTRATION.md).
- A `fit` conjecture is **not** a certificate and must never be reported as one.
