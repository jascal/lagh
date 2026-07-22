# lagh MCP server — `recover` / `verify` / `fit`

The tool-shape from [`docs/DIRECTION_TOOLSHAPE.md`](../../docs/DIRECTION_TOOLSHAPE.md),
made runnable. An LLM parses a problem, samples the black-box oracle, and calls these
three tools; lagh returns a **certified law or a reasoned abstention** (never a
confident-wrong answer). The composite is provably non-degrading: lagh only ever adds
certified answers or abstains.

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
| **`verify(X, y, form, sigma)`** | bounded | check a *declared* form | `Certificate` \| `Abstention` |
| **`fit(X, y, sigma)`** | **UNBOUNDED** | best-guess + identifiability diagnosis | `Conjecture[]` + `Diagnosis` — **no `certified` field** |

- **Certificate strength** (`recover`/`verify`): `pinned` = this exact law, no rival
  within the noise; `consistent` = a *declared irrational* (`x_0**E`) fits but the
  constant is not identifiable from finite data. Only `pinned` is the strong claim.
- **The wall is structural.** `fit`'s result has no `certified` key — a guarantee
  cannot be read off it by accident. `fit` is a scout; to make a conjecture
  trustworthy, feed it to `recover` or `verify`.

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
- Under noise, `recover` certifies the exact **structure** or abstains; coefficients
  are noise-limited (`RNOISE_STUDY.md`: structure-or-abstain holds to ~1% relative).
- A `fit` conjecture is **not** a certificate and must never be reported as one.
