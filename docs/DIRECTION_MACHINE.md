# Direction: the research loop as a verified state machine (not a prompt)

Status: **proposal + working prototype** (2026-07-22). Supersedes the prompt-driven
orchestration in `deploy/hermes/` for the research loop. Companion to
`DIRECTION_CHARACTERIZATION.md` (the ladder + routing this machine executes).

## The category error we were making

The Hermes composite drove the loop with an LLM and *prompted* it into the right control
flow: "call discover FIRST," "don't abstain without discover," "at most one research
move," "use `lagh.recover` not `lab.discover`," "route on `research.move`." Every one of
those is a **transition rule** written as English and handed to a non-deterministic agent.
The failures were all the agent not honoring the prose: M3 wandered, GPT under-delegated,
Grok re-called discover → the 150s timeouts. We kept patching the prompt because a prompt
isn't an executable.

The tell: **we had already designed the state machine.** `characterize() → research_move`
*is* its transition function. We built the transition table, then spent days teaching it
to a model via prose. Don't teach it — **execute it.**

## The architecture

`machine/research.orca.md` — the control flow as a verified Orca state machine.
`machine/driver.py` — a deterministic driver that binds each state to its lagh call.

```
discovering ──DISCOVERED_LAW─────▶ done_certified   (proved)
     │
     └─DISCOVER_ABSTAINED─▶ routing
                              ├─MOVE_STOP────────▶ done_characterized  (empirical)
                              ├─MOVE_ACQUIRE─▶ acquiring ──RECOVERED──▶ done_certified
                              │                      └─RECOVER_ABSTAINED▶ done_characterized
                              └─MOVE_PROPOSE▶ proposing ─FORM_PROPOSED▶ verifying ─VERIFIED─▶ done_certified
                                                 │                          └─VERIFY_REFUTED▶ done_characterized
                                                 └─NO_FORM──────────────────────────────────▶ done_characterized
```

**The machine drives; the LLM is a bounded proposer** called at exactly one state
(`proposing`). Everything mechanical — the delegation discipline, the stop rule, the
no-loop guarantee, the routing — is deterministic control flow the model cannot violate,
because the model isn't holding the wheel. This is "**bold proposer, sound checker**" made
architectural instead of aspirational.

## What is now GUARANTEED (topology, not prose)

`orca verify` (and `tests/test_machine.py`) check these against the spec:

- **Zero-wrong is structural.** `done_certified` has NO inbound transition except
  `DISCOVERED_LAW` / `RECOVERED` / `VERIFIED` — all three are lagh *certificate* events. So
  "a proved outcome ⟹ lagh certified it" is a graph invariant. The abstain/refute events
  reach ONLY `done_characterized`. No prompt can break this.
- **It always terminates in the degradation ladder.** The graph is a DAG; every path ends
  in `done_certified` (proved) or `done_characterized` (empirical hedge).
- **At most ONE research move, by construction.** `acquiring`/`proposing` lead only to
  terminals — there is no path back to `routing`. The thrash-to-timeout is structurally
  impossible; no "stop rule in a prompt" needed.
- **`lab.discover` runs at most once.** `discovering` is the initial state and never
  re-entered.

## The driver, and graceful degradation

`machine/driver.py` binds states to computation: `discovering`→`recover(box_search)`,
`acquiring`→ hand-sample + `recover` (data path, fast — never a second box-search),
`proposing`→ the injected `propose_fn` (the ONE LLM call), `verifying`→ `verify`. The LLM
is `propose_fn(ctx, box) -> form|None`, injected. **With `propose_fn=None` the whole loop
runs on lagh alone** and every abstain routes to `report_and_stop` — so the composite
*provably* never does worse than the tool, and the LLM is a pure optional enhancement at
one point. True model-agnosticism: the control flow is byte-identical across models; only
proposal quality varies.

## Validation (2026-07-22)

End-to-end, **no LLM** (`propose_fn=None`), real NewtonBench cells:

- `m0_gravity/easy/v0` → **4.7s → done_certified / proved**, `6.674e-5*x_0*x_1/x_2**(3/2)`.
- `m4_snell_law/easy/v0` → **86.9s → done_characterized / empirical**, class `unresolved` —
  the `report_and_stop` path, **bounded, no thrash**. The entire timeout saga is gone: the
  machine *cannot* loop.

Topology + drive tests: `tests/test_machine.py` (6, all pass). Structural parse verified
against the runtime; `orca verify` (MCP/CLI) is the formal gate for the `## properties`.

## Relationship to the Hermes prototype

The Hermes agent approach was **informative thinking-out-loud in code**: it surfaced the
characterization ladder, the routing table, the report_and_stop terminal, and every failure
mode. That logic is now debugged and packaged as a *verifiable executable* instead of a
fragile prompt. Hermes is no longer the orchestrator; if used at all, it is just an LLM API
client for `propose_fn`.

## Next

1. Wire `propose_fn` to an LLM API call (Grok/GPT/MiniMax) — a *narrow* prompt: "given
   these samples and this characterization, propose ONE sympy form to verify, or say none."
2. Run the **machine composite** on NewtonBench-dev: expect the full-108 result (67/108,
   0 CW) with the abstains now *fast and clean* (no timeouts), and any `declare_and_verify`
   gains the proposer earns — all under a structurally-guaranteed zero-wrong.
3. Package `machine/` as a first-class lagh entry point (a `lagh-research` script / MCP
   verb) so the verified loop is the product surface, not the bench harness.
