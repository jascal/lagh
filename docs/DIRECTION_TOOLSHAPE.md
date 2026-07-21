# Direction (registered 2026-07-21): lagh as a certified tool an LLM orchestrates

**Status:** strategic architecture. Reframes the benchmark-win target and interacts with
`STRATEGY.md` (the blind benchmark) and `DIRECTION_SIGNIFICANCE.md` (α as the routing signal).
Recorded before building so the framing is not reverse-justified.

## The problem it solves

Most high-profile benchmarks are **LLM-shaped**: unstructured natural-language / multimodal
problems where interpretation *is* the task (NewtonBench is itself an LLM agent probing a black
box). `lagh` cannot read text, so *"lagh alone beats benchmark X"* is impossible for exactly the
benchmarks worth winning. The reframe: **lagh is a certified sub-solver an LLM orchestrates.** The
LLM parses the problem, sets up the queryable oracle, and calls lagh; lagh returns a certified law
or abstains; the LLM integrates the result.

## The composite guarantee (the whole point)

    composite accuracy = LLM_alone (non-lagh subset) + certified (lagh subset)

Because lagh **abstains rather than fabricates** (zero-wrong, 0/174+), it can only *add* correct
answers, never subtract:

> **composite ≥ LLM_alone, provably.** Certified tool-use that never degrades the base model and
> sometimes boosts it.

This is the exact inverse of wyly's founding disaster — a bounded expert that cut a host LLM
0.683 → 0.267 *because it never abstained* (`wyly/DIAGNOSIS.md` §3.3). `DIAGNOSIS.md` argued the
composition should work; the old count-table expert couldn't abstain, so it couldn't. lagh's
design — zero-wrong + first-class machine-readable abstention — is precisely what makes the
composition net-positive. The abstention is not a limitation here; it is the safety property that
licenses the guarantee.

## It fixes the three recorded serving blockers (`wyly/ASSETS.md`)

- **B1** (citation is a string, not a derivation): lagh's certificate + derivation tree *is* the
  citation — the recovered law, its domain, its tier/class provenance, its discriminating queries.
- **B2** (abstention string-sniffed at the hub, `a.find("isn't covered")`): lagh emits a structured
  `Abstain{domain|structural|noise|surrogate|numerical|range}` — the discrete refusal signal B2
  needed.
- **B3** (ranking cannot compare experts; flat constants 0.9 > 0.6): the certificate's **α**
  (`DIRECTION_SIGNIFICANCE.md`) is the principled cross-expert ranking. `certified-at-α=10⁻⁵⁶`
  strictly outranks an LLM's unquantified guess. The significance and tool-shape directions
  *combine*: α is the hub's routing/ranking function.

## Two concrete tool shapes (compatible, not competing)

1. **MCP server + skill.** lagh exposes `discover` / `certify` as MCP tools; a skill tells the LLM
   *when* a subproblem is law-discovery-shaped and how to set up the oracle. Portable to any
   MCP-speaking agent. This is the modern, host-agnostic form.
2. **claymore spoke (`sgiandubh`-shaped).** lagh served as a bounded OpenAI-compatible expert;
   `claymore` hub fans the query, drops abstainers (the bound *is* the router), ranks survivors by
   α, falls back to a hub LLM. The workspace-native form; the routing logic is already designed.

MCP is the transport; claymore is the routing doctrine. A first build can be the MCP server (small)
and let any agent, including a claymore hub, consume it.

## Honest risks (recorded up front)

- **Framing honesty.** The win is *LLM + lagh*, not lagh alone. Never claim lagh beats GPT-5 by
  itself; claim lagh makes an LLM agent strictly better on the certifiable subset. The composite
  guarantee is the honest, strong claim.
- **The oracle-setup surface (the real danger).** If the LLM misinterprets the problem and sets up
  the *wrong* oracle, lagh will faithfully certify a law for a question nobody asked — a
  technically-sound certificate answering the wrong thing. This is the "understand before you
  derive" trap that produced the MMLU-chemistry collapse. Mitigation: the certificate states
  exactly what oracle it answered; the LLM (and any eval) can check the oracle matches the problem
  before trusting the certified law. lagh's soundness is *relative to the oracle given* — that
  boundary must be explicit in every served answer.
- **Latency / query budget** on real black-box benchmarks with rate-limited oracles: the active
  loop's budget ledger already meters this; the acquisition policy caps it.

## Interaction with the current plan

- The NewtonBench-dev sweep still matters — it establishes **which law families lagh certifies**,
  i.e. the size of the "boost" subset the composite provides. Better lagh coverage → bigger
  provable boost. So the readiness bar (`STRATEGY.md`) is unchanged; it now also sizes the tool's
  value.
- The **blind benchmark** should be chosen with this shape in mind: an LLM-agent SR/discovery
  benchmark (LLM-SRBench, or NewtonBench's own agent protocol on the *sealed* blind analog) where
  the entry is *"LLM agent + lagh tool vs. LLM agent alone,"* and the deliverable is the composite
  guarantee demonstrated: strict non-degradation + measured boost.

## The one-sentence version

Stop trying to make a text-blind exact solver beat text-shaped benchmarks; make it the certified,
abstaining, α-ranked tool that provably upgrades an LLM agent — which is the claymore/sgiandubh
thesis the workspace already holds, finally paired with an expert disciplined enough (zero-wrong +
structured abstention) to make the composition safe.

## Generalization: this is the workspace's "LLM-as-scientist" architecture (i-orca et al.)

The tool-shape is not lagh-specific. It is the pattern the whole Orca/decompilation workspace
already embodies: **an LLM-as-scientist orchestrating a federation of certified bounded experts,
each built on the untrusted-proposer / sound-checker discipline, each tool-served, each
epistemically honest about what "verified" means.**

`i-orca` is the proof-domain exemplar and already has every piece:
- LLM **produces** proofs "at the register LLMs naturally produce" (Markdown tables) — the
  productive, untrusted proposer;
- a **cheap static `verify`** (structural linter, decidable, *strictly weaker*) — a fast pre-gate;
- the **Isabelle kernel `check`** — the sound checker; *"truth is delegated to Isabelle's kernel."*
- an **MCP server** already shipped (`pip install -e ".[mcp]"`);
- and the exact epistemic honesty lagh insists on: *green verify = skeleton well-formed ≠ proof
  true; only the kernel certifies truth* — the proof-domain twin of lagh's *certified-over-a-stated-
  domain ≠ proved-for-the-world*.

The same shape recurs across the workspace, each a sound checker for an LLM's untrusted proposals:

| tool | proposer (untrusted) | sound checker (certifies) | domain |
|---|---|---|---|
| **i-orca** | LLM proof in Markdown | Isabelle kernel | mathematical proof |
| **lagh** | curriculum + LLM law forms | exhaustive residual / exact check | law discovery |
| **rosetta** | mined circuits | Soufflé `equiv.dl` | circuit equivalence |
| **ergo** | — | verified rule core (fixpoint) | deduction |
| **fieldrun** | native decompiler | emitted certified Datalog | model → program |
| **n-orca** | Markdown NN spec | verify → runnable `nn.Module` | architecture |

**The single discipline underneath all of them** is `LEARNER.md` §2 / `equiv.dl`'s corollary /
Isabelle's kernel: *the generator can be untrusted and productive; a sound cheap checker makes it
safe. You do not need a correct proposer, only a productive one plus a sound checker.* The LLM is
that productive-but-untrusted scientist; the tools are the checkers.

**What generalizes from lagh's two directions to the whole federation:**
- The **composite guarantee** (`composite ≥ LLM_alone`) holds for *any* checker that abstains
  rather than fabricates — i-orca abstains (no kernel proof), rosetta abstains (no clean
  certificate), lagh abstains (no law certifies). Each can only add certified answers to the LLM,
  never subtract. The federation is provably non-degrading by construction.
- The **α / significance** ranking is the common currency across tools: a kernel-checked proof, an
  `equiv.dl` certificate over `|D|`, a lagh law at significance α — all are *quantified* trust
  signals `claymore` can rank against each other and against an LLM's unquantified guess. This is
  the principled multi-expert `B3` ranking, now federation-wide.

So the endpoint is not "lagh beats a benchmark" but **"the LLM-as-scientist, given a federation of
certified abstaining tools (lagh + i-orca + rosetta + …) served over MCP and routed by `claymore`
on quantified trust, provably matches-or-beats the bare LLM on any benchmark, and beats it on the
subset any tool can certify."** That is the workspace's actual thesis, and lagh + the tool-shape is
one spoke of it — the one that turns law-discovery subproblems into certified, non-degrading boosts.
