---
name: law-discovery
description: Discover an exact physical law from a queryable oracle using the certified lagh tools (recover/verify/fit) and lab (problems/sample/discover). On abstain, lagh returns a structural CHARACTERIZATION + a research move; follow it. Use for any task that asks you to find y = f(x) by sampling a black box. Written to be reliable regardless of how strong the model is.
version: 1.1.0
author: lagh
---

# Discovering an exact law — the reliable recipe

You must return the EXACT law the certified tools confirm, or abstain honestly (with a
hedged characterization if lagh gives you one). The true law is **often a variant** of the
textbook one — NEVER assume a formula; work from samples. Only a tool's certificate counts.

**The frame:** you are a *bold proposer*, lagh is a *sound checker*. Every form you invent
is run through `verify`/`recover`; a wrong guess is rejected, never certified. So explore
**boldly** — the honesty constraint is at the OUTPUT (report only what certified, plus any
labeled hedge), not at the exploration.

## Step 1 — DELEGATE first (MANDATORY — your very first tool call)
Your FIRST action for every problem MUST be `lab.discover(problem)`. lagh runs its own
adaptive sampling loop — far better than hand-picking points — and returns a certificate or
a reasoned abstention.
- **Certified → you are done.** Report its `law` verbatim with `via: "discover"` and STOP.
- Only if it **abstains** do you continue.

**HARD RULE — you may NOT abstain without having called `lab.discover` first.** Reporting
`{"status":"abstained","via":"none"}` is ONLY legal after `lab.discover` has actually run.
"I'm not sure" is never a reason to abstain — it is the reason to call `lab.discover`.
Skipping the tool and abstaining is a FAILURE, not honesty.

## Step 2 — Read the diagnosis and follow the research move
On abstain, lagh returns a **`characterization`** (a GUESS about the function's structure —
NOT a certificate) with a **`research.move`**. Do exactly what the move says, then re-check:

| `research.move` | what you do |
|-----------------|-------------|
| **`report_and_stop`** | An exact law is impossible/unreachable here (e.g. an irrational exponent, or no structure). **Report the characterization as a hedge and abstain. Do NOT sample again.** |
| **`acquire_divergent`** | Rivals fit the sampled range. Sample a **WIDER / asymptotic** regime (widen the box 10×, favor extremes) where the candidate forms separate, then `lagh.recover`. |
| **`acquire_more_data`** | Sample **≥20 more rows**, each input **independently at random** across a **2–10× wider** box (log-uniform if `hi/lo>100`), then `lagh.recover`. |
| **`declare_and_verify`** | The shape looks non-algebraic (trig / inverse-trig / saturating). **Declare a specific form** from prior knowledge and `lagh.verify(X, y, "<form>")`. A `consistent` result = fits but a constant isn't identifiable — keep that hedge. |

**STOP RULE — at most 2 research moves.** If two moves do not produce a certificate, switch
to `report_and_stop`: report the characterization and abstain. Do NOT loop
discover→recover→discover; it rarely changes the verdict and wastes time.

## Honesty — non-negotiable
- Report ONLY what a tool **certified** (`recover`/`verify` → `pinned`/`consistent`). If
  nothing certifies after the steps above, abstain — but include the characterization as a
  **hedge** ("could not certify; the samples look like <class>"), never as a law.
- A `lagh.fit` output or a `characterization` is a **guess**, never a result.
- `pinned` = the exact law, uniquely identified. `consistent` = fits, but a constant is not
  pinned. Report whichever the tool returned; never upgrade it.
- Attribute: `[discover]`/`[recover]` = certified by lagh; `[verify]` = your hypothesis
  checked; `[characterization]` = lagh's uncertified structural guess; `[me]` = your own.

## Final output — exactly one line, raw text
No markdown, no bold/italics. Copy the law VERBATIM from the tool with every `*` and `**`
operator intact (e.g. `3*x_0**2`, never `3x_0^2` or `3x_0*2`). On abstain, put lagh's
characterization `class` (or null) in `characterization`:

RESULT: {"status": "certified"|"consistent"|"abstained", "law": "<sympy expr in x_0.. or null>", "via": "discover"|"recover"|"verify"|"none", "characterization": "<class or null>"}
