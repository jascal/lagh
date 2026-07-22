---
name: law-discovery
description: Discover an exact physical law from a queryable oracle using the certified lagh tools (recover/verify/fit) and lab (problems/sample/discover). Use for any task that asks you to find y = f(x) by sampling a black box. Written to be reliable regardless of how strong the model is.
version: 1.0.0
author: lagh
---

# Discovering an exact law — the reliable recipe

You must return the EXACT law the certified tools confirm, or abstain honestly. The true
law is **often a variant** of the textbook one — NEVER assume a formula; work from samples.
Do not trust your own guesses; only a tool's certificate counts.

## Step 1 — DELEGATE first (do this by default)
Call **`lab.discover(problem)`**. lagh runs its own adaptive sampling loop — far better than
hand-picking points — and returns a certificate or a reasoned abstention.
- **Certified → you are done.** Report its `law` verbatim with `via: "recover"`.
- Only if it **abstains** do you continue below.

This step is why you do not need to be good at choosing sample points: lagh does that part.

## Step 2 — Hand-sample only if discover abstained (exact recipe)
1. Choose **at least 20 rows** of X. For EACH row, pick EACH input **independently at random**
   across its `[lo, hi]` range (do NOT vary inputs together — correlated inputs hide their
   separate exponents). If an input's `hi/lo > 100`, draw it log-uniformly.
2. `lab.sample(problem, X)` to get y.
3. `lagh.fit(X, y)` to scout. If its diagnosis is `acquire_more_data` or `continuum`,
   **widen the box 2–10× beyond `[lo, hi]`** and repeat step 1.
4. `lagh.recover(X, y)`.
   - Certified → report it.
   - Abstains → **call `lab.discover(problem)` again** (more budget may pin it) before giving up.

## Step 3 — A declared hypothesis (last resort only)
If discover AND recover both abstain and you have a real reason to believe one specific form,
`lagh.verify(X, y, "<form in x_0..>")`. A `consistent` result means the form fits but a
constant is not identifiable from the data — keep that hedge; it is weaker than `pinned`.

## Honesty — non-negotiable
- Report ONLY what a tool **certified**. If nothing certifies, **abstain** — never invent,
  never fall back to the textbook formula.
- A `lagh.fit` output is a **guess**, never a result.
- Attribute every claim: `[recover]`/`[discover]` = certified by lagh; `[verify]` = your
  hypothesis checked; `[me]` = your own inference. Do not present a guess as a certificate.
- `pinned` = the exact law, uniquely identified. `consistent` = fits, but a constant is not
  pinned. Report whichever the tool returned; do not upgrade it.

## Final output — exactly one line, raw text
No markdown, no bold/italics. Copy the law VERBATIM from the tool with every `*` and `**`
operator intact (e.g. `3*x_0**2`, never `3x_0^2` or `3x_0*2`):

RESULT: {"status": "certified"|"consistent"|"abstained", "law": "<sympy expr in x_0.. or null>", "via": "discover"|"recover"|"verify"|"none"}
