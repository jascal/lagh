# Stretch goal: make a refusal choose the next measurement

**Goal in one sentence.** Turn a structural abstain into the specific measurement
that would resolve it, and beat the fixed acquisition strategies already in the
code on how many measurements that takes — at unchanged soundness.

## Why this one

Of everything on lagh's strategic list, this is the capability that is hard to
copy. Checking residuals and auditing pipelines are valuable and commoditisable.
Saying *"your current data cannot separate these two explanations, and here is the
measurement that would"* is not, and C3 is the existence proof: a refusal that
decoded into a drag measurement agreeing with an independent route to 2%.

lagh already has every part except the decision.

| piece | where | what it does |
|---|---|---|
| the rival classes | `certify.coherent` (certify.py:1151) | returns the materially-different laws that all certify |
| where rivals differ | `certify.max_divergence` (certify.py:345) | max disagreement of two laws over a probe box |
| disagreement-scored queries | `acquisition.run_active` (acquisition.py:114) | already scores candidate points by `w_disagreement` |

What is missing is that **nothing uses the refusal's own content to choose where
to measure next.**

## The baselines already exist, and they are fixed

- `mcp.core.recover` answers a structural abstain with
  `suggested_box = [lo/10, hi*10]` (core.py:226).
- `acquisition._box_ladder` (acquisition.py:274-278) tries five fixed rungs:
  as-declared, expand-x10, shift-up-decade, shift-down-decade, expand-x100.

Neither looks at what actually certified. Those are the numbers to beat.

## The apex target: `rational-d1`

It is the **single uncertified cell of 36** in the reach audit
(`experiments/results/test_selection_reach.json`), and has been since July.
`experiments/reach/audit.py:82` defines it as `(2x + 1)/(x + 3)`.

`MUNTZ_ARBITRATION.md` already diagnosed it: the dense channel proposes two
fractional-power twins, dof 34 and 45 of n = 80, **neither of them the truth**,
and it registered that such twins *"diverge only outside the sampled box"* and
that the abstain is permanent **under that box's evidence**. Acquisition outside
the box is therefore exactly the honest route that registration leaves open.

Resolving it by measurement — not by widening the hypothesis class — would
convert a documented permanent abstain into a certificate. That is the apex; it
is not required for the run to be a success.

## What to build

From a structural abstain, take the surviving classes, choose the query that
maximises their disagreement per unit cost, acquire there, and re-certify.
Manufacture a twin set on a substrate you own — the reach ladder and
LawSystemBench give known truth and an oracle you can query — so that provenance
is not the bottleneck it was in the tweezers campaign.

## The constraint that makes or breaks it

**Choosing the next box from the data you will certify on is the C2 drive-scale
defect wearing new clothes**: the estimator returns its own assumption, and the
error looks like agreement. Declare the split before you look. If the box cannot
be chosen without touching the certification sample, stop and say so — that is a
result, not a failure to be engineered around.

Second trap: acquiring outside the box **changes the domain the certificate
claims**. A resolved twin must state its new domain, never inherit the old one.

## Register before the scored run

Pick numbers from a pilot: how many of N twin pairs resolve, and in how many
queries against each fixed baseline. **Predict a fraction that stays
unresolved** — the Müntz registration says some twins are genuinely
indistinguishable, so a method that resolves all of them is suspicious rather
than impressive.

## Both outcomes are results

- **Fewer measurements to resolution at unchanged soundness** is the first
  demonstrated reach gain in this project.
- **No gain, with a measured account of why**, retires an idea that has been
  assumed to work since the acquisition module was written.

Zero confident-wrong across every arm is the invariant, and it is checkable here
because the twins have known truth.

## Preconditions, both small

1. **Restate the headline claim with its real scope.** The falsifiability pass
   showed `check()` itself accepted four kinds of invalid evidence — empty
   domain, NaN observation, NaN input row under a constant law, NaN band — and
   the zero-wrong record was protected by callers pre-filtering, not by the
   checker refusing. "Never a confident wrong answer" needs its domain attached.
   One honest paragraph in `README.md` and `docs/INSTRUMENT_REPORT.md`, backed by
   the 44 committed witnesses. Not a documentation programme.
2. **Every new gate ships with its falsifying witness**, per
   `docs/FALSIFIABILITY_AUDIT.md`. And remember its own lesson: a witness that is
   too specific makes a non-fix look like a fix.

## House rules (unchanged)

Register predictions with kill criteria before running. Tag every claim
`proved` / `empirical` / `open`. Commit artifacts and keep them reproducible.
Run tests one file per process. No paid proposer. No publication push.

Before crediting any test, exhibit the input that makes it fail.
