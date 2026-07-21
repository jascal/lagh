# Direction (registered 2026-07-21): box-search on abstain (Policy v3)

**Status:** capability, registered before build. Motivated by the measured finding that
**box degeneracy is the most prevalent failure shape** on NewtonBench-dev (BE overflow, snell
degrees, decay/malus/harmonic regime issues). Couples tightly with `DIRECTION_SIGNIFICANCE.md`.

## The finding

Most abstentions in the dev sweep are not "no law exists" — they are "the declared box is in a
degenerate regime": arguments driven to overflow (BE `1/(e^{ω/T}−1)`), to a constant
(BE `1/(e^{Cω/T}+1)` with `Cω/T≈0` → flat 0.5), or to a saturated edge (snell angle-in-degrees
box → all outputs ≈90°). The law is in-class; the *window* hides it. And in the blind/general case
the correct box is **unknown**, so a single fixed declaration is a structural failure point.

## The capability

On abstain, **search the box** before concluding no law exists. Extends the existing adaptive
ranging (which only *contracts* when signal falls below the floor) to also **expand** and
**relocate**:

- a bounded ladder of `K` box transforms: ×10 / ×100 expansion per axis, decade shifts of the
  center, and (existing) contraction — chosen to move the input arguments through O(1) regimes;
- retry discovery on each; stop at the first that certifies **and passes the guards below**;
- if none in `K` certifies, abstain — still first-class, now with "searched K boxes" recorded.

## The soundness guards (non-negotiable — this is why it needs registration)

Box-search is a multiple-testing procedure; unguarded it manufactures spurious certifications
(the same shape as the confident-wrong the min-domain guard just fixed). Two guards make it sound:

1. **Held-out-box certification.** A law found on box `B` must ALSO certify on a *fresh* box `B'`
   sampled from `B`'s regime but disjoint from the points used to find it — the box-level fit/cert
   split. A cherry-picked spurious fit does not survive an independent box.
2. **Multiple-testing accounting.** `K` boxes multiply the effective search: `α ≤ K·|H|·q^h`. The
   certificate records `K` (boxes tried) so the significance is honest. The min-domain guard
   (≥8 valid points) still applies per box.

With both, box-broadening is the general-case robustness the blind benchmark needs; without them
it is p-hacking. They are the *same coin* as the significance direction — box-search is precisely
why `α` must account for the search size, not just the per-fit chance.

## Interaction with the plan

- Directly lifts the current R-cap gaps (decay/malus/harmonic/BE/snell), which are box-shaped.
- Becomes standard policy (v3), so the blind benchmark is not exposed to a hand-guessed box.
- The `K` in every certificate feeds the eventual significance layer; the held-out box is a second
  independent-domain check on top of the held-out points.

## Honest limits

- It cannot rescue a law genuinely outside C1–C6 (snell's nested arc-trig may stay abstained even
  with the right box — a capability gap, not a box gap; box-search will *reveal* which is which by
  finding a good regime and *still* abstaining).
- Unbounded search would defeat the significance accounting; `K` is capped and recorded.
