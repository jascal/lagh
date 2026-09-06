# Refusal acquisition pilot findings

Status: **empirical** for recorded runs; **open** for an efficiency gain.

P1's sixteen outcomes are retained under experiments/results/acquisition_p1.
Code and parameters are hashed in every artifact; every query plan was flushed
to its plan ledger before querying the oracle. New gate witnesses are in
tests/test_refusal_acquisition.py (11 passing cases, mocked discovery).

| Initial case | Initial structural? | Guided | Fixed x10 | Matched ladder | Native ladder |
|---|---|---:|---:|---:|---:|
| rational-a, seed10 | no: certified approximant |480 wrong|480 wrong|480 wrong|480 wrong|
| rational-b, seed10 | yes |560 exact|560 exact|640 exact|598 exact|
| rational-c, seed10 | yes |560 exact|560 exact|640 exact|598 exact|
| restricted-a, seed10 | no: same approximant |480 wrong|480 wrong|480 wrong|480 wrong|

Counts include400 initial observations and80 fresh final observations. Native
counts include its own60-row holdout as well. The nonstructural rows did not run
acquisition at all: P1's initial-certificate branch accepted the existing
candidate after an in-box final check. They are outside the intended structural
refusal workflow but remain soundness failures of the implemented P1 wrapper.
P1 therefore cannot establish unchanged soundness. Its restricted-control
prediction failed: no twin refusal existed on that draw. No scored claim follows.

**Empirical:** the seed10 approximant is structurally wrong and diverges from
known rational truth by about7e126 on the independent [.005,300] probe, despite
passing the original-domain residual checks. This is a concrete counterexample
to turning local residual coverage into an exact-form claim. No grammar or
tolerance is amended to hide it.

**Empirical:** guided chose x=300 in both structural cases. It saved one batch
against the matched ladder, but tied x10. Native found the exact rational law
using fewer initial fitting rows on the original domain; it did not need a box
expansion. Its extra holdouts account for the raw38-query difference. These
results do not establish that refusal content is superior to simple expansion.

## Split classes versus surviving explanations

P1 recorded full-design checks of its retained representatives. Rational-b's
two representatives have0 and2 misses on all400 initial rows; rational-c's four
have1,0,2,2 misses. Only one in each set remains fully compatible with the
available observations. Neither law is the rational truth.

Reconstructing P0's historical apex draw and checking both stored expressions
on all400 initial rows gives2 and3 misses respectively. Both fail before a new
measurement. They certify on their internal certification split, not on every
observation already available. The sticky structural refusal is real, but
calling these two globally surviving explanations would be incorrect.

**Open implication:** disagreement among these expressions may still be a useful
search heuristic; a measurement chosen from it cannot be credited with first
falsifying a law that old observations already falsify. The strict surviving-
rival workflow needs at least two representatives that pass all design rows.
The apex remains a diagnostic of reach ordering; it is not quietly dropped.

## P2 survivor screen

**Empirical:** only3/16 seeds qualify, missing the registered at-least4 prediction.
All16 records are retained under acquisition_p2 (6400 design observations).
Seed30 has3 fully compatible representatives; seeds32 and35 have2 each. Eight
other draws already certify structurally wrong approximants before acquisition.
The remaining five are structural but have fewer than2 full-data survivors.

These three qualified, frozen design states are a manufactured twin bank for
follow-up experiments, not a random sample of successful discovery tasks. Their
follow-up and final observations have not been sampled by P2. The ordinary
engine's initial false-exact behavior remains a finding, not a repaired guarantee.

## P3 corrected pilot

**Empirical:** guided, fixed x10 and the matched ladder each resolve3/6 problem
conditions, with all three restricted conditions unresolved. Successful totals
are560,560,640 respectively. No corrected terminal certificate is wrong across
the24 arms. The guided restricted seed30 case spends80 additional measurements
before stopping; seeds32 and35 have no above-band informative probe and stop at400.

The native ladder returns the exact truth in its raw result, but those original-
box measurements leave older rival explanations compatible. All six native
results therefore fail the COMMON twin-resolution criterion and are reported
unresolved at518 total queries, with the raw law/verdict retained. This is not
a claim that native cannot recover the correct formula; it does. The comparison
distinguishes correct output from measured elimination of known competitors.

**Empirical apex diagnostic:** the pinned P1 heuristic recovers exact
(2*x+1)/(x+3) at560 queries, tying fixed x10. Matched ladder takes640; native
takes598 including its own holdout and the common final guard. All four laws
are exactly equivalent to truth. Guided's new finite domain extends to x=300;
this does not amend the35/36 passive reach artifact. The old rivals already
failed original observations, so this measures escape from reach ordering,
not the first experimental falsification of those rivals.

P3 gate counter-inputs are replayed against pinned P1 in
acquisition_gate_witnesses.json: four new gate cases fail before repair; all15
cases pass afterward. The old unit fixture(x,2x) was itself not a compatible
pair on exact x data; positive tests now use(x,x+1e-16*x**8), compatible in the
initial box but separated at300. The former fixture remains as a negative case.
