# machine LawResearch

> The certified law-discovery research loop as a verifiable state machine
> (docs/DIRECTION_CHARACTERIZATION.md, docs/DIRECTION_MACHINE.md). The topology IS the
> guarantee: `done_certified` is reachable ONLY via a lagh certificate event
> (DISCOVERED_LAW / RECOVERED / VERIFIED), so "a certified output means lagh certified it"
> is a graph invariant, not a prompt. The machine drives; the LLM is a bounded proposer
> called at exactly one state (proposing). At most ONE research move by construction (the
> DAG has no path back to routing), so the loop cannot thrash.

## context

| Field | Type | Default |
|-------|------|---------|
| problem | string | "" |
| certified | bool | false |
| law | string | "" |
| strength | string | "" |
| characterization | string | "" |
| move | string | "" |
| outcome | string | "" |

## events

- DISCOVERED_LAW
- DISCOVER_ABSTAINED
- MOVE_STOP
- MOVE_ACQUIRE
- MOVE_PROPOSE
- RECOVERED
- RECOVER_ABSTAINED
- FORM_PROPOSED
- NO_FORM
- VERIFIED
- VERIFY_REFUTED

## state discovering [initial] "lab.discover: lagh runs its own adaptive box-search loop"
> on DISCOVERED_LAW -> done_certified
> on DISCOVER_ABSTAINED -> routing

## state routing "read the characterization's research move (set by lagh, not the model)"
> on MOVE_STOP -> done_characterized
> on MOVE_ACQUIRE -> acquiring
> on MOVE_PROPOSE -> proposing

## state acquiring "hand-sample a wider/divergent regime, then lagh.recover (data path, fast)"
> on RECOVERED -> done_certified
> on RECOVER_ABSTAINED -> done_characterized

## state proposing "the ONE LLM call: propose a declared form to check"
> on FORM_PROPOSED -> verifying
> on NO_FORM -> done_characterized

## state verifying "lagh.verify the declared form (sound checker)"
> on VERIFIED -> done_certified
> on VERIFY_REFUTED -> done_characterized

## state done_certified [final] "outcome=proved: emit the certified law + strength"

## state done_characterized [final] "outcome=empirical: emit the uncertified structural hedge (an 'unresolved' class is the degenerate 'open' rung)"

## transitions

| Source | Event | Guard | Target | Action |
|--------|-------|-------|--------|--------|
| discovering | DISCOVERED_LAW | | done_certified | emit_law |
| discovering | DISCOVER_ABSTAINED | | routing | record_characterization |
| routing | MOVE_STOP | | done_characterized | emit_hedge |
| routing | MOVE_ACQUIRE | | acquiring | |
| routing | MOVE_PROPOSE | | proposing | |
| acquiring | RECOVERED | | done_certified | emit_law |
| acquiring | RECOVER_ABSTAINED | | done_characterized | emit_hedge |
| proposing | FORM_PROPOSED | | verifying | |
| proposing | NO_FORM | | done_characterized | emit_hedge |
| verifying | VERIFIED | | done_certified | emit_law |
| verifying | VERIFY_REFUTED | | done_characterized | emit_hedge |

## actions

| Name | Signature |
|------|-----------|
| emit_law | `(ctx, event) -> Context` |
| record_characterization | `(ctx, event) -> Context` |
| emit_hedge | `(ctx, event) -> Context` |

## properties

- reachable: done_certified
- reachable: done_characterized

<!-- SAFETY INVARIANT (verify with `orca verify`): done_certified has NO inbound
transition except on DISCOVERED_LAW / RECOVERED / VERIFIED -- all three are lagh
certificate events. So "a proved outcome implies lagh certified it" is structural. The
abstain events (DISCOVER_ABSTAINED / RECOVER_ABSTAINED / NO_FORM / VERIFY_REFUTED) can
reach ONLY done_characterized. This is the zero-wrong guarantee as graph topology. -->

