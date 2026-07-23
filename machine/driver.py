"""Deterministic driver for the LawResearch Orca machine (machine/research.orca.md).

The verified state machine owns the CONTROL FLOW; this driver owns the COMPUTATION. For
each state it runs the corresponding lagh call, derives the semantic event, and lets the
machine route. The topology guarantees (terminates, ends in the degradation ladder, at
most one research move, `done_certified` reachable ONLY via a certificate event) come from
the spec + `orca verify`, not from this code.

The LLM is a BOUNDED PROPOSER, injected as `propose_fn` and called at exactly ONE state
(`proposing`). With `propose_fn=None` the whole loop runs on lagh alone and every abstain
routes to `report_and_stop` -- so the composite provably never does worse than the tool.

Requires the orca runtime (`pip install orca-runtime-python`) so the loop is driven by the
verified spec, not a hand-rolled mirror of it.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from lagh.mcp.core import recover, verify

_SPEC = Path(__file__).with_name("research.orca.md")


def _pl(p) -> dict:
    return p if isinstance(p, dict) else (getattr(p, "payload", {}) or {})


def _load_machine():
    from orca_runtime_python import OrcaMachine, parse_orca_auto
    parsed = parse_orca_auto(_SPEC.read_text())
    mdef = parsed.machines[0] if hasattr(parsed, "machines") else parsed
    finals = {s.name for s in mdef.states if s.is_final}
    m = OrcaMachine(mdef, context={"problem": "", "certified": False, "law": "",
                                   "strength": "", "characterization": "", "move": "",
                                   "form": "", "outcome": ""})
    # bookkeeping actions only record the driver's result; the FLOW is the verified topology
    m.register_action("emit_law", lambda ctx, p: {
        "certified": True, "law": _pl(p).get("law", ""),
        "strength": _pl(p).get("strength", ""), "outcome": "proved"})
    m.register_action("record_characterization", lambda ctx, p: {
        "characterization": _pl(p).get("class", ""), "move": _pl(p).get("move", "")})
    m.register_action("emit_hedge", lambda ctx, p: {"certified": False, "outcome": "empirical"})
    return m, finals


def _sample(oracle, box, n=200, widen=1.0, seed=0):
    lo = np.asarray(box[0], float) / widen
    hi = np.asarray(box[1], float) * widen
    rng = np.random.default_rng(seed)
    X = np.exp(rng.uniform(np.log(np.maximum(lo, 1e-12)), np.log(np.maximum(hi, 1e-12)),
                           (n, len(lo))))
    return X, np.asarray(oracle(X), float)


def _step(state, ctx, oracle, box, propose_fn, seed):
    """Map the current state to its lagh/LLM computation; return (event_name, payload)."""
    if state == "discovering":
        r = recover(oracle=oracle, box=box, box_search=True, seed=seed)
        if r.get("certified"):
            return "DISCOVERED_LAW", {"law": r["law"], "strength": r.get("strength", "")}
        ch = r.get("characterization", {}) or {}
        return "DISCOVER_ABSTAINED", {"class": ch.get("class", "unresolved"),
                                      "move": ch.get("research", {}).get("move", "report_and_stop")}
    if state == "routing":
        move = ctx.get("move", "report_and_stop")
        if move in ("acquire_divergent", "acquire_more_data"):
            return "MOVE_ACQUIRE", {}
        if move == "declare_and_verify":
            return "MOVE_PROPOSE", {}
        return "MOVE_STOP", {}                                # report_and_stop (the common case)
    if state == "acquiring":
        widen = 10.0 if ctx.get("move") == "acquire_divergent" else 4.0
        X, y = _sample(oracle, box, widen=widen, seed=seed + 1)
        r = recover(X.tolist(), y.tolist())                  # DATA path -- fast, no 2nd box-search
        if r.get("certified"):
            return "RECOVERED", {"law": r["law"], "strength": r.get("strength", "")}
        return "RECOVER_ABSTAINED", {}
    if state == "proposing":
        form = propose_fn(ctx, box) if propose_fn else None   # the ONE bounded LLM call
        return ("FORM_PROPOSED", {"form": form}) if form else ("NO_FORM", {})
    if state == "verifying":
        X, y = _sample(oracle, box, seed=seed + 2)
        r = verify(X.tolist(), y.tolist(), ctx.get("form", ""))
        if r.get("certified"):
            return "VERIFIED", {"law": r["law"], "strength": r.get("strength", "")}
        return "VERIFY_REFUTED", {}
    raise RuntimeError(f"no step for state {state!r}")


async def run(problem: str, oracle, box, *, propose_fn=None, seed: int = 0) -> dict:
    """Drive the verified LawResearch machine to a terminal. Returns
    {outcome: 'proved'|'empirical', law, strength, characterization, final_state}.
    `propose_fn(ctx, box) -> form_str|None` is the injected LLM proposer (None -> lagh alone)."""
    m, finals = _load_machine()
    m.context["problem"] = problem
    await m.start()
    guard = 0
    while str(m.state) not in finals and guard < 8:          # 8 = spec's longest path; a backstop
        guard += 1
        name, payload = _step(str(m.state), m.context, oracle, box, propose_fn, seed)
        if name == "FORM_PROPOSED":                          # carry the form for verify (context)
            m.context["form"] = payload.get("form", "")
        await m.send(name, payload)
    c = m.context
    return {"outcome": c.get("outcome", "empirical"), "law": c.get("law", ""),
            "strength": c.get("strength", ""), "characterization": c.get("characterization", ""),
            "final_state": str(m.state)}


def run_sync(problem: str, oracle, box, *, propose_fn=None, seed: int = 0) -> dict:
    import asyncio
    return asyncio.run(run(problem, oracle, box, propose_fn=propose_fn, seed=seed))
