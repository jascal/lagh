"""The LawResearch machine's guarantees are STRUCTURAL -- verified here against the spec
itself, not against a model's behavior. Needs the `[machine]` extra (orca runtime)."""
import asyncio
from pathlib import Path

import pytest

orp = pytest.importorskip("orca_runtime_python")

SPEC = Path(__file__).resolve().parents[1] / "machine" / "research.orca.md"


def _mdef():
    parsed = orp.parse_orca_auto(SPEC.read_text())
    return parsed.machines[0] if hasattr(parsed, "machines") else parsed


def test_topology_is_well_formed():
    m = _mdef()
    states = {s.name: s for s in m.states}
    inits = [s.name for s in m.states if s.is_initial]
    finals = {s.name for s in m.states if s.is_final}
    assert inits == ["discovering"]                       # exactly one initial
    assert finals == {"done_certified", "done_characterized"}
    srcs = {t.source for t in m.transitions}
    # no deadlock: every non-final state has an outgoing transition
    assert all(s in srcs for s in states if s not in finals)
    # finals are terminal: no outgoing
    assert not (finals & srcs)
    # deterministic: no (source, event) pair appears twice
    pairs = [(t.source, t.event) for t in m.transitions]
    assert len(pairs) == len(set(pairs))


def test_safety_invariant_certified_only_via_certificate_events():
    # THE guarantee: a proved outcome is reachable ONLY through a lagh certificate event.
    m = _mdef()
    into_certified = {t.event for t in m.transitions if t.target == "done_certified"}
    assert into_certified == {"DISCOVERED_LAW", "RECOVERED", "VERIFIED"}
    # every abstain/refute event can reach ONLY done_characterized, never done_certified
    abstain_events = {"DISCOVER_ABSTAINED", "RECOVER_ABSTAINED", "NO_FORM", "VERIFY_REFUTED"}
    for t in m.transitions:
        if t.event in abstain_events:
            assert t.target != "done_certified"


def _drive(events):
    """Drive the machine through a sequence of (event, payload) and return the outcome."""
    m = _mdef()
    finals = {s.name for s in m.states if s.is_final}
    mach = orp.OrcaMachine(m, context={"outcome": "", "law": "", "move": ""})
    _pl = lambda p: p if isinstance(p, dict) else (getattr(p, "payload", {}) or {})  # noqa: E731
    mach.register_action("emit_law", lambda c, p: {"outcome": "proved", "law": _pl(p).get("law", "")})
    mach.register_action("record_characterization", lambda c, p: {"move": _pl(p).get("move", "")})
    mach.register_action("emit_hedge", lambda c, p: {"outcome": "empirical"})

    async def go():
        await mach.start()
        for name, payload in events:
            await mach.send(name, payload)
        return str(mach.state), mach.context

    state, ctx = asyncio.run(go())
    return state, ctx, finals


def test_certify_path_reaches_proved():
    state, ctx, finals = _drive([("DISCOVERED_LAW", {"law": "2*x_0"})])
    assert state == "done_certified" and state in finals
    assert ctx["outcome"] == "proved" and ctx["law"] == "2*x_0"


def test_report_and_stop_reaches_empirical_in_two_steps():
    state, ctx, finals = _drive([("DISCOVER_ABSTAINED", {"move": "report_and_stop"}),
                                 ("MOVE_STOP", {})])
    assert state == "done_characterized" and state in finals
    assert ctx["outcome"] == "empirical"


def test_acquire_then_recover_reaches_proved():
    # the gain path: abstain -> acquire -> recover certifies
    state, ctx, finals = _drive([("DISCOVER_ABSTAINED", {"move": "acquire_divergent"}),
                                 ("MOVE_ACQUIRE", {}), ("RECOVERED", {"law": "sqrt(x_0)"})])
    assert state == "done_certified" and ctx["outcome"] == "proved"


def test_propose_then_refute_reaches_empirical():
    state, ctx, finals = _drive([("DISCOVER_ABSTAINED", {"move": "declare_and_verify"}),
                                 ("MOVE_PROPOSE", {}), ("FORM_PROPOSED", {"form": "sin(x_0)"}),
                                 ("VERIFY_REFUTED", {})])
    assert state == "done_characterized" and ctx["outcome"] == "empirical"
