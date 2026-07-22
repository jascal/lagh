"""The lab is source-pluggable, never leaks the truth, and the server respects LAB_SOURCE."""

import numpy as np
import pytest

from lagh.lab import available, get_source


def test_proxy_source_cards_never_leak_the_truth():
    src = get_source("proxy")
    cards = src.problems()
    assert {c["id"] for c in cards} >= {"orbit", "steps", "exotic"}
    blob = str(cards)
    for bad in ("fn", "truth", "lambda", "6.674e-11", "cos(theta)"):
        assert bad not in blob, f"card leaked {bad!r}"


def test_proxy_sample_returns_true_oracle_values():
    src = get_source("proxy")
    y = src.sample("orbit", [[2.0, 3.0, 1.5], [4.0, 5.0, 2.0]])
    exp = [6.674e-11 * 2 * 3 / 1.5 ** 2, 6.674e-11 * 4 * 5 / 2.0 ** 2]
    assert np.allclose(y, exp)


def test_registry_and_unknown_source():
    assert set(available()) >= {"proxy", "newtonbench"}
    with pytest.raises(ValueError, match="unknown LAB_SOURCE"):
        get_source("nope")


def test_newtonbench_source_same_contract():
    # the REAL dev benchmark behind the same interface (needs the newtonbench adapter)
    nb = pytest.importorskip("lagh.adapters.newtonbench")
    src = get_source("newtonbench")
    cards = src.problems()
    assert len(cards) == 108 and all("/" in c["id"] for c in cards)
    c0 = next(c for c in cards if c["id"].startswith("m0_gravity/easy"))
    y = src.sample(c0["id"], [c0["suggested_domain"][0]])  # sample the lo corner
    assert len(y) == 1 and np.isfinite(y[0])


def test_server_respects_lab_source_and_hides_truth(monkeypatch):
    import sys
    import types
    reg = {}

    class Stub:
        def __init__(self, *a, **k):
            self.instructions = k.get("instructions", "")
        def tool(self):
            def d(f):
                reg[f.__name__] = f
                return f
            return d
        def run(self):
            pass

    fake = types.ModuleType("mcp.server.fastmcp"); fake.FastMCP = Stub
    monkeypatch.setitem(sys.modules, "mcp", types.ModuleType("mcp"))
    monkeypatch.setitem(sys.modules, "mcp.server", types.ModuleType("mcp.server"))
    monkeypatch.setitem(sys.modules, "mcp.server.fastmcp", fake)
    monkeypatch.setenv("LAB_SOURCE", "proxy")

    from lagh.lab.server import build_server
    build_server()
    # safe surface: query tools + the delegate-to-lagh `discover`; never a ground-truth tool
    assert set(reg) == {"problems", "sample", "discover"} and "truth" not in reg
    out = reg["sample"]("orbit", [[2.0, 3.0, 1.5]])
    assert out["source"] == "proxy" and "y" in out
