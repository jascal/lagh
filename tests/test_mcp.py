"""In-process tests for the three MCP tool acts (no SDK / no wire).

Validates the contract in docs/DIRECTION_TOOLSHAPE.md: recover/verify are bounded
(cert-or-abstain, with a pinned/consistent strength) and fit is unbounded (conjectures
+ diagnosis, and structurally CANNOT carry a certificate -- no `certified` key)."""

import numpy as np

from lagh.mcp import fit, recover, verify


def _data(f, dim, lo, hi, n=120, sigma=0.0, seed=0):
    rng = np.random.default_rng(seed)
    X = np.exp(rng.uniform(np.log(lo), np.log(hi), (n, dim)))
    y = f(X)
    if sigma:
        y = y * (1 + sigma * rng.standard_normal(len(y)))
    return X, y


def test_recover_certifies_a_clean_rational_law():
    X, y = _data(lambda X: 3 * X[:, 0] ** 2, 1, 0.5, 4.0)
    r = recover(X, y)
    assert r["certified"] is True and r["tag"] == "proved"
    assert r["strength"] == "pinned"


def test_recover_abstains_returns_structured_reason_not_a_guess():
    # a jagged non-smooth target: recover must abstain with a reason, never fabricate
    X, y = _data(lambda X: np.round(X[:, 0] * 7) % 3, 1, 0.5, 4.0)
    r = recover(X, y)
    assert r["certified"] is False and r["tag"] == "open"
    assert "abstain" in r


def test_verify_confirms_a_correct_declared_form_pinned():
    X, y = _data(lambda X: 2.5 * X[:, 0] * X[:, 1], 2, 0.5, 3.0)
    r = verify(X, y, "x_0*x_1")            # scale refit to 5/2
    assert r["certified"] is True
    assert r["strength"] == "pinned"


def test_verify_refutes_a_wrong_declared_form():
    X, y = _data(lambda X: 2.5 * X[:, 0] * X[:, 1], 2, 0.5, 3.0)
    r = verify(X, y, "x_0**2*x_1")         # wrong structure
    assert r["certified"] is False and r["tag"] == "open"


def test_verify_declared_irrational_is_consistent_never_pinned():
    e = np.e
    X, y = _data(lambda X: 2.0 * X[:, 0] ** e, 1, 1.0, 3.0)
    r = verify(X, y, "x_0**E")             # declared irrational exponent
    assert r["certified"] is True
    assert r["strength"] == "consistent"   # fits, but e is not identifiable


def test_fit_has_no_certified_field_and_carries_a_diagnosis():
    X, y = _data(lambda X: 3 * X[:, 0] ** 2, 1, 0.5, 4.0)
    r = fit(X, y)
    assert "certified" not in r            # the structural wall
    assert r["tag"] == "exploratory"
    assert "diagnosis" in r and "next_action" in r
    assert isinstance(r["conjectures"], list)


def test_fit_flags_a_continuum_exponent_as_not_pinning():
    # x^e -- a free exponent that does not pin to a small rational
    e = np.e
    X, y = _data(lambda X: 2.0 * X[:, 0] ** e, 1, 1.0, 4.0)
    r = fit(X, y)
    assert r["diagnosis"]["kind"] == "continuum"
    assert r["next_action"] in ("declare_and_verify", "acquire_more_data")


def test_server_registers_the_three_acts_and_delegates_to_core(monkeypatch):
    """Prove the transport wiring without the real SDK: stub FastMCP, build the
    server, and confirm exactly recover/verify/fit are registered and each delegates
    to core (so `fit` still returns no `certified` field through the transport)."""
    import sys, types

    registered = {}

    class _StubServer:
        def __init__(self, *a, **k):
            pass

        def tool(self):
            def deco(fn):
                registered[fn.__name__] = fn
                return fn
            return deco

        def run(self):
            pass

    fake = types.ModuleType("mcp.server.fastmcp")
    fake.FastMCP = _StubServer
    monkeypatch.setitem(sys.modules, "mcp", types.ModuleType("mcp"))
    monkeypatch.setitem(sys.modules, "mcp.server", types.ModuleType("mcp.server"))
    monkeypatch.setitem(sys.modules, "mcp.server.fastmcp", fake)

    from lagh.mcp.server import build_server
    build_server()
    assert set(registered) == {"recover", "verify", "fit"}

    import numpy as np
    X = np.exp(np.random.default_rng(0).uniform(0, 1, (120, 1)))
    out = registered["fit"](X.tolist(), (3 * X[:, 0] ** 2).tolist())
    assert "certified" not in out and out["tag"] == "exploratory"


def test_recover_active_drives_the_oracle_and_reports_acquisition():
    """Active mode: give recover a live oracle + box; lagh runs the acquisition loop
    (ranging + adaptive queries) and returns the certificate WITH provenance."""
    oracle = lambda X: 3 * X[:, 0] ** 2                        # noqa: E731
    r = recover(oracle=oracle, box=[[0.5], [4.0]])
    assert r["certified"] is True and r["strength"] == "pinned"
    assert r["acquisition"]["mode"] == "active"
    assert r["acquisition"]["queries_used"] > 0                # it actually queried


def test_recover_box_search_reports_boxes_and_heldout_guard():
    oracle = lambda X: 2.0 * X[:, 0] * X[:, 1]                 # noqa: E731
    r = recover(oracle=oracle, box=[[0.5, 0.5], [3.0, 3.0]], box_search=True)
    assert r["certified"] is True
    assert r["acquisition"]["mode"] == "box-search"
    assert "boxes_tried" in r["acquisition"] and "heldout_box_ok" in r["acquisition"]


def test_recover_data_abstain_offers_a_broadened_box_for_the_caller_loop():
    import numpy as np
    X, y = _data(lambda X: np.round(X[:, 0] * 7) % 3, 1, 0.5, 4.0)
    r = recover(X, y)
    assert r["certified"] is False
    # abstain now carries a hedged characterization + a research-move pointer (the middle
    # rung of the degradation ladder), alongside the broadened box for the caller loop.
    assert "suggested_box" in r and "characterization" in r
    ch = r["characterization"]
    assert ch["certified"] is False and "law" not in ch          # never a certificate
    assert r["next_action"] == ch["research"]["move"]
    assert r["next_action"] in ("acquire_more_data", "acquire_divergent",
                                "declare_and_verify", "report_and_stop")
