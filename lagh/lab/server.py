"""MCP server for the problem lab: the queryable oracle the model discovers against.

Exposes `problems` (cards) and `sample` (query the hidden oracle) over a fixed contract;
the benchmark BEHIND that contract is chosen by `LAB_SOURCE` (default `proxy`; also
`newtonbench`). It never exposes ground-truth forms. Thin FastMCP wrapper; SDK imported
lazily. Run: `LAB_SOURCE=proxy python -m lagh.lab.server`.
"""

from __future__ import annotations

import os

from .sources import get_source


def build_server():
    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as e:                                   # noqa: BLE001
        raise SystemExit("the MCP SDK is not installed. `pip install lagh[mcp]`.") from e

    src = get_source(os.environ.get("LAB_SOURCE"))
    server = FastMCP(f"lab", instructions=(
        f"A lab of hidden law-discovery problems (source: {src.name}). `problems` lists "
        "them (inputs, a suggested domain, a hint); `sample` queries a problem's oracle "
        "at points you choose. You set up the sampling, then hand the (X, y) to the lagh "
        "tools (`fit` to scout, `recover` to certify). The true law is never revealed -- "
        "discover it, or abstain honestly like lagh does."))

    @server.tool()
    def problems() -> list[dict]:
        """List the lab's problems: each has an `id`, `inputs`, a `suggested_domain`
        ([lo, hi] per input), and a one-line `hint`. No answers -- you discover them."""
        return src.problems()

    @server.tool()
    def sample(problem: str, X: list[list[float]]) -> dict:
        """Query a problem's hidden oracle. `problem` is an id from `problems`; `X` is a
        list of input rows (each row has n_inputs values). Returns the oracle's `y` per
        row. Sample enough points over a wide-enough range to pin the law (lagh needs
        >= 8 total and a non-degenerate box). You may sample OUTSIDE the suggested
        domain -- that's how you broaden when `fit` says acquire_more_data."""
        try:
            y = src.sample(problem, X)
        except Exception as e:                                # noqa: BLE001
            return {"error": str(e)}
        return {"problem": problem, "source": src.name, "n": len(y), "y": y}

    @server.tool()
    def discover(problem: str) -> dict:
        """DELEGATE the full discovery to lagh. lagh runs its OWN active-acquisition loop
        -- adaptive ranging, budgeted multi-objective queries, and a broaden-the-box
        ladder on abstain -- directly against the problem's oracle, and returns a
        certificate (certified law + strength + acquisition provenance) or a reasoned
        abstention. Use this when your own sample+recover abstains: lagh samples the
        oracle far more effectively than a manual pass, so it often certifies where a
        hand-picked set of points cannot. `problem` is an id from `problems`."""
        import numpy as np
        from lagh.mcp.core import recover as _recover
        cards = {c["id"]: c for c in src.problems()}
        if problem not in cards:
            return {"error": f"unknown problem {problem!r}"}
        lo, hi = cards[problem]["suggested_domain"]

        def oracle(X):
            return np.asarray(src.sample(problem, np.asarray(X, float).tolist()), float)

        r = _recover(oracle=oracle, box=[lo, hi], box_search=True)
        r["problem"] = problem
        return r

    return server


def main():
    build_server().run()


if __name__ == "__main__":
    main()
