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

    return server


def main():
    build_server().run()


if __name__ == "__main__":
    main()
