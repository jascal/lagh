"""MCP transport for the three lagh acts (docs/DIRECTION_TOOLSHAPE.md).

Thin wrapper over `core.py`. The MCP SDK is imported here (not in `__init__`), so
`import lagh.mcp` never requires it; only running the server does. Install with
`pip install lagh[mcp]`, then `python -m lagh.mcp.server` (or the `lagh-mcp` script).

The tool docstrings below ARE the contract the calling model reads. They are worded
to make the bounded/unbounded distinction unmissable: `recover`/`verify` return a
certificate-or-abstention; `fit` returns conjectures that are explicitly NOT
certificates. The disjoint return shapes enforce the same wall structurally -- a
`fit` result has no `certified` field to misread.
"""

from __future__ import annotations

from . import core


def build_server():
    """Construct the FastMCP server. Imported lazily so the SDK is optional."""
    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as e:                                   # noqa: BLE001
        raise SystemExit(
            "the MCP SDK is not installed. `pip install lagh[mcp]` "
            "(or `pip install mcp`) to run the server.") from e

    server = FastMCP("lagh", instructions=(
        "lagh turns law-discovery subproblems into CERTIFIED results or explicit "
        "abstentions. Sample your black-box oracle, pass the points, and: use "
        "`recover` to discover an exact law, `verify` to check a form you already "
        "suspect, `fit` to scout the data before you commit. `recover`/`verify` are "
        "bounded (certificate or a reasoned abstention -- never a guess). `fit` is a "
        "SCOUT: its output is a conjecture with a diagnosis, NOT a certificate -- to "
        "turn a conjecture into a guarantee, feed it to `recover` or `verify`."))

    @server.tool()
    def recover(X: list[list[float]], y: list[float], sigma: float = 0.0) -> dict:
        """BOUNDED. Discover an exact closed-form law from data points (X, y).

        Returns a CERTIFICATE (`tag:"proved"`, `certified:true`, the `law`, and a
        `strength` of `pinned`) exhaustively checked over the stated finite domain,
        OR an ABSTENTION (`tag:"open"`, `certified:false`, a machine-readable
        `abstain` reason). It never fabricates: no law that fits is a refusal, not a
        wrong answer. `sigma` = relative measurement noise (0 for a clean oracle);
        under noise it certifies the exact STRUCTURE or abstains (coefficients are
        noise-limited). The certificate is over the domain you sampled, not a claim
        about the world -- check the `bounds` match your question.

        ACTIVE ACQUISITION over the wire is CALLER-DRIVEN: on a thin/under-determined
        abstain this returns `next_action:"acquire"` + a broadened `suggested_box` --
        re-sample your oracle over that box and call `recover` again. (In-process
        callers can instead hand lagh the oracle directly and let it run the full
        adaptive loop; that path can't cross the JSON wire.)
        """
        return core.recover(X, y, sigma=sigma)

    @server.tool()
    def verify(X: list[list[float]], y: list[float], form: str,
               sigma: float = 0.0) -> dict:
        """BOUNDED. Check a form YOU declare against data (X, y).

        `form` is a sympy expression in `x_0..x_{d-1}` (e.g. `"x_0*x_1/x_2**2"`); its
        overall scale is refit, then it is checked over the domain. Returns a
        CERTIFICATE or an ABSTENTION (a wrong form is REFUTED, not tolerated). The
        `strength` field is load-bearing: a rational form certifies `pinned` (this
        law, no rival within the noise); a form carrying a declared irrational (e.g.
        `x_0**E`) can only certify `consistent` -- it fits, but the constant is not
        identifiable from finite data. Use this for hypothesis-checking, including
        irrational/known-constant forms you have an a-priori reason to believe.
        """
        return core.verify(X, y, form, sigma=sigma)

    @server.tool()
    def fit(X: list[list[float]], y: list[float], sigma: float = 0.0) -> dict:
        """UNBOUNDED SCOUT -- returns a CONJECTURE, NOT a certificate.

        Best-guess forms plus an identifiability DIAGNOSIS (`pinned` / `continuum` /
        `under_determined`) and a `next_action` pointer that names your next move.
        Its value is the diagnosis, not the guess: it tells you whether a clean
        rational is there to `recover`, whether a free exponent is a continuum
        reaching for a named constant (→ `declare_and_verify`), or whether forms tie
        and you need to `acquire_more_data`. There is deliberately NO `certified`
        field here -- nothing in this result may be reported as a guarantee. To make
        any conjecture trustworthy, pass it to `recover` or `verify`.
        """
        return core.fit(X, y, sigma=sigma)

    return server


def main():
    build_server().run()


if __name__ == "__main__":
    main()
