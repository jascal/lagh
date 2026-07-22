"""lagh MCP tool-shape: the three acts (docs/DIRECTION_TOOLSHAPE.md).

`recover` / `verify` (bounded, cert-or-abstain) and `fit` (unbounded scout). The
core logic is transport-free (importable and testable without the MCP SDK); the
server wrapper lives in `.server` and is imported lazily so `import lagh.mcp` never
requires the SDK.
"""

from .core import fit, recover, verify

__all__ = ["recover", "verify", "fit"]
