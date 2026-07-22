"""A pluggable problem lab for exercising the LLM + lagh loop.

The lab server speaks a fixed MCP contract -- `problems` (cards: inputs, domain, hint)
and `sample` (query the hidden oracle) -- so the LLM loop and scoring are
benchmark-agnostic. The benchmark behind it is a swappable `Source` (`LAB_SOURCE`):

    proxy        6 self-contained numpy problems -- runs anywhere, spans the honesty
                 spectrum (recoverable / needs-broader-box / non-smooth / irrational)
    newtonbench  the real NewtonBench-dev cells (needs the newtonbench adapter)

Adding a benchmark = adding a `Source` in sources.py; nothing else changes. Ground-truth
forms live in the sources for OFFLINE scoring only and are never exposed over MCP.
"""

from .sources import Source, available, get_source

__all__ = ["get_source", "available", "Source"]
