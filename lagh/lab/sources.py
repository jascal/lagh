"""Pluggable problem sources: swap the benchmark behind ONE interface.

The lab server speaks a fixed contract -- `problems()` (cards) and `sample(pid, X)`
(query the hidden oracle) -- so the LLM loop and scoring are benchmark-agnostic. A
`Source` binds that contract to a specific benchmark; pick one with `LAB_SOURCE`
(default `proxy`). Adding a new benchmark = adding a `Source`, nothing else changes.

    proxy        6 self-contained problems (pure numpy) -- runs anywhere, no deps
    newtonbench  the real NewtonBench-dev cells -- needs the newtonbench adapter (openai)
"""

from __future__ import annotations

import os

from . import problems as _proxy


class Source:
    """A benchmark bound to the lab contract. `truth` is OFFLINE scoring only."""

    name = "?"

    def problems(self) -> list[dict]:      # public cards, no answers
        raise NotImplementedError

    def sample(self, pid: str, X) -> list[float]:
        raise NotImplementedError

    def truth(self, pid: str) -> str:      # never exposed over MCP
        return "(hidden)"


class ProxySource(Source):
    name = "proxy"

    def problems(self):
        return [_proxy.card(p) for p in _proxy.PROBLEMS]

    def sample(self, pid, X):
        return _proxy.sample(pid, X)

    def truth(self, pid):
        return _proxy.truth(pid)


class NewtonBenchSource(Source):
    """The real NewtonBench-dev cells behind the same interface. Ids are
    'module/difficulty/version'. Requires the newtonbench adapter (imports openai), so
    it only loads where those deps are present -- the proxy source runs anywhere."""

    name = "newtonbench"

    def __init__(self, difficulties=("easy", "medium", "hard")):
        from ..adapters.newtonbench import MODULES, available_versions
        self._M = MODULES
        self._ids = [f"{m}/{d}/{v}" for m in MODULES for d in difficulties
                     for v in available_versions(m, d)]

    def problems(self):
        out = []
        for pid in self._ids:
            m, d, _ = pid.split("/")
            inp, lo, hi = self._M[m]
            out.append({"id": pid, "inputs": list(inp),
                        "suggested_domain": [list(map(float, lo)), list(map(float, hi))],
                        "hint": f"{m.split('_', 1)[1].replace('_', ' ')} ({d})",
                        "n_inputs": len(inp)})
        return out

    def sample(self, pid, X):
        import numpy as np
        from ..adapters.newtonbench import make_oracle
        m, d, v = pid.split("/")
        orc = make_oracle(m, v, d)
        Xa = np.asarray(X, float)
        Xa = Xa[:, None] if Xa.ndim == 1 else Xa
        return [float(t) for t in np.asarray(orc(Xa), float).ravel()]

    def truth(self, pid):
        # best-effort docstring formula (host-side scoring); else score by re-querying
        try:
            import glob
            import re
            m, d, v = pid.split("/")
            for f in glob.glob(f"*NewtonBench*/modules/{m}/laws.py") or \
                    glob.glob(f"../*NewtonBench*/modules/{m}/laws.py"):
                src = open(f).read()
                mt = re.search(rf'_ground_truth_law_{d}_{v}\(.*?"""(.*?)"""', src, re.S)
                if mt:
                    return mt.group(1).strip().splitlines()[-1].strip()
        except Exception:                                     # noqa: BLE001
            pass
        return f"{pid} (score by re-querying the oracle)"


_REGISTRY = {ProxySource.name: ProxySource, NewtonBenchSource.name: NewtonBenchSource}


def available() -> list[str]:
    return list(_REGISTRY)


def get_source(name: str | None = None) -> Source:
    name = name or os.environ.get("LAB_SOURCE", "proxy")
    if name not in _REGISTRY:
        raise ValueError(f"unknown LAB_SOURCE {name!r}; known: {available()}")
    return _REGISTRY[name]()
