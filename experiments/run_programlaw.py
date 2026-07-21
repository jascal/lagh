"""Run the registered program-law targets F1-F10 (docs/TESTBED_PROGRAMLAW.md).

Each function is a black-box oracle; its source (== mathematical definition) is used
ONLY for the dense-grid reference at scoring time, never to guide recovery. The
candidate list and per-function domain boxes are frozen in the registration.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np

from lagh.acquisition import run_active

# (id, callable, arity, box_lo, box_hi) -- frozen per registration
FUNCS = [
    ("F1_hypot", lambda X: np.hypot(X[:, 0], X[:, 1]), 2, [0.1, 0.1], [100, 100]),
    ("F2_log2", lambda X: np.log2(X[:, 0]), 1, [0.1], [100]),
    ("F3_gamma", lambda X: np.vectorize(math.gamma)(X[:, 0]), 1, [0.5], [6.0]),
    ("F4_erf", lambda X: np.vectorize(math.erf)(X[:, 0]), 1, [0.1], [3.0]),
    ("F5_atan2", lambda X: np.arctan2(X[:, 1], X[:, 0]), 2, [0.1, 0.1], [10, 10]),
    ("F6_logaddexp", lambda X: np.logaddexp(X[:, 0], X[:, 1]), 2, [0.1, 0.1], [5, 5]),
    ("F7_expm1", lambda X: np.expm1(X[:, 0]), 1, [0.1], [3.0]),
    ("F8_copysign", lambda X: np.copysign(X[:, 0], X[:, 1] - 5), 2, [0.1, 0.1], [10, 10]),
    ("F9_remainder", lambda X: np.remainder(X[:, 0], X[:, 1]), 2, [0.1, 0.1], [10, 10]),
    ("F10_hypot3", lambda X: np.sqrt(X[:, 0]**2 + X[:, 1]**2 + X[:, 2]**2), 3,
     [0.1]*3, [100]*3),
]


def dense_ref_ok(expr, oracle, lo, hi, dim, tol=1e-9, n=200):
    if expr is None:
        return None
    import sympy as sp
    rng = np.random.default_rng(777)
    X = np.exp(rng.uniform(np.log(lo), np.log(hi), (n, dim)))
    y = oracle(X)
    syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
    try:
        got = sp.lambdify(syms, expr, "numpy")(*X.T)
        got = np.broadcast_to(np.asarray(got, float), y.shape)
    except Exception:                                     # noqa: BLE001
        return False
    ok = np.isfinite(got) & np.isfinite(y) & (np.abs(y) > 1e-12)
    if ok.sum() < n // 2:
        return False
    return bool(np.max(np.abs(got[ok] - y[ok]) / np.abs(y[ok])) < tol)


def main() -> int:
    out = Path("experiments/results"); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for fid, fn, dim, lo, hi in FUNCS:
        t0 = time.time()
        r = run_active(fn, lo, hi, seed=1)
        expr = r.result.expr
        correct = dense_ref_ok(expr, fn, np.array(lo), np.array(hi), dim)
        rec = {"func": fid, "arity": dim, "certified": r.result.certificate.certified,
               "abstain": r.result.certificate.abstain,
               "law": str(expr) if expr is not None else None,
               "dense_ref_correct": correct,
               "confident_wrong": bool(r.result.certificate.certified
                                       and correct is False),
               "queries": r.queries_used, "tier": r.result.tier,
               "seconds": round(time.time() - t0, 1)}
        rows.append(rec)
        flag = ("CERT " + ("ok" if correct else "**WRONG**")) \
            if r.result.certificate.certified else f"abstain[{r.result.certificate.abstain}]"
        print(f"{fid:14s} {flag:16s} tier={r.result.tier} q={r.queries_used} "
              f"({rec['seconds']}s)", flush=True)
        if rec["law"]:
            print(f"     law: {rec['law'][:80]}", flush=True)

    (out / "programlaw.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    rec_ok = sum(r["certified"] and r["dense_ref_correct"] for r in rows)
    ab = sum(r["abstain"] is not None for r in rows)
    cw = sum(r["confident_wrong"] for r in rows)
    print(f"\nrecovered {rec_ok} | abstained {ab} | confident-wrong {cw} (must be 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
