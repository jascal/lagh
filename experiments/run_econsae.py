"""Run the registered econ-sae targets E1-E4 through lagh (docs/TESTBED_ECONSAE.md).

Active acquisition per target; dense-grid reference scoring; zero-wrong invariant
checked. Ground truth is the simulator itself -- never a formula -- so scoring is honest
even where no closed form exists.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import sympy as sp

from lagh.acquisition import run_active
from lagh.adapters.econsae import TARGETS, make_oracle


def dense_reference_ok(expr, oracle, box, dim, tol=1e-3, n=200) -> bool | None:
    """Correct iff the recovered law matches a fresh dense grid to `tol` relative.
    Returns None if expr is None (abstained)."""
    if expr is None:
        return None
    lo, hi = np.array(box[0]), np.array(box[1])
    rng = np.random.default_rng(12345)
    X = np.exp(rng.uniform(np.log(lo), np.log(hi), (n, dim)))
    y = oracle(X)
    syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
    got = sp.lambdify(syms, expr, "numpy")(*X.T)
    got = np.broadcast_to(np.asarray(got, float), y.shape)
    ok = np.isfinite(got) & np.isfinite(y) & (np.abs(y) > 1e-12)
    if ok.sum() < n // 2:
        return False
    return bool(np.max(np.abs(got[ok] - y[ok]) / np.abs(y[ok])) < tol)


def main() -> int:
    out = Path("experiments/results"); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for tid, spec in TARGETS.items():
        dim = len(spec["inputs"])
        oracle = make_oracle(tid)
        t0 = time.time()
        r = run_active(oracle, spec["box"][0], spec["box"][1], seed=1)
        cert = r.result.certificate
        law = str(r.result.expr) if r.result.expr is not None else None
        correct = dense_reference_ok(r.result.expr, oracle, spec["box"], dim)
        rec = {"target": tid, "aggregate": spec["aggregate"],
               "certified": cert.certified, "abstain": cert.abstain,
               "law": law, "dense_ref_correct": correct,
               "confident_wrong": bool(cert.certified and correct is False),
               "queries": r.queries_used, "tier": r.result.tier,
               "box_initial": r.box_initial.tolist(),
               "box_final": np.array(r.box_final).tolist(),
               "ranging_contractions": len(r.ranging_trajectory) - 1,
               "seconds": round(time.time() - t0, 1)}
        rows.append(rec)
        status = ("CERT " + ("OK" if correct else "**WRONG**")) if cert.certified \
            else f"abstain[{cert.abstain}]"
        print(f"{tid} {spec['aggregate']:16s} {status:16s} "
              f"q={r.queries_used} tier={r.result.tier} "
              f"contract={rec['ranging_contractions']} ({rec['seconds']}s)",
              flush=True)
        if law:
            print(f"     law: {law[:80]}", flush=True)

    (out / "econsae.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    cw = sum(r["confident_wrong"] for r in rows)
    cert = sum(r["certified"] for r in rows)
    print(f"\ncertified {cert}/{len(rows)} | confident-wrong {cw} "
          f"(invariant: must be 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
