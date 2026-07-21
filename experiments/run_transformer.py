"""Run the registered transformer-idiom targets (docs/TESTBED_TRANSFORMER.md).

Train grokked + undertrained checkpoints, recover the idiom via C6 from argmax queries,
score against the true modulus and an extended-range check. Weights are never read.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

from lagh import discover
from tiny_transformer import oracle_fn, train

P = 11
SLICES = [0, 3, 7]
T_MAX = 48
EXT = list(range(49, 89))


def recover_slice(model, b0):
    oracle = oracle_fn(model, b0)
    a = np.arange(1, T_MAX + 1, dtype=float)[:, None]
    y = oracle(a).astype(float)
    i = np.random.default_rng(0).permutation(T_MAX)
    f, s = 30, 39
    r = discover(a[i[:f]], y[i[:f]], a[i[f:s]], y[i[f:s]], a[i[s:]], y[i[s:]])
    if not r.certificate.certified:
        return {"certified": False, "abstain": r.certificate.abstain}
    period = getattr(r.expr, "period", None)
    # extended-range check vs the model itself (never a formula)
    ext_ok = all(int(r.expr(t)) == int(oracle(np.array([[t]]))[0]) for t in EXT)
    return {"certified": True, "period": period, "note": r.certificate.notes[0],
            "period_correct": period == P, "ext_ref_ok": ext_ok,
            "idiom_correct": bool(period == P and ext_ok)}


def main() -> int:
    out = Path("experiments/results"); out.mkdir(parents=True, exist_ok=True)
    rows = []
    checkpoints = {"grokked": train(P, 4000, seed=0),
                   "undertrained": train(P, 50, seed=0)}
    for name, model in checkpoints.items():
        acc = model._train_acc
        for b0 in SLICES:
            res = recover_slice(model, b0)
            rec = {"checkpoint": name, "exact_acc": round(acc, 4),
                   "b0": b0, **res,
                   "confident_wrong": bool(res.get("certified")
                                           and res.get("idiom_correct") is False)}
            rows.append(rec)
            if res["certified"]:
                flag = "CERT " + ("idiom ok" if res.get("idiom_correct") else "**WRONG**")
                extra = f"period={res.get('period')}"
            else:
                flag = f"abstain[{res['abstain']}]"
                extra = ""
            print(f"{name:12s} acc={acc:.3f} b0={b0}: {flag:18s} {extra}", flush=True)

    (out / "transformer.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    gk = [r for r in rows if r["checkpoint"] == "grokked"]
    ut = [r for r in rows if r["checkpoint"] == "undertrained"]
    print(f"\ngrokked: {sum(r.get('idiom_correct', False) for r in gk)}/{len(gk)} idiom recovered")
    print(f"undertrained: {sum(not r['certified'] for r in ut)}/{len(ut)} abstained")
    print(f"confident-wrong: {sum(r['confident_wrong'] for r in rows)} (must be 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
