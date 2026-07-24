"""R-noise re-confirmation on the post-capability instrument, PASSIVE regime
(docs/RNOISE_STUDY.md protocol, corrected structural scoring).

Relative Gaussian noise at three levels, correctly-declared sigma, one fixed n=250
loguniform dataset per cell (the passive-sweep seeds), discover_passive. Each
certified law is scored against the CLEAN-data recovery (the 87-cell truth set from
newtonbench_passive.jsonl / loguniform):

- STRUCTURAL match: strip each additive term's multiplicative coefficient, keep
  exponents/atoms; a certified law with the wrong structure is a STRUCTURAL
  confident-wrong (the R-noise gate metric -- must be 0 at 60/40 dB).
- GROSS numeric wrong: certified law's predictions off the CLEAN oracle values by
  > 3x the injected noise anywhere on the dataset (the no-fabrication metric).
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import numpy as np, sympy as sp
from lagh.passive import discover_passive
from lagh.adapters.newtonbench import MODULES, available_versions, make_oracle
from experiments.run_newtonbench_passive import make_dataset, N_POINTS  # noqa: F401

LEVELS = [0.001, 0.01, 0.1]


def structure(expr_str):
    """Frozen-coefficient structural signature: the set of additive terms with
    multiplicative numeric coefficients stripped."""
    try:
        e = sp.expand(sp.sympify(expr_str))
    except Exception:                                          # noqa: BLE001
        return None
    sig = set()
    terms = e.as_ordered_terms() if e.is_Add else [e]
    for t in terms:
        c, rest = t.as_coeff_Mul()
        sig.add(sp.srepr(sp.expand(rest)))
    return frozenset(sig)


def main():
    truth = {}
    for line in open("experiments/results/newtonbench_passive.jsonl"):
        r = json.loads(line)
        if r["sampling"] == "loguniform" and r["certified"] and r["correct"]:
            truth[(r["module"], r["difficulty"], r["version"])] = r["law"]
    out = Path("experiments/results/rnoise_passive.jsonl")
    rows = []
    for sigma in LEVELS:
        for mi, (module, (inputs, lo, hi)) in enumerate(MODULES.items()):
            dim = len(inputs)
            for di, diff in enumerate(["easy", "medium", "hard"]):
                for vi, v in enumerate(available_versions(module, diff)):
                    key = (module, diff, v)
                    oracle = make_oracle(module, v, diff)
                    seed = 100000 + 1000 * mi + 100 * di + vi     # passive-sweep seeds
                    X, y = make_dataset(oracle, lo, hi, dim, "loguniform", seed)
                    rng = np.random.default_rng(seed + 7)
                    yn = y * (1.0 + sigma * rng.standard_normal(len(y)))
                    t0 = time.time()
                    r = discover_passive(X, yn, sigma=sigma, n_resplits=1)
                    rec = {"sigma": sigma, "module": module, "difficulty": diff,
                           "version": v, "certified": bool(r.certified),
                           "law": str(r.result.expr) if r.certified else None,
                           "has_truth": key in truth,
                           "seconds": round(time.time() - t0, 1)}
                    if r.certified:
                        m = np.isfinite(y) & np.isfinite(yn)
                        syms = [sp.Symbol(f"x_{i}") for i in range(dim)]
                        try:
                            f = sp.lambdify(syms, r.result.expr, "numpy")
                            pred = np.broadcast_to(
                                np.asarray(f(*X[m].T), float), y[m].shape)
                            relerr = np.abs(pred - y[m]) / np.maximum(np.abs(y[m]), 1e-30)
                            rec["max_relerr_vs_clean"] = float(np.nanmax(relerr))
                            rec["gross_wrong"] = bool(np.nanmax(relerr) > 3 * sigma)
                        except Exception:                          # noqa: BLE001
                            rec["gross_wrong"] = None
                        if key in truth:
                            rec["structural_match"] = \
                                structure(rec["law"]) == structure(truth[key])
                    rows.append(rec)
                    mk = ("ok" if rec.get("structural_match") else
                          "STRUCT-WRONG" if rec.get("structural_match") is False else
                          ("cert" if rec["certified"] else "-"))
                    print(f"s={sigma:<6} {module:22s} {diff:6s} {v} {mk:12s} "
                          f"{rec['seconds']}s", flush=True)
        out.write_text("\n".join(json.dumps(r) for r in rows) + "\n")

    print("\n== R-NOISE PASSIVE SUMMARY ==")
    for sigma in LEVELS:
        sub = [r for r in rows if r["sigma"] == sigma]
        witht = [r for r in sub if r["has_truth"]]
        cert = [r for r in witht if r["certified"]]
        sm = [r for r in cert if r.get("structural_match")]
        sw = [r for r in cert if r.get("structural_match") is False]
        gross = [r for r in sub if r.get("gross_wrong")]
        print(f"sigma={sigma}: truth-cells={len(witht)} certified={len(cert)} "
              f"structural-ok={len(sm)} STRUCTURAL-CW={len(sw)} "
              f"gross-wrong(>3x noise)={len(gross)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
