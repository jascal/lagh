"""Full NewtonBench-dev sweep: all difficulties x modules x versions, coverage matrix.
DEV MEASUREMENT (STRATEGY.md) -- drives capability toward the readiness bar, no win claim."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import numpy as np, sympy as sp
from lagh.acquisition import run_active
from lagh.adapters.newtonbench import MODULES, available_versions, make_oracle

def dense_ok(expr, oracle, lo, hi, dim, tol=1e-6, n=200):
    """Scored on the ORACLE's own domain: where the oracle answers (finite y), the law
    must cover >=90% of it and match to tol; where the oracle declines, the law is not
    penalized. The old ok>=n/2 rule falsely scored a correct restricted-domain law
    (snell acos, oracle NaN on most of the box) as WRONG for lack of raw coverage."""
    if expr is None: return None
    rng=np.random.default_rng(4321); lo,hi=np.array(lo),np.array(hi)
    X=np.exp(rng.uniform(np.log(np.maximum(lo,1e-6)),np.log(hi),(n,dim))); y=oracle(X)
    syms=[sp.Symbol(f"x_{i}") for i in range(dim)]
    try: got=np.broadcast_to(np.asarray(sp.lambdify(syms,expr,"numpy")(*X.T),float),y.shape)
    except Exception: return False
    val=np.isfinite(y)&(np.abs(y)>1e-9)
    if val.sum()<20: return None            # oracle too sparse here to score at all
    cov=np.isfinite(got[val])
    if cov.mean()<0.9: return False         # law undefined where the oracle answers
    ok=val.copy(); ok[val]=cov
    return bool(np.max(np.abs(got[ok]-y[ok])/np.abs(y[ok]))<tol)

def main():
    out=Path("experiments/results/newtonbench_all.jsonl"); out.parent.mkdir(parents=True,exist_ok=True)
    diffs=["easy","medium","hard"]; rows=[]
    for module,(inputs,lo,hi) in MODULES.items():
        dim=len(inputs)
        for diff in diffs:
            for v in available_versions(module,diff):
                oracle=make_oracle(module,v,diff); t0=time.time()
                r=run_active(oracle,lo,hi,seed=1)
                correct=dense_ok(r.result.expr,oracle,lo,hi,dim)
                rec={"module":module,"difficulty":diff,"version":v,
                     "certified":r.result.certificate.certified,"correct":correct,
                     "abstain":r.result.certificate.abstain,
                     "law":str(r.result.expr) if r.result.expr is not None else None,
                     "confident_wrong":bool(r.result.certificate.certified and correct is False),
                     "seconds":round(time.time()-t0,1)}
                rows.append(rec)
                mk="ok" if (rec["certified"] and correct) else ("WRONG" if rec["confident_wrong"] else "-")
                print(f"{module:22s} {diff:6s} {v} {mk:5s} {rec['seconds']}s",flush=True)
    out.write_text("\n".join(json.dumps(r) for r in rows)+"\n")
    # coverage matrix: module x difficulty (any version recovered)
    print("\nCOVERAGE (modules recovered per difficulty, any version):")
    for diff in diffs:
        mods={}
        for r in rows:
            if r["difficulty"]==diff: mods.setdefault(r["module"],[]).append(r["certified"] and r["correct"])
        rec=sum(any(v) for v in mods.values())
        print(f"  {diff:6s}: {rec}/{len(mods)} modules")
    cw=sum(r["confident_wrong"] for r in rows)
    print(f"confident-wrong: {cw} (must be 0) | total tasks {len(rows)}")
    print("readiness R-cap: >=10/12 on >=2 difficulties")
    return 0
if __name__=="__main__": raise SystemExit(main())
