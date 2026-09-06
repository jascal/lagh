"""Pilot-frozen selection and independent evaluation on an Itô reach ladder."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp

from lagh.ito_selection import choose_test
from lagh.ito import build_rows, certify_drift, ItoBand
from lagh.certify import check
from experiments.stochastic.generator import ou_paths, gbm_paths


def paths(case, seed):
    if case=='OU':
        return ou_paths(theta=1.,b=.5,T=256,dt=.002,n_traj=8,seed=seed)
    return gbm_paths(mu=.015 if case=='GBM' else 0.,b=.2,T=128,dt=.002,n_traj=8,seed=seed)


def main():
    out=[]
    for case in ('OU','GBM','null'):
        for pilot_seed in range(100,104):
            t,pilot=paths(case,pilot_seed)
            choice=choose_test(t,pilot)
            # Creation of certification data occurs ONLY AFTER the choice.
            ct,cert=paths(case,pilot_seed+100)
            cert_hash=hashlib.sha256(ct.tobytes()+cert.tobytes()).hexdigest()
            assert choice.pilot_sha256!=cert_hash
            result={}
            for label,f in [('baseline','x**2/2'),('selected',choice.function)]:
                rows=build_rows(ct,cert,('1','x','x**2'),fs=(f,),
                                half=4000 if case=='OU' else 8000)
                verdict=certify_drift(rows,delta=.05,seed=0)
                coefficient=-1. if case=='OU' else .015 if case=='GBM' else 0.
                truth=sp.Float(coefficient)*sp.Symbol('x_1')
                band=ItoBand(rows.qv,rows.qv_se,rows.corr_se,rows.quad,len(rows.y),
                             delta=.05,y=rows.y,feat_names=list(rows.names))
                checked=check(truth,band.syms,rows.A,rows.y,band)
                components=verdict.get('partial',{}).get('components',{})
                excluded=[]
                for name,value in zip(('1','x','x**2'),(0.,coefficient,0.)):
                    c=components.get(name,{})
                    lo,hi=c.get('lo'),c.get('hi')
                    if (lo is not None and lo>value) or (hi is not None and hi<value):
                        excluded.append(name)
                result[label]=dict(function=f,certified=verdict['certified'],
                    abstain=verdict.get('abstain'),partial=verdict.get('partial'),
                    true_law_misses=checked['nmiss'],true_law_uncovered=checked['nuncov'],
                    truth_exclusions=excluded,n_rows=len(rows.y),
                    median_signal_to_band=verdict.get('median_signal_to_band'))
            out.append(dict(case=case,pilot_seed=pilot_seed,certification_seed=pilot_seed+100,
                            pilot=asdict(choice),certification_sha256=cert_hash,results=result))
            print(case,pilot_seed,choice.function,
                  [(k,v['certified'],v['true_law_misses'],v['truth_exclusions']) for k,v in result.items()],flush=True)
    def rounded(v):
        if isinstance(v,dict): return {k:rounded(a) for k,a in v.items()}
        if isinstance(v,(list,tuple)): return [rounded(a) for a in v]
        if isinstance(v,float): return float(format(v,'.12g')) if np.isfinite(v) else None
        return v
    Path('experiments/results/test_selection.json').write_text(json.dumps(rounded(
        dict(evidence='empirical',runs=out,delta_per_run=.05,
             claim='Pilot-only proposal; empirical reach, no new coverage theorem.')),indent=2)+'\n')


if __name__=='__main__':
    main()
