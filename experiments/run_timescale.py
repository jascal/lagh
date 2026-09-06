"""Registered synthetic estimator comparison; never opens raw tweezers data."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tests'))
from test_instrument import bead, ou, DT, FS, RD_UM_PER_V, KBT, G_BULK
from lagh.instrument import acf_timescale, covariance_timescale, retention


def estimate(x, filtered):
    a = acf_timescale(x, DT, min_lag=8, max_lag=200)['theta']
    b = covariance_timescale(x, DT)['theta']
    try:
        p = 2*np.pi*retention(x, FS, fit_range=(100,2300),
                             alpha=.451 if filtered else 1., f_diode=11684.).fc_fit
    except (ValueError, RuntimeError):
        p = None
    return {'legacy': a, 'gls': b, 'psd': p}


def main(initial=False):
    records=[]
    for n in (32768,131072,524288):
        for filtered in (False,True):
            for seed in list(range(20,28 if initial else 32))+list(range(40,44) if initial else range(60,64)):
                theta=.1194e-3/G_BULK
                if filtered:
                    _, x, truth=bead(n=n, seed=seed, driven=False)
                    theta=truth['theta']
                else:
                    x=ou(n=n, seed=seed, theta=theta, b2=2*KBT/G_BULK*1e18)[0]/(RD_UM_PER_V*1e3)
                for contaminated in (False,True):
                    v=x.copy()
                    if contaminated:
                        v+=2*x.std()*np.sin(2*np.pi*7*np.arange(n)*DT)
                    full=estimate(v,filtered)
                    halves=[estimate(a,filtered) for a in np.array_split(v,2)]
                    records.append(dict(n=n,filtered=filtered,contaminated=contaminated,
                                        seed=seed,split='calibration' if seed<(40 if initial else 60) else 'evaluation',
                                        truth_theta=theta,full=full,halves=halves))
            print(n,filtered,'done',flush=True)
    summaries=[]
    for n in (32768,131072,524288):
        for filtered in (False,True):
            for contaminated in (False,True):
                subset=[r for r in records if (r['n'],r['filtered'],r['contaminated'])==(n,filtered,contaminated)]
                for method in ('legacy','gls','psd'):
                    def metrics(r):
                        vals=[r['full'][method]]+[h[method] for h in r['halves']]
                        if any(v is None or v<=0 or not np.isfinite(v) for v in vals):
                            return [float('inf')]*3
                        truth= max(abs(v/r['truth_theta']-1) for v in vals)
                        half=max(vals[1:])/min(vals[1:])-1
                        psd=r['full']['psd']
                        disagree=abs(vals[0]/psd-1) if psd and psd>0 else float('inf')
                        return truth,half,disagree
                    cal=[max(metrics(r)[:2]) for r in subset if r['split']=='calibration']
                    # Frozen rank correction is registered; no estimator setting changed.
                    rank=int(np.ceil((len(cal)+1)*.9))
                    tau=sorted(cal)[rank-1] if rank<=len(cal) else float('inf')
                    ev=np.array([metrics(r) for r in subset if r['split']=='evaluation'])
                    summaries.append(dict(n=n,filtered=filtered,contaminated=contaminated,
                      method=method,evaluation_max_truth_error=float(ev[:,0].max()),
                      evaluation_max_half_disagreement=float(ev[:,1].max()),
                      evaluation_max_psd_disagreement=float(ev[:,2].max()),
                      calibration_empirical_max=max(cal), calibrated_tolerance=tau,
                      quantile_rank=rank, calibration_count=len(cal),
                      within_eight_percent=bool(np.all(ev[:,:2]<.08))))
    def finite(v):
        if isinstance(v,dict): return {k:finite(x) for k,x in v.items()}
        if isinstance(v,list): return [finite(x) for x in v]
        if isinstance(v,float): return float(format(v,'.12g')) if np.isfinite(v) else None
        return v
    out=dict(evidence='empirical',records=records,summary=summaries,
             note=('Null tolerance means no finite 90% bound at registered pilot size; no classifier threshold changed.' if initial else
                   'Per-regime calibration rank ceil(.9*(m+1)); null means unbounded. Synthetic exchangeability only; no joint classifier tolerance changed.'))
    Path('experiments/results/timescale_initial.json' if initial else 'experiments/results/timescale.json').write_text(json.dumps(finite(out),indent=2,allow_nan=False)+'\n')
    for s in summaries:
        if s['method']=='gls': print(s,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--initial',action='store_true')
    main(parser.parse_args().initial)
