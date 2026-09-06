"""Wrong-law witnesses through consumers; includes the rejected aliasing design."""
from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import sympy as sp
from lagh import engine, passive, acquisition
from lagh.certify import check, Certificate


def mutation(legacy_aliasing=False):
    y=np.array([1.,3.]); eps=np.array([.1,.1]); ids=np.array([8,2])
    r=check(sp.S.One,[],np.empty((2,0)),y,eps,row_indices=ids)
    if legacy_aliasing:
        # Explicit reconstruction of the rejected intermediate implementation;
        # never a production option or an empirical test credited as a pass.
        r.observed, r.epsilon, r.row_indices = y, eps, ids
    y[:]=1; eps[:]=100; ids[:]=[0,1]
    return dict(input='check 1 vs [1,3] at .1, ids[8,2]; mutate y=1, eps=100, ids[0,1]',
                verdict=dict(r),measurement_after_mutation=r.measurement())


def main():
    x=sp.Symbol('x_0')
    X=np.linspace(1,2,100)[:,None]; y=2*X[:,0]
    proposed=SimpleNamespace(expr=x,complexity=1,channel='witness')
    with patch.object(engine,'_tier_candidates',return_value=[proposed]):
        result=engine.discover(X[:60],y[:60],X[60:80],y[60:80],X[80:],y[80:],max_tier=1)
    def false_certificate(*args,**kwargs):
        return engine.Result(Certificate(True,0,0,20,[(1,2)],str(x)),x,1,1)
    with patch.object(passive,'discover',side_effect=false_certificate):
        p=passive.discover_passive(X,y,n_resplits=1)
    a=SimpleNamespace(result=false_certificate(),box_final=[[1],[2]])
    accepted=acquisition._heldout_box_ok(lambda X:2*X[:,0],a,1e-12,0)
    out=dict(evidence='empirical',input='only proposed law x_0; observations 2*x_0 on [1,2]',
        engine=asdict(result.certificate),passive=asdict(p.result.certificate),
        acquisition=dict(heldout_passed=accepted,measurement=a.result.certificate.measurement),
        mutation=mutation())
    Path('experiments/results/residual_consumers.json').write_text(json.dumps(out,indent=2)+'\n')
    Path('experiments/results/residual_mutation_baseline.json').write_text(json.dumps(mutation(True),indent=2)+'\n')
    print('engine',result.certificate.measurement is not None,
          'passive',p.result.certificate.measurement is not None,
          'acquisition',a.result.certificate.measurement is not None)


if __name__=='__main__':
    main()
