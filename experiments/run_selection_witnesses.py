"""Exhibit selector and serialization counter-inputs before crediting tests."""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from lagh.ito_selection import choose_test
from lagh.certify import check
from experiments.compare_artifacts import differences


def main():
    out=[]
    for label,t,p in [
        ('constant',np.arange(10.),np.ones((2,10))),
        ('nonfinite',np.arange(10.),np.full((2,10),np.nan)),
        ('irregular',np.arange(10.)**2,np.tile(np.arange(10.),(2,1)))]:
        try:
            choose_test(t,p)
            out.append(dict(input=label,rejected=False))
        except ValueError as e:
            out.append(dict(input=label,rejected=True,reason=str(e)))
    x=sp.Symbol('x')
    out.append(dict(input='wrong drift +x on noiseless y=-x, x=1..10, eps=.01',
        result=dict(check(x,[x],np.arange(1.,11.)[:,None],-np.arange(1.,11.),.01))))
    out.append(dict(input='unit float + 1e-5',differences=differences({'x':1.},{'x':1.00001})))
    out.append(dict(input='unit float + 1 ULP',differences=differences({'x':1.},{'x':float(np.nextafter(1.,2.))})))
    Path('experiments/results/selection_artifact_witnesses.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
