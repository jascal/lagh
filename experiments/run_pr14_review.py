"""Review counter-inputs; run --baseline before changes, then without it."""
from types import SimpleNamespace
from unittest.mock import patch
import argparse
import json
from pathlib import Path
import numpy as np
import sympy as sp
from lagh.certify import check, Certificate
from lagh import engine
from lagh.mcp import core
from experiments.compare_artifacts import differences
from lagh.instrument import faxen_height, gate_record, covariance_timescale


def run():
    x=sp.Symbol('x_0')
    values=np.linspace(1,2,5)
    out=[]
    def observe(name, fn):
        try:
            result=fn()
        except Exception as exc:
            result={'exception':type(exc).__name__, 'message':str(exc)}
        out.append({'input':name,'observed':result})
    for label, X in [('1D',values),('scalar',np.array(1.)),('3D',values[:,None,None]),
                     ('wrong columns',np.ones((5,2))),('non-numeric',np.full((5,1),'bad'))]:
        observe('check X='+label,lambda X=X:dict(check(x,[x],X,2*values,1e-6)))
    for ids in [[0,1],[-1,0,1,2,3],[0,0,1,2,3]]:
        observe('check row_indices='+repr(ids),lambda ids=ids:dict(check(x,[x],values[:,None],2*values,1e-6,row_indices=ids)))
    for F in [1.,np.nextafter(1.,2.),1.+1e-12,1.+1e-9,1.5]:
        for se in [None,0.,.1,float('nan'),float('inf'),-.1]:
            observe(f'faxen F={F!r}, radius=1, se={se!r}',lambda F=F,se=se:faxen_height(F,1.,se))
    observe('gate_record ratios=1, no variance',lambda:gate_record({'a':{'ratio':1.},'b':{'ratio':1.}}))
    observe('GLS zeros(100)',lambda:covariance_timescale(np.zeros(100),.001))
    xdata=np.linspace(1,2,30)[:,None]
    checked=check(x,[x],xdata,2*xdata[:,0],.01)
    cert=Certificate(False,30,0,30,[(1,2)],str(x),abstain='structural')
    cert.measurement=checked.measurement(domain='fixture')['measurement']
    result=engine.Result(cert,None,1,1)
    with patch.object(core,'discover_passive',return_value=SimpleNamespace(result=result)):
        observe('recover passive failed x against 2x, retained measurement',
                lambda: core.recover(xdata,2*xdata[:,0],max_tier=1))
    active=SimpleNamespace(result=result,box_final=[[1],[2]],queries_used=30,
                           ledger=SimpleNamespace(spent=30),ranging_trajectory=[])
    with patch.object(core,'run_active',return_value=active), patch.object(core,'_characterize_oracle',return_value=None):
        observe('recover active failed x against 2x, retained measurement',
                lambda:core.recover(oracle=lambda X:2*X[:,0],box=[[1],[2]]))
    observe('verify exact 2x with NaN band',lambda:core.verify(xdata,2*xdata[:,0],'x_0',floor_abs=float('nan')))
    v=1.234567890125
    neighbor=np.nextafter(v,np.inf if format(v,'.12g')=='1.23456789012' else -np.inf)
    a,b=(float(format(z,'.12g')) for z in (v,neighbor))
    observe(f'one-ULP neighbors {v!r},{neighbor!r}, rounded to {a!r},{b!r}; raw profile',
            lambda:differences(a,b))
    observe('same rounded pair; 12-significant-digit profile',lambda:differences(a,b,significant_digits=12))
    observe('rounded profile, unit-sized value perturbed1e-5',lambda:differences(a,a*1.00001,significant_digits=12))
    return {'evidence':'empirical','witnesses':out}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',action='store_true')
    parser.add_argument('--historical-worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.baseline:
        import io, os, shutil, subprocess, sys, tarfile, tempfile
        archive=subprocess.check_output(['git','archive','150c617','lagh','experiments/compare_artifacts.py'])
        with tempfile.TemporaryDirectory() as tmp:
            with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
                tar.extractall(tmp,filter='data')
            output=Path(tmp)/'experiments/results'
            output.mkdir(parents=True)
            subprocess.run([sys.executable,str(Path(__file__).resolve()),'--historical-worker'],
                           cwd=tmp,env=dict(os.environ,PYTHONPATH=tmp),check=True)
            shutil.copyfile(output/'pr14_review_baseline.json',
                            'experiments/results/pr14_review_baseline.json')
        raise SystemExit(0)
    target='pr14_review_baseline.json' if args.historical_worker else 'pr14_review.json'
    # Nonfinite inputs are named textually. Nonfinite outputs, if any, remain
    # explicit strings so the baseline is valid JSON, never sanitized into passes.
    def safe(v):
        if isinstance(v,dict): return {k:safe(x) for k,x in v.items()}
        if isinstance(v,list): return [safe(x) for x in v]
        if isinstance(v,float) and not np.isfinite(v): return repr(v)
        return v
    Path('experiments/results',target).write_text(json.dumps(safe(run()),indent=2)+'\n')
    print(target)
