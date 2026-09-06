"""Registered apex diagnostic using P1's historical refusal-model heuristic.

The old models are already contradicted by original data. This is an engine
reach-ordering experiment, explicitly not a surviving-twin separation result.
"""
from dataclasses import asdict
import importlib.util
import json
import subprocess
import sys
import zlib

import numpy as np
import sympy as sp

from experiments.run_refusal_comparison import ARMS, ROOT, score
from lagh.certify import check,epsilon
from lagh.passive import discover_passive

PIN = '5f64d9a'


def main():
    source = subprocess.run(['git','show',f'{PIN}:lagh/refusal_acquisition.py'],cwd=ROOT,
                            capture_output=True,text=True,check=True).stdout
    name='lagh._registered_p1_acquisition'
    spec=importlib.util.spec_from_loader(name,loader=None)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    exec(compile(source,f'{PIN}:lagh/refusal_acquisition.py','exec'),module.__dict__)
    seed=zlib.crc32(b'rational-d1')
    X=np.random.default_rng(seed).uniform(.5,3.,(400,1))
    y=(2*X[:,0]+1)/(X[:,0]+3)
    x=sp.Symbol('x_0'); truth=(2*x+1)/(x+3)
    initial=discover_passive(X,y,sigma=0.,seed=0)
    witnesses=[dict(check(e,[x],X,y,epsilon(y))) for e in initial.result.rivals]
    if initial.certified or len(witnesses)!=2 or any(w['certified'] for w in witnesses):
        raise RuntimeError('historical initial state changed; prediction needs rescoring')
    directory=ROOT/'experiments/results/acquisition_apex'
    directory.mkdir(exist_ok=True)
    for arm in ARMS:
        path=directory/f'{arm}.json'
        if path.exists():
            print(f'{arm}: completed artifact retained',flush=True)
            continue
        plans=directory/f'{arm}.plans.jsonl'
        if plans.exists():
            raise RuntimeError('unfinished apex ledger; inspect before restarting')
        count=0
        def oracle(points):
            nonlocal count
            count+=len(points)
            return (2*points[:,0]+1)/(points[:,0]+3)
        policy=module.RefusalPolicy(((.005,),(300.,)))
        with plans.open('x') as stream:
            def record(plan):
                stream.write(json.dumps(plan,allow_nan=False)+'\n');stream.flush()
                print(f'{arm}: {plan["role"]}, next{plan["count"]}, spent{plan["query_start"]}',flush=True)
            out=module.run_refusal_acquisition(oracle,X,y,initial_box=[[.5],[3.]],
                                               policy=policy,seed=99,strategy=arm,
                                               initial_result=initial,on_plan=record)
        if out.queries !=400+count:
            raise RuntimeError('oracle ledger mismatch')
        truth_score=score(out.result.expr,truth) if out.certified else score(None,truth)
        data={'tag':'empirical','scope':'apex reach-ordering diagnostic; old rivals already refuted',
              'initial_data_seed':seed,'followup_seed':99,'policy':asdict(policy),
              'pinned_runner':PIN,'arm':arm,'initial_rival_checks':witnesses,
              'certified':out.certified,'law':str(out.result.expr) if out.certified else None,
              'queries':out.queries,'actual_new_queries':count,'domain':out.domain,
              'final_passed':out.final_passed,'truth_score':truth_score,
              'wrong_certified':out.certified and not truth_score['exact_truth'],
              'history':out.history,'notes':out.result.certificate.notes}
        path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
        print(f'{arm}: certified={out.certified}, queries={out.queries}, '
              f'wrong={data["wrong_certified"]}',flush=True)


if __name__=='__main__':
    main()
