"""Numerical JSON reproducibility contract; never changes scientific verdicts."""
import argparse
import json
import math
from pathlib import Path

RTOL, ATOL = 1e-12, 1e-15


def differences(a, b, path='$'):
    if type(a) is not type(b):
        return [path+': type differs']
    if isinstance(a, dict):
        if a.keys()!=b.keys():
            return [path+': keys differ']
        return [d for k in a for d in differences(a[k],b[k],path+'.'+k)]
    if isinstance(a, list):
        if len(a)!=len(b):
            return [path+': length differs']
        return [d for k,(x,y) in enumerate(zip(a,b)) for d in differences(x,y,f'{path}[{k}]')]
    if isinstance(a,float):
        ok=math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=RTOL,abs_tol=ATOL)
    else:
        ok=a==b
    return [] if ok else [path+': value differs']


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('original',type=Path)
    parser.add_argument('rerun',type=Path)
    args=parser.parse_args()
    found=differences(json.loads(args.original.read_text()),json.loads(args.rerun.read_text()))
    print(json.dumps(dict(rtol=RTOL,atol=ATOL,differences=found),indent=2))
    raise SystemExit(bool(found))
