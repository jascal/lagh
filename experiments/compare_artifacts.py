"""Numerical JSON reproducibility contract; never changes scientific verdicts."""
import argparse
import json
import math
from pathlib import Path

RTOL, ATOL = 1e-12, 1e-15


def differences(a, b, path='$', *, significant_digits=None):
    if significant_digits is not None and not 1 <= significant_digits <= 17:
        raise ValueError('significant_digits must be in [1,17]')
    if type(a) is not type(b):
        return [path+': type differs']
    if isinstance(a, dict):
        if a.keys()!=b.keys():
            return [path+': keys differ']
        return [d for k in a for d in differences(a[k],b[k],path+'.'+k,significant_digits=significant_digits)]
    if isinstance(a, list):
        if len(a)!=len(b):
            return [path+': length differs']
        return [d for k,(x,y) in enumerate(zip(a,b)) for d in differences(x,y,f'{path}[{k}]',significant_digits=significant_digits)]
    if isinstance(a,float):
        # Each p-significant-digit rounding contributes at most
        # .5*10**(1-p) relative error. Allow both storage errors in
        # addition to the raw numerical tolerance; never use this profile
        # implicitly for an unrounded artifact such as C4.
        quantum = 0. if significant_digits is None else 10.**(1-significant_digits)
        ok=math.isfinite(a) and math.isfinite(b) and math.isclose(
            a,b,rel_tol=RTOL+quantum,abs_tol=ATOL)
    else:
        ok=a==b
    return [] if ok else [path+': value differs']


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('original',type=Path)
    parser.add_argument('rerun',type=Path)
    parser.add_argument('--significant-digits',type=int,choices=range(1,18))
    args=parser.parse_args()
    found=differences(json.loads(args.original.read_text()),json.loads(args.rerun.read_text()),significant_digits=args.significant_digits)
    print(json.dumps(dict(raw_rtol=RTOL,atol=ATOL,significant_digits=args.significant_digits,differences=found),indent=2))
    raise SystemExit(bool(found))
