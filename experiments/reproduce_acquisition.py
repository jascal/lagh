"""Re-execute one scored condition in a fresh directory and compare all four arms.

Copies no completed output into the temporary run: this checks reproduction,
not the runner's skip-if-already-complete path. Model code must match S1's hash.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import scipy
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    with tempfile.TemporaryDirectory(prefix='lagh-acquisition-reproduce-') as folder:
        target=Path(folder)
        shutil.copytree(ROOT/'lagh',target/'lagh',ignore=shutil.ignore_patterns('__pycache__'))
        exp=target/'experiments';exp.mkdir()
        for source in (ROOT/'experiments').glob('*.py'):
            shutil.copyfile(source,exp/source.name)
        shutil.copyfile(ROOT/'experiments/acquisition_scored_protocol.json',
                        exp/'acquisition_scored_protocol.json')
        results=exp/'results';results.mkdir()
        shutil.copytree(ROOT/'experiments/results/acquisition_p2',results/'acquisition_p2')
        env=dict(os.environ,PYTHONPATH=str(target),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        run=subprocess.run([sys.executable,'-m','experiments.run_survivor_comparison',
                            '--stage','scored','--design-seed','30','--domain','wide',
                            '--followup-seed','100'],cwd=target,env=env,
                           capture_output=True,text=True)
        if run.returncode:
            raise RuntimeError(run.stdout+run.stderr)
        generated=results/'acquisition_scored'
        comparisons=[]
        for new in sorted(generated.iterdir()):
            old=ROOT/'experiments/results/acquisition_scored'/new.name
            comparisons.append({'file':new.name,'original_sha256':digest(old),
                                'reproduced_sha256':digest(new),'byte_identical':digest(old)==digest(new)})
        if len(comparisons)!=8 or not all(c['byte_identical'] for c in comparisons):
            raise RuntimeError(json.dumps(comparisons,indent=2))
    sources=[]
    fixed=[ROOT/'lagh/acquisition.py',ROOT/'lagh/certify.py',
           *sorted((ROOT/'lagh/classes').glob('*.py'))]
    for path in fixed:
        rel=str(path.relative_to(ROOT))
        old=subprocess.run(['git','show',f'e810b98:{rel}'],cwd=ROOT,
                            capture_output=True,check=True).stdout
        same=old==path.read_bytes()
        if not same:
            raise RuntimeError(f'fixed baseline or hypothesis class changed: {rel}')
        sources.append({'file':rel,'unchanged_from_base':same,'sha256':digest(path)})
    data={'tag':'empirical','scope':'one fresh scored condition, all four arms and plan ledgers',
          'comparisons':comparisons,'unchanged_sources':sources,
          'environment':{'python':platform.python_version(),'numpy':np.__version__,
                         'scipy':scipy.__version__,'sympy':sp.__version__,
                         'machine':platform.machine(),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}}
    out=ROOT/'experiments/results/acquisition_reproduction.json'
    out.write_text(json.dumps(data,indent=2)+'\n')
    print('8/8 artifacts byte-identical after actual fresh execution; baseline and class sources unchanged')


if __name__=='__main__':
    main()
