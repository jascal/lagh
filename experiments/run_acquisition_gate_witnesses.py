"""Before/after P3 gates: actual failing inputs, one test file per process."""
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
BASELINE = '5f64d9a'
TEST = 'test_refusal_acquisition.py'


def replay(old):
    with tempfile.TemporaryDirectory(prefix='lagh-acquisition-gates-') as tmp:
        target = Path(tmp)
        if old:
            archive = subprocess.run(['git','archive',BASELINE,'lagh'],cwd=ROOT,
                                     capture_output=True,check=True).stdout
            with tarfile.open(fileobj=io.BytesIO(archive)) as source:
                source.extractall(target,filter='data')
        else:
            shutil.copytree(ROOT/'lagh',target/'lagh',
                            ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copyfile(ROOT/'tests'/TEST,target/TEST)
        env = dict(os.environ,PYTHONPATH=str(target),OPENBLAS_NUM_THREADS='1',
                   OMP_NUM_THREADS='1')
        run = subprocess.run([sys.executable,'-m','pytest','-q',TEST,'--junitxml=out.xml'],
                             cwd=target,env=env,capture_output=True,text=True)
        if not (target/'out.xml').exists():
            raise RuntimeError(run.stdout+run.stderr)
        rows = {}
        for case in ET.parse(target/'out.xml').iter('testcase'):
            if case.find('error') is not None or case.find('skipped') is not None:
                raise RuntimeError(run.stdout+run.stderr)
            failed = case.find('failure')
            rows[case.attrib['name']] = {'passed':failed is None,
                                        'failure':None if failed is None else failed.attrib.get('message')}
        return {'exit_code':run.returncode,'cases':rows}


def main():
    old,new = replay(True),replay(False)
    expected = {'test_initial_certificate_is_not_a_refusal_resolution',
                'test_already_refuted_model_cannot_choose_the_query',
                'test_tiny_disagreement_is_not_a_useful_query',
                'test_later_proposal_list_cannot_erase_a_live_twin'}
    actual = {k for k,v in old['cases'].items() if not v['passed']}
    if (old['exit_code'] != 1 or actual != expected or new['exit_code'] != 0
            or len(new['cases']) != 15 or not all(v['passed'] for v in new['cases'].values())):
        raise RuntimeError(f'witness prediction failed: {old}, {new}')
    data={'tag':'empirical','scope':'mocked discovery gate witnesses; not scientific reach',
          'baseline':subprocess.run(['git','rev-parse',BASELINE],cwd=ROOT,
                                    capture_output=True,text=True,check=True).stdout.strip(),
          'before':old,'after':new}
    out=ROOT/'experiments/results/acquisition_gate_witnesses.json'
    out.write_text(json.dumps(data,indent=2)+'\n')
    print('4 specified inputs fail before; all15 cases pass after',flush=True)


if __name__=='__main__':
    main()
