"""Replay transport counter-inputs against pinned preparation and current code.

One test file per process. These mocked-split unit witnesses measure retention,
not acquisition reach. No scientific oracle is queried.
"""
from __future__ import annotations

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
BASELINE = 'cef9fee'
TEST = 'tests/test_refusal_transport.py'
OUT = ROOT / 'experiments/results/refusal_transport_witnesses.json'


def replay(baseline: bool) -> dict:
    with tempfile.TemporaryDirectory(prefix='lagh-refusal-transport-') as folder:
        target = Path(folder)
        if baseline:
            archive = subprocess.run(['git', 'archive', BASELINE, 'lagh'],
                                     cwd=ROOT, check=True, capture_output=True).stdout
            with tarfile.open(fileobj=io.BytesIO(archive)) as source:
                source.extractall(target, filter='data')
        else:
            shutil.copytree(ROOT / 'lagh', target / 'lagh',
                            ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copyfile(ROOT / TEST, target / 'test_refusal_transport.py')
        env = dict(os.environ, PYTHONPATH=str(target), OPENBLAS_NUM_THREADS='1',
                   OMP_NUM_THREADS='1')
        run = subprocess.run([sys.executable, '-m', 'pytest', '-q',
                              'test_refusal_transport.py', '--junitxml=results.xml'],
                             cwd=target, env=env, capture_output=True, text=True)
        xml = target / 'results.xml'
        if not xml.exists():
            raise RuntimeError(run.stdout + run.stderr)
        cases = {}
        for item in ET.parse(xml).iter('testcase'):
            failed = item.find('failure')
            if item.find('error') is not None or item.find('skipped') is not None:
                raise RuntimeError(run.stdout + run.stderr)
            cases[item.attrib['name']] = {
                'passed': failed is None,
                'failure': None if failed is None else failed.attrib.get('message'),
            }
        return {'exit_code': run.returncode, 'cases': cases}


def main():
    before = replay(True)
    after = replay(False)
    if (len(before['cases']) != 5 or before['exit_code'] != 1
            or any(c['passed'] for c in before['cases'].values())
            or after['exit_code'] != 0 or len(after['cases']) != 5
            or not all(c['passed'] for c in after['cases'].values())):
        raise RuntimeError(f'Counter-input prediction failed: {before}, {after}')
    revision = subprocess.run(['git', 'rev-parse', BASELINE], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    OUT.write_text(json.dumps({'tag': 'empirical', 'baseline_revision': revision,
                               'scope': 'mocked-split transport only; no acquisition reach',
                               'before': before, 'after': after}, indent=2) + '\n')
    print(f'5/5 fail on baseline, 5/5 pass after repair: {OUT.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
