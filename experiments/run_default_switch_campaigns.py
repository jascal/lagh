"""Registered default-switch study, CAMPAIGN arm (docs/ESCALATION_REGISTRATION.md, A3).

Runs every real-data campaign script twice, each in its own clean git worktree
of the same commit:

  marginal  the worktree as committed (the shipped default)
  joint     the same worktree with the two defaults flipped by the registered
            one-line edit: `coefficient_gate: str = "marginal"` ->
            `coefficient_gate: str = "joint"` in lagh/engine.py and
            lagh/passive.py, nothing else

Both arms run the scripts with the same parallelism, because `recover` carries a
wall-clock time budget and load must not differ between the arms. Every result
file a script writes (git-detected) is copied to
experiments/results/default_switch/campaigns/<arm>/ together with its stdout,
and `--compare` lists every certified/law/alpha difference between the arms.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'experiments/results/default_switch/campaigns'
SCRIPTS = ['gaia/run_c0.py', 'gaia/run_p1.py', 'gaia/run_p2.py', 'gaia/run_p3.py',
           'solar/run_case_study.py', 'macro/run_case_study.py',
           'exoplanet/run_c0.py', 'exoplanet/run_c1.py', 'exoplanet/run_c2.py',
           'exoplanet/run_c5.py', 'exoplanet/run_ph2.py',
           'materials/run_c0.py', 'materials/run_c1.py', 'materials/run_c2.py']
FLIP_FROM = 'coefficient_gate: str = "marginal"'
ARMS = ('marginal', 'joint')   # A3; A4 runs ('marginal', 'joint_modulo') via --gate


def worktree(arm: str, base: Path) -> Path:
    wt = base / f'lagh-{arm}'
    if not wt.exists():
        subprocess.run(['git', 'worktree', 'add', '-q', '--detach', str(wt), 'HEAD'],
                       cwd=ROOT, check=True)
        (wt / '.venv').symlink_to(ROOT / '.venv')
        # every arm's worktree has its default set EXPLICITLY to the arm's gate
        # (the committed default changed in PR #17, so no arm is "unflipped")
        import re
        for f in ('lagh/engine.py', 'lagh/passive.py'):
            p = wt / f
            src = p.read_text()
            new_src, n = re.subn(r'coefficient_gate: str = "[a-z_]+"',
                                 f'coefficient_gate: str = "{arm}"', src)
            if n != 1:
                raise RuntimeError(f'registered flip does not apply to {f}')
            p.write_text(new_src)
    return wt


def run_one(arm: str, wt: Path, script: str) -> None:
    dest = OUT / arm / script.replace('/', '__').removesuffix('.py')
    if (dest / 'stdout.txt').exists():
        print(f'{arm} {script}: retained', flush=True)
        return
    env = dict(os.environ, PYTHONPATH=str(wt), OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1')
    before = set(subprocess.run(['git', 'status', '--porcelain', '-uall', 'experiments'],
                                cwd=wt, capture_output=True, text=True).stdout.splitlines())
    t0 = time.time()
    # CLOCK_BOOTTIME advances through a system suspend, CLOCK_MONOTONIC does
    # not: their drift over the run is the time the machine slept (a sleep
    # also eats `recover`'s wall-clock budget, so such a run is not usable)
    s0 = time.clock_gettime(time.CLOCK_BOOTTIME) - time.monotonic()
    proc = subprocess.run([str(wt / '.venv/bin/python'), f'experiments/{script}'],
                          cwd=wt, env=env, capture_output=True, text=True)
    seconds = round(time.time() - t0, 1)
    slept = round(time.clock_gettime(time.CLOCK_BOOTTIME) - time.monotonic() - s0, 1)
    after = subprocess.run(['git', 'status', '--porcelain', '-uall', 'experiments'],
                           cwd=wt, capture_output=True, text=True).stdout.splitlines()
    dest.mkdir(parents=True, exist_ok=True)
    for line in sorted(set(after) - before):
        path = line[3:]
        if path.endswith('.json') or path.endswith('.jsonl'):
            shutil.copy(wt / path, dest / Path(path).name)
    (dest / 'timing.json').write_text(json.dumps({'seconds': seconds,
                                                  'suspended_seconds': slept,
                                                  'exit': proc.returncode}) + '\n')
    (dest / 'stdout.txt').write_text(proc.stdout + f'\n[exit {proc.returncode}]\n'
                                     + proc.stderr[-4000:])
    print(f'{arm} {script}: exit {proc.returncode} {seconds}s slept={slept}s', flush=True)


def _walk(obj, path=''):
    if isinstance(obj, dict):
        if 'certified' in obj or 'law' in obj:
            yield path, {k: obj.get(k) for k in ('certified', 'law', 'alpha_log10',
                                                 'abstain')}
        for k, v in obj.items():
            yield from _walk(v, f'{path}/{k}')
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk(v, f'{path}[{i}]')


def collect(base: Path) -> None:
    """Copy each worktree's FINAL result files (every results file git sees as
    changed) to <arm>/final/. Per-script capture is not used for comparison:
    scripts run concurrently in one worktree, so a git-status diff taken around
    one script also sees files its neighbours wrote (measured on the first
    rerun: 41 spurious 'differences')."""
    for arm in ARMS:
        wt = base / f'lagh-{arm}'
        lines = subprocess.run(['git', 'status', '--porcelain', '-uall',
                                'experiments/results'], cwd=wt, capture_output=True,
                               text=True).stdout.splitlines()
        dest = OUT / arm / 'final'
        dest.mkdir(parents=True, exist_ok=True)
        for line in lines:
            path = line[3:]
            if path.endswith('.json') or path.endswith('.jsonl'):
                shutil.copy(wt / path, dest / Path(path).name)


def compare() -> None:
    first, second = ARMS
    report = {'arms': list(ARMS), 'differences': [], f'certified_{first}': 0,
              f'certified_{second}': 0, 'files': []}
    for f in sorted((OUT / first / 'final').glob('*.json')):
        other = OUT / second / 'final' / f.name
        a = dict(_walk(json.loads(f.read_text())))
        b = dict(_walk(json.loads(other.read_text()))) if other.exists() else {}
        report['files'].append(f.name)
        report[f'certified_{first}'] += sum(v['certified'] is True for v in a.values())
        report[f'certified_{second}'] += sum(v['certified'] is True for v in b.values())
        for key in sorted(set(a) | set(b)):
            if a.get(key) != b.get(key):
                report['differences'].append({'file': f.name, 'key': key,
                                              first: a.get(key),
                                              second: b.get(key)})
    timing = {}
    for arm in ARMS:
        for t in (OUT / arm).glob('*/timing.json'):
            timing.setdefault(t.parent.name, {})[arm] = json.loads(t.read_text())
    report['timing'] = timing
    (OUT / 'compare.json').write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'timing'}, indent=1,
                     default=str))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', required=True, help='directory for the two worktrees')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--compare', action='store_true')
    ap.add_argument('--gate', default='joint',
                    help='second arm: the default the worktree is flipped to')
    ap.add_argument('--first', default='marginal',
                    help='first arm (A6 compares joint_modulo with joint_quotient)')
    ap.add_argument('--out', help='results directory (default: the A3 location)')
    ap.add_argument('--concurrent', action='store_true',
                    help='run both arms at once, so they share the same load (A4)')
    a = ap.parse_args()
    global ARMS, OUT
    ARMS = (a.first, a.gate)
    if a.out:
        OUT = ROOT / a.out
    if a.compare:
        collect(Path(a.base))
        compare()
        return
    base = Path(a.base)
    trees = {arm: worktree(arm, base) for arm in ARMS}
    if a.concurrent:
        with ThreadPoolExecutor(2 * a.jobs) as pool:
            # interleave the arms so each script runs beside its twin
            work = [(arm, s) for s in SCRIPTS for arm in ARMS]
            list(pool.map(lambda w: run_one(w[0], trees[w[0]], w[1]), work))
        return
    for arm in ARMS:
        with ThreadPoolExecutor(a.jobs) as pool:
            list(pool.map(lambda s: run_one(arm, trees[arm], s), SCRIPTS))


if __name__ == '__main__':
    main()
