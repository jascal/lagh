"""Registered oracle-owned acquisition pilot; P0 before numerical scored promises."""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import zlib

import numpy as np

from lagh.certify import free_dof
from lagh.measurement_design import choose_measurement
from lagh.passive import discover_passive


OUT = Path('experiments/results/acquisition_pilot_p0.json')


def main():
    seed = zlib.crc32(b'rational-d1')
    X = np.random.default_rng(seed).uniform(.5, 3., (400, 1))
    y = (2*X[:, 0] + 1)/(X[:, 0] + 3)
    print('P0: historical apex draw, 400 design observations', flush=True)
    result = discover_passive(X, y, sigma=0., seed=0)
    r = result.result
    choice = choose_measurement(r.rivals, np.geomspace(.005, 300., 257)[:, None],
                                bounds=[[.005], [300.]])
    data = {
        'tag': 'empirical', 'stage': 'pilot P0, design only',
        'seed': seed, 'design_queries': len(X), 'final_queries': 0,
        'initial_bounds': [[float(X.min())], [float(X.max())]],
        'certified': result.certified, 'law': str(r.expr),
        'abstain': r.certificate.abstain, 'notes': r.certificate.notes,
        'rivals': [{'expr': str(e), 'dof': free_dof(e)} for e in r.rivals],
        'choice': asdict(choice),
        'P0': not result.certified and len(r.rivals) >= 2,
        'P0b': choice.point is not None and not .5 <= choice.point[0] <= 3.,
    }
    OUT.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    print(f"P0={data['P0']}, P0b={data['P0b']}, "
          f"rivals={len(r.rivals)}, point={choice.point}", flush=True)


if __name__ == '__main__':
    main()
