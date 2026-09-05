import json

import numpy as np

from lagh.mcp.core import verify
from lagh.refusal import residual_measurement


def test_refusal_preserves_signed_rows_and_domain():
    x = np.linspace(1, 2, 100)
    y = 2*x + .1*x*x
    x[0] = np.nan
    r = verify(x, y, 'x_0', floor_abs=1e-6)
    assert not r['certified'] and r['abstain'] == 'structural'
    m = r['measurement']
    indices = np.array(m['row_indices'])
    assert 0 not in indices and len(indices) == 20
    assert m['attribution'] == 'unresolved'
    assert m['domain'] == 'certification rows of supplied dataset'
    residual = np.array(m['residual'])
    # Every retained discrepancy exceeds the declared band exactly when the gate
    # counts a miss, and the partial envelope contains every per-row interval.
    assert np.any(np.abs(residual) > m['epsilon'])
    envelope = r['partial']['components']['row_discrepancy_envelope']
    assert np.all(residual - m['epsilon'] >= envelope['lo'])
    assert np.all(residual + m['epsilon'] <= envelope['hi'])
    json.dumps(r, allow_nan=False)


def test_full_domain_and_clean_control():
    x = np.linspace(1, 2, 100)
    y = 2*x
    a, b = np.random.default_rng(0).permutation(100)[:2]
    y[a] += .1
    y[b] -= .1*x[a]/x[b]
    r = verify(x, y, 'x_0', floor_abs=1e-6)
    assert not r['certified']
    assert r['measurement']['domain'] == 'all finite rows of supplied dataset'
    assert r['measurement']['row_indices'] == list(range(100))
    clean = verify(x, 2*x, 'x_0', floor_abs=1e-6)
    assert clean['certified'] and clean['law'] == '2*x_0'
    assert 'measurement' not in clean


def test_nonfinite_measurement_is_omitted():
    assert residual_measurement([1e308], [-1e308], [1], [0], domain='row') == {}
    assert residual_measurement([1], None, [1], [0], domain='row') == {}
