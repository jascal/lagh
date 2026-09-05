import json

import numpy as np
import sympy as sp

from lagh.base import eval_expr
from lagh.mcp.core import verify
from lagh.refusal import residual_measurement


def measure(y, pred, eps, indices):
    return residual_measurement(y, pred, eps, indices, domain='checked rows',
                                candidate='2*x_0')


def test_public_residual_corresponds_to_original_rows_and_candidate():
    x = np.linspace(1, 2, 100)
    y = 2*x + .1*x*x
    x[0] = np.nan
    r = verify(x, y, 'x_0', floor_abs=1e-6)
    assert not r['certified'] and r['abstain'] == 'structural'
    assert 'partial' not in r and 'law' not in r
    m = r['measurement']
    indices = np.array(m['row_indices'])
    assert 0 not in indices and len(indices) == 20
    pred = eval_expr(sp.sympify(m['candidate']), [sp.Symbol('x_0')], x[indices, None])
    np.testing.assert_allclose(m['residual'], y[indices] - pred, rtol=0, atol=0)
    assert m['attribution'] == 'unresolved' and m['evidence'] == 'empirical'
    assert m['domain'] == 'certification rows of supplied dataset'
    assert m['n_exceeding'] == 20 and m['max_band_excess'] > 0
    json.dumps(r, allow_nan=False)


def test_full_domain_and_clean_control():
    x = np.linspace(1, 2, 100)
    y = 2*x
    # Frozen fixture: RNG(0)'s original first two fitting rows, not recomputed.
    a, b = 82, 36
    y[a] += .1
    y[b] -= .1*x[a]/x[b]
    r = verify(x, y, 'x_0', floor_abs=1e-6)
    assert not r['certified'] and 'partial' not in r
    m = r['measurement']
    assert m['domain'] == 'all finite rows of supplied dataset'
    assert m['n_checked'] == 100 and m['n_exceeding'] == 2
    assert set(m['row_indices'][:2]) == {a, b}
    clean = verify(x, 2*x, 'x_0', floor_abs=1e-6)
    assert clean['certified'] and clean['law'] == '2*x_0'
    assert 'measurement' not in clean


def test_known_signed_residuals_bands_and_excess():
    m = measure([3., -4., 0.], [1., -1., 0.], [.5, 1., .1], [7, 9, 12])['measurement']
    assert m['row_indices'] == [9, 7, 12]
    assert m['residual'] == [-3., 2., 0.]
    assert m['epsilon'] == [1., .5, .1]
    assert m['max_band_excess'] == 2. and m['n_exceeding'] == 2


def test_undefined_rows_preserve_finite_neighbors():
    m = measure([1., 2., 3., 4.], [0., np.nan, 0., 0.],
                [.1, .1, np.inf, -.1], [7, 8, 9, 10])['measurement']
    assert m['row_indices'] == [7] and m['residual'] == [1.]
    assert m['n_invalid'] == 3 and m['n_measurable'] == 1
    json.dumps(m, allow_nan=False)


def test_bad_alignment_and_unsupported_band_explain_omission():
    for indices in ([7], [1, 1, 2], [-1, 2, 3], [1., 2., 3.]):
        assert 'measurement_omitted' in measure([1, 2, 3], [0, 0, 0], .1, indices)
    assert 'measurement_omitted' in measure([1], [0], lambda _: .1, [0])
    assert 'measurement_omitted' in measure([1], None, .1, [0])
    assert 'measurement_omitted' in measure([1e308], [-1e308], .1, [0])


def test_payload_is_capped_without_losing_miss_counts():
    y = np.zeros(20000)
    y[-100:] = 2
    m = measure(y, np.zeros_like(y), .1, np.arange(len(y)))['measurement']
    assert len(m['residual']) == 64 and m['n_exceeding'] == 100
    assert m['n_elided'] == 19936 and m['n_invalid'] == 0
    assert min(m['row_indices']) >= 19900
    assert m['max_band_excess'] == 1.9
    assert len(json.dumps(m)) < 6000


def test_exact_band_boundary_is_not_a_miss():
    m = measure([1., -1.], [0., 0.], 1., [0, 1])['measurement']
    assert m['n_exceeding'] == 0 and m['max_band_excess'] == 0
