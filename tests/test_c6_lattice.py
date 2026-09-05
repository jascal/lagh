"""The evaluation interface a recovered QuasiPoly needs to reach public consumers
(jascal/lagh#6)."""

import numpy as np
import sympy as sp

from lagh.base import eval_expr
from lagh.classes.c6_quasipoly import recover_integer
from lagh.passive import discover_passive

X0 = sp.Symbol("x_0")


def _qp(n=60, offset=100000.0):
    X = np.arange(1., n + 1.)[:, None]
    return X, offset + X[:, 0] ** 2 + (X[:, 0] % 2) * 100



def test_eval_expr_evaluates_a_recovered_quasipoly():
    t = np.arange(1, 25)
    q = recover_integer(t, t ** 2 + t % 2).quasipoly
    v = eval_expr(q, [X0], t[:, None].astype(float))
    assert np.array_equal(v, (t ** 2 + t % 2).astype(float))
    v2 = eval_expr(q, [X0], np.array([[2.0], [2.5]]))   # off-lattice: undefined, not rounded
    assert v2[0] == 4.0 and np.isnan(v2[1])


def test_passive_discovery_passes_a_quasipoly_through_the_full_data_gate():
    X, y = _qp()
    r = discover_passive(X, y, n_resplits=1, max_tier=6)
    assert r.certified and r.result.tier == 6
    assert r.full_check_passed is True
    assert str(r.result.expr).startswith("quasipoly")
