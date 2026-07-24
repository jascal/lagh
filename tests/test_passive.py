"""Passive-data mode (docs/DIRECTION_PASSIVE.md): fixed dataset, no oracle."""

import numpy as np
import sympy as sp
import pytest

import lagh.passive as passive_mod
from lagh.passive import PassiveResult, discover_passive


def _data(fn, dim, n=250, lo=0.5, hi=10.0, seed=7):
    rng = np.random.default_rng(seed)
    X = np.exp(rng.uniform(np.log(lo), np.log(hi), (n, dim)))
    return X, fn(X)


def test_recovers_polynomial_from_fixed_data():
    X, y = _data(lambda X: 3 * X[:, 0] ** 2 + 2 * X[:, 1], dim=2)
    r = discover_passive(X, y)
    assert r.certified
    x0, x1 = sp.symbols("x_0 x_1")
    assert sp.simplify(r.result.expr - (3 * x0**2 + 2 * x1)) == 0


def test_recovers_power_law_from_fixed_data():
    X, y = _data(lambda X: 6.674e-5 * X[:, 0] * X[:, 1] / X[:, 2] ** 2, dim=3)
    r = discover_passive(X, y)
    assert r.certified


def test_honest_abstain_out_of_class():
    # irrational exponent: no exact closed form exists -- must abstain, never guess
    X, y = _data(lambda X: X[:, 0] ** np.e, dim=1)
    r = discover_passive(X, y)
    assert not r.certified
    assert r.result.certificate.abstain


def test_nonfinite_rows_dropped_not_fatal():
    X, y = _data(lambda X: 5 * X[:, 0] ** 3, dim=1)
    y = y.copy()
    y[::17] = np.nan          # saturated cells, as a passive dataset may contain
    r = discover_passive(X, y)
    assert r.certified


def test_full_data_gate_demotes_split_artifact(monkeypatch):
    """A law that certifies on its split but fails any point of the full dataset must
    come back DEMOTED to abstain -- the gate that keeps K re-splits sound."""
    from lagh.certify import Certificate
    from lagh.engine import Result

    bogus = sp.sympify("x_0**2")            # wrong law for the data below

    def fake_discover(*a, **k):
        cert = Certificate(True, 0, 0, 50, [(0.5, 10.0)], str(bogus))
        return Result(cert, bogus, 1, 1)

    monkeypatch.setattr(passive_mod, "discover", fake_discover)
    X, y = _data(lambda X: X[:, 0] ** 3, dim=1)
    r = discover_passive(X, y, n_resplits=2)
    assert isinstance(r, PassiveResult)
    assert not r.certified
    assert not r.result.certificate.certified
    assert r.full_check_passed is False
    assert any("full-data" in n for n in r.result.certificate.notes)
