"""Two-track submission policy (docs/DIRECTION_OUTPUT_POLICY.md)."""

import numpy as np

from lagh.submit import submission


def _data(fn, dim, n=200, seed=11):
    rng = np.random.default_rng(seed)
    X = np.exp(rng.uniform(np.log(0.5), np.log(10), (n, dim)))
    return X, fn(X)


def test_in_class_law_goes_track_a():
    X, y = _data(lambda X: 3 * X[:, 0] ** 2 / X[:, 1], dim=2)
    s = submission(X, y)
    assert s["track"] == "certified" and s["tag"] == "proved"
    assert s["expr"]


def test_out_of_class_goes_track_b_labeled():
    X, y = _data(lambda X: X[:, 0] ** np.e, dim=1)   # irrational exponent: no exact form
    s = submission(X, y)
    assert s["track"] == "conjecture" and s["tag"] == "empirical"
    assert s["expr"]                                 # still submits SOMETHING, labeled
    assert "NOT certified" in s["detail"]


def test_unusable_data_abstains():
    X = np.ones((30, 1))
    y = np.full(30, np.nan)
    s = submission(X, y)
    assert s["track"] == "abstain" and s["expr"] is None
