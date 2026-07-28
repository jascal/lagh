"""SRBench entry for lagh (certified symbolic law discovery).

Contract objects: `est` (sklearn-compatible regressor), `model(est, X)`
(sympy-compatible string), `eval_kwargs`.
"""
from lagh.sklearn import LaghRegressor

est = LaghRegressor(max_time=3600)


def model(est, X=None):                                        # noqa: ARG001
    return est.model()


eval_kwargs = {}
