"""Scikit-learn compatible wrapper (the SRBench contract).

SRBench expects an estimator with fit/predict, a controllable `max_time` and
`random_state`, and a `model()` string in sympy syntax. lagh's two-track
output maps as: the conjecture track always answers (their harness requires
a model on every problem); certificate status rides along as metadata
(`track_`, `tag_`, `alpha_log10_`) — the abstention discipline becomes a
reported column rather than a refusal.
"""
from __future__ import annotations

import signal

import numpy as np


class LaghRegressor:
    """Certified law discovery as a sklearn-style regressor.

    Parameters
    ----------
    sigma : declared relative noise scale (0 = clean/exact data).
    max_time : soft wall-clock budget in seconds (SIGALRM guard; on timeout
        the affine-OLS fallback conjecture is returned).
    random_state : passive re-split seed (discovery itself is deterministic).
    """

    def __init__(self, sigma: float = 0.0, max_time: int = 3600,
                 random_state: int = 0):
        self.sigma = sigma
        self.max_time = max_time
        self.random_state = random_state

    def get_params(self, deep=True):                           # noqa: ARG002
        return {"sigma": self.sigma, "max_time": self.max_time,
                "random_state": self.random_state}

    def set_params(self, **params):
        for k, v in params.items():
            setattr(self, k, v)
        return self

    def _affine_fallback(self, X, y):
        A = np.column_stack([X, np.ones(len(X))])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        terms = [f"({float(c[j])!r})*x_{j}" for j in range(X.shape[1])]
        return " + ".join(terms + [f"({float(c[-1])!r})"])

    def fit(self, X, y):
        X = np.asarray(X, float)
        y = np.asarray(y, float).ravel()
        self.n_features_in_ = X.shape[1]
        from .submit import submission

        result = {}

        def _run():
            result.update(submission(X, y, sigma=self.sigma,
                                     seed=int(self.random_state)))

        if self.max_time and hasattr(signal, "SIGALRM"):
            def _timeout(signum, frame):                       # noqa: ARG001
                raise TimeoutError

            old = signal.signal(signal.SIGALRM, _timeout)
            signal.alarm(int(self.max_time))
            try:
                _run()
            except TimeoutError:
                result.clear()
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old)
        else:
            _run()

        if not result or not result.get("expr"):
            self.expr_str_ = self._affine_fallback(X, y)
            self.track_ = "conjecture"
            self.tag_ = "empirical"
            self.alpha_log10_ = None
            self.detail_ = "timeout/abstain -> affine-OLS fallback"
        else:
            self.expr_str_ = str(result["expr"])
            self.track_ = result.get("track")
            self.tag_ = result.get("tag")
            self.alpha_log10_ = result.get("alpha_log10")
            self.detail_ = result.get("detail", "")
        import sympy as sp
        self._expr = sp.sympify(self.expr_str_)
        syms = [sp.Symbol(f"x_{j}") for j in range(self.n_features_in_)]
        self._fn = sp.lambdify(syms, self._expr, "numpy")
        return self

    def predict(self, X):
        X = np.asarray(X, float)
        with np.errstate(all="ignore"):
            out = self._fn(*[X[:, j] for j in range(self.n_features_in_)])
        out = np.broadcast_to(np.asarray(out, float), (len(X),)).astype(float)
        return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)

    def model(self):
        return self.expr_str_


def model(est, X=None):                                        # noqa: ARG001
    """SRBench hook: sympy-compatible string of the final model."""
    return est.model()


est = LaghRegressor()
