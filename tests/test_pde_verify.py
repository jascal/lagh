"""The PDE verify track (experiments/pde/verify.py, docs/CASE_STUDY_PDE_C2.md).

A weak-form certificate is a claim about patch integrals; these check the
machinery that turns it into a claim about DYNAMICS, and in particular that the
band declares every error the forecast introduces -- solver, parameter interval,
and the initial condition's own noise.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from experiments.pde.verify import (ic_noise_bound, integrate,  # noqa: E402
                                    solver_bound, verify)

X = np.linspace(0.0, 2 * np.pi, 128, endpoint=False)
T = np.linspace(0.0, 0.3, 11)


def heat(nu=0.1):
    return np.exp(-nu * T)[None, :] * np.sin(X)[:, None]


def test_forecast_reproduces_an_exact_solution():
    u = heat()
    got = integrate(u[:, 0], X, T, {"u_xx": 0.1}, rtol=1e-10)
    assert got is not None
    assert np.max(np.abs(got - u)) < 1e-9


def test_solver_bound_is_scalar_and_covers_the_ladder():
    u = heat()
    y, b = solver_bound(u[:, 0], X, T, {"u_xx": 0.1}, rtol=1e-9)
    assert isinstance(b, float) and b > 0
    assert np.max(np.abs(y - u)) <= b          # the declared bound covers


def test_exact_law_verifies_and_a_wrong_one_does_not():
    u = heat()
    ok = verify(u, u[:, 0], X, T, {"u_xx": (0.1, 0.1)}, rtol=1e-9)
    assert ok["verified"] and ok["n_outside"] == 0
    bad = verify(u, u[:, 0], X, T, {"u_xx": (0.11, 0.11)}, rtol=1e-9)
    assert not bad["verified"]


def test_interval_certificate_forecasts_an_envelope():
    """The certificate claims a family of laws, so the forecast is a family: a
    wider interval must widen the band, never move the verdict by luck."""
    u = heat()
    narrow = verify(u, u[:, 0], X, T, {"u_xx": (0.0999, 0.1001)}, rtol=1e-9)
    wide = verify(u, u[:, 0], X, T, {"u_xx": (0.099, 0.101)}, rtol=1e-9)
    assert wide["envelope_width_med"] > 5 * narrow["envelope_width_med"]
    assert narrow["verified"] and wide["verified"]


def test_ic_noise_is_declared_not_assumed_away():
    """A forecast started from measured data carries that data's noise for the
    whole window (measured: 2.8 sigma for advection, which is non-dissipative,
    and enough to fail every noisy rung when the band omitted it)."""
    u = heat()
    assert ic_noise_bound(u[:, 0], X, T, {"u_xx": 0.1}, 0.0) == 0.0
    b1 = ic_noise_bound(u[:, 0], X, T, {"u_xx": 0.1}, 1e-6, rtol=1e-9)
    b2 = ic_noise_bound(u[:, 0], X, T, {"u_xx": 0.1}, 1e-5, rtol=1e-9)
    assert b1 > 0 and 5 < b2 / b1 < 20          # linear in sigma
    # advection carries it further than heat: nothing damps it
    adv = ic_noise_bound(u[:, 0], X, T, {"u_x": -0.7}, 1e-5, rtol=1e-9)
    assert adv > b2
