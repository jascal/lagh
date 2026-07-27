"""Synthetic validation of the deterministic astronomer (dev-only orbits)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.gravitybench.integrator import G_SI, TwoBody, make_circularish  # noqa: E402
from experiments.gravitybench import astronomer as ast  # noqa: E402


def _run(tb, maxtime):
    tb.run(maxtime)
    return ast.plan_and_observe(tb.observe, maxtime, budget=100, per_request=10)


def test_period_and_masses():
    tb, P_true = make_circularish(m1=2e30, m2=1e30, a=1.5e11, ecc=0.3)
    obs, _ = _run(tb, 4.0 * P_true)
    P = ast.task_period(obs)
    assert abs(P - P_true) / P_true < 0.02, (P, P_true)
    m1, m2 = ast.task_masses(obs, P)
    assert abs(m1 - 2e30) / 2e30 < 0.10, m1
    assert abs(m2 - 1e30) / 1e30 < 0.10, m2


def test_gravity_exponent_newtonian():
    tb, P_true = make_circularish(ecc=0.4)
    obs, _ = _run(tb, 4.0 * P_true)
    slope, cert = ast.task_gravity_exponent(obs)
    assert slope is not None and abs(slope - (-2.0)) < 0.15, slope


def test_gravity_exponent_modified():
    tb, P_kepler = make_circularish(ecc=0.35, mod_gravity_exponent=-2.5)
    obs, _ = _run(tb, 4.0 * P_kepler)     # period estimate approximate; fine
    slope, cert = ast.task_gravity_exponent(obs)
    assert slope is not None and abs(slope - (-2.5)) < 0.2, slope


def test_drag_tau():
    tb, P_true = make_circularish(ecc=0.3, drag_tau=None)
    tb2, _ = make_circularish(ecc=0.3)
    tau_true = 8.0 * P_true
    tb3, _ = make_circularish(ecc=0.3, drag_tau=tau_true)
    obs, _ = _run(tb3, 5.0 * P_true)
    tau = ast.task_drag_tau(obs)
    assert tau is not None and abs(tau - tau_true) / tau_true < 0.5, (tau, tau_true)


if __name__ == "__main__":
    test_period_and_masses(); print("period+masses ok")
    test_gravity_exponent_newtonian(); print("newton exponent ok")
    test_gravity_exponent_modified(); print("modified exponent ok")
    test_drag_tau(); print("drag ok")
