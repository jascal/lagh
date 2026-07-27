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




def test_twin_end_to_end():
    """system_id -> Twin -> answers vs analytic ground truth on a dev orbit."""
    from experiments.gravitybench.twin import Twin, system_id
    m1t, m2t, a_t, e_t = 2e30, 1e30, 1.5e11, 0.3
    tb, P_true = make_circularish(m1=m1t, m2=m2t, a=a_t, ecc=e_t)
    obs, _ = _run(tb, 4.0 * P_true)
    st = system_id(obs)
    assert abs(st["q"] - 2.0) < 0.1, st["q"]
    assert abs(st["M"] - 3e30) / 3e30 < 0.05, st["M"]
    assert st["p"] == -2.0 and st["tau"] is None
    tw = Twin(st, 4.0 * P_true)
    # twin validated against the observations it was fitted from
    assert tw.validate(obs) < 0.05, tw.validate(obs)
    checks = {
        "period": (P_true, 0.03),
        "total_mass": (3e30, 0.06),
        "mass_star1": (m1t, 0.1),
        "eccentricity": (e_t, 0.1),
        "semi_major_axis": (a_t, 0.05),
        "apoastron": (a_t * (1 + e_t), 0.05),
        "periastron": (a_t * (1 - e_t), 0.05),
        "reduced_mass": (m1t * m2t / 3e30, 0.1),
    }
    for task, (truth, tol) in checks.items():
        got = tw.answer(task)
        assert got is not None and abs(got - truth) / abs(truth) < tol, \
            (task, got, truth)
    assert tw.answer("is_bound") is True
    assert tw.answer("kepler_3rd_law") is True
    # max velocity of star2 at periastron (analytic vis-viva, COM frame)
    import numpy as np
    G = 6.67430e-11
    vrel_peri = np.sqrt(G * 3e30 * (2 / (a_t * (1 - e_t)) - 1 / a_t))
    v2_peri = vrel_peri * m1t / 3e30
    got = tw.answer("max_velocity_star2")
    assert abs(got - v2_peri) / v2_peri < 0.08, (got, v2_peri)


def test_twin_modified_gravity():
    from experiments.gravitybench.twin import Twin, system_id
    tb, P_k = make_circularish(ecc=0.35, mod_gravity_exponent=-2.5)
    obs, _ = _run(tb, 4.0 * P_k)
    st = system_id(obs)
    assert abs(st["alpha"] - 0.5) < 0.2, st["alpha"]
    tw = Twin(st, 4.0 * P_k)
    assert tw.answer("kepler_3rd_law") is False
    assert abs(tw.answer("modified_gravity_power_law") - 0.5) < 0.2


if __name__ == "__main__":
    test_period_and_masses(); print("period+masses ok")
    test_gravity_exponent_newtonian(); print("newton exponent ok")
    test_gravity_exponent_modified(); print("modified exponent ok")
    test_drag_tau(); print("drag ok")
    test_twin_end_to_end(); print("twin ok")
    test_twin_modified_gravity(); print("twin modgrav ok")
