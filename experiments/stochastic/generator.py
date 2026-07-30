"""Level 0 systems for the stochastic suite (docs/DIRECTION_STOCHASTIC.md).

Three systems and two nulls, all seeded, and every one of them simulated EXACTLY
where an exact transition law exists. That is deliberate: Level 0 measures the
INSTRUMENT's coverage, and an Euler-Maruyama discretization bias would be
indistinguishable in the results from a mis-calibrated band. OU and geometric
Brownian motion both have exact transitions, so neither needs a scheme.

Systems
  ou            dX = -theta X dt + b dW          (linear drift, additive noise)
  gbm           dS = mu S dt + b S dW            (linear drift, multiplicative)
  ode_obs       dX = -theta X dt, observed with sigma_obs (NO process noise)

Nulls
  pure_noise    dX = dW                          -- no drift; must certify nothing
  coarse_ou     ou sampled at dt >> 1/theta      -- drift unidentifiable AT this
                                                    dt; must abstain with a
                                                    RESOLUTION reason, not a wide
                                                    interval

`ode_obs` is the rung the curriculum renamed: at sigma_process = 0 an SDE is an
ODE, so it tests existing machinery rather than stochastic discovery. It is kept
because it is the one Level 0 case where process and measurement noise are
separable BY CONSTRUCTION, and it predicts something specific: realized quadratic
variation cannot tell them apart -- observational noise contributes 2 n sigma_obs^2
to [X], which DIVERGES as dt -> 0. So the same estimator that makes the martingale
band honest for a real SDE reports a spurious diffusion here, and the size of that
spurious value is a measurement worth having (S6, and
docs/DIRECTION_ERROR_PROVENANCE.md).
"""
from __future__ import annotations

import numpy as np

from lagh.stochcheck import Consumer, Declaration, Task, component


def ou_paths(*, theta: float, b: float, T: float, dt: float, n_traj: int,
             x0: float = 1.0, seed: int = 0):
    """EXACT OU transition: X_{t+dt} = X_t e^{-theta dt} + b sqrt((1-e^{-2 theta
    dt})/(2 theta)) Z. No discretization bias of any order."""
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    t = np.arange(n) * dt
    dec = np.exp(-theta * dt)
    sd = b * np.sqrt((1.0 - dec ** 2) / (2.0 * theta))
    X = np.empty((n_traj, n))
    X[:, 0] = x0
    Z = rng.standard_normal((n_traj, n - 1))
    for k in range(1, n):
        X[:, k] = dec * X[:, k - 1] + sd * Z[:, k - 1]
    return t, X


def gbm_paths(*, mu: float, b: float, T: float, dt: float, n_traj: int,
              s0: float = 1.0, seed: int = 0):
    """EXACT GBM: S_t = s0 exp((mu - b^2/2) t + b W_t)."""
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    t = np.arange(n) * dt
    dW = rng.standard_normal((n_traj, n - 1)) * np.sqrt(dt)
    W = np.concatenate([np.zeros((n_traj, 1)), np.cumsum(dW, axis=1)], axis=1)
    return t, s0 * np.exp((mu - 0.5 * b ** 2) * t[None, :] + b * W)


def ode_obs_paths(*, theta: float, sigma_obs: float, T: float, dt: float,
                  n_traj: int, x0: float = 1.0, seed: int = 0):
    """Deterministic decay, observed with independent Gaussian error. The
    trajectories differ ONLY in their observation noise, which is the point: the
    process is the same one every time."""
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    t = np.arange(n) * dt
    truth = x0 * np.exp(-theta * t)
    return t, truth[None, :] + sigma_obs * rng.standard_normal((n_traj, n))


def brownian_paths(*, b: float, T: float, dt: float, n_traj: int,
                   seed: int = 0):
    """dX = b dW: the null with no drift at all."""
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    t = np.arange(n) * dt
    dW = rng.standard_normal((n_traj, n - 1)) * np.sqrt(dt) * b
    return t, np.concatenate([np.zeros((n_traj, 1)), np.cumsum(dW, axis=1)],
                             axis=1)


# The declared library for every Level 0 task. Registered here, before any run:
# the claim's vocabulary is these four terms and nothing else, so a task's `truth`
# lists all four INCLUDING the zeros (docs/STOCHASTIC_CHECKER.md §2).
LIBRARY = ("1", "x", "x**2", "x**3")


def _task(task_id, system, truth, expectation, sampling, *, level=0, null=False,
          null_reason="", declarations=(), invariants=()):
    return Task(task_id=task_id, level=level, system=system, state_dim=1,
                truth={component("drift", g): truth.get(g, 0.0)
                       for g in LIBRARY},
                expectation={component("drift", g): expectation.get(g, "interval")
                             for g in LIBRARY},
                sampling=sampling, declarations=tuple(declarations),
                invariants=tuple(invariants), null=null, null_reason=null_reason)


# ---------------------------------------------------------------- the L0 suite
# Registered expectations, per component, BEFORE the run. The reasoning, so it can
# be judged rather than taken on trust:
#
#   drift:x       the only true term in OU/GBM/ode_obs. An interval, not `exact`:
#                 a martingale band cannot pin a rational, and claiming it could
#                 is what RNOISE_STUDY.md's parametric gate exists to refuse.
#   drift:1       true coefficient 0. `interval` -- a band around zero is exactly
#                 what the data supports, and it is honest.
#   drift:x**2,
#   drift:x**3    true coefficient 0, and expected `interval` for the same reason.
#                 If they come back UNCONSTRAINED that is a reach loss to report,
#                 not a failure.
#   nulls         every component `abstain`: no law exists (pure_noise) or none is
#                 identifiable at the given dt (coarse_ou).

def level0_tasks(*, theta=1.0, b=1.4, mu=0.8, b_gbm=0.02, sigma_obs=1e-3,
                 T=640.0, dt=1e-3, n_traj=8, seed=0) -> list:
    """The defaults are MEASURED, not guessed (2026-07-29), and each one is a
    consequence of the identifiability arithmetic in `lagh/ito.py`:

    * `T`/`theta` satisfy theta*L > 2*kappa^2 at L = T/10, the threshold below
      which an OU drift is vacuous at ANY sampling rate. At the doc's original
      T = 40 every OU verdict was ABSTAIN[noise], correctly.
    * `b = 1.4` at theta = 1 puts the stationary spread at ~1.0 rather than 0.23.
      A narrow visited range does not separate {1, x, x^2, x^3} -- the joint bound
      on the x coefficient was [-24, +5.7] at spread 0.13 -- which is the
      "unexplored state space" null appearing on a SYSTEM task.
    * `b_gbm = 0.02` is where GBM's x term becomes RESOLVED (its joint bound
      excludes zero). The bound's width is proportional to b, so this is a
      declared operating point on a measured curve, not a special value.
    """
    lin, zero = "interval", "interval"
    return [
        _task("L0-ou", "ou", {"x": -theta},
              {"x": lin, "1": zero, "x**2": zero, "x**3": zero},
              {"dt": dt, "T": T, "n_traj": n_traj, "seed": seed,
               "sigma_obs": 0.0, "exact_transition": True},
              declarations=(Declaration(Consumer.DIFFUSION_QV, b ** 2,
                                        provenance="declared",
                                        note="b^2 is constant for OU, so the "
                                             "reference state is every state"),)),
        _task("L0-gbm", "gbm", {"x": mu},
              {"x": lin, "1": zero, "x**2": zero, "x**3": zero},
              {"dt": 1e-4, "T": 6.0, "n_traj": n_traj, "seed": seed + 1,
               "sigma_obs": 0.0, "exact_transition": True},
              declarations=(Declaration(Consumer.DIFFUSION_QV, b_gbm ** 2,
                                        provenance="declared",
                                        note="reference state x = 1: b^2 x^2 is "
                                             "b^2 there"),)),
        _task("L0-ode-obs", "ode_obs", {"x": -theta},
              {"x": lin, "1": zero, "x**2": zero, "x**3": zero},
              {"dt": dt, "T": T, "n_traj": n_traj, "seed": seed + 2,
               "sigma_obs": sigma_obs, "exact_transition": True},
              declarations=(Declaration(Consumer.OBSERVATION, sigma_obs,
                                        provenance="declared"),)),
        _task("L0-null-noise", "pure_noise", {},
              {g: "abstain" for g in LIBRARY},
              {"dt": dt, "T": T, "n_traj": n_traj, "seed": seed + 3,
               "sigma_obs": 0.0, "exact_transition": True},
              null=True, null_reason="noise",
              declarations=(Declaration(Consumer.DIFFUSION_QV, b ** 2,
                                        provenance="declared"),)),
        _task("L0-null-coarse", "coarse_ou", {},
              {g: "abstain" for g in LIBRARY},
              {"dt": 4.0 / theta, "T": T * 20, "n_traj": n_traj, "seed": seed + 4,
               "sigma_obs": 0.0, "exact_transition": True,
               "note": "sampled at 4/theta -- beyond the correlation time"},
              null=True, null_reason="resolution"),
    ]


def paths_for(task: Task, **over):
    """Simulate the trajectories a task describes. Returns (t, paths, truth_info)."""
    s = dict(task.sampling)
    s.update(over)
    T, dt, n_traj, seed = s["T"], s["dt"], s["n_traj"], s["seed"]
    if task.system == "ou":
        theta = -task.truth[component("drift", "x")]
        b = float(next(d.value for d in task.declarations
                       if d.consumer is Consumer.DIFFUSION_QV)) ** 0.5
        t, X = ou_paths(theta=theta, b=b, T=T, dt=dt, n_traj=n_traj, seed=seed)
        return t, X, {"b": b, "theta": theta}
    if task.system == "coarse_ou":
        # the null's truth is all zeros (no identifiable law at this dt), so its
        # generating parameters live here rather than in `truth`
        theta, b = 1.5, 0.4
        t, X = ou_paths(theta=theta, b=b, T=T, dt=dt, n_traj=n_traj, seed=seed)
        return t, X, {"b": b, "theta": theta,
                      "generated_from": "ou at dt = 4/theta"}
    if task.system == "gbm":
        mu = task.truth[component("drift", "x")]
        b = float(next(d.value for d in task.declarations
                       if d.consumer is Consumer.DIFFUSION_QV)) ** 0.5
        t, X = gbm_paths(mu=mu, b=b, T=T, dt=dt, n_traj=n_traj, seed=seed)
        return t, X, {"b": b, "mu": mu}
    if task.system == "ode_obs":
        theta = -task.truth[component("drift", "x")]
        t, X = ode_obs_paths(theta=theta, sigma_obs=s["sigma_obs"], T=T, dt=dt,
                             n_traj=n_traj, seed=seed)
        return t, X, {"theta": theta, "sigma_obs": s["sigma_obs"]}
    if task.system == "pure_noise":
        b = float(next(d.value for d in task.declarations
                       if d.consumer is Consumer.DIFFUSION_QV)) ** 0.5
        t, X = brownian_paths(b=b, T=T, dt=dt, n_traj=n_traj, seed=seed)
        return t, X, {"b": b}
    raise ValueError(f"unknown system {task.system!r}")
