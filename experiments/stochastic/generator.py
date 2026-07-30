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


def sde_path(drift, vol, x0, *, T: float, dt: float, seed: int = 0,
             substeps: int = 8):
    """Substepped Euler-Maruyama, OBSERVED on the `dt` grid.

    For a system with no exact transition. The substeps are not a nicety: the
    scheme's own weak error is O(dt_sim), and a run must not confuse the
    GENERATOR's bias with the instrument's band. Simulating at dt/substeps and
    observing on the dt grid puts the scheme's error `substeps` times below the
    quadrature bound the instrument declares for that grid.
    """
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    ds = dt / max(int(substeps), 1)
    sq = np.sqrt(ds)
    x = np.asarray(x0, float).copy()
    out = np.empty((len(x), n))
    out[:, 0] = x
    for k in range(1, n):
        for _ in range(substeps):
            x = x + drift(x) * ds + vol(x) * sq * rng.standard_normal(len(x))
        out[:, k] = x
    return np.arange(n) * dt, out


def double_well_paths(*, theta: float = 1.0, b: float = 1.4, T: float = 400.0,
                      dt: float = 1e-3, n_traj: int = 8, seed: int = 0,
                      multiplicative: bool = False, substeps: int = 8):
    """dX = theta(X - X^3) dt + b dW, or b X dW when `multiplicative`.

    The multiplicative variant is kept because of what it does rather than what it
    is: b(x) = b x vanishes at the origin, so the origin is unreachable and every
    trajectory is TRAPPED in one well (measured: zero well-crossings over T = 400
    across 8 trajectories, against ~21000 for the additive case). The registered
    third null -- trajectories that do not visit enough of state space -- arriving
    as a property of a system.
    """
    x0 = np.array([-1.0, 1.0] * (n_traj // 2) + [1.0] * (n_traj % 2))
    return sde_path(lambda x: theta * (x - x ** 3),
                    (lambda x: b * x) if multiplicative
                    else (lambda x: b * np.ones_like(x)),
                    x0, T=T, dt=dt, seed=seed, substeps=substeps)


def cir_paths(*, theta: float, m: float, b: float, T: float, dt: float,
              n_traj: int, seed: int = 0, x0: float | None = None):
    """EXACT CIR transition: dX = theta(m - X)dt + b sqrt(X) dW, so b^2(x) = b^2 x.

    The Level 1 system with a genuinely STATE-DEPENDENT diffusion, and exact rather
    than schemed so a verdict is the instrument's. The transition is a scaled
    noncentral chi-square:

        X_{t+dt} = ncx2(df = 4 theta m / b^2, nc = c X_t e^{-theta dt}) / c,
        c = 4 theta / (b^2 (1 - e^{-theta dt}))

    The 4 in `c` is load-bearing and was wrong once: with a 2 there the stationary
    mean comes out at 2m instead of m. `tests/test_ito.py` checks both stationary
    moments (mean m, variance m b^2 / 2 theta) for exactly that reason -- an exact
    sampler with a wrong constant is a silent generator bug, and the instrument
    would have been blamed for it.

    The Feller condition 2 theta m >= b^2 keeps X strictly positive; above it the
    boundary is attainable and the sampler still works, but sqrt(X) can reach 0.
    """
    from scipy.stats import ncx2
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    c = 4.0 * theta / (b ** 2 * (1.0 - np.exp(-theta * dt)))
    d = 4.0 * theta * m / b ** 2
    X = np.empty((n_traj, n))
    X[:, 0] = m if x0 is None else x0
    dec = np.exp(-theta * dt)
    for k in range(1, n):
        X[:, k] = ncx2.rvs(df=d, nc=c * X[:, k - 1] * dec, size=n_traj,
                           random_state=rng) / c
    return np.arange(n) * dt, X


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


# ------------------------------------------------------------- the L1 suite
# The diffusion library, registered alongside the drift library above. Both are
# part of every Level 1 task's declared vocabulary, so a task's truth lists every
# term of both -- zeros included.
DIFF_LIBRARY = ("1", "x", "x**2")


def _l1_task(task_id, system, drift, diff, expect_drift, expect_diff, sampling,
             *, null=False, null_reason="", declarations=()):
    truth = {component("drift", g): drift.get(g, 0.0) for g in LIBRARY}
    truth.update({component("diffusion", h): diff.get(h, 0.0)
                  for h in DIFF_LIBRARY})
    exp = {component("drift", g): expect_drift.get(g, "interval")
           for g in LIBRARY}
    exp.update({component("diffusion", h): expect_diff.get(h, "interval")
                for h in DIFF_LIBRARY})
    return Task(task_id=task_id, level=1, system=system, state_dim=1,
                truth=truth, expectation=exp, sampling=sampling,
                declarations=tuple(declarations), null=null,
                null_reason=null_reason)


def level1_tasks(*, seed: int = 0) -> list:
    """Four systems, drift AND diffusion scored per component.

    **The expectations here are CALIBRATED, not blind.** Level 0 and the two Level 1
    increments probed these systems first, so the registered expectation per
    component reflects what the instrument was measured to manage rather than a
    prediction made in ignorance. That is stated because it changes what the scored
    table is evidence FOR: it tests zero-confident-wrong (S5) and exercises the
    frozen checker on the diffusion path, and it is NOT a blind test of the
    expectations. A genuinely blind Level 1 would need systems nobody has probed.
    """
    return [
        _l1_task("L1-dw-additive", "double_well",
                 {"x": 1.0, "x**3": -1.0}, {"1": 1.4 ** 2},
                 # the wells sit where the drift vanishes, so the drift is expected
                 # to stay an interval; the constant diffusion should resolve
                 {"x": "interval", "x**3": "interval"},
                 {"1": "interval", "x": "interval", "x**2": "interval"},
                 {"dt": 1e-3, "T": 400.0, "n_traj": 8, "seed": seed,
                  "sigma_obs": 0.0, "substeps": 8, "b": 1.4, "theta": 1.0}),
        _l1_task("L1-dw-multiplicative", "double_well_mult",
                 {"x": 1.0, "x**3": -1.0}, {"x**2": 0.7 ** 2},
                 # b(x) = b x vanishes at the origin: the state space DISCONNECTS
                 # and each trajectory sees one well only
                 {g: "abstain" for g in LIBRARY},
                 {"1": "interval", "x": "interval", "x**2": "interval"},
                 {"dt": 1e-3, "T": 400.0, "n_traj": 8, "seed": seed + 1,
                  "sigma_obs": 0.0, "substeps": 8, "b": 0.7, "theta": 1.0}),
        _l1_task("L1-cir", "cir",
                 {"1": 1.0, "x": -1.0}, {"x": 1.0 ** 2},
                 {"1": "interval", "x": "interval"},
                 {"1": "interval", "x": "interval", "x**2": "interval"},
                 {"dt": 1e-3, "T": 200.0, "n_traj": 6, "seed": seed + 2,
                  "sigma_obs": 0.0, "b": 1.0, "theta": 1.0, "m": 1.0}),
        _l1_task("L1-gbm-mult", "gbm",
                 {"x": 0.8}, {"x**2": 0.2 ** 2},
                 {"x": "interval"},
                 {"1": "interval", "x": "interval", "x**2": "interval"},
                 {"dt": 1e-4, "T": 6.0, "n_traj": 8, "seed": seed + 3,
                  "sigma_obs": 0.0, "b": 0.2, "mu": 0.8}),
    ]


def l1_paths(task: Task):
    """Simulate the trajectories a Level 1 task describes."""
    s = task.sampling
    if task.system == "double_well":
        return double_well_paths(theta=s["theta"], b=s["b"], T=s["T"],
                                 dt=s["dt"], n_traj=s["n_traj"], seed=s["seed"],
                                 substeps=s["substeps"])
    if task.system == "double_well_mult":
        return double_well_paths(theta=s["theta"], b=s["b"], T=s["T"],
                                 dt=s["dt"], n_traj=s["n_traj"], seed=s["seed"],
                                 multiplicative=True, substeps=s["substeps"])
    if task.system == "cir":
        return cir_paths(theta=s["theta"], m=s["m"], b=s["b"], T=s["T"],
                         dt=s["dt"], n_traj=s["n_traj"], seed=s["seed"])
    if task.system == "gbm":
        return gbm_paths(mu=s["mu"], b=s["b"], T=s["T"], dt=s["dt"],
                         n_traj=s["n_traj"], seed=s["seed"])
    raise ValueError(f"unknown Level 1 system {task.system!r}")


def vdp_paths(*, mu: float = 1.0, b: float = 0.5, T: float = 200.0,
              dt: float = 1e-3, n_traj: int = 4, seed: int = 0,
              substeps: int = 8):
    """Van der Pol with additive noise on y ONLY:
        dx = y dt,   dy = (mu(1 - x^2) y - x) dt + b dW.

    The x component carries NO noise, which is why this system is here rather than
    only for being 2-D: its x equation is a deterministic weak-form identity, and
    measuring how much more determinable that makes it is the point.
    """
    rng = np.random.default_rng(seed)
    n = int(round(T / dt)) + 1
    ds = dt / max(int(substeps), 1)
    sq = np.sqrt(ds)
    x = rng.uniform(-2, 2, n_traj)
    y = rng.uniform(-2, 2, n_traj)
    X = np.empty((n_traj, n))
    Y = np.empty((n_traj, n))
    X[:, 0], Y[:, 0] = x, y
    for k in range(1, n):
        for _ in range(substeps):
            xn = x + y * ds
            y = y + (mu * (1 - x ** 2) * y - x) * ds \
                + b * sq * rng.standard_normal(n_traj)
            x = xn
        X[:, k], Y[:, k] = x, y
    return np.arange(n) * dt, X, Y


# Van der Pol's per-component drift libraries, registered here with the systems.
VDP_LIBRARIES = {"x": ("1", "x", "y", "x**2", "x*y", "y**2"),
                 "y": ("1", "x", "y", "x**2", "x*y", "x**2*y")}
VDP_FIELDS = ("x", "y")


def vdp_task(*, mu: float = 1.0, b: float = 0.5, seed: int = 10) -> Task:
    """Van der Pol as a scored task, using the frozen `part[index]:term` convention.

    Component INDEX is what makes a per-component drift expressible: `drift[0]:y` is
    the y term of the x component's drift. The convention was in the interface from
    the freeze and this is the first producer to need it.

    Expectations are CALIBRATED (measured per equation first), like the rest of
    Level 1: the noise-free x equation certifies, the driven y equation is vacuous at
    this configuration.
    """
    truth, exp = {}, {}
    for i, fld in enumerate(VDP_FIELDS):
        for g in VDP_LIBRARIES[fld]:
            truth[component("drift", g, i)] = 0.0
            # the driven component is measured vacuous here; the noise-free one is
            # measured determinable
            exp[component("drift", g, i)] = "interval" if i == 0 else "abstain"
    truth[component("drift", "y", 0)] = 1.0                     # dx = y dt
    truth[component("drift", "y", 1)] = mu                      # dy = mu y ...
    truth[component("drift", "x**2*y", 1)] = -mu                #      - mu x^2 y ...
    truth[component("drift", "x", 1)] = -1.0                    #      - x
    for i in range(2):
        for h in ("1", "x", "x**2"):
            truth[component("diffusion", h, i)] = 0.0
            exp[component("diffusion", h, i)] = "interval"
    truth[component("diffusion", "1", 1)] = b ** 2              # (bb^T)_yy = b^2
    return Task(task_id="L1-vdp", level=1, system="vdp", state_dim=2,
                truth=truth, expectation=exp,
                sampling={"dt": 1e-3, "T": 200.0, "n_traj": 4, "seed": seed,
                          "sigma_obs": 0.0, "substeps": 8, "mu": mu, "b": b},
                declarations=(Declaration(Consumer.DIFFUSION_QV, b ** 2,
                                          provenance="declared",
                                          note="(bb^T)_yy is the constant b^2; the "
                                               "x component carries no noise"),))
