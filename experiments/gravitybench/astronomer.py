"""The deterministic lagh astronomer (H2a): observation planning + quantity
playbook, no LLM anywhere. Developed against experiments/gravitybench/integrator
ONLY; the sealed benchmark is touched once, at the registered read.

Playbook entries all trace to fits; lagh certificates (with alpha) attach where
the corresponding law certifies. Best-fit values are ALWAYS submitted (numeric-
accuracy scoring; the two-track lesson).
"""
from __future__ import annotations

import numpy as np

G_SI = 6.67430e-11


# ------------------------------------------------------------------ planning

def plan_and_observe(observe, maxtime, budget=100, per_request=10):
    """Fixed policy (registered): coarse scan -> period estimate -> two dense
    phase-covering passes + three close-spaced triplets for accelerations.
    `observe(times)->dict of columns` is the only interface used."""
    n_coarse = min(40, budget // 2)
    t_coarse = np.linspace(0.0, maxtime, n_coarse)
    obs = _merge(None, observe(t_coarse))
    P = estimate_period(obs)
    remaining = budget - n_coarse
    if P is None or not np.isfinite(P):
        # unbound / under-covered: uniform fill PLUS strict-cadence acceleration
        # triplets across the span (the mass/force-law estimators need them --
        # measured battery failure on the flyby case)
        n_trip = 9
        t_extra = np.linspace(0.0, maxtime, remaining - n_trip + 2)[1:-1]
        obs = _merge(obs, observe(t_extra))
        dt = maxtime / 2000.0
        trips = []
        for phase in (0.2, 0.5, 0.8):
            tc = phase * maxtime
            trips.append((tc - dt, tc, tc + dt))
        flat = np.clip(np.asarray(trips).ravel(), 0, maxtime)
        obs = _merge(obs, observe(flat))
        obs["triplets"] = np.asarray(trips)
        return obs, None
    # dense pass over ~1.5 periods + acceleration triplets at three phases
    n_trip = 9
    n_dense = remaining - n_trip
    t0 = 0.25 * maxtime
    t_dense = t0 + np.linspace(0.0, min(1.5 * P, 0.6 * maxtime), max(n_dense, 4))
    obs = _merge(obs, observe(np.clip(t_dense, 0, maxtime)))
    dt = P / 400.0
    trips = []
    for phase in (0.15, 0.45, 0.8):
        tc = t0 + phase * P
        trips.append((tc - dt, tc, tc + dt))
    flat = np.clip(np.asarray(trips).ravel(), 0, maxtime)
    obs = _merge(obs, observe(flat))
    obs["triplets"] = np.asarray(trips)
    return obs, P


def _merge(a, b):
    """Merge observation frames; the 'triplets' metadata key is carried through
    untouched (it is planner bookkeeping, not a data column)."""
    meta = a.pop("triplets", None) if isinstance(a, dict) else None
    b = {k: np.asarray(v, float) for k, v in b.items() if k != "triplets"}
    if a is None:
        out = b
    else:
        out = {k: np.concatenate([a[k], b[k]]) for k in a}
    order = np.argsort(out["time"])
    dedup = {k: v[order] for k, v in out.items()}
    _, uniq = np.unique(dedup["time"], return_index=True)
    out = {k: v[uniq] for k, v in dedup.items()}
    if meta is not None:
        out["triplets"] = meta
    return out


# ------------------------------------------------------------------ features

def separation(obs):
    d = np.stack([obs[f"star2_{ax}"] - obs[f"star1_{ax}"] for ax in "xyz"], 1)
    return np.sqrt((d ** 2).sum(1))


def estimate_period(obs):
    """Angle-based primary estimator: unwrap the in-plane angle of the relative
    vector; P = 2*pi * span / total swept angle. Works for ANY bound orbit
    including near-circular (where separation barely modulates and phase-
    dispersion is blind -- measured battery failure). Returns None when less
    than ~0.7 revolutions swept (unbound / too-short window)."""
    t = obs["time"]
    if len(t) < 8:
        return None
    d = np.stack([obs[f"star2_{ax}"] - obs[f"star1_{ax}"] for ax in "xyz"], 1)
    # orbital plane via the two principal axes of the relative positions
    dc = d - d.mean(0)
    _, _, Vt = np.linalg.svd(dc, full_matrices=False)
    u, v = dc @ Vt[0], dc @ Vt[1]
    theta = np.unwrap(np.arctan2(v, u))
    swept = abs(theta[-1] - theta[0])
    if swept < 0.7 * 2 * np.pi:
        return None
    span = t[-1] - t[0]
    return float(2 * np.pi * span / swept)


def refine_period(obs, P0):
    t, r = obs["time"], separation(obs)
    r = (r - r.mean()) / (r.std() + 1e-300)
    best = (np.inf, P0)
    for P in np.linspace(0.9 * P0, 1.1 * P0, 4001):
        ph = np.mod(t, P) / P
        o = np.argsort(ph)
        rough = float(np.sum(np.diff(r[o]) ** 2))
        if rough < best[0]:
            best = (rough, P)
    return best[1]


# ------------------------------------------------------------------ playbook

def task_period(obs):
    P = estimate_period(obs)
    return refine_period(obs, P) if P else None


def task_mass_ratio(obs):
    """m1/m2 from amplitudes about the (drift-corrected) center of motion:
    m1*A1 = m2*A2. Robust to proper motion via linear detrend per axis."""
    amps = {}
    t = obs["time"]
    A = np.ones((len(t), 2)); A[:, 1] = t
    for s in ("star1", "star2"):
        amp2 = 0.0
        for ax in "xyz":
            v = obs[f"{s}_{ax}"]
            c, *_ = np.linalg.lstsq(A, v, rcond=None)
            res = v - A @ c
            amp2 += float(res.var())
        amps[s] = np.sqrt(amp2)
    if amps["star1"] <= 0:
        return None
    return amps["star2"] / amps["star1"]        # m1/m2 = A2/A1


def task_total_mass(obs, P):
    """Kepler III on the relative orbit: a from (r_max+r_min)/2 over a full
    period of dense coverage; M = 4 pi^2 a^3 / (G P^2). Newtonian tasks only."""
    r = separation(obs)
    a = 0.5 * (np.max(r) + np.min(r))
    return 4 * np.pi ** 2 * a ** 3 / (G_SI * P ** 2)


def task_masses(obs, P):
    q = task_mass_ratio(obs)          # m1/m2
    M = task_total_mass(obs, P)
    if q is None or M is None:
        return None, None
    m2 = M / (1 + q)
    m1 = M - m2
    return m1, m2


def _accels(obs):
    """Second differences at the planner's EXPLICIT triplets (metadata carried on
    the obs dict). Inferring cadence from the merged series was twice-measured
    fragile: dense-pass rows slipped under absolute caps (0.59 mass error), and
    accidental near-coincidences between passes defined a bogus finest-cadence
    under self-scaling (near-circular short-window nan). When no metadata exists
    (the full-table variant), fall back to uniform-cadence scanning."""
    t = obs["time"]
    d = np.stack([obs[f"star2_{ax}"] - obs[f"star1_{ax}"] for ax in "xyz"], 1)
    rs, accs = [], []
    trips = obs.get("triplets")
    if trips is not None:
        for ta, tb, tc in np.asarray(trips):
            ia = int(np.argmin(np.abs(t - ta)))
            ib = int(np.argmin(np.abs(t - tb)))
            ic = int(np.argmin(np.abs(t - tc)))
            dt1, dt2 = t[ib] - t[ia], t[ic] - t[ib]
            if dt1 <= 0 or dt2 <= 0 or abs(dt1 - dt2) > 0.05 * max(dt1, dt2):
                continue                 # clipped/collided triplet: skip
            acc = (d[ic] - 2 * d[ib] + d[ia]) / (dt1 * dt2)
            rs.append(float(np.sqrt(d[ib] @ d[ib])))
            accs.append(float(np.sqrt(acc @ acc)))
        return np.asarray(rs), np.asarray(accs)
    # full-table fallback: uniform cadence -- every interior point qualifies
    for i in range(1, len(t) - 1):
        dt1, dt2 = t[i] - t[i - 1], t[i + 1] - t[i]
        if dt1 <= 0 or dt2 <= 0 or abs(dt1 - dt2) > 0.02 * max(dt1, dt2):
            continue
        acc = (d[i + 1] - 2 * d[i] + d[i - 1]) / (dt1 * dt2)
        rs.append(float(np.sqrt(d[i] @ d[i])))
        accs.append(float(np.sqrt(acc @ acc)))
    return np.asarray(rs), np.asarray(accs)


def task_gravity_exponent(obs):
    """F ~ r^p: lagh C3 log-log recovery on |acc_rel| vs r -- certificate+alpha
    when it certifies; the raw slope as best-fit either way."""
    rs, accs = _accels(obs)
    if len(rs) < 3:
        return None, None
    L = np.column_stack([np.ones(len(rs)), np.log(rs)])
    c, *_ = np.linalg.lstsq(L, np.log(accs), rcond=None)
    slope = float(c[1])
    cert = None
    try:
        from lagh.mcp.core import recover
        out = recover(rs.reshape(-1, 1).tolist(), accs.tolist(), sigma=1e-3)
        if out.get("certified"):
            cert = out
    except Exception:                                          # noqa: BLE001
        pass
    return slope, cert


def task_drag_tau(obs):
    """Linear drag: the relative-orbit amplitude decays ~ exp(-t/tau) (momentum
    damping); tau from the log-envelope slope of |r - r_trend|."""
    t = obs["time"]
    r = separation(obs)
    A = np.ones((len(t), 2)); A[:, 1] = t
    c, *_ = np.linalg.lstsq(A, r, rcond=None)
    env = np.abs(r - A @ c)
    # PEAK sequence only: the raw envelope includes near-zero crossings whose
    # log-scatter swamps the slope's standard error, so even a true 34% decay
    # failed the significance gate (measured). Local maxima decay cleanly as
    # e^{-2t/tau}, one per half-orbit.
    peaks = [i for i in range(1, len(env) - 1)
             if env[i] >= env[i - 1] and env[i] >= env[i + 1]
             and env[i] > 0.05 * np.max(env)]
    if len(peaks) < 4:
        return None
    t = t[peaks]
    env = env[peaks]
    m = np.ones(len(t), bool)
    L = np.column_stack([np.ones(int(m.sum())), t[m]])
    yfit = np.log(env[m])
    ce, *_ = np.linalg.lstsq(L, yfit, rcond=None)
    lam = float(ce[1])
    # virial energy balance: dE/dt = -2K/tau with <K> = -E gives |E| ~ e^{+2t/tau},
    # a ~ 1/|E| ~ e^{-2t/tau} -- the envelope decays at rate 2/tau (measured
    # factor-2 bias before this)
    if lam >= 0:
        return None
    # SIGNIFICANCE gate (supersedes two magnitude windows, both measured wrong:
    # 3*span rejected true 34%-decay drag; 10*span readmitted flat-envelope noise
    # slopes at short windows, which silently amputated the twin's backward
    # branch). A claimed decay must be statistically significant against the
    # envelope fit's own scatter: |lam| >= 3 * stderr(lam).
    resid = yfit - L @ ce
    n = len(yfit)
    tbar = L[:, 1].mean()
    denom = float(np.sum((L[:, 1] - tbar) ** 2))
    if n <= 2 or denom <= 0:
        return None
    se_lam = float(np.sqrt(np.sum(resid ** 2) / (n - 2) / denom))
    if abs(lam) < 3.0 * se_lam:
        return None
    return -2.0 / lam
