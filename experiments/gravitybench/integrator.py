"""Self-built two-body integrator for H2a DEVELOPMENT ONLY
(docs/BLIND_READ_REGISTRATION_GRAVITYBENCH.md section 2.3).

Velocity-Verlet on the two-body problem with the benchmark's physics variants:
optional linear drag (dv/dt -= v/tau) and modified gravity (F ~ r^p instead of
r^-2). Emits the same observable schema as the sealed benchmark: (time,
star1_x/y/z, star2_x/y/z). The astronomer is developed and validated against
THIS; GravityBench scenarios stay untouched until the one-shot read.
"""
from __future__ import annotations

import numpy as np

G_SI = 6.67430e-11


class TwoBody:
    def __init__(self, m1, m2, r1, r2, v1, v2, *, G=G_SI, drag_tau=None,
                 mod_gravity_exponent=None, softening=0.0):
        self.m1, self.m2 = float(m1), float(m2)
        self.r1 = np.asarray(r1, float).copy()
        self.r2 = np.asarray(r2, float).copy()
        self.v1 = np.asarray(v1, float).copy()
        self.v2 = np.asarray(v2, float).copy()
        self.G = G
        self.drag_tau = drag_tau
        self.p = -2.0 if mod_gravity_exponent is None else float(mod_gravity_exponent)
        self.soft = softening

    def _acc(self, r1, r2, v1, v2):
        d = r2 - r1
        r = np.sqrt(d @ d + self.soft ** 2)
        # F = G m1 m2 * r^p * rhat  (p = -2 is Newton)
        fmag = self.G * self.m1 * self.m2 * r ** self.p
        f = fmag * d / r
        a1 = f / self.m1
        a2 = -f / self.m2
        if self.drag_tau is not None:
            a1 = a1 - v1 / self.drag_tau
            a2 = a2 - v2 / self.drag_tau
        return a1, a2

    def run(self, maxtime, n_steps=200_000):
        """Integrate and return dense arrays (t, r1(t), r2(t)) for interpolation."""
        dt = float(maxtime) / n_steps
        t = 0.0
        r1, r2, v1, v2 = self.r1.copy(), self.r2.copy(), self.v1.copy(), self.v2.copy()
        a1, a2 = self._acc(r1, r2, v1, v2)
        T = np.empty(n_steps + 1)
        R1 = np.empty((n_steps + 1, 3))
        R2 = np.empty((n_steps + 1, 3))
        T[0], R1[0], R2[0] = t, r1, r2
        for i in range(1, n_steps + 1):
            v1h = v1 + 0.5 * dt * a1
            v2h = v2 + 0.5 * dt * a2
            r1 = r1 + dt * v1h
            r2 = r2 + dt * v2h
            a1n, a2n = self._acc(r1, r2, v1h, v2h)
            v1 = v1h + 0.5 * dt * a1n
            v2 = v2h + 0.5 * dt * a2n
            a1, a2 = a1n, a2n
            t += dt
            T[i], R1[i], R2[i] = t, r1, r2
        self._T, self._R1, self._R2 = T, R1, R2
        return T, R1, R2

    def observe(self, times):
        """The benchmark-shaped observation call: positions at requested times
        (linear interpolation on the dense grid -- adequate at n_steps 2e5)."""
        times = np.asarray(times, float)
        out = {"time": times}
        for name, R in (("star1", self._R1), ("star2", self._R2)):
            for j, ax in enumerate("xyz"):
                out[f"{name}_{ax}"] = np.interp(times, self._T, R[:, j])
        return out


def make_circularish(m1=2e30, m2=1e30, a=1.5e11, ecc=0.3, seed=0, **kw):
    """A bound eccentric binary, COM at rest: standard dev scenario."""
    rng = np.random.default_rng(seed)
    M = m1 + m2
    # apoapsis start: r = a(1+e), v = sqrt(G M (1-e)/(a(1+e))) tangential
    r = a * (1 + ecc)
    v = np.sqrt(G_SI * M * (1 - ecc) / (a * (1 + ecc)))
    # split about COM
    r1 = np.array([-m2 / M * r, 0, 0])
    r2 = np.array([m1 / M * r, 0, 0])
    v1 = np.array([0, -m2 / M * v, 0])
    v2 = np.array([0, m1 / M * v, 0])
    period = 2 * np.pi * np.sqrt(a ** 3 / (G_SI * M))
    tb = TwoBody(m1, m2, r1, r2, v1, v2, **kw)
    return tb, period
