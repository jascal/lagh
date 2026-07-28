"""Exact analytic PDE solutions for the weak-form dev campaign C0
(docs/CASE_STUDY_PDE_DEV.md).

No solver, no external data: every field here satisfies its PDE to machine
precision, so a failed certification is the instrument's verdict and never the
integrator's error. That is the whole point of running C0 before PDEBench.
"""
from __future__ import annotations

import numpy as np

L = 2 * np.pi


def grid(nx=257, nt=81, tmax=1.0, x0=0.0, xL=L):
    return np.linspace(x0, xL, nx), np.linspace(0.0, tmax, nt)


def heat(x, t, nu=0.1, modes=(1, 2), amps=(1.0, 0.4), phases=(0.0, 0.7)):
    """u_t = nu u_xx, exactly: a truncated Fourier series."""
    u = np.zeros((len(x), len(t)))
    for k, a, ph in zip(modes, amps, phases):
        u += a * np.exp(-nu * k ** 2 * t)[None, :] * np.sin(k * x[:, None] + ph)
    return u


def advection(x, t, c=0.7, modes=(1, 2), amps=(1.0, 0.3), phases=(0.0, 1.1)):
    """u_t = -c u_x, exactly: any profile translated."""
    z = x[:, None] - c * t[None, :]
    u = np.zeros_like(z)
    for k, a, ph in zip(modes, amps, phases):
        u += a * np.sin(k * z + ph)
    return u


def burgers_coleHopf(x, t, nu=0.2, modes=(1, 2), amps=(0.3, 0.15),
                     phases=(0.0, 0.5)):
    """u_t = nu u_xx - u u_x, exactly, via Cole-Hopf: u = -2 nu phi_x / phi with
    phi any positive solution of the heat equation. Gives as many DISTINCT
    initial conditions as we want, all exact -- which is what the multi-solution
    identifiability predictions need."""
    phi = np.ones((len(x), len(t)))
    dphi = np.zeros((len(x), len(t)))
    for k, a, ph in zip(modes, amps, phases):
        e = np.exp(-nu * k ** 2 * t)[None, :]
        phi += a * e * np.cos(k * x[:, None] + ph)
        dphi += -a * k * e * np.sin(k * x[:, None] + ph)
    assert phi.min() > 0.1, "Cole-Hopf potential must stay positive"
    return -2.0 * nu * dphi / phi


def burgers_wave(x, t, nu=0.2, a=1.0, c=0.5):
    """The traveling wave u = c - a tanh(a z / (2 nu)), z = x - c t: an exact
    Burgers solution on which u_t = -c u_x holds IDENTICALLY. The single-solution
    degeneracy case (N5)."""
    z = x[:, None] - c * t[None, :]
    return c - a * np.tanh(a * z / (2 * nu))


def kdv_soliton(x, t, k=0.5, x0=0.0):
    """u_t = -6 u u_x - u_xxx, exactly: 12k^2 sech^2(k(x - 4k^2 t)). Also a
    traveling wave -- the second N5 case."""
    z = x[:, None] - 4 * k ** 2 * t[None, :] - x0
    return 12 * k ** 2 / np.cosh(k * z) ** 2


def smooth_random(x, t, seed=0, kmax=3):
    """A smooth field solving no library PDE: band-limited noise in x, an
    unrelated smooth time dependence. The hard null (N1b) -- 'smooth' must not
    be enough to certify."""
    rng = np.random.default_rng(seed)
    u = np.zeros((len(x), len(t)))
    for k in range(1, kmax + 1):
        a, b = rng.normal(0, 1, 2)
        w = rng.uniform(0.5, 2.0)
        u += a * np.sin(k * x[:, None] + b) * (1.0 + 0.5 * np.sin(w * t[None, :]))
    return u


def ic_family(fn, n_ic, **kw):
    """n_ic distinct exact solutions of the same PDE (different initial data):
    the multi-solution evidence the single-trajectory finding says is required."""
    out = []
    rng = np.random.default_rng(11)
    for i in range(n_ic):
        modes = tuple(sorted(rng.choice([1, 2, 3], size=2, replace=False).tolist()))
        amps = tuple(np.round(rng.uniform(0.2, 1.0, 2), 3).tolist())
        phases = tuple(np.round(rng.uniform(0, 2 * np.pi, 2), 3).tolist())
        out.append(dict(modes=modes, amps=amps, phases=phases, **kw))
    return out
