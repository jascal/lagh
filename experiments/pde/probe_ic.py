"""SCOPING PROBE (not a campaign): can an initial condition be certified
mode-by-mode, and does the resolution curve behave as the physics says?

Skips the weak-form factory entirely. The law is KNOWN, the domain is periodic,
so the forward operator is diagonal in Fourier space: mode k of the initial
condition reaches the observation time multiplied by d_k (exp(-nu k^2 T) for
heat, a unit-modulus phase for advection). The question is what the DECLARED
band then permits for each mode -- measured with the production interval search
(`certify.parameter_interval`), not a bespoke one, so the numbers preview what a
real state certificate would say.

Prediction under test: half-width(k) ~ eps / d_k, i.e. exp(+nu k^2 T) * sigma
for heat (the exact inverse of the forward decay) and FLAT in k for advection.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.certify import epsilon, parameter_interval           # noqa: E402

NX, T = 256, 0.3
X = np.linspace(0.0, 2 * np.pi, NX, endpoint=False)
A = sp.Symbol("x_0")


def decay(kind, k, nu=0.1, c=0.7):
    if kind == "heat":
        return np.exp(-nu * k ** 2 * T)
    return 1.0                       # advection: unitary, |d_k| = 1


def field_at_T(kind, coeffs, nu=0.1, c=0.7):
    """u(x, T) for an initial condition given as cosine-mode amplitudes."""
    u = np.zeros(NX)
    for k, a in coeffs.items():
        if kind == "heat":
            u += a * decay("heat", k, nu) * np.cos(k * X)
        else:                        # translation by c*T
            u += a * np.cos(k * (X - c * T))
    return u


def mode_interval(kind, k, coeffs, sigma, nu=0.1, c=0.7):
    """The interval of amplitudes for mode k that still certifies against the
    observation at time T, everything else held -- the production predicate."""
    y = field_at_T(kind, coeffs, nu, c)
    if sigma > 0:
        y = y + np.random.default_rng(4).normal(0, sigma, NX)
    others = {kk: v for kk, v in coeffs.items() if kk != k}
    base = field_at_T(kind, others, nu, c)
    if kind == "heat":
        col = decay("heat", k, nu) * np.cos(k * X)
    else:
        col = np.cos(k * (X - c * T))
    resid = y - base
    # sigma here is ABSOLUTE field noise, so it must be declared with a unit
    # proportional term -- epsilon's default treats sigma as relative to |y|
    eps = epsilon(y, sigma=sigma, prop=np.ones_like(y), floor_abs=0.0)
    a0 = coeffs[k]
    iv = parameter_interval(sp.Float(a0) * A, [A], col.reshape(-1, 1), resid,
                            eps, sp.Float(a0), max_rel=1e6)
    return iv


def main():
    coeffs = {1: 1.0, 2: 0.6, 3: 0.4, 4: 0.3, 5: 0.25, 6: 0.2, 8: 0.15,
              10: 0.12, 12: 0.1}
    for kind in ("advection", "heat"):
        print(f"\n=== {kind}, T = {T}, N = {NX}")
        for sigma in (1e-6, 1e-4):
            print(f"  sigma = {sigma:g}")
            for k in sorted(coeffs):
                iv = mode_interval(kind, k, coeffs, sigma)
                d = decay(kind, k)
                if iv is None:
                    print(f"    k={k:3d} d_k={d:.3e}  UNDETERMINED (abstain)")
                    continue
                hw = 0.5 * (iv[1] - iv[0])
                print(f"    k={k:3d} d_k={d:.3e}  half-width {hw:.3e}  "
                      f"hw/sigma {hw/sigma:8.2f}  hw*d_k/sigma {hw*d/sigma:6.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
