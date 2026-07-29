"""SCOPING PROBE (not a campaign): does a COUPLED system of PDEs certify
equation-by-equation over shared patch rows?

Uses an exactly-solvable linear system, so any failure is the instrument's and
never an integrator's:

    u_t = a u_xx + b v
    v_t = c v_xx + d u

which diagonalizes per Fourier mode (2x2 matrix exponential), giving exact
fields for as many initial conditions as we want.

The probe deliberately needs NO factory extension: for a linear system every
feature is a single-field term of u or of v, so the existing single-field
`build()` run on each field, on the SAME patches, produces every column. What it
tests is the part that is genuinely new -- two targets sharing one row set, each
equation's features spanning both fields, and whether the coupling costs
identifiability.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.engine import discover                              # noqa: E402
from lagh.weakform import (PatchEpsilon, build,               # noqa: E402
                           multiscale_patches)

A, B, C, D = 0.1, 0.5, 0.05, -0.3          # the true constants
NX, NT, TMAX = 257, 81, 1.0
SCALES = [(24, 12), (32, 16), (40, 20)]
TERMS = ["u_t", "u_xx", "u_x", "u", "1"]


def exact_pair(seed):
    """Exact (u, v) for one initial condition, via the per-mode 2x2 propagator."""
    x = np.linspace(0.0, 2 * np.pi, NX)
    t = np.linspace(0.0, TMAX, NT)
    rng = np.random.default_rng(seed)
    modes = [1, 2, 3]
    amps = rng.uniform(0.2, 1.0, (2, len(modes)))
    # ONE phase per mode, shared by both fields: the 2x2 propagator acts on the
    # amplitude pair of the SAME basis function, so independent phases would
    # make the coupling term B*v spatially orthogonal to u and the system would
    # not hold at all (measured: the truth then misses its own band by 0.85)
    ph = rng.uniform(0, 2 * np.pi, len(modes))
    u = np.zeros((NX, NT))
    v = np.zeros((NX, NT))
    for j, k in enumerate(modes):
        M = np.array([[-A * k ** 2, B], [D, -C * k ** 2]])
        for i, ti in enumerate(t):
            e = expm(M * ti)
            a0 = np.array([amps[0, j], amps[1, j]])
            at = e @ a0
            u[:, i] += at[0] * np.cos(k * x + ph[j])
            v[:, i] += at[1] * np.cos(k * x + ph[j])
    return u, v, x, t


def system_rows(sigma, n_ic=4, seed0=0):
    """Weak-form rows for both fields on the SAME patches, normalized."""
    rng = np.random.default_rng(99)
    cols, dets, grams, sol = [], [], [], []
    names = None
    for i in range(n_ic):
        u, v, x, t = exact_pair(seed0 + i)
        if sigma > 0:
            u = u + rng.normal(0, sigma, u.shape)
            v = v + rng.normal(0, sigma, v.shape)
        pa = multiscale_patches(x, t, SCALES, n_x=4, n_t=3)
        su = build(u, x, t, TERMS, pa, p=16, sigma=sigma).normalize(by="1")
        sv = build(v, x, t, TERMS, pa, p=16, sigma=sigma).normalize(by="1")
        if len(su.A) == 0 or len(su.A) != len(sv.A):
            continue
        names = ([f"u:{n}" for n in su.names] + [f"v:{n}" for n in sv.names])
        cols.append(np.hstack([su.A, sv.A]))
        dets.append(np.hstack([su.quad + su.roundoff, sv.quad + sv.roundoff]))
        g = np.zeros((len(su.A), len(names), len(names)))
        nu_ = su.gram.shape[1]
        g[:, :nu_, :nu_] = su.gram          # independent field noise -> block
        g[:, nu_:, nu_:] = sv.gram          # diagonal, no cross terms
        grams.append(g)
        sol.append(np.full(len(su.A), i))
    return (np.vstack(cols), np.vstack(dets), np.concatenate(grams),
            np.concatenate(sol), names)


def run_equation(target, sigma, A_all, det, gram, sol, names):
    j = names.index(target)
    feat = [n for n in names if n not in (target, "u:u_t", "v:u_t")]
    cols = [names.index(n) for n in feat]
    m = PatchEpsilon(names, target, A_all[:, j], A_all[:, cols], det, gram,
                     sigma=sigma, floor_abs=0.0, coeff_max=2.0)
    held = np.unique(sol)[-1]
    tr = np.random.default_rng(0).permutation(np.where(sol != held)[0])
    ce = np.where(sol == held)[0]
    a = int(0.75 * len(tr))
    X, y = A_all[:, cols], A_all[:, j]
    r = discover(X[tr[:a]], y[tr[:a]], X[tr[a:]], y[tr[a:]], X[ce], y[ce],
                 sigma=sigma, eps_model=m.subset(ce), max_tier=3,
                 declared_basis=True, band_sel=m.subset(tr[a:])(None))
    law = str(r.expr) if r.certificate.certified else ""
    for i, nm in enumerate(feat):
        law = law.replace(f"x_{i}", f"[{nm}]")
    return r.certificate.certified, r.certificate.abstain, law


def main():
    for sigma in (0.0, 1e-6, 1e-4):
        A_all, det, gram, sol, names = system_rows(sigma)
        print(f"\nsigma = {sigma:g}   rows {len(A_all)}  features {len(names)}")
        for tgt, want in (("u:u_t", "0.1*[u:u_xx] + 0.5*[v:u]"),
                          ("v:u_t", "0.05*[v:u_xx] - 0.3*[u:u]")):
            ok, ab, law = run_equation(tgt, sigma, A_all, det, gram, sol, names)
            print(f"  {tgt:6s} -> {'CERT   ' if ok else 'ABSTAIN'} "
                  f"{str(ab or ''):12s} {law[:70]}")
            print(f"         truth: {want}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
