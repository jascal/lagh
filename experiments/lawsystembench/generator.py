"""LawSystemBench v1 generator (docs/LAWSYSTEMBENCH.md): five seeded families,
8 draws each, clean + sigma_rep tiers -> 80 problems with exact ground truth.

Observables include rate columns (SINDy convention); columns are shuffled and
neutrally renamed c0..cK so roles are genuinely unlabeled. Ground truth stores
the equations (in original names AND the neutral mapping), shared constants,
and invariants.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

OUT = Path("experiments/results/lawsystembench_v1.jsonl")
N_DRAWS = 8
N_SAMPLES = 400
SIGMA_REP = 1e-6          # the noisy tier: float32-ish representation noise


def _traj(deriv, x0, t_end, n_steps=40000):
    """RK4 on dx/dt = deriv(x); returns dense states + derivatives."""
    x = np.asarray(x0, float)
    dt = t_end / n_steps
    X, D = [x.copy()], [np.asarray(deriv(x))]
    for _ in range(n_steps):
        k1 = np.asarray(deriv(x))
        k2 = np.asarray(deriv(x + 0.5 * dt * k1))
        k3 = np.asarray(deriv(x + 0.5 * dt * k2))
        k4 = np.asarray(deriv(x + dt * k3))
        x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        X.append(x.copy()); D.append(np.asarray(deriv(x)))
    return np.asarray(X), np.asarray(D)


def make_family(fam, rng):
    if fam == "F1_chain":
        k1, k2 = rng.uniform(0.3, 2.0), rng.uniform(0.3, 2.0)
        names = ["A", "B", "C"]
        deriv = lambda s: np.array([-k1 * s[0], k1 * s[0] - k2 * s[1], k2 * s[1]])
        x0 = [rng.uniform(1, 5), rng.uniform(0.1, 1), rng.uniform(0.0, 0.5) + 0.1]
        gt_eqs = {"dA_dt": f"-{k1}*A", "dB_dt": f"{k1}*A - {k2}*B",
                  "dC_dt": f"{k2}*B"}
        inv = ["A + B + C"]
        shared = {"k1": ["dA_dt", "dB_dt"], "k2": ["dB_dt", "dC_dt"]}
        t_end = 3.0 / min(k1, k2)
    elif fam == "F2_sir":
        beta, gamma = rng.uniform(0.2, 1.0), rng.uniform(0.05, 0.5)
        names = ["S", "I", "R"]
        deriv = lambda s: np.array([-beta * s[0] * s[1],
                                    beta * s[0] * s[1] - gamma * s[1],
                                    gamma * s[1]])
        x0 = [rng.uniform(0.7, 0.95), rng.uniform(0.01, 0.1), 0.05]
        gt_eqs = {"dS_dt": f"-{beta}*S*I", "dI_dt": f"{beta}*S*I - {gamma}*I",
                  "dR_dt": f"{gamma}*I"}
        inv = ["S + I + R"]
        shared = {"beta": ["dS_dt", "dI_dt"], "gamma": ["dI_dt", "dR_dt"]}
        t_end = 12.0 / gamma
    elif fam == "F3_lv":
        a, b, c, d = (rng.uniform(0.4, 1.4), rng.uniform(0.1, 0.6),
                      rng.uniform(0.4, 1.4), rng.uniform(0.1, 0.6))
        names = ["x", "y"]
        deriv = lambda s: np.array([a * s[0] - b * s[0] * s[1],
                                    -c * s[1] + d * s[0] * s[1]])
        x0 = [rng.uniform(1, 3), rng.uniform(1, 3)]
        gt_eqs = {"dx_dt": f"{a}*x - {b}*x*y", "dy_dt": f"-{c}*y + {d}*x*y"}
        inv = [f"{c}*log(x) - {d}*x + {a}*log(y) - {b}*y"]
        shared = {"b": ["dx_dt"], "d": ["dy_dt"]}
        t_end = 20.0 / a
    elif fam == "F4_osc":
        m1, m2 = rng.uniform(0.5, 2.0), rng.uniform(0.5, 2.0)
        k1, k2, kc = (rng.uniform(0.5, 3.0), rng.uniform(0.5, 3.0),
                      rng.uniform(0.2, 1.5))
        names = ["x1", "v1", "x2", "v2"]
        deriv = lambda s: np.array([
            s[1], (-(k1 + kc) * s[0] + kc * s[2]) / m1,
            s[3], (-(k2 + kc) * s[2] + kc * s[0]) / m2])
        x0 = [rng.uniform(-1, 1), 0.0, rng.uniform(-1, 1), 0.0]
        gt_eqs = {"dx1_dt": "v1",
                  "dv1_dt": f"(-({k1 + kc})*x1 + {kc}*x2)/{m1}",
                  "dx2_dt": "v2",
                  "dv2_dt": f"(-({k2 + kc})*x2 + {kc}*x1)/{m2}"}
        inv = [f"{m1}*v1**2/2 + {m2}*v2**2/2 + {k1}*x1**2/2 + {k2}*x2**2/2 "
               f"+ {kc}*(x1 - x2)**2/2"]
        shared = {"kc": ["dv1_dt", "dv2_dt"]}
        t_end = 12.0
    elif fam == "F5_kepler":
        GM = rng.uniform(0.5, 3.0)
        names = ["x", "y", "vx", "vy"]

        def deriv(s):
            r3 = (s[0] ** 2 + s[1] ** 2) ** 1.5
            return np.array([s[2], s[3], -GM * s[0] / r3, -GM * s[1] / r3])
        r0 = rng.uniform(0.8, 1.6)
        v0 = np.sqrt(GM / r0) * rng.uniform(0.85, 1.1)
        x0 = [r0, 0.0, 0.0, v0]
        gt_eqs = {"dx_dt": "vx", "dy_dt": "vy",
                  "dvx_dt": f"-{GM}*x/(x**2 + y**2)**(3/2)",
                  "dvy_dt": f"-{GM}*y/(x**2 + y**2)**(3/2)"}
        inv = ["x*vy - y*vx",
               f"(vx**2 + vy**2)/2 - {GM}/sqrt(x**2 + y**2)"]
        shared = {"GM": ["dvx_dt", "dvy_dt"]}
        t_end = 3 * 2 * np.pi * np.sqrt(r0 ** 3 / GM)
    else:
        raise KeyError(fam)
    return names, deriv, x0, t_end, gt_eqs, inv, shared


def main():
    rows = []
    rng_master = np.random.default_rng(20260727)
    for fam in ("F1_chain", "F2_sir", "F3_lv", "F4_osc", "F5_kepler"):
        for draw in range(N_DRAWS):
            rng = np.random.default_rng(rng_master.integers(2 ** 32))
            names, deriv, x0, t_end, gt_eqs, inv, shared = make_family(fam, rng)
            X, D = _traj(deriv, x0, t_end)
            idx = np.sort(rng.choice(len(X), N_SAMPLES, replace=False))
            cols = {}
            for j, nm in enumerate(names):
                cols[nm] = X[idx, j]
                cols[f"d{nm}_dt"] = D[idx, j]
            for tier, sig in (("clean", 0.0), ("noisy", SIGMA_REP)):
                data = {k: (v * (1 + sig * rng.standard_normal(len(v)))
                            if sig else v).tolist() for k, v in cols.items()}
                # neutral renaming: shuffle and mask the column names
                orig = sorted(data)
                perm = rng.permutation(len(orig))
                mapping = {orig[p]: f"c{i}" for i, p in enumerate(perm)}
                rows.append({
                    "id": f"{fam}_d{draw}_{tier}", "family": fam, "tier": tier,
                    "sigma": sig,
                    "columns": {mapping[k]: data[k] for k in orig},
                    "mapping": mapping,           # gt-only: scorer uses it
                    "gt_equations": gt_eqs, "gt_invariants": inv,
                    "gt_shared": shared})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"{len(rows)} problems -> {OUT}")


if __name__ == "__main__":
    main()
