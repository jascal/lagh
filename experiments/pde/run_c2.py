"""PDE dev campaign C2 -- the verify track (registered: docs/CASE_STUDY_PDE_C2.md).

Takes the certificates C1b actually produced (experiments/results/pde_c1.json),
integrates each certified law forward from an initial condition NO stage of the
pipeline has seen, and asks whether the field stays inside the declared band.
"""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import fields as F                      # noqa: E402
from experiments.pde.verify import verify                    # noqa: E402

IN = Path("experiments/results/pde_c1.json")
OUT = Path("experiments/results/pde_c2.json")
NX, TMAX, NT = 256, 0.4, 21
RTOL = 1e-9
WRONG_FACTOR = 10.0            # V2: a law this far outside its own interval

TERM_OF = {"heat": "u_xx", "advection": "u_x"}


def fresh_solution(name, sigma, seed=101):
    """A 5th initial condition: not fitted, not certified against, not the
    held-out solution either. On a PERIODIC grid, which the patch grid was not."""
    x = np.linspace(0.0, 2 * np.pi, NX, endpoint=False)
    t = np.linspace(0.0, TMAX, NT)
    kw = F.ic_family(None, 5)[4]
    if name == "heat":
        u = F.heat(x, t, nu=0.1, **kw)
    elif name == "advection":
        u = F.advection(x, t, c=0.7, **kw)
    else:
        amps = np.array(kw["amps"])
        kw["amps"] = tuple(0.5 * amps / amps.sum())
        u = F.burgers_coleHopf(x, t, nu=0.2, **kw)
    clean = u.copy()
    if sigma > 0:
        u = u + np.random.default_rng(seed).normal(0, sigma, u.shape)
    return u, clean, x, t


def intervals_of(entry):
    """The certificate's parameter intervals, or degenerate ones from the point
    coefficients when the rung reported an exact value."""
    coeffs = entry.get("coefficients") or {}
    note = entry.get("interval_note") or ""
    got = re.findall(r"\[([-\d.e+]+), ([-\d.e+]+)\]", note)
    if not got:
        return {k: (v, v) for k, v in coeffs.items()}
    # notes list intervals in the same order as the law's atoms, smallest |v|
    # first; match each to the coefficient it brackets
    out = {}
    for lo, hi in got:
        lo, hi = float(lo), float(hi)
        for k, v in coeffs.items():
            if lo <= v <= hi and k not in out:
                out[k] = (lo, hi)
                break
    for k, v in coeffs.items():
        out.setdefault(k, (v, v))
    return out


def main():
    src = json.loads(IN.read_text())
    res = {}
    for key, entry in src.items():
        name = key.split("_sigma")[0]
        if name not in ("heat", "advection", "burgers") or not entry.get("certified"):
            continue
        sigma = float(key.split("_sigma")[1])
        ivs = intervals_of(entry)
        u, clean, x, t = fresh_solution(name, sigma)
        t0 = time.time()
        r = verify(u, u[:, 0], x, t, ivs, sigma=sigma, rtol=RTOL,
                   u_clean=clean)
        r["intervals"] = {k: list(v) for k, v in ivs.items()}
        r["seconds"] = round(time.time() - t0, 1)
        res[key] = r
        print(f"{key:24s} law={'OK ' if r['verified'] else 'FAIL':4s} "
              f"data={'OK ' if r.get('data_verified') else 'FAIL':4s} "
              f"law_err={r.get('law_max_err', 0):.1e} "
              f"outside={r.get('n_outside')}/{r.get('n_points')} "
              f"env={r.get('envelope_width_med', 0):.1e} "
              f"ic_noise={r.get('ic_noise_bound', 0):.1e} {r['seconds']}s",
              flush=True)

        # V2: the same rung with every coefficient pushed WRONG_FACTOR interval
        # half-widths off centre must fail
        wrong = {}
        for k, (lo, hi) in ivs.items():
            c = 0.5 * (lo + hi)
            hw = max(0.5 * (hi - lo), abs(c) * 1e-6)
            wrong[k] = (c + WRONG_FACTOR * hw, c + WRONG_FACTOR * hw)
        rw = verify(u, u[:, 0], x, t, wrong, sigma=sigma, rtol=RTOL,
                    u_clean=clean)
        res[key + "_WRONG"] = {"verified": rw["verified"],
                               "n_outside": rw.get("n_outside"),
                               "worst_excess": rw.get("worst_excess"),
                               "coefficients": {k: v[0] for k, v in wrong.items()}}
        print(f"{'  ^ 10x-wrong control':24s} "
              f"{'VERIFIED (BAD)' if rw['verified'] else 'FAILED (good)':14s} "
              f"law_outside={rw.get('law_n_outside')}/{rw.get('n_points')}",
              flush=True)
    OUT.write_text(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
