"""PDE dev campaign C1 -- the noise ladder (registered: docs/CASE_STUDY_PDE_C1.md).

Same analytic solutions as C0, now with declared field noise, certified through
the errors-in-variables epsilon (per-candidate band from the patch Gram) and
reported with interval parameters where an exact value is not what the data
determines.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import fields as F                          # noqa: E402
from experiments.pde.run_c0 import (advection_solutions,         # noqa: E402
                                    burgers_solutions, heat_solutions)
from lagh.engine import discover                                 # noqa: E402
from lagh.weakform import PatchEpsilon, build, make_patches      # noqa: E402

OUT = Path("experiments/results/pde_c1.json")
FAMILY = dict(nx_half=24, nt_half=12, n_x=6, n_t=4)
P_BUMP = 16
SIGMAS = [0.0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]

# system -> (library, target, the true law in column terms)
SYSTEMS = {
    "heat":      (["u_t", "u_xx", "u*u_x", "u_x", "u", "1"], {"u_xx": 0.1}),
    "advection": (["u_t", "u_x", "u_xx", "u*u_x", "u", "1"], {"u_x": -0.7}),
    "burgers":   (["u_t", "u_xx", "u*u_x", "u_x", "u", "1"],
                  {"u_xx": 0.2, "u*u_x": -1.0}),
}


def solutions(name, n_ic=4):
    return {"heat": heat_solutions, "advection": advection_solutions,
            "burgers": burgers_solutions}[name](n_ic)


def assemble(sols, terms, sigma, seed=0):
    """Noisy weak system pooled over solutions + the errors-in-variables model."""
    rng = np.random.default_rng(seed)
    A, det, gram, sol, rejected = [], [], [], [], 0
    for i, (u, x, t) in enumerate(sols):
        un = u + (rng.normal(0, sigma, u.shape) if sigma > 0 else 0.0)
        s = build(un, x, t, terms, make_patches(x, t, **FAMILY), p=P_BUMP,
                  sigma=sigma)
        rejected += s.rejected
        if len(s.A) == 0:
            continue
        A.append(s.A); det.append(s.quad + s.roundoff); gram.append(s.gram)
        sol.append(np.full(len(s.A), i))
    if not A:
        return None
    A = np.vstack(A)
    j = terms.index("u_t")
    cols = [k for k in range(len(terms)) if k != j]
    m = PatchEpsilon(terms, "u_t", A[:, j], A[:, cols], np.vstack(det),
                     np.concatenate(gram), sigma=sigma, floor_abs=0.0,
                     coeff_max=2.0)
    return A[:, cols], A[:, j], m, np.concatenate(sol), rejected


def score(expr, feat, truth):
    """(support_ok, coefficients) for the certified law against the true law."""
    if expr is None:
        return None, {}
    e = sp.expand(expr)
    syms = [sp.Symbol(f"x_{i}") for i in range(len(feat))]
    got = {}
    for i, nm in enumerate(feat):
        c = e.coeff(syms[i])
        if c.free_symbols:
            return False, {}
        if abs(float(c)) > 0:
            got[nm] = float(c)
    const = float(e - sum(sp.Float(v) * syms[feat.index(k)]
                          for k, v in got.items()))
    if abs(const) > 0:
        got["1"] = const
    return set(got) == set(truth), got


def run(name, sigma, seed=0):
    terms, truth = SYSTEMS[name]
    got = assemble(solutions(name), terms, sigma, seed=seed)
    if got is None:
        return {"certified": False, "abstain": "all-patches-rejected"}
    X, y, m, sol, rejected = got
    feat = [n for n in terms if n != "u_t"]
    held = np.unique(sol)[-1]
    tr = np.random.default_rng(seed).permutation(np.where(sol != held)[0])
    ce = np.where(sol == held)[0]
    a = int(0.75 * len(tr))
    t0 = time.time()
    r = discover(X[tr[:a]], y[tr[:a]], X[tr[a:]], y[tr[a:]], X[ce], y[ce],
                 sigma=sigma, eps_model=m.subset(ce), max_tier=3,
                 declared_basis=True, band_sel=m.subset(tr[a:])(None))
    c = r.certificate
    sup_ok, coeffs = score(r.expr if c.certified else None, feat, truth)
    law = str(r.expr) if c.certified else ""
    for i, nm in enumerate(feat):
        law = law.replace(f"x_{i}", f"[{nm}]")
    iv = [n for n in c.notes
          if "interval-parameter" in str(n) or "parameter precision" in str(n)]
    # does every certified coefficient sit within (interval, or exactly) truth?
    err = {k: (coeffs.get(k, 0.0) - v) / abs(v) for k, v in truth.items()}
    return {"certified": bool(c.certified), "law": law, "abstain": c.abstain,
            "alpha_log10": c.alpha_log10, "support_ok": sup_ok,
            "coefficients": coeffs, "rel_error": err,
            "interval_note": (str(iv[0]) if iv else None),
            "n_patches": int(len(y)), "n_cert_patches": int(len(ce)),
            "patches_rejected": int(rejected), "tier": r.tier,
            "notes": [str(n)[:200] for n in c.notes][:2],
            "seconds": round(time.time() - t0, 1)}


def run_null(kind, sigma, seed=0):
    terms = SYSTEMS["heat"][0]
    x, t = F.grid()
    rng = np.random.default_rng(7)
    if kind == "iid":
        sols = [(rng.normal(0, 1, (len(x), len(t))), x, t)]
    else:
        sols = [(F.smooth_random(x, t, seed=s), x, t) for s in range(4)]
    got = assemble(sols, terms, sigma, seed=seed)
    if got is None:
        return {"certified": False, "abstain": "all-patches-rejected"}
    X, y, m, sol, rejected = got
    if len(np.unique(sol)) < 2:
        return {"certified": False, "abstain": "single-solution",
                "patches_rejected": int(rejected)}
    held = np.unique(sol)[-1]
    tr = np.random.default_rng(seed).permutation(np.where(sol != held)[0])
    ce = np.where(sol == held)[0]
    a = int(0.75 * len(tr))
    r = discover(X[tr[:a]], y[tr[:a]], X[tr[a:]], y[tr[a:]], X[ce], y[ce],
                 sigma=sigma, eps_model=m.subset(ce), max_tier=3,
                 declared_basis=True, band_sel=m.subset(tr[a:])(None))
    return {"certified": bool(r.certificate.certified),
            "abstain": r.certificate.abstain,
            "law": str(r.expr) if r.certificate.certified else "",
            "patches_rejected": int(rejected)}


def main(only=None):
    res = {}
    for name in SYSTEMS:
        if only and only != name:
            continue
        for sigma in SIGMAS:
            k = f"{name}_sigma{sigma:g}"
            res[k] = run(name, sigma)
            v = res[k]
            print(f"{k:24s} {'CERT' if v.get('certified') else 'ABSTAIN':8s} "
                  f"{str(v.get('abstain') or ''):12s} sup_ok={v.get('support_ok')} "
                  f"{str(v.get('coefficients'))[:60]}", flush=True)
    for kind in ("iid", "smooth"):
        for sigma in (0.0, 1e-6, 1e-3):
            k = f"null_{kind}_sigma{sigma:g}"
            res[k] = run_null(kind, sigma)
            print(f"{k:24s} {'CERT' if res[k].get('certified') else 'ABSTAIN':8s} "
                  f"{str(res[k].get('abstain') or ''):12s}", flush=True)
    # M6: single traveling wave under noise
    xw = np.linspace(-8.0, 8.0, 257)
    _, t = F.grid()
    for sigma in (0.0, 1e-4):
        got = assemble([(F.burgers_wave(xw, t), xw, t)], SYSTEMS["heat"][0],
                       sigma)
        k = f"single_wave_sigma{sigma:g}"
        res[k] = {"certified": False, "abstain": "single-solution",
                  "n_patches": 0 if got is None else int(len(got[1]))}
        print(f"{k:24s} ABSTAIN  single-solution", flush=True)
    OUT.write_text(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else None))
