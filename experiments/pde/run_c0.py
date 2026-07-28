"""PDE dev campaign C0 (registered: docs/CASE_STUDY_PDE_DEV.md).

Weak-form certification against EXACT analytic solutions -- no solver, no
downloads. Scores the registered predictions N1-N7: nulls certify nothing,
heat/advection/Burgers certify from pooled multi-IC patches, a single traveling
wave must refuse, and the verdict must not depend on the patch family.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import fields as F                      # noqa: E402
from lagh.engine import discover                             # noqa: E402
from lagh.weakform import build, make_patches                # noqa: E402

OUT = Path("experiments/results/pde_c0.json")

LIB = ["u_t", "u_xx", "u*u_x", "u_x", "u", "u_xxx", "1"]
TARGET = "u_t"
# every declared band assumes |coefficient| <= this; coeff_audit() checks it
COEFF_MAX = 2.0

# The bump exponent p sets how many derivatives of phi vanish at the support
# edge, and so how fast the patch quadrature converges for the HIGH-derivative
# weights. Measured on the heat field: the u_xxx bound falls 1.5e-6 -> 7.5e-11
# going p = 8 -> 16 at the same patch size, and the declared band is set by the
# roughest term in the library, so p = 16 is the default.
FAMILIES = {
    "A": dict(nx_half=24, nt_half=12, n_x=6, n_t=4, p=16),
    "B": dict(nx_half=32, nt_half=16, n_x=4, n_t=3, p=16),
    "C": dict(nx_half=24, nt_half=12, n_x=6, n_t=4, p=12),
}


def coeff_audit(expr, n_feat, coeff_max):
    """The band was assembled assuming every coefficient is bounded by
    `coeff_max` (weakform.declared_epsilon). VERIFY it against the law that
    actually certified -- an assumption stated and then checked, never assumed
    silently. Returns (ok, max|c|, linear?)."""
    if expr is None:
        return True, 0.0, True
    syms = [sp.Symbol(f"x_{i}") for i in range(n_feat)]
    e = sp.expand(expr)
    coeffs, rest = [], e
    for s in syms:
        c = e.coeff(s)
        if c.free_symbols:
            return False, float("inf"), False      # not linear in the columns
        coeffs.append(abs(float(c)))
        rest = rest - c * s
    if sp.simplify(rest).free_symbols:
        return False, float("inf"), False
    m = max(coeffs) if coeffs else 0.0
    return m <= coeff_max, m, True


def pooled(solutions, family="A", terms=None):
    """Stack patch rows from several exact solutions into one weak system,
    remembering WHICH solution each row came from."""
    terms = terms or LIB
    fam = dict(FAMILIES[family])
    p = fam.pop("p")
    A, band, sol, names, rejected, kept = [], [], [], None, 0, 0
    for i, (u, x, t) in enumerate(solutions):
        pa = make_patches(x, t, **fam)
        s = build(u, x, t, terms, pa, p=p)
        rejected += s.rejected
        if len(s.A) == 0:
            continue
        names = s.names
        A.append(s.A)
        band.append(s.declared_epsilon(TARGET, coeff_max=COEFF_MAX))
        sol.append(np.full(len(s.A), i))
        kept += len(s.A)
    if not A:
        return (None, None, None, terms, rejected, 0)
    return (np.vstack(A), np.concatenate(band), np.concatenate(sol),
            names, rejected, kept)


def run(solutions, family="A", terms=None, seed=0, max_tier=3, holdout=True):
    """`holdout`: certify on patches from a HELD-OUT SOLUTION, never on rows of
    the solutions that were fitted (DIRECTION_PDE.md piece 3, the multi-solution
    coherence probe).

    A PDE is a claim about a family of solutions; on one solution the on-shell
    relation (a traveling wave satisfies u_t = -c u_x identically) is a simpler
    law that fits the same data, and MEASURED: with a row split, the KdV soliton
    certifies `u_t = -u_x` -- true of that field, not the KdV equation. Holding
    out a whole solution makes a single-solution PDE claim structurally
    unattemptable, which is the honest verdict, and it is the probe the
    constrained-input closure uses in dynamics form."""
    A, band, sol, names, rejected, kept = pooled(solutions, family, terms)
    if A is None:
        return {"certified": False, "abstain": "all-patches-rejected",
                "n_patches": 0, "patches_rejected": int(rejected),
                "notes": ["every patch failed the resolution/aliasing gate: "
                          "the grid does not represent this field"]}
    n_sol = len(np.unique(sol))
    if holdout and n_sol < 2:
        return {"certified": False, "abstain": "single-solution",
                "n_patches": int(len(A)), "patches_rejected": int(rejected),
                "notes": ["a PDE claim needs patches from a HELD-OUT solution; "
                          "one solution supports only an on-shell statement "
                          "about that solution"]}
    j = names.index(TARGET)
    keep = [k for k in range(len(names)) if k != j]
    X, y = A[:, keep], A[:, j]
    feat = [names[k] for k in keep]
    rng = np.random.default_rng(seed)
    if holdout:
        held = np.unique(sol)[-1]
        tr = np.where(sol != held)[0]
        ce = np.where(sol == held)[0]
        tr = rng.permutation(tr)
        a = int(0.75 * len(tr))
        idx = np.concatenate([tr, rng.permutation(ce)])
        b = len(tr)
    else:
        idx = rng.permutation(len(y))
        a, b = int(0.6 * len(y)), int(0.8 * len(y))
    t0 = time.time()
    r = discover(X[idx[:a]], y[idx[:a]], X[idx[a:b]], y[idx[a:b]],
                 X[idx[b:]], y[idx[b:]], hard_cert=band[idx[b:]],
                 max_tier=max_tier)
    c = r.certificate
    law = str(r.expr) if c.certified else ""
    for i, nm in enumerate(feat):                     # x_0 -> the term it is
        law = law.replace(f"x_{i}", f"[{nm}]")
    ok, cmax, linear = coeff_audit(r.expr if c.certified else None, len(feat),
                                   COEFF_MAX)
    return {"certified": bool(c.certified and ok), "law": law,
            "coeff_max_declared": COEFF_MAX, "coeff_max_certified": cmax,
            "coeff_audit_ok": bool(ok), "law_linear_in_columns": bool(linear),
            "abstain": c.abstain if (c.certified and ok) else
            (c.abstain or "coeff-bound-violated"),
            "alpha_log10": c.alpha_log10, "n_patches": int(len(y)),
            "n_solutions": int(n_sol), "holdout": bool(holdout),
            "n_cert_patches": int(len(y) - b),
            "patches_rejected": int(rejected), "features": feat,
            "median_signal_to_band": float(np.median(np.abs(y) / band)),
            "tier": r.tier, "notes": [str(n)[:200] for n in c.notes][:3],
            "seconds": round(time.time() - t0, 1)}


def heat_solutions(n_ic=4, nu=0.1):
    x, t = F.grid()
    return [(F.heat(x, t, nu=nu, **kw), x, t)
            for kw in F.ic_family(F.heat, n_ic)]


def advection_solutions(n_ic=4, c=0.7):
    x, t = F.grid()
    return [(F.advection(x, t, c=c, **kw), x, t)
            for kw in F.ic_family(F.advection, n_ic)]


def burgers_solutions(n_ic=4, nu=0.2):
    x, t = F.grid()
    out = []
    for kw in F.ic_family(F.burgers_coleHopf, n_ic):
        amps = np.array(kw["amps"]); kw["amps"] = tuple(0.5 * amps / amps.sum())
        out.append((F.burgers_coleHopf(x, t, nu=nu, **kw), x, t))
    return out


def main():
    res = {}
    x, t = F.grid()

    # N1a: i.i.d. random field. N1b: smooth field solving no library PDE.
    rng = np.random.default_rng(0)
    res["N1a_iid_null"] = run([(rng.normal(0, 1, (len(x), len(t))), x, t)])
    res["N1b_smooth_null"] = run([(F.smooth_random(x, t, seed=s), x, t)
                                  for s in range(4)])

    # N2/N3/N4: multi-IC recovery of the true law
    res["N2_heat_multiIC"] = run(heat_solutions())
    res["N3_advection_multiIC"] = run(advection_solutions())
    res["N4_burgers_multiIC"] = run(burgers_solutions())

    # N5: single traveling-wave solutions -- must NOT certify a general PDE.
    # Recorded BOTH ways: the row split is what a naive pipeline does, the
    # solution holdout is what this instrument does.
    xw = np.linspace(-8.0, 8.0, 257)
    xk = np.linspace(-20.0, 20.0, 257)
    kdv_lib = ["u_t", "u_xxx", "u*u_x", "u_x", "u", "u_xx", "1"]
    res["N5a_burgers_single_wave"] = run([(F.burgers_wave(xw, t), xw, t)])
    res["N5a_burgers_single_wave_rowsplit"] = run(
        [(F.burgers_wave(xw, t), xw, t)], holdout=False)
    res["N5b_kdv_single_soliton"] = run([(F.kdv_soliton(xk, t), xk, t)],
                                        terms=kdv_lib)
    res["N5b_kdv_single_soliton_rowsplit"] = run(
        [(F.kdv_soliton(xk, t), xk, t)], terms=kdv_lib, holdout=False)
    # ...and the same PDE with several distinct solutions must still work
    res["N5c_burgers_multiIC_control"] = res["N4_burgers_multiIC"]

    # N6: the verdict must not depend on the patch family
    for fam in ("B", "C"):
        res[f"N6_heat_family_{fam}"] = run(heat_solutions(), family=fam)
        res[f"N6_burgers_family_{fam}"] = run(burgers_solutions(), family=fam)

    OUT.write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        print(f"{k:28s} {'CERT' if v.get('certified') else 'ABSTAIN':8s} "
              f"{str(v.get('abstain') or ''):12s} n={v.get('n_patches')} "
              f"s/b={v.get('median_signal_to_band', 0):.1e} "
              f"{str(v.get('law'))[:70]}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
