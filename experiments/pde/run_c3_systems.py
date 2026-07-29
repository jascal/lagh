"""PDE dev campaign C3 -- systems of PDEs (registered:
docs/CASE_STUDY_PDE_C3.md, direction docs/DIRECTION_PDE_SYSTEMS.md).

One row set, one target per equation, features spanning every field. The
vocabulary of each stage is a REGISTERED list, never a generated cross-product:
cross terms multiply fast and |H| enters alpha directly.

Scores the registered predictions Y1-Y5.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from experiments.pde import systems_fields as S                   # noqa: E402
from lagh.pdesystem import (agreement, assemble, conjoin,          # noqa: E402
                            discover_equation, truth_check, weakest)
from lagh.weakform import (Term, field_terms, multiscale_patches,  # noqa: E402
                           multiscale_patches_nd)

OUT = Path("experiments/results/pde_c3.json")
SCALES = [(24, 12), (32, 16), (40, 20)]
FAMILY = dict(n_x=4, n_t=3)
P_BUMP = 16
ONE = Term("1", 0, 0, "1")


def patch_fn(scales=None, family=None):
    sc = scales or SCALES
    fam = family or FAMILY
    return lambda coords: multiscale_patches(coords[0], coords[1], sc, **fam)


# --------------------------------------------------------------------------
# The registered vocabularies and the true systems
# --------------------------------------------------------------------------

LINEAR_TERMS = (field_terms("u", ["u_t", "u_xx", "u_x", "u"])
                + field_terms("v", ["u_t", "u_xx", "u_x", "u"]) + [ONE])
LINEAR_TRUTH = {"u:u_t": {"u:u_xx": 0.1, "v:u": 0.5},
                "v:u_t": {"v:u_xx": 0.05, "u:u": -0.3}}

FHN_TERMS = (field_terms("u", ["u_t", "u_xx", "u_x", "u", "u^3"])
             + field_terms("v", ["u_t", "u_xx", "u_x", "u"])
             + [Term("u*v", 0, 0, "u*v"), ONE])
FHN_TRUTH = {
    "u:u_t": {"u:u_xx": S.FHN_D, "u:u": 1.0, "u:u^3": -1.0 / 3.0,
              "v:u": -1.0, "1": S.FHN_I},
    "v:u_t": {"u:u": S.FHN_EPS, "v:u": -S.FHN_EPS * S.FHN_B,
              "1": S.FHN_EPS * S.FHN_A}}

BRU_TERMS = (field_terms("u", ["u_t", "u_xx", "u_x", "u", "u^2"])
             + field_terms("v", ["u_t", "u_xx", "u_x", "u"])
             + [Term("u*v", 0, 0, "u*v"), Term("u^2*v", 0, 0, "u**2*v"), ONE])
BRU_TRUTH = {
    "u:u_t": {"u:u_xx": S.BRU_DU, "u:u": -(S.BRU_B + 1.0), "u^2*v": 1.0,
              "1": S.BRU_A},
    "v:u_t": {"v:u_xx": S.BRU_DV, "u:u": S.BRU_B, "u^2*v": -1.0}}

SW_TERMS = [
    Term("h_t", 0, 1, "h"), Term("(hu)_t", 0, 1, "h*u"),
    Term("(hu)_x", 1, 0, "h*u"), Term("(hu^2)_x", 1, 0, "h*u**2"),
    Term("(h^2)_x", 1, 0, "h**2"), Term("h_x", 1, 0, "h"),
    Term("u_x", 1, 0, "u"), Term("h_xx", 2, 0, "h"), Term("u_xx", 2, 0, "u"),
    Term("h", 0, 0, "h"), Term("u", 0, 0, "u"), Term("hu", 0, 0, "h*u"), ONE]
SW_TRUTH = {"h_t": {"(hu)_x": -1.0},
            "(hu)_t": {"(hu^2)_x": -1.0, "(h^2)_x": -0.5 * S.SW_G}}

# Stage 4: 2-D space. Four fields are handed over precisely because three of
# them are tied to the fourth by MACHINE-EXACT relations -- psi_xx + psi_yy = -w
# and u_x + v_y = 0 -- which is the constrained-input situation the engine
# closed for algebraic constraints, appearing here as physics.
NS_TERMS = [
    Term("w_t", alpha=(0, 0, 1), gexpr="w"),
    Term("w_xx", alpha=(2, 0, 0), gexpr="w"),
    Term("w_yy", alpha=(0, 2, 0), gexpr="w"),
    Term("(uw)_x", alpha=(1, 0, 0), gexpr="u*w"),
    Term("(vw)_y", alpha=(0, 1, 0), gexpr="v*w"),
    Term("psi_xx", alpha=(2, 0, 0), gexpr="psi"),
    Term("psi_yy", alpha=(0, 2, 0), gexpr="psi"),
    Term("w", alpha=(0, 0, 0), gexpr="w"),
    Term("u_x", alpha=(1, 0, 0), gexpr="u"),
    Term("v_y", alpha=(0, 1, 0), gexpr="v"),
    Term("1", alpha=(0, 0, 0), gexpr="1")]
NS_TRUTH = {"w_t": {"w_xx": S.NS_NU, "w_yy": S.NS_NU,
                    "(uw)_x": -1.0, "(vw)_y": -1.0}}


# --------------------------------------------------------------------------

def score(eq, truth, rows=None, sigma=0.0):
    """Y1 per equation: is the SUPPORT the true one, does every reported
    interval contain the truth -- and, when the support differs, does the
    certified law still AGREE with the truth everywhere the certificate applies?

    The last question is not a softener. A certificate is a domain claim, so a
    law that agrees with the truth on its whole certified domain is not wrong;
    it is under-determined, and the instrument is supposed to say so (the
    engine's constrained-input path does, by marking the certificate
    domain-restricted). What would be confident-wrong is a law that DISAGREES
    with the truth by more than its own declared band.
    """
    if not eq.get("certified"):
        return {"support_ok": None, "intervals_contain_truth": None,
                "confident_wrong": False}
    got = eq.get("coefficients") or {}
    sup = set(got) == set(truth)
    ivs = eq.get("intervals") or {}
    contain, worst = True, 0.0
    for k, v in truth.items():
        iv = ivs.get(k)
        if iv is None:
            contain = False
            continue
        lo, hi = iv
        pad = 1e-12 * max(1.0, abs(v))
        if not (lo - pad <= v <= hi + pad):
            contain = False
        worst = max(worst, abs(got.get(k, 0.0) - v) / max(abs(v), 1e-30))
    out = {"support_ok": bool(sup),
           "intervals_contain_truth": bool(contain),
           "worst_rel_error": float(worst)}
    dr = any("domain-restricted certificate" in str(n)
             for n in eq.get("notes", []))
    out["domain_restricted"] = bool(dr)
    if rows is not None:
        out.update(agreement(rows, eq["target"], got, truth, sigma=sigma))
    agrees = out.get("agrees_on_certified_domain", sup and contain)
    out["confident_wrong"] = bool(not agrees)
    # a support that differs while agreeing on the whole certified domain is an
    # under-determination; it is only a SILENT one when the certificate did not
    # flag its own domain restriction
    out["degenerate_support_unflagged"] = bool(agrees and not sup and not dr)
    return out


def run_stage(name, terms, truth, solutions, sigma, *, field_err=0.0, seed=0,
              scales=None, family=None, max_tier=3, holdout=True,
              patch_fn_override=None):
    t0 = time.time()
    rows = assemble(solutions, terms,
                    patch_fn_override or patch_fn(scales, family), sigma=sigma,
                    field_err=field_err, p=P_BUMP)
    if rows is None:
        return {"stage": name, "certified": False,
                "abstain": "all-patches-rejected", "sigma": sigma}
    eqs = []
    for target in truth:
        # the discipline the scoping probe bought: ask whether the TRUTH
        # certifies against its own band BEFORE reading any abstain as a finding
        tc = truth_check(rows, target, truth[target], sigma=sigma)
        eq = discover_equation(rows, target, sigma=sigma, seed=seed,
                               max_tier=max_tier, holdout=holdout)
        eq["truth"] = tc
        eq["score"] = score(eq, truth[target], rows=rows, sigma=sigma)
        eqs.append(eq)
    cert = conjoin(eqs)
    wt, wa = weakest(eqs)
    return {"stage": name, "sigma": sigma, "field_err": field_err,
            "n_rows": int(len(rows.A)), "n_terms": len(rows.names),
            "n_solutions": rows.n_solutions, "rejected": int(rows.rejected),
            "vocabulary": rows.names,
            "equations": eqs,
            "system_certified": bool(all(e.get("certified") for e in eqs)),
            "alpha_log10_total": cert.alpha_log10_total,
            "weakest_equation": wt, "weakest_alpha_log10": wa,
            "shared_constants": cert.shared,
            "seconds": round(time.time() - t0, 1)}


def linear_solutions(n_ic=4, n_modes=3, seed0=0):
    return [S.linear_pair(seed0 + i, modes=tuple(range(1, n_modes + 1)))
            for i in range(n_ic)]


def solved_solutions(fn, n_ic=4, n_modes=2, seed0=0):
    """Solved stages: fields plus ONE declared field error for the family (the
    largest of the per-solution ladder bounds -- a uniform declaration)."""
    out, err = [], 0.0
    for i in range(n_ic):
        f, c, e = fn(seed0 + i, n_modes=n_modes)
        if f is None:
            continue
        out.append((f, c))
        err = max(err, e)
    return out, err


def main(only=None):
    res = {}
    sigmas = [0.0, 1e-6, 1e-4]

    # ---- stage 1: linear coupled (exact fields, no solver error at all)
    if only in (None, "linear"):
        for sigma in sigmas:
            k = f"S1_linear_sigma{sigma:g}"
            res[k] = run_stage("linear", LINEAR_TERMS, LINEAR_TRUTH,
                               linear_solutions(), sigma)
            report(k, res[k])

    # ---- stage 2a: FitzHugh-Nagumo (nonlinear in one field, linear coupling)
    if only in (None, "fhn"):
        sols, ferr = solved_solutions(S.fitzhugh_nagumo)
        for sigma in sigmas:
            k = f"S2a_fhn_sigma{sigma:g}"
            res[k] = run_stage("fhn", FHN_TERMS, FHN_TRUTH, sols, sigma,
                               field_err=ferr)
            report(k, res[k])

    # ---- stage 2b: Brusselator (the first CROSS term, u^2 v)
    if only in (None, "brusselator"):
        sols, ferr = solved_solutions(S.brusselator)
        for sigma in sigmas:
            k = f"S2b_brusselator_sigma{sigma:g}"
            res[k] = run_stage("brusselator", BRU_TERMS, BRU_TRUTH, sols,
                               sigma, field_err=ferr)
            report(k, res[k])

    # ---- stage 3: shallow water (conservation form with a real flux)
    if only in (None, "shallow"):
        sols, ferr = solved_solutions(S.shallow_water)
        for sigma in sigmas:
            k = f"S3_shallow_sigma{sigma:g}"
            res[k] = run_stage("shallow_water", SW_TERMS, SW_TRUTH, sols,
                               sigma, field_err=ferr)
            report(k, res[k])

    # ---- stage 4: Navier-Stokes vorticity-streamfunction, 2-D space
    if only in (None, "ns"):
        sols, ferr = [], 0.0
        for i in range(4):
            f, c, e = S.ns_vorticity(i)
            if f is None:
                continue
            sols.append((f, c))
            ferr = max(ferr, e)
        pf = (lambda coords: multiscale_patches_nd(
            coords, [(12, 12, 8), (16, 16, 8)], (4, 4, 2)))
        for sigma in (0.0, 1e-6):
            k = f"S4_ns_sigma{sigma:g}"
            res[k] = run_stage("ns_vorticity", NS_TERMS, NS_TRUTH, sols, sigma,
                               field_err=ferr, patch_fn_override=pf)
            report(k, res[k])

    # ---- Y3: single-solution data must refuse, for systems too
    if only in (None, "single"):
        res["Y3_linear_single_solution"] = run_stage(
            "linear-single", LINEAR_TERMS, LINEAR_TRUTH,
            linear_solutions(n_ic=1), 0.0)
        report("Y3_linear_single_solution", res["Y3_linear_single_solution"])
        res["Y3_linear_single_rowsplit"] = run_stage(
            "linear-single-rowsplit", LINEAR_TERMS, LINEAR_TRUTH,
            linear_solutions(n_ic=1), 0.0, holdout=False)
        report("Y3_linear_single_rowsplit", res["Y3_linear_single_rowsplit"])

    # ---- nulls: two smooth fields solving no system
    if only in (None, "null"):
        x = np.linspace(0.0, 2 * np.pi, 257)
        t = np.linspace(0.0, 1.0, 81)
        sols = [(S.smooth_random_pair(x, t, seed=s), (x, t)) for s in range(4)]
        for sigma in (0.0, 1e-4):
            k = f"N_smooth_pair_sigma{sigma:g}"
            res[k] = run_stage("null-smooth", LINEAR_TERMS, LINEAR_TRUTH,
                               sols, sigma)
            report(k, res[k])
        rng = np.random.default_rng(0)
        sols = [({"u": rng.normal(0, 1, (len(x), len(t))),
                  "v": rng.normal(0, 1, (len(x), len(t)))}, (x, t))
                for _ in range(4)]
        res["N_iid_pair"] = run_stage("null-iid", LINEAR_TERMS, LINEAR_TRUTH,
                                      sols, 0.0)
        report("N_iid_pair", res["N_iid_pair"])

    # ---- Y4: how much spectral richness does each stage need?
    if only in (None, "modes"):
        for n_modes in (1, 2, 3):
            k = f"Y4_linear_modes{n_modes}"
            res[k] = run_stage("linear", LINEAR_TERMS, LINEAR_TRUTH,
                               linear_solutions(n_modes=n_modes), 0.0)
            report(k, res[k])
            sols, ferr = solved_solutions(S.brusselator, n_modes=n_modes)
            k = f"Y4_brusselator_modes{n_modes}"
            res[k] = run_stage("brusselator", BRU_TERMS, BRU_TRUTH, sols, 0.0,
                               field_err=ferr)
            report(k, res[k])

    prev = json.loads(OUT.read_text()) if OUT.exists() and only else {}
    prev.update(res)
    OUT.write_text(json.dumps(prev, indent=1))
    return 0


def report(k, v):
    if "equations" not in v:
        print(f"{k:30s} {v.get('abstain')}", flush=True)
        return
    print(f"{k:30s} rows={v['n_rows']:4d} "
          f"{'SYSTEM-CERT' if v['system_certified'] else 'partial    '} "
          f"alpha<=1e{v['alpha_log10_total']:.0f}" if v["alpha_log10_total"]
          else f"{k:30s} rows={v['n_rows']:4d} no-alpha", flush=True)
    for e in v["equations"]:
        sc = e.get("score", {})
        tr = e.get("truth", {})
        print(f"   {e['target']:12s} "
              f"{'CERT   ' if e.get('certified') else 'ABSTAIN'} "
              f"{str(e.get('abstain') or ''):12s} "
              f"truth/band={tr.get('truth_max_ratio', float('nan')):.2e} "
              f"sup={sc.get('support_ok')} iv={sc.get('intervals_contain_truth')} "
              f"{str(e.get('law', ''))[:70]}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else None))
