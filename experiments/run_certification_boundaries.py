"""BND1: bounded, append-only certification-path study. See registration.

Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1. No oracle/API/network calls.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import time

import numpy as np
import sympy as sp

from lagh.base import eval_expr
from lagh.certify import check, epsilon, free_atoms, parameter_interval
from lagh.engine import discover

ROOT = Path(__file__).resolve().parents[1]
REGIMES = {
    "clean": (0., 1e-12), "noise6": (1e-6, 1e-12),
    "noise3": (1e-3, 1e-12), "floor6": (0., 1e-6),
    "floor3": (0., 1e-3),
}
FAMILIES = ("affine", "rational", "power", "small-term", "mono3",
            "small-term3", "affine3", "sphere3", "basis-independent",
            "basis-near", "basis-exact", "null")
S = list(sp.symbols("x_0:3"))
TRUTHS = {
    "affine": 13*S[0]/10 + sp.Rational(7, 10),
    "rational": (3*S[0]+2)/(S[0]+4),
    "power": sp.Rational(13, 10)*S[0]**sp.Rational(3, 2),
    "small-term": S[0] + S[0]**2/100000,
    "mono3": sp.Rational(13, 10)*S[0]*sp.sqrt(S[1])/S[2],
    "small-term3": sp.Rational(13, 10)*S[0]*sp.sqrt(S[1])/S[2] + S[0]**2/100000,
    "affine3": 13*S[0]/10 - 7*S[1]/10 + 3*S[2]/10,
    "sphere3": 13*S[0]/10 - 7*S[1]/10 + 3*S[2]/10,
    "basis-independent": 13*S[0]/10 - 7*S[1]/10,
    "basis-near": 13*S[0]/10 - 7*S[1]/10,
    "basis-exact": 13*S[0]/10 - 7*S[1]/10,
}


def digest(data):
    return hashlib.sha256(np.ascontiguousarray(data, dtype="<f8").tobytes()).hexdigest()


def source_hashes():
    paths = sorted((ROOT / "lagh").rglob("*.py")) + [Path(__file__),
        ROOT / "docs/CERTIFICATION_BOUNDARIES_REGISTRATION.md"]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def data(family, regime, seed, *, n=240, salt=0, wide=False):
    rng = np.random.default_rng(np.random.SeedSequence([seed, salt]))
    dim = 3 if family.endswith("3") else 2 if family.startswith("basis-") else 1
    X = rng.uniform(*((.25, 6.) if wide else (.5, 3.)), (n, dim))
    if family == "sphere3":
        X = rng.normal(size=(n, 3))
        X /= np.linalg.norm(X, axis=1, keepdims=True)
    if family in ("basis-near", "basis-exact"):
        X[:, 1] = X[:, 0] + (1e-6 * rng.uniform(-1, 1, n)
                              if family == "basis-near" else 0)
    truth = TRUTHS.get(family)
    y0 = eval_expr(truth, S[:dim], X) if truth is not None else rng.uniform(.5, 3., n)
    sigma, floor = REGIMES[regime]
    y = y0.copy()
    if sigma:
        y += rng.normal(size=n) * sigma * np.abs(y0)
    elif regime.startswith("floor"):
        y += floor / 4 * np.sin(np.arange(n) * np.sqrt(2))
    return X, y, y0, truth


class Limit(BaseException):
    pass


def alarm(_signum, _frame):
    raise Limit("registered time limit")


def finite_json(value):
    if isinstance(value, dict):
        return {str(k): finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite_json(v) for v in value]
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, sp.Basic):
        return str(value)
    return value


def interval_diagnostics(expr, syms, X, y, eps, truth):
    """Coordinate slices only; preserve limits rather than claim joint coverage."""
    atoms = free_atoms(expr)
    if len(atoms) > 10:
        return {"omitted": "more than 10 numeric atoms", "n_atoms": len(atoms)}
    slices = []
    for a in atoms:
        iv = parameter_interval(expr, syms, X, y, eps, a, iters=20)
        row = {"atom": str(a), "value": float(a), "bounds": iv}
        if iv is not None:
            row["endpoints_feasible"] = [check(expr.xreplace({a: sp.Float(v)}),
                syms, X, y, eps)["certified"] for v in iv]
            # Match this numeric coordinate to the known generator only if
            # replacing it alone makes the two expressions symbolically equal.
            matching = [b for b in free_atoms(truth) if
                        sp.expand(expr.xreplace({a: b}) - truth) == 0] if truth is not None else []
            if len(matching) == 1:
                v = float(matching[0])
                row.update(truth_value=v, contains_truth=iv[0] <= v <= iv[1])
        slices.append(row)
    bounded = [(a, r["bounds"]) for a, r in zip(atoms, slices) if r["bounds"]]
    corners = None
    if 1 < len(bounded) <= 5:
        corners = [check(expr.xreplace({a: sp.Float(iv[i]) for (a, iv), i
                                       in zip(bounded, bits)}), syms, X, y, eps)["certified"]
                   for bits in itertools.product((0, 1), repeat=len(bounded))]
    return {"kind": "coordinate_slices_other_atoms_fixed", "slices": slices,
            "simultaneous_corners_feasible": corners,
            "truth_inclusion_scope": "only uniquely matched single-coordinate replacements"}


def score(expr, family, regime, seed, X, y, truth, cert_idx, out):
    syms = S[:X.shape[1]]
    sigma, floor = REGIMES[regime]
    out["all_observed_rows_fit"] = check(expr, syms, X, y,
                         epsilon(y, sigma=sigma, floor_abs=floor))["certified"]
    for label, wide, salt in (("fresh", False, 1), ("extended", True, 2)):
        Z, _, z, _ = data(family, regime, seed, n=512, salt=salt, wide=wide)
        p = eval_expr(expr, syms, Z)
        out[label + "_relative_error"] = (float(np.max(np.abs(p-z)) / max(np.max(np.abs(z)), 1e-300))
            if p is not None and np.all(np.isfinite(p)) else None)
        out[label + "_truth_band_misses"] = (int(np.sum(np.abs(p-z) >
            epsilon(z, sigma=sigma, floor_abs=floor))) if p is not None else 512)
    if truth is not None:
        # Domain-restricted scoring: canonicalize the exact linear variety.
        t, e = truth, expr
        if family == "basis-exact":
            t, e = t.subs(S[1], S[0]), e.subs(S[1], S[0])
        delta = sp.cancel(e-t)
        out["symbolically_equal"] = bool(delta == 0)
        if t.is_polynomial(*syms) and e.is_polynomial(*syms):
            tp, ep = sp.Poly(t, *syms), sp.Poly(e, *syms)
            out["same_polynomial_support"] = set(tp.monoms()) == set(ep.monoms())
            out["max_coefficient_error"] = max(abs(float(tp.coeff_monomial(m)-ep.coeff_monomial(m)))
                                               for m in set(tp.monoms()) | set(ep.monoms()))
        out["generator_mismatch"] = bool(
            out.get("same_polynomial_support") is False or
            out.get("max_coefficient_error", 0) > 1e-8 or
            out["extended_relative_error"] is None or out["extended_relative_error"] > 1e-8)
    out["interval_diagnostics"] = interval_diagnostics(expr, syms, X[cert_idx],
        y[cert_idx], epsilon(y[cert_idx], sigma=sigma, floor_abs=floor), truth)
    return out


def run_case(args):
    family, regime, seed, folder = args
    name = f"{family}__{regime}__{seed}"
    path = Path(folder) / (name + ".json")
    if path.exists():
        return name, "retained"
    X, y, y0, truth = data(family, regime, seed)
    sigma, floor = REGIMES[regime]
    a, b = 144, 192
    idx = np.arange(b, len(X))
    rec = {"case": name, "family": family, "regime": regime, "seed": seed,
           "input_hash": digest(np.column_stack([X, y])), "truth": str(truth),
           "sigma": sigma, "floor_abs": floor,
           "truth_within_observed_band": bool(np.all(np.abs(y-y0) <= epsilon(y, sigma=sigma, floor_abs=floor))),
           "truth_within_cert_band": bool(np.all(np.abs(y[idx]-y0[idx]) <= epsilon(y[idx], sigma=sigma, floor_abs=floor))),
           "status": None}
    old = signal.signal(signal.SIGALRM, alarm)
    start = time.monotonic()
    signal.setitimer(signal.ITIMER_REAL, 20)
    try:
        r = discover(X[:a], y[:a], X[a:b], y[a:b], X[b:], y[b:], sigma=sigma,
                     floor_abs=floor, max_tier=3,
                     declared_basis=family.startswith("basis-"),
                     linear_basis=family.startswith("basis-"))
        rec.update(status="certified" if r.certificate.certified else "refused",
                   certificate=asdict(r.certificate), law=str(r.expr),
                   tier=r.tier, n_candidates=r.n_candidates,
                   route="prepass" if "CAP-S cheap pre-pass" in r.certificate.notes else "tier-loop")
    except Limit:
        rec["status"] = "timeout"
    except Exception as exc:
        rec.update(status="error", error=f"{type(exc).__name__}: {exc}")
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        rec["discovery_seconds"] = time.monotonic()-start
    if rec["status"] == "certified":
        start = time.monotonic()
        signal.setitimer(signal.ITIMER_REAL, 15)
        rec["score"] = {}
        try:
            score(r.expr, family, regime, seed, X, y, truth, idx, rec["score"])
        except Limit:
            rec["score_omitted"] = "15-second diagnostic limit"
        except Exception as exc:
            rec["score_omitted"] = f"{type(exc).__name__}: {exc}"
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            rec["scoring_seconds"] = time.monotonic()-start
    signal.signal(signal.SIGALRM, old)
    path.write_text(json.dumps(finite_json(rec), indent=2, allow_nan=False)+"\n")
    return name, rec["status"]


def summary(folder, seeds):
    expected = {f"{f}__{r}__{s}" for f in FAMILIES for r in REGIMES for s in seeds}
    rows = [json.loads((folder / (name+".json")).read_text()) for name in sorted(expected)
            if (folder / (name+".json")).exists()]
    out = {"expected": len(expected), "present": len(rows), "by_regime": {}}
    for regime in REGIMES:
        group = [r for r in rows if r["regime"] == regime]
        certs = [r for r in group if r["status"] == "certified"]
        out["by_regime"][regime] = {
            "n": len(group),
            **{k: sum(r["status"] == k for r in group) for k in ("certified", "refused", "timeout", "error")},
            "generator_mismatch": [r["case"] for r in certs if r.get("score", {}).get("generator_mismatch")],
            "missing_scores": [r["case"] for r in certs if "score_omitted" in r],
            "null_certificates": [r["case"] for r in certs if r["family"] == "null"],
            "median_discovery_seconds": float(np.median([r["discovery_seconds"] for r in group])) if group else None,
            "refusal_reasons": {reason: sum(r.get("certificate", {}).get("abstain") == reason for r in group)
                for reason in sorted({r["certificate"]["abstain"] for r in group if r["status"] == "refused"})},
        }
    (folder / "summary.json").write_text(json.dumps(out, indent=2)+"\n")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--seeds", type=int, nargs=2, default=[710, 711])
    ap.add_argument("--workers", type=int, choices=range(1, 5), default=4)
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()
    if args.summarize:
        print(json.dumps(summary(args.output, args.seeds), indent=2))
        return
    if os.environ.get("OPENBLAS_NUM_THREADS") != "1" or os.environ.get("OMP_NUM_THREADS") != "1":
        raise SystemExit("set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 before launch")
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = {"sources": source_hashes(), "seeds": args.seeds,
                "python": platform.python_version(), "numpy": np.__version__,
                "sympy": sp.__version__, "platform": platform.platform(),
                "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()}
    dest = args.output / "manifest.json"
    if dest.exists() and json.loads(dest.read_text()) != manifest:
        raise SystemExit("source/environment differs from retained manifest; choose a new arm")
    dest.write_text(json.dumps(manifest, indent=2)+"\n")
    jobs = [(f, r, s, str(args.output)) for f in FAMILIES for r in REGIMES for s in args.seeds]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in as_completed([pool.submit(run_case, job) for job in jobs]):
            print(*result.result(), flush=True)
    print(json.dumps(summary(args.output, args.seeds), indent=2))


if __name__ == "__main__":
    main()
