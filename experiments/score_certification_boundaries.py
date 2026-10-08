"""Audit BND1 completeness, original evidence, claim scope and interval endpoints.

No discovery, new observations for selection, or mutation of arm artifacts.
Generator mismatch counts are retained even after claim wording is corrected.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import sympy as sp

from experiments.run_certification_boundaries import FAMILIES, REGIMES, S, data, digest
from lagh.certify import check, epsilon


def read_arm(folder, seeds, *, require_scope):
    expected = {f"{f}__{r}__{s}" for f in FAMILIES for r in REGIMES for s in seeds}
    paths = {p.stem: p for p in folder.glob("*__*.json")}
    if set(paths) != expected:
        raise ValueError(f"incomplete/extra bank in {folder}: missing={sorted(expected-set(paths))}, extra={sorted(set(paths)-expected)}")
    rows = {k: json.loads(p.read_text()) for k, p in paths.items()}
    checks = {"input_hashes": 0, "certified_row_rechecks": 0,
              "interval_endpoint_rechecks": 0}
    issues = []
    for name, r in rows.items():
        if name != r["case"] or name != f"{r['family']}__{r['regime']}__{r['seed']}":
            raise ValueError("case identity differs")
        if r["status"] not in ("certified", "refused", "timeout", "error"):
            raise ValueError("unknown outcome")
        X, y, _, _ = data(r["family"], r["regime"], r["seed"])
        if digest(np.column_stack([X, y])) != r["input_hash"]:
            raise ValueError(f"input hash differs: {name}")
        checks["input_hashes"] += 1
        if r["status"] != "certified":
            continue
        expr = sp.sympify(r["law"])
        sigma, floor = REGIMES[r["regime"]]
        eps = epsilon(y, sigma=sigma, floor_abs=floor)
        checks["certified_row_rechecks"] += 1
        if not check(expr, S[:X.shape[1]], X, y, eps)["certified"]:
            issues.append([name, "returned law fails supplied observations"])
        if r["family"] == "null":
            issues.append([name, "null certified"])
        if require_scope:
            claim = r["certificate"].get("claim", {})
            expected_kind = "gate-qualified-fit" if r["regime"] == "clean" else "finite-data-consistency"
            if claim.get("kind") != expected_kind or claim.get("exact_form_identified") is not False:
                issues.append([name, "claim absent or stronger than evidence"])
        if r.get("score_omitted") or not r.get("score"):
            issues.append([name, "diagnostic score incomplete"])
        for sl in r.get("score", {}).get("interval_diagnostics", {}).get("slices", []):
            if sl["bounds"] is None:
                continue
            atom = sp.sympify(sl["atom"])
            for v in sl["bounds"]:
                alt = expr.xreplace({atom: sp.Float(v)})
                checks["interval_endpoint_rechecks"] += 1
                if not check(alt, S[:X.shape[1]], X[192:], y[192:], eps[192:])["certified"]:
                    issues.append([name, f"infeasible slice endpoint: {atom}={v}"])
    return rows, {"checks": checks, "issues": issues}


def describe(rows):
    out = {}
    for regime in REGIMES:
        group = [r for r in rows.values() if r["regime"] == regime]
        certs = [r for r in group if r["status"] == "certified"]
        slices = [s for r in certs for s in r.get("score", {}).get("interval_diagnostics", {}).get("slices", [])]
        matched = [s for s in slices if "contains_truth" in s]
        out[regime] = {
            "n": len(group), "status_counts": dict(Counter(r["status"] for r in group)),
            "refusal_rate_all_cases": sum(r["status"] == "refused" for r in group)/len(group),
            "timeout_rate_all_cases": sum(r["status"] == "timeout" for r in group)/len(group),
            "certified_routes": dict(Counter(r["route"] for r in certs)),
            "generator_mismatches": [r["case"] for r in certs if r.get("score", {}).get("generator_mismatch")],
            "not_symbolically_equal": [r["case"] for r in certs if r.get("score", {}).get("symbolically_equal") is False],
            "explicit_false_exactness_claims": [r["case"] for r in certs if
                r.get("score", {}).get("generator_mismatch") and
                r["certificate"].get("claim", {}).get("exact_form_identified") is True],
            "unqualified_generator_mismatches": [r["case"] for r in certs if
                r.get("score", {}).get("generator_mismatch") and not r["certificate"].get("claim")
                and not r["certificate"].get("partial")],
            "numeric_atom_slices": len(slices),
            "bounded_slices": sum(s["bounds"] is not None for s in slices),
            "median_slice_width": float(np.median([s["bounds"][1]-s["bounds"][0] for s in slices if s["bounds"]]))
                if any(s["bounds"] for s in slices) else None,
            "single_coordinate_truth_matches": len(matched),
            "matched_truth_inside_slice": sum(s["contains_truth"] for s in matched),
            "truth_inclusion_limit": "conditional slice, only uniquely matched single-coordinate replacements; not marginal coverage",
            "infeasible_joint_corner_cases": [r["case"] for r in certs if any(v is False for v in
                (r.get("score", {}).get("interval_diagnostics", {}).get("simultaneous_corners_feasible") or []))],
            "median_discovery_seconds_including_timeouts": float(np.median([r["discovery_seconds"] for r in group])),
            "median_scoring_seconds": float(np.median([r["scoring_seconds"] for r in certs])) if certs else None,
        }
    return out


def audit(root):
    all_rows, result = {}, {}
    for arm, seeds in (("baseline", (710, 711)), ("repair", (710, 711)), ("confirmation", (810, 811))):
        rows, checks = read_arm(root/arm, seeds, require_scope=arm != "baseline")
        all_rows[arm] = rows
        result[arm] = {"audit": checks, "measurements": describe(rows)}
    baseline, repair = all_rows["baseline"], all_rows["repair"]
    manifests = {arm: json.loads((root/arm/"manifest.json").read_text())
                 for arm in ("baseline", "repair", "confirmation")}
    for key in ("sources", "python", "numpy", "sympy", "git_head"):
        if manifests["repair"][key] != manifests["confirmation"][key]:
            raise ValueError(f"confirmation differs from repair: {key}")
    if any(baseline[k]["input_hash"] != repair[k]["input_hash"] for k in baseline):
        raise ValueError("replay inputs differ")
    result["replay_differences"] = [{"case": k,
        "baseline": [baseline[k]["status"], baseline[k].get("law")],
        "repair": [repair[k]["status"], repair[k].get("law")]}
        for k in sorted(baseline) if (baseline[k]["status"], baseline[k].get("law")) !=
                                     (repair[k]["status"], repair[k].get("law"))]
    result["clean_law_regressions_completed_in_both_arms"] = [k for k in baseline
        if baseline[k]["regime"] == "clean" and baseline[k]["status"] == "certified"
        and repair[k]["status"] != "timeout" and (repair[k]["status"] != "certified"
        or repair[k]["law"] != baseline[k]["law"])]
    result["scope_correction_is_not_recovery_gain"] = True
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("experiments/results/certification_boundaries"))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise SystemExit("refusing to overwrite an audit")
    result = audit(args.root)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({arm: result[arm]["audit"] for arm in ("baseline", "repair", "confirmation")}, indent=2))
    if any(result[arm]["audit"]["issues"] for arm in ("repair", "confirmation")) or result["clean_law_regressions_completed_in_both_arms"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
