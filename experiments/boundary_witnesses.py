"""BND1 mechanism witnesses; controlled proposals are NOT discovery successes."""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import sympy as sp

from lagh.base import Candidate, eval_expr
from lagh.certify import check, epsilon, float_pinned, joint_pinned, parameter_interval
from lagh.engine import discover

S = list(sp.symbols("x_0:3"))


def prepass_fixture():
    rng = np.random.default_rng(713)
    X = rng.uniform(.5, 3., (240, 3))
    X[:, 1] = X[:, 0] + rng.uniform(-1e-7, 1e-7, len(X))
    expr = sp.Float(np.pi)*S[0] + sp.Float(np.sqrt(2))*S[1]
    y = eval_expr(expr, S, X)
    return X, y, expr


def prepass_result(gate="joint_quotient"):
    X, y, expr = prepass_fixture()
    with ExitStack() as ctx:
        ctx.enter_context(patch("lagh.classes.c3_powerlaw.candidates",
                                side_effect=lambda _: [Candidate(expr, 3, "controlled-witness")]))
        ctx.enter_context(patch("lagh.classes.c8_angular.candidates", return_value=[]))
        ctx.enter_context(patch("lagh.classes.c9_genmonomial.candidates", return_value=[]))
        ctx.enter_context(patch("lagh.classes.c5_transforms.transforms", return_value=[]))
        ctx.enter_context(patch("lagh.engine._tier_candidates", return_value=[]))
        return discover(X[:144], y[:144], X[144:192], y[144:192],
                        X[192:], y[192:], max_tier=3, coefficient_gate=gate)


def interval_fixture():
    rng = np.random.default_rng(712)
    X = rng.uniform(.5, 3., (80, 2))
    X[:, 1] = X[:, 0] + rng.uniform(-.01, .01, len(X))
    a, b = sp.Float(1.3), sp.Float(.7)
    expr = a*S[0]+b*S[1]
    y = eval_expr(expr, S[:2], X)
    eps = np.full(len(X), .01)
    return X, y, eps, expr, (a, b)


def measure():
    X, y, eps, expr, (a, b) = interval_fixture()
    ia = parameter_interval(expr, S[:2], X, y, eps, a)
    ib = parameter_interval(expr, S[:2], X, y, eps, b)
    alternative = sp.Float(1.4)*S[0]+sp.Float(.6)*S[1]
    corner = expr.xreplace({a: sp.Float(ia[1]), b: sp.Float(ib[1])})
    Xp, yp, ep = prepass_fixture()
    xp, yc = Xp[192:], yp[192:]
    marginal, snapped = float_pinned(ep, S, xp, yc, epsilon(yc))
    results = {gate: prepass_result(gate) for gate in ("marginal", "joint_quotient")}
    z = np.linspace(.5, 3., 80)[:, None]
    yz = 1.3*z[:, 0]
    return {
        "coordinate_slices": {
            "intervals": [ia, ib],
            "each_upper_endpoint_fits": [check(expr.xreplace({atom: sp.Float(iv[1])}),
                S[:2], X, y, eps)["certified"] for atom, iv in ((a, ia), (b, ib))],
            "simultaneous_upper_corner_fits": check(corner, S[:2], X, y, eps)["certified"],
            "alternative_vector_fits": check(alternative, S[:2], X, y, eps)["certified"],
            "alternative_coefficients_inside_slices": [ia[0] <= 1.4 <= ia[1], ib[0] <= .6 <= ib[1]],
        },
        "controlled_prepass": {
            "marginal_gate_passes": marginal,
            "joint_gate_passes_snapped_candidate": joint_pinned(snapped, S, xp, yc, epsilon(yc)),
            "arms": {g: {"certified": r.certificate.certified, "law": str(r.expr),
                         "notes": r.certificate.notes} for g, r in results.items()},
        },
        "indistinguishable_coefficients": {
            regime: [check(sp.Float(c)*S[0], S[:1], z, yz, band)["certified"]
                     for c in (1.3, 1.300001)]
            for regime, band in (("noise", epsilon(yz, sigma=1e-3)),
                                  ("floor", epsilon(yz, floor_abs=1e-3)))},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise SystemExit("refusing to overwrite witness evidence")
    result = measure()
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
