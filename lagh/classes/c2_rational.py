"""C2: rational laws P/Q without nonlinear optimisation (implicit-SINDy), plus the
pure-term-denominator pass.

Both normalisation gaps are lessons: Q's constant normalised to 1 makes constant-free
denominators (y = P/x_1) unreachable, so a second pass solves y*d = P for single
terms d. Both passes are STLSQ-guided with a numeric prefilter before any symbolic
construction -- prefiltering which candidates to CHECK cannot weaken the checker.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import sympy as sp

from ..base import Candidate, design_matrix, lstsq, snap_all, to_expr
from ..engine_util import stlsq_supports

TIER = 2
PREFILTER_REL = 1e-6


def candidates(ctx) -> list[Candidate]:
    terms, X_fit, y_fit = ctx.terms, ctx.X_fit, ctx.y_fit
    X_sel, y_sel = ctx.X_sel, ctx.y_sel
    M_tr, M_va = design_matrix(terms, X_fit), design_matrix(terms, X_sel)
    y_tr = np.asarray(y_fit, float).ravel()
    y_va = np.asarray(y_sel, float).ravel()
    yscale = float(np.sqrt(np.mean(y_va**2))) + 1e-300
    out: list[Candidate] = []

    denom_idx = [i for i, t in enumerate(terms)
                 if t.name != "1" and ("*" not in t.name or "(" in t.name)]

    # pure-term denominators: y*d = P, plain lstsq, STLSQ-guided numerators
    for j in denom_idx:
        d_tr, d_va = M_tr[:, j], M_va[:, j]
        if np.abs(d_tr).min() < 1e-9 or np.abs(d_va).min() < 1e-9:
            continue
        target = y_tr * d_tr
        sups = stlsq_supports(M_tr, target) | {(k,) for k in range(len(terms))}
        # STLSQ cannot isolate sparse pairs in a heavily collinear feature library
        # (measured: it returned one size-52 support while the exact pair passed
        # every gate). Add targeted pairs: every {term, constant}, plus all pairs
        # among the best dozen singles by fit residual.
        one = next((i for i, t in enumerate(terms) if t.name == "1"), None)
        if one is not None:
            sups |= {(k, one) for k in range(len(terms)) if k != one}
        singles = []
        for k in range(len(terms)):
            c1 = lstsq(M_tr[:, [k]], target)
            if c1 is None:
                continue
            r = target - M_tr[:, [k]] @ c1
            singles.append((float(r @ r), k))
        singles.sort()
        top = [k for _, k in singles[:12]]
        sups |= {tuple(sorted((a, b))) for ai, a in enumerate(top)
                 for b in top[ai + 1:]}
        for sup in sups:
            if len(sup) > 4:
                continue
            cols = list(sup)
            c = lstsq(M_tr[:, cols], target)
            if c is None:
                continue
            pred = (M_va[:, cols] @ c) / d_va
            if not np.all(np.isfinite(pred)):
                continue
            vr = float(np.sqrt(np.mean((pred - y_va) ** 2)))
            if vr > PREFILTER_REL * yscale:
                continue
            expr = to_expr([terms[i] for i in cols], snap_all(c)) / terms[j].sympy()
            out.append(Candidate(expr=expr, complexity=int(sp.count_ops(expr)),
                                 channel="c2-pure", val_residual=vr))

    # implicit-SINDy: y*(1 + sum q T) = sum p S -- linear in [p, q].
    # The budget counts ATTEMPTS, not constructed candidates: with a prefilter that
    # gates construction, a constructed-count budget never triggers and enumeration
    # runs exhaustive over the feature-enlarged library (measured: suite timeout).
    # Exhaustive numerators are limited to size <= 2; size-3 comes from STLSQ.
    attempts = 0
    MAX_ATTEMPTS = 40000
    p_pool3 = stlsq_supports(M_tr, y_tr)
    for q_size in (1, 2):
        for q_sup in combinations(denom_idx, q_size):
            Tq_tr = M_tr[:, list(q_sup)] * y_tr[:, None]
            Tq_va = M_va[:, list(q_sup)] * y_va[:, None]
            p_iter = list(combinations(range(len(terms)), 1)) + \
                     list(combinations(range(len(terms)), 2)) + \
                     [s_ for s_ in p_pool3 if len(s_) == 3]
            for p_sup in p_iter:
                    p_size = len(p_sup)
                    attempts += 1
                    if attempts > MAX_ATTEMPTS:
                        return out
                    A = np.column_stack([M_tr[:, list(p_sup)], -Tq_tr])
                    sol = lstsq(A, y_tr)
                    if sol is None:
                        continue
                    p_c, q_c = sol[:p_size], sol[p_size:]
                    Q_va = 1.0 + M_va[:, list(q_sup)] @ q_c
                    if np.abs(Q_va).min() < 1e-3:
                        continue
                    pred = (M_va[:, list(p_sup)] @ p_c) / Q_va
                    if not np.all(np.isfinite(pred)):
                        continue
                    vr = float(np.sqrt(np.mean((pred - y_va) ** 2)))
                    if vr > PREFILTER_REL * yscale:
                        continue
                    num = to_expr([terms[i] for i in p_sup], snap_all(p_c))
                    den = sp.Integer(1) + to_expr([terms[i] for i in q_sup],
                                                  snap_all(q_c))
                    if den == 0:
                        continue
                    expr = num / den
                    out.append(Candidate(expr=expr,
                                         complexity=int(sp.count_ops(expr)),
                                         channel="c2-implicit", val_residual=vr))
    return out
