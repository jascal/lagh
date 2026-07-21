"""C4: inner-scaled transcendental features T(b*u) with amplitude-product terms.

The inner scale b is judged by the residual of the FULL linear block (feature times
low-degree monomials AND ratio columns 1/x_i -- ratio amplitudes like exp(u)/x_0 fail
the gate without them), grid + golden-section refined, then snapped to the SMALLEST
rational the fit supports. STLSQ composes additively, so the products f, x_i*f,
x_i*x_j*f enter the library explicitly.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np

from ..base import Term, finite_guard, lstsq, snap_small

TIER = 4
B_GRID = np.concatenate([-np.logspace(-2, 2, 13), np.logspace(-2, 2, 13)])
TOP_K = 8
TOP_PAIR = 6   # feature x feature products among the best few


def _cores(dim: int):
    out = [(f"x_{i}", lambda X, i=i: X[:, i]) for i in range(dim)]
    for i, j in combinations(range(dim), 2):
        out.append((f"x_{i}*x_{j}", lambda X, i=i, j=j: X[:, i] * X[:, j]))
        out.append((f"x_{i}/x_{j}", lambda X, i=i, j=j: X[:, i] / X[:, j]))
        out.append((f"x_{j}/x_{i}", lambda X, i=i, j=j: X[:, j] / X[:, i]))
    return out


def _block(X: np.ndarray) -> np.ndarray:
    n, dim = X.shape
    cols = [np.ones(n)] + [X[:, i] for i in range(dim)]
    cols += [X[:, i] * X[:, j] for i, j in combinations(range(dim), 2)]
    with np.errstate(all="ignore"):
        for i in range(dim):
            inv = 1.0 / X[:, i]
            if np.all(np.isfinite(inv)):
                cols.append(inv)
    return np.column_stack(cols)


def _resid(feat: np.ndarray, M: np.ndarray, y: np.ndarray) -> float:
    A = np.column_stack([M * feat[:, None], M])
    c = lstsq(A, y)
    if c is None:
        return np.inf
    r = y - A @ c
    return float(r @ r)


def terms(dim: int, X_fit: np.ndarray, y_fit: np.ndarray, X_cert=None) -> list[Term]:
    X, y = X_fit, np.asarray(y_fit, float).ravel()
    M = _block(X)
    base = float(y @ y) + 1e-300
    found = []
    for T_name, T in (("exp", np.exp), ("sin", np.sin), ("cos", np.cos)):
        for u_name, u_fn in _cores(dim):
            with np.errstate(all="ignore"):
                u = u_fn(X)
            if not np.all(np.isfinite(u)):
                continue
            best_b, best_r = None, np.inf
            for b in B_GRID:
                with np.errstate(all="ignore"):
                    f = T(b * u)
                if not np.all(np.isfinite(f)) or np.std(f) < 1e-14:
                    continue
                r = _resid(f, M, y)
                if r < best_r:
                    best_b, best_r = float(b), r
            if best_b is None or best_r > 0.5 * base:
                continue
            lo, hi = sorted((best_b / 3, best_b * 3))
            for _ in range(40):
                m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
                with np.errstate(all="ignore"):
                    if _resid(T(m1 * u), M, y) < _resid(T(m2 * u), M, y):
                        hi = m2
                    else:
                        lo = m1
            b_ref = 0.5 * (lo + hi)
            with np.errstate(all="ignore"):
                if _resid(T(b_ref * u), M, y) < best_r:
                    best_b = b_ref
                    best_r = _resid(T(b_ref * u), M, y)
            bq = snap_small(best_b,
                            lambda bb: _resid(T(bb * u), M, y), best_r)
            bf = float(bq)
            name = f"{T_name}({bq.numerator}/{bq.denominator}*({u_name}))" \
                if bq.denominator != 1 else f"{T_name}({bq.numerator}*({u_name}))"
            fn = (lambda X_, T=T, bf=bf, u_fn=u_fn: T(bf * u_fn(X_)))
            found.append((best_r, T_name, Term(name, fn, finite_guard(fn), 3)))
    found.sort(key=lambda t: t[0])
    feats = [t for _, _, t in found[:TOP_K]]
    out = list(feats)
    for ft in feats:
        for i in range(dim):
            fn = (lambda X_, i=i, ft=ft: X_[:, i] * ft.fn(X_))
            out.append(Term(f"x_{i}*{ft.name}", fn, ft.guard, ft.complexity + 1))
        for i, j in combinations(range(dim), 2):
            fn = (lambda X_, i=i, j=j, ft=ft: X_[:, i] * X_[:, j] * ft.fn(X_))
            out.append(Term(f"x_{i}*x_{j}*{ft.name}", fn, ft.guard,
                            ft.complexity + 2))
    # FEATURE x FEATURE products: reaches damped oscillations exp(-b*t)*cos(w*t)
    # and the "envelope x carrier" family. The trap is that neither factor fits y
    # ALONE (each is modulated by the other), so a global top-K ranking misses
    # them. Guarantee coverage instead: take the best few features of EACH
    # transcendental family and pair ACROSS families, so exp x cos is always
    # formed regardless of single-feature residual.
    from collections import defaultdict
    by_family = defaultdict(list)
    # per-family cap scales DOWN with arity: the pair count is O(pool^2) and the
    # base library already grows with dim, so a fixed pool explodes discover on
    # 3-input modules (measured: snell hung). Keep coverage (>=2 per family) while
    # bounding the pool.
    per_family = max(2, TOP_PAIR - (dim - 1) * 2)
    for _, fam, t in found:
        if len(by_family[fam]) < per_family:
            by_family[fam].append(t)
    pool = [t for fam in by_family for t in by_family[fam]]
    seen = set()
    for a in range(len(pool)):
        for b in range(a, len(pool)):
            fa, fb = pool[a], pool[b]
            key = tuple(sorted((fa.name, fb.name)))
            if key in seen:
                continue
            seen.add(key)
            fn = (lambda X_, fa=fa, fb=fb: fa.fn(X_) * fb.fn(X_))
            out.append(Term(f"({fa.name})*({fb.name})", fn,
                            (lambda X_, fa=fa, fb=fb:
                             fa.guard(X_) and fb.guard(X_)),
                            fa.complexity + fb.complexity))
    return out
