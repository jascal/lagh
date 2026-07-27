"""The system discoverer (H3b, docs/LAWSYSTEMBENCH.md): certified discovery of
equation SYSTEMS -- roles unlabeled, shared constants, conserved quantities.

Three passes over a dict of named observable columns:
1. role/equation discovery: each column attempted as a certified function of
   the others (discover_passive, unchanged verdict machinery);
2. invariant discovery: nullspace of the centered term-library design matrix,
   sparsified, then CERTIFIED for constancy with minimality and alpha -- the
   inherently system-level product no per-equation run can emit;
3. cross-equation coherence: constants appearing in several certified
   equations must agree; the SystemCertificate conjoins everything with a
   union-bound alpha.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp

from .certify import epsilon, significance_log10
from .passive import discover_passive

_NULL_TOL = 1e-6          # singular-direction threshold (relative to largest)
_SPARSE_DROP = 0.02       # coefficient sparsification threshold (relative)


@dataclass
class SystemCertificate:
    equations: dict = field(default_factory=dict)     # target_col -> result dict
    invariants: list = field(default_factory=list)    # {expr, value, alpha_log10}
    roles: dict = field(default_factory=dict)         # col -> 'dependent'|'free'
    shared: list = field(default_factory=list)        # consistency verdicts
    alpha_log10_total: float | None = None


def _library(cols, data, n):
    """Registered invariant-search library over ALL columns: 1, x, x^2, x_i*x_j,
    ln x (guarded positive). Small on purpose -- |H| enters alpha."""
    names, mats = [], []
    for c in cols:
        v = data[c]
        names.append(c); mats.append(v)
        names.append(f"{c}**2"); mats.append(v ** 2)
        if np.all(v > 0):
            names.append(f"log({c})"); mats.append(np.log(v))
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            names.append(f"{a}*{b}"); mats.append(data[a] * data[b])
    return names, np.column_stack(mats)


def discover_invariants(data, *, sigma=0.0, max_invariants=10):
    """Exhaustive SMALL-SUBSET scan (supersedes global nullspace: the null
    space of a closed system is high-dimensional and an arbitrary basis vector
    is a many-term mixture, not the physical invariant -- measured 29-term
    output). For each term subset of size 2-4 (plus greedy size-5 extensions),
    the minimum-variance combination is checked for constancy and certified.
    Minimality is by construction (smallest subset first); each accepted
    invariant's span is excluded from later, larger candidates."""
    cols = sorted(data.keys())
    n = len(data[cols[0]])
    names, M = _library(cols, data, n)
    scale = np.sqrt(np.mean(M ** 2, axis=0)) + 1e-300
    Mn = (M - M.mean(0)) / scale
    T = M.shape[1]
    from itertools import combinations as _comb

    def combo(idx):
        sub = Mn[:, list(idx)]
        cov = sub.T @ sub / n
        try:
            w, V = np.linalg.eigh(cov)
        except np.linalg.LinAlgError:
            return None, np.inf
        return V[:, 0], float(max(w[0], 0.0) / (np.trace(cov) / len(idx) + 1e-300))

    accepted = []          # list of (set(idx), dict)
    cand_sizes = [list(_comb(range(T), k)) for k in (2, 3, 4)]
    near = []              # for greedy size-5
    for subsets in cand_sizes:
        for idx in subsets:
            iset = set(idx)
            # skip if it contains an accepted invariant's support (not minimal)
            if any(a[0] <= iset for a in accepted):
                continue
            v, ratio = combo(idx)
            if v is None:
                continue
            if ratio > 1e-14:
                if ratio < 1e-6:
                    near.append((ratio, idx))
                continue
            c = np.zeros(T)
            c[list(idx)] = v / scale[list(idx)]
            vals = M @ c
            mu = float(np.mean(vals))
            spread = float(np.max(np.abs(vals - mu)))
            vscale = float(np.mean(np.abs(vals))) + abs(mu) + 1e-300
            if spread > max(1e-7 * vscale, 4 * sigma * vscale):
                continue
            lead = c[list(idx)[0]]
            c = c / lead
            mu2 = float(np.mean(M @ c))
            expr = sum(sp.Float(ci) * sp.sympify(nm)
                       for ci, nm in zip(c, names) if ci != 0.0)
            alpha = significance_log10(expr, M @ c,
                                       epsilon(M @ c, sigma=sigma), T ** 4)
            accepted.append((iset, {"expr": str(expr), "value": mu2,
                                    "alpha_log10": alpha,
                                    "n_terms": len(idx)}))
            if len(accepted) >= max_invariants:
                return [a[1] for a in accepted]
    # greedy size-5: extend the best near-constant size-4 subsets
    near.sort()
    for _, idx in near[:200]:
        if len(idx) != 4:
            continue
        for extra in range(T):
            if extra in idx:
                continue
            iset = set(idx) | {extra}
            if any(a[0] <= iset for a in accepted):
                continue
            v, ratio = combo(tuple(sorted(iset)))
            if v is None or ratio > 1e-14:
                continue
            ids = tuple(sorted(iset))
            c = np.zeros(T)
            c[list(ids)] = v / scale[list(ids)]
            vals = M @ c
            mu = float(np.mean(vals))
            vscale = float(np.mean(np.abs(vals))) + abs(mu) + 1e-300
            if float(np.max(np.abs(vals - mu))) > max(1e-7 * vscale,
                                                      4 * sigma * vscale):
                continue
            c = c / c[ids[0]]
            expr = sum(sp.Float(ci) * sp.sympify(nm)
                       for ci, nm in zip(c, names) if ci != 0.0)
            alpha = significance_log10(expr, M @ c,
                                       epsilon(M @ c, sigma=sigma), T ** 5)
            accepted.append((iset, {"expr": str(expr),
                                    "value": float(np.mean(M @ c)),
                                    "alpha_log10": alpha, "n_terms": 5}))
            if len(accepted) >= max_invariants:
                return [a[1] for a in accepted]
    return [a[1] for a in accepted]


def _numeric_constants(expr_str):
    e = sp.sympify(expr_str)
    out = []
    for a in e.atoms(sp.Number):
        v = abs(float(a))
        if v not in (0.0, 1.0) and 1e-12 < v < 1e12:
            out.append(v)
    return out


def discover_system(data: dict, *, sigma: float = 0.0, seed: int = 0,
                    per_target_s: float = 60.0, max_tier: int = 4
                    ) -> SystemCertificate:
    """data: {column_name: 1-D array}, roles UNLABELED. Per-target wall-clock
    budget: an abstaining column must not starve the rest of the system pass
    (a system problem runs |columns| discoveries)."""
    import signal as _sig

    class _TO(Exception):
        pass

    def _h(s_, f_):
        raise _TO()

    cols = sorted(data.keys())
    n = len(data[cols[0]])
    cert = SystemCertificate()
    prev = _sig.getsignal(_sig.SIGALRM)
    # 1. per-column certified discovery
    for target in cols:
        others = [c for c in cols if c != target]
        X = np.column_stack([data[c] for c in others])
        try:
            _sig.signal(_sig.SIGALRM, _h)
            _sig.setitimer(_sig.ITIMER_REAL, per_target_s)
            r = discover_passive(X, np.asarray(data[target], float), sigma=sigma,
                                 seed=seed, n_resplits=1, max_tier=max_tier)
        except _TO:
            cert.roles[target] = "free"
            continue
        finally:
            _sig.setitimer(_sig.ITIMER_REAL, 0)
            _sig.signal(_sig.SIGALRM, prev)
        if r.certified:
            # rename x_i back to column names
            sub = {sp.Symbol(f"x_{i}"): sp.Symbol(others[i])
                   for i in range(len(others))}
            expr = sp.sympify(str(r.result.expr)).xreplace(sub)
            cert.equations[target] = {
                "expr": str(expr),
                "alpha_log10": r.result.certificate.alpha_log10}
            cert.roles[target] = "dependent"
        else:
            cert.roles[target] = "free"
    # 2. invariants (system-level)
    cert.invariants = discover_invariants(data, sigma=sigma)
    # 3. cross-equation shared-constant coherence
    consts = {t: _numeric_constants(e["expr"]) for t, e in cert.equations.items()}
    targets = list(consts)
    for i, a in enumerate(targets):
        for b in targets[i + 1:]:
            for va in consts[a]:
                for vb in consts[b]:
                    if abs(va - vb) <= 1e-6 * max(va, vb):
                        cert.shared.append(
                            {"between": [a, b], "value": va,
                             "consistent": True})
    alphas = [e["alpha_log10"] for e in cert.equations.values()
              if e["alpha_log10"] is not None]
    alphas += [iv["alpha_log10"] for iv in cert.invariants
               if iv["alpha_log10"] is not None]
    if alphas:
        # union bound in log space: log10(sum 10^a) ~ max + log10(count)
        cert.alpha_log10_total = float(max(alphas) + np.log10(len(alphas)))
    return cert
