"""Systems of PDEs over one weak-form row set (docs/DIRECTION_PDE_SYSTEMS.md,
registration docs/CASE_STUDY_PDE_C3.md).

A coupled system is NOT n independent discovery problems. What makes it one
problem is that every equation is certified over the SAME patch rows, from the
same fields, under one declared error model -- and that the conjoined claim's
alpha is a union bound over the equations, dominated by the weakest of them.

Everything the scalar arc built transfers unchanged: the by-parts identity does
not care whether g is u**2/2 or u*v; with independent per-field noise the
sensitivity vectors concatenate, so the Gram is block-diagonal over fields and
the band is still one quadratic form a'Ga (weakform.PatchEpsilon); the solution
holdout, multi-scale pooling, row normalization and interval parameters are all
field-agnostic. This module is the part a per-equation loop cannot do:

  * one row set, one target per equation, features spanning every field,
  * time-derivative columns are TARGETS and never features (a system of
    evolution equations; using another equation's u_t as an input would be a
    different -- differential-algebraic -- claim),
  * parameter intervals measured per equation on its own certification rows,
  * a conjoined SystemCertificate: cross-equation constant agreement plus the
    union-bound alpha (lagh.systems).

The library growth is the thing to watch: cross terms multiply fast and |H|
enters alpha directly, so the vocabulary is a REGISTERED list per curriculum
stage, never a generated cross-product. That discipline lives in the runner's
declaration, and this module reports the vocabulary it was handed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp

from .certify import band, free_atoms, parameter_interval
from .engine import discover
from .systems import SystemCertificate, union_alpha_log10
from .weakform import LIBRARY, PatchEpsilon, build_nd


@dataclass
class SystemRows:
    """Weak-form rows pooled over solutions, for a system of fields."""
    A: np.ndarray                 # (n_rows, n_terms)
    names: list                   # term names, column order of A
    terms: list                   # the Term objects (vocabulary, as registered)
    det: np.ndarray               # (n_rows, n_terms) deterministic error
    gram: np.ndarray              # (n_rows, n_terms, n_terms) noise Gram
    sol: np.ndarray               # (n_rows,) which solution each row came from
    rejected: int = 0
    ndim: int = 2
    fields: tuple = ()
    field_err: float = 0.0

    @property
    def n_solutions(self) -> int:
        return int(len(np.unique(self.sol)))

    def targets(self) -> list:
        """The evolution targets: every term carrying a time derivative."""
        return [tm.name for tm in self.terms if tm.multi(self.ndim)[-1] > 0]

    def features(self) -> list:
        """Everything that is not a target -- the candidate vocabulary."""
        return [tm.name for tm in self.terms if tm.multi(self.ndim)[-1] == 0]


def assemble(solutions, terms, patch_fn, *, sigma: float = 0.0,
             field_err: float = 0.0, p: int = 16,
             normalize_by: str | None = "1") -> SystemRows | None:
    """Pool weak-form rows over several SOLUTIONS of the same system.

    `solutions` is a list of (fields, coords): fields is {name: array}, coords is
    the tuple of axis coordinate vectors (space..., time last). `patch_fn(coords)`
    returns that solution's patch family -- one family per solution, because
    different solutions may live on different grids.

    A field produced by a reference solver carries the solver's own declared
    error; pass it as `field_err` and it enters deterministically through the L1
    sensitivity (weakform.WeakSystem.det). Declared, never assumed away: stage 1
    of the curriculum is exactly solvable so that this term is provably zero, and
    every later stage has to state it.
    """
    A, det, gram, sol, rejected = [], [], [], [], 0
    names = None
    for i, (fields, coords) in enumerate(solutions):
        s = build_nd(fields, coords, terms, patch_fn(coords), p=p, sigma=sigma)
        rejected += s.rejected
        if len(s.A) == 0:
            continue
        d = s.det(field_err)
        if normalize_by is not None:
            j = s.names.index(normalize_by)
            w = s.A[:, j]
            keep = [k for k in range(len(s.names)) if k != j]
            A.append((s.A / w[:, None])[:, keep])
            det.append((d / np.abs(w)[:, None])[:, keep])
            gram.append((s.gram / (w ** 2)[:, None, None])[:, keep][:, :, keep])
            names = [s.names[k] for k in keep]
        else:
            A.append(s.A)
            det.append(d)
            gram.append(s.gram)
            names = list(s.names)
        sol.append(np.full(len(s.A), i))
    if not A:
        return None
    T = [LIBRARY[tm] if isinstance(tm, str) else tm for tm in terms]
    keep_terms = [tm for tm in T if tm.name in names]
    return SystemRows(A=np.vstack(A), names=names, terms=keep_terms,
                      det=np.vstack(det), gram=np.concatenate(gram),
                      sol=np.concatenate(sol), rejected=rejected,
                      ndim=len(solutions[0][1]),
                      fields=tuple(sorted(solutions[0][0])),
                      field_err=field_err)


def linear_coefficients(expr, feat) -> tuple:
    """(is_linear, {term: coefficient}) reading a certified law in the term
    vocabulary. A constant term is reported under the name '1'."""
    if expr is None:
        return False, {}
    e = sp.expand(sp.sympify(expr))
    syms = [sp.Symbol(f"x_{i}") for i in range(len(feat))]
    got, rest = {}, e
    for i, s in enumerate(syms):
        c = e.coeff(s)
        if c.free_symbols:
            return False, {}
        if float(c) != 0.0:
            got[feat[i]] = float(c)
        rest = rest - c * s
    rest = sp.simplify(rest)
    if rest.free_symbols:
        return False, {}
    if float(rest) != 0.0:
        got["1"] = float(rest)
    return True, got


def intervals_for(expr, syms, X, y, eps, coeffs, feat) -> dict:
    """{term: (lo, hi)} -- what the declared band actually determines about each
    coefficient, measured by bisection on the certification predicate itself
    (certify.parameter_interval), not by a covariance.

    An atom the search cannot bound within `max_rel` is reported as UNDETERMINED
    (None) rather than as its point value; the runner must not quietly turn that
    into a number."""
    out = {}
    atoms = free_atoms(expr)
    for a in atoms:
        iv = parameter_interval(expr, syms, X, y, eps, a)
        v = float(a)
        # attach the atom to the term whose coefficient it is
        for nm, c in coeffs.items():
            if nm in out:
                continue
            if abs(c - v) <= 1e-12 * max(1.0, abs(v)) or \
                    abs(abs(c) - abs(v)) <= 1e-12 * max(1.0, abs(v)):
                s = 1.0 if abs(c - v) <= abs(c + v) else -1.0
                out[nm] = None if iv is None else (
                    (iv[0], iv[1]) if s > 0 else (-iv[1], -iv[0]))
                break
    for nm, c in coeffs.items():
        out.setdefault(nm, (c, c))          # integers/±1: exact by construction
    return out


def truth_check(rows: SystemRows, target: str, truth: dict, *,
                sigma: float = 0.0, coeff_max: float = 2.0) -> dict:
    """Does the TRUE law sit inside its own declared band on these rows?

    Run this BEFORE reading any abstain as a finding. Measured the hard way in
    the system scoping probe: the first run abstained structurally at every
    sigma including 0, which looks exactly like an identifiability result -- and
    was a construction bug (independent random phases for u and v made the
    coupling term spatially orthogonal to u, so the system did not hold at all;
    the truth missed its own band by 0.85 against a band of 4e-10). A null
    result that is really a construction bug is indistinguishable from a finding
    unless the truth is checked against its own band first.

    Returns the worst |residual| / eps over the rows; <= 1 means the fields and
    the vocabulary really do satisfy the claimed system.
    """
    feat = [n for n in rows.names if n != target]
    j = rows.names.index(target)
    cols = [rows.names.index(n) for n in feat]
    X, y = rows.A[:, cols], rows.A[:, j]
    m = PatchEpsilon(rows.names, target, y, X, rows.det, rows.gram,
                     sigma=sigma, floor_abs=0.0, coeff_max=coeff_max,
                     feat_names=feat)
    syms = [sp.Symbol(f"x_{i}") for i in range(len(feat))]
    expr = sp.S.Zero
    for nm, c in truth.items():
        expr = expr + (sp.Float(c) if nm == "1"
                       else sp.Float(c) * syms[feat.index(nm)])
    pred = np.zeros(len(y))
    for nm, c in truth.items():
        if nm != "1":
            pred = pred + c * X[:, feat.index(nm)]
        else:
            pred = pred + c
    eps = band(m, expr)
    ratio = np.abs(pred - y) / eps
    return {"truth_max_ratio": float(np.max(ratio)),
            "truth_median_ratio": float(np.median(ratio)),
            "truth_certifies": bool(np.max(ratio) <= 1.0),
            "median_band": float(np.median(eps)),
            "median_target": float(np.median(np.abs(y)))}


def agreement(rows: SystemRows, target: str, got: dict, truth: dict, *,
              sigma: float = 0.0, coeff_max: float = 2.0) -> dict:
    """Do the certified law and the true law AGREE where the certificate applies?

    Support equality is the wrong test on its own. With spectrally poor initial
    data a library column can be an exact multiple of another ACROSS THE WHOLE
    solution family -- single-mode data makes u_xx = -u identically, so the
    solution holdout cannot break it -- and the honest verdict there is a
    domain-restricted certificate, not a wrong one: the certified law and the
    truth are the same function everywhere the domain claim applies. What would
    be wrong is a law that DISAGREES with the truth on the certified rows by
    more than the declared band. These are different questions and this returns
    both."""
    feat = [n for n in rows.names if n != target]
    j = rows.names.index(target)
    cols = [rows.names.index(n) for n in feat]
    X, y = rows.A[:, cols], rows.A[:, j]
    m = PatchEpsilon(rows.names, target, y, X, rows.det, rows.gram,
                     sigma=sigma, floor_abs=0.0, coeff_max=coeff_max,
                     feat_names=feat)
    syms = [sp.Symbol(f"x_{i}") for i in range(len(feat))]

    def evaluate(c):
        v = np.zeros(len(y))
        e = sp.S.Zero
        for nm, co in c.items():
            if nm == "1":
                v = v + co
                e = e + sp.Float(co)
            else:
                v = v + co * X[:, feat.index(nm)]
                e = e + sp.Float(co) * syms[feat.index(nm)]
        return v, e

    vg, eg = evaluate(got)
    vt, _ = evaluate(truth)
    eps = band(m, eg)
    d = np.abs(vg - vt) / eps
    return {"max_disagreement_over_band": float(np.max(d)),
            "agrees_on_certified_domain": bool(np.max(d) <= 1.0)}


def discover_equation(rows: SystemRows, target: str, *, sigma: float = 0.0,
                      seed: int = 0, max_tier: int = 3,
                      features: list | None = None,
                      holdout: bool = True, coeff_max: float = 2.0) -> dict:
    """One equation of the system: `target` as a certified function of the
    feature columns, certified on rows from a HELD-OUT SOLUTION.

    The holdout is not a convenience. On one solution the field satisfies
    on-shell relations (a traveling wave has u_t = -c u_x identically), so rival
    laws differing by a multiple of an on-shell-zero quantity are
    indistinguishable -- measured in C0, where a row split let the KdV soliton
    certify u_t = -u_x. For a SYSTEM the exposure is larger, not smaller: a
    coupled pair on one trajectory can satisfy relations holding in both
    equations at once.
    """
    feat = list(features) if features is not None else rows.features()
    if target in feat:
        feat.remove(target)
    j = rows.names.index(target)
    cols = [rows.names.index(n) for n in feat]
    X, y = rows.A[:, cols], rows.A[:, j]
    m = PatchEpsilon(rows.names, target, y, X, rows.det, rows.gram,
                     sigma=sigma, floor_abs=0.0, coeff_max=coeff_max,
                     feat_names=feat)
    n_sol = rows.n_solutions
    if holdout and n_sol < 2:
        return {"target": target, "certified": False,
                "abstain": "single-solution", "n_rows": int(len(y)),
                "features": feat,
                "notes": ["a PDE claim needs patches from a HELD-OUT solution; "
                          "one solution supports only an on-shell statement "
                          "about that solution"]}
    rng = np.random.default_rng(seed)
    if holdout:
        held = np.unique(rows.sol)[-1]
        tr = rng.permutation(np.where(rows.sol != held)[0])
        ce = np.where(rows.sol == held)[0]
    else:
        idx = rng.permutation(len(y))
        b = int(0.8 * len(y))
        tr, ce = idx[:b], idx[b:]
    a = int(0.75 * len(tr))
    r = discover(X[tr[:a]], y[tr[:a]], X[tr[a:]], y[tr[a:]], X[ce], y[ce],
                 sigma=sigma, eps_model=m.subset(ce), max_tier=max_tier,
                 declared_basis=True, linear_basis=True,
                 band_sel=m.subset(tr[a:])(None))
    c = r.certificate
    out = {"target": target, "certified": bool(c.certified),
           "abstain": c.abstain, "alpha_log10": c.alpha_log10,
           "n_rows": int(len(y)), "n_cert_rows": int(len(ce)),
           "n_solutions": n_sol, "features": feat, "tier": r.tier,
           "notes": [str(n)[:220] for n in c.notes][:3]}
    if not c.certified:
        return out
    syms = [sp.Symbol(f"x_{i}") for i in range(len(feat))]
    linear, coeffs = linear_coefficients(r.expr, feat)
    eps_ce = m.subset(ce)
    out["law"] = _readable(str(r.expr), feat)
    out["expr"] = str(r.expr)
    out["linear_in_columns"] = bool(linear)
    out["coefficients"] = coeffs
    out["coeff_max_declared"] = coeff_max
    out["coeff_max_certified"] = max([abs(v) for v in coeffs.values()] or [0.0])
    out["intervals"] = {k: (None if v is None else [float(v[0]), float(v[1])])
                        for k, v in intervals_for(r.expr, syms, X[ce], y[ce],
                                                  eps_ce, coeffs, feat).items()}
    out["median_signal_to_band"] = float(np.median(
        np.abs(y[ce]) / band(eps_ce, r.expr)))
    return out


def _readable(law: str, feat) -> str:
    for i, nm in sorted(enumerate(feat), key=lambda kv: -kv[0]):
        law = law.replace(f"x_{i}", f"[{nm}]")
    return law


def conjoin(equations: list, *, invariants=None,
            rel_tol: float = 1e-6) -> SystemCertificate:
    """The system-level certificate: the equations conjoined, their shared
    constants checked for agreement, and ONE union-bound alpha.

    The union bound is the honest conjunction: the chance that ANY of the
    conjoined claims is a coincidence is at most the sum of the per-claim
    bounds, so the total is dominated by the WEAKEST equation. Reporting the
    strongest equation's alpha for the system, or a product of them, would be
    claiming the conjunction is stronger than its parts.
    """
    cert = SystemCertificate()
    for e in equations:
        if e.get("certified"):
            cert.equations[e["target"]] = {
                "expr": e.get("law", e.get("expr", "")),
                "alpha_log10": e.get("alpha_log10")}
            cert.roles[e["target"]] = "dependent"
        else:
            cert.roles[e["target"]] = "free"
    cert.invariants = list(invariants or [])
    consts = {e["target"]: [abs(v) for v in (e.get("coefficients") or {}).values()
                            if abs(v) not in (0.0, 1.0)]
              for e in equations if e.get("certified")}
    tg = list(consts)
    for i, a in enumerate(tg):
        for b in tg[i + 1:]:
            for va in consts[a]:
                for vb in consts[b]:
                    if abs(va - vb) <= rel_tol * max(va, vb):
                        cert.shared.append({"between": [a, b], "value": va,
                                            "consistent": True})
    alphas = [e["alpha_log10"] for e in equations
              if e.get("certified") and e.get("alpha_log10") is not None]
    alphas += [iv["alpha_log10"] for iv in cert.invariants
               if iv.get("alpha_log10") is not None]
    cert.alpha_log10_total = union_alpha_log10(alphas)
    return cert


def weakest(equations: list):
    """(target, alpha) of the equation that dominates the union bound -- the
    quantity Y2 says must be reported instead of the strongest one."""
    scored = [(e["alpha_log10"], e["target"]) for e in equations
              if e.get("certified") and e.get("alpha_log10") is not None]
    if not scored:
        return None, None
    a, t = max(scored)
    return t, a
