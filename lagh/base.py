"""Shared datatypes: Term, Candidate, hardening (C2').

Extracted from wyly (docs/DISCOVERER.md P1). Coefficients are exact rationals with
bounded denominator; inner scales snap to SMALL rationals when the fit supports it --
a refined parameter carries optimizer error that can never certify at machine
precision, and a genuinely non-nice parameter should fail honestly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, Sequence

import numpy as np
import sympy as sp

MAX_DENOM = 10**6


@dataclass(frozen=True)
class Term:
    """One typed basis term: name, evaluator, domain guard, parsimony cost."""

    name: str
    fn: Callable[[np.ndarray], np.ndarray]
    guard: Callable[[np.ndarray], bool]
    complexity: int

    def __call__(self, X: np.ndarray) -> np.ndarray:
        return self.fn(X)

    def sympy(self) -> sp.Expr:
        return sp.sympify(self.name)


@dataclass
class Candidate:
    """One proposed law. Untrusted by construction; the checker decides."""

    expr: sp.Expr
    complexity: int
    channel: str
    val_residual: float = 0.0


ALWAYS = lambda X: True                                       # noqa: E731


def finite_guard(fn: Callable) -> Callable[[np.ndarray], bool]:
    def g(X: np.ndarray) -> bool:
        with np.errstate(all="ignore"):
            v = fn(X)
        return bool(np.all(np.isfinite(v)))
    return g


def snap(c: float, max_denom: int = MAX_DENOM) -> Fraction:
    """Smallest-cap rational that reproduces the fitted float. The 10^6 cap alone was
    a measured reach bug: a multi-term law carrying a small physical constant
    (G=6.674e-5 needs denominator 5e7) mis-snapped by ~1e-8 relative and failed
    certification -- so gravity/coulomb sums were unreachable by every snapped
    channel. Escalation triggers only when the small cap cannot reproduce the float
    to 1e-11 relative, so optimizer noise still snaps to nice rationals (1-3e-13
    -> 1); the checker (plus the coefficient gate on winners) still decides."""
    if not np.isfinite(c):
        raise ValueError(f"cannot snap non-finite coefficient {c}")
    f = Fraction(float(c))
    scale = max(abs(float(c)), 1e-300)
    for md in (max_denom, 10**9, 10**12):
        s = f.limit_denominator(md)
        if abs(float(s) - float(c)) <= 1e-11 * scale:
            return s
    return f            # exact dyadic representation of the float


def snap_all(coeffs: Sequence[float]) -> list[Fraction]:
    return [snap(c) for c in coeffs]


def snap_small(b: float, residual_of: Callable[[float], float],
               best_r: float, denominators=(1, 2, 3, 4, 6, 12)) -> Fraction:
    """Smallest-denominator rational the fit supports, else full precision."""
    for den in denominators:
        cand = Fraction(b).limit_denominator(den)
        if cand == 0:
            continue
        if residual_of(float(cand)) <= 1.05 * best_r + 1e-300:
            return cand
    return Fraction(b).limit_denominator(MAX_DENOM)


def to_expr(terms: Sequence[Term], coeffs: Sequence[Fraction]) -> sp.Expr:
    expr: sp.Expr = sp.Integer(0)
    for t, c in zip(terms, coeffs):
        if c == 0:
            continue
        expr = expr + sp.Rational(c.numerator, c.denominator) * t.sympy()
    return expr


def design_matrix(terms: Sequence[Term], X: np.ndarray) -> np.ndarray:
    with np.errstate(all="ignore"):
        return (np.column_stack([t(X) for t in terms])
                if terms else np.empty((len(X), 0)))


def admissible(terms: Sequence[Term], *Xs: np.ndarray) -> list[Term]:
    """Terms whose guard holds on EVERY provided point set (fit and certification:
    a law undefined anywhere in its own certification domain must never be proposed)."""
    keep = []
    for t in terms:
        try:
            if all(t.guard(X) for X in Xs):
                keep.append(t)
        except Exception:                                     # noqa: BLE001
            continue
    return keep


def lstsq(M: np.ndarray, y: np.ndarray):
    if M.shape[1] == 0:
        return None
    # A non-finite design matrix (NaN/Inf from e.g. arcsin out-of-domain or log<=0 on the
    # inverse-trig class) makes LAPACK's gelsd spam `DLASCL illegal value` on stderr and
    # grind -- guard BEFORE the call. Semantics are unchanged (such fits return None
    # anyway via the finiteness check below); this just short-circuits cheaply and quietly.
    if not (np.isfinite(M).all() and np.isfinite(y).all()):
        return None
    try:
        c, *_ = np.linalg.lstsq(M, y, rcond=None)
    except np.linalg.LinAlgError:
        return None
    return c if np.all(np.isfinite(c)) else None


def lambdify(syms, expr):
    return sp.lambdify(syms, expr, "numpy")


def eval_expr(expr, syms, X: np.ndarray) -> np.ndarray | None:
    """Numeric evaluation, never raising: None means 'undefined somewhere'."""
    if expr.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
        return None
    try:
        f = lambdify(syms, expr)
        with np.errstate(all="ignore"):
            v = np.broadcast_to(np.asarray(f(*X.T), float), (len(X),)).copy()
    except Exception:                                         # noqa: BLE001
        return None
    return v
