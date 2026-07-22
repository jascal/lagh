"""Core logic for the three lagh tools (docs/DIRECTION_TOOLSHAPE.md).

Pure Python, no MCP-SDK dependency, so it is testable in-process. `server.py` is a
thin wrapper that registers these over the MCP transport.

The three acts, named for what they DO (Orca convention -- never `certify`):

    recover(X, y)         bounded  -- discover an exact law     -> Certificate | Abstention
    verify(X, y, form)    bounded  -- check a DECLARED form      -> Certificate | Abstention
    fit(X, y)             UNBOUNDED-- best-guess + diagnosis      -> Conjecture[] + Diagnosis

The wall is structural: `fit`'s result has NO `certified` field -- a guarantee cannot
be read off it by accident. `recover`/`verify` return a `certified` bool and a
`strength` of `pinned` (rational, no rival within the noise) or `consistent` (a declared
irrational fits, but the constant is not identifiable). See the doc for the full rationale.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import sympy as sp

from ..base import eval_expr, lstsq
from ..certify import (Abstain, coherent, epsilon, pinned, sample_box)
from ..engine import discover

# nameable constants a free-fit exponent might be reaching for (fit's diagnosis)
_NAMED = [("e", float(np.e)), ("pi", float(np.pi)), ("sqrt2", float(np.sqrt(2))),
          ("1/e", 1.0 / float(np.e)), ("ln2", float(np.log(2))),
          ("golden", (1 + 5 ** 0.5) / 2)]


def _syms(dim: int) -> list[sp.Symbol]:
    return [sp.Symbol(f"x_{i}") for i in range(dim)]


def _prep(X, y):
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    y = np.asarray(y, float).ravel()
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    return X[m], y[m]


def _split(X, y, seed=0):
    n = len(X)
    i = np.random.default_rng(seed).permutation(n)
    a, b = int(0.6 * n), int(0.8 * n)
    return (X[i[:a]], y[i[:a]], X[i[a:b]], y[i[a:b]], X[i[b:]], y[i[b:]])


def _has_irrational(expr) -> bool:
    """A declared irrational the data cannot pin: e/pi/golden constants, or a Float
    that is not a clean bounded rational."""
    if expr.has(sp.E, sp.pi, sp.GoldenRatio, sp.EulerGamma):
        return True
    for f in expr.atoms(sp.Float):
        if Fraction(float(f)).limit_denominator(10 ** 6) != Fraction(float(f)):
            # a Float that only approximates -- treat as un-pinnable
            if abs(float(f) - round(float(f))) > 1e-9:
                return True
    return False


def _strength(expr, syms, X_cert, y_cert, eps, sigma) -> str:
    """`consistent` if the form carries a declared irrational (never identifiable);
    otherwise `pinned` (rational structure, and -- under noise -- no neighbour rival)."""
    if _has_irrational(expr):
        return "consistent"
    P = sample_box(X_cert, extend=0.5)
    yscale = float(np.sqrt(np.mean(y_cert ** 2))) + 1e-300
    return "pinned" if pinned(expr, syms, X_cert, y_cert, eps, P, yscale, sigma) \
        else "not-pinned"


# --------------------------------------------------------------------------- recover

def recover(X, y, *, sigma: float = 0.0, max_tier: int = 7) -> dict:
    """Bounded. Discover an exact law from (X, y). Returns a Certificate or an
    Abstention -- never a bare guess."""
    X, y = _prep(X, y)
    if len(X) < 8:
        return {"tag": "open", "tool": "recover", "certified": False,
                "abstain": Abstain.RANGE.value,
                "note": f"only {len(X)} finite points; too thin to certify"}
    dim = X.shape[1]
    r = discover(*_split(X, y), sigma=float(sigma), max_tier=max_tier)
    c = r.certificate
    if not c.certified:
        return {"tag": "open", "tool": "recover", "certified": False,
                "abstain": c.abstain, "domain_size": c.domain_size,
                "note": "; ".join(map(str, c.notes)) if c.notes else ""}
    syms = _syms(dim)
    eps = epsilon(y, sigma=float(sigma))
    return {"tag": "proved", "tool": "recover", "certified": True,
            "law": str(r.expr), "strength": _strength(r.expr, syms, X, y, eps, sigma),
            "domain_size": c.domain_size, "tier": r.tier,
            "bounds": [[float(X[:, j].min()), float(X[:, j].max())] for j in range(dim)],
            "note": "certified over the stated finite domain, not proved for the world"}


# ---------------------------------------------------------------------------- verify

def verify(X, y, form: str, *, sigma: float = 0.0) -> dict:
    """Bounded. Check a caller-DECLARED form. The form is a sympy expression in
    x_0..x_{d-1}; its single overall scale is refit, then it is checked over the
    domain. A rational form can certify `pinned`; a declared irrational only
    `consistent`."""
    X, y = _prep(X, y)
    dim = X.shape[1]
    syms = _syms(dim)
    try:
        expr = sp.sympify(form)
    except Exception as e:                                     # noqa: BLE001
        return {"tag": "open", "tool": "verify", "certified": False,
                "abstain": "malformed-form", "note": str(e)[:200]}
    if len(X) < 8:
        return {"tag": "open", "tool": "verify", "certified": False,
                "abstain": Abstain.RANGE.value,
                "note": f"only {len(X)} finite points; too thin to certify"}
    Xf, yf, Xs, ys, Xc, yc = _split(X, y)
    base = eval_expr(expr, syms, Xf)
    if base is None or not np.all(np.isfinite(base)):
        return {"tag": "open", "tool": "verify", "certified": False,
                "abstain": Abstain.NUMERICAL.value,
                "note": "declared form does not evaluate finitely on the domain"}
    d2 = float(np.dot(base, base))
    alpha = float(np.dot(base, yf) / d2) if d2 > 0 else 1.0
    # snap the refit scale to a rational when clean (keeps the exact claim exact)
    a_snap = Fraction(alpha).limit_denominator(10 ** 6)
    scaled = sp.Rational(a_snap.numerator, a_snap.denominator) * expr
    eps = epsilon(yc, sigma=float(sigma))
    pred = eval_expr(scaled, syms, Xc)
    if pred is None or not np.all(np.isfinite(pred)):
        return {"tag": "open", "tool": "verify", "certified": False,
                "abstain": Abstain.NUMERICAL.value, "note": "form diverges on cert split"}
    miss = int(np.sum(np.abs(pred - yc) > eps))
    if miss:
        return {"tag": "open", "tool": "verify", "certified": False,
                "abstain": Abstain.STRUCTURAL.value,
                "note": f"declared form refuted: {miss}/{len(yc)} points exceed eps"}
    strength = _strength(scaled, syms, Xc, yc, eps, sigma)
    if strength == "not-pinned":
        return {"tag": "open", "tool": "verify", "certified": False,
                "abstain": Abstain.PARAMETRIC.value,
                "note": "form fits but a neighbour-rational fits within the noise too"}
    return {"tag": "proved", "tool": "verify", "certified": True,
            "law": str(scaled), "strength": strength, "domain_size": len(Xc),
            "note": ("consistent: fits within eps, but the irrational constant is not "
                     "identifiable from the data" if strength == "consistent"
                     else "pinned: this exact form, no rival within the noise")}


# ------------------------------------------------------------------------------- fit

def fit(X, y, *, sigma: float = 0.0, top: int = 5) -> dict:
    """UNBOUNDED scout. Best-guess conjectures + an identifiability diagnosis + a
    next_action pointer. NO `certified` field -- a guarantee cannot be read off this.
    Its real product is the diagnosis (what pins, what does not, what it would take)."""
    X, y = _prep(X, y)
    dim = X.shape[1]
    syms = _syms(dim)
    out: dict = {"tag": "exploratory", "tool": "fit",
                 "note": "conjectures, NOT certificates; to certify call recover/verify"}
    conj: list[dict] = []
    diagnosis = {"kind": "unknown", "detail": ""}
    next_action = "recover"

    # 1) free-exponent power-law probe (surfaces continuum / irrational exponents)
    if np.all(X > 0) and np.all(y > 0):
        L = np.column_stack([np.ones(len(X)), np.log(X)])
        c = lstsq(L, np.log(y))
        if c is not None:
            exps = [float(a) for a in c[1:]]
            coeff = float(np.exp(c[0]))
            # how far is each exponent from a small rational?
            rat = [Fraction(a).limit_denominator(12) for a in exps]
            gaps = [abs(float(rat[i]) - exps[i]) for i in range(dim)]
            worst = max(gaps) if gaps else 0.0
            expr = sp.Float(coeff)
            for i, a in enumerate(exps):
                expr = expr * syms[i] ** sp.Rational(rat[i].numerator, rat[i].denominator)
            conj.append({"form": str(expr), "residual": worst,
                         "raw_exponents": [round(a, 5) for a in exps],
                         "snapped": [f"{r.numerator}/{r.denominator}" for r in rat]})
            if worst < 1e-4:
                diagnosis = {"kind": "pinned",
                             "detail": f"power-law exponents pin to {conj[-1]['snapped']}"}
                next_action = "recover"
            else:
                # a free exponent is NOT pinning to a small rational -> continuum;
                # is it reaching for a nameable constant?
                near = []
                for a in exps:
                    for nm, v in _NAMED:
                        if abs(a - v) < 5e-3:
                            near.append(nm)
                diagnosis = {"kind": "continuum",
                             "detail": (f"a free exponent {[round(a,4) for a in exps]} does not "
                                        f"pin to a small rational"
                                        + (f"; nearest constant(s): {near}" if near else ""))}
                next_action = "declare_and_verify" if near else "acquire_more_data"

    # 2) bounded-grammar rivalry read: how many materially-different forms fit loosely?
    try:
        Xf, yf, Xs, ys, Xc, yc = _split(X, y)
        r = discover(Xf, yf, Xs, ys, Xc, yc, sigma=max(float(sigma), 1e-6))
        if r.expr is not None:
            conj.insert(0, {"form": str(r.expr), "residual": 0.0,
                            "source": "bounded-grammar best certifiable form"})
            if diagnosis["kind"] == "unknown":
                diagnosis = {"kind": "pinned", "detail": "a bounded exact law certifies"}
                next_action = "recover"
        elif r.certificate.abstain == Abstain.STRUCTURAL.value and diagnosis["kind"] == "unknown":
            diagnosis = {"kind": "under_determined",
                         "detail": "multiple materially-different forms fit; not separable here"}
            next_action = "acquire_more_data"
    except Exception:                                          # noqa: BLE001
        pass

    out["conjectures"] = conj[:top]
    out["diagnosis"] = diagnosis
    out["next_action"] = next_action
    return out
