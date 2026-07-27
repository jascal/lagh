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

from ..acquisition import run_active, run_active_boxsearch
from ..base import eval_expr, lstsq
from ..certify import (Abstain, coherent, epsilon, pinned, sample_box)
from ..characterize import characterize
from ..engine import discover
from ..passive import discover_passive

# nameable constants a free-fit exponent might be reaching for (fit's diagnosis)
_NAMED = [("e", float(np.e)), ("pi", float(np.pi)), ("sqrt2", float(np.sqrt(2))),
          ("1/e", 1.0 / float(np.e)), ("ln2", float(np.log(2))),
          ("golden", (1 + 5 ** 0.5) / 2)]


def _syms(dim: int) -> list[sp.Symbol]:
    return [sp.Symbol(f"x_{i}") for i in range(dim)]


def _characterize_oracle(oracle, box_final, abstain_reason, seed=0):
    """Characterize the abstain of an ACTIVE (oracle) recover: the acquired points are not
    kept on ActiveResult, so re-sample box_final (log-uniform, as acquisition does) and
    diagnose. Isolated and best-effort -- any failure just omits the characterization."""
    try:
        box = np.asarray(box_final, float)
        lo, hi = box[0], box[1]
        rng = np.random.default_rng(seed + 7)
        Xc = np.exp(rng.uniform(np.log(np.maximum(lo, 1e-12)), np.log(np.maximum(hi, 1e-12)),
                                (200, len(lo))))
        yc = np.asarray(oracle(Xc), float)
        return characterize(Xc, yc, abstain_reason=abstain_reason)
    except Exception:                                              # noqa: BLE001
        return None


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
    """A DECLARED irrational constant the data cannot pin -- the symbolic e/pi/golden
    that a declared form like x**E introduces. A fitted Float coefficient (a physical
    constant like 6.674e-5) is NOT this: it is an identified value, so it must not make
    the law 'consistent' -- that would mislabel every constant-carrying law."""
    return bool(expr.has(sp.E, sp.pi, sp.GoldenRatio, sp.EulerGamma))


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

def recover(X=None, y=None, *, oracle=None, box=None, sigma: float = 0.0,
            floor_abs: float = 1e-12,
            max_tier: int = 7, budget: int = 200, box_search: bool = False,
            seed: int = 0, time_budget_s: float | None = 45.0) -> dict:
    """Bounded. Discover an exact law. Two modes:

    * **active** (`oracle` + `box` given, in-process only): lagh DRIVES the oracle --
      adaptive ranging, budget-metered multi-objective queries, per-round
      micro-predictions -- via `run_active` (or `run_active_boxsearch` on abstain if
      `box_search=True`, the "broaden-the-box" ladder with a held-out-box guard). This
      is the real acquisition loop; the certificate comes with `acquisition` provenance
      (queries used, budget spent, the box it settled on).
    * **data** (`X`, `y` given, the MCP-wire default): one-shot certification of the
      points you provide. On abstain it returns `next_action:"acquire"` + a broadened
      `suggested_box` so the CALLER can run the acquisition loop (JSON can't carry a
      live oracle across the wire).

    Returns a Certificate or an Abstention -- never a bare guess."""
    # ---- active-acquisition path: lagh drives a live oracle ----
    if oracle is not None:
        if not box or len(box) != 2:
            return {"tag": "open", "tool": "recover", "certified": False,
                    "abstain": "bad-request",
                    "note": "active recover needs box=[lo_vector, hi_vector]"}
        lo, hi = np.asarray(box[0], float), np.asarray(box[1], float)
        sig = float(sigma) if sigma and sigma > 0 else None
        if box_search:
            bs = run_active_boxsearch(oracle, lo, hi, budget=budget, seed=seed,
                                      time_budget_s=time_budget_s)
            active = bs.active
        else:
            active = run_active(oracle, lo, hi, budget=budget, sigma_declared=sig, seed=seed,
                                time_budget_s=time_budget_s)
        r = active.result
        c = r.certificate
        bf = np.asarray(active.box_final, float)
        acq = {"mode": "box-search" if box_search else "active",
               "queries_used": int(active.queries_used),
               "budget_spent": int(active.ledger.spent),
               "box_final": bf.tolist(), "ranging_steps": len(active.ranging_trajectory)}
        if box_search:
            acq["boxes_tried"] = bs.boxes_tried
            acq["heldout_box_ok"] = bs.heldout_box_ok
        if not c.certified:
            ch = _characterize_oracle(oracle, bf, c.abstain, seed)
            out = {"tag": "open", "tool": "recover", "certified": False,
                   "abstain": c.abstain, "domain_size": c.domain_size,
                   "acquisition": acq,
                   "note": "; ".join(map(str, c.notes)) if c.notes else ""}
            if ch is not None:                    # middle rung: a hedged diagnosis, not a law
                out["characterization"] = ch
                out["next_action"] = ch["research"]["move"]
            return out
        # the parametric gate already ran inside discover() -> certified ⇒ pinned
        strength = "consistent" if _has_irrational(r.expr) else "pinned"
        return {"tag": "proved", "tool": "recover", "certified": True,
                "law": str(r.expr), "strength": strength, "domain_size": c.domain_size,
                "tier": r.tier, "bounds": bf.tolist(), "acquisition": acq,
                "note": "certified over the actively-acquired domain, not proved for the world"}

    # ---- data-only path (MCP wire): the PASSIVE regime (docs/DIRECTION_PASSIVE.md):
    # K deterministic re-splits + the full-data exhaustive gate, one code path shared
    # with benchmark submission track A.
    X, y = _prep(X, y)
    if len(X) < 8:
        return {"tag": "open", "tool": "recover", "certified": False,
                "abstain": Abstain.RANGE.value,
                "note": f"only {len(X)} finite points; too thin to certify"}
    dim = X.shape[1]
    syms = _syms(dim)
    if len(X) < 15:
        # TINY-DATA mode (measured need: Kepler III from 8 planets): the
        # 60/20/20 split machinery collapses below ~15 points. Certify
        # EXHAUSTIVELY on all n points with fit=sel=cert; the selection
        # exposure this permits is exactly what alpha = |H|*q^h quantifies
        # (h is dof-discounted), and the certificate says so.
        r = discover(X, y, X, y, X, y, sigma=float(sigma),
                     floor_abs=float(floor_abs), max_tier=max_tier)
        r.certificate.notes.append(
            "tiny-data mode: exhaustive certification on all points; "
            "selection exposure bounded by the stated alpha")
        c = r.certificate
    else:
        pr = discover_passive(X, y, sigma=float(sigma),
                              floor_abs=float(floor_abs), max_tier=max_tier,
                              seed=seed)
        r = pr.result
        c = r.certificate
    if not c.certified:
        lo, hi = X.min(axis=0), X.max(axis=0)
        ch = characterize(X, y, sigma=float(sigma), abstain_reason=c.abstain)
        return {"tag": "open", "tool": "recover", "certified": False,
                "abstain": c.abstain, "domain_size": c.domain_size,
                "next_action": ch["research"]["move"],
                "characterization": ch,
                "suggested_box": [(lo / 10).tolist(), (hi * 10).tolist()],
                "note": ((("; ".join(map(str, c.notes)) + " | ") if c.notes else "")
                         + "see characterization.research for the next move")}
    eps = epsilon(y, sigma=float(sigma), floor_abs=float(floor_abs))
    return {"tag": "proved", "tool": "recover", "certified": True,
            "law": str(r.expr), "strength": _strength(r.expr, syms, X, y, eps, sigma),
            "alpha_log10": c.alpha_log10, "n_hypotheses": c.n_hypotheses,
            "domain_size": c.domain_size, "tier": r.tier,
            "bounds": [[float(X[:, j].min()), float(X[:, j].max())] for j in range(dim)],
            "note": "certified over the stated finite domain, not proved for the world"}


# ---------------------------------------------------------------------------- verify

def verify(X, y, form: str, *, sigma: float = 0.0,
           floor_abs: float = 1e-12) -> dict:
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
    eps = epsilon(yc, sigma=float(sigma), floor_abs=float(floor_abs))
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

    # last-resort AFFINE conjecture (always available -- the power-law probe
    # requires positive data and real-world series cross zero; measured: Okun
    # on raw macro data abstained with NO conjecture)
    try:
        A = np.column_stack([np.ones(len(X))] + [X[:, j] for j in range(dim)])
        c_aff = lstsq(A, y)
        if c_aff is not None:
            e_aff = sp.Float(float(c_aff[0])) + sum(
                sp.Float(float(c_aff[j + 1])) * syms[j] for j in range(dim))
            rr = float(np.sqrt(np.mean((A @ c_aff - y) ** 2)))
            conj.append({"form": str(e_aff), "residual": rr,
                         "source": "affine-OLS fallback"})
    except Exception:                                          # noqa: BLE001
        pass
    out["conjectures"] = conj[:top]
    out["diagnosis"] = diagnosis
    out["next_action"] = next_action
    return out
