"""Characterization-on-abstain: a SOUND, UNCERTIFIED structural diagnosis of y=f(X).

The middle rung of the degradation ladder (docs/DIRECTION_CHARACTERIZATION.md). When
recover/discover abstains, a bare reason enum wastes the caller -- it gives up or
thrashes. This produces grounded, measured observations about the function (trend,
power-law slope, whether that slope is a clean rational or an irrational reaching for a
named constant) plus a single research-move pointer for the caller's next step.

INVARIANT: the return is `empirical`, never a certificate. `certified` is always False and
there is no `law` key. It only reports what the samples MEASURE; it never snaps to a
specific form and calls it the answer. It cannot produce a confident-wrong -- there is no
certificate to be wrong about.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np

from .base import lstsq

# nameable constants a free exponent might be reaching for -- an irrational exponent that
# lands on one of these is the "wedge": no exact RATIONAL closed form exists.
_NAMED = [("e", float(np.e)), ("pi", float(np.pi)), ("sqrt2", float(np.sqrt(2))),
          ("1/e", 1.0 / float(np.e)), ("ln2", float(np.log(2))),
          ("golden", (1 + 5 ** 0.5) / 2)]

_RAT_TOL = 1e-3      # exponent within this of a small rational => "pins"
_NAMED_TOL = 5e-3    # non-pinning exponent within this of a named constant => wedge hint
_FIT_CLEAN = 0.05    # log-fit residual below this => a clean fit of that shape
_FIT_POOR = 0.10     # log-fit residual above this => that shape does NOT explain the data
_STRUCT_RESID = 0.03  # residual below this => a near-monomial worth a research move; above
                      # it, sampling won't help -> report_and_stop (keeps the timeout bounded).
                      # Tight on purpose: the only observed gain (underdamped, residual ~0.01)
                      # is a near-monomial; hopeless cells (snell 0.05-0.5, BE 0.08) sit above.


def _prep(X, y):
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    y = np.asarray(y, float).ravel()
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    return X[m], y[m]


def _trend(x, y):
    """Sign of rank-correlation -- a trend read robust to the other inputs' spread.
    Sound as 'y tends to rise/fall with this input', NOT a monotonicity certificate."""
    if len(x) < 8:
        return "none"
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.sqrt(np.dot(rx, rx) * np.dot(ry, ry)))
    if denom <= 0:
        return "none"
    rho = float(np.dot(rx, ry)) / denom
    return "increasing" if rho > 0.3 else "decreasing" if rho < -0.3 else "none"


def _power_law(X, y):
    """log-log least-squares exponents + the rationality flag (the wedge detector).
    Only valid on strictly-positive data; returns None otherwise."""
    if not (np.all(X > 0) and np.all(y > 0)) or len(X) < 8:
        return None
    dim = X.shape[1]
    L = np.column_stack([np.ones(len(X)), np.log(X)])
    c = lstsq(L, np.log(y))
    if c is None:
        return None
    exps = [float(a) for a in c[1:]]
    resid = float(np.sqrt(np.mean((L @ c - np.log(y)) ** 2)))   # pure-monomial fit residual
    rats = [Fraction(a).limit_denominator(12) for a in exps]
    gaps = [abs(float(rats[i]) - exps[i]) for i in range(dim)]
    pins = all(g < _RAT_TOL for g in gaps)
    irr = None
    # An irrational-exponent claim is only meaningful when the data IS a clean monomial:
    # for a compound/multi-term law the log-log "exponent" is a meaningless summary that
    # coincidentally lands near a named constant (measured false positives near ln2/sqrt2
    # on hard multi-term cells). Require a low fit residual before trusting the hint --
    # under-claim (miss a compound irrational) rather than mislead. A genuine clean
    # irrational monomial (x**e) has residual ~0 and is still caught.
    if not pins and resid < _FIT_CLEAN:
        for i, a in enumerate(exps):
            if gaps[i] < _RAT_TOL:
                continue
            for nm, v in _NAMED:
                if abs(a - v) < _NAMED_TOL:
                    irr = {"input": i, "value": round(a, 5), "nearest": nm}
                    break
            if irr:
                break
    return {"exponents": [round(a, 5) for a in exps],
            "snapped": [f"{r.numerator}/{r.denominator}" for r in rats],
            "pins_to_rational": bool(pins), "fit_residual": round(resid, 5),
            "irrational_hint": irr}


def _exponential(X, y):
    """log-LINEAR fit (log y vs x, not log x): flags an exp(k*x) dependence. Returns the
    best single-input log-linear residual for the class read."""
    if not np.all(y > 0) or len(X) < 8:
        return None
    best = None
    for j in range(X.shape[1]):
        A = np.column_stack([np.ones(len(X)), X[:, j]])
        c = lstsq(A, np.log(y))
        if c is None:
            continue
        resid = float(np.sqrt(np.mean((A @ c - np.log(y)) ** 2)))
        if best is None or resid < best[1]:
            best = (j, resid, float(c[1]))
    if best is None:
        return None
    return {"input": best[0], "fit_residual": round(best[1], 5), "rate": round(best[2], 5)}


def _has_structure(pl, ex):
    """Did lagh find something EXPLOITABLE -- a NEAR-MONOMIAL (residual < _STRUCT_RESID) or a
    clean exponential? If not, more sampling won't help, so the caller must stop, not thrash.
    This discriminator keeps the timeout bounded: the only observed gain candidate
    (underdamped, residual ~0.01) has structure; hopeless out-of-grammar cells (snell ~0.05-0.5,
    BE ~0.08) do not. Kept tight ON PURPOSE -- a false 'has structure' costs a research loop."""
    return bool((pl and pl["fit_residual"] < _STRUCT_RESID)
                or (ex and ex["fit_residual"] < _STRUCT_RESID))


def _research_move(cls, pl, ex, abstain_reason):
    """The playbook pointer: one next move for the caller. `report_and_stop` is the terminal
    signal that ends thrashing -- it is the DEFAULT for any cell without exploitable
    structure, so hopeless cells abstain fast instead of researching to the timeout. Only a
    cell where lagh actually found near-monomial structure earns a (bounded) research move,
    and that move uses the FAST data-path lagh.recover -- never a second slow lab.discover."""
    if cls == "irrational-power":
        return ("report_and_stop",
                "an exact rational law is impossible here; report the characterization as "
                "a hedge and abstain -- do NOT keep sampling")
    if cls == "non-algebraic":
        # bounded output + no power-law/exp fit -> a trig/inverse-trig shot is worth ONE
        # propose+verify. In the machine this is structurally bounded (proposing->verifying
        # ->terminal, no loop); the sound checker rejects a wrong form, so it is safe.
        return ("declare_and_verify",
                "declare a trig / inverse-trig (or saturating) form and lagh.verify it ONCE; "
                "if it does not certify, report_and_stop")
    if not _has_structure(pl, ex):
        return ("report_and_stop",
                "no exploitable structure in the samples -- more sampling is unlikely to "
                "help; report the characterization as a hedge and abstain. (You MAY try ONE "
                "declared form via lagh.verify if you have a strong prior, then stop.)")
    # near-monomial structure WAS found -> ONE bounded research move is worth it
    if abstain_reason == "structural":
        return ("acquire_divergent",
                "near-monomial but rivals fit the sampled range; hand-sample a WIDER / "
                "asymptotic regime where they separate, then call lagh.recover(X, y) on those "
                "points (NOT lab.discover -- that re-runs the slow loop)")
    return ("acquire_more_data",
            "near-monomial that did not pin; hand-sample independent-random points across a "
            "wider box, then call lagh.recover(X, y) (NOT lab.discover)")


def _classify(pl, ex, bounded):
    """Conservative class synthesis from the probes. Order = most-specific first."""
    if pl and pl["irrational_hint"]:
        h = pl["irrational_hint"]
        return ("irrational-power",
                f"power-law exponent on input x_{h['input']} is {h['value']} (irrational, "
                f"~{h['nearest']}); no exact RATIONAL closed form exists")
    if pl and pl["pins_to_rational"] and pl["fit_residual"] < _FIT_CLEAN:
        return ("power-law", f"clean monomial: exponents pin to {pl['snapped']}")
    if ex and ex["fit_residual"] < _FIT_CLEAN and (pl is None or ex["fit_residual"] < pl["fit_residual"]):
        return ("exponential",
                f"log(y) is ~linear in x_{ex['input']} (rate {ex['rate']}): an exp dependence")
    if bounded and (pl is None or pl["fit_residual"] > _FIT_POOR):
        return ("non-algebraic",
                "output stays bounded and no clean power-law/exponential fits -- a periodic "
                "/ inverse-trig / saturating component is likely")
    if pl and pl["fit_residual"] > _FIT_POOR:
        return ("additive-or-mixed",
                "a single monomial does not fit (large log-log residual) -- likely a "
                "multi-term / additive law")
    return ("unresolved", "no clean algebraic structure detected in the samples")


def characterize(X, y, *, sigma: float = 0.0, abstain_reason: str | None = None) -> dict:
    """Sound, UNCERTIFIED structural diagnosis of y=f(X). See the module docstring for the
    invariant. Never raises for data reasons -- returns an `unresolved` characterization
    instead, so it can be attached to any abstain without risk."""
    try:
        X, y = _prep(X, y)
        if len(X) < 8:
            return {"tag": "open", "certified": False, "kind": "characterization",
                    "class": "unresolved", "why": "too few finite points to characterize",
                    "research": {"move": "acquire_more_data",
                                 "detail": "sample more points and widen the box"},
                    "note": "not a certificate"}
        dim = X.shape[1]
        ymin, ymax = float(np.min(y)), float(np.max(y))
        sign = "positive" if ymin > 0 else "negative" if ymax < 0 else "mixed"
        yscale = float(np.median(np.abs(y))) + 1e-300
        bounded = float(np.max(np.abs(y))) / yscale < 20.0
        shape = {"trend": [_trend(X[:, j], y) for j in range(dim)],
                 "sign": sign, "bounded": bool(bounded)}
        pl = _power_law(X, y)
        ex = _exponential(X, y)
        cls, why = _classify(pl, ex, bounded)
        move, detail = _research_move(cls, pl, ex, abstain_reason)
        out = {"tag": "empirical", "certified": False, "kind": "characterization",
               "shape": shape, "class": cls, "why": why,
               "research": {"move": move, "detail": detail},
               "note": "a GUESS about the function's structure, NOT a certificate; report "
                       "it as a hedge, never as a law"}
        if pl:
            out["power_law"] = pl
        if ex:
            out["exponential"] = ex
        return out
    except Exception as e:                                          # noqa: BLE001
        # characterization must NEVER break the abstain it decorates
        return {"tag": "open", "certified": False, "kind": "characterization",
                "class": "unresolved", "why": f"characterization failed: {str(e)[:80]}",
                "research": {"move": "report_and_stop", "detail": "no diagnosis available"},
                "note": "not a certificate"}
