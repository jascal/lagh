"""Error provenance: what KIND of error is in this field?
(docs/DIRECTION_ERROR_PROVENANCE.md, registered predictions P1-P6.)

The band has two channels and they are not interchangeable:

  * unstructured/stochastic  ->  KAPPA * sigma * sqrt(a'Ga)   (L2, cancels
    across one realization),
  * structured/deterministic ->  field_err * ||w g'||_1       (L1, nothing
    cancels).

L2 on a deterministic error UNDER-declares and admits impostors; L1 on a
stochastic one OVER-declares by roughly sqrt(n) and loses laws that are there.
Until now that choice was made by hand, once per dataset, by someone who might be
wrong -- and it is why every PDEBench certificate is +-17% rather than +-1e-4.

This module measures the STRUCTURE and reports it. It never returns a
certificate: the verdict is `empirical`, it names what was measured, and
`undetermined` is a first-class answer rather than a default. The usual causes --
an instrument for the unstructured case, a pipeline (a solver, an assimilation
system, a photometric reduction) for the structured one -- are reported as the
likely reading, never as the finding. Nothing here can distinguish a physical
process from a sufficiently good simulation of one, and a certificate never
needed it to.

The sharp diagnostic is the MODIFIED EQUATION. A scheme integrating
u_t + beta u_x = 0 does not solve that equation; it solves
u_t + beta u_x = c2 u_xx + c3 u_xxx + ..., whose leading coefficients are the
discretization's signature (c2 dissipation, c3 dispersion). So the residual of a
stated law, regressed on the next derivatives in the hierarchy, IS that
signature -- and those derivatives are already weak-form columns, so it costs
one lstsq. Measured on PDEBench advection: u_xxx explains 84% of the residual
variance at c3 = 1.09e-7 while u_xx explains 1.5%.

Explaining a residual is NOT certifying a term, and the report says so.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# A higher-derivative term is taken to EXPLAIN the residual at this fraction of
# its variance. Chosen an order below the measured advection signature (0.84) and
# well above the terms that do not explain it there (0.015, 0.009).
EXPLAINS = 0.5
# ...and the residual counts as GROWING when its late-window level exceeds its
# early-window level by this factor. Storage rounding is flat in t; a scheme's
# accumulated error is not (measured: 0 -> 2.75e-2 over the advection window).
GROWTH = 3.0


@dataclass
class ErrorReport:
    verdict: str                       # structured-deterministic | unstructured
                                       # -stochastic | both | undetermined
    likely_reading: str = ""
    recommended_channel: str = ""      # "L1 (field_err)" | "L2 (sigma)" | both
    evidence: dict = field(default_factory=dict)
    modified_equation: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    certified: bool = False            # INVARIANT: never a certificate

    def one_line(self) -> str:
        return (f"{self.verdict} -> {self.recommended_channel}"
                f"  ({self.likely_reading})")


def modified_equation(residual, columns: dict) -> dict:
    """Regress a stated law's residual on the next derivatives in the hierarchy.

    `columns` maps a term name to its weak-form column. Returns each term's
    variance-explained and coefficient, plus the best single explanation and the
    joint fit. A term explaining most of the residual is the discretization's
    leading truncation term -- which is a DIAGNOSIS, not a certificate: it says
    the residual has that shape, not that the term is separably certifiable
    (measured on PDEBench advection, where it is not)."""
    r = np.asarray(residual, float).ravel()
    rr = float(r @ r)
    out = {"residual_l2": float(np.sqrt(rr)), "terms": {}}
    if rr <= 0 or not columns:
        return out
    for name, col in columns.items():
        c = np.asarray(col, float).ravel()
        cc = float(c @ c)
        if cc <= 0:
            continue
        coef = float(c @ r / cc)
        expl = 1.0 - float(np.sum((r - coef * c) ** 2)) / rr
        out["terms"][name] = {"coefficient": coef,
                              "variance_explained": float(expl)}
    if out["terms"]:
        best = max(out["terms"], key=lambda k: out["terms"][k]["variance_explained"])
        out["best"] = best
        out["best_variance_explained"] = out["terms"][best]["variance_explained"]
        out["explains"] = bool(out["best_variance_explained"] >= EXPLAINS)
    return out


def required_declarations(resid, det_l1, gram_a, *, kappa=4.0, q=1.0) -> dict:
    """How large a declaration each CHANNEL would need to cover this residual.

    Analytic rather than a bisection scan: the L1 channel enters the band as
    field_err * ||w g'||_1 and the L2 channel as kappa * sigma * sqrt(a'Ga), so
    the declaration each would need is the residual divided by its own
    multiplier, taken at quantile `q` over the rows.

    These are the numbers a run would otherwise take by hand. They are DERIVED
    FROM THE RESIDUAL BEING BANDED, which is circular if used to certify the law
    that produced them -- see `characterize_rows`, which returns them labelled so
    the circularity cannot travel silently.
    """
    r = np.abs(np.asarray(resid, float).ravel())
    out = {}
    d = np.asarray(det_l1, float).ravel()
    m = d > 0
    if m.any():
        out["field_err"] = float(np.quantile(r[m] / d[m], q))
    g = np.sqrt(np.maximum(np.asarray(gram_a, float).ravel(), 0.0))
    m = g > 0
    if m.any():
        out["sigma"] = float(np.quantile(r[m] / (kappa * g[m]), q))
    return out


def characterize_rows(rows, target, *, extensions=(), sigma=0.0,
                      library=None, stated_law=None, q=1.0):
    """Characterize the error in a weak-form row set, WITHOUT being told the law.

    The residual is taken from the best fit over the REGISTERED library (or from
    `stated_law` when one is known, which is the pipeline-decode case), and the
    modified-equation regression runs against EXTENSION columns held out of that
    fit -- a higher derivative left in the fitting library would simply be
    absorbed, and the residual would look structureless by construction.

    Returns (report, required) where `required` carries the declaration each
    channel would need AND its provenance. The provenance is the point: a
    magnitude derived from the residual of the very law it will then band is
    circular, and a run that certifies at such a declaration must say so.
    """
    import numpy as _np
    names = list(rows.names)
    j = names.index(target)
    ext = [e for e in extensions if e in names]
    lib = [n for n in (library if library is not None else rows.features())
           if n != target and n not in ext]
    X = rows.A[:, [names.index(n) for n in lib]]
    y = rows.A[:, j]
    if stated_law:
        coef = _np.array([stated_law.get(n, 0.0) for n in lib], float)
        provenance = "residual of the STATED law (pipeline decode)"
    else:
        coef, *_ = _np.linalg.lstsq(X, y, rcond=None)
        provenance = ("residual of the best fit over the registered library: a "
                      "LOWER BOUND only, and NOT usable as a declaration -- the "
                      "fit absorbs the very error being measured")
    resid = y - X @ coef
    cols = {e: rows.A[:, names.index(e)] for e in ext}
    rep = characterize(resid, columns=cols or None,
                       field_scale=float(_np.median(_np.abs(y))))
    a = _np.zeros((len(y), len(names)))
    a[:, j] = 1.0
    for n, c in zip(lib, coef):
        a[:, names.index(n)] = -c
    det_l1 = (rows.field_l1[:, j]
              + _np.abs(rows.field_l1[:, [names.index(n) for n in lib]]
                        @ _np.abs(coef))) if rows.field_l1 is not None else None
    gram_a = _np.einsum("nk,nkl,nl->n", a, rows.gram, a)
    req = required_declarations(resid, det_l1, gram_a, q=q) \
        if det_l1 is not None else {}
    req["provenance"] = provenance
    req["from_stated_law"] = bool(stated_law)
    req["usable_as_declaration"] = bool(stated_law)
    req["fitted_coefficients"] = dict(zip(lib, [float(c) for c in coef]))
    req["recommended"] = ("field_err" if rep.verdict == "structured-deterministic"
                          else "sigma" if rep.verdict == "unstructured-stochastic"
                          else "field_err (conservative: the kind is unknown)")
    if not stated_law:
        # MEASURED, and it is the dangerous direction. On PDEBench advection --
        # where the solver error is independently known to be 2.75e-2 -- the
        # free fit moves beta to -0.70003 to soak up the dispersion's mean
        # effect, leaving a residual that is both smaller and less structured:
        # u_xxx's share of it falls from 84% (stated law) to 13.6%, the verdict
        # flips to unstructured-stochastic, and the derived field_err comes out
        # at 2.6e-6 -- FOUR ORDERS too small. Auto-setting a band this way would
        # have made every certificate four orders too tight, which is how
        # confident-wrongs are manufactured. The numbers stay, labelled; the
        # automation does not happen.
        rep.notes.insert(0,
            "NOT A DECLARATION: with no stated law the residual comes from a fit "
            "that absorbs the error being measured, so both the verdict and the "
            "magnitude are biased toward 'small and unstructured'. Measured on "
            "PDEBench advection: 13.6% vs 84% explained, and a magnitude 4 "
            "orders below the independently measured solver error. Use an "
            "independent reference (an exact solution, a separate solve, or "
            "replicates), or declare it and report the sensitivity.")
        rep.verdict = "undetermined"
        rep.recommended_channel = ("undetermined -- the free-fit residual "
                                   "cannot decide it")
    return rep, req


def characterize(residual, *, times=None, columns=None, replicates=None,
                 field_scale=None) -> ErrorReport:
    """Measure the error's STRUCTURE from a stated law's residual.

    `times` -- the time coordinate of each residual row, enabling the growth
    test. `columns` -- higher-derivative weak-form columns for the
    modified-equation regression. `replicates` -- residuals of repeated
    measurements of the SAME state, the one direct way to see a stochastic
    component.

    With none of those the honest answer is `undetermined`, and it is returned:
    a single trajectory with no replicates and no reference cannot distinguish a
    deterministic error that happens to look irregular from a stochastic one.
    """
    r = np.asarray(residual, float).ravel()
    rep = ErrorReport(verdict="undetermined")
    ev = rep.evidence
    ev["n"] = int(r.size)
    ev["residual_median"] = float(np.median(np.abs(r)))
    if field_scale:
        ev["residual_relative"] = float(np.median(np.abs(r)) / field_scale)

    structured = []
    # (a) does it GROW? storage rounding cannot; an accumulated scheme error does
    if times is not None and len(np.unique(times)) >= 4:
        t = np.asarray(times, float).ravel()
        lo, hi = np.quantile(t, 0.25), np.quantile(t, 0.75)
        early = np.abs(r[t <= lo])
        late = np.abs(r[t >= hi])
        if early.size and late.size and np.median(early) > 0:
            g = float(np.median(late) / np.median(early))
            ev["growth_late_over_early"] = g
            if g >= GROWTH:
                structured.append(f"residual grows {g:.1f}x across the window")

    # (b) is it EXPLAINED by a higher-derivative term? the modified equation
    if columns:
        me = modified_equation(r, columns)
        rep.modified_equation = me
        if me.get("explains"):
            structured.append(
                f"{me['best_variance_explained']:.0%} of the residual variance "
                f"is explained by {me['best']} (coefficient "
                f"{me['terms'][me['best']]['coefficient']:.3g}) -- a modified-"
                f"equation signature")

    # (c) do REPLICATES disagree? that part cannot be deterministic
    if replicates is not None:
        reps = np.asarray(replicates, float)
        spread = float(np.median(np.std(reps, axis=0)))
        ev["replicate_spread"] = spread
        ev["stochastic_fraction"] = (float(spread / np.median(np.abs(r)))
                                     if np.median(np.abs(r)) > 0 else None)

    has_stochastic = (ev.get("stochastic_fraction") or 0) > 0.5
    if structured and has_stochastic:
        rep.verdict = "both"
        rep.recommended_channel = "L1 (field_err) AND L2 (sigma)"
    elif structured:
        rep.verdict = "structured-deterministic"
        rep.recommended_channel = "L1 (field_err)"
        rep.likely_reading = ("a pipeline -- a solver, an assimilation system, "
                              "a reduction -- rather than an instrument")
    elif replicates is not None:
        rep.verdict = "unstructured-stochastic"
        rep.recommended_channel = "L2 (sigma)"
        rep.likely_reading = "an instrument rather than a pipeline"
    elif columns or times is not None:
        # the tests ran and found no structure. That is evidence FOR the
        # stochastic reading, but without replicates it is not conclusive: a
        # deterministic error orthogonal to every column tested looks identical.
        rep.verdict = "unstructured-stochastic"
        rep.recommended_channel = "L2 (sigma)"
        rep.likely_reading = "an instrument rather than a pipeline"
        rep.notes.append(
            "no structure found by the tests that ran; without replicates this "
            "is evidence for the stochastic reading and not a demonstration -- "
            "a deterministic error orthogonal to every column tested would look "
            "the same")
    else:
        rep.notes.append(
            "no growth test, no columns and no replicates: nothing was measured "
            "that could separate the two, so the answer is UNDETERMINED rather "
            "than a default")
        rep.recommended_channel = ("L1 (field_err) -- the conservative choice "
                                   "while the kind is unknown")
    rep.notes.extend(structured)
    rep.notes.append("empirical, never a certificate: this reports measured "
                     "STRUCTURE; instrument and pipeline are the usual causes "
                     "and are the likely reading, not the finding")
    return rep
