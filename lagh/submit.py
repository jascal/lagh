"""Two-track benchmark submission (docs/DIRECTION_OUTPUT_POLICY.md).

Per problem: submit a law that passed the finite-data checks when one exists
(track A, tag `proved`, with its explicit `claim` scope); otherwise the best available LABELED
conjecture (track B, tag `empirical`, never upgraded); otherwise an explicit abstain.
Accuracy-scored benchmarks see both tracks; the certified partition is reported
separately and cannot be redrawn after scoring. Certification does not establish
exact recovery of the generating form or coefficients.
"""

from __future__ import annotations

import numpy as np

from .mcp.core import fit
from .passive import discover_passive


def submission(X, y, *, sigma: float = 0.0, seed: int = 0) -> dict:
    """One benchmark problem -> one submission dict:
    {track: 'certified'|'conjecture'|'abstain', expr, tag, detail}."""
    r = discover_passive(X, y, sigma=sigma, seed=seed)
    if r.certified:
        c = r.result.certificate
        return {"track": "certified", "expr": str(r.result.expr), "tag": "proved",
                "alpha_log10": c.alpha_log10,
                "claim": dict(c.claim),
                "detail": "machine-checked consistency over the dataset's "
                          "finite domain under the supplied error model"
                          + (f"; chance-fit significance alpha <= 1e{c.alpha_log10:.0f}"
                             if c.alpha_log10 is not None else "")}
    # Track B, in the pre-registered order: fit scout, then the engine's best
    # non-certifying candidate, then nothing.
    abstain = r.result.certificate.abstain or "structural"
    Xa = np.atleast_2d(np.asarray(X, float))
    ya = np.asarray(y, float).ravel()
    n_finite = int((np.isfinite(ya) & np.all(np.isfinite(Xa), axis=1)).sum())
    if n_finite < 10:      # a conjecture from (almost) no data is noise, not a guess
        return {"track": "abstain", "expr": None, "tag": "open",
                "detail": f"only {n_finite} finite points ({abstain})"}
    try:
        f = fit(np.asarray(X, float).tolist(), np.asarray(y, float).ravel().tolist(),
                sigma=sigma)
        conj = (f.get("conjectures") or [{}])[0].get("form")
    except Exception:                                          # noqa: BLE001
        conj = None
    if not conj and r.result.certificate.law:
        conj = r.result.certificate.law
    if conj:
        return {"track": "conjecture", "expr": conj, "tag": "empirical",
                "detail": f"lagh abstained ({abstain}); best labeled conjecture, "
                          "NOT certified"}
    return {"track": "abstain", "expr": None, "tag": "open",
            "detail": f"no certificate and no conjecture ({abstain})"}
