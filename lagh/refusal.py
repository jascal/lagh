"""Measurements retained from failed checks; no causal attribution or new gates."""
import numpy as np

from .certify import determination, domain_qualifier


def residual_measurement(y, pred, eps, row_indices, *, domain):
    """Signed y-pred and declared epsilon envelopes on the checked rows only.

    Envelopes are conditional on the supplied error model and the fitted candidate;
    they are not confidence intervals for an omitted physical parameter. No refit.
    Return no diagnostic if finite arithmetic cannot represent the measurement.
    """
    if pred is None:
        return {}
    y, pred, eps = (np.asarray(v, float).ravel() for v in (y, pred, eps))
    if not (y.shape == pred.shape == eps.shape) or not len(y):
        return {}
    with np.errstate(over="ignore", invalid="ignore"):
        residual = y - pred
        lo, hi = residual - eps, residual + eps
    if np.any(eps < 0) or not np.all(np.isfinite([residual, lo, hi])):
        return {}
    partial = determination(
        [("row_discrepancy_envelope", float(lo.min()), float(hi.max()))],
        status="refuted-form", qualifier=domain_qualifier(domain),
        note="envelope over checked rows, conditional on declared epsilon; "
             "not a common offset estimate or a parameter confidence interval")
    return {"measurement": {
        "tag": "empirical", "quantity": "observation-minus-candidate",
        "row_indices": [int(i) for i in row_indices],
        "residual": residual.tolist(), "epsilon": eps.tolist(),
        "domain": domain, "attribution": "unresolved",
        "extra_observable": "independent instrument calibration or a controlled "
                            "intervention separating the proposed physical causes",
        "uncertainty": "declared epsilon envelope; no additional coverage claim"},
        "partial": partial}
