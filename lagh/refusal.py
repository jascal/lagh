"""Bounded measurements from failed checks; no attribution or determination."""
import numpy as np

MAX_RESIDUAL_ROWS = 64


def residual_measurement(y, pred, eps, row_indices, *, domain, candidate):
    """Report a bounded sample and aggregate counts on already checked rows.

    Accept evaluated predictions and scalar/row epsilon only. Callable bands must
    be resolved by the caller for the checked candidate. Invalid input reports an
    omission, never a determination. Undefined rows do not discard finite peers.
    """
    def omitted(reason):
        return {"measurement_omitted": reason}

    if pred is None:
        return omitted("prediction unavailable")
    if callable(eps):
        return omitted("band must be evaluated for the checked candidate")
    try:
        y, pred = (np.asarray(v, float).ravel() for v in (y, pred))
        eps = np.broadcast_to(np.asarray(eps, float), y.shape)
        indices = np.asarray(row_indices)
        if (y.shape != pred.shape or indices.shape != y.shape
                or indices.dtype.kind not in "iu" or np.any(indices < 0)
                or len(np.unique(indices)) != len(indices)):
            return omitted("invalid row alignment")
        if not len(y):
            return omitted("no checked rows")
    except (TypeError, ValueError, OverflowError):
        return omitted("invalid arrays or band shape")
    with np.errstate(over="ignore", invalid="ignore"):
        residual = y - pred
        excess = np.abs(residual) - eps
    valid = np.isfinite(residual) & np.isfinite(eps) & (eps >= 0) & np.isfinite(excess)
    available = np.flatnonzero(valid)
    # Stable order breaks ties by checked-row order, with misses before agreement.
    chosen = available[np.argsort(-excess[available], kind="stable")[:MAX_RESIDUAL_ROWS]]
    invalid_count = int((~valid).sum())
    if not len(available):
        return {**omitted("no finite residuals with valid bands"),
                "measurement_invalid_rows": invalid_count}
    return {"measurement": {
        "evidence": "empirical", "quantity": "observation-minus-candidate",
        "candidate": str(candidate), "domain": str(domain),
        "row_indices": [int(i) for i in indices[chosen]],
        "residual": residual[chosen].tolist(), "epsilon": eps[chosen].tolist(),
        "n_checked": len(y), "n_measurable": len(available),
        "n_exceeding": int(np.sum(excess[available] > 0)),
        "max_band_excess": float(excess[available].max()),
        "n_invalid": invalid_count, "n_elided": len(available) - len(chosen),
        "selection": "up to 64 rows, descending absolute residual minus epsilon",
        "attribution": "unresolved",
        "extra_observable": "independent instrument calibration or a controlled "
                            "intervention separating the proposed physical causes",
        "uncertainty": "declared epsilon envelope; no additional coverage claim"}}
