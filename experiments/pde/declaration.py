"""What declaration does a stated law actually need, in the WEAK FORM?

The band a weak-form certificate consumes is a bound on a LOCAL violation over
one patch. That is not the same quantity as a solver's POINTWISE deviation from
an exact solution accumulated over a whole trajectory, and PDEBench proved the
difference is not academic: declaring the second where the first was wanted
over-declared by ~3900x and cost three and a half orders of interval width
(`docs/DIRECTION_ERROR_PROVENANCE.md`, and the correction at the top of
`docs/CASE_STUDY_PDEBENCH.md`).

So this module answers the question directly rather than by analogy: given rows
and a stated law, bisect the declared field error to the point where the truth
exactly meets its own band. That number is what the weak form requires; anything
larger is an over-declaration, which is safe and costs width.

It is also the measurement the STRATEGY rule "scan the declaration" asks for --
if a certificate survives three decades below the number declared, the
declaration is not what is binding it and the gap is worth finding before the
result is written up.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.pdesystem import truth_check                             # noqa: E402

DECADES = [0.0] + [10.0 ** e for e in range(-10, 1)]


def required_declaration(base, target, law, sigma, *, grid=None, iters=40):
    """The declared field error at which `law` exactly meets its band.

    Returns (exact, grid_point, scan). `grid_point` is the coarse decadal answer
    the scan alone would have given, kept because it is what a decadal scan
    reports and because the gap between the two is the quantization this module
    exists to remove: a decadal grid rounds every requirement up to the next
    decade, so a ratio of two of them can only ever be a power of ten.

    `exact` is 0.0 when the law holds on sigma alone (nothing to declare) and
    None when it never holds on the grid.

    The bisection is assumption-free: `rebanded` is affine in field_err, so the
    truth's ratio to its band is monotone decreasing in it.
    """
    scan, grid_point, below = [], None, 0.0
    for fe in (grid if grid is not None else DECADES):
        tc = truth_check(base.rebanded(fe), target, law, sigma=sigma)
        holds = bool(tc["truth_certifies"] and not tc["vacuous"])
        scan.append({"field_err": fe, "truth_over_band": tc["truth_max_ratio"],
                     "signal_over_band": tc["signal_to_band"],
                     "holds": holds, "vacuous": bool(tc["vacuous"])})
        if grid_point is None:
            if holds:
                grid_point = fe
            else:
                below = fe
    if grid_point is None:
        return None, None, scan
    lo, hi = below, grid_point
    for _ in range(iters):
        mid = (lo * hi) ** 0.5 if lo > 0 else hi / 2.0
        tc = truth_check(base.rebanded(mid), target, law, sigma=sigma)
        if tc["truth_certifies"] and not tc["vacuous"]:
            hi = mid
        else:
            lo = mid
    return (hi if grid_point else 0.0), grid_point, scan
