"""The math-class registry. The curriculum is COMPLETE at instrument registration
(docs/DISCOVERER.md 7): per-target class choices do not exist, hence cannot be abused.

A class contributes any of:
    terms(dim, X_fit, y_fit, X_cert) -> list[Term]      library terms (guarded)
    candidates(ctx) -> list[Candidate]                  direct law proposals
    transforms(y)   -> list[(name, t(y), inverse)]      target transforms

The engine escalates tier by tier, only when the certifying set is empty.
"""

from . import c1_polynomial, c2_rational, c3_powerlaw, c4_inner, c5_transforms

CURRICULUM = [
    (1, c1_polynomial),
    (2, c2_rational),
    (3, c3_powerlaw),
    (4, c4_inner),
    (5, c5_transforms),
]
