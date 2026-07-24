"""The math-class registry. The curriculum is COMPLETE at instrument registration
(docs/DISCOVERER.md 7): per-target class choices do not exist, hence cannot be abused.

A class contributes any of:
    terms(dim, X_fit, y_fit, X_cert) -> list[Term]      library terms (guarded)
    candidates(ctx) -> list[Candidate]                  direct law proposals
    transforms(y)   -> list[(name, t(y), inverse)]      target transforms

The engine escalates tier by tier, only when the certifying set is empty.
"""

from . import (c1_polynomial, c2_rational, c3_powerlaw, c4_inner,
               c5_transforms, c6_quasipoly, c8_angular, c9_genmonomial)

# C1-C5 are FLOAT candidate-list classes driven by the engine's proposal loop.
# C6 (quasi-polynomial, exact integer) is NOT a candidate list -- it is a full
# exact-arithmetic recovery the engine escalates into when C1-C5 fail on an
# integer-lattice target. It is registered here for curriculum completeness; the
# engine invokes c6_quasipoly directly.
CURRICULUM = [
    (1, c1_polynomial),
    (2, c2_rational),
    (3, c3_powerlaw),
    (4, c4_inner),
    (4, c9_genmonomial),
    (5, c5_transforms),
    (5, c8_angular),
]
EXACT_TIERS = [(6, c6_quasipoly)]
