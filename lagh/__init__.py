"""lagh -- a general-purpose certified law discoverer.

Design: wyly's docs/DISCOVERER.md. Product invariant: never a confident wrong
answer -- a symbolic law with an exhaustive certificate over a stated domain, or a
machine-readable refusal saying why.
"""

from .certify import Abstain, Certificate, epsilon, TAU
from .engine import Result, discover

__all__ = ["discover", "Result", "Certificate", "Abstain", "epsilon", "TAU"]
__version__ = "0.1.0"
