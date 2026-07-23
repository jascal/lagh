"""characterize(): sound structural diagnosis, and the hard invariant that it can NEVER
be read as a certificate (no `certified: True`, no `law`). See docs/DIRECTION_CHARACTERIZATION.md."""
import numpy as np

from lagh.characterize import characterize


def _X(dim, n=200, lo=0.5, hi=8.0, seed=0):
    rng = np.random.default_rng(seed)
    return np.exp(rng.uniform(np.log(lo), np.log(hi), (n, dim)))


def _never_a_certificate(ch):
    assert ch["certified"] is False
    assert ch["tag"] in ("empirical", "open")
    assert "law" not in ch                       # the structural wall: no form to mistake for a result
    assert ch["kind"] == "characterization"
    assert ch["research"]["move"] in ("acquire_more_data", "acquire_divergent",
                                      "declare_and_verify", "report_and_stop")


def test_clean_power_law_is_classified_and_pins():
    X = _X(2)
    y = 3.0 * X[:, 0] ** 2 * X[:, 1] ** 1  # rational monomial
    ch = characterize(X, y)
    _never_a_certificate(ch)
    assert ch["class"] == "power-law"
    assert ch["power_law"]["pins_to_rational"] is True


def test_irrational_exponent_is_the_wedge_not_a_law():
    # y = x0^e : the NewtonBench wedge -- must be flagged irrational, NEVER given a law
    X = _X(1)
    y = X[:, 0] ** np.e
    ch = characterize(X, y)
    _never_a_certificate(ch)
    assert ch["class"] == "irrational-power"
    assert ch["power_law"]["irrational_hint"]["nearest"] == "e"
    assert ch["power_law"]["pins_to_rational"] is False
    # the wedge is terminal: tell the caller to stop, not thrash
    assert ch["research"]["move"] == "report_and_stop"


def test_exponential_dependence_detected():
    X = _X(1, lo=0.5, hi=5.0)
    y = np.exp(-0.7 * X[:, 0])
    ch = characterize(X, y)
    _never_a_certificate(ch)
    assert ch["class"] == "exponential"


def test_bounded_nonalgebraic_routes_to_declare_and_verify():
    # bounded, non-monotone, no clean power-law -> a trig/inverse-trig-shaped hedge
    X = _X(1, lo=0.1, hi=3.0)
    y = np.sin(X[:, 0]) + 2.0                    # strictly positive, bounded
    ch = characterize(X, y)
    _never_a_certificate(ch)
    assert ch["class"] in ("non-algebraic", "unresolved")
    assert ch["class"] != "power-law"


def test_trend_and_sign_are_reported():
    X = _X(1)
    y = 5.0 * X[:, 0] ** 2                       # increasing, positive
    ch = characterize(X, y)
    assert ch["shape"]["trend"][0] == "increasing"
    assert ch["shape"]["sign"] == "positive"


def test_abstain_reason_shapes_the_research_move():
    X = _X(2)
    y = 3.0 * X[:, 0] * X[:, 1]
    # a structural abstain on an algebraic-looking function -> acquire where rivals diverge
    ch = characterize(X, y, abstain_reason="structural")
    assert ch["research"]["move"] in ("acquire_divergent", "acquire_more_data")


def test_too_few_points_is_graceful():
    X = _X(1, n=4)
    y = X[:, 0] ** 2
    ch = characterize(X, y)
    _never_a_certificate(ch)
    assert ch["class"] == "unresolved"
