"""Significance arbitration at the Müntz boundary (MUNTZ_ARBITRATION.md).

The rival shapes here are the MEASURED ones from the reach audit, not
hand-imagined whales: the rivals arbitration may dismiss are 60+ term
interpolations whose free parameters consume the whole certification sample
(h/n ~ 0.01), and the rivals it must NOT dismiss retain roughly half the
sample as held-out evidence (rational-d1: dof 34 vs 45 of n = 80).
"""
import numpy as np
import sympy as sp

from lagh.base import Candidate
from lagh.certify import arbitrate_significance, epsilon, free_dof


def _c(s, cx):
    return Candidate(expr=sp.sympify(s), complexity=cx, channel="linear")


def _classes(*cands):
    return [(c.expr, [c]) for c in cands]


def _dense(n_terms, seed):
    """A dense fractional-power fit: 2 free numbers per term (coefficient +
    rational exponent), the shape the linear channel proposes at the Müntz
    boundary."""
    x = sp.Symbol("x_0")
    rng = np.random.default_rng(seed)
    e = sp.Integer(0)
    for k, c in enumerate(rng.uniform(-1, 1, n_terms)):
        e += sp.Float(round(float(c), 6)) * x ** sp.Rational(k + 3, k + 2)
    return e


Y = np.random.default_rng(0).uniform(1, 10, 80)
EPS = epsilon(Y, sigma=0.0)


def test_muntz_contest_resolves_to_low_dof():
    """The measured winning contest: a 2-dof law against an interpolation with
    no held-out evidence left (mixed-4term-d2 / sparse6-d2)."""
    truth = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    whale = Candidate(expr=_dense(40, 1), complexity=113, channel="linear")
    assert free_dof(whale.expr) >= 0.9 * len(Y)
    r = arbitrate_significance(_classes(truth, whale), Y, EPS, 500)
    assert r is not None
    assert r[0][1][0].expr == truth.expr
    assert "domain claim" in r[1] and "evidence-starved" in r[1]


def test_constrained_rival_blocks_arbitration():
    """The rational-d1 counterexample that amended the rule: the margin is
    astronomical (the rival carries far more dof), but the rival still retains
    ~half the sample as held-out evidence -- a law the data constrains, so the
    contest is genuine and the abstain must stand."""
    truth = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    rival = Candidate(expr=_dense(16, 2), complexity=98, channel="linear")
    h = 1 - free_dof(rival.expr) / len(Y)
    assert 0.3 < h < 0.7                      # genuinely constrained, as measured
    assert arbitrate_significance(_classes(truth, rival), Y, EPS, 500) is None


def test_marginal_contest_stays_abstained():
    a = _c("x_0**3 - 2*x_0**2 + x_0 - 1", 6)
    b = _c("3*x_0**2 - 1", 4)
    assert arbitrate_significance(_classes(a, b), Y, EPS, 500) is None


def test_single_class_is_noop():
    a = _c("x_0 + 1", 2)
    assert arbitrate_significance(_classes(a), Y, EPS, 500) is None


def test_two_whales_no_winner():
    w1 = Candidate(expr=_dense(40, 3), complexity=113, channel="linear")
    w2 = Candidate(expr=_dense(39, 4), complexity=110, channel="linear")
    assert arbitrate_significance(_classes(w1, w2), Y, EPS, 500) is None


def test_coherence_stops_once_two_classes_retain_evidence():
    """EXACT early exit, not an approximation: arbitration can crown a winner
    only when every defeated rival is evidence-starved, so two classes that both
    retain held-out evidence settle a structural abstain no matter what else
    would have clustered. Measured motivation: a loose declared band produced 662
    classes on PDEBench CFD and 1553 s of pairwise clustering to reach the
    abstain the second class had already decided."""
    import numpy as np
    import sympy as sp

    from lagh.certify import coherent

    class C:
        def __init__(self, e):
            self.expr, self.complexity = e, 2

    x = sp.Symbol("x_0")
    # materially different, each with ONE free parameter, so h/n is ~1
    cands = [C(sp.Float(c) * x) for c in (1.0, 2.0, 3.0, 4.0, 5.0)]
    P = np.linspace(1.0, 2.0, 60).reshape(-1, 1)
    full = coherent(cands, [x], P, 1.0)
    assert len(full) == 5                       # all distinct without the hint
    early = coherent(cands, [x], P, 1.0, n_evidence=50)
    assert len(early) == 2                      # stopped at the settled verdict
    # ...and the verdict both would produce is the same: structural
    assert len(full) > 1 and len(early) > 1


def test_support_proposal_is_bounded_by_the_certification_split():
    """A support with dof >= n_cert has h = 0, so alpha_log10 = 0 and the
    significance gate demotes it however well it fits. Proposing it can only
    inflate the certifying set that coherence must cluster."""
    from lagh.engine import _basis_supports

    full, reached = _basis_supports(8)
    assert reached == 8 and max(len(s) for s in full) == 8
    capped, reached_c = _basis_supports(8, max_size=4)
    assert reached_c == 4
    assert max(len(s) for s in capped) == 4
    assert capped < full                        # a strict subset, nothing new


def test_partial_determination_reports_what_every_consistent_law_agrees_on():
    """A structural abstain discards the part that WAS determined. Measured on
    PDEBench CFD: 662 classes certified, the verdict reported nothing about the
    coefficients, and the truth check knew the stated law sat at 0.003 of its
    band. The invariant content is a claim about vocabulary+data+band, not about
    nature, so it cannot weaken zero-confident-wrong."""
    import sympy as sp

    from lagh.certify import invariant_content

    class C:
        def __init__(self, e):
            self.expr, self.complexity = e, 2

    x0, x1, x2 = sp.symbols("x_0 x_1 x_2")
    # every consistent law: uses x_0 with a coefficient in [1.9, 2.1], never
    # uses x_2, and uses x_1 only sometimes
    cands = [C(sp.Float(2.0) * x0), C(sp.Float(1.9) * x0 + sp.Float(0.5) * x1),
             C(sp.Float(2.1) * x0 - sp.Float(0.3) * x1)]
    got = invariant_content(cands, [x0, x1, x2], names=["u_xx", "u_x", "u_xxx"])
    assert got["n_certifying_read"] == 3
    assert got["required_terms"] == ["u_xx"]        # in all of them
    assert got["excluded_terms"] == ["u_xxx"]       # in none
    c = got["coefficients"]
    assert c["u_xx"]["lo"] == 1.9 and c["u_xx"]["hi"] == 2.1
    assert not c["u_x"]["always_present"] and not c["u_x"]["never_present"]
    assert "NOT a certificate" in got["claim"]


def test_partial_determination_skips_candidates_it_cannot_read_linearly():
    import sympy as sp

    from lagh.certify import invariant_content

    class C:
        def __init__(self, e):
            self.expr, self.complexity = e, 2

    x0, x1 = sp.symbols("x_0 x_1")
    cands = [C(sp.Float(2.0) * x0), C(x0 * x1), C(sp.sin(x0))]
    got = invariant_content(cands, [x0, x1])
    assert got["n_certifying_read"] == 1            # the nonlinear ones skipped


def test_one_vocabulary_for_partial_determination():
    """Five mechanisms stated partial determination in five encodings. A
    certified verdict and an under-determined one must be readable by one
    consumer, or a per-component checker cannot score either."""
    from lagh.certify import determination

    d = determination([("a", 2.0, 2.0), ("b", -1.0, 0.5), ("c", None, None),
                       ("e", 0.2, 0.9)], status="certified")
    assert d["exact"] == ["a"] and d["unconstrained"] == ["c"]
    assert set(d["interval"]) == {"b", "e"}
    # RESOLVED is the threshold-free 'is it there at all' test: an interval that
    # straddles zero is bounded but does not establish presence
    assert d["components"]["b"]["resolved"] is False      # [-1, 0.5] spans zero
    assert d["components"]["e"]["resolved"] is True       # [0.2, 0.9] does not
    assert d["components"]["a"]["resolved"] is True
    assert d["n_resolved"] == 2
    # the status names WHAT produced it: a range from one certified law and a
    # range over a certifying SET are different claims
    assert d["status"] == "certified"


def test_a_domain_qualifier_scopes_the_whole_record_and_blocks_composition():
    """DOMAIN is the fifth dimension and the one that is not a component: a
    restriction qualifies every entry at once. Registered 2026-07-29 as the
    variable-coefficient route -- certify where the coefficient is locally
    constant and report where. The rule that earns it: two records established
    on DIFFERENT regions must not conjoin, because the conjunction is defined
    only where both were, and this layer cannot intersect two predicates."""
    from lagh.certify import (conjoin_determination, determination,
                              domain_qualifier)

    hi = domain_qualifier("a == a_hi", coverage=0.61)
    lo = domain_qualifier("a == a_lo", coverage=0.39)
    d_hi = determination([("beta", 0.0999, 0.1001)], status="certified",
                         qualifier=hi)
    d_lo = determination([("beta", 0.048, 0.052)], status="certified",
                         qualifier=lo)
    assert d_hi["qualifier"]["predicate"] == "a == a_hi"
    assert determination([("beta", 0.1, 0.1)], status="certified") \
        .get("qualifier") is None                      # unqualified = everywhere

    refused = conjoin_determination([d_hi, d_lo])
    assert refused["status"] == "refused"
    assert refused["domains"] == ["a == a_hi", "a == a_lo"]

    # same domain: the conjunction is the INTERSECTION, and it keeps the scope
    tighter = determination([("beta", 0.09995, 0.10005)], status="certified",
                            qualifier=hi)
    both = conjoin_determination([d_hi, tighter])
    assert both["components"]["beta"]["lo"] == 0.09995
    assert both["components"]["beta"]["hi"] == 0.10005
    assert both["qualifier"]["predicate"] == "a == a_hi"
    assert "contradiction" not in both

    # an EMPTY intersection is a finding about the inputs, not something to drop
    clash = conjoin_determination(
        [d_hi, determination([("beta", 0.2, 0.3)], status="certified",
                             qualifier=hi)])
    assert clash["contradiction"] == ["beta"]
