"""The PDE system driver (lagh/pdesystem.py, docs/CASE_STUDY_PDE_C3.md).

What a per-equation loop cannot do: one row set with one target per equation,
time-derivative columns kept OUT of the features, a conjoined certificate whose
alpha is a union bound dominated by the weakest equation, and the truth check
that must run before any abstain is read as a finding.
"""
import numpy as np

from lagh.pdesystem import (agreement, assemble, conjoin, discover_equation,
                            truth_check, weakest)
from lagh.systems import union_alpha_log10
from lagh.weakform import field_terms, multiscale_patches
from lagh.weakform import Term

A, B, C, D = 0.1, 0.5, 0.05, -0.3
TERMS = (field_terms("u", ["u_t", "u_xx", "u_x", "u"])
         + field_terms("v", ["u_t", "u_xx", "u_x", "u"])
         + [Term("1", 0, 0, "1")])
TRUTH = {"u:u_t": {"u:u_xx": A, "v:u": B},
         "v:u_t": {"v:u_xx": C, "u:u": D}}


def exact_pair(seed, nx=193, nt=65, tmax=1.0, modes=(1, 2, 3)):
    """The stage-1 system, exactly: a 2x2 propagator per Fourier mode, with ONE
    phase per mode shared by both fields (independent phases make the coupling
    term orthogonal and the system does not hold at all)."""
    from scipy.linalg import expm
    x = np.linspace(0.0, 2 * np.pi, nx)
    t = np.linspace(0.0, tmax, nt)
    rng = np.random.default_rng(seed)
    amps = rng.uniform(0.2, 1.0, (2, len(modes)))
    ph = rng.uniform(0, 2 * np.pi, len(modes))
    u = np.zeros((nx, nt))
    v = np.zeros((nx, nt))
    for j, k in enumerate(modes):
        M = np.array([[-A * k ** 2, B], [D, -C * k ** 2]])
        for i, ti in enumerate(t):
            at = expm(M * ti) @ np.array([amps[0, j], amps[1, j]])
            u[:, i] += at[0] * np.cos(k * x + ph[j])
            v[:, i] += at[1] * np.cos(k * x + ph[j])
    return {"u": u, "v": v}, (x, t)


def rows(n_sol=3):
    sols = [exact_pair(s) for s in range(n_sol)]
    pf = (lambda coords: multiscale_patches(coords[0], coords[1],
                                            [(20, 10), (28, 14)], n_x=4, n_t=3))
    return assemble(sols, TERMS, pf, sigma=0.0, p=16)


def test_time_derivatives_are_targets_and_never_features():
    r = rows()
    assert set(r.targets()) == {"u:u_t", "v:u_t"}
    assert "u:u_t" not in r.features() and "v:u_t" not in r.features()


def test_the_truth_sits_inside_its_own_band():
    """Run BEFORE reading any abstain as a finding -- the discipline the system
    scoping probe bought the hard way."""
    r = rows()
    for target, truth in TRUTH.items():
        tc = truth_check(r, target, truth)
        assert tc["truth_certifies"], (target, tc)


def test_the_truth_check_reports_vacuity_rather_than_a_green_light():
    """When the band swallows the target, EVERY law sits inside it -- so the
    truth doing so is not evidence. Measured on PDEBench's reaction-diffusion,
    whose field is frozen over most of its record: the check returned True while
    the target column was 1.9e-14 against a band of 37."""
    r = rows()
    j = r.names.index("u:u_t")
    r.det = r.det + 1e3 * np.max(np.abs(r.A))      # a band that swallows all
    tc = truth_check(r, "u:u_t", TRUTH["u:u_t"])
    assert tc["vacuous"]
    assert not tc["truth_certifies"]               # not a green light
    assert tc["signal_to_band"] < 1.0
    assert "VACUOUS" in tc["note"]
    del j


def test_both_equations_certify_the_true_support_over_shared_rows():
    r = rows()
    eqs = [discover_equation(r, t_, sigma=0.0, max_tier=3) for t_ in TRUTH]
    for eq in eqs:
        assert eq["certified"], (eq["target"], eq.get("abstain"))
        assert set(eq["coefficients"]) == set(TRUTH[eq["target"]])
        for k, v in TRUTH[eq["target"]].items():
            lo, hi = eq["intervals"][k]
            assert lo <= v <= hi


def test_a_single_solution_refuses():
    """A coupled system on one trajectory supports only an on-shell statement
    about that trajectory."""
    r = rows(n_sol=1)
    eq = discover_equation(r, "u:u_t", sigma=0.0)
    assert not eq["certified"]
    assert eq["abstain"] == "single-solution"


def test_the_conjoined_alpha_is_dominated_by_the_weakest_equation():
    eqs = [{"target": "a", "certified": True, "alpha_log10": -300.0,
            "coefficients": {"x": 0.5}},
           {"target": "b", "certified": True, "alpha_log10": -50.0,
            "coefficients": {"y": 0.5}}]
    cert = conjoin(eqs)
    assert cert.alpha_log10_total == union_alpha_log10([-300.0, -50.0])
    assert cert.alpha_log10_total > -50.0          # the union bound is WEAKER
    assert weakest(eqs) == ("b", -50.0)
    # a constant shared by two equations is reported as consistent
    assert cert.shared and cert.shared[0]["value"] == 0.5


def test_agreement_separates_a_wrong_law_from_an_under_determined_one():
    """A law that agrees with the truth everywhere its domain claim applies is
    under-determined, not wrong; one that disagrees beyond its band is wrong."""
    r = rows()
    same = agreement(r, "u:u_t", TRUTH["u:u_t"], TRUTH["u:u_t"])
    assert same["agrees_on_certified_domain"]
    off = dict(TRUTH["u:u_t"])
    off["u:u_xx"] = A * 1.5
    bad = agreement(r, "u:u_t", off, TRUTH["u:u_t"])
    assert not bad["agrees_on_certified_domain"]
