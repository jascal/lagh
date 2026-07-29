"""State certificates (lagh/statecert.py, docs/CASE_STUDY_PDE_C4.md).

A state certificate rests on one identity -- the boundary term a one-sided-in-
time test function keeps -- and on one honesty property: a mode the observations
do not determine must come back UNDETERMINED rather than as a number. These
check both, plus the resolution bound that says how many modes may be claimed at
all.
"""
import numpy as np

from lagh.statecert import (MIN_HELDOUT, assemble_state, backpropagate,
                            certify_state, fourier_basis, truth_check_state)
from lagh.weakform import bump_derivatives, ic_columns, onesided_patches

NX, NT = 256, 161
COEFFS = {1: 1.0, 2: 0.6, 3: 0.4}


def grid():
    return (np.linspace(0.0, 2 * np.pi, NX, endpoint=False),
            np.linspace(0.0, 0.6, NT))


def heat_field(x, t, nu=0.1, coeffs=COEFFS):
    u = np.zeros((len(x), len(t)))
    for k, a in coeffs.items():
        u += a * np.exp(-nu * k ** 2 * t)[None, :] * np.cos(k * x)[:, None]
    return u


def state_patches(x, t, t0_index=0, n_x=12):
    out = []
    for nxh, cells in ((10, 40), (14, 56)):
        out += onesided_patches(x, t, nx_half=nxh, nt_cells=cells, n_x=n_x,
                                t0_index=t0_index)
    return out


def test_one_sided_window_is_one_at_the_initial_time():
    """The whole mechanism: phi(., t0) != 0, so the by-parts identity in time
    KEEPS its boundary term, and that term is the initial condition."""
    one = bump_derivatives(8, 0, onesided=True)[0]
    assert np.isclose(one(np.array([0.0]))[0], 1.0)
    assert np.isclose(one(np.array([1.0]))[0], 0.0)     # vanishes at the far end
    assert np.isclose(one(np.array([-0.5]))[0], 0.0)    # and behind t0
    two = bump_derivatives(8, 0)[0]
    assert np.isclose(two(np.array([-0.5]))[0], two(np.array([0.5]))[0])


def test_ic_columns_are_analytic_and_carry_no_data_error():
    x, t = grid()
    pa = state_patches(x, t)
    labels, fns = fourier_basis(3)
    B, err = ic_columns(fns, x, pa, p=16)
    assert B.shape == (len(pa), len(labels))
    # the columns are integrals of two KNOWN smooth functions: their error is a
    # quadrature error at the float floor, not a measurement error
    assert np.max(err) < 1e-12 * max(1.0, np.max(np.abs(B)))


def test_a_clean_state_certifies_and_contains_the_truth():
    x, t = grid()
    u = heat_field(x, t)
    labels, fns = fourier_basis(4)
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": 0.1},
                                     state_patches(x, t), fns, p=16, sigma=1e-6)
    cert = certify_state(B, y, eps, labels, window=(0.0, float(t[-1])),
                         info=info)
    assert cert.certified, cert.abstain
    for k, a in COEFFS.items():
        lo, hi = cert.modes[f"cos{k}"]["interval"]
        assert lo <= a <= hi                      # zero confident-wrong
    # a mode that is not in the field must be BOUNDED but not RESOLVED: its
    # interval straddles zero, and reporting it as a recovered mode would be
    # the whole failure mode this distinction exists to prevent
    assert cert.modes["cos4"]["determined"]
    assert not cert.modes["cos4"]["resolved"]
    assert cert.modes["cos1"]["resolved"]


def test_the_resolution_bound_refuses_before_it_interpolates():
    """dof must stay below the number of independent patch equations, with
    margin -- the 35-term-interpolation lesson, stated up front."""
    x, t = grid()
    u = heat_field(x, t)
    labels, fns = fourier_basis(40)               # 81 unknowns
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": 0.1},
                                     state_patches(x, t, n_x=6), fns, p=16,
                                     sigma=1e-6)
    cert = certify_state(B, y, eps, labels, info=info)
    assert not cert.certified
    assert cert.abstain == "resolution"
    assert cert.heldout < MIN_HELDOUT


def test_a_basis_that_cannot_represent_the_state_refuses():
    """A state the declared basis cannot express must REFUSE, and the truth
    check must show why -- the basis, not the instrument."""
    x, t = grid()
    u = heat_field(x, t, coeffs={1: 1.0, 7: 0.5})   # mode 7 is outside the basis
    labels, fns = fourier_basis(3)
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": 0.1},
                                     state_patches(x, t), fns, p=16, sigma=1e-8)
    cert = certify_state(B, y, eps, labels, info=info)
    assert not cert.certified
    assert cert.abstain == "no-state-explains-the-observations"
    a_true = np.zeros(len(labels))
    a_true[labels.index("cos1")] = 1.0
    assert truth_check_state(B, y, eps, a_true)["truth_max_ratio"] > 1.0


def test_backpropagation_widens_by_the_inverse_decay():
    """Carrying a certificate of the state at t0 back to t = 0 multiplies each
    mode's interval by exp(+nu k^2 t0): the ill-posedness is in the propagation,
    not in the reading."""
    x, t = grid()
    nu, i0 = 0.1, 40
    u = heat_field(x, t, nu=nu)
    labels, fns = fourier_basis(4)
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": nu},
                                     state_patches(x, t, t0_index=i0), fns,
                                     p=16, sigma=1e-6)
    cert = certify_state(B, y, eps, labels, info=info)
    assert cert.certified, cert.abstain
    bp = backpropagate(cert, "heat", float(t[i0]), nu=nu)
    for k in (1, 2, 3):
        gain = np.exp(-nu * k ** 2 * float(t[i0]))
        got = bp["modes"][f"cos{k}"]["half_width"]
        want = cert.modes[f"cos{k}"]["half_width"] / gain
        assert np.isclose(got, want, rtol=1e-9)
        # ...and the back-propagated interval still contains the true amplitude
        lo, hi = bp["modes"][f"cos{k}"]["interval"]
        assert lo <= COEFFS[k] <= hi
    assert bp["k_cut"] >= 3


def test_joint_intervals_are_at_least_as_wide_as_conditional_ones():
    """The reported interval is the projection of the whole feasible set; the
    conditional bisection holds every other mode fixed and must be narrower."""
    x, t = grid()
    u = heat_field(x, t)
    labels, fns = fourier_basis(3)
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": 0.1},
                                     state_patches(x, t), fns, p=16, sigma=1e-5)
    cert = certify_state(B, y, eps, labels, info=info)
    assert cert.certified, cert.abstain
    for m in cert.modes.values():
        if m["interval"] is None or m["conditional_interval"] is None:
            continue
        joint = m["interval"][1] - m["interval"][0]
        cond = m["conditional_interval"][1] - m["conditional_interval"][0]
        assert joint >= cond * (1 - 1e-6)


def test_the_partial_record_says_the_same_thing_as_the_modes_dict():
    """MODE is one of the five places lagh states partial determination, and
    since `certify.determination` exists they must not drift: a mode the joint
    projection leaves unbounded is `unconstrained` there, a bounded one is an
    `interval` with the same bounds, and `resolved` (the interval excludes zero,
    so the certificate can tell the mode is present at all) must agree
    component-wise. Two encodings of one claim is how the vocabulary was lost
    the first time."""
    x, t = grid()
    u = heat_field(x, t)
    labels, fns = fourier_basis(3)
    B, y, eps, info = assemble_state(u, x, t, {"u_xx": 0.1},
                                     state_patches(x, t), fns, p=16, sigma=1e-5)
    cert = certify_state(B, y, eps, labels, info=info)
    assert cert.certified, cert.abstain
    assert cert.partial["status"] == "state"
    assert set(cert.partial["components"]) == set(cert.modes)
    assert sorted(cert.partial["unconstrained"]) == sorted(cert.undetermined)
    for lab, m in cert.modes.items():
        rec = cert.partial["components"][lab]
        assert rec["resolved"] == m["resolved"]
        if m["interval"] is None:
            assert rec["kind"] == "unconstrained"
            assert not m["determined"]
        else:
            assert rec["kind"] in ("interval", "exact")
            assert rec["lo"] == m["interval"][0]
            assert rec["hi"] == m["interval"][1]
