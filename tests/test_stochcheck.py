"""The frozen stochastic checker, and the coverage factor it rests on.

Step 1 of `docs/DIRECTION_STOCHASTIC.md`. Two halves: the coverage-factor answer
(what kappa means when the residual is INTRINSIC, and how the union bound over
patches enters), and the checker interface, which is scored against and is
therefore frozen -- `test_the_interface_is_frozen` is the lock.
"""
import numpy as np
import pytest

from lagh.certify import (KAPPA, Abstain, conjoin_determination, coverage_budget,
                          coverage_factor, determination, domain_qualifier)
from lagh.stochcheck import (ABSTAIN_REASONS, CONSUMER_QUANTITY, EXPECTATIONS,
                             PARTS, PROVENANCE, RECORD_STATUSES,
                             SUBMISSION_KINDS, Consumer, Coverage, Declaration,
                             Outcome, Submission, Task, TaskScore, component,
                             parse_component, rank, rank_key, score_task,
                             suite_totals, validate_submission, validate_task)

# ------------------------------------------------------------ coverage factor


def test_the_coverage_factor_inverts_a_union_bound_over_rows():
    """kappa(n, delta) = sqrt(2 ln(2n/delta)): the exponential martingale
    inequality per row, union-bounded over an exhaustive check, inverted for a
    declared false-abstain budget."""
    assert coverage_factor(300, 0.05) == pytest.approx(4.3342, abs=1e-3)
    assert coverage_factor(300, 0.01) == pytest.approx(4.6909, abs=1e-3)
    # monotone the right ways: more rows need more band, a tighter budget too
    ks = [coverage_factor(n, 0.05) for n in (30, 100, 300, 1000, 10000)]
    assert ks == sorted(ks)
    assert coverage_factor(300, 0.01) > coverage_factor(300, 0.05)
    with pytest.raises(ValueError):
        coverage_factor(300, 1.5)


def test_the_existing_kappa_is_the_same_number_for_a_stated_reason():
    """The continuity claimed when option (1) was registered, made numerical:
    kappa(300, 0.05) = 4.33 against the deterministic constant 4. What the
    constant did NOT have is the budget, and it is worth knowing -- 20% under the
    martingale bound at 300 rows, and no coverage at all past ~1491."""
    assert abs(coverage_factor(300, 0.05) - KAPPA) / KAPPA < 0.09
    assert coverage_budget(4, 300) == pytest.approx(0.2013, abs=1e-3)
    assert coverage_budget(4, 300, tail="gaussian") == pytest.approx(0.0190,
                                                                    abs=1e-3)
    # the bound is the price of assuming nothing about the residual's law
    assert (coverage_budget(4, 300)
            > 10 * coverage_budget(4, 300, tail="gaussian"))
    assert coverage_budget(4, 1490) < 1.0 <= coverage_budget(4, 1492)
    with pytest.raises(ValueError):
        coverage_budget(4, 300, tail="student")


def test_patch_count_is_a_winning_resource():
    """S7. kappa costs sqrt(log n) while alpha's exponent gains linearly in the
    held-out count, so exhaustive certification over MORE patches is cheap. This
    is the opposite of the usual multiple-comparisons intuition and it is the
    reason the stochastic regime is workable at all."""
    k_small, k_big = coverage_factor(300, 0.05), coverage_factor(10000, 0.05)
    assert k_big / k_small < 1.25          # 33x the rows, under 25% more band
    # ...and 33x the rows in the exponent of alpha <= |H| q^h
    assert (10000 - 4) / (300 - 4) > 30


def test_the_martingale_bound_holds_on_simulated_stochastic_integrals():
    """The inequality being registered, measured rather than trusted.

    M = sum b(X_k) sqrt(dt) Z_k with b STATE-DEPENDENT and driven by the same
    increments, so <M> is random and M is not a scaled Gaussian. Counting only
    the inequality's own event {|M| >= kappa sqrt(V) and <M> <= V}, the observed
    rate must sit under 2 exp(-kappa^2/2).
    """
    rng = np.random.default_rng(11)
    n_trials, n_steps, dt = 6000, 200, 1.0 / 200
    Z = rng.standard_normal((n_trials, n_steps))
    X = np.cumsum(0.3 * np.sqrt(dt) * Z, axis=1)
    b = 1.0 + 0.5 * X ** 2
    M = np.sum(b * np.sqrt(dt) * Z, axis=1)
    qv = np.sum(b ** 2 * dt, axis=1)
    V = 1.1                                     # a declared bound on <M>
    inside = qv <= V
    # the conditioning genuinely bites (it drops 12% of the paths), so the test
    # exercises the inequality's actual event rather than an unconditional tail
    assert 0.5 < inside.mean() < 0.95
    rates = []
    for kappa in (1.5, 2.0, 2.5, 3.0):
        rate = float(np.mean((np.abs(M) >= kappa * np.sqrt(V)) & inside))
        rates.append(rate)
        assert rate <= 2 * np.exp(-0.5 * kappa ** 2)
    assert rates[0] > 0 and rates[1] > 0        # the assertion above can fail
    # ...and the bound is roughly a decade loose, as it is at kappa = 4: that is
    # the price of assuming nothing about the law of the residual
    assert 2 * np.exp(-0.5 * 2.0 ** 2) / rates[1] > 5


# ---------------------------------------------------------------- the freeze

def test_the_interface_is_frozen():
    """It is scored against, so it cannot drift after tasks exist. This test IS
    the freeze: every tuple and every field name below is the interface."""
    assert PARTS == ("drift", "diffusion", "jump", "switch")
    assert EXPECTATIONS == ("exact", "interval", "abstain")
    assert SUBMISSION_KINDS == ("answer", "abstain")
    assert PROVENANCE == ("measured", "declared", "residual-derived")
    assert [o.value for o in Outcome] == [
        "exact", "covered", "abstained-correctly", "missed", "confident-wrong"]
    assert [c.value for c in Consumer] == ["drift-band", "diffusion-qv",
                                           "observation"]
    # the three stochastic abstain reasons live in the ONE abstain vocabulary
    for r in ("resolution", "coverage", "exploration"):
        assert r in ABSTAIN_REASONS
    assert set(RECORD_STATUSES) >= {"certified", "state", "conjoined",
                                    "structural-abstain", "structural"}
    fields = lambda c: tuple(c.__dataclass_fields__)          # noqa: E731
    assert fields(Task) == ("task_id", "level", "system", "state_dim", "truth",
                            "expectation", "sampling", "declarations",
                            "invariants", "tol_rel", "null", "null_reason")
    assert fields(Submission) == ("task_id", "kind", "record", "abstain",
                                  "reason_detail", "declarations", "coverage",
                                  "alpha_log10", "alpha_kind", "invariants",
                                  "law", "submission_id")
    assert fields(Coverage) == ("kappa", "delta", "n_rows", "n_disjoint",
                                "qv_provenance")
    assert fields(Declaration) == ("consumer", "value", "provenance",
                                   "quantity", "note")
    # the score sheet's shape is part of the interface too: a consumer reads it
    assert fields(TaskScore) == (
        "task_id", "components", "n_confident_wrong", "n_exact", "n_covered",
        "n_abstained_correctly", "n_missed", "n_resolved", "n_informative",
        "reach", "abstention", "invariants", "certificate_valid",
        "declaration_audit", "refusals", "exceeded_expectation", "domains",
        "notes")
    row = next(iter(score_task(_ou_task(), [_answer([("drift:x", -1.52, -1.48)])]
                               ).components.values()))
    assert set(row) == {"component", "domain", "truth", "in_library", "kind",
                        "lo", "hi", "resolved", "informative", "expectation",
                        "outcome", "note", "voided"}
    assert set(suite_totals([]).keys()) == {
        "n_tasks", "confident_wrong", "exact", "covered", "abstained_correctly",
        "missed", "resolved", "informative", "invariants",
        "invalid_certificates", "refused", "exceeded_expectation"}
    # each consumer states the QUANTITY it bounds, and no two are the same
    assert len(CONSUMER_QUANTITY) == len(Consumer)
    assert len(set(CONSUMER_QUANTITY.values())) == len(Consumer)


def test_component_names_round_trip_and_reject_junk():
    assert component("drift", "x") == "drift:x"
    assert component("diffusion", "x0*x1", 1) == "diffusion[1]:x0*x1"
    assert parse_component("diffusion[1]:x0*x1") == ("diffusion", 1, "x0*x1")
    assert parse_component("drift:x**3") == ("drift", None, "x**3")
    with pytest.raises(ValueError):
        component("dispersion", "x")
    with pytest.raises(ValueError):
        parse_component("drift-x")
    with pytest.raises(ValueError):
        parse_component("velocity:x")


# ------------------------------------------------------------------- fixtures

def _cov(**kw):
    d = dict(kappa=4.4, delta=0.05, n_rows=300, n_disjoint=60,
             qv_provenance="measured")
    d.update(kw)
    return Coverage(**d)


def _ou_task(**kw):
    """Level 0 OU: dX = -theta X dt + sigma dW, with a cubic term the generator
    does NOT have (the zero coefficient a null is made of)."""
    d = dict(
        task_id="L0-ou", level=0, system="ou", state_dim=1,
        truth={"drift:x": -1.5, "drift:x**3": 0.0, "diffusion:1": 0.4},
        expectation={"drift:x": "interval", "drift:x**3": "interval",
                     "diffusion:1": "abstain"},
        sampling={"dt": 1e-3, "T": 20.0, "n_traj": 8, "seed": 0},
        declarations=(Declaration(Consumer.OBSERVATION, 1e-3,
                                  provenance="declared"),),
        tol_rel=1e-9)
    d.update(kw)
    return Task(**d)


def _answer(entries, **kw):
    d = dict(task_id="L0-ou", kind="answer",
             record=determination(entries, status="certified"),
             coverage=_cov(), submission_id="a1")
    d.update(kw)
    return Submission(**d)


# ------------------------------------------------------------------- scoring

def test_interval_coverage_of_a_true_coefficient_is_the_scored_thing():
    """The registered scoring is per component and three-regime, which is why
    partial determination was a prerequisite: an interval that CONTAINS the truth
    is the good outcome and it is read straight off a `determination` record."""
    s = score_task(_ou_task(), [_answer([("drift:x", -1.52, -1.48),
                                        ("drift:x**3", -0.01, 0.01),
                                        ("diffusion:1", None, None)])])
    got = {k: v["outcome"] for k, v in s.components.items()}
    assert got["drift:x@all"] == "covered"
    assert got["drift:x**3@all"] == "covered"       # bounded, and around zero
    assert got["diffusion:1@all"] == "abstained-correctly"
    assert s.n_confident_wrong == 0
    assert s.n_covered == 2 and s.n_abstained_correctly == 1
    # RESOLVED is reach, not correctness: the linear drift is established as
    # present, the cubic one is honestly bounded without claiming presence
    assert s.n_resolved == 1
    assert s.reach == ("drift:x", "drift:x**3")


def test_a_certified_interval_that_excludes_the_truth_is_a_confident_wrong():
    s = score_task(_ou_task(), [_answer([("drift:x", -1.2, -1.1),
                                        ("drift:x**3", 0.0, 0.0),
                                        ("diffusion:1", None, None)])])
    assert s.components["drift:x@all"]["outcome"] == "confident-wrong"
    assert s.n_confident_wrong == 1
    assert "EXCLUDES" in s.components["drift:x@all"]["note"]


def test_an_exact_claim_is_scored_exactly():
    t = _ou_task(expectation={"drift:x": "exact", "drift:x**3": "interval",
                              "diffusion:1": "abstain"})
    ok = score_task(t, [_answer([("drift:x", -1.5, -1.5),
                                 ("drift:x**3", 0.0, 0.0),
                                 ("diffusion:1", None, None)])])
    assert ok.n_exact == 2 and ok.n_confident_wrong == 0
    wrong = score_task(t, [_answer([("drift:x", -1.5000001, -1.5000001)])])
    assert wrong.components["drift:x@all"]["outcome"] == "confident-wrong"


def test_silence_is_not_abstention():
    """The registered format says abstention is an EXPLICIT token with a
    structured reason. A component nobody mentioned is a MISS."""
    s = score_task(_ou_task(), [_answer([("drift:x", -1.52, -1.48)])])
    row = s.components["diffusion:1@-"]
    assert row["outcome"] == "missed"
    assert "silence is not abstention" in row["note"]
    assert s.n_abstained_correctly == 0 and s.n_missed == 2


def test_unconstrained_where_the_data_should_have_determined_it_is_a_miss():
    """A reach loss, never a wrong answer -- the distinction the whole ranking
    turns on."""
    s = score_task(_ou_task(), [_answer([("drift:x", None, None),
                                        ("drift:x**3", None, None),
                                        ("diffusion:1", None, None)])])
    assert s.n_missed == 2 and s.n_abstained_correctly == 1
    assert s.n_confident_wrong == 0


def test_an_explicit_abstention_scores_its_reason_against_the_registered_one():
    """S4: the Delta-t-unidentifiable drift must abstain with a RESOLUTION
    reason, not a wide interval."""
    t = _ou_task(task_id="L0-null-dt", null=True,
                 truth={"drift:x": 0.0, "diffusion:1": 0.0},
                 expectation={"drift:x": "abstain", "diffusion:1": "abstain"},
                 null_reason=Abstain.RESOLUTION.value)
    right = score_task(t, [Submission(task_id="L0-null-dt", kind="abstain",
                                      abstain="resolution",
                                      reason_detail="drift unidentifiable at dt")])
    assert right.abstention["reason_correct"] is True
    assert right.n_abstained_correctly == 2 and right.n_confident_wrong == 0
    assert right.abstention["accepted"] is True
    wrong_reason = score_task(t, [Submission(task_id="L0-null-dt",
                                             kind="abstain", abstain="noise")])
    assert wrong_reason.abstention["reason_correct"] is False
    assert wrong_reason.n_abstained_correctly == 2   # abstaining was still right
    # a reason OUTSIDE the frozen vocabulary is refused, and then the abstention
    # earns nothing -- "offered" alone must not read as credit
    junk = score_task(t, [Submission(task_id="L0-null-dt", kind="abstain",
                                     abstain="it looked hard")])
    assert junk.refusals and junk.abstention["accepted"] is False
    assert junk.n_abstained_correctly == 0 and junk.n_missed == 2


def test_an_abstention_may_carry_partial_determination():
    """lagh's primary output shape: "ABSTAIN[structural], and here is what every
    consistent law agrees on". Dropping the record on an abstain submission left
    the interface unable to score its own producer (found on Itô Level 0). The
    token now covers only the components the record does not mention."""
    t = _ou_task()
    sub = Submission(
        task_id="L0-ou", kind="abstain", abstain="structural",
        reason_detail="materially different laws certify",
        record=determination([("drift:x", -1.6, -1.4)], status="structural"),
        coverage=_cov(), submission_id="partial")
    assert validate_submission(t, sub) == []
    s = score_task(t, [sub])
    # the record's component is SCORED, not dropped
    assert s.components["drift:x@all"]["outcome"] == "covered"
    # ...and the token speaks for the rest, exactly once each
    assert s.components["diffusion:1@abstain"]["outcome"] == "abstained-correctly"
    assert "drift:x@abstain" not in s.components
    assert s.n_covered == 1 and s.n_abstained_correctly == 1
    assert s.abstention["reason"] == "structural"
    # interval content without a coverage statement is refused
    naked = Submission(task_id="L0-ou", kind="abstain", abstain="structural",
                       record=determination([("drift:x", -1.6, -1.4)],
                                            status="structural"))
    assert any("coverage statement" in r
               for r in validate_submission(t, naked))


def test_a_null_punishes_only_a_claim_of_presence():
    """Pure noise: the truth is all-zero. An interval AROUND zero is honest and
    scores covered; a resolved claim is a confident-wrong. The null needs no
    special case in the scorer, which is the sign the vocabulary is right."""
    t = _ou_task(task_id="L0-null-noise", null=True,
                 truth={"drift:x": 0.0, "diffusion:1": 0.0},
                 expectation={"drift:x": "abstain", "diffusion:1": "abstain"})
    honest = score_task(t, [_answer([("drift:x", -0.05, 0.05),
                                     ("diffusion:1", None, None)],
                                    task_id="L0-null-noise")])
    assert honest.n_confident_wrong == 0 and honest.n_covered == 1
    assert honest.exceeded_expectation == ["drift:x"]     # a finding, not a fine
    claimed = score_task(t, [_answer([("drift:x", 0.2, 0.4)],
                                     task_id="L0-null-noise")])
    assert claimed.n_confident_wrong == 1


def test_a_term_outside_the_declared_library_may_not_be_claimed_present():
    s = score_task(_ou_task(), [_answer([("drift:x", -1.52, -1.48),
                                        ("drift:sin(x)", 0.3, 0.5)])])
    assert s.components["drift:sin(x)@all"]["outcome"] == "confident-wrong"
    assert s.components["drift:sin(x)@all"]["in_library"] is False


def test_a_spelling_difference_cannot_manufacture_a_confident_wrong():
    """A false confident-wrong in the CHECKER would be worse than one in the
    instrument, so terms are compared canonically."""
    s = score_task(_ou_task(), [_answer([("drift:x*x*x", -0.01, 0.01)])])
    assert s.components["drift:x**3@all"]["outcome"] == "covered"
    assert s.components["drift:x**3@all"]["in_library"] is True


# ----------------------------------------------------------------- refusals

def test_a_run_that_cannot_state_its_coverage_refuses():
    """Registered with the certificate-kind decision, and enforced rather than
    exhorted."""
    bad = validate_submission(_ou_task(),
                              _answer([("drift:x", -1.52, -1.48)], coverage=None))
    assert any("cannot STATE its coverage" in r for r in bad)


def test_a_kappa_that_does_not_meet_its_own_budget_refuses():
    """This is what makes the coverage-factor answer load-bearing instead of
    decorative: kappa, delta and n_rows are checked against each other."""
    bad = validate_submission(
        _ou_task(), _answer([("drift:x", -1.52, -1.48)],
                            coverage=_cov(kappa=3.0)))
    assert any("does not meet the stated budget" in r for r in bad)
    assert validate_submission(
        _ou_task(), _answer([("drift:x", -1.52, -1.48)],
                            coverage=_cov(kappa=4.34))) == []


def test_one_magnitude_may_not_do_two_jobs():
    """The PDEBench field_err mistake, encoded as a refusal: three consumers want
    'noise intensity' in three different units."""
    bad = validate_submission(_ou_task(), _answer(
        [("drift:x", -1.52, -1.48)],
        declarations=(Declaration(Consumer.OBSERVATION, 2.75e-2),
                      Declaration(Consumer.DIFFUSION_QV, 2.75e-2))))
    assert any("SAME magnitude" in r for r in bad)
    twice = validate_submission(_ou_task(), _answer(
        [("drift:x", -1.52, -1.48)],
        declarations=(Declaration(Consumer.OBSERVATION, 1e-3),
                      Declaration(Consumer.OBSERVATION, 2e-3))))
    assert any("declared twice" in r for r in twice)
    # equal magnitudes are the mistake's SIGNATURE, not a proof of it: sigma_obs
    # and b^2 can coincide, so an attested coincidence is accepted rather than
    # falsely refused
    attested = validate_submission(_ou_task(), _answer(
        [("drift:x", -1.52, -1.48)],
        declarations=(Declaration(Consumer.OBSERVATION, 0.1,
                                  note="instrument precision, independently "
                                       "measured from replicates"),
                      Declaration(Consumer.DIFFUSION_QV, 0.1,
                                  note="b^2 at the reference state; equal to "
                                       "sigma_obs by coincidence"))))
    assert attested == []


def test_a_declaration_that_misstates_what_it_bounds_refuses():
    bad = validate_submission(_ou_task(), _answer(
        [("drift:x", -1.52, -1.48)],
        declarations=(Declaration(Consumer.OBSERVATION, 1e-3,
                                  quantity="the diffusion rate b^2"),)))
    assert any("states it bounds" in r for r in bad)


def test_a_band_derived_from_its_own_residual_refuses():
    """`errormodel.characterize_rows` returns those numbers LABELLED so they
    cannot travel; in a scored setting the label becomes a refusal."""
    bad = validate_submission(
        _ou_task(), _answer([("drift:x", -1.52, -1.48)],
                            coverage=_cov(qv_provenance="residual-derived")))
    assert any("circular" in r for r in bad)


def test_alpha_without_the_independence_discount_refuses():
    bad = validate_submission(
        _ou_task(), _answer([("drift:x", -1.52, -1.48)], alpha_log10=-40.0,
                            coverage=_cov(n_disjoint=None)))
    assert any("n_disjoint" in r for r in bad)


def test_a_refusal_voids_credit_but_not_exposure():
    """Otherwise 'forget the coverage line' would be a way to launder a wrong
    answer into a non-answer."""
    s = score_task(_ou_task(), [_answer([("drift:x", -1.2, -1.1),
                                        ("drift:x**3", -0.01, 0.01)],
                                       coverage=None)])
    assert s.refusals and s.n_confident_wrong == 1
    assert s.n_covered == 0                      # the good rows earn nothing
    assert s.components["drift:x**3@all#refused:a1"]["voided"] is True
    assert s.components["drift:x@all#refused:a1"]["voided"] is False


def test_a_refused_submission_cannot_overwrite_an_accepted_one():
    """Otherwise 'one good answer plus one refused answer' would lose the good
    answer's credit to the voided row that shares its key."""
    good = _answer([("drift:x", -1.52, -1.48)], submission_id="good")
    refused = _answer([("drift:x", -1.52, -1.48)], submission_id="bad",
                      coverage=None)
    s = score_task(_ou_task(), [good, refused])
    assert s.components["drift:x@all"]["voided"] is False
    assert s.components["drift:x@all#refused:bad"]["voided"] is True
    assert s.n_covered == 1 and s.n_confident_wrong == 0


def test_a_task_missing_a_registered_expectation_is_not_registrable():
    t = _ou_task(expectation={"drift:x": "interval"})
    assert any("no registered expectation" in p for p in validate_task(t))
    with pytest.raises(ValueError):
        score_task(t, [])
    assert any("nonzero true coefficient" in p
               for p in validate_task(_ou_task(null=True)))


# ------------------------------------------------------------------- domains

def test_a_domain_qualified_answer_is_scored_on_its_own_domain():
    """Darcy is the precedent: beta = 0.1000 in one conductivity phase and
    nowhere else. Narrowness is honest, so it is reported, never penalised."""
    hi = domain_qualifier("x > 0", coverage=0.55)
    lo = domain_qualifier("x <= 0", coverage=0.45)
    t = _ou_task()
    subs = [
        Submission("L0-ou", "answer", submission_id="hi", coverage=_cov(),
                   record=determination([("drift:x", -1.52, -1.48)],
                                        status="certified", qualifier=hi)),
        Submission("L0-ou", "answer", submission_id="lo", coverage=_cov(),
                   record=determination([("drift:x", -1.55, -1.45)],
                                        status="certified", qualifier=lo)),
    ]
    s = score_task(t, subs)
    assert s.components["drift:x@x > 0"]["outcome"] == "covered"
    assert s.components["drift:x@x <= 0"]["outcome"] == "covered"
    assert s.domains["x > 0"]["coverage"] == 0.55
    assert s.n_confident_wrong == 0
    # they were never conjoined across domains, so no refusal note was needed
    assert not any("would not conjoin" in n for n in s.notes)


def test_same_domain_records_are_conjoined_and_cannot_lose_a_covered_truth():
    """Intersection is sound in the direction that matters: if the truth is in
    both inputs it is in the intersection, so conjoining can only ever expose a
    claim that was already wrong."""
    q = domain_qualifier("x > 0", coverage=0.5)
    t = _ou_task()
    both = [Submission("L0-ou", "answer", submission_id=f"s{i}", coverage=_cov(),
                       record=determination([("drift:x", lo, hi)],
                                            status="certified", qualifier=q))
            for i, (lo, hi) in enumerate([(-1.6, -1.4), (-1.52, -1.48)])]
    s = score_task(t, both)
    assert s.components["drift:x@x > 0"]["outcome"] == "covered"
    assert (s.components["drift:x@x > 0"]["lo"],
            s.components["drift:x@x > 0"]["hi"]) == (-1.52, -1.48)
    clash = [both[0],
             Submission("L0-ou", "answer", submission_id="c", coverage=_cov(),
                        record=determination([("drift:x", -1.2, -1.1)],
                                             status="certified", qualifier=q))]
    sc = score_task(t, clash)
    assert any("CONTRADICT" in n for n in sc.notes)
    assert sc.n_confident_wrong == 1
    # the conjoined record is what certify.conjoin_determination reports
    assert conjoin_determination([both[0].record, clash[1].record]
                                 ).get("contradiction") == ["drift:x"]


# --------------------------------------------------------- audit and re-check

def test_the_declaration_audit_reports_the_ratio_and_flags_a_decade():
    """'Scan the declaration', made routine. PDEBench over-declared 3900x and it
    took a session to notice; here it is a column in the score sheet."""
    t = _ou_task()
    over = score_task(t, [_answer([("drift:x", -1.52, -1.48)], declarations=(
        Declaration(Consumer.OBSERVATION, 3.9),))])
    row = over.declaration_audit["observation"]
    assert row["true"] == 1e-3 and row["ratio"] == pytest.approx(3900.0)
    assert "OVER-declared" in row["flag"]
    under = score_task(t, [_answer([("drift:x", -1.52, -1.48)], declarations=(
        Declaration(Consumer.OBSERVATION, 1e-6),))])
    assert "DANGEROUS" in under.declaration_audit["observation"]["flag"]
    # the drift band depends on a patch family only the submitter knows, so the
    # task states no reference and the audit says so rather than inventing one
    none = score_task(t, [_answer([("drift:x", -1.52, -1.48)], declarations=(
        Declaration(Consumer.DRIFT_BAND, 1e-4, provenance="measured"),))])
    dr = none.declaration_audit["drift-band"]
    assert dr["true"] is None and dr["ratio"] is None
    assert "not auditable" in dr["note"] and dr["provenance"] == "measured"


def test_the_independent_checker_is_its_own_axis():
    """A claim that does not survive re-verification at the TASK's declarations is
    not promoted to a confident-wrong, and is not ignored either."""
    sub = _answer([("drift:x", -1.52, -1.48)])
    s = score_task(_ou_task(), [sub],
                   recheck=lambda x: {"valid": False, "note": "row 12 outside eps"})
    assert s.certificate_valid is False
    assert s.n_confident_wrong == 0 and s.n_covered == 1
    assert any("REJECTED" in n for n in s.notes)
    ok = score_task(_ou_task(), [sub], recheck=lambda x: {"valid": True})
    assert ok.certificate_valid is True
    assert score_task(_ou_task(), [sub]).certificate_valid is None


def test_invariants_are_scored_up_to_affine_reparametrization():
    """An invariant is defined only up to scale and offset, and the samples must
    span more than one trajectory or every constant matches every invariant."""
    t = _ou_task(invariants=("x**2 + y**2",))
    xs = np.linspace(0.5, 2.0, 12)
    vals = {"x**2 + y**2": xs, "3*(x**2 + y**2) - 7": 3 * xs - 7,
            "x*y": np.linspace(-1, 1, 12) ** 2}
    ev = vals.get
    hit = score_task(t, [_answer([("drift:x", -1.52, -1.48)],
                                 invariants=("3*(x**2 + y**2) - 7",))],
                     invariant_eval=ev)
    assert hit.invariants["n_hits"] == 1
    miss = score_task(t, [_answer([("drift:x", -1.52, -1.48)],
                                 invariants=("x*y",))], invariant_eval=ev)
    assert miss.invariants["n_hits"] == 0
    # with no hook, only symbolic equality is available and the report says so
    bare = score_task(t, [_answer([("drift:x", -1.52, -1.48)],
                                  invariants=("3*(x**2 + y**2) - 7",))])
    assert bare.invariants["n_hits"] == 0
    assert any("symbolic equality only" in n for n in bare.invariants["notes"])


# ------------------------------------------------------------------- ranking

def test_any_confident_wrong_dominates_the_ranking():
    """No amount of reach trades against one. This is the product invariant
    expressed as a sort key."""
    t = _ou_task()
    careful = score_task(t, [_answer([("drift:x", -1.6, -1.4),
                                      ("drift:x**3", None, None),
                                      ("diffusion:1", None, None)])])
    reachy = score_task(t, [_answer([("drift:x", -1.52, -1.48),
                                     ("drift:x**3", 0.0, 0.0),
                                     ("diffusion:1", 0.3, 0.5)])])
    assert reachy.n_covered + reachy.n_exact > careful.n_covered
    assert reachy.n_confident_wrong == 0
    wrong = score_task(t, [_answer([("drift:x", -1.52, -1.48),
                                    ("drift:x**3", 0.4, 0.6),
                                    ("diffusion:1", 0.3, 0.5)])])
    order = [name for name, _ in rank({"careful": [careful], "reachy": [reachy],
                                       "wrong": [wrong]})]
    assert order[-1] == "wrong"
    assert order == ["reachy", "careful", "wrong"]
    # an invalid certificate ranks below a correct abstention but above a wrong
    tot = suite_totals([reachy])
    assert rank_key(tot)[0] == 0
    assert rank_key(suite_totals([wrong]))[0] == 1


def test_a_vacuous_interval_cannot_top_the_ranking():
    """The hole a coverage-only key leaves: [-1e9, 1e9] everywhere is covered on
    everything with zero confident-wrongs. lagh's own vacuity doctrine, applied
    per component -- honest, and worth nothing."""
    t = _ou_task()
    vacuous = score_task(t, [_answer([("drift:x", -1e9, 1e9),
                                      ("drift:x**3", -1e9, 1e9),
                                      ("diffusion:1", -1e9, 1e9)])])
    assert vacuous.n_covered == 3 and vacuous.n_confident_wrong == 0
    assert vacuous.n_informative == 0
    assert all(not r["informative"] for r in vacuous.components.values())
    assert "VACUOUS" in vacuous.components["drift:x@all"]["note"]
    tight = score_task(t, [_answer([("drift:x", -1.52, -1.48),
                                    ("drift:x**3", -0.01, 0.01),
                                    ("diffusion:1", None, None)])])
    assert tight.n_informative == 2
    # the vacuous entrant has MORE covered components and still ranks lower
    assert suite_totals([vacuous])["covered"] > suite_totals([tight])["covered"]
    order = [n for n, _ in rank({"vacuous": [vacuous], "tight": [tight]})]
    assert order == ["tight", "vacuous"]
    # ...and a wrong answer still ranks below the vacuous one
    wrong = score_task(t, [_answer([("drift:x", -1.2, -1.1)])])
    assert [n for n, _ in rank({"vacuous": [vacuous], "wrong": [wrong]})] == \
        ["vacuous", "wrong"]


def test_determinations_come_before_exactness_in_the_key():
    """The one ordering choice inside the key that the registration left open:
    DETERMINATIONS first, exactness as the tie-break. Three exact claims must not
    beat two exact plus three tight intervals -- five components determined is
    more than three, and 'exact recovery' is a tie-break among equals."""
    fewer = {"confident_wrong": 0, "invalid_certificates": 0, "informative": 3,
             "exact": 3, "covered": 0, "abstained_correctly": 0,
             "invariants": 0, "missed": 0, "refused": 0}
    more = dict(fewer, informative=5, exact=2, covered=3)
    assert rank_key(more) < rank_key(fewer)
    tie = dict(fewer, informative=5, exact=5, covered=0)
    assert rank_key(tie) < rank_key(more)        # same determinations, more exact
