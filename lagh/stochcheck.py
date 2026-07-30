"""The FROZEN checker for the stochastic suite (docs/STOCHASTIC_CHECKER.md).

Step 1 of `docs/DIRECTION_STOCHASTIC.md`. It is scored against, so it is frozen
BEFORE any Level 0 task exists -- a checker that moves after the tasks arrive is
a scoring key someone can fit.

Three commitments, each of them the cash-in of something measured elsewhere:

1. **Component claims travel ONLY in a `certify.determination` record.** That is
   why partial determination was a prerequisite rather than cleanup: the
   registered scoring is per component and three-regime (exact / interval /
   abstain), and `certified: bool` plus an abstain enum cannot express it.
   Interval COVERAGE of a true coefficient is then a containment test on a record
   this program already emits, from `engine`, `pdesystem` and `statecert` alike.

2. **A declared magnitude names its CONSUMER.** PDEBench fed one `field_err` to a
   weak-form band and a pointwise forecast check and each got the other's number
   (docs/DIRECTION_ERROR_PROVENANCE.md). Here at least three consumers want
   "noise intensity" and they want three different QUANTITIES in three different
   units, so the interface makes the mapping explicit and refuses a submission
   that gives one magnitude two jobs.

3. **A run that cannot state its coverage refuses.** Registered with the
   certificate-kind decision. It is enforced, not exhorted: an answer whose
   stated kappa does not meet its own stated budget at its own row count
   (`certify.coverage_factor`) is REFUSED.

The scoring vocabulary is five outcomes and one of them dominates:

    exact / covered / abstained-correctly / missed / CONFIDENT-WRONG

`missed` is a reach loss -- the instrument left something on the table. Only
`confident-wrong` is a failure of the product invariant, and it dominates the
ranking, so no amount of reach buys one off.

Coverage alone would be gameable -- [-1e9, 1e9] on every component is covered
everywhere with zero confident-wrongs -- so `certify.vacuous`'s doctrine is
applied per component: an interval that excludes nothing is honest, and worth
nothing (`_informative`, and `rank_key` puts determinations ahead of coverage).

Three rules worth stating up front because they are the ones that could have been
got wrong in the convenient direction:

* **Silence is not abstention.** A component the record never mentions scores as
  `missed`, never as a correct abstention. The registered format says abstention
  is an EXPLICIT token with a structured reason.
* **A refusal does not launder a wrong claim.** A submission refused on
  procedure (no coverage, a circular declaration) still has its component claims
  checked for confident-wrongs; what it loses is credit, not exposure.
* **An expectation is a prediction, not a ceiling.** Beating one is recorded as a
  finding about the expectation; it earns no penalty, and it shields nothing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np
import sympy as sp

from .certify import (MACHINE_REL, TAU, Abstain, conjoin_determination,
                      coverage_factor)

# ---------------------------------------------------------------- vocabularies
# Every tuple below is part of the FROZEN interface. `tests/test_stochcheck.py`
# asserts them verbatim; changing one is changing the scoring key.

PARTS = ("drift", "diffusion", "jump", "switch")
EXPECTATIONS = ("exact", "interval", "abstain")
SUBMISSION_KINDS = ("answer", "abstain")
# `determination(status=)` values a submission may carry: what PRODUCED the
# record. A range from one certified law and a range over a certifying SET are
# different claims (certify.determination) and the checker keeps them labelled.
RECORD_STATUSES = ("certified", "state", "conjoined",
                   "structural-abstain") + tuple(a.value for a in Abstain)
# The abstention reasons a stochastic run may give. The three added for this
# layer (resolution / coverage / exploration) live in `certify.Abstain` so the
# vocabulary stays enumerable by a consumer.
ABSTAIN_REASONS = tuple(a.value for a in Abstain)
PROVENANCE = ("measured", "declared", "residual-derived")
# Over/under-declaration flags for the audit. Under-declaring is the DANGEROUS
# direction (it admits impostors); over-declaring merely loses laws -- but it
# loses them by 3900x if nobody looks, so both are reported and a decade is
# flagged either way ("scan the declaration").
DECLARATION_FLAG = 10.0


class Consumer(str, Enum):
    """WHO consumes a declared noise magnitude. One job each, and the units are
    the argument that they cannot share a number."""
    DRIFT_BAND = "drift-band"
    DIFFUSION_QV = "diffusion-qv"
    OBSERVATION = "observation"


CONSUMER_QUANTITY = {
    Consumer.DRIFT_BAND: (
        "an upper bound on the martingale variance of ONE weak-form row, "
        "<M_p> = int phi_p^2 f'(X)^2 b(X)^2 dt; units (row target)^2; consumed "
        "as the band eps_p = kappa*sqrt(<M_p>)"),
    Consumer.DIFFUSION_QV: (
        "the diffusion RATE b^2(x) at the task's stated reference state; units "
        "[state]^2/time; consumed as realized quadratic variation per unit time"),
    Consumer.OBSERVATION: (
        "the measurement error on the observed state, sigma_obs; units [state]; "
        "consumed in the L2 channel, where it cancels within one realization"),
}


class Outcome(str, Enum):
    EXACT = "exact"                        # exact claim, matches the truth
    COVERED = "covered"                    # interval contains the truth
    ABSTAINED_CORRECTLY = "abstained-correctly"
    MISSED = "missed"                      # reach loss, not a wrong answer
    CONFIDENT_WRONG = "confident-wrong"    # dominates the ranking


# --------------------------------------------------------------- the interface

@dataclass(frozen=True)
class Declaration:
    """One declared magnitude, bound to the consumer it belongs to.

    `quantity` is optional and is a CHECK, not an input: when a submitter states
    what they believe the number bounds, the checker compares it against the
    frozen text for that consumer. A submitter who thinks their sigma bounds
    something else is then caught here rather than four orders of magnitude
    later.
    """
    consumer: Consumer
    value: float
    provenance: str = "declared"
    quantity: str = ""
    note: str = ""


@dataclass(frozen=True)
class Coverage:
    """The stochastic certificate's coverage statement.

    `kappa`/`delta`/`n_rows` are the coverage-factor triple: the band is
    kappa*sqrt(<M_p>) per row, exhaustive over n_rows, at a false-abstain budget
    delta. `n_disjoint` is the number of patches with DISJOINT support -- the
    honest held-out count for alpha, because q^h assumes independence and
    overlapping patches share the driving path. `qv_provenance` says where the
    band's SCALE came from; it may not come from the candidate.
    """
    kappa: float
    delta: float
    n_rows: int
    n_disjoint: int | None = None
    qv_provenance: str = "measured"


@dataclass(frozen=True)
class Task:
    """One stochastic system as a scoreable task.

    `truth` is the generating law, keyed by component name (see `component`),
    INCLUDING the zero coefficients -- the declared library is part of the task,
    and a term whose true coefficient is zero is what a null is made of.

    `expectation` is a REGISTERED PREDICTION per component, from `EXPECTATIONS`:
    what a calibrated instrument should manage at this sampling regime. It is not
    a ceiling. Beating it (certifying where it said abstain) is recorded as
    `exceeded_expectation` and is a finding about the expectation, never a
    penalty; and it never shields -- an answer that beats the expectation and
    misses the truth is still a confident-wrong.

    `declarations` are the task's OWN reference magnitudes, per consumer, and are
    what makes the submission's declarations auditable. A consumer the task
    cannot pin (the drift band depends on a patch family only the submitter
    knows) is simply absent, and the audit says so instead of inventing a ratio.
    """
    task_id: str
    level: int
    system: str
    state_dim: int
    truth: dict
    expectation: dict
    sampling: dict = field(default_factory=dict)
    declarations: tuple = ()
    invariants: tuple = ()
    tol_rel: float = 1e-9
    null: bool = False
    null_reason: str = ""


@dataclass(frozen=True)
class Submission:
    """One answer or one explicit abstention, for one task, on one domain.

    A domain-qualified answer is a legitimate partial answer (Darcy is the live
    precedent: beta = 0.1000 in one conductivity phase and nowhere else), so
    several submissions per task are allowed and each is scored ON ITS OWN
    DOMAIN. Records sharing a domain are conjoined first, which is where a
    same-domain contradiction surfaces.
    """
    task_id: str
    kind: str
    record: dict | None = None          # a certify.determination record
    abstain: str | None = None          # a reason from ABSTAIN_REASONS
    reason_detail: str = ""
    declarations: tuple = ()
    coverage: Coverage | None = None
    alpha_log10: float | None = None
    alpha_kind: str = ""
    invariants: tuple = ()
    law: str = ""                       # human-readable; the RECORD is scored
    submission_id: str = ""


@dataclass
class TaskScore:
    task_id: str
    components: dict = field(default_factory=dict)
    n_confident_wrong: int = 0
    n_exact: int = 0
    n_covered: int = 0
    n_abstained_correctly: int = 0
    n_missed: int = 0
    n_resolved: int = 0                 # reach: right AND established as present
    n_informative: int = 0              # right AND not vacuous at the task scale
    reach: tuple = ()                   # distinct components determined somewhere
    abstention: dict = field(default_factory=dict)
    invariants: dict = field(default_factory=dict)
    certificate_valid: bool | None = None
    declaration_audit: dict = field(default_factory=dict)
    refusals: list = field(default_factory=list)
    exceeded_expectation: list = field(default_factory=list)
    domains: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def one_line(self) -> str:
        if self.n_confident_wrong:
            return (f"{self.task_id}: {self.n_confident_wrong} CONFIDENT-WRONG "
                    f"(dominates) | exact {self.n_exact} covered {self.n_covered}")
        return (f"{self.task_id}: exact {self.n_exact} covered {self.n_covered} "
                f"abstained-ok {self.n_abstained_correctly} missed "
                f"{self.n_missed} | 0 confident-wrong")


# ------------------------------------------------------------ component naming

def component(part: str, term: str, index: int | None = None) -> str:
    """The frozen component name: `drift:x`, `diffusion[1]:x0*x1`.

    A function rather than a convention in prose, so a task author and a
    submitter cannot spell it differently.
    """
    if part not in PARTS:
        raise ValueError(f"part must be one of {PARTS}, got {part!r}")
    head = part if index is None else f"{part}[{int(index)}]"
    return f"{head}:{term}"


def parse_component(name: str) -> tuple:
    """(part, index, term). Raises on anything the convention cannot read."""
    head, _, term = str(name).partition(":")
    if not term:
        raise ValueError(f"component {name!r} has no ':<term>'")
    idx = None
    if head.endswith("]") and "[" in head:
        head, _, rest = head.partition("[")
        idx = int(rest[:-1])
    if head not in PARTS:
        raise ValueError(f"component {name!r} names no part in {PARTS}")
    return head, idx, term


def _canon(name: str) -> str:
    """Canonical component key, with the TERM normalized through sympy.

    `x*x` and `x**2` are the same column, and a spelling difference must not
    manufacture a confident-wrong in the CHECKER -- which would be worse than one
    in the instrument. Unsympifiable terms fall back to their raw text.
    """
    part, idx, term = parse_component(name)
    try:
        t = sp.srepr(sp.sympify(term))
    except Exception:                                          # noqa: BLE001
        t = term.strip()
    return f"{part}[{idx}]:{t}"


# ----------------------------------------------------------------- validation

def validate_task(task: Task) -> list:
    """Authoring errors in a task. A task that fails this is not registered."""
    bad = []
    if task.level not in (0, 1, 2, 3):
        bad.append(f"level {task.level} outside the registered curriculum 0-3")
    if not task.truth:
        bad.append("no truth: a task must declare its library, zeros included")
    for name in task.truth:
        try:
            parse_component(name)
        except ValueError as e:
            bad.append(f"truth key: {e}")
    missing = [n for n in task.truth if n not in task.expectation]
    if missing:
        bad.append(f"no registered expectation for {sorted(missing)[:4]} -- an "
                   "expectation is a prediction and must exist before any run")
    for n, e in task.expectation.items():
        if e not in EXPECTATIONS:
            bad.append(f"expectation {e!r} for {n} outside {EXPECTATIONS}")
    if task.null and task.null_reason and task.null_reason not in ABSTAIN_REASONS:
        bad.append(f"null_reason {task.null_reason!r} outside {ABSTAIN_REASONS}")
    if task.null and any(v for v in task.truth.values()):
        bad.append("a null task declares a nonzero true coefficient")
    for d in task.declarations:
        if not isinstance(d, Declaration):
            bad.append("task declarations must be Declaration instances")
    return bad


def validate_submission(task: Task, sub: Submission) -> list:
    """Why the checker refuses this submission. Empty means accepted.

    A refusal voids CREDIT, not exposure: `score_task` still checks the refused
    record's claims for confident-wrongs.
    """
    bad = []
    if sub.task_id != task.task_id:
        bad.append(f"task_id {sub.task_id!r} != {task.task_id!r}")
    if sub.kind not in SUBMISSION_KINDS:
        bad.append(f"kind {sub.kind!r} outside {SUBMISSION_KINDS}")
    if sub.kind == "abstain":
        if sub.abstain not in ABSTAIN_REASONS:
            bad.append(f"abstention reason {sub.abstain!r} outside the frozen "
                       f"vocabulary {ABSTAIN_REASONS}")
        # an abstention carrying INTERVAL content is making an interval claim, and
        # any interval claim states its coverage
        if (isinstance(sub.record, dict) and sub.record.get("components")
                and sub.coverage is None):
            bad.append("this abstention carries partial-determination components, "
                       "which are interval claims at a declared band, so it needs "
                       "a coverage statement like any other interval claim")
    if sub.kind == "answer":
        if not isinstance(sub.record, dict) or not sub.record:
            bad.append("an answer carries its component claims in a "
                       "certify.determination record; there is no other channel")
        else:
            st = sub.record.get("status")
            if st not in RECORD_STATUSES:
                bad.append(f"record status {st!r} outside {RECORD_STATUSES}")
            for n in sub.record.get("components", {}):
                try:
                    parse_component(n)
                except ValueError as e:
                    bad.append(f"record component: {e}")
        c = sub.coverage
        if c is None:
            bad.append("no coverage: a run that cannot STATE its coverage "
                       "refuses (registered with the certificate-kind decision)")
        else:
            if not (0.0 < c.delta < 1.0):
                bad.append(f"delta {c.delta} is not a probability")
            elif c.kappa + 1e-9 < coverage_factor(c.n_rows, c.delta):
                bad.append(
                    f"kappa {c.kappa:.4g} does not meet the stated budget: an "
                    f"exhaustive check over {c.n_rows} rows at delta="
                    f"{c.delta:g} needs kappa >= "
                    f"{coverage_factor(c.n_rows, c.delta):.4g}")
            if c.qv_provenance not in PROVENANCE:
                bad.append(f"qv_provenance {c.qv_provenance!r} outside {PROVENANCE}")
            elif c.qv_provenance == "residual-derived":
                bad.append("the band's scale is derived from the residual of the "
                           "law it bands: circular (errormodel returns those "
                           "numbers labelled precisely so they do not travel)")
            if sub.alpha_log10 is not None and not c.n_disjoint:
                bad.append("alpha is reported without n_disjoint: q^h assumes "
                           "independent held-out rows and overlapping patches "
                           "share the driving path, so the discount must be stated")
            if c.n_disjoint and c.n_disjoint > c.n_rows:
                bad.append(f"n_disjoint {c.n_disjoint} exceeds n_rows {c.n_rows}")
    seen = {}
    for d in sub.declarations:
        if not isinstance(d, Declaration):
            bad.append("declarations must be Declaration instances")
            continue
        if d.consumer in seen:
            bad.append(f"consumer {d.consumer.value} declared twice")
        seen[d.consumer] = d
        if d.provenance not in PROVENANCE:
            bad.append(f"provenance {d.provenance!r} outside {PROVENANCE}")
        elif d.provenance == "residual-derived":
            bad.append(f"the {d.consumer.value} magnitude is derived from the "
                       "residual of the law it will band: circular")
        if d.quantity and d.quantity.strip() != CONSUMER_QUANTITY[d.consumer].strip():
            bad.append(f"the {d.consumer.value} declaration states it bounds "
                       f"{d.quantity!r}; that consumer bounds "
                       f"{CONSUMER_QUANTITY[d.consumer]!r}")
    # Equal magnitudes across consumers is the SIGNATURE of the PDEBench mistake
    # (one field_err reused for two jobs), but it is a signature and not a proof:
    # sigma_obs and b^2 can coincide. So it refuses unless BOTH declarations carry
    # a note -- an explicit attestation that the coincidence is real. The friction
    # stays where it belongs and the check cannot be tripped silently.
    ds = list(seen.values())
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            a, b = ds[i], ds[j]
            if a.value == b.value and a.value != 0 and not (a.note and b.note):
                bad.append(
                    f"{a.consumer.value} and {b.consumer.value} were given the "
                    f"SAME magnitude {a.value:g}; they bound quantities in "
                    "different units, so one of them likely has the other's "
                    "number (the PDEBench field_err mistake). If the coincidence "
                    "is real, say so in both declarations' `note`")
    return bad


# -------------------------------------------------------------------- scoring

def _score_one(truth_v, expect, rec, tol_rel):
    """(outcome, note) for one component claim against one true coefficient."""
    kind = rec.get("kind")
    lo, hi = rec.get("lo"), rec.get("hi")
    if kind == "unconstrained":
        if expect == "abstain":
            return Outcome.ABSTAINED_CORRECTLY, "unconstrained, as expected"
        return Outcome.MISSED, ("reported unconstrained where the expectation "
                                f"was {expect!r}: reach loss, not a wrong answer")
    if lo is None or hi is None:
        return Outcome.MISSED, f"kind {kind!r} with no bounds"
    if kind == "exact":
        if abs(lo - truth_v) <= tol_rel * max(abs(truth_v), 1.0):
            return Outcome.EXACT, "exact claim matches the truth"
        return Outcome.CONFIDENT_WRONG, (f"exact claim {lo!r} against truth "
                                         f"{truth_v!r}")
    if kind == "interval":
        slack = MACHINE_REL * max(abs(truth_v), abs(lo), abs(hi))
        if lo - slack <= truth_v <= hi + slack:
            return Outcome.COVERED, f"interval [{lo:g}, {hi:g}] covers the truth"
        return Outcome.CONFIDENT_WRONG, (f"certified interval [{lo:g}, {hi:g}] "
                                         f"EXCLUDES the truth {truth_v!r}")
    return Outcome.MISSED, f"unreadable component kind {kind!r}"


def _informative(rec, scale) -> bool:
    """Is this component claim worth anything -- lagh's VACUITY test, per
    component.

    `certify.vacuous` refuses a certificate when the ZERO law also certifies:
    the band swallowed the signal and the certificate proves nothing. An interval
    has the same failure mode and it is the hole a ranking on COVERAGE alone
    leaves open -- submit [-1e9, 1e9] for every component and you are covered
    everywhere with zero confident-wrongs. So an interval counts as informative
    only when its span is under 2*TAU of the task's own coefficient scale (TAU is
    the material-difference scale the whole engine uses), and an exact claim
    always is. Reported per row and counted, never scored as WRONG: a vacuous
    interval is honest, it is just worth nothing.
    """
    if rec.get("kind") == "exact":
        return True
    if rec.get("kind") != "interval":
        return False
    lo, hi = rec.get("lo"), rec.get("hi")
    if lo is None or hi is None:
        return False
    return bool((hi - lo) <= 2.0 * TAU * scale)


def _rows(record, truth, expect, tol_rel, domain, names=None, scale=1.0, *,
          voided=False):
    """Score every component a record mentions. Returns {key: row}.

    `names` maps a canonical key back to the TASK's spelling, so every row -- and
    therefore `reach` -- reports one name per component however the submitter
    spelled it. `scale` is the task's coefficient scale, for the vacuity test.
    """
    out = {}
    names = names or {}
    for name, rec in record.get("components", {}).items():
        try:
            key = _canon(name)
        except ValueError:
            continue
        known = key in truth
        tv = truth.get(key, 0.0)
        label = names.get(key, name)
        if known:
            outcome, note = _score_one(tv, expect.get(key, "interval"), rec,
                                       tol_rel)
        elif rec.get("resolved"):
            outcome = Outcome.CONFIDENT_WRONG
            note = ("a term OUTSIDE the task's declared library is claimed "
                    "present (resolved); the generator does not have it")
        else:
            outcome = Outcome.MISSED
            note = ("outside the declared library, and not claimed present -- "
                    "no penalty beyond this note")
        info = _informative(rec, scale)
        if outcome is Outcome.COVERED and not info:
            note += ("; but the interval is VACUOUS at this task's coefficient "
                     "scale -- it excludes nothing, so it earns coverage and no "
                     "determination")
        out[f"{label}@{domain or 'all'}"] = {
            "component": label, "domain": domain or None, "truth": tv,
            "in_library": known, "kind": rec.get("kind"),
            "lo": rec.get("lo"), "hi": rec.get("hi"),
            "resolved": bool(rec.get("resolved")), "informative": info,
            "expectation": expect.get(key), "outcome": outcome.value,
            "note": note,
            # a refused submission keeps its EXPOSURE and loses its credit
            "voided": bool(voided and outcome is not Outcome.CONFIDENT_WRONG),
        }
    return out


def _audit(task: Task, subs) -> dict:
    """Declared vs true magnitude, per consumer. The 'scan the declaration' step,
    made routine instead of discovered a session later."""
    ref = {d.consumer: d for d in task.declarations}
    out = {}
    for sub in subs:
        for d in sub.declarations:
            row = {"declared": float(d.value), "provenance": d.provenance,
                   "bounds": CONSUMER_QUANTITY[d.consumer], "true": None,
                   "ratio": None, "flag": None}
            t = ref.get(d.consumer)
            if t is None:
                row["note"] = ("the task states no reference for this consumer, "
                               "so the magnitude is not auditable here; its "
                               "provenance still is")
            else:
                row["true"] = float(t.value)
                if t.value:
                    r = float(d.value) / float(t.value)
                    row["ratio"] = r
                    if r >= DECLARATION_FLAG:
                        row["flag"] = (f"OVER-declared {r:.3g}x -- the band is "
                                       "needlessly wide and laws will be lost")
                    elif r <= 1.0 / DECLARATION_FLAG:
                        row["flag"] = (f"UNDER-declared {r:.3g}x -- the "
                                       "DANGEROUS direction: a band this tight "
                                       "admits impostors")
            prev = out.get(d.consumer.value)
            if prev is not None and prev["declared"] != row["declared"]:
                # two submissions to one task declared the same consumer
                # differently. Legal per submission, and worth surfacing: the
                # audit is not free to silently keep whichever came last.
                row["conflict"] = (f"another submission declared "
                                   f"{prev['declared']:g} for this consumer")
            out[d.consumer.value] = row
    return out


def _invariant_hits(claimed, true, invariant_eval, tol_rel):
    """Which claimed invariants match a true one, up to affine reparametrization.

    An invariant is defined only up to scale and offset, so equality is the wrong
    test. With `invariant_eval` the check is numeric: the claimed quantity must be
    a NONCONSTANT affine function of a true one over the task's samples. Without
    it, only normalized symbolic equality is available and the report says so.

    The hook's samples must span MORE THAN ONE TRAJECTORY. Within one trajectory
    every invariant is constant and every constant matches every invariant.
    """
    hits, notes = [], []
    if invariant_eval is None:
        for c in claimed:
            for t in true:
                try:
                    if sp.simplify(sp.sympify(c) - sp.sympify(t)) == 0:
                        hits.append((c, t, "symbolic equality"))
                        break
                except Exception:                              # noqa: BLE001
                    continue
        notes.append("no invariant_eval hook: scored by symbolic equality only, "
                     "so an invariant right up to scale or offset reads as a miss")
        return hits, notes
    for c in claimed:
        v = invariant_eval(c)
        if v is None:
            notes.append(f"{c!r} did not evaluate")
            continue
        v = np.asarray(v, float).ravel()
        for t in true:
            w = invariant_eval(t)
            if w is None:
                continue
            w = np.asarray(w, float).ravel()
            if v.size != w.size or v.size < 4 or np.ptp(w) == 0:
                notes.append("samples do not span more than one trajectory (the "
                             "true invariant is constant on them), so no "
                             "invariant claim can be decided")
                continue
            A = np.column_stack([np.ones_like(w), w])
            coef, *_ = np.linalg.lstsq(A, v, rcond=None)
            resid = float(np.max(np.abs(v - A @ coef)))
            if abs(coef[1]) > 0 and resid <= tol_rel * max(np.ptp(v), 1e-300):
                hits.append((c, t, f"affine in the true invariant "
                                   f"(slope {coef[1]:.6g}, max resid {resid:.3g})"))
                break
    return hits, notes


def score_task(task: Task, submissions, *, recheck=None, invariant_eval=None
               ) -> TaskScore:
    """Score every submission for one task.

    `recheck(submission) -> {"valid": bool, "note": str, "alpha_log10": float}`
    is the INDEPENDENT checker: the harness holds the task's data and re-verifies
    the claim at the TASK's declarations rather than the submitter's. An invalid
    certificate is a procedural failure scored on its own axis -- it is not
    silently promoted to a confident-wrong, and it is not ignored either.

    `invariant_eval(expr_str) -> array | None` evaluates a claimed invariant on
    the task's own samples (see `_invariant_hits`). Both hooks are optional, and
    what could not be checked is reported rather than assumed.
    """
    problems = validate_task(task)
    if problems:
        raise ValueError(f"task {task.task_id} is not registrable: {problems}")
    truth = {_canon(k): float(v) for k, v in task.truth.items()}
    expect = {_canon(k): v for k, v in task.expectation.items()}
    names = {_canon(k): k for k in task.truth}
    # the task's own coefficient scale, for the per-component vacuity test. All
    # zero (a null) leaves no scale in the truth, so the unit scale is used and
    # said so in the note the rows carry.
    scale = max([abs(v) for v in truth.values()] or [0.0]) or 1.0
    score = TaskScore(task_id=task.task_id)

    subs = list(submissions)
    groups: dict = {}
    for sub in subs:
        bad = validate_submission(task, sub)
        if bad:
            score.refusals.append({"submission": sub.submission_id or None,
                                   "kind": sub.kind, "reasons": bad})
        dom = ""
        if isinstance(sub.record, dict):
            dom = (sub.record.get("qualifier") or {}).get("predicate", "") or ""
        groups.setdefault(dom, []).append((sub, bool(bad)))

    mentioned: set = set()
    blanket_abstain = None
    for dom, members in sorted(groups.items()):
        # AN ABSTENTION MAY CARRY PARTIAL DETERMINATION, and its record is scored.
        # "ABSTAIN[structural], and here is what every consistent law agrees on" is
        # lagh's primary output shape -- the whole point of step 0 -- and dropping
        # the record on an abstain submission made the interface unable to score
        # its own producer (measured on Itô Level 0). The token then covers only
        # the components the record does not mention.
        accepted = [s.record for s, bad in members
                    if not bad and isinstance(s.record, dict) and s.record]
        if len(accepted) > 1:
            rec = conjoin_determination(accepted)
            if rec.get("status") == "refused":
                score.notes.append(
                    f"domain {dom or 'all'}: {len(accepted)} records would not "
                    f"conjoin ({rec.get('refusal')})")
                rec = None
            elif rec.get("contradiction"):
                score.notes.append(
                    f"domain {dom or 'all'}: same-domain records CONTRADICT on "
                    f"{rec['contradiction']} -- they cannot all be right, and "
                    "the intersection can only exclude a truth an input already did")
        elif accepted:
            rec = accepted[0]
        else:
            rec = None
        if rec is not None:
            rows = _rows(rec, truth, expect, task.tol_rel, dom, names, scale)
            score.components.update(rows)
            mentioned |= {_canon(r["component"]) for r in rows.values()
                          if r["in_library"]}
            q = rec.get("qualifier") or {}
            if dom:
                score.domains[dom] = {"coverage": q.get("coverage"),
                                      "n_components": len(rows)}
        for s, bad in members:
            if bad and isinstance(s.record, dict) and s.record:
                # suffixed keys: a refused submission must not overwrite an
                # ACCEPTED row for the same component and domain, which would
                # turn "one good answer plus one refused one" into a loss of the
                # good one's credit
                tag = s.submission_id or "refused"
                score.components.update({
                    f"{k}#refused:{tag}": v for k, v in
                    _rows(s.record, truth, expect, task.tol_rel, dom, names,
                          scale, voided=True).items()})
            if s.kind == "abstain" and not bad and blanket_abstain is None:
                blanket_abstain = s

    if blanket_abstain is not None:
        for key, tv in truth.items():
            if key in mentioned:
                continue           # the abstention's own record already spoke
            e = expect.get(key, "interval")
            ok = e == "abstain" or task.null
            score.components[f"{names[key]}@abstain"] = {
                "component": names[key], "domain": None, "truth": tv,
                "in_library": True, "kind": "abstain-token", "lo": None,
                "hi": None, "resolved": False, "informative": False,
                "expectation": e,
                "outcome": (Outcome.ABSTAINED_CORRECTLY if ok
                            else Outcome.MISSED).value,
                "note": ("explicit abstention with reason "
                         f"{blanket_abstain.abstain!r}"), "voided": False}
        mentioned |= set(truth)

    # SILENCE IS NOT ABSTENTION: anything no accepted record and no abstention
    # token ever mentioned is a miss, and says why it is not scored as one.
    for key, tv in truth.items():
        if key in mentioned:
            continue
        score.components[f"{names[key]}@-"] = {
            "component": names[key], "domain": None, "truth": tv,
            "in_library": True,
            "kind": None, "lo": None, "hi": None, "resolved": False,
            "informative": False,
            "expectation": expect.get(key), "outcome": Outcome.MISSED.value,
            "note": ("never mentioned by any ACCEPTED submission (a refused one "
                     "does not count); silence is not abstention, so it is a "
                     "miss rather than a correct abstention"), "voided": False}

    counts = {o: 0 for o in Outcome}
    reach = set()
    for row in score.components.values():
        o = Outcome(row["outcome"])
        if o is Outcome.CONFIDENT_WRONG or not row["voided"]:
            counts[o] += 1
        if o in (Outcome.EXACT, Outcome.COVERED) and not row["voided"]:
            reach.add(row["component"])
            if row["resolved"] and row["truth"] != 0.0:
                score.n_resolved += 1
            if row["informative"]:
                score.n_informative += 1
        if (row["expectation"] == "abstain"
                and o in (Outcome.EXACT, Outcome.COVERED) and not row["voided"]):
            score.exceeded_expectation.append(row["component"])
    score.n_confident_wrong = counts[Outcome.CONFIDENT_WRONG]
    score.n_exact = counts[Outcome.EXACT]
    score.n_covered = counts[Outcome.COVERED]
    score.n_abstained_correctly = counts[Outcome.ABSTAINED_CORRECTLY]
    score.n_missed = counts[Outcome.MISSED]
    score.reach = tuple(sorted(reach))

    tokens = [s for s in subs if s.kind == "abstain"]
    if tokens:
        t = blanket_abstain or tokens[0]
        correct = None
        if task.null_reason:
            correct = bool(t.abstain == task.null_reason)
        score.abstention = {
            # `accepted` matters: an abstention whose reason is outside the frozen
            # vocabulary was REFUSED, so it earned no abstained-correctly rows,
            # and "offered" alone would read as credit it did not get
            "offered": True, "accepted": blanket_abstain is not None,
            "reason": t.abstain, "detail": t.reason_detail,
            "registered_reason": task.null_reason or None,
            "reason_correct": correct,
            "note": ("S4: a Delta-t-unidentifiable drift must abstain with a "
                     "RESOLUTION reason, not a wide interval")}
    else:
        score.abstention = {"offered": False}

    claimed_inv = [i for s in subs for i in s.invariants]
    hits, inv_notes = _invariant_hits(claimed_inv, list(task.invariants),
                                      invariant_eval, task.tol_rel)
    score.invariants = {"claimed": claimed_inv, "true": list(task.invariants),
                        "hits": hits, "n_hits": len(hits), "notes": inv_notes}

    score.declaration_audit = _audit(task, subs)

    if recheck is not None:
        verdicts = []
        for s in subs:
            if s.kind != "answer":
                continue
            r = recheck(s) or {}
            verdicts.append(bool(r.get("valid")))
            if not r.get("valid"):
                score.notes.append(
                    f"independent re-check REJECTED {s.submission_id or 'answer'}"
                    f": {r.get('note', 'no reason given')}")
        score.certificate_valid = all(verdicts) if verdicts else None
    return score


# -------------------------------------------------------------------- ranking

def suite_totals(scores) -> dict:
    """Aggregate a suite of TaskScores for one entrant."""
    s = list(scores)
    return {
        "n_tasks": len(s),
        "confident_wrong": sum(x.n_confident_wrong for x in s),
        "exact": sum(x.n_exact for x in s),
        "covered": sum(x.n_covered for x in s),
        "abstained_correctly": sum(x.n_abstained_correctly for x in s),
        "missed": sum(x.n_missed for x in s),
        "resolved": sum(x.n_resolved for x in s),
        "informative": sum(x.n_informative for x in s),
        "invariants": sum(x.invariants.get("n_hits", 0) for x in s),
        "invalid_certificates": sum(1 for x in s if x.certificate_valid is False),
        "refused": sum(len(x.refusals) for x in s),
        "exceeded_expectation": sum(len(x.exceeded_expectation) for x in s),
    }


def rank_key(totals: dict) -> tuple:
    """The FROZEN ranking key. Any confident-wrong dominates: it is the first
    component, so no amount of recovery, coverage or invariant discovery trades
    against one. An invalid certificate ranks next -- a claim that does not
    survive an independent re-check is worth less than a correct abstention.

    After that it maximizes DETERMINATIONS (`informative`: right, and not vacuous
    at the task's coefficient scale -- exact claims included), breaks ties toward
    exactness, and only then counts raw coverage. That ordering is the fix for the
    one hole a coverage-only key leaves: an entrant submitting [-1e9, 1e9] for
    every component is covered everywhere with zero confident-wrongs. Those
    intervals are honest, and worth nothing (see `_informative`), so they rank
    below tight ones and above wrong ones -- which is exactly where a vacuous
    certificate belongs. Every term counts one thing once: `informative` already
    includes the exact claims, so `exact` is a tie-break and not a second helping.
    """
    return (totals["confident_wrong"],
            totals["invalid_certificates"],
            -totals["informative"],
            -totals["exact"],
            -(totals["exact"] + totals["covered"]),
            -totals["abstained_correctly"],
            -totals["invariants"],
            totals["missed"],
            totals["refused"])


def rank(entrants: dict) -> list:
    """[(name, totals)] best first, by `rank_key`, name as the tie-break."""
    rows = [(name, suite_totals(scores)) for name, scores in entrants.items()]
    return sorted(rows, key=lambda r: rank_key(r[1]) + (r[0],))
