# The frozen checker for the stochastic suite

**Registered 2026-07-29. Step 1 of [`DIRECTION_STOCHASTIC.md`](./DIRECTION_STOCHASTIC.md),
frozen before any Level 0 task exists.** The implementation is
`lagh/stochcheck.py`; `tests/test_stochcheck.py::test_the_interface_is_frozen`
is the lock, asserting every vocabulary tuple and every dataclass field list
verbatim. A checker that moves after the tasks arrive is a scoring key someone
can fit, so the order is deliberate: interface, then tasks, then results.

Two things had to be settled first. The certificate-kind decision was registered
in `DIRECTION_STOCHASTIC.md` (coverage-factor calibration). The open question it
left — **what κ means when the dominant term is intrinsic, and how a union bound
over n patches enters α** — is answered in `certify.coverage_factor` and
summarized in §1, because the answer is what makes the interface checkable rather
than declarative.

## 1. κ for an intrinsic term, and the union bound

Under the true law a weak-form row's residual is not an error — it is the system:

    M_p  =  ∫ φ_p f'(X) b(X) dW,        ⟨M_p⟩ = ∫ φ_p² f'(X)² b(X)² dt

a mean-zero martingale, unbounded, so no finite band holds for every row almost
surely. What does hold is the exponential martingale inequality: for any V > 0,

    P( |M_p| ≥ κ√V  **and**  ⟨M_p⟩ ≤ V )  ≤  2 exp(−κ²/2)

which is pathwise and assumes no Gaussianity, no CLT and no independence.
Exhaustiveness over n rows costs a union bound, and inverting for a declared
false-abstain budget δ gives the whole answer:

> **κ(n, δ) = √(2 ln(2n/δ))**

Four consequences, and each of them lands somewhere in the interface.

**(a) The continuity claimed for option (1) is numerical, not rhetorical.**
κ(300, 0.05) = **4.334** against the existing constant `KAPPA = 4` — the same
number to within 8%. What the constant never had is the budget: κ = 4 over 300
rows buys δ = 0.201 under the martingale bound, or 0.019 under an exact Gaussian
tail (`certify.coverage_budget`). The bound is ~10.6× looser at κ = 4, and that
gap is the price of assuming nothing about the residual's law. Past
n = e⁸/2 ≈ 1491 rows the martingale reading of κ = 4 exceeds 1 and states **no
coverage at all** — worth knowing, since weak-form runs routinely carry
thousands of patches.

**(b) κ becomes a function of the run and must be reported.** It is no longer a
module constant. A `Submission` carries `Coverage(kappa, delta, n_rows, …)` and
the checker refuses an answer whose κ does not meet its own stated budget at its
own row count. This is what makes the answer load-bearing instead of decorative:
the three numbers are checked against each other.

**(c) The band's SCALE is measured, not declared — and never by the candidate.**
The inequality needs an upper bound V on ⟨M_p⟩, and ⟨M_p⟩ depends on b, which is
one of the objects being discovered. A candidate allowed to set its own band could
widen it at will: the loose-ε failure mode, now self-serving. So V comes from
**realized quadratic variation** — a functional of the data — with only the φ²f'²
weights coming from the test function. `Coverage.qv_provenance` records this, and
`residual-derived` is refused as circular — the label
`errormodel.characterize_rows` already attaches to those numbers, hardened into a
refusal because this setting is scored. The residual self-declaration measured
there came out four orders too small; in a benchmark that is a way to manufacture
a certificate.

**(d) The union bound enters α through the DISJOINT patch count, and the trade is
favourable.** α ≤ |H|·∏q_k with q_k = 2ε_k/R assumes independent held-out rows
(`DIRECTION_SIGNIFICANCE.md`, "honest subtleties"), and overlapping patches share
the driving path. The honest held-out count is therefore the number of patches
with **disjoint support**, which is computable from the patch geometry and must be
stated: an answer reporting α without `n_disjoint` is refused. Then note the two
opposite dependences on n:

| | scaling in n | 300 → 10 000 rows |
|---|---|---|
| κ (band width) | √(2 ln n) | 4.33 → 5.08 (+17%) |
| h = n − dof (α's exponent) | linear in n | ×33 |
| log(1/q) per row | shrinks as log(1/κ) | −3% at q ≈ 0.01 |

The two right-hand columns do not cancel: log α ≈ h·log q, so widening the band
costs a few percent of each row's contribution while the row count multiplies by
33 — net, log α improves by ~32×.
**Patch count is a strictly winning resource in the stochastic regime** — the
opposite of the usual multiple-comparisons intuition, and the reason the regime is
workable at all. Registered as **S7** in `DIRECTION_STOCHASTIC.md`. The
independence discount itself is validated empirically by the null suite, the same
way `MUNTZ_ARBITRATION.md` P3 was closed at 0/200.

**What κ does not do.** It does not make the certificate pathwise. A stochastic
certificate at coverage 1 − δ is a different object from a deterministic one and
must never be compared against a deterministic α unlabelled — the discipline
`statecert.alpha_kind` already uses. Option (3) (a separate statistical class)
stays deferred; its *labelling* does not, and rides on
`Submission.alpha_kind` + `Coverage`.

The inequality is measured, not trusted:
`test_the_martingale_bound_holds_on_simulated_stochastic_integrals` simulates
state-dependent stochastic integrals (so ⟨M⟩ is random and M is not a scaled
Gaussian) and checks the observed exceedance rate against the bound at four κ.

## 2. Component claims travel only in a determination record

This is why partial determination was step 0 rather than cleanup. The registered
scoring is per component and three-regime, and `certified: bool` plus an abstain
enum cannot express it. A submission's coefficient claims live in a
`certify.determination` record and nowhere else — the same records
`engine`, `pdesystem` and `statecert` already emit, so lagh's own answers are
submissions by construction, and **interval COVERAGE of a true coefficient is a
containment test on a record that already exists**.

Component names are frozen as `part[index]:term` (`drift:x`,
`diffusion[1]:x0*x1`) and built by `stochcheck.component` so a task author and a
submitter cannot spell them differently. Terms are compared **canonically**
(through sympy): `x*x*x` and `x**3` are one column, because a false
confident-wrong in the *checker* would be worse than one in the instrument.

A task declares its whole library, **zeros included** — a term whose true
coefficient is zero is what a null is made of, and it lets nulls be scored by the
same rules as everything else with no special case. That the vocabulary needs no
null branch is the sign it is the right one.

## 3. Expectations are predictions, not a ceiling

Per component, a task registers an `expectation` from `exact | interval |
abstain`: what a calibrated instrument should manage at that sampling regime.
`validate_task` refuses a task that omits one, on the LawSystemBench rule that a
prediction registered after the run is not one.

* Beating an expectation (certifying where it said abstain) is recorded as
  `exceeded_expectation` — a finding about the expectation, never a penalty.
* An expectation **never shields**. An answer that beats the expectation and
  misses the truth is still a confident-wrong.

## 4. A declared magnitude names its consumer

The PDEBench lesson (`DIRECTION_ERROR_PROVENANCE.md`), encoded as validation
rather than a comment. Three consumers want "noise intensity" and they want three
different quantities in three different units:

| consumer | bounds | units |
|---|---|---|
| `drift-band` | ⟨M_p⟩ = ∫φ²f'²b²dt for one row | (row target)² |
| `diffusion-qv` | the diffusion rate b²(x) at the task's stated reference state | [state]²/time |
| `observation` | the measurement error on the observed state, σ_obs | [state] |

The units are the argument that no number can serve two of them, so the checker
refuses a submission that declares one consumer twice, states a `quantity` text
that does not match the frozen one for its consumer, or uses `residual-derived`
provenance. It also refuses **equal magnitudes across two consumers** — the
signature of the original mistake — but that is a signature and not a proof, since
σ_obs and b² can coincide, so the refusal is waived when both declarations carry a
`note` attesting the coincidence. The friction stays where it belongs and cannot be
tripped silently.

**The audit is the routine version of "scan the declaration."** For every consumer
the task can pin, the score sheet carries `declared / true` and flags a decade in
either direction — over-declaration as *laws will be lost*, under-declaration as
*the dangerous direction*. PDEBench over-declared 3900× and it took a session to
notice; here it is a column. Where the task cannot state a reference — the drift
band depends on a patch family only the submitter knows — the audit says the
magnitude is **not auditable** and checks the provenance instead of inventing a
ratio.

Level 2 is what makes the `observation` row meaningful: it constructs a known
process-vs-measurement decomposition, which exists nowhere else in the repo
(S6).

## 5. The five outcomes, and the one that dominates

    exact  /  covered  /  abstained-correctly  /  missed  /  CONFIDENT-WRONG

`missed` is a reach loss — the instrument left something on the table. Only
`confident-wrong` breaks the product invariant. Per component:

| claim | truth | outcome |
|---|---|---|
| `exact` v | \|v − truth\| ≤ tol·max(\|truth\|,1) | **exact** |
| `exact` v | otherwise | **confident-wrong** |
| `interval` [lo,hi] | truth inside (machine slack only) | **covered** |
| `interval` [lo,hi] | truth outside | **confident-wrong** |
| `unconstrained` | expectation was `abstain` | **abstained-correctly** |
| `unconstrained` | expectation was determinable | **missed** |
| out-of-library term, claimed `resolved` | — | **confident-wrong** |
| out-of-library term, not resolved | — | noted, no penalty |

`resolved` (the interval excludes zero, i.e. the term is established as *present*)
is scored as **reach**, not as correctness: an interval straddling zero around a
true nonzero coefficient is honest and partial, and reads `covered` with
`resolved: false`.

### What an `interval` component MEANS (clarified 2026-07-29, Level 0)

`interval` scoring reads the record as **"the true coefficient is in [lo, hi]"** —
that is what interval COVERAGE means, and it was implicit until a real producer
made it matter. Two things that look like such a bound are not one, and neither may
be submitted as an interval component:

* **A conditional interval.** `certify.parameter_interval` bisects the
  certification predicate with every *other* coefficient held. "This law still
  certifies for c_x in [lo, hi]" is true and useful and says nothing about where
  the truth is, because the true law may differ in the other coefficients too.
* **A range over a certifying set.** `certify.invariant_content` ranges over the
  laws the SEARCH found. Measured on Itô Level 0: the true drift −5x certified
  every row inside its band while that range reported the coefficient in
  [−0.66, 0]. Read as a bound, that is a confident-wrong.

`certify.admissible_interval` is the bound that holds: min/max of each coefficient
by LP over the consistent polytope {c : |y − Ac| ≤ ε}. If the truth is in the
vocabulary and the band covers, the truth is in the polytope, hence in the
interval. `lagh/ito.py` emits records from that and reports the conditional
intervals separately, labelled.

The interface shape is unchanged by this — it is a statement of what the existing
`interval` kind always meant, so the freeze holds.

### Vacuity, per component

A key that rewards coverage alone has one hole, and it is a large one: submit
`[-1e9, 1e9]` for every component and you are **covered everywhere with zero
confident-wrongs**. lagh already has the answer to this — `certify.vacuous`
refuses a certificate when the ZERO law also certifies, because then the band
swallowed the signal and the certificate proves nothing. An interval has the same
failure mode, so the same test is applied per component: an interval is
`informative` only when its span is under 2·`TAU` of the task's own coefficient
scale (`TAU` = 0.05, the material-difference scale the engine already uses); an
exact claim always is.

A vacuous interval is **not wrong** — it is honest, and it is worth nothing. So it
is never scored `confident-wrong`; it earns `covered`, is flagged in its row's
note, and is excluded from `n_informative`, which ranks **ahead** of raw coverage
(§8). A vacuous entrant therefore lands below a tight one and above a wrong one,
which is exactly where a vacuous certificate belongs.

Two rules chosen against convenience:

* **Silence is not abstention.** A component no accepted submission ever
  mentions scores `missed`, never `abstained-correctly`. The registered format
  says abstention is an explicit token with a structured reason.
* **A refusal voids credit, not exposure.** A submission refused on procedure
  (no coverage line, a circular declaration) still has its claims checked for
  confident-wrongs; its correct rows are marked `voided` and earn nothing.
  Otherwise "forget the coverage line" would launder a wrong answer into a
  non-answer.

### An abstention may carry partial determination

Added 2026-07-29 after Level 0, and the interface's field shapes did not change —
`Submission` already had both `abstain` and `record`. What changed is that the
record on an abstain submission is now **scored** rather than dropped.

"ABSTAIN[structural], and here is what every consistent law agrees on" is lagh's
primary output shape — it is the whole point of step 0 — and silently discarding
that record left the checker unable to score its own producer. The rule: the
reason is scored for abstention correctness, the record's components are scored for
coverage, and the token speaks only for the components the record does **not**
mention, so nothing is counted twice. An abstention carrying interval content is
making an interval claim, so it must state its coverage like any other.

Abstention reasons come from the one `certify.Abstain` vocabulary, which gains
three members for this layer: **`resolution`** (the sampling rate cannot separate
the quantity at this Δt — the S4 reason), **`coverage`** (the run cannot state its
coverage, so it refuses), **`exploration`** (the trajectories never visit enough
of state space). They live in the shared enum rather than as loose strings because
the alternative already happened: `pdesystem`'s `single-solution` is a bare string
no consumer can enumerate.

## 6. Domains, and several answers to one task

A domain-qualified answer is a legitimate *partial* answer — Darcy is the live
precedent, β = 0.1000 in one conductivity phase and nowhere else — so several
submissions per task are allowed and each is scored **on its own domain**, with
`coverage` reported. Narrowness is honest and is never scored as wrong.

Records sharing a domain are conjoined by `certify.conjoin_determination` before
scoring, which is where a same-domain contradiction surfaces as a finding. The
direction that matters is sound: if the truth lies in both inputs it lies in the
intersection, so conjoining can only ever expose a claim that was already wrong.
Records on **different** domains are never conjoined — that refusal is the whole
point of the qualifier.

## 7. Certificate validity is its own axis

`score_task(..., recheck=…)` takes the independent checker: the harness holds the
task's data and re-verifies the claim at the **task's** declarations rather than
the submitter's. A claim that fails re-verification is **not** promoted to a
confident-wrong (that would conflate a procedural failure with a false claim) and
is not ignored either — it is its own field, and ranks immediately after
confident-wrongs.

## 8. The frozen ranking key

```
(confident_wrong, invalid_certificates,
 −informative, −exact, −(exact + covered),
 −abstained_correctly, −invariants,
 missed, refused, name)
```

Ascending, best first. `confident_wrong` is first, so no amount of recovery,
coverage or invariant discovery trades against one. An invalid certificate ranks
below a correct abstention. Then the key maximizes **determinations**
(`informative` — right and not vacuous, exact claims included), breaks ties toward
exactness, and only then counts raw coverage: coverage alone is gameable by a wide
enough interval (§5), and a claim that excludes nothing must not outrank one that
excludes something. Each term counts one thing once — `informative` already
contains the exact claims, so `exact` is a tie-break, not a second helping. Invariants are a bonus, scored **up to affine
reparametrization** (an invariant is defined only up to scale and offset) via an
optional evaluation hook — and the hook's samples must span more than one
trajectory, or every constant matches every invariant.

## What this does NOT settle

* **STRATUM**, the fifth partial-determination dimension, remains a campaign
  convention rather than a record, and is not in the interface.
* **Level 3** (generator/Fokker–Planck-only identifiability) is deferred, so no
  outcome describes a claim about a stationary density.
* The **jump / switch** parts are in `PARTS` and reachable by the same component
  vocabulary, but Level 2 is where they are first exercised; nothing here says how
  a switching structure is named beyond the `part[index]:term` convention.
* κ is registered and its inequality is checked in simulation; **it is not yet
  calibrated on a discovery run.** That is Level 0, step 2, and it is what comes
  next.
