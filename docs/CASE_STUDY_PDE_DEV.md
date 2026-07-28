# PDE dev campaign C0 — registration

**Registered 2026-07-28, before any implementation or run.** Design and
rationale: `DIRECTION_PDE.md`. This document fixes what C0 is allowed to
claim, what it must refuse, and what result would falsify the approach — in
advance, as every campaign in this program does.

## What C0 certifies over

**The weak form, and only the weak form.** A PDE claim is a relation among
derivative fields; derivatives estimated from sampled data by finite
differences amplify noise like σ/hᵏ plus a truncation error that is neither
measured nor declared, and this program does not certify against numbers whose
error was not declared (the σ_rep/float32 lesson). So every candidate term is
integrated against a smooth compactly-supported test function φ and every
derivative is moved onto φ by parts:

    ∫∫ φ · ∂^α(g(u)) dx dt  =  (−1)^{|α|} ∫∫ (∂^α φ) · g(u) dx dt

The right-hand side touches only the RAW field u through a pointwise g, and
weights it by an ANALYTIC known function. A "data point" is one patch integral;
the certification domain is the patch family, not the grid.

**Library restriction (stated, not hidden):** the term library is
{∂^α(g(u))} — divergence-form terms only. `u·u_xx` and friends cannot be moved
onto φ and are OUT of C0's reach by construction. This is the
conservation-law-compatible library; it covers heat, advection, Burgers, KdV,
and is exactly the restriction that makes conservation-law detection natural
later.

## The declared error model for a patch

Per patch p, the residual band ε_p is assembled from parts that are each either
computable or declared — no term is assumed:

1. **Float summation** R_p = ε_mach · Σ|contributions| — rigorous and computed
   per patch (the quadrature sum is a cancelling sum; this bounds its roundoff).
2. **Quadrature truncation** Q_p — MEASURED on a three-level resolution ladder
   (h, 2h, 4h) per patch, with the Richardson estimate at h as the declared
   value and the observed convergence order recorded. A patch whose ladder does
   not converge is **rejected as unresolved** — it is not certified against and
   not silently kept.
3. **Machine floor** — the existing `MACHINE_REL·|y|` term.
4. **Noise** — σ‖w‖₂ through the known linear functional, with the coverage
   factor the ε model already applies to declared σ.

C0 runs the σ = 0 case only, where parts 1–3 are the whole model. Part 4 is
deliberately deferred: the design matrix columns are themselves noisy
functionals of u, which is an errors-in-variables problem the current ε model
(exact X, banded y) does not express. **C1 is the noise ladder plus the EIV
extension**; claiming a noisy weak-form certificate under the current model
would be exactly the kind of undeclared-error claim this arc exists to avoid.

## Registered predictions

Data: analytic solutions only — no solver error, no downloads (heat via a
truncated Fourier series, advection via translation, Burgers via the exact
tanh traveling wave and Cole–Hopf, KdV via the exact soliton). Patches: tensor
bump test functions ψ(s) = (1−s²)^p, multiple scales and locations.

- **N1 — nulls certify nothing.** Weak-form systems built from (a) i.i.d.
  random fields and (b) smooth random fields not solving any library PDE
  certify **0 laws** out of 50 trials each.
- **N2 — heat certifies multi-IC.** Patches pooled from ≥ 3 distinct initial
  conditions of u_t = ν u_xx certify that law with ν pinned to the exact
  rational, α well below 10⁻⁶.
- **N3 — advection certifies multi-IC.** Same for u_t = −c u_x.
- **N4 — Burgers certifies multi-IC.** Pooled patches from ≥ 3 ICs of
  u_t = ν u_xx − u u_x certify that law.
- **N5 — the single-solution refusal.** On a SINGLE traveling-wave solution
  (Burgers tanh wave, or the KdV soliton) the true PDE must NOT be certified as
  a general law: on that solution u_t ≡ −c u_x holds identically, so the u_t
  column is exactly collinear with the u_x column and rival PDEs differing by
  multiples of that on-shell zero are indistinguishable. Two verdicts are
  acceptable — (a) a structural abstain, or (b) a certificate explicitly
  domain-restricted to the on-shell constraint. **I predict (a)**: the closed
  constrained-input machinery detects polynomial constraints among the INPUT
  columns, and this constraint binds the TARGET column to an input, which that
  detector does not see. Certifying an unrestricted PDE here is a
  form-overclaim and would falsify the approach; certifying a DIFFERENT PDE is
  a confident-wrong.
- **N6 — patch-family dependence is visible, not silent.** Changing the patch
  family (scales, count, locations) must not change WHICH law certifies. If it
  does, the certificate was a patch artifact and the arc stops until the
  dependence is understood.
- **N7 — zero confident-wrong throughout**, as everywhere in this program: no
  certified law that is not the truth of the field it was built from.

## What would end the arc

A certified weak-form law that is not the truth (N7), or a law that certifies
only for one patch family (N6). Both are checkable inside C0 with no external
data, which is why C0 is entirely analytic and runs before PDEBench is touched.

## Results

Run 2026-07-28 (`experiments/pde/run_c0.py`, results
`experiments/results/pde_c0.json`; factory `lagh/weakform.py`, tests
`tests/test_weakform.py`). Patch family A: 49×25-point patches, bump exponent
p = 16, 24 patches per solution, 4 solutions per PDE.

| prediction | verdict | measured |
|---|---|---|
| N1a i.i.d. null | **met** | every patch rejected at the aliasing gate — the null never forms a weak system |
| N1b smooth null | **met** | 96 patches, "no law certifies through tier 3" — structural abstain |
| N2 heat multi-IC | **met** | certified `u_t = [u_xx]/10` — ν = 1/10 exact, α ≤ 10⁻¹⁶³·⁹ |
| N3 advection multi-IC | **met** | certified `u_t = −7[u_x]/10` — c = 7/10 exact, α ≤ 10⁻¹⁷⁴·⁵ |
| N4 Burgers multi-IC | **met** | certified `u_t = [u_xx]/5 − [u·u_x]` — ν = 1/5 exact, α ≤ 10⁻¹⁶⁹·⁹ |
| N5 single-solution refusal | **FALSIFIED as registered, then closed** | see below |
| N6 patch-family independence | **met** | identical laws across families A/B/C (α 10⁻¹⁰⁰·⁴…10⁻¹⁴⁵·³); the band moves by 3 orders, the law does not |
| N7 zero confident-wrong | **met**, with the scope point below stated precisely | |

Every certification above is on patches from a **held-out solution**: fit and
select on three initial conditions, certify on a fourth the fit never saw.

**Deviation from the registration, stated:** N1 was registered as 50 trials of
each null; C0 ran one i.i.d. field (24 patches, every one rejected at the
resolution gate) and one 4-solution smooth null (96 patches, structural
abstain). The verdicts are unambiguous but the trial count is smaller than
registered, so N1 is scored "met" at 1/50th the registered evidence — the full
50-trial null sweep is the first item of C1, not a thing to quietly call done.

### N5: the single-solution result, and what it forced

Registered prediction: a single traveling wave abstains. What a naive row split
actually does, measured:

- **KdV single soliton → CERTIFIED `u_t = −[u_x]`, α ≤ 10⁻³⁷·⁸.** That is the
  on-shell traveling-wave relation (the soliton moves at speed 4k² = 1), which
  is TRUE of that field and is *not* the KdV equation. Diagnosis: the cheap
  single-term pre-pass finds the simplest law that fits, it certifies alone, and
  the true 2-term KdV combination never enters the contest at all — the same
  shape as the Müntz finding in `MUNTZ_ARBITRATION.md` (what goes wrong is not
  a bad choice between rivals, it is the truth being absent from the contest).
- **Burgers traveling wave → structural abstain** ("2 materially different
  classes certify at tier 1"), the registered mechanism, on the same data shape.
  So the registered prediction was right about one field and wrong about the
  other; which one you get is a property of the library and the pre-pass, not of
  the physics. That is not a basis for a claim.

The fix is the piece `DIRECTION_PDE.md` already registered (multi-solution
coherence), sharpened by this measurement into a gate: **certification patches
must come from a solution the fit never saw.** A PDE is a claim about a family
of solutions; one solution supports only an on-shell statement about that
solution. With the gate, single-solution data returns `abstain: single-solution`
structurally — not by luck of the library — and multi-solution certification
becomes strictly stronger evidence than the row split it replaces. Both verdicts
are recorded side by side in the results file (`*_rowsplit` keys) because the
contrast IS the finding.

On N7: the row-split KdV certificate was a true statement about its data, not a
false one — a scope overclaim, not a confident-wrong. The holdout gate removes
the scope, which is the only honest way to hold the distinction.

### Measured facts about the factory (each one cost a wrong first attempt)

- **The bump exponent sets the band.** p controls how many derivatives of φ
  vanish at the support edge. On the heat field the `u_xxx` bound falls
  1.5×10⁻⁶ → 7.5×10⁻¹¹ from p = 8 to p = 16 at fixed patch size.
- **The band is set by the roughest term in the library.** At p = 8, adding
  `u_xxx` widens every patch's declared band ~1000× — so the declared band is a
  property of the LIBRARY, not of the data, and must be reported with it.
- **The refinement ladder alone can be fooled.** A field at ~4 points per
  wavelength aliases to a different low frequency at each coarse level, so the
  three integrals agree smoothly and the ladder reports convergence order 17 on
  garbage. The direct test is spectral: energy above half-Nyquist in the
  windowed patch field. Measured separation — resolved fields 7×10⁻⁷…1.2×10⁻⁶
  (the window's own tail, the floor), aliased field and white noise 0.3…0.8.
- **Richardson with the observed coarse-pair order UNDER-declares the band.**
  The (4h, 2h) pair reported order ~10 where the true h-level ratio was ~2⁷; the
  resulting band was ~7× too small and rejected the true advection law. The
  order credited to the bound is capped at 4; the observed order is recorded.
- **The band's own assumption is checked, not assumed.** ε is assembled before
  the coefficients are known, so it bounds them by a declared `coeff_max` = 2;
  `coeff_audit` then verifies the certified law against it. Measured: every
  certified law is linear in the columns with max |c| ≤ 1.0 — the assumption
  held, and a violation would demote the certificate rather than pass silently.
- Deterministic bounds enter ε with **coefficient 1** (`certify.epsilon(hard=)`),
  not through the κ/λ coverage factors — those exist to turn a stochastic scale
  into a band, and multiplying a computed bound by 4 loosens ε in the direction
  that admits impostors.

### Open, and deliberately not claimed

σ > 0 (needs the errors-in-variables ε model — the design matrix columns are
themselves noisy functionals of u), variable coefficients, systems and 2D,
non-divergence-form terms, and every claim against external data. PDEBench is
untouched: C0 is analytic on purpose, so a failure here is the instrument's and
never the integrator's.
