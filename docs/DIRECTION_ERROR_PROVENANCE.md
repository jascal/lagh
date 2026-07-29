# Direction: error provenance — is this a measurement or a simulation?

**Scoped 2026-07-29** (user direction), out of the PDEBench dev pass, which
supplied both the motivating failure and most of the diagnostics.

Handed a field, this program should produce the best certified symbolic law it
can whether the field came from an instrument or from a solver — and it should
say WHICH it is looking at, because the two carry different errors and the
difference is not cosmetic. It decides which band channel is correct, and
getting that wrong fails in both directions.

## The crux: the error kind chooses the channel

| | observation error | simulation error |
|---|---|---|
| nature | stochastic, one realization | deterministic, reproducible |
| in time | flat | **grows** (PDEBench advection: 0 → 2.75e-2 by t = 2) |
| in wavenumber | white, or an instrument spectrum | **structured** (β_eff 0.70001 → 0.70035 with k) |
| vs the field's own content | independent | **scales with it** (8.9e-4 at k_max = 4 → 2.0e-2 at k_max = 34) |
| across conserved quantities | uniform | **asymmetric** — the scheme's invariants (PDEBench CFD: mass 1e-4, momentum 1e-2) |
| correct band channel | `KAPPA·σ·sqrt(a'Ga)` — **L2, cancels across the realization** | `field_err·‖w·g'‖₁` — **L1, nothing cancels** |

Using L2 for a deterministic error **under-declares** and admits impostors;
using L1 for a stochastic error **over-declares** by roughly √n and loses laws
that are really there. Every PDEBench certificate in `CASE_STUDY_PDEBENCH.md` is
±17% rather than ±1e-4 because the L1 channel was chosen — correctly, and **by
hand**. That hand is what this direction removes.

## The sharp instrument: recover the modified equation

A scheme integrating `u_t + βu_x = 0` does not solve that equation. It solves its
MODIFIED equation,

    u_t + β u_x = c₂ u_xx + c₃ u_xxx + ...

whose leading coefficients are the discretization's signature: c₂ is numerical
dissipation, c₃ numerical dispersion. So the residual of a stated law, regressed
on the next derivatives in the hierarchy, IS that signature — and it is measured
with machinery this program already has, since those derivatives are just more
weak-form columns.

**Measured 2026-07-29 on PDEBench advection (β = 0.7), and it works:** the
residual of `u_t + 0.7u_x = 0` is explained at **84% of its variance by `u_xxx`
alone**, coefficient **1.09e-7**, against 1.5% for `u_xx` and 0.9% for `u_xxxx`.
Purely dispersive, no dissipation — which independently reproduces the spectral
measurement (amplitudes preserved to 1.0000, phase speed rising with k). The two
agree quantitatively: a `c₃u_xxx` term gives phase speed `β + c₃k²`, predicting
`1.088e-7 × (2π·14)² = 8.4e-4` at mode 14 against a spectrally measured 3.4e-4.
Two unrelated diagnostics, same scheme property, same order.

So the honest description of that file is `u_t = −0.7u_x + 1.09e-7 u_xxx`: the
solver's equation, not advection.

### ...and its measured limit, which is the interesting part

Certifying the modified equation does NOT separate it from the physical law on
this data. Scanning the declaration:

| declared field error | plain `u_t = −0.7u_x` | modified (+ c₃u_xxx) | verdict |
|---|---|---|---|
| 1e-4 | 0.071 | 0.038 | both certify |
| 1e-5 | 0.71 | 0.38 | both certify |
| 3e-6 | 2.35 | 1.25 | neither |

The two residual ratios differ by a constant ~1.85×, so they pass and fail
together and no band separates them. The mechanism is structural: the weak form
integrates `u_xxx` against `φ_xxx`, which a smooth bump window damps hard, so the
term's contribution to a patch integral is small — while the 16% of the residual
it does NOT explain is comparable. Parsimony then correctly drops the term
(measured: the engine certifies `u_t = −0.70003u_x` alone at every band where
anything certifies).

What that buys anyway is large: moving from an opaque L1 declaration (0.0275) to
the modified-equation vocabulary at 5e-4 tightened β from **±17% to ±0.3%**, and
at 3e-5 to **±0.02%** with α from 1e-16 to 1e-87. The diagnosis pays even when
the diagnosed term is not itself certifiable.

## What the verdict is a claim ABOUT

The labels `observation` and `simulation` are operational, not ontological, and
the characterizer must say so in its own output. What is measured is STRUCTURE —
deterministic reproducibility, growth in time, monotone drift in wavenumber,
asymmetry across conserved quantities — and "simulation error" is shorthand for
that structure, never a claim about where the numbers came from. Nothing in this
program can distinguish a physical process from a sufficiently good simulation of
one, and it does not need to: a certificate claims that over a stated domain,
within a declared band, a relation holds. That claim is untouched by what the
substrate turns out to be.

The distinction is not idle, though, because the same diagnostic is used in
physics for exactly the substrate question. The measurement made here —
fitting `β_eff(k) = β + c₃k²` and attributing the drift to a `c₃u_xxx` term — is
formally the modified-equation fit that Lorentz-invariance-violation searches run
against gamma-ray burst arrival times (`v(E) = c(1 ± E/E_QG)`), where a
wavenumber-dependent propagation speed would be evidence of a discrete
substrate rather than of a discrete SOLVER. Same logic, same fit, different
subject; the experimental bounds there are severe (E_QG ≳ 10¹⁹ GeV), which is
what makes the question empirical rather than metaphysical.

So the discipline is the one this arc has already been bitten by twice: the
mundane explanation is checked against its own band FIRST. The c₃ = 1.09e-7
measured here is a finite-volume scheme's truncation error, fully accounted for
by the scheme, and reading anything else into it would be the interpretive form
of reading a construction bug as a finding.

## What to build

1. **`characterize_error(fields, coords, law=None, replicates=None)`** — measures
   the diagnostics in the table above and returns an `empirical` verdict, never a
   certificate: `structured-deterministic` / `unstructured-stochastic` / `both` /
   **`undetermined`**, with the evidence and the recommended channel and
   magnitude. The verdict names the STRUCTURE it measured; `simulation` and
   `observation` are the usual causes of each and are reported as the likely
   reading, not as the finding.
2. **UNDETERMINED must be first-class.** A single trajectory with no replicates
   and no reference cannot distinguish a deterministic error that happens to look
   irregular from a stochastic one. Saying so is the correct answer and is the
   same discipline as the state certificate's undetermined modes.
3. **The modified-equation regression** as a reportable diagnostic, with its
   variance-explained and the coefficient — plus the honest note that explaining
   the residual is not the same as certifying the term.
4. **Wire the verdict to the band**: the characterizer's output selects the
   channel rather than the analyst selecting it.

## Naming

**Observational error** and **pipeline error** are the pair to use. "Pipeline"
because a solver is one, and so is a data-assimilation system, a photometric
reduction, or a counting-and-binning step — the program already says "pipeline
decode" for exactly this in `DIRECTION_PDE.md`. Where the pipeline is
specifically a discretization, **scheme error** is the narrower term. The
measured verdict stays `structured-deterministic` / `unstructured-stochastic`,
because that is what is measured; the two names above are the usual causes.

## A learned operator (FNO) as an instrument for this

Registered here rather than in `PDEBENCH_READINESS.md`, because error provenance
turns out to be the better use for one than equation proposal is. Four roles, in
order of value, and none of them ever touches the band:

1. **A reference where none exists.** The measurement needs something to compare
   against: advection had a closed form, Burgers an independent high-accuracy
   solve, and CFD had NEITHER -- which is why that run could only report the
   declaration required rather than measure the error. A learned operator
   supplies a reference for families this program cannot integrate.
2. **Resolution extrapolation → the discretization error itself.** An FNO is
   parameterized in mode space and is resolution-transferable, so training across
   the resolutions a benchmark ships and extrapolating the learned operator to
   h → 0 estimates the scheme's error directly. Richardson extrapolation in
   operator space, no closed form needed. This is the role that would unblock
   the CFD case.
3. **Symbol probing → the modified equation, better conditioned.** Perturbing a
   trained operator with single Fourier modes at small amplitude and reading the
   linearized response recovers lambda-hat(k) directly, giving beta_eff(k) or
   nu_eff(k). The same quantity the weak-form residual regression measured
   (c3 = 1.09e-7), but fitted over all 10^4 samples rather than 4.
4. **The deterministic/stochastic split, i.e. sigma WITHOUT replicates.** This is
   the case the characterizer would otherwise have to call `undetermined`. A
   learned operator is a deterministic function of its input, so what it predicts
   ON HELD-OUT SAMPLES is reproducible-given-input (pipeline), and the
   irreducible residual is what no deterministic operator can recover
   (observational). Held-out is load-bearing -- enough capacity memorizes noise.

**Limits, which matter more than the promises:**

* It cannot detect error shared by the whole training distribution. If every
  sample carries the same dispersion, the operator learns it as SIGNAL. It
  separates reproducible-given-input from not, never true from false -- and the
  advection dispersion measured here would be invisible to it for exactly that
  reason.
* It carries no declarable bound of its own, so it may inform which CHANNEL to
  use and never the magnitude that goes into it.
* On shocks it rings, so its residual there is its own artifact rather than the
  data's -- which is the regime it would otherwise be most wanted for.

## Registered predictions

- **P1.** On our own C1/C2 fields with declared σ (exact analytic solutions plus
  known noise), the characterizer returns `observation`: no higher-derivative
  term explains more than 10% of the residual variance, and the residual is flat
  in t.
- **P2.** On PDEBench advection it returns `simulation`, names dispersion, and
  recovers c₃ within a factor of 2 of 1.09e-7 across independent sample sets.
- **P3.** On a single trajectory with no replicates and no stated law it returns
  `undetermined` rather than guessing.
- **P4.** Zero confident-wrong in the DANGEROUS direction: it never returns
  `observation` for a field whose error is deterministic. The reverse error
  (calling observation error `simulation`) is merely conservative — it widens the
  band — and is permitted.
- **P5.** On PDEBench CFD the modified-equation terms differ between the
  continuity and momentum equations, consistent with the measured 100× asymmetry
  in the declaration each needs.
- **P6 (learned operator).** Trained on a family with declared σ and held out by
  sample, an operator's irreducible residual estimates that σ within a factor of
  2 — recovering the noise scale with no replicates. On PDEBench advection the
  same procedure returns a residual floor far BELOW the measured 2.75e-2
  deviation, because the dispersion is common to the training distribution and is
  learned as signal. Both halves are the prediction: it finds observational error
  and is blind to distribution-wide pipeline error, and that asymmetry is exactly
  why it informs the channel and never the magnitude.

## Why this matters beyond PDEBench

Every external data source this program has touched or plans to touch falls on
one side or the other, and several are mixed: Gaia and the exoplanet archive are
instruments with pipelines, ERA5 and PDEBench are model output, traffic data is
an instrument with a discretization (counting and binning). The two-strata
treatment the traffic case study already plans — certifiable conservation
against a conjectured closure — is the same idea one level down. An instrument
that can say which kind of error it is looking at can pick its own channel, and
an instrument that cannot has to be told, once per dataset, by someone who might
be wrong.
