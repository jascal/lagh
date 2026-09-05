# Refusals as measurements — first audit

2026-09-05. Tag: **empirical** for the source audit and executed fixtures;
**open** for proposed attribution experiments. No new proved claims.

This audits producers in `lagh/`, including the public tool and wrappers. A
reason is not a sufficient statistic: `structural` can mean failed residuals,
multiple fitting classes, or an exhausted proposal search. Consumers must keep
producer and gate identity. A refusal alone is not evidence of randomness or
of a physical cause. Some gates compute counts rather than physical quantities.

| `certify.Abstain` | What the live gate computes / retains | Extra observable and attribution limit |
|---|---|---|
| domain | Reserved: no enum-based producer found | A query location plus the established cover would measure membership, not cause. No measurement to ship from this enum today. |
| structural | `engine.discover`: certifying class count, probe disagreements, or no candidate; partial invariant content already survives some branches. `mcp.core.verify`: per-row residual against epsilon. | Rival predictions and independent observations where they separate distinguish declared rivals; absent proposals remain open. A residual needs independently calibrated instrument response or a controlled intervention before naming a physical cause. |
| noise | `engine`: zero-law feasibility or log10 significance versus -6; `verify` likewise; state certificates also enforce vacuity/significance. Alpha already survives significance demotion. | Replicate observations with independent error calibration can separate measurement scatter from reproducible signal. An alpha deficit measures evidence accounting, not noise amplitude. |
| surrogate | Reserved: no enum-based producer found | Would require direct observations independent of the surrogate to attribute surrogate error. No current emitted measurement. |
| numerical | `verify`: nonfinite form evaluation or refit scale | Valid-domain observations and independently established units can distinguish domain misuse from a numerical representation failure. No finite residual is available on undefined rows. |
| range | Split/finite-point counts (`engine`, `passive`, `verify`); output range versus floor (`acquisition`); integer or Ito row availability (`quasipoly`, `ito`) | New valid rows or a wider sampled range test insufficient evidence; replicate precision estimates test a floor explanation. Counts cannot identify a physical parameter. |
| parametric | `engine`, `verify`: feasibility after coefficient perturbation; parameter intervals on applicable engine paths | Independent parameter calibration or observations at inputs separating the perturbed predictions can break degeneracy. A perturbation pass/fail alone is not an uncertainty interval. |
| heldout | `acquisition`: fresh-box exhaustive check reduced to bool, then demotion; `mcp` forwards failed holdout | Repeated independent acquisition with recorded residuals and calibration can distinguish sampling variation from a persistent mismatch. Current boolean cannot quantify a physical discrepancy. Calling every failure a box-selection artifact overattributes its cause. |
| resolution | `ito`: surviving row counts, rejected-row records from quadrature/observation-noise gates; `statecert`: row count minus basis degrees of freedom below the evidence floor | A second sampling rate with independently measured instrument response separates discretization loss from response loss. Row survival alone does not estimate diffusion or drift. |
| coverage | Reserved in the enum/checker vocabulary; no direct enum-based production path found | Independent quadratic-variation bound and declared tail model would support a coverage statement. The coverage factor is already computable; missing premises are not a physical measurement. |
| exploration | Reserved in the enum/checker vocabulary; no direct enum-based production path found | Independent trajectories/initial conditions can test exploration; a missing exploration declaration has no measured mixing time. |

## Vocabulary outside the enum

`pdesystem` emits `single-solution`; `ito` emits `single-trajectory` (counts
and held-out-family availability, requiring a distinct solution/trajectory).
`statecert` emits `amplitude-bound-violated` and
`no-state-explains-the-observations` (declared amplitude feasibility and LP
feasibility; independent amplitude/initial-state observations or an instrument
calibration are needed before attributing failure to the dynamics).
The MCP layer also emits `bad-request` and `malformed-form`; these are input
validation outcomes, not measurements of the system. `stochcheck` accepts the
enumerated refusal vocabulary as submissions and scores it, rather than
measuring all the corresponding quantities itself. Instrument axis/record gates
have their own statuses and measured ratios; C3 is already their attribution
experiment, not an additional `Abstain` emitter.

No vocabulary migration is bundled here: changing reason strings would affect
consumers and does not recover a measurement by itself.

## Shipped measurement and registered result

`lagh.refusal.residual_measurement`, used at both exhaustive residual failures
in public `verify`, retains y minus the fitted candidate and epsilon on each
checked row. Original supplied-row indices survive finite-row filtering. The
partial determination record reports the envelope over all per-row discrepancy
intervals, explicitly conditional on the declared error model. This is neither
a fitted common offset nor a confidence interval for a physical parameter.
The candidate remains uncertified. Cause is `unresolved`; no new candidate is
fit and no gate threshold changes. Nonfinite measurement arithmetic omits the
diagnostic while preserving the refusal.

Artifact: `experiments/results/refusal_measurements.json`, reproduced by
`.venv/bin/python -m experiments.run_refusal_measurements`.

Registered P1 and P2 met: omitted-quadratic case refuses on 20 certification
rows; the training-only mismatch refuses on all 100 finite supplied rows.
P3 control retains `2*x_0`; validation results are recorded in the registration.
P4 attribution stays unresolved. Observational equivalence of additive
instrument error and an additive omitted contribution is an argument about the
specified observations, not a machine-checked theorem or an empirical discovery
of either cause. No reach gain or new physical attribution is claimed.

## Remaining work

Other live gates still discard statistics; this first patch does not finish
instrument-wide integration. Highest-value next registered steps are retaining
fresh-box holdout residuals and coefficient perturbation witnesses, followed by
independent-observable attribution experiments. Optimal held-out weak tests,
per-case approximant arbitration and LawSystemBench v1.3 remain **open**.
The blanket dense-channel restriction remains in force.
