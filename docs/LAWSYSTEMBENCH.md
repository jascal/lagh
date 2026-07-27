# LawSystemBench — certified discovery of equation SYSTEMS (H3)

**Registered 2026-07-27, spec before generator, predictions before baseline.**
The benchmark class that doesn't exist (`ROADMAP.md` H3): every current SR
benchmark assumes one equation per dataset, pre-labeled roles, explicit form.
Real theories are coupled systems with shared constants and conserved
quantities. LawSystemBench scores exactly that.

## Design principles

1. **Roles unlabeled.** A problem is a set of NAMED OBSERVABLE COLUMNS (state
   variables and their rates). Nothing marks which columns are "targets" —
   discovering the dependency structure is part of the task.
2. **Systems, not equations.** Ground truth is a SET of simultaneous
   equations, with CONSTANTS SHARED across equations where the physics shares
   them, plus CONSERVED QUANTITIES (invariants) that are functions of the
   observables, constant along every trajectory.
3. **Exact ground truth, honest scoring.** Every instance is generated from a
   seeded simulator committed alongside; scoring is mechanical: per-equation
   structural accuracy (the judge-v4 comparator), invariant recovery,
   shared-constant consistency, role identification.
4. **Two difficulty axes**: number of coupled equations (2–4) and noise
   (clean float64 vs declared σ_rep).

## v1 families (self-generated, seeded; econ-sae identities registered as the
v1.1 external family)

| family | system | shared constants | invariants |
|---|---|---|---|
| F1 chain | A→B→C decay: dA=−k₁A; dB=k₁A−k₂B; dC=k₂B | k₁ (2 eqs), k₂ (2 eqs) | A+B+C |
| F2 SIR | dS=−βSI; dI=βSI−γI; dR=γI | β (2), γ (2) | S+I+R |
| F3 Lotka–Volterra | dx=ax−bxy; dy=−cy+dxy | b,d cross-appear | c·ln x−d·x+a·ln y−b·y |
| F4 coupled oscillators | dv₁=−(k₁+k_c)/m₁·x₁+k_c/m₁·x₂; dv₂=… | k_c (2 eqs) | total energy |
| F5 planar central force | dvx=−GM·x/r³; dvy=−GM·y/r³ | GM (2 eqs) | L=x·vy−y·vx and E |

Per family: 8 seeded parameter draws × clean + one σ_rep tier → **80 problems
v1**. Observables include the rate columns (SINDy convention — recovery is
algebraic, not ODE-solving). **Amendment (2026-07-27, after 15 baseline
problems, before any results were published):** rate columns are LABELED as
rates (`r0..rN`, states `c0..cN`) with the rate↔state assignment shuffled and
hidden — the original fully-unlabeled design let the discoverer certify
true-but-degenerate alternates (`dB = −dA − dC` is exactly true in a closed
system), which structure scoring cannot distinguish from misses. Equations are
sought as rate = f(states); WHICH rate belongs to WHICH state remains part of
the discovery task. The 15 pre-amendment rows were discarded.

## The system discoverer (`lagh/systems.py`)

- **Role + equation discovery**: for every column, attempt certified passive
  discovery of it as a function of the others; the certified set defines the
  dependency structure.
- **Invariant discovery**: nullspace method — build the registered term
  library over ALL columns, center the design matrix, take smallest singular
  directions, sparsify, then CERTIFY constancy (every sample within ε of the
  constant) with minimality and α. This is inherently system-level: no
  per-equation run can produce it.
- **Cross-equation coherence**: constants appearing in multiple certified
  equations must agree (exact rationals: equality; floats: within ε); the
  system certificate conjoins equations + invariants + consistency, with a
  combined α (union bound).

## Registered predictions (before the baseline run)

- **P1**: ≥ 80% of clean problems get ALL equations structurally correct
  (these are C1/C2-class forms with the CAP-P/CAP-B features present).
- **P2**: linear invariants (F1, F2) recovered on ≥ 90% of clean problems;
  the nonlinear LV invariant on ≥ 50% (log terms exist in the library);
  energy/L on ≥ 50%.
- **P3**: zero structurally-wrong certified equations (the invariant carries).
- **P4 (the kill-criterion measurement)**: system-level machinery must add
  value beyond independent per-equation runs — quantified as (a) invariants
  recovered (impossible per-equation), (b) shared-constant consistency
  verdicts issued. If both are empty across the sweep, H3b reduces to
  bookkeeping and the benchmark stands alone.

## Results (v1.2 baseline, 2026-07-27/28 — 80 problems, judge-v4 comparator,
level-set invariant scoring)

| family (clean) | equations | wrong | invariants |
|---|---|---|---|
| F1 chain | 16/24 | 0 | **8/8** |
| F2 SIR | 17/24 | 0 | **8/8** |
| F3 Lotka–Volterra | 11/16 | 0 | **8/8 (the 4-term log-invariant, every draw)** |
| F4 coupled oscillators | **32/32** | 0 | 0/8 (energy: 5-term search gap) |
| F5 Kepler | 16/32 | 0 | 8/16 (L every draw; E needs 1/√(x²+y²), not in library) |
| **clean total** | **92/128 (72%)** | **0** | **32/48** |

Noisy tier (σ_rep=1e-6): certified equations 0 BY DESIGN — the significance
boundary makes dense linear channels empirical-only under declared noise; the
registered refinement is α-justified certification for SMALL supports (a
2-term support over ~30 terms has |H|≈450 — the boundary was drawn for dense
classes). Noisy invariants 0/48: constancy-threshold scaling under noise is a
registered v1.3 diagnostic.

**Prediction scoring:**
- **P1 FAILED** (all-equations-correct 12/40 vs ≥80%): residuals enumerated —
  F1/F2 middle equations (2-term, remaining zero-crossing/snap interactions),
  Kepler force laws (composite-base `(x²+y²)^{3/2}` — honest grammar boundary).
- **P2 PARTIAL, exceeded where it counts**: linear invariants 16/16 (bar 90%),
  LV log-invariant 8/8 (bar 50%), angular momentum 8/8; energies 0 (F4's
  5-term greedy path and Kepler's missing library term — both named).
- **P3 MET: zero structurally-wrong certified equations across all 80** —
  restored by the multi-trajectory redesign (the single-trajectory wrongs were
  identifiability artifacts, not instrument fabrications).
- **P4 MET — the kill criterion resolves in favor of H3b**: 32 invariant
  recoveries (unreachable by any per-equation run) + shared-constant
  consistency verdicts issued. System-level machinery earns its existence.

**The central finding (the benchmark's contribution):** single-trajectory
system data lies on a curve, under-determining BOTH the vector field and the
invariant level sets; lagh's coherence machinery *detected this as structural
abstention* before it was understood, and the multi-trajectory redesign
(varied ICs + run-id + level-set scoring) is what makes system discovery
well-posed. Registered v1.3: small-support α certification under noise, F4/F5
invariant-library extensions, the 2-term equation residual, econ-sae
accounting-identity family.
