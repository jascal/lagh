# BLIND READ REGISTRATION — Gravity-Bench-v1 (Horizon 2a)

**Registered 2026-07-27, BEFORE the dataset/scenarios are downloaded or
inspected.** Second blind read of the program; executes under the discipline
battle-tested on LLM-SRBench (`BLIND_READ_REPORT.md`), including its lessons.

## 1. Frozen public facts (metadata only: arXiv 2501.18411 / ICML 2025,
gravitybench.github.io, fetched 2026-07-27)

- Environment benchmark: agent-as-astronomer over simulated binary-star
  systems; scenarios include standard orbits, unbound systems, **modified
  gravity (power-law exponents), drag forces**, proper motion, varied units —
  including OOD physics that deviates from the real world.
- Interface: the agent REQUESTS observations (times of its choosing) under a
  budget (~100 points); a full-observation variant exists. Tasks ask for
  quantities: orbital period, stellar masses, gravity-modification exponent,
  drag coefficient, etc.
- Scoring: per-task correctness within problem-specific accuracy thresholds
  against expert solutions.
- **Frozen baselines:** o4-mini-high **74%** (full data) / **49%** (budgeted
  100-observation); other frontier models lower. Code:
  github.com/NolanKoblischke/GravityBench; data: HF GravityBench/GravityBench.

## 2. The entry: a DETERMINISTIC lagh-powered astronomer (no LLM)

The agent is the Orca-machine pattern with a hand-authored playbook instead of
a proposer — zero LLM calls (cost directive honored; the fit is natural since
the tasks are quantitative, not open-ended):

1. **Observation planning** = lagh's active-acquisition policy adapted to time
   sampling: coarse uniform pass → period detection on the trajectory (the
   oscillation IS the signal) → refined sampling at the detected timescale
   (phase coverage), all within budget. Registered as a FIXED policy, no
   per-task tuning.
2. **Quantity extraction playbook**, hand-authored ONCE against the paper's
   task-type list (metadata), never against instances: period from certified
   periodic fits; masses via Kepler/center-of-mass relations on fitted orbital
   elements; gravity-modification exponent via lagh C3 recovery on the fitted
   force law; drag coefficient via decay-rate fits (C4/C9 family). Every
   quantity traces to a lagh fit; certificates + α attached where the fit
   certifies, best-fit values submitted regardless (the benchmark scores
   numeric accuracy — the two-track lesson: abstention earns nothing, labels
   travel anyway).
3. **Development happens ONLY on self-built synthetic orbits** (our own
   two-body integrator, written for this purpose): the GravityBench scenarios
   stay sealed until the one-shot read. Their REPO CODE (interface, task list)
   is protocol and may be read; their scenario data/parameters may not.

## 3. Declared-noise rule (the float32 lesson, now standing policy)

At read time, BEFORE discovery: inspect array dtypes and declare
σ_rep = eps(dtype) as representation noise; simulation-integrator tolerance,
if documented in their code, may also be declared. No other per-task noise
tuning.

## 4. Accounting (pre-registered)

Report: overall % correct (their thresholds) for BOTH budgeted and
full-observation variants vs the frozen baselines; per-task-type breakdown;
observation-budget usage; the certified-fit fraction (how many submitted
quantities trace to a certificate, with α) vs best-fit-only; and any task type
the playbook cannot express (reported as such, scored 0 — no improvisation at
read time).

## 5. Walk-away and honesty commitments

- If repo-code reading reveals a scoring/interface regime the deterministic
  agent cannot express (e.g. free-text physics explanations dominate the
  score), walk away BEFORE the read; the benchmark stays sealed.
- One read per variant (budgeted, full). Crash-fixes logged; scoring-affecting
  changes prohibited once scenarios are opened.
- The report states the result vs the frozen 74%/49% — win, lose, or mixed.
