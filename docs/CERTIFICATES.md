# Public certificate release 1

Every certificate the real-data program has produced, in one place: the law,
its significance bound α (chance-fit probability under the null, log₁₀), the
domain it is certified over, the frozen data artifact, and the command that
reproduces it. Certificates are exhaustive per-point checks under the
declared error model — every point within ε, or no certificate. Everything
not listed here that the case studies examined ended as a labeled conjecture
or an explicit abstention (the doctrine: certified / conjectured / abstained,
never guessed). Engine as of the loose-ε parsimony and CAP-T reach closures
(2026-07-27; exact commits in git history).

Conventions: `x_0, x_1` are the stated inputs in order; α ≤ 10^(value shown);
rationals are exact (the coefficient gate demands the data pin them).

## Gaia DR3 (live archive, frozen ADQL snapshots in `experiments/gaia/data/`)

Reproduce: `.venv/bin/python experiments/gaia/run_c0.py`
(study: `docs/CASE_STUDY_GAIA_C0.md`; snapshot `c0_cc027870a063.csv`,
400 sources, floor 5×10⁻⁶ mag from measured column precision)

| Law | Inputs | α ≤ | Note |
|---|---|---|---|
| `mag_G = 25.687367… − (5/2)·x_0` (zp exact rational 12660019/492850) | x_0 = log₁₀ flux_G | 10⁻⁴⁶¹ | slope is the EXACT rational 5/2; zp matches the DR3 G Vega zeropoint 25.6874 |
| `mag_BP = 25.338542… − (5/2)·x_0` | x_0 = log₁₀ flux_BP | 10⁻⁴⁶⁶ | |
| `mag_RP = 24.747895… − (5/2)·x_0` | x_0 = log₁₀ flux_RP | 10⁻⁴⁵⁸ | |
| `bp_rp − mag_BP + mag_RP = 0` | invariant over the three columns | 10⁻⁴⁵⁸⁶ | exact unit coefficients; α under the term-scale null |

What did NOT certify, correctly: the color–magnitude relation (conjecture);
Stefan–Boltzmann over FLAME point-estimate columns (structural abstain — the
published marginal medians do not compose; `docs/CASE_STUDY_GAIA_P1.md`);
every kinematic relation (Strömberg stays a labeled conjecture, k ≈ 109 km/s).

## Solar system (NSSDC fact-sheet tables, committed CSVs)

Reproduce: `.venv/bin/python experiments/solar/run_case_study.py`
(study: `docs/CASE_STUDY_SOLAR.md`; declared σ: 2×10⁻³ planets, 1×10⁻³ moons,
3×10⁻³ Saturn secondary with J₂ envelope)

| Law | Inputs | α ≤ | Note |
|---|---|---|---|
| `P_yr = 0.99994·x_0^(3/2)` | x_0 = a (AU), 8 planets | 10⁻¹⁴ | Kepler III; tiny-data mode, α bounds the selection exposure |
| `P_d = 6.4577×10⁻⁹·x_0^(3/2)` | x_0 = a (km), Galilean moons | 10⁻¹⁷ | |
| `P_d = 1.1812×10⁻⁸·x_0^(3/2)` | x_0 = a (km), Saturn moons, Hyperion excluded, J₂ envelope declared | 10⁻¹⁵ | the primary run honestly abstained; residuals decode into Saturn's J₂ + the 4:3 Titan resonance |

## Macroeconomy (FRED quarterly US accounts, committed snapshot)

Reproduce: `.venv/bin/python experiments/macro/run_case_study.py`
(study: `docs/CASE_STUDY_MACRO.md`; σ declared from vintage-revision scale)

| Law | Inputs | α ≤ | Note |
|---|---|---|---|
| `C + I + G + NX − Y = 0` | national-accounts columns, invariant | 10⁻⁸⁵⁷ | the GDP accounting identity, recovered blind from raw columns; five scaled variants certify alongside (α from 10⁻³⁸⁸); bound re-baselined under term-scale α (logged in the case study) |

## Gaia DR3 Phase 2 (`docs/CASE_STUDY_GAIA_P2.md`)

Reproduce: `.venv/bin/python experiments/gaia/run_p2.py`
(snapshot `c4rr_bf9b07ed638a.csv`, 800 RRab stars; floor 4.2×10⁻⁷ dex by the
v2 propagated-quantization procedure)

| Law | Inputs | α ≤ | Note |
|---|---|---|---|
| `[Fe/H] ≈ 3.03368 − 20.12750·x_0 + 1.36843·x_1 + 6.2700001·x_0·x_1 − 0.7200000·x_1²` (exact rationals in the certificate) | x_0 = pf (d), x_1 = φ31_G | 10⁻¹⁰⁰⁹ | the DR3 RRab photometric-metallicity pipeline formula, recovered OUTRIGHT after the CAP-T reach closure; the cross/quadratic coefficients match Nemec et al. (2013) exactly (6.27, −0.72); linear-term provenance left open (a φ-offset composition matches only to ~10⁻⁶). First certified via the verify track (pinned), then re-earned by `recover` — both surfaces agree to 8×10⁻⁸ on the snapshot |

What did NOT certify, correctly: any cross-catalog Kepler law over the
astrometric binaries (mass-dependent family — structural abstain); the
Leavitt law (labeled conjecture, slope −2.21 mag/dex).

## Gaia DR3 Phase 3 (`docs/CASE_STUDY_GAIA_P3.md`)

Reproduce: `.venv/bin/python experiments/gaia/run_p3.py`
(snapshot `c6_9fe3e2b6bd35.csv`, 400 sources; open-discovery sweep, 134
cells, certificates only where the pipeline defines the column)

| Law | Inputs | α ≤ | Note |
|---|---|---|---|
| `mag_G = 25.687367… − (978839/901544)·ln(flux_G)` | x_0 = raw flux_G | 10⁻⁴⁶⁰ | the C0 anchor re-found blind by the open sweep; 978839/901544 = 5/(2 ln 10) to 7 digits |
| `parallax_over_error = C·parallax/parallax_error`, C = 1 to ~10⁻⁸ | ratio triple | 10⁻⁶²¹ | σ_rep = 3×10⁻⁸ declared (two float32 roundings compose) |
| `flux_over_error_G = C·flux_G/flux_error_G`, C = 1 to 5×10⁻¹⁰ | ratio triple | 10⁻⁵⁵⁰ | same declaration |

What did NOT certify, correctly: the IAU galactic-frame rotation row — a
structural abstain that decodes into the unit-sphere constraint on the
direction-cosine inputs (every rival class equals the IAU row modulo the
constraint; registered open issue: constrained-input coherence); all
kinematics (Oort A = 15.8 in band, B = −14.9 an honest band miss, both
conjectures); 131 of 134 open-sweep cells (heavy abstention as designed).
The phase also converted a caught vacuous certificate (a 35-term
interpolation of 12 points with α ≤ 1) into an engine-wide significance
gate: certification now requires α ≤ 10⁻⁶.
