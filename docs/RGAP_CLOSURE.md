# R-gap closure: every curriculum gap bounded with a one-line cause

**Registered 2026-07-22.** R-gap (`STRATEGY.md`) requires *every known curriculum gap either
closed or explicitly bounded with a one-line cause*. Of 108 NewtonBench-dev cells, **66
recovered** (closed); the **42 non-recovered** are each bounded below, grouped into four
buckets. **Zero cells are unexplained.**

## Buckets

| bucket | count | meaning |
|---|---|---|
| **OUT-OF-CLASS** | 11 | an exact-rational instrument fundamentally cannot — correct-abstain, the honest ceiling |
| **UNITS / IRRATIONAL-IN-PRACTICE** | 9 | the benchmark's degree units inject π/180; effectively out-of-class for this oracle |
| **REGISTERED-CAPABILITY** | 14 | scoped, buildable — registered in NEWTONBENCH_GAP_PLAN.md |
| **REACH-BOUNDED** | 8 | exact closed form exists; a stated new tier/capability would be required |

Only the **11 OUT-OF-CLASS** cells are a permanent exact-instrument limit (the irrational-`^e`
wedge). The **9 units** cells are benchmark-specific (degrees). The remaining
**14+8 = 22** are the *addressable* frontier — each named, none silent.

## Per-cell bounding (all 42)

### OUT-OF-CLASS (11)

| cell | cause |
|---|---|
| `m11_heat_transfer hard v0` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m11_heat_transfer hard v1` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m11_heat_transfer hard v2` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m1_coulomb_force hard v0` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m1_coulomb_force hard v2` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m5_radioactive_decay hard v0` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m5_radioactive_decay hard v1` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m5_radioactive_decay hard v2` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m7_malus_law hard v1` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m7_malus_law hard v2` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |
| `m7_malus_law medium v2` | irrational exponent (^e) — no exact-rational form exists; correct-abstain (the wedge) |

### UNITS / IRRATIONAL-IN-PRACTICE (9)

| cell | cause |
|---|---|
| `m4_snell_law easy v0` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law easy v1` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law easy v2` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law hard v0` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law hard v1` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law hard v2` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law medium v0` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law medium v1` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |
| `m4_snell_law medium v2` | inverse-trig (acos/asin/atan) in DEGREES → π/180 factors make it irrational-in-practice; radian form would need an arcsin/atan tier |

### REGISTERED-CAPABILITY (14)

| cell | cause |
|---|---|
| `m0_gravity hard v0` | high-degree monomial product / (sum)^n / r² — CAP-C, unbuilt |
| `m0_gravity hard v2` | high-degree monomial product / (sum)^n / r² — CAP-C, unbuilt |
| `m10_be_distribution easy v1` | 1/(e^u−1) power-law inner — CAP-E; reverted as numerically unsound at scale, needs a log1p-stable transform (structure recoverable, coeff a float like G) |
| `m10_be_distribution easy v2` | 1/(e^u−1) power-law inner — CAP-E; reverted as numerically unsound at scale, needs a log1p-stable transform (structure recoverable, coeff a float like G) |
| `m10_be_distribution hard v1` | 1/(e^u−1) power-law inner — CAP-E; reverted as numerically unsound at scale, needs a log1p-stable transform (structure recoverable, coeff a float like G) |
| `m10_be_distribution hard v2` | 1/(e^u−1) power-law inner — CAP-E; reverted as numerically unsound at scale, needs a log1p-stable transform (structure recoverable, coeff a float like G) |
| `m10_be_distribution medium v1` | 1/(e^u−1) power-law inner — CAP-E; reverted as numerically unsound at scale, needs a log1p-stable transform (structure recoverable, coeff a float like G) |
| `m10_be_distribution medium v2` | 1/(e^u−1) power-law inner — CAP-E; reverted as numerically unsound at scale, needs a log1p-stable transform (structure recoverable, coeff a float like G) |
| `m1_coulomb_force hard v1` | high-degree monomial product / (sum)^n / r² — CAP-C, unbuilt |
| `m8_sound_speed hard v0` | transcendental×power (e^γ, ln γ) — CAP-F, unbuilt |
| `m8_sound_speed hard v1` | transcendental×power (e^γ, ln γ) — CAP-F, unbuilt |
| `m9_hooke_law hard v1` | multi-term fractional-power SUM — CAP-A part 2; deferred additive C1 fractional features (inflation risk) |
| `m9_hooke_law hard v2` | multi-term fractional-power SUM — CAP-A part 2; deferred additive C1 fractional features (inflation risk) |
| `m9_hooke_law medium v2` | multi-term fractional-power SUM — CAP-A part 2; deferred additive C1 fractional features (inflation risk) |

### REACH-BOUNDED (8)

| cell | cause |
|---|---|
| `m2_magnetic_force hard v0` | (I1+I2)^{3/2} — rational power of a SUM (exact form exists), not a monomial; needs a composite-fractional-power capability |
| `m5_radioactive_decay easy v2` | exp(−(λt)^p) fractional inner-transcendental — exact form, beyond current C4 |
| `m5_radioactive_decay medium v2` | exp(−(λt)^p) fractional inner-transcendental — exact form, beyond current C4 |
| `m6_underdamped_harmonic easy v0` | fractional-exponent rational (m^1.3) or degree-5 numerator — exact form, beyond current reach |
| `m6_underdamped_harmonic hard v0` | fractional-exponent rational (m^1.3) or degree-5 numerator — exact form, beyond current reach |
| `m6_underdamped_harmonic hard v1` | fractional-exponent rational (m^1.3) or degree-5 numerator — exact form, beyond current reach |
| `m6_underdamped_harmonic hard v2` | fractional-exponent rational (m^1.3) or degree-5 numerator — exact form, beyond current reach |
| `m6_underdamped_harmonic medium v2` | fractional-exponent rational (m^1.3) or degree-5 numerator — exact form, beyond current reach |

## R-gap: MET

Every non-recovered cell carries a one-line cause; nothing is silently unexplained. R-gap is satisfied.

## Readiness bar — COMPLETE

| gate | status |
|---|---|
| **R-cap** | ✅ MET (easy 11/12, medium 11/12) |
| **R-zero** | ✅ MET (0/108 confident-wrong, clean) |
| **R-noise** | ◐ SETTLED (`RNOISE_STUDY.md`) — exact-structure-or-abstain to ~1% noise; registrable for SNR ≥ ~40 dB |
| **R-gap** | ✅ MET (this doc) |

`R-cap ∧ R-zero ∧ R-noise ∧ R-gap` now hold/settled on NewtonBench-dev. The reserved **blind read is
ELIGIBLE** — but it stays a deliberate, user-authorized one-shot: choose the sealed benchmark, freeze
its published SOTA into a registration, restrict to SNR ≥ 40 dB for the structural guarantee, then read
**once**. Never automatic from this doc.
