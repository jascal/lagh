"""Solar case study (docs/CASE_STUDY_SOLAR.md — predictions frozen).

P1: Kepler III certified three times (planets, Jovian moons, Saturnian moons)
    from real fact-sheet tables at declared sigma.
P2: per-system constants agree within / differ across systems.
P3: sunspot cycle -- the astronomer's period estimator as labeled conjecture.
P4: Waldmeier effect as labeled conjecture.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.mcp.core import recover  # noqa: E402

D = Path(__file__).parent / "data"
OUT = Path("experiments/results/solar_case_study.json")


def load_system(fname):
    rows = [l.split(",") for l in (D / fname).read_text().splitlines()
            if l and not l.startswith("#")]
    a = np.array([float(r[1]) for r in rows])
    P = np.array([float(r[2]) for r in rows])
    return a, P


def kepler(fname, sigma):
    a, P = load_system(fname)
    r = recover(a.reshape(-1, 1).tolist(), P.tolist(), sigma=sigma)
    out = {"certified": r.get("certified", False), "law": r.get("law"),
           "alpha_log10": r.get("alpha_log10"),
           "abstain": r.get("abstain"), "n": len(a)}
    if r.get("certified"):
        import sympy as sp
        e = sp.sympify(r["law"])
        x = sp.Symbol("x_0")
        # exponent check: P = C * a^p -> p from the expression
        p_exp = sp.degree(sp.log(e).expand(force=True), sp.log(x)) \
            if False else None
        out["exponent_is_3_2"] = "x_0**1.5" in r["law"].replace("(3/2)", "1.5") \
            or "3/2" in r["law"]
        # constant: C = P/a^1.5 empirical
        out["C_fit"] = float(np.mean(P / a ** 1.5))
    return out


def sunspots():
    raw = (D / "SN_m_tot_V2.0.csv").read_text().splitlines()
    t, sn = [], []
    for line in raw:
        parts = [p.strip() for p in line.split(";")]
        if len(parts) >= 4 and float(parts[3]) >= 0:
            t.append(float(parts[2]))
            sn.append(float(parts[3]))
    t, sn = np.array(t), np.array(sn)
    # certification attempt on the raw series as function of time: must abstain
    r = recover(t.reshape(-1, 1).tolist(), sn.tolist(), sigma=0.0)
    # period estimate: phase-dispersion on the smoothed series (13-month box)
    k = 13
    sm = np.convolve(sn, np.ones(k) / k, mode="valid")
    ts = t[k // 2: k // 2 + len(sm)]
    best = (np.inf, None)
    z = (sm - sm.mean()) / (sm.std() + 1e-300)
    for P in np.linspace(6.0, 16.0, 4001):
        ph = np.mod(ts, P) / P
        o = np.argsort(ph)
        rough = float(np.sum(np.diff(z[o]) ** 2))
        if rough < best[0]:
            best = (rough, P)
    # Waldmeier: per-cycle rise time vs amplitude (minima-to-maxima on smoothed)
    mins = [i for i in range(60, len(sm) - 60)
            if sm[i] == sm[i - 60: i + 60].min()]
    dedup = []
    for i in mins:
        if not dedup or ts[i] - ts[dedup[-1]] > 6.0:
            dedup.append(i)
    rises, amps = [], []
    for i0, i1 in zip(dedup, dedup[1:]):
        seg = sm[i0:i1]
        j = int(np.argmax(seg))
        rises.append(ts[i0 + j] - ts[i0])
        amps.append(float(seg[j]))
    corr = float(np.corrcoef(rises, amps)[0, 1]) if len(rises) > 3 else None
    return {"raw_series_certified": r.get("certified", False),
            "raw_abstain": r.get("abstain"),
            "cycle_period_conjecture_yr": best[1],
            "n_cycles": len(rises),
            "waldmeier_corr_rise_vs_amp": corr}


def main():
    out = {}
    out["P1_planets"] = kepler("planets.csv", sigma=2e-3)
    out["P1_jupiter_moons"] = kepler("jupiter_moons.csv", sigma=1e-3)
    out["P1_saturn_moons"] = kepler("saturn_moons.csv", sigma=1e-3)
    # secondary: exclude the known 4:3 Titan-resonant moon (Hyperion) and
    # declare the J2-oblateness envelope -- the physics the primary refusal
    # detected (measured: deviations track J2 theory to ~1e-4 per moon)
    rows = [l.split(",") for l in (D / "saturn_moons.csv").read_text().splitlines()
            if l and not l.startswith("#") and not l.startswith("Hyperion")]
    a_s = np.array([float(r[1]) for r in rows])
    P_s = np.array([float(r[2]) for r in rows])
    r_s = recover(a_s.reshape(-1, 1).tolist(), P_s.tolist(), sigma=3e-3)
    out["P1_saturn_secondary_no_resonance_J2_envelope"] = {
        "certified": r_s.get("certified"), "law": r_s.get("law"),
        "alpha_log10": r_s.get("alpha_log10")}
    # P2: within-system constant stability (split halves) vs across-system
    import itertools
    consts = {}
    for name, f in (("planets", "planets.csv"), ("jup", "jupiter_moons.csv"),
                    ("sat", "saturn_moons.csv")):
        a, P = load_system(f)
        consts[name] = [float(np.mean(P[:len(P)//2] / a[:len(a)//2] ** 1.5)),
                        float(np.mean(P[len(P)//2:] / a[len(a)//2:] ** 1.5))]
    out["P2_constants"] = consts
    out["P2_within_agree"] = {k: abs(v[0] - v[1]) / v[0] < 5e-3
                              for k, v in consts.items()}
    out["P2_across_differ"] = all(
        abs(consts[a_][0] - consts[b_][0]) / consts[a_][0] > 0.1
        for a_, b_ in itertools.combinations(consts, 2))
    out["P3_P4_sunspots"] = sunspots()
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
