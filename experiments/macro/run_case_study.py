"""Raw-US-macro case study (docs/CASE_STUDY_MACRO.md — predictions frozen).

P1: invariant scan over NOMINAL {Y,C,I,G,NX} -> the expenditure identity,
    certified with alpha at declared BEA-rounding noise.
P2: the same scan on CHAINED-dollar components -> must NOT certify (chain
    residual is real).
P3/P4: Okun and Phillips through the two-track submission -> honest labeled
    conjectures, never certificates.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from lagh.submit import submission  # noqa: E402
from lagh.systems import discover_invariants  # noqa: E402

D = Path(__file__).parent / "data"
OUT = Path("experiments/results/macro_case_study.json")


def q(series, start="1990-01-01", end="2019-12-31", agg=None):
    df = pd.read_csv(D / f"{series}.csv", parse_dates=[0])
    df.columns = ["date", "v"]
    df["v"] = pd.to_numeric(df["v"], errors="coerce")
    df = df[(df.date >= start) & (df.date <= end)].dropna()
    if agg:                                   # monthly -> quarterly mean
        df = df.set_index("date").resample("QS").mean().reset_index()
    return df.set_index("date")["v"]


def main():
    out = {}
    # ---- P1: nominal expenditure identity ----
    nom = {k: q(s).to_numpy() for k, s in
           [("Y", "GDP"), ("C", "PCEC"), ("I", "GPDI"), ("G", "GCE"),
            ("NX", "NETEXP")]}
    n = min(len(v) for v in nom.values())
    nom = {k: v[:n] for k, v in nom.items()}
    sigma_rep = 0.05 / float(np.mean(np.abs(nom["Y"])))    # $0.05B half-ulp
    invs = discover_invariants(nom, sigma=max(sigma_rep, 1e-6))
    out["P1_nominal_invariants"] = [
        {"expr": iv["expr"], "value": iv["value"],
         "alpha_log10": iv["alpha_log10"], "n_terms": iv["n_terms"]}
        for iv in invs]
    out["P1_sigma_declared"] = max(sigma_rep, 1e-6)

    # ---- P2: chained-dollar components must NOT certify the identity ----
    real = {k: q(s).to_numpy() for k, s in
            [("Y", "GDPC1"), ("C", "PCECC96"), ("I", "GPDIC1"),
             ("G", "GCEC1")]}
    nr = min(len(v) for v in real.values())
    real = {k: v[:nr] for k, v in real.items()}
    # NX_real derived as residual would beg the question; test the 4-var scan:
    # is there any certified linear relation among chained Y,C,I,G? (There
    # should NOT be, at rounding noise -- the chain residual is material.)
    invs_r = discover_invariants(real, sigma=1e-5)
    lin_certified = [iv for iv in invs_r
                     if "**" not in iv["expr"] and "log" not in iv["expr"]
                     and "*" not in iv["expr"].replace("**", "")]
    out["P2_chained_linear_certified"] = [iv["expr"] for iv in lin_certified]

    # ---- P3: Okun (du vs real GDP growth), quarterly ----
    gdp_r = q("GDPC1")
    u = q("UNRATE", agg=True)
    growth = 100 * np.diff(np.log(gdp_r.to_numpy()))
    du = np.diff(u.to_numpy()[: len(gdp_r)])
    m = min(len(growth), len(du))
    sub = submission(growth[:m].reshape(-1, 1), du[:m], sigma=0.0)
    slope = np.polyfit(growth[:m], du[:m], 1)
    out["P3_okun"] = {"track": sub["track"], "tag": sub["tag"],
                      "expr": sub.get("expr"),
                      "ols_slope": float(slope[0]),
                      "corr": float(np.corrcoef(growth[:m], du[:m])[0, 1])}

    # ---- P4: Phillips (inflation vs unemployment level) ----
    cpi = q("CPIAUCSL", agg=True)
    infl = 400 * np.diff(np.log(cpi.to_numpy()))
    ul = u.to_numpy()[1: len(infl) + 1]
    m2 = min(len(infl), len(ul))
    sub2 = submission(ul[:m2].reshape(-1, 1), infl[:m2], sigma=0.0)
    sl2 = np.polyfit(ul[:m2], infl[:m2], 1)
    out["P4_phillips"] = {"track": sub2["track"], "tag": sub2["tag"],
                          "expr": sub2.get("expr"),
                          "ols_slope": float(sl2[0]),
                          "corr": float(np.corrcoef(ul[:m2], infl[:m2])[0, 1])}

    OUT.write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
