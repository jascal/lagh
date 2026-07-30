"""The digital twin (H2a): fit the full system once, answer every task from it.

system_id(obs) estimates {q, M, m1, m2, p (gravity exponent), tau (drag),
epoch state, COM motion}; Twin re-integrates OUR integrator with the fitted
parameters and answers all 50 registered task types from the dense fitted
trajectory in the barycenter frame. The twin is validated by predicting
held-out observations (predictions-as-harness): its max relative prediction
error is reported with every answer set.
"""
from __future__ import annotations

import numpy as np

from experiments.gravitybench.integrator import G_SI, TwoBody
from experiments.gravitybench import astronomer as ast

DAY = 86400.0

# The twin ABSTAINS above this reconstruction error. Registered 2026-07-29 from
# the SYNTHETIC BATTERY ALONE, before the benchmark effect was measured: across
# 43 dev cases the 42 whose worst answer error is <= 15% all validate at
# <= 0.0231, and the one that fails (drag, 17.3%) validates at 1.7562 -- a 76x
# gap with nothing in it. 0.05 sits 2.2x above the worst passing dev case and
# 35x below the failing one, and is the bar `test_twin_end_to_end` already used
# for the same quantity.
#
# It gates EVERY task, not only the trajectory-derived ones. The narrower
# reading is defensible -- masses and the exponent come from fits rather than
# from integrating, so a bad trajectory is only indirect evidence against them
# -- and it was rejected: this pipeline's claim is "every answer traceable to a
# fitted, PREDICTION-VALIDATED model", and a twin that cannot reproduce the
# observations it was fitted from has not validated anything it emits. The
# price is answers that would have been right anyway (measured: 7 of 13 gated
# budget instances), which is the price of not being right by luck.
TWIN_VALIDATION_MAX = 0.05


# ---------------------------------------------------------------- system id

def _com_ratio_scan(obs):
    """q = m1/m2 chosen so COM_q(t) = (q r1 + r2)/(q+1) is maximally linear in t.
    Works bound/unbound/proper-motion. Golden-section on log q."""
    t = obs["time"]
    A = np.column_stack([np.ones(len(t)), t])
    R1 = np.stack([obs[f"star1_{ax}"] for ax in "xyz"], 1)
    R2 = np.stack([obs[f"star2_{ax}"] for ax in "xyz"], 1)

    def resid(logq):
        q = np.exp(logq)
        com = (q * R1 + R2) / (q + 1.0)
        r = 0.0
        for j in range(3):
            c, *_ = np.linalg.lstsq(A, com[:, j], rcond=None)
            r += float(np.sum((com[:, j] - A @ c) ** 2))
        return r

    lo, hi = np.log(1e-3), np.log(1e3)
    for _ in range(80):
        m1 = lo + 0.382 * (hi - lo)
        m2 = lo + 0.618 * (hi - lo)
        if resid(m1) < resid(m2):
            hi = m2
        else:
            lo = m1
    q = float(np.exp(0.5 * (lo + hi)))
    com = (q * R1 + R2) / (q + 1.0)
    drift = []
    for j in range(3):
        c, *_ = np.linalg.lstsq(A, com[:, j], rcond=None)
        drift.append(c)
    drift = np.asarray(drift)          # (3, 2): [offset, velocity] per axis
    return q, drift


def system_id(obs):
    q, drift = _com_ratio_scan(obs)
    # force law from acceleration triplets: log|a_rel| = log(G M) + p log r
    rs, accs = ast._accels(obs)
    if len(rs) >= 3:
        L = np.column_stack([np.ones(len(rs)), np.log(rs)])
        c, *_ = np.linalg.lstsq(L, np.log(accs), rcond=None)
        p = float(c[1])
        GM = float(np.exp(c[0]))
    else:
        p, GM = -2.0, np.nan
    # snap to Newtonian when within tolerance (alpha definition: r^-(2+alpha))
    p_used = -2.0 if abs(p + 2.0) < 0.05 else p
    # ...and RE-FIT THE INTERCEPT AT THE SNAPPED EXPONENT. Keeping exp(c[0])
    # after moving the slope leaves GM describing a force law nobody is using:
    # the intercept was fitted against p, so switching to p_used biases it by
    # r_mid^(p - p_used). Measured on the dev orbit -- p_raw = -1.999828 snapped
    # to -2, r_mid = 1.5e11, bias 0.444%, and the total mass came back 0.445%
    # low. That is small enough to pass every direct mass tolerance and large
    # enough to fail the twin: a 0.44% mass error is a 0.22% period error, which
    # over the four-period validation window is a 2*pi*4*0.0022 = 5.6% phase
    # drift (11% at periastron on an e = 0.3 orbit). The twin's whole premise is
    # that it can be integrated forward, so a bias the direct answers survive
    # is one the trajectory does not.
    if p_used != p and len(rs) >= 3:
        GM = float(np.exp(np.mean(np.log(accs) - p_used * np.log(rs))))
    M = GM / G_SI
    m2 = M / (1.0 + q)
    m1 = M - m2
    # drag from envelope decay of the separation
    tau = ast.task_drag_tau(obs)
    lam_ok = tau is not None and tau > 0
    # epoch: center of a strict-cadence triplet nearest the median time
    t = obs["time"]
    R1 = np.stack([obs[f"star1_{ax}"] for ax in "xyz"], 1)
    R2 = np.stack([obs[f"star2_{ax}"] for ax in "xyz"], 1)
    best = None
    trips = obs.get("triplets")
    if trips is not None:
        for ta, tb_, tc in np.asarray(trips):
            ia = int(np.argmin(np.abs(t - ta)))
            ib = int(np.argmin(np.abs(t - tb_)))
            ic = int(np.argmin(np.abs(t - tc)))
            dt1, dt2 = t[ib] - t[ia], t[ic] - t[ib]
            if dt1 <= 0 or dt2 <= 0 or abs(dt1 - dt2) > 0.05 * max(dt1, dt2):
                continue
            score = abs(t[ib] - np.median(t))
            if best is None or score < best[0]:
                v1 = (R1[ic] - R1[ia]) / (dt1 + dt2)
                v2 = (R2[ic] - R2[ia]) / (dt1 + dt2)
                best = (score, t[ib], R1[ib], R2[ib], v1, v2)
    else:
        for i in range(1, len(t) - 1):
            dt1, dt2 = t[i] - t[i - 1], t[i + 1] - t[i]
            if dt1 <= 0 or dt2 <= 0 or abs(dt1 - dt2) > 0.02 * max(dt1, dt2):
                continue
            score = abs(t[i] - np.median(t))
            if best is None or score < best[0]:
                v1 = (R1[i + 1] - R1[i - 1]) / (dt1 + dt2)
                v2 = (R2[i + 1] - R2[i - 1]) / (dt1 + dt2)
                best = (score, t[i], R1[i], R2[i], v1, v2)
    if best is None:
        raise RuntimeError("no strict-cadence triplet available for epoch state")
    _, t0, r1_0, r2_0, v1_0, v2_0 = best
    # per-triplet specific orbital energies (for the direct is_bound verdict)
    energies = []
    trips = obs.get("triplets")
    if trips is not None:
        for ta, tb_, tc in np.asarray(trips):
            ia = int(np.argmin(np.abs(t - ta)))
            ib = int(np.argmin(np.abs(t - tb_)))
            ic = int(np.argmin(np.abs(t - tc)))
            dt1, dt2 = t[ib] - t[ia], t[ic] - t[ib]
            if dt1 <= 0 or dt2 <= 0:
                continue
            drel_a = R2[ia] - R1[ia]
            drel_c = R2[ic] - R1[ic]
            vrel = (drel_c - drel_a) / (dt1 + dt2)
            rmid = float(np.sqrt(((R2[ib] - R1[ib]) ** 2).sum()))
            energies.append(0.5 * float((vrel ** 2).sum()) - G_SI * M / rmid)
    return {"q": q, "M": M, "m1": m1, "m2": m2, "p": p_used, "p_raw": p,
            "triplet_energies": np.asarray(energies),
            "alpha": -(p + 2.0), "tau": tau if lam_ok else None,
            "drift": drift, "epoch": {"t0": float(t0), "r1": r1_0, "r2": r2_0,
                                      "v1": v1_0, "v2": v2_0}}


# ---------------------------------------------------------------------- twin

class Twin:
    def __init__(self, state, maxtime, n_steps=200_000):
        self.s = state
        e = state["epoch"]
        self.t0 = e["t0"]
        self.maxtime = float(maxtime)
        mod = None if state["p"] == -2.0 else state["p"]
        # forward from epoch
        fwd = TwoBody(state["m1"], state["m2"], e["r1"], e["r2"],
                      e["v1"], e["v2"],
                      drag_tau=state["tau"], mod_gravity_exponent=mod)
        n_f = max(int(n_steps * (maxtime - self.t0) / maxtime), 1000)
        fwd.run(maxtime - self.t0, n_steps=n_f)
        # backward from epoch (time-reversal: negate velocities, drag sign flips
        # is NOT exact under drag -- backward segment is skipped when tau is set)
        if state["tau"] is None:
            bwd = TwoBody(state["m1"], state["m2"], e["r1"], e["r2"],
                          -e["v1"], -e["v2"], mod_gravity_exponent=mod)
            n_b = max(int(n_steps * self.t0 / maxtime), 1000)
            bwd.run(self.t0, n_steps=n_b)
            T = np.concatenate([self.t0 - bwd._T[::-1], self.t0 + fwd._T[1:]])
            R1 = np.vstack([bwd._R1[::-1], fwd._R1[1:]])
            R2 = np.vstack([bwd._R2[::-1], fwd._R2[1:]])
            V1 = np.vstack([-bwd._V1[::-1], fwd._V1[1:]])
            V2 = np.vstack([-bwd._V2[::-1], fwd._V2[1:]])
        else:
            T = self.t0 + fwd._T
            R1, R2 = fwd._R1, fwd._R2
            V1, V2 = fwd._V1, fwd._V2
        self.T, self.R1, self.R2, self.V1, self.V2 = T, R1, R2, V1, V2
        # barycenter frame
        q = state["q"]
        com = (q * R1 + R2) / (q + 1.0)
        vcom = (q * V1 + V2) / (q + 1.0)
        self.R1c, self.R2c = R1 - com, R2 - com
        self.V1c, self.V2c = V1 - vcom, V2 - vcom
        self.rel = R2 - R1
        self.vrel = V2 - V1
        self.sep = np.sqrt((self.rel ** 2).sum(1))

    def gated_answer(self, task, obs, *, max_validation=TWIN_VALIDATION_MAX):
        """The answer, or an ABSTAIN when the twin failed its own validation.

        Returns (answer, validation, refusal). `answer` is None exactly when
        `refusal` is set, so a caller cannot read a number without having been
        told the model behind it does not reproduce the data.

        This is the gate the instrument was already computing and not using:
        `validate` was recorded next to every answer in the read journal and
        never consulted. Two wrong answers in that journal
        (max_velocity_star1 and max_angular_velocity_star1 on
        9p6_M_3p1_M_Proper_Motion2) came from a twin validating at 0.62 -- 62%
        reconstruction error -- and a third, earlier, from one at 4.09.
        """
        v = self.validate(obs)
        if not np.isfinite(v) or v > max_validation:
            return None, v, (f"twin validation {v:.4g} > {max_validation:g}: "
                             "the fitted model does not reproduce the "
                             "observations it was fitted from, so no answer "
                             "read off it is validated")
        return self.answer(task), v, None

    # ---- validation (predictions-as-harness) ----
    def validate(self, obs_holdout):
        errs = []
        for s, R in (("star1", (self.T, self.R1)), ("star2", (self.T, self.R2))):
            for j, ax in enumerate("xyz"):
                pred = np.interp(obs_holdout["time"], R[0], R[1][:, j])
                scale = np.max(np.abs(obs_holdout[f"{s}_{ax}"])) + 1e-300
                errs.append(np.max(np.abs(pred - obs_holdout[f"{s}_{ax}"]) / scale))
        return float(np.max(errs))

    # ---- one-period window for "over the orbit" quantities ----
    def _period_slice(self):
        P = self.period()
        if P is None:
            return slice(0, len(self.T))
        i0 = np.searchsorted(self.T, self.t0)
        i1 = np.searchsorted(self.T, min(self.t0 + P, self.T[-1]))
        return slice(max(i0, 0), max(i1, i0 + 100))

    def period(self):
        # angle method on the twin's dense relative trajectory (works from ~1
        # revolution; separation minima needed two periastra and returned None
        # on 1.5-period windows -- measured battery failure)
        dc = self.rel - self.rel.mean(0)
        _, _, Vt = np.linalg.svd(dc[:: max(1, len(dc) // 2000)],
                                 full_matrices=False)
        u, v = dc @ Vt[0], dc @ Vt[1]
        theta = np.unwrap(np.arctan2(v, u))
        swept = theta - theta[0]
        total = abs(swept[-1])
        if total < 2 * np.pi:
            return None
        # SAME-ANGLE CROSSING: 2*pi*span/swept is exact only over integer
        # revolutions (measured 17% bias on 1.5 eccentric revs); the time for
        # theta to advance exactly 2*pi is exact for any eccentricity
        sgn = np.sign(swept[-1])
        target = sgn * 2 * np.pi
        idx = np.searchsorted(sgn * swept, sgn * target)
        if idx >= len(self.T):
            return None
        f = (target - swept[idx - 1]) / (swept[idx] - swept[idx - 1] + 1e-300)
        t_cross = self.T[idx - 1] + f * (self.T[idx] - self.T[idx - 1])
        return float(t_cross - self.T[0])

    def is_bound(self):
        # energy sign DIRECTLY from observed triplets (stored by system_id):
        # specific orbital energy 0.5 v_rel^2 - GM/r per triplet, majority sign.
        # Twin-midpoint energy inherited integration error on flybys (measured).
        tv = self.s.get("triplet_energies")
        if tv is not None and len(tv):
            return bool(np.median(tv) < 0)
        mu = self.s["m1"] * self.s["m2"] / self.s["M"]
        i = len(self.T) // 2
        v2 = float((self.vrel[i] ** 2).sum())
        K = 0.5 * mu * v2
        U = -G_SI * self.s["m1"] * self.s["m2"] / self.sep[i]
        return bool(K + U < 0)

    # ---- the dispatch ----
    def answer(self, task):
        s = self.s
        sl = self._period_slice()
        P = self.period()

        def star(i):
            return (self.R1c, self.V1c, s["m1"]) if i == 1 else \
                   (self.R2c, self.V2c, s["m2"])

        def speeds(i):
            _, V, _ = star(i)
            return np.sqrt((V[sl] ** 2).sum(1))

        def accs(i):
            _, V, _ = star(i)
            dt = np.gradient(self.T[sl])
            return np.sqrt(((np.gradient(V[sl], axis=0) / dt[:, None]) ** 2).sum(1))

        def angvel(i):
            R, V, _ = star(i)
            r2 = (R[sl] ** 2).sum(1)
            Lz = np.cross(R[sl], V[sl])
            return np.sqrt((Lz ** 2).sum(1)) / np.maximum(r2, 1e-300)

        if task == "period":
            return P
        if task == "mass_ratio":
            return s["q"]
        if task == "total_mass":
            return s["M"]
        if task == "mass_star1":
            return s["m1"]
        if task == "mass_star2":
            return s["m2"]
        if task == "mass_largest_star":
            return max(s["m1"], s["m2"])
        if task == "reduced_mass":
            return s["m1"] * s["m2"] / s["M"]
        if task == "modified_gravity_power_law":
            return s["alpha"]
        if task == "linear_drag":
            return s["tau"]
        if task == "is_bound":
            return self.is_bound()
        if task == "kepler_3rd_law":
            return bool(abs(s["p"] + 2.0) < 0.05)
        if task == "virial_theorem":
            if not self.is_bound() or s["tau"] is not None or s["p"] != -2.0:
                return False
            return True
        if task in ("apoastron", "periastron"):
            r = self.sep[sl]
            return float(np.max(r)) if task == "apoastron" else float(np.min(r))
        if task == "eccentricity":
            r = self.sep[sl]
            ra, rp = float(np.max(r)), float(np.min(r))
            return (ra - rp) / (ra + rp)
        if task in ("semi_major_axis", "semi_minor_axis"):
            r = self.sep[sl]
            a = 0.5 * (np.max(r) + np.min(r))
            if task == "semi_major_axis":
                return float(a)
            e = (np.max(r) - np.min(r)) / (np.max(r) + np.min(r))
            return float(a * np.sqrt(1 - e ** 2))
        if task.startswith("semi_major_axis_star") or \
                task.startswith("semi_minor_axis_star"):
            i = int(task[-1])
            R, _, _ = star(i)
            d = np.sqrt((R[sl] ** 2).sum(1))
            ra, rp = float(np.max(d)), float(np.min(d))
            a = 0.5 * (ra + rp)
            if "major" in task:
                return a
            e = (ra - rp) / (ra + rp)
            return a * np.sqrt(1 - e ** 2)
        if task.startswith("max_velocity_star"):
            return float(np.max(speeds(int(task[-1]))))
        if task.startswith("min_velocity_star"):
            return float(np.min(speeds(int(task[-1]))))
        if task.startswith("max_momentum_star"):
            i = int(task[-1])
            return float(np.max(speeds(i)) * star(i)[2])
        if task.startswith("min_momentum_star"):
            i = int(task[-1])
            return float(np.min(speeds(i)) * star(i)[2])
        if task.startswith("max_acceleration_star"):
            return float(np.max(accs(int(task[-1]))))
        if task.startswith("min_acceleration_star"):
            return float(np.min(accs(int(task[-1]))))
        if task.startswith("max_angular_velocity_star"):
            return float(np.max(angvel(int(task[-1]))))
        if task.startswith("min_angular_velocity_star"):
            return float(np.min(angvel(int(task[-1]))))
        if task.startswith("avg_distance_COM_star"):
            i = int(task[-1])
            R, _, _ = star(i)
            return float(np.mean(np.sqrt((R[sl] ** 2).sum(1))))
        if task.startswith("orbital_area_star"):
            i = int(task[-1])
            R, _, _ = star(i)
            x, y = R[sl][:, 0], R[sl][:, 1]     # planar orbits (z ~ 0)
            return float(0.5 * np.abs(np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))))
        if task in ("area_swept_over_time_apo", "area_swept_over_time_peri"):
            r = self.sep[sl]
            i = int(np.argmax(r)) if task.endswith("apo") else int(np.argmin(r))
            cr = np.cross(self.rel[sl][i], self.vrel[sl][i])
            return float(0.5 * np.sqrt((cr ** 2).sum()))
        if task == "specific_angular_momentum":
            i = len(self.T) // 2
            cr = np.cross(self.rel[i], self.vrel[i])
            return float(np.sqrt((cr ** 2).sum()))
        if task == "K+U":
            i = len(self.T) // 2
            K = 0.5 * s["m1"] * float((self.V1c[i] ** 2).sum()) + \
                0.5 * s["m2"] * float((self.V2c[i] ** 2).sum())
            U = -G_SI * s["m1"] * s["m2"] / float(self.sep[i])
            return K + U
        if task == "roche_lobe_radius":
            q = s["q"]
            r = self.sep[sl]
            a = 0.5 * (np.max(r) + np.min(r))
            return float(a * 0.49 * q ** (2 / 3) /
                         (0.6 * q ** (2 / 3) + np.log(1 + q ** (1 / 3))))
        if task == "multiply_mass_period":
            return (P / (21.0 * DAY)) ** 2 if P else None
        if task == "time_fraction_acceleraton_below_mean":
            a = accs(1)
            return float(np.mean(a < np.mean(a)))
        if task.startswith("travel_time_orbital_"):
            frac = 0.2 if "20per" in task else 0.7
            R, _, _ = star(1)
            r = self.sep[sl]
            i_peri = int(np.argmin(r))
            seg = R[sl][i_peri:]
            d = np.sqrt((np.diff(seg, axis=0) ** 2).sum(1))
            arc = np.concatenate([[0.0], np.cumsum(d)])
            # total path over one period from pericenter
            if P is None:
                return None
            t_seg = self.T[sl][i_peri:] - self.T[sl][i_peri]
            j_end = np.searchsorted(t_seg, P)
            total = arc[min(j_end, len(arc) - 1)]
            j = np.searchsorted(arc, frac * total)
            return float(t_seg[min(j, len(t_seg) - 1)])
        raise KeyError(f"no playbook entry for task {task!r}")
