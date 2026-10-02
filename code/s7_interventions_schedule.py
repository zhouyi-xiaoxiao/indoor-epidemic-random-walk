#!/usr/bin/env python
"""s7_interventions_schedule.py -- checks of the duty-cycle rule and of the time-average rule (Section 7).

(A) Duty cycle (office scene, 1 m kernel, N = 100, D0 = 1 and 10).  The room is occupied during the first
    fraction f = 8/24 of every day (movement and transmission on); otherwise people are frozen and nothing is
    transmitted; recovery runs all the time.  Proposition s7_interventions:prop:duty shows that the mean-field
    threshold is R0(gamma/f) = 1 exactly.  As a COUNT of secondary cases the substitution gamma -> gamma/f is
    an approximation; its error depends on the phase of the day at which the index case becomes infectious.
    Here the secondary cases of the index (only the index transmits) are counted
      (i)  for an index that becomes infectious at the start of a room period  (as in
           code/bsc_sim/experiments/11_dwell_time.py, which used one batch of 20,000 runs, seed 1100, and
           found 0.918; the re-check found 0.944 +- 0.006 with a second simulator), and
      (ii) for an index that becomes infectious at a time uniformly distributed over the room period
           (the situation of every secondary case), implemented by shifting the schedule by 20 equally spaced
           phases.
    200,000 runs each, in batches of 10,000 with distinct seeds.  The exact expected in-room infectious
    times are computed in closed form for comparison.

(B) Time average (metro scene, D0 = 10, N = 310 fixed, 1 m kernel).  This part's experiment
    (code/bsc_sim/experiments/09_metro_time.py; arms: constant peak values / hourly schedule / constants
    equal to the time averages) was run at D0 = 1 and 100 only; the same three arms are run here at D0 = 10
    with the same schedule definition, 2000 epidemics per arm.

(C) Deterministic check of Proposition s7_interventions:prop:duty: R0^(f)(D0) = f R0(f D0), and the one-day
    propagator of the linearised mean-field model equals exp(f J_f).

(D) Exact pair-level expectation with the on/off schedule.  The counts of (A)
    are Monte Carlo estimates of a quantity that the pair formalism of Section 5 gives exactly.  Let
    A = gamma + H - M(+)M be the pair operator (symmetric positive definite), u(z) the probability that the
    given person escapes infection when the pair starts in state z, c = 1 - exp(-gamma (1-f)) and
    d = exp(-gamma (1-f)).  During a room period of length f the pair evolves and the hazard acts; during the
    rest of the day nothing moves and only recovery proceeds.  For an index case that becomes infectious at
    the start of a room period,
        u0 = gamma A^-1 (1 - e^{-A f} 1) + e^{-A f} (c 1 + d u0),
    a fixed-point equation with contraction factor <= exp(-gamma), solved by iteration; for an index that
    becomes infectious at a time uniform over the room period the same expression is averaged over the
    remaining room time r in (0, f):
        u_bar = gamma [A^-1 1 - (1/f) A^-1 A^-1 (1 - e^{-A f} 1)] + (1/f) A^-1 (v - e^{-A f} v),  v = c 1 + d u0.
    The expected number of secondary cases is (N - 1) * mean(1 - u).  Office (f = 8/24) and classroom
    (f = 1/24), D0 = 1 and 10, 1 m kernel.

(E) Larger Monte Carlo samples of the two office counts of (A), with other seeds:
    2,000,000 runs each at D0 = 1 and 800,000 each at D0 = 10, compared with the exact values of (D).

Output: ../data/s7_interventions_schedule.json
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))          # repository root
sys.path.insert(0, os.path.join(ROOT, "code", "bsc_sim", "src"))
from bsc_sim import theory as th                                       # noqa: E402
from bsc_sim.scenes import load_scene                                  # noqa: E402
from bsc_sim.sim import mean_ci, simulate, summarize                   # noqa: E402

BETA, GAMMA = 0.5, 0.14
ELL_C = 1.0
PATH = os.path.join(HERE, "..", "data", "s7_interventions_schedule.json")


def in_room_time_start(gamma, f):
    """Expected in-room infectious time of a case that becomes infectious at the start of a room period."""
    return (1 - np.exp(-gamma * f)) / gamma / (1 - np.exp(-gamma))


def in_room_time_uniform(gamma, f):
    """Same, for a case that becomes infectious at a time uniform over the room period [0, f)."""
    a = (1 - (1 - np.exp(-gamma * f)) / (gamma * f)) / gamma
    b = (np.exp(gamma * f) - 1) / (gamma * f) * (1 - np.exp(-gamma * f)) / gamma \
        * np.exp(-gamma) / (1 - np.exp(-gamma))
    return a + b


def duty_cycle(out):
    f = 8.0 / 24.0
    n_batch, n_phase = 10000, 20
    for D0 in (1.0, 10.0):
        key = f"duty|office|D0={D0:g}"
        if key in out:
            continue
        t0 = time.time()
        sc = load_scene("office", D0=D0)
        P = th.contact_kernel(sc, ELL_C / sc.a_m)
        frozen = np.zeros(sc.M)
        row = dict(scene="office", D0=D0, f=f, N=sc.N,
                   R0_continuous=th.R0_dense(sc, BETA, GAMMA, P=P),
                   R0_duty=th.R0_dense(sc, BETA, GAMMA / f, P=P),
                   wellmixed_duty=BETA * sc.q.mean() * f / GAMMA)
        p = th.pair_infection_probability(sc, BETA, GAMMA / f, P=P)
        row["R1bar_duty"] = float((sc.N - 1) * p.mean())
        p = th.pair_infection_probability(sc, BETA, GAMMA, P=P)
        row["R1bar_continuous"] = float((sc.N - 1) * p.mean())
        # (i) index infectious from the start of a room period
        sched = dict(period=1.0, segments=[(f, sc.D, 1.0), (1.0, frozen, 0.0)])
        off = []
        for b in range(n_phase):
            r = simulate(sc, BETA, GAMMA, n_rep=n_batch, seed=73000 + b + int(100 * D0), P=P, gmax=0,
                         schedule=sched, t_max=400.0, dt_out=50.0)
            off.append(r.index_offspring)
        batch_means = [float(o.mean()) for o in off]
        m, lo, hi = mean_ci(np.concatenate(off))
        row.update(count_start=m, count_start_lo=lo, count_start_hi=hi, count_start_se=(hi - m) / 1.96,
                   count_start_batch_means=batch_means, n_start=n_batch * n_phase)
        # (ii) index infectious from a uniformly distributed time within the room period:
        # shift the schedule so that at t = 0 a time s of the room period has already elapsed.
        off = []
        for b in range(n_phase):
            s = (b + 0.5) * f / n_phase
            sched_s = dict(period=1.0, segments=[(f - s, sc.D, 1.0), (1.0 - s, frozen, 0.0), (1.0, sc.D, 1.0)])
            r = simulate(sc, BETA, GAMMA, n_rep=n_batch, seed=74000 + b + int(100 * D0), P=P, gmax=0,
                         schedule=sched_s, t_max=400.0, dt_out=50.0)
            off.append(r.index_offspring)
        m, lo, hi = mean_ci(np.concatenate(off))
        row.update(count_uniform_phase=m, count_uniform_phase_lo=lo, count_uniform_phase_hi=hi,
                   count_uniform_phase_se=(hi - m) / 1.96,
                   count_uniform_phase_by_phase=[float(o.mean()) for o in off], n_uniform_phase=n_batch * n_phase)
        row.update(in_room_time_rule=f / GAMMA, in_room_time_start=in_room_time_start(GAMMA, f),
                   in_room_time_uniform=in_room_time_uniform(GAMMA, f), seconds=time.time() - t0)
        out[key] = row
        json.dump(out, open(PATH, "w"), indent=1, default=float)
        print(key, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()
                    if not isinstance(v, list)}, flush=True)


def duty_cycle_exact(out, tol=1e-12):
    """(D) exact expected number of secondary cases of a lone index case under the on/off schedule."""
    import scipy.sparse.linalg as spla
    for name, f in (("office", 8.0 / 24.0), ("classroom", 1.0 / 24.0)):
        for D0 in (1.0, 10.0):
            key = f"duty_exact|{name}|D0={D0:g}"
            if key in out:
                continue
            t0 = time.time()
            sc = load_scene(name, D0=D0)
            P = th.contact_kernel(sc, ELL_C / sc.a_m)
            A, h, N, M = th._pair_operator(sc, BETA, GAMMA, None, None, None, P, None)
            assert abs(A - A.T).max() < 1e-12
            one = np.ones(M * M)
            Af = (-f) * A.tocsc()
            ex = lambda v: spla.expm_multiply(Af, v)                 # e^{-A f} v
            c, d = 1 - np.exp(-GAMMA * (1 - f)), np.exp(-GAMMA * (1 - f))
            e1 = ex(one)
            g = th._cg(A, one - e1, rtol=1e-12)                     # A^-1 (1 - e^{-Af} 1)
            a = GAMMA * g
            u0 = one.copy()
            for it in range(2000):
                new = a + c * e1 + d * ex(u0)
                err = float(np.abs(new - u0).max())
                u0 = new
                if err < tol:
                    break
            v = c * one + d * u0
            w1 = th._cg(A, one, rtol=1e-12)
            ubar = GAMMA * (w1 - th._cg(A, g, rtol=1e-12) / f) + th._cg(A, v - ex(v), rtol=1e-12) / f
            # consistency: with f -> the whole day the fixed point is the continuous pair solve
            row = dict(scene=name, D0=D0, f=f, N=N, pair_states=int(M * M), iterations=it + 1, fixed_point_residual=err,
                       exact_start=float((N - 1) * (1 - u0).mean()), exact_uniform_phase=float((N - 1) * (1 - ubar).mean()),
                       rule_gamma_over_f=float((N - 1) * th.pair_infection_probability(sc, BETA, GAMMA / f, P=P).mean()),
                       seconds=time.time() - t0)
            row["start_over_rule_pct"] = 100 * (row["exact_start"] / row["rule_gamma_over_f"] - 1)
            row["uniform_over_rule_pct"] = 100 * (row["exact_uniform_phase"] / row["rule_gamma_over_f"] - 1)
            mc = out.get(f"duty|{name}|D0={D0:g}")
            if mc is not None:
                row["z_mc_start"] = (mc["count_start"] - row["exact_start"]) / mc["count_start_se"]
                row["z_mc_uniform_phase"] = (mc["count_uniform_phase"] - row["exact_uniform_phase"]) / mc["count_uniform_phase_se"]
            out[key] = row
            json.dump(out, open(PATH, "w"), indent=1, default=float)
            print(key, {k: (round(v_, 5) if isinstance(v_, float) else v_) for k, v_ in row.items()}, flush=True)
    # control of the method: f = 1 (room always occupied) must reproduce the continuous pair solve
    key = "duty_exact|control_f1|office|D0=1"
    if key not in out:
        sc = load_scene("office", D0=1.0)
        P = th.contact_kernel(sc, ELL_C / sc.a_m)
        A, h, N, M = th._pair_operator(sc, BETA, GAMMA, None, None, None, P, None)
        one = np.ones(M * M)
        f = 1.0
        ex = lambda v: spla.expm_multiply((-f) * A.tocsc(), v)
        e1 = ex(one)
        a = GAMMA * th._cg(A, one - e1, rtol=1e-12)
        u0 = one.copy()
        for it in range(2000):
            new = a + ex(u0)                                         # c = 0, d = 1
            err = float(np.abs(new - u0).max())
            u0 = new
            if err < tol:
                break
        ref = float((N - 1) * th.pair_infection_probability(sc, BETA, GAMMA, P=P).mean())
        out[key] = dict(exact_start=float((N - 1) * (1 - u0).mean()), continuous_pair_solve=ref,
                        abs_diff=abs(float((N - 1) * (1 - u0).mean()) - ref))
        json.dump(out, open(PATH, "w"), indent=1, default=float)
        print(key, out[key], flush=True)


def identity(out):
    """Deterministic check of Proposition s7_interventions:prop:duty(c): R0^(f)(D0) = f R0(f D0), and of the
    exact linearised solution I(k) = exp(k f J_f) I(0) against a direct product of matrix exponentials."""
    import scipy.linalg as sla
    f = 8.0 / 24.0
    res = {}
    for name in ("office", "classroom"):
        for D0 in (1.0, 10.0):
            sc = load_scene(name, D0=D0)
            P = th.contact_kernel(sc, ELL_C / sc.a_m)
            lhs = th.R0_dense(sc, BETA, GAMMA / f, P=P)
            sc2 = load_scene(name, D0=f * D0)
            rhs = f * th.R0_dense(sc2, BETA, GAMMA, P=th.contact_kernel(sc2, ELL_C / sc2.a_m))
            # one-day propagator: exp((1-f)(-gamma)) exp(f (M + F - gamma)) versus exp(f J_f)
            L = th.generator(sc).toarray()
            F = (P.T @ np.diag(BETA * sc.q)) if not hasattr(P, "toarray") else (P.T.toarray() @ np.diag(BETA * sc.q))
            B = L + F
            day = np.exp(-GAMMA * (1 - f)) * sla.expm(f * (B - GAMMA * np.eye(sc.M)))
            Jf = B - (GAMMA / f) * np.eye(sc.M)
            direct = sla.expm(f * Jf)
            growth = float(np.log(np.max(np.abs(np.linalg.eigvals(day)))))
            res[f"{name}|D0={D0:g}"] = dict(R0_duty=lhs, f_times_R0_at_fD0=rhs, abs_diff=abs(lhs - rhs),
                                            propagator_max_abs_diff=float(np.abs(day - direct).max()),
                                            daily_log_growth=growth,
                                            f_times_lambda1_at_gamma_over_f=float(
                                                f * th.growth_rate(sc, BETA, GAMMA / f, P=P)))
    out["identity"] = res
    json.dump(out, open(PATH, "w"), indent=1, default=float)
    for k, v in res.items():
        print("identity", k, {kk: round(vv, 12) for kk, vv in v.items()}, flush=True)


# ---- (B) metro schedule, copied from code/bsc_sim/experiments/09_metro_time.py ----------------------
RHO_BASE, RHO_PEAK, SIGMA_H = 1.5, 3.1, 1.0


def rho_of_hour(h):
    h = np.asarray(h, float)
    r = RHO_BASE + np.zeros_like(h)
    for pk in (8.0, 18.0):
        d = np.minimum(np.abs(h - pk), 24 - np.abs(h - pk))
        r = r + RHO_PEAK * np.exp(-d ** 2 / (2 * SIGMA_H ** 2))
    return r


def metro(out, D0=10.0, n_rep=2000):
    hours = np.arange(24) + 0.5
    rho_h = rho_of_hour(hours)
    m_h = np.where(rho_h >= 3.0, 1.0, 0.8)
    sc = load_scene("metro", D0=D0)
    P = th.contact_kernel(sc, ELL_C / sc.a_m)
    zones = sc.meta["zones"]
    coef = np.array([zones[z]["D_coef"] for z in sc.site_zone])
    expo = np.array([zones[z]["D_exp"] for z in sc.site_zone])
    D_h = [D0 * coef / r ** expo for r in rho_h]
    D_mean = np.mean(D_h, axis=0)
    m_mean = float(m_h.mean())
    arms = {
        "peak": dict(D=sc.D, schedule=None, m=1.0),
        "varying": dict(D=D_h[0], schedule=dict(period=1.0, segments=[((k + 1) / 24.0, D_h[k], float(m_h[k]))
                                                                      for k in range(24)]), m=None),
        "mean": dict(D=D_mean, schedule=dict(period=1.0, segments=[(1.0, D_mean, m_mean)]), m=m_mean),
    }
    for arm, cfg in arms.items():
        key = f"metro|D0={D0:g}|{arm}"
        if key in out:
            continue
        t0 = time.time()
        res = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=900, P=P, D=cfg["D"], schedule=cfg["schedule"],
                       t_max=300.0, dt_out=0.25)
        row = dict(D0=D0, arm=arm, n_rep=n_rep, N=sc.N, mean_D_over_D0=float(np.mean(cfg["D"]) / D0),
                   m_mean=m_mean)
        if cfg["m"] is not None:
            L = th.generator(sc, cfg["D"])
            row["R0_ngm"] = th.R0_dense(sc, BETA * cfg["m"], GAMMA, P=P, L=L)
            p = th.pair_infection_probability(sc, BETA * cfg["m"], GAMMA, P=P, L=L)
            row["R1_pair_uniform"] = float((sc.N - 1) * p.mean())
        row.update(summarize(res))
        row["seconds"] = time.time() - t0
        out[key] = row
        json.dump(out, open(PATH, "w"), indent=1, default=float)
        print(key, {k: round(row[k], 4) for k in ("p_major", "attack_major_mean", "peak_prev_major_mean",
                                                    "peak_time_major_mean")}, flush=True)


def duty_cycle_large(out):
    """(E) Larger Monte Carlo samples of the counts of (A), with other seeds.
    The 200,000-run samples of (A) lie 2.0 and 2.5 standard errors above the exact values of (D) at D0 = 1.
    Fresh samples: office, D0 = 1: 2,000,000 runs for each of the two cases; D0 = 10: 800,000 runs each.
    Checkpointed per batch."""
    f = 8.0 / 24.0
    plan = {1.0: (40, 50000), 10.0: (40, 20000)}
    for D0, (n_phase, n_batch) in plan.items():
        key = f"duty_large|office|D0={D0:g}"
        row = out.get(key, dict(scene="office", D0=D0, f=f, n_batch=n_batch, n_phase=n_phase, start_sum=0.0,
                                start_sumsq=0.0, start_n=0, start_batches=[], uni_sum=0.0, uni_sumsq=0.0,
                                uni_n=0, uni_batches=[], seconds=0.0))
        if row.get("done"):
            continue
        sc = load_scene("office", D0=D0)
        P = th.contact_kernel(sc, ELL_C / sc.a_m)
        frozen = np.zeros(sc.M)
        sched = dict(period=1.0, segments=[(f, sc.D, 1.0), (1.0, frozen, 0.0)])
        while len(row["start_batches"]) < n_phase:
            b = len(row["start_batches"])
            t0 = time.time()
            r = simulate(sc, BETA, GAMMA, n_rep=n_batch, seed=830000 + b + int(1000 * D0), P=P, gmax=0,
                         schedule=sched, t_max=400.0, dt_out=50.0)
            x = np.asarray(r.index_offspring, float)
            row["start_sum"] += float(x.sum()); row["start_sumsq"] += float((x * x).sum()); row["start_n"] += len(x)
            row["start_batches"].append(float(x.mean()))
            row["seconds"] += time.time() - t0
            out[key] = row
            json.dump(out, open(PATH, "w"), indent=1, default=float)
        while len(row["uni_batches"]) < n_phase:
            b = len(row["uni_batches"])
            s = (b + 0.5) * f / n_phase
            sched_s = dict(period=1.0, segments=[(f - s, sc.D, 1.0), (1.0 - s, frozen, 0.0), (1.0, sc.D, 1.0)])
            t0 = time.time()
            r = simulate(sc, BETA, GAMMA, n_rep=n_batch, seed=840000 + b + int(1000 * D0), P=P, gmax=0,
                         schedule=sched_s, t_max=400.0, dt_out=50.0)
            x = np.asarray(r.index_offspring, float)
            row["uni_sum"] += float(x.sum()); row["uni_sumsq"] += float((x * x).sum()); row["uni_n"] += len(x)
            row["uni_batches"].append(float(x.mean()))
            row["seconds"] += time.time() - t0
            out[key] = row
            json.dump(out, open(PATH, "w"), indent=1, default=float)
        for tag in ("start", "uni"):
            n = row[f"{tag}_n"]
            m = row[f"{tag}_sum"] / n
            var = (row[f"{tag}_sumsq"] - n * m * m) / (n - 1)
            row[f"count_{tag}"] = m
            row[f"count_{tag}_se"] = float(np.sqrt(var / n))
        ex = out.get(f"duty_exact|office|D0={D0:g}")
        if ex is not None:
            row["exact_start"] = ex["exact_start"]
            row["exact_uniform_phase"] = ex["exact_uniform_phase"]
            row["z_start"] = (row["count_start"] - ex["exact_start"]) / row["count_start_se"]
            row["z_uni"] = (row["count_uni"] - ex["exact_uniform_phase"]) / row["count_uni_se"]
        row["done"] = True
        out[key] = row
        json.dump(out, open(PATH, "w"), indent=1, default=float)
        print(key, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()
                    if not isinstance(v, list)}, flush=True)


if __name__ == "__main__":
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    duty_cycle(out)
    metro(out)
    identity(out)
    duty_cycle_exact(out)
    duty_cycle_large(out)
