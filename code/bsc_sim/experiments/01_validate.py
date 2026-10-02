#!/usr/bin/env python
"""Validation of the simulator (all tests write to data/01_validation.json).

A. C core vs pure-Python reference implementation (same model, independent
   code and RNG) on a small heterogeneous lattice, for the same-cell kernel,
   a radius-1 kernel and the per-cell rule.
B. Well-mixed limit: one cell, N people -> exact final-size distribution of
   the Markovian SIR model (dynamic programming on the embedded jump chain).
C. Exposure identity: the time-integrated potential infection rate of an
   index case born at x has expectation (1^T K)_x; with the index drawn from
   the Perron vector its expectation is R0 = rho(K).
D. Pair theory: mean number of secondary cases of the index case when only the
   index transmits (gmax = 0) vs the exact pair-lattice solve.
E. Per-cell rule (rate beta q n_I / n_total in the same cell): sub-critical
   at the occupancies of the scenes.
F. Mean-field limit: with many people per cell and fast mixing the stochastic
   epidemic curve converges to the lattice reaction-diffusion ODE.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.refsim import simulate_one                # noqa: E402
from bsc_sim.scenes import Scene, load_scene, uniform_room  # noqa: E402
from bsc_sim.sim import mean_ci, simulate              # noqa: E402

BETA, GAMMA = 0.5, 0.14
OUT = os.path.join(HERE, "..", "data", "01_validation.json")
results = {}


def small_scene():
    """4 x 3 lattice with one barrier, two q zones and two D zones."""
    zone = np.array([list("HHLL"), list("HBLL"), list("HHLL")])
    access = zone != "B"
    D = np.where(zone == "H", 0.4, 1.5) * access
    q = np.where(zone == "H", 2.0, 0.6) * access
    return Scene(name="small", nx=4, ny=3, a_m=1.0, zone=zone, access=access,
                 D_grid=D, q_grid=q, N=10, meta={})


def two_sample(a, b):
    """z-score of the difference of two sample means."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    return float((a.mean() - b.mean()) / se) if se > 0 else 0.0


def chi2_hist(a, b, nbins):
    """Chi-square statistic comparing two samples of integers in 0..nbins-1
    (bins with fewer than 10 pooled counts are merged)."""
    ha = np.bincount(a, minlength=nbins).astype(float)
    hb = np.bincount(b, minlength=nbins).astype(float)
    # merge sparse bins
    A, B, ca, cb = [], [], 0.0, 0.0
    for x, y in zip(ha, hb):
        ca += x
        cb += y
        if ca + cb >= 10:
            A.append(ca)
            B.append(cb)
            ca = cb = 0.0
    if ca + cb > 0 and A:
        A[-1] += ca
        B[-1] += cb
    A, B = np.array(A), np.array(B)
    na, nb = A.sum(), B.sum()
    stat = np.sum((A / na - B / nb) ** 2 / ((A + B) / (na * nb)))
    return float(stat), int(len(A) - 1)


# ---------------------------------------------------------------- A
def test_A():
    sc = small_scene()
    out = []
    n_ref, n_c = 1500, 20000
    for label, kw in [("same-cell kernel", dict(kernel_radius=0.0)),
                      ("radius-1 kernel", dict(kernel_radius=1.0)),
                      ("per-cell rule", dict(mode="percell")),
                      ("same-cell, gmax=0", dict(kernel_radius=0.0, gmax=0))]:
        rng = np.random.default_rng(100)
        t0 = time.time()
        ref = [simulate_one(sc, BETA, GAMMA, rng, t_max=500.0, **kw)
               for _ in range(n_ref)]
        t_ref = time.time() - t0
        fs_ref = np.array([r["final_size"] for r in ref])
        off_ref = np.array([r["index_offspring"] for r in ref])
        t0 = time.time()
        res = simulate(sc, BETA, GAMMA, n_rep=n_c, seed=7, t_max=500.0, **kw)
        t_c = time.time() - t0
        chi, dof = chi2_hist(fs_ref, res.final_size, sc.N + 1)
        row = dict(case=label,
                   ref_final_mean=float(fs_ref.mean()),
                   c_final_mean=float(res.final_size.mean()),
                   z_final=two_sample(fs_ref, res.final_size),
                   ref_offspring_mean=float(off_ref.mean()),
                   c_offspring_mean=float(res.index_offspring.mean()),
                   z_offspring=two_sample(off_ref, res.index_offspring),
                   chi2=chi, dof=dof, n_ref=n_ref, n_c=n_c,
                   sec_per_run_python=t_ref / n_ref, sec_per_run_c=t_c / n_c)
        print("A", row)
        out.append(row)
    results["A_c_vs_python"] = out


# ---------------------------------------------------------------- B
def exact_final_size(N, lam, gamma):
    """Final-size distribution of Markovian SIR with pair hazard lam, one
    initial infective, by dynamic programming over (s, i)."""
    prob = np.zeros((N + 1, N + 2))
    prob[N - 1, 1] = 1.0
    out = np.zeros(N + 1)
    for s in range(N - 1, -1, -1):
        for i in range(N - s, 0, -1):
            p = prob[s, i]
            if p == 0.0:
                continue
            pinf = lam * s / (lam * s + gamma)
            if s > 0:
                prob[s - 1, i + 1] += p * pinf
            if i == 1:
                out[N - s] += p * (1 - pinf)
            else:
                prob[s, i - 1] += p * (1 - pinf)
    return out


def test_B():
    N = 40
    sc = uniform_room(1, 1, N, D=1.0, q=1.0)
    res = simulate(sc, BETA, GAMMA, n_rep=40000, seed=11, t_max=2000.0)
    lam = BETA / (N - 1)
    pex = exact_final_size(N, lam, GAMMA)
    emp = np.bincount(res.final_size, minlength=N + 1) / len(res.final_size)
    n = len(res.final_size)
    # chi-square against the exact law (merge bins with expected < 5)
    E, O, ce, co = [], [], 0.0, 0.0
    for k in range(1, N + 1):
        ce += pex[k] * n
        co += emp[k] * n
        if ce >= 5:
            E.append(ce)
            O.append(co)
            ce = co = 0.0
    if ce > 0:
        E[-1] += ce
        O[-1] += co
    E, O = np.array(E), np.array(O)
    chi = float(np.sum((O - E) ** 2 / E))
    row = dict(N=N, n_rep=n, exact_mean=float(np.arange(N + 1) @ pex),
               sim_mean=float(res.final_size.mean()),
               sim_se=float(res.final_size.std(ddof=1) / np.sqrt(n)),
               exact_p_minor_le3=float(pex[:4].sum()),
               sim_p_minor_le3=float(emp[:4].sum()),
               chi2=chi, dof=int(len(E) - 1))
    print("B", row)
    results["B_wellmixed_exact"] = row
    np.savez_compressed(os.path.join(HERE, "..", "data", "01_B_finalsize.npz"),
                        exact=pex, empirical=emp)


# ---------------------------------------------------------------- C
def test_C():
    out = []
    for name in ("office", "classroom"):
        sc = load_scene(name)
        rad = 1.0 / sc.a_m
        P = th.contact_kernel(sc, rad if rad >= 1 else 0.0)
        R0, v, u = th.R0_dense(sc, BETA, GAMMA, P=P, return_vectors=True)
        col = th.offspring_mean_by_birth_cell(sc, BETA, GAMMA, P=P)
        n = 40000
        # index drawn from the Perron vector -> E[exposure] = R0
        res = simulate(sc, BETA, GAMMA, n_rep=n, seed=21, P=P, gmax=0,
                       index_site=v, t_max=400.0)
        m, lo, hi = mean_ci(res.index_exposure)
        row = dict(scene=name, kernel_cells=float(P.getnnz() / sc.M),
                   R0_ngm=R0, exposure_perron_mean=m, lo=lo, hi=hi, n_rep=n)
        # uniform index -> E[exposure] = beta <q> / gamma
        res = simulate(sc, BETA, GAMMA, n_rep=n, seed=22, P=P, gmax=0,
                       t_max=400.0)
        m, lo, hi = mean_ci(res.index_exposure)
        row.update(wellmixed=BETA * sc.q.mean() / GAMMA,
                   exposure_uniform_mean=m, lo_u=lo, hi_u=hi)
        # fixed cells: the hottest and the coldest column sum
        for tag, x in (("max", int(np.argmax(col))), ("min", int(np.argmin(col)))):
            res = simulate(sc, BETA, GAMMA, n_rep=n // 4, seed=23, P=P, gmax=0,
                           index_site=x, t_max=400.0)
            m, lo, hi = mean_ci(res.index_exposure)
            row[f"cell_{tag}_theory"] = float(col[x])
            row[f"cell_{tag}_mean"] = m
            row[f"cell_{tag}_lo"] = lo
            row[f"cell_{tag}_hi"] = hi
        print("C", row)
        out.append(row)
    results["C_exposure_identity"] = out


# ---------------------------------------------------------------- D
def test_D():
    out = []
    n = 20000
    for name, D0 in [("office", 0.1), ("office", 1.0), ("office", 10.0),
                     ("office", 100.0), ("supermarket", 1.0),
                     ("classroom", 1.0), ("metro", 1.0)]:
        sc = load_scene(name, D0=D0)
        for ell in ((0.0, 1.0) if sc.a_m <= 1.0 else (0.0,)):
            P = th.contact_kernel(sc, ell / sc.a_m)
            p = th.pair_infection_probability(sc, BETA, GAMMA, P=P)
            R1 = (sc.N - 1) * p.mean()
            res = simulate(sc, BETA, GAMMA, n_rep=n, seed=31, P=P, gmax=0,
                           t_max=600.0)
            m, lo, hi = mean_ci(res.index_offspring)
            row = dict(scene=name, D0=D0, contact_range_m=ell,
                       kernel_cells=float(P.getnnz() / sc.M),
                       pair_theory=float(R1), sim_mean=m, lo=lo, hi=hi,
                       n_rep=n, meanfield=BETA * sc.q.mean() / GAMMA)
            print("D", row)
            out.append(row)
    results["D_pair_theory"] = out


# ---------------------------------------------------------------- E
def test_E():
    out = []
    n = 20000
    cases = [("office (scene, N=100)", load_scene("office")),
             ("uniform 20x20, q=1, D=1, N=100", uniform_room(20, 20, 100))]
    for label, sc in cases:
        res = simulate(sc, BETA, GAMMA, n_rep=n, seed=41, mode="percell",
                       t_max=400.0)
        m, lo, hi = mean_ci(res.index_offspring)
        k = int(res.major(0.1).sum())
        row = dict(case=label, index_offspring_mean=m, lo=lo, hi=hi,
                   n_major=k, n_rep=n,
                   mean_final_size=float(res.final_size.mean()),
                   max_final_size=int(res.final_size.max()))
        print("E", row)
        out.append(row)
    results["E_percell_rule"] = out


# ---------------------------------------------------------------- F
def test_F():
    out = []
    curves = {}
    t = np.arange(0, 80.25, 0.25)
    curves["t"] = t
    for N, D0 in [(100, 1.0), (1000, 10.0), (5000, 100.0)]:
        sc = load_scene("office", N=N, D0=D0)
        I0 = max(1, N // 100)                 # 1 % initial infectives
        n_rep = 200 if N <= 1000 else 40
        res = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=51, N=N,
                       n_index=I0, t_max=80.0, dt_out=0.25)
        ode = th.meanfield_ode(sc, BETA, GAMMA, t, I0=I0, N=N)
        Im = res.I.mean(axis=0) / N
        curves[f"N{N}_sim_mean"] = Im
        curves[f"N{N}_sim_q05"] = np.quantile(res.I / N, 0.05, axis=0)
        curves[f"N{N}_sim_q95"] = np.quantile(res.I / N, 0.95, axis=0)
        curves[f"N{N}_ode"] = ode["I"] / N
        row = dict(N=N, D0=D0, I0=I0, n_rep=n_rep,
                   occupancy=N / sc.M,
                   ode_peak=float(ode["I"].max() / N),
                   ode_peak_time=float(t[ode["I"].argmax()]),
                   sim_peak_of_mean=float(Im.max()),
                   sim_peak_time_of_mean=float(t[Im.argmax()]),
                   ode_final=float(1 - ode["S"][-1] / N),
                   sim_final_mean=float(1 - res.S[:, -1].mean() / N),
                   max_abs_diff_curve=float(np.max(np.abs(Im - ode["I"] / N))),
                   R0_ngm=th.R0_dense(sc, BETA, GAMMA))
        print("F", row)
        out.append(row)
    results["F_meanfield_limit"] = out
    np.savez_compressed(os.path.join(HERE, "..", "data", "01_F_curves.npz"),
                        **curves)


# ---------------------------------------------------------------- G
def test_G():
    """Piecewise-constant schedules.
    (i)  three identical segments per day must reproduce the unscheduled
         model (same law, different random streams);
    (ii) on/off schedule (in the room for the first third of each day): the
         index exposure must equal beta<q> times the expected in-room part of
         an Exp(gamma) infectious period that starts at phase 0."""
    sc = load_scene("office")
    n = 20000
    base = simulate(sc, BETA, GAMMA, n_rep=n, seed=61, t_max=400.0, dt_out=1.0)
    sched = dict(period=1.0, segments=[(0.2, sc.D, 1.0), (0.7, sc.D, 1.0),
                                       (1.0, sc.D, 1.0)])
    seg = simulate(sc, BETA, GAMMA, n_rep=n, seed=62, t_max=400.0, dt_out=1.0,
                   schedule=sched)
    row = dict(n_rep=n,
               final_mean_noschedule=float(base.final_size.mean()),
               final_mean_schedule=float(seg.final_size.mean()),
               z_final=two_sample(base.final_size, seg.final_size),
               offspring_noschedule=float(base.index_offspring.mean()),
               offspring_schedule=float(seg.index_offspring.mean()),
               z_offspring=two_sample(base.index_offspring, seg.index_offspring),
               p_major_noschedule=float(base.major().mean()),
               p_major_schedule=float(seg.major().mean()))
    w = 1.0 / 3.0
    onoff = dict(period=1.0, segments=[(w, sc.D, 1.0), (1.0, np.zeros(sc.M), 0.0)])
    r = simulate(sc, BETA, GAMMA, n_rep=40000, seed=63, gmax=0, t_max=600.0,
                 dt_out=100.0, schedule=onoff)
    inroom = (1 - np.exp(-GAMMA * w)) / GAMMA / (1 - np.exp(-GAMMA))
    m, lo, hi = mean_ci(r.index_exposure)
    row.update(onoff_exposure_theory=float(BETA * sc.q.mean() * inroom),
               onoff_exposure_sim=m, onoff_lo=lo, onoff_hi=hi,
               inroom_days_theory=float(inroom))
    print("G", row)
    results["G_schedules"] = row


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "ABCDEFG"
    for k in which:
        t0 = time.time()
        globals()["test_" + k]()
        print(f"--- test {k} done in {time.time() - t0:.1f} s", flush=True)
        merged = json.load(open(OUT)) if os.path.exists(OUT) else {}
        merged.update(results)
        json.dump(merged, open(OUT, "w"), indent=1, default=float)
