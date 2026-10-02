#!/usr/bin/env python
"""Fair benchmark: spectral / linear-algebra evaluation vs Monte Carlo at the
SAME accuracy target, for the SAME quantity.

A comparison of one scalar formula with full epidemic simulations says
nothing about efficiency.  Here each row compares two ways of
computing one and the same number, and Monte Carlo is given its best estimator
and is charged only for the replicates it needs to reach the stated accuracy
(relative standard error 1 %, or absolute standard error 0.01 for a
probability).  All times are single-core CPU times (time.process_time),
medians of repeated measurements; the stochastic simulator is the C core.

Tasks (office scene, D0 = 1 and 10)
  T1  mean-field R0 = rho(K)
        spectral: sparse LU of (gamma I - L) + power iteration
        MC      : generation-wise particle power method (walkers carry weight
                  beta q integrated along their path), best case for MC
  T2  mean-field secondary cases of a case born at a given cell, (1^T K)_x
        spectral: one sparse solve;  MC: single-walker path integral
  T3  discrete-people secondary cases R1 (uniform index)
        theory  : one CG solve on the pair lattice
        MC      : individual-based runs with only the index transmitting
  T4  probability of a major outbreak
        theory  : spatial branching fixed point (mean-field)
        MC      : full individual-based epidemics
  T5  attack rate of a major outbreak
        theory  : mean-field lattice ODE;  MC: full epidemics

Output data/10_benchmark.json.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.refsim import simulate_one                # noqa: E402
from bsc_sim.scenes import load_scene                  # noqa: E402
from bsc_sim.sim import simulate, summarize            # noqa: E402

BETA, GAMMA = 0.5, 0.14
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "10_benchmark.json")


def cpu_median(f, repeat=7):
    ts, val = [], None
    for _ in range(repeat):
        t0 = time.process_time()
        val = f()
        ts.append(time.process_time() - t0)
    return float(np.median(ts)), val


def walker_paths(sc, starts, rng, ptr, idx, W, Wtot, cum):
    """Path integrals of beta*q along single random walks killed at rate gamma.
    Returns (total weight per walker, per-cell weight summed over walkers)."""
    n = len(starts)
    x = starts.copy()
    life = rng.exponential(1.0 / GAMMA, n)
    wsum = np.zeros(n)
    cell = np.zeros(sc.M)
    alive = np.arange(n)
    bq = BETA * sc.q
    while len(alive):
        xa = x[alive]
        hold = rng.exponential(1.0, len(alive)) / Wtot[xa]
        dt = np.minimum(hold, life[alive])
        w = bq[xa] * dt
        wsum[alive] += w
        np.add.at(cell, xa, w)
        life[alive] -= dt
        moving = hold < life[alive] + dt          # hop happened before death
        # choose the neighbour for the walkers that hop
        mv = alive[moving]
        xm = x[mv]
        u = rng.random(len(mv)) * Wtot[xm]
        # up to 4 neighbours: cumulative search
        k = (cum[xm] < u[:, None]).sum(axis=1)
        x[mv] = idx[ptr[xm] + np.minimum(k, (ptr[xm + 1] - ptr[xm]) - 1)]
        alive = mv
    return wsum, cell


def main():
    out = {}
    for D0 in (1.0, 10.0):
        sc = load_scene("office", D0=D0)
        L = th.generator(sc)
        ptr, idx, W = th.hop_rates(sc)
        Wtot = np.add.reduceat(W, ptr[:-1])
        cum = np.full((sc.M, 4), np.inf)
        for s in range(sc.M):
            c = np.cumsum(W[ptr[s]:ptr[s + 1]])
            cum[s, :len(c)] = c
        rng = np.random.default_rng(10)
        res = {}

        # ---------------- T1: R0 = rho(K)
        R_dense = th.R0_dense(sc, BETA, GAMMA)
        t_sparse, R_sparse = cpu_median(lambda: th.R0_sparse(sc, BETA, GAMMA, L=L))
        t_dense, _ = cpu_median(lambda: th.R0_dense(sc, BETA, GAMMA, L=L))
        # particle power method: n walkers per generation; the burn-in is set
        # from the subdominant eigenvalue so that the start-up bias is < 0.1 %
        ev = np.sort(np.abs(np.linalg.eigvals(th.ngm(sc, BETA, GAMMA, L=L))))[::-1]
        ratio = float(ev[1] / ev[0])
        n_w, G = 20000, 30
        burn = int(np.ceil(np.log(1e-3) / np.log(ratio)))
        t0 = time.process_time()
        dist = np.ones(sc.M) / sc.M
        ks = []
        for g in range(burn + G):
            starts = rng.choice(sc.M, size=n_w, p=dist)
            wsum, cell = walker_paths(sc, starts, rng, ptr, idx, W, Wtot, cum)
            if g >= burn:
                ks.append(wsum.mean())
            dist = cell / cell.sum()
        t_mc = time.process_time() - t0
        ks = np.array(ks)
        est = ks.mean()
        se = ks.std(ddof=1) / np.sqrt(len(ks))     # generations nearly independent
        rel = se / est
        # time to reach 1 % relative SE (scales as 1/n)
        t_mc_1pct = t_mc * (rel / 0.01) ** 2
        res["T1_R0_meanfield"] = dict(
            value_dense=R_dense, value_sparse=R_sparse,
            abs_err_sparse=abs(R_sparse - R_dense),
            cpu_sparse_s=t_sparse, cpu_dense_s=t_dense,
            mc_estimate=float(est), mc_rel_se=float(rel), mc_cpu_s=t_mc,
            subdominant_ratio=ratio, mc_burn_in_generations=burn,
            mc_walkers=n_w * (burn + G), mc_cpu_for_1pct_s=float(t_mc_1pct),
            mc_cpu_for_0p1pct_s=float(t_mc * (rel / 0.001) ** 2),
            speedup_at_1pct=float(t_mc_1pct / t_sparse),
            speedup_at_0p1pct=float(t_mc * (rel / 0.001) ** 2 / t_sparse))

        # ---------------- T2: (1^T K)_x for the hottest cell
        col = th.offspring_mean_by_birth_cell(sc, BETA, GAMMA, L=L)
        x = int(np.argmax(col))
        A = (GAMMA * sp.identity(sc.M, format="csc") - L.tocsc())

        def one_solve():
            return spla.spsolve(A.T.tocsc(), BETA * sc.q)[x]
        t_solve, val = cpu_median(one_solve)
        n_w = 200000
        t0 = time.process_time()
        wsum, _ = walker_paths(sc, np.full(n_w, x), rng, ptr, idx, W, Wtot, cum)
        t_mc = time.process_time() - t0
        rel = wsum.std(ddof=1) / np.sqrt(n_w) / wsum.mean()
        res["T2_column_sum"] = dict(
            cell=x, value_solve=float(val), cpu_solve_s=t_solve,
            mc_estimate=float(wsum.mean()), mc_rel_se=float(rel), mc_cpu_s=t_mc,
            mc_cv_single_walker=float(wsum.std(ddof=1) / wsum.mean()),
            mc_cpu_for_1pct_s=float(t_mc * (rel / 0.01) ** 2),
            speedup_at_1pct=float(t_mc * (rel / 0.01) ** 2 / t_solve))

        # ---------------- T3: discrete-people R1, uniform index
        t_pair, p = cpu_median(
            lambda: th.pair_infection_probability(sc, BETA, GAMMA, L=L), repeat=5)
        R1 = float((sc.N - 1) * p.mean())
        n = 40000
        t0 = time.process_time()
        r = simulate(sc, BETA, GAMMA, n_rep=n, seed=1010, gmax=0, t_max=600.0,
                     dt_out=50.0)
        t_mc = time.process_time() - t0
        off = r.index_offspring
        rel = off.std(ddof=1) / np.sqrt(n) / off.mean()
        res["T3_R1_discrete"] = dict(
            value_pair_solve=R1, cpu_pair_s=t_pair,
            mc_estimate=float(off.mean()), mc_rel_se=float(rel), mc_cpu_s=t_mc,
            mc_runs=n, mc_cv_single_run=float(off.std(ddof=1) / off.mean()),
            mc_runs_for_1pct=float(n * (rel / 0.01) ** 2),
            mc_cpu_for_1pct_s=float(t_mc * (rel / 0.01) ** 2),
            speedup_at_1pct=float(t_mc * (rel / 0.01) ** 2 / t_pair),
            theory_minus_mc_in_se=float((R1 - off.mean()) / (rel * off.mean())))

        # ---------------- T4/T5: outbreak probability and attack rate
        t_br, ext = cpu_median(
            lambda: th.extinction_probability(sc, BETA, GAMMA, L=L), repeat=5)
        p_br = float(1 - ext.mean())
        tt = np.arange(0, 400.5, 0.5)
        t_ode, ode = cpu_median(
            lambda: th.meanfield_ode(sc, BETA, GAMMA, tt, L=L), repeat=3)
        n = 8000
        t0 = time.process_time()
        full = simulate(sc, BETA, GAMMA, n_rep=n, seed=1011, t_max=400.0,
                        dt_out=0.5)
        t_mc = time.process_time() - t0
        s = summarize(full)
        pm = s["p_major"]
        se_p = np.sqrt(pm * (1 - pm) / n)
        ar = full.attack_rate[full.major()]
        rel_ar = ar.std(ddof=1) / np.sqrt(len(ar)) / ar.mean()
        res["T4_p_major"] = dict(
            value_branching_meanfield=p_br, cpu_branching_s=t_br,
            mc_estimate=pm, mc_se=float(se_p), mc_runs=n, mc_cpu_s=t_mc,
            mc_cpu_per_run_s=t_mc / n,
            mc_runs_for_se_0p01=float(pm * (1 - pm) / 0.01 ** 2),
            mc_cpu_for_se_0p01_s=float(pm * (1 - pm) / 0.01 ** 2 * t_mc / n),
            theory_bias=float(p_br - pm),
            speedup_at_se_0p01=float(pm * (1 - pm) / 0.01 ** 2 * t_mc / n / t_br))
        res["T5_attack_rate"] = dict(
            value_ode=float(1 - ode["S"][-1] / sc.N), cpu_ode_s=t_ode,
            mc_estimate=float(ar.mean()), mc_rel_se=float(rel_ar),
            mc_runs_for_1pct=float(n * (rel_ar / 0.01) ** 2),
            mc_cpu_for_1pct_s=float(t_mc * (rel_ar / 0.01) ** 2),
            theory_bias=float(1 - ode["S"][-1] / sc.N - ar.mean()),
            speedup_at_1pct=float(t_mc * (rel_ar / 0.01) ** 2 / t_ode))

        # ---------------- context: the readable pure-Python simulator
        rng2 = np.random.default_rng(3)
        t0 = time.process_time()
        nref = 3 if D0 == 1.0 else 1
        evs = [simulate_one(sc, BETA, GAMMA, rng2, t_max=400.0) for _ in range(nref)]
        n_ev = sum(e["n_events"] for e in evs)
        res["python_reference_cpu_per_event_s"] = (time.process_time() - t0) / max(n_ev, 1)
        res["python_reference_cpu_per_run_s"] = (
            res["python_reference_cpu_per_event_s"] * float(full.n_events.mean()))
        res["c_core_cpu_per_run_s"] = t_mc / n
        res["events_per_run"] = float(full.n_events.mean())
        out[f"D0={D0:g}"] = res
        json.dump(out, open(PATH, "w"), indent=1, default=float)
        for k, v in res.items():
            print(D0, k, v, flush=True)


if __name__ == "__main__":
    main()
