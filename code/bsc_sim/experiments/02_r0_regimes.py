#!/usr/bin/env python
"""Reproduction numbers: mean-field NGM vs pair-level theory vs counted
secondary cases in the individual-based simulation.

Part 1 (finite size, homogeneous room): L x L rooms with uniform q and D at
the office occupancy.
Part 2 (mobility, office scene): R as a function of D0 from 0.03 to 1000
lattice^2/day.

Writes data/02_finite_size.json and data/02_mobility.json (checkpointed).
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import load_scene, uniform_room    # noqa: E402
from bsc_sim.sim import mean_ci, simulate              # noqa: E402

BETA, GAMMA = 0.5, 0.14
DATA = os.path.join(HERE, "..", "data")


def load(path):
    return json.load(open(path)) if os.path.exists(path) else {}


def dump(path, obj):
    json.dump(obj, open(path, "w"), indent=1, default=float)


def perron(K):
    w, V = np.linalg.eig(K)
    k = int(np.argmax(w.real))
    v = np.abs(V[:, k].real)
    return float(w[k].real), v / v.sum()


# ------------------------------------------------------------------ part 1
def finite_size():
    path = os.path.join(DATA, "02_finite_size.json")
    out = load(path)
    occ = 99 / 260            # other people per cell, as in the office lattice
    for D in (1.0, 0.51, 10.0):
        for Lside in (4, 5, 7, 10, 14, 20, 28, 40):
            key = f"D{D}_L{Lside}"
            if key in out:
                continue
            N = int(round(occ * Lside * Lside)) + 1
            sc = uniform_room(Lside, Lside, N, D=D, q=1.0)
            t0 = time.time()
            Rcf, g0 = th.pair_R0_homogeneous(Lside, Lside, N, BETA, GAMMA, D)
            Rper, g0p = th.pair_R0_homogeneous(Lside, Lside, N, BETA, GAMMA, D,
                                               "periodic")
            Rex = np.nan
            if Lside <= 28:
                p = th.pair_infection_probability(sc, BETA, GAMMA)
                Rex = float((N - 1) * p.mean())
            n = 20000
            res = simulate(sc, BETA, GAMMA, n_rep=n, seed=200 + Lside, gmax=0,
                           t_max=500.0)
            m, lo, hi = mean_ci(res.index_offspring)
            me, loe, hie = mean_ci(res.index_exposure)
            out[key] = dict(
                D=D, L=Lside, N=N, n_rep=n,
                R0_meanfield=BETA / GAMMA,
                R0_ngm=th.R0_dense(sc, BETA, GAMMA),
                pair_closed_form=Rcf, g0=g0, pair_periodic=Rper,
                pair_exact=Rex, sim_offspring=m, sim_lo=lo, sim_hi=hi,
                sim_exposure=me, sim_exposure_lo=loe, sim_exposure_hi=hie,
                wellmixed_finiteN=(BETA / GAMMA) / (1 + BETA / (GAMMA * (N - 1))),
                seconds=time.time() - t0)
            print(key, out[key], flush=True)
            dump(path, out)
    # infinite-lattice limit of g0 at the same density (numerical k-integral)
    for D in (1.0, 0.51, 10.0):
        k = (np.arange(4000) + 0.5) * np.pi / 4000
        mu = 2 * D * (1 - np.cos(k))[:, None] + 2 * D * (1 - np.cos(k))[None, :]
        g_inf = float(np.mean(1.0 / (GAMMA + 2 * mu)))
        out[f"D{D}_Linf"] = dict(D=D, g0=g_inf,
                                 pair_closed_form=(BETA / GAMMA) /
                                 (1 + BETA / occ * g_inf))
    dump(path, out)


# ------------------------------------------------------------------ part 2
def mobility():
    path = os.path.join(DATA, "02_mobility.json")
    out = load(path)
    for D0 in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0):
        key = f"D0_{D0}"
        if key in out:
            continue
        sc = load_scene("office", D0=D0)
        t0 = time.time()
        R0, v, u = th.R0_dense(sc, BETA, GAMMA, return_vectors=True)
        lo_b, hi_b = th.R0_bounds(sc, BETA, GAMMA)
        K1 = th.pair_ngm(sc, BETA, GAMMA)
        R1, v1 = perron(K1)
        R1_uniform = float(K1.sum(axis=0).mean())
        n = 20000 if D0 <= 30 else (6000 if D0 <= 100 else 2000)
        # index case born according to the pair-level Perron vector
        res = simulate(sc, BETA, GAMMA, n_rep=n, seed=300, gmax=0,
                       index_site=v1, t_max=600.0)
        m, lo, hi = mean_ci(res.index_offspring)
        # uniform index case
        res_u = simulate(sc, BETA, GAMMA, n_rep=n, seed=301, gmax=0,
                         t_max=600.0)
        mu_, lou, hiu = mean_ci(res_u.index_offspring)
        # second generation: index + its direct cases transmit (gmax = 1)
        res2 = simulate(sc, BETA, GAMMA, n_rep=n, seed=302, gmax=1,
                        index_site=v1, t_max=900.0)
        g = res2.generation_sizes(2)
        ratio = g[:, 2].sum() / g[:, 1].sum()
        # bootstrap CI for the ratio of means
        rng = np.random.default_rng(5)
        bs = []
        for _ in range(400):
            ii = rng.integers(0, n, n)
            bs.append(g[ii, 2].sum() / g[ii, 1].sum())
        out[key] = dict(
            D0=D0, n_rep=n, R0_ngm=R0, lower_bound=lo_b, upper_bound=hi_b,
            R1_pair_ngm=R1, R1_pair_uniform_index=R1_uniform,
            sim_offspring_perron=m, sim_perron_lo=lo, sim_perron_hi=hi,
            sim_offspring_uniform=mu_, sim_uniform_lo=lou, sim_uniform_hi=hiu,
            sim_gen2_over_gen1=float(ratio),
            sim_gen2_lo=float(np.quantile(bs, 0.025)),
            sim_gen2_hi=float(np.quantile(bs, 0.975)),
            seconds=time.time() - t0)
        print(key, out[key], flush=True)
        dump(path, out)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "12"
    if "1" in which:
        finite_size()
    if "2" in which:
        mobility()
