#!/usr/bin/env python
"""Spatial heterogeneity in q at fixed mean <q> raises R0.

Two-zone office on the 20 x 13 lattice:
  workspace  70 % of the area, D = 0.3 D0, q = q_W
  corridor   30 % of the area, D = 2.0 D0, q = q_C   (central band, 6 columns)
The area mean <q> = 1.2 is held fixed while the contrast c varies:
  q_C = 1.2 (1 - c),  q_W = 1.2 (1 + c * 0.3/0.7);   c = 0.5833 is the
reference case (q_W = 1.5, q_C = 0.5).  The mobility field D(r) is the same for every c,
so the only thing that changes is how unevenly q is distributed.

For each c we compute the mean-field R0 (NGM), the pair-level R1 for N = 100
discrete people, and we count secondary cases in the simulation; we also run
full epidemics.  Outputs data/05_heterogeneity.json.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import Scene                       # noqa: E402
from bsc_sim.sim import mean_ci, simulate, summarize   # noqa: E402

BETA, GAMMA = 0.5, 0.14
QBAR = 1.2
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "05_heterogeneity.json")


def two_zone(c: float, D0: float = 1.0, N: int = 100,
             variant: str = "zoneD") -> Scene:
    """variant 'zoneD': D = 0.3 (workspace) / 2.0 (corridor);
    variant 'uniformD': the same q field but D = 0.81 D0 everywhere (the area
    mean), which removes the correlation between q and mobility."""
    nx, ny = 20, 13
    zone = np.full((ny, nx), "W")
    zone[:, 7:13] = "C"                       # 6 of 20 columns = 30 %
    access = np.ones((ny, nx), dtype=bool)
    qW = QBAR * (1 + c * 0.3 / 0.7)
    qC = QBAR * (1 - c)
    D = np.where(zone == "W", 0.3, 2.0) * D0
    if variant == "uniformD":
        D = np.full((ny, nx), 0.7 * 0.3 + 0.3 * 2.0) * D0
    q = np.where(zone == "W", qW, qC)
    return Scene(name=f"twozone_c{c:g}", nx=nx, ny=ny, a_m=1.5, zone=zone,
                 access=access, D_grid=D, q_grid=q, N=N, meta=dict(c=c, D0=D0))


def perron(K):
    w, V = np.linalg.eig(K)
    k = int(np.argmax(w.real))
    v = np.abs(V[:, k].real)
    return float(w[k].real), v / v.sum()


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    c_ref = 1 - 0.5 / QBAR
    cs = [0.0, 0.2, 0.4, c_ref, 0.8, 1.0]
    # (a) mean-field R0 on a fine grid of c and D0 (cheap)
    if "meanfield_grid" not in out:
        grid = {}
        for D0 in (0.1, 1.0, 10.0, 100.0):
            grid[str(D0)] = []
            for c in np.linspace(0, 1, 26):
                sc = two_zone(float(c), D0)
                grid[str(D0)].append(dict(
                    c=float(c), R0=th.R0_dense(sc, BETA, GAMMA),
                    qW=float(sc.q.max()), qC=float(sc.q.min())))
        out["meanfield_grid"] = grid
        json.dump(out, open(PATH, "w"), indent=1, default=float)
    # (b) simulation
    for D0 in (1.0, 10.0):
        for N in (100, 1000):
            for c in cs:
                key = f"D0={D0:g}|N={N}|c={c:.4f}"
                if key in out:
                    continue
                t0 = time.time()
                sc = two_zone(c, D0, N)
                R0, v, _ = th.R0_dense(sc, BETA, GAMMA, return_vectors=True)
                lo, hi = th.R0_bounds(sc, BETA, GAMMA)
                K1 = th.pair_ngm(sc, BETA, GAMMA)
                R1, v1 = perron(K1)
                n = 20000 if N == 100 else 4000
                # index born according to the pair-level Perron vector
                r = simulate(sc, BETA, GAMMA, n_rep=n, seed=500, gmax=0,
                             index_site=v1, t_max=600.0, dt_out=5.0)
                m, mlo, mhi = mean_ci(r.index_offspring)
                # uniformly placed index: expectation beta<q>/gamma in mean field
                ru = simulate(sc, BETA, GAMMA, n_rep=n, seed=501, gmax=0,
                              t_max=600.0, dt_out=5.0)
                mu, ulo, uhi = mean_ci(ru.index_offspring)
                # generation ratio with the first generation also infectious
                r2 = simulate(sc, BETA, GAMMA, n_rep=n, seed=502, gmax=1,
                              index_site=v1, t_max=900.0, dt_out=5.0)
                g = r2.generation_sizes(2)
                rng = np.random.default_rng(1)
                bs = [g[ii, 2].sum() / g[ii, 1].sum()
                      for ii in (rng.integers(0, n, n) for _ in range(300))]
                # full epidemics, uniform index
                nf = 2000 if N == 100 else 300
                full = simulate(sc, BETA, GAMMA, n_rep=nf, seed=503,
                                t_max=600.0, dt_out=0.5)
                s = summarize(full)
                out[key] = dict(
                    D0=D0, N=N, c=c, qW=float(sc.q.max()), qC=float(sc.q.min()),
                    R0_ngm=R0, wellmixed=lo, frozen=hi,
                    R1_pair_ngm=R1, R1_pair_uniform=float(K1.sum(0).mean()),
                    sim_R_perron=m, sim_R_perron_lo=mlo, sim_R_perron_hi=mhi,
                    sim_R_uniform=mu, sim_R_uniform_lo=ulo, sim_R_uniform_hi=uhi,
                    sim_gen2_over_gen1=float(g[:, 2].sum() / g[:, 1].sum()),
                    sim_gen2_lo=float(np.quantile(bs, 0.025)),
                    sim_gen2_hi=float(np.quantile(bs, 0.975)),
                    n_rep=n, n_full=nf, p_major=s["p_major"],
                    p_major_lo=s["p_major_lo"], p_major_hi=s["p_major_hi"],
                    attack_major=s.get("attack_major_mean"),
                    attack_major_lo=s.get("attack_major_lo"),
                    attack_major_hi=s.get("attack_major_hi"),
                    attack_all=s["attack_all_mean"],
                    peak_prev_major=s.get("peak_prev_major_mean"),
                    peak_time_major=s.get("peak_time_major_mean"),
                    seconds=time.time() - t0)
                json.dump(out, open(PATH, "w"), indent=1, default=float)
                print(key, f"R0={R0:.3f} R1={R1:.3f} sim={m:.3f} [{mlo:.3f},{mhi:.3f}] "
                      f"gen2/gen1={out[key]['sim_gen2_over_gen1']:.3f} "
                      f"Pmaj={s['p_major']:.3f} AR={s['attack_all_mean']:.3f} "
                      f"[{time.time() - t0:.0f} s]", flush=True)


def extra():
    """Controls that isolate the mechanism (N = 100, D0 = 1 and 10):
    (i) 'uniformD' - same q contrast, uniform mobility;
    (ii) negative contrast - the hot zone is the mobile corridor
         (q_C > q_W), still with <q> = 1.2.
    Written to data/05_heterogeneity_extra.json."""
    PATH = os.path.join(DATA, "05_heterogeneity_extra.json")
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    c_ref = 1 - 0.5 / QBAR
    jobs = [("uniformD", c) for c in (0.0, 0.2, 0.4, c_ref, 0.8, 1.0)] + \
           [("zoneD", c) for c in (-0.5, -1.0, -1.5, -2.0)]
    for D0 in (1.0, 10.0):
        for variant, c in jobs:
            key = f"extra|{variant}|D0={D0:g}|N=100|c={c:.4f}"
            if key in out:
                continue
            t0 = time.time()
            sc = two_zone(c, D0, 100, variant)
            R0 = th.R0_dense(sc, BETA, GAMMA)
            lo, hi = th.R0_bounds(sc, BETA, GAMMA)
            K1 = th.pair_ngm(sc, BETA, GAMMA)
            R1, v1 = perron(K1)
            n = 20000
            r = simulate(sc, BETA, GAMMA, n_rep=n, seed=510, gmax=0,
                         index_site=v1, t_max=600.0, dt_out=5.0)
            m, mlo, mhi = mean_ci(r.index_offspring)
            full = simulate(sc, BETA, GAMMA, n_rep=2000, seed=511, t_max=600.0,
                            dt_out=0.5)
            s = summarize(full)
            out[key] = dict(variant=variant, D0=D0, N=100, c=c,
                            qW=float(sc.q_grid[0, 0]), qC=float(sc.q_grid[0, 8]),
                            R0_ngm=R0, wellmixed=lo, frozen=hi, R1_pair_ngm=R1,
                            R1_pair_uniform=float(K1.sum(0).mean()),
                            sim_R_perron=m, sim_R_perron_lo=mlo,
                            sim_R_perron_hi=mhi, p_major=s["p_major"],
                            p_major_lo=s["p_major_lo"], p_major_hi=s["p_major_hi"],
                            attack_all=s["attack_all_mean"],
                            attack_major=s.get("attack_major_mean"),
                            seconds=time.time() - t0)
            json.dump(out, open(PATH, "w"), indent=1, default=float)
            print(key, f"R0={R0:.3f} R1={R1:.3f} sim={m:.3f} Pmaj={s['p_major']:.3f} "
                  f"AR={s['attack_all_mean']:.3f}", flush=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "main"
    if which == "main":
        main()
    else:
        extra()
