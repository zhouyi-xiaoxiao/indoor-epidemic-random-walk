#!/usr/bin/env python
"""Occupancy dependence (office scene, D0 = 1 and 10 lattice^2/day).

Scene differences and the effect of capacity limits are often attributed to
crowd density, but a force of infection beta q I/N is frequency-dependent,
so in the mean-field model density cancels out of R0.  This experiment shows
what occupancy does under two explicit assumptions:

  FD  frequency-dependent calibration: pair hazard beta q / rho_bar(N).
      Mean-field R0 is independent of N; discrete people still feel N through
      pair saturation (few distinct contacts at low occupancy).
  DD  density-dependent: pair hazard fixed at its N = 100 value.  Mean-field
      R0 is proportional to (N - 1).

Outputs data/07_occupancy.json.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import load_scene                  # noqa: E402
from bsc_sim.sim import simulate, summarize            # noqa: E402

BETA, GAMMA = 0.5, 0.14
N_REF = 100
NS = [20, 35, 50, 75, 100, 150, 200, 300, 400]
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "07_occupancy.json")


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    for D0 in (1.0, 10.0):
        for N in NS:
            for cal in ("FD", "DD"):
                key = f"D0={D0:g}|N={N}|{cal}"
                if key in out:
                    continue
                t0 = time.time()
                sc = load_scene("office", N=N, D0=D0)
                rho_ref = None if cal == "FD" else (N_REF - 1) / sc.M
                scale = 1.0 if cal == "FD" else (N - 1) / (N_REF - 1)
                R0 = th.R0_dense(sc, BETA, GAMMA) * scale
                p = th.pair_infection_probability(sc, BETA, GAMMA, N=N,
                                                  rho_ref=rho_ref)
                R1 = float((N - 1) * p.mean())
                n_rep = 2000
                res = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=700, N=N,
                               rho_ref=rho_ref, t_max=500.0, dt_out=0.5)
                s = summarize(res)
                r0 = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=701, N=N,
                              rho_ref=rho_ref, gmax=0, t_max=500.0, dt_out=5.0)
                off = r0.index_offspring
                # mean-field branching and ODE predictions
                qeff = sc.q * scale
                ext = th.extinction_probability(sc, BETA, GAMMA, q=qeff)
                t = np.arange(0, 500.5, 0.5)
                ode = th.meanfield_ode(sc, BETA, GAMMA, t, I0=1.0, N=N,
                                       rho_ref=rho_ref)
                row = dict(D0=D0, N=N, calibration=cal, occupancy=N / sc.M,
                           density_per_m2=N / sc.area_m2, R0_ngm=R0,
                           R1_pair_uniform=R1,
                           R_index_alone=float(off.mean()),
                           R_index_alone_se=float(off.std(ddof=1) / np.sqrt(n_rep)),
                           p_major_branching=float(1 - ext.mean()),
                           ode_attack=float(1 - ode["S"][-1] / N),
                           ode_peak_prev=float(ode["I"].max() / N),
                           ode_peak_time=float(t[ode["I"].argmax()]),
                           seconds=time.time() - t0)
                row.update(s)
                out[key] = row
                json.dump(out, open(PATH, "w"), indent=1, default=float)
                print(key, f"R0={R0:.2f} R1={R1:.2f} Rsim={row['R_index_alone']:.2f} "
                      f"Pmaj={s['p_major']:.3f} AR|maj={s.get('attack_major_mean', float('nan')):.3f} "
                      f"peak={s.get('peak_prev_major_mean', float('nan')):.3f} "
                      f"tpk={s.get('peak_time_major_mean', float('nan')):.1f} "
                      f"[{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
