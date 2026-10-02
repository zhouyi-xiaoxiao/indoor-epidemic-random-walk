#!/usr/bin/env python
"""s6_scenes_ode_pointseed.py -- mean-field lattice ODE started from ONE infective in ONE cell.

The scene table of the simulation part (data/bsc_sim/03_scenes.json) quotes the mean-field
ODE started from one infective spread uniformly over all cells.  A simulated outbreak starts in one
cell, so part of the difference between the simulated and the ODE peak is the initial condition and
not the mean-field approximation.  This script integrates the same ODE (code/bsc_sim/src/bsc_sim/
theory.py, meanfield_ode; eq. s2_model:eq:meanfield of the article) with the infective placed in a
single cell, for a regular subsample of seed cells (every STEP-th walkable cell), and records the
mean, minimum and maximum over seed cells of the peak prevalence, the peak day and the attack rate.

Scenes: office, supermarket, classroom, metro; D0 = 1, 10, 100 cell^2/day; 1 m contact kernel;
beta = 0.5/day, gamma = 0.14/day; output grid 0.5 day up to day 400 (as in 03_scenes.py).

Output: data/s6_scenes_ode_pointseed.json
Run:    python code/s6_scenes_ode_pointseed.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                                   # repository root
sys.path.insert(0, str(ROOT / "code" / "bsc_sim" / "src"))
from bsc_sim import theory as th                         # noqa: E402
from bsc_sim.scenes import SCENE_NAMES, load_scene       # noqa: E402

BETA, GAMMA, ELL_C = 0.5, 0.14, 1.0
T_MAX, DT_OUT = 400.0, 0.5
STEP = 3                                                 # every third walkable cell is used as a seed
OUT = HERE.parent / "data" / "s6_scenes_ode_pointseed.json"


def main():
    out = json.load(open(OUT)) if OUT.exists() else {}
    t = np.arange(0.0, T_MAX + DT_OUT, DT_OUT)
    for D0 in (1.0, 10.0, 100.0):
        for name in SCENE_NAMES:
            key = f"{name}|D0={D0:g}|range1m"
            if key in out:
                continue
            t0 = time.time()
            sc = load_scene(name, D0=D0)
            P = th.contact_kernel(sc, ELL_C / sc.a_m)
            uni = th.meanfield_ode(sc, BETA, GAMMA, t, I0=1.0, P=P)
            seeds = list(range(0, sc.M, STEP))
            peak, day, att = [], [], []
            for s in seeds:
                o = th.meanfield_ode(sc, BETA, GAMMA, t, I0=1.0, I0_site=s, P=P)
                peak.append(float(o["I"].max() / sc.N))
                day.append(float(t[o["I"].argmax()]))
                att.append(float(1 - o["S"][-1] / sc.N))
            peak, day, att = map(np.array, (peak, day, att))
            out[key] = dict(
                scene=name, D0=D0, kernel="range1m", N=sc.N, cells=sc.M, n_seeds=len(seeds), step=STEP,
                uniform_seed=dict(peak=float(uni["I"].max() / sc.N), day=float(t[uni["I"].argmax()]),
                                  attack=float(1 - uni["S"][-1] / sc.N)),
                point_seed=dict(peak_mean=float(peak.mean()), peak_min=float(peak.min()),
                                peak_max=float(peak.max()), day_mean=float(day.mean()),
                                day_min=float(day.min()), day_max=float(day.max()),
                                attack_mean=float(att.mean()), attack_min=float(att.min()),
                                attack_max=float(att.max())),
                seconds=time.time() - t0)
            json.dump(out, open(OUT, "w"), indent=1)
            print(key, out[key]["uniform_seed"], out[key]["point_seed"], f"[{time.time() - t0:.0f} s]",
                  flush=True)


if __name__ == "__main__":
    main()
