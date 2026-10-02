#!/usr/bin/env python
"""s7_interventions_barriers.py -- barrier experiment with identical starts at D = 1 and D = 10.

Purpose (Section 7, Table s7_interventions:tab:barriers).  The simulation part ran the barrier scan only at
D = 1 cell^2/day (code/bsc_sim/experiments/06_obstacles.py: 8 layouts x 500 runs) and the re-check
re-ran it with 40 layouts x 1000 runs, again at D = 1 (code/bsc_sim/verify/v08b_obstacles_many.py).
The article reports every simulation result at the default mobility AND at a higher one, so this script
repeats the re-check's design with the released simulator (code/bsc_sim/src/bsc_sim) at D = 1 and D = 10:

* room 20 x 13 cells, N = 100 people, uniform q = 1.2 (beta q / gamma = 4.29), uniform mobility D, same-cell
  contact kernel;
* nested random barrier cells at nominal densities 0, 0.10, 0.30 (pockets cut off from the main walkable
  component are closed, so the realised density at 0.30 is larger); the centre cell is never a barrier;
* identical starts: in every arm of a layout the people start on the same cells (drawn uniformly from the
  walkable set of the densest arm), the index case is the same person on the centre cell, and each
  replicate uses the same random seed;
* two calibrations of the pair hazard: FD (frequency dependent, hazard beta q / rho_bar with rho_bar
  recomputed for the remaining floor) and DD (density dependent, hazard frozen at its barrier-free value);
* statistics at the LAYOUT level: every quantity is first averaged over the replicates of one layout and
  the standard error is taken across layouts (the paired interval of this part treated runs as independent).

Usage:  python s7_interventions_barriers.py 1      (D = 1,  40 layouts x 1000 replicates)
        python s7_interventions_barriers.py 10     (D = 10, 40 layouts x  500 replicates)
Output: ../data/s7_interventions_barriers.json  (key "D=1" / "D=10"; merged across invocations)
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
from bsc_sim.scenes import keep_largest_component, uniform_room        # noqa: E402
from bsc_sim.sim import simulate                                       # noqa: E402

BETA, GAMMA = 0.5, 0.14
NX, NY, N, Q = 20, 13, 100, 1.2
CENTRE = (10, 6)
FRACS = [0.0, 0.10, 0.30]
N_LAYOUT = 40
N_REP = {1.0: 1000, 10.0: 500}
DT_OUT = 0.5
T_MAX = {1.0: 600.0, 10.0: 300.0}
PATH = os.path.join(HERE, "..", "data", "s7_interventions_barriers.json")
SEED0 = 71000                                   # layout seeds 71000, 71001, ... (part: 6000+, re-check: 5001+)


def nested_access(seed):
    """Nested walkable masks for the densities in FRACS (centre cell never a barrier)."""
    rng = np.random.default_rng(seed)
    cells = [i for i in range(NX * NY) if i != CENTRE[1] * NX + CENTRE[0]]
    perm = rng.permutation(cells)
    acc = {}
    for f in FRACS:
        b = np.zeros(NX * NY, dtype=bool)
        b[perm[:int(round(f * NX * NY))]] = True
        acc[f] = keep_largest_component(~b.reshape(NY, NX))
    return acc, rng


def valid_layouts():
    """First N_LAYOUT seeds whose densest arm keeps the centre cell and whose walkable sets are nested."""
    lay, seed = 0, SEED0
    while lay < N_LAYOUT:
        acc, rng = nested_access(seed)
        seed += 1
        dense = acc[FRACS[-1]]
        if not dense[CENTRE[1], CENTRE[0]]:
            continue
        if any((dense & ~acc[f]).any() for f in FRACS):
            continue
        yield lay, seed - 1, acc, rng
        lay += 1


def main(D):
    n_rep, t_max = N_REP[D], T_MAX[D]
    rho0 = (N - 1) / (NX * NY)                      # barrier-free occupancy (DD reference)
    per = {(f, cal): {k: [] for k in ("attack", "p_major", "attack_major", "peak_major", "peak_day_major",
                                       "R0", "R1_centre", "R1_uniform", "realised", "unfinished")}
           for f in FRACS for cal in ("FD", "DD")}
    seeds = []
    t0 = time.time()
    for lay, seed, acc, rng in valid_layouts():
        seeds.append(seed)
        ys, xs = np.nonzero(acc[FRACS[-1]])
        pick = rng.integers(0, len(xs), size=(n_rep, N))
        sx, sy = xs[pick], ys[pick]
        sx[:, 0], sy[:, 0] = CENTRE
        for f in FRACS:
            sc = uniform_room(NX, NY, N, D=D, q=Q, barriers=~acc[f])
            init = sc.site_index[sy, sx]
            assert (init >= 0).all()
            c_idx = int(sc.site_index[CENTRE[1], CENTRE[0]])
            for cal in ("FD", "DD"):
                rho_ref = None if cal == "FD" else rho0
                res = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=72000 + lay, init_pos=init, t_max=t_max,
                               dt_out=DT_OUT, rho_ref=rho_ref)
                ar = res.attack_rate
                maj = ar >= 0.1
                scale = 1.0 if cal == "FD" else ((N - 1) / sc.M) / rho0
                p = th.pair_infection_probability(sc, BETA, GAMMA, rho_ref=rho_ref)
                d = per[(f, cal)]
                d["attack"].append(float(ar.mean()))
                d["p_major"].append(float(maj.mean()))
                d["attack_major"].append(float(ar[maj].mean()))
                d["peak_major"].append(float((res.peak_I[maj] / N).mean()))
                d["peak_day_major"].append(float(res.peak_time[maj].mean()))
                d["R0"].append(float(th.R0_dense(sc, BETA, GAMMA) * scale))
                d["R1_centre"].append(float((N - 1) * p[c_idx].mean()))
                d["R1_uniform"].append(float((N - 1) * p.mean()))
                d["realised"].append(float(1 - sc.M / (NX * NY)))
                d["unfinished"].append(int(np.isnan(res.t_ext).sum()))
        print(f"D={D:g} layout {lay} done [{time.time() - t0:.0f} s]", flush=True)

    def ms(v):
        v = np.asarray(v, float)
        return [float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))]

    out = dict(D=D, n_layout=N_LAYOUT, n_rep=n_rep, t_max=t_max, dt_out=DT_OUT, layout_seeds=seeds,
               design="nested barriers, identical starts and seeds across arms; mean and SE across layouts",
               arms={}, paired={}, per_layout={})
    for (f, cal), d in per.items():
        key = f"rhoB={f:.2f}|{cal}"
        out["arms"][key] = {k: ms(v) for k, v in d.items() if k != "unfinished"}
        out["arms"][key]["unfinished_total"] = int(np.sum(d["unfinished"]))
        out["per_layout"][key] = {k: [float(x) for x in d[k]]
                                  for k in ("attack", "p_major", "peak_major", "peak_day_major")}
        if f > 0:
            base = per[(0.0, cal)]
            out["paired"][key] = {k: ms(np.asarray(d[k]) - np.asarray(base[k]))
                                  for k in ("attack", "p_major", "attack_major", "peak_major",
                                            "peak_day_major", "R1_centre", "R1_uniform", "R0")}
    out["seconds"] = time.time() - t0
    allout = json.load(open(PATH)) if os.path.exists(PATH) else {}
    allout[f"D={D:g}"] = out
    json.dump(allout, open(PATH, "w"), indent=1)
    for k, v in out["arms"].items():
        print(k, {kk: (round(vv[0], 4), round(vv[1], 4)) for kk, vv in v.items() if isinstance(vv, list)})
    for k, v in out["paired"].items():
        print("diff", k, {kk: (round(vv[0], 4), round(vv[1], 4)) for kk, vv in v.items()})


if __name__ == "__main__":
    main(float(sys.argv[1]))
