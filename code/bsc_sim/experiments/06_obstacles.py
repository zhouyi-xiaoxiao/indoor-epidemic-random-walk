#!/usr/bin/env python
"""Obstacle experiment with IDENTICAL initial conditions across arms.

Room: 20 x 13 lattice, N = 100, uniform D = D0 = 1 lattice^2/day and uniform
q = 1.2 (R0 = beta q / gamma = 4.29).  Random barrier
cells at density rho_B in {0, 0.05, ..., 0.30}.

Design
* Nested barriers: for one layout the barrier set at a lower density is a
  subset of the set at a higher density.
* The 100 people start on cells that are walkable in EVERY arm (the walkable
  set of the densest arm), the index case is the same person at the same
  cell (the room centre), and each replicate uses the same
  random-number seed in every arm.  So all arms of a replicate start from the
  same microscopic state; they differ only by the barriers.
* 8 independent barrier layouts x 500 replicates per arm.
* Two calibrations of the pair hazard: frequency-dependent (FD, the default
  calibration: R0 does not depend on occupancy) and density-dependent (DD: the pair
  hazard is fixed at its barrier-free value, so removing floor area raises the
  crowding and R0).

Outputs data/06_obstacles.json, data/06_obstacles_runs.npz (per-run
outcomes), data/06_obstacle_curves.npz, data/06_snapshots.npz.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import keep_largest_component, uniform_room  # noqa: E402
from bsc_sim.sim import simulate, summarize            # noqa: E402

BETA, GAMMA = 0.5, 0.14
NX, NY, N, Q, D = 20, 13, 100, 1.2, 1.0
CENTRE = (10, 6)
FRACS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
N_LAYOUT, N_REP = 8, 500
T_MAX, DT_OUT = 400.0, 0.5
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "06_obstacles.json")


def nested_masks(seed):
    """Nested barrier masks that never cover the centre cell; pockets cut off
    from the main walkable component are closed too."""
    rng = np.random.default_rng(seed)
    cells = [i for i in range(NX * NY) if i != CENTRE[1] * NX + CENTRE[0]]
    perm = rng.permutation(cells)
    masks = {}
    for f in FRACS:
        k = int(round(f * NX * NY))
        b = np.zeros(NX * NY, dtype=bool)
        b[perm[:k]] = True
        acc = keep_largest_component(~b.reshape(NY, NX))
        masks[f] = ~acc
    return masks


def valid_masks(lay):
    """First nested family (seeds 6000+lay, +100, ...) whose densest arm keeps
    the centre cell inside the main walkable component."""
    for k in range(50):
        masks = nested_masks(6000 + lay + 100 * k)
        if (~masks[FRACS[-1]])[CENTRE[1], CENTRE[0]]:
            return masks
    raise RuntimeError("no valid layout")


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    runs = {}
    curves = {}
    rho0 = (N - 1) / (NX * NY)              # barrier-free occupancy (DD reference)
    for lay in range(N_LAYOUT):
        masks = valid_masks(lay)
        dense_acc = ~masks[FRACS[-1]]
        for f in FRACS:                      # nestedness check
            assert not (dense_acc & masks[f]).any(), "walkable sets not nested"
        rng = np.random.default_rng(6100 + lay)
        ys, xs = np.nonzero(dense_acc)
        pick = rng.integers(0, len(xs), size=(N_REP, N))
        start_xy = np.stack([xs[pick], ys[pick]], axis=-1)      # (rep, N, 2)
        start_xy[:, 0, 0], start_xy[:, 0, 1] = CENTRE           # index case
        for f in FRACS:
            for cal in ("FD", "DD"):
                key = f"layout{lay}|rhoB={f:.2f}|{cal}"
                sc = uniform_room(NX, NY, N, D=D, q=Q, barriers=masks[f])
                init = sc.site_index[start_xy[..., 1], start_xy[..., 0]]
                assert (init >= 0).all()
                t0 = time.time()
                res = simulate(sc, BETA, GAMMA, n_rep=N_REP, seed=6200 + lay,
                               init_pos=init, t_max=T_MAX, dt_out=DT_OUT,
                               rho_ref=None if cal == "FD" else rho0)
                s = summarize(res)
                scale = 1.0 if cal == "FD" else ((N - 1) / sc.M) / rho0
                R0 = th.R0_dense(sc, BETA, GAMMA) * scale
                p = th.pair_infection_probability(
                    sc, BETA, GAMMA, rho_ref=None if cal == "FD" else rho0)
                c_idx = int(sc.site_index[CENTRE[1], CENTRE[0]])
                r0 = simulate(sc, BETA, GAMMA, n_rep=N_REP, seed=6300 + lay,
                              init_pos=init, gmax=0, t_max=T_MAX, dt_out=5.0,
                              rho_ref=None if cal == "FD" else rho0)
                row = dict(layout=lay, rhoB_nominal=f,
                           rhoB_realised=float(1 - sc.M / (NX * NY)),
                           cells=sc.M, calibration=cal, R0_ngm=R0,
                           R1_pair_centre=float((N - 1) * p[c_idx].mean()),
                           R1_pair_uniform=float((N - 1) * p.mean()),
                           R_index_alone=float(r0.index_offspring.mean()),
                           seconds=time.time() - t0)
                row.update(s)
                out[key] = row
                runs[key + "|attack"] = res.attack_rate.astype(np.float32)
                runs[key + "|peak"] = (res.peak_I / N).astype(np.float32)
                runs[key + "|peak_time"] = res.peak_time.astype(np.float32)
                runs[key + "|R_alone"] = r0.index_offspring.astype(np.int16)
                if f in (0.0, 0.10) and cal == "FD":
                    curves[key + "|I"] = res.I.astype(np.int16)
                    curves["t"] = res.t
                print(key, f"cells={sc.M} R0={R0:.2f} R1c={row['R1_pair_centre']:.2f} "
                      f"Rsim={row['R_index_alone']:.2f} Pmaj={s['p_major']:.3f} "
                      f"AR={s['attack_all_mean']:.3f}", flush=True)
        json.dump(out, open(PATH, "w"), indent=1, default=float)
    np.savez_compressed(os.path.join(DATA, "06_obstacles_runs.npz"), **runs)
    np.savez_compressed(os.path.join(DATA, "06_obstacle_curves.npz"), **curves)

    # ---- snapshots for the spatial figure: layout 0, rho_B = 0 and 0.10, FD
    masks = valid_masks(0)
    dense_acc = ~masks[FRACS[-1]]
    rng = np.random.default_rng(6100)
    ys, xs = np.nonzero(dense_acc)
    n_snap_rep = 4000
    pick = rng.integers(0, len(xs), size=(n_snap_rep, N))
    start_xy = np.stack([xs[pick], ys[pick]], axis=-1)
    start_xy[:, 0, 0], start_xy[:, 0, 1] = CENTRE
    snap = {}
    for f in (0.0, 0.10):
        sc = uniform_room(NX, NY, N, D=D, q=Q, barriers=masks[f])
        init = sc.site_index[start_xy[..., 1], start_xy[..., 0]]
        res = simulate(sc, BETA, GAMMA, n_rep=n_snap_rep, seed=6400,
                       init_pos=init, t_max=60.0, dt_out=1.0,
                       snap_times=(5, 10, 15, 20))
        dens = res.snap.mean(axis=0)                     # (4, M) mean infectives per cell
        grid = np.stack([sc.to_grid(d) for d in dens])
        snap[f"rhoB{f:.2f}_density"] = grid
        snap[f"rhoB{f:.2f}_barriers"] = masks[f]
        # one single realisation as well (replicate 0)
        snap[f"rhoB{f:.2f}_single"] = np.stack(
            [sc.to_grid(d) for d in res.snap[0].astype(float)])
        # cumulative infection sites of replicate-averaged first 20 days
        snap[f"rhoB{f:.2f}_meanI"] = res.I.mean(axis=0)
    snap["times"] = np.array([5, 10, 15, 20])
    np.savez_compressed(os.path.join(DATA, "06_snapshots.npz"), **snap)


if __name__ == "__main__":
    main()
