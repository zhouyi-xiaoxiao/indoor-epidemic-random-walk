#!/usr/bin/env python
"""Interventions, computed from the model.

Scenes: office and classroom, D0 = 1 and 10 lattice^2/day, contact range 1 m.
Arms
  baseline
  ventilation       q -> 0.5 q everywhere            (air changes x 2)
  hotspot vent.     q -> 0.5 q on the 25 % of cells with the largest
                    NGM elasticity (reproductive value x birth density)
  coldspot vent.    q -> 0.5 q on the 25 % of cells with the smallest one
  masks             beta -> 0.3 beta                 (assumed mask factor)
  capacity 50 %     N -> N/2, under FD and under DD calibration
  partitions        barrier cells added on 15 % of the lattice (random cells of
                    the main zone, floor kept connected), FD and DD
  combination       ventilation + capacity 50 % + partitions, FD and DD

For each arm: mean-field R0, pair-level R1, counted secondary cases, and full
stochastic epidemics.  Output data/08_interventions.json.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import is_connected, load_scene    # noqa: E402
from bsc_sim.sim import simulate, summarize            # noqa: E402

BETA, GAMMA = 0.5, 0.14
ELL_C = 1.0
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "08_interventions.json")
MAIN_ZONE = {"office": ("W",), "classroom": ("d", "D")}


def add_partitions(sc, frac=0.15, seed=8):
    """Close frac * (lattice cells) cells of the main zone, keeping the floor
    connected.  Returns a new Scene."""
    rng = np.random.default_rng(seed)
    target = int(round(frac * sc.nx * sc.ny))
    access = sc.access.copy()
    cand = [(y, x) for y in range(sc.ny) for x in range(sc.nx)
            if sc.zone[y, x] in MAIN_ZONE[sc.name]]
    order = rng.permutation(len(cand))
    added = 0
    for k in order:
        if added == target:
            break
        y, x = cand[k]
        access[y, x] = False
        if is_connected(access):
            added += 1
        else:
            access[y, x] = True
    zone = sc.zone.copy()
    zone[~access & sc.access] = "X"
    return sc.copy_with(access=access, zone=zone,
                        D_grid=np.where(access, sc.D_grid, 0.0),
                        q_grid=np.where(access, sc.q_grid, 0.0)), added


def evaluate(sc, beta, q, N, rho_ref, scale, n_rep, seed):
    P = th.contact_kernel(sc, ELL_C / sc.a_m)
    R0 = th.R0_dense(sc, beta, GAMMA, P=P, q=q) * scale
    p = th.pair_infection_probability(sc, beta, GAMMA, N=N, q=q, P=P,
                                      rho_ref=rho_ref)
    res = simulate(sc, beta, GAMMA, n_rep=n_rep, seed=seed, N=N, q=q, P=P,
                   rho_ref=rho_ref, t_max=500.0, dt_out=0.5)
    r0 = simulate(sc, beta, GAMMA, n_rep=n_rep, seed=seed + 1, N=N, q=q, P=P,
                  rho_ref=rho_ref, gmax=0, t_max=500.0, dt_out=5.0)
    row = dict(R0_ngm=R0, R1_pair_uniform=float((N - 1) * p.mean()),
               R_index_alone=float(r0.index_offspring.mean()),
               R_index_alone_se=float(r0.index_offspring.std(ddof=1) / np.sqrt(n_rep)),
               N=N, cells=sc.M)
    row.update(summarize(res))
    return row


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    n_rep = 2000
    for name in ("office", "classroom"):
        for D0 in (1.0, 10.0):
            base = load_scene(name, D0=D0)
            N0 = base.N
            part, n_added = add_partitions(base)
            P = th.contact_kernel(base, ELL_C / base.a_m)
            R0, v, u = th.R0_dense(base, BETA, GAMMA, P=P, return_vectors=True)
            elast = u * v / (u @ v)                 # NGM elasticity per cell
            order = np.argsort(-elast)
            k = int(round(0.25 * base.M))
            q_hot = base.q.copy()
            q_hot[order[:k]] *= 0.5
            q_cold = base.q.copy()
            q_cold[order[-k:]] *= 0.5
            # the same targeting rule, but with the pair-level matrix K1
            K1 = th.pair_ngm(base, BETA, GAMMA, P=P)
            w1, vl1, vr1 = sla.eig(K1, left=True, right=True)
            k1 = int(np.argmax(w1.real))
            el1 = np.abs(vl1[:, k1].real) * np.abs(vr1[:, k1].real)
            order1 = np.argsort(-el1)
            q_hot1 = base.q.copy()
            q_hot1[order1[:k]] *= 0.5
            rho_ref0 = (N0 - 1) / base.M            # DD reference (baseline room)
            arms = [
                ("baseline", base, BETA, base.q, N0, "FD"),
                ("ventilation q x 0.5", base, BETA, 0.5 * base.q, N0, "FD"),
                ("hotspot ventilation (25 % of cells)", base, BETA, q_hot, N0, "FD"),
                ("coldspot ventilation (25 % of cells)", base, BETA, q_cold, N0, "FD"),
                ("pair-level hotspot ventilation (25 % of cells)", base, BETA,
                 q_hot1, N0, "FD"),
                ("masks beta x 0.3", base, 0.3 * BETA, base.q, N0, "FD"),
                ("capacity 50 %", base, BETA, base.q, N0 // 2, "FD"),
                ("capacity 50 %", base, BETA, base.q, N0 // 2, "DD"),
                (f"partitions (+{n_added} barrier cells)", part, BETA, part.q, N0, "FD"),
                (f"partitions (+{n_added} barrier cells)", part, BETA, part.q, N0, "DD"),
                ("combination", part, BETA, 0.5 * part.q, N0 // 2, "FD"),
                ("combination", part, BETA, 0.5 * part.q, N0 // 2, "DD"),
            ]
            for label, sc, beta, q, N, cal in arms:
                key = f"{name}|D0={D0:g}|{label}|{cal}"
                if key in out:
                    continue
                t0 = time.time()
                rho_ref = None if cal == "FD" else rho_ref0
                scale = 1.0 if cal == "FD" else ((N - 1) / sc.M) / rho_ref0
                row = evaluate(sc, beta, q, N, rho_ref, scale, n_rep, seed=800)
                row.update(scene=name, D0=D0, arm=label, calibration=cal,
                           seconds=time.time() - t0)
                out[key] = row
                json.dump(out, open(PATH, "w"), indent=1, default=float)
                print(key, f"R0={row['R0_ngm']:.2f} R1={row['R1_pair_uniform']:.2f} "
                      f"Rsim={row['R_index_alone']:.2f} Pmaj={row['p_major']:.3f} "
                      f"AR={row['attack_all_mean']:.3f} "
                      f"peak|maj={row.get('peak_prev_major_mean', float('nan')):.3f}",
                      flush=True)


if __name__ == "__main__":
    main()
