#!/usr/bin/env python
"""Dwell time: what the closed-room reproduction numbers mean for scenes where
people stay minutes or hours.

Every closed-room R0 assumes that an infective spends its whole
infectious period (1/gamma = 7 days) in the room with the same closed group.
The scene files give dwell times of 8 h (office), 50-75 min (classroom), 27 min
(supermarket) and 20 min (metro).  Two consistent corrections:

(1) Recurrent closed cohort (office colleagues, a class): the same group is
    in the room for a fraction w of every day, and nothing happens in the
    room otherwise.  For a daily cycle much shorter than the infectious period
    this is equivalent to continuous occupancy with the recovery clock running
    1/w times faster in room time:
        R0(w) = rho( P^T diag(beta q) (gamma/w I - L)^{-1} ),
    and the same substitution gamma -> gamma/w in the pair-level formula.
    Checked against a simulation with an explicit daily on/off schedule.

(2) One visit (supermarket, metro; also one office day or one lecture): the
    expected number of people one infective infects during a single stay of
    length T, counted directly in the simulation (only the index transmits,
    run stopped at T) and compared with the mean-field value beta <q> T and
    the pair-level value.

Output data/11_dwell_time.json.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import SCENE_NAMES, load_scene     # noqa: E402
from bsc_sim.sim import mean_ci, simulate, summarize   # noqa: E402

BETA, GAMMA = 0.5, 0.14
ELL_C = 1.0
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "11_dwell_time.json")
# hours per day in the room (dwell times of the scene files)
HOURS = {"office": 8.0, "classroom": 1.0, "supermarket": 27 / 60,
         "metro": 2 * 20 / 60}
VISIT_H = {"office": 8.0, "classroom": 1.0, "supermarket": 27 / 60,
           "metro": 20 / 60}


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    for D0 in (1.0, 10.0):
        for name in SCENE_NAMES:
            key = f"{name}|D0={D0:g}"
            if key in out:
                continue
            t0 = time.time()
            sc = load_scene(name, D0=D0)
            P = th.contact_kernel(sc, ELL_C / sc.a_m)
            w = HOURS[name] / 24.0
            row = dict(scene=name, D0=D0, hours_per_day=HOURS[name], w=w)
            # ---- (1) duty-cycle reproduction numbers
            row["R0_continuous"] = th.R0_dense(sc, BETA, GAMMA, P=P)
            row["R0_duty_cycle"] = th.R0_dense(sc, BETA, GAMMA / w, P=P)
            row["wellmixed_duty_cycle"] = BETA * sc.q.mean() * w / GAMMA
            p = th.pair_infection_probability(sc, BETA, GAMMA / w, P=P)
            row["R1_pair_duty_cycle"] = float((sc.N - 1) * p.mean())
            if name in ("office", "classroom"):
                # explicit daily schedule: in the room for the first w of each
                # day (movement + transmission), frozen and no transmission
                # for the rest; recovery runs all the time.
                frozen = np.zeros(sc.M)
                sched = dict(period=1.0,
                             segments=[(w, sc.D, 1.0), (1.0, frozen, 0.0)])
                n = 20000
                r = simulate(sc, BETA, GAMMA, n_rep=n, seed=1100, P=P, gmax=0,
                             schedule=sched, t_max=400.0, dt_out=50.0)
                m, lo, hi = mean_ci(r.index_offspring)
                row.update(sim_R_schedule=m, sim_R_schedule_lo=lo,
                           sim_R_schedule_hi=hi, n_rep_R=n)
                full = simulate(sc, BETA, GAMMA, n_rep=2000, seed=1101, P=P,
                                schedule=sched, t_max=1500.0, dt_out=1.0)
                s = summarize(full)
                ext = th.extinction_probability(sc, BETA, GAMMA / w, P=P)
                row.update(p_major=s["p_major"], p_major_lo=s["p_major_lo"],
                           p_major_hi=s["p_major_hi"],
                           p_major_branching=float(1 - ext.mean()),
                           attack_all=s["attack_all_mean"],
                           attack_major=s.get("attack_major_mean"),
                           peak_prev_major=s.get("peak_prev_major_mean"),
                           peak_time_major=s.get("peak_time_major_mean"),
                           n_unfinished=s["n_unfinished"])
            # ---- (2) one visit
            T = VISIT_H[name] / 24.0
            n, chunk = 200000, 20000
            offs, expo = [], []
            for c in range(n // chunk):
                r = simulate(sc, BETA, GAMMA, n_rep=chunk, seed=1102 + c, P=P,
                             gmax=0, t_max=T, dt_out=T)
                offs.append(r.index_offspring)
                expo.append(r.index_exposure)
            m, lo, hi = mean_ci(np.concatenate(offs))
            me, loe, hie = mean_ci(np.concatenate(expo))
            row.update(visit_hours=VISIT_H[name],
                       per_visit_meanfield=BETA * sc.q.mean() * T
                       * (1 - np.exp(-GAMMA * T)) / (GAMMA * T),
                       per_visit_sim=m, per_visit_sim_lo=lo,
                       per_visit_sim_hi=hi, per_visit_exposure_sim=me,
                       per_visit_exposure_lo=loe, per_visit_exposure_hi=hie,
                       n_rep_visit=n, seconds=time.time() - t0)
            out[key] = row
            json.dump(out, open(PATH, "w"), indent=1, default=float)
            print(key, {k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in row.items()}, flush=True)


if __name__ == "__main__":
    main()
