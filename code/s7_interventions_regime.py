#!/usr/bin/env python
"""s7_interventions_regime.py -- what remains of the effects of Section 6.4 (Supplementary Sections S5.6-S5.11) at higher mobility.

The experiments of the simulations were run at D0 = 1 and 10 cell^2/day.  Section 2 of the article
estimates 2.4e3 cell^2/day for people who walk 5 % of the time (a = 1.5 m).  This script evaluates, at
D0 = 1, 10, 100 and 2400,

 (det)  the mean-field R0 and the pair-level numbers (R1-bar: exact expected number of secondary cases of a
        uniformly placed lone index case; rho(K1) for the two-zone room) for
          A. the two-zone room of Supplementary Section S5.6 at contrast c = 0 and c = 0.5833 (N = 100, same-cell kernel, FD);
          B. the barrier room of Supplementary Section S5.7 at nominal barrier density 0, 0.10, 0.30 (FD and DD), mean over
             the 40 layouts of s7_interventions_barriers.py (same seeds);
          C. the office scene: baseline, ventilation (q/2), capacity 50 % (FD and DD);
 (sim)  full epidemics at D0 = 100 only (the cost of the exact simulator grows in proportion to D0):
          A. two-zone room, c = 0 and c = 0.5833, 2000 epidemics each, uniform index;
          B. barrier room, FD, densities 0, 0.10, 0.30, identical starts, the first 20 layouts x 250 epidemics.

Usage:  python s7_interventions_regime.py det      (about 3 minutes)
        python s7_interventions_regime.py sim      (about 12 minutes)
Output: ../data/s7_interventions_regime.json
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "code", "bsc_sim", "src"))
sys.path.insert(0, HERE)
from bsc_sim import theory as th                                       # noqa: E402
from bsc_sim.scenes import Scene, load_scene, uniform_room             # noqa: E402
from bsc_sim.sim import simulate, summarize                            # noqa: E402
import s7_interventions_barriers as bar                                # noqa: E402

BETA, GAMMA = 0.5, 0.14
QBAR = 1.2
C_EXAMPLE = 1 - 0.5 / QBAR
D_LIST = (1.0, 10.0, 100.0, 2400.0)
PATH = os.path.join(HERE, "..", "data", "s7_interventions_regime.json")


def two_zone(c, D0, N=100):
    """Two-zone room of code/bsc_sim/experiments/05_heterogeneity.py (two-zone variant with zone-dependent mobility)."""
    nx, ny = 20, 13
    zone = np.full((ny, nx), "W")
    zone[:, 7:13] = "C"
    access = np.ones((ny, nx), dtype=bool)
    q = np.where(zone == "W", QBAR * (1 + c * 0.3 / 0.7), QBAR * (1 - c))
    D = np.where(zone == "W", 0.3, 2.0) * D0
    return Scene(name=f"twozone_c{c:g}", nx=nx, ny=ny, a_m=1.5, zone=zone, access=access, D_grid=D, q_grid=q,
                 N=N, meta=dict(c=c, D0=D0))


def save(out):
    json.dump(out, open(PATH, "w"), indent=1, default=float)


def det(out):
    res = out.setdefault("det", {})
    t0 = time.time()
    for D0 in D_LIST:
        # A. two-zone room
        for c in (0.0, C_EXAMPLE):
            key = f"twozone|D0={D0:g}|c={c:.4f}"
            if key in res:
                continue
            sc = two_zone(c, D0)
            p = th.pair_infection_probability(sc, BETA, GAMMA)
            K1 = th.pair_ngm(sc, BETA, GAMMA)
            res[key] = dict(R0=th.R0_dense(sc, BETA, GAMMA), R1bar=float((sc.N - 1) * p.mean()),
                            rho_K1=float(np.max(np.linalg.eigvals(K1).real)))
            save(out)
            print(key, res[key], f"[{time.time() - t0:.0f} s]", flush=True)
        # B. barrier room, mean over the 40 layouts
        key = f"barriers|D0={D0:g}"
        if key not in res:
            rho0 = (bar.N - 1) / (bar.NX * bar.NY)
            acc_rows = {(f, cal): [] for f in bar.FRACS for cal in ("FD", "DD")}
            for lay, seed, acc, rng in bar.valid_layouts():
                for f in bar.FRACS:
                    sc = uniform_room(bar.NX, bar.NY, bar.N, D=D0, q=bar.Q, barriers=~acc[f])
                    c_idx = int(sc.site_index[bar.CENTRE[1], bar.CENTRE[0]])
                    for cal in ("FD", "DD"):
                        rho_ref = None if cal == "FD" else rho0
                        scale = 1.0 if cal == "FD" else ((bar.N - 1) / sc.M) / rho0
                        p = th.pair_infection_probability(sc, BETA, GAMMA, rho_ref=rho_ref)
                        acc_rows[(f, cal)].append((th.R0_dense(sc, BETA, GAMMA) * scale,
                                                   (bar.N - 1) * p.mean(), (bar.N - 1) * p[c_idx].mean()))
            d = {}
            for (f, cal), rows in acc_rows.items():
                a = np.array(rows)
                d[f"rhoB={f:.2f}|{cal}"] = dict(R0=float(a[:, 0].mean()), R1bar=float(a[:, 1].mean()),
                                                R1bar_se=float(a[:, 1].std(ddof=1) / np.sqrt(len(a))),
                                                R1_centre=float(a[:, 2].mean()), n_layout=len(a))
            res[key] = d
            save(out)
            print(key, {k: (round(v["R0"], 3), round(v["R1bar"], 3)) for k, v in d.items()},
                  f"[{time.time() - t0:.0f} s]", flush=True)
        # C. office scene
        key = f"office|D0={D0:g}"
        if key not in res:
            sc = load_scene("office", D0=D0)
            P = th.contact_kernel(sc, 1.0 / sc.a_m)
            rho_ref0 = (sc.N - 1) / sc.M
            R0 = th.R0_dense(sc, BETA, GAMMA, P=P)
            lo, hi = th.R0_bounds(sc, BETA, GAMMA)
            d = dict(R0=R0, wellmixed=lo, frozen=hi)
            d["R1bar_baseline"] = float((sc.N - 1) * th.pair_infection_probability(sc, BETA, GAMMA, P=P).mean())
            d["R1bar_ventilation"] = float((sc.N - 1) * th.pair_infection_probability(
                sc, BETA, GAMMA, P=P, q=0.5 * sc.q).mean())
            d["R1bar_capacity50_FD"] = float(49 * th.pair_infection_probability(sc, BETA, GAMMA, N=50, P=P).mean())
            d["R1bar_capacity50_DD"] = float(49 * th.pair_infection_probability(
                sc, BETA, GAMMA, N=50, P=P, rho_ref=rho_ref0).mean())
            d["R0_capacity50_DD"] = R0 * (49 / sc.M) / rho_ref0
            res[key] = d
            save(out)
            print(key, {k: round(v, 4) for k, v in d.items()}, f"[{time.time() - t0:.0f} s]", flush=True)


def sim(out, D0=100.0):
    res = out.setdefault("sim", {})
    t0 = time.time()
    for c in (0.0, C_EXAMPLE):
        key = f"twozone|D0={D0:g}|c={c:.4f}"
        if key in res:
            continue
        sc = two_zone(c, D0)
        full = simulate(sc, BETA, GAMMA, n_rep=2000, seed=75000, t_max=300.0, dt_out=0.5)
        s = summarize(full)
        res[key] = {k: s[k] for k in ("n_rep", "p_major", "p_major_lo", "p_major_hi", "attack_all_mean",
                                      "attack_all_sd", "attack_major_mean", "peak_prev_major_mean",
                                      "peak_time_major_mean", "n_unfinished")}
        save(out)
        print(key, res[key], f"[{time.time() - t0:.0f} s]", flush=True)
    key = f"barriers|D0={D0:g}|FD"
    if key not in res:
        n_lay, n_rep = 20, 250
        per = {f: dict(attack=[], p_major=[], peak_major=[], peak_day_major=[]) for f in bar.FRACS}
        for lay, seed, acc, rng in bar.valid_layouts():
            if lay >= n_lay:
                break
            ys, xs = np.nonzero(acc[bar.FRACS[-1]])
            pick = rng.integers(0, len(xs), size=(n_rep, bar.N))
            sx, sy = xs[pick], ys[pick]
            sx[:, 0], sy[:, 0] = bar.CENTRE
            for f in bar.FRACS:
                sc = uniform_room(bar.NX, bar.NY, bar.N, D=D0, q=bar.Q, barriers=~acc[f])
                r = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=76000 + lay, init_pos=sc.site_index[sy, sx],
                             t_max=300.0, dt_out=0.5)
                ar = r.attack_rate
                maj = ar >= 0.1
                per[f]["attack"].append(float(ar.mean()))
                per[f]["p_major"].append(float(maj.mean()))
                per[f]["peak_major"].append(float((r.peak_I[maj] / bar.N).mean()))
                per[f]["peak_day_major"].append(float(r.peak_time[maj].mean()))
            print(f"sim barriers D={D0:g} layout {lay} [{time.time() - t0:.0f} s]", flush=True)

        def ms(v):
            v = np.asarray(v, float)
            return [float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))]

        d = dict(n_layout=n_lay, n_rep=n_rep, arms={}, paired={})
        for f in bar.FRACS:
            d["arms"][f"rhoB={f:.2f}"] = {k: ms(v) for k, v in per[f].items()}
            if f > 0:
                d["paired"][f"rhoB={f:.2f}"] = {k: ms(np.asarray(per[f][k]) - np.asarray(per[0.0][k]))
                                                for k in per[f]}
        res[key] = d
        save(out)
        print(key, d["arms"], d["paired"], flush=True)


if __name__ == "__main__":
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    which = sys.argv[1] if len(sys.argv) > 1 else "det"
    if which == "det":
        det(out)
    else:
        sim(out)
