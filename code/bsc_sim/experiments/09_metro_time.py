#!/usr/bin/env python
"""Time-dependent parameters in the metro scene.

The model allows D(r, t) = D(r) f(t) and q(r, t) = q(r) g(t) with a
fixed layout and a fixed population.  Here the crowd-density signal is

    rho(t) = rho_base + rho_peak * sum_{peaks at 08:00, 18:00}
                           exp(-(t - t_peak)^2 / (2 sigma^2)),   period 1 day,

with rho_base = 1.5, rho_peak = 3.1 persons/m^2 (peak = crush load 4.6) and
sigma = 1 h.  It drives the mobility through the density rule of the scene and the
transmission multiplier m(t) = 1.0 in peak hours (rho >= 3), 0.8 otherwise
(an assumption).  The day is cut into 24 one-hour segments; inside a segment
all rates are constant, so the event-driven simulation stays exact.

Arms (N = 310 people in every arm, contact range 1 m)
  peak     constant parameters at the peak (the default metro scene)
  varying  hourly schedule
  mean     constant parameters equal to the time averages of D(r,t) and m(t)

Outputs data/09_metro_time.json, data/09_metro_curves.npz.
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
ELL_C = 1.0
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "09_metro_time.json")
RHO_BASE, RHO_PEAK, SIGMA_H = 1.5, 3.1, 1.0
AREA = 67.5


def rho_of_hour(h):
    h = np.asarray(h, float)
    r = RHO_BASE + np.zeros_like(h)
    for pk in (8.0, 18.0):
        d = np.minimum(np.abs(h - pk), 24 - np.abs(h - pk))
        r = r + RHO_PEAK * np.exp(-d ** 2 / (2 * SIGMA_H ** 2))
    return r


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    hours = np.arange(24) + 0.5
    rho_h = rho_of_hour(hours)
    m_h = np.where(rho_h >= 3.0, 1.0, 0.8)
    curves = {"hours": hours, "rho": rho_h, "m": m_h}
    n_rep = 2000
    for D0 in (1.0, 100.0):
        sc = load_scene("metro", D0=D0)              # N = 310, rho = 4.59
        P = th.contact_kernel(sc, ELL_C / sc.a_m)
        # exact rho for each hour (not rounded to whole people)
        zones = sc.meta["zones"]
        coef = np.array([zones[z]["D_coef"] for z in sc.site_zone])
        expo = np.array([zones[z]["D_exp"] for z in sc.site_zone])
        D_h = [D0 * coef / r ** expo for r in rho_h]
        D_mean = np.mean(D_h, axis=0)
        m_mean = float(m_h.mean())
        arms = {
            "peak": dict(D=sc.D, schedule=None, m=1.0),
            "varying": dict(D=D_h[0], schedule=dict(
                period=1.0,
                segments=[((k + 1) / 24.0, D_h[k], float(m_h[k]))
                          for k in range(24)]), m=None),
            "mean": dict(D=D_mean, schedule=dict(
                period=1.0, segments=[(1.0, D_mean, m_mean)]), m=m_mean),
        }
        for arm, cfg in arms.items():
            key = f"D0={D0:g}|{arm}"
            t0 = time.time()
            res = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=900, P=P,
                           D=cfg["D"], schedule=cfg["schedule"], t_max=300.0,
                           dt_out=0.25)
            s = summarize(res)
            row = dict(D0=D0, arm=arm, n_rep=n_rep)
            if cfg["m"] is not None:
                L = th.generator(sc, cfg["D"])
                row["R0_ngm"] = th.R0_dense(sc, BETA * cfg["m"], GAMMA, P=P, L=L)
                p = th.pair_infection_probability(sc, BETA * cfg["m"], GAMMA,
                                                  P=P, L=L)
                row["R1_pair_uniform"] = float((sc.N - 1) * p.mean())
            row.update(s)
            row["seconds"] = time.time() - t0
            out[key] = row
            maj = res.major()
            curves[f"{key}|t"] = res.t
            curves[f"{key}|mean_all"] = (res.I / sc.N).mean(axis=0)
            if maj.sum():
                curves[f"{key}|mean_major"] = (res.I[maj] / sc.N).mean(axis=0)
                curves[f"{key}|q05"] = np.quantile(res.I[maj] / sc.N, 0.05, axis=0)
                curves[f"{key}|q95"] = np.quantile(res.I[maj] / sc.N, 0.95, axis=0)
            json.dump(out, open(PATH, "w"), indent=1, default=float)
            print(key, f"R0={row.get('R0_ngm', float('nan')):.2f} "
                  f"Pmaj={s['p_major']:.3f} AR|maj={s.get('attack_major_mean', float('nan')):.3f} "
                  f"peak={s.get('peak_prev_major_mean', float('nan')):.3f} "
                  f"tpk={s.get('peak_time_major_mean', float('nan')):.2f} "
                  f"[{time.time() - t0:.0f} s]", flush=True)
    out["schedule"] = dict(hours=hours.tolist(), rho=rho_h.tolist(),
                           m=m_h.tolist(), m_mean=float(m_h.mean()),
                           rho_mean=float(rho_h.mean()))
    json.dump(out, open(PATH, "w"), indent=1, default=float)
    np.savez_compressed(os.path.join(DATA, "09_metro_curves.npz"), **curves)


if __name__ == "__main__":
    main()
