#!/usr/bin/env python
"""Why is the simulated outbreak probability below the mean-field branching
value?  Three predictions for the probability of a major outbreak, per scene:

  (1) mean-field spatial branching process (Poisson contacts along the walk)
  (2) Galton-Watson process with the *simulated* offspring distribution of a
      lone index case (discrete people, pair saturation and over-dispersion
      included; clustering of later generations still ignored)
  (3) full simulation (data/03_scenes.json)

Output data/03b_outbreak_probability.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import SCENE_NAMES, load_scene     # noqa: E402
from bsc_sim.sim import simulate                       # noqa: E402

BETA, GAMMA, ELL_C = 0.5, 0.14, 1.0
DATA = os.path.join(HERE, "..", "data")
S = json.load(open(os.path.join(DATA, "03_scenes.json")))
out = {}
for D0 in (1.0, 10.0):
    for name in SCENE_NAMES:
        sc = load_scene(name, D0=D0)
        P = th.contact_kernel(sc, ELL_C / sc.a_m)
        n = 40000
        r = simulate(sc, BETA, GAMMA, n_rep=n, seed=3300, P=P, gmax=0,
                     t_max=600.0, dt_out=100.0)
        off = r.index_offspring
        pk = np.bincount(off) / n
        k = np.arange(len(pk))
        s = 0.0
        for _ in range(100000):                 # smallest root of s = G(s)
            new = float(np.sum(pk * s ** k))
            if abs(new - s) < 1e-13:
                break
            s = new
        row = S[f"{name}|D0={D0:g}|range1m"]
        out[f"{name}|D0={D0:g}"] = dict(
            scene=name, D0=D0, n_rep=n, offspring_mean=float(off.mean()),
            offspring_var=float(off.var(ddof=1)),
            var_over_mean=float(off.var(ddof=1) / off.mean()),
            p_zero=float(pk[0]),
            p_major_gw_empirical=1.0 - s,
            p_major_poisson_same_mean=None,
            p_major_branching_meanfield=row["p_major_branching"],
            p_major_sim=row["p_major"], p_major_sim_lo=row["p_major_lo"],
            p_major_sim_hi=row["p_major_hi"])
        # Poisson offspring with the same mean, for reference
        m = off.mean()
        s2 = 0.0
        for _ in range(100000):
            new = np.exp(-m * (1 - s2))
            if abs(new - s2) < 1e-13:
                break
            s2 = new
        out[f"{name}|D0={D0:g}"]["p_major_poisson_same_mean"] = float(1 - s2)
        print(out[f"{name}|D0={D0:g}"], flush=True)
        json.dump(out, open(os.path.join(DATA, "03b_outbreak_probability.json"), "w"),
                  indent=1)
