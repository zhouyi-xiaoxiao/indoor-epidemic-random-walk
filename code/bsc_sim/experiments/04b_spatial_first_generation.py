#!/usr/bin/env python
"""Addendum to 04: cell-by-cell test of the pair-level matrix K1 for the
first generation with only the index case transmitting (no competition), where
K1 1/|Omega| is exact.  Reads the predictions stored by 04 and writes
data/04b_first_generation.json."""
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
out = {}
rng = np.random.default_rng(44)
for name in SCENE_NAMES:
    for D0 in (1.0, 10.0):
        sc = load_scene(name, D0=D0)
        P = th.contact_kernel(sc, ELL_C / sc.a_m)
        z = np.load(os.path.join(DATA, f"04_spatial_{name}_D{D0:g}.npz"))
        pp, pm = z["pred_pair_gen1"], z["pred_mf_gen1"]
        # per-cell mean and variance of the number of cases per run.  Cases of
        # one run are clustered around the index case, so cell counts are
        # over-dispersed relative to a multinomial; the test below uses the
        # run-to-run variance instead.
        n_rep, chunk = 100000, 20000
        s1 = np.zeros(sc.M)
        s2 = np.zeros(sc.M)
        tot_run = []
        for c in range(n_rep // chunk):
            r = simulate(sc, BETA, GAMMA, n_rep=chunk, seed=4400 + c, P=P,
                         gmax=0, t_max=600.0, dt_out=100.0)
            rr, aa = np.nonzero(r.gen == 1)
            per = np.zeros((chunk, sc.M), dtype=np.int16)
            np.add.at(per, (rr, r.site_inf[rr, aa]), 1)
            s1 += per.sum(axis=0)
            s2 += (per.astype(np.int64) ** 2).sum(axis=0)
            tot_run.append(per.sum(axis=1))
        tot = int(s1.sum())
        mean = s1 / n_rep                       # cases per run in each cell
        var = s2 / n_rep - mean ** 2
        se = np.sqrt(var / n_rep)
        cases_per_run = tot / n_rep
        emp = s1 / tot
        # compare the *shape*: observed share vs predicted share, with the
        # run-level standard error of the share (delta method, total fixed)
        se_share = se / cases_per_run
        zscore = (emp - pp) / se_share
        chi2 = float(np.sum(zscore ** 2))
        ceil = float(np.mean([np.corrcoef(pp, pp + se_share * rng.standard_normal(sc.M))[0, 1]
                              for _ in range(400)]))
        out[f"{name}|D0={D0:g}"] = dict(
            scene=name, D0=D0, n_rep=n_rep, n_cases=tot,
            r_pair=float(np.corrcoef(pp, emp)[0, 1]),
            r_meanfield=float(np.corrcoef(pm, emp)[0, 1]),
            r_expected_if_exact=ceil, chi2_pair=chi2, dof=sc.M - 1,
            overdispersion=float(np.mean(var / mean)),
            max_abs_z=float(np.abs(zscore).max()),
            cases_per_run=cases_per_run)
        print(out[f"{name}|D0={D0:g}"], flush=True)
        json.dump(out, open(os.path.join(DATA, "04b_first_generation.json"), "w"),
                  indent=1)
