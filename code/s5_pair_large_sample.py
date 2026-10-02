"""s5_pair_large_sample.py -- large-sample counts of the secondary cases of a lone, uniformly placed index
case in the office scene.

The same quantity, R_hat for a uniform index case with only the index transmitting, is estimated in several
places of the article from samples of 2,000 to 40,000 runs (tables of Sections 5 and 6 and of the Supplementary Material), and its
expectation is known exactly (R1_bar, pair solve).  One of the 20,000-run batches (mobility table, D0 = 10,
seed 301) lies 3.5 standard errors below the exact value.  This script draws fresh, larger samples with
other seeds so that the headline numbers of the article can be quoted from one precise estimate:

    D0 = 1   : 20 batches x 20,000 runs
    D0 = 10  : 20 batches x 20,000 runs
    D0 = 100 : 20 batches x  5,000 runs

and re-runs the stored batch of the mobility table (seed 301, 20,000 runs, D0 = 10) to confirm that it is
reproduced bit for bit (a chance-low batch, not a defect).

    python s5_pair_large_sample.py            # about 15 minutes on one core, < 1 GB; checkpointed per batch
Output: ../data/s5_pair_large_sample.json
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
from bsc_sim.scenes import load_scene                                  # noqa: E402
from bsc_sim.sim import simulate                                       # noqa: E402

BETA, GAMMA = 0.5, 0.14
PATH = os.path.join(HERE, "..", "data", "s5_pair_large_sample.json")
PLAN = {1.0: (20, 20000), 10.0: (20, 20000), 100.0: (20, 5000)}
SEED0 = 910000


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    for D0, (n_batch, n_rep) in PLAN.items():
        key = f"office|D0={D0:g}"
        row = out.get(key, dict(D0=D0, n_batch=n_batch, n_rep_per_batch=n_rep, batch_means=[], batch_vars=[],
                                seeds=[], seconds=0.0))
        sc = load_scene("office", D0=D0)
        if "R1bar_exact" not in row:
            row["R1bar_exact"] = float((sc.N - 1) * th.pair_infection_probability(sc, BETA, GAMMA).mean())
            row["R0_ngm"] = float(th.R0_dense(sc, BETA, GAMMA))
        while len(row["batch_means"]) < n_batch:
            b = len(row["batch_means"])
            seed = SEED0 + int(1000 * D0) + b
            t0 = time.time()
            r = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=seed, gmax=0, t_max=600.0)
            x = np.asarray(r.index_offspring, float)
            row["batch_means"].append(float(x.mean()))
            row["batch_vars"].append(float(x.var(ddof=1)))
            row["seeds"].append(seed)
            row["seconds"] += time.time() - t0
            n = len(row["batch_means"]) * n_rep
            m = float(np.mean(row["batch_means"]))
            se = float(np.sqrt(np.mean(row["batch_vars"]) / n))
            row.update(n_runs=n, mean=m, se=se, z_vs_exact=(m - row["R1bar_exact"]) / se,
                       batch_sd_observed=float(np.std(row["batch_means"], ddof=1)) if b else None,
                       batch_sd_expected=float(np.sqrt(np.mean(row["batch_vars"]) / n_rep)),
                       batch_min=float(min(row["batch_means"])), batch_max=float(max(row["batch_means"])))
            out[key] = row
            json.dump(out, open(PATH, "w"), indent=1)
            print(key, b + 1, n_batch, round(m, 4), round(se, 4), flush=True)
    key = "rerun_mobility_table_seed301|D0=10"
    if key not in out:
        sc = load_scene("office", D0=10.0)
        r = simulate(sc, BETA, GAMMA, n_rep=20000, seed=301, gmax=0, t_max=600.0)
        x = np.asarray(r.index_offspring, float)
        stored = json.load(open(os.path.join(ROOT, "code", "bsc_sim", "data", "02_mobility.json")))["D0_10.0"]
        out[key] = dict(mean=float(x.mean()), se=float(x.std(ddof=1) / np.sqrt(len(x))),
                        stored_mean=stored["sim_offspring_uniform"],
                        identical_to_stored=bool(float(x.mean()) == stored["sim_offspring_uniform"]),
                        R1bar_exact=out["office|D0=10"]["R1bar_exact"])
        out[key]["z_vs_exact"] = (out[key]["mean"] - out[key]["R1bar_exact"]) / out[key]["se"]
        json.dump(out, open(PATH, "w"), indent=1)
        print(key, out[key], flush=True)


if __name__ == "__main__":
    main()
