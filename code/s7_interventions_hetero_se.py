"""s7_interventions_hetero_se.py -- uncertainties of the whole-epidemic outcomes of the two-zone room.

The stored results of code/bsc_sim/experiments/05_heterogeneity.py give the probability of a major
outbreak with its Wilson interval but the mean attack rate without a standard error, and for N = 1000 they
rest on 300 epidemics per point, too few to resolve the dependence on the contrast c.  This script

  (A) re-runs every stored batch of full epidemics from its seed (503 for the main runs, 511 for the
      controls), asserts that the stored probability of a major outbreak and mean attack rate are reproduced
      exactly, and records the standard deviation of the attack rate, hence its standard error;
  (B) for N = 1000 adds 2,000 further epidemics per point with another seed (the six contrasts at D0 = 1 and
      D0 = 10) and reports the 2,000-run sample next to the stored 300-run one.

    python s7_interventions_hetero_se.py            # about 25 minutes on one core, < 1 GB; checkpointed
Output: ../data/s7_interventions_hetero_se.json
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))          # repository root
SIM = os.path.join(ROOT, "code", "bsc_sim")
sys.path.insert(0, os.path.join(SIM, "src"))
from bsc_sim.sim import simulate, summarize, wilson                    # noqa: E402

spec = importlib.util.spec_from_file_location("het", os.path.join(SIM, "experiments", "05_heterogeneity.py"))
het = importlib.util.module_from_spec(spec)
spec.loader.exec_module(het)

BETA, GAMMA = 0.5, 0.14
PATH = os.path.join(HERE, "..", "data", "s7_interventions_hetero_se.json")
MAIN = json.load(open(os.path.join(SIM, "data", "05_heterogeneity.json")))
EXTRA = json.load(open(os.path.join(SIM, "data", "05_heterogeneity_extra.json")))
N_NEW, SEED_NEW = 2000, 5503


def record(s):
    n = s["n_rep"]
    return dict(n=n, p_major=s["p_major"], p_major_lo=s["p_major_lo"], p_major_hi=s["p_major_hi"],
                p_major_se=float(np.sqrt(s["p_major"] * (1 - s["p_major"]) / n)),
                attack_all=s["attack_all_mean"], attack_all_sd=s["attack_all_sd"],
                attack_all_se=float(s["attack_all_sd"] / np.sqrt(n)),
                attack_major=s.get("attack_major_mean"), n_major=s["n_major"])


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}

    def save():
        json.dump(out, open(PATH, "w"), indent=1, default=float)

    # (A) stored batches, main runs
    for key, row in MAIN.items():
        if key == "meanfield_grid" or ("stored|" + key) in out:
            continue
        t0 = time.time()
        sc = het.two_zone(row["c"], row["D0"], row["N"])
        full = simulate(sc, BETA, GAMMA, n_rep=row["n_full"], seed=503, t_max=600.0, dt_out=0.5)
        s = summarize(full)
        assert s["p_major"] == row["p_major"] and abs(s["attack_all_mean"] - row["attack_all"]) < 1e-12, key
        out["stored|" + key] = dict(record(s), reproduces_stored=True, seconds=time.time() - t0)
        save()
        print("stored", key, round(s["p_major"], 4), round(s["attack_all_mean"], 4), flush=True)
    # (A) stored batches, controls
    for key, row in EXTRA.items():
        if ("stored|" + key) in out:
            continue
        t0 = time.time()
        sc = het.two_zone(row["c"], row["D0"], 100, row["variant"])
        full = simulate(sc, BETA, GAMMA, n_rep=2000, seed=511, t_max=600.0, dt_out=0.5)
        s = summarize(full)
        assert s["p_major"] == row["p_major"] and abs(s["attack_all_mean"] - row["attack_all"]) < 1e-12, key
        out["stored|" + key] = dict(record(s), reproduces_stored=True, seconds=time.time() - t0)
        save()
        print("stored", key, round(s["p_major"], 4), round(s["attack_all_mean"], 4), flush=True)
    # (B) larger sample at N = 1000
    for D0 in (1.0, 10.0):
        for key, row in MAIN.items():
            if key == "meanfield_grid" or row["N"] != 1000 or row["D0"] != D0 or ("new|" + key) in out:
                continue
            t0 = time.time()
            sc = het.two_zone(row["c"], row["D0"], 1000)
            full = simulate(sc, BETA, GAMMA, n_rep=N_NEW, seed=SEED_NEW, t_max=600.0, dt_out=0.5)
            s = summarize(full)
            rec = record(s)
            old = out["stored|" + key]
            rec["z_p_major_new_minus_stored"] = (rec["p_major"] - old["p_major"]) / float(
                np.hypot(rec["p_major_se"], old["p_major_se"]))
            rec["z_attack_new_minus_stored"] = (rec["attack_all"] - old["attack_all"]) / float(
                np.hypot(rec["attack_all_se"], old["attack_all_se"]))
            out["new|" + key] = dict(rec, seed=SEED_NEW, c=row["c"], D0=D0, N=1000, seconds=time.time() - t0)
            save()
            print("new", key, round(s["p_major"], 4), round(s["attack_all_mean"], 4), round(time.time() - t0), flush=True)
    save()


if __name__ == "__main__":
    main()
