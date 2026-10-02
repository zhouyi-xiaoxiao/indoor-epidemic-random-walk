#!/usr/bin/env python
"""Where do infections happen?  Theory vs simulation.

The lowest diffusion mode |phi_0|^2 is constant for reflecting walls and
predicts nothing.  The predictors are

  generation-g birth-cell distribution   proportional to  K^g 1   (uniform index)
  asymptotic distribution of birth cells  = right Perron vector of K

with K the next-generation matrix: the mean-field K = P^T diag(beta q) G for
many people per cell, or the pair-level K1 for discrete people.

For each scene we count, over many simulated outbreaks started from a
uniformly placed index case, the cell in which every case of generation 1, 2
and 3 was infected, and compare the three maps with K^g 1 and K1^g 1.

Outputs data/04_spatial_<scene>.npz and data/04_spatial.json.
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
from bsc_sim.sim import simulate                       # noqa: E402

BETA, GAMMA = 0.5, 0.14
ELL_C = 1.0
DATA = os.path.join(HERE, "..", "data")


def pearson_ci(x, y, counts_total, rng, n_boot=1000):
    """Pearson r between predicted and observed cell probabilities, with a
    parametric bootstrap interval (multinomial resampling of the observed
    counts)."""
    r = float(np.corrcoef(x, y)[0, 1])
    p = y / y.sum()
    bs = []
    for _ in range(n_boot):
        c = rng.multinomial(counts_total, p)
        bs.append(np.corrcoef(x, c)[0, 1])
    return r, float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def main():
    path = os.path.join(DATA, "04_spatial.json")
    out = json.load(open(path)) if os.path.exists(path) else {}
    rng = np.random.default_rng(404)
    for name in SCENE_NAMES:
        for D0 in (1.0, 10.0):
            key = f"{name}|D0={D0:g}"
            if key in out:
                continue
            t0 = time.time()
            sc = load_scene(name, D0=D0)
            P = th.contact_kernel(sc, ELL_C / sc.a_m)
            K = th.ngm(sc, BETA, GAMMA, P=P)
            K1 = th.pair_ngm(sc, BETA, GAMMA, P=P)
            R0, v, _ = th.R0_dense(sc, BETA, GAMMA, P=P, return_vectors=True)
            n_rep = 40000
            res = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=4000, P=P,
                           gmax=2, t_max=600.0, dt_out=5.0)
            row = dict(scene=name, D0=D0, n_rep=n_rep, R0_ngm=R0)
            save = dict(perron_mf=v, q=sc.q)
            pred_mf = np.ones(sc.M) / sc.M
            pred_pair = np.ones(sc.M) / sc.M
            for g in (1, 2, 3):
                pred_mf = K @ pred_mf
                pred_pair = K1 @ pred_pair
                sel = res.gen == g
                counts = np.bincount(res.site_inf[sel], minlength=sc.M)
                tot = int(counts.sum())
                pm = pred_mf / pred_mf.sum()
                pp = pred_pair / pred_pair.sum()
                r_mf = pearson_ci(pm, counts.astype(float), tot, rng)
                r_pr = pearson_ci(pp, counts.astype(float), tot, rng)
                # correlation ceiling: prediction = truth, multinomial noise only
                emp = counts / tot
                ceil = np.mean([np.corrcoef(emp, rng.multinomial(tot, emp))[0, 1]
                                for _ in range(200)])
                # goodness of fit of the pair-level prediction
                chi2 = float(np.sum((counts - tot * pp) ** 2 / (tot * pp)))
                row[f"gen{g}"] = dict(
                    n_cases=tot, r_meanfield=r_mf[0], r_meanfield_lo=r_mf[1],
                    r_meanfield_hi=r_mf[2], r_pair=r_pr[0], r_pair_lo=r_pr[1],
                    r_pair_hi=r_pr[2], r_noise_ceiling=float(ceil),
                    chi2_pair=chi2, dof=sc.M - 1,
                    mean_cases_per_run=tot / n_rep,
                    pred_cases_per_run_mf=float(pred_mf.sum()),
                    pred_cases_per_run_pair=float(pred_pair.sum()))
                save[f"counts_gen{g}"] = counts
                save[f"pred_mf_gen{g}"] = pm
                save[f"pred_pair_gen{g}"] = pp
            np.savez_compressed(
                os.path.join(DATA, f"04_spatial_{name}_D{D0:g}.npz"), **save)
            row["seconds"] = time.time() - t0
            out[key] = row
            json.dump(out, open(path, "w"), indent=1, default=float)
            print(key, {g: (round(row[g]["r_meanfield"], 3),
                            round(row[g]["r_pair"], 3),
                            round(row[g]["r_noise_ceiling"], 3))
                        for g in ("gen1", "gen2", "gen3")},
                  f"[{row['seconds']:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
