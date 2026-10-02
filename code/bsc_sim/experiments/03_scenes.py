#!/usr/bin/env python
"""The four scenes: theory and stochastic simulation.

For every scene, mobility level D0 in {1, 10, 100} lattice^2/day and contact
kernel (physical range 1 m; plus the strict same-cell rule where it differs)
this script computes

  theory      R0 (NGM), its bounds, early growth rate, branching-process
              major-outbreak probability, pair-level R1, mean-field ODE curve
  simulation  n_rep exact stochastic epidemics from one random index case:
              major-outbreak probability, attack rate, peak prevalence and
              time, duration, counted secondary cases of the index case
  per-cell    the same scene under the per-cell rule (beta q n_I / n_total,
              same cell)

Outputs: data/03_scenes.json (summary), data/03_scenes_table.csv,
data/03_curves_<scene>_<config>.npz (epidemic curves and per-run outcomes).
"""
import csv
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from bsc_sim import theory as th                       # noqa: E402
from bsc_sim.scenes import SCENE_NAMES, load_scene     # noqa: E402
from bsc_sim.sim import simulate, summarize            # noqa: E402

BETA, GAMMA = 0.5, 0.14
ELL_C = 1.0                    # physical contact range, metres
DATA = os.path.join(HERE, "..", "data")
PATH = os.path.join(DATA, "03_scenes.json")
T_MAX, DT_OUT = 400.0, 0.5
THRESH = 0.10                  # attack-rate threshold defining a major outbreak


def perron_radius(K):
    return float(np.max(np.linalg.eigvals(K).real))


def run_config(name, D0, kernel, n_rep, seed, with_pair_ngm):
    sc = load_scene(name, D0=D0)
    rad = (ELL_C / sc.a_m) if kernel == "range1m" else 0.0
    P = th.contact_kernel(sc, rad)
    row = dict(scene=name, D0=D0, kernel=kernel, N=sc.N, cells=sc.M,
               occupancy=sc.N / sc.M, a_m=sc.a_m,
               kernel_cells=float(P.getnnz() / sc.M),
               q_mean=float(sc.q.mean()), q_max=float(sc.q.max()),
               D_min=float(sc.D.min()), D_max=float(sc.D.max()))
    # ---- theory
    t0 = time.time()
    R0, v, u = th.R0_dense(sc, BETA, GAMMA, P=P, return_vectors=True)
    lo, hi = th.R0_bounds(sc, BETA, GAMMA)
    s = th.extinction_probability(sc, BETA, GAMMA, P=P)
    row.update(R0_ngm=R0, R0_wellmixed=lo, R0_frozen=hi,
               growth_rate=th.growth_rate(sc, BETA, GAMMA, P=P),
               p_major_branching=float(1 - s.mean()),
               p_major_homog=float(1 - 1 / R0) if R0 > 1 else 0.0)
    p = th.pair_infection_probability(sc, BETA, GAMMA, P=P)
    row["R1_pair_uniform_index"] = float((sc.N - 1) * p.mean())
    if with_pair_ngm:
        row["R1_pair_ngm"] = perron_radius(th.pair_ngm(sc, BETA, GAMMA, P=P))
    t_ode = np.arange(0.0, T_MAX + DT_OUT, DT_OUT)
    ode = th.meanfield_ode(sc, BETA, GAMMA, t_ode, I0=1.0, P=P)
    row.update(ode_attack=float(1 - ode["S"][-1] / sc.N),
               ode_peak_prev=float(ode["I"].max() / sc.N),
               ode_peak_time=float(t_ode[ode["I"].argmax()]),
               wellmixed_attack_at_R0=th.wellmixed_final_size(R0),
               theory_seconds=time.time() - t0)
    # ---- simulation, pair-hazard rule
    t0 = time.time()
    res = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=seed, P=P, t_max=T_MAX,
                   dt_out=DT_OUT)
    row.update(summarize(res, THRESH))
    row["sim_seconds"] = time.time() - t0
    # second-generation reproduction number counted from the infection tree
    g = res.generation_sizes(3)
    row["gen1_mean"] = float(g[:, 1].mean())
    row["gen2_over_gen1"] = float(g[:, 2].sum() / max(g[:, 1].sum(), 1))
    # index case alone transmitting (no competition) -> compare with R1
    r0 = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=seed + 1, P=P, gmax=0,
                  t_max=T_MAX)
    off = r0.index_offspring
    row["R_index_alone_mean"] = float(off.mean())
    row["R_index_alone_se"] = float(off.std(ddof=1) / np.sqrt(len(off)))
    row["R_index_alone_var"] = float(off.var(ddof=1))
    # ---- simulation, per-cell rule (same cell, n_I / n_total)
    rt = simulate(sc, BETA, GAMMA, n_rep=n_rep, seed=seed + 2, mode="percell",
                  t_max=T_MAX, dt_out=DT_OUT)
    st = summarize(rt, THRESH)
    row.update(percell_rule_index_offspring=st["index_offspring_mean"],
               percell_rule_p_major=st["p_major"],
               percell_rule_attack_all=st["attack_all_mean"])
    # ---- curves for plotting
    maj = res.major(THRESH)
    Ifrac = res.I / sc.N
    curves = dict(t=res.t, mean_all=Ifrac.mean(axis=0),
                  ode=ode["I"] / sc.N,
                  attack=res.attack_rate, peak_prev=res.peak_I / sc.N,
                  peak_time=res.peak_time, duration=res.t_ext,
                  index_offspring=res.index_offspring, major=maj)
    if maj.sum() > 0:
        curves.update(mean_major=Ifrac[maj].mean(axis=0),
                      q05_major=np.quantile(Ifrac[maj], 0.05, axis=0),
                      q50_major=np.quantile(Ifrac[maj], 0.50, axis=0),
                      q95_major=np.quantile(Ifrac[maj], 0.95, axis=0))
    np.savez_compressed(
        os.path.join(DATA, f"03_curves_{name}_D{D0:g}_{kernel}.npz"), **curves)
    return row


def main():
    out = json.load(open(PATH)) if os.path.exists(PATH) else {}
    for D0 in (1.0, 10.0, 100.0):
        for name in SCENE_NAMES:
            sc = load_scene(name)
            kernels = ["range1m"] + (["samecell"] if sc.a_m <= ELL_C else [])
            for kernel in kernels:
                key = f"{name}|D0={D0:g}|{kernel}"
                if key in out:
                    continue
                n_rep = 2000
                t0 = time.time()
                row = run_config(name, D0, kernel, n_rep, seed=1000,
                                 with_pair_ngm=(D0 == 1.0))
                out[key] = row
                json.dump(out, open(PATH, "w"), indent=1, default=float)
                print(f"{key}: R0={row['R0_ngm']:.2f} R1={row['R1_pair_uniform_index']:.2f} "
                      f"Rsim={row['R_index_alone_mean']:.2f} Pmaj={row['p_major']:.3f} "
                      f"(branch {row['p_major_branching']:.3f}) AR|maj="
                      f"{row.get('attack_major_mean', float('nan')):.3f} "
                      f"[{time.time() - t0:.0f} s]", flush=True)
    # flat CSV
    keys = sorted({k for r in out.values() for k in r})
    with open(os.path.join(DATA, "03_scenes_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in out.values():
            w.writerow(r)


if __name__ == "__main__":
    main()
