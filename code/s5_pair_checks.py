"""s5_pair_checks.py -- numbers quoted in Section 5 (pair-level theory) that are not stored verbatim in
the research directories, recomputed from the scene files and the stored results.

What it does
  1. finite_size   : reads data/bsc_sim/02_finite_size.json; accuracy of the reflecting-wall
                     closed form against the exact pair solve; number of simulated points outside their
                     95 % interval; rows of Table s5_pair:tab:finite.
  2. mobility      : reads data/bsc_sim/02_mobility.json; rows of Table s5_pair:tab:mobility and
                     the ratios quoted in the text.
  3. torus         : separate check of Proposition s5_pair:prop:torus on rectangular periodic lattices
                     (dense solve of the pair equation against the closed form), limits D -> 0 and
                     D -> infinity, monotonicity in D.
  4. infinite      : elliptic-integral form of g0 against a k-space quadrature; mobility at which
                     R1/R0 reaches 0.5, 0.9, 0.99 on the infinite lattice (validity map).
  5. saturation    : entrywise check 0 <= K1 <= K, R1(x) <= n(x), rho(K1) <= R0 (Proposition
                     s5_pair:prop:saturation) on the office (same-cell kernel), the classroom (1 m kernel)
                     and the metro carriage (same-cell kernel), D0 = 1; gives rho(K1) and R1-bar of the
                     metro carriage with the same-cell kernel (the "not a threshold" example).
  5b. wellmixed    : generation-2/generation-1 ratio in a well-mixed Markovian SIR population of 100 with
                     the design of the last column of Table s5_pair:tab:mobility (200 000 runs, fixed seed).
  6. encounter     : encounter factor of the per-cell simulation rule and the annealed predictions; also copies
                     the re-checks' counts used in the text (finite size, mobility, per-cell rule).

Output: ../data/s5_pair_checks.json and ../data/s5_pair_tables.tex (table bodies, for checking the
typeset tables against the data).

Run:  python s5_pair_checks.py      (about 10 seconds on one core, < 1 GB)

Units: lengths in lattice cells, time in days; beta = 0.5/day, gamma = 0.14/day.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.special import ellipk

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                       # repository root
sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))

from bsc_sim import theory as th             # noqa: E402
from bsc_sim.scenes import load_scene        # noqa: E402

BETA, GAMMA = 0.5, 0.14
SIM = ROOT / "code/bsc_sim"
OUT = HERE.parent / "data"
OUT.mkdir(exist_ok=True)
RES: dict = {"beta": BETA, "gamma": GAMMA}
TEX: list[str] = []


# --------------------------------------------------------------------------- 1. finite size
def finite_size():
    d = json.load(open(SIM / "data/02_finite_size.json"))
    rows = {k: v for k, v in d.items() if "Linf" not in k}
    dev = {}
    for k, v in rows.items():
        if v["pair_exact"] == v["pair_exact"]:          # not NaN
            dev[k] = v["pair_closed_form"] / v["pair_exact"] - 1.0
    worst = max(dev, key=lambda k: abs(dev[k]))
    by_D = {}
    for D in (0.51, 1.0, 10.0):
        ks = [k for k in dev if k.startswith(f"D{D}_")]
        by_D[str(D)] = max(abs(dev[k]) for k in ks)
    exact_out = [k for k, v in rows.items() if v["pair_exact"] == v["pair_exact"]
                 and not (v["sim_lo"] <= v["pair_exact"] <= v["sim_hi"])]
    closed_out = [k for k, v in rows.items()
                  if not (v["sim_lo"] <= v["pair_closed_form"] <= v["sim_hi"])]
    RES["finite_size"] = dict(
        n_points=len(rows), n_exact=len(dev),
        max_abs_rel_dev_closed_vs_exact=abs(dev[worst]), worst_point=worst,
        max_abs_rel_dev_by_D=by_D,
        rel_dev=dev,
        exact_outside_sim_ci=exact_out, closed_outside_sim_ci=closed_out,
        occupancy_used={k: (v["N"] - 1) / v["L"] ** 2 for k, v in rows.items()},
        ratio_to_meanfield={k: v["pair_closed_form"] / (BETA / GAMMA) for k, v in d.items()},
    )
    TEX.append("% Table s5_pair:tab:finite  (src: data/bsc_sim/02_finite_size.json)")
    for D in (1.0, 10.0):
        TEX.append(f"% D = {D:g}")
        for L in (4, 10, 20, 40):
            v = d[f"D{D}_L{L}"]
            ex = "--" if v["pair_exact"] != v["pair_exact"] else f"{v['pair_exact']:.3f}"
            TEX.append(f"{D:g} & " * (L == 4) + "   & " * (L != 4) +
                       f"{L} & {v['N']} & {v['R0_ngm']:.3f} & "
                       f"{v['pair_closed_form']:.3f} & {ex} & {v['sim_offspring']:.3f} "
                       f"({v['sim_lo']:.3f}--{v['sim_hi']:.3f}) \\\\")
        v = d[f"D{D}_Linf"]
        TEX.append(f"   & $\\infty$ & -- & {BETA / GAMMA:.3f} & "
                   f"{v['pair_closed_form']:.3f} & -- & -- \\\\")


# --------------------------------------------------------------------------- 2. mobility
def mobility():
    d = json.load(open(SIM / "data/02_mobility.json"))
    rows = sorted(d.values(), key=lambda r: r["D0"])
    out = {}
    TEX.append("% Table s5_pair:tab:mobility  (src: data/bsc_sim/02_mobility.json)")
    for r in rows:
        out[f"D0={r['D0']:g}"] = dict(
            R1bar_over_R0=r["R1_pair_uniform_index"] / r["R0_ngm"],
            rhoK1_over_R0=r["R1_pair_ngm"] / r["R0_ngm"],
            gen_ratio_over_rhoK1=r["sim_gen2_over_gen1"] / r["R1_pair_ngm"],
            R0_over_counted_perron=r["R0_ngm"] / r["sim_offspring_perron"],
            perron_inside_ci=bool(r["sim_perron_lo"] <= r["R1_pair_ngm"] <= r["sim_perron_hi"]),
            uniform_inside_ci=bool(r["sim_uniform_lo"] <= r["R1_pair_uniform_index"] <= r["sim_uniform_hi"]),
        )
        TEX.append(f"{r['D0']:g} & {r['R0_ngm']:.3f} & {r['R1_pair_ngm']:.3f} & "
                   f"{r['sim_offspring_perron']:.3f} ({r['sim_perron_lo']:.3f}--{r['sim_perron_hi']:.3f}) & "
                   f"{r['R1_pair_uniform_index']:.3f} & "
                   f"{r['sim_offspring_uniform']:.3f} ({r['sim_uniform_lo']:.3f}--{r['sim_uniform_hi']:.3f}) & "
                   f"{r['sim_gen2_over_gen1']:.3f} ({r['sim_gen2_lo']:.3f}--{r['sim_gen2_hi']:.3f}) "
                   f"& {r['n_rep']} \\\\")
    office = load_scene("office")
    bq = BETA * office.q.mean()
    out["office_wellmixed_finite_population"] = bq / GAMMA / (1 + bq / (GAMMA * (office.N - 1)))
    out["office_beta_meanq_over_gamma"] = bq / GAMMA
    RES["mobility"] = out


# --------------------------------------------------------------------------- 3. torus
def torus_generator(L1, L2, D):
    n = L1 * L2
    idx = lambda i, j: (i % L1) * L2 + (j % L2)          # noqa: E731
    A = np.zeros((n, n))
    for i in range(L1):
        for j in range(L2):
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                A[idx(i + di, j + dj), idx(i, j)] += D
    return A - np.diag(A.sum(axis=0))


def torus_closed(L1, L2, N, bq, D):
    k1 = 2 * np.pi * np.arange(L1) / L1
    k2 = 2 * np.pi * np.arange(L2) / L2
    mu = 2 * D * (2 - np.cos(k1)[:, None] - np.cos(k2)[None, :])
    g0 = float(np.mean(1.0 / (GAMMA + 2 * mu)))
    rho = (N - 1) / (L1 * L2)
    return bq / GAMMA / (1 + bq / rho * g0), g0


def torus():
    out = {}
    for (L1, L2, N, D) in ((5, 5, 11, 1.0), (6, 4, 10, 1.0), (7, 3, 9, 0.3), (6, 5, 13, 10.0)):
        n = L1 * L2
        M = torus_generator(L1, L2, D)
        I = np.eye(n)
        M2 = np.kron(M, I) + np.kron(I, M)               # state (x, y) -> x*n + y
        h = (BETA / ((N - 1) / n)) * I.ravel()           # hazard on the diagonal x = y
        u = GAMMA * np.linalg.solve(GAMMA * np.eye(n * n) + np.diag(h) - M2.T, np.ones(n * n))
        p = 1.0 - u.reshape(n, n)
        R1x = (N - 1) / n * p.sum(axis=1)
        Rc, g0 = torus_closed(L1, L2, N, BETA, D)
        out[f"L{L1}x{L2}_N{N}_D{D:g}"] = dict(R1_exact_min=float(R1x.min()), R1_exact_max=float(R1x.max()),
                                              R1_closed=Rc, g0=g0,
                                              abs_err=float(np.max(np.abs(R1x - Rc))))
    # limits and monotonicity on a 20 x 13 torus with N = 100
    L1, L2, N = 20, 13, 100
    Ds = np.geomspace(1e-6, 1e8, 141)
    R = np.array([torus_closed(L1, L2, N, BETA, D)[0] for D in Ds])
    h = BETA / ((N - 1) / (L1 * L2))
    out["limits_20x13_N100"] = dict(
        R1_at_D_small=float(R[0]), limit_D_to_0=(N - 1) / (L1 * L2) * h / (h + GAMMA),
        R1_at_D_large=float(R[-1]), limit_D_to_inf=BETA / GAMMA / (1 + BETA / (GAMMA * (N - 1))),
        increasing_in_D=bool(np.all(np.diff(R) > 0)),
    )
    RES["torus"] = out


# --------------------------------------------------------------------------- 4. infinite lattice
def g0_inf(D):
    k = 8 * D / (GAMMA + 8 * D)
    return (2 / np.pi) * ellipk(k ** 2) / (GAMMA + 8 * D)   # scipy takes the parameter m = k^2


def infinite():
    out = {}
    rho = 99 / 260
    for D in (0.51, 1.0, 10.0, 100.0, 2.4e3, 4.9e4):
        m = 4000
        kk = (np.arange(m) + 0.5) * np.pi / m
        mu = 2 * D * (2 - np.cos(kk)[:, None] - np.cos(kk)[None, :])
        gnum = float(np.mean(1.0 / (GAMMA + 2 * mu))) if D <= 10 else float("nan")
        g = float(g0_inf(D))
        out[f"D={D:g}"] = dict(g0_elliptic=g, g0_quadrature=gnum,
                               g0_large_D_asymptote=float(np.log(64 * D / GAMMA) / (8 * np.pi * D)),
                               R1_inf=BETA / GAMMA / (1 + BETA / rho * g),
                               ratio=1 / (1 + BETA / rho * g))
    for bq in (0.5, 1.25):
        for occ, tag in ((rho, "rho=99/260"), (1.0, "rho=1"), (10.0, "rho=10")):
            for level in (0.5, 0.9, 0.99):
                f = lambda lD: 1 / (1 + bq / occ * g0_inf(10 ** lD)) - level   # noqa: E731
                try:
                    out[f"D_for_ratio_{level}_bq={bq}_{tag}"] = float(10 ** brentq(f, -6, 9))
                except ValueError:
                    out[f"D_for_ratio_{level}_bq={bq}_{tag}"] = None
    RES["infinite"] = out


# --------------------------------------------------------------------------- 5. saturation, K1 <= K
def saturation():
    out = {}
    cases = (("office", 0.0, "same-cell (= 1 m)"), ("classroom", 1.0, "1 m"), ("metro", 0.0, "same-cell"))
    for name, radius_m, label in cases:
        t0 = time.time()
        scn = load_scene(name, D0=1.0)
        P = None if radius_m == 0 else th.contact_kernel(scn, th.kernel_radius_cells(scn, radius_m))
        K = th.ngm(scn, BETA, GAMMA, P=P)
        K1 = th.pair_ngm(scn, BETA, GAMMA, P=P)
        ev = np.linalg.eigvals(K1)
        rho1 = float(np.max(ev.real))
        R0 = float(np.max(np.linalg.eigvals(K).real))
        R1x = K1.sum(axis=0)
        nx = K.sum(axis=0)
        p = th.pair_infection_probability(scn, BETA, GAMMA, P=P)
        R1x_direct = (scn.N - 1) * p.mean(axis=1)
        out[name] = dict(
            kernel=label, N=scn.N, cells=scn.M, rho_bar=(scn.N - 1) / scn.M,
            R0=R0, rho_K1=rho1, R1bar=float(R1x.mean()), R1_min=float(R1x.min()), R1_max=float(R1x.max()),
            meanfield_offspring_mean=float(nx.mean()),
            max_K1_minus_K=float((K1 - K).max()), min_K1=float(K1.min()),
            max_R1x_over_nx=float((R1x / nx).max()),
            max_abs_colsum_minus_direct=float(np.max(np.abs(R1x - R1x_direct))),
            seconds=time.time() - t0,
        )
        print(name, out[name], flush=True)
    # mean time between hops of a passenger in the metro carriage at D0 = 1
    metro = load_scene("metro", D0=1.0)
    exit_rate = -th.generator(metro).diagonal()
    out["metro"]["mean_exit_rate_per_day"] = float(exit_rate.mean())
    out["metro"]["mean_time_between_hops_days"] = float(1.0 / exit_rate.mean())
    out["metro"]["median_cell_time_between_hops_days"] = float(np.median(1.0 / exit_rate))
    # counted values for the metro carriage, same-cell kernel (re-check's second simulator)
    v = json.load(open(SIM / "verify/v05_scenes.json"))["metro|1.0|same"]
    out["metro"]["check_counts"] = dict(n=v["n"], p_major=v["p_major"], p_hi=v["p_hi"],
                                           gen2_over_gen1=v["gen2_over_gen1"],
                                           R_alone=v["R_alone"], R_alone_se=v["R_alone_se"])
    # the same case in the simulation part (data/bsc_sim/03_scenes.json)
    t = json.load(open(SIM / "data/03_scenes.json"))["metro|D0=1|samecell"]
    out["metro"]["first_counts"] = {k: t[k] for k in (
        "n_rep", "n_major", "p_major_hi", "attack_all_mean", "gen2_over_gen1", "R1_pair_ngm",
        "R1_pair_uniform_index", "R_index_alone_mean", "R_index_alone_se", "R0_ngm")}
    o = json.load(open(SIM / "data/03_scenes.json"))["office|D0=1|range1m"]
    out["office"]["full_epidemic_gen2_over_gen1"] = dict(
        first=o["gen2_over_gen1"], first_n=o["n_rep"],
        check=json.load(open(SIM / "verify/v05_scenes.json"))["office|1.0|1m"]["gen2_over_gen1"])
    RES["saturation"] = out


# --------------------------------------------------------------------------- 5b. well-mixed generation ratio
def wellmixed_generation_ratio(N=100, bq=0.5 * 1.3, n_rep=200000, seed=20261001):
    """Well-mixed Markovian SIR with pair rate bq/(N-1), one index case; the index and its direct cases
    transmit, second-generation cases do not (the design of the last column of Table s5_pair:tab:mobility).
    Returns sum(generation 2)/sum(generation 1) and the mean first generation.  Embedded jump chain."""
    rng = np.random.default_rng(seed)
    b = bq / (N - 1)
    S = np.full(n_rep, N - 1)
    I0 = np.ones(n_rep, dtype=np.int64)
    I1 = np.zeros(n_rep, dtype=np.int64)
    g1 = np.zeros(n_rep, dtype=np.int64)
    g2 = np.zeros(n_rep, dtype=np.int64)
    active = np.arange(n_rep)
    while len(active):
        s_, i0, i1 = S[active], I0[active], I1[active]
        r = np.stack([b * s_ * i0, b * s_ * i1, GAMMA * i0, GAMMA * i1], axis=1)
        tot = r.sum(axis=1)
        u = rng.random(len(active)) * tot
        c = np.cumsum(r, axis=1)
        ev = (u[:, None] >= c).sum(axis=1)
        a0, a1, a2, a3 = (active[ev == k] for k in range(4))
        S[a0] -= 1; I1[a0] += 1; g1[a0] += 1
        S[a1] -= 1; g2[a1] += 1
        I0[a2] -= 1
        I1[a3] -= 1
        active = active[(I0[active] + I1[active]) > 0]
    ratio = g2.sum() / g1.sum()
    bs = []
    for _ in range(200):
        ii = rng.integers(0, n_rep, n_rep)
        bs.append(g2[ii].sum() / g1[ii].sum())
    pair_value = (bq / GAMMA) / (1 + bq / (GAMMA * (N - 1)))
    RES["wellmixed_generation_ratio"] = dict(
        N=N, beta_q=bq, n_rep=n_rep, gen1_mean=float(g1.mean()),
        gen1_se=float(g1.std(ddof=1) / np.sqrt(n_rep)), pair_level_value=pair_value,
        gen2_over_gen1=float(ratio), lo=float(np.quantile(bs, 0.025)), hi=float(np.quantile(bs, 0.975)),
        ratio_over_pair_value=float(ratio / pair_value))


# --------------------------------------------------------------------------- 6. encounter factor
def c_factor(N, pi):
    return 1 - (1 - (1 - pi) ** N) / (N * pi)


def encounter():
    office = load_scene("office")
    out = {}
    c_off = c_factor(100, 1 / office.M)
    c_400 = c_factor(100, 1 / 400)
    c_260 = c_factor(100, 1 / 260)
    R0_off = th.R0_dense(office, BETA, GAMMA)
    out["office"] = dict(c=c_off, cells=office.M, N=100,
                         prob_index_alone_in_cell=(1 - 1 / office.M) ** 99,
                         annealed_uniform_index=c_off * BETA * office.q.mean() / GAMMA,
                         annealed_spectral=c_off * R0_off)
    out["uniform_20x20"] = dict(c=c_400, annealed=c_400 * BETA / GAMMA,
                                prob_index_alone_in_cell=(1 - 1 / 400) ** 99)
    out["uniform_20x13_q1.4"] = dict(c=c_260, annealed=c_260 * BETA * 1.4 / GAMMA)
    e = json.load(open(SIM / "data/01_validation.json"))["E_percell_rule"]
    out["counted"] = {r["case"]: dict(mean=r["index_offspring_mean"], lo=r["lo"], hi=r["hi"],
                                      n_rep=r["n_rep"],
                                      mean_final_size=r["mean_final_size"],
                                      max_final_size=r["max_final_size"]) for r in e}
    v3 = json.load(open(SIM / "verify/v03_percell_rule.json"))
    out["counted_check"] = dict(office=dict(mean=v3["office"]["off"], se=v3["office"]["se"]),
                                   uniform_20x20=dict(mean=v3["uniform20"]["off"], se=v3["uniform20"]["se"]))
    v4 = json.load(open(SIM / "verify/v04_finite.json"))
    out_f = {k: dict(exact=v["exact"], sim=v["sim"], se=v["se"], n=v["n"],
                     z=(v["sim"] - v["exact"]) / v["se"]) for k, v in v4.items()
             if "sim" in v and v["exact"] == v["exact"]}
    RES["finite_size_check"] = out_f
    v6 = json.load(open(SIM / "verify/v06_mobility.json"))
    RES["mobility_check"] = {k: dict(sim_uniform=v["sim_uniform"], se=v["sim_uniform_se"],
                                        R1_uniform=v["R1_uniform"], n=v["n"]) for k, v in v6.items()}
    out["office"]["counted_over_annealed"] = e[0]["index_offspring_mean"] / out["office"]["annealed_uniform_index"]
    out["uniform_20x20"]["counted_over_annealed"] = e[1]["index_offspring_mean"] / out["uniform_20x20"]["annealed"]
    # re-check's runs of the theory part (homogeneous 20 x 13 room, q = 1.4, 100 people)
    out["check_20x13"] = {"D0=10": dict(mean=0.825, se=0.013, n=9000, src="code/bsc_theory/verify/v4_ibm_big.log"),
                             "D0=1": dict(mean=0.754, se=0.022, src="code/bsc_theory/verify/v3_dynamics.log"),
                             "D0=0.1": dict(mean=0.545, se=0.017, src="code/bsc_theory/verify/v3_dynamics.log")}
    for k, v in out["check_20x13"].items():
        v["ratio_to_annealed"] = v["mean"] / out["uniform_20x13_q1.4"]["annealed"]
    RES["encounter"] = out


if __name__ == "__main__":
    t0 = time.time()
    finite_size()
    mobility()
    torus()
    infinite()
    encounter()
    saturation()
    wellmixed_generation_ratio()
    RES["runtime_s"] = time.time() - t0
    json.dump(RES, open(OUT / "s5_pair_checks.json", "w"), indent=1, default=float)
    open(OUT / "s5_pair_tables.tex", "w").write("\n".join(TEX) + "\n")
    print("\n".join(TEX))
    print(json.dumps({k: v for k, v in RES.items() if k != "finite_size"}, indent=1, default=float))
    print(json.dumps({k: v for k, v in RES["finite_size"].items() if k not in ("rel_dev", "occupancy_used", "ratio_to_meanfield")}, indent=1, default=float))
