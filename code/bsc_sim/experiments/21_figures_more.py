"""Further figure definitions (loaded by 20_figures.py, which injects
``register``, ``jload``, ``have`` and ``panel_label``)."""
import os

import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.colors import ListedColormap

from bsc_sim.plotting import C, tr
from bsc_sim.scenes import SCENE_NAMES, load_scene

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
BETA, GAMMA = 0.5, 0.14
SCOL = {"office": C["blue"], "supermarket": C["green"],
        "classroom": C["orange"], "metro": C["red"]}


def scene_title(name, lang):
    sc = load_scene(name)
    return sc.meta["title_en"] if lang == "en" else sc.meta["title_zh"]


def npz(name):
    return np.load(os.path.join(DATA, name), allow_pickle=True)


# =========================================================================
# Validation
# =========================================================================
@register("validation", "figA_validation", "01_validation.json",
          "01_B_finalsize.npz", "01_F_curves.npz")
def fig_validation(lang):
    v = jload("01_validation.json")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
    # (a) exact final-size law, well-mixed limit
    ax = axes[0]
    b = npz("01_B_finalsize.npz")
    k = np.arange(len(b["exact"]))
    ax.bar(k, b["empirical"], width=0.9, color=C["sky"],
           label=tr(lang, "simulation (40 000 runs)", "模拟（40 000次）"))
    ax.plot(k, b["exact"], color=C["black"], lw=1.0, marker=".", ms=3,
            label=tr(lang, "exact Markov SIR", "精确马尔可夫SIR"))
    ax.set_yscale("log")
    ax.set_ylim(2e-4, 0.5)
    ax.set_xlabel(tr(lang, "final size (one cell, $N$ = 40)",
                     "最终感染人数（单格点，$N$ = 40）"))
    ax.set_ylabel(tr(lang, "probability", "概率"))
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=7,
              borderaxespad=0.2)
    panel_label(ax, "(a)", x=-0.3, y=1.02)
    # (b) pair theory vs counted secondary cases
    ax = axes[1]
    rows = v["D_pair_theory"]
    x = np.array([r["pair_theory"] for r in rows])
    y = np.array([r["sim_mean"] for r in rows])
    lo = np.array([r["lo"] for r in rows])
    hi = np.array([r["hi"] for r in rows])
    mf = np.array([r["meanfield"] for r in rows])
    ax.plot([0, 6.2], [0, 6.2], color=C["grey"], lw=0.8)
    ax.errorbar(x, y, yerr=[y - lo, hi - y], ls="none", marker="o", ms=4,
                mfc="white", color=C["blue"], capsize=2,
                label=tr(lang, "against pair level $\\bar{R}_1$", "对比配对层次 $\\bar{R}_1$"))
    ax.plot(mf, y, ls="none", marker="s", ms=3.5, color=C["red"], alpha=0.8,
            label=tr(lang, "against mean field $\\beta\\langle q\\rangle/\\gamma$",
                     "对比平均场 $\\beta\\langle q\\rangle/\\gamma$"))
    ax.set_xlim(0, 9.5)
    ax.set_ylim(0, 6.2)
    ax.set_xlabel(tr(lang, "predicted secondary cases", "理论预测的二代病例数"))
    ax.set_ylabel(tr(lang, "counted in simulation", "模拟中计数值"))
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=7,
              borderaxespad=0.2)
    panel_label(ax, "(b)", x=-0.3, y=1.02)
    # (c) convergence to the mean-field ODE
    ax = axes[2]
    c = npz("01_F_curves.npz")
    t = c["t"]
    for N, col, D0 in ((100, C["orange"], 1), (1000, C["green"], 10),
                       (5000, C["blue"], 100)):
        ax.plot(t, c[f"N{N}_sim_mean"], color=col,
                label=tr(lang, f"$N$={N}, $D_0$={D0}", f"$N$={N}，$D_0$={D0}"))
    ax.plot(t, c["N5000_ode"], color=C["black"], ls="--", lw=1.0,
            label=tr(lang, "mean-field ODE", "平均场ODE"))
    ax.set_xlabel(tr(lang, "time (days)", "时间（天）"))
    ax.set_ylabel(tr(lang, "mean prevalence $I/N$", "平均感染比例 $I/N$"))
    ax.set_xlim(0, 60)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.05, 1.0), fontsize=6.6,
              ncol=2, columnspacing=0.7, handlelength=1.5, borderaxespad=0.2)
    panel_label(ax, "(c)", x=-0.3, y=1.02)
    fig.tight_layout()
    return fig


# =========================================================================
# R versus mobility (office)
# =========================================================================
@register("mobility", "fig5_5_R_vs_mobility", "02_mobility.json")
def fig_mobility(lang):
    d = jload("02_mobility.json")
    rows = sorted(d.values(), key=lambda r: r["D0"])
    D0 = np.array([r["D0"] for r in rows])
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    ax.fill_between(D0, rows[0]["lower_bound"], rows[0]["upper_bound"],
                    color=C["grey"], alpha=0.18, lw=0,
                    label=tr(lang, "bounds $\\beta\\langle q\\rangle/\\gamma$ – $\\beta q_{\\max}/\\gamma$",
                             "界 $\\beta\\langle q\\rangle/\\gamma$ – $\\beta q_{\\max}/\\gamma$"))
    ax.plot(D0, [r["R0_ngm"] for r in rows], color=C["black"],
            label=tr(lang, "mean-field $R_0=\\rho(K)$", "平均场 $R_0=\\rho(K)$"))
    ax.plot(D0, [r["R1_pair_ngm"] for r in rows], color=C["blue"],
            label=tr(lang, "pair level $\\rho(K_1)$, $N$ = 100",
                     "配对层次 $\\rho(K_1)$，$N$ = 100"))
    ax.plot(D0, [r["R1_pair_uniform_index"] for r in rows], color=C["blue"],
            ls="--", lw=1.1,
            label=tr(lang, "pair level $\\bar{R}_1$ (uniform index)",
                     "配对层次 $\\bar{R}_1$（均匀放置的指示病例）"))
    y = np.array([r["sim_offspring_perron"] for r in rows])
    lo = np.array([r["sim_perron_lo"] for r in rows])
    hi = np.array([r["sim_perron_hi"] for r in rows])
    ax.errorbar(D0, y, yerr=[y - lo, hi - y], ls="none", marker="o", ms=4,
                mfc="white", color=C["orange"], capsize=2,
                label=tr(lang, "simulation: $\\widehat{R}$, Perron-weighted index",
                         "模拟：$\\widehat{R}$（按Perron向量放置指示病例）"))
    y = np.array([r["sim_gen2_over_gen1"] for r in rows])
    lo = np.array([r["sim_gen2_lo"] for r in rows])
    hi = np.array([r["sim_gen2_hi"] for r in rows])
    ax.errorbar(D0, y, yerr=[y - lo, hi - y], ls="none", marker="^", ms=4,
                color=C["green"], capsize=2,
                label=tr(lang, "simulation: generation 2 / generation 1",
                         "模拟：第2代/第1代"))
    ax.axhline(1.0, color=C["red"], lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_ylim(0, 5.8)
    ax.set_xlabel(tr(lang, "mobility scale $D_0$ (cell$^2$ day$^{-1}$)",
                     "迁移率尺度 $D_0$（格$^2$ 天$^{-1}$）"))
    ax.set_ylabel(tr(lang, "reproduction number", "再生数"))
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2),
              ncol=2)
    fig.tight_layout()
    return fig


# =========================================================================
# =========================================================================
def _curve_file(name, D0):
    return f"03_curves_{name}_D{D0:g}_range1m.npz"


@register("scene_curves", "fig5_3_scene_epidemic_curves", "03_scenes.json",
          *[_curve_file(n, d) for n in SCENE_NAMES for d in (1, 10)])
def fig_scene_curves(lang):
    S = jload("03_scenes.json")
    fig, axes = plt.subplots(2, 4, figsize=(7.6, 4.2), sharey=True)
    for i, D0 in enumerate((1, 10)):
        for j, name in enumerate(SCENE_NAMES):
            ax = axes[i, j]
            c = npz(_curve_file(name, D0))
            row = S[f"{name}|D0={D0:g}|range1m"]
            t = c["t"]
            col = SCOL[name]
            if "mean_major" in c:
                ax.fill_between(t, c["q05_major"], c["q95_major"], color=col,
                                alpha=0.22, lw=0)
                ax.plot(t, c["mean_major"], color=col)
                k = int(np.argmax(c["mean_major"]))
                ax.plot(t[k], c["mean_major"][k], marker="o", ms=3.5, color=col)
            ax.plot(t, c["mean_all"], color=col, lw=0.9, ls=":")
            ax.plot(t, c["ode"], color=C["grey"], lw=1.0, ls="--")
            ax.set_xlim(0, 160 if D0 == 1 else 80)
            ax.set_ylim(0, 0.62)
            if i == 0:
                ax.set_title(scene_title(name, lang), fontsize=9)
            if i == 1:
                ax.set_xlabel(tr(lang, "time (days)", "时间（天）"))
            ax.text(0.97, 0.95,
                    f"$D_0$ = {D0}\n$P_{{\\mathrm{{maj}}}}$ = {row['p_major']:.2f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
        axes[i, 0].set_ylabel(tr(lang, "prevalence $I/N$", "感染比例 $I/N$"))
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h = [Line2D([], [], color="k", label=tr(lang, "mean of major outbreaks",
                                             "大规模暴发的均值")),
         Patch(color="k", alpha=0.22, label=tr(lang, "5–95 % band (major)",
                                               "5–95%分位带（大规模暴发）")),
         Line2D([], [], color="k", ls=":", lw=0.9,
                label=tr(lang, "mean of all runs", "全部模拟的均值")),
         Line2D([], [], color=C["grey"], ls="--",
                label=tr(lang, "mean-field ODE", "平均场ODE"))]
    fig.legend(handles=h, loc="lower center", ncol=4, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return fig


@register("scene_attack", "fig5_4_attack_rate_distributions", "03_scenes.json",
          *[_curve_file(n, d) for n in SCENE_NAMES for d in (1, 10)])
def fig_scene_attack(lang):
    fig, axes = plt.subplots(2, 4, figsize=(7.6, 3.6), sharex=True, sharey=True)
    bins = np.linspace(0, 1, 26)
    for i, D0 in enumerate((1, 10)):
        for j, name in enumerate(SCENE_NAMES):
            ax = axes[i, j]
            c = npz(_curve_file(name, D0))
            ar = c["attack"]
            ax.hist(ar, bins=bins, weights=np.ones(len(ar)) / len(ar),
                    color=SCOL[name], alpha=0.85)
            ax.axvline(0.10, color=C["black"], lw=0.8, ls="--")
            ax.set_yscale("log")
            ax.set_ylim(3e-4, 1)
            if i == 0:
                ax.set_title(scene_title(name, lang), fontsize=9)
            if i == 1:
                ax.set_xlabel(tr(lang, "attack rate", "罹患率"))
            ax.text(0.5, 0.93, f"$D_0$ = {D0}", transform=ax.transAxes,
                    ha="center", va="top", fontsize=7.5)
        axes[i, 0].set_ylabel(tr(lang, "fraction of runs", "模拟次数占比"))
    fig.tight_layout()
    return fig


# =========================================================================
# =========================================================================
def _spatial(lang, name, D0, gen=2):
    J = jload("04_spatial.json")[f"{name}|D0={D0:g}"][f"gen{gen}"]
    z = npz(f"04_spatial_{name}_D{D0:g}.npz")
    sc = load_scene(name, D0=D0)
    M = sc.M
    pred = z[f"pred_pair_gen{gen}"] * M
    pmf = z[f"pred_mf_gen{gen}"] * M
    obs = z[f"counts_gen{gen}"] / z[f"counts_gen{gen}"].sum() * M
    wide = sc.nx / sc.ny > 4
    if wide:
        fig = plt.figure(figsize=(7.4, 5.4))
        gs = fig.add_gridspec(5, 1, height_ratios=[0.5, 0.5, 0.07, 0.35, 1.5],
                              hspace=0.5)
        axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        sub0 = gs[2].subgridspec(1, 3, width_ratios=[0.3, 1, 0.3])
        cax = fig.add_subplot(sub0[0, 1])
        sub = gs[4].subgridspec(1, 3, width_ratios=[0.6, 1, 0.6])
        axc = fig.add_subplot(sub[0, 1])
    else:
        fig = plt.figure(figsize=(7.6, 2.9))
        gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.22, 0.8], wspace=0.1)
        axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
        axc = fig.add_subplot(gs[0, 3])
    vmax = max(pred.max(), obs.max())
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#bbbbbb")
    for ax, v, ttl in ((axa, pred, tr(lang, "(a) theory (pair level)",
                                      "(a) 理论（配对层次）")),
                       (axb, obs, tr(lang, "(b) simulation", "(b) 模拟"))):
        im = ax.imshow(np.ma.masked_invalid(sc.to_grid(v)), cmap=cmap, vmin=0,
                       vmax=vmax, origin="upper", interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp_ in ax.spines.values():
            sp_.set_visible(True)
        ax.set_title(ttl, fontsize=9)
    if wide:
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    else:
        cb = fig.colorbar(im, ax=[axa, axb], orientation="horizontal",
                          fraction=0.07, pad=0.08, aspect=40)
    cb.set_label(tr(lang, "relative infection frequency (1 = uniform)",
                    "相对感染频率（1 = 均匀）"), fontsize=8)
    lim = 1.05 * max(vmax, pmf.max())
    axc.plot([0, lim], [0, lim], color=C["grey"], lw=0.8)
    axc.plot(pmf, obs, ls="none", marker="s", ms=2.5, color=C["red"],
             alpha=0.55,
             label=tr(lang, f"mean field, $r$ = {J['r_meanfield']:.2f}",
                      f"平均场，$r$ = {J['r_meanfield']:.2f}"))
    axc.plot(pred, obs, ls="none", marker="o", ms=2.8, color=C["blue"],
             alpha=0.75,
             label=tr(lang, f"pair level, $r$ = {J['r_pair']:.2f}",
                      f"配对层次，$r$ = {J['r_pair']:.2f}"))
    axc.set_xlim(0, lim)
    axc.set_ylim(0, lim)
    axc.set_xlabel(tr(lang, "predicted", "理论预测"))
    axc.set_ylabel(tr(lang, "simulated", "模拟结果"))
    axc.set_title(tr(lang, "(c) cell by cell", "(c) 逐格点比较"), fontsize=9)
    axc.legend(fontsize=7, loc="upper left", handletextpad=0.2)
    return fig


for _n in SCENE_NAMES:
    for _d in (1, 10):
        _key = f"spatial_{_n}_D{_d}"
        _fname = ("fig5_2_spatial_pattern" if (_n == "office" and _d == 1)
                  else f"figS_spatial_pattern_{_n}_D{_d}")
        register(_key, _fname, "04_spatial.json",
                 f"04_spatial_{_n}_D{_d:g}.npz")(
            lambda lang, n=_n, d=_d: _spatial(lang, n, d))


# =========================================================================
# Heterogeneity raises R0
# =========================================================================
@register("heterogeneity", "fig4_heterogeneity", "05_heterogeneity.json",
          "05_heterogeneity_extra.json")
def fig_heterogeneity(lang):
    d = jload("05_heterogeneity.json")
    e = jload("05_heterogeneity_extra.json")
    c_th = 1 - 0.5 / 1.2
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.5))
    # ---------------- (a) mean field
    ax = axes[0]
    g = d["meanfield_grid"]
    cs = np.array([r["c"] for r in g["1.0"]])
    qW = np.array([r["qW"] for r in g["1.0"]])
    ax.fill_between(cs, BETA * 1.2 / GAMMA, BETA * qW / GAMMA, color=C["grey"],
                    alpha=0.18, lw=0)
    for D0, col in (("0.1", C["red"]), ("1.0", C["orange"]),
                    ("10.0", C["green"]), ("100.0", C["blue"])):
        ax.plot(cs, [r["R0"] for r in g[D0]], color=col,
                label=f"$D_0$ = {float(D0):g}")
    ax.axvline(c_th, color=C["black"], lw=0.5, ls="--")
    ax.set_ylim(1.5, 6.4)
    ax.set_xlabel(tr(lang, "contrast $c$ ($\\langle q\\rangle$ fixed)",
                     "反差 $c$（$\\langle q\\rangle$ 不变）"))
    ax.set_ylabel(tr(lang, "mean-field $R_0=\\rho(K)$", "平均场 $R_0=\\rho(K)$"))
    ax.set_title(tr(lang, "mean field", "平均场"), fontsize=8.5)
    ax.legend(fontsize=6.3, loc="upper left")
    panel_label(ax, "(a)", x=-0.22)

    def series(src, prefix, D0, N):
        rows = [v for k, v in src.items() if k.startswith(prefix)
                and v.get("D0") == D0 and v.get("N") == N]
        return sorted(rows, key=lambda r: r["c"])

    main_1 = series(d, "D0=1|N=100|", 1.0, 100)
    neg_1 = series(e, "extra|zoneD|D0=1|", 1.0, 100)
    uni_1 = series(e, "extra|uniformD|D0=1|", 1.0, 100)
    big_1 = series(d, "D0=1|N=1000|", 1.0, 1000)
    sets = [
        (sorted(neg_1 + main_1, key=lambda r: r["c"]), C["orange"], "o",
         tr(lang, "$N$=100, two-zone mobility", "$N$=100，两区扩散")),
        (uni_1, C["purple"], "^",
         tr(lang, "$N$=100, uniform mobility", "$N$=100，均匀扩散")),
        (big_1, C["blue"], "s",
         tr(lang, "$N$=1000, two-zone mobility", "$N$=1000，两区扩散")),
    ]
    # ---------------- (b) discrete people: reproduction number
    ax = axes[1]
    for rows, col, mk, lab in sets:
        if not rows:
            continue
        c = np.array([r["c"] for r in rows])
        ax.plot(c, [r["R1_pair_ngm"] for r in rows], color=col, lw=1.1)
        y = np.array([r["sim_R_perron"] for r in rows])
        lo = np.array([r["sim_R_perron_lo"] for r in rows])
        hi = np.array([r["sim_R_perron_hi"] for r in rows])
        ax.errorbar(c, y, yerr=[y - lo, hi - y], ls="none", marker=mk, ms=3.5,
                    mfc="white", color=col, capsize=1.5, label=lab)
    ax.axvline(c_th, color=C["black"], lw=0.5, ls="--")
    ax.axvline(0, color=C["grey"], lw=0.5)
    ax.set_xlabel(tr(lang, "contrast $c$  ($c<0$: corridor hotter)",
                     "反差 $c$（$c<0$：走廊更高）"))
    ax.set_ylabel(tr(lang, "secondary cases (lines: $\\rho(K_1)$)",
                     "二代病例数（线：$\\rho(K_1)$）"))
    ax.set_title(tr(lang, "discrete people, $D_0$ = 1", "离散个体，$D_0$ = 1"),
                 fontsize=8.5)
    ax.set_ylim(0, 6)
    ax.legend(fontsize=6.0, loc="lower left")
    panel_label(ax, "(b)", x=-0.2)
    # ---------------- (c) epidemic outcome
    ax = axes[2]
    for rows, col, mk, lab in sets:
        if not rows:
            continue
        c = np.array([r["c"] for r in rows])
        ax.plot(c, [r["attack_all"] for r in rows], marker=mk, ms=3.5,
                mfc="white", color=col, lw=1.1)
    ax.axvline(c_th, color=C["black"], lw=0.5, ls="--")
    ax.axvline(0, color=C["grey"], lw=0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel(tr(lang, "contrast $c$", "反差 $c$"))
    ax.set_ylabel(tr(lang, "mean attack rate (all runs)", "平均罹患率（全部模拟）"))
    ax.set_title(tr(lang, "epidemic outcome, $D_0$ = 1", "疫情结局，$D_0$ = 1"),
                 fontsize=8.5)
    panel_label(ax, "(c)", x=-0.2)
    fig.tight_layout()
    return fig


# =========================================================================
# Obstacles
# =========================================================================
FRACS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]


def _obst_agg(d, cal, field):
    """Mean over layouts and 95 % CI half-width (t, n-1 dof) per density."""
    from scipy import stats
    m, h = [], []
    for f in FRACS:
        v = np.array([r[field] for k, r in d.items()
                      if r["rhoB_nominal"] == f and r["calibration"] == cal
                      and r.get(field) is not None], float)
        v = v[np.isfinite(v)]
        m.append(v.mean())
        h.append(stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
                 if len(v) > 1 else 0.0)
    return np.array(m), np.array(h)


@register("obstacle_scan", "fig6_3_obstacle_scan", "06_obstacles.json")
def fig_obstacle_scan(lang):
    d = jload("06_obstacles.json")
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2))
    x = np.array(FRACS)
    lab = {"FD": "FD", "DD": "DD"}
    lab_long = {"FD": tr(lang, "FD (frequency-dependent)", "FD（频率依赖）"),
                "DD": tr(lang, "DD (density-dependent)", "DD（密度依赖）")}
    col = {"FD": C["blue"], "DD": C["red"]}
    ax = axes[0, 0]
    for cal in ("FD", "DD"):
        m, _ = _obst_agg(d, cal, "R0_ngm")
        ax.plot(x, m, color=col[cal], ls="--", lw=1.0,
                label=tr(lang, f"mean-field $R_0$, {lab[cal]}",
                         f"平均场 $R_0$，{lab[cal]}"))
        m, _ = _obst_agg(d, cal, "R1_pair_centre")
        ax.plot(x, m, color=col[cal], lw=1.2,
                label=tr(lang, f"pair level $R_1(x_c)$, {lab[cal]}",
                         f"配对层次 $R_1(x_c)$，{lab[cal]}"))
        m, h = _obst_agg(d, cal, "R_index_alone")
        ax.errorbar(x, m, yerr=h, ls="none", marker="o", ms=3.5, mfc="white",
                    color=col[cal], capsize=2)
    ax.set_ylabel(tr(lang, "reproduction number", "再生数"))
    ax.set_ylim(0, 9.5)
    ax.legend(fontsize=6.8, loc="upper left", ncol=2, columnspacing=1.0,
              handlelength=1.8)
    panel_label(ax, "(a)", x=0.0, y=1.02)
    for ax, field, yl, pl in (
            (axes[0, 1], "peak_prev_major_mean",
             tr(lang, "peak prevalence (major)", "峰值感染比例（大规模暴发）"), "(b)"),
            (axes[1, 0], "peak_time_major_mean",
             tr(lang, "peak day (major)", "峰值时间（天，大规模暴发）"), "(c)"),
            (axes[1, 1], "attack_all_mean",
             tr(lang, "mean attack rate (all runs)", "平均罹患率（全部模拟）"), "(d)")):
        for cal in ("FD", "DD"):
            m, h = _obst_agg(d, cal, field)
            ax.errorbar(x, m, yerr=h, marker="o", ms=3.5, color=col[cal],
                        capsize=2, lw=1.1, label=lab_long[cal])
        ax.set_ylabel(yl)
        ax.set_ylim(0, None)
        panel_label(ax, pl, x=0.0, y=1.02)
    axes[0, 1].legend(fontsize=7, loc="lower left")
    for ax in axes[1]:
        ax.set_xlabel(tr(lang, "nominal barrier density $\\rho_B$", "名义障碍物密度 $\\rho_B$"))
    fig.tight_layout()
    return fig


@register("obstacle_curves", "fig6_1_obstacle_curves", "06_obstacle_curves.npz",
          "06_obstacles_runs.npz")
def fig_obstacle_curves(lang):
    z = npz("06_obstacle_curves.npz")
    runs = npz("06_obstacles_runs.npz")
    t = z["t"]
    N = 100
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9),
                             gridspec_kw=dict(width_ratios=[1.6, 1]))
    for f, col, lab in ((0.0, C["blue"], tr(lang, "no barriers", "无障碍物")),
                        (0.10, C["red"], tr(lang, "10 % barriers", "10%障碍物"))):
        keys = [k for k in z.files if k.endswith("|I") and f"rhoB={f:.2f}|FD" in k]
        I = np.concatenate([z[k] for k in keys])
        att = np.concatenate([runs[k[:-2] + "|attack"] for k in keys])
        maj = att >= 0.10
        ax = axes[0]
        ax.plot(t, I[maj].mean(axis=0) / N, color=col,
                label=tr(lang, f"{lab}: major outbreaks ({maj.mean():.0%} of runs)",
                         f"{lab}：大规模暴发（占{maj.mean():.0%}）"))
        q05, q95 = np.quantile(I[maj] / N, [0.05, 0.95], axis=0)
        ax.fill_between(t, q05, q95, color=col, alpha=0.13, lw=0)
        ax.plot(t, I.mean(axis=0) / N, color=col, ls=":", lw=1.1,
                label=tr(lang, f"{lab}: all runs", f"{lab}：全部模拟"))
        axes[1].plot(t, I.mean(axis=0), color=col, label=lab)
    ax = axes[0]
    ax.set_xlim(0, 110)
    ax.set_ylim(0, None)
    ax.set_xlabel(tr(lang, "time (days)", "时间（天）"))
    ax.set_ylabel(tr(lang, "prevalence $I/N$", "感染比例 $I/N$"))
    ax.legend(fontsize=7)
    panel_label(ax, "(a)", x=-0.12)
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.set_xlabel(tr(lang, "time (days)", "时间（天）"))
    ax.set_ylabel(tr(lang, "mean number of infectives (all runs)",
                     "平均感染者人数（全部模拟）"))
    ax.set_title(tr(lang, "identical start: $I(0)$ = 1", "相同初始条件：$I(0)$ = 1"),
                 fontsize=8)
    ax.legend(fontsize=7)
    panel_label(ax, "(b)", x=-0.2)
    fig.tight_layout()
    return fig


@register("obstacle_snapshots", "fig6_2_spatial_spread", "06_snapshots.npz")
def fig_obstacle_snapshots(lang):
    z = npz("06_snapshots.npz")
    times = z["times"]
    fig, axes = plt.subplots(2, 4, figsize=(7.6, 3.0))
    vmax = max(np.nanmax(z["rhoB0.00_density"]), np.nanmax(z["rhoB0.10_density"]))
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#9a9a9a")
    for i, f in enumerate((0.0, 0.10)):
        dens = z[f"rhoB{f:.2f}_density"]
        for j in range(4):
            ax = axes[i, j]
            im = ax.imshow(np.ma.masked_invalid(dens[j]), cmap=cmap, vmin=0,
                           vmax=vmax, origin="upper", interpolation="nearest")
            ax.plot([10], [6], marker="+", color="white", ms=5, mew=0.8)
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(True)
            if i == 0:
                ax.set_title(tr(lang, f"$t$ = {times[j]} d", f"$t$ = {times[j]} 天"),
                             fontsize=9)
        axes[i, 0].set_ylabel(tr(lang, "no barriers" if f == 0 else "10 % barriers",
                                 "无障碍物" if f == 0 else "10%障碍物"), fontsize=8.5)
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label(tr(lang, "mean infectives per cell", "每格点平均感染者数"),
                 fontsize=8)
    return fig


# =========================================================================
# Occupancy
# =========================================================================
@register("occupancy", "fig6_4_occupancy", "07_occupancy.json")
def fig_occupancy(lang):
    from matplotlib.lines import Line2D
    d = jload("07_occupancy.json")
    fig, axes = plt.subplots(2, 3, figsize=(7.6, 5.0), sharex=True)
    lab = {"FD": tr(lang, "frequency-dependent", "频率依赖"),
           "DD": tr(lang, "density-dependent", "密度依赖")}
    col = {"FD": C["blue"], "DD": C["red"]}
    for i, D0 in enumerate((1, 10)):
        for cal in ("FD", "DD"):
            rows = sorted([r for k, r in d.items()
                           if k.startswith(f"D0={D0}|") and k.endswith(cal)],
                          key=lambda r: r["N"])
            if not rows:
                continue
            N = np.array([r["N"] for r in rows])
            ax = axes[i, 0]
            ax.plot(N, [r["R0_ngm"] for r in rows], color=col[cal], ls="--",
                    lw=1.0)
            ax.plot(N, [r["R1_pair_uniform"] for r in rows], color=col[cal], lw=1.2)
            y = np.array([r["R_index_alone"] for r in rows])
            e = 1.96 * np.array([r["R_index_alone_se"] for r in rows])
            ax.errorbar(N, y, yerr=e, ls="none", marker="o", ms=3.5, mfc="white",
                        color=col[cal], capsize=1.5)
            # the 10 % threshold is meaningless for very small groups
            big = N >= 50
            ax = axes[i, 1]
            ax.plot(N[big], np.array([r["p_major_branching"] for r in rows])[big],
                    color=col[cal], ls="--", lw=1.0)
            y = np.array([r["p_major"] for r in rows])
            lo = np.array([r["p_major_lo"] for r in rows])
            hi = np.array([r["p_major_hi"] for r in rows])
            ax.errorbar(N[big], y[big], yerr=[(y - lo)[big], (hi - y)[big]],
                        ls="none", marker="o", ms=3.5, color=col[cal], capsize=1.5)
            ax = axes[i, 2]
            ax.plot(N[big], np.array([r["ode_peak_prev"] for r in rows])[big],
                    color=col[cal], ls="--", lw=1.0)
            y = np.array([r.get("peak_prev_major_mean", np.nan) for r in rows], float)
            lo = np.array([r.get("peak_prev_major_lo", np.nan) for r in rows], float)
            hi = np.array([r.get("peak_prev_major_hi", np.nan) for r in rows], float)
            ax.errorbar(N[big], y[big], yerr=[(y - lo)[big], (hi - y)[big]],
                        ls="none", marker="o", ms=3.5, color=col[cal], capsize=1.5)
        axes[i, 0].set_ylabel(tr(lang, f"$D_0$ = {D0}\nreproduction number",
                                 f"$D_0$ = {D0}\n再生数"))
        axes[i, 0].set_ylim(0, 12)
        axes[i, 1].set_ylabel(tr(lang, "$P$(major outbreak)", "大规模暴发概率"))
        axes[i, 1].set_ylim(0, 1)
        axes[i, 2].set_ylabel(tr(lang, "peak prevalence (major)", "峰值感染比例（大规模暴发）"))
        axes[i, 2].set_ylim(0, 0.85)
    for ax in axes[1]:
        ax.set_xlabel(tr(lang, "people in the office $N$", "办公室人数 $N$"))
    for ax, s_ in zip(axes.ravel(), "abcdef"):
        panel_label(ax, f"({s_})", x=-0.22)
    h = [Line2D([], [], color=col["FD"], lw=2, label=lab["FD"]),
         Line2D([], [], color=col["DD"], lw=2, label=lab["DD"]),
         Line2D([], [], color="k", ls="--", lw=1.0,
                label=tr(lang, "mean field ($R_0$ / branching / ODE)",
                         "平均场（$R_0$/分支过程/ODE）")),
         Line2D([], [], color="k", lw=1.2,
                label=tr(lang, "pair level $\\bar{R}_1$", "配对层次 $\\bar{R}_1$")),
         Line2D([], [], color="k", ls="none", marker="o", ms=3.5,
                label=tr(lang, "simulation (95 % CI)", "模拟（95%置信区间）"))]
    fig.legend(handles=h, loc="lower center", ncol=5, fontsize=7,
               bbox_to_anchor=(0.5, -0.005), columnspacing=1.0)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    return fig


# =========================================================================
# Metro, time-dependent parameters
# =========================================================================
@register("metro_time", "fig5_6_metro_time_dependence", "09_metro_time.json",
          "09_metro_curves.npz")
def fig_metro_time(lang):
    J = jload("09_metro_time.json")
    z = npz("09_metro_curves.npz")
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.8))
    ax = axes[0]
    h = np.arange(25)
    rho = np.r_[z["rho"], z["rho"][-1]]
    ax.step(h, rho, where="post", color=C["blue"])
    ax.set_xlabel(tr(lang, "hour of day", "一天中的时刻（时）"))
    ax.set_ylabel(tr(lang, "crowd density $\\rho(t)$ (persons/m$^2$)",
                     "人员密度 $\\rho(t)$（人/m$^2$）"), color=C["blue"])
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_ylim(0, 5)
    ax2 = ax.twinx()
    ax2.spines["right"].set_visible(True)
    ax2.step(h, np.r_[z["m"], z["m"][-1]], where="post", color=C["red"], lw=1.1)
    ax2.set_ylim(0, 1.25)
    ax2.set_ylabel(tr(lang, "multiplier $m(t)$", "调制因子 $m(t)$"), color=C["red"])
    panel_label(ax, "(a)", x=-0.2)
    names = {"peak": tr(lang, "constant, peak values", "恒定（高峰值）"),
             "varying": tr(lang, "time-varying", "随时间变化"),
             "mean": tr(lang, "constant, time averages", "恒定（时间平均值）")}
    cols = {"peak": C["blue"], "varying": C["red"], "mean": C["black"]}
    lss = {"peak": "-", "varying": "-", "mean": "--"}
    for ax, D0, pl in ((axes[1], 1, "(b)"), (axes[2], 100, "(c)")):
        for arm in ("peak", "varying", "mean"):
            key = f"D0={D0}|{arm}"
            if f"{key}|mean_major" not in z.files:
                continue
            t = z[f"{key}|t"]
            ax.plot(t, z[f"{key}|mean_major"], color=cols[arm], ls=lss[arm],
                    label=names[arm])
            if arm != "mean":
                ax.fill_between(t, z[f"{key}|q05"], z[f"{key}|q95"],
                                color=cols[arm], alpha=0.12, lw=0)
        ax.set_xlim(0, 90)
        ax.set_ylim(0, None)
        ax.set_xlabel(tr(lang, "time (days)", "时间（天）"))
        ax.set_ylabel(tr(lang, "prevalence $I/N$ (major)", "感染比例 $I/N$（大规模暴发）"))
        ax.set_title(f"$D_0$ = {D0}", fontsize=8.5)
        panel_label(ax, pl, x=-0.2)
    axes[1].legend(fontsize=6.8)
    fig.tight_layout()
    return fig


# =========================================================================
# Benchmark
# =========================================================================
@register("benchmark", "fig4_benchmark_theory_vs_mc", "10_benchmark.json")
def fig_benchmark(lang):
    J = jload("10_benchmark.json")
    tasks = [
        ("T1_R0_meanfield", "cpu_sparse_s", "mc_cpu_for_1pct_s",
         tr(lang, "mean-field $R_0$", "平均场 $R_0$")),
        ("T2_column_sum", "cpu_solve_s", "mc_cpu_for_1pct_s",
         tr(lang, "cases per cell $(1^TK)_x$", "分格点病例数 $(1^TK)_x$")),
        ("T3_R1_discrete", "cpu_pair_s", "mc_cpu_for_1pct_s",
         tr(lang, "pair-level $R_1$", "配对层次 $R_1$")),
        ("T4_p_major", "cpu_branching_s", "mc_cpu_for_se_0p01_s",
         tr(lang, "$P$(major outbreak)", "大规模暴发概率")),
        ("T5_attack_rate", "cpu_ode_s", "mc_cpu_for_1pct_s",
         tr(lang, "attack rate", "罹患率")),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), sharey=True)
    for ax, key in zip(axes, ("D0=1", "D0=10")):
        R = J[key]
        x = np.arange(len(tasks))
        th_t = [R[t][a] for t, a, b, _ in tasks]
        mc_t = [R[t][b] for t, a, b, _ in tasks]
        ax.bar(x - 0.2, th_t, 0.4, color=C["blue"],
               label=tr(lang, "linear algebra / ODE", "线性代数/ODE"))
        ax.bar(x + 0.2, mc_t, 0.4, color=C["orange"],
               label=tr(lang, "Monte Carlo, same accuracy (1 %)",
                        "蒙特卡罗，相同精度（1%）"))
        biased = [False, False, False,
                  abs(R["T4_p_major"]["theory_bias"]) > 0.01,
                  abs(R["T5_attack_rate"]["theory_bias"])
                  > 0.01 * R["T5_attack_rate"]["mc_estimate"]]
        for xi, a, b, bad in zip(x, th_t, mc_t, biased):
            r_ = b / a
            lab_ = f"×{r_:,.0f}" if r_ >= 10 else f"×{r_:.2g}"
            ax.text(xi, max(a, b) * 1.6, lab_ + ("*" if bad else ""),
                    ha="center", fontsize=7)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels([t[3] for t in tasks], rotation=25, ha="right",
                           fontsize=7.5)
        ax.set_title(tr(lang, f"office, $D_0$ = {key[3:]}", f"办公室，$D_0$ = {key[3:]}"),
                     fontsize=9)
        ax.set_ylim(1e-4, 3e3)
    axes[0].set_ylabel(tr(lang, "single-core CPU time (s)", "单核CPU时间（秒）"))
    axes[0].legend(fontsize=7, loc="upper left")
    fig.text(0.5, 0.0, tr(lang,
             "* the deterministic (mean-field) value is off by more than the accuracy target",
             "* 确定性（平均场）结果的偏差超过精度目标"), ha="center", fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig


# =========================================================================
# When is the mean-field NGM accurate?  R1 / R0 for a homogeneous floor
# =========================================================================
@register("validity_map", "fig4_meanfield_validity_map")
def fig_validity_map(lang):
    """Closed form on the infinite lattice, same-cell rule:
    R1/R0 = 1 / (1 + (beta q / rho) g0),  g0 = (2/pi) K(k) / (gamma + 8D),
    k = 8D / (gamma + 8D)  (lattice Green's function of the relative walk)."""
    from scipy.special import ellipk
    D = np.geomspace(0.01, 1000, 240)
    rho = np.geomspace(0.05, 30, 200)
    k = 8 * D / (GAMMA + 8 * D)
    g0 = (2 / np.pi) * ellipk(k ** 2) / (GAMMA + 8 * D)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2), sharey=True)
    for ax, bq, lab in zip(axes, (0.5, 1.25), "ab"):
        ratio = 1.0 / (1.0 + bq / rho[:, None] * g0[None, :])
        cs = ax.contourf(D, rho, ratio, levels=np.linspace(0, 1, 11),
                         cmap="viridis")
        cl = ax.contour(D, rho, ratio, levels=[0.5, 0.9], colors="white",
                        linewidths=1.0, linestyles=["--", "-"])
        ax.clabel(cl, fmt="%.1f", fontsize=7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(tr(lang, "mobility $D$ (cell$^2$ day$^{-1}$)", "迁移率 $D$（格$^2$ 天$^{-1}$）"))
        ax.set_title(tr(lang, f"({lab}) $\\beta q$ = {bq:g} day$^{{-1}}$", f"({lab}) $\\beta q$ = {bq:g} 天$^{{-1}}$"),
                     fontsize=9)
        for name, occ, dlo, dhi in (("office", 0.43, 0.3, 1.5),
                                    ("supermarket", 0.54, 0.1, 2.0)):
            ax.plot([dlo, dhi], [occ, occ], color="white", lw=2.2,
                    solid_capstyle="butt")
            ax.plot([dlo, dhi], [occ, occ], color=SCOL[name], lw=1.2)
            ax.text(dhi * 1.3, occ, scene_title(name, lang), color="white",
                    fontsize=6.5, va="center")
    axes[0].set_ylabel(tr(lang, "occupancy $\\bar\\rho$ (people per cell)",
                          "占据率 $\\bar\\rho$（人/格点）"))
    cb = fig.colorbar(cs, ax=axes, pad=0.02)
    cb.set_label(tr(lang, "$R_1/R_0$ (pair level / mean field)",
                    "$R_1/R_0$（配对层次/平均场）"), fontsize=8)
    return fig
