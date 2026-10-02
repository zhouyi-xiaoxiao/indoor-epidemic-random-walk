#!/usr/bin/env python
"""s7_interventions_fig_hetero.py -- Figure s7_interventions:fig:hetero (two-zone room, contrast in q).

Redraws code/bsc_sim/figures/fig4_heterogeneity_* from the same data
(data/bsc_sim/05_heterogeneity.json and 05_heterogeneity_extra.json; no new simulation) with two
changes that follow from the re-check (notes/VERIFICATION.md of the repository, Section 4.2, further point 5):
  * the counted secondary cases are shown for BOTH index laws (index placed with the Perron vector of the
    pair-level matrix K1, and index placed uniformly), because the rise with contrast at N = 1000 is seen
    only for the Perron-weighted index;
  * the results at D0 = 10 are shown next to those at D0 = 1.
Labels use the article's notation (R0, rho(K1), R1-bar).
Panel (d) shows 95 % intervals of the mean attack rate, and for N = 1000 the 2,000-run
samples of ../data/s7_interventions_hetero_se.json (s7_interventions_hetero_se.py) in place of the 300-run ones.

Output: ../figures/s7_interventions_hetero_en.pdf/.png and ..._zh.pdf/.png (bilingual via figstyle).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figstyle import C, bilingual, tr  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROOT = HERE.parents[0]
MAIN = json.loads((ROOT / "data/bsc_sim/05_heterogeneity.json").read_text())
EXTRA = json.loads((ROOT / "data/bsc_sim/05_heterogeneity_extra.json").read_text())
SE = json.loads((HERE.parent / "data" / "s7_interventions_hetero_se.json").read_text())
BETA, GAMMA, QBAR = 0.5, 0.14, 1.2
C_EX = 1 - 0.5 / QBAR                        # contrast of the two-zone example (q_W = 1.5, q_C = 0.5)


def rows_main(D0, N):
    r = [v for k, v in MAIN.items() if k != "meanfield_grid" and v["D0"] == D0 and v["N"] == N]
    return sorted(r, key=lambda v: v["c"])


def rows_extra(D0, variant, with_main=True):
    r = [v for v in EXTRA.values() if v["D0"] == D0 and v["variant"] == variant]
    if variant == "zoneD" and with_main:      # negative contrasts continue the N = 100 series of MAIN
        r = r + rows_main(D0, 100)
    return sorted(r, key=lambda v: v["c"])


def panel(ax, letter, title):
    ax.set_title(f"({letter}) {title}", loc="left")


def secondary_panel(ax, D0, lang):
    cfg = [("zone", 100, C["orange"], "o"), ("uniformD", 100, C["purple"], "^"), ("zone", 1000, C["blue"], "s")]
    for variant, N, col, mk in cfg:
        if variant == "uniformD":
            r = rows_extra(D0, "uniformD")
        elif N == 100:
            r = rows_extra(D0, "zoneD")
        else:
            r = rows_main(D0, N)
        c = np.array([v["c"] for v in r])
        # Perron-weighted index: points with 95 % CI, line = rho(K1)
        ax.plot(c, [v["R1_pair_ngm"] for v in r], color=col, lw=1.3)
        y = np.array([v["sim_R_perron"] for v in r])
        err = np.array([[v["sim_R_perron"] - v["sim_R_perron_lo"], v["sim_R_perron_hi"] - v["sim_R_perron"]]
                        for v in r]).T
        ax.errorbar(c, y, yerr=err, fmt=mk, color=col, ms=4, capsize=1.5, lw=0.8)
        # uniformly placed index (not simulated for the controls): open markers, dashed line = R1-bar
        ru = [v for v in r if "sim_R_uniform" in v]
        if ru:
            cu = np.array([v["c"] for v in ru])
            ax.plot(cu, [v["R1_pair_uniform"] for v in ru], color=col, lw=1.0, ls="--")
            yu = np.array([v["sim_R_uniform"] for v in ru])
            erru = np.array([[v["sim_R_uniform"] - v["sim_R_uniform_lo"],
                              v["sim_R_uniform_hi"] - v["sim_R_uniform"]] for v in ru]).T
            ax.errorbar(cu, yu, yerr=erru, fmt=mk, color=col, mfc="white", ms=4, capsize=1.5, lw=0.8)
    ax.axvline(0.0, color=C["grey"], lw=0.7)
    ax.axvline(C_EX, color=C["black"], lw=0.7, ls=":")
    ax.set_xlim(-2.15, 1.15)
    ax.set_ylim(0, 6)
    ax.set_xlabel(tr(lang, r"contrast $c$  ($c<0$: corridor hotter)", r"对比度 $c$（$c<0$：走廊更高）"))
    ax.set_ylabel(tr(lang, "secondary cases of the index case", "指示病例的二代病例数"))


def draw(lang):
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.0))

    # (a) mean-field R0 against the contrast
    ax = axs[0, 0]
    grid = MAIN["meanfield_grid"]
    cg = np.array([g["c"] for g in grid["1.0"]])
    lo = np.full_like(cg, BETA * QBAR / GAMMA)
    hi = np.array([BETA * g["qW"] / GAMMA for g in grid["1.0"]])
    ax.fill_between(cg, lo, hi, color="0.9", zorder=0)
    for D0, col in (("0.1", C["red"]), ("1.0", C["orange"]), ("10.0", C["green"]), ("100.0", C["blue"])):
        ax.plot(cg, [g["R0"] for g in grid[D0]], color=col, lw=1.6, label=rf"$D_0={float(D0):g}$")
    ax.axvline(C_EX, color=C["black"], lw=0.7, ls=":")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(1.5, 6.4)
    ax.set_xlabel(tr(lang, r"contrast $c$  (area mean $\langle q\rangle=1.2$ fixed)",
                     r"对比度 $c$（面积平均 $\langle q\rangle=1.2$ 固定）"))
    ax.set_ylabel(tr(lang, r"mean-field $R_0=\rho(\mathbf{K})$", r"平均场 $R_0=\rho(\mathbf{K})$"))
    ax.legend(loc="upper left", fontsize=7)
    panel(ax, "a", tr(lang, "mean field", "平均场"))

    # (b), (c) counted secondary cases
    secondary_panel(axs[0, 1], 1.0, lang)
    panel(axs[0, 1], "b", tr(lang, r"discrete people, $D_0=1$", r"离散个体，$D_0=1$"))
    secondary_panel(axs[1, 0], 10.0, lang)
    panel(axs[1, 0], "c", tr(lang, r"discrete people, $D_0=10$", r"离散个体，$D_0=10$"))
    handles = [
        Line2D([], [], color=C["orange"], marker="o", ms=4, lw=1.3,
               label=tr(lang, r"$N=100$, zone mobility", r"$N=100$，分区迁移率")),
        Line2D([], [], color=C["purple"], marker="^", ms=4, lw=1.3,
               label=tr(lang, r"$N=100$, uniform mobility", r"$N=100$，均匀迁移率")),
        Line2D([], [], color=C["blue"], marker="s", ms=4, lw=1.3,
               label=tr(lang, r"$N=1000$, zone mobility", r"$N=1000$，分区迁移率")),
        Line2D([], [], color="0.3", marker="o", ms=4, lw=1.3,
               label=tr(lang, r"index placed with the Perron vector of $\mathbf{K}_1$; line $\rho(\mathbf{K}_1)$",
                        r"指示病例按 $\mathbf{K}_1$ 的 Perron 向量放置；曲线 $\rho(\mathbf{K}_1)$")),
        Line2D([], [], color="0.3", marker="o", mfc="white", ms=4, lw=1.0, ls="--",
               label=tr(lang, r"index placed uniformly; line $\bar R_1$", r"指示病例均匀放置；曲线 $\bar R_1$")),
    ]
    axs[1, 0].legend(handles=handles, loc="lower left", fontsize=6.3, handlelength=2.6)

    # (d) mean attack rate over all runs, both mobilities
    ax = axs[1, 1]
    for D0, ls, mfc_white in ((1.0, "-", False), (10.0, "--", True)):
        for variant, N, col, mk in (("zone", 100, C["orange"], "o"), ("uniformD", 100, C["purple"], "^"),
                                    ("zone", 1000, C["blue"], "s")):
            if variant == "uniformD":
                r = rows_extra(D0, "uniformD")
            elif N == 100:
                r = rows_extra(D0, "zoneD")
            else:
                r = rows_main(D0, N)
            cs = [v["c"] for v in r]
            if N == 1000:      # 2,000 further epidemics per point (s7_interventions_hetero_se.py), not the 300 stored
                rec = [SE[f"new|D0={D0:g}|N=1000|c={c:.4f}"] for c in cs]
            else:
                pre = "stored|" if variant == "zone" else "stored|extra|uniformD|"
                rec = [SE[(pre if c >= 0 or variant != "zone" else "stored|extra|zoneD|") + f"D0={D0:g}|N=100|c={c:.4f}"]
                       for c in cs]
            ax.errorbar(cs, [x["attack_all"] for x in rec], yerr=[1.96 * x["attack_all_se"] for x in rec], ls=ls,
                        marker=mk, color=col, ms=4, lw=1.1, capsize=1.5, elinewidth=0.8,
                        mfc="white" if mfc_white else col)
    ax.axvline(0.0, color=C["grey"], lw=0.7)
    ax.axvline(C_EX, color=C["black"], lw=0.7, ls=":")
    ax.set_xlim(-2.15, 1.15)
    ax.set_ylim(0, 1)
    ax.set_xlabel(tr(lang, r"contrast $c$", r"对比度 $c$"))
    ax.set_ylabel(tr(lang, "mean attack rate (all runs)", "平均罹患率（全部模拟）"))
    ax.legend(handles=[Line2D([], [], color="0.3", marker="o", ms=4, lw=1.1, label=r"$D_0=1$"),
                       Line2D([], [], color="0.3", marker="o", mfc="white", ms=4, lw=1.1, ls="--",
                              label=r"$D_0=10$")], loc="lower left", fontsize=7)
    panel(ax, "d", tr(lang, "epidemic outcome", "流行结局"))
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    print(bilingual(draw, "s7_interventions_hetero"))
