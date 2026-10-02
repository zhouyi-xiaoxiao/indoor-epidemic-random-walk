"""fig_finite_size.py -- Fig. "s5_pair:fig:finite" of the article (finite-size effect for discrete people).

Homogeneous L x L rooms with reflecting walls: mean-field R0, the pair-level closed form and exact pair solve,
the infinite-lattice value and the counted secondary cases of simulation.  No number is computed here.

Input : data/bsc_sim/02_finite_size.json   (written by code/bsc_sim/experiments/02_r0_regimes.py)
Output: ../figures/fig5_1_finite_size_en.pdf/.png and ..._zh.pdf/.png

    python fig_finite_size.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figstyle import C, bilingual, tr  # noqa: E402  (sets the Agg backend)

import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402

ROOT = HERE.parents[0]                       # repository root
D = json.loads((ROOT / "code" / "bsc_sim" / "data" / "02_finite_size.json").read_text())
BETA, GAMMA = 0.5, 0.14


def fig_finite_size(lang):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, Dm, lab in zip(axes, (1.0, 10.0), "ab"):
        rows = sorted([v for k, v in D.items() if k.startswith(f"D{Dm}_L") and "Linf" not in k],
                      key=lambda r: r["L"])
        Ls = np.array([r["L"] for r in rows])
        ax.axhline(BETA / GAMMA, color=C["black"], lw=1.2,
                   label=tr(lang, "mean field: $R_0=\\beta q/\\gamma$", "平均场：$R_0=\\beta q/\\gamma$"))
        ax.plot(Ls, [r["pair_closed_form"] for r in rows], color=C["blue"],
                label=tr(lang, "pair level, closed form for $\\bar{R}_1$", "配对层次：$\\bar{R}_1$ 的闭式"))
        ex = [(r["L"], r["pair_exact"]) for r in rows if r["pair_exact"] == r["pair_exact"]]
        ax.plot(*zip(*ex), ls="none", marker="x", ms=5, color=C["blue"],
                label=tr(lang, "pair level, exact $\\bar{R}_1$", "配对层次：精确的 $\\bar{R}_1$"))
        y = np.array([r["sim_offspring"] for r in rows])
        lo = np.array([r["sim_lo"] for r in rows])
        hi = np.array([r["sim_hi"] for r in rows])
        ax.errorbar(Ls, y, yerr=[y - lo, hi - y], ls="none", marker="o", ms=4, mfc="white",
                    color=C["orange"], capsize=2,
                    label=tr(lang, "simulation: counted $\\widehat{R}$", "模拟：计数值 $\\widehat{R}$"))
        ax.axhline(D[f"D{Dm}_Linf"]["pair_closed_form"], color=C["blue"], ls=":", lw=1.0,
                   label=tr(lang, "closed form, $L\\to\\infty$", "闭式，$L\\to\\infty$"))
        ax.set_xscale("log")
        ax.set_xticks([4, 5, 7, 10, 14, 20, 28, 40])
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xlim(3.5, 46)
        ax.set_ylim(0, 4.0)
        ax.set_xlabel(tr(lang, "room side $L$ (cells)", "房间边长 $L$（格点）"))
        ax.set_title(tr(lang, f"({lab}) $D$ = {Dm:g} cell$^2$ day$^{{-1}}$",
                        f"({lab}) $D$ = {Dm:g} 格$^2$ 天$^{{-1}}$"))
    axes[0].set_ylabel(tr(lang, "reproduction number", "再生数"))
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, -0.02), columnspacing=1.2)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    return fig


if __name__ == "__main__":
    print(bilingual(fig_finite_size, "fig5_1_finite_size"))
