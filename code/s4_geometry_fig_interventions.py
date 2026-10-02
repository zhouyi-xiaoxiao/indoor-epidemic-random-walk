"""s4_geometry_fig_interventions.py -- Fig. "s4_geometry:fig:interventions" of the article.

Mean-field worked examples on the canonical office scene (code/bsc_sim/scenes/office.json):
  (a) corridor mobility multiplied by a factor (symmetric rule and a departure-site rule),
  (b) walls on workstation-workstation edges and 39 workstation cells turned into barrier cells, at
      D0 = 1, 10 and 100 cell^2/day,
  (c) targeted ventilation (q -> q/2 on a fraction of the floor) at D0 = 1 cell^2/day,
  (d) the same at D0 = 100 cell^2/day.

No number is computed here.

Input : ../data/plan_theory_checks.json        (written by plan_theory_checks.py)
        ../data/plan_corridor_sweep.json       (written by plan_fig_theory.py; 25 multipliers)
        ../data/s4_geometry_barriers_leak.json (written by s4_geometry_barriers_leak.py; panel (b))
Output: ../figures/fig_s4_interventions_en.pdf/.png and ..._zh.pdf/.png

    python s4_geometry_fig_interventions.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figstyle import C, bilingual, tr  # noqa: E402  (sets the Agg backend)

import matplotlib.pyplot as plt  # noqa: E402

DATA = HERE.parent / "data"
RES = json.loads((DATA / "plan_theory_checks.json").read_text())
CORR = json.loads((DATA / "plan_corridor_sweep.json").read_text())
BL = json.loads((DATA / "s4_geometry_barriers_leak.json").read_text())
OFF = RES["office"]


def panel(ax, letter, title):
    ax.set_title(f"({letter}) {title}", loc="left")


def draw(lang):
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.6))

    # ---------------------------------------------------------------- (a) corridor mobility
    ax = axs[0, 0]
    ax.semilogx(CORR["mult"], CORR["symmetric"], color=C["blue"], lw=1.6,
                label=tr(lang, "symmetric rule", "对称跳跃规则"))
    ax.semilogx(CORR["mult"], CORR["departure"], color=C["orange"], lw=1.4, ls="--",
                label=tr(lang, "departure-site rule", "出发点规则"))
    ax.axhline(OFF["upper_bound"], color=C["black"], ls=":", lw=0.8)
    ax.text(9.5, OFF["upper_bound"] + 0.008, r"$\beta q_{\max}/\gamma$", ha="right", va="bottom", fontsize=7.5)
    ax.axvline(1.0, color=C["grey"], lw=0.8, ls="-.")
    ax.set_xlabel(tr(lang, "corridor mobility multiplier", "走廊迁移率倍数"))
    ax.set_ylabel(r"$R_0$")
    ax.set_ylim(4.95, 5.42)
    ax.legend(loc="lower left", fontsize=7)
    panel(ax, "a", tr(lang, "corridor mobility, $D_0=1$", "走廊迁移率，$D_0=1$"))

    # ---------------------------------------------------------------- (b) walls and barrier cells
    # Walls keep the cell set (edge-wise monotonicity applies: R0 cannot fall); barrier cells remove cells and
    # with them their q, so the sign of the change depends on the mobility scale.
    ax = axs[0, 1]
    A = BL["A_partitions"]
    keys = ("walls_10pct", "walls_30pct", "walls_50pct", "walls_100pct")
    xs = [0] + [100 * A["walls"][k]["mean_closed_edges"] / A["n_WW_edges"] for k in keys]
    bc = A["barrier_cells_39_workstations"]
    for D0, col, mk in (("1", C["blue"], "o"), ("10", C["orange"], "s"), ("100", C["green"], "^")):
        base = A["baseline_R0"][f"D0={D0}"]
        ys = [0] + [A["walls"][k][f"D0={D0}"]["change_pct"] for k in keys]
        es = [0] + [100 * A["walls"][k][f"D0={D0}"]["sd"] / base for k in keys]
        ax.errorbar(xs, ys, yerr=es, fmt=mk + "-", color=col, ms=3.5, capsize=2, lw=1.2,
                    label=tr(lang, f"walls, $D_0={D0}$", f"隔断墙，$D_0={D0}$"))
        ax.axhline(bc[f"D0={D0}"]["change_pct"], color=col, lw=1.0, ls="--")
        ax.text(61, bc[f"D0={D0}"]["change_pct"] + 0.04, tr(lang, f"39 barrier cells, $D_0={D0}$",
                                                           f"39 个障碍格，$D_0={D0}$"),
                ha="right", va="bottom", fontsize=6.5, color=col)
    ax.axhline(0, color=C["black"], lw=0.8)
    ax.set_xlim(-2, 62)
    ax.set_ylim(-1.75, 1.05)
    ax.set_xlabel(tr(lang, "workstation-workstation edges closed (%)", "被封闭的工位间连边（%）"))
    ax.set_ylabel(tr(lang, r"change of $R_0$ (%)", r"$R_0$ 的变化（%）"))
    ax.legend(loc="upper left", fontsize=6.5, ncol=3, columnspacing=0.8, handlelength=1.4, borderaxespad=0.2)
    panel(ax, "b", tr(lang, "walls and barrier cells", "隔断墙与障碍格"))

    # ---------------------------------------------------------------- (c), (d) targeted ventilation
    def targeting(ax, key, letter, legend):
        T = OFF["targeted_ventilation"][key]
        rows = T["rows"]
        b = T["baseline"]
        fr = [0] + [100 * r["fraction"] for r in rows]
        ax.plot(fr, [b] + [r["one_shot"] for r in rows], "-.", color=C["grey"], lw=1.1, marker="d", ms=3,
                label=tr(lang, "baseline ranking, not recomputed", "按基线弹性排序、不重算"))
        ax.errorbar(fr, [b] + [r["random_mean"] for r in rows], yerr=[0] + [r["random_sd"] for r in rows], fmt="^:",
                    color=C["orange"], ms=3.5, lw=1.0, capsize=2, label=tr(lang, "random cells", "随机选格"))
        ax.errorbar(fr, [b] + [r["highest_q_mean"] for r in rows], yerr=[0] + [r["highest_q_sd"] for r in rows], fmt="s--",
                    color=C["green"], ms=3.5, lw=1.0, capsize=2, label=tr(lang, "highest $q$ first", "$q$ 最大者优先"))
        ax.plot(fr, [b] + [r["adaptive"] for r in rows], "o-", color=C["blue"], lw=1.6, ms=3,
                label=tr(lang, "elasticity-guided, recomputed", "按弹性分布、逐步重算"))
        ax.set_xlabel(tr(lang, r"area with $q\to q/2$ (%)", r"$q\to q/2$ 的面积占比（%）"))
        ax.set_ylabel(r"$R_0$")
        ax.set_ylim(1.45, 5.3)
        ax.set_xlim(-3, 103)
        if legend:
            ax.legend(loc="lower left", fontsize=7)
        d0 = key.split("=")[1]
        panel(ax, letter, tr(lang, f"targeted ventilation, $D_0={d0}$", f"定向通风，$D_0={d0}$"))

    targeting(axs[1, 0], "D0=1", "c", True)
    targeting(axs[1, 1], "D0=100", "d", False)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    print(bilingual(draw, "fig_s4_interventions"), flush=True)
