#!/usr/bin/env python
"""s6_scenes_fig_spatial.py -- Fig. s6_scenes:fig:spatial (file fig5_2_spatial_pattern_{en,zh}).

Where infections happen in the office scene at D0 = 1 cell^2/day: mean-field prediction, pair-level
prediction and exact simulation (40 000 outbreaks, uniformly placed index case, generations 0-2
infectious).  Regenerated from the data of the simulations; replaces the copy of
code/bsc_sim/figures/fig5_2_spatial_pattern_* (whose legend overlapped the points) and adds the
mean-field map and the breakdown by cell class: the hottest cells are workstation cells bordering a
corridor, not the corridor cells (notes/VERIFICATION.md of the repository, Section 4.2, further point 3).

Inputs  data/bsc_sim/04_spatial_office_D1.npz   counts and predictions per cell, generations 1-3
        data/bsc_sim/04_spatial.json            correlations quoted in the legend
        data/s6_scenes_numbers.json          cell classes (code/s6_scenes_tables.py)
Output  figures/fig5_2_spatial_pattern_en.pdf|png and ..._zh.pdf|png

Run:    python code/s6_scenes_tables.py      (first)
        python code/s6_scenes_fig_spatial.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
SIM = ROOT / "code" / "bsc_sim"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(SIM / "src"))
import matplotlib.pyplot as plt                          # noqa: E402
from matplotlib.colors import Normalize                  # noqa: E402
from matplotlib.lines import Line2D                      # noqa: E402

from figstyle import C, bilingual, tr                    # noqa: E402
from bsc_sim.scenes import load_scene                    # noqa: E402

NAME, D0, GEN = "office", 1, 2
Z = np.load(SIM / "data" / f"04_spatial_{NAME}_D{D0}.npz")
JS = json.load(open(SIM / "data" / "04_spatial.json"))[f"{NAME}|D0={D0}"]
NUM = json.load(open(HERE.parent / "data" / "s6_scenes_numbers.json"))["office_cell_classes"][f"D0={D0}"]
SC = load_scene(NAME, D0=float(D0))
INK, MUTED = "#222222", "#6b6b6b"


def zone_edges(ax):
    """Thin lines on the edges between cells of different zones (barriers included)."""
    zone = SC.zone
    ny, nx = zone.shape
    for y in range(ny):
        for x in range(nx):
            if x + 1 < nx and zone[y, x] != zone[y, x + 1]:
                ax.plot([x + 0.5, x + 0.5], [y - 0.5, y + 0.5], color=INK, lw=0.55, solid_capstyle="butt")
            if y + 1 < ny and zone[y, x] != zone[y + 1, x]:
                ax.plot([x - 0.5, x + 0.5], [y + 0.5, y + 0.5], color=INK, lw=0.55, solid_capstyle="butt")


def draw(lang):
    M = SC.M
    pair = Z[f"pred_pair_gen{GEN}"] * M
    mf = Z[f"pred_mf_gen{GEN}"] * M
    obs = Z[f"counts_gen{GEN}"] / Z[f"counts_gen{GEN}"].sum() * M
    J = JS[f"gen{GEN}"]

    fig = plt.figure(figsize=(7.4, 5.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.78, 1.25], hspace=0.28)
    top = gs[0].subgridspec(1, 4, width_ratios=[1, 1, 1, 0.045], wspace=0.08)
    bot = gs[1].subgridspec(1, 2, width_ratios=[0.60, 1.0], wspace=0.30)

    # ---- maps: diverging scale centred on the uniform value 1 (warm = above, cool = below)
    cmap = plt.get_cmap("BrBG_r").copy()
    cmap.set_bad("#8c8c8c")
    norm = Normalize(vmin=0.2, vmax=1.8)
    panels = ((mf, tr(lang, "(a) mean field", "(a) 平均场")),
              (pair, tr(lang, "(b) pair level", "(b) 配对层次")),
              (obs, tr(lang, "(c) simulation", "(c) 模拟")))
    for k, (v, ttl) in enumerate(panels):
        ax = fig.add_subplot(top[0, k])
        im = ax.imshow(np.ma.masked_invalid(SC.to_grid(v)), cmap=cmap, norm=norm, origin="upper",
                       interpolation="nearest")
        zone_edges(ax)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_linewidth(0.6)
        ax.set_title(ttl, fontsize=9, loc="left")
    cax = fig.add_subplot(top[0, 3])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(tr(lang, "relative infection\nfrequency (1 = uniform)", "相对感染频率\n（1 = 均匀）"), fontsize=8)
    cb.outline.set_linewidth(0.5)

    # ---- (d) cell by cell
    ax = fig.add_subplot(bot[0, 0])
    lim = 1.8
    ax.plot([0, lim], [0, lim], color=MUTED, lw=0.8)
    ax.plot(mf, obs, ls="none", marker="s", ms=3.0, mfc=C["red"], mec="white", mew=0.3, alpha=0.8)
    ax.plot(pair, obs, ls="none", marker="o", ms=3.2, mfc=C["blue"], mec="white", mew=0.3, alpha=0.85)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1.0, 1.5])
    ax.set_yticks([0, 0.5, 1.0, 1.5])
    ax.set_xlabel(tr(lang, "predicted relative frequency", "预测的相对频率"))
    ax.set_ylabel(tr(lang, "simulated relative frequency", "模拟的相对频率"))
    ax.set_title(tr(lang, "(d) cell by cell, generation 2", "(d) 逐格比较（第 2 代）"), fontsize=9, loc="left")
    h = [Line2D([], [], ls="none", marker="s", ms=4.5, mfc=C["red"], mec="white", mew=0.3,
                label=tr(lang, f"mean field, $r$ = {J['r_meanfield']:.2f}",
                         f"平均场，$r$ = {J['r_meanfield']:.2f}").replace("-", "−")),
         Line2D([], [], ls="none", marker="o", ms=4.5, mfc=C["blue"], mec="white", mew=0.3,
                label=tr(lang, f"pair level, $r$ = {J['r_pair']:.2f}", f"配对层次，$r$ = {J['r_pair']:.2f}"))]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1, fontsize=7.5,
              handletextpad=0.3, borderaxespad=0.0)

    # ---- (e) by cell class, generations 1-3
    ax = fig.add_subplot(bot[0, 1])
    labels = dict(
        W_border=tr(lang, "workstations\nby a corridor", "工位\n（邻走廊）"),
        W_interior=tr(lang, "other work-\nstations", "工位\n（其余）"),
        C=tr(lang, "corridor", "走廊"), M=tr(lang, "meeting\nroom", "会议室"), K=tr(lang, "kitchen", "茶水间"))
    order = ["W_border", "W_interior", "C", "M", "K"]
    dx = 0.24
    ax.axhline(1.0, color=MUTED, lw=0.8)
    for i, key in enumerate(order):
        e = NUM["classes"][key]
        xs = i + dx * (np.arange(3) - 1)
        ax.plot(xs, e["meanfield"], color=C["red"], lw=1.0, ls=":", marker="s", ms=4, mec="white", mew=0.4)
        ax.plot(xs, e["pair"], color=C["blue"], lw=1.0, ls="--", marker="o", ms=4, mec="white", mew=0.4)
        ax.plot(xs, e["sim"], color=INK, lw=1.4, marker="D", ms=4.2, mec="white", mew=0.4)
        for x, g in zip(xs, (1, 2, 3)):
            ax.text(x, 0.335, str(g), ha="center", va="bottom", fontsize=6.5, color=MUTED)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([labels[k] for k in order], fontsize=6.8)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_xlim(-0.55, len(order) - 0.45)
    ax.set_ylim(0.32, 1.52)
    ax.set_ylabel(tr(lang, "share of infections / share of floor", "感染占比 / 面积占比"))
    ax.set_title(tr(lang, "(e) by cell class, generations 1, 2, 3", "(e) 按格点类别（第 1、2、3 代）"),
                 fontsize=9, loc="left")
    h = [Line2D([], [], color=INK, lw=1.4, marker="D", ms=4.2, mec="white", mew=0.4,
                label=tr(lang, "simulation", "模拟")),
         Line2D([], [], color=C["blue"], lw=1.0, ls="--", marker="o", ms=4, mec="white", mew=0.4,
                label=tr(lang, "pair level", "配对层次")),
         Line2D([], [], color=C["red"], lw=1.0, ls=":", marker="s", ms=4, mec="white", mew=0.4,
                label=tr(lang, "mean field", "平均场"))]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=3, fontsize=7.5,
              handletextpad=0.4, columnspacing=1.4, borderaxespad=0.0)
    return fig


if __name__ == "__main__":
    print(bilingual(draw, "fig5_2_spatial_pattern"))
