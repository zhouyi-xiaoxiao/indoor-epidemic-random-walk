#!/usr/bin/env python
"""All figures of the simulator part, each in English (_en) and Simplified
Chinese (_zh), as vector PDF and PNG.  Reads only files in data/.

    python experiments/20_figures.py            # all figures whose data exist
    python experiments/20_figures.py scenes finite_size
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import matplotlib                                         # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt                           # noqa: E402
import matplotlib.ticker                                  # noqa: E402
from matplotlib.colors import ListedColormap              # noqa: E402
from matplotlib.patches import Patch                      # noqa: E402

from bsc_sim.plotting import C, bilingual, tr             # noqa: E402
from bsc_sim.scenes import SCENE_NAMES, load_scene        # noqa: E402

DATA = os.path.join(HERE, "..", "data")
BETA, GAMMA = 0.5, 0.14


def jload(name):
    return json.load(open(os.path.join(DATA, name)))


def have(*names):
    return all(os.path.exists(os.path.join(DATA, n)) for n in names)


def panel_label(ax, s, x=-0.14, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontweight="bold", fontsize=10,
            va="bottom", ha="left")


# =========================================================================
# Scene maps
# =========================================================================
ZONE_COLOURS = {
    "office": {"W": C["sky"], "C": C["yellow"], "M": C["red"], "K": C["green"],
               "B": "#333333"},
    "supermarket": {"S": C["sky"], "A": C["yellow"], "F": C["green"],
                    "Q": C["red"], "E": C["purple"], "B": "#333333"},
    "classroom": {"d": C["sky"], "D": C["blue"], "P": C["orange"],
                  "I": C["yellow"], "O": "#333333"},
    "metro": {"T": C["blue"], "N": C["sky"], "G": C["yellow"], "R": C["red"]},
}


def fig_scenes(lang):
    fig = plt.figure(figsize=(7.4, 9.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.36], hspace=1.0,
                          wspace=0.22)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[2, :])]
    for ax, name, lab in zip(axes, SCENE_NAMES, "abcd"):
        sc = load_scene(name)
        zones = list(ZONE_COLOURS[name].keys())
        code = np.zeros((sc.ny, sc.nx), dtype=int)
        for k, z in enumerate(zones):
            code[sc.zone == z] = k
        cmap = ListedColormap([ZONE_COLOURS[name][z] for z in zones])
        ax.imshow(code, cmap=cmap, vmin=-0.5, vmax=len(zones) - 0.5,
                  origin="upper", interpolation="nearest",
                  extent=(0, sc.nx * sc.a_m, sc.ny * sc.a_m, 0))
        ax.set_xticks(np.arange(0, sc.nx + 1) * sc.a_m, minor=True)
        ax.set_yticks(np.arange(0, sc.ny + 1) * sc.a_m, minor=True)
        ax.grid(which="minor", color="white", lw=0.3)
        ax.tick_params(which="minor", length=0)
        ax.set_aspect("equal")
        ax.set_xlabel(tr(lang, "x (m)", "x（米）"), labelpad=1)
        ax.set_ylabel(tr(lang, "y (m)", "y（米）"), labelpad=1)
        for s in ("top", "right"):
            ax.spines[s].set_visible(True)
        title = sc.meta["title_en"] if lang == "en" else sc.meta["title_zh"]
        ax.set_title(
            f"({lab}) {title}: {sc.nx}×{sc.ny}, $a$ = {sc.a_m:g} m, $N$ = {sc.N}",
            fontsize=8.5)
        handles = []
        zp = sc.meta["zones"]
        rho = sc.meta["rho_per_m2"]
        for z in zones:
            p = zp[z]
            nm = p["name_en"] if lang == "en" else p["name_zh"]
            if p.get("barrier"):
                txt = nm
            else:
                d = p["D"] if "D" in p else p["D_coef"] / rho ** p["D_exp"]
                txt = f"{nm}: $D/D_0$ = {d:.2g}, $q$ = {p['q']:g}"
            handles.append(Patch(facecolor=ZONE_COLOURS[name][z],
                                 edgecolor="k", lw=0.3, label=txt))
        if name == "classroom":
            ax.legend(handles=handles, loc="center left",
                      bbox_to_anchor=(1.12, 0.5), fontsize=7.5)
        elif name == "metro":
            ax.legend(handles=handles, loc="upper center", ncol=2,
                      bbox_to_anchor=(0.5, -0.75), fontsize=7.5)
        else:
            ax.legend(handles=handles, loc="upper center", ncol=1,
                      bbox_to_anchor=(0.5, -0.24), fontsize=7.2,
                      handlelength=1.2, labelspacing=0.25)
    return fig


# =========================================================================
# Finite-size figure
# =========================================================================
def fig_finite_size(lang):
    d = jload("02_finite_size.json")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, D, lab in zip(axes, (1.0, 10.0), "ab"):
        rows = sorted([v for k, v in d.items()
                       if k.startswith(f"D{D}_L") and "Linf" not in k],
                      key=lambda r: r["L"])
        Ls = np.array([r["L"] for r in rows])
        Lf = np.geomspace(3.5, 60, 200)
        ax.axhline(BETA / GAMMA, color=C["black"], lw=1.2,
                   label=tr(lang, "mean field: $R_0=\\beta q/\\gamma$",
                            "平均场：$R_0=\\beta q/\\gamma$"))
        ax.plot(Ls, [r["pair_closed_form"] for r in rows], color=C["blue"],
                label=tr(lang, "pair level, closed form for $\\bar{R}_1$",
                         "配对层次：$\\bar{R}_1$ 的闭式"))
        ex = [(r["L"], r["pair_exact"]) for r in rows
              if r["pair_exact"] == r["pair_exact"]]
        ax.plot(*zip(*ex), ls="none", marker="x", ms=5, color=C["blue"],
                label=tr(lang, "pair level, exact $\\bar{R}_1$", "配对层次：精确的 $\\bar{R}_1$"))
        y = np.array([r["sim_offspring"] for r in rows])
        lo = np.array([r["sim_lo"] for r in rows])
        hi = np.array([r["sim_hi"] for r in rows])
        ax.errorbar(Ls, y, yerr=[y - lo, hi - y], ls="none", marker="o", ms=4,
                    mfc="white", color=C["orange"], capsize=2,
                    label=tr(lang, "simulation: counted $\\widehat{R}$",
                             "模拟：计数值 $\\widehat{R}$"))
        ax.axhline(d[f"D{D}_Linf"]["pair_closed_form"], color=C["blue"],
                   ls=":", lw=1.0,
                   label=tr(lang, "closed form, $L\\to\\infty$",
                            "闭式，$L\\to\\infty$"))
        ax.set_xscale("log")
        ax.set_xticks([4, 5, 7, 10, 14, 20, 28, 40])
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xlim(3.5, 46)
        ax.set_ylim(0, 4.0)
        ax.set_xlabel(tr(lang, "room side $L$ (cells)", "房间边长 $L$（格点）"))
        ax.set_title(tr(lang, f"({lab}) $D$ = {D:g} cell$^2$ day$^{{-1}}$",
                        f"({lab}) $D$ = {D:g} 格$^2$ 天$^{{-1}}$"))
    axes[0].set_ylabel(tr(lang, "reproduction number", "再生数"))
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.02), columnspacing=1.2)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    return fig


FIGS = {
    "scenes": (fig_scenes, "fig3_scene_maps", ()),
    "finite_size": (fig_finite_size, "fig5_1_finite_size",
                    ("02_finite_size.json",)),
}


def register(key, name, *needs):
    def deco(f):
        FIGS[key] = (f, name, needs)
        return f
    return deco


if __name__ == "__main__":
    # later figure definitions live in 21_figures_more.py to keep files short
    import importlib.util
    more = os.path.join(HERE, "21_figures_more.py")
    if os.path.exists(more):
        spec = importlib.util.spec_from_file_location("figs_more", more)
        mod = importlib.util.module_from_spec(spec)
        mod.register = register
        mod.jload, mod.have, mod.panel_label = jload, have, panel_label
        spec.loader.exec_module(mod)
    which = sys.argv[1:] or list(FIGS)
    for key in which:
        f, name, needs = FIGS[key]
        if not have(*needs):
            print(f"skip {key}: missing data {needs}")
            continue
        paths = bilingual(f, name)
        print("wrote", [os.path.basename(p) for p in paths])
