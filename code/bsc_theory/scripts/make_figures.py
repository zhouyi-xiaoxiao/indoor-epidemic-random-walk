"""make_figures.py -- all figures of the theory part, each in English (_en) and Simplified Chinese (_zh),
as vector PDF and PNG, written to ../figures/.   Data come from ../results/*.json (produced by the other scripts);
only the room-size curves of Fig. 2 are computed here (sparse eigen-solves, a few seconds).

    python make_figures.py            # all figures
    python make_figures.py 3 5        # only figures 3 and 5
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.colors import ListedColormap, LogNorm

sys.path.insert(0, str(Path(__file__).parent))
import ngm_lattice as ng

ROOT = Path(__file__).resolve().parent.parent
RESD, FIGD = ROOT / "results", ROOT / "figures"
FIGD.mkdir(exist_ok=True)
BETA, GAMMA = 0.5, 0.14
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00", purple="#CC79A7", sky="#56B4E9", grey="#666666", black="#000000")


def style(lang):
    plt.rcParams.update({
        "font.family": ["Songti SC", "STIXGeneral"] if lang == "zh" else ["STIXGeneral"],
        "mathtext.fontset": "stix", "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42,
        "axes.unicode_minus": False, "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 100, "savefig.dpi": 300, "legend.frameon": False,
    })


def save(fig, name, lang):
    for ext in ("pdf", "png"):
        fig.savefig(FIGD / f"{name}_{lang}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote", name, lang, flush=True)


def T(lang, en, zh):
    return zh if lang == "zh" else en


def panel(ax, letter):
    t = ax.get_title()
    ax.set_title("")
    ax.set_title(f"({letter}) {t}", loc="left")


WE = json.loads((RESD / "worked_examples.json").read_text())


# ======================================================================================= Fig 1
def fig1(lang):
    style(lang)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.25), sharey=False)
    E1, E2 = WE["E1_two_zone"], WE["E2_office"]
    D = np.array(E1["Dsweep"])
    ax = axs[0]
    names = {"side": T(lang, "side corridor", "单侧走廊"), "aisles": T(lang, "parallel aisles", "平行过道"),
             "grid": T(lang, "aisle grid", "网格过道"), "scatter": T(lang, "scattered cells", "随机散布")}
    cols = {"side": C["blue"], "aisles": C["green"], "grid": C["orange"], "scatter": C["purple"]}
    for k in ("side", "aisles", "grid", "scatter"):
        ax.semilogx(D, E1[k]["symmetric"]["sweep"], color=cols[k], lw=1.6, label=names[k])
    lo, hi = E1["side"]["symmetric"]["lower"], E1["side"]["symmetric"]["upper"]
    ax.axhline(hi, color=C["black"], ls="--", lw=0.9)
    ax.axhline(lo, color=C["black"], ls="--", lw=0.9)
    ax.text(D[-1], hi + 0.05, r"$\beta q_{\max}/\gamma$", ha="right", va="bottom", fontsize=8)
    ax.text(D[-1], lo - 0.07, r"$\beta\langle q\rangle_\pi/\gamma$", ha="right", va="top", fontsize=8)
    ax.axvline(1.0, color=C["grey"], lw=0.8, ls="-.")
    ax.text(1.15, 3.1, T(lang, "default\n$D_0=1$", "默认取值\n$D_0=1$"), color=C["grey"], fontsize=7.5, va="center")
    ax.axvspan(2.4e3, D[-1], color="0.9", zorder=0)
    ax.text(1.4e4, 3.1, T(lang, "walking\n$\\geq$5% of\nthe time", "行走时间\n占比$\\geq$5%"), fontsize=7.5, ha="center", va="center", color="0.3")
    ax.set_ylim(1.7, 5.75)
    ax.set_xlim(D[0], D[-1])
    ax.set_xlabel(T(lang, r"mobility scale $D_0$ (lattice$^2$/day)", r"迁移率尺度 $D_0$（格点$^2$/天）"))
    ax.set_ylabel(T(lang, r"basic reproduction number $R_0$", r"基本再生数 $R_0$"))
    ax.set_title(T(lang, "Two-zone office (workspace 69%, corridor 31%)", "两区域办公室（工位区 69%，走廊 31%）"))
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.57))
    panel(ax, "a")

    ax = axs[1]
    D2 = np.array(E2["Dsweep"])
    ax.semilogx(D2, E2["O1"]["symmetric"]["sweep"], color=C["blue"], lw=1.6, label=T(lang, "layout O1, symmetric rates (3-26)", "布局 O1，对称跳跃率 (3-26)"))
    ax.semilogx(D2, E2["O2"]["symmetric"]["sweep"], color=C["green"], lw=1.6, label=T(lang, "layout O2, symmetric rates", "布局 O2，对称跳跃率"))
    ax.semilogx(D2, E2["O1"]["departure"]["sweep"], color=C["orange"], lw=1.6, label=T(lang, "layout O1, departure-site rates (3-25)", "布局 O1，出发点跳跃率 (3-25)"))
    o1 = dict(E2["O1"]["symmetric"])
    o1["largeD_C"], o1["smallD_muZ"] = E2["O1"]["largeD_C"], E2["O1"]["smallD_muZ"]
    lo, hi = o1["lower"], o1["upper"]
    # asymptotes for O1
    Dl = D2[D2 >= 0.6]
    ax.semilogx(Dl, lo + BETA * o1["largeD_C"] / (o1["mean_q_pi"] * Dl), color=C["blue"], lw=0.9, ls=":")
    Ds = D2[D2 <= 2.2]
    ax.semilogx(Ds, hi * (1 - Ds * o1["smallD_muZ"] / GAMMA), color=C["blue"], lw=0.9, ls=":",
                label=T(lang, "asymptotes (Thm 3) for O1", "O1 的渐近式（定理 3）"))
    for yv, lab in ((hi, r"$\beta q_{\max}/\gamma=5.36$"), (lo, r"$\beta\langle q\rangle/\gamma=4.64$")):
        ax.axhline(yv, color=C["black"], ls="--", lw=0.9)
    ax.text(D2[-1], hi + 0.03, r"$\beta q_{\max}/\gamma=5.36$", ha="right", va="bottom", fontsize=8)
    ax.text(D2[-1], lo - 0.03, r"$\beta\langle q\rangle/\gamma=4.64$", ha="right", va="top", fontsize=8)
    ax.axhline(E2["O1"]["departure"]["lower"], color=C["orange"], ls="--", lw=0.8)
    ax.axvline(1.0, color=C["grey"], lw=0.8, ls="-.")
    ax.axvspan(2.4e3, D2[-1], color="0.9", zorder=0)
    ax.set_ylim(3.3, 5.6)
    ax.set_xlim(D2[0], D2[-1])
    ax.set_xlabel(T(lang, r"mobility scale $D_0$ (lattice$^2$/day)", r"迁移率尺度 $D_0$（格点$^2$/天）"))
    ax.set_title(T(lang, "Full office scene", "完整办公室场景"))
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.115), fontsize=7)
    panel(ax, "b")
    fig.tight_layout()
    save(fig, "fig1_R0_vs_mobility", lang)


# ======================================================================================= Fig 2
_F2 = {}


def _fig2_data():
    if _F2:
        return _F2
    Ls = np.array([3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 30, 40])
    door, wall, refl = [], [], []
    for L in Ls:
        lat = ng.Lattice(np.ones((L, L), bool))
        M = ng.generator_symmetric(lat, np.full(lat.n, 0.5))
        refl.append(ng.R0_dense(M, np.ones(lat.n), BETA, GAMMA) * GAMMA / BETA if L <= 20 else ng.R0_sparse(M, np.ones(lat.n), BETA, GAMMA) * GAMMA / BETA)
        for store, cells in ((door, [lat.index[0, L // 2]]), (wall, list(lat.index[0, :]))):
            kap = np.zeros(lat.n)
            kap[cells] = 0.5
            mu = ng.mu_leak(M, kap) if lat.n <= 500 else -ng.growth_rate_sparse(M, np.zeros(lat.n), 0.0, 0.0, kappa=kap)
            store.append(GAMMA / (GAMMA + mu))
    L1 = np.array([3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 30, 40, 60, 100])
    r1 = []
    for L in L1:
        lat = ng.Lattice(np.ones((1, L), bool))
        M = ng.generator_symmetric(lat, np.ones(L))
        kap = np.zeros(L)
        kap[0] = kap[-1] = 1.0
        r1.append(ng.R0_dense(M, np.ones(L), BETA, GAMMA, kappa=kap))
    _F2.update(Ls=Ls, door=np.array(door), wall=np.array(wall), refl=np.array(refl), L1=L1, r1=np.array(r1))
    (RESD / "fig2_room_size_data.json").write_text(json.dumps({k: v.tolist() for k, v in _F2.items()}, indent=1))
    return _F2


def fig2(lang):
    style(lang)
    d = _fig2_data()
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.1))
    ax = axs[0]
    Lc = np.linspace(2.5, 42, 300)
    ax.plot(d["Ls"], d["refl"], "o", color=C["blue"], ms=4.5, label=T(lang, "reflecting walls (exact, all $D$)", "反射壁（精确，任意 $D$）"))
    ax.plot(d["Ls"], d["door"], "s-", color=C["green"], ms=3.5, lw=1.2, label=T(lang, "one door cell (leak)", "一个门格点（泄漏）"))
    ax.plot(d["Ls"], d["wall"], "^-", color=C["orange"], ms=3.5, lw=1.2, label=T(lang, "one fully open side", "一侧完全敞开"))
    ax.set_xscale("log")
    ax.set_xticks([3, 5, 10, 20, 40])
    ax.set_xticklabels(["3", "5", "10", "20", "40"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0.15, 1.04)
    ax.set_xlabel(T(lang, r"room side $L$ (cells)", r"房间边长 $L$（格点）"))
    ax.set_ylabel(T(lang, r"$R_0\,/\,(\beta q/\gamma)$", r"$R_0\,/\,(\beta q/\gamma)$"))
    ax.set_title(T(lang, r"Square room, uniform $q$, $D=0.5$, $\gamma=0.14$", r"方形房间，均匀 $q$，$D=0.5$，$\gamma=0.14$"))
    ax.legend(loc="lower right")
    panel(ax, "a")

    ax = axs[1]
    Lc = np.linspace(2.5, 110, 400)
    ax.plot(d["L1"], d["r1"], "o", color=C["blue"], ms=4.5, label=T(lang, "NGM, corridor with open ends", "NGM：两端开放的走廊"))
    ax.plot(Lc, BETA / (GAMMA + 2 * (1 - np.cos(np.pi / (Lc + 1)))), color=C["blue"], lw=1.1,
            label=r"$\beta q/[\gamma+2w(1-\cos\frac{\pi}{L+1})]$")
    ax.plot(Lc, BETA / (GAMMA + np.pi**2 / Lc**2), color=C["red"], ls=":", lw=1.4, label=r"$\beta q/(\gamma+\pi^2D/L^2)$")
    ax.axhline(BETA / GAMMA, color=C["black"], ls="--", lw=0.9)
    ax.text(110, BETA / GAMMA + 0.04, T(lang, r"closed room: $\beta q/\gamma$", r"封闭房间：$\beta q/\gamma$"), ha="right", va="bottom", fontsize=8)
    ax.set_xscale("log")
    ax.set_xticks([3, 10, 30, 100])
    ax.set_xticklabels(["3", "10", "30", "100"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0, 3.95)
    ax.set_xlabel(T(lang, r"corridor length $L$ (cells)", r"走廊长度 $L$（格点）"))
    ax.set_ylabel(T(lang, r"room reproduction number $R_0^{\rm room}$", r"房间再生数 $R_0^{\rm room}$"))
    ax.set_title(T(lang, r"1-D corridor, absorbing ends, $D=w=1$", r"一维走廊，两端吸收，$D=w=1$"))
    ax.legend(loc="lower right")
    panel(ax, "b")
    fig.tight_layout()
    save(fig, "fig2_room_size", lang)


# ======================================================================================= Fig 3
def fig3(lang):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    style(lang)
    plt.rcParams.update({"axes.spines.left": False, "axes.spines.bottom": False})
    m = WE["E2_office"]["O1"]["maps"]
    zones = np.array(m["zones"])
    n = WE["E2_office"]["O1"]["symmetric"]["n"]
    fig, axs = plt.subplots(2, 3, figsize=(7.6, 4.2))
    code = {"W": 0, "C": 1, "M": 2, "K": 3, "B": 4}
    zi = np.vectorize(code.get)(zones)
    cmapz = ListedColormap(["#c6dbef", "#f7f7f7", "#fdae6b", "#a1d99b", "#525252"])
    ax = axs[0, 0]
    ax.imshow(zi, origin="lower", cmap=cmapz, vmin=-0.5, vmax=4.5, interpolation="nearest")
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.06)
    cax.axis("off")
    labs = [T(lang, "desks W", "工位 W"), T(lang, "corridor C", "走廊 C"), T(lang, "meeting room M", "会议室 M"),
            T(lang, "kitchen K", "茶水间 K"), T(lang, "obstacle B", "障碍物 B")]
    handles = [plt.Rectangle((0, 0), 1, 1, fc=cmapz(i), ec="0.4", lw=0.4) for i in range(5)]
    ax.legend(handles, labs, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=3, fontsize=6.5, handlelength=1, columnspacing=0.8)
    ax.set_title(T(lang, "(a) layout O1 (20$\\times$13 cells)", "(a) 布局 O1（20$\\times$13 格点）"), loc="left")

    def show(ax, arr, title, norm=None, cmap="viridis", cbl=None, vmin=None, vmax=None):
        a = np.array(arr, float)
        cm = plt.get_cmap(cmap).copy()
        cm.set_bad("#525252")
        im = ax.imshow(np.ma.masked_invalid(a), origin="lower", cmap=cm, norm=norm, interpolation="nearest",
                       vmin=None if norm else vmin, vmax=None if norm else vmax)
        cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.06)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=7)
        if cbl:
            cb.set_label(cbl, fontsize=7.5)
        ax.set_title(title, loc="left")

    show(axs[0, 1], m["q"], T(lang, r"(b) infection efficiency $q(\mathbf{r})$", r"(b) 传染效率 $q(\mathbf{r})$"), cmap="YlOrRd", vmin=0.7, vmax=1.6)
    lab = T(lang, r"relative density $n\,u_{\mathbf{r}}$", r"相对密度 $n\,u_{\mathbf{r}}$")
    phi0 = np.array(m["phi0"], float) ** 2 * n
    show(axs[0, 2], phi0, T(lang, r"(c) lowest mode $n|\phi_0|^2$ (constant)", r"(c) 最低模 $n|\phi_0|^2$（常数）"), norm=LogNorm(0.03, 30))
    for ax, key, d0, letter in ((axs[1, 0], "u", 1, "d"), (axs[1, 1], "u_D10", 10, "e"), (axs[1, 2], "u_D100", 100, "f")):
        u = np.array(m[key], float) * n
        show(ax, u, T(lang, rf"({letter}) Perron vector $u$, $D_0={d0}$", rf"({letter}) Perron 向量 $u$，$D_0={d0}$"), norm=LogNorm(0.03, 30), cbl=lab)
    for ax in axs.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout(h_pad=1.6)
    save(fig, "fig3_hotspot_maps", lang)


# ======================================================================================= Fig 4
def fig4(lang):
    style(lang)
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
    E4, E3, E5 = WE["E4_corridor_D"], WE["E3_partitions"], WE["E5_targeting"]
    ax = axs[0]
    c = np.array(E4["c"])
    ax.semilogx(c, E4["O1"]["symmetric"]["R0"], color=C["blue"], lw=1.6, label=T(lang, "O1, symmetric (3-26)", "O1，对称 (3-26)"))
    ax.semilogx(c, E4["O1"]["departure"]["R0"], color=C["orange"], lw=1.6, label=T(lang, "O1, departure-site (3-25)", "O1，出发点 (3-25)"))
    ax.semilogx(c, E4["two_zone_aisles"]["symmetric"]["R0"], color=C["blue"], lw=1.2, ls="--", label=T(lang, "two-zone aisles, symmetric", "两区过道，对称"))
    ax.semilogx(c, E4["two_zone_aisles"]["departure"]["R0"], color=C["orange"], lw=1.2, ls="--", label=T(lang, "two-zone aisles, departure", "两区过道，出发点"))
    ax.axvline(1, color=C["grey"], lw=0.8, ls="-.")
    ax.set_xlabel(T(lang, r"corridor mobility multiplier", r"走廊迁移率倍数"))
    ax.set_ylabel(r"$R_0$")
    ax.set_title(T(lang, "'Strategy 2': slower corridors", "“策略 2”：降低走廊迁移率"))
    ax.legend(loc="lower right", fontsize=6.3)
    ax.annotate("", xy=(0.02, 5.33), xytext=(0.6, 5.33), arrowprops=dict(arrowstyle="->", color=C["grey"], lw=0.8))
    ax.text(0.11, 5.38, T(lang, "restriction", "限制流动"), fontsize=7, color=C["grey"], ha="center")
    ax.set_ylim(1.9, 5.5)
    panel(ax, "a")

    ax = axs[1]
    base = E3["R0_base"]
    ps = sorted(float(k) for k in E3["edge_partitions"])
    mean = np.array([E3["edge_partitions"][str(p) if str(p) in E3["edge_partitions"] else repr(p)]["mean"] for p in ps])
    sd = np.array([E3["edge_partitions"][str(p) if str(p) in E3["edge_partitions"] else repr(p)]["sd"] for p in ps])
    ncl = np.array([E3["edge_partitions"][str(p) if str(p) in E3["edge_partitions"] else repr(p)]["mean_closed_edges"] for p in ps])
    frac = ncl / E3["n_WW_edges"]
    ax.errorbar(np.concatenate([[0], frac]) * 100, np.concatenate([[0], 100 * (mean / base - 1)]),
                yerr=np.concatenate([[0], 100 * sd / base]), fmt="o-", color=C["blue"], ms=3.5, lw=1.2, capsize=2,
                label=T(lang, "walls on desk–desk edges", "工位间隔板（封闭边）"))
    ob = E3["obstacle_cells_plus15pct"]
    ax.errorbar([25], [100 * (ob["mean"] / base - 1)], yerr=[100 * ob["sd"] / base], fmt="s", color=C["green"], ms=4.5, capsize=2,
                label=T(lang, "39 desk cells (25%) $\\to$ obstacles", "39 个工位格点（25%）改为障碍物"))
    ax.axhline(0, color=C["black"], lw=0.7)
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel(T(lang, "desk–desk edges closed (%)", "被封闭的工位间连接（%）"))
    ax.set_ylabel(T(lang, r"change of $R_0$ (%)", r"$R_0$ 的变化（%）"))
    ax.set_title(T(lang, "Partitions (layout O1, $D_0=1$)", "隔板（布局 O1，$D_0=1$）"))
    ax.legend(loc="upper left", fontsize=6.5)
    panel(ax, "b")

    ax = axs[2]
    e5 = E5["D0=1"]
    f = np.array(e5["fracs"]) * 100
    g = np.array(e5["greedy_curve"])
    ax.plot(np.concatenate([[0], g[:, 0] * 100]), np.concatenate([[e5["R0_base"]], g[:, 1]]), "-", color=C["blue"], lw=1.6,
            label=T(lang, "hotspot-guided (adaptive elasticity)", "热点引导（自适应弹性）"))
    ax.plot(np.concatenate([[0], f]), np.concatenate([[e5["R0_base"]], e5["by_q"]]), "s--", color=C["green"], ms=3, lw=1.1,
            label=T(lang, "highest $q$ first", "按 $q$ 从高到低"))
    rm, rs = np.array(e5["random_mean"]), np.array(e5["random_sd"])
    ax.errorbar(np.concatenate([[0], f]), np.concatenate([[e5["R0_base"]], rm]), yerr=np.concatenate([[0], rs]), fmt="^:", color=C["orange"],
                ms=3, lw=1.1, capsize=2, label=T(lang, "random cells", "随机选取格点"))
    ax.set_xlabel(T(lang, "area with improved ventilation (%)", "加强通风的面积占比（%）"))
    ax.set_ylabel(r"$R_0$")
    ax.set_title(T(lang, r"Targeted ventilation ($q\to q/2$)", r"定点通风（$q\to q/2$）"))
    ax.legend(loc="upper right", fontsize=6.5)
    ax.set_ylim(2.3, 5.9)
    panel(ax, "c")
    fig.tight_layout()
    save(fig, "fig4_interventions", lang)


# ======================================================================================= Fig 5
def fig5(lang):
    style(lang)
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
    fam = [T(lang, "symmetric", "对称"), T(lang, "departure-site", "出发点"), T(lang, "non-reversible", "不可逆"), T(lang, "dense random", "稠密随机")]
    cols = [C["blue"], C["orange"], C["green"], C["purple"]]
    pts = np.load(RESD / "T1_points.npy")
    ax = axs[0]
    ax.axhspan(0, 1e3, xmin=0.5, color="0.93", zorder=0)
    ax.axhspan(-1e3, 0, xmax=0.5, color="0.93", zorder=0)
    for k in range(4):
        s = pts[:, 2] == k
        ax.plot(np.log10(pts[s, 0]), np.arcsinh(pts[s, 1] * 5) / np.log(10), ".", color=cols[k], ms=3.5, label=fam[k])
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-1.3, 1.3)
    ax.set_xlabel(r"$\log_{10}R_0$")
    ax.set_ylabel(T(lang, r"$\lambda_1/\gamma$ (asinh scale)", r"$\lambda_1/\gamma$（asinh 标度）"))
    yt = np.array([-1, -0.3, 0, 0.3, 1, 3, 10, 30])
    ax.set_yticks(np.arcsinh(yt * 5) / np.log(10))
    ax.set_yticklabels([f"{v:g}" for v in yt])
    ax.set_ylim(np.arcsinh(-1.1 * 5) / np.log(10), np.arcsinh(60 * 5) / np.log(10))
    ax.set_title(T(lang, "Thm 1: sign$(R_0-1)$=sign$(\\lambda_1)$", "定理 1：sign$(R_0-1)$=sign$(\\lambda_1)$"))
    ax.legend(loc="upper left", fontsize=6.5, handletextpad=0.1)
    panel(ax, "a")

    ax = axs[1]
    pos = np.load(RESD / "T2_positions.npy")
    bins = np.linspace(0, 1, 26)
    bottom = np.zeros(len(bins) - 1)
    for k in range(4):
        h, _ = np.histogram(pos[pos[:, 0] == k, 1], bins)
        ax.bar(bins[:-1], h, width=np.diff(bins), align="edge", bottom=bottom, color=cols[k], label=fam[k], lw=0)
        bottom += h
    ax.set_xlim(-0.05, 1.05)
    ax.axvline(0, color="k", ls="--", lw=0.9)
    ax.axvline(1, color="k", ls="--", lw=0.9)
    ax.set_xlabel(r"$(R_0-\beta\langle q\rangle_\pi/\gamma)\,/\,(\beta q_{\max}/\gamma-\beta\langle q\rangle_\pi/\gamma)$", fontsize=7.5)
    ax.set_ylabel(T(lang, "number of random instances", "随机算例数"))
    ax.set_title(T(lang, "Thm 2: 1996 instances, all inside", "定理 2：1996 个算例均在界内"))
    panel(ax, "b")

    ax = axs[2]
    for name, col in (("O1", C["blue"]), ("O2", C["green"])):
        r = dict(WE["E2_office"][name])
        r["R0"] = r["symmetric"]["R0"]
        N = sorted(int(k) for k in r["ritz"])
        v = [r["ritz"][str(k)] for k in N]
        ax.semilogx(N, v, "o-", color=col, ms=3.5, lw=1.2, label=T(lang, f"Galerkin, layout {name}", f"Galerkin，布局 {name}"))
        ax.axhline(r["R0"], color=col, ls="--", lw=0.8)
        ax.plot([1.0], [r["diagonal_modal_max"]], marker="*", color=col, ms=9, mec="k", mew=0.4, ls="none")
    ax.plot([], [], marker="*", color="0.6", ms=9, mec="k", mew=0.4, ls="none", label=T(lang, "largest diagonal term (lower bound)", "最大对角项（下界）"))
    ax.plot([], [], color="0.4", ls="--", lw=0.8, label=T(lang, r"exact $R_0=\rho(K)$", r"精确值 $R_0=\rho(K)$"))
    ax.set_ylim(3.4, 5.25)
    ax.set_xlabel(T(lang, "number of diffusion modes $N$", "扩散本征模数 $N$"))
    ax.set_ylabel(r"$R_0^{(N)}$")
    ax.set_title(T(lang, "Mode truncation: from below", "模式截断：自下收敛"))
    ax.legend(loc="center right", fontsize=6.5, bbox_to_anchor=(1.0, 0.42))
    panel(ax, "c")
    fig.tight_layout()
    save(fig, "fig5_theorem_checks", lang)


# ======================================================================================= Fig 6
def fig6(lang):
    style(lang)
    NL = json.loads((RESD / "nonlinear_check.json").read_text())
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
    ax = axs[0]
    cols = [C["blue"], C["sky"], C["orange"], C["red"], C["purple"]]
    for (k, r), col in zip(sorted(NL["A_threshold"].items(), key=lambda kv: float(kv[0])), cols):
        if "ts" not in r:
            continue
        ax.semilogy(r["ts"], r["prevalence"], color=col, lw=1.4, label=rf"$R_0={float(k):.2f}$")
    ax.set_xlim(0, 1500)
    ax.set_ylim(1e-9, 1)
    ax.set_xlabel(T(lang, "time (days)", "时间（天）"))
    ax.set_ylabel(T(lang, "prevalence $I/N$", "感染者比例 $I/N$"))
    ax.set_title(T(lang, r"Threshold at $R_0=1$ ($\beta$ rescaled)", r"阈值 $R_0=1$（缩放 $\beta$）"))
    ax.legend(loc="upper right", fontsize=6.5)
    panel(ax, "a")

    ax = axs[1]
    B = NL["B_heterogeneity"]
    if "ts" in B:
        ax.semilogy(B["ts"], B["prevalence_heterogeneous"], color=C["red"], lw=1.6,
                    label=T(lang, rf"heterogeneous $q$: $R_0={B['R0']:.3f}$", rf"异质 $q$：$R_0={B['R0']:.3f}$"))
        ax.semilogy(B["ts"], B["prevalence_homogenised"], color=C["blue"], lw=1.6,
                    label=T(lang, r"same room, $q\equiv\langle q\rangle$: $R_0=0.900$", r"同一房间，$q\equiv\langle q\rangle$：$R_0=0.900$"))
    ax.set_ylim(1e-9, 1)
    ax.set_xlim(0, 3000)
    ax.set_xlabel(T(lang, "time (days)", "时间（天）"))
    ax.set_ylabel(T(lang, "prevalence $I/N$", "感染者比例 $I/N$"))
    ax.set_title(T(lang, r"Heterogeneity raises $R_0$", r"空间异质性抬高 $R_0$"))
    ax.text(0.97, 0.58, T(lang, f"final attack fraction\n{100 * B['final_attack_heterogeneous']:.1f}% vs {100 * B['final_attack_homogenised']:.3f}%",
                          f"最终感染比例\n{100 * B['final_attack_heterogeneous']:.1f}% 对 {100 * B['final_attack_homogenised']:.3f}%"),
            transform=ax.transAxes, ha="right", fontsize=7.5)
    ax.legend(loc="upper right", fontsize=6.5)
    panel(ax, "b")

    ax = axs[2]
    cc = {"1.0": C["red"], "10.0": C["orange"], "100.0": C["blue"]}
    for d0, r in NL["C_hotspot"].items():
        for seed, ls in (("far desk", "-"), ("meeting room", "--")):
            if seed not in r:
                continue
            s = r[seed]
            ts, tv, prev = np.array(s["ts"]), np.array(s["tv"]), np.array(s["prevalence"])
            ok = prev < 1e-2
            ax.semilogy(ts[ok], np.maximum(tv[ok], 1e-6), color=cc[d0], ls=ls, lw=1.4)
    for d0, col in cc.items():
        ax.plot([], [], color=col, lw=1.4, label=rf"$D_0={float(d0):g}$")
    ax.plot([], [], color="0.3", ls="-", lw=1.2, label=T(lang, "seed: far desk", "初始病例：远端工位"))
    ax.plot([], [], color="0.3", ls="--", lw=1.2, label=T(lang, "seed: meeting room", "初始病例：会议室"))
    ax.set_ylim(1e-5, 2)
    ax.set_xlabel(T(lang, "time (days)", "时间（天）"))
    ax.set_ylabel(T(lang, r"TV distance of $I(t)/|I(t)|$ to $u$", r"$I(t)/|I(t)|$ 与 $u$ 的全变差距离"))
    ax.set_title(T(lang, "Approach to hotspot map $u$", "向热点图 $u$ 的收敛"))
    ax.set_xlim(-2, 118)
    ax.set_xticks([0, 20, 40, 60])
    ax.legend(loc="center right", fontsize=6.3, bbox_to_anchor=(1.02, 0.60))
    panel(ax, "c")
    fig.tight_layout()
    save(fig, "fig6_nonlinear_checks", lang)


# ======================================================================================= Fig 7
def fig7(lang):
    style(lang)
    IB = json.loads((RESD / "ibm_check.json").read_text())
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.9))
    for ax, rule, title in ((axs[0], "local", T(lang, r"per-cell rule: rate $\beta q\,n_I/n_{\rm tot}$ (local count)", r"逐格规则：感染率 $\beta q\,n_I/n_{\rm tot}$（局部人数）")),
                            (axs[1], "mean", T(lang, r"rate $\beta q\,n_I/(N\pi_{\mathbf{r}})$ (mean occupancy)", r"感染率 $\beta q\,n_I/(N\pi_{\mathbf{r}})$（平均占据数）"))):
        rows = [r for r in IB.values() if r["rule"] == rule]
        labels, means, ses, preds = [], [], [], []
        for r in rows:
            labels.append(f"$N$={r['N']}\n$D_0$={r['D0']:g}")
            means.append(r["mean"])
            ses.append(r["se"])
            preds.append(r["annealed_prediction"])
        x = np.arange(len(rows))
        ax.errorbar(x, means, yerr=2 * np.array(ses), fmt="o", color=C["blue"], ms=4, capsize=2.5,
                    label=T(lang, "simulated offspring number ($\\pm2$ s.e.)", "模拟的一代子代数（$\\pm2$ 标准误）"))
        ax.plot(x, preds, "_", color=C["red"], ms=14, mew=1.8, label=T(lang, "annealed NGM prediction", "退火近似 NGM 预测"))
        ax.axhline(rows[0]["mean_field"], color="k", ls="--", lw=0.9)
        ax.text(-0.42, rows[0]["mean_field"] + 0.07, r"$\beta q/\gamma=5$", ha="left", va="bottom", fontsize=8)
        ax.axhline(1.0, color=C["grey"], ls=":", lw=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylim(0, 5.9)
        ax.set_xlim(-0.5, len(rows) - 0.5)
        ax.set_title(title)
        ax.set_ylabel(T(lang, "secondary cases of one index case", "单个指示病例的二代病例数"))
        if rule == "local":
            ax.legend(loc="center left", fontsize=6.8, bbox_to_anchor=(0.0, 0.33))
        else:
            ax.legend(loc="lower right", fontsize=6.8, bbox_to_anchor=(1.0, 0.2))
    panel(axs[0], "a")
    panel(axs[1], "b")
    fig.tight_layout()
    save(fig, "fig7_ibm_encounter", lang)


FIGS = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6, 7: fig7}
if __name__ == "__main__":
    which = [int(a) for a in sys.argv[1:]] or list(FIGS)
    for k in which:
        for lang in ("en", "zh"):
            try:
                FIGS[k](lang)
            except FileNotFoundError as e:
                print(f"fig{k}: missing data ({e}); skipped")
