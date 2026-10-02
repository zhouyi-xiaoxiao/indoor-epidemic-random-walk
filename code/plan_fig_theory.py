"""plan_fig_theory.py -- mean-field (next-generation) figures of the article, drawn for the canonical
scenes (code/bsc_sim/scenes/*.json).  Each figure is written to ../figures as <name>_en.pdf/png and
<name>_zh.pdf/png.

Needs ../data/plan_theory_checks.json and ../data/plan_theory_sweeps.npz (run plan_theory_checks.py first)
and data/bsc_theory/fig2_room_size_data.json.

    python plan_fig_theory.py            # all five figures (about 1 minute)

Figures
    fig_r0_mobility            R0 against the mobility scale D0: office scene (bounds, expansions, two
                               alternative layouts) and all four scenes relative to the well-mixed value
    fig_room_size_leak         closed room vs leaky rooms; corridor with absorbing ends
    fig_hotspot_maps           office scene: zones, q, prevalence map u, incidence map w, elasticity map e
    fig_meanfield_interventions  corridor mobility, partitions as movement barriers, targeted ventilation
    fig_theorem_checks         random-instance checks of the theorems; modal truncation on the office scene
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figstyle import C, bilingual, tr  # noqa: E402  (sets the Agg backend)

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
from matplotlib.colors import ListedColormap, LogNorm  # noqa: E402

ROOT = HERE.parents[0]
DATA = HERE.parent / "data"
RES = json.loads((DATA / "plan_theory_checks.json").read_text())
NPZ = np.load(DATA / "plan_theory_sweeps.npz")
ROOM = json.loads((ROOT / "data/bsc_theory/fig2_room_size_data.json").read_text())
BETA, GAMMA = RES["beta"], RES["gamma"]
OFF = RES["office"]
D_WALK = 2.4e3          # cell^2/day: walking 5 % of the time (see Section 2 of the article)


def panel(ax, letter, title):
    ax.set_title(f"({letter}) {title}", loc="left")


# ------------------------------------------------------------------------------------------------
def fig_r0_mobility(lang):
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2))
    ax = axs[0]
    D = NPZ["office_Dgrid"]
    lo, hi = OFF["lower_bound"], OFF["upper_bound"]
    ax.axvspan(D_WALK, D[-1], color="0.9", zorder=0)
    ax.semilogx(D, NPZ["layout_O1_R0_sweep"], color=C["grey"], lw=0.9, ls="-",
                label=tr(lang, "two other layouts, same zone fractions", "分区比例相同的另两种布局"))
    ax.semilogx(D, NPZ["layout_O2_R0_sweep"], color=C["grey"], lw=0.9, ls="-")
    ax.semilogx(D, NPZ["office_R0_sweep"], color=C["blue"], lw=1.8, label=tr(lang, r"office scene, $R_0=\rho(K)$", r"办公室场景，$R_0=\rho(K)$"))
    Dl = D[D >= 0.45]
    ax.semilogx(Dl, lo + OFF["largeD_coefficient_betaC_over_qbar"] / Dl, color=C["red"], lw=1.0, ls=":",
                label=tr(lang, "fast- and slow-mixing expansions", "快混合与慢混合展开"))
    Ds = D[D <= 2.5]
    ax.semilogx(Ds, hi * (1 - Ds * OFF["smallD_slope_muZ0_over_gamma"]), color=C["red"], lw=1.0, ls=":")
    for yv in (lo, hi):
        ax.axhline(yv, color=C["black"], ls="--", lw=0.8)
    ax.text(D[0] * 1.5, hi + 0.02, rf"$\beta q_{{\max}}/\gamma={hi:.2f}$", ha="left", va="bottom", fontsize=8)
    ax.text(D[0] * 1.5, lo - 0.03, rf"$\beta\langle q\rangle/\gamma={lo:.2f}$", ha="left", va="top", fontsize=8)
    ax.axvline(1.0, color=C["grey"], lw=0.8, ls="-.")
    ax.text(1.25, 4.49, tr(lang, "$D_0=1$\n(default)", "$D_0=1$\n（默认值）"), color="0.35", fontsize=7.5, va="center", ha="left")
    ax.text(1.5e4, 4.3, tr(lang, "walking\n$\\geq$5% of\nthe time", "行走时间\n占比$\\geq$5%"), fontsize=7.5, ha="center", va="center", color="0.3")
    ax.set_ylim(4.05, 5.5)
    ax.set_xlim(D[0], D[-1])
    ax.set_xlabel(tr(lang, r"mobility scale $D_0$ (cell$^2$ day$^{-1}$)", r"迁移率尺度 $D_0$（格$^2$ 天$^{-1}$）"))
    ax.set_ylabel(tr(lang, r"basic reproduction number $R_0$", r"基本再生数 $R_0$"))
    ax.legend(loc="lower left", fontsize=7)
    panel(ax, "a", tr(lang, "office scene", "办公室场景"))

    ax = axs[1]
    Dg = NPZ["scenes_Dgrid"]
    names = dict(office=tr(lang, "office", "办公室"), supermarket=tr(lang, "supermarket", "超市"),
                 classroom=tr(lang, "classroom", "教室"), metro=tr(lang, "metro carriage", "地铁车厢"))
    cols = dict(office=C["blue"], supermarket=C["green"], classroom=C["orange"], metro=C["red"])
    ax.axvspan(D_WALK, Dg[-1], color="0.9", zorder=0)
    for k in ("office", "supermarket", "classroom", "metro"):
        low = RES["scenes"][k]["lower"]
        ax.semilogx(Dg, NPZ[f"{k}_R0_sweep_1m"] / low, color=cols[k], lw=1.6, label=names[k])
        if k in ("classroom", "metro"):
            ax.semilogx(Dg, NPZ[f"{k}_R0_sweep_samecell"] / low, color=cols[k], lw=1.0, ls="--")
    ax.plot([], [], color="0.3", lw=1.0, ls="--", label=tr(lang, "same-cell kernel", "同格接触核"))
    ax.axhline(1.0, color=C["black"], ls="--", lw=0.8)
    ax.axvline(1.0, color=C["grey"], lw=0.8, ls="-.")
    ax.set_xlim(Dg[0], Dg[-1])
    ax.set_ylim(0.95, 2.3)
    ax.set_xlabel(tr(lang, r"mobility scale $D_0$ (cell$^2$ day$^{-1}$)", r"迁移率尺度 $D_0$（格$^2$ 天$^{-1}$）"))
    ax.set_ylabel(tr(lang, r"$R_0\,/\,(\beta\langle q\rangle/\gamma)$", r"$R_0\,/\,(\beta\langle q\rangle/\gamma)$"))
    ax.legend(loc="upper right", fontsize=7)
    panel(ax, "b", tr(lang, "four scenes, relative to well mixed", "四个场景，相对于均匀混合值"))
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------
def fig_room_size_leak(lang):
    d = {k: np.array(v) for k, v in ROOM.items()}
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.1))
    ax = axs[0]
    Lc = np.linspace(2.5, 42, 300)
    ax.plot(d["Ls"], d["refl"], "o", color=C["blue"], ms=4.5, label=tr(lang, "closed room, reflecting walls (any $D$)", "封闭房间，反射壁（任意 $D$）"))
    ax.plot(d["Ls"], d["door"], "s-", color=C["green"], ms=3.5, lw=1.2, label=tr(lang, "one door cell (leak)", "一个门格（泄漏）"))
    ax.plot(d["Ls"], d["wall"], "^-", color=C["orange"], ms=3.5, lw=1.2, label=tr(lang, "one side fully open", "一侧完全敞开"))
    ax.set_xscale("log")
    ax.set_xticks([3, 5, 10, 20, 40])
    ax.set_xticklabels(["3", "5", "10", "20", "40"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0.15, 1.05)
    ax.set_xlabel(tr(lang, r"room side $L$ (cells)", r"房间边长 $L$（格）"))
    ax.set_ylabel(r"$R_0\,/\,(\beta q/\gamma)$")
    ax.legend(loc="lower right", fontsize=7)
    panel(ax, "a", tr(lang, r"square room, uniform $q$, $D=0.5$", r"方形房间，均匀 $q$，$D=0.5$"))

    ax = axs[1]
    Lc = np.linspace(2.5, 110, 400)
    ax.plot(d["L1"], d["r1"], "o", color=C["blue"], ms=4.5, label=tr(lang, r"$\rho(K)$, corridor with absorbing ends", r"$\rho(K)$：两端吸收的走廊"))
    ax.plot(Lc, BETA / (GAMMA + 2 * (1 - np.cos(np.pi / (Lc + 1)))), color=C["blue"], lw=1.1,
            label=r"$\beta q/[\gamma+2w(1-\cos\frac{\pi}{L+1})]$")
    ax.plot(Lc, BETA / (GAMMA + np.pi**2 / Lc**2), color=C["red"], ls=":", lw=1.4, label=r"$\beta q/(\gamma+\pi^2D/L^2)$")
    ax.axhline(BETA / GAMMA, color=C["black"], ls="--", lw=0.8)
    ax.text(110, BETA / GAMMA + 0.04, tr(lang, r"closed room: $\beta q/\gamma$", r"封闭房间：$\beta q/\gamma$"), ha="right", va="bottom", fontsize=8)
    ax.set_xscale("log")
    ax.set_xticks([3, 10, 30, 100])
    ax.set_xticklabels(["3", "10", "30", "100"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0, 3.95)
    ax.set_xlabel(tr(lang, r"corridor length $L$ (cells)", r"走廊长度 $L$（格）"))
    ax.set_ylabel(tr(lang, r"room reproduction number $R_0^{\rm room}$", r"房间再生数 $R_0^{\rm room}$"))
    ax.legend(loc="lower right", fontsize=7)
    panel(ax, "b", tr(lang, r"corridor, absorbing ends, $w=D=1$, $q=1$", r"走廊，两端吸收，$w=D=1$，$q=1$"))
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------
def fig_hotspot_maps(lang):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    plt.rcParams.update({"axes.spines.left": False, "axes.spines.bottom": False})
    zi = NPZ["office_zone_code"]
    fig, axs = plt.subplots(2, 3, figsize=(7.6, 4.3))
    zcols = ["#56B4E9", "#F0E442", "#D55E00", "#009E73", "#333333"]
    znames = [tr(lang, "workstations", "工位区"), tr(lang, "corridor", "走廊"), tr(lang, "meeting room", "会议室"),
              tr(lang, "kitchen", "茶水间"), tr(lang, "barrier", "障碍物")]
    barrier = np.ma.masked_where(zi != 4, np.ones_like(zi))

    def hide(ax):
        ax.set_xticks([])
        ax.set_yticks([])

    ax = axs[0, 0]
    ax.imshow(zi, cmap=ListedColormap(zcols), vmin=-0.5, vmax=4.5, interpolation="nearest")
    hide(ax)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=zcols[i], ec="0.3", lw=0.4) for i in range(5)]
    ax.legend(handles, znames, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=3, fontsize=6.5, handlelength=1.0, columnspacing=0.8)
    panel(ax, "a", tr(lang, "zones (20 x 13 cells)", "分区（20 x 13 格）"))

    def heat(ax, arr, letter, title, cmap, norm=None, label=None, vmin=None, vmax=None):
        im = ax.imshow(np.ma.masked_invalid(arr), cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.imshow(barrier, cmap=ListedColormap(["#333333"]), interpolation="nearest")
        hide(ax)
        cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.05)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=7)
        if label:
            cb.set_label(label, fontsize=7.5)
        panel(ax, letter, title)

    heat(axs[0, 1], NPZ["office_q"], "b", tr(lang, r"infection efficiency $q$", r"感染效率 $q$"), "YlOrRd", vmin=0.7, vmax=1.6)
    norm = LogNorm(vmin=0.03, vmax=30)
    lab = tr(lang, "relative to uniform", "相对于均匀分布")
    heat(axs[0, 2], NPZ["office_u_D1"], "c", tr(lang, r"prevalence map $n\,u$, $D_0=1$", r"现患分布 $n\,u$，$D_0=1$"), "viridis", norm=norm, label=lab)
    heat(axs[1, 0], NPZ["office_u_D100"], "d", tr(lang, r"prevalence map $n\,u$, $D_0=100$", r"现患分布 $n\,u$，$D_0=100$"), "viridis", norm=norm, label=lab)
    heat(axs[1, 1], NPZ["office_w_D1"], "e", tr(lang, r"incidence map $n\,w$, $D_0=1$", r"发病分布 $n\,w$，$D_0=1$"), "viridis", norm=norm, label=lab)
    heat(axs[1, 2], NPZ["office_e_D1"], "f", tr(lang, r"elasticity map $n\,e$, $D_0=1$", r"弹性分布 $n\,e$，$D_0=1$"), "viridis", norm=norm, label=lab)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------
_CORR = {}


def _corridor_sweep():
    """R0 of the office scene when the corridor mobility is multiplied by m (symmetric rule and the
    departure-site reading).  Computed here (a few seconds) and cached in ../data/plan_corridor_sweep.json."""
    if _CORR:
        return _CORR
    cache = DATA / "plan_corridor_sweep.json"
    if cache.exists():
        _CORR.update(json.loads(cache.read_text()))
        return _CORR
    sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))
    sys.path.insert(0, str(ROOT / "code/bsc_theory/scripts"))
    from bsc_sim import scenes as sc
    from bsc_sim import theory as th
    import ngm_lattice as ng
    import plan_theory_checks as pc
    s1 = sc.load_scene("office", D0=1.0)
    lat = ng.Lattice(s1.access)
    mult = np.logspace(-2, 1, 25)
    sym, dep = [], []
    for m in mult:
        D = s1.D.copy()
        D[s1.site_zone == "C"] *= m
        sym.append(pc.r0_only(th.generator(s1, D).toarray(), s1.q))
        dep.append(pc.r0_only(ng.generator_departure(lat, D).toarray(), s1.q))
    _CORR.update(mult=mult.tolist(), symmetric=sym, departure=dep)
    cache.write_text(json.dumps(_CORR, indent=1))
    return _CORR


def fig_meanfield_interventions(lang):
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
    ax = axs[0]
    cs = _corridor_sweep()
    ax.semilogx(cs["mult"], cs["symmetric"], color=C["blue"], lw=1.6, label=tr(lang, "symmetric rule", "对称跳跃规则"))
    ax.semilogx(cs["mult"], cs["departure"], color=C["orange"], lw=1.4, ls="--", label=tr(lang, "departure-site rule", "出发点跳跃规则"))
    ax.axvline(1.0, color=C["grey"], lw=0.8, ls="-.")
    ax.set_xlabel(tr(lang, "corridor mobility multiplier", "走廊迁移率倍数"))
    ax.set_ylabel(r"$R_0$")
    ax.set_ylim(4.9, 5.45)
    ax.legend(loc="lower left", fontsize=6.5)
    panel(ax, "a", tr(lang, "slower corridors", "降低走廊迁移率"))

    ax = axs[1]
    P = OFF["partitions"]
    base = OFF["D0=1"]["R0"]
    xs = [0] + [100 * P[k]["mean_closed_edges"] / P["n_WW_edges"] for k in ("walls_10pct", "walls_30pct", "walls_50pct", "walls_100pct")]
    ys = [0] + [P[k]["change_pct"] for k in ("walls_10pct", "walls_30pct", "walls_50pct", "walls_100pct")]
    es = [0] + [100 * P[k]["sd"] / base for k in ("walls_10pct", "walls_30pct", "walls_50pct", "walls_100pct")]
    ax.errorbar(xs, ys, yerr=es, fmt="o-", color=C["blue"], ms=4, capsize=2, label=tr(lang, "walls between workstation cells", "工位格之间加隔断"))
    ob = P["obstacles_39_workstations"]
    ax.errorbar([25], [ob["change_pct"]], yerr=[100 * ob["sd"] / base], fmt="s", color=C["green"], ms=5, capsize=2,
                label=tr(lang, "39 workstation cells made barriers", "39 个工位格改为障碍"))
    ax.axhline(0, color=C["black"], lw=0.8)
    ax.set_ylim(-0.05, 0.16)
    ax.set_xlabel(tr(lang, "workstation-workstation edges closed (%)", "被封闭的工位间连边（%）"))
    ax.set_ylabel(tr(lang, r"change of $R_0$ (%)", r"$R_0$ 的变化（%）"))
    ax.legend(loc="upper left", fontsize=6.5)
    panel(ax, "b", tr(lang, "partitions, $D_0=1$", "隔断，$D_0=1$"))

    ax = axs[2]
    T = OFF["targeted_ventilation"]["D0=1"]
    rows = T["rows"]
    fr = [0] + [100 * r["fraction"] for r in rows]
    b = T["baseline"]
    ax.plot(fr, [b] + [r["adaptive"] for r in rows], "-", color=C["blue"], lw=1.6, label=tr(lang, "elasticity-guided, recomputed", "按弹性分布、逐步重算"))
    ax.errorbar(fr, [b] + [r["highest_q_mean"] for r in rows], yerr=[0] + [r["highest_q_sd"] for r in rows], fmt="s--", color=C["green"], ms=3.5, lw=1.0, capsize=2,
                label=tr(lang, "highest $q$ first", "$q$ 最大者优先"))
    ax.errorbar(fr, [b] + [r["random_mean"] for r in rows], yerr=[0] + [r["random_sd"] for r in rows], fmt="^:", color=C["orange"], ms=3.5, lw=1.0, capsize=2,
                label=tr(lang, "random cells", "随机选格"))
    ax.plot(fr, [b] + [r["one_shot"] for r in rows], "-.", color=C["grey"], lw=1.0, label=tr(lang, "baseline ranking, not recomputed", "按基线排序、不重算"))
    ax.set_xlabel(tr(lang, r"area with $q\to q/2$ (%)", r"$q\to q/2$ 的面积占比（%）"))
    ax.set_ylabel(r"$R_0$")
    ax.set_ylim(2.4, 6.6)
    ax.legend(loc="upper right", fontsize=6.5)
    panel(ax, "c", tr(lang, "targeted ventilation, $D_0=1$", "定向通风，$D_0=1$"))
    fig.tight_layout()
    return fig

# ------------------------------------------------------------------------------------------------
def fig_theorem_checks(lang):
    """Random-instance checks of the threshold theorem and of the bounds (data of the theory part:
    data/bsc_theory/T1_points.npy, T2_positions.npy) and convergence of the modal (Galerkin)
    formula on the office scene (data/plan_theory_checks.json)."""
    RESD = ROOT / "code/bsc_theory/results"
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
    fam = [tr(lang, "symmetric", "对称"), tr(lang, "departure-site", "出发点规则"), tr(lang, "non-reversible", "不可逆"), tr(lang, "dense random", "稠密随机")]
    cols = [C["blue"], C["orange"], C["green"], C["purple"]]
    pts = np.load(RESD / "T1_points.npy")
    ax = axs[0]
    ax.axhspan(0, 1e3, xmin=0.5, color="0.93", zorder=0)
    ax.axhspan(-1e3, 0, xmax=0.5, color="0.93", zorder=0)
    for k in range(4):
        m = pts[:, 2] == k
        ax.plot(np.log10(pts[m, 0]), np.arcsinh(pts[m, 1] * 5) / np.log(10), ".", color=cols[k], ms=3.5, label=fam[k])
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(-2.6, 2.6)
    yt = np.array([-1, -0.3, 0, 0.3, 1, 3, 10, 30])
    ax.set_yticks(np.arcsinh(yt * 5) / np.log(10))
    ax.set_yticklabels([f"{v:g}" for v in yt])
    ax.set_ylim(np.arcsinh(-1.1 * 5) / np.log(10), np.arcsinh(60 * 5) / np.log(10))
    ax.set_xlabel(r"$\log_{10}R_0$")
    ax.set_ylabel(tr(lang, r"$\lambda_1/\gamma$ (asinh scale)", r"$\lambda_1/\gamma$（asinh 标度）"))
    ax.legend(loc="upper left", fontsize=6.5, handletextpad=0.1)
    panel(ax, "a", tr(lang, r"sign of $\lambda_1$ and of $R_0-1$", r"$\lambda_1$ 与 $R_0-1$ 的符号"))

    ax = axs[1]
    pos = np.load(RESD / "T2_positions.npy")
    bins = np.linspace(0, 1, 26)
    bottom = np.zeros(len(bins) - 1)
    for k in range(4):
        h, _ = np.histogram(pos[pos[:, 0] == k, 1], bins)
        ax.bar(bins[:-1], h, width=np.diff(bins), align="edge", bottom=bottom, color=cols[k], lw=0)
        bottom += h
    ax.set_xlim(-0.05, 1.05)
    ax.axvline(0, color="k", ls="--", lw=0.9)
    ax.axvline(1, color="k", ls="--", lw=0.9)
    ax.set_xlabel(r"$(R_0-\beta\langle q\rangle_\pi/\gamma)\,/\,(\beta q_{\max}/\gamma-\beta\langle q\rangle_\pi/\gamma)$", fontsize=7)
    ax.set_ylabel(tr(lang, "number of random instances", "随机算例数"))
    panel(ax, "b", tr(lang, f"bounds: {len(pos)} instances", f"上下界：{len(pos)} 个算例"))

    ax = axs[2]
    g = OFF["galerkin_R0_by_modes"]
    Ns = sorted(int(k) for k in g)
    ax.semilogx(Ns, [g[str(k)] for k in Ns], "o-", color=C["blue"], ms=4, label=tr(lang, "modal formula, $m$ modes", "模态公式，$m$ 个模态"))
    ax.axhline(OFF["D0=1"]["R0"], color=C["blue"], ls="--", lw=0.9, label=tr(lang, r"$R_0=\rho(K)$", r"$R_0=\rho(K)$"))
    ax.plot([1], [OFF["diagonal_modal_max"]], "*", color=C["red"], ms=9, label=tr(lang, "largest diagonal term (lower bound)", "最大对角项（下界）"))
    ax.set_xlabel(tr(lang, "number of modes $m$", "模态数 $m$"))
    ax.set_ylabel(r"$R_0^{(m)}$")
    ax.set_ylim(4.15, 5.2)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.08), fontsize=6, handlelength=1.5)
    panel(ax, "c", tr(lang, "office scene, $D_0=1$", "办公室场景，$D_0=1$"))
    fig.tight_layout()
    return fig


FIGS = dict(fig_r0_mobility=fig_r0_mobility, fig_room_size_leak=fig_room_size_leak,
            fig_hotspot_maps=fig_hotspot_maps, fig_meanfield_interventions=fig_meanfield_interventions,
            fig_theorem_checks=fig_theorem_checks)

if __name__ == "__main__":
    want = sys.argv[1:] or list(FIGS)
    for name in want:
        print(bilingual(FIGS[name], name), flush=True)
