"""Figures of the validation section, each in English (_en) and Simplified
Chinese (_zh), vector PDF + PNG, drawn from the stored results only."""
import numpy as np

import _common as C
from valmod import plotting as PL
from valmod.plotting import C as COL, tr

plt = PL.plt
MC = {"M0": COL["grey"], "M1": COL["orange"], "M2": COL["blue"]}
EV = ["T1", "T2", "T5", "C1"]


def evname(lang, ev):
    return {"T1": tr(lang, "T1 Zhejiang bus", "T1 浙江大巴"),
            "T2": tr(lang, "T2 Hunan coach", "T2 湖南长途客车"),
            "T5": tr(lang, "T5 flight VN54", "T5 航班VN54"),
            "C1": tr(lang, "C1 Marin classroom", "C1 马林县教室"),
            "T4": tr(lang, "T4 trains", "T4 高铁"),
            "R1": tr(lang, "R1 restaurant", "R1 广州餐厅"),
            "T3": tr(lang, "T3 minibus", "T3 小巴"),
            "H1": tr(lang, "H1 choir", "H1 合唱团"),
            "O2": tr(lang, "O2 office", "O2 办公室"),
            "R2": tr(lang, "R2 restaurant", "R2 全州餐厅"),
            "C5": tr(lang, "C5 fitness", "C5 健身课")}[ev]


def mname(lang, m):
    return {"M0": tr(lang, "M0 well mixed", "M0 均匀混合"),
            "M1": tr(lang, "M1 (same-site)", "M1（同格点传播）"),
            "M2": tr(lang, "M2 (diffusing aerosol)", "M2（扩散气溶胶）")}[m]


# ---------------------------------------------------------------------------
def fig_calibration(lang):
    t4 = C.load_json("01_cal_T4_primary.json")
    r1 = C.load_json("01_cal_R1_primary.json")
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 6.6))
    c = t4["cells"]
    x = np.arange(len(c["dr"]))
    pct, lo, hi = np.array(c["pct"]), np.array(c["lo"]), np.array(c["hi"])
    used = np.array(t4["mask"])
    a = ax[0, 0]
    lo_plot = np.maximum(lo, 0.012)
    a.errorbar(x[used], pct[used], yerr=[pct[used] - lo_plot[used], hi[used] - pct[used]], fmt="o", ms=3.5,
               color="k", lw=0.8, capsize=1.5, label=tr(lang, "observed (95% CI)", "观测值 (95% CI)"))
    a.errorbar(x[~used], pct[~used], yerr=[pct[~used] - lo_plot[~used], hi[~used] - pct[~used]], fmt="o", ms=3.5,
               mfc="w", color="k", lw=0.8, capsize=1.5,
               label=tr(lang, "adjacent seat (excluded)", "相邻座位 (不参与拟合)"))
    for m in ("M0", "M1", "M2"):
        a.plot(x, 100 * np.array(t4[m]["p"]), "-" if m != "M1" else "--", color=MC[m], label=mname(lang, m))
    a.set_yscale("log")
    a.set_xticks(x)
    a.set_xticklabels([f"{r},{cc}" for r, cc in zip(c["dr"], c["dc"])], rotation=90, fontsize=6.5)
    a.set_xlabel(tr(lang, "seat offset (rows apart, seats apart)", "座位偏移 (相隔排数, 相隔座位数)"))
    a.set_ylabel(tr(lang, "attack rate (%)", "感染率 (%)"))
    a.set_title(tr(lang, "(a) Train cohort: calibration fit", "(a) 高铁队列：标定拟合"))
    a.legend(fontsize=6.8, ncol=2)
    a = ax[0, 1]
    logD = np.array(t4["logD"])
    a.plot(logD, t4["M1"]["profile_logL"], "--", color=MC["M1"], label=mname(lang, "M1"))
    a.plot(logD, t4["M2"]["marg_logL"], color=MC["M2"], label=mname(lang, "M2"))
    a.axhline(t4["M0"]["logL"], color=MC["M0"], label=mname(lang, "M0"))
    a.set_ylim(-135, -65)
    a.set_xlabel(tr(lang, r"$\log_{10}$ mixing coefficient $D$ (m$^2$/h)", r"混合系数 $D$ 的 $\log_{10}$ (m$^2$/h)"))
    a.set_ylabel(tr(lang, "profile log-likelihood", "剖面对数似然"))
    a.set_title(tr(lang, "(b) Train cohort: likelihood of $D$", "(b) 高铁队列：$D$ 的似然"))
    a.legend(fontsize=7)
    a = ax[1, 0]
    names = r1["tables"]
    n = np.array(r1["n"])
    grp = {"TB": [names.index("TB")], "TC": [names.index("TC")], "T18": [names.index("T18")],
           tr(lang, "13 remote\ntables", "13张远桌"): [i for i, t in enumerate(names) if t not in ("TB", "TC", "T18")]}
    xs = np.arange(len(grp))
    wbar = 0.25
    for j, m in enumerate(("M0", "M1", "M2")):
        p = np.array(r1[m]["p"])
        vals = [np.sum(p[idx] * n[idx]) / np.sum(n[idx]) for idx in grp.values()]
        a.bar(xs + (j - 1) * wbar, np.array(vals) * 100, wbar, color=MC[m], label=mname(lang, m))
    obs_lo = [100 * 1 / 4, 100 * 1 / 7, 0, 0]
    obs_hi = [100 * 3 / 4, 100 * 2 / 7, 0, 0]
    for xi, l, h in zip(xs, obs_lo, obs_hi):
        a.plot([xi - 0.4, xi + 0.4], [h, h], color="k", lw=1.2)
        a.plot([xi - 0.4, xi + 0.4], [l, l], color="k", lw=1.2, ls=":")
    a.plot([], [], color="k", lw=1.2, label=tr(lang, "observed, upper bound", "观测上界"))
    a.plot([], [], color="k", lw=1.2, ls=":", label=tr(lang, "observed, lower bound", "观测下界"))
    a.set_xticks(xs)
    a.set_xticklabels(list(grp.keys()))
    a.set_ylabel(tr(lang, "attack rate (%)", "感染率 (%)"))
    a.set_title(tr(lang, "(c) Guangzhou restaurant: fit per table group", "(c) 广州餐厅：各桌组拟合"))
    a.legend(fontsize=6.8)
    a = ax[1, 1]
    logD = np.array(r1["logD"])
    a.plot(logD, r1["M1"]["profile_logL"], "--", color=MC["M1"], label=mname(lang, "M1"))
    a.plot(logD, r1["M2"]["profile_logL"], color=MC["M2"], label=mname(lang, "M2"))
    a.axhline(r1["M0"]["logL"], color=MC["M0"], label=mname(lang, "M0"))
    a.set_xlabel(tr(lang, r"$\log_{10}$ mixing coefficient $D$ (m$^2$/h)", r"混合系数 $D$ 的 $\log_{10}$ (m$^2$/h)"))
    a.set_ylabel(tr(lang, "profile log-likelihood", "剖面对数似然"))
    a.set_title(tr(lang, "(d) Guangzhou restaurant: likelihood of $D$", "(d) 广州餐厅：$D$ 的似然"))
    a.legend(fontsize=7)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def fig_holdout_pmf(lang):
    sh = C.load_json("03_holdout_shape.json")
    pred = C.load_json("02_predictive.json")
    ex = C.load_json("02b_m1_exact.json")["events"]
    fig, axs = plt.subplots(2, 2, figsize=(9.2, 6.4))
    for a, ev in zip(axs.ravel(), EV):
        e = pred["events"][ev]
        r = sh["events"][ev]
        kmax = min(e["n_near"], e["K"])
        k = np.arange(kmax + 1)
        a.bar(k, np.array(e["M0"]["pmf"])[:kmax + 1], color=MC["M0"], alpha=0.45, width=0.9, label=mname(lang, "M0"))
        unatt = False
        if r["M1_uses_exact_simulator"]:
            pge = float(np.mean([q["p_total_ge_K"] for q in ex[ev]["rows"]]))
            unatt = pge < 0.01
            p1 = np.array(ex[ev]["pmf_sim"])
            lab1 = mname(lang, "M1")
        else:
            p1 = np.array(e["M1"]["pmf"])
            lab1 = mname(lang, "M1")
        if not unatt:
            a.plot(k, p1[:kmax + 1], "s--", ms=3, color=MC["M1"], label=lab1)
        a.plot(k, np.array(e["M2"]["pmf"])[:kmax + 1], "o-", ms=3, color=MC["M2"], label=mname(lang, "M2"))
        ko = r["obs"]["near_k"]
        a.axvline(ko, color="k", ls=":", lw=1.3, label=tr(lang, "observed", "观测值"))
        lo = max(0, e["K"] - e["n_far"])
        a.set_xlim(lo - 0.7, kmax + 0.7)
        txt = "\n".join(
            (f"M1: " + tr(lang, "total unattainable", "总数不可达")) if (m == "M1" and unatt)
            else f"{m}: p = {r[m]['p_two_sided']:.3g}" for m in ("M0", "M1", "M2"))
        a.text(0.03, 0.97, txt, transform=a.transAxes, va="top", fontsize=7.5)
        a.set_title(f"{evname(lang, ev)}: {e['K']} " + tr(lang, "cases", "例") +
                    f" ({e['n_near']} " + tr(lang, "near", "近区") + f", {e['n_far']} " + tr(lang, "far", "远区") + ")")
        a.set_xlabel(tr(lang, "number of cases in the near stratum", "近区病例数"))
        a.set_ylabel(tr(lang, "predictive probability", "预测概率"))
    h, l = axs[1, 0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


# ---------------------------------------------------------------------------
def fig_forest(lang):
    sh = C.load_json("03_holdout_shape.json")
    ex = C.load_json("02b_m1_exact.json")["events"]
    fig, a = plt.subplots(figsize=(7.6, 5.6))
    y = 0
    yt, yl = [], []
    items = [(ev, sh["events"][ev], evname(lang, ev)) for ev in EV]
    for key, r in sh["secondary"].items():
        lab = {"T1_rows6-10": tr(lang, "T1 bus, rows 6-10 (secondary)", "T1 大巴 第6-10排 (次要)"),
               "T2_driver_side": tr(lang, "T2 coach, driver side (secondary)", "T2 客车 司机侧 (次要)")}[key]
        items.append((key, r, lab))
    for key, r, lab in items:
        ro = r["M0"]["rr_obs"]
        se = r["M0"]["se_log_rr_obs"]
        a.errorbar(ro, y, xerr=[[ro - ro * np.exp(-1.96 * se)], [ro * np.exp(1.96 * se) - ro]], fmt="D", color="k",
                   ms=5, capsize=2, label=tr(lang, "observed (95% CI)", "观测值 (95% CI)") if y == 0 else None)
        for j, m in enumerate(("M0", "M1", "M2")):
            x = r[m]
            yy = y - 0.2 * (j + 1)
            if m == "M1" and key in ex and np.mean([q["p_total_ge_K"] for q in ex[key]["rows"]]) < 0.01:
                a.text(1.05, yy, tr(lang, "M1: observed total unattainable", "M1：观测总数不可达"), color=MC[m],
                       fontsize=7, va="center")
                continue
            lo, hi = x["rr_pred_90"]
            med = min(max(x["rr_pred"], lo), hi)
            a.plot([lo, hi], [yy, yy], color=MC[m], lw=2.2,
                   label=(mname(lang, m) + tr(lang, " (90% predictive interval)", " (90%预测区间)")) if y == 0 else None)
            a.plot(med, yy, "|", color=MC[m], ms=9, mew=2)
        yt.append(y - 0.3)
        yl.append(lab)
        y -= 1
    a.axvline(1, color="k", lw=0.6, ls=":")
    a.set_xscale("log")
    a.set_yticks(yt)
    a.set_yticklabels(yl)
    a.set_xlabel(tr(lang, "risk ratio, near vs far stratum", "近区与远区的风险比"))
    h, l = a.get_legend_handles_labels()
    for m in ("M1",):
        if mname(lang, m) + tr(lang, " (90% predictive interval)", " (90%预测区间)") not in l:
            h.append(plt.Line2D([], [], color=MC[m], lw=2.2))
            l.append(mname(lang, m) + tr(lang, " (90% predictive interval)", " (90%预测区间)"))
    a.legend(h, l, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def fig_scores(lang):
    sh = C.load_json("03_holdout_shape.json")
    fig, a = plt.subplots(figsize=(6.6, 3.6))
    labs = [evname(lang, ev) for ev in EV] + [tr(lang, "pooled", "合计")]
    x = np.arange(5)
    for j, m in enumerate(("M1", "M2")):
        v = [sh["events"][ev][m]["delta_vs_M0"] for ev in EV] + [sh["pooled"][m]["delta"]]
        a.bar(x + (j - 0.5) * 0.36, v, 0.36, color=MC[m], label=mname(lang, m))
    a.axhline(0, color="k", lw=0.7)
    a.set_xticks(x)
    a.set_xticklabels(labs, fontsize=7.5)
    a.set_ylabel(tr(lang, "log-score difference against M0", "相对于M0的对数评分差"))
    a.text(0.02, 0.04, tr(lang, f"pooled: M2 p = {sh['pooled']['M2']['p_under_M0']:.2g} (exact, under M0)",
                          f"合计：M2 p = {sh['pooled']['M2']['p_under_M0']:.2g} (M0下精确检验)"),
           transform=a.transAxes, fontsize=7.5)
    a.legend(fontsize=7.5)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def fig_callcentre(lang):
    o1 = C.load_json("06_callcentre.json")
    fig, a = plt.subplots(figsize=(6.4, 3.6))
    bins = np.linspace(-3.5, 14.5, 46)
    for m in ("M0", "M1", "M2"):
        a.hist(o1[m]["z_sample"], bins=bins, density=True, histtype="step", color=MC[m], lw=1.5,
               label=mname(lang, m))
    a.axvline(o1["observed"]["z"], color="k", ls=":", lw=1.4,
              label=tr(lang, f"observed z = {o1['observed']['z']:.2f}", f"观测 z = {o1['observed']['z']:.2f}"))
    a.set_xlabel(tr(lang, "join-count z (case-case adjacent desk pairs, north wing)",
                    "邻接计数 z (北翼相邻工位的病例-病例对)"))
    a.set_ylabel(tr(lang, "predictive density", "预测密度"))
    a.legend(fontsize=7.5)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def fig_m1_ceiling(lang):
    ex = C.load_json("02b_m1_exact.json")["events"]
    fig, a = plt.subplots(figsize=(6.2, 3.6))
    x = np.arange(4)
    obs = [ex[ev]["K"] for ev in EV]
    med, lo, hi = [], [], []
    for ev in EV:
        t = np.array([q["mean_total_sim"] for q in ex[ev]["rows"]])
        med.append(np.median(t))
        lo.append(t.min())
        hi.append(t.max())
    a.bar(x - 0.2, obs, 0.4, color="k", alpha=0.75, label=tr(lang, "observed secondary cases", "观测二代病例数"))
    a.bar(x + 0.2, med, 0.4, color=MC["M1"],
          label=tr(lang, "M1, exact simulator: mean total at the tuned beta,\nor at the largest beta when the total cannot be reached\n(median and range over posterior draws)",
                   "M1 精确模拟：调定β下的平均总数；\n无法达到时取最大β下的平均总数\n(后验抽样的中位数与范围)"))
    a.errorbar(x + 0.2, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt="none", color="k", capsize=3, lw=0.9)
    a.set_xticks(x)
    a.set_xticklabels([evname(lang, ev) for ev in EV], fontsize=7.5)
    a.set_ylabel(tr(lang, "number of secondary cases", "二代病例数"))
    a.legend(fontsize=7)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def fig_level(lang):
    lv = C.load_json("04_level_controls.json")
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.8))
    a = ax[0]
    evs = ["T2", "T3", "C5", "O2", "R2", "H1"]
    for i, ev in enumerate(evs):
        x = lv["implied"][ev]
        a.errorbar(x["E_implied"], i, xerr=[[x["E_implied"] - x["E_range_all"][0]], [x["E_range_all"][1] - x["E_implied"]]],
                   fmt="o", color=COL["blue"], capsize=2, ms=4)
    et = lv["E_T4"]["q"]
    a.errorbar(et[1], -1, xerr=[[et[1] - et[0]], [et[2] - et[1]]], fmt="s", color=COL["green"], capsize=2, ms=4)
    a.set_yticks(range(-1, len(evs)))
    a.set_yticklabels([tr(lang, "T4 trains (cohort mean)", "T4 高铁 (队列均值)")] + [evname(lang, e) for e in evs])
    a.set_xscale("log")
    a.set_xlabel(tr(lang, "implied emission rate (quanta/h at 0.5 m$^3$/h)", "反推的释放率 (quanta/h, 呼吸率0.5 m$^3$/h)"))
    a.set_title(tr(lang, "(a) Infectiousness has to be fitted per event", "(a) 传染强度需逐事件拟合"))
    a = ax[1]
    rows = [(tr(lang, "T3 minibus from T2 (M0)", "由T2预测T3 (M0)"), lv["tied_T2_T3"]["M0"], 2),
            (tr(lang, "T3 minibus from T2 (M1)", "由T2预测T3 (M1)"), lv["tied_T2_T3"]["M1"], 2),
            (tr(lang, "T3 minibus from T2 (M2)", "由T2预测T3 (M2)"), lv["tied_T2_T3"]["M2"], 2),
            (tr(lang, "C3 classrooms from T4 (M2, masks)", "由T4预测C3 (M2, 口罩)"), lv["C3"]["M2_masks_0.35"], 5),
            (tr(lang, "C3 classrooms from T4 (M2, no masks)", "由T4预测C3 (M2, 无口罩)"), lv["C3"]["M2_no_masks"], 5),
            (tr(lang, "C3 classrooms from T4 (M1, masks)", "由T4预测C3 (M1, 口罩)"), lv["C3"]["M1_masks_0.35"], 5)]
    for i, (lab, r, ko) in enumerate(rows):
        a.plot(r["pi90"], [i, i], color=COL["blue"], lw=2.2)
        a.plot(r["mean"], i, "|", color=COL["blue"], ms=9, mew=2)
        a.plot(ko, i, "D", color="k", ms=4.5)
    a.set_yticks(range(len(rows)))
    a.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    a.set_xlabel(tr(lang, "number of cases (bar: 90% predictive interval; diamond: observed)",
                    "病例数 (线段：90%预测区间；菱形：观测值)"))
    a.set_title(tr(lang, "(b) Level checks with a defined reference", "(b) 有明确参照的水平检验"))
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def fig_heterogeneity(lang):
    lo = C.load_json("07_loo.json")
    logD = np.array(lo["logD"])
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.7), sharey=True)
    cols = [COL["black"], COL["red"], COL["green"], COL["purple"]]
    for a, cls, ttl in ((ax[0], "vehicle", tr(lang, "(a) Vehicle cabins, M2", "(a) 交通工具舱室, M2")),
                        (ax[1], "room", tr(lang, "(b) Seated rooms, M2", "(b) 就座房间, M2"))):
        for (ev, cur), col in zip(lo["classes"][cls]["M2"]["curves"].items(), cols):
            cur = np.array(cur)
            a.plot(logD, cur - cur.max(), color=col, label=evname(lang, ev))
        a.set_ylim(-12, 0.5)
        a.axhline(-1.92, color="k", lw=0.5, ls=":")
        a.set_xlabel(tr(lang, r"$\log_{10}$ mixing coefficient $D_{\rm air}$ (m$^2$/h)",
                        r"混合系数 $D_{\rm air}$ 的 $\log_{10}$ (m$^2$/h)"))
        a.set_title(ttl)
        a.legend(fontsize=7.5)
    ax[0].set_ylabel(tr(lang, "log-likelihood relative to each event's maximum", "相对于各事件最大值的对数似然"))
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
def fig_sensitivity(lang):
    se = C.load_json("09_sensitivity.json")["variants"]
    tags = list(se.keys())
    P = np.array([[se[t]["events"][ev]["M2"]["p"] for ev in EV] for t in tags])
    R = np.array([[se[t]["events"][ev]["M2"]["rr_pred"] for ev in EV] for t in tags])
    fig, a = plt.subplots(figsize=(6.4, 4.8))
    a.imshow(np.where(P < 0.05, 0.0, 1.0), cmap="RdYlGn", vmin=-0.6, vmax=1.6, aspect="auto")
    for i in range(len(tags)):
        for j in range(4):
            a.text(j, i, f"RR {R[i, j]:.1f}\np = {P[i, j]:.2g}", ha="center", va="center", fontsize=6.5)
    a.set_xticks(range(4))
    a.set_xticklabels([evname(lang, ev) for ev in EV], fontsize=7.5)
    a.set_yticks(range(len(tags)))
    a.set_yticklabels([t.replace("primary (closed form for M1)", tr(lang, "primary", "主分析")) for t in tags], fontsize=7.5)
    a.set_title(tr(lang, "M2: predicted risk ratio and predictive p under each variant",
                   "M2：各敏感性方案下的预测风险比与预测 p 值"))
    for s in a.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    for fn, name in ((fig_calibration, "figV1_calibration"), (fig_holdout_pmf, "figV2_holdout_shape"),
                     (fig_forest, "figV3_risk_ratios"), (fig_scores, "figV4_scores"),
                     (fig_callcentre, "figV5_callcentre"), (fig_m1_ceiling, "figV6_m1_ceiling"),
                     (fig_level, "figV7_level"), (fig_heterogeneity, "figV8_mixing_heterogeneity"),
                     (fig_sensitivity, "figV9_sensitivity")):
        print(PL.bilingual(fn, name), flush=True)
