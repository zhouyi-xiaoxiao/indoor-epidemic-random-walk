#!/usr/bin/env python
"""v2_figures.py -- figures of the second set of tests (Sections 7.5-7.9; Supplementary Sections S9 and S10).

Every figure is drawn from stored results (nothing is recomputed) and is written with English and Chinese
labels (<name>_en.pdf, <name>_zh.pdf) by figstyle.bilingual:

    v2_rr_forest           near/far risk ratios of the eleven further outbreaks, with pooled estimates
    v2_holdout_strata      observed and expected cases by distance stratum, five held-out events
    v2_summary_scores      log-score differences (model minus well mixed, model minus constant ratio), every
                           event of both sets
    v2_contacts_epidemic   SIR outcomes on recorded contacts against model-generated contacts, six settings
    v2_contacts_structure  model / observed ratio of contact statistics, six settings (Supplementary Section S9)
    v2_tracer_kernels      tracer measurements and mixing coefficients
    v2_tracer_predictions  cases in the near stratum under the tracer-fixed kernel and its competitors
    v2_zone                zone-level test: mixing shares, held-out scores by school, dispersion, floors
    v2_D_profiles          log-likelihood of the mixing coefficient for each of the eleven events (Supplementary Section S9)

Run:  python code/v2_figures.py [name ...]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
sys.path.insert(0, str(HERE))
from figstyle import C, bilingual, plt, tr  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

A1 = "code/bsc_validation2/a1_more_outbreaks/"
A2 = "code/bsc_validation2/a2_contact_mobility/"
A3 = "code/bsc_validation2/a3_tracer_physics/"
A4 = "code/bsc_validation2/a4_zone_level/"
# results, tracer and school tables, and re-checks (written out in full: the release build maps each
# of these directories to its place in the repository)
A1O = "data/bsc_validation2/a1_more_outbreaks/"
A1V = "code/checks/outbreaks2/"
A2O = "data/bsc_validation2/a2_contact_mobility/"
A2V = "code/checks/contacts/"
A2SUM = "code/checks/contacts/check_summary.json"
A3R = "data/bsc_validation2/a3_tracer_physics/"
A3D = "data/tracer/"
A4R = "data/bsc_validation2/a4_zone_level/"
A4D = "data/schools/"
A4V = "code/checks/zones/"
R1 = "code/bsc_validation/"


def _type(size: float = 8.5):
    """Raise the default type of a figure that is printed smaller than it is drawn (call at the top of a figure
    function; figstyle.setup resets the defaults before the next figure)."""
    plt.rcParams.update({"font.size": size, "axes.titlesize": size + 0.5, "axes.labelsize": size,
                         "legend.fontsize": size, "xtick.labelsize": size, "ytick.labelsize": size})


def load(rel):
    return json.loads((ROOT / rel).read_text())


def read_csv(rel):
    import csv
    with open(ROOT / rel, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


COL = dict(obs="#222222", M0=C["grey"], CRR=C["orange"], M2=C["blue"], P=C["green"], M1=C["red"])

EV_EN = {"F1": "Hong Kong–Beijing (SARS-CoV-1)", "F2": "Chicago–Honolulu (tuberculosis)",
         "F3": "Los Angeles–Auckland (influenza)", "F4": "Cancun–Birmingham (influenza)",
         "F5": "Flight to Naha", "F6": "Sydney–Perth", "F7": "Dubai–Auckland (EK448)", "F8": "Tel Aviv–Frankfurt",
         "W1": "Ward 8A, students (SARS-CoV-1)", "W2": "Ward 8A, inpatients (SARS-CoV-1)",
         "P1": "Meat-processing line"}
EV_ZH = {"F1": "香港–北京（SARS-CoV-1）", "F2": "芝加哥–檀香山（结核）", "F3": "洛杉矶–奥克兰（流感）",
         "F4": "坎昆–伯明翰（流感）", "F5": "飞往那霸的航班", "F6": "悉尼–珀斯", "F7": "迪拜–奥克兰（EK448）",
         "F8": "特拉维夫–法兰克福", "W1": "8A 病房医学生（SARS-CoV-1）", "W2": "8A 病房住院病人（SARS-CoV-1）",
         "P1": "肉类加工线"}
ORDER = ["F1", "F2", "F3", "F4", "W1", "W2", "F5", "F6", "F7", "F8", "P1"]
HOLD = ["F5", "F6", "F7", "F8", "P1"]


def ev_label(lang, e):
    return (EV_EN if lang == "en" else EV_ZH)[e]


# ------------------------------------------------------------------------------------------ risk ratios
def rr_forest(lang):
    des = load(A1O + "07_descriptive.json")
    pooled = [("all 11 events", tr(lang, "All eleven events", "全部十一起")),
              ("pathogen: SARS-CoV-2", tr(lang, "SARS-CoV-2 (5 events)", "SARS-CoV-2（5 起）")),
              ("pathogen: SARS-CoV-1", tr(lang, "SARS-CoV-1 (3 events)", "SARS-CoV-1（3 起）")),
              ("pathogen: influenza A(H1N1)pdm09", tr(lang, "Influenza A(H1N1)pdm09 (2 events)", "甲型 H1N1 流感（2 起）"))]
    _type(8)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))    # printed at 0.82 text width: type of about 8 pt
    n = len(ORDER)
    ys = list(range(n + len(pooled) + 1, len(pooled) + 1, -1))
    for y, e in zip(ys, ORDER):
        r = des["events"][e]
        hold = r["role"] == "hold-out"
        ax.plot([r["lo"], min(r["hi"], 400)], [y, y], color=COL["obs"], lw=1)
        ax.plot(r["RR"], y, "s", color=COL["obs"], mfc=COL["obs"] if hold else "white", ms=5.5, zorder=3)
        for m, mk, dy in (("CRR", "D", 0.24), ("M2", "o", -0.24)):
            v = r["pred_RR"][m]
            if v is not None and math.isfinite(v):
                ax.plot(min(v, 400), y + dy, mk, color=COL[m], ms=4, zorder=2)
    yp = list(range(len(pooled), 0, -1))
    for y, (key, lab) in zip(yp, pooled):
        p = des["pooled"][key]
        ax.plot([p["lo_r"], p["hi_r"]], [y, y], color="k", lw=2.2)
        ax.plot(p["RR_random"], y, "D", color="k", ms=6)
    ax.axvline(1, color="#999999", lw=0.8, ls=":")
    ax.axhline(len(pooled) + 1, color="#cccccc", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.15, 500)
    ax.set_yticks(ys + yp)
    ax.set_yticklabels([ev_label(lang, e) for e in ORDER] + [lab for _, lab in pooled], fontsize=8)
    ax.set_xlabel(tr(lang, "risk ratio, near / far (log scale)", "近/远风险比（对数坐标）"))
    ax.legend(handles=[
        Line2D([], [], marker="s", color=COL["obs"], ls="-", lw=1, label=tr(lang, "observed, held out", "观测值（留出）")),
        Line2D([], [], marker="s", color=COL["obs"], mfc="white", ls="-", lw=1,
               label=tr(lang, "observed, calibration", "观测值（标定）")),
        Line2D([], [], marker="D", color="k", ls="-", lw=2, label=tr(lang, "pooled, random effects", "合并（随机效应）")),
        Line2D([], [], marker="D", color=COL["CRR"], ls="", ms=4, label=tr(lang, "constant risk ratio", "恒定风险比")),
        Line2D([], [], marker="o", color=COL["M2"], ls="", ms=4, label="M2")],
        fontsize=7.5, loc="lower right", frameon=True, framealpha=0.9, edgecolor="none")
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------- hold-out strata
def holdout_strata(lang):
    H = load(A1O + "10_corrected.json")["primary_corrected"]["events"]
    binl = {4: ["0–1", "2", "3–5", "6+"], 3: ["0–1", "2", "3+"]}
    fig, axs = plt.subplots(1, 5, figsize=(8.2, 3.0))     # printed at text width: type of about 7 pt
    for ax, e in zip(axs, HOLD):
        r = H[e]
        nb = len(r["sizes"])
        x = np.arange(nb)
        w = 0.2
        ax.bar(x - 1.5 * w, r["obs"], w, color=COL["obs"], label=tr(lang, "observed", "观测"))
        for k, (m, lab) in enumerate((("M0", tr(lang, "well mixed (M0)", "均匀混合（M0）")),
                                      ("CRR", tr(lang, "constant risk ratio", "恒定风险比")),
                                      ("M2", "M2"))):
            ax.bar(x + (k - 0.5) * w, r[m]["expected"], w, color=COL[m], label=lab)
        lab = ["0–4", "4–8", "8–12", ">12"] if e == "P1" else binl[nb]
        ax.set_xticks(x)
        ax.set_xticklabels([f"{a}\n({s})" for a, s in zip(lab, r["sizes"])], fontsize=8.5)
        ax.tick_params(axis="y", labelsize=8.5)
        ax.set_xlabel(tr(lang, "distance from the index\ncase (m)", "与指示病例的距离（m）") if e == "P1"
                      else tr(lang, "rows from the\nnearest source", "距最近传染源的排数"), fontsize=8.5)
        ax.set_title(ev_label(lang, e), fontsize=9)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    axs[0].set_ylabel(tr(lang, "cases in the stratum", "该层病例数"))
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=8.5, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(w_pad=0.4, rect=(0, 0.08, 1, 1))
    return fig


# --------------------------------------------------------------------------------------- summary figure
DIG = 0.11                                       # vertical offset of the digit markers in panel (b)


def summary_scores(lang):
    N = load("data/v2_numbers.json")
    first = N["first.per_event"]["value"]
    ev3 = load(A3R + "04_evaluation.json")
    ex3 = load(A3R + "E2_ext_evaluation.json")
    cor = load(A1O + "10_corrected.json")["primary_corrected"]["events"]
    rows = [("h", tr(lang, "First dataset (kernel calibrated on trains / restaurant)", "第一组（核由列车队列、餐厅标定）")),
            ("T1", tr(lang, "Zhejiang bus", "浙江大巴")), ("T2", tr(lang, "Hunan coach", "湖南客车")),
            ("T5", tr(lang, "Flight VN54", "VN54 航班")), ("C1", tr(lang, "Marin classroom", "马林县教室")),
            ("R1", tr(lang, "Guangzhou restaurant", "广州餐厅")),
            ("h", tr(lang, "Second dataset (kernel calibrated on pre-2020 outbreaks)", "第二组（核由 2020 年前的疫情标定）")),
            ("F5", ev_label(lang, "F5")), ("F6", ev_label(lang, "F6")), ("F7", ev_label(lang, "F7")),
            ("F8", ev_label(lang, "F8")), ("P1", ev_label(lang, "P1"))]
    ny = len(rows)
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 5.0), sharey=True)    # printed at text width: type of about 7 pt
    for ax in axs:
        ax.axvline(0, color="#666666", lw=0.8)
    yt, yl = [], []
    for i, (e, lab) in enumerate(rows):
        y = ny - 1 - i
        yt.append(y)
        yl.append(lab)
        if e == "h":
            for ax in axs:
                ax.axhspan(y - 0.5, y + 0.5, color="#eeeeee", zorder=0)
            continue
        a, b = axs
        if e in first:                                                # outbreak-calibrated M2, first set
            f = first[e]
            a.plot(f["d0"], y + 0.17, "o", color=COL["M2"], ms=6)
            b.plot([f["dC2"], f["dC3"]], [y + 0.17, y + 0.17], color=COL["M2"], lw=1.2)
            b.plot(0.5 * (f["dC2"] + f["dC3"]), y + 0.17, "o", color=COL["M2"], ms=3.5)
            b.plot(f["dC2"], y + 0.17 + DIG, marker="$2$", color=COL["M2"], ms=5.5)    # digits above / below the
            b.plot(f["dC3"], y + 0.17 - DIG, marker="$3$", color=COL["M2"], ms=5.5)    # line, so that they cannot overprint
        if e in cor:                                                  # outbreak-calibrated M2, second set
            r = cor[e]
            a.plot(r["M2"]["logp"] - r["M0"]["logp"], y + 0.17, "o", color=COL["M2"], ms=6)
            b.plot(r["M2"]["logp"] - r["CRR"]["logp"], y + 0.17, "o", color=COL["M2"], ms=6)
        src = ev3["events"].get(e)
        if src is not None:                                           # tracer-fixed kernel, first set
            m = src["models"]
            d0 = m["P"]["logscore"] - m["M0"]["logscore"]
            d2 = m["P"]["logscore"] - m["CRR2"]["logscore"]
            d3 = m["P"]["logscore"] - m["CRR3"]["logscore"]
            a.plot(d0, y - 0.17, "s", color=COL["P"], ms=5.5)
            b.plot([d2, d3], [y - 0.17, y - 0.17], color=COL["P"], lw=1.2)
            b.plot(0.5 * (d2 + d3), y - 0.17, "s", color=COL["P"], ms=3.5)
            b.plot(d2, y - 0.17 + DIG, marker="$2$", color=COL["P"], ms=5.5)
            b.plot(d3, y - 0.17 - DIG, marker="$3$", color=COL["P"], ms=5.5)
        if e in ex3["events"]:                                        # tracer-fixed kernel, exploratory
            m = ex3["events"][e]["models"]
            d0 = m["P"]["logp"] - m["M0"]["logp"]
            d2 = m["P"]["logp"] - m["CRR2"]["logp"]
            d3 = m["P"]["logp"] - m["CRR3"]["logp"]
            a.plot(d0, y - 0.17, "s", color=COL["P"], mfc="white", ms=5.5)
            b.plot([d2, d3], [y - 0.17, y - 0.17], color=COL["P"], lw=1.2, ls=":")
            b.plot(0.5 * (d2 + d3), y - 0.17, "s", color=COL["P"], mfc="white", ms=3.5)
            b.plot(d2, y - 0.17 + DIG, marker="$2$", color=COL["P"], ms=5.5, alpha=0.6)
            b.plot(d3, y - 0.17 - DIG, marker="$3$", color=COL["P"], ms=5.5, alpha=0.6)
    axs[0].set_yticks(yt)
    axs[0].set_yticklabels(yl, fontsize=8.5)
    for t, (e, _) in zip(axs[0].get_yticklabels(), rows):
        if e == "h":
            t.set_fontweight("bold")
    axs[0].set_ylim(-0.6, ny - 0.4)
    axs[0].set_xlabel(tr(lang, "log score:\nmodel minus well mixed", "对数得分：\n模型减均匀混合"))
    axs[1].set_xlabel(tr(lang, "log score: model\nminus constant ratio", "对数得分：\n模型减恒定风险比"))
    axs[0].set_title(tr(lang, "(a) against well mixed", "(a) 对比均匀混合"), loc="left")
    axs[1].set_title(tr(lang, "(b) against a constant ratio", "(b) 对比恒定风险比"), loc="left")
    fig.legend(handles=[
        # filled column by column: row 1 holds M2 and the digits, row 2 the two tracer-fixed entries
        Line2D([], [], marker="o", color=COL["M2"], ls="", label=tr(lang, "M2, mixing coefficient fitted to outbreaks", "M2，混合系数由疫情拟合")),
        Line2D([], [], marker="s", color=COL["P"], ls="", label=tr(lang, "mixing coefficient fixed by tracer measurements", "混合系数由示踪测量确定")),
        Line2D([], [], marker="$2$", color="k", ls="-", lw=1.0, label=tr(lang, "against fixed risk ratios 2 (above) and 3 (below)", "对比固定风险比 2（上）与 3（下）")),
        Line2D([], [], marker="s", color=COL["P"], mfc="white", ls="", label=tr(lang, "tracer-fixed kernel, exploratory (outcome-aware)", "示踪确定的核，探索性（已知结果）"))],
        fontsize=8.5, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(w_pad=1.0, rect=(0, 0.11, 1, 1))
    return fig


# -------------------------------------------------------------------------------------- contact records
DS = ["InVS13", "LyonSchool", "LH10", "InVS15", "Thiers13", "SFHH"]
STAGE = {"InVS13": "s1", "LyonSchool": "s1", "LH10": "s1", "InVS15": "s3a", "Thiers13": "s3a", "SFHH": "s3a"}
EXTRA = {"RWclus": ("s2b", "s3b"), "RWstickZ": ("s2d", "s3c"), "RWstickO": ("s2d", "s3d")}


def _contact_models(lang):
    return [("WM", "s", C["grey"], tr(lang, "well mixed", "均匀混合")),
            ("BLOCK", "D", C["orange"], tr(lang, "constant within/between-group ratio", "恒定的组内/组间比")),
            ("RW0", "o", C["red"], tr(lang, "lattice walk, homogeneous", "格点随机游走（均匀）")),
            ("RWhet", "v", C["purple"], tr(lang, "lattice walk, two zones", "格点随机游走（两区）")),
            ("RWclus", "P", C["sky"], tr(lang, "home-anchored walk (post hoc)", "锚定于座位的游走（事后）")),
            ("RWstickO", "*", C["blue"], tr(lang, "home-anchored, slow own seat (post hoc)", "锚定且本座位慢（事后）"))]


def _contact_file(ds, model):
    if model in EXTRA:
        tag = EXTRA[model][0] if STAGE[ds] == "s1" else EXTRA[model][1]
    else:
        tag = STAGE[ds]
    return load(A2O + f"{tag}_{ds}.json")


def contacts_epidemic(lang):
    models = _contact_models(lang)
    quant = [("R_index", tr(lang, "secondary cases of the index case", "首例的二代病例数")),
             ("P_major", tr(lang, "probability of a major outbreak", "发生大规模暴发的概率")),
             ("attack", tr(lang, "mean attack rate", "平均罹患率"))]
    _type(8.5)
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.9))    # printed at text width: type of about 7 pt
    for ax, (q, lab) in zip(axs, quant):
        hi = 0
        for model, mk, col, _ in models:
            for ds in DS:
                d = _contact_file(ds, model)
                conf = STAGE[ds] != "s1"
                for sc, v in d["epi"].items():
                    x, y = v["REAL"][q], v[model][q]
                    hi = max(hi, x, y)
                    ax.plot(x, y, mk, color=col, mfc=col if conf else "white", ms=4.2, mew=0.9, alpha=0.9)
        # static network of training-day pair rates (re-check; exploratory)
        for f in ("v4_static_InVS13_LH10_InVS15.json", "v4_static_LyonSchool_Thiers13_SFHH.json"):
            for ds, d in load(A2V + f).items():
                for sc, v in d["epi"].items():
                    x, y = v["real"][q], v["STATIC"][q]
                    hi = max(hi, x, y)
                    ax.plot(x, y, "x", color="k", ms=4.5, mew=1.0)
        hi *= 1.06
        xs = np.linspace(0, hi, 50)
        if q == "R_index":
            ax.fill_between(xs, 0.9 * xs, xs / 0.9, color=C["green"], alpha=0.18, lw=0)
        else:
            ax.fill_between(xs, xs - 0.05, xs + 0.05, color=C["green"], alpha=0.18, lw=0)
        ax.plot([0, hi], [0, hi], color="k", lw=0.8)
        ax.set_xlim(0, hi)
        ax.set_ylim(0, hi)
        ax.set_xlabel(lab + tr(lang, ",\nrecorded contacts", "\n（实测接触）"))
        ax.set_ylabel(lab + tr(lang, ",\nmodel contacts", "\n（模型接触）"))
    handles = [Line2D([], [], marker=mk, color=col, ls="", label=lab) for _, mk, col, lab in models]
    handles.append(Line2D([], [], marker="x", color="k", ls="", label=tr(lang, "static network of training-day pair rates (exploratory)",
                                                                       "训练日配对接触率的静态网络（探索性）")))
    handles.append(plt.Rectangle((0, 0), 1, 1, color=C["green"], alpha=0.18, label=tr(lang, "pre-specified tolerance", "预先设定的容差")))
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.25, 1, 1))
    return fig


def contacts_structure(lang):
    models = _contact_models(lang)
    stats = [("S0_total", tr(lang, "total contact time", "接触总时间")),
             ("S1_mean_dur", tr(lang, "mean contact duration", "平均接触时长")),
             ("S1_timefrac_ge15", tr(lang, "share of contact time in contacts ≥ 5 min", "≥5 分钟接触所占时间份额")),
             ("S2_recurrence", tr(lang, "recurrence of pairs", "配对重复率")),
             ("S3_deg_mean", tr(lang, "partners per person-day", "每人每日接触对象数")),
             ("S3_deg_cv", tr(lang, "CV of partners", "接触对象数的变异系数")),
             ("S4_strength_cv", tr(lang, "CV of individual contact time", "个体接触时间的变异系数")),
             ("S5_persistence", tr(lang, "pairs repeated next day", "次日重复的配对比例")),
             ("S6_within_frac", tr(lang, "within-group share", "组内接触份额"))]
    names = {"InVS13": tr(lang, "office 2013", "办公楼 2013"), "LyonSchool": tr(lang, "primary school", "小学"),
             "LH10": tr(lang, "hospital ward", "医院病房"), "InVS15": tr(lang, "office 2015", "办公楼 2015"),
             "Thiers13": tr(lang, "high school", "高中"), "SFHH": tr(lang, "conference", "学术会议")}
    _type(8.5)
    fig, axs2 = plt.subplots(2, 3, figsize=(7.2, 6.8), sharey=True)    # printed at text width: type of about 7.5 pt
    axs = axs2.ravel()
    for ax, ds in zip(axs, DS):
        ax.axvspan(0.8, 1.25, color=C["green"], alpha=0.18, lw=0)
        ax.axvline(1, color="k", lw=0.8)
        for j, (model, mk, col, _) in enumerate(models):
            d = _contact_file(ds, model)
            for i, (k, _) in enumerate(stats):
                o, m = d["obs"].get(k), d["models"][model]["stats"].get(k)
                if o is None or m is None or not o or (isinstance(o, float) and math.isnan(o)):
                    continue
                y = len(stats) - 1 - i + (2.5 - j) * 0.12
                ax.plot(min(max(m / o, 0.02), 6), y, mk, color=col, ms=3.6, mew=0.8)
        ax.set_xscale("log")
        ax.set_xlim(0.018, 7)
        ax.set_title(names[ds])
        ax.set_xlabel(tr(lang, "model / observed", "模型 / 观测"))
    for row in axs2:
        row[0].set_yticks(range(len(stats)))
        row[0].set_yticklabels([s[1] for s in stats][::-1])
    fig.legend(handles=[Line2D([], [], marker=mk, color=col, ls="", label=lab) for _, mk, col, lab in models],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.11, 1, 1), w_pad=0.4)
    return fig


# ----------------------------------------------------------------------------------- tracer measurements
def tracer_kernels(lang):
    ker = load(A3R + "01_kernels.json")
    des = load(A3R + "05_descriptive.json")
    _type(8.5)
    fig, ax = plt.subplots(2, 2, figsize=(7.6, 7.0))    # printed at text width: type of about 7 pt
    lab_meas = tr(lang, "decay length of the\nmeasured coefficient", "实测系数对应的衰减长度")
    lab_fit = tr(lang, "decay length of the coefficient\nfitted to outbreaks (25 m$^2$ h$^{-1}$)",
                 "疫情拟合系数（25 m$^2$ h$^{-1}$）\n对应的衰减长度")
    # (a) coach: ethane tracer and tracer-validated CFD by seat row
    a = ax[0, 0]
    seats = read_csv(A3D + "ou2022_B1_seats.csv")
    cfd = [(float(s["row"]), float(s["s4_cfd"])) for s in seats if s["s4_cfd"] not in ("", "nan") and s["row"]]
    mx = max(v for _, v in cfd)
    a.scatter([r for r, _ in cfd], [v / mx for _, v in cfd], s=10, color=C["grey"], alpha=0.6, lw=0,
              label=tr(lang, "tracer-validated CFD, all seats", "经示踪验证的 CFD（全部座位）"))
    ms = [(float(s["row"]), float(s["s3_measured"])) for s in seats if s["s3_measured"] not in ("", "nan")]
    a.scatter([r for r, _ in ms], [v for _, v in ms], s=30, color=C["blue"], edgecolor="white", lw=0.8, zorder=4,
              label=tr(lang, "measured (ethane)", "实测（乙烷）"))
    pitch, kap = 0.877, 6.5                    # row pitch of the coach lattice (m); air change rate of the test (1/h)
    rr = np.linspace(1, 13, 100)
    for D, col, lab in ((ker["coach"]["D_hat"], C["blue"], lab_meas), (25.1, C["red"], lab_fit)):
        a.plot(rr, np.exp(-np.abs(rr - 12) * pitch / math.sqrt(D / kap)), color=col, lw=1.4, label=lab)
    a.set_yscale("log")
    a.set_ylim(0.004, 60)
    a.set_xlabel(tr(lang, "seat row (index case in row 12)", "座位排号（指示病例在第 12 排）"))
    a.set_ylabel(tr(lang, "tracer concentration (normalised)", "示踪物浓度（归一化）"))
    a.set_title(tr(lang, "(a) Coach", "(a) 长途客车"), loc="left")
    a.legend(fontsize=7.8, loc="upper left")
    # (b) cabin: time-integrated counts by row, two releases
    a = ax[0, 1]
    tidy = read_csv(A3D + "kinahan_tidy.csv")
    for rel, col, mk, dx, lab in (("5L", C["blue"], "o", -0.12, tr(lang, "release in seat 5L", "5L 座释放")),
                                  ("5A", C["purple"], "s", 0.12, tr(lang, "release in seat 5A", "5A 座释放"))):
        acc = {}
        for r in tidy:
            if (r["airframe"] == "777" and r["section"] == "FWD" and r["release"] == rel and r["kind"] == "B"
                    and r["mask"] == "No" and r["row"] not in ("", "nan") and r["count"] not in ("", "nan")):
                acc.setdefault(r["sensor"], []).append(float(r["count"]))
        xs, ys = [], []
        for sensor, v in acc.items():
            m = sum(v) / len(v)
            if m > 0 and sensor != rel:
                xs.append(float("".join(ch for ch in sensor if ch.isdigit())) - 5 + dx)
                ys.append(m)
        a.scatter(xs, ys, s=16, color=col, marker=mk, edgecolor="white", lw=0.6, label=lab, alpha=0.9)
    a.set_yscale("log")
    a.set_xlabel(tr(lang, "rows from the release seat (negative: forward)", "与释放座位相隔排数（负值为向前）"))
    a.set_ylabel(tr(lang, "time-integrated particle count", "时间积分粒子计数"))
    a.set_title(tr(lang, "(b) Wide-body forward cabin", "(b) 宽体客机前舱"), loc="left")
    a.set_ylim(50, 3e5)
    a.legend(fontsize=7.8, loc="upper center", ncol=2)
    # (c) rail carriage: salt aerosol along the saloon (middle release)
    a = ax[1, 0]
    pts = [p for p in read_csv(A3D + "woodward2022_fig9_points.csv")
           if p["panel"] == "B_middle" and p["near_source"] == "False"]
    a.scatter([float(p["x_over_L"]) for p in pts], [float(p["c_norm"]) for p in pts], s=18, color=C["blue"],
              edgecolor="white", lw=0.6, zorder=3, label=tr(lang, "measured (salt aerosol)", "实测（盐气溶胶）"))
    Ls, x0, kap = 19.85, -0.135, 13.0          # saloon length (m), release position, air change rate (1/h)
    xx = np.linspace(-0.5, 0.5, 200)
    top = max(float(p["c_norm"]) for p in pts)
    for D, col, lab in ((ker["train"]["B_middle"]["D_hat"], C["blue"], lab_meas), (25.1, C["red"], lab_fit)):
        a.plot(xx, top * np.exp(-np.abs(xx - x0) * Ls / math.sqrt(D / kap)), color=col, lw=1.4, label=lab)
    a.axvline(x0, color="#999999", lw=0.8)
    a.set_yscale("log")
    a.set_xlabel(tr(lang, "position along the saloon, $x/L$", "沿客室位置 $x/L$"))
    a.set_ylabel(tr(lang, "tracer concentration (normalised)", "示踪物浓度（归一化）"))
    a.set_title(tr(lang, "(c) Rail carriage", "(c) 铁路车厢"), loc="left")
    a.set_ylim(0.004, 80)
    a.legend(fontsize=7.8, loc="upper left")
    # (d) coefficients: measured against fitted to outbreaks
    a = ax[1, 1]
    rows = [des["D_table"][i] for i in (0, 1, 2, 4, 5)]
    encl = [tr(lang, "coach", "长途客车"), tr(lang, "aircraft cabin", "客机客舱"), tr(lang, "rail carriage", "铁路车厢"),
            tr(lang, "classroom", "教室"), tr(lang, "open-plan office", "开放式办公室")]
    cor = load(A1O + "10_corrected.json")["primary_corrected"]["post"]["air|M2"]
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i
        a.plot([r["lo"], r["hi"]], [y, y], color=C["blue"], lw=2, alpha=0.45)
        a.plot(r["D_phys"], y, "o", color=C["blue"], ms=6, zorder=4)
        a.plot(r["D_calibrated"], y, "D", color=C["red"], ms=5, zorder=4)
        if r["D_outbreak_alone"]:
            a.plot(r["D_outbreak_alone"], y, "^", color=C["orange"], ms=5, zorder=4)
    a.plot([cor["q05"], cor["q95"]], [len(rows) - 2 - 0.25] * 2, color=C["red"], lw=1.2)
    a.plot(cor["map"]["D"], len(rows) - 2 - 0.25, "D", color=C["red"], mfc="white", ms=5, zorder=4)
    a.set_xscale("log")
    a.set_yticks(range(len(rows)))
    a.set_yticklabels(encl[::-1])
    a.set_ylim(-0.6, len(rows) - 0.4)
    a.set_xlim(0.03, 3e5)
    a.set_xlabel(tr(lang, "mixing coefficient (m$^2$ h$^{-1}$, log scale)", "混合系数（m$^2$ h$^{-1}$，对数坐标）"))
    a.set_title(tr(lang, "(d) Measured and fitted coefficients", "(d) 实测与拟合的系数"), loc="left")
    a.legend(handles=[
        Line2D([], [], marker="o", color=C["blue"], ls="-", lw=2, label=tr(lang, "measured by tracer", "示踪实测")),
        Line2D([], [], marker="D", color=C["red"], ls="", label=tr(lang, "fitted to outbreaks, first test", "由疫情拟合（第一次检验）")),
        Line2D([], [], marker="D", color=C["red"], mfc="white", ls="-", lw=1.2,
               label=tr(lang, "fitted to four pre-2020 flights", "由 2020 年前四次航班拟合")),
        Line2D([], [], marker="^", color=C["orange"], ls="", label=tr(lang, "preferred by the outbreak alone", "单个疫情自身偏好的值"))],
        fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1)
    fig.tight_layout()
    return fig


def tracer_predictions(lang):
    ev = load(A3R + "04_evaluation.json")
    evs = [("T2", tr(lang, "Hunan coach\n(rear / front)", "湖南客车\n（后部 / 前部）")),
           ("T1", tr(lang, "Zhejiang bus\n(rows 5–11 / rest)", "浙江大巴\n（5–11 排 / 其余）")),
           ("C1", tr(lang, "Marin classroom\n(front / back)", "马林县教室\n（前排 / 后排）")),
           ("T5", tr(lang, "Flight VN54\n(≤2 seats / >2)", "VN54 航班\n（≤2 座 / >2 座）")),
           ("R1", tr(lang, "Restaurant, 2 cases\n(adjacent / remote tables)", "餐厅，2 例\n（邻桌 / 远桌）"))]
    keys = [("M0", COL["M0"], "o", tr(lang, "well mixed", "均匀混合")),
            ("P", COL["P"], "s", tr(lang, "tracer-fixed kernel", "示踪确定的核")),
            ("CRR2", COL["CRR"], "<", tr(lang, "constant ratio 2", "恒定风险比 2")),
            ("CRR3", COL["CRR"], ">", tr(lang, "constant ratio 3", "恒定风险比 3")),
            ("M2cal_mode", COL["M2"], "o", tr(lang, "M2 fitted to outbreaks", "由疫情拟合的 M2"))]
    fig, axs = plt.subplots(1, 5, figsize=(8.4, 3.1))     # printed at text width: type of about 7 pt
    for j, (e, title) in enumerate(evs):
        a = axs[j]
        d = ev["events"][e]
        for i, (k, col, mk, _) in enumerate(keys):
            if k not in d["models"]:
                continue
            r = d["models"][k]
            y = len(keys) - 1 - i
            a.plot(r["int95"], [y, y], color=col, lw=2.2, solid_capstyle="round", alpha=0.55)
            a.plot(r["mean_near"], y, mk, color=col, ms=5.5, zorder=4)
        a.axvline(d["k_obs"], color="k", lw=1.2)
        a.set_ylim(-0.6, len(keys) - 0.4)
        a.set_yticks(range(len(keys)))
        a.set_yticklabels([k[3] for k in keys][::-1] if j == 0 else [""] * len(keys), fontsize=8.5)
        a.tick_params(axis="x", labelsize=8.5)
        a.set_title(title, fontsize=8.5)
        a.set_xlim(-0.5, min(d["n_near"], d["K"]) + 0.5)
        a.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    fig.supxlabel(tr(lang, "number of cases in the near stratum (vertical line: observed)",
                     "近区病例数（竖线为观测值）"), fontsize=9)
    fig.tight_layout(w_pad=0.4)
    return fig


# ------------------------------------------------------------------------------------------ zone level
def zone(lang):
    R = load(A4R + "matsumoto_primary.json")
    P = R["primary"]
    S = load(A4R + "seoul_floors.json")
    lyon = load(A4D + "lyon_mixing.json")
    shares = lyon.get("shares_time") or lyon.get("w") or lyon.get("shares")
    if isinstance(shares, dict):
        shares = [shares[k] for k in list(shares)[:3]]
    _type(8.5)
    fig, ax = plt.subplots(2, 2, figsize=(7.6, 7.0))    # printed at text width: type of about 7 pt
    a = ax[0, 0]
    x = np.arange(3)
    w = np.array(P["w_hat"])
    ci = np.array(R["w_hat_boot_ci"])
    a.bar(x - 0.19, shares, 0.36, color=C["blue"], label=tr(lang, "contact durations measured\nin another school", "在另一所学校实测的接触时长"))
    a.bar(x + 0.19, w, 0.36, color=C["orange"], label=tr(lang, "fitted to the epidemic\n(95 % bootstrap over schools)", "由疫情拟合（按学校自助法 95% 区间）"))
    a.errorbar(x + 0.19, w, yerr=[w - ci[:, 0], ci[:, 1] - w], fmt="none", ecolor="k", capsize=3, lw=1)
    for i in range(3):
        a.text(x[i] - 0.19, shares[i] + 0.015, f"{shares[i]:.2f}", ha="center", fontsize=8)
        a.text(x[i] + 0.19, ci[i, 1] + 0.015, f"{w[i]:.2f}", ha="center", fontsize=8)
    a.set_xticks(x)
    a.set_xticklabels([tr(lang, "own\nclass", "本班"), tr(lang, "other classes,\nsame grade", "同年级\n其他班"),
                       tr(lang, "rest of\nschool", "学校\n其余部分")])
    a.set_ylabel(tr(lang, "share of within-school transmission", "校内传播的份额"))
    a.set_ylim(0, 1.4)
    a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.legend(fontsize=7.8, loc="upper right")
    a.set_title("(a)", loc="left")
    b = ax[0, 1]
    dH, dX, dC = (np.array(P[f"D_Z_minus_{k}_per_school"]) for k in "HXC")
    o = np.argsort(dH)
    xs = np.arange(len(o))
    b.axhline(0, color="#999999", lw=0.8)
    b.plot(xs, dH[o], "o-", color=C["blue"], ms=3.5, lw=0.9,
           label=tr(lang, f"against a well-mixed school (total {dH.sum():+.0f})", f"对比全校均匀混合（合计 {dH.sum():+.0f}）"))
    b.plot(xs, dX[o], "s-", color=C["orange"], ms=3.5, lw=0.9,
           label=tr(lang, f"against independent classes (total {dX.sum():+.0f})", f"对比各班相互独立（合计 {dX.sum():+.0f}）"))
    b.plot(xs, dC[o], "^-", color=C["grey"], ms=3.5, lw=0.9,
           label=tr(lang, f"against one fitted class/school ratio (total {dC.sum():+.0f})", f"对比拟合的单一班/校比值（合计 {dC.sum():+.0f}）"))
    b.set_xlabel(tr(lang, "school (sorted)", "学校（已排序）"))
    b.set_ylabel(tr(lang, "held-out log-likelihood difference", "留出对数似然之差"))
    b.set_xticks([])
    b.set_ylim(top=max(dH.max(), dX.max(), dC.max()) * 1.45)
    b.legend(fontsize=7.8, loc="upper left")
    b.set_title("(b)", loc="left")
    c = ax[1, 0]
    cols = [C["blue"], C["grey"], C["orange"]]
    names = [tr(lang, "zone model", "分区模型"), tr(lang, "school well mixed", "全校均匀混合"),
             tr(lang, "classes independent", "各班相互独立")]
    for i, m in enumerate("ZHX"):
        lo, md, hi = P[f"disp_pred_{m}"]
        c.plot([lo, hi], [i, i], color=cols[i], lw=6, solid_capstyle="butt")
        c.plot(md, i, "|", color="k", ms=11)
        c.text(hi + 0.06, i, f"{lo:.2f}–{hi:.2f}", va="center", fontsize=8)
    c.axvline(P["disp_obs"], color="k", ls="--", lw=1)
    c.set_yticks(range(3))
    c.set_yticklabels(names)
    c.set_ylim(-0.6, 2.6)
    c.invert_yaxis()
    c.set_xlim(0.5, 4.8)
    c.set_xlabel(tr(lang, "between-class heterogeneity of attack rates\n(dashed: observed)", "班间罹患率异质性（虚线为观测值）"))
    c.set_title("(c)", loc="left")
    d = ax[1, 1]
    lo, hi = S["pred_count_interval95"]
    md, wm = S["pred_count_median"], S["well_mixed_expected"]
    d.plot([lo, hi], [0, 0], color=C["blue"], lw=6, solid_capstyle="butt")
    d.plot(md, 0, "|", color="k", ms=11)
    d.text(hi * 1.25, 0, f"{lo}–{hi}", va="center", fontsize=8)
    d.plot(wm, 1, "o", color=C["grey"], ms=6)
    d.text(wm * 1.3, 1, f"{wm:.0f}", va="center", fontsize=8)
    d.plot(0, 2, "o", color=C["orange"], ms=6)
    d.text(0.25, 2, "0", va="center", fontsize=8)
    d.axvline(S["obs_other"], color="k", ls="--", lw=1)
    d.set_xscale("symlog", linthresh=1)
    d.set_xlim(-0.3, 300)
    d.set_xticks([0, 1, 3, 10, 30, 100])
    d.set_xticklabels(["0", "1", "3", "10", "30", "100"])
    d.set_yticks(range(3))
    d.set_yticklabels([tr(lang, "floors with a shared\nlift and lobby zone", "楼层 + 共用电梯/大厅"),
                       tr(lang, "building well mixed", "整栋楼均匀混合"), tr(lang, "floors independent", "各楼层相互独立")])
    d.set_ylim(-0.6, 2.6)
    d.invert_yaxis()
    d.set_xlabel(tr(lang, "cases among the 927 occupants of other\nfloors (dashed: observed)", "其他楼层 927 人中的病例数（虚线为观测值）"))
    d.set_title("(d)", loc="left")
    fig.tight_layout(w_pad=2.0)
    return fig


# ---------------------------------------------------------------------------------- likelihood profiles
def D_profiles(lang):
    cor = load(A1O + "10_corrected.json")["preferred_D_corrected"]
    x = 10 ** np.round(np.linspace(-2.0, 4.0, 61), 6)          # grid of the analysis (log10 D from -2 to 4, step 0.1)
    cm = plt.get_cmap("tab10")
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 6.2), sharey=True)    # printed at text width: type of about 8 pt
    for ax, cls, title in zip(axs, ("air", "room"), (tr(lang, "(a) aircraft cabins", "(a) 客机客舱"),
                                                      tr(lang, "(b) rooms", "(b) 房间"))):
        tot = 0
        for k, (e, cur) in enumerate(cor[cls]["curves"].items()):
            cur = np.array(cur)
            ax.plot(x, cur - cur[-1], color=cm(k), lw=1.2, ls="-" if e in HOLD else "--", label=ev_label(lang, e))
            tot = tot + cur - cur[-1]
        if cls == "air":
            ax.plot(x, tot, color="k", lw=2.0, label=tr(lang, "sum (common coefficient)", "合计（共同系数）"))
        ax.axhline(0, color="#999999", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_ylim(-12, 24)
        ax.set_xlabel(tr(lang, "mixing coefficient $D_{\\mathrm{air}}$ (m$^2$ h$^{-1}$)", "混合系数 $D_{\\mathrm{air}}$（m$^2$ h$^{-1}$）"))
        ax.set_title(title, loc="left")
        ax.legend(fontsize=8, ncol=1, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
        ax.set_ylabel(tr(lang, "log-likelihood relative to\n$D_{\\mathrm{air}}=10^{4}$ m$^2$ h$^{-1}$",
                         "相对于 $D_{\\mathrm{air}}=10^{4}$ m$^2$ h$^{-1}$\n的对数似然"))
    fig.tight_layout()
    return fig


FIGS = {"v2_rr_forest": rr_forest, "v2_holdout_strata": holdout_strata, "v2_summary_scores": summary_scores,
        "v2_contacts_epidemic": contacts_epidemic, "v2_contacts_structure": contacts_structure,
        "v2_tracer_kernels": tracer_kernels, "v2_tracer_predictions": tracer_predictions, "v2_zone": zone,
        "v2_D_profiles": D_profiles}

if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGS):
        print(bilingual(FIGS[name], name))
