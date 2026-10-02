"""appB_dataset_fig_attack.py -- Fig. appB_dataset:fig:attack (bilingual).

Attack rates of the 23 records of the outbreak dataset that have an exposed denominator, with exact
(Clopper-Pearson) 95 % intervals, grouped as in Table appB_dataset:tab:records.

Regenerated from fig. 1 of code/bsc_outbreaks/scripts/make_figures.py with the shared style of
code/figstyle.py: no in-figure title, white background, English and Chinese labels, records in the order
of the table, and the grade column carries the re-check's caveat (A* = graded A in the dataset but not
meeting every stated grade-A criterion).

Input : ../data/appB_dataset_numbers.json   (written by appB_dataset_table.py from
                                             data/bsc_outbreaks/outbreaks.csv; run that script first)
Output: ../figures/appB_dataset_attack_{en,zh}.pdf (+ .png)
        ../data/appB_dataset_fig_attack.json      (the rows drawn)

    python appB_dataset_table.py && python appB_dataset_fig_attack.py
"""
from __future__ import annotations

import json
from pathlib import Path

from figstyle import C, bilingual, plt, tr

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
numbers = json.loads((DATA / "appB_dataset_numbers.json").read_text())
recs = {r["id"]: r for r in numbers["records"] if r.get("n")}
assert len(recs) == 23

NAME = {  # id: (English, Chinese)
    "O1": ("Seoul call centre, 11th floor", "首尔呼叫中心（11 层）"),
    "O2": ("Zurich open-plan office team", "苏黎世开放式办公室团队"),
    "O3": ("England public-facing office", "英格兰对外服务办公室"),
    "O4": ("German meat plant, fixed stations", "德国肉类加工厂（固定工位）"),
    "O5": ("Office building (tuberculosis)", "办公楼（结核病）"),
    "S1": ("Liaocheng supermarket, staff", "聊城超市员工"),
    "S2": ("Tianjin department store, staff", "天津百货商场员工"),
    "S4": ("Massachusetts grocery, staff", "马萨诸塞州食品店员工"),
    "C1": ("Marin County classroom", "马林县小学教室"),
    "C2": ("Jerusalem high school, students", "耶路撒冷中学学生"),
    "C3": ("Utah primary classrooms (masked)", "犹他州小学教室（戴口罩）"),
    "C5": ("Cheonan fitness dance classes", "天安健身舞蹈班"),
    "T1": ("Zhejiang bus", "浙江大巴"),
    "T2": ("Hunan coach", "湖南长途客车"),
    "T3": ("Hunan minibus", "湖南小巴"),
    "T4": ("China high-speed trains", "中国高铁"),
    "T5": ("Flight VN54, business cabin", "VN54 航班公务舱"),
    "R1": ("Guangzhou restaurant (other tables)", "广州餐厅（其他餐桌）"),
    "R2": ("Jeonju restaurant", "全州餐厅"),
    "H1": ("Skagit Valley choir", "斯卡吉特合唱团排练"),
    "H2": ("Sydney church services", "悉尼教堂礼拜"),
    "X1": ("Diamond Princess (cruise ship)", "钻石公主号邮轮"),
    "X2": ("Singapore dormitories", "新加坡外籍劳工宿舍"),
}
GROUPS = [("O", "Offices and workplaces", "办公室与工作场所"),
          ("S", "Supermarket and retail", "超市与零售"),
          ("C", "Classrooms and class-like rooms", "教室及类似场所"),
          ("T", "Vehicle cabins (metro stand-ins)", "交通工具车厢（代替地铁车厢）"),
          ("R", "Restaurants", "餐厅"),
          ("H", "Choir and church", "合唱与教堂"),
          ("X", "Outside the scope of the model", "超出模型范围")]


def dur_label(r, lang):
    T = r["exposure_duration_h"]
    if T is None:
        return tr(lang, "multi-day", "多日")
    if r["id"] == "C5":
        return tr(lang, "50 min per class", "每节课 50 分钟")
    if r["id"] == "T4":
        return tr(lang, "mean 2.1 h", "平均 2.1 小时")
    if T < 2:
        return tr(lang, f"{T * 60:.0f} min", f"{T * 60:.0f} 分钟")
    return tr(lang, f"{T:.3g} h", f"{T:.3g} 小时")


def grade_label(r):
    return r["grade"] + ("" if r["grade_criteria_met"] else "*")


rows = []
for prefix, en, zh in GROUPS:
    rows.append(("header", (en, zh)))
    for sid in sorted(k for k in recs if k[0] == prefix):
        rows.append(("rec", recs[sid]))

figdata = [dict(id=r["id"], k=r["k"], n=r["n"], attack_rate_pct=100 * r["attack_rate"],
                ci95_pct=[100 * v for v in r["ci95"]], alt=r["alt"], grade=grade_label(r),
                exposure=dur_label(r, "en")) for kind, r in rows if kind == "rec"]
(DATA / "appB_dataset_fig_attack.json").write_text(json.dumps(figdata, indent=1))


def draw(lang: str):
    n = len(rows)
    # fixed layout (inches): labels 2.35 | axes 2.5 | cases, exposure, grade columns 2.15
    W, H = 7.0, 0.19 * n + 0.6
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([2.35 / W, 0.45 / H, 2.5 / W, 1 - 0.6 / H])
    yax = ax.get_yaxis_transform()
    ink, muted = "#222222", "#666666"
    y = n
    for kind, item in rows:
        y -= 1
        if kind == "header":
            ax.text(-0.93, y, tr(lang, *item), fontsize=8, fontweight="bold", color=ink, va="center", ha="left",
                    transform=yax)
            continue
        r = item
        ar, (lo, hi) = 100 * r["attack_rate"], [100 * v for v in r["ci95"]]
        ax.plot([lo, hi], [y, y], color=C["blue"], lw=1.3, solid_capstyle="round", zorder=2)
        if r["alt"] is not None:
            ax.plot(100 * r["alt"][0] / r["alt"][1], y, "o", ms=5, mfc="white", mec=C["blue"], mew=1.1, zorder=3)
        ax.plot(ar, y, "o", ms=5, color=C["blue"], mec="white", mew=0.8, zorder=4, clip_on=False)
        ax.text(-0.03, y, f"{r['id']}  " + tr(lang, *NAME[r["id"]]), fontsize=7.5, color=ink, va="center",
                ha="right", transform=yax)
        dag = "†" if r["includes_index"] else ""
        ax.text(1.04, y, f"{r['k']:,}/{r['n']:,}{dag}", fontsize=7, color=ink, va="center", ha="left", transform=yax)
        ax.text(1.40, y, dur_label(r, lang), fontsize=7, color=ink, va="center", ha="left", transform=yax)
        ax.text(1.77, y, grade_label(r), fontsize=7, color=ink, va="center", ha="left", transform=yax)
    for x, en, zh in ((1.04, "cases/exposed", "病例/暴露人数"), (1.40, "exposure", "暴露时长"), (1.77, "grade", "等级")):
        ax.text(x, n - 0.25, tr(lang, en, zh), fontsize=7, color=muted, ha="left", va="bottom", transform=yax)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(tr(lang, "attack rate among exposed persons (%), exact 95% interval",
                     "暴露人群罹患率（%），精确 95% 置信区间"))
    ax.xaxis.grid(True, color="#dddddd", lw=0.6)
    ax.set_axisbelow(True)
    return fig


if __name__ == "__main__":
    print("\n".join(bilingual(draw, "appB_dataset_attack")))
