"""s8_outbreaks_fig_callcentre.py -- Fig. s8_outbreaks:fig:callcentre (bilingual).

Seoul call centre, 11th floor (record O1 of the outbreak dataset):
  (a) seat map with case and non-case desks.  THIRD-PARTY DATA: the desk positions and case flags were
      digitised from Fig. 2 of Park et al. (2020), Emerg. Infect. Dis. 26(8), doi:10.3201/eid2608.201274,
      by code/bsc_outbreaks/scripts/digitize_callcentre.py.  The figure has no scale bar: coordinates are
      in desk pitches.  84 of the 94 reported case seats appear in the published plan.
  (b) join-count test inside the north wing: number of adjacent desk pairs (centres within 1.25 desk pitches)
      that are both case desks, against its distribution under random relabelling of the 79 cases among the
      137 desks (20,000 permutations stored by code/bsc_outbreaks/scripts/analyze_dataset.py).

  (c) probability, under model M2, of a join-count z as small as observed, as a function of the mixing
      coefficient (quadrature over a grid of 25 values in the re-check), with the
      90 % posterior interval of the restaurant-calibrated coefficient, the values integrated over that
      posterior, and the values for the well-mixed model and for the tracer-based coefficient of Section 7.8.
      This panel replaces the first analysis with 40 posterior draws, which had read as borderline.

Regenerated from panels a and c of fig. 4 of code/bsc_outbreaks/scripts/make_figures.py with the shared
style of code/figstyle.py (no in-figure titles, white background, English and Chinese labels).

Inputs : data/bsc_outbreaks/tables/park2020_callcentre_seats_digitized.csv
         data/bsc_outbreaks/tables/callcentre_joincount_null.npy
         data/validation_checks/callcentre_quadrature.json        (panel c)
         data/s8_outbreaks_numbers.json (calibration.restaurant.M2) (panel c)
         data/bsc_validation2/a3_tracer_physics/04_evaluation.json, 05_descriptive.json (panel c)
Outputs: ../figures/s8_outbreaks_callcentre_{en,zh}.pdf (+ .png)
         ../data/s8_outbreaks_callcentre_figdata.json   (the numbers printed in the figure)

    python s8_outbreaks_fig_callcentre.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from figstyle import C, bilingual, plt, tr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
OB = ROOT / "code/bsc_outbreaks/data"

seats = pd.read_csv(OB / "park2020_callcentre_seats_digitized.csv")
null = np.load(OB / "callcentre_joincount_null.npy")

nw = seats[seats.region == "north_wing"]
sw = seats[seats.region == "south_wing"]
kn, nn, ks, ns = int(nw.case.sum()), len(nw), int(sw.case.sum()), len(sw)

# join count in the north wing (same definition as analyze_dataset.py)
xy = nw[["x_pitch_units", "y_pitch_units"]].to_numpy()
d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
adj = (d > 0) & (d <= 1.25)
c = nw.case.to_numpy().astype(bool)
iu = np.triu_indices(len(nw), 1)
obs = int((c[iu[0]] & c[iu[1]] & adj[iu]).sum())
z = (obs - null.mean()) / null.std(ddof=1)
p = (np.sum(null >= obs) + 1) / (null.size + 1)

quad = json.loads((ROOT / "data/validation_checks/callcentre_quadrature.json").read_text())
post = json.loads((ROOT / "data/s8_outbreaks_numbers.json").read_text())["calibration"]["restaurant"]["M2"]
tr_ev = json.loads((ROOT / "data/bsc_validation2/a3_tracer_physics/04_evaluation.json").read_text())["events"]["O1"]["models"]
tr_D = json.loads((ROOT / "data/bsc_validation2/a3_tracer_physics/05_descriptive.json").read_text())["D_table"][5]

figdata = dict(north=[kn, nn], south=[ks, ns], case_seats_shown=int(seats.case.sum()), adjacent_pairs=int(adj[iu].sum()),
               observed_case_case_pairs=obs, null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)), z=float(z),
               p_one_sided=float(p), n_perm=int(null.size),
               M2_P_z_le_obs_equal_weight=quad["equal_weight"], M2_P_z_le_obs_pooled_weight=quad["pooled_weight"],
               restaurant_posterior_q05_q95=[post["q05"], post["q95"]],
               wellmixed_P_z_le_obs=tr_ev["M0"]["P_le_obs_pooled"], tracer_P_z_le_obs=tr_ev["P"]["P_le_obs_pooled"],
               tracer_D_office=[tr_D["lo"], tr_D["D_phys"], tr_D["hi"]],
               source="seat map digitised from Fig. 2 of Park et al. (2020), doi:10.3201/eid2608.201274")
(HERE.parent / "data" / "s8_outbreaks_callcentre_figdata.json").write_text(json.dumps(figdata, indent=1))
print(json.dumps(figdata, indent=1))


def draw(lang: str):
    fig = plt.figure(figsize=(9.4, 3.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.0, 1.05], wspace=0.36)

    # ---- (a) seat map -----------------------------------------------------------------------------
    a = fig.add_subplot(gs[0, 0])
    non, cas = seats[seats.case == 0], seats[seats.case == 1]
    a.plot(non.x_pitch_units, -non.y_pitch_units, "s", ms=3.3, mfc="white", mec=C["grey"], mew=0.6,
           label=tr(lang, "non-case desk", "非病例工位"))
    a.plot(cas.x_pitch_units, -cas.y_pitch_units, "s", ms=3.3, color=C["blue"], mec="white", mew=0.3,
           label=tr(lang, "case desk", "病例工位"))
    a.add_patch(plt.Rectangle((5.6, -19.2), 15.6, 9.6, fc="#eeeeee", ec=C["grey"], lw=0.5))
    a.text(13.4, -14.4, tr(lang, "service core\n(lifts, stairs, toilets)", "核心筒\n（电梯、楼梯、卫生间）"),
           ha="center", va="center", fontsize=7, color="#444444")
    a.text(27.3, -14.2, tr(lang, "east offices (not digitised)", "东侧办公室（未数字化）"), ha="center", va="center",
           fontsize=6.5, color=C["grey"], rotation=90)
    a.text(0.3, -9.6, tr(lang, f"north wing: {kn}/{nn} desks ({100 * kn / nn:.1f}%)",
                          f"北翼：{kn}/{nn} 个工位（{100 * kn / nn:.1f}%）"), fontsize=7.5, ha="left", va="center")
    a.text(7.0, -28.3, tr(lang, f"south wing: {ks}/{ns} desks ({100 * ks / ns:.1f}%)",
                          f"南翼：{ks}/{ns} 个工位（{100 * ks / ns:.1f}%）"), fontsize=7.5, ha="left", va="center")
    a.set_aspect("equal")
    a.set_xlim(-0.5, 28.6)
    a.set_ylim(-29.6, 2.6)
    a.set_xlabel(tr(lang, "x (desk pitches)", "x（工位间距）"))
    a.set_ylabel(tr(lang, "y (desk pitches)", "y（工位间距）"))
    a.set_yticks([-25, -20, -15, -10, -5, 0])
    a.set_yticklabels(["25", "20", "15", "10", "5", "0"])
    a.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=7, handletextpad=0.2, columnspacing=1.2)
    a.text(-0.16, 1.02, "(a)", transform=a.transAxes, fontsize=10, va="bottom", ha="left")

    # ---- (b) join-count permutation distribution -----------------------------------------------------
    b = fig.add_subplot(gs[0, 1])
    bins = np.arange(null.min() - 0.5, max(null.max(), obs) + 1.5)
    b.hist(null, bins=bins, color="#b0b0b0", lw=0)
    b.axvline(obs, color=C["blue"], lw=1.6)
    top = b.get_ylim()[1] * 1.25
    b.set_xlim(obs - 34, max(null.max(), obs) + 2)
    b.set_ylim(0, top)
    b.text(obs - 0.8, top * 0.98,
           tr(lang, f"observed: {obs} pairs\n", f"观测值：{obs} 对\n") + f"$z = {z:.2f}$\n$p = {p:.2f}$",
           fontsize=7.5, va="top", ha="right")
    b.text(0.97, 0.80, tr(lang, f"random relabelling:\n{null.mean():.1f} $\\pm$ {null.std(ddof=1):.1f}",
                          f"随机重标记：\n{null.mean():.1f} $\\pm$ {null.std(ddof=1):.1f}"),
           transform=b.transAxes, fontsize=7.5, va="top", ha="right", color="#444444")
    b.set_xlabel(tr(lang, "adjacent desk pairs that are both cases\n(north wing)", "两个工位均为病例的相邻工位对数\n（北翼）"))
    b.set_ylabel(tr(lang, f"permutations (of {null.size:,})", f"置换次数（共 {null.size:,} 次）"))
    b.text(-0.2, 1.02, "(b)", transform=b.transAxes, fontsize=10, va="bottom", ha="left")

    # ---- (c) probability of the observed join count under M2, by mixing coefficient ---------------------
    cax = fig.add_subplot(gs[0, 2])
    D = 10 ** np.array(quad["grid_log10D"])
    cax.axvspan(post["q05"], post["q95"], color=C["orange"], alpha=0.18, lw=0)
    cax.plot(D, quad["P_z_le_obs_given_D"], "o-", color=C["blue"], ms=3, lw=1.1)
    cax.axhline(0.025, color="k", lw=0.8, ls="--")
    cax.axhline(tr_ev["M0"]["P_le_obs_pooled"], color=C["grey"], lw=1.0, ls=":")
    cax.plot([tr_D["lo"], tr_D["hi"]], [tr_ev["P"]["P_le_obs_pooled"]] * 2, color=C["green"], lw=2.2, alpha=0.6)
    cax.plot(tr_D["D_phys"], tr_ev["P"]["P_le_obs_pooled"], "s", color=C["green"], ms=5)
    cax.set_xscale("log")
    cax.set_ylim(-0.012, 0.36)
    cax.set_xlabel(tr(lang, "mixing coefficient $D_{\\mathrm{air}}$ (m$^2$ h$^{-1}$)", "混合系数 $D_{\\mathrm{air}}$（m$^2$ h$^{-1}$）"))
    cax.set_ylabel(tr(lang, "$P(z\\leq z_{\\mathrm{obs}})$ under M2", "M2 下的 $P(z\\leq z_{\\mathrm{obs}})$"))
    cax.text(0.03, 0.97, tr(lang, f"restaurant posterior (90 %, shaded):\nintegrated {quad['equal_weight']:.3f} / {quad['pooled_weight']:.3f}",
                            f"餐厅后验（90%，阴影）：\n积分值 {quad['equal_weight']:.3f} / {quad['pooled_weight']:.3f}"),
             transform=cax.transAxes, fontsize=7.2, va="top", ha="left")
    cax.text(D[0], 0.040, tr(lang, "limit of the central 95 % interval", "中心 95% 区间的下限"), fontsize=6.8, va="bottom", ha="left")
    cax.text(D[0], tr_ev["M0"]["P_le_obs_pooled"] + 0.006, tr(lang, "well mixed", "均匀混合"), fontsize=6.8, va="bottom",
             ha="left", color="#444444")
    cax.text(tr_D["lo"] * 0.8, tr_ev["P"]["P_le_obs_pooled"], tr(lang, "tracer-based\ncoefficient", "示踪确定的\n系数"), fontsize=6.8,
             va="center", ha="right", color=C["green"])
    cax.text(-0.22, 1.02, "(c)", transform=cax.transAxes, fontsize=10, va="bottom", ha="left")
    return fig


if __name__ == "__main__":
    print(bilingual(draw, "s8_outbreaks_callcentre"))
