"""Additional checks of the outbreak dataset, and figures of the re-check.

B. Heterogeneity of the near/far risk ratios between single events (Cochran Q): is 'three regimes' a
   statistical finding?
C. Discriminating power of the decision rule of the protocol proposal (code/bsc_outbreaks/PROTOCOL_PROPOSAL.md).
Figures: figV2 (call centre), figV3 (train matrix), figV4 (near/far forest + protocol band).
Reads results/stats_check.json and results/callcentre_check.json written by v02 and v03.
"""
import json, math, csv
from pathlib import Path
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = Path(__file__).resolve().parents[1]
RES = V / "results"; FIG = V / "figures"; FIG.mkdir(exist_ok=True)
S = json.load(open(RES / "stats_check.json")); C = json.load(open(RES / "callcentre_check.json"))
out = {}

def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)

# ------------------------------------------------------------------ B. heterogeneity of near/far RR
EV = {r["id"]: r for r in S["risk_ratios"]}
use = ["C1", "T1a", "T2", "T5"]          # one contrast per single-index event with non-zero cells
y = np.array([math.log(EV[i]["rr"]) for i in use])
se = np.array([(math.log(EV[i]["hi"]) - math.log(EV[i]["lo"])) / (2 * 1.959964) for i in use])
w = 1 / se ** 2; pooled = float((w * y).sum() / w.sum()); Q = float((w * (y - pooled) ** 2).sum())
out["heterogeneity"] = dict(events=use, logRR=y.tolist(), se=se.tolist(), pooled_RR=math.exp(pooled),
    pooled_CI=[math.exp(pooled - 1.959964 / math.sqrt(w.sum())), math.exp(pooled + 1.959964 / math.sqrt(w.sum()))],
    Q=Q, dof=len(use) - 1, p=float(stats.chi2.sf(Q, len(use) - 1)),
    I2=max(0.0, (Q - (len(use) - 1)) / Q) if Q > 0 else 0.0)
pw = {}
for a in range(len(use)):
    for b in range(a + 1, len(use)):
        z = (y[a] - y[b]) / math.sqrt(se[a] ** 2 + se[b] ** 2)
        pw[f"{use[a]} vs {use[b]}"] = dict(z=float(z), p=float(2 * stats.norm.sf(abs(z))))
out["heterogeneity"]["pairwise"] = pw
# exact (conditional) homogeneity is not needed; also give Breslow-Day-free check: bus vs classroom
out["heterogeneity"]["note"] = "inverse-variance weights on log RR (Katz SE); T1 uses the rows 5-11 split"

# ------------------------------------------------------------------ C. power of the decision rule
lo_all = max(EV[i]["lo"] for i in use); hi_all = min(EV[i]["hi"] for i in use)
out["protocol_power"] = dict(
    constant_RR_band_passing_all_four_single_event_shape_tests=[lo_all, hi_all],
    wellmixed_RR1_z={i: float(-math.log(EV[i]["rr"]) / s) for i, s in zip(use, se)},
    wellmixed_fails=[i for i, s in zip(use, se) if abs(math.log(EV[i]["rr"]) / s) > 1.959964])
g = S["hazard_per_hour"]["gsd"]; f90 = g ** 1.6449
ex = {}
for ar_med in (0.05, 0.2, 0.5):
    H = -math.log(1 - ar_med); ex[str(ar_med)] = [1 - math.exp(-H / f90), 1 - math.exp(-H * f90)]
out["protocol_power"]["level_test"] = dict(gsd=g, factor_90pct=f90, predictive_AR_interval_for_median_AR=ex,
    note="if event infectiousness is log-normal with the GSD observed across the nine events, a 90% predictive interval "
         "for one held-out event spans this multiplicative factor either side of the median hazard")
# how many of the nine observed events would fall inside a 90% interval centred on the geometric mean hazard?
hz = S["hazard_per_hour"]["values"]; gm = S["hazard_per_hour"]["geometric_mean"]
out["protocol_power"]["events_inside_90pct_band_around_GM"] = sum(1 for v in hz.values() if gm / f90 <= v <= gm * f90)

json.dump(out, open(RES / "extra_checks.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))

# ================================================================== figures
BLUE, ORANGE, AQUA, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#d9d8d2"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "font.family": "DejaVu Sans",
                     "axes.titlesize": 9.5, "axes.titleweight": "bold", "legend.frameon": False, "pdf.fonttype": 42})
def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight"); fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight"); plt.close(fig)

# ---- figV2: call centre
seat = list(csv.DictReader(open(RES / "callcentre_manual_seatmap.csv")))
null = np.load(RES / "callcentre_joincount_null_rookdiag.npy")
fig = plt.figure(figsize=(9.2, 6.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], hspace=0.38, wspace=0.28)
axs = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
ax = axs[0]
for d in seat:
    x, yv, c = float(d["x"]), float(d["y"]), int(d["case"])
    wdt, hgt = (34, 15) if d["col"] == "H" else (16, 28)
    ax.add_patch(plt.Rectangle((x - wdt / 2, yv - hgt / 2), wdt, hgt, facecolor=BLUE if c else "white", edgecolor=MUTED, lw=0.6))
ax.set_xlim(0, 890); ax.set_ylim(262, 30); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("a  North wing: 79 case desks of 137 (re-check's reading of Park et al. 2020, Fig. 2)", loc="left")
ax.text(0, 268, "filled = seat of a confirmed case; positions in pixels of the published 900-px figure (no scale bar)", fontsize=7.5, color=MUTED, va="top")
ax = axs[1]
ax.hist(null, bins=np.arange(null.min() - 0.5, null.max() + 1.5, 2), color=GRID, edgecolor="white", lw=0.4)
ob = C["joincount"]["rook+diag"]["observed"]
ax.axvline(ob, color=BLUE, lw=2); ax.text(ob - 1, ax.get_ylim()[1] * 1.03, f"observed {ob}", color=INK, ha="right", fontsize=8)
for eps, col in (("1.0", ORANGE),):
    m, s = C["power_contagion_alternative"]["by_eps"][eps]["mean_joins"], C["power_contagion_alternative"]["by_eps"][eps]["sd"]
    ax.axvspan(m - s, m + s, color=ORANGE, alpha=0.18, lw=0); ax.axvline(m, color=ORANGE, lw=1.5, ls="--")
    ax.text(m - 1, ax.get_ylim()[1] * 0.62, "neighbour-driven\nalternative (ε = 1)\nmean ± SD", color=INK, fontsize=7.5, ha="right")
ax.set_xlabel("Case–case adjacent desk pairs (286 pairs)"); ax.set_ylabel("Permutations (of 20,000)")
ax.set_title("b  Join count vs random relabelling", loc="left")
ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
ax = axs[2]
eps = sorted(C["power_contagion_alternative"]["by_eps"], key=float)
pw_ = [C["power_contagion_alternative"]["by_eps"][e]["power_at_alpha05"] for e in eps]
share = [1 / (1 + float(e)) for e in eps]
ax.plot([float(e) for e in eps], pw_, "-o", color=BLUE, lw=2, ms=6)
for e, p in zip(eps, pw_): ax.annotate(f"{p:.2f}", (float(e), p), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=7.5, color=INK)
ax.set_xscale("log"); ax.set_ylim(0, 1.12); ax.set_xlabel("Background weight ε (each infected neighbour has weight 1)")
ax.set_ylabel("Power, one-sided test at α = 0.05"); ax.yaxis.grid(True, color=GRID, lw=0.6)
ax.set_title("c  Power against neighbour-driven spread", loc="left")
save(fig, "figV2_callcentre")

# ---- figV3: train matrix
tm = S["train_matrix"]; cells = tm["cells"]
fig, axs = plt.subplots(2, 1, figsize=(7.4, 6.6), gridspec_kw={"height_ratios": [1.1, 1], "hspace": 0.5})
ax = axs[0]
M = np.full((4, 6), np.nan); N = np.full((4, 6), np.nan); K = np.full((4, 6), np.nan)
for c in cells:
    col = c["printed_col"] + 1 if c["row"] == 0 else c["printed_col"]
    M[c["row"], col] = c["pct"]; N[c["row"], col] = c["n"]; K[c["row"], col] = c["k"]
im = ax.imshow(np.log10(M), cmap="Blues", aspect="auto", vmin=-1.6, vmax=0.6)
for r in range(4):
    for c in range(6):
        if np.isnan(M[r, c]): ax.text(c, r, "index", ha="center", va="center", color=MUTED, fontsize=8); continue
        ax.text(c, r, f"{M[r, c]:.2f}%\n{int(K[r, c])}/{int(N[r, c])}", ha="center", va="center", fontsize=7.6,
                color="white" if M[r, c] > 1 else INK)
ax.set_xticks(range(6)); ax.set_yticks(range(4)); ax.set_xlabel("Columns apart (aisle counted as a seat)"); ax.set_ylabel("Rows apart")
ax.set_title("a  Hu et al. Table 1: attack rate by seat offset, with (k, n) recovered from the printed CIs", loc="left")
for s in ax.spines.values(): s.set_visible(False)
ax = axs[1]
x = np.arange(6); pm = [tm["printed_column_means"][str(c)] for c in range(6)]
a_ = [tm["column_rates_printed_layout"][str(c)] for c in range(6)]; b_ = [tm["column_rates_row0_shifted"][str(c)] for c in range(6)]
ax.bar(x - 0.27, pm, 0.25, color=INK, label="printed column mean"); ax.bar(x, b_, 0.25, color=BLUE, label="pooled, same-row cells at columns 1–5")
ax.bar(x + 0.27, a_, 0.25, color=ORANGE, label="pooled, same-row cells at columns 0–4 (as on publisher page)")
ax.set_xticks(x); ax.set_xlabel("Columns apart"); ax.set_ylabel("Attack rate, %"); ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_ylim(0, 1.0); ax.legend(fontsize=7.6, loc="upper right"); ax.set_title("b  Which placement of the same-row cells reproduces the printed column means", loc="left")
save(fig, "figV3_train_matrix")

# ---- figV4: forest of near/far RR + protocol band
lab = {"O1": "Seoul call centre: north vs south wing", "C2": "Jerusalem school: junior vs senior wing", "T5": "Flight VN54: <2 vs >2 seats",
       "C1": "Marin classroom: front 2 vs back 3 rows", "T1b": "Zhejiang bus: rows 6–10 vs others", "T1a": "Zhejiang bus: rows 5–11 vs others",
       "T2": "Hunan coach: rear vs front rows", "T2side": "Hunan coach: driver side vs opposite"}
ordr = ["O1", "C2", "T5", "C1", "T1b", "T1a", "T2", "T2side"]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
ax.axvspan(lo_all, hi_all, color=ORANGE, alpha=0.16, lw=0)
ax.axvline(1, color=MUTED, lw=0.8)
for i, k in enumerate(ordr[::-1]):
    r = EV[k]; single = k in use
    ax.plot([r["lo"], r["hi"]], [i, i], color=INK if single else MUTED, lw=1.3)
    ax.plot(r["rr"], i, "o" if single else "s", color=BLUE if single else "white", mec=BLUE, ms=6, zorder=3)
    ax.text(60, i, f"{r['near']} vs {r['far']}", va="center", fontsize=7.6, color=MUTED)
ax.set_xscale("log"); ax.set_xlim(0.2, 50); ax.set_yticks(range(len(ordr))); ax.set_yticklabels([lab[k] for k in ordr[::-1]], color=INK)
ax.set_xticks([0.25, 0.5, 1, 2, 4, 8, 16, 32]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8", "16", "32"])
ax.set_xlabel("Risk ratio near/far (Katz 95% CI)")
ax.text(math.sqrt(lo_all * hi_all), len(ordr) - 0.45, f"any constant RR in {lo_all:.2f}–{hi_all:.2f}\npasses all four single-event tests",
        ha="center", va="bottom", fontsize=7.6, color=INK)
ax.set_ylim(-0.6, len(ordr) + 0.4)
ax.set_title("Near/far contrasts: filled = the four single-index events used in the protocol's shape test", loc="left")
save(fig, "figV4_near_far")
print("figures written")
