"""Model-free analysis of the outbreak dataset (no model simulation is run here).

Computes and saves to data/analysis_results.json:
  1. near/far risk ratios with exact conditional tests for every stratified outbreak;
  2. Seoul call centre: wing contrast, join-count clustering statistic inside the north wing
     with a permutation null (seed fixed), epidemic growth rate from the digitized onset curve;
  3. per-event cumulative hazard and hazard per hour (level heterogeneity between events);
  4. upper bounds for the negative controls.
Figures are produced by make_figures.py from the JSON written here.
"""
import json, math, pathlib, sys

import numpy as np
import pandas as pd
from scipy import stats, optimize

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

SEED = 20261001
rng = np.random.default_rng(SEED)
res = {"seed": SEED}

ob = json.loads((ROOT / "outbreaks.json").read_text())["records"]
byid = {r["id"]: r for r in ob}
strata = pd.read_csv(ROOT / "data" / "spatial_strata.csv")


# ---------------------------------------------------------------- 1. risk ratios
def rr_ci(k1, n1, k0, n0):
    """Risk ratio near/far with a log-normal (Katz) 95% CI; 0.5 continuity correction if a zero cell."""
    a, b, c, d = k1, n1 - k1, k0, n0 - k0
    cc = 0.5 if min(a, b, c, d) == 0 else 0.0
    p1, p0 = (a + cc) / (n1 + 2 * cc), (c + cc) / (n0 + 2 * cc)
    rr = p1 / p0
    se = math.sqrt(1 / (a + cc) - 1 / (n1 + 2 * cc) + 1 / (c + cc) - 1 / (n0 + 2 * cc))
    return rr, rr * math.exp(-1.96 * se), rr * math.exp(1.96 * se), cc > 0


pairs = [
    ("O1_seoul_callcentre", "north wing desks (digitized)", "south wing desks (digitized)", "wing with most cases vs other wing (index position unknown)"),
    ("C1_marin_classroom", "front two rows (nearest teacher)", "back three rows", "rows nearest the index vs farther rows"),
    ("T1_zhejiang_bus", "rows 5-11 (within ~2 m of index row 8)", "rows 1-4 and 12-15", "within ~2 m (3 rows) vs beyond"),
    ("T1_zhejiang_bus", "rows 6-10 (within 2 rows)", "rows 1-5 and 11-15", "within 2 rows vs beyond (alternative split)"),
    ("T2_hunan_coach", "rear seats rows 8-13 (index zone)", "front seats rows 1-7", "index half vs other half of coach"),
    ("T5_vn54_flight", "business class, '≤2 seats away' (within ~2 m)", "business class, '>2 seats away'", "within ~2 m vs beyond, same cabin"),
    ("R1_guangzhou_restaurant", "immediate neighbouring tables (TB, TC, T18)", "remote tables (T04-T17)", "adjacent tables vs remote tables"),
    ("R1_guangzhou_restaurant", "ABC air-conditioning zone, non-A tables (TB, TC)", "other air-conditioning zones", "same airflow zone vs other zones"),
    ("C2_jerusalem_highschool", "grades 7-9 (junior wing)", "grades 10-12 (senior wing)", "building wing with index classes vs other wing"),
    ("S1_liaocheng_supermarket", "employees (shift-long exposure)", "customers (short visits)", "long-dwell staff vs short-dwell customers"),
]
rr_rows = []
for oid, near, far, desc in pairs:
    a = strata[(strata.outbreak_id == oid) & (strata.stratum == near)].iloc[0]
    b = strata[(strata.outbreak_id == oid) & (strata.stratum == far)].iloc[0]
    k1, n1, k0, n0 = int(a.n_infected), int(a.n_exposed), int(b.n_infected), int(b.n_exposed)
    rr, lo, hi, cc = rr_ci(k1, n1, k0, n0)
    p = stats.fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])[1]
    rr_rows.append(dict(outbreak_id=oid, contrast=desc, near=near, far=far, near_k=k1, near_n=n1, far_k=k0, far_n=n0,
                        near_ar=k1 / n1, far_ar=k0 / n0, risk_ratio=rr, rr_ci_low=lo, rr_ci_high=hi,
                        continuity_corrected=bool(cc), fisher_p=p))
res["risk_ratios"] = rr_rows
pd.DataFrame(rr_rows).to_csv(ROOT / "data" / "near_far_risk_ratios.csv", index=False)

# ---------------------------------------------------------------- 2. call centre
seats = pd.read_csv(ROOT / "data" / "park2020_callcentre_seats_digitized.csv")
cc = {}
nw = seats[seats.region == "north_wing"].reset_index(drop=True)
sw = seats[seats.region == "south_wing"]
cc["north"] = [int(nw.case.sum()), len(nw)]
cc["south"] = [int(sw.case.sum()), len(sw)]
cc["fisher_p_north_vs_south"] = float(stats.fisher_exact([[cc["north"][0], cc["north"][1] - cc["north"][0]],
                                                           [cc["south"][0], cc["south"][1] - cc["south"][0]]])[1])
cc["cases_with_seat_in_figure"] = int(seats.case.sum())
cc["cases_reported_floor"] = 94
# sensitivity: allocate the 10 unlocated cases (i) all to north wing desks, (ii) proportionally to non-digitized seats
cc["north_ar_if_all_10_unlocated_in_north"] = (cc["north"][0] + 10) / cc["north"][1]
cc["digitized_desks"] = int(len(nw) + len(sw))
cc["north_share_of_digitized_desks"] = len(nw) / (len(nw) + len(sw))
cc["north_share_of_216_employees"] = len(nw) / 216

# join-count statistic within the north wing: number of adjacent desk pairs that are both cases.
xy = nw[["x_pitch_units", "y_pitch_units"]].to_numpy()
d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
adj = (d > 0) & (d <= 1.25)          # side-by-side, face-to-face or front/back neighbours (<= 1.25 desk pitches)
iu = np.triu_indices(len(nw), 1)
pair_mask = adj[iu]
case = nw.case.to_numpy()


def joins(c):
    return int((c[iu[0]] & c[iu[1]] & pair_mask).sum())


obs = joins(case.astype(bool))
NPERM = 20000
null = np.empty(NPERM, dtype=int)
cb = case.astype(bool)
for i in range(NPERM):
    null[i] = joins(rng.permutation(cb))
cc["joincount"] = dict(adjacent_pairs=int(pair_mask.sum()), observed_case_case_pairs=obs, null_mean=float(null.mean()),
                       null_sd=float(null.std(ddof=1)), z=float((obs - null.mean()) / null.std(ddof=1)),
                       p_one_sided_greater=float((np.sum(null >= obs) + 1) / (NPERM + 1)), n_perm=NPERM,
                       neighbour_radius_pitch=1.25)
np.save(ROOT / "data" / "callcentre_joincount_null.npy", null)

# column-level heterogeneity in the north wing (desk blocks): chi-square of cases by desk column block
nw = nw.assign(block=np.round(nw.x_pitch_units / 2.5).astype(int))
tab = nw.groupby("block").case.agg(["sum", "count"])
chi2, p_chi, dof, _ = stats.chi2_contingency(np.vstack([tab["sum"], tab["count"] - tab["sum"]]))
cc["north_wing_blocks"] = dict(n_blocks=int(len(tab)), cases=tab["sum"].tolist(), desks=tab["count"].tolist(),
                               chi2=float(chi2), dof=int(dof), p=float(p_chi))
# west-east gradient within the north wing
rho, p_rho = stats.spearmanr(nw.x_pitch_units, nw.case)
cc["north_wing_x_gradient"] = dict(spearman_rho=float(rho), p=float(p_rho))

# epidemic growth rate on the 11th floor from the digitized onset curve (Poisson log-linear fit, growth phase)
epi = pd.read_csv(ROOT / "data" / "park2020_callcentre_epicurve_digitized.csv", parse_dates=["onset_date"])
full = pd.DataFrame({"onset_date": pd.date_range("2020-02-25", "2020-03-20")}).merge(epi, how="left").fillna(0)
full["t"] = (full.onset_date - pd.Timestamp("2020-02-25")).dt.days


def poisson_growth(df):
    t, y = df.t.to_numpy(float), df.floor11.to_numpy(float)

    def nll(p):
        mu = np.exp(p[0] + p[1] * t)
        return float(np.sum(mu - y * np.log(mu)))

    o = optimize.minimize(nll, [0.0, 0.2], method="BFGS")
    r = o.x[1]
    se = math.sqrt(o.hess_inv[1, 1])
    return r, se


growth = {}
for end in ("2020-03-06", "2020-03-07", "2020-03-09"):
    sub = full[full.onset_date <= end]
    r, se = poisson_growth(sub)
    growth[f"2020-02-25..{end}"] = dict(r_per_day=r, se=se, ci95=[r - 1.96 * se, r + 1.96 * se], doubling_time_d=math.log(2) / r,
                                        n_cases=int(sub.floor11.sum()))
cc["onset_growth"] = growth
cc["onset_total_floor11"] = int(epi.floor11.sum())
cc["onset_peak_window"] = "2020-03-06..2020-03-09 (58 of 86 dated 11th-floor onsets)"
cc["onsets_mar6_to_mar9"] = int(epi[(epi.onset_date >= "2020-03-06") & (epi.onset_date <= "2020-03-09")].floor11.sum())
res["callcentre"] = cc

# ---------------------------------------------------------------- 3. level heterogeneity
lev = []
for r in ob:
    if r.get("hazard_per_hour") is not None:
        lev.append(dict(id=r["id"], T_h=r["exposure_duration_h"], AR=r["attack_rate"], H=r["cumulative_hazard"], h_per_hour=r["hazard_per_hour"]))
res["single_event_levels"] = lev
# range over SARS-CoV-2 single-index events only (O5 is tuberculosis)
hs = np.array([x["h_per_hour"] for x in lev if byid[x["id"]]["pathogen"].startswith("SARS-CoV-2")])
res["hazard_per_hour_range"] = dict(min=float(hs.min()), max=float(hs.max()), ratio=float(hs.max() / hs.min()), n=len(hs),
                                    geometric_mean=float(np.exp(np.log(hs).mean())), gsd=float(np.exp(np.log(hs).std(ddof=1))),
                                    scope="SARS-CoV-2 single-index events")

# upper bounds for negative controls (exact one-sided 95%)
negs = [("S1 customers", 0, 8224), ("T1 unexposed bus", 0, 60), ("R1 remote tables", 0, 63), ("R1 other zones", 0, 68), ("T5 premium economy", 0, 35),
        ("C3 Utah classroom contacts", 5, 728), ("O1 floors 7-9", 1, 595)]
res["negative_controls"] = [dict(name=nm, k=k, n=n, ar=k / n, upper95_one_sided=float(stats.beta.ppf(0.95, k + 1, n - k))) for nm, k, n in negs]

(ROOT / "data" / "analysis_results.json").write_text(json.dumps(res, indent=1))
print(json.dumps({k: res[k] for k in ("risk_ratios", "callcentre", "hazard_per_hour_range", "negative_controls")}, indent=1))
