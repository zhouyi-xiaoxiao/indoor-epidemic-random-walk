"""s8_outbreaks_numbers.py -- collect, from the stored results of the research directories, every number quoted in
Section 7 that is not recomputed by s8_outbreaks_checks.py, and write them to one file so that the text can be
checked line by line.  Nothing is computed here except rounding-free extraction and three ratios.

Inputs (read only):
    data/bsc_validation/{01_cal_T4_primary,01_cal_R1_primary,03_holdout_shape,02b_m1_exact,
                                  04_level_controls,06_callcentre,07_loo,09_sensitivity,10_decision}.json
    data/bsc_outbreaks/tables/analysis_results.json
    data/validation_checks/callcentre_quadrature.json
Output: ../data/s8_outbreaks_numbers.json

    python s8_outbreaks_numbers.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
VAL = ROOT / "code/bsc_validation/data"
OB = ROOT / "code/bsc_outbreaks/data"
ART = HERE.parent / "data"


def load(p):
    return json.loads(Path(p).read_text())


out: dict = {}

# ---- Seoul call centre and per-hour hazards (model-free) ----------------------------------------------
a = load(OB / "analysis_results.json")
out["callcentre_desks"] = dict(digitized_desks=a["callcentre"]["digitized_desks"], north=a["callcentre"]["north"],
                               south=a["callcentre"]["south"], joincount=a["callcentre"]["joincount"])
out["hazard_per_hour_range"] = a["hazard_per_hour_range"]

# ---- calibration: mode and 90 % posterior interval of the mixing coefficient (flat prior on log10 D) ----
def post_summary(logD, ll):
    logD, ll = np.asarray(logD), np.asarray(ll)
    w = np.exp(ll - ll.max())
    w /= w.sum()
    c = np.cumsum(w)
    return dict(mode=float(10 ** logD[ll.argmax()]), q05=float(10 ** np.interp(0.05, c, logD)),
                q95=float(10 ** np.interp(0.95, c, logD)), max_logL=float(ll.max()))


t4 = load(VAL / "01_cal_T4_primary.json")
r1 = load(VAL / "01_cal_R1_primary.json")
ll2 = np.asarray(t4["M2"]["profile_logL"])          # grid of D x ventilation (removal-rate) nodes
m = ll2.max()
with np.errstate(divide="ignore"):
    ll2m = m + np.log(np.exp(ll2 - m).mean(axis=1))  # marginal over the ventilation prior (equal-weight nodes)
iD, ik = np.unravel_index(ll2.argmax(), ll2.shape)   # joint maximum: the value reported in Table V1 of this part
out["calibration"] = dict(
    trains=dict(n_cells=t4["n_cells"], M0_logL=t4["M0"]["logL"], logL_saturated=t4["logL_saturated"],
                M0_deviance=2 * (t4["logL_saturated"] - t4["M0"]["logL"]),
                M1=post_summary(t4["logD"], t4["M1"]["profile_logL"]), M2=post_summary(t4["logD"], ll2m)),
    restaurant=dict(n_tables=len(r1["tables"]), M0_logL=r1["M0"]["logL"],
                    M1=post_summary(r1["logD"], r1["M1"]["profile_logL"]),
                    M2=post_summary(r1["logD"], r1["M2"]["profile_logL"])))
out["calibration"]["trains"]["M2"]["note"] = ("mode, q05, q95, max_logL refer to the posterior marginalised over the "
                                              "ventilation prior; joint_* to the maximum over D and the removal-rate nodes")
out["calibration"]["trains"]["M2"]["joint_max_logL"] = float(ll2.max())
out["calibration"]["trains"]["M2"]["joint_max_D"] = float(10 ** np.asarray(t4["logD"])[iD])
out["calibration"]["trains"]["M2"]["joint_max_removal_rate"] = float(t4["kappa_nodes"][ik])
out["calibration"]["trains"]["M1"]["deviance"] = 2 * (t4["logL_saturated"] - out["calibration"]["trains"]["M1"]["max_logL"])
out["calibration"]["trains"]["M2"]["deviance"] = 2 * (t4["logL_saturated"] - float(ll2.max()))
out["calibration"]["trains"]["M1"]["gain_over_M0"] = out["calibration"]["trains"]["M1"]["max_logL"] - t4["M0"]["logL"]
out["calibration"]["trains"]["M2"]["gain_over_M0"] = float(ll2.max()) - t4["M0"]["logL"]

# ---- hold-out shape tests -------------------------------------------------------------------------------
h = load(VAL / "03_holdout_shape.json")
out["holdout"] = {ev: {mm: {q: h["events"][ev][mm][q] for q in ("p_two_sided", "rr_pred", "k_pred_mean", "log_score")}
                       for mm in ("M0", "M1", "M2")} | dict(obs=h["events"][ev]["obs"])
                  for ev in ("T1", "T2", "T5", "C1")}
out["holdout_pooled"] = {mm: {q: h["pooled"][mm][q] for q in ("delta", "p_under_M0", "per_event", "failed")}
                         for mm in ("M1", "M2")}

# ---- M1 reach ceiling -----------------------------------------------------------------------------------
ex = load(VAL / "02b_m1_exact.json")["events"]
out["m1_ceiling"] = {}
for ev, e in ex.items():
    rows = e["rows"]
    cannot = [r["mean_total_sim"] for r in rows if not r["reached_K"]]
    out["m1_ceiling"][ev] = dict(K=e["K"], draws=len(rows), draws_reaching=int(sum(r["reached_K"] for r in rows)),
                                 mean_total_cannot_median=float(np.median(cannot)) if cannot else None,
                                 mean_total_cannot_range=[float(min(cannot)), float(max(cannot))] if cannot else None,
                                 P_total_ge_K=float(np.mean([r["p_total_ge_K"] for r in rows])),
                                 Dp_of_first_row=rows[0]["D"])

# ---- levels ---------------------------------------------------------------------------------------------
lv = load(VAL / "04_level_controls.json")
out["levels"] = {k: {q: v[q] for q in ("k", "n", "T_h", "kappa_median", "E_implied", "E_ci_binomial") if q in v}
                 | ({"E_implied_alt": v["E_implied_alt"]} if "E_implied_alt" in v else {})
                 for k, v in lv["implied"].items()}
out["levels"]["T4_mean_emission_q05_q50_q95"] = lv["E_T4"]["q"]
out["levels"]["spread_max_over_min"] = lv["implied"]["H1"]["E_implied"] / lv["implied"]["T2"]["E_implied"]
out["levels"]["spread_confirmed_choir_only"] = lv["implied"]["R2"]["E_implied"] / lv["implied"]["T2"]["E_implied"]
out["level_checks"] = dict(tied_T2_T3_M2={q: lv["tied_T2_T3"]["M2"][q] for q in ("mean", "pi90", "inside")},
                           C3={k: {q: lv["C3"][k][q] for q in ("mean", "pi90", "inside")}
                               for k in ("M2_masks_0.35", "M2_masks_0.2", "M1_masks_0.35")},
                           T5_premium_economy=dict(bound_pct=lv["T5_PE"]["bound_pct"], M2=lv["T5_PE"]["M2"]))

# ---- call centre ----------------------------------------------------------------------------------------
cc = load(VAL / "06_callcentre.json")
out["callcentre_first_analysis_40_draws"] = {mm: dict(z_quantiles=cc[mm]["q"],  # 2.5, 5, 50, 95, 97.5 %
                                                      P_pooled=cc[mm]["P_z_le_obs"],
                                                      P_equal_weight=cc[mm]["P_z_le_obs_equal_weight"])
                                             for mm in ("M0", "M1", "M2")}
q = load(ART / "validation_checks/callcentre_quadrature.json")
out["callcentre_quadrature_M2"] = dict(equal_weight=q["equal_weight"], pooled_weight=q["pooled_weight"],
                                       prior_upper_limit_sensitivity=q["prior_upper_limit_sensitivity"])

# ---- sensitivity variants: in how many is A1 passed, and in how many does the coach fail? ---------------------
sv = load(VAL / "09_sensitivity.json")["variants"]
rows = {}
for tag, v in sv.items():
    pool = v["M2"]
    rows[tag] = dict(T1_p=v["events"]["T1"]["M2"]["p"], T2_p=v["events"]["T2"]["M2"]["p"],
                     T5_p=v["events"]["T5"]["M2"]["p"], C1_p=v["events"]["C1"]["M2"]["p"],
                     pooled_delta=pool["delta"], pooled_p=pool["p_under_M0"], A1=pool["A1"], failed=pool["failed"])
out["sensitivity_M2"] = dict(
    variants=rows, n_variants=len(rows),
    n_coach_fails=int(sum(r["T2_p"] < 0.05 for r in rows.values())),
    n_pooled_test_passes=int(sum(bool(r["A1"]) for r in rows.values())))

(ART / "s8_outbreaks_numbers.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
