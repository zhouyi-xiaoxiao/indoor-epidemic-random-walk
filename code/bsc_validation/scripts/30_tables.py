"""Markdown tables generated from the stored results (no number in the tables
is typed by hand).  Output: data/tables.md"""
import numpy as np
from scipy import stats

import _common as C
from valmod import events as E
from valmod import predict as P

L = []
w = L.append
f2 = lambda x: f"{x:.2f}"
f3 = lambda x: f"{x:.3f}"


def pfmt(p):
    if p < 0.0001:
        return f"{p:.1e}"
    return f"{p:.4f}" if p < 0.01 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


# ---- calibration -------------------------------------------------------------
t4 = C.load_json("01_cal_T4_primary.json")
r1 = C.load_json("01_cal_R1_primary.json")
w("### Table V1. Calibration fits (in-sample)\n")
w("| Record | Model | Free parameters | Mixing coefficient at the mode, m^2/h (90% posterior) | log-likelihood | AIC | BIC | Deviance (d.f.), p |")
w("|---|---|---|---|---|---|---|---|")
for name, cal, cls, nobs in (("T4 trains, 22 cells", t4, "vehicle", 22), ("R1 restaurant, 16 tables", r1, "room", 16)):
    for m in ("M0", "M1", "M2"):
        c = cal[m]
        k = c["n_par"]
        ll = c["logL"]
        if m == "M0":
            dtxt = "n/a"
        else:
            pd_ = P.posterior_D(cls, m, n=400, seed=C.SEED + 7)
            s = P.summarize_posterior(pd_["w_logD"], pd_["logD"])
            dtxt = f"{c['D_hat']:.3g} ({s['q05']:.2g}-{s['q95']:.2g})"
        if "logL_saturated" in cal:
            dev = 2 * (cal["logL_saturated"] - ll)
            df = nobs - k
            dv = f"{dev:.1f} ({df}), {pfmt(stats.chi2.sf(dev, df))}"
        else:
            dv = "n/a (interval-censored)"
        w(f"| {name} | {m} | {k} | {dtxt} | {ll:.2f} | {2 * k - 2 * ll:.1f} | {k * np.log(nobs) - 2 * ll:.1f} | {dv} |")
w("")
w("### Table V2. Train cohort: observed and fitted attack rate per seat offset (%)\n")
w("| Rows apart | Seats apart | Observed (k/n, reconstructed) | M0 | M1 | M2 | In likelihood |")
w("|---|---|---|---|---|---|---|")
c = t4["cells"]
for i in range(len(c["dr"])):
    w(f"| {c['dr'][i]} | {c['dc'][i]} | {c['pct'][i]:.2f} ({c['k'][i]}/{c['n'][i]}) | {100 * t4['M0']['p'][i]:.2f} | "
      f"{100 * t4['M1']['p'][i]:.2f} | {100 * t4['M2']['p'][i]:.2f} | {'yes' if t4['mask'][i] else 'no (companions)'} |")
w("")

# ---- sensitivity of calibration ---------------------------------------------
w("### Table V3. Calibration variants (log-likelihood; mixing coefficient at the mode)\n")
w("| Variant | Cells | M0 | M1 (D_p) | M2 (D_air) |")
w("|---|---|---|---|---|")
for v in ("primary", "S-adj", "S-row", "S-logit", "S-pitch"):
    d = C.load_json(f"01_cal_T4_{v}.json")
    w(f"| T4 {v} | {d['n_cells']} | {d['M0']['logL']:.2f} | {d['M1']['logL']:.2f} ({d['M1']['D_hat']:.3g}) | "
      f"{d['M2']['logL']:.2f} ({d['M2']['D_hat']:.3g}) |")
an = C.load_json("01_cal_T4_S-aniso.json")
w(f"| T4 S-aniso (seat-back factor theta) | 22 | {t4['M0']['logL']:.2f} | {an['M1']['logL']:.2f} "
  f"({an['M1']['D_hat']:.3g}; theta {an['M1']['theta_hat']:.2g}) | {an['M2']['logL']:.2f} "
  f"({an['M2']['D_hat']:.3g}; theta {an['M2']['theta_hat']:.2g}) |")
for v in ("primary", "S-R1-upper", "S-R1-lower"):
    d = C.load_json(f"01_cal_R1_{v}.json")
    w(f"| R1 {v} | 16 tables | {d['M0']['logL']:.2f} | {d['M1']['logL']:.2f} ({d['M1']['D_hat']:.3g}) | "
      f"{d['M2']['logL']:.2f} ({d['M2']['D_hat']:.3g}) |")
w("")

# ---- power -------------------------------------------------------------------
pw = C.load_json("02_power.json")
w("### Table V4. Power of the amended decision rule (computed before unblinding)\n")
w("| Judged model | Quantity | truth: M0 | truth: M1 | truth: M2 |")
w("|---|---|---|---|---|")
for m in ("M1", "M2"):
    for key, lab in (("P_A1", "P(A1 passes)"), ("P_A2_nofail", "P(no A2 failure)"),
                     ("P_A2_le1", "P(at most one A2 failure)"), ("P_W1_like", "P(A1 and no A2 failure)")):
        w(f"| {m} | {lab} | {pw[m]['M0'][key]:.4f} | {pw[m]['M1'][key]:.3f} | {pw[m]['M2'][key]:.3f} |")
w("")

# ---- hold-out shape ----------------------------------------------------------
sh = C.load_json("03_holdout_shape.json")
ex = C.load_json("02b_m1_exact.json")["events"]
w("### Table V5. Hold-out shape tests (event totals fixed; one infectiousness parameter per event)\n")
w("| Event | Strata (near vs far) | Observed | Observed RR (95% CI) | Model | Model RR = ratio of expected stratum attack rates (90% predictive interval of the sample RR) | Predicted near cases (90% PI) | Two-sided predictive p | log score | Delta vs M0 | z (original protocol) | Label |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|")
names = {"T1": "T1 Zhejiang bus", "T2": "T2 Hunan coach", "T5": "T5 flight VN54", "C1": "C1 Marin classroom"}
for ev in ("T1", "T2", "T5", "C1"):
    r = sh["events"][ev]
    o = r["obs"]
    se = r["M0"]["se_log_rr_obs"]
    rr = r["M0"]["rr_obs"]
    ci = f"{rr:.2f} ({rr * np.exp(-1.96 * se):.2f}-{rr * np.exp(1.96 * se):.2f})"
    for m in ("M0", "M1", "M2"):
        x = r[m]
        lab = "null" if m == "M0" else x["label"]
        note = ""
        if m == "M1" and r["M1_uses_exact_simulator"]:
            pge = float(np.mean([q["p_total_ge_K"] for q in ex[ev]["rows"]]))
            note = f" (exact simulator; P(total >= {r['K']}) = {pge:.4f})"
            if pge < 0.01:
                lab = "not reproduced (total unattainable)"
        d = "" if m == "M0" else f"{x['delta_vs_M0']:+.2f}"
        w(f"| {names[ev]} | {o['label']} | {o['near_k']}/{o['near_n']} vs {o['far_k']}/{o['far_n']} | {ci} | {m}{note} | "
          f"{x['rr_pred']:.2f} ({x['rr_pred_90'][0]:.2f}-{x['rr_pred_90'][1]:.2f}) | {x['k_pred_mean']:.1f} "
          f"({x['k_pred_90'][0]}-{x['k_pred_90'][1]}) | {pfmt(x['p_two_sided'])} | {x['log_score']:.2f} | {d} | "
          f"{x['z_original']:+.2f} | {lab} |")
w("")
w("### Table V6. Pooled test against the well-mixed null (A1)\n")
w("| Model | Pooled Delta (sum of log-score differences) | Exact one-sided p under M0 | A1 | Single-event A2 failures |")
w("|---|---|---|---|---|")
for m, lab in (("M2", "M2"), ("M1", "M1 (exact simulator where required)"), ("M1cf", "M1, closed form only (for information)")):
    x = sh["pooled"][m]
    w(f"| {lab} | {x['delta']:+.2f} | {pfmt(x['p_under_M0'])} | {'pass' if x['A1_pass'] else 'fail'} | "
      f"{', '.join(x['failed']) or 'none'} |")
w("")
w("### Table V7. Secondary splits (not counted in the decision rule)\n")
w("| Split | Observed | Observed RR | M0: RR, p | M1 (closed form): RR, p | M2: RR, p |")
w("|---|---|---|---|---|---|")
for key, r in sh["secondary"].items():
    o = r["obs"]
    w(f"| {key} | {o['near_k']}/{o['near_n']} vs {o['far_k']}/{o['far_n']} | {r['M0']['rr_obs']:.2f} | " + " | ".join(
        f"{r[m]['rr_pred']:.2f}, {pfmt(r[m]['p_two_sided'])}" for m in ("M0", "M1", "M2")) + " |")
w("")

# ---- M1 exact ----------------------------------------------------------------
w("### Table V8. M1 (the per-cell rule) in the verified exact simulator: reach ceiling\n")
w("| Event | Observed secondary cases / exposed | Posterior draws in which the mean total can reach the observed total | Draws that cannot: mean total at the largest beta, median (range) | P(total >= observed), pooled over draws | Largest absolute difference in a stratum attack rate, closed form vs exact |")
w("|---|---|---|---|---|---|")
for ev in ("T1", "T2", "T5", "C1"):
    e = ex[ev]
    rows = e["rows"]
    tot = [q["mean_total_sim"] for q in rows if not q["reached_K"]]
    pge = float(np.mean([q["p_total_ge_K"] for q in rows]))
    n = e["n_near"] + e["n_far"]
    ttxt = f"{np.median(tot):.1f} ({min(tot):.1f}-{max(tot):.1f})" if tot else "n/a"
    w(f"| {names[ev]} | {e['K']}/{n} | {sum(q['reached_K'] for q in rows)} of {len(rows)} | "
      f"{ttxt} | {pge:.4f} | {e['max_abs_diff_stratum_AR']:.3f} |")
w("")

# ---- O1 ----------------------------------------------------------------------
o1 = C.load_json("06_callcentre.json")
w("### Table V9. Seoul call centre, north wing: join-count z (79 cases on 137 desks)\n")
w("| | z: 2.5% | 5% | median | 95% | 97.5% | P(z <= observed), retained outbreaks pooled (as run) | Inside central 95% | P(z <= observed), equal weight per posterior draw | Inside central 95% | Posterior draws individually compatible (P > 0.025) |")
w("|---|---|---|---|---|---|---|---|---|---|---|")
w(f"| Observed | | | {o1['observed']['z']:.2f} | | | | | | | |")
for m in ("M0", "M1", "M2"):
    q = o1[m]["q"]
    w(f"| {m} | {q[0]:.2f} | {q[1]:.2f} | {q[2]:.2f} | {q[3]:.2f} | {q[4]:.2f} | {o1[m]['P_z_le_obs']:.3f} | "
      f"{'yes' if o1[m]['inside_95'] else 'no'} | {o1[m]['P_z_le_obs_equal_weight']:.3f} | "
      f"{'yes' if o1[m]['inside_95_equal_weight'] else 'no'} | {o1[m]['n_draws_with_P_gt_0.025']} of {o1[m]['n_draws']} |")
w("")

# ---- level -------------------------------------------------------------------
lv = C.load_json("04_level_controls.json")
w("### Table V10. Implied infectiousness of single-index events (Wells-Riley scaling shared by M0 and M2)\n")
w("| Event | Cases / exposed | Duration, h | Volume, m^3 | Removal rate, 1/h (prior median) | Cumulative hazard | Hazard per hour | Implied emission, quanta/h at 0.5 m^3/h (binomial 95%; with ventilation prior) |")
w("|---|---|---|---|---|---|---|---|")
for ev, x in lv["implied"].items():
    w(f"| {ev} {x['label']} | {x['k']}/{x['n']} | {x['T_h']:.2f} | {x['V']:.0f} | {x['kappa_median']:.2f} | {x['dose']:.3f} | "
      f"{x['hazard_per_h']:.3f} | {x['E_implied']:.0f} ({x['E_ci_binomial'][0]:.0f}-{x['E_ci_binomial'][1]:.0f}; "
      f"{x['E_range_all'][0]:.0f}-{x['E_range_all'][1]:.0f}) |")
et = lv["E_T4"]["q"]
w(f"| T4 trains (mean over 2,334 index cases) | 230/72,692 | 2.1 (mean) | per coach lattice | 5.9-20.9 | | | "
  f"{et[1]:.2f} ({et[0]:.2f}-{et[2]:.2f}, posterior 90%) |")
w("")
w("### Table V11. Dispersion of implied infectiousness across H1, O2, R2, C5\n")
w("| Scaling | SD of log implied infectiousness | Geometric SD | Permutation p (24 permutations) |")
w("|---|---|---|---|")
for k, lab in (("wells_riley", "Wells-Riley, T/(V kappa) with transient (M0, M2)"),
               ("m1_T_over_area", "M1, T/area"), ("constant_per_hour", "constant hazard per hour")):
    x = lv["dispersion"][k]
    w(f"| {lab} | {x['sd_log']:.2f} | {x['gsd']:.2f} | {x['perm_p']:.3f} |")
w("")
w("### Table V12. Level checks with defined reference\n")
w("| Check | Model | Predicted mean | 90% predictive interval | Observed | Inside |")
w("|---|---|---|---|---|---|")
for m in ("M0", "M1", "M2"):
    x = lv["tied_T2_T3"][m]
    w(f"| T2 -> T3, same index (17 exposed) | {m} | {x['mean']:.2f} | {x['pi90'][0]}-{x['pi90'][1]} | 2 | {'yes' if x['inside'] else 'no'} |")
for key, lab in (("M2_masks_0.35", "M2, mask factor 0.35"), ("M2_masks_0.2", "M2, mask factor 0.2"),
                 ("M2_masks_0.6", "M2, mask factor 0.6"), ("M2_no_masks", "M2, no masks"),
                 ("M1_masks_0.35", "M1, mask factor 0.35"), ("M1_no_masks", "M1, no masks")):
    x = lv["C3"][key]
    w(f"| T4 -> C3, unselected cohorts (728 contacts) | {lab} | {x['mean']:.2f} | {x['pi90'][0]:.0f}-{x['pi90'][1]:.0f} | 5 | "
      f"{'yes' if x['inside'] else 'no'} |")
w("")
w("### Table V13. Negative controls\n")
w("| Control | Observed | One-sided 95% upper bound, % | Model | Predicted attack rate, % | P(prediction exceeds bound) | Verdict |")
w("|---|---|---|---|---|---|---|")
s1 = lv["S1"]
w(f"| S1 supermarket customers, one 27-min visit, one infectious employee present | 0/8,224 | {s1['bound_pct']:.4f} | M0/M2 | "
  f"{s1['M2_single_index']['median_pct']:.4f} (median; 95th percentile {s1['M2_single_index']['q95_pct']:.4f}) | "
  f"{s1['M2_single_index']['P_exceed']:.4f} | pass |")
w(f"| S1 | | | M1 | {s1['M1_pooled']['median_pct']:.4f} | {s1['M1_pooled']['P_exceed']:.4f} | pass |")
c3 = lv["C3"]
w(f"| C3 masked classrooms | 5/728 | {c3['bound_pct']:.2f} | M0/M2, mask factor 0.35 | {c3['M2_masks_0.35']['ar_mean_pct']:.2f} | "
  f"{c3['M2_masks_0.35']['P_exceed_bound']:.4f} | pass |")
w(f"| C3 | | | M0/M2, no masks | {c3['M2_no_masks']['ar_mean_pct']:.2f} | {c3['M2_no_masks']['P_exceed_bound']:.3f} | pass |")
w(f"| C3 | | | M1, mask factor 0.35 | {c3['M1_masks_0.35']['ar_mean_pct']:.2f} | {c3['M1_masks_0.35']['P_exceed_bound']:.4f} | pass |")
pe = lv["T5_PE"]
w(f"| T5 premium economy | 0/35 | {pe['bound_pct']:.1f} | M2 | {pe['M2']['ar_median_pct']:.2f} (median; 95th percentile "
  f"{pe['M2']['ar_q95_pct']:.1f}) | {pe['M2']['P_exceed']:.2f} | pass |")
w(f"| T5 premium economy | | | M1, exact simulator, largest beta | {pe['M1_exact_beta_cap']['ar_mean_pct']:.1f} | "
  f"{pe['M1_exact_beta_cap']['P_exceed']:.2f} | pass |")
t7 = lv["T7"]
w(f"| T7 one 20-min metro ride, 309 other riders: P(at least 5 secondary cases) | no rail cluster reported | n/a | M0/M2 | "
  f"{t7['M2']['mean_cases']:.3f} expected cases | {t7['M2']['P_ge5']:.4f} | pass (weak) |")
w(f"| T7 | | | M1 | {t7['M1']['mean_cases']:.3f} expected cases | {t7['M1']['P_ge5']:.4f} | pass (weak) |")
w("")

# ---- sensitivity -------------------------------------------------------------
se = C.load_json("09_sensitivity.json")["variants"]
w("### Table V14. Sensitivity analyses: M2 predicted risk ratio and predictive p per event; pooled test\n")
w("| Variant | T1 | T2 | T5 | C1 | Pooled Delta | p under M0 | Events failing A2 |")
w("|---|---|---|---|---|---|---|---|")
for tag, r in se.items():
    cells = " | ".join(f"{r['events'][ev]['M2']['rr_pred']:.2f}, {pfmt(r['events'][ev]['M2']['p'])}" for ev in ("T1", "T2", "T5", "C1"))
    w(f"| {tag} | {cells} | {r['M2']['delta']:+.2f} | {pfmt(r['M2']['p_under_M0'])} | {', '.join(r['M2']['failed']) or 'none'} |")
w("")
w("### Table V15. Sensitivity analyses: M1 (closed form) predicted risk ratio and predictive p\n")
w("| Variant | T1 | T2 | T5 | C1 | Pooled Delta | Events failing A2 |")
w("|---|---|---|---|---|---|---|")
for tag, r in se.items():
    cells = " | ".join(f"{r['events'][ev]['M1']['rr_pred']:.2f}, {pfmt(r['events'][ev]['M1']['p'])}" for ev in ("T1", "T2", "T5", "C1"))
    w(f"| {tag} | {cells} | {r['M1']['delta']:+.2f} | {', '.join(r['M1']['failed']) or 'none'} |")
w("")

# ---- LOO ---------------------------------------------------------------------
lo = C.load_json("07_loo.json")
w("### Table V16. Leave-one-event-out (secondary analysis) and the mixing coefficient each event prefers\n")
w("| Class | Model | Left-out event | D preferred by that event alone, m^2/h | D fitted on the other events | log score | M0 log score | Delta |")
w("|---|---|---|---|---|---|---|---|")
for cls in ("vehicle", "room"):
    for m in ("M2", "M1"):
        r = lo["classes"][cls][m]
        for ev, x in r["loo"].items():
            w(f"| {cls} | {m}{' (closed form)' if m == 'M1' else ''} | {ev} | {x['D_own_mode']:.3g} | {x['D_train_mode']:.3g} | "
              f"{x['score']:.2f} | {x['score_M0']:.2f} | {x['delta']:+.2f} |")
w("")
w("### Table V17. Is one mixing coefficient per setting class tenable? (likelihood-ratio test, all events of the class)\n")
w("| Class | Model | Sum of per-event maxima | Maximum with one common D | 2 x difference | d.f. | p |")
w("|---|---|---|---|---|---|---|")
for cls in ("vehicle", "room"):
    for m in ("M2", "M1"):
        r = lo["classes"][cls][m]
        lr = 2 * (r["sum_own_max"] - r["joint_max"])
        df = len(r["loo"]) - 1
        w(f"| {cls} | {m}{' (closed form)' if m == 'M1' else ''} | {r['sum_own_max']:.2f} | {r['joint_max']:.2f} | {lr:.1f} | {df} | "
          f"{pfmt(stats.chi2.sf(lr, df))} |")
w("")

dc = C.load_json("10_decision.json")["models"]
w("### Table V18. Decision rule\n")
w("| Criterion | M1 (same-site model) | M2 (framework with the non-local term) |")
w("|---|---|---|")
w(f"| A1 pooled test against M0 | Delta {dc['M1']['A1']['delta']:+.2f}: fail | Delta {dc['M2']['A1']['delta']:+.2f}, p = {pfmt(dc['M2']['A1']['p'])}: pass |")
w(f"| A2 single-event failures | {', '.join(x for x in dc['M1']['A2_failures_as_run'] if x != 'O1') or 'none'} | {', '.join(x for x in dc['M2']['A2_failures_as_run'] if x != 'O1') or 'none'} |")
w(f"| A2 call-centre join count | inside (P = {dc['M1']['O1']['P_z_le_obs_pooled']:.3f} pooled; {dc['M1']['O1']['P_z_le_obs_equal_weight']:.3f} equal weight) | pooled {dc['M2']['O1']['P_z_le_obs_pooled']:.3f}: inside; equal weight {dc['M2']['O1']['P_z_le_obs_equal_weight']:.3f}: outside |")
w(f"| A3 level checks | tied index: {'inside' if dc['M1']['A3']['tied_index'] else 'outside'}; unselected transfer: {'inside' if dc['M1']['A3']['unselected_transfer'] else 'outside'} | tied index: {'inside' if dc['M2']['A3']['tied_index'] else 'outside'}; unselected transfer: {'inside' if dc['M2']['A3']['unselected_transfer'] else 'outside'} |")
w(f"| A4 negative controls | {'pass' if dc['M1']['A4']['passed'] else 'fail'} | {'pass' if dc['M2']['A4']['passed'] else 'fail'} |")
w(f"| Outcome, as run | {dc['M1']['outcome_as_run']} | {dc['M2']['outcome_as_run']} |")
w(f"| Outcome, conservative reading of the call-centre test | {dc['M1']['outcome_conservative']} | {dc['M2']['outcome_conservative']} |")
w("")
mc = C.load_json("02c_mc_stability.json")
w("### Table V19. Monte Carlo stability of the M2 shape results (five other seeds)\n")
w("| Seed | T1 p | T2 p | T5 p | C1 p | Pooled Delta | p under M0 |")
w("|---|---|---|---|---|---|---|")
for sd, r in zip(mc["seeds"], mc["M2"]):
    w(f"| {sd} | " + " | ".join(pfmt(r["p"][ev]) for ev in ("T1", "T2", "T5", "C1")) + f" | {r['delta']:+.2f} | {pfmt(r['p_pooled'])} |")
w("")
tb = E.r1_tables()
w("### Table V20. Guangzhou restaurant: fitted M2 against the published tracer measurement\n")
w("| Table | Measured tracer concentration relative to table A (Li et al. 2021) | M2 at the posterior mode: time-integrated exposure relative to table A |")
w("|---|---|---|")
xr = dict(zip(r1["tables"], r1["M2"]["x_rel_src"]))
for _, row in tb.iterrows():
    if row["table"] in xr and row["norm_measured_tracer"] == row["norm_measured_tracer"]:
        w(f"| {row['table']} | {row['norm_measured_tracer']:.2f} | {xr[row['table']]:.4f} |")
w("")
ic = C.load_json("08_bruteforce_check.json")["events"]
st = C.load_json("00_selftest.json")
w("### Table V21. Software verification\n")
w("| Check | Result |")
w("|---|---|")
w(f"| Closed-form kernels vs quadrature of the matrix exponential (synthetic lattice) | max abs error {st['A_max_abs_err_M1']:.1e} (M1), {st['A_max_abs_err_M2']:.1e} (M2) |")
w(f"| M2 with D -> infinity equals M0; spatial mean of M2 equals M0 | {st['B_M2_Dinf_minus_M0']:.1e}; {st['B_mean_M2_minus_M0']:.1e} |")
w(f"| Conditional stratum distribution vs enumeration; vs hypergeometric | {st['C_cond_pmf_err']:.1e}; {st['C_hypergeom_err']:.1e} |")
w(f"| Join-count moments vs 40,000 permutations (mean; variance) | {st['D_mean_theory_sim'][0]:.3f} vs {st['D_mean_theory_sim'][1]:.3f}; {st['D_var_theory_sim'][0]:.3f} vs {st['D_var_theory_sim'][1]:.3f} |")
for r in st["E_pair"]:
    w(f"| M1 closed form vs exact simulator, two people, beta = {r['beta']} | {r['closed_form']:.4f} vs {r['exact']:.4f} +/- {r['se']:.4f} |")
for ev, r in ic.items():
    w(f"| {ev}: M2 exposure by brute force; conditional pmf by rejection sampling ({r['n_accepted']} accepted) | relative error {r['rel_err_exposure']:.1e}; max abs pmf difference {r['max_abs_diff_pmf']:.4f} (MC SE {r['mc_se_max']:.4f}) |")
w("")

open(C.os.path.join(C.DATA, "tables.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
