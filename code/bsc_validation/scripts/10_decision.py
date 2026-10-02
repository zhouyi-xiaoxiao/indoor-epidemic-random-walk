"""Apply the pre-registered decision rule (section 8) to the stored results."""
import numpy as np

import _common as C

sh = C.load_json("03_holdout_shape.json")
o1 = C.load_json("06_callcentre.json")
lv = C.load_json("04_level_controls.json")
ex = C.load_json("02b_m1_exact.json")["events"]
pw = C.load_json("02_power.json")
EV = ["T1", "T2", "T5", "C1"]
out = {"stamp": C.stamp(), "models": {}}
for m in ("M1", "M2"):
    r = {}
    r["A1"] = dict(delta=sh["pooled"][m]["delta"], p=sh["pooled"][m]["p_under_M0"], passed=sh["pooled"][m]["A1_pass"],
                   power_if_true=pw[m][m]["P_A1"])
    fails, labels = [], {}
    for ev in EV:
        x = sh["events"][ev][m]
        lab = x["label"]
        bad = x["p_two_sided"] < 0.05
        if m == "M1" and sh["events"][ev]["M1_uses_exact_simulator"]:
            pge = float(np.mean([q["p_total_ge_K"] for q in ex[ev]["rows"]]))
            if pge < 0.01:
                bad, lab = True, "not reproduced (total unattainable)"
        labels[ev] = dict(p=x["p_two_sided"], label=lab)
        if bad:
            fails.append(ev)
    r["A2_single_events"] = labels
    r["O1"] = dict(P_z_le_obs_pooled=o1[m]["P_z_le_obs"], inside_95_pooled=o1[m]["inside_95"],
                   P_z_le_obs_equal_weight=o1[m]["P_z_le_obs_equal_weight"],
                   inside_95_equal_weight=o1[m]["inside_95_equal_weight"])
    r["A2_failures_as_run"] = fails + ([] if o1[m]["inside_95"] else ["O1"])
    r["A2_failures_conservative"] = fails + ([] if (o1[m]["inside_95"] and o1[m]["inside_95_equal_weight"]) else ["O1"])
    tied = lv["tied_T2_T3"][m]["inside"]
    c3 = lv["C3"][f"{m}_masks_0.35"]["inside"]
    r["A3"] = dict(tied_index=tied, unselected_transfer=c3, passed=bool(tied and c3))
    pe = lv["T5_PE"]["M2"]["P_exceed"] if m == "M2" else lv["T5_PE"]["M1_exact_beta_cap"]["P_exceed"]
    ctrl = dict(S1=lv["S1"]["M2_single_index" if m == "M2" else "M1_pooled"]["P_exceed"],
                C3=lv["C3"][f"{m}_masks_0.35"]["P_exceed_bound"], T5_PE=pe, T7=lv["T7"][m]["P_ge5"])
    r["A4"] = dict(P_exceed=ctrl, passed=bool(all(v <= 0.5 for v in ctrl.values())))

    def verdict(nf):
        if not r["A4"]["passed"]:
            return "W5"
        if nf >= 2:
            return "W4"
        if not r["A1"]["passed"]:
            return "W3"
        if nf == 0 and r["A3"]["passed"] and r["A1"]["power_if_true"] >= 0.8:
            return "W1"
        return "W2"
    r["outcome_as_run"] = verdict(len(r["A2_failures_as_run"]))
    r["outcome_conservative"] = verdict(len(r["A2_failures_conservative"]))
    out["models"][m] = r
    print(m, "A1", r["A1"], "\n   A2 fails (as run)", r["A2_failures_as_run"], "| conservative", r["A2_failures_conservative"],
          "\n   A3", r["A3"], "A4", r["A4"]["passed"], "->", r["outcome_as_run"], "/", r["outcome_conservative"])
C.save_json("10_decision.json", out)
