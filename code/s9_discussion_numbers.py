#!/usr/bin/env python3
"""s9_discussion_numbers.py -- collect every number quoted in the discussion, sections/m8_discussion.tex (Section 8)
(plus a few neighbouring values of the same files, kept so that the context of each quoted number is visible).

The discussion section introduces no new computation.  This script only READS the results files of the
research directories and of the article-level scripts and writes the values quoted in the section, with the path and key of
each, to  data/s9_discussion_numbers.json , so that every figure in the text can be traced and
re-checked with one command:

    python code/s9_discussion_numbers.py

Paths are relative to the repository root.  Nothing is simulated, fitted or rounded here; rounding
happens only in the text.  No figure is produced (the section has none), so there is no _en/_zh variant.
"""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # repository root
OUT = ROOT / "data/s9_discussion_numbers.json"


def load(rel):
    with open(ROOT / rel) as fh:
        return json.load(fh)


out = {}


def put(name, value, src, key):
    out[name] = {"value": value, "src": src, "key": key}


# ---------------------------------------------------------------- mean-field theory, office scene
src = "data/plan_theory_checks.json"
office = load(src)["office"]
put("office_R0_lower_bound", office["lower_bound"], src, "office.lower_bound")
put("office_R0_upper_bound", office["upper_bound"], src, "office.upper_bound")
put("office_R0_D0_1", office["D0=1"]["R0"], src, "office.D0=1.R0")
for row in office["mobility_table"]:
    if row["D0"] in (1.0, 10.0, 100.0, 2400.0, 48700.0):
        put(f"office_excess_pct_D0_{row['D0']:g}", row["excess_pct"], src,
            f"office.mobility_table[D0={row['D0']:g}].excess_pct")
put("office_meeting_room_exit_time_days", office["zone_bound_meeting_room"]["exit_time_days"], src,
    "office.zone_bound_meeting_room.exit_time_days")
for k, v in office["corridor_multiplier_symmetric"].items():
    put(f"office_R0_corridor_{k}", v, src, f"office.corridor_multiplier_symmetric.{k}")
for k in ("walls_10pct", "walls_100pct"):
    put(f"office_partitions_{k}_change_pct", office["partitions"][k]["change_pct"], src,
        f"office.partitions.{k}.change_pct")

src = "data/bsc_theory/worked_examples.json"
phys = load(src)["E9_physical_D"]
for k, v in phys.items():
    if k.startswith("walking") or k.startswith("continuous"):
        put("physical_D0_" + ("walk5pct" if k.startswith("walking") else "continuous"),
            v["D0_lattice2_per_day"], src, f"E9_physical_D[{k}].D0_lattice2_per_day")

# ---------------------------------------------------------------- scenes (exact simulation)
src = "data/bsc_sim/03_scenes.json"
sc = load(src)
r = sc["office|D0=1|range1m"]
put("office_counted_R", r["R_index_alone_mean"], src, "office|D0=1|range1m.R_index_alone_mean")
put("office_counted_R_se", r["R_index_alone_se"], src, "office|D0=1|range1m.R_index_alone_se")
put("office_R1bar", r["R1_pair_uniform_index"], src, "office|D0=1|range1m.R1_pair_uniform_index")
put("office_p_major", r["p_major"], src, "office|D0=1|range1m.p_major")
put("office_p_branching", r["p_major_branching"], src, "office|D0=1|range1m.p_major_branching")
for key in ("metro|D0=1|range1m", "metro|D0=1|samecell", "classroom|D0=1|range1m",
            "classroom|D0=1|samecell"):
    put("p_major_" + key, sc[key]["p_major"], src, key + ".p_major")
    put("n_rep_" + key, sc[key]["n_rep"], src, key + ".n_rep")
put("metro_samecell_R0", sc["metro|D0=1|samecell"]["R0_ngm"], src, "metro|D0=1|samecell.R0_ngm")
put("metro_samecell_R1bar", sc["metro|D0=1|samecell"]["R1_pair_uniform_index"], src,
    "metro|D0=1|samecell.R1_pair_uniform_index")
put("metro_samecell_rhoK1", sc["metro|D0=1|samecell"]["R1_pair_ngm"], src,
    "metro|D0=1|samecell.R1_pair_ngm")
src = "code/bsc_sim/verify/v05_scenes.json"
v05 = load(src)["metro|1.0|same"]
put("metro_samecell_check_runs", v05["n"], src, "metro|1.0|same.n")
put("metro_samecell_check_p_major", v05["p_major"], src, "metro|1.0|same.p_major")

# ---------------------------------------------------------------- dwell time / duty cycle
src = "data/bsc_sim/11_dwell_time.json"
dw = load(src)
for d0 in ("1", "10"):
    r = dw[f"office|D0={d0}"]
    for k in ("R0_continuous", "R0_duty_cycle", "wellmixed_duty_cycle", "R1_pair_duty_cycle",
              "sim_R_schedule", "sim_R_schedule_lo", "sim_R_schedule_hi", "p_major", "per_visit_meanfield",
              "per_visit_sim"):
        put(f"office_8h_D0_{d0}_{k}", r[k], src, f"office|D0={d0}.{k}")
for scn in ("supermarket", "classroom", "metro"):
    r = dw[f"{scn}|D0=1"]
    put(f"{scn}_per_visit_meanfield", r["per_visit_meanfield"], src, f"{scn}|D0=1.per_visit_meanfield")
    put(f"{scn}_visit_hours", r["visit_hours"], src, f"{scn}|D0=1.visit_hours")
# pooled duty-cycle count of the re-check (numbers of the re-check of the simulations, data/checks/)
src = "data/checks/simulation.json"
pooled = load(src)["values"]["duty_cycle_office_D0_1"]["count_explicit_schedule_pooled"]
put("office_8h_D0_1_counted_check_pooled", [pooled["mean"], pooled["se"]], src,
    "values.duty_cycle_office_D0_1.count_explicit_schedule_pooled (mean, se)")

# ---------------------------------------------------------------- interventions
src = "data/bsc_sim/08_interventions.json"
iv = load(src)
for d0 in ("1", "10"):
    base = iv[f"office|D0={d0}|baseline|FD"]
    part = iv[f"office|D0={d0}|partitions (+39 barrier cells)|FD"]
    capF = iv[f"office|D0={d0}|capacity 50 %|FD"]
    capD = iv[f"office|D0={d0}|capacity 50 %|DD"]
    vent = iv[f"office|D0={d0}|ventilation q x 0.5|FD"]
    mask = iv[f"office|D0={d0}|masks beta x 0.3|FD"]
    put(f"iv_D0_{d0}_baseline_attack_all", base["attack_all_mean"], src, f"office|D0={d0}|baseline|FD.attack_all_mean")
    put(f"iv_D0_{d0}_partitions_attack_all", part["attack_all_mean"], src,
        f"office|D0={d0}|partitions (+39 barrier cells)|FD.attack_all_mean")
    put(f"iv_D0_{d0}_partitions_attack_change_pct",
        100 * (part["attack_all_mean"] / base["attack_all_mean"] - 1), src, "derived: ratio of the two above")
    put(f"iv_D0_{d0}_baseline_R1bar", base["R1_pair_uniform"], src, f"office|D0={d0}|baseline|FD.R1_pair_uniform")
    put(f"iv_D0_{d0}_capacity_FD_R1bar", capF["R1_pair_uniform"], src, f"office|D0={d0}|capacity 50 %|FD.R1_pair_uniform")
    put(f"iv_D0_{d0}_capacity_DD_R1bar", capD["R1_pair_uniform"], src, f"office|D0={d0}|capacity 50 %|DD.R1_pair_uniform")
    put(f"iv_D0_{d0}_capacity_FD_R0", capF["R0_ngm"], src, f"office|D0={d0}|capacity 50 %|FD.R0_ngm")
    put(f"iv_D0_{d0}_capacity_DD_R0", capD["R0_ngm"], src, f"office|D0={d0}|capacity 50 %|DD.R0_ngm")
    put(f"iv_D0_{d0}_ventilation_R0", vent["R0_ngm"], src, f"office|D0={d0}|ventilation q x 0.5|FD.R0_ngm")
    put(f"iv_D0_{d0}_masks_R0", mask["R0_ngm"], src, f"office|D0={d0}|masks beta x 0.3|FD.R0_ngm")

src = "data/bsc_sim/07_occupancy.json"
occ = load(src)
for n in (50, 100, 400):
    r = occ[f"D0=1|N={n}|FD"]
    put(f"occupancy_D0_1_N_{n}_FD_R1bar", r["R1_pair_uniform"], src, f"D0=1|N={n}|FD.R1_pair_uniform")
    put(f"occupancy_D0_1_N_{n}_FD_R0", r["R0_ngm"], src, f"D0=1|N={n}|FD.R0_ngm")

src = "code/bsc_sim/verify/v08b_obstacles_many.json"
ob = load(src)
put("barriers_attack_0", ob["0.0"]["ar"], src, "0.0.ar [mean, SE over layouts]")
put("barriers_attack_10pct", ob["0.1"]["ar"], src, "0.1.ar")
put("barriers_diff_10pct_attack", ob["diff_0.1_ar"], src, "diff_0.1_ar")
put("barriers_diff_10pct_peak_day", ob["diff_0.1_pkt"], src, "diff_0.1_pkt")

src = "data/bsc_sim/05_heterogeneity.json"
het = load(src)
for d0 in ("1", "10"):
    a = het[f"D0={d0}|N=100|c=0.0000"]
    b = het[f"D0={d0}|N=100|c=0.5833"]
    for k in ("R0_ngm", "sim_R_perron", "p_major", "attack_all"):
        put(f"het_D0_{d0}_{k}", [a[k], b[k]], src, f"D0={d0}|N=100|c=0.0000 and c=0.5833: {k}")

# ---------------------------------------------------------------- outbreak validation
src = "data/bsc_validation/03_holdout_shape.json"
ho = load(src)["events"]
for ev in ("T1", "T2", "T5", "C1"):
    put(f"holdout_{ev}_M2_p", ho[ev]["M2"]["p_two_sided"], src, f"events.{ev}.M2.p_two_sided")
    put(f"holdout_{ev}_M0_p", ho[ev]["M0"]["p_two_sided"], src, f"events.{ev}.M0.p_two_sided")
    put(f"holdout_{ev}_M2_rr_pred", ho[ev]["M2"]["rr_pred"], src, f"events.{ev}.M2.rr_pred")
    put(f"holdout_{ev}_rr_obs", ho[ev]["M2"]["rr_obs"], src, f"events.{ev}.M2.rr_obs")

src = "data/validation_checks/callcentre_quadrature.json"
cq = load(src)
put("callcentre_quadrature", {k: v for k, v in cq.items() if not isinstance(v, list) and k != "description"},
    src, "(scalar fields)")

vv = ROOT / "data/validation_checks"
t = (vv / "v11_breakdown.txt").read_text()
m = re.search(r"at posterior mode D=([0-9.]+): pmf\(obs\)=([0-9.]+), two-sided p=([0-9.]+)", t)
put("classroom_at_restaurant_mode", {"D_m2_per_h": float(m.group(1)), "p": float(m.group(3))},
    "data/validation_checks/v11_breakdown.txt", "line 'at posterior mode'")
m = re.search(r"train-calibrated D_air \(exploratory, post hoc\): pmf\(obs\)=([0-9.]+) p=([0-9.]+) E\[k_near\]=([0-9.]+)", t)
put("classroom_with_train_coefficient", {"p": float(m.group(2)), "expected_near_cases": float(m.group(3))},
    "data/validation_checks/v11_breakdown.txt", "line 'C1 with the train-calibrated D_air'")
m = re.search(r"T2 at train posterior mode D=20: pmf\(obs\)=([0-9.]+) p=([0-9.]+)", t)
put("coach_at_train_mode_p", float(m.group(2)),
    "data/validation_checks/v11_breakdown.txt", "line 'T2 at train posterior mode'")

t = (vv / "v10_c1_loo.txt").read_text()
pref = dict(re.findall(r"^(T[125]) D preferred alone ([0-9.eE+-]+)", t, flags=re.M))
put("preferred_D_air_alone_m2_per_h", {k: float(v) for k, v in pref.items()},
    "data/validation_checks/v10_c1_loo.txt", "lines '<event> D preferred alone'")
m = re.search(r"T4 preferred ([0-9.]+)", t)
put("train_preferred_D_air", float(m.group(1)), "data/validation_checks/v10_c1_loo.txt", "T4 preferred")
m = re.search(r"LR ([0-9.]+) p\(3 df\) ([0-9.]+)", t)
put("common_D_test", {"LR": float(m.group(1)), "df": 3, "p": float(m.group(2))},
    "data/validation_checks/v10_c1_loo.txt", "line 'sum of maxima ... LR'")

t = (vv / "v05_generic.txt").read_text()
gen = re.findall(r"constant RR=([0-9.]+): pooled delta ([+-][0-9.]+).*failures (\[.*\])", t)
put("generic_constant_RR_model", [{"RR": float(a), "pooled_delta": float(b), "failures": c} for a, b, c in gen],
    "data/validation_checks/v05_generic.txt", "lines 'constant RR='")
m = re.search(r"M2 \(re-check re-implementation\): pooled \(np\.float64\(([0-9.]+)\)", t)
put("M2_pooled_delta_check", float(m.group(1)), "data/validation_checks/v05_generic.txt", "first line")

t = (vv / "v12_pooled_variants.txt").read_text()
m = re.search(r"without T5: \(([-0-9.]+),", t)
put("pooled_delta_without_flight", float(m.group(1)), "data/validation_checks/v12_pooled_variants.txt",
    "line 'without T5'")

# M1 reach ceiling and implied emission: generated tables of the validation part
src = "data/bsc_validation/tables.md"
t = (ROOT / src).read_text()
m = re.search(r"\| T1 Zhejiang bus \| 23/67 \| 0 of 24 \| ([0-9.]+) \(([0-9.]+)-([0-9.]+)\)", t)
put("M1_ceiling_bus_mean_total", [float(m.group(1)), float(m.group(2)), float(m.group(3))], src, "Table V8, row T1")
m = re.search(r"\| H1 Skagit choir \| 52/60 .*\| ([0-9]+) \(", t)
put("emission_choir_quanta_per_h", float(m.group(1)), src, "Table V10, row H1")
m = re.search(r"\| T2 Hunan coach \| 7/46 .*\| ([0-9]+) \(", t)
put("emission_coach_quanta_per_h", float(m.group(1)), src, "Table V10, row T2")

# restaurant: tracer against fitted exposure (Table V20)
rows = re.findall(r"^\| (T[A-Z0-9]+) \| ([0-9.]+) \| ([0-9.]+) \|$", (ROOT / "data/bsc_validation/tables.md").read_text(),
                  flags=re.M)
put("restaurant_tracer_vs_M2", [{"table": a, "tracer_rel": float(b), "M2_exposure_rel": float(c)} for a, b, c in rows],
    "data/bsc_validation/tables.md", "Table V20")

# dataset size
src = "data/bsc_outbreaks/outbreaks.json"
ob = load(src)
recs = ob if isinstance(ob, list) else ob.get("outbreaks", ob.get("records", ob))
put("dataset_n_records", len(recs), src, "number of records")

# ---------------------------------------------------------------- pair-level and schedule values
src = "data/s5_pair_large_sample.json"
ls = load(src)
for d0 in ("1", "10", "100"):
    r = ls[f"office|D0={d0}"]
    put(f"office_D0_{d0}_counted_large_sample", [r["mean"], r["se"], r["n_runs"]], src, f"office|D0={d0}.mean, se, n_runs")
    put(f"office_D0_{d0}_R1bar_exact", r["R1bar_exact"], src, f"office|D0={d0}.R1bar_exact")
src = "data/s7_interventions_schedule.json"
sch = load(src)
for d0 in ("1", "10"):
    r = sch[f"duty_exact|office|D0={d0}"]
    put(f"office_8h_D0_{d0}_exact_start", r["exact_start"], src, f"duty_exact|office|D0={d0}.exact_start")
    put(f"office_8h_D0_{d0}_exact_uniform_phase", r["exact_uniform_phase"], src,
        f"duty_exact|office|D0={d0}.exact_uniform_phase")
src = "data/s2_model_scenes.json"
s2 = load(src)
put("office_meeting_room_residence_time_days", s2["office_regime"]["meeting_room_1_over_mu_Z_days"], src,
    "office_regime.meeting_room_1_over_mu_Z_days")
put("office_meeting_room_mean_exit_time_uniform_start_days",
    s2["office_regime"]["meeting_room_mean_exit_time_uniform_start_days"], src,
    "office_regime.meeting_room_mean_exit_time_uniform_start_days")
put("walking_mobility_range_cell2_per_day",
    [min(v["walking_5pct_D0_cell2_per_day"] for v in s2["scenes"].values()),
     max(v["walking_D0_cell2_per_day"] for v in s2["scenes"].values())], src,
    "scenes.*.walking_5pct_D0_cell2_per_day (min), walking_D0_cell2_per_day (max)")
src = "data/s4_geometry_barriers_leak.json"
bl = load(src)["A_partitions"]
put("office_barrier_cells_change_pct", [bl["barrier_cells_39_workstations"][f"D0={d}"]["change_pct"] for d in ("1", "10", "100")],
    src, "A_partitions.barrier_cells_39_workstations.D0=*.change_pct")
put("office_walls_change_pct_D0_1", [bl["walls"][k]["D0=1"]["change_pct"] for k in bl["walls"]], src,
    "A_partitions.walls.*.D0=1.change_pct")
src = "data/s5_pair_checks.json"
fs = load(src)["finite_size"]
put("closed_form_max_rel_deviation_reflecting_walls", fs["max_abs_rel_dev_closed_vs_exact"], src,
    "finite_size.max_abs_rel_dev_closed_vs_exact")

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fh:
    json.dump(out, fh, indent=1)
print(f"wrote {OUT.relative_to(ROOT)} ({len(out)} entries)")
for k, v in out.items():
    print(f"{k:48s} {json.dumps(v['value'])[:110]}")
