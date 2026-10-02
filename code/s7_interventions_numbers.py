#!/usr/bin/env python
"""s7_interventions_numbers.py -- every number quoted on interventions (sections/m6_simulation.tex, supp_sim.tex),
collected from the results files, plus the LaTeX rows of its tables.

No simulation is run here.  Sources (paths relative to the repository root):
  data/bsc_sim/05_heterogeneity.json, 05_heterogeneity_extra.json   (two-zone room)
  code/bsc_sim/verify/v09_het.json                                      (re-check, two-zone room)
  data/bsc_sim/06_obstacles.json                                   (barrier scan, 8 layouts x 500)
  code/bsc_sim/verify/v08b_obstacles_many.json                          (re-check, 40 layouts x 1000)
  data/s7_interventions_barriers.json                           (this article, D = 1 and 10)
  data/bsc_sim/07_occupancy.json, 08_interventions.json
  data/bsc_sim/09_metro_time.json, 11_dwell_time.json
  code/bsc_sim/verify/v11.json, v14_schedule.json                       (re-check)
  data/s7_interventions_schedule.json                           (this article)
  data/s7_interventions_hetero_se.json                          (this article; standard errors)
  data/s7_interventions_regime.json                             (this article)
  data/bsc_sim/10_benchmark.json, code/bsc_sim/verify/v12_benchmark_rerun.json, v13_bench.json

Output: ../data/s7_interventions_numbers.json and ../data/s7_interventions_tables.tex (the table rows).

    python s7_interventions_numbers.py            # write both files
    python s7_interventions_numbers.py --check    # also verify that every generated table row occurs verbatim
                                                  # in ../sections/m6_simulation.tex or supp_sim.tex
                                                  # (tables of Section 6 and of Supplementary Section S5)
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
SIM = ROOT / "code/bsc_sim"
ART = HERE.parent / "data"
TEXFILES = [HERE.parent / "sections" / "m6_simulation.tex", HERE.parent / "sections" / "supp_sim.tex"]


def load(p):
    return json.loads(Path(p).read_text())


HET = load(SIM / "data/05_heterogeneity.json")
HETX = load(SIM / "data/05_heterogeneity_extra.json")
VHET = load(SIM / "verify/v09_het.json")
OBS = load(SIM / "data/06_obstacles.json")
VOBS = load(SIM / "verify/v08b_obstacles_many.json")
BAR = load(ART / "s7_interventions_barriers.json")
OCC = load(SIM / "data/07_occupancy.json")
INT = load(SIM / "data/08_interventions.json")
MET = load(SIM / "data/09_metro_time.json")
DWL = load(SIM / "data/11_dwell_time.json")
V11 = load(SIM / "verify/v11.json")
V14 = load(SIM / "verify/v14_schedule.json")
SCH = load(ART / "s7_interventions_schedule.json")
HSE = load(ART / "s7_interventions_hetero_se.json")      # standard errors; 2,000-run samples for N = 1000
REG = load(ART / "s7_interventions_regime.json")
BEN = load(SIM / "data/10_benchmark.json")
V12 = load(SIM / "verify/v12_benchmark_rerun.json")
V13 = load(SIM / "verify/v13_bench.json")

NUM: dict = {}
TEX: list[str] = []          # table rows (checked against the section)
NOTE: list[str] = []         # comments


def f(x, n=2):
    return f"{x:.{n}f}"


def pm(m, s, n=3):
    return f"${m:.{n}f}\\pm{s:.{n}f}$"


def spm(m, s, n=3):
    sign = "+" if round(m, n) >= 0 else "-"
    return f"${sign}{abs(m):.{n}f}\\pm{s:.{n}f}$"


def pct(x, n=1):
    sign = "+" if x >= 0 else "-"
    return f"${sign}{abs(100 * x):.{n}f}$"


# ------------------------------------------------------------------------------------------------ 1
def hetero():
    out = {}
    NOTE.append("% ---- Table s7_interventions:tab:hetero")
    c58 = "0.5833"

    def row(D0, N, mob, v, vv=None, first=True, se=None):
        head = [str(D0), str(N), mob] if first else ["", "", ""]
        cval = "0.58" if abs(v["c"] - 0.5833) < 1e-3 else (f"${v['c']:.2f}$" if v["c"] < 0 else f(v["c"], 2))
        cells = head + [cval, f(v["R0_ngm"], 2), f(v["R1_pair_ngm"], 2),
                        f"{v['sim_R_perron']:.2f} ({v['sim_R_perron_lo']:.2f}--{v['sim_R_perron_hi']:.2f})"]
        if "sim_R_uniform" in v:
            cells.append(f"{v['sim_R_uniform']:.2f} ({v['sim_R_uniform_lo']:.2f}--{v['sim_R_uniform_hi']:.2f})")
        elif vv is not None:
            cells.append(f"{vv['sim_uniform']:.2f}$^{{\\dagger}}$")
        else:
            cells.append("--")
        # outbreak probability and mean attack rate with standard errors: the stored
        # 2,000-run batches for N = 100, and 2,000 further epidemics (other seed) for N = 1000, where the
        # stored batches have 300 epidemics only (s7_interventions_hetero_se.py)
        cells += [pm(se["p_major"], se["p_major_se"]), pm(se["attack_all"], se["attack_all_se"])]
        TEX.append(" & ".join(cells) + r" \\")

    plan = {1: [(100, "zones", ("0.0000", c58, "1.0000")), (1000, "zones", ("0.0000", c58, "1.0000")),
                (100, "uniform", ("0.0000", c58)), (100, "zones", ("-2.0000",))],
            10: [(100, "zones", ("0.0000", c58, "1.0000")), (1000, "zones", ("0.0000", c58)),
                 (100, "uniform", (c58,)), (100, "zones", ("-2.0000",))]}
    for D0, blocks in plan.items():
        for N, mob, cs in blocks:
            for i, c in enumerate(cs):
                if mob == "uniform":
                    v = HETX[f"extra|uniformD|D0={D0}|N=100|c={c}"]
                    vv = VHET.get(f"uniformD|D0={D0}.0|N=100|c={float(c):.3f}")
                elif float(c) < 0:
                    v = HETX[f"extra|zoneD|D0={D0}|N=100|c={c}"]
                    vv = VHET.get(f"zoneD|D0={D0}.0|N=100|c={float(c):.3f}")
                else:
                    v, vv = HET[f"D0={D0}|N={N}|c={c}"], None
                out[f"D0={D0}|N={N}|{mob}|c={c}"] = {k: v[k] for k in v if k != "seconds"}
                if N == 1000:
                    se = HSE[f"new|D0={D0}|N=1000|c={c}"]
                elif mob == "uniform":
                    se = HSE[f"stored|extra|uniformD|D0={D0}|N=100|c={c}"]
                elif float(c) < 0:
                    se = HSE[f"stored|extra|zoneD|D0={D0}|N=100|c={c}"]
                else:
                    se = HSE[f"stored|D0={D0}|N=100|c={c}"]
                if N != 1000:
                    assert abs(se["p_major"] - v["p_major"]) < 1e-12 and abs(se["attack_all"] - v["attack_all"]) < 1e-12
                out[f"D0={D0}|N={N}|{mob}|c={c}"]["epidemics_se"] = {k: se[k] for k in (
                    "n", "p_major", "p_major_se", "attack_all", "attack_all_se")}
                row(D0, N, mob, v, vv, first=(i == 0), se=se)
    a, b = HET["D0=1|N=100|c=0.0000"], HET["D0=1|N=100|c=0.5833"]
    out["changes_D1_N100"] = dict(
        R0_rel=b["R0_ngm"] / a["R0_ngm"] - 1, perron_rel=b["sim_R_perron"] / a["sim_R_perron"] - 1,
        uniform_rel=b["sim_R_uniform"] / a["sim_R_uniform"] - 1,
        gen_ratio_rel=b["sim_gen2_over_gen1"] / a["sim_gen2_over_gen1"] - 1,
        attack_rel=b["attack_all"] / a["attack_all"] - 1, R0_over_beta_gamma=b["R0_ngm"] / (0.5 / 0.14))
    a, b = HET["D0=10|N=100|c=0.0000"], HET["D0=10|N=100|c=0.5833"]
    out["changes_D10_N100"] = dict(
        R0_rel=b["R0_ngm"] / a["R0_ngm"] - 1, perron_rel=b["sim_R_perron"] / a["sim_R_perron"] - 1,
        uniform_rel=b["sim_R_uniform"] / a["sim_R_uniform"] - 1, attack_rel=b["attack_all"] / a["attack_all"] - 1,
        p_major_rel=b["p_major"] / a["p_major"] - 1)
    a, b = HET["D0=1|N=1000|c=0.0000"], HET["D0=1|N=1000|c=0.5833"]
    out["changes_D1_N1000"] = dict(attack_rel=b["attack_all"] / a["attack_all"] - 1)
    out["check"] = VHET
    NUM["hetero"] = out


# ------------------------------------------------------------------------------------------------ 2
def barriers():
    out = {"first_8x500": {}, "check_40x1000": VOBS, "article_40_layouts": {}}
    for cal in ("FD", "DD"):
        for fr in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30):
            rows = [OBS[f"layout{l}|rhoB={fr:.2f}|{cal}"] for l in range(8)]
            d = {}
            for k in ("rhoB_realised", "R0_ngm", "R1_pair_centre", "R1_pair_uniform", "R_index_alone", "p_major",
                      "attack_all_mean", "attack_major_mean", "peak_prev_major_mean", "peak_time_major_mean"):
                v = np.array([r[k] for r in rows])
                # mean, standard error and 95 % half-width (t quantile, 7 d.f.) across the 8 layouts
                d[k] = [float(v.mean()), float(v.std(ddof=1) / np.sqrt(8)),
                        float(2.365 * v.std(ddof=1) / np.sqrt(8))]
            out["first_8x500"][f"rhoB={fr:.2f}|{cal}"] = d
    NOTE.append("% ---- Table s7_interventions:tab:barriers")
    for Dk in ("D=1", "D=10"):
        b = BAR[Dk]
        out["article_40_layouts"][Dk] = dict(n_layout=b["n_layout"], n_rep=b["n_rep"], arms=b["arms"],
                                             paired=b["paired"])
        for cal in ("FD", "DD"):
            for fr in ("0.00", "0.10", "0.30"):
                if cal == "DD" and fr == "0.00":
                    continue                       # identical to the FD row
                a = b["arms"][f"rhoB={fr}|{cal}"]
                cells = [Dk.split("=")[1], cal, fr, f(a["R0"][0], 2), f(a["R1_centre"][0], 2),
                         f(a["p_major"][0], 3), f(a["attack"][0], 3), f(100 * a["peak_major"][0], 1),
                         f(a["peak_day_major"][0], 1)]
                if fr == "0.00":
                    cells += ["--", "--", "--"]
                else:
                    p = b["paired"][f"rhoB={fr}|{cal}"]
                    cells += [spm(*p["attack"]), spm(100 * p["peak_major"][0], 100 * p["peak_major"][1], 1),
                              spm(*p["peak_day_major"], n=2)]
                TEX.append(" & ".join(cells) + r" \\")
        if Dk == "D=1":
            v = VOBS
            for fr in ("0.1", "0.3"):
                TEX.append(" & ".join(["1", "FD$^{\\ddagger}$", f"{float(fr):.2f}", "4.29", "--",
                                       f(v[fr]["pm"][0], 3), f(v[fr]["ar"][0], 3), f(100 * v[fr]["pk"][0], 1),
                                       f(v[fr]["pkt"][0], 1), spm(*v[f"diff_{fr}_ar"]),
                                       spm(100 * v[f"diff_{fr}_pk"][0], 100 * v[f"diff_{fr}_pk"][1], 1),
                                       spm(*v[f"diff_{fr}_pkt"], n=2)]) + r" \\")
    NUM["barriers"] = out


# ------------------------------------------------------------------------------------------------ 3
def occupancy():
    out = {}
    for D0 in (1, 10):
        for N in (20, 50, 75, 100, 150, 200, 400):
            for cal in ("FD", "DD"):
                v = OCC[f"D0={D0}|N={N}|{cal}"]
                out[f"D0={D0}|N={N}|{cal}"] = {k: v[k] for k in (
                    "R0_ngm", "R1_pair_uniform", "R_index_alone", "R_index_alone_se", "p_major", "p_major_lo",
                    "p_major_hi", "p_major_branching", "attack_all_mean", "attack_major_mean",
                    "peak_prev_major_mean", "peak_time_major_mean")}
        for cal in ("FD", "DD"):
            lo, hi = OCC[f"D0={D0}|N=75|{cal}"], OCC[f"D0={D0}|N=150|{cal}"]
            out[f"elasticity|D0={D0}|{cal}"] = dict(
                pair=float(np.log(hi["R1_pair_uniform"] / lo["R1_pair_uniform"]) / np.log(2)),
                counted=float(np.log(hi["R_index_alone"] / lo["R_index_alone"]) / np.log(2)),
                meanfield=float(np.log(hi["R0_ngm"] / lo["R0_ngm"]) / np.log(2)))
    NUM["occupancy"] = out


# ------------------------------------------------------------------------------------------------ 4
def interventions():
    out = {}
    NOTE.append("% ---- Table s7_interventions:tab:interventions")
    arms = [("baseline", "FD", "baseline"), ("ventilation q x 0.5", "FD", r"ventilation, $q\to q/2$"),
            ("masks beta x 0.3", "FD", r"masks, $\beta\to0.3\beta$"),
            ("capacity 50 %", "FD", r"capacity 50\,\%"), ("capacity 50 %", "DD", r"capacity 50\,\%"),
            ("partitions (+39 barrier cells)", "FD", "partitions"),
            ("partitions (+39 barrier cells)", "DD", "partitions"),
            ("combination", "FD", "combination"), ("combination", "DD", "combination")]
    for D0 in (1, 10):
        for arm, cal, lab in arms:
            v = INT[f"office|D0={D0}|{arm}|{cal}"]
            out[f"office|D0={D0}|{arm}|{cal}"] = {k: v.get(k) for k in (
                "N", "R0_ngm", "R1_pair_uniform", "R_index_alone", "R_index_alone_se", "p_major", "p_major_lo",
                "p_major_hi", "attack_all_mean", "attack_all_sd", "attack_major_mean", "peak_prev_major_mean",
                "peak_time_major_mean", "n_rep")}
            TEX.append(" & ".join([str(D0), lab, cal, str(v["N"]), f(v["R0_ngm"], 2),
                                   f(v["R1_pair_uniform"], 2),
                                   pm(v["R_index_alone"], v["R_index_alone_se"], 2), f(v["p_major"], 3),
                                   f(v["attack_all_mean"], 3), f(100 * v["peak_prev_major_mean"], 1),
                                   f(v["peak_time_major_mean"], 1)]) + r" \\")
        for arm in ("hotspot ventilation (25 % of cells)", "coldspot ventilation (25 % of cells)",
                    "pair-level hotspot ventilation (25 % of cells)"):
            v = INT[f"office|D0={D0}|{arm}|FD"]
            out[f"office|D0={D0}|{arm}|FD"] = {k: v.get(k) for k in (
                "R0_ngm", "R1_pair_uniform", "R_index_alone", "p_major", "attack_all_mean", "attack_all_sd",
                "n_rep")}
            out[f"office|D0={D0}|{arm}|FD"]["attack_all_se"] = v["attack_all_sd"] / np.sqrt(v["n_rep"])
    for D0 in (1, 10):
        for arm, cal in (("baseline", "FD"), ("combination", "FD"), ("combination", "DD")):
            v = INT[f"classroom|D0={D0}|{arm}|{cal}"]
            out[f"classroom|D0={D0}|{arm}|{cal}"] = {k: v.get(k) for k in (
                "N", "R0_ngm", "R1_pair_uniform", "R_index_alone", "R_index_alone_se", "p_major",
                "attack_all_mean")}
    NUM["interventions"] = out


# ------------------------------------------------------------------------------------------------ 5
def dwell():
    out = {"first": {}, "article": {}, "identity": SCH.get("identity"),
           "check": {k: V11[k] for k in V11 if k.startswith(("office_duty", "visit"))},
           "check_schedule": V14, "metro": {}}
    NOTE.append("% ---- Table s7_interventions:tab:dwell (a) duty cycle")
    names = dict(office="office", classroom="classroom", supermarket="supermarket", metro="metro carriage")
    for scene in ("office", "classroom"):
        for D0 in (1, 10):
            v = DWL[f"{scene}|D0={D0}"]
            out["first"][f"{scene}|D0={D0}"] = v
            cells = [f"{names[scene]}, {v['hours_per_day']:g}", str(D0), f(v["R0_continuous"], 2),
                     f(v["R0_duty_cycle"], 3), f(v["R1_pair_duty_cycle"], 3)]
            key = f"duty|{scene}|D0={D0}"
            if key in SCH:
                s = SCH[key]
                out["article"][key] = {k: s[k] for k in s if not isinstance(s[k], list)}
                out["article"][key]["count_start_over_rule"] = s["count_start"] / s["R1bar_duty"]
                out["article"][key]["count_uniform_over_rule"] = s["count_uniform_phase"] / s["R1bar_duty"]
            # exact pair-level expectations with the on/off schedule replace the counts
            ex = SCH[f"duty_exact|{scene}|D0={D0}"]
            out["article"][f"duty_exact|{scene}|D0={D0}"] = ex
            assert abs(ex["rule_gamma_over_f"] - v["R1_pair_duty_cycle"]) < 1e-6
            cells += [f(ex["exact_start"], 3), f(ex["exact_uniform_phase"], 3)]
            lg = SCH.get(f"duty_large|{scene}|D0={D0}")
            if lg is not None:
                out["article"][f"duty_large|{scene}|D0={D0}"] = {k: lg[k] for k in lg if not isinstance(lg[k], list)}
            se = (v["sim_R_schedule_hi"] - v["sim_R_schedule"]) / 1.96
            out["article"][f"count_single_batch|{scene}|D0={D0}"] = dict(mean=v["sim_R_schedule"], se=se)
            cells.append(f(v["p_major"], 3) if round(v["p_major"], 3) != 0.001 or v["p_major"] >= 0.001
                         else f"{v['p_major']:.4f}")
            TEX.append(" & ".join(cells) + r" \\")
    for scene in ("supermarket", "metro"):
        v = DWL[f"{scene}|D0=1"]
        out["first"][f"{scene}|D0=1"] = v
        out["first"][f"{scene}|D0=10"] = DWL[f"{scene}|D0=10"]
        TEX.append(" & ".join([f"{names[scene]}, {v['hours_per_day']:.2f}", "1", f(v["R0_continuous"], 2),
                               f"\\multicolumn{{2}}{{c}}{{{v['wellmixed_duty_cycle']:.3f}}}", "--", "--", "--"])
                   + r" \\")
    NOTE.append("% ---- Table s7_interventions:tab:dwell (b) one visit")
    vis = dict(office="office day", classroom="lecture", metro="metro ride", supermarket="supermarket visit")
    for scene in ("office", "classroom", "metro", "supermarket"):
        a, b = DWL[f"{scene}|D0=1"], DWL[f"{scene}|D0=10"]
        se1 = (a["per_visit_sim_hi"] - a["per_visit_sim"]) / 1.96
        se10 = (b["per_visit_sim_hi"] - b["per_visit_sim"]) / 1.96
        TEX.append(" & ".join([vis[scene], f"{a['visit_hours']:.2f}", f(a["per_visit_meanfield"], 4),
                               pm(a["per_visit_sim"], se1, 4), pm(b["per_visit_sim"], se10, 4)]) + r" \\")
        out["first"][f"{scene}|exposure_rel_diff"] = [a["per_visit_exposure_sim"] / a["per_visit_meanfield"] - 1,
                                                      b["per_visit_exposure_sim"] / b["per_visit_meanfield"] - 1]
    NOTE.append("% ---- Table s7_interventions:tab:schedule (metro)")
    arms = dict(peak="peak", varying="schedule", mean="average")
    for D0 in (1, 10, 100):
        for arm in ("peak", "varying", "mean"):
            v = MET.get(f"D0={D0}|{arm}") or SCH.get(f"metro|D0={D0}|{arm}")
            out["metro"][f"D0={D0}|{arm}"] = {k: v.get(k) for k in (
                "R0_ngm", "R1_pair_uniform", "p_major", "p_major_lo", "p_major_hi", "attack_major_mean",
                "peak_prev_major_mean", "peak_prev_major_lo", "peak_prev_major_hi", "peak_time_major_mean",
                "peak_time_major_lo", "peak_time_major_hi", "n_rep")}
            TEX.append(" & ".join([
                str(D0), arms[arm], f(v["R0_ngm"], 2) if "R0_ngm" in v else "--",
                f"{v['p_major']:.3f} ({v['p_major_lo']:.3f}--{v['p_major_hi']:.3f})",
                f(100 * v["attack_major_mean"], 1),
                f"{100 * v['peak_prev_major_mean']:.1f} ({100 * v['peak_prev_major_lo']:.1f}--{100 * v['peak_prev_major_hi']:.1f})",
                f"{v['peak_time_major_mean']:.1f} ({v['peak_time_major_lo']:.1f}--{v['peak_time_major_hi']:.1f})"])
                + r" \\")
    out["metro_schedule"] = {k: MET["schedule"][k] for k in ("m_mean", "rho_mean")}
    NUM["dwell"] = out


# ------------------------------------------------------------------------------------------------ 6
def regime():
    """Relative changes (per cent) of the deterministic numbers at D0 = 1, 10, 100, 2400, and the simulated
    attack rates where available."""
    det, sim = REG["det"], REG.get("sim", {})
    out = {"det": det, "sim": sim, "rows": {}}
    NOTE.append("% ---- Table s7_interventions:tab:regime")
    Ds = ("1", "10", "100", "2400")
    c58 = "0.5833"

    def rel(fun):
        vals = [fun(D) for D in Ds]
        return vals

    rows = [
        ("two-zone contrast $c=0.58$", r"$\Rz$",
         rel(lambda D: det[f"twozone|D0={D}|c={c58}"]["R0"] / det[f"twozone|D0={D}|c=0.0000"]["R0"] - 1), 1),
        ("", r"$\Ronebar$",
         rel(lambda D: det[f"twozone|D0={D}|c={c58}"]["R1bar"] / det[f"twozone|D0={D}|c=0.0000"]["R1bar"] - 1), 1),
        ("barriers, $\\rho_B=0.30$, FD", r"$\Ronebar$",
         rel(lambda D: det[f"barriers|D0={D}"]["rhoB=0.30|FD"]["R1bar"] / det[f"barriers|D0={D}"]["rhoB=0.00|FD"]["R1bar"] - 1), 1),
        ("barriers, $\\rho_B=0.30$, DD", r"$\Ronebar$",
         rel(lambda D: det[f"barriers|D0={D}"]["rhoB=0.30|DD"]["R1bar"] / det[f"barriers|D0={D}"]["rhoB=0.00|DD"]["R1bar"] - 1), 0),
        ("capacity 50\\,\\%, FD (office)", r"$\Ronebar$",
         rel(lambda D: det[f"office|D0={D}"]["R1bar_capacity50_FD"] / det[f"office|D0={D}"]["R1bar_baseline"] - 1), 1),
        ("capacity 50\\,\\%, DD (office)", r"$\Ronebar$",
         rel(lambda D: det[f"office|D0={D}"]["R1bar_capacity50_DD"] / det[f"office|D0={D}"]["R1bar_baseline"] - 1), 1),
        ("ventilation, $q\\to q/2$ (office)", r"$\Ronebar$",
         rel(lambda D: det[f"office|D0={D}"]["R1bar_ventilation"] / det[f"office|D0={D}"]["R1bar_baseline"] - 1), 1),
    ]
    for lab, qty, vals, nd in rows:
        out["rows"][f"{lab}|{qty}"] = vals
        cells = [lab, qty]
        for v in vals:
            n = nd if abs(100 * v) >= 0.095 else 2
            cells.append(pct(v, n))
        TEX.append(" & ".join(cells) + r" \\")
    ratio = [det[f"office|D0={D}"]["R1bar_baseline"] / det[f"office|D0={D}"]["R0"] for D in Ds]
    out["rows"]["office R1bar/R0"] = ratio
    TEX.append(" & ".join(["office scene, baseline", r"$\Ronebar/\Rz$"] + [f(v, 2) for v in ratio]) + r" \\")
    NUM["regime"] = out


# ------------------------------------------------------------------------------------------------ 7
def cost():
    out = {"first": BEN, "check_rerun_of_first_script": V12, "check_compiled": V13}
    NOTE.append("% ---- Table s7_interventions:tab:cost")
    b1, b10, r1, r10, c1, c10 = BEN["D0=1"], BEN["D0=10"], V12["D0=1"], V12["D0=10"], V13["1.0"], V13["10.0"]

    def ms(x):
        return f"{1e3 * x:.1f}\\,ms" if x < 0.01 else (f"{1e3 * x:.0f}\\,ms" if x < 0.2 else f"{x:.2f}\\,s")

    TEX.append(" & ".join([
        r"$\Rz=\srad(\ngm)$", "sparse LU, power iteration", ms(b1["T1_R0_meanfield"]["cpu_sparse_s"]),
        ms(b10["T1_R0_meanfield"]["cpu_sparse_s"]),
        f"{b1['T1_R0_meanfield']['speedup_at_1pct']:.0f}; {r1['T1_R0_meanfield']['speedup_at_1pct']:.0f}; {c1['speedup_1pct']:.0f}",
        f"{b10['T1_R0_meanfield']['speedup_at_1pct']:.0f}; {r10['T1_R0_meanfield']['speedup_at_1pct']:.0f}; {c10['speedup_1pct']:.0f}",
        "none"]) + r" \\")
    TEX.append(" & ".join([
        r"$\noff(x)$, one cell", "one sparse solve", f"{1e3 * b1['T2_column_sum']['cpu_solve_s']:.2f}\\,ms",
        f"{1e3 * b10['T2_column_sum']['cpu_solve_s']:.2f}\\,ms",
        f"{b1['T2_column_sum']['speedup_at_1pct']:.0f}; {r1['T2_column_sum']['speedup_at_1pct']:.0f}",
        f"{b10['T2_column_sum']['speedup_at_1pct']:.0f}; {r10['T2_column_sum']['speedup_at_1pct']:.0f}",
        "none"]) + r" \\")
    TEX.append(" & ".join([
        r"$\Ronebar$", "CG on the pair lattice", ms(b1["T3_R1_discrete"]["cpu_pair_s"]),
        ms(b10["T3_R1_discrete"]["cpu_pair_s"]),
        f"{b1['T3_R1_discrete']['speedup_at_1pct']:.0f}; {r1['T3_R1_discrete']['speedup_at_1pct']:.0f}",
        f"{b10['T3_R1_discrete']['speedup_at_1pct']:.0f}; {r10['T3_R1_discrete']['speedup_at_1pct']:.0f}",
        "none"]) + r" \\")
    TEX.append(" & ".join([
        r"$P(\text{major})$", "mean-field branching", ms(b1["T4_p_major"]["cpu_branching_s"]),
        ms(b10["T4_p_major"]["cpu_branching_s"]), "--", "--",
        f"${b1['T4_p_major']['theory_bias']:+.2f}$; ${b10['T4_p_major']['theory_bias']:+.2f}$"]) + r" \\")
    TEX.append(" & ".join([
        "attack rate, major", r"lattice equations \eqref{s2_model:eq:meanfield}",
        ms(b1["T5_attack_rate"]["cpu_ode_s"]), ms(b10["T5_attack_rate"]["cpu_ode_s"]), "--", "--",
        f"${b1['T5_attack_rate']['theory_bias']:+.2f}$; ${b10['T5_attack_rate']['theory_bias']:+.2f}$"]) + r" \\")
    for Dk, b, r, c in (("D0=1", b1, r1, c1), ("D0=10", b10, r10, c10)):
        NOTE.append(f"% {Dk}: speed-up at 0.1 %: {b['T1_R0_meanfield']['speedup_at_0p1pct']:.0f} / "
                    f"{r['T1_R0_meanfield']['speedup_at_0p1pct']:.0f} / {c['speedup_01pct']:.0f}; abs err sparse "
                    f"{b['T1_R0_meanfield']['abs_err_sparse']:.1e}; pair solve theory-MC in SE "
                    f"{b['T3_R1_discrete']['theory_minus_mc_in_se']:.2f}; MC CPU for 1 % attack "
                    f"{b['T5_attack_rate']['mc_cpu_for_1pct_s']:.2f} s ({b['T5_attack_rate']['mc_runs_for_1pct']:.0f} runs); "
                    f"C core {b['c_core_cpu_per_run_s'] * 1e3:.2f} ms/run, {b['events_per_run']:.0f} events")
    NUM["cost"] = out


if __name__ == "__main__":
    hetero()
    barriers()
    occupancy()
    interventions()
    dwell()
    regime()
    cost()
    (ART / "s7_interventions_numbers.json").write_text(json.dumps(NUM, indent=1, default=float))
    (ART / "s7_interventions_tables.tex").write_text("\n".join(NOTE + TEX) + "\n")
    print("\n".join(TEX))
    print("\n".join(NOTE))
    if "--check" in sys.argv:
        tex = " ".join(" ".join(p.read_text().split()) for p in TEXFILES)
        tex = re.sub(r"\\(eq)?ref\{[MS]-", lambda m: "\\" + (m.group(1) or "") + "ref{", tex)   # xr-hyper prefixes
        missing = [r for r in TEX if " ".join(r.split()) not in tex]
        print(f"\n{len(TEX) - len(missing)} of {len(TEX)} generated table rows found verbatim in "
              + " + ".join(p.name for p in TEXFILES))
        for r in missing:
            print("MISSING:", r)
        sys.exit(1 if missing else 0)
