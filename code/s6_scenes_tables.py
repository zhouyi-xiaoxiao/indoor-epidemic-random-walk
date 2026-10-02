#!/usr/bin/env python
"""s6_scenes_tables.py -- numbers and LaTeX table bodies of Section 6 and Supplementary Section S5.

Nothing is simulated here.  The script reads the result files of the simulation part and of its
re-check, derives the quantities quoted in the prose of sections/m6_simulation.tex and sections/supp_sim.tex, and
writes the three table bodies, so that no number in the section is typed by hand.

Inputs (relative to the repository root)
  data/bsc_sim/03_scenes.json                 scene runs (2 000 epidemics per row)
  data/bsc_sim/03_curves_<scene>_D<D0>_range1m.npz   per-run outcomes and mean curves
  data/bsc_sim/03b_outbreak_probability.json  offspring law of a lone index (40 000 runs)
  data/bsc_sim/04_spatial.json, 04_spatial_office_D{1,10}.npz   infection maps (40 000 outbreaks)
  data/bsc_sim/04b_first_generation.json      first generation alone (100 000 runs)
  code/bsc_sim/verify/v05_scenes.json, v07_hotspot.json   re-check's second simulator
  data/s6_scenes_ode_pointseed.json        (code/s6_scenes_ode_pointseed.py)

Outputs
  data/s6_scenes_numbers.json     every derived number quoted in the section
  data/s6_scenes_tables.tex       bodies of Tables s6_scenes:tab:scenes, :tab:pmajor, :tab:hotspot

Run:  python code/s6_scenes_tables.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                                   # repository root
SIM = ROOT / "code" / "bsc_sim"
sys.path.insert(0, str(SIM / "src"))
from bsc_sim.scenes import load_scene                    # noqa: E402

DATA = HERE.parent / "data"
SCENES = ["office", "supermarket", "classroom", "metro"]
TITLE = dict(office="Office", supermarket="Supermarket", classroom="Classroom", metro="Metro carriage")
D0S = (1, 10, 100)


def J(p):
    return json.load(open(p))


S = J(SIM / "data" / "03_scenes.json")
B = J(SIM / "data" / "03b_outbreak_probability.json")
SP = J(SIM / "data" / "04_spatial.json")
FG = J(SIM / "data" / "04b_first_generation.json")
V5 = J(SIM / "verify" / "v05_scenes.json")
V7 = J(SIM / "verify" / "v07_hotspot.json")
ODE_PT = DATA / "s6_scenes_ode_pointseed.json"
OP = J(ODE_PT) if ODE_PT.exists() else {}

num = {}          # derived numbers
tex = []          # LaTeX lines


def f(x, d=2):
    return f"{x:.{d}f}"


def pct(x, d=1):
    return f"{100 * x:.{d}f}"


def neg(s):
    """typographic minus for LaTeX text mode numbers"""
    return s.replace("-", "$-$")


# =====================================================================================================
# Table s6_scenes:tab:scenes
# =====================================================================================================
def scene_rows():
    out = []
    ncol = 12
    for name in SCENES:
        kernels = ["range1m"] + (["samecell"] if f"{name}|D0=1|samecell" in S else [])
        for kern in kernels:
            r1 = S[f"{name}|D0=1|{kern}"]
            if kern == "samecell":
                klabel = "same-cell kernel"
            elif abs(r1["kernel_cells"] - 1.0) < 1e-9:
                klabel = "1\\,m kernel (identical to the same-cell kernel)"
            else:
                klabel = f"1\\,m kernel ({r1['kernel_cells']:.1f} cells on average)"
            head = (f"\\emph{{{TITLE[name]}}}, $\\Npop={r1['N']}$, $\\ncell={r1['cells']}$, {klabel}; "
                    f"$\\beta\\langle q\\rangle/\\gamma={r1['R0_wellmixed']:.2f}$, "
                    f"$\\beta\\qmax/\\gamma={r1['R0_frozen']:.2f}$")
            out.append("\\addlinespace")
            out.append(f"\\multicolumn{{{ncol}}}{{@{{}}l}}{{{head}}}\\\\")
            for D0 in D0S:
                r = S[f"{name}|D0={D0}|{kern}"]
                cells = [str(D0), f(r["R0_ngm"]), f(r["R1_pair_uniform_index"]),
                         f"{r['R_index_alone_mean']:.2f}\\,$\\pm$\\,{r['R_index_alone_se']:.2f}",
                         f"{r['p_major']:.3f} ({r['p_major_lo']:.3f}--{r['p_major_hi']:.3f})",
                         f"{r['p_major_branching']:.3f}"]
                if r["n_major"] > 0:
                    cells += [pct(r["attack_major_mean"]), pct(r["attack_all_mean"]),
                              pct(r["peak_prev_major_mean"]), f(r["peak_time_major_mean"], 1)]
                else:
                    cells += ["--", pct(r["attack_all_mean"]), "--", "--"]
                cells.append(f"{pct(r['ode_attack'])} / {pct(r['ode_peak_prev'])} / {r['ode_peak_time']:.1f}")
                k = f"{name}|D0={D0}|range1m"
                if kern == "range1m" and k in OP:
                    p = OP[k]["point_seed"]
                    cells.append(f"{pct(p['peak_mean'])} / {p['day_mean']:.1f}")
                else:
                    cells.append("")
                out.append(" & ".join(cells) + "\\\\")
    return out


tex.append("% ---- body of Table s6_scenes:tab:scenes "
           "(src: data/bsc_sim/03_scenes.json; last column data/s6_scenes_ode_pointseed.json)")
tex += scene_rows()
tex.append("")

# ---- derived numbers for the prose ------------------------------------------------------------------
sc = {}
for name in SCENES:
    for D0 in D0S:
        r = S[f"{name}|D0={D0}|range1m"]
        d = dict(R0=r["R0_ngm"], R1bar=r["R1_pair_uniform_index"], Rhat=r["R_index_alone_mean"],
                 Rhat_se=r["R_index_alone_se"], p_major=r["p_major"], p_branch=r["p_major_branching"],
                 branch_minus_sim=r["p_major_branching"] - r["p_major"],
                 die_out=1 - r["p_major"], n_major=r["n_major"], n_rep=r["n_rep"],
                 attack_major=r["attack_major_mean"], attack_major_lo=r["attack_major_lo"],
                 attack_major_hi=r["attack_major_hi"], attack_major_sd=r["attack_major_sd"],
                 attack_all=r["attack_all_mean"],
                 ode_attack=r["ode_attack"], ode_minus_sim_attack_points=100 * (r["ode_attack"] - r["attack_major_mean"]),
                 sim_attack_below_ode_percent=100 * (1 - r["attack_major_mean"] / r["ode_attack"]),
                 peak=r["peak_prev_major_mean"], peak_day=r["peak_time_major_mean"],
                 ode_peak=r["ode_peak_prev"], ode_day=r["ode_peak_time"],
                 duration_major=r["duration_major_mean"],
                 index_offspring_full_epidemic=r["index_offspring_mean"],
                 gen2_over_gen1=r["gen2_over_gen1"], rho_K1=r.get("R1_pair_ngm"),
                 z_Rhat_vs_R1bar=(r["R_index_alone_mean"] - r["R1_pair_uniform_index"]) / r["R_index_alone_se"],
                 Rhat_over_R0=r["R_index_alone_mean"] / r["R0_ngm"],
                 people_within_reach=r["occupancy"] * r["kernel_cells"],
                 percell_rule_index_cases=r["percell_rule_index_offspring"], percell_rule_p_major=r["percell_rule_p_major"],
                 percell_rule_attack_all=r["percell_rule_attack_all"])
        c = np.load(SIM / "data" / f"03_curves_{name}_D{D0}_range1m.npz")
        ar = c["attack"]
        d["P_attack_ge_5pct"] = float((ar >= 0.05 - 1e-12).mean())
        d["P_attack_ge_10pct"] = float((ar >= 0.10 - 1e-12).mean())
        d["P_attack_ge_20pct"] = float((ar >= 0.20 - 1e-12).mean())
        d["frac_runs_attack_10_to_50pct"] = float(((ar >= 0.10 - 1e-12) & (ar < 0.5)).mean())
        if "mean_major" in c:
            k = int(np.argmax(c["mean_major"]))
            d["mean_curve_major_peak"] = float(c["mean_major"][k])
            d["mean_curve_major_peak_day"] = float(c["t"][k])
        sc[f"{name}|D0={D0}"] = d
    for D0 in D0S:
        k = f"{name}|D0={D0}|samecell"
        if k in S:
            r = S[k]
            sc[f"{name}|D0={D0}|samecell"] = dict(
                R0=r["R0_ngm"], R1bar=r["R1_pair_uniform_index"], Rhat=r["R_index_alone_mean"],
                Rhat_se=r["R_index_alone_se"], p_major=r["p_major"], p_major_lo=r["p_major_lo"],
                p_major_hi=r["p_major_hi"], n_major=r["n_major"], n_rep=r["n_rep"], p_branch=r["p_major_branching"],
                attack_all=r["attack_all_mean"], gen2_over_gen1=r["gen2_over_gen1"], rho_K1=r.get("R1_pair_ngm"),
                z_Rhat_vs_R1bar=(r["R_index_alone_mean"] - r["R1_pair_uniform_index"]) / r["R_index_alone_se"])
num["scenes"] = sc

zs = {k: v["z_Rhat_vs_R1bar"] for k, v in sc.items()}
kmax = max(zs, key=lambda k: abs(zs[k]))
num["Rhat_vs_R1bar"] = dict(n_rows=len(zs), largest_abs_z=abs(zs[kmax]), row=kmax,
                            n_rows_abs_z_gt_2=int(sum(abs(z) > 2 for z in zs.values())))

summ = {}
for D0 in D0S:
    rows = [sc[f"{n}|D0={D0}"] for n in SCENES]
    summ[f"D0={D0}"] = dict(
        branch_minus_sim_min=min(r["branch_minus_sim"] for r in rows),
        branch_minus_sim_max=max(r["branch_minus_sim"] for r in rows),
        ode_minus_sim_attack_points_min=min(r["ode_minus_sim_attack_points"] for r in rows),
        ode_minus_sim_attack_points_max=max(r["ode_minus_sim_attack_points"] for r in rows),
        sim_attack_below_ode_percent_min=min(r["sim_attack_below_ode_percent"] for r in rows),
        sim_attack_below_ode_percent_max=max(r["sim_attack_below_ode_percent"] for r in rows),
        percell_rule_index_cases_min=min(r["percell_rule_index_cases"] for r in rows),
        percell_rule_index_cases_max=max(r["percell_rule_index_cases"] for r in rows),
        percell_rule_attack_all_min=min(r["percell_rule_attack_all"] for r in rows),
        percell_rule_attack_all_max=max(r["percell_rule_attack_all"] for r in rows),
        percell_rule_p_major_min=min(r["percell_rule_p_major"] for r in rows),
        percell_rule_p_major_max=max(r["percell_rule_p_major"] for r in rows),
        gen2_over_gen1={n: sc[f"{n}|D0={D0}"]["gen2_over_gen1"] for n in SCENES},
        frac_runs_attack_10_to_50pct={n: sc[f"{n}|D0={D0}"]["frac_runs_attack_10_to_50pct"] for n in SCENES})
num["summary_by_D0"] = summ

# re-check's second simulator (D0 = 1 and 10, 1 m kernel)
ver = {}
for name in SCENES:
    for D0 in (1, 10):
        k = f"{name}|{D0:.1f}|1m"
        if k in V5:
            v = V5[k]
            pt, nt = sc[f"{name}|D0={D0}"]["p_major"], sc[f"{name}|D0={D0}"]["n_rep"]
            se = np.sqrt(pt * (1 - pt) / nt + v["p_major"] * (1 - v["p_major"]) / v["n"])
            ver[f"{name}|D0={D0}"] = dict(n=v["n"], p_major=v["p_major"], p_lo=v["p_lo"], p_hi=v["p_hi"],
                                           z_first_minus_check=float((pt - v["p_major"]) / se),
                                           gen2_over_gen1=v["gen2_over_gen1"],
                                           peak_exact_minus_grid_points=100 * (v["peak_exact_major"]
                                                                               - v["peak_grid_major"]),
                                           attack_major=v["attack_major"], peak_grid_major=v["peak_grid_major"],
                                           peak_exact_major=v["peak_exact_major"],
                                           peak_day_grid_major=v["peak_day_grid_major"], R_alone=v["R_alone"],
                                           R_alone_se=v["R_alone_se"],
                                           first_p_major=sc[f"{name}|D0={D0}"]["p_major"])
for k in ("metro|1.0|same", "classroom|1.0|same"):
    v = V5[k]
    ver[k] = dict(n=v["n"], p_major=v["p_major"], p_hi=v["p_hi"], attack_all=v["attack_all"],
                  gen2_over_gen1=v["gen2_over_gen1"], R_alone=v["R_alone"], R_alone_se=v["R_alone_se"])
num["check_scenes"] = ver

if OP:
    num["ode_point_seed"] = {k: dict(uniform=v["uniform_seed"], point=v["point_seed"], n_seeds=v["n_seeds"])
                             for k, v in OP.items()}


# =====================================================================================================
# Table s6_scenes:tab:pmajor
# =====================================================================================================
tex.append("% ---- body of Table s6_scenes:tab:pmajor (src: data/bsc_sim/03b_outbreak_probability.json)")
pm = {}
for D0 in (1, 10):
    tex.append("\\addlinespace")
    for name in SCENES:
        b = B[f"{name}|D0={D0}"]
        geo = 1 - 1 / b["offspring_mean"]
        tex.append(" & ".join([
            TITLE[name] if name != "metro" else "Metro", str(D0), f(b["offspring_mean"]), f(b["offspring_var"], 1),
            f(b["var_over_mean"]), f(b["p_zero"], 3), f(b["p_major_poisson_same_mean"], 3),
            f(b["p_major_branching_meanfield"], 3), f(b["p_major_gw_empirical"], 3),
            f"{b['p_major_sim']:.3f} ({b['p_major_sim_lo']:.3f}--{b['p_major_sim_hi']:.3f})"]) + "\\\\")
        pm[f"{name}|D0={D0}"] = dict(
            mean=b["offspring_mean"], var=b["offspring_var"], var_over_mean=b["var_over_mean"],
            geometric_var_over_mean=1 + b["offspring_mean"], p_major_geometric_same_mean=geo,
            gw_minus_sim=b["p_major_gw_empirical"] - b["p_major_sim"],
            branching_minus_sim=b["p_major_branching_meanfield"] - b["p_major_sim"],
            poisson_minus_sim=b["p_major_poisson_same_mean"] - b["p_major_sim"])
tex.append("")
for name in ("office", "supermarket", "classroom"):
    b = B[f"{name}|D0=1"]
    p5 = sc[f"{name}|D0=1"]["P_attack_ge_5pct"]
    pm[f"{name}|D0=1"]["gw_minus_sim_5pct"] = b["p_major_gw_empirical"] - p5
    pm[f"{name}|D0=1"]["branching_minus_sim_5pct"] = b["p_major_branching_meanfield"] - p5
num["pmajor"] = pm
num["threshold_sensitivity_D0_10_max_change"] = max(
    abs(sc[f"{n}|D0=10"]["P_attack_ge_5pct"] - sc[f"{n}|D0=10"]["P_attack_ge_20pct"]) for n in SCENES)
num["pmajor_ranges"] = dict(
    var_over_mean_min=min(v["var_over_mean"] for v in pm.values()),
    var_over_mean_max=max(v["var_over_mean"] for v in pm.values()),
    geometric_var_over_mean_min=min(v["geometric_var_over_mean"] for v in pm.values()),
    geometric_var_over_mean_max=max(v["geometric_var_over_mean"] for v in pm.values()),
    gw_minus_sim_D0_1_three_scenes=[pm[f"{n}|D0=1"]["gw_minus_sim"] for n in ("office", "supermarket", "classroom")],
    gw_minus_sim_other_rows_max_abs=max(abs(pm[k]["gw_minus_sim"]) for k in pm
                                        if k not in ("office|D0=1", "supermarket|D0=1", "classroom|D0=1")))


# =====================================================================================================
# Table s6_scenes:tab:hotspot  (a) correlations, (b) office cell classes
# =====================================================================================================
tex.append("% ---- body of Table s6_scenes:tab:hotspot, panel (a) "
           "(src: data/bsc_sim/04_spatial.json, 04b_first_generation.json)")
hs = {}
for D0 in (1, 10):
    tex.append("\\addlinespace")
    for name in SCENES:
        r = SP[f"{name}|D0={D0}"]
        g = [r[f"gen{i}"] for i in (1, 2, 3)]
        fg = FG[f"{name}|D0={D0}"]
        tex.append(" & ".join([
            TITLE[name] if name != "metro" else "Metro", str(D0),
            neg(" / ".join(f(x["r_meanfield"]) for x in g)),
            neg(" / ".join(f(x["r_pair"]) for x in g)),
            neg(f(fg["r_meanfield"], 3)), f(fg["r_pair"], 3), f(fg["r_expected_if_exact"], 3),
            " / ".join([f(g[1]["mean_cases_per_run"], 1), f(g[1]["pred_cases_per_run_pair"], 1),
                        f(g[1]["pred_cases_per_run_mf"], 1)])]) + "\\\\")
        hs[f"{name}|D0={D0}"] = dict(
            r_meanfield=[x["r_meanfield"] for x in g], r_pair=[x["r_pair"] for x in g],
            noise_ceiling_multinomial=[x["r_noise_ceiling"] for x in g],
            cases_per_run_sim=[x["mean_cases_per_run"] for x in g],
            cases_per_run_pair=[x["pred_cases_per_run_pair"] for x in g],
            cases_per_run_meanfield=[x["pred_cases_per_run_mf"] for x in g],
            first_gen_alone=dict(r_pair=fg["r_pair"], r_meanfield=fg["r_meanfield"],
                                 r_expected_if_exact=fg["r_expected_if_exact"], chi2=fg["chi2_pair"],
                                 dof=fg["dof"], chi2_over_dof=fg["chi2_pair"] / fg["dof"],
                                 p_value=float(chi2_dist.sf(fg["chi2_pair"], fg["dof"])),
                                 n_rep=fg["n_rep"]))
tex.append("")
num["hotspot_correlations"] = hs
_p = {k: v["first_gen_alone"]["p_value"] for k, v in hs.items()}
num["first_generation_chi2"] = dict(
    p_values=_p, smallest=min(_p, key=_p.get),
    chi2_over_dof_min=min(v["first_gen_alone"]["chi2_over_dof"] for v in hs.values()),
    chi2_over_dof_max=max(v["first_gen_alone"]["chi2_over_dof"] for v in hs.values()))
num["hotspot_pair_r_range_D0_1"] = [min(min(hs[f"{n}|D0=1"]["r_pair"]) for n in SCENES),
                                    max(max(hs[f"{n}|D0=1"]["r_pair"]) for n in SCENES)]
num["hotspot_meanfield_r_range_D0_1"] = [min(min(hs[f"{n}|D0=1"]["r_meanfield"]) for n in SCENES),
                                         max(max(hs[f"{n}|D0=1"]["r_meanfield"]) for n in SCENES)]
num["check_hotspot"] = V7


def office_classes(D0):
    """Relative infection rate (share of infections / share of walkable floor) by cell class."""
    s = load_scene("office", D0=float(D0))
    z = np.load(SIM / "data" / f"04_spatial_office_D{D0}.npz")
    zone, xy, idx = s.site_zone, s.site_xy, s.site_index

    def nbrs(i):
        x, y = xy[i]
        o = []
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < s.nx and 0 <= yy < s.ny and idx[yy, xx] >= 0:
                o.append(idx[yy, xx])
        return o

    border = np.array([zone[i] == "W" and any(zone[j] == "C" for j in nbrs(i)) for i in range(s.M)])
    classes = [("W_border", "Workstations bordering a corridor", border),
               ("W_interior", "Other workstations", (zone == "W") & ~border),
               ("C", "Corridor", zone == "C"), ("M", "Meeting room", zone == "M"),
               ("K", "Kitchen", zone == "K")]
    res = dict(n_cells=s.M, classes={})
    for key, label, m in classes:
        e = dict(label=label, cells=int(m.sum()), floor_share=float(m.mean()), sim=[], pair=[], meanfield=[])
        for g in (1, 2, 3):
            c = z[f"counts_gen{g}"]
            e["sim"].append(float(c[m].sum() / c.sum() / m.mean()))
            e["pair"].append(float(z[f"pred_pair_gen{g}"][m].sum() / m.mean()))
            e["meanfield"].append(float(z[f"pred_mf_gen{g}"][m].sum() / m.mean()))
        res["classes"][key] = e
    res["top15"] = {}
    for g in (1, 2, 3):
        c = z[f"counts_gen{g}"]
        top = np.argsort(-c)[:15]
        res["top15"][f"gen{g}"] = dict(zones="".join(zone[top]), n_border_workstations=int(border[top].sum()),
                                       hottest_cell_relative_rate=float(c.max() / c.sum() * s.M),
                                       coldest_cell_relative_rate=float(c.min() / c.sum() * s.M))
    return res, border


tex.append("% ---- body of Table s6_scenes:tab:hotspot, panel (b) "
           "(src: data/s6_scenes_numbers.json: office_cell_classes; recomputed from "
           "data/bsc_sim/04_spatial_office_D{1,10}.npz)")
oc = {}
for D0 in (1, 10):
    res, _ = office_classes(D0)
    oc[f"D0={D0}"] = res
    tex.append("\\addlinespace")
    for key, e in res["classes"].items():
        tex.append(" & ".join([
            e["label"], str(D0), str(e["cells"]), pct(e["floor_share"]),
            " / ".join(f(x) for x in e["sim"]), " / ".join(f(x) for x in e["pair"]),
            " / ".join(f(x) for x in e["meanfield"])]) + "\\\\")
num["office_cell_classes"] = oc
num["office_first_generation_zone_shares_check"] = {k: V7[k] for k in ("zone_W", "zone_C", "zone_M", "zone_K")}

# ranking of the scenes
num["ranking_D0_1"] = dict(
    by_R0=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=1"]["R0"]),
    by_Rhat=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=1"]["Rhat"]),
    by_p_major=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=1"]["p_major"]),
    by_attack_all=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=1"]["attack_all"]),
    by_attack_major=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=1"]["attack_major"]))
num["ranking_D0_10"] = dict(
    by_R0=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=10"]["R0"]),
    by_Rhat=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=10"]["Rhat"]),
    by_p_major=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=10"]["p_major"]),
    by_attack_all=sorted(SCENES, key=lambda n: -sc[f"{n}|D0=10"]["attack_all"]))

DATA.mkdir(exist_ok=True)
json.dump(num, open(DATA / "s6_scenes_numbers.json", "w"), indent=1)
open(DATA / "s6_scenes_tables.tex", "w").write("\n".join(tex) + "\n")
print("wrote", DATA / "s6_scenes_numbers.json", "and", DATA / "s6_scenes_tables.tex")
