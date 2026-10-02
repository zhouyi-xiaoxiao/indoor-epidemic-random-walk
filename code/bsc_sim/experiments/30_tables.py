#!/usr/bin/env python
"""Print Markdown tables from the JSON results.

    python experiments/30_tables.py scenes | mobility | spatial | hetero |
                                   obstacles | occupancy | interventions |
                                   metro | benchmark | dwell
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
SCENES = ["office", "supermarket", "classroom", "metro"]


def J(name):
    return json.load(open(os.path.join(DATA, name)))


def f(x, n=2):
    if x is None or (isinstance(x, float) and x != x):
        return "—"
    return f"{x:.{n}f}"


def ci(r, key, n=2):
    if r.get(key + "_mean") is None:
        return "—"
    return f"{f(r[key + '_mean'], n)} ({f(r[key + '_lo'], n)}–{f(r[key + '_hi'], n)})"


def scenes():
    S = J("03_scenes.json")
    for D0 in (1, 10, 100):
        print(f"\n**D0 = {D0} cells²/day** (2 000 epidemics per row, one uniformly placed index case)\n")
        print("| scene | kernel | N / cells | mean-field R0 (bounds) | pair R1 | counted, index alone (±SE) | P(major) sim (95 % CI) | P(major) branching | attack rate, major (95 % CI) | attack, all runs | peak prevalence, major | peak day, major | duration, major (d) | ODE: attack / peak / day |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for name in SCENES:
            for kern in ("range1m", "samecell"):
                k = f"{name}|D0={D0}|{kern}"
                if k not in S:
                    continue
                r = S[k]
                kn = {"range1m": "1 m", "samecell": "same cell"}[kern]
                if kern == "range1m" and r["kernel_cells"] == 1.0:
                    kn = "1 m (= same cell)"
                print(f"| {name} | {kn} | {r['N']} / {r['cells']} | "
                      f"{f(r['R0_ngm'])} ({f(r['R0_wellmixed'])}–{f(r['R0_frozen'])}) | "
                      f"{f(r['R1_pair_uniform_index'])} | "
                      f"{f(r['R_index_alone_mean'])} ± {f(r['R_index_alone_se'])} | "
                      f"{f(r['p_major'], 3)} ({f(r['p_major_lo'], 3)}–{f(r['p_major_hi'], 3)}) | "
                      f"{f(r['p_major_branching'], 3)} | {ci(r, 'attack_major', 3)} | "
                      f"{f(r['attack_all_mean'], 3)} | {ci(r, 'peak_prev_major', 3)} | "
                      f"{ci(r, 'peak_time_major', 1)} | {ci(r, 'duration_major', 1)} | "
                      f"{f(r['ode_attack'], 3)} / {f(r['ode_peak_prev'], 3)} / {f(r['ode_peak_time'], 1)} |")
    print("\n**Second-generation growth and the per-cell rule (D0 = 1)**\n")
    print("| scene | kernel | pair-level ρ(K1) | counted gen-2 / gen-1 (full epidemic) | growth rate λ of mean field (1/d) | per-cell rule: index secondary cases | per-cell rule: P(attack ≥ 10 %) | per-cell rule: mean attack |")
    print("|---|---|---|---|---|---|---|---|")
    for name in SCENES:
        for kern in ("range1m", "samecell"):
            k = f"{name}|D0=1|{kern}"
            if k not in S:
                continue
            r = S[k]
            print(f"| {name} | {kern} | {f(r.get('R1_pair_ngm'))} | {f(r['gen2_over_gen1'])} | "
                  f"{f(r['growth_rate'], 3)} | {f(r['percell_rule_index_offspring'])} | "
                  f"{f(r['percell_rule_p_major'], 3)} | {f(r['percell_rule_attack_all'], 3)} |")
    print("\n**Spread of single outcomes (major outbreaks, D0 = 1 and 10, 1 m kernel): mean ± SD [5 %, 95 %]**\n")
    print("| scene | D0 | n major / n | attack rate | peak prevalence | peak day |")
    print("|---|---|---|---|---|---|")
    for D0 in (1, 10):
        for name in SCENES:
            r = S[f"{name}|D0={D0}|range1m"]

            def sp(key, n):
                return (f"{f(r[key + '_mean'], n)} ± {f(r[key + '_sd'], n)} "
                        f"[{f(r[key + '_q05'], n)}, {f(r[key + '_q95'], n)}]")
            print(f"| {name} | {D0} | {r['n_major']} / {r['n_rep']} | {sp('attack_major', 3)} | "
                  f"{sp('peak_prev_major', 3)} | {sp('peak_time_major', 1)} |")


def mobility():
    d = J("02_mobility.json")
    print("| D0 | mean-field R0 | pair ρ(K1) | counted, Perron index (95 % CI) | pair, uniform index | counted, uniform index (95 % CI) | counted gen-2 / gen-1 (95 % CI) |")
    print("|---|---|---|---|---|---|---|")
    for r in sorted(d.values(), key=lambda r: r["D0"]):
        print(f"| {r['D0']:g} | {f(r['R0_ngm'], 3)} | {f(r['R1_pair_ngm'], 3)} | "
              f"{f(r['sim_offspring_perron'], 3)} ({f(r['sim_perron_lo'], 3)}–{f(r['sim_perron_hi'], 3)}) | "
              f"{f(r['R1_pair_uniform_index'], 3)} | "
              f"{f(r['sim_offspring_uniform'], 3)} ({f(r['sim_uniform_lo'], 3)}–{f(r['sim_uniform_hi'], 3)}) | "
              f"{f(r['sim_gen2_over_gen1'], 3)} ({f(r['sim_gen2_lo'], 3)}–{f(r['sim_gen2_hi'], 3)}) |")


def spatial():
    """Pearson r between predicted and simulated cell probabilities, with the
    Fisher-z 95 % interval over the M cells, and the ceiling that sampling
    noise alone would allow for a perfect prediction."""
    d = J("04_spatial.json")

    def fisher(x, y):
        r = np.corrcoef(x, y)[0, 1]
        z, se = np.arctanh(r), 1 / np.sqrt(len(x) - 3)
        return r, np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    print("| scene | D0 | generation | cases counted | r, mean-field K (95 % CI) | r, pair-level K1 (95 % CI) | r ceiling (sampling noise only) | cases per run: sim / pair / mean field |")
    print("|---|---|---|---|---|---|---|---|")
    for k, r in d.items():
        z = np.load(os.path.join(DATA, f"04_spatial_{r['scene']}_D{r['D0']:g}.npz"))
        for g in (1, 2, 3):
            q = r[f"gen{g}"]
            c = z[f"counts_gen{g}"].astype(float)
            rm = fisher(z[f"pred_mf_gen{g}"], c)
            rp = fisher(z[f"pred_pair_gen{g}"], c)
            print(f"| {r['scene']} | {r['D0']:g} | {g} | {q['n_cases']} | "
                  f"{f(rm[0], 3)} ({f(rm[1], 3)}–{f(rm[2], 3)}) | "
                  f"{f(rp[0], 3)} ({f(rp[1], 3)}–{f(rp[2], 3)}) | "
                  f"{f(q['r_noise_ceiling'], 3)} | {f(q['mean_cases_per_run'])} / "
                  f"{f(q['pred_cases_per_run_pair'])} / {f(q['pred_cases_per_run_mf'])} |")


def hetero():
    d = J("05_heterogeneity.json")
    print("**Two-zone office with zone-dependent mobility (workspace D = 0.3 D0, corridor D = 2 D0)**\n")
    print("| D0 | N | contrast c | q_W / q_C | mean-field R0 | bounds | pair ρ(K1) | counted, Perron index (95 % CI) | counted, uniform index | gen-2 / gen-1 | P(major) | attack, all runs |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = [v for k, v in d.items() if k != "meanfield_grid"]
    for r in sorted(rows, key=lambda r: (r["D0"], r["N"], r["c"])):
        print(f"| {r['D0']:g} | {r['N']} | {f(r['c'], 3)} | {f(r['qW'], 2)} / {f(r['qC'], 2)} | "
              f"{f(r['R0_ngm'], 3)} | {f(r['wellmixed'])}–{f(r['frozen'])} | {f(r['R1_pair_ngm'], 3)} | "
              f"{f(r['sim_R_perron'], 3)} ({f(r['sim_R_perron_lo'], 3)}–{f(r['sim_R_perron_hi'], 3)}) | "
              f"{f(r['sim_R_uniform'], 3)} | {f(r['sim_gen2_over_gen1'], 3)} | "
              f"{f(r['p_major'], 3)} ({f(r['p_major_lo'], 3)}–{f(r['p_major_hi'], 3)}) | {f(r['attack_all'], 3)} |")
    e = J("05_heterogeneity_extra.json")
    print("\n**Controls (N = 100): uniform mobility, and the hot zone moved to the mobile corridor (c < 0)**\n")
    print("| variant | D0 | contrast c | q_W / q_C | mean-field R0 | pair ρ(K1) | counted, Perron index (95 % CI) | P(major) (95 % CI) | attack, all runs |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in sorted(e.values(), key=lambda r: (r["D0"], r["variant"], r["c"])):
        vn = {"uniformD": "uniform D = 0.81 D0", "zoneD": "zone D, hot corridor"}[r["variant"]]
        print(f"| {vn} | {r['D0']:g} | {f(r['c'], 3)} | {f(r['qW'], 2)} / {f(r['qC'], 2)} | {f(r['R0_ngm'], 3)} | "
              f"{f(r['R1_pair_ngm'], 3)} | {f(r['sim_R_perron'], 3)} ({f(r['sim_R_perron_lo'], 3)}–{f(r['sim_R_perron_hi'], 3)}) | "
              f"{f(r['p_major'], 3)} ({f(r['p_major_lo'], 3)}–{f(r['p_major_hi'], 3)}) | {f(r['attack_all'], 3)} |")


def obstacles():
    from scipy import stats
    d = J("06_obstacles.json")
    fr = sorted({r["rhoB_nominal"] for r in d.values()})
    for cal in ("FD", "DD"):
        print(f"\n**{cal} calibration** — mean over 8 barrier layouts ± 95 % CI half-width across layouts (500 runs per layout and arm)\n")
        print("| ρ_B nominal | realised | mean-field R0 | pair R1 (index at centre) | counted, index alone | P(major) | attack, all runs | attack, major | peak prevalence, major | peak day, major |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for fv in fr:
            rows = [r for r in d.values() if r["rhoB_nominal"] == fv and r["calibration"] == cal]

            def ag(key, n=3):
                v = np.array([r[key] for r in rows if r.get(key) is not None], float)
                v = v[np.isfinite(v)]
                if len(v) < 2:
                    return f(v.mean(), n) if len(v) else "—"
                h = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
                return f"{f(v.mean(), n)} ± {f(h, n)}"
            print(f"| {fv:.2f} | {ag('rhoB_realised')} | {ag('R0_ngm', 2)} | {ag('R1_pair_centre', 2)} | "
                  f"{ag('R_index_alone', 2)} | {ag('p_major')} | {ag('attack_all_mean')} | "
                  f"{ag('attack_major_mean')} | {ag('peak_prev_major_mean')} | {ag('peak_time_major_mean', 1)} |")
    # paired differences against the barrier-free arm
    z = np.load(os.path.join(DATA, "06_obstacles_runs.npz"))
    print("\n**Paired differences against ρ_B = 0 (same people, same start, same seed), FD, pooled over 8 layouts × 500 replicates: mean ± 95 % CI**\n")
    print("| ρ_B | Δ attack rate (all runs) | Δ peak prevalence (all runs) | Δ secondary cases of index (alone) |")
    print("|---|---|---|---|")
    for fv in fr[1:]:
        out = []
        for key in ("attack", "peak", "R_alone"):
            diff = np.concatenate([
                z[f"layout{lay}|rhoB={fv:.2f}|FD|{key}"].astype(float)
                - z[f"layout{lay}|rhoB=0.00|FD|{key}"].astype(float)
                for lay in range(8)])
            m = diff.mean()
            h = 1.96 * diff.std(ddof=1) / np.sqrt(len(diff))
            out.append(f"{m:+.3f} ± {h:.3f}")
        print(f"| {fv:.2f} | {out[0]} | {out[1]} | {out[2]} |")


def occupancy():
    d = J("07_occupancy.json")
    for D0 in (1, 10):
        print(f"\n**D0 = {D0}**\n")
        print("| N | people/m² | calibration | mean-field R0 | pair R1 | counted, index alone (±SE) | P(major) (95 % CI) | branching | attack, all runs | attack, major | peak prevalence, major | peak day, major |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        rows = [r for k, r in d.items() if k.startswith(f"D0={D0}|")]
        for r in sorted(rows, key=lambda r: (r["N"], r["calibration"] == "DD")):
            print(f"| {r['N']} | {f(r['density_per_m2'], 3)} | {r['calibration']} | {f(r['R0_ngm'])} | "
                  f"{f(r['R1_pair_uniform'])} | {f(r['R_index_alone'])} ± {f(r['R_index_alone_se'])} | "
                  f"{f(r['p_major'], 3)} ({f(r['p_major_lo'], 3)}–{f(r['p_major_hi'], 3)}) | "
                  f"{f(r['p_major_branching'], 3)} | {f(r['attack_all_mean'], 3)} | {ci(r, 'attack_major', 3)} | "
                  f"{ci(r, 'peak_prev_major', 3)} | {ci(r, 'peak_time_major', 1)} |")


def interventions():
    d = J("08_interventions.json")
    for name in ("office", "classroom"):
        for D0 in (1, 10):
            print(f"\n**{name}, D0 = {D0}** (2 000 epidemics per arm)\n")
            print("| arm | calibration | N | mean-field R0 | pair R1 | counted, index alone (±SE) | P(major) (95 % CI) | attack, all runs | attack, major | peak prevalence, major | peak day, major |")
            print("|---|---|---|---|---|---|---|---|---|---|---|")
            for k, r in d.items():
                if r["scene"] != name or r["D0"] != D0:
                    continue
                print(f"| {r['arm']} | {r['calibration']} | {r['N']} | {f(r['R0_ngm'])} | "
                      f"{f(r['R1_pair_uniform'])} | {f(r['R_index_alone'])} ± {f(r['R_index_alone_se'])} | "
                      f"{f(r['p_major'], 3)} ({f(r['p_major_lo'], 3)}–{f(r['p_major_hi'], 3)}) | "
                      f"{f(r['attack_all_mean'], 3)} | {ci(r, 'attack_major', 3)} | "
                      f"{ci(r, 'peak_prev_major', 3)} | {ci(r, 'peak_time_major', 1)} |")


def metro():
    d = J("09_metro_time.json")
    print("| D0 | arm | mean-field R0 | pair R1 | P(major) (95 % CI) | attack, major | peak prevalence, major | peak day, major |")
    print("|---|---|---|---|---|---|---|---|")
    for k, r in d.items():
        if k == "schedule":
            continue
        print(f"| {r['D0']:g} | {r['arm']} | {f(r.get('R0_ngm'))} | {f(r.get('R1_pair_uniform'))} | "
              f"{f(r['p_major'], 3)} ({f(r['p_major_lo'], 3)}–{f(r['p_major_hi'], 3)}) | "
              f"{ci(r, 'attack_major', 3)} | {ci(r, 'peak_prev_major', 3)} | {ci(r, 'peak_time_major', 1)} |")
    s = d["schedule"]
    print(f"\nschedule: mean m = {s['m_mean']:.3f}, mean rho = {s['rho_mean']:.2f} persons/m²")


def benchmark():
    d = J("10_benchmark.json")
    for key, R in d.items():
        print(f"\n**office, {key}** — single-core CPU seconds\n")
        print("| quantity | deterministic method | value | CPU (s) | Monte Carlo estimator | MC value (rel. SE of the run) | MC CPU for the target | speed-up at equal accuracy | deterministic − MC |")
        print("|---|---|---|---|---|---|---|---|---|")
        t = R["T1_R0_meanfield"]
        print(f"| mean-field R0 = ρ(K) | sparse LU + power iteration | {t['value_sparse']:.4f} | {t['cpu_sparse_s']:.2e} | particle power method ({t['mc_burn_in_generations']} burn-in generations) | {t['mc_estimate']:.3f} ({100 * t['mc_rel_se']:.2f} %) | {t['mc_cpu_for_1pct_s']:.2e} (1 %) | {t['speedup_at_1pct']:.0f}× (1 %); {t['speedup_at_0p1pct']:.0f}× (0.1 %) | {t['value_sparse'] - t['mc_estimate']:+.3f} |")
        t = R["T2_column_sum"]
        print(f"| cases of a case born at one cell, (1ᵀK)_x | one sparse solve | {t['value_solve']:.4f} | {t['cpu_solve_s']:.2e} | single-walker path integral | {t['mc_estimate']:.3f} ({100 * t['mc_rel_se']:.2f} %) | {t['mc_cpu_for_1pct_s']:.2e} (1 %) | {t['speedup_at_1pct']:.0f}× | {t['value_solve'] - t['mc_estimate']:+.3f} |")
        t = R["T3_R1_discrete"]
        print(f"| discrete-people R1 (uniform index) | CG on the pair lattice | {t['value_pair_solve']:.4f} | {t['cpu_pair_s']:.2e} | individual-based runs, index alone | {t['mc_estimate']:.3f} ({100 * t['mc_rel_se']:.2f} %) | {t['mc_cpu_for_1pct_s']:.2e} (1 %, {t['mc_runs_for_1pct']:.0f} runs) | {t['speedup_at_1pct']:.0f}× | {t['value_pair_solve'] - t['mc_estimate']:+.3f} |")
        t = R["T4_p_major"]
        print(f"| P(major outbreak) | mean-field branching fixed point | {t['value_branching_meanfield']:.3f} | {t['cpu_branching_s']:.2e} | full epidemics | {t['mc_estimate']:.3f} (SE {t['mc_se']:.3f}) | {t['mc_cpu_for_se_0p01_s']:.2e} (SE 0.01, {t['mc_runs_for_se_0p01']:.0f} runs) | {t['speedup_at_se_0p01']:.0f}× | {t['theory_bias']:+.3f} |")
        t = R["T5_attack_rate"]
        print(f"| attack rate of a major outbreak | mean-field lattice ODE | {t['value_ode']:.3f} | {t['cpu_ode_s']:.2e} | full epidemics | {t['mc_estimate']:.3f} ({100 * t['mc_rel_se']:.2f} %) | {t['mc_cpu_for_1pct_s']:.2e} (1 %, {t['mc_runs_for_1pct']:.0f} runs) | {t['speedup_at_1pct']:.2g}× | {t['theory_bias']:+.3f} |")
        print(f"\nfull epidemic: {R['events_per_run']:.0f} events, C core {R['c_core_cpu_per_run_s']:.2e} s per run; "
              f"pure-Python reference {1e3 * R['python_reference_cpu_per_event_s']:.2f} ms per event, i.e. {R['python_reference_cpu_per_run_s']:.1f} s per run "
              f"({R['python_reference_cpu_per_run_s'] / R['c_core_cpu_per_run_s']:.0f}× slower).")


def dwell():
    d = J("11_dwell_time.json")
    print("| scene | D0 | hours/day (w) | R0 continuous | R0 duty cycle | β⟨q⟩w/γ | pair R1 duty cycle | counted with explicit schedule (95 % CI) | P(major) schedule sim | branching | attack, major |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for k, r in d.items():
        sim = "—"
        if r.get("sim_R_schedule") is not None:
            sim = f"{f(r['sim_R_schedule'], 3)} ({f(r['sim_R_schedule_lo'], 3)}–{f(r['sim_R_schedule_hi'], 3)})"
        pm = "—" if r.get("p_major") is None else f"{f(r['p_major'], 3)} ({f(r['p_major_lo'], 3)}–{f(r['p_major_hi'], 3)})"
        print(f"| {r['scene']} | {r['D0']:g} | {r['hours_per_day']:.2f} ({r['w']:.4f}) | {f(r['R0_continuous'])} | "
              f"{f(r['R0_duty_cycle'], 3)} | {f(r['wellmixed_duty_cycle'], 3)} | {f(r['R1_pair_duty_cycle'], 3)} | "
              f"{sim} | {pm} | {f(r.get('p_major_branching'), 3)} | {f(r.get('attack_major'), 3)} |")
    print("\n| scene | D0 | one visit (h) | infections per visit: mean field | exposure counted (95 % CI) | infections counted (95 % CI) |")
    print("|---|---|---|---|---|---|")
    for k, r in d.items():
        print(f"| {r['scene']} | {r['D0']:g} | {r['visit_hours']:.2f} | {r['per_visit_meanfield']:.4f} | "
              f"{r['per_visit_exposure_sim']:.4f} ({r['per_visit_exposure_lo']:.4f}–{r['per_visit_exposure_hi']:.4f}) | "
              f"{r['per_visit_sim']:.4f} ({r['per_visit_sim_lo']:.4f}–{r['per_visit_sim_hi']:.4f}) |")


if __name__ == "__main__":
    for name in sys.argv[1:]:
        globals()[name]()
