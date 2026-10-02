"""s2_model_scene_table.py -- numbers of Section 2 ("Model, scenes and parameter regime").

Recomputes, from the canonical scene files code/bsc_sim/scenes/*.json and the generator and
contact kernel of code/bsc_sim/src/bsc_sim/theory.py:

  * Table s2_model:tab:scenes -- lattice, spacing a, walkable cells n, people N, rho_bar = (N-1)/n,
    zone table (cells, D/D0, q), area means of q and of D/D0, beta<q>/gamma, beta q_max/gamma,
    mean number of cells in the 1 m contact kernel, mean hop rate per person at D0 = 1;
  * Section s2_model:sec:regime -- what D0 = 1 cell^2/day means in metres, the Green-Kubo walking
    mobility D = v^2 tau / 2 in cell^2/day for each lattice spacing, R0 of every scene (1 m kernel)
    at D0 = 1, 10, 100 and at the walking values, and the exit rate of the office meeting room.

Cross-checks against ../data/plan_theory_checks.json are stored under "crosscheck".

Outputs: ../data/s2_model_scenes.json and ../data/s2_model_table_rows.tex (rows of the table, so that
the numbers in the article are not retyped by hand).

Run:  python s2_model_scene_table.py        (a few seconds, < 200 MB)

Units: lengths in lattice cells, time in days, D0 in cell^2/day; beta = 0.5/day, gamma = 0.14/day.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                       # repository root
sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))

from bsc_sim import scenes as sc             # noqa: E402
from bsc_sim import theory as th             # noqa: E402

BETA, GAMMA = 0.5, 0.14
V_WALK, TAU = 1.3, 1.5                       # assumed walking speed (m/s) and persistence time (s)
WALK_FRACTION = 0.05                         # fraction of the time spent walking (assumption)
DAY = 86400.0
OUT = HERE.parent / "data"
OUT.mkdir(exist_ok=True)


def r0(scene, P, D0):
    """rho(P^T diag(beta q)(gamma - D0 L0)^-1) for the scene's zone multipliers scaled by D0."""
    L = th.generator(scene, D0 * scene.D).toarray()
    n = scene.M
    K = (BETA * scene.q)[:, None] * np.linalg.inv(GAMMA * np.eye(n) - L)
    if P is not None:
        K = P.T @ K
    return float(np.max(np.linalg.eigvals(K).real))


def main():
    res: dict = {"beta": BETA, "gamma": GAMMA, "beta_over_gamma": BETA / GAMMA, "scenes": {}}

    # ---- walking mobility (Green-Kubo, two dimensions): D = v^2 tau / 2 -------------------------
    D_walk_m2_s = V_WALK ** 2 * TAU / 2.0
    D_walk_m2_day = D_walk_m2_s * DAY
    res["walking"] = {
        "v_m_per_s": V_WALK, "tau_s": TAU, "walk_fraction": WALK_FRACTION,
        "D_m2_per_s": D_walk_m2_s, "D_m2_per_day": D_walk_m2_day,
    }

    rows_a, rows_b = [], []
    for name in sc.SCENE_NAMES:
        s = sc.load_scene(name)              # D0 = 1: s.D are the zone multipliers d_x
        n, N, a = s.M, s.N, s.a_m
        zones = s.meta["zones"]
        rho = s.meta["rho_per_m2"]
        zt = {}
        for z, zp in zones.items():
            cells = int((s.zone == z).sum())
            if zp.get("barrier", False):
                zt[z] = {"name": zp["name_en"], "cells": cells, "barrier": True}
                continue
            d = zp["D"] if "D" in zp else zp["D_coef"] / rho ** zp["D_exp"]
            zt[z] = {"name": zp["name_en"], "cells": cells, "d": float(d), "q": float(zp["q"])}
        P1 = th.contact_kernel(s, th.kernel_radius_cells(s, 1.0))
        kernel_sizes = np.diff(P1.indptr)
        ptr, idx, W = th.hop_rates(s)
        out_rate = np.add.reduceat(W, ptr[:-1])            # total hop rate of a person in each cell
        Dw = D_walk_m2_day / a ** 2                        # walking mobility in cell^2/day
        r = {
            "lattice": [s.nx, s.ny], "a_m": a, "floor_m": [s.nx * a, s.ny * a],
            "cells_total": s.nx * s.ny, "cells_walkable": n, "cells_barrier": s.nx * s.ny - n,
            "N": N, "people_per_cell": N / n, "rho_bar": (N - 1) / n,
            "people_per_m2_walkable": N / (n * a * a), "density_used_for_metro_rule_per_m2": rho,
            "zones": zt,
            "q_mean": float(s.q.mean()), "q_min": float(s.q.min()), "q_max": float(s.q.max()),
            "d_mean": float(s.D.mean()), "d_min": float(s.D.min()), "d_max": float(s.D.max()),
            "R0_wellmixed": BETA * float(s.q.mean()) / GAMMA,
            "R0_frozen_samecell": BETA * float(s.q.max()) / GAMMA,
            "kernel1m_cells_mean": float(kernel_sizes.mean()),
            "kernel1m_cells_max": int(kernel_sizes.max()),
            "kernel1m_is_samecell": bool(kernel_sizes.max() == 1),
            "hop_rate_mean_per_day_D0=1": float(out_rate.mean()),
            "hop_rate_min_per_day_D0=1": float(out_rate.min()),
            "hop_rate_max_per_day_D0=1": float(out_rate.max()),
            "D0_cell2_per_day_in_m2_per_day": a * a,
            "rms_displacement_m_per_day_at_D0=1": float(np.sqrt(4.0) * a),
            "walking_D0_cell2_per_day": Dw,
            "walking_5pct_D0_cell2_per_day": WALK_FRACTION * Dw,
            "D0=1_below_walking_factor": Dw, "D0=1_below_5pct_walking_factor": WALK_FRACTION * Dw,
        }
        r["R0_1m"] = {}
        for lab, D0 in [("D0=1", 1.0), ("D0=10", 10.0), ("D0=100", 100.0),
                        ("walking_5pct", WALK_FRACTION * Dw), ("walking", Dw)]:
            val = r0(s, None if r["kernel1m_is_samecell"] else P1, D0)
            r["R0_1m"][lab] = {"D0": D0, "R0": val,
                               "excess_over_wellmixed_pct": 100 * (val / r["R0_wellmixed"] - 1)}
        # same-cell kernel where it differs from the 1 m kernel (classroom, metro carriage); these rows are part
        # of Table s2_model:tab:regime
        if not r["kernel1m_is_samecell"]:
            r["R0_samecell"] = {}
            for lab, D0 in [("D0=1", 1.0), ("D0=10", 10.0), ("D0=100", 100.0),
                            ("walking_5pct", WALK_FRACTION * Dw), ("walking", Dw)]:
                val = r0(s, None, D0)
                r["R0_samecell"][lab] = {"D0": D0, "R0": val,
                                         "excess_over_wellmixed_pct": 100 * (val / r["R0_wellmixed"] - 1)}
        res["scenes"][name] = r

        # ---- table rows -------------------------------------------------------------------------
        title = {"office": "Office", "supermarket": "Supermarket", "classroom": "Classroom",
                 "metro": "Metro carriage"}[name]
        rows_a.append(
            f"{title} & ${s.nx}\\times{s.ny}$ & {a:g} & {n} & {N} & {r['rho_bar']:.3f} & "
            f"{r['d_mean']:.3f} & {r['q_mean']:.3f} & {r['R0_wellmixed']:.3f} & "
            f"{r['R0_frozen_samecell']:.3f} & {r['kernel1m_cells_mean']:.1f} \\\\")
        parts = []
        for z, e in zt.items():
            if e.get("barrier"):
                parts.append(f"{e['name']} ({e['cells']}; barrier)")
            else:
                dstr = f"{e['d']:.2g}" if e["d"] < 0.02 else f"{e['d']:g}" if "D" in zones[z] else f"{e['d']:.2g}"
                parts.append(f"{e['name']} ({e['cells']}; {dstr}; {e['q']:g})")
        rows_b.append(f"{title} & " + ", ".join(parts) + " \\\\")

    # ---- office: what D0 = 1 means ---------------------------------------------------------------
    off = sc.load_scene("office")
    L = th.generator(off).toarray()
    Z = np.where(off.site_zone == "M")[0]
    LZ = L[np.ix_(Z, Z)]
    muZ = -float(np.max(np.linalg.eigvals(LZ).real))
    t_exit_uniform = float(np.linalg.solve(-LZ.T, np.ones(len(Z))).mean())   # mean exit time, uniform start
    W_ = np.where(off.site_zone == "W")[0]
    res["office_regime"] = {
        "D0_m2_per_day": off.a_m ** 2,
        "rms_m_per_day_D0": float(np.sqrt(4.0 * 1.0) * off.a_m),
        "rms_m_per_day_workstation_d=0.3": float(np.sqrt(4.0 * 0.3) * off.a_m),
        "rms_m_per_infectious_period_D0": float(np.sqrt(4.0 / GAMMA) * off.a_m),
        "meeting_room_cells": int(len(Z)),
        "meeting_room_mu_Z_per_day": muZ,
        "meeting_room_1_over_mu_Z_days": 1.0 / muZ,
        "meeting_room_mean_exit_time_uniform_start_days": t_exit_uniform,
        "infectious_period_days": 1.0 / GAMMA,
        "workstation_cells": int(len(W_)),
    }

    # ---- per-cell rule: how often is an infective alone in its cell? --------------------------------
    # the other N-1 walkers are independent and uniform: P(alone) = (1 - 1/n)^(N-1)
    res["alone_probability"] = {
        "office_N100_n234": {"people_per_cell": 100 / off.M, "p_alone": (1 - 1 / off.M) ** 99},
        "uniform_20x20_N100": {"people_per_cell": 100 / 400, "p_alone": (1 - 1 / 400) ** 99},
    }

    # ---- cross-check against plan_theory_checks.json ----------------------------------------------
    plan = json.load(open(OUT / "plan_theory_checks.json"))
    cc = {}
    for name in sc.SCENE_NAMES:
        p = plan["scenes"][name]
        r = res["scenes"][name]
        cc[name] = {
            "q_mean_diff": r["q_mean"] - p["q_mean"],
            "d_mean_diff": r["d_mean"] - p["D_mean_over_D0"],
            "kernel_cells_diff": r["kernel1m_cells_mean"] - p["kernel_cells_mean"],
            "R0_D0=1_diff": r["R0_1m"]["D0=1"]["R0"] - p["R0_1m_D0=1"],
            "R0_D0=10_diff": r["R0_1m"]["D0=10"]["R0"] - p["R0_1m_D0=10"],
            "R0_D0=100_diff": r["R0_1m"]["D0=100"]["R0"] - p["R0_1m_D0=100"],
        }
    cc["office_mu_Z_diff"] = res["office_regime"]["meeting_room_mu_Z_per_day"] - \
        plan["office"]["zone_bound_meeting_room"]["mu_Z"]
    res["crosscheck"] = cc
    worst = max(abs(v) for d in cc.values() for v in (d.values() if isinstance(d, dict) else [d]))
    res["crosscheck_max_abs_diff"] = worst
    assert worst < 1e-8, worst

    # ---- rows of Table s2_model:tab:regime (excess of R0 over beta<q>/gamma, per cent) -----------
    def sci(x):
        e = int(np.floor(np.log10(x)))
        return f"${x / 10 ** e:.1f}\\times10^{{{e}}}$"

    rows_c = []
    for name in sc.SCENE_NAMES:
        r = res["scenes"][name]
        title = {"office": "Office", "supermarket": "Supermarket", "classroom": "Classroom",
                 "metro": "Metro carriage"}[name]
        for kind, t in (("1m", r["R0_1m"]), ("same", r.get("R0_samecell"))):
            if t is None:
                continue
            kern = ("same cell" if r["kernel1m_is_samecell"] else "1\\,m") if kind == "1m" else "same cell"
            first = f"{title} & {kern} & {r['R0_wellmixed']:.3f}" if kind == "1m" else f" & {kern} & "
            rows_c.append(
                f"{first} & {t['D0=1']['R0']:.3f} ({t['D0=1']['excess_over_wellmixed_pct']:.1f}) & "
                f"{t['D0=10']['R0']:.3f} ({t['D0=10']['excess_over_wellmixed_pct']:.2f}) & "
                f"{t['D0=100']['R0']:.3f} ({t['D0=100']['excess_over_wellmixed_pct']:.2f}) & "
                f"{sci(t['walking_5pct']['D0'])} & {t['walking_5pct']['excess_over_wellmixed_pct']:.4f} & "
                f"{sci(t['walking']['D0'])} & {t['walking']['excess_over_wellmixed_pct']:.5f} \\\\")

    with open(OUT / "s2_model_scenes.json", "w") as f:
        json.dump(res, f, indent=1)
    with open(OUT / "s2_model_table_rows.tex", "w") as f:
        f.write("% generated by code/s2_model_scene_table.py -- do not edit\n% panel (a)\n")
        f.write("\n".join(rows_a) + "\n% panel (b): zone (cells; D/D0; q)\n" + "\n".join(rows_b) + "\n")
        f.write("% Table s2_model:tab:regime: scene, kernel, beta<q>/gamma, R0 (excess %) at D0 = 1, 10, 100, "
                "D0 and excess % at 5 % walking, D0 and excess % at continuous walking\n")
        f.write("\n".join(rows_c) + "\n")
    print("\n".join(rows_c))

    # ---- console summary --------------------------------------------------------------------------
    print("\n".join(rows_a))
    print("\n".join(rows_b))
    for name, r in res["scenes"].items():
        print(name, {k: (round(v["D0"], 1), round(v["R0"], 5), round(v["excess_over_wellmixed_pct"], 5))
                     for k, v in r["R0_1m"].items()})
        print("   hop rate mean/min/max per day at D0=1:", r["hop_rate_mean_per_day_D0=1"],
              r["hop_rate_min_per_day_D0=1"], r["hop_rate_max_per_day_D0=1"],
              "| rho_bar", r["rho_bar"], "| people/m2", r["people_per_m2_walkable"])
    print(json.dumps(res["walking"], indent=1))
    print(json.dumps(res["office_regime"], indent=1))
    print(json.dumps(res["alone_probability"], indent=1))
    print("cross-check max |diff| =", worst)


if __name__ == "__main__":
    main()
