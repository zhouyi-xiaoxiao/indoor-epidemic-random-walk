#!/usr/bin/env python
"""appA_numerics_checks.py -- numbers of Supplementary Section S11 (numerical verification) that are not stored verbatim
in a results file of the research directories.

Run (from anywhere):
    python appA_numerics_checks.py            # all parts
    python appA_numerics_checks.py derived    # one part: derived | adversarial | simulator

Output: ../data/appA_numerics_checks.json   (one top-level key per part; parts are merged, so a part can be
re-run without repeating the others).

Parts
  derived      z-scores and chi-square p-values computed from the stored results of the research directories
               (data/bsc_sim/01_validation.json, code/bsc_sim/verify/v02_simcheck.json), and the
               counts of random instances per generator family used in the theorem table
               (data/bsc_theory/verify_theorems.json, T2_positions.npy).  No new simulation.
  adversarial  optimiser-driven searches for a counter-example to (i) the lower bound
               beta<q>_pi/gamma <= R0 and (ii) the monotone decrease of R0 in the mobility scale, for
               non-reversible dense generators.  Every optimised instance is re-evaluated in 60-digit
               arithmetic (mpmath) with a generator whose column sums are exactly zero.  This repeats, with
               other seeds and with one wider range of rates, the searches of the re-check
               (code/bsc_theory/verify/v1b_adv.py, v1d_adv3_exact.py).
  simulator    re-run of the re-check's second simulator (code/bsc_sim/verify/vsim.c through
               vlib.py; shares no code with the simulator of the article) on (a) the one-cell model, whose
               final-size law is known exactly, 4 x 10^6 + 2 x 10^7 runs with each of the two simulators, and (b) the exposure identity in the office
               scene (uniform and Perron-distributed index case).

Seeds are fixed; the script is deterministic.  Chinese edition: this script produces numbers only (no figure).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2dist

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                       # repository root
OUT = HERE.parent / "data" / "appA_numerics_checks.json"
THEORY = ROOT / "code" / "bsc_theory"
SIM = ROOT / "code" / "bsc_sim"


def save(key, value):
    merged = json.loads(OUT.read_text()) if OUT.exists() else {}
    merged[key] = value
    OUT.write_text(json.dumps(merged, indent=1, default=float))


# ====================================================================================================
# part 1: derived statistics
# ====================================================================================================
def part_derived():
    out = {}
    val = json.loads((SIM / "data" / "01_validation.json").read_text())
    z95 = 1.96

    # --- test A: C core vs pure-Python reference
    A = []
    for r in val["A_c_vs_python"]:
        A.append(dict(case=r["case"], z_final=r["z_final"], z_offspring=r["z_offspring"], chi2=r["chi2"],
                      dof=r["dof"], p_chi2=float(chi2dist.sf(r["chi2"], r["dof"])),
                      ms_per_run_python=1e3 * r["sec_per_run_python"], ms_per_run_c=1e3 * r["sec_per_run_c"],
                      speed_ratio=r["sec_per_run_python"] / r["sec_per_run_c"]))
    out["A"] = dict(rows=A, max_abs_z=max(max(abs(a["z_final"]), abs(a["z_offspring"])) for a in A),
                    min_p_chi2=min(a["p_chi2"] for a in A),
                    speed_ratio_range=[min(a["speed_ratio"] for a in A), max(a["speed_ratio"] for a in A)])

    # --- test B: exact final-size law of the one-cell model
    B = val["B_wellmixed_exact"]
    out["B"] = dict(z_mean=(B["sim_mean"] - B["exact_mean"]) / B["sim_se"], chi2=B["chi2"], dof=B["dof"],
                    p_chi2=float(chi2dist.sf(B["chi2"], B["dof"])))

    # --- test C: exposure identity (z-scores from the stored 95 % intervals)
    C = []
    for r in val["C_exposure_identity"]:
        for tag, th, m, lo, hi in (("perron", r["R0_ngm"], r["exposure_perron_mean"], r["lo"], r["hi"]),
                                   ("uniform", r["wellmixed"], r["exposure_uniform_mean"], r["lo_u"], r["hi_u"]),
                                   ("cell_max", r["cell_max_theory"], r["cell_max_mean"], r["cell_max_lo"], r["cell_max_hi"]),
                                   ("cell_min", r["cell_min_theory"], r["cell_min_mean"], r["cell_min_lo"], r["cell_min_hi"])):
            se = (hi - lo) / (2 * z95)
            C.append(dict(scene=r["scene"], index=tag, theory=th, sim=m, se=se, z=(m - th) / se,
                          inside_95=bool(lo <= th <= hi)))
    out["C"] = dict(rows=C, n=len(C), max_abs_z=max(abs(c["z"]) for c in C),
                    n_outside_95=sum(not c["inside_95"] for c in C),
                    chi2_sum_z2=float(sum(c["z"] ** 2 for c in C)),
                    p_sum_z2=float(chi2dist.sf(sum(c["z"] ** 2 for c in C), len(C))))

    # --- test D: pair theory vs counted secondary cases
    D = []
    for r in val["D_pair_theory"]:
        se = (r["hi"] - r["lo"]) / (2 * z95)
        D.append(dict(scene=r["scene"], D0=r["D0"], kernel_cells=r["kernel_cells"], theory=r["pair_theory"],
                      sim=r["sim_mean"], se=se, z=(r["sim_mean"] - r["pair_theory"]) / se,
                      rel_diff=r["sim_mean"] / r["pair_theory"] - 1,
                      inside_95=bool(r["lo"] <= r["pair_theory"] <= r["hi"]),
                      ratio_to_meanfield=r["pair_theory"] / r["meanfield"]))
    out["D"] = dict(rows=D, n=len(D), max_abs_z=max(abs(d["z"]) for d in D),
                    max_abs_rel_diff=max(abs(d["rel_diff"]) for d in D),
                    n_outside_95=sum(not d["inside_95"] for d in D),
                    chi2_sum_z2=float(sum(d["z"] ** 2 for d in D)),
                    p_sum_z2=float(chi2dist.sf(sum(d["z"] ** 2 for d in D), len(D))))

    # --- test G: schedules
    G = val["G_schedules"]
    se = (G["onoff_hi"] - G["onoff_lo"]) / (2 * z95)
    out["G"] = dict(z_final=G["z_final"], z_offspring=G["z_offspring"],
                    onoff_z=(G["onoff_exposure_sim"] - G["onoff_exposure_theory"]) / se, onoff_se=se)

    # --- the re-check's stored checks of the second simulator
    v = json.loads((SIM / "verify" / "v02_simcheck.json").read_text())
    V = {}
    w = v["wellmixed"]
    V["wellmixed_40k"] = dict(z_mean=(w["sim"] - w["exact"]) / w["se"], chi2=w["chi2"], dof=w["dof"],
                              p_chi2=float(chi2dist.sf(w["chi2"], w["dof"])))
    for k in ("small_rad0.0", "small_rad1.0"):
        V[k] = dict(z=(v[k]["sim"] - v[k]["R1"]) / v[k]["se"])
    for k in ("small_expo0.0", "small_expo1.0"):
        V[k] = dict(z=(v[k]["sim"] - v[k]["theory"]) / v[k]["se"])
    out["check_v02"] = V
    v4 = json.loads((SIM / "verify" / "v04_finite.json").read_text().replace("NaN", "null"))
    zs = {k: (r["sim"] - r["exact"]) / r["se"] for k, r in v4.items() if r.get("exact") is not None and "sim" in r}
    out["check_v04"] = dict(z=zs, n=len(zs), max_abs_z=max(abs(x) for x in zs.values()))

    # --- theorem table: instances per generator family
    th = json.loads((THEORY / "results" / "verify_theorems.json").read_text())
    pos = np.load(THEORY / "results" / "T2_positions.npy")
    out["theorems"] = dict(T1_kinds=th["T1"]["kinds"], T2_nonconstant_q=int(len(pos)),
                           T2_nonconstant_by_family=[int((pos[:, 0] == k).sum()) for k in range(4)],
                           T2_max_position=float(pos[:, 1].max()), T2_min_position=float(pos[:, 1].min()),
                           corridor_max_abs_err=max(abs(r["mu"] - r["lattice_exact"]) for r in th["T4d_1d"].values()),
                           size_scan_max_abs_dev=max(abs(r["R0"] - 0.5 / 0.14) for r in th["T4a_size_scan"].values()),
                           runtime_s=th["runtime_s"])
    # --- software environment and stored timings (for the reproducibility paragraph)
    import platform
    import subprocess

    import matplotlib
    import mpmath
    import scipy
    try:
        cc = subprocess.run(["cc", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
        cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True).stdout.strip()
    except Exception:
        cc, cpu = "unknown", "unknown"
    out["environment"] = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                              matplotlib=matplotlib.__version__, mpmath=mpmath.__version__, cc=cc, cpu=cpu)
    bench = json.loads((SIM / "data" / "10_benchmark.json").read_text())["D0=1"]
    out["timing"] = dict(c_core_ms_per_office_epidemic=1e3 * bench["c_core_cpu_per_run_s"],
                         events_per_office_epidemic=bench["events_per_run"],
                         python_reference_s_per_office_epidemic=bench["python_reference_cpu_per_run_s"],
                         python_reference_ms_per_event=1e3 * bench["python_reference_cpu_per_event_s"],
                         verify_theorems_runtime_s=th["runtime_s"])
    save("derived", out)
    print(json.dumps({k: v for k, v in out.items() if k in ("environment", "timing", "theorems")}, indent=1, default=float))


# ====================================================================================================
# part 2: adversarial searches, re-evaluated in exact-conservative 60-digit arithmetic
# ====================================================================================================
def part_adversarial():
    import mpmath as mp
    mp.mp.dps = 60
    GAM = "0.3"

    def R0_float(M, q, g=0.3):
        n = len(q)
        K = q[:, None] * np.linalg.inv(g * np.eye(n) - M)
        return float(np.max(np.abs(np.linalg.eigvals(K))))

    def stat_solve(M):
        n = M.shape[0]
        A = M.copy()
        A[-1, :] = 1.0
        r = np.zeros(n)
        r[-1] = 1.0
        return np.linalg.solve(A, r)

    def unpack(z, n, clip):
        z = np.clip(z, -clip, clip)
        W = np.exp(z[:n * n]).reshape(n, n)
        np.fill_diagonal(W, 0.0)
        q = np.exp(z[n * n:n * n + n])
        return W, q, z

    def exact_generator(W):
        """mp matrix with the float off-diagonal rates of W and a diagonal equal to minus the EXACT column sums."""
        n = W.shape[0]
        M = mp.matrix(n, n)
        for j in range(n):
            s = mp.mpf(0)
            for i in range(n):
                if i != j:
                    M[i, j] = mp.mpf(float(W[i, j]))
                    s += M[i, j]
            M[j, j] = -s
        return M

    def R0_exact(Mm, q, D=1):
        n = len(q)
        V = mp.mpf(GAM) * mp.eye(n) - mp.mpf(D) * Mm
        Vi = V ** -1
        K = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                K[i, j] = mp.mpf(float(q[i])) * Vi[i, j]
        return max(abs(e) for e in mp.eig(K, left=False, right=False))

    def stat_exact(Mm):
        n = Mm.rows
        A = Mm.copy()
        for j in range(n):
            A[n - 1, j] = 1
        return mp.lu_solve(A, mp.matrix([0] * (n - 1) + [1]))

    out = {}
    t0 = time.time()
    # ---- (i) lower bound  R0 / (<q>_pi / gamma) - 1 >= 0   (beta = 1)
    rng = np.random.default_rng(20261002)
    lower = {}
    for clip, nstart in ((12.0, 120),):
        worst_f, worst_e, rows = np.inf, mp.mpf("inf"), []
        for trial in range(nstart):
            n = int(rng.integers(3, 6))

            def f(z, n=n, clip=clip):
                W, q, _ = unpack(z, n, clip)
                M = W - np.diag(W.sum(axis=0))
                pi = stat_solve(M)
                return R0_float(M, q) / ((pi @ q) / 0.3) - 1.0
            res = minimize(f, rng.normal(0, 1, n * n + n), method="Nelder-Mead",
                           options=dict(maxiter=3000, xatol=1e-7, fatol=1e-15))
            W, q, _ = unpack(res.x, n, clip)
            Mm = exact_generator(W)
            pim = stat_exact(Mm)
            lo = sum(pim[i] * mp.mpf(float(q[i])) for i in range(n)) / mp.mpf(GAM)
            ex = R0_exact(Mm, q) / lo - 1
            worst_f = min(worst_f, res.fun)
            worst_e = min(worst_e, ex)
            rows.append((float(res.fun), float(ex), float(np.ptp(q) / q.max())))
        rows = np.array(rows)
        lower[f"clip{clip:g}"] = dict(n_starts=nstart, min_float64=float(worst_f), min_exact=mp.nstr(worst_e, 6),
                                      min_exact_float=float(worst_e), n_exact_negative=int((rows[:, 1] < 0).sum()),
                                      median_exact=float(np.median(rows[:, 1])),
                                      n_float64_negative=int((rows[:, 0] < 0).sum()))
        print("lower bound", lower[f"clip{clip:g}"], f"[{time.time() - t0:.0f}s]", flush=True)
    out["lower_bound"] = lower

    # ---- (ii) monotonicity  R0(D2) - R0(D1) <= 0 for D2 > D1   (beta = 1)
    rng = np.random.default_rng(20261003)
    mono = {}
    for clip in (4.0, 8.0, 15.0):
        best_e, best_f, rows, nfail = mp.mpf("-inf"), -np.inf, [], 0
        for trial in range(40):
            n = int(rng.integers(3, 6))

            def f(z, n=n, clip=clip):
                W, q, zz = unpack(z, n, clip)
                M0 = W - np.diag(W.sum(axis=0))
                D1 = np.exp(zz[-2])
                D2 = D1 * (1 + np.exp(zz[-1]))
                return -(R0_float(D2 * M0, q) - R0_float(D1 * M0, q))
            res = minimize(f, rng.normal(0, 1.5, n * n + n + 2), method="Nelder-Mead",
                           options=dict(maxiter=2500, xatol=1e-7, fatol=1e-15))
            W, q, zz = unpack(res.x, n, clip)
            D1 = float(np.exp(zz[-2]))
            D2 = D1 * (1 + float(np.exp(zz[-1])))
            Mm = exact_generator(W)
            try:
                r1, r2 = R0_exact(Mm, q, D1), R0_exact(Mm, q, D2)
            except Exception as e:            # mp.eig did not converge: retry with more digits
                mp.mp.dps = 120
                try:
                    r1, r2 = R0_exact(Mm, q, D1), R0_exact(Mm, q, D2)
                except Exception:
                    nfail += 1
                    mp.mp.dps = 60
                    continue
                mp.mp.dps = 60
            ex = r2 - r1
            if ex > best_e:
                best_e = ex
            best_f = max(best_f, -res.fun)
            rows.append((float(-res.fun), float(ex), float(ex / r1), D1, D2))
        rows = np.array(rows)
        mono[f"clip{clip:g}"] = dict(n_starts=40, n_evaluated=int(len(rows)), n_eig_failures=nfail,
                                     max_float64=float(best_f), max_exact=mp.nstr(best_e, 6),
                                     max_exact_float=float(best_e), max_exact_relative=float(rows[:, 2].max()),
                                     n_exact_positive=int((rows[:, 1] > 0).sum()),
                                     n_float64_positive=int((rows[:, 0] > 0).sum()),
                                     rate_range=[float(np.exp(-clip)), float(np.exp(clip))])
        print("monotonicity", f"clip{clip:g}", mono[f"clip{clip:g}"], f"[{time.time() - t0:.0f}s]", flush=True)
    out["monotonicity"] = mono
    out["digits"] = 60
    out["runtime_s"] = time.time() - t0
    save("adversarial", out)


# ====================================================================================================
# part 3: the re-check's second simulator
# ====================================================================================================
def exact_final_size(N, lam, gamma):
    """Final-size law of the Markovian SIR epidemic with pair hazard lam and one initial infective
    (dynamic programming on the embedded jump chain; own implementation)."""
    prob = np.zeros((N + 1, N + 2))
    prob[N - 1, 1] = 1.0
    out = np.zeros(N + 1)
    for s in range(N - 1, -1, -1):
        for i in range(N - s, 0, -1):
            p = prob[s, i]
            if p == 0.0:
                continue
            pinf = lam * s / (lam * s + gamma)
            if s > 0:
                prob[s - 1, i + 1] += p * pinf
            if i == 1:
                out[N - s] += p * (1 - pinf)
            else:
                prob[s, i - 1] += p * (1 - pinf)
    return out


def chi2_against(pex, counts, n):
    E = pex * n
    chi, dof, ce, co = 0.0, -1, 0.0, 0.0
    for k in range(1, len(pex)):
        ce += E[k]
        co += counts[k]
        if ce >= 5:
            chi += (co - ce) ** 2 / ce
            dof += 1
            ce = co = 0.0
    return float(chi), int(dof)


def part_simulator():
    sys.path.insert(0, str(SIM / "verify"))
    import vlib as V                       # the re-check's code; imports nothing from the article's simulator
    B, G = V.BETA, V.GAMMA
    out = {}
    t0 = time.time()
    # ---- (a) one cell, N = 40: exact final-size law; both simulators.
    #      Stage 1 (planned): 20 batches of 200 000 runs.  Stage 2 (replication, added after stage 1 showed a
    #      2.6-standard-error deviation of the mean for the article's simulator): 50 batches of 400 000 runs
    #      with other seeds.  Both stages are reported; nothing is discarded.
    sys.path.insert(0, str(SIM / "src"))
    from bsc_sim.scenes import uniform_room            # the article's simulator (C core gillespie.c)
    from bsc_sim.sim import simulate
    from scipy.stats import kstest
    N = 40
    pex = exact_final_size(N, B / (N - 1), G)
    kk = np.arange(N + 1)
    exact_mean = float(kk @ pex)
    exact_sd = float(np.sqrt(kk ** 2 @ pex - exact_mean ** 2))
    room_v = V.uniform(1, 1, N)
    room_a = uniform_room(1, 1, N, D=1.0, q=1.0)
    stages = dict(stage1=dict(nb=20, n=200000,
                              second=lambda b: V.sim(room_v, 200000, (b + 1) * 1000003)["final"],
                              article=lambda b: simulate(room_a, B, G, n_rep=200000, seed=7001 + b, t_max=2000.0,
                                                         dt_out=1000.0).final_size),
                  stage2=dict(nb=50, n=400000,
                              second=lambda b: V.sim(room_v, 400000, (b + 101) * 1000003)["final"],
                              article=lambda b: simulate(room_a, B, G, n_rep=400000, seed=90000 + b, t_max=2000.0,
                                                         dt_out=1000.0).final_size))

    def summary(counts):
        ntot = int(counts.sum())
        chi, dof = chi2_against(pex, counts, ntot)
        m = float(kk @ counts / ntot)
        sd = float(np.sqrt(kk ** 2 @ counts / ntot - m ** 2))
        pm, pe = float(counts[:4].sum() / ntot), float(pex[:4].sum())
        return dict(n=ntot, mean=m, se=sd / np.sqrt(ntot), z=(m - exact_mean) / (sd / np.sqrt(ntot)), chi2=chi,
                    dof=dof, p=float(chi2dist.sf(chi, dof)), p_minor_le3=pm, exact_p_minor_le3=pe,
                    z_minor=(pm - pe) / np.sqrt(pe * (1 - pe) / ntot))

    wm = dict(N=N, exact_mean=exact_mean, exact_sd=exact_sd)
    for name in ("article", "second"):
        total = np.zeros(N + 1)
        wm[name] = {}
        for stage, cfg in stages.items():
            rows, pooled = [], np.zeros(N + 1)
            for bidx in range(cfg["nb"]):
                fs = np.asarray(cfg[name](bidx))
                counts = np.bincount(fs, minlength=N + 1)
                pooled += counts
                r = summary(counts)
                rows.append(dict(batch=bidx, mean=r["mean"], z=r["z"], chi2=r["chi2"], p=r["p"]))
            total += pooled
            ps = np.array([r["p"] for r in rows])
            wm[name][stage] = dict(n_batches=cfg["nb"], n_per_batch=cfg["n"], pooled=summary(pooled),
                                   batch_mean_range=[min(r["mean"] for r in rows), max(r["mean"] for r in rows)],
                                   batch_p_range=[float(ps.min()), float(ps.max())],
                                   n_batches_p_below_0p05=int((ps < 0.05).sum()),
                                   ks_p_uniformity_of_batch_p=float(kstest(ps, "uniform").pvalue),
                                   sum_batch_chi2=float(sum(r["chi2"] for r in rows)), sum_batch_dof=39 * len(rows),
                                   p_sum_batch_chi2=float(chi2dist.sf(sum(r["chi2"] for r in rows), 39 * len(rows))),
                                   batches=rows)
            print(name, stage, {k: v for k, v in wm[name][stage].items() if k != "batches"},
                  f"[{time.time() - t0:.0f}s]", flush=True)
        wm[name]["combined"] = summary(total)
        print(name, "combined", wm[name]["combined"], flush=True)
    out["wellmixed"] = wm
    save("simulator", out)

    # ---- (b) exposure identity, office scene, D0 = 1, same-cell kernel (a = 1.5 m > 1 m)
    office = V.load("office", D0=1.0)
    K = V.ngm(office)
    w_, v_ = np.linalg.eig(K)
    k = int(np.argmax(w_.real))
    R0 = float(w_[k].real)
    perron = np.abs(v_[:, k].real)
    perron /= perron.sum()
    wellmixed = float(B * office.q.mean() / G)
    res = {"R0": R0, "wellmixed": wellmixed, "cells": int(office.M), "N": int(office.N)}
    n_u, n_p = 400000, 200000
    s = V.sim(office, n_u, 20261004, gmax=0, move_all=True)
    m, se = V.mci(s["expo"])
    res["uniform_index"] = dict(n=n_u, theory=wellmixed, sim=m, se=se, z=(m - wellmixed) / se)
    print("office uniform", res["uniform_index"], f"[{time.time() - t0:.0f}s]", flush=True)
    rng = np.random.default_rng(20261005)
    idx = rng.choice(office.M, size=n_p, p=perron).astype(np.int32)
    s = V.sim(office, n_p, 20261006, gmax=0, move_all=True, index_pos=idx)
    m, se = V.mci(s["expo"])
    res["perron_index"] = dict(n=n_p, theory=R0, sim=m, se=se, z=(m - R0) / se)
    print("office perron", res["perron_index"], f"[{time.time() - t0:.0f}s]", flush=True)
    out["office_exposure"] = res
    out["runtime_s"] = time.time() - t0
    save("simulator", out)


if __name__ == "__main__":
    parts = sys.argv[1:] or ["derived", "adversarial", "simulator"]
    for p in parts:
        {"derived": part_derived, "adversarial": part_adversarial, "simulator": part_simulator}[p]()
