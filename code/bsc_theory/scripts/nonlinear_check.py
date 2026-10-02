"""nonlinear_check.py -- checks of the theory against the full NONLINEAR lattice SIR mean-field model

    dS_r/dt = (M S)_r - beta q_r S_r I_r / N_r
    dI_r/dt = (M I)_r + beta q_r S_r I_r / N_r - gamma I_r
    dR_r/dt = (M R)_r + gamma I_r ,        N_r = S_r + I_r + R_r

(A) threshold: with beta rescaled so that R0 = 0.9 / 1.1 the seed dies out / grows (sign equivalence, Thm 1);
(B) heterogeneity raises R0: a room whose well-mixed estimate beta<q>/gamma is below 1 still has an outbreak
    when the NGM R0 is above 1 (Thm 2);
(C) hotspot map: the early spatial profile of I converges to the Perron vector u (Thm 5), and the early growth
    rate equals lambda_1.
Outputs ../results/nonlinear_check.json (+ .log, + arrays for the figures).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).parent))
import layouts as lay
import ngm_lattice as ng

OUT = Path(__file__).resolve().parent.parent / "results"
LOG, RES = [], {}


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def run(M, q, beta, gamma, I0, T, n_out=400, rtol=1e-9, atol=1e-14):
    n = M.shape[0]
    Md = M.toarray()
    Nstar = np.full(n, 1.0)  # uniform stationary occupancy (symmetric rates), 1 person per cell (units irrelevant)

    def rhs(t, y):
        S, I, R = y[:n], y[n:2 * n], y[2 * n:3 * n]
        N = S + I + R
        inf = beta * q * S * I / N
        return np.concatenate([Md @ S - inf, Md @ I + inf - gamma * I, Md @ R + gamma * I, inf])

    # 4th block: cumulative incidence by cell of infection, C_r(t) = int_0^t beta q_r S_r I_r / N_r
    y0 = np.concatenate([Nstar - I0, I0, np.zeros(n), np.zeros(n)])
    ts = np.linspace(0, T, n_out)
    sol = solve_ivp(rhs, (0, T), y0, method="LSODA", t_eval=ts, rtol=rtol, atol=atol)
    S, I, R = sol.y[:n], sol.y[n:2 * n], sol.y[2 * n:3 * n]
    run.cum_incidence = sol.y[3 * n:]
    return ts, S, I, R


def build(z, Dmap, qmap, D0=1.0):
    mask, Dg, qg = lay.fields(z, Dmap, qmap)
    lat = ng.Lattice(mask)
    return lat, ng.generator_symmetric(lat, lat.vec(Dg) * D0), lat.vec(qg)


GAMMA = 0.14
# ------------------------------------------------------------------ (A) threshold
log("== (A) threshold in the nonlinear model (office O1, D0 = 1) ==")
lat, M, q = build(lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE)
R0_unit = ng.R0_dense(M, q, 1.0, GAMMA)  # R0 is linear in beta
A = {}
seed = np.full(lat.n, 1e-6)  # spatially uniform small seed
for target in (0.8, 0.95, 1.05, 1.2, 2.0):
    beta = target / R0_unit
    ts, S, I, R = run(M, q, beta, GAMMA, seed.copy(), T=3000.0, n_out=600)
    tot_I = I.sum(axis=0)
    attack = (R[:, -1].sum() + I[:, -1].sum()) / lat.n
    lam = ng.growth_rate_dense(M, q, beta, GAMMA)[0]
    peak = tot_I.max() / lat.n
    A[target] = dict(beta=beta, lambda1=lam, final_attack=attack, peak_prevalence=peak, grew=bool(peak > 2 * 1e-6),
                     ts=ts.tolist(), prevalence=(tot_I / lat.n).tolist())
    log(f"  R0 = {target:4.2f} (beta = {beta:.5f}, lambda1 = {lam:+.5f}): peak prevalence {peak:.3e}, final attack fraction {attack:.3e}")
RES["A_threshold"] = A

# ------------------------------------------------------------------ (B) heterogeneity raises R0
log("== (B) well-mixed estimate < 1 < NGM R0  (two-zone 'side' office, D0 = 1) ==")
lat2, M2, q2 = build(lay.two_zone("side"), lay.D_TWO, lay.Q_TWO)
qbar = q2.mean()
R0_unit2 = ng.R0_dense(M2, q2, 1.0, GAMMA)
beta = 0.9 * GAMMA / qbar  # well-mixed estimate beta<q>/gamma = 0.9
R0_true = beta * R0_unit2
B = dict(beta=beta, well_mixed=0.9, R0=R0_true)
seed = np.full(lat2.n, 1e-6)
ts, S, I, R = run(M2, q2, beta, GAMMA, seed.copy(), T=6000.0, n_out=1200)
att_het = (R[:, -1].sum() + I[:, -1].sum()) / lat2.n
B["ts"] = ts.tolist()
B["prevalence_heterogeneous"] = (I.sum(axis=0) / lat2.n).tolist()
B["cum_incidence_map"] = lat2.grid(run.cum_incidence[:, -1]).tolist()
zvec = lay.two_zone("side")[lat2.mask]
B["final_attack_heterogeneous"] = att_het
cum = run.cum_incidence[:, -1]
B["share_of_infections_in_workspace"] = float(cum[zvec == "W"].sum() / cum.sum())
B["area_share_workspace"] = float((zvec == "W").mean())
B["infections_per_cell_workspace"] = float(cum[zvec == "W"].mean())
B["infections_per_cell_corridor"] = float(cum[zvec == "C"].mean())
# same room, q replaced by its mean
ts, S, I, R = run(M2, np.full(lat2.n, qbar), beta, GAMMA, seed.copy(), T=6000.0, n_out=1200)
att_hom = (R[:, -1].sum() + I[:, -1].sum()) / lat2.n
B["prevalence_homogenised"] = (I.sum(axis=0) / lat2.n).tolist()
B["final_attack_homogenised"] = att_hom
log(f"  beta<q>/gamma = 0.900, NGM R0 = {R0_true:.4f}: final attack fraction {att_het:.4f}; "
    f"{100 * B['share_of_infections_in_workspace']:.1f}% of infections occur in the workspace ({100 * B['area_share_workspace']:.1f}% of the area); "
    f"same room with q = <q>: attack {att_hom:.2e}")
RES["B_heterogeneity"] = B

# ------------------------------------------------------------------ (C) hotspot profile
log("== (C) early spatial profile vs Perron vector u (office O1) ==")
C = {}
BETA = 0.5
for D0 in (1.0, 10.0, 100.0):
    lat, M, q = build(lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE, D0=D0)
    lam, u, l, _ = ng.growth_rate_dense(M, q, BETA, GAMMA)
    ev = np.sort(np.linalg.eigvals(M.toarray() + np.diag(BETA * q - GAMMA)).real)
    gap = lam - ev[-2]
    C[D0] = dict(lambda1=lam, gap=gap)
    for seed_name, cell in (("far desk", (2, 2)), ("meeting room", (10, 17))):
        I0 = np.zeros(lat.n)
        I0[lat.index[cell]] = 1e-12
        T = 60.0
        ts, S, I, R = run(M, q, BETA, GAMMA, I0, T=T, n_out=601, rtol=1e-10, atol=1e-30)
        tot = I.sum(axis=0)
        prof = I / tot
        tv = 0.5 * np.abs(prof - u[:, None]).sum(axis=0)  # total-variation distance to u
        win = (tot / lat.n > 1e-9) & (tot / lat.n < 1e-4)
        slope = np.polyfit(ts[win], np.log(tot[win]), 1)[0]
        i4 = int(np.argmax(tot / lat.n > 1e-4))
        i2 = int(np.argmax(tot / lat.n > 1e-2))
        C[D0][seed_name] = dict(fitted_growth=slope, t_prev_1e4=float(ts[i4]), tv_at_prev_1e4=float(tv[i4]),
                                t_prev_1e2=float(ts[i2]), tv_at_prev_1e2=float(tv[i2]), tv_initial=float(tv[0]),
                                ts=ts.tolist(), tv=tv.tolist(), prevalence=(tot / lat.n).tolist())
        log(f"  D0 = {D0:5.0f}, seed in {seed_name:12s}: lambda1 = {lam:.4f}, gap {gap:.4f}/day; growth rate fitted for prevalence in "
            f"[1e-9,1e-4]: {slope:.4f}; TV(profile,u): {tv[0]:.3f} at t=0 -> {tv[i4]:.4f} at prevalence 1e-4 (t = {ts[i4]:.1f} d) "
            f"-> {tv[i2]:.4f} at prevalence 1e-2 (t = {ts[i2]:.1f} d)")
    if D0 == 10.0:
        np.save(OUT / "nonlinear_profile_D10.npy", np.stack([lat.grid(prof[:, i4]), lat.grid(u)]))
RES["C_hotspot"] = C


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


(OUT / "nonlinear_check.json").write_text(json.dumps(_clean(RES), indent=1))
(OUT / "nonlinear_check.log").write_text("\n".join(LOG) + "\n")
