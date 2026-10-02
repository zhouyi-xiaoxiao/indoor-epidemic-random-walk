"""worked_examples.py -- worked examples of the next-generation-matrix theory.

Outputs ../results/worked_examples.json (+ .log).  All numbers are computed, none assumed.
Defaults: beta = 0.5 /day, gamma = 0.14 /day, D0 = 1 lattice^2/day, 20 x 13 lattice (a = 1.5 m).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).parent))
import layouts as lay
import ngm_lattice as ng

OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(exist_ok=True)
BETA, GAMMA, D0 = 0.5, 0.14, 1.0
LOG, RES = [], {}
t0 = time.time()


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def build(z, Dmap, qmap, D0=D0, model="symmetric", edge_scale=None):
    mask, Dg, qg = lay.fields(z, Dmap, qmap)
    lat = ng.Lattice(mask)
    D = lat.vec(Dg) * D0
    q = lat.vec(qg)
    if model == "symmetric":
        M = ng.generator_symmetric(lat, D, edge_scale=edge_scale)
        pi = np.full(lat.n, 1.0 / lat.n)
    else:
        M = ng.generator_departure(lat, D, edge_scale=edge_scale)
        pi = (1.0 / D) / (1.0 / D).sum()
    return lat, M, q, pi, D


def final_size(R0):
    return brentq(lambda z: 1 - z - np.exp(-R0 * z), 1e-9, 1 - 1e-15) if R0 > 1 else 0.0


# =============================================================================== E1 two-zone office
log("== E1. Two-zone office (W: D=0.3, q=1.5, ~70%; C: D=2.0, q=0.5, ~30%) ==")
Dsweep = np.logspace(-3, 5, 81)
E1 = {}
for kind in ("side", "aisles", "grid", "scatter"):
    z = lay.two_zone(kind)
    cnt = lay.counts(z)
    rec = dict(counts=cnt)
    for model in ("symmetric", "departure"):
        lat, M, q, pi, D = build(z, lay.D_TWO, lay.Q_TWO, model=model)
        lo, hi = ng.bounds(q, BETA, GAMMA, pi)
        R0 = ng.R0_dense(M, q, BETA, GAMMA)
        lam = ng.growth_rate_dense(M, q, BETA, GAMMA)[0]
        sweep = [ng.R0_dense(M * d, q, BETA, GAMMA) for d in Dsweep]
        rec[model] = dict(R0=R0, R0_over_beta_gamma=R0 * GAMMA / BETA, lower=lo, upper=hi, lambda1=lam,
                          mean_q_pi=float(pi @ q), sweep=sweep)
        C_, qb_ = ng.large_D_coefficient(M, q, pi)
        Zm_ = np.nonzero(q >= q.max() - 1e-12)[0]
        rec[model]["largeD_C"] = C_
        rec[model]["smallD_muZ"] = ng.mu_zone(M, Zm_)
        log(f"  {kind:8s} {model:9s}: <q>_pi = {pi @ q:.4f}; bounds [{lo:.3f}, {hi:.3f}]; R0(D0=1) = {R0:.4f} "
            f"= {R0 * GAMMA / BETA:.4f} beta/gamma; lambda1 = {lam:.4f}/day")
    E1[kind] = rec
E1["Dsweep"] = Dsweep.tolist()
RES["E1_two_zone"] = E1

# =============================================================================== E2 full office
log("== E2. Full office scene (zone counts W156 C39 M26 K13 B26) ==")
E2 = {}
for name, zfun in (("O1", lay.office_O1), ("O2", lay.office_O2)):
    z = zfun()
    rec = dict(counts=lay.counts(z))
    for model in ("symmetric", "departure"):
        lat, M, q, pi, D = build(z, lay.D_OFFICE, lay.Q_OFFICE, model=model)
        assert lat.n_components() == 1
        lo, hi = ng.bounds(q, BETA, GAMMA, pi)
        R0, w, zl, K = ng.R0_dense(M, q, BETA, GAMMA, vectors=True)
        lam, u, l, _ = ng.growth_rate_dense(M, q, BETA, GAMMA)
        sweep = [ng.R0_dense(M * d, q, BETA, GAMMA) for d in Dsweep]
        rec[model] = dict(n=lat.n, R0=R0, lower=lo, upper=hi, lambda1=lam, mean_q_pi=float(pi @ q), sweep=sweep,
                          doubling_time=float(np.log(2) / lam))
        log(f"  {name} {model:9s}: n = {lat.n}; <q>_pi = {pi @ q:.4f}; bounds [{lo:.3f}, {hi:.3f}]; R0(D0=1) = {R0:.4f}; "
            f"lambda1 = {lam:.4f}/day (doubling {np.log(2) / lam:.2f} d)")
        if model == "symmetric":
            ritz = {}
            for nm in (1, 2, 3, 5, 10, 20, 50, 100, lat.n):
                Rr, diag_terms, mu, phi = ng.ritz_R0(M, q, BETA, GAMMA, nm)
                ritz[nm] = Rr
            modal_max = float(diag_terms.max())
            rec["ritz"] = ritz
            rec["diagonal_modal_max"] = modal_max
            rec["mu1"] = float(mu[1])
            log(f"     largest diagonal modal term max_n beta<phi_n|q|phi_n>/(gamma+mu_n) = {modal_max:.4f} "
                f"(attained at n = {int(np.argmax(diag_terms))}); mu_1 = {mu[1]:.5f}/day")
            log("     Galerkin (Ritz) R0 with N slowest modes: " + ", ".join(f"N={k}: {v:.4f}" for k, v in ritz.items()))
            # zone bounds
            zones = {}
            zvec = z[lat.mask]
            for zn in ("M", "W", "K"):
                Z = np.nonzero(zvec == zn)[0]
                muZ = ng.mu_zone(M, Z)
                zones[zn] = dict(mu=muZ, bound=BETA * q[Z].min() / (GAMMA + muZ))
                log(f"     zone bound, Z = all '{zn}' cells: mu_Z = {muZ:.5f}/day -> R0 >= {zones[zn]['bound']:.4f}")
            rec["zone_bounds"] = zones
            C, qbar = ng.large_D_coefficient(M, q, pi)
            Zmax = np.nonzero(q >= q.max() - 1e-12)[0]
            rec["largeD_C"] = C
            rec["smallD_muZ"] = ng.mu_zone(M, Zmax)
            log(f"     large-D expansion: R0 ~ {BETA * qbar / GAMMA:.4f} + {BETA * C / qbar:.4f}/D0;  "
                f"small-D: R0 ~ {hi:.4f} (1 - {rec['smallD_muZ'] / GAMMA:.4f} D0)")
            rec["maps"] = dict(q=lat.grid(q).tolist(), u=lat.grid(u).tolist(), w=lat.grid(w).tolist(),
                               z=lat.grid(zl).tolist(), e=lat.grid(ng.elasticity(w, zl)).tolist(),
                               phi0=lat.grid(np.full(lat.n, 1 / np.sqrt(lat.n))).tolist(), zones=z.tolist())
            # hotspot maps at faster mobility for the figure
            for dd in (10.0, 100.0):
                R0d, wd, zd, _ = ng.R0_dense(M * dd, q, BETA, GAMMA, vectors=True)
                lamd, ud, _, _ = ng.growth_rate_dense(M * dd, q, BETA, GAMMA)
                rec["maps"][f"u_D{int(dd)}"] = lat.grid(ud).tolist()
                rec["maps"][f"e_D{int(dd)}"] = lat.grid(ng.elasticity(wd, zd)).tolist()
                rec[f"R0_D{int(dd)}"] = R0d
            # share of elasticity by zone
            e = ng.elasticity(w, zl)
            rec["elasticity_by_zone"] = {zn: float(e[zvec == zn].sum()) for zn in "WCMK"}
            rec["area_share_by_zone"] = {zn: float((zvec == zn).mean()) for zn in "WCMK"}
            log(f"     elasticity share by zone {rec['elasticity_by_zone']} vs area share {rec['area_share_by_zone']}")
    E2[name] = rec
E2["Dsweep"] = Dsweep.tolist()
RES["E2_office"] = E2

# =============================================================================== E3 partitions
log("== E3. Partitions ==")
rng = np.random.default_rng(42)
z = lay.office_O1()
lat, M, q, pi, D = build(z, lay.D_OFFICE, lay.Q_OFFICE)
zvec = z[lat.mask]
R_base = ng.R0_dense(M, q, BETA, GAMMA)
isWW = (zvec[lat.edges[:, 0]] == "W") & (zvec[lat.edges[:, 1]] == "W")
E3 = dict(R0_base=R_base, n_WW_edges=int(isWW.sum()))
log(f"  baseline O1: R0 = {R_base:.4f}; {isWW.sum()} desk-desk edges")
# (i) partition walls on a random fraction p of desk-desk edges (cells kept; connectivity enforced)
edge_part = {}
def random_walls(p, rng):
    """close a fraction p of the desk-desk edges one at a time, never disconnecting the room"""
    sc = np.ones(len(lat.edges))
    cand = rng.permutation(np.nonzero(isWW)[0])
    target_n = int(round(p * len(cand)))
    closed = 0
    for e_ in cand:
        if closed >= target_n:
            break
        sc[e_] = 0.0
        if lat.n_components(sc) != 1:
            sc[e_] = 1.0
        else:
            closed += 1
    return sc, closed


for p in (0.1, 0.2, 0.3, 0.4, 0.5, 1.0):  # 1.0 = close as many as possible while keeping the room connected
    vals, ncl = [], []
    for _ in range(20):
        sc, closed = random_walls(p, rng)
        Mp = ng.generator_symmetric(lat, D, edge_scale=sc)
        vals.append(ng.R0_dense(Mp, q, BETA, GAMMA))
        ncl.append(closed)
    edge_part[p] = dict(mean=float(np.mean(vals)), sd=float(np.std(vals)), min=float(np.min(vals)), max=float(np.max(vals)),
                        n=len(vals), mean_closed_edges=float(np.mean(ncl)))
    log(f"  (i) walls on {int(p * 100)}% of desk-desk edges ({np.mean(ncl):.0f} edges closed, room kept connected): "
        f"R0 = {np.mean(vals):.4f} +- {np.std(vals):.4f} (min {np.min(vals):.4f}, n = {len(vals)}) -> change {100 * (np.mean(vals) / R_base - 1):+.2f}%")
E3["edge_partitions"] = edge_part
# (iii) 15% more obstacle cells (39 desk cells become obstacles)
vals, qbars = [], []
tries = 0
Widx = np.nonzero(zvec == "W")[0]
while len(vals) < 40 and tries < 4000:
    tries += 1
    rm = rng.choice(Widx, 39, replace=False)
    z2 = z.copy()
    cells = lat.cells[rm]
    z2[cells[:, 0], cells[:, 1]] = "B"
    lat2, M2, q2, pi2, D2 = build(z2, lay.D_OFFICE, lay.Q_OFFICE)
    if lat2.n_components() != 1:
        continue
    vals.append(ng.R0_dense(M2, q2, BETA, GAMMA))
    qbars.append(float(pi2 @ q2))
E3["obstacle_cells_plus15pct"] = dict(mean=float(np.mean(vals)), sd=float(np.std(vals)), min=float(np.min(vals)),
                                      max=float(np.max(vals)), n=len(vals), lower_bound=BETA * float(np.mean(qbars)) / GAMMA)
log(f"  (iii) 39 desk cells -> obstacles (+15% of 260): R0 = {np.mean(vals):.4f} +- {np.std(vals):.4f} "
    f"[{np.min(vals):.4f}, {np.max(vals):.4f}] (n = {len(vals)}); change {100 * (np.mean(vals) / R_base - 1):+.2f}%; "
    f"new lower bound beta<q>/gamma = {BETA * np.mean(qbars) / GAMMA:.4f}")
# (iv) partitions modelled as droplet barriers: they cut q at the desks by a factor c (movement unchanged)
desk_q = {}
for c in (1.0, 0.75, 0.5, 0.25, 0.0):
    q2 = q.copy()
    q2[zvec == "W"] *= c
    desk_q[c] = ng.R0_dense(M, q2, BETA, GAMMA)
E3["desk_q_factor_scan"] = desk_q
log("  (iv) partitions as droplet barriers at desks (q_W -> c q_W): " + ", ".join(f"c={c}: R0={v:.4f}" for c, v in desk_q.items()))
log(f"       acting on desks alone at D0 = 1, the meeting room (q = 1.5) keeps R0 >= {desk_q[0.0]:.3f}")
# same scan at fast mixing (D0 = 100)
lat_f, M_f, q_f, pi_f, D_f = build(z, lay.D_OFFICE, lay.Q_OFFICE, D0=100.0)
R_f = ng.R0_dense(M_f, q_f, BETA, GAMMA)
desk_q_fast = {}
for c in (1.0, 0.75, 0.5, 0.25, 0.0):
    q2 = q_f.copy()
    q2[zvec == "W"] *= c
    desk_q_fast[c] = ng.R0_dense(M_f, q2, BETA, GAMMA)
E3["desk_q_factor_scan_D100"] = desk_q_fast
log("       at D0 = 100 (fast mixing): " + ", ".join(f"c={c}: R0={v:.4f}" for c, v in desk_q_fast.items())
)
# (v) ventilation (q -> 0.5 q everywhere) + 50% capacity + partitions
R_comb = 0.5 * edge_part[0.3]["mean"]
E3["combination_vent_capacity_partitions"] = R_comb
log(f"  (v) ventilation q->0.5q + 50% capacity + partitions (30% edges): R0 = 0.5 x {edge_part[0.3]['mean']:.4f} = {R_comb:.4f} "
    f"(capacity has no effect under frequency-dependent incidence)")
log(f"      masks (beta->0.3beta) + ventilation (q->0.5q): R0 = 0.15 x {R_base:.4f} = {0.15 * R_base:.4f}")
E3["masks_plus_ventilation"] = 0.15 * R_base
RES["E3_partitions"] = E3

# =============================================================================== E4 corridor-D strategy
log("== E4. 'Strategy 2': lower D in the connecting corridors ==")
cgrid = np.logspace(-2, 1, 31)
E4 = dict(c=cgrid.tolist())
for name, z, Dmap, qmap in (("O1", lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE), ("two_zone_side", lay.two_zone("side"), lay.D_TWO, lay.Q_TWO),
                            ("two_zone_aisles", lay.two_zone("aisles"), lay.D_TWO, lay.Q_TWO)):
    rec = {}
    for model in ("symmetric", "departure"):
        vals, los = [], []
        for c in cgrid:
            Dm = dict(Dmap)
            Dm["C"] = Dmap["C"] * c
            lat, M, q, pi, D = build(z, Dm, qmap, model=model)
            vals.append(ng.R0_dense(M, q, BETA, GAMMA))
            los.append(BETA * float(pi @ q) / GAMMA)
        rec[model] = dict(R0=vals, lower=los)
        i1 = int(np.argmin(np.abs(cgrid - 1)))
        i01 = int(np.argmin(np.abs(cgrid - 0.1)))
        log(f"  {name:16s} {model:9s}: R0 at corridor-D x1 = {vals[i1]:.4f}, x0.1 = {vals[i01]:.4f}, x0.01 = {vals[0]:.4f}, x10 = {vals[-1]:.4f}"
            f"   (monotone in c: {'decreasing' if np.all(np.diff(vals) <= 1e-12) else 'increasing' if np.all(np.diff(vals) >= -1e-12) else 'non-monotone'})")
    E4[name] = rec
RES["E4_corridor_D"] = E4

# =============================================================================== E5 targeted ventilation
log("== E5. Targeted ventilation guided by the hotspot (elasticity) map, office O1 ==")
E5 = {}
z = lay.office_O1()
fracs = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]
for dd in (1.0, 100.0):
    lat, M, q, pi, D = build(z, lay.D_OFFICE, lay.Q_OFFICE, D0=dd)
    R0b, w, zl, _ = ng.R0_dense(M, q, BETA, GAMMA, vectors=True)
    e = ng.elasticity(w, zl)
    rec = dict(R0_base=R0b, fracs=fracs, by_elasticity=[], by_q=[], random_mean=[], random_sd=[], greedy=[])
    order_e = np.argsort(-e)
    order_q = np.lexsort((np.random.default_rng(0).random(lat.n), -q))
    rr = np.random.default_rng(5)
    for f in fracs:
        k = int(round(f * lat.n))
        for key, order in (("by_elasticity", order_e), ("by_q", order_q)):
            q2 = q.copy()
            q2[order[:k]] *= 0.5
            rec[key].append(ng.R0_dense(M, q2, BETA, GAMMA))
        rv = []
        for _ in range(30):
            q2 = q.copy()
            q2[rr.choice(lat.n, k, replace=False)] *= 0.5
            rv.append(ng.R0_dense(M, q2, BETA, GAMMA))
        rec["random_mean"].append(float(np.mean(rv)))
        rec["random_sd"].append(float(np.std(rv)))
    # adaptive greedy: treat 5% of cells at a time, recomputing the elasticity map
    q2 = q.copy()
    treated = np.zeros(lat.n, bool)
    step = int(round(0.05 * lat.n))
    gre = {}
    kdone = 0
    while kdone < lat.n:
        R0c, wc, zc, _ = ng.R0_dense(M, q2, BETA, GAMMA, vectors=True)
        ec = ng.elasticity(wc, zc)
        ec[treated] = -1
        pick = np.argsort(-ec)[: min(step, lat.n - kdone)]
        q2[pick] *= 0.5
        treated[pick] = True
        kdone += len(pick)
        gre[kdone / lat.n] = ng.R0_dense(M, q2, BETA, GAMMA)
    rec["greedy_curve"] = [[float(k), float(v)] for k, v in gre.items()]
    E5[f"D0={dd:g}"] = rec
    log(f"  D0 = {dd:g}: baseline R0 = {R0b:.4f}")
    gk = np.array([k for k, _ in rec["greedy_curve"]])
    gv = np.array([v for _, v in rec["greedy_curve"]])
    for i, f in enumerate(fracs):
        gi = int(np.argmin(np.abs(gk - f)))
        log(f"     treat {int(f * 100):3d}% of cells (q -> 0.5q): adaptive elasticity-greedy {gv[gi]:.4f} (at {100 * gk[gi]:.0f}%) | "
            f"one-shot elasticity rank {rec['by_elasticity'][i]:.4f} | q-ranked {rec['by_q'][i]:.4f} | "
            f"random {rec['random_mean'][i]:.4f} +- {rec['random_sd'][i]:.4f}")
RES["E5_targeting"] = E5

# =============================================================================== E6 room size
log("== E6. Room size: reflecting walls vs a door (leak) ==")
E6 = {}
log("  square rooms (D = 0.5, gamma = 0.14, uniform q): ratio R0/R0_inf")
for L in (3, 13, 20, 40):
    lat = ng.Lattice(np.ones((L, L), bool))
    Dv = np.full(lat.n, 0.5)
    M = ng.generator_symmetric(lat, Dv)
    R_refl = ng.R0_sparse(M, np.ones(lat.n), BETA, GAMMA) / (BETA / GAMMA)
    # one door cell in the middle of a wall: a person on the door cell leaves at the hopping rate w = D
    kap = np.zeros(lat.n)
    kap[lat.index[0, L // 2]] = 0.5
    mu1 = -ng.growth_rate_sparse(M, np.zeros(lat.n), 0.0, 0.0, kappa=kap) if L > 3 else ng.mu_leak(M, kap)
    # one whole wall open
    kap2 = np.zeros(lat.n)
    kap2[lat.index[0, :]] = 0.5
    mu2 = -ng.growth_rate_sparse(M, np.zeros(lat.n), 0.0, 0.0, kappa=kap2) if L > 3 else ng.mu_leak(M, kap2)
    E6[L] = dict(reflecting=R_refl,
                 one_door_mu=mu1, one_door_ratio=GAMMA / (GAMMA + mu1), open_wall_mu=mu2, open_wall_ratio=GAMMA / (GAMMA + mu2),
                 open_wall_mu_continuum=0.5 * np.pi**2 / (4 * L**2))
    log(f"   L = {L:2d}: reflecting {R_refl:.6f} | "
        f"one door cell: mu = {mu1:.5f}, ratio {GAMMA / (GAMMA + mu1):.4f} | one open wall: mu = {mu2:.5f} "
        f"(continuum pi^2 D/(4L^2) = {0.5 * np.pi**2 / (4 * L**2):.5f}), ratio {GAMMA / (GAMMA + mu2):.4f}")
# narrow-escape scaling of the one-door eigenvalue
ne = {}
for L in (8, 16, 32, 64, 128):
    lat = ng.Lattice(np.ones((L, L), bool))
    M = ng.generator_symmetric(lat, np.ones(lat.n))
    kap = np.zeros(lat.n)
    kap[lat.index[0, L // 2]] = 1.0
    mu = -ng.growth_rate_sparse(M, np.zeros(lat.n), 0.0, 0.0, kappa=kap)
    kinf = np.zeros(lat.n)
    kinf[lat.index[0, L // 2]] = 1e6
    mu_inf = -ng.growth_rate_sparse(M, np.zeros(lat.n), 0.0, 0.0, kappa=kinf)
    c_fit = L / np.exp(np.pi / (mu_inf * L * L))      # mu_abs = pi D / (L^2 ln(L/c))
    series = 1.0 / (1.0 / mu_inf + L * L / 1.0)          # 1/mu ~ 1/mu_abs + 1/<kappa>_pi
    ne[L] = dict(mu=mu, mu_absorbing=mu_inf, mu_times_L2=mu * L * L, mu_abs_times_L2=mu_inf * L * L,
                 c_fit=c_fit, series_approx=series, series_rel_err=series / mu - 1)
    log(f"   one door, D = 1, L = {L:3d}: absorbing door cell L^2 mu_abs = {mu_inf * L * L:.4f} = pi/ln(L/c), c = {c_fit:.3f}; "
        f"door exit rate kappa = w: mu = {mu:.4e}; series formula (1/mu_abs + 1/<kappa>_pi)^-1 = {series:.4e} "
        f"(rel. err {100 * (series / mu - 1):+.2f}%)")
E6["narrow_escape"] = ne
RES["E6_room_size"] = E6

# uniform leak = dwell time: R_visit = beta q/(gamma + 1/T)
log("  uniform exit rate 1/T (dwell time T): R = beta <q>/(gamma + 1/T); dwell times and mean q of the scenes")
dw = {}
for scene, T_min, qbar in (("office", 8 * 60, 1.2), ("supermarket", 27, 1.2), ("classroom", 60, 1.5), ("metro", 20, 2.5)):
    T = T_min / 1440.0
    Rv = BETA * qbar / (GAMMA + 1 / T)
    dw[scene] = dict(T_min=T_min, qbar=qbar, R_per_visit=Rv, closed_room=BETA * qbar / GAMMA)
    log(f"   {scene:11s}: T = {T_min:4d} min, q = {qbar}: R per visit = {Rv:.4f}; closed-room value beta q/gamma = {BETA * qbar / GAMMA:.3f}")
RES["E6_dwell"] = dw

# =============================================================================== E7 call centre
log("== E7. Uniform q = 1.68 (call-centre example) ==")
R_cc = BETA * 1.68 / GAMMA
E7 = dict(R0=R_cc, final_size_R0=final_size(R_cc))
# confirm numerically on a lattice with D = 0.2 and 10% random obstacles
rngc = np.random.default_rng(3)
while True:
    mask = rngc.random((13, 20)) > 0.1
    latc = ng.Lattice(mask)
    if latc.n_components() == 1:
        break
Mc = ng.generator_symmetric(latc, np.full(latc.n, 0.2))
E7["R0_lattice_check"] = ng.R0_dense(Mc, np.full(latc.n, 1.68), BETA, GAMMA)
log(f"  R0 = beta q/gamma = {R_cc:.4f} for every layout and D (lattice check with D=0.2, 10% obstacles: {E7['R0_lattice_check']:.6f}); "
    f"SIR final size at R0=6: {E7['final_size_R0']:.4f}")
RES["E7_call_centre"] = E7

# =============================================================================== E9 physical mobility
log("== E9. How large is D? (Green-Kubo: D = v^2 tau / 2 in 2-D) ==")
E9 = {}
lat, M, q, pi, D = build(lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE)
lo, hi = ng.bounds(q, BETA, GAMMA, pi)
for label, Dphys in (("default D0 = 1 lattice^2/day", 1.0), ("D0 = 100", 100.0),
                     ("walking 5% of the time: 0.05 * v^2 tau/2 (v=1.3 m/s, tau=1.5 s)", 0.05 * 1.3**2 * 1.5 / 2 * 86400 / 1.5**2),
                     ("continuous walking: v^2 tau/2", 1.3**2 * 1.5 / 2 * 86400 / 1.5**2)):
    R = ng.R0_dense(M * Dphys, q, BETA, GAMMA)
    E9[label] = dict(D0_lattice2_per_day=Dphys, R0=R, excess_over_mean_pct=100 * (R / lo - 1), position=(R - lo) / (hi - lo))
    log(f"  {label}: D0 = {Dphys:.3g} lattice^2/day -> R0 = {R:.4f} ({100 * (R / lo - 1):+.3f}% above beta<q>/gamma = {lo:.4f})")
RES["E9_physical_D"] = E9

RES["runtime_s"] = time.time() - t0
log(f"runtime {RES['runtime_s']:.1f} s")


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


(OUT / "worked_examples.json").write_text(json.dumps(_clean(RES), indent=1))
(OUT / "worked_examples.log").write_text("\n".join(LOG) + "\n")
