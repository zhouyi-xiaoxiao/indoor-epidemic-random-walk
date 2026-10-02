"""plan_theory_checks.py -- recomputation of the mean-field (next-generation) quantities
on the CANONICAL scenes of the article (code/bsc_sim/scenes/*.json), so that the theory
sections and the simulation sections quote numbers for the same floor plans.

The theory analysis (code/bsc_theory) used its own office layouts O1/O2; the simulations
(code/bsc_sim) used scenes/office.json.  This script evaluates the quantities of the theory
analysis on scenes/office.json and adds further checks (part 4 below).

Outputs: ../data/plan_theory_checks.json (all numbers), ../data/plan_theory_sweeps.npz (curves, maps).

Run:  python plan_theory_checks.py          (about 2-4 minutes on one core, < 1 GB)

Conventions: lengths in lattice cells, time in days; D0 in cell^2/day; hop rate between
neighbouring cells = harmonic mean of the two cell values of D (symmetric rule); beta = 0.5/day,
gamma = 0.14/day.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.special as spsp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                       # repository root
sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))
sys.path.insert(0, str(ROOT / "code/bsc_theory/scripts"))

from bsc_sim import scenes as sc             # noqa: E402
from bsc_sim import theory as th             # noqa: E402
import layouts as lay                        # noqa: E402
import ngm_lattice as ng                     # noqa: E402

BETA, GAMMA = 0.5, 0.14
OUT = HERE.parent / "data"
OUT.mkdir(exist_ok=True)
RES: dict = {"beta": BETA, "gamma": GAMMA}
NPZ: dict = {}


def dense_L(scene, D=None):
    return th.generator(scene, D).toarray()


def r0_and_vectors(Ld, q, P=None):
    """R0 = rho(K), K = P^T diag(beta q)(gamma - L)^-1; right (w) and left (z) Perron vectors."""
    n = Ld.shape[0]
    Vinv = np.linalg.inv(GAMMA * np.eye(n) - Ld)
    K = (BETA * q)[:, None] * Vinv
    if P is not None:
        K = P.T @ K
    ev, vl, vr = sla.eig(K, left=True, right=True)
    k = int(np.argmax(ev.real))
    w = np.abs(vr[:, k].real)
    z = np.abs(vl[:, k].real)
    return float(ev[k].real), w / w.sum(), z / z.sum(), K


def r0_only(Ld, q, P=None, gamma=GAMMA):
    n = Ld.shape[0]
    Vinv = np.linalg.inv(gamma * np.eye(n) - Ld)
    K = (BETA * q)[:, None] * Vinv
    if P is not None:
        K = P.T @ K
    return float(np.max(np.linalg.eigvals(K).real))


def growth(Ld, q, P=None):
    """lambda_1 = s(L + F - gamma), Perron vector u (sum 1), spectral gap to the next eigenvalue."""
    F = np.diag(BETA * q)
    if P is not None:
        F = P.T @ F
    J = Ld + F - GAMMA * np.eye(Ld.shape[0])
    ev, vr = np.linalg.eig(J)
    order = np.argsort(-ev.real)
    lam1 = float(ev[order[0]].real)
    gap = float(lam1 - ev[order[1]].real)
    u = np.abs(vr[:, order[0]].real)
    return lam1, u / u.sum(), gap, J


def zone_shares(vec, zones):
    return {z: float(vec[zones == z].sum()) for z in sorted(set(zones.tolist()))}


# ------------------------------------------------------------------------------------------------
# 1. canonical office scene
# ------------------------------------------------------------------------------------------------
def office_block():
    o = {}
    s1 = sc.load_scene("office", D0=1.0)
    n = s1.M
    q = s1.q
    zones = s1.site_zone
    o["cells"] = int(n)
    o["zone_counts"] = {z: int((zones == z).sum()) for z in sorted(set(zones.tolist()))}
    o["zone_counts"]["B"] = int((~s1.access).sum())
    o["q_mean"] = float(q.mean())
    o["q_max"] = float(q.max())
    o["lower_bound"] = BETA * float(q.mean()) / GAMMA
    o["upper_bound"] = BETA * float(q.max()) / GAMMA
    o["lower_bound_qmin"] = BETA * float(q.min()) / GAMMA
    L1 = dense_L(s1)                      # generator at D0 = 1; L(D0) = D0 * L1 (D enters linearly)
    assert np.allclose(L1, L1.T) and np.allclose(L1.sum(axis=0), 0)

    # --- point values at D0 = 1, 10, 100
    for D0 in (1.0, 10.0, 100.0):
        Ld = D0 * L1
        R0, w, z, K = r0_and_vectors(Ld, q)
        lam1, u, gap, J = growth(Ld, q)
        e = z * w / (z @ w)
        o[f"D0={D0:g}"] = dict(
            R0=R0, lambda1=lam1, doubling_time=float(np.log(2) / lam1), gap=gap,
            excess_over_wellmixed_pct=100 * (R0 / o["lower_bound"] - 1),
            elasticity_share=zone_shares(e, zones),
            incidence_share_w=zone_shares(w, zones),
            prevalence_share_u=zone_shares(u, zones),
            area_share={zz: float((zones == zz).mean()) for zz in sorted(set(zones.tolist()))},
            offspring_min=float(K.sum(axis=0).min()), offspring_max=float(K.sum(axis=0).max()),
            offspring_uniform_mean=float(K.sum(axis=0).mean()),
            first_generation_share_uniform_index=zone_shares(K.sum(axis=1) / K.sum(), zones),
        )
        NPZ[f"office_u_D{D0:g}"] = s1.to_grid(u * n)
        NPZ[f"office_w_D{D0:g}"] = s1.to_grid(w * n)
        NPZ[f"office_e_D{D0:g}"] = s1.to_grid(e * n)
        NPZ[f"office_offspring_D{D0:g}"] = s1.to_grid(K.sum(axis=0))
    NPZ["office_q"] = s1.to_grid(q)
    NPZ["office_D"] = s1.to_grid(s1.D)
    NPZ["office_zone_code"] = np.vectorize({"W": 0, "C": 1, "M": 2, "K": 3, "B": 4}.get)(s1.zone).astype(float)

    # --- elasticity of R0 with respect to gamma (central difference), D0 = 1, 10
    for D0 in (1.0, 10.0):
        h = 1e-4
        rp = r0_only(D0 * L1, q, gamma=GAMMA * (1 + h))
        rm = r0_only(D0 * L1, q, gamma=GAMMA * (1 - h))
        o[f"D0={D0:g}"]["elasticity_gamma"] = float((np.log(rp) - np.log(rm)) / (np.log(1 + h) - np.log(1 - h)))

    # --- mobility sweep
    Dgrid = np.logspace(-3, 5, 81)
    sweep = np.array([r0_only(D0 * L1, q) for D0 in Dgrid])
    NPZ["office_Dgrid"] = Dgrid
    NPZ["office_R0_sweep"] = sweep
    table_D = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 2.4e3, 4.87e4]
    o["mobility_table"] = [dict(D0=float(D0), R0=r0_only(D0 * L1, q),
                                excess_pct=100 * (r0_only(D0 * L1, q) / o["lower_bound"] - 1)) for D0 in table_D]

    # --- asymptotic coefficients: R0 ~ lower + beta*C/(<q> D0);  R0 ~ upper*(1 - D0*muZ0/gamma)
    C, qbar = ng.large_D_coefficient(L1, q, pi=np.full(n, 1.0 / n))
    Zmax = np.nonzero(np.isclose(q, q.max()))[0]
    muZ0 = ng.mu_zone(L1, Zmax)
    o["largeD_C"] = float(C)
    o["largeD_coefficient_betaC_over_qbar"] = float(BETA * C / qbar)
    o["smallD_muZ0"] = float(muZ0)
    o["smallD_slope_muZ0_over_gamma"] = float(muZ0 / GAMMA)

    # --- zone bound for the meeting room (Z = cells of zone M), D0 = 1
    ZM = np.nonzero(zones == "M")[0]
    muM = ng.mu_zone(L1, ZM)
    o["zone_bound_meeting_room"] = dict(mu_Z=float(muM), exit_time_days=float(1 / muM),
                                        q_Z=float(q[ZM].min()), bound=float(BETA * q[ZM].min() / (GAMMA + muM)),
                                        cells=int(len(ZM)))
    ZW = np.nonzero(zones == "W")[0]
    muW = ng.mu_zone(L1, ZW)
    o["zone_bound_workstations"] = dict(mu_Z=float(muW), exit_time_days=float(1 / muW),
                                        q_Z=float(q[ZW].min()), bound=float(BETA * q[ZW].min() / (GAMMA + muW)),
                                        cells=int(len(ZW)))

    # --- Galerkin / modal formula (symmetric generator, uniform pi), D0 = 1
    modes = [1, 2, 3, 5, 10, 20, 50, n]
    gal = {}
    for N in modes:
        Rr, diag, mu, phi = ng.ritz_R0(L1, q, BETA, GAMMA, N)
        gal[str(N)] = float(Rr)
    _, diag_all, mu_all, _ = ng.ritz_R0(L1, q, BETA, GAMMA, n)
    o["galerkin_R0_by_modes"] = gal
    o["diagonal_modal_max"] = float(diag_all.max())
    o["diagonal_modal_argmax"] = int(np.argmax(diag_all))
    o["mu_1"] = float(mu_all[1])

    # --- transient: linearised dynamics from a point seed, time to grow 10^4-fold, TV distance to u
    tr = {}
    far = int(s1.site_index[11, 1])          # workstation cell far from the meeting room
    mr = int(ZM[len(ZM) // 2])
    for D0 in (1.0, 10.0, 100.0):
        lam1, u, gap, J = growth(D0 * L1, q)
        ev, V = np.linalg.eigh(0.5 * (J + J.T))
        for name, seed in (("far_workstation", far), ("meeting_room", mr)):
            c0 = V.T @ np.eye(n)[seed]
            ts = np.linspace(0, 80, 1601)
            tot = np.array([(V @ (np.exp(ev * t) * c0)).sum() for t in ts])
            k = int(np.argmax(tot >= 1e4)) if (tot >= 1e4).any() else len(ts) - 1
            It = V @ (np.exp(ev * ts[k]) * c0)
            tv = 0.5 * np.abs(It / It.sum() - u).sum()
            k2 = int(np.argmax(tot >= 1e2))
            rate = float((np.log(tot[k]) - np.log(tot[k2])) / (ts[k] - ts[k2]))
            tr[f"D0={D0:g}|{name}"] = dict(t_1e4=float(ts[k]), tv_to_u=float(tv), fitted_rate_1e2_to_1e4=rate,
                                           lambda1=lam1, gap=gap)
    o["transient_linear"] = tr
    o["transient_seed_cells_yx"] = {"far_workstation": [11, 1], "meeting_room": [int(s1.site_xy[mr][1]), int(s1.site_xy[mr][0])]}

    # --- corridor mobility multiplier (symmetric rule) and the departure-site reading
    corr = {}
    for mult in (1.0, 0.3, 0.1, 0.01):
        D = s1.D.copy()
        D[zones == "C"] *= mult
        corr[f"x{mult:g}"] = r0_only(dense_L(s1, D), q)
    o["corridor_multiplier_symmetric"] = corr
    lat = ng.Lattice(s1.access)
    assert np.array_equal(lat.cells[:, ::-1], s1.site_xy)          # same row-major ordering
    dep = {}
    for mult in (1.0, 0.1, 0.01):
        D = s1.D.copy()
        D[zones == "C"] *= mult
        Md = ng.generator_departure(lat, D).toarray()
        pi = ng.stationary(Md)
        dep[f"x{mult:g}"] = dict(R0=r0_only(Md, q), q_mean_pi=float(pi @ q), lower=BETA * float(pi @ q) / GAMMA)
    o["corridor_multiplier_departure_reading"] = dep

    # --- partitions as movement barriers (symmetric rule): walls on workstation-workstation edges
    rng = np.random.default_rng(20261001)
    Msym = ng.generator_symmetric(lat, s1.D)
    assert np.allclose(Msym.toarray(), L1)
    ww = np.nonzero((zones[lat.edges[:, 0]] == "W") & (zones[lat.edges[:, 1]] == "W"))[0]
    base = o["D0=1"]["R0"]
    part = {"n_WW_edges": int(len(ww))}
    def random_walls(frac):
        """close a fraction of the workstation-workstation edges one at a time, never disconnecting the room"""
        scale = np.ones(len(lat.edges))
        cand = rng.permutation(ww)
        target = int(round(frac * len(cand)))
        closed = 0
        for e_ in cand:
            if closed >= target:
                break
            scale[e_] = 0.0
            if lat.n_components(scale) != 1:
                scale[e_] = 1.0
            else:
                closed += 1
        return scale, closed

    for frac in (0.1, 0.3, 0.5, 1.0):     # 1.0 = as many walls as connectivity allows
        vals, ncl = [], []
        for _ in range(20):
            scale, closed = random_walls(frac)
            vals.append(r0_only(ng.generator_symmetric(lat, s1.D, edge_scale=scale).toarray(), q))
            ncl.append(closed)
        part[f"walls_{int(frac * 100)}pct"] = dict(mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)),
                                                 min=float(np.min(vals)), n_draws=len(vals),
                                                 mean_closed_edges=float(np.mean(ncl)),
                                                 change_pct=100 * (float(np.mean(vals)) / base - 1))
    # obstacles: 39 workstation cells removed (15 % of the 260 lattice cells), room kept connected
    vals = []
    widx = np.nonzero(zones == "W")[0]
    tries = 0
    while len(vals) < 40 and tries < 4000:
        tries += 1
        rm = rng.choice(widx, 39, replace=False)
        acc = s1.access.copy()
        acc[s1.site_xy[rm, 1], s1.site_xy[rm, 0]] = False
        if not sc.is_connected(acc):
            continue
        keep = np.ones(n, bool)
        keep[rm] = False
        s2 = s1.copy_with(access=acc, D_grid=np.where(acc, s1.D_grid, 0.0), q_grid=np.where(acc, s1.q_grid, 0.0))
        vals.append(r0_only(dense_L(s2), s2.q))
    part["obstacles_39_workstations"] = dict(mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)),
                                             min=float(np.min(vals)), n_draws=len(vals), tries=tries,
                                             change_pct=100 * (float(np.mean(vals)) / base - 1))
    o["partitions"] = part

    # --- targeted ventilation: q -> q/2 on a fraction of the cells, exactly equal treated area
    tv = {}
    for D0 in (1.0, 100.0):
        Ld = D0 * L1
        fr = [0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00]
        rows = []
        # adaptive elasticity-guided: treat 6 cells per step (2.6 % of the floor), recompute, stop at k
        qa = q.copy()
        treated = np.zeros(n, bool)
        adaptive = {0: r0_only(Ld, qa)}
        while treated.sum() < n:
            R0, w, z, _ = r0_and_vectors(Ld, qa)
            e = z * w
            e[treated] = -1
            for idx in np.argsort(-e)[: min(6, n - int(treated.sum()))]:
                treated[idx] = True
                qa[idx] = q[idx] / 2
                adaptive[int(treated.sum())] = None
            adaptive[int(treated.sum())] = r0_only(Ld, qa)
        # fill intermediate counts exactly: redo with truncation at each target count
        def adaptive_at(k):
            qa = q.copy()
            treated = np.zeros(n, bool)
            while treated.sum() < k:
                R0, w, z, _ = r0_and_vectors(Ld, qa)
                e = z * w
                e[treated] = -1
                take = min(6, k - int(treated.sum()))
                for idx in np.argsort(-e)[:take]:
                    treated[idx] = True
                    qa[idx] = q[idx] / 2
            return r0_only(Ld, qa)

        R0b, wb, zb, _ = r0_and_vectors(Ld, q)
        eb = zb * wb
        for f in fr:
            k = int(round(f * n))
            # one-shot ranking by the baseline elasticity
            q1 = q.copy()
            q1[np.argsort(-eb)[:k]] /= 2
            # highest q first, ties broken at random (30 draws)
            hq = []
            for _ in range(30):
                key = q + 1e-6 * rng.random(n)
                q2 = q.copy()
                q2[np.argsort(-key)[:k]] /= 2
                hq.append(r0_only(Ld, q2))
            rd = []
            for _ in range(30):
                q3 = q.copy()
                q3[rng.choice(n, k, replace=False)] /= 2
                rd.append(r0_only(Ld, q3))
            rows.append(dict(fraction=f, cells=k, adaptive=adaptive_at(k), one_shot=r0_only(Ld, q1),
                             highest_q_mean=float(np.mean(hq)), highest_q_sd=float(np.std(hq, ddof=1)),
                             random_mean=float(np.mean(rd)), random_sd=float(np.std(rd, ddof=1))))
        tv[f"D0={D0:g}"] = dict(baseline=R0b, rows=rows)
    o["targeted_ventilation"] = tv
    RES["office"] = o
    return s1, L1


# ------------------------------------------------------------------------------------------------
# 2. all four scenes: R0 against mobility (kernel = 1 m contact range, as in the simulations)
# ------------------------------------------------------------------------------------------------
def scenes_block():
    out = {}
    Dgrid = np.logspace(-3, 5, 41)
    NPZ["scenes_Dgrid"] = Dgrid
    for name in sc.SCENE_NAMES:
        s1 = sc.load_scene(name, D0=1.0)
        P1 = th.contact_kernel(s1, th.kernel_radius_cells(s1, 1.0))
        P1d = None if np.allclose(P1.toarray(), np.eye(s1.M)) else P1.toarray()
        L1 = dense_L(s1)
        q = s1.q
        rec = dict(cells=int(s1.M), N=int(s1.N), a_m=float(s1.a_m), q_mean=float(q.mean()), q_max=float(q.max()),
                   q_min=float(q.min()), lower=BETA * float(q.mean()) / GAMMA, upper=BETA * float(q.max()) / GAMMA,
                   lower_qmin=BETA * float(q.min()) / GAMMA, D_mean_over_D0=float(s1.D.mean()),
                   kernel_cells_mean=float(1.0 if P1d is None else (P1d > 0).sum(axis=1).mean()))
        for D0 in (1.0, 10.0, 100.0, 2.4e3):
            rec[f"R0_1m_D0={D0:g}"] = r0_only(D0 * L1, q, P1d)
            rec[f"R0_samecell_D0={D0:g}"] = r0_only(D0 * L1, q)
            rec[f"lambda1_1m_D0={D0:g}"] = growth(D0 * L1, q, P1d)[0]
        # gamma-elasticity with the 1 m kernel
        for D0 in (1.0, 10.0):
            h = 1e-4
            rp = r0_only(D0 * L1, q, P1d, gamma=GAMMA * (1 + h))
            rm = r0_only(D0 * L1, q, P1d, gamma=GAMMA * (1 - h))
            rec[f"elasticity_gamma_1m_D0={D0:g}"] = float((np.log(rp) - np.log(rm)) / (np.log(1 + h) - np.log(1 - h)))
        NPZ[f"{name}_R0_sweep_1m"] = np.array([r0_only(D0 * L1, q, P1d) for D0 in Dgrid])
        NPZ[f"{name}_R0_sweep_samecell"] = np.array([r0_only(D0 * L1, q) for D0 in Dgrid])
        out[name] = rec
        print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
    RES["scenes"] = out


# ------------------------------------------------------------------------------------------------
# 3. alternative office layouts with the same zone fractions (layouts O1, O2 of the theory analysis)
# ------------------------------------------------------------------------------------------------
def layouts_block():
    out = {}
    Dgrid = np.logspace(-3, 5, 81)
    for nm, z in (("O1", lay.office_O1()), ("O2", lay.office_O2())):
        mask, D, q = lay.fields(z, lay.D_OFFICE, lay.Q_OFFICE)
        lat = ng.Lattice(mask)
        L1 = ng.generator_symmetric(lat, lat.vec(D)).toarray()
        qv = lat.vec(q)
        out[nm] = dict(R0_D0_1=r0_only(L1, qv), R0_D0_100=r0_only(100 * L1, qv), q_mean=float(qv.mean()), cells=int(lat.n))
        NPZ[f"layout_{nm}_R0_sweep"] = np.array([r0_only(D0 * L1, qv) for D0 in Dgrid])
    RES["alternative_layouts"] = out


# ------------------------------------------------------------------------------------------------
# 4. further checks: contact-kernel counter-example, pair-level closed forms, uniform q with a kernel
# ------------------------------------------------------------------------------------------------
def checks_block():
    c = {}
    # (i) lower bound beta<q>/gamma fails for a contact kernel: 3 cells in a row, radius-1 kernel
    P = np.array([[0.5, 0.5, 0.0], [1 / 3, 1 / 3, 1 / 3], [0.0, 0.5, 0.5]])
    Lp = np.array([[-1.0, 1.0, 0.0], [1.0, -2.0, 1.0], [0.0, 1.0, -1.0]])
    for qv in ([1.0, 0.0, 1.0], [1.0, 0.2, 1.0]):
        qv = np.array(qv)
        rows = []
        for D in (0.0, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e3):
            K = P.T @ np.diag(BETA * qv) @ np.linalg.inv(GAMMA * np.eye(3) - D * Lp)
            rows.append(dict(D=D, R0_over_beta_gamma=float(np.max(np.linalg.eigvals(K).real)) * GAMMA / BETA))
        c[f"kernel_counterexample_q={qv.tolist()}"] = dict(mean_q=float(qv.mean()), min_q=float(qv.min()), rows=rows)
    # (ii) infinite square lattice: g0 = (2/pi) K(k) / (gamma + 8D), k = 8D/(gamma+8D); scipy ellipk takes m = k^2
    for D in (0.51, 1.0, 10.0):
        k = 8 * D / (GAMMA + 8 * D)
        g0 = (2 / np.pi) * spsp.ellipk(k * k) / (GAMMA + 8 * D)
        # numerical k-space integral of 1/(gamma + 2*mu(k)), mu = 2D(2 - cos kx - cos ky)
        m = 4000
        kk = (np.arange(m) + 0.5) * np.pi / m
        cx = np.cos(kk)
        g_num = float(np.mean(1.0 / (GAMMA + 4 * D * (2 - cx[:, None] - cx[None, :]))))
        rho = 99 / 260
        bq = 0.5
        c[f"g0_infinite_D={D:g}"] = dict(g0_elliptic=float(g0), g0_numeric=g_num,
                                         R1_inf=float((bq / GAMMA) / (1 + bq / rho * g0)), rho_bar=rho, beta_q=bq)
    # (iii) torus closed form against a direct pair solve on a small torus (L = 5, N = 11)
    Lt = 5
    n = Lt * Lt
    A = np.zeros((n, n))
    D = 1.0
    for x in range(Lt):
        for y in range(Lt):
            i = x * Lt + y
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                j = ((x + dx) % Lt) * Lt + (y + dy) % Lt
                A[j, i] += D
    A -= np.diag(A.sum(axis=0))
    Npeople = 11
    rho = (Npeople - 1) / n
    h = 0.5 / rho
    I = np.eye(n)
    L2 = np.kron(A, I) + np.kron(I, A)
    H = np.zeros(n * n)
    H[np.arange(n) * n + np.arange(n)] = h
    u = np.linalg.solve(GAMMA * np.eye(n * n) + np.diag(H) - L2, GAMMA * np.ones(n * n))
    p = 1 - u.reshape(n, n)
    R1_exact = (Npeople - 1) * p.mean(axis=1)
    f = 2 * np.pi * np.arange(Lt) / Lt
    mu = 2 * D * (1 - np.cos(f))[:, None] + 2 * D * (1 - np.cos(f))[None, :]
    g0 = float(np.mean(1 / (GAMMA + 2 * mu)))
    c["torus_L5_N11"] = dict(R1_exact_min=float(R1_exact.min()), R1_exact_max=float(R1_exact.max()),
                             R1_closed=float((0.5 / GAMMA) / (1 + 0.5 / rho * g0)), g0=g0)
    # (iv) uniform q: R0 = beta q/gamma on the office floor plan for any D field and the 1 m kernel of the classroom
    s = sc.load_scene("classroom", D0=1.0)
    P1 = th.contact_kernel(s, th.kernel_radius_cells(s, 1.0)).toarray()
    c["uniform_q_classroom_kernel"] = dict(R0=r0_only(dense_L(s), np.full(s.M, 1.3), P1), expected=BETA * 1.3 / GAMMA)
    RES["checks"] = c


if __name__ == "__main__":
    office_block()
    print("office done", flush=True)
    scenes_block()
    layouts_block()
    checks_block()
    (OUT / "plan_theory_checks.json").write_text(json.dumps(RES, indent=1))
    np.savez_compressed(OUT / "plan_theory_sweeps.npz", **NPZ)
    print(json.dumps({k: v for k, v in RES["office"].items() if k not in ("mobility_table",)}, indent=1)[:6000])
    print(json.dumps(RES["alternative_layouts"], indent=1))
    print(json.dumps(RES["checks"], indent=1))
