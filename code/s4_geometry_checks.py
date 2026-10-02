"""s4_geometry_checks.py -- numbers and re-checks for Section 4 ("Room size, leaks and hotspot maps").

Everything quoted in sections/m4_geometry.tex and sections/supp_geometry.tex that is not read directly from
    data/plan_theory_checks.json                                   (canonical office scene), or
    data/bsc_theory/worked_examples.json               (room size, call centre), or
    data/bsc_theory/fig2_room_size_data.json           (curves of Fig. room-size/leak)
is computed here and written to data/s4_geometry_checks.json.  The script also re-derives, with its own
generator code, the room-size and leak numbers of those files (block A) and checks the statements of the
section on the canonical office scene (block C).

    python s4_geometry_checks.py          # about 20 s, < 1 GB

Blocks
  A  square rooms with uniform q (own implementation): closed room, one door cell, one open side,
     corridor with absorbing ends (closed form, continuum ratio, critical length);
     one-door (narrow-escape) fit.
  B  closed room with uniform q and a general contact kernel, random barriers and mobility field.
  C  canonical office scene: identities of the hotspot theorem (w ~ Q xi, elasticity formula against
     finite differences, e = q*l*xi normalised), the slow-mixing bound on the elasticity outside
     argmax q used in the proof of part (iv), and derived numbers quoted in the text.
  D  derived numbers from data/plan_theory_checks.json (ratios, differences) quoted in the text/tables.

Conventions: lengths in cells, time in days, beta = 0.5/day, gamma = 0.14/day; column index = source cell.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                                   # repository root
DATA = HERE.parent / "data"
BETA, GAMMA = 0.5, 0.14
OUT: dict = {"beta": BETA, "gamma": GAMMA}


# ------------------------------------------------------------------------------------------------
# lattice helpers (own implementation, independent of the research directories)
# ------------------------------------------------------------------------------------------------
def grid_generator(mask: np.ndarray, D: np.ndarray) -> sp.csr_matrix:
    """Symmetric (harmonic-mean) generator on the walkable cells of ``mask``; M[y, x] = hop rate x -> y."""
    ny, nx = mask.shape
    idx = -np.ones(mask.shape, int)
    idx[mask] = np.arange(mask.sum())
    rows, cols, vals = [], [], []
    for (dy, dx) in ((0, 1), (1, 0)):
        a = idx[: ny - dy, : nx - dx].ravel()
        b = idx[dy:, dx:].ravel()
        ok = (a >= 0) & (b >= 0)
        a, b = a[ok], b[ok]
        Da, Db = D.ravel()[np.flatnonzero(mask.ravel())][a], D.ravel()[np.flatnonzero(mask.ravel())][b]
        w = 2 * Da * Db / (Da + Db)
        rows += [a, b]
        cols += [b, a]
        vals += [w, w]
    n = int(mask.sum())
    A = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()
    return (A - sp.diags(np.asarray(A.sum(axis=0)).ravel())).tocsr(), idx


def top_eig_sym(A: sp.spmatrix) -> float:
    """Largest eigenvalue of a symmetric matrix (dense for small sizes, shift-invert Lanczos otherwise)."""
    n = A.shape[0]
    if n <= 1700:
        return float(np.linalg.eigvalsh(A.toarray())[-1])
    return float(spla.eigsh(A.tocsc(), k=1, sigma=1e-9, which="LM", return_eigenvectors=False)[0])


def connected(mask: np.ndarray) -> bool:
    ys, xs = np.nonzero(mask)
    seen = np.zeros(mask.shape, bool)
    stack = [(ys[0], xs[0])]
    seen[ys[0], xs[0]] = True
    while stack:
        y, x = stack.pop()
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            yy, xx = y + dy, x + dx
            if 0 <= yy < mask.shape[0] and 0 <= xx < mask.shape[1] and mask[yy, xx] and not seen[yy, xx]:
                seen[yy, xx] = True
                stack.append((yy, xx))
    return bool(seen[mask].all())


# ------------------------------------------------------------------------------------------------
# A. room size and leaks, uniform q
# ------------------------------------------------------------------------------------------------
def block_A():
    A = {}
    D = 0.5
    rooms = {}
    for L in (3, 13, 20, 40):
        mask = np.ones((L, L), bool)
        M, idx = grid_generator(mask, np.full((L, L), D))
        n = L * L
        # closed room: rho(beta q (gamma - M)^-1) = beta q / (gamma - s(M)); s(M) computed, not assumed
        sM = top_eig_sym(M)
        closed_ratio = GAMMA / (GAMMA - sM)
        # one door cell in the middle of a wall row, exit rate = hop rate D
        kap = np.zeros(n)
        kap[idx[0, L // 2]] = D
        mu_door = -top_eig_sym(M - sp.diags(kap))
        # one whole side open
        kap2 = np.zeros(n)
        kap2[idx[0, :]] = D
        mu_side = -top_eig_sym(M - sp.diags(kap2))
        rooms[str(L)] = dict(
            closed_ratio=closed_ratio, s_M=sM,
            door_mu=mu_door, door_ratio=GAMMA / (GAMMA + mu_door), door_residence_days=1 / mu_door,
            door_kappa_mean=D / n,
            side_mu=mu_side, side_ratio=GAMMA / (GAMMA + mu_side), side_residence_days=1 / mu_side,
            side_kappa_mean=D * L / n,
        )
    A["square_rooms_D0.5"] = rooms
    # mu_kappa is proportional to the mobility scale when kappa = D at the door cells (M_kappa = D * const)
    res = [r[k] for r in rooms.values() for k in ("door_residence_days", "side_residence_days")]
    A["residence_days_range_D0.5"] = [min(res), max(res)]
    A["max_mu_D0.5"] = max(max(r["door_mu"], r["side_mu"]) for r in rooms.values())
    A["residence_range_at_D2400"] = dict(factor=2.4e3 / D, min_minutes=min(res) / (2.4e3 / D) * 1440,
                                         max_days=max(res) / (2.4e3 / D))
    A["beta_q_over_gamma_q1"] = BETA / GAMMA

    # corridor of L cells, hop rate w = 1, a hop beyond either end leaves the room
    cor = {}
    for L in (3, 4, 5, 13, 100):
        T = sp.diags([np.full(L - 1, 1.0), np.full(L, -2.0), np.full(L - 1, 1.0)], [-1, 0, 1])
        mu_num = -float(np.linalg.eigvalsh(T.toarray())[-1])
        mu_closed = 2 * (1 - np.cos(np.pi / (L + 1)))
        cor[str(L)] = dict(mu_numeric=mu_num, mu_closed_form=mu_closed, mu_continuum=np.pi**2 / L**2,
                           ratio_to_continuum=mu_closed / (np.pi**2 / L**2),
                           R_room=BETA / (GAMMA + mu_closed), R_continuum=BETA / (GAMMA + np.pi**2 / L**2))
    A["corridor_w1_q1"] = cor
    # critical length: R_room = 1  <=>  mu = beta q - gamma
    A["corridor_critical_L_lattice"] = float(np.pi / np.arccos(1 - (BETA - GAMMA) / 2) - 1)
    A["corridor_critical_L_continuum"] = float(np.pi / np.sqrt(BETA - GAMMA))

    # one-door square room, D = 1: absorbing door cell and door with exit rate w = D
    ne = {}
    for L in (8, 16, 32, 64, 128):
        mask = np.ones((L, L), bool)
        M, idx = grid_generator(mask, np.ones((L, L)))
        n = L * L
        door = idx[0, L // 2]
        keep = np.setdiff1d(np.arange(n), [door])
        mu_abs = -top_eig_sym(M[keep][:, keep])              # principal submatrix: door cell absorbing
        kap = np.zeros(n)
        kap[door] = 1.0
        mu_k = -top_eig_sym(M - sp.diags(kap))
        series = 1 / (1 / mu_abs + n / 1.0)                  # 1/mu ~ 1/mu_abs + 1/<kappa>_pi, <kappa>_pi = 1/n
        ne[str(L)] = dict(mu_abs=mu_abs, mu_abs_L2=mu_abs * n, c_fit=L / np.exp(np.pi / (mu_abs * n)),
                          mu_kappa=mu_k, series=series, series_rel_err_pct=100 * (series / mu_k - 1))
    A["one_door_D1"] = ne
    OUT["A_room_size_and_leaks"] = A

    # compare with the theory part
    ref = json.loads((ROOT / "data/bsc_theory/worked_examples.json").read_text())["E6_room_size"]
    dev = 0.0
    for L in ("3", "13", "20", "40"):
        dev = max(dev, abs(rooms[L]["door_ratio"] - ref[L]["one_door_ratio"]),
                  abs(rooms[L]["side_ratio"] - ref[L]["open_wall_ratio"]),
                  abs(rooms[L]["closed_ratio"] - ref[L]["reflecting"]))
    for L in ("8", "16", "32", "64", "128"):
        dev = max(dev, abs(ne[L]["c_fit"] - ref["narrow_escape"][L]["c_fit"]),
                  abs(ne[L]["series_rel_err_pct"] / 100 - ref["narrow_escape"][L]["series_rel_err"]))
    A["max_abs_deviation_from_worked_examples_json"] = dev
    assert dev < 2e-4, dev


# ------------------------------------------------------------------------------------------------
# B. closed room, uniform q, general kernel / barriers / mobility field
# ------------------------------------------------------------------------------------------------
def block_B():
    rng = np.random.default_rng(404)
    worst = 0.0
    worst_lam = 0.0
    worst_off = 0.0
    n_inst = 0
    while n_inst < 40:
        ny, nx = rng.integers(4, 10, size=2)
        mask = rng.random((ny, nx)) > 0.15
        if mask.sum() < 6 or not connected(mask):
            continue
        Dfield = 10 ** rng.uniform(-2, 2, size=(ny, nx))
        M, _ = grid_generator(mask, Dfield)
        M = M.toarray()
        n = M.shape[0]
        if n_inst % 2:                                        # non-symmetric irreducible generator as well
            W = np.where(M > 0, M * 10 ** rng.uniform(-1, 1, size=M.shape), 0.0)
            M = W - np.diag(W.sum(axis=0))
        P = rng.random((n, n)) * (rng.random((n, n)) < 0.3) + np.eye(n) * rng.random(n)
        P[P.sum(axis=1) == 0, 0] = 1.0
        P = P / P.sum(axis=1, keepdims=True)                  # row-stochastic contact kernel
        q0 = rng.uniform(0.3, 3.0)
        F = BETA * P.T * q0                                   # F = beta P^T Q with Q = q0 * identity
        K = F @ np.linalg.inv(GAMMA * np.eye(n) - M)
        R0 = float(np.max(np.abs(np.linalg.eigvals(K))))
        lam = float(np.max(np.linalg.eigvals(M + F - GAMMA * np.eye(n)).real))
        worst = max(worst, abs(R0 / (BETA * q0 / GAMMA) - 1))
        worst_lam = max(worst_lam, abs(lam - (BETA * q0 - GAMMA)))
        worst_off = max(worst_off, float(np.max(np.abs(K.sum(axis=0) / (BETA * q0 / GAMMA) - 1))))
        n_inst += 1
    OUT["B_closed_room_general_kernel"] = dict(instances=n_inst, max_rel_err_R0=worst, max_abs_err_lambda1=worst_lam,
                                               max_rel_err_offspring=worst_off)
    assert worst < 1e-10 and worst_lam < 1e-10


# ------------------------------------------------------------------------------------------------
# C. canonical office scene: hotspot identities and the slow-mixing bound
# ------------------------------------------------------------------------------------------------
def perron(A, left=False):
    ev, vec = np.linalg.eig(A.T if left else A)
    k = int(np.argmax(ev.real))
    v = np.abs(vec[:, k].real)
    return float(ev[k].real), v / v.sum()


def block_C():
    sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))
    from bsc_sim import scenes as sc
    from bsc_sim import theory as th
    s1 = sc.load_scene("office", D0=1.0)
    L1 = th.generator(s1).toarray()
    q = s1.q
    n = s1.M
    zones = s1.site_zone
    C = {}
    Z = np.isclose(q, q.max())
    c_rate = float(np.max(-np.diag(L1)))                      # max_x |M0_xx|
    C["max_exit_rate_M0"] = c_rate
    for D0 in (1e-3, 1e-2, 1e-1, 1.0, 100.0):
        M = D0 * L1
        V = GAMMA * np.eye(n) - M
        K = (BETA * q)[:, None] * np.linalg.inv(V)
        R0, w = perron(K)
        _, z = perron(K, left=True)
        e = z * w / (z @ w)
        # xi: beta Q xi = R0 (gamma - M) xi ; l: left null vector of A = M + beta Q / R0 - gamma
        Amat = M + np.diag(BETA * q / R0) - GAMMA * np.eye(n)
        sA, xi = perron(Amat)
        _, ell = perron(Amat, left=True)
        p = ell * xi / (ell @ xi)
        e2 = q * p / (q @ p)
        wq = q * xi / (q @ xi)
        rec = dict(R0=R0, s_A=sA, max_abs_e_minus_qlxi=float(np.abs(e - e2).max()),
                   max_abs_w_minus_Qxi=float(np.abs(w - wq).max()), sum_e=float(e.sum()),
                   e_outside_argmax_q=float(e[~Z].sum()))
        # slow-mixing bound: p_x <= D0 * c / (gamma (1 - q_x/q_max)) for x outside argmax q
        bound = D0 * c_rate / (GAMMA * (1 - q[~Z] / q.max()))
        rec["bound_p_holds"] = bool(np.all(p[~Z] <= bound + 1e-12))
        rec["max_p_over_bound"] = float(np.max(p[~Z] / bound))
        if D0 == 1.0:
            # finite-difference check of the elasticity on 12 cells
            rng = np.random.default_rng(7)
            errs = []
            for x in rng.choice(n, 12, replace=False):
                h = 1e-5
                qp, qm = q.copy(), q.copy()
                qp[x] *= 1 + h
                qm[x] *= 1 - h
                Rp = perron((BETA * qp)[:, None] * np.linalg.inv(V))[0]
                Rm = perron((BETA * qm)[:, None] * np.linalg.inv(V))[0]
                errs.append(abs((np.log(Rp) - np.log(Rm)) / (np.log(1 + h) - np.log(1 - h)) - e[x]))
            rec["elasticity_fd_max_abs_err"] = float(max(errs))
            # first generation of a uniformly placed index case: K pi = beta q pi / gamma
            rec["first_gen_identity_err"] = float(np.abs(K @ np.full(n, 1 / n) - BETA * q / (GAMMA * n)).max())
        C[f"D0={D0:g}"] = rec
        assert rec["bound_p_holds"] and rec["max_abs_e_minus_qlxi"] < 1e-8
    # --- why the targeting strategies differ at D0 = 1 (zone composition of the treated cells and of the
    #     elasticity after treatment); same rules as code/plan_theory_checks.py
    def zshare(vec):
        return {zz: float(vec[zones == zz].sum()) for zz in ("W", "C", "M", "K")}

    def zcount(sel):
        return {zz: int((zones[sel] == zz).sum()) for zz in ("W", "C", "M", "K")}

    def elast(M, qv):
        K = (BETA * qv)[:, None] * np.linalg.inv(GAMMA * np.eye(n) - M)
        R0, w = perron(K)
        _, z = perron(K, left=True)
        return R0, z * w / (z @ w)

    M1 = L1
    R0b, eb = elast(M1, q)
    T = {}
    # (i) recomputed-elasticity rule, 6 cells per step, stopped at exactly k cells
    for frac in (0.10, 0.50):
        k = int(round(frac * n))
        qa = q.copy()
        treated = np.zeros(n, bool)
        while treated.sum() < k:
            _, e = elast(M1, qa)
            e = e.copy()
            e[treated] = -1
            for i in np.argsort(-e)[: min(6, k - int(treated.sum()))]:
                treated[i] = True
                qa[i] = q[i] / 2
        R, e_after = elast(M1, qa)
        T[f"recomputed_{int(100 * frac)}pct"] = dict(R0=R, treated_by_zone=zcount(treated), elasticity_share_after=zshare(e_after))
    # (ii) baseline ranking, not recomputed
    for frac in (0.20, 0.50):
        k = int(round(frac * n))
        sel = np.argsort(-eb)[:k]
        q1 = q.copy()
        q1[sel] /= 2
        R, e_after = elast(M1, q1)
        T[f"not_recomputed_{int(100 * frac)}pct"] = dict(R0=R, treated_by_zone=zcount(sel), elasticity_share_after=zshare(e_after))
    # (iii) highest q first at 75 %: all meeting-room cells and 150 of the 156 workstation cells
    rng = np.random.default_rng(11)
    Rs, shares = [], []
    for _ in range(10):
        key = q + 1e-6 * rng.random(n)
        sel = np.argsort(-key)[: int(round(0.75 * n))]
        q2 = q.copy()
        q2[sel] /= 2
        R, e_after = elast(M1, q2)
        Rs.append(R)
        shares.append(zshare(e_after))
    T["highest_q_75pct"] = dict(R0_mean=float(np.mean(Rs)), treated_by_zone=zcount(sel),
                                elasticity_share_after_mean={zz: float(np.mean([s_[zz] for s_ in shares])) for zz in ("W", "C", "M", "K")})
    # zone bound contributed by the untreated kitchen (Proposition zone)
    ZK = np.flatnonzero(zones == "K")
    muK = -float(np.max(np.linalg.eigvalsh(M1[np.ix_(ZK, ZK)])))
    T["kitchen_zone"] = dict(cells=int(len(ZK)), q=float(q[ZK].min()), mu_Z=muK, bound=BETA * float(q[ZK].min()) / (GAMMA + muK))
    C["targeting_explanations_D0=1"] = T

    # number of lattice edges that join the meeting room (the kitchen) to the rest of the floor
    for zz, name in (("M", "meeting_room"), ("K", "kitchen")):
        inside = zones == zz
        C[f"{name}_boundary_edges"] = int(((L1 > 0) & inside[None, :] & ~inside[:, None]).sum())
    # limits of the meeting-room shares as D0 -> infinity (Theorem hotspots (iv))
    C["meeting_room_area_share"] = float((zones == "M").mean())
    C["meeting_room_share_limit_e_and_w"] = float(q[zones == "M"].sum() / q.sum())
    OUT["C_office_hotspot_identities"] = C


# ------------------------------------------------------------------------------------------------
# D. derived numbers from data/plan_theory_checks.json
# ------------------------------------------------------------------------------------------------
def block_D():
    P = json.loads((DATA / "plan_theory_checks.json").read_text())
    off = P["office"]
    d1 = off["D0=1"]
    D = {}
    D["mixing_time_1_over_gap_days"] = {k: 1 / off[k]["gap"] for k in ("D0=1", "D0=10", "D0=100")}
    D["linear_phase_ln100_over_lambda1_days"] = {k: float(np.log(100) / off[k]["lambda1"]) for k in ("D0=1", "D0=10", "D0=100")}
    zw = off["zone_bound_workstations"]
    D["workstation_zone_growth_rate"] = BETA * zw["q_Z"] - GAMMA - zw["mu_Z"]
    zm = off["zone_bound_meeting_room"]
    D["meeting_zone_growth_rate"] = BETA * zm["q_Z"] - GAMMA - zm["mu_Z"]
    D["uniform_q_half"] = d1["R0"] / 2
    D["uniform_mask0.3_and_q_half"] = d1["R0"] * 0.15
    D["partition_values"] = {k: dict(mean=v["mean"], sd=v["sd"], min=v["min"], closed=v.get("mean_closed_edges"),
                                     change_pct=v["change_pct"], above_baseline=bool(v["min"] > d1["R0"]))
                             for k, v in off["partitions"].items() if isinstance(v, dict)}
    cs = off["corridor_multiplier_symmetric"]
    D["corridor_change_pct"] = {k: 100 * (v / cs["x1"] - 1) for k, v in cs.items()}
    tv = {}
    for key in ("D0=1", "D0=100"):
        rows = off["targeted_ventilation"][key]["rows"]
        tv[key] = {f"{int(round(100 * r['fraction']))}%": dict(
            cells=r["cells"], recomputed=r["adaptive"], highest_q=r["highest_q_mean"], highest_q_sd=r["highest_q_sd"],
            random=r["random_mean"], random_sd=r["random_sd"], not_recomputed=r["one_shot"],
            gain_recomputed_vs_random=r["random_mean"] - r["adaptive"],
            spread=max(r["adaptive"], r["highest_q_mean"], r["random_mean"], r["one_shot"])
            - min(r["adaptive"], r["highest_q_mean"], r["random_mean"], r["one_shot"])) for r in rows}
    D["targeted_ventilation"] = tv
    D["mobility_excess_pct"] = {f"{r['D0']:g}": r["excess_pct"] for r in off["mobility_table"]}
    cc = json.loads((ROOT / "data/bsc_theory/worked_examples.json").read_text())["E7_call_centre"]
    D["call_centre"] = dict(R0=cc["R0"], lattice_check=cc["R0_lattice_check"])
    OUT["D_derived"] = D


if __name__ == "__main__":
    block_A()
    block_B()
    block_C()
    block_D()
    (DATA / "s4_geometry_checks.json").write_text(json.dumps(OUT, indent=1))
    print(json.dumps(OUT, indent=1))
