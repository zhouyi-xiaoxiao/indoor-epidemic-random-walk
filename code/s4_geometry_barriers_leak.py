"""s4_geometry_barriers_leak.py -- two checks of Section 4.

(A) Partitions on the canonical office scene at three mobility scales (D0 = 1, 10, 100 cell^2/day).
    Theorem "edges" (edge-wise monotonicity) covers operations that KEEP the cell set: lowering edge rates
    and closing edges (walls).  Turning cells into barrier cells REMOVES cells, and with them their q, from
    the walkable set; the theorem does not apply and the sign of the change of R0 depends on the mobility.
      walls         : a fraction of the 242 workstation-workstation edges is closed (room kept connected);
      barrier cells : 39 workstation cells (15 % of the 260 lattice cells) are removed (room kept connected).
    The random draws are those of plan_theory_checks.py (same seed, same order of calls; that script
    evaluated D0 = 1 only), re-evaluated at the three mobility scales (the generator is linear in D0); the
    D0 = 1 means are asserted to reproduce data/plan_theory_checks.json.

(B) Leaky room with heterogeneous q (Theorem "leaky", part (4)): the substitution gamma -> gamma + mu_leak
    is exact for uniform q or for a uniform leak, and is NOT exact for heterogeneous q with a localised
    door.  Office scene with one door cell at the end of a corridor arm; compared quantities:
      R_room    = rho(beta Q (gamma - M + diag kappa)^-1)                       (exact)
      R_subst   = rho(beta Q ((gamma + mu_leak) - M)^-1)                         (substitution)
      bounds    = beta <q>_nu / (gamma + mu_leak) <= R_room <= beta q_max / (gamma + mu_leak)
    and, as a control, a uniform leak kappa = 1/T, for which R_room = R_subst to rounding.

    python s4_geometry_barriers_leak.py        # about 1 minute, < 1 GB
Output: ../data/s4_geometry_barriers_leak.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                                   # repository root
sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))
sys.path.insert(0, str(ROOT / "code/bsc_theory/scripts"))

from bsc_sim import scenes as sc                         # noqa: E402
from bsc_sim import theory as th                         # noqa: E402
import ngm_lattice as ng                                 # noqa: E402

BETA, GAMMA = 0.5, 0.14
DATA = HERE.parent / "data"
D0S = (1.0, 10.0, 100.0)


def r0(L, q, gamma=GAMMA, kappa=None):
    """rho(beta Q (gamma - L + diag kappa)^-1) for a symmetric generator L (dense)."""
    n = L.shape[0]
    V = gamma * np.eye(n) - L
    if kappa is not None:
        V = V + np.diag(kappa)
    K = (BETA * q)[:, None] * np.linalg.inv(V)
    return float(np.max(np.linalg.eigvals(K).real))


def block_A():
    """Same random draws as plan_theory_checks.py (seed 20261001, same order of calls), evaluated at three
    mobility scales.  The D0 = 1 values are asserted to reproduce data/plan_theory_checks.json."""
    plan = json.loads((DATA / "plan_theory_checks.json").read_text())["office"]["partitions"]
    s1 = sc.load_scene("office", D0=1.0)
    n, q, zones = s1.M, s1.q, s1.site_zone
    L1 = th.generator(s1).toarray()
    assert np.allclose(L1, L1.T) and np.allclose(L1.sum(axis=0), 0)
    lat = ng.Lattice(s1.access)
    assert np.allclose(ng.generator_symmetric(lat, s1.D).toarray(), L1)
    ww = np.nonzero((zones[lat.edges[:, 0]] == "W") & (zones[lat.edges[:, 1]] == "W"))[0]
    base = {D0: r0(D0 * L1, q) for D0 in D0S}
    out = {"n_cells": int(n), "n_WW_edges": int(len(ww)), "baseline_R0": {f"D0={D0:g}": base[D0] for D0 in D0S},
           "wellmixed_baseline": BETA * float(q.mean()) / GAMMA, "seed": 20261001,
           "draws": "identical to code/plan_theory_checks.py (office/partitions)"}
    rng = np.random.default_rng(20261001)

    def stats(vals, D0):
        v = np.array(vals)
        return dict(mean=float(v.mean()), sd=float(v.std(ddof=1)), min=float(v.min()), max=float(v.max()),
                    change_pct=100 * (float(v.mean()) / base[D0] - 1),
                    min_change_pct=100 * (float(v.min()) / base[D0] - 1),
                    max_change_pct=100 * (float(v.max()) / base[D0] - 1),
                    n_below_baseline=int((v < base[D0] - 1e-12).sum()),
                    n_above_baseline=int((v > base[D0] + 1e-12).sum()))

    # ---- walls: close edges, keep the cell set (covered by the edge-wise monotonicity theorem)
    def random_walls(frac):
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

    walls = {}
    for frac in (0.1, 0.3, 0.5, 1.0):
        vals = {D0: [] for D0 in D0S}
        ncl = []
        for _ in range(20):
            scale, closed = random_walls(frac)
            Lw = ng.generator_symmetric(lat, s1.D, edge_scale=scale).toarray()
            ncl.append(closed)
            for D0 in D0S:
                vals[D0].append(r0(D0 * Lw, q))
        key = f"walls_{int(frac * 100)}pct"
        rec = {"n_draws": 20, "mean_closed_edges": float(np.mean(ncl))}
        for D0 in D0S:
            rec[f"D0={D0:g}"] = stats(vals[D0], D0)
        assert abs(rec["D0=1"]["mean"] - plan[key]["mean"]) < 1e-9, key          # same draws as plan_theory_checks.py
        assert abs(rec["mean_closed_edges"] - plan[key]["mean_closed_edges"]) < 1e-12
        walls[key] = rec
    out["walls"] = walls
    out["walls_all_draws_above_baseline"] = {f"D0={D0:g}": bool(all(walls[k][f"D0={D0:g}"]["n_below_baseline"] == 0
                                                                     for k in walls)) for D0 in D0S}

    # ---- barrier cells: 39 workstation cells removed (changes the cell set and the mean of q)
    widx = np.nonzero(zones == "W")[0]
    vals = {D0: [] for D0 in D0S}
    wm = []
    tries = 0
    while len(wm) < 40 and tries < 4000:
        tries += 1
        rm = rng.choice(widx, 39, replace=False)
        acc = s1.access.copy()
        acc[s1.site_xy[rm, 1], s1.site_xy[rm, 0]] = False
        if not sc.is_connected(acc):
            continue
        s2 = s1.copy_with(access=acc, D_grid=np.where(acc, s1.D_grid, 0.0), q_grid=np.where(acc, s1.q_grid, 0.0))
        L2 = th.generator(s2).toarray()
        wm.append(BETA * float(s2.q.mean()) / GAMMA)
        for D0 in D0S:
            vals[D0].append(r0(D0 * L2, s2.q))
    rec = {"n_draws": len(wm), "tries": tries, "cells_removed": 39, "q_removed": 1.4,
           "q_mean_before": float(q.mean()), "q_mean_after": float(np.mean(wm)) * GAMMA / BETA,
           "wellmixed_after": float(np.mean(wm)), "wellmixed_after_spread": float(np.ptp(wm)),
           "wellmixed_change_pct": 100 * (float(np.mean(wm)) / out["wellmixed_baseline"] - 1)}
    for D0 in D0S:
        rec[f"D0={D0:g}"] = stats(vals[D0], D0)
    po = plan["obstacles_39_workstations"]
    assert abs(rec["D0=1"]["mean"] - po["mean"]) < 1e-9 and rec["tries"] == po["tries"]
    out["barrier_cells_39_workstations"] = rec
    return out


def block_B():
    s1 = sc.load_scene("office", D0=1.0)
    n, q, zones = s1.M, s1.q, s1.site_zone
    L1 = th.generator(s1).toarray()
    # door cell: the corridor cell on the lattice boundary that is farthest from the meeting room
    cor = np.flatnonzero(zones == "C")
    xy = s1.site_xy
    on_edge = (xy[cor, 0] == 0) | (xy[cor, 0] == s1.nx - 1) | (xy[cor, 1] == 0) | (xy[cor, 1] == s1.ny - 1)
    cand = cor[on_edge]
    mroom = xy[zones == "M"].mean(axis=0)
    door = int(cand[np.argmax(np.abs(xy[cand] - mroom).sum(axis=1))])
    out = {"door_cell_xy": [int(xy[door, 0]), int(xy[door, 1])], "door_zone": str(zones[door]), "cases": {}}
    for D0, kap in ((1.0, 1.0), (10.0, 100.0), (10.0, 1.0), (1.0, 100.0)):
        L = D0 * L1
        kappa = np.zeros(n)
        kappa[door] = kap
        Mk = L - np.diag(kappa)
        ev, vec = np.linalg.eigh(Mk)                       # symmetric: left = right Perron vector
        mu = -float(ev[-1])
        u = np.abs(vec[:, -1])
        nu = u * u / (u @ u)
        R_room = r0(L, q, kappa=kappa)
        R_sub = r0(L, q, gamma=GAMMA + mu)
        lo = BETA * float(nu @ q) / (GAMMA + mu)
        hi = BETA * float(q.max()) / (GAMMA + mu)
        out["cases"][f"D0={D0:g}|kappa={kap:g}"] = dict(
            mu_leak=mu, kappa_mean_pi=kap / n, R_closed=r0(L, q), R_room=R_room, R_substitution=R_sub,
            substitution_minus_exact=R_sub - R_room, substitution_rel_err_pct=100 * (R_sub / R_room - 1),
            q_mean_nu=float(nu @ q), lower_bound=lo, upper_bound=hi,
            bounds_hold=bool(lo <= R_room + 1e-12 and R_room <= hi + 1e-12),
            substitution_below_lower_bound=bool(R_sub < lo))
        assert out["cases"][f"D0={D0:g}|kappa={kap:g}"]["bounds_hold"]
    # control: uniform leak kappa = 1/T (Corollary "cases" (a)): substitution exact
    T = 8.0 / 24.0
    ctrl = {}
    for D0 in (1.0, 10.0):
        L = D0 * L1
        ctrl[f"D0={D0:g}"] = dict(R_room=r0(L, q, kappa=np.full(n, 1 / T)), R_substitution=r0(L, q, gamma=GAMMA + 1 / T))
        assert abs(ctrl[f"D0={D0:g}"]["R_room"] - ctrl[f"D0={D0:g}"]["R_substitution"]) < 1e-10
    out["uniform_leak_T_8h"] = ctrl
    # random check of part (4) on non-reversible instances (left and right Perron vectors differ)
    rng = np.random.default_rng(5)
    worst_lo = worst_hi = 0.0
    n_sub_below = 0
    for _ in range(1500):
        m = int(rng.integers(3, 9))
        W = rng.lognormal(0, 1, size=(m, m)) * (rng.random((m, m)) < 0.7)
        np.fill_diagonal(W, 0)
        W = W + np.diag(np.ones(m - 1), 1) * 0.1 + np.diag(np.ones(m - 1), -1) * 0.1   # irreducible
        np.fill_diagonal(W, 0)
        M = W - np.diag(W.sum(axis=0))
        qq = rng.uniform(0.2, 2.0, m)
        kappa = rng.uniform(0, 2, m) * (rng.random(m) < 0.5)
        if not kappa.any():
            kappa[0] = 1.0
        g = rng.uniform(0.05, 1.0)
        Mk = M - np.diag(kappa)
        ev, vr = np.linalg.eig(Mk)
        k = int(np.argmax(ev.real))
        mu = -float(ev[k].real)
        ur = np.abs(vr[:, k].real)
        evl, vl = np.linalg.eig(Mk.T)
        ul = np.abs(vl[:, int(np.argmax(evl.real))].real)
        nu = ul * ur / (ul @ ur)
        K = (BETA * qq)[:, None] * np.linalg.inv(g * np.eye(m) - Mk)
        R = float(np.max(np.linalg.eigvals(K).real))
        lo, hi = BETA * float(nu @ qq) / (g + mu), BETA * float(qq.max()) / (g + mu)
        worst_lo = max(worst_lo, (lo - R) / R)
        worst_hi = max(worst_hi, (R - hi) / R)
        Ks = (BETA * qq)[:, None] * np.linalg.inv((g + mu) * np.eye(m) - M)
        n_sub_below += int(float(np.max(np.linalg.eigvals(Ks).real)) < lo)
    out["random_nonreversible"] = dict(instances=1500, max_rel_violation_lower=worst_lo, max_rel_violation_upper=worst_hi,
                                       substitution_below_lower_bound_in=n_sub_below)
    assert worst_lo < 1e-9 and worst_hi < 1e-9
    return out


if __name__ == "__main__":
    OUT = {"beta": BETA, "gamma": GAMMA, "A_partitions": block_A(), "B_leak_heterogeneous_q": block_B()}
    (DATA / "s4_geometry_barriers_leak.json").write_text(json.dumps(OUT, indent=1))
    print(json.dumps(OUT, indent=1))
