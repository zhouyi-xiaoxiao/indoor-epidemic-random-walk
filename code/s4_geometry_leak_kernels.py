"""s4_geometry_leak_kernels.py -- counter-examples for the leaky-room results with contact kernels other than
the same-cell kernel (Remark after Corollary 4.3 of the article).

Theorem 4.2(3)-(4) (leak formula R_room = beta q0 / (gamma + mu_leak) for uniform q, and the bounds
beta <q>_nu / (gamma + mu_leak) <= R_room <= beta q_max / (gamma + mu_leak)) is proved for the same-cell kernel,
F = beta Q.  With a row-stochastic kernel P, F = beta P^T Q, R_room := rho(beta P^T Q V_kappa^{-1}).  This script
evaluates three small instances in which the formula or a bound fails:

  (a) uniform q = 1, P uniform over a cell and its neighbours, three cells in a row, hop rate w = 1/day, a door
      with exit rate kappa = 1/day at one end: R_room differs from beta / (gamma + mu_leak);
  (b) the same corridor and door, P with rows (0,0,1), (0,1,0), (0,0,1): R_room exceeds the upper bound
      beta q_max / (gamma + mu_leak);
  (c) the counter-example of Proposition 3.15(iv) (q = (1,0,1), P uniform over a cell and its neighbours) with
      hop rate 0.01/day and a small door: R_room lies below the lower bound beta <q>_nu / (gamma + mu_leak),
      which tends to beta <q> / gamma as kappa -> 0 while R_room tends to R0 of the closed room.

  (d) office scene, one door cell at the end of a corridor arm (the door of s4_geometry_barriers_leak.py):
      mu_leak and R_room / R0 against the mobility scale, (i) for a door with a fixed exit rate kappa = 1/day,
      where mu_leak is bounded by <kappa>_pi (Theorem 4.2(2)) and saturates, and (ii) for kappa equal to the
      hop rate of the door cell (kappa proportional to D0, the convention of Table 5), where mu_leak is
      proportional to D0.

beta = 0.5/day and gamma = 0.14/day, the default values.  Pure linear algebra; a few seconds.
Output: ../data/s4_geometry_leak_kernels.json
"""
from __future__ import annotations

import json
from pathlib import Path

import sys

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "data" / "s4_geometry_leak_kernels.json"
BETA, GAMMA = 0.5, 0.14


def generator(n: int, w: float) -> np.ndarray:
    """Generator of a walk on n cells in a row with hop rate w between neighbours (reflecting ends)."""
    G = np.zeros((n, n))
    for x in range(n - 1):
        G[x + 1, x] += w
        G[x, x + 1] += w
    G -= np.diag(G.sum(axis=0))
    return G


def leaky(q, P, w, kappa):
    q, P, kappa = np.asarray(q, float), np.asarray(P, float), np.asarray(kappa, float)
    n = len(q)
    Gk = generator(n, w) - np.diag(kappa)
    V = GAMMA * np.eye(n) - Gk
    R = float(max(abs(np.linalg.eigvals(BETA * P.T @ np.diag(q) @ np.linalg.inv(V)))))
    ev, U = np.linalg.eig(Gk)
    i = int(np.argmax(ev.real))
    mu = float(-ev[i].real)
    u = np.abs(U[:, i].real)
    evl, Ul = np.linalg.eig(Gk.T)
    l = np.abs(Ul[:, int(np.argmax(evl.real))].real)
    nu = l * u / (l @ u)
    return dict(R_room=R, mu_leak=mu, nu=nu.tolist(),
                formula_uniform_q=float(BETA * q.max() / (GAMMA + mu)) if np.ptp(q) == 0 else None,
                upper_bound=float(BETA * q.max() / (GAMMA + mu)),
                lower_bound=float(BETA * (nu @ q) / (GAMMA + mu)))


nbr = np.array([[1 / 2, 1 / 2, 0], [1 / 3, 1 / 3, 1 / 3], [0, 1 / 2, 1 / 2]])
res = {"parameters": dict(beta=BETA, gamma=GAMMA, cells=3,
                          kernel_convention="F = beta P^T Q, P row-stochastic; V_kappa = gamma I - (G - diag kappa)")}
a = leaky([1, 1, 1], nbr, 1.0, [1.0, 0, 0])
a["relative_difference_from_formula"] = a["R_room"] / a["formula_uniform_q"] - 1
res["a_uniform_q_neighbour_kernel_end_door"] = a
b = leaky([1, 1, 1], [[0, 0, 1], [0, 1, 0], [0, 0, 1]], 1.0, [1.0, 0, 0])
b["excess_over_upper_bound"] = b["R_room"] / b["upper_bound"] - 1
res["b_uniform_q_skew_kernel_end_door"] = b
c = {}
for where, kap in (("end door", [None, 0, 0]), ("middle door", [0, None, 0])):
    for k in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
        kv = [k if v is None else v for v in kap]
        r = leaky([1, 0, 1], nbr, 0.01, kv)
        c[f"{where}|kappa={k:g}"] = dict(R_room=r["R_room"], lower_bound=r["lower_bound"], mu_leak=r["mu_leak"],
                                         bound_holds=r["R_room"] >= r["lower_bound"] - 1e-12)
res["c_prop_3_15_iv_kernel_q101_D0_0.01"] = c
res["c_closed_room_limit"] = dict(beta_mean_q_over_gamma=BETA * (2 / 3) / GAMMA,
                                  R0_closed_formula=BETA / GAMMA * (GAMMA + 4 * 0.01) / (2 * GAMMA + 6 * 0.01))

# (d) office scene
ROOT = HERE.parents[0]
sys.path.insert(0, str(ROOT / "code/bsc_sim/src"))
from bsc_sim import scenes as sc            # noqa: E402
from bsc_sim import theory as th            # noqa: E402

s1 = sc.load_scene("office", D0=1.0)
n, q, zones = s1.M, s1.q, s1.site_zone
L1 = th.generator(s1).toarray()
cor = np.flatnonzero(zones == "C")
xy = s1.site_xy
on_edge = (xy[cor, 0] == 0) | (xy[cor, 0] == s1.nx - 1) | (xy[cor, 1] == 0) | (xy[cor, 1] == s1.ny - 1)
cand = cor[on_edge]
mroom = xy[zones == "M"].mean(axis=0)
door = int(cand[np.argmax(np.abs(xy[cand] - mroom).sum(axis=1))])
hop1 = float(np.max(np.delete(L1[:, door], door)))      # hop rate of the door cell to a neighbour at D0 = 1


def r_spec(L, kappa):
    V = GAMMA * np.eye(n) - L + np.diag(kappa)
    return float(np.max(np.linalg.eigvals((BETA * q)[:, None] * np.linalg.inv(V)).real))


d = {"door_cell_xy": [int(xy[door, 0]), int(xy[door, 1])], "n_cells": int(n), "door_hop_rate_at_D0_1": hop1,
     "kappa_mean_pi_for_kappa_1": 1.0 / n, "fixed_kappa_1": {}, "kappa_equal_door_hop_rate": {}}
for D0 in (1.0, 10.0, 100.0, 2400.0):
    L = D0 * L1
    R0 = r_spec(L, np.zeros(n))
    for name, kap in (("fixed_kappa_1", 1.0), ("kappa_equal_door_hop_rate", hop1 * D0)):
        kappa = np.zeros(n)
        kappa[door] = kap
        mu = -float(np.linalg.eigvalsh(L - np.diag(kappa))[-1])
        d[name][f"D0={D0:g}"] = dict(kappa=kap, mu_leak=mu, mu_over_D0=mu / D0, R0_closed=R0,
                                     R_room=r_spec(L, kappa), R_room_over_R0=r_spec(L, kappa) / R0)
res["d_office_door_against_mobility"] = d
OUT.write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
