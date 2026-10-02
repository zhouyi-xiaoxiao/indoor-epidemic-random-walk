"""ibm_check.py -- why the per-cell individual-based rule is sub-critical, and the encounter-weighted NGM.

Per-cell rule: a susceptible sharing a cell with n_I infectives among n_tot individuals is infected at
rate beta q n_I / n_tot  (local frequency dependence with the INSTANTANEOUS local count).  One infective sharing its
cell with k susceptibles therefore infects at total rate beta q k/(k+1): zero when alone.

Supplementary Section S4.3 of the article: if the other individuals are at stationarity and independent (annealed approximation, valid when
partners are exchanged quickly), the first-generation offspring number of an index case is obtained from the NGM with
q replaced by   q~_r = q_r c_r ,   c_r = E[K/(K+1)],  K ~ Binomial(N-1, pi_r)
              = 1 - (1 - (1 - pi_r)^N) / (N pi_r).
This script measures the first-generation offspring number directly (secondary cases are made non-infectious, so every
infection is caused by the index case) and compares with (i) the mean-field value beta q/gamma, (ii) the annealed
formula, for the per-cell rule and for the 'mean-occupancy' rule  beta q n_I /(N pi_r).

Homogeneous open room 20 x 13, uniform D0 (per-neighbour hop rate, lattice units), uniform q = 1.4, beta = 0.5,
gamma = 0.14.  Output ../results/ibm_check.json (+ .log).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent.parent / "results"
NX, NY = 20, 13
BETA, GAMMA, Q = 0.5, 0.14, 1.4
LOG, RES = [], {}


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def c_factor(N, p):
    return 1.0 - (1.0 - (1.0 - p) ** N) / (N * p)


def offspring(N, D0, nrep, rule, seed, dt=None):
    """first-generation offspring numbers of one index case among N individuals (nrep independent rooms)."""
    rng = np.random.default_rng(seed)
    n = NX * NY
    if dt is None:
        dt = min(0.005, 0.025 / D0)  # total hop probability per step 4 D0 dt <= 0.1
    x = rng.integers(0, NX, (nrep, N))
    y = rng.integers(0, NY, (nrep, N))  # stationary (uniform) start; individual 0 is the index case
    sus = np.ones((nrep, N), bool)
    sus[:, 0] = False
    life = rng.exponential(1 / GAMMA, nrep)
    count = np.zeros(nrep)
    alive = np.arange(nrep)
    t = 0.0
    dxs = np.array([1, -1, 0, 0])
    dys = np.array([0, 0, 1, -1])
    while len(alive):
        m = len(alive)
        # movement: each individual attempts each of the 4 directions with prob D0*dt; blocked by walls -> stays
        r = rng.random((m, N))
        d = np.floor(r / (D0 * dt)).astype(int)  # 0..3 with prob D0 dt each, >= 4 means stay
        mv = d < 4
        dd = np.where(mv, d, 0)
        nx_ = x[alive] + np.where(mv, dxs[dd], 0)
        ny_ = y[alive] + np.where(mv, dys[dd], 0)
        ok = (nx_ >= 0) & (nx_ < NX) & (ny_ >= 0) & (ny_ < NY)
        xa = np.where(ok, nx_, x[alive])
        ya = np.where(ok, ny_, y[alive])
        x[alive], y[alive] = xa, ya
        # infection by the index case
        same = (xa == xa[:, [0]]) & (ya == ya[:, [0]])
        ntot = same.sum(axis=1)  # includes the index case
        if rule == "local":
            rate = BETA * Q / ntot
        else:  # 'mean': normalise by the stationary mean occupancy N/n
            rate = np.full(m, BETA * Q / (N / n))
        pinf = 1.0 - np.exp(-rate * dt)
        sa = sus[alive]
        hit = same & sa & (rng.random((m, N)) < pinf[:, None])
        count[alive] += hit.sum(axis=1)
        sa &= ~hit
        sus[alive] = sa
        t += dt
        alive = alive[life[alive] > t]
    return count


runs = [
    # (N, D0, nrep)
    (100, 0.1, 4000), (100, 1.0, 4000), (100, 10.0, 4000),
    (400, 10.0, 2000), (1600, 10.0, 800),
]
t0 = time.time()
mf = BETA * Q / GAMMA
log(f"mean-field (continuum) value beta q/gamma = {mf:.4f};  room {NX}x{NY} = {NX * NY} cells, uniform q = {Q}")
for N, D0, nrep in runs:
    p = 1.0 / (NX * NY)
    for rule in ("local", "mean"):
        cnt = offspring(N, D0, nrep, rule, seed=1000 + N + int(10 * D0) + (0 if rule == "local" else 7))
        mean, se = cnt.mean(), cnt.std(ddof=1) / np.sqrt(nrep)
        pred = mf * c_factor(N, p) if rule == "local" else mf * (N - 1) / N
        key = f"N={N},D0={D0:g},rule={rule}"
        RES[key] = dict(N=N, D0=D0, nrep=nrep, rule=rule, mean=float(mean), se=float(se), annealed_prediction=float(pred),
                        mean_field=mf, occupancy=N * p, p_zero_offspring=float((cnt == 0).mean()))
        log(f"  N = {N:4d} (occupancy {N * p:.2f}/cell), D0 = {D0:5.1f}, rule = {rule:5s}: offspring mean = {mean:.3f} +- {se:.3f} "
            f"(annealed prediction {pred:.3f}; mean-field {mf:.3f}); P(no secondary case) = {(cnt == 0).mean():.3f}  [{time.time() - t0:.0f} s]")
        (OUT / "ibm_check.json").write_text(json.dumps(RES, indent=1))  # checkpoint
(OUT / "ibm_check.log").write_text("\n".join(LOG) + "\n")
