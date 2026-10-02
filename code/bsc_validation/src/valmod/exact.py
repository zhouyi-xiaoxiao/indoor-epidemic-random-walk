"""Bridge to the verified exact simulator (../bsc_sim) for model M1.

Runs the per-cell rule (same site, beta q n_I / n_total) at event
scale: only the index case is infectious (gmax = 0), nobody recovers
(gamma = 0), everybody starts on a seat and performs the lattice walk with
diffusion coefficient D (m^2/h) on the rectangular lattice of the event.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_SRC = os.path.normpath(os.path.join(HERE, "..", "..", "..", "bsc_sim", "src"))
if SIM_SRC not in sys.path:
    sys.path.insert(0, SIM_SRC)

from bsc_sim import sim as _sim            # noqa: E402
from bsc_sim.scenes import Scene           # noqa: E402

from .lattice import Lattice               # noqa: E402


def _scene(lat: Lattice, N: int) -> Scene:
    ny, nx = lat.ny, lat.nx
    return Scene(name="event", nx=nx, ny=ny, a_m=1.0,
                 zone=np.full((ny, nx), "W"), access=np.ones((ny, nx), bool),
                 D_grid=np.ones((ny, nx)), q_grid=np.ones((ny, nx)), N=N)


def _aniso_hop_rates(lat: Lattice, D: float):
    wx = D / lat.ax ** 2
    wy = D / lat.ay ** 2 * lat.theta_y
    ptr, idx, W = [0], [], []
    for s in range(lat.M):
        x, y = s % lat.nx, s // lat.nx
        for dx, dy, w in ((1, 0, wx), (-1, 0, wx), (0, 1, wy), (0, -1, wy)):
            xx, yy = x + dx, y + dy
            if 0 <= xx < lat.nx and 0 <= yy < lat.ny:
                idx.append(yy * lat.nx + xx)
                W.append(w)
        ptr.append(len(idx))
    return (np.array(ptr, np.int32), np.array(idx, np.int32), np.array(W, float))


def simulate_m1(lat: Lattice, src: int, rec, beta: float, D: float, T: float,
                n_rep: int = 4000, seed: int = 0) -> np.ndarray:
    """Infection frequency of every susceptible in ``rec`` over one exposure
    segment of duration T hours.  Returns (p_j, n_rep)."""
    rec = np.asarray(rec, dtype=np.int32)
    N = len(rec) + 1
    sc = _scene(lat, N)
    init = np.tile(np.concatenate([[src], rec]).astype(np.int32), (n_rep, 1))
    rates = _aniso_hop_rates(lat, D)
    orig = _sim.hop_rates
    _sim.hop_rates = lambda scene, Darr=None: rates
    try:
        res = _sim.simulate(sc, beta=beta, gamma=0.0, n_rep=n_rep, seed=seed,
                            t_max=T, dt_out=T, mode="percell", init_pos=init,
                            gmax=0)
    finally:
        _sim.hop_rates = orig
    inf = (res.t_inf[:, 1:] >= 0)
    return inf.mean(axis=0), inf
