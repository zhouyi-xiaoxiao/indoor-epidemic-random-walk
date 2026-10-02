"""Event-scale exposure kernels on a rectangular lattice with reflecting walls.

Lattice: nx sites across (pitch ax, metres), ny sites along (pitch ay).  A
continuous-time nearest-neighbour walk with hop rates D/ax^2 (across) and
D/ay^2 (along), D in m^2/h, has the propagator

    G(r, s | r0) = sum_k phi_k(r) phi_k(r0) exp(-D mu_k s),

    phi_k(r) = c_kx cos(pi kx (x + 1/2) / nx) * c_ky cos(pi ky (y + 1/2) / ny),
    mu_k     = 2 (1 - cos(pi kx / nx)) / ax^2 + 2 (1 - cos(pi ky / ny)) / ay^2,

(the cosine-mode expansion).  The three exposure kernels
of the pre-registration are sums over the same modes with different weights:

    M0  X = (1/M) [T/kappa - (1 - e^{-kappa T}) / kappa^2]           (k = 0 only)
    M1  X = int_0^T G(r, 2t | r0) dt       -> w_k = (1 - e^{-2 D mu T}) / (2 D mu)
    M2  X = int_0^T (T - s) e^{-kappa s} G(r, s | r0) ds
                                           -> w_k = (L T - 1 + e^{-L T}) / L^2,
                                              L = kappa + D mu

Infection probability of susceptible j:  p_j = 1 - exp(-eps * X_j).
Units: time in hours, D in m^2/h, kappa in 1/h.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Lattice:
    nx: int
    ny: int
    ax: float = 1.0
    ay: float = 1.0
    H: float = 2.5                      # ceiling height (m), level only
    theta_y: float = 1.0                # factor on the along-axis hop rate (S-aniso)
    Phi: np.ndarray = field(init=False, repr=False)
    mu: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        nx, ny = self.nx, self.ny
        kx = np.arange(nx)
        ky = np.arange(ny)
        x = np.arange(nx)
        y = np.arange(ny)
        cx = np.sqrt(np.where(kx == 0, 1.0, 2.0) / nx)
        cy = np.sqrt(np.where(ky == 0, 1.0, 2.0) / ny)
        Px = cx[None, :] * np.cos(np.pi * np.outer(x + 0.5, kx) / nx)   # (x, kx)
        Py = cy[None, :] * np.cos(np.pi * np.outer(y + 0.5, ky) / ny)   # (y, ky)
        mux = 2.0 * (1.0 - np.cos(np.pi * kx / nx)) / self.ax ** 2
        muy = 2.0 * (1.0 - np.cos(np.pi * ky / ny)) / self.ay ** 2 * self.theta_y
        # site index s = y * nx + x ; mode index k = ky * nx + kx
        self.Phi = (Py[:, None, :, None] * Px[None, :, None, :]).reshape(nx * ny, nx * ny)
        self.mu = (muy[:, None] + mux[None, :]).reshape(nx * ny)
        self.mu[0] = 0.0

    @property
    def M(self) -> int:
        return self.nx * self.ny

    @property
    def cell_volume(self) -> float:
        return self.ax * self.ay * self.H

    def site(self, x, y):
        return np.asarray(y) * self.nx + np.asarray(x)

    def xy_m(self, s):
        s = np.asarray(s)
        return (s % self.nx + 0.5) * self.ax, (s // self.nx + 0.5) * self.ay


# ---------------------------------------------------------------------------
# mode weights
# ---------------------------------------------------------------------------
def _w_m2(L, T):
    """(L T - 1 + exp(-L T)) / L^2, stable for small L T."""
    L = np.asarray(L, float)
    x = L * T
    out = np.empty_like(x)
    small = x < 1e-4
    xs = x[small]
    out[small] = T * T * (0.5 - xs / 6.0 + xs * xs / 24.0)
    xl = x[~small]
    out[~small] = (xl - 1.0 + np.exp(-xl)) / L[~small] ** 2
    return out


def _w_m1(lam, T):
    """(1 - exp(-2 lam T)) / (2 lam), -> T as lam -> 0."""
    lam = np.asarray(lam, float)
    x = 2.0 * lam * T
    out = np.empty_like(x)
    small = x < 1e-8
    out[small] = T * (1.0 - x[small] / 2.0)
    out[~small] = -np.expm1(-x[~small]) / (2.0 * lam[~small])
    return out


def weights(lat: Lattice, model: str, D: float, kappa: float, T: float):
    """Mode weights w_k for one exposure segment of duration T (hours)."""
    if model == "M0":
        w = np.zeros(lat.M)
        w[0] = _w_m2(np.array([kappa]), T)[0]
        return w
    if model == "M1":
        return _w_m1(D * lat.mu, T)
    if model == "M2":
        return _w_m2(kappa + D * lat.mu, T)
    raise ValueError(model)


def exposure(lat: Lattice, model: str, D: float, kappa: float, segments,
             src: int, rec) -> np.ndarray:
    """Exposure X_j of receivers ``rec`` (site indices) to one source at ``src``.

    ``segments`` is a list of durations (hours); the air (M2) and the seating
    (M1) are reset at the start of every segment.
    """
    rec = np.asarray(rec)
    w = np.zeros(lat.M)
    for T in np.atleast_1d(segments):
        w += weights(lat, model, D, kappa, float(T))
    return lat.Phi[rec] @ (w * lat.Phi[src])


def exposure_matrix(lat: Lattice, model: str, D: float, kappa: float, segments,
                    sites=None) -> np.ndarray:
    """Symmetric exposure matrix X[i, j] between all pairs of ``sites``."""
    w = np.zeros(lat.M)
    for T in np.atleast_1d(segments):
        w += weights(lat, model, D, kappa, float(T))
    P = lat.Phi if sites is None else lat.Phi[np.asarray(sites)]
    return (P * w) @ P.T


def steady_weights(lat: Lattice, model: str, D: float, kappa: float, T: float = 1.0):
    """Per-unit-time weights for long exposures (M0/M2: steady concentration)."""
    if model == "M0":
        w = np.zeros(lat.M)
        w[0] = 1.0 / kappa
        return w
    if model == "M2":
        return 1.0 / (kappa + D * lat.mu)
    return _w_m1(D * lat.mu, T)


# ---------------------------------------------------------------------------
# brute-force reference (used by the self-test only)
# ---------------------------------------------------------------------------
def generator(lat: Lattice, D: float) -> np.ndarray:
    """Dense generator Q[r, r0] of the lattice walk (columns sum to zero)."""
    M = lat.M
    Q = np.zeros((M, M))
    wx = D / lat.ax ** 2
    wy = D / lat.ay ** 2 * lat.theta_y
    for y in range(lat.ny):
        for x in range(lat.nx):
            s = y * lat.nx + x
            for dx, dy, w in ((1, 0, wx), (-1, 0, wx), (0, 1, wy), (0, -1, wy)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < lat.nx and 0 <= yy < lat.ny:
                    t = yy * lat.nx + xx
                    Q[t, s] += w
                    Q[s, s] -= w
    return Q
