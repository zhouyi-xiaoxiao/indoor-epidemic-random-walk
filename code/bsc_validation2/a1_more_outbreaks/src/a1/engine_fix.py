"""Round-off-safe M2 exposure (post-hoc numerical correction, POSTHOC_LOG P5).

The frozen spectral kernel (lattice.py) sums cosine modes; its absolute error is ~1e-16 times the
largest exposure, so exposures smaller than that (far sites at small D) are round-off noise of either
sign.  The lattice walk is separable, and each 1-D reflecting propagator is a sum of positive terms
(method of images, modified Bessel functions):

    g(x, s | x0) = sum_m e^{-2ws} [ I_{|x - x0 - 2nm|}(2ws) + I_{|x + x0 + 1 - 2nm|}(2ws) ],  w = D / a^2.

X = int_0^T (T - s) e^{-kappa s} gx gy ds is then a quadrature of a positive integrand (Gauss-Legendre in
log s), accurate component-wise.  The spectral value is kept wherever it is resolved (X > 1e-8 max X);
the image value replaces it elsewhere.  The M0 limit (xbar) is unchanged.
"""
from __future__ import annotations
import numpy as np
from scipy.special import ive
from . import lattice as L
from . import engine as G

NQ = 240
_gl = np.polynomial.legendre.leggauss(NQ)
MIMG = np.arange(-3, 4)
THRESH = 1e-8


def _g1d(n, x0, w, s):
    """(n, len(s)) positive 1-D reflecting-lattice propagator from site x0."""
    x = np.arange(n)
    z = 2.0 * w * s                                               # (q,)
    o1 = np.abs(x[:, None] - x0 - 2 * n * MIMG[None, :])          # (n, m)
    o2 = np.abs(x[:, None] + x0 + 1 - 2 * n * MIMG[None, :])
    return ive(o1[:, :, None], z[None, None, :]).sum(1) + ive(o2[:, :, None], z[None, None, :]).sum(1)


def x_images(d, D):
    lat, kappa = d["lat"], d["kappa"]
    X = np.zeros(lat.M)
    for T in d["segs"]:
        lo, hi = np.log(T * 1e-9), np.log(T)
        ls = 0.5 * (hi - lo) * _gl[0] + 0.5 * (hi + lo)
        s = np.exp(ls)
        wq = 0.5 * (hi - lo) * _gl[1] * s * (T - s) * np.exp(-kappa * s)
        for src in np.atleast_1d(d["src"]):
            x0, y0 = int(src % lat.nx), int(src // lat.nx)
            gx = _g1d(lat.nx, x0, D / lat.ax ** 2, s)
            gy = _g1d(lat.ny, y0, D / lat.ay ** 2 * lat.theta_y, s)
            X += np.einsum("yq,xq,q->yx", gy, gx, wq).reshape(-1)
    return X


def x_m2(d, D, model="M2"):
    Xs, xbar = G._x_m2_spectral(d, D, model)
    if model != "M2":
        return Xs, xbar
    bad = Xs < THRESH * Xs.max()
    if bad.any():
        Xi = x_images(d, D)[d["rec"]]
        Xs = np.where(bad, Xi, Xs)
    return np.clip(Xs, 1e-300, None), xbar


def install():
    if not hasattr(G, "_x_m2_spectral"):
        G._x_m2_spectral = G.x_m2
        G.x_m2 = x_m2
