"""Adversarial search for an increase of R0(D) with nonreversible generators; candidates re-evaluated in 60-digit arithmetic
with an EXACT zero-column-sum generator (diagonal = minus exact sum of the float off-diagonals)."""
import numpy as np, mpmath as mp
from scipy.optimize import minimize
rng = np.random.default_rng(2026)
mp.mp.dps = 60
def R0f(M, q, b, g):
    n = len(q); K = b * q[:, None] * np.linalg.inv(g * np.eye(n) - M); return float(np.max(np.abs(np.linalg.eigvals(K))))
def R0mp_exact(W, q, D, g='0.3'):
    n = len(q); Wm = mp.matrix(W.tolist())
    M = mp.matrix(n, n)
    for j in range(n):
        s = mp.mpf(0)
        for i in range(n):
            if i != j: M[i, j] = Wm[i, j]; s += Wm[i, j]
        M[j, j] = -s
    V = mp.mpf(g) * mp.eye(n) - mp.mpf(D) * M; Vi = V ** -1
    K = mp.matrix(n, n)
    for i in range(n):
        for j in range(n): K[i, j] = mp.mpf(float(q[i])) * Vi[i, j]
    return max(abs(e) for e in mp.eig(K, left=False, right=False))
def unpack(z, n, clip):
    z = np.clip(z, -clip, clip); W = np.exp(z[:n * n]).reshape(n, n); np.fill_diagonal(W, 0); q = np.exp(z[n * n:n * n + n]); return W, q, z
for clip in (4.0, 8.0, 15.0):
    best = None
    for trial in range(40):
        n = int(rng.integers(3, 6))
        def f(z):
            W, q, zz = unpack(z, n, clip); M0 = W - np.diag(W.sum(axis=0)); D1 = np.exp(zz[-2]); D2 = D1 * (1 + np.exp(zz[-1]))
            return -(R0f(D2 * M0, q, 1.0, 0.3) - R0f(D1 * M0, q, 1.0, 0.3))
        res = minimize(f, rng.normal(0, 1.5, n * n + n + 2), method='Nelder-Mead', options=dict(maxiter=2500, xatol=1e-7, fatol=1e-15))
        W, q, zz = unpack(res.x, n, clip); D1 = np.exp(zz[-2]); D2 = D1 * (1 + np.exp(zz[-1]))
        ex = R0mp_exact(W, q, D2) - R0mp_exact(W, q, D1)
        if best is None or ex > best[0]: best = (ex, -res.fun, D1, D2)
    print(f'clip {clip}: max over 40 optimised instances of EXACT-generator 60-digit R0(D2)-R0(D1) = {mp.nstr(best[0], 8)} (float64 objective there {best[1]:.3e}; D1={best[2]:.3g}, D2={best[3]:.3g})', flush=True)
