import numpy as np, mpmath as mp, sys
from scipy.optimize import minimize
rng = np.random.default_rng(99)
def R0f(M, q, b, g):
    n = len(q); K = b * q[:, None] * np.linalg.inv(g * np.eye(n) - M); return float(np.max(np.abs(np.linalg.eigvals(K))))
def R0mp(M, q, b, g, D):
    n = len(q); Mm = mp.matrix(M.tolist()) * mp.mpf(D); V = mp.mpf(g) * mp.eye(n) - Mm; Vi = V ** -1
    K = mp.matrix(n, n)
    for i in range(n):
        for j in range(n): K[i, j] = mp.mpf(b) * mp.mpf(float(q[i])) * Vi[i, j]
    ev = mp.eig(K, left=False, right=False); return max(abs(e) for e in ev)
def unpack(z, n, clip):
    z = np.clip(z, -clip, clip)
    W = np.exp(z[:n * n]).reshape(n, n); np.fill_diagonal(W, 0); M = W - np.diag(W.sum(axis=0)); q = np.exp(z[n * n:n * n + n]); return M, q, z
mp.mp.dps = 50
for clip in (6.0, 15.0, 40.0):
    worst = -1e9; worst_mp = None
    for trial in range(60):
        n = int(rng.integers(3, 6))
        def f(z):
            M0, q, zz = unpack(z, n, clip); D1 = np.exp(zz[-2]); D2 = D1 * (1 + np.exp(zz[-1]))
            return -(R0f(D2 * M0, q, 1.0, 0.3) - R0f(D1 * M0, q, 1.0, 0.3))
        res = minimize(f, rng.normal(0, 1.5, n * n + n + 2), method='Nelder-Mead', options=dict(maxiter=2500, xatol=1e-7, fatol=1e-15))
        if -res.fun > worst:
            worst = -res.fun
            M0, q, zz = unpack(res.x, n, clip); D1 = np.exp(zz[-2]); D2 = D1 * (1 + np.exp(zz[-1]))
            worst_mp = (R0mp(M0, q, 1.0, 0.3, D2) - R0mp(M0, q, 1.0, 0.3, D1), np.exp(zz[:n*n]).min(), np.exp(zz[:n*n]).max(), D1, D2, np.linalg.cond(0.3*np.eye(n)-D2*M0))
    print(f'clip {clip}: float64 worst increase {worst:.3e}; same instance in 50-digit arithmetic: {mp.nstr(worst_mp[0], 8)}; rates [{worst_mp[1]:.2e},{worst_mp[2]:.2e}] D1={worst_mp[3]:.2e} D2={worst_mp[4]:.2e} cond(V)={worst_mp[5]:.2e}', flush=True)
