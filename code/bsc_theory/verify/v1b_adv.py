import numpy as np, mpmath as mp
from scipy.optimize import minimize
rng = np.random.default_rng(424243)
def stat_eig(M):
    w, v = np.linalg.eig(M); k = np.argmin(np.abs(w)); p = np.real(v[:, k]); return p / p.sum()
def stat_solve(M):
    n = M.shape[0]; A = M.copy(); A[-1, :] = 1; r = np.zeros(n); r[-1] = 1; return np.linalg.solve(A, r)
def R0f(M, q, b, g):
    n = len(q); K = b * q[:, None] * np.linalg.inv(g * np.eye(n) - M); return float(np.max(np.abs(np.linalg.eigvals(K))))
def unpack(z, n):
    W = np.exp(z[:n * n]).reshape(n, n); np.fill_diagonal(W, 0); M = W - np.diag(W.sum(axis=0)); q = np.exp(z[n * n:n * n + n]); return M, q
worst = (1e9, None)
for trial in range(60):
    n = int(rng.integers(3, 6))
    def f(z, stat=stat_eig):
        M, q = unpack(z, n); pi = stat(M); return R0f(M, q, 1.0, 0.3) / ((pi @ q) / 0.3) - 1.0
    res = minimize(f, rng.normal(0, 1, n * n + n), method='Nelder-Mead', options=dict(maxiter=1500, xatol=1e-6, fatol=1e-14))
    if res.fun < worst[0]: worst = (res.fun, (res.x.copy(), n))
    if res.fun < -1e-9:
        M, q = unpack(res.x, n); pe = stat_eig(M); ps = stat_solve(M)
        print(f'trial {trial} n={n}: f(eig-stat)={res.fun:.3e}; pi_eig min {pe.min():.2e} |M pi_eig| {np.abs(M@pe).max():.2e}; pi_solve min {ps.min():.2e} |M pi_solve| {np.abs(M@ps).max():.2e}; f(solve-stat)={f(res.x, stat_solve):.3e}; rate range [{np.exp(res.x[:n*n]).min():.2e},{np.exp(res.x[:n*n]).max():.2e}] cond eigs {np.sort(np.abs(np.linalg.eigvals(M)))[:2]}')
        # high precision check
        mp.mp.dps = 60
        Mm = mp.matrix(M.tolist()); A = Mm.copy()
        for j in range(n): A[n - 1, j] = 1
        rhs = mp.matrix([0] * (n - 1) + [1]); pim = mp.lu_solve(A, rhs)
        V = mp.mpf('0.3') * mp.eye(n) - Mm; Vi = V ** -1
        K = mp.matrix(n, n)
        for i in range(n):
            for j in range(n): K[i, j] = mp.mpf(float(q[i])) * Vi[i, j]
        ev = mp.eig(K, left=False, right=False); R = max(abs(e) for e in ev)
        lo = sum(pim[i] * mp.mpf(float(q[i])) for i in range(n)) / mp.mpf('0.3')
        print('    mpmath: R0/lower - 1 =', mp.nstr(R / lo - 1, 10), ' min pi =', mp.nstr(min(pim), 5))
print('worst with eig-stat', worst[0])
# redo adversarial search with robust stationary solver
worst2 = 1e9
for trial in range(120):
    n = int(rng.integers(3, 6))
    def f2(z):
        M, q = unpack(np.clip(z, -12, 12), n); pi = stat_solve(M); return R0f(M, q, 1.0, 0.3) / ((pi @ q) / 0.3) - 1.0
    res = minimize(f2, rng.normal(0, 1, n * n + n), method='Nelder-Mead', options=dict(maxiter=3000, xatol=1e-7, fatol=1e-15))
    worst2 = min(worst2, res.fun)
print('robust adversarial min of R0/lower-1 (120 starts):', worst2)
