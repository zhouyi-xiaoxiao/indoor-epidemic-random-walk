import sys, numpy as np
sys.path.insert(0, '../scripts'); import layouts as lay
B_, G_ = 0.5, 0.14
def build(z, Dmap, qmap, variant):
    ny, nx = z.shape; mask = z != 'B'; idx = -np.ones(z.shape, int); idx[mask] = np.arange(mask.sum()); n = int(mask.sum())
    D = np.zeros(n); q = np.zeros(n); nb = [[] for _ in range(n)]
    for y in range(ny):
        for x in range(nx):
            if mask[y, x]:
                i = idx[y, x]; D[i] = Dmap[z[y, x]]; q[i] = qmap[z[y, x]]
                for dy, dx in ((0, 1), (1, 0), (0, -1), (-1, 0)):
                    yy, xx = y + dy, x + dx
                    if 0 <= yy < ny and 0 <= xx < nx and mask[yy, xx]: nb[i].append(idx[yy, xx])
    M = np.zeros((n, n))
    for i in range(n):
        for j in nb[i]:
            if variant == 'D per neighbour': r = D[i]
            elif variant == 'D/4 per neighbour (T=1/4, rejection at walls)': r = D[i] / 4
            elif variant == 'D/deg per neighbour (T uniform over accessible)': r = D[i] / len(nb[i])
            M[j, i] += r
    M -= np.diag(M.sum(axis=0)); return M, q
def R0(M, q):
    n = len(q); return float(np.max(np.abs(np.linalg.eigvals(B_ * q[:, None] * np.linalg.inv(G_ * np.eye(n) - M)))))
def stat(M):
    n = M.shape[0]; A = M.copy(); A[-1, :] = 1; r = np.zeros(n); r[-1] = 1; return np.linalg.solve(A, r)
for name, z, Dm, qm in (('O1', lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE), ('two-zone aisles', lay.two_zone('aisles'), lay.D_TWO, lay.Q_TWO)):
    for v in ('D per neighbour', 'D/4 per neighbour (T=1/4, rejection at walls)', 'D/deg per neighbour (T uniform over accessible)'):
        M, q = build(z, Dm, qm, v); pi = stat(M); print(f'{name:16s} {v:50s}: <q>_pi={pi@q:.4f} lower {B_*(pi@q)/G_:.3f} R0={R0(M,q):.4f}')
