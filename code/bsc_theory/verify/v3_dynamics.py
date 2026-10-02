"""Checks with separately written code: (B) nonlinear two-zone outbreak below well-mixed threshold; (A) threshold; (C) transient TV distance (linear, eigen-decomposition);
   IBM first-generation offspring numbers (own event-driven simulation with small dt)."""
import sys, numpy as np, time
from scipy.integrate import solve_ivp
sys.path.insert(0, '../scripts'); import layouts as lay
sys.path.insert(0, '.'); 
B_, G_ = 0.5, 0.14
def P(*a): print(*a, flush=True)
def build(z, Dmap, qmap, D0=1.0):
    ny, nx = z.shape; mask = z != 'B'; idx = -np.ones(z.shape, int); idx[mask] = np.arange(mask.sum()); n = int(mask.sum())
    D = np.zeros(n); q = np.zeros(n); zone = np.empty(n, dtype='<U1'); M = np.zeros((n, n))
    for y in range(ny):
        for x in range(nx):
            if mask[y, x]: D[idx[y, x]] = Dmap[z[y, x]] * D0; q[idx[y, x]] = qmap[z[y, x]]; zone[idx[y, x]] = z[y, x]
    for y in range(ny):
        for x in range(nx):
            if not mask[y, x]: continue
            for dy, dx in ((0, 1), (1, 0)):
                yy, xx = y + dy, x + dx
                if yy < ny and xx < nx and mask[yy, xx]:
                    i, j = idx[y, x], idx[yy, xx]; w = 2 * D[i] * D[j] / (D[i] + D[j]); M[j, i] += w; M[i, j] += w
    M -= np.diag(M.sum(axis=0)); return M, q, zone, idx
def R0(M, q, b, g=G_):
    n = len(q); return float(np.max(np.abs(np.linalg.eigvals(b * q[:, None] * np.linalg.inv(g * np.eye(n) - M)))))
def ode(M, q, b, I0, T):
    n = len(q)
    def rhs(t, y):
        S, I, R = y[:n], y[n:2*n], y[2*n:3*n]; inf = b * q * S * I / (S + I + R)
        return np.concatenate([M @ S - inf, M @ I + inf - G_ * I, M @ R + G_ * I, inf])
    y0 = np.concatenate([1 - I0, I0, np.zeros(n), np.zeros(n)])
    sol = solve_ivp(rhs, (0, T), y0, method='LSODA', rtol=1e-9, atol=1e-13, t_eval=[T]); y = sol.y[:, -1]
    return (y[n:2*n].sum() + y[2*n:3*n].sum()) / n, y[3*n:], y[n:2*n].sum() / n
t0 = time.time()
# (B)
M, q, zone, idx = build(lay.two_zone('side'), lay.D_TWO, lay.Q_TWO); n = len(q)
b = 0.9 * G_ / q.mean(); r = R0(M, q, b); P(f'(B) beta={b:.5f}: well-mixed {b*q.mean()/G_:.3f}, NGM R0={r:.4f}')
att, cum, Iend = ode(M, q, b, np.full(n, 1e-6), 6000.0); P(f'    heterogeneous: final attack {att:.4f} (I at end {Iend:.2e}); share of infections in W {cum[zone=="W"].sum()/cum.sum():.4f} (area {np.mean(zone=="W"):.3f})')
att2, _, _ = ode(M, np.full(n, q.mean()), b, np.full(n, 1e-6), 6000.0); P(f'    homogenised q: final attack {att2:.3e}')
# (A)
M, q, zone, idx = build(lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE); n = len(q); ru = R0(M, q, 1.0)
for tgt in (0.8, 0.95, 1.05, 1.2, 2.0):
    att, _, _ = ode(M, q, tgt / ru, np.full(n, 1e-6), 3000.0); P(f'(A) R0={tgt}: final attack {att:.3e}')
# (C) linear transient via symmetric eigen-decomposition; TV distance to u when prevalence (sum I / n) crosses 1e-4, 1e-2 from a 1e-12 seed
for D0 in (1.0, 10.0, 100.0):
    M, q, zone, idx = build(lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE, D0); n = len(q)
    J = M + np.diag(B_ * q - G_); w, U = np.linalg.eigh(J); lam = w[-1]; u = np.abs(U[:, -1]); u /= u.sum()
    for nm, cell in (('far desk', (2, 2)), ('meeting room', (10, 17))):
        I0 = np.zeros(n); I0[idx[cell]] = 1e-12; c = U.T @ I0
        def prof(t): It = U @ (np.exp(w * t) * c); return It
        ts = np.linspace(0, 70, 7001); tot = np.array([prof(t).sum() / n for t in ts])
        out = []
        for lev in (1e-4, 1e-2):
            i = int(np.argmax(tot > lev)); It = prof(ts[i]); out.append((ts[i], 0.5 * np.abs(It / It.sum() - u).sum()))
        win = (tot > 1e-9) & (tot < 1e-4); slope = np.polyfit(ts[win], np.log(tot[win]), 1)[0]
        P(f'(C) D0={D0:g} seed {nm}: lam1={lam:.4f} gap={w[-1]-w[-2]:.4f}; fitted rate {slope:.4f}; TV at 1e-4: {out[0][1]:.4f} (t={out[0][0]:.1f}); TV at 1e-2 (linear approx): {out[1][1]:.4f} (t={out[1][0]:.1f})')
P(f'ODE part {time.time()-t0:.0f}s')
# IBM
def c_factor(N, p): return 1 - (1 - (1 - p) ** N) / (N * p)
def ibm(N, D0, nrep, rule, seed, NX=20, NY=13, Q=1.4):
    rng = np.random.default_rng(seed); n = NX * NY; dt = min(0.004, 0.02 / D0)
    x = rng.integers(0, NX, (nrep, N)); y = rng.integers(0, NY, (nrep, N)); sus = np.ones((nrep, N), bool); sus[:, 0] = False
    life = rng.exponential(1 / G_, nrep); cnt = np.zeros(nrep); alive = np.arange(nrep); t = 0.0
    while len(alive):
        m = len(alive); xa = x[alive]; ya = y[alive]
        r = rng.random((m, N)); p = D0 * dt
        dx = np.where(r < p, 1, np.where(r < 2 * p, -1, 0)); dy = np.where((r >= 2 * p) & (r < 3 * p), 1, np.where((r >= 3 * p) & (r < 4 * p), -1, 0))
        xn = xa + dx; yn = ya + dy; ok = (xn >= 0) & (xn < NX) & (yn >= 0) & (yn < NY); xa = np.where(ok, xn, xa); ya = np.where(ok, yn, ya)
        x[alive] = xa; y[alive] = ya
        same = (xa == xa[:, :1]) & (ya == ya[:, :1]); ntot = same.sum(axis=1)
        rate = B_ * Q / ntot if rule == 'local' else np.full(m, B_ * Q * n / N)
        sa = sus[alive]; hit = same & sa & (rng.random((m, N)) < (1 - np.exp(-rate * dt))[:, None]); cnt[alive] += hit.sum(axis=1); sus[alive] = sa & ~hit
        t += dt; alive = alive[life[alive] > t]
    return cnt
for N, D0, nrep, rule in ((100, 10.0, 2500, 'local'), (100, 1.0, 2500, 'local'), (100, 0.1, 2500, 'local'), (100, 10.0, 1500, 'mean'), (400, 10.0, 1200, 'local')):
    c = ibm(N, D0, nrep, rule, seed=5 + N + int(D0 * 10)); pred = (B_ * 1.4 / G_) * (c_factor(N, 1 / 260) if rule == 'local' else (N - 1) / N)
    P(f'IBM N={N} D0={D0} rule={rule}: offspring {c.mean():.3f} +- {c.std(ddof=1)/np.sqrt(nrep):.3f}; annealed prediction {pred:.3f}  [{time.time()-t0:.0f}s]')
