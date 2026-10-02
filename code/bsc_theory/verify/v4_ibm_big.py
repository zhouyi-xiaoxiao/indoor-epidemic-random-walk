import numpy as np, time, sys
B_, G_ = 0.5, 0.14
def c_factor(N, p): return 1 - (1 - (1 - p) ** N) / (N * p)
def ibm(N, D0, nrep, seed, NX=20, NY=13, Q=1.4, dtfac=1.0):
    rng = np.random.default_rng(seed); n = NX * NY; dt = dtfac * min(0.004, 0.02 / D0)
    x = rng.integers(0, NX, (nrep, N)); y = rng.integers(0, NY, (nrep, N)); sus = np.ones((nrep, N), bool); sus[:, 0] = False
    life = rng.exponential(1 / G_, nrep); cnt = np.zeros(nrep); expo = np.zeros(nrep); alive = np.arange(nrep); t = 0.0
    while len(alive):
        m = len(alive); xa = x[alive]; ya = y[alive]
        r = rng.random((m, N)); p = D0 * dt
        dx = np.where(r < p, 1, np.where(r < 2 * p, -1, 0)); dy = np.where((r >= 2 * p) & (r < 3 * p), 1, np.where((r >= 3 * p) & (r < 4 * p), -1, 0))
        xn = xa + dx; yn = ya + dy; ok = (xn >= 0) & (xn < NX) & (yn >= 0) & (yn < NY); xa = np.where(ok, xn, xa); ya = np.where(ok, yn, ya)
        x[alive] = xa; y[alive] = ya
        same = (xa == xa[:, :1]) & (ya == ya[:, :1]); ntot = same.sum(axis=1)
        rate = B_ * Q / ntot
        expo[alive] += rate * (ntot - 1) * dt        # integrated total hazard = what the annealed formula predicts in expectation
        sa = sus[alive]; hit = same & sa & (rng.random((m, N)) < (1 - np.exp(-rate * dt))[:, None]); cnt[alive] += hit.sum(axis=1); sus[alive] = sa & ~hit
        t += dt; alive = alive[life[alive] > t]
    return cnt, expo
t0 = time.time()
for N, D0, nrep in ((100, 10.0, 9000), (100, 100.0, 2500)):
    c, e = ibm(N, D0, nrep, seed=77 + int(D0)); pred = (B_ * 1.4 / G_) * c_factor(N, 1 / 260)
    print(f'IBM local N={N} D0={D0}: offspring {c.mean():.3f} +- {c.std(ddof=1)/np.sqrt(nrep):.3f}; integrated hazard {e.mean():.3f} +- {e.std(ddof=1)/np.sqrt(nrep):.3f}; annealed prediction {pred:.3f} [{time.time()-t0:.0f}s]', flush=True)
