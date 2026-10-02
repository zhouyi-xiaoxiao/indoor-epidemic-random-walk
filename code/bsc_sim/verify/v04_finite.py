"""Finite-size: closed form vs the exact pair solve vs the simulated counts of the re-check (gmax=0)."""
import json, time, numpy as np, vlib as V
B, G = V.BETA, V.GAMMA
out = {}
occ = 99 / 260
def closed(L, N, D, bc='refl'):
    f = np.pi if bc == 'refl' else 2 * np.pi
    m1 = 2 * D * (1 - np.cos(f * np.arange(L) / L))
    mu = m1[:, None] + m1[None, :]
    g0 = np.mean(1 / (G + 2 * mu))
    return (B / G) / (1 + B / ((N - 1) / L ** 2) * g0), g0
for D, L, n in [(1.0, 4, 200000), (1.0, 10, 100000), (1.0, 20, 60000), (10.0, 10, 100000), (10.0, 14, 60000),
                (0.51, 20, 60000), (0.51, 5, 100000), (1.0, 40, 30000)]:
    N = int(round(occ * L * L)) + 1
    r = V.uniform(L, L, N, D=D)
    t0 = time.time()
    cf, g0 = closed(L, N, D)
    ex = np.nan
    if L <= 20:
        p = V.pair_solve(r)
        ex = (N - 1) * p.mean()
    s = V.sim(r, n, 31337 + L + int(10 * D), gmax=0)
    m, se = V.mci(s['off'])
    out[f'D{D}_L{L}'] = dict(N=N, closed=cf, exact=ex, sim=m, se=se, n=n)
    print(f'D={D} L={L} N={N}: closed {cf:.4f} exact {ex:.4f} sim {m:.4f}+-{se:.4f} z_exact={(m-ex)/se if ex==ex else (m-cf)/se:.2f} [{time.time()-t0:.0f}s]', flush=True)
    json.dump(out, open('v04_finite.json', 'w'), indent=1, default=float)
# periodic lattice: closed form should be EXACT -> test with exact pair solve on a torus
def torus_gen(L, D):
    M = L * L
    Lm = np.zeros((M, M))
    for x in range(L):
        for y in range(L):
            k = x * L + y
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                t = ((x + dx) % L) * L + (y + dy) % L
                Lm[t, k] += D; Lm[k, k] -= D
    return Lm
for L, D in ((6, 1.0), (8, 0.51)):
    N = int(round(occ * L * L)) + 1
    r = V.uniform(L, L, N, D=D)
    p = V.pair_solve(r, Lm=torus_gen(L, D))
    cf, g0 = closed(L, N, D, 'per')
    print(f'torus L={L} D={D}: closed {cf:.8f} exact {(N-1)*p.mean():.8f}')
    out[f'torus_L{L}_D{D}'] = dict(closed=cf, exact=(N - 1) * p.mean())
json.dump(out, open('v04_finite.json', 'w'), indent=1, default=float)
