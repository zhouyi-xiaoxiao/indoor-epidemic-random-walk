"""Deterministic checks: NGM values, uniform-q identity, bounds, closed forms."""
import json, numpy as np, time
from scipy.special import ellipk
import vlib as V
B, G = V.BETA, V.GAMMA
out = {}
# 1. uniform q : R0 = beta q/gamma for any L, D, kernel, barriers
rng = np.random.default_rng(0)
rows = []
for L in (4, 7, 13, 20, 40):
    for D in (0.51, 1.0, 10.0):
        r = V.uniform(L, L, 10, D=D)
        K = V.ngm(r)
        rows.append((L, D, 'I', V.rho(K), float(np.abs(K.sum(0) - B / G).max())))
        P = V.kernel(r, 1.5)
        K = V.ngm(r, P)
        rows.append((L, D, 'r1.5', V.rho(K), float(np.abs(K.sum(0) - B / G).max())))
# heterogeneous D + barriers, uniform q
bar = rng.random((13, 20)) < 0.15
r = V.uniform(20, 13, 10, barriers=bar)
# keep largest comp not enforced: R0 identity holds regardless of connectivity
r.D = rng.uniform(0.05, 3.0, r.M)
K = V.ngm(r, V.kernel(r, 2.0))
rows.append(('het D+barriers', 'rand', 'r2', V.rho(K), float(np.abs(K.sum(0) - B / G).max())))
out['uniform_q'] = rows
print('uniform q: max |R0 - beta/gamma| =', max(abs(x[3] - B / G) for x in rows))
# 2. scene R0, bounds
sc = {}
for name in ('office', 'supermarket', 'classroom', 'metro'):
    for D0 in (1.0, 10.0, 100.0):
        r = V.load(name, D0)
        rad = 1.0 / r.a_m
        for kn, P in (('1m', V.kernel(r, rad) if rad >= 1 else None), ('same', None)):
            K = V.ngm(r, P)
            sc[f'{name}|{D0}|{kn}'] = (V.rho(K), B * r.q.mean() / G, B * r.q.max() / G, r.M, r.N)
            print(name, D0, kn, 'R0=%.3f  bounds %.3f %.3f  M=%d' % sc[f'{name}|{D0}|{kn}'][:4])
out['scenes'] = sc
# area-mean q and D
for name in ('office', 'supermarket', 'classroom', 'metro'):
    r = V.load(name, 1.0)
    print(name, 'mean q %.3f mean D %.4f' % (r.q.mean(), r.D.mean()), {z: (float(r.D[r.zone == z][0])) for z in set(r.zone)})
# office mobility
for D0 in (0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000):
    r = V.load('office', D0)
    print('office D0', D0, 'R0 %.3f' % V.rho(V.ngm(r)))
# 3. lower bound counter-example for kernel NGM
acc = np.ones((1, 3), bool)
r = V.Room(acc, [[1e-6] * 3], [[1.0, 0.0, 1.0]], 5)
P = V.kernel(r, 1.0)
print('counterexample 3 cells q=(1,0,1), D->0, radius-1 kernel: R0/(beta/gamma)=', V.rho(V.ngm(r, P)) / (B / G), ' <q>=', r.q.mean())
r = V.Room(acc, [[1e-6] * 3], [[1.0, 0.2, 1.0]], 5)
print('  q=(1,0.2,1):', V.rho(V.ngm(r, P)) / (B / G), ' <q>=', r.q.mean())
# does the lower bound hold in the four scenes at low D0 with the 1 m kernel?
for name in ('classroom', 'metro'):
    for D0 in (1e-3, 0.03, 1.0):
        r = V.load(name, D0)
        K = V.ngm(r, V.kernel(r, 1.0 / r.a_m))
        print(name, D0, 'R0 %.3f  lower %.3f upper %.3f' % (V.rho(K), B * r.q.mean() / G, B * r.q.max() / G))
# 4. closed form g0 : elliptic
for D in (0.51, 1.0, 10.0):
    k = 8 * D / (G + 8 * D)
    g_ell = 2 / np.pi * ellipk(k * k) / (G + 8 * D)
    kk = (np.arange(6000) + 0.5) * np.pi / 6000
    mu = 2 * D * (1 - np.cos(kk))[:, None] + 2 * D * (1 - np.cos(kk))[None, :]
    g_num = np.mean(1 / (G + 2 * mu))
    occ = 99 / 260
    print('D', D, 'g0 elliptic %.6f numeric %.6f  R1_inf %.4f' % (g_ell, g_num, (B / G) / (1 + B / occ * g_ell)))
json.dump(out, open('v01_det.json', 'w'), indent=1, default=float)
