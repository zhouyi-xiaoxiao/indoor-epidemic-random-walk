"""Validate the second simulator against exact results, then test the key claims."""
import json, numpy as np, time, sys
import vlib as V
B, G = V.BETA, V.GAMMA
out = {}
# (a) one cell, N=40: exact final-size law of Markovian SIR with pair hazard beta/(N-1)
def exact_fs(N, lam, g):
    # state (s,i) ; prob of final size
    from functools import lru_cache
    import sys
    sys.setrecursionlimit(10000)
    P = {(N - 1, 1): 1.0}
    res = np.zeros(N + 1)
    for tot in range(N, 0, -1):            # s+i decreasing or same..., do BFS by layers of r=N-s-i then s
        pass
    # iterative: order by number removed r then by s descending
    probs = np.zeros((N + 1, N + 1))
    probs[N - 1, 1] = 1
    for r in range(0, N + 1):
        for s in range(N - r, -1, -1):
            i = N - r - s
            if i <= 0: continue
            p = probs[s, i]
            if p == 0: continue
            a = lam * s * i; b = g * i
            if s > 0: probs[s - 1, i + 1] += p * a / (a + b)
            if i == 1: res[N - s] += p * b / (a + b)
            else: probs[s, i - 1] += p * b / (a + b)
    return res
N = 40
r = V.uniform(1, 1, N)
t0 = time.time()
s = V.sim(r, 40000, 12345)
pex = exact_fs(N, B / (N - 1), G)
emp = np.bincount(s['final'], minlength=N + 1)
m, se = V.mci(s['final'])
E = pex * 40000
# chi2 merging bins with E<5
chi = 0; dof = -1; ce = co = 0
for k in range(1, N + 1):
    ce += E[k]; co += emp[k]
    if ce >= 5: chi += (co - ce) ** 2 / ce; dof += 1; ce = co = 0
print('well-mixed N=40: exact mean %.3f sim %.3f +- %.3f ; chi2 %.1f / %d dof ; sum pex %.6f ; %.1fs' % (np.arange(N + 1) @ pex, m, se, chi, dof, pex.sum(), time.time() - t0))
out['wellmixed'] = dict(exact=float(np.arange(N + 1) @ pex), sim=m, se=se, chi2=chi, dof=dof)
# (b) pair theory vs the simulated counts of the re-check (gmax=0), small heterogeneous lattice -> exact LU solve
acc = np.array([[1, 1, 1, 1], [1, 0, 1, 1], [1, 1, 1, 1]], bool)
Dg = [[0.4, 0.4, 1.5, 1.5]] * 3; qg = [[2.0, 2.0, 0.6, 0.6]] * 3
r = V.Room(acc, Dg, qg, 10)
for rad in (0.0, 1.0):
    P = V.kernel(r, rad) if rad > 0 else None
    p = V.pair_solve(r, P)
    R1 = (r.N - 1) * p.mean()
    s = V.sim(r, 200000, 777 + int(rad), P=P, gmax=0)
    m, se = V.mci(s['off'])
    print('small lattice rad', rad, 'pair R1 %.4f  sim %.4f +- %.4f  z=%.2f' % (R1, m, se, (m - R1) / se))
    out[f'small_rad{rad}'] = dict(R1=R1, sim=m, se=se)
    # exposure identity with move_all
    s = V.sim(r, 200000, 778 + int(rad), P=P, gmax=0, move_all=True)
    K = V.ngm(r, P)
    m, se = V.mci(s['expo'])
    print('   exposure %.4f +- %.4f vs beta<q>/gamma %.4f (mean col sum %.4f)' % (m, se, B * r.q.mean() / G, K.sum(0).mean()))
    out[f'small_expo{rad}'] = dict(theory=float(K.sum(0).mean()), sim=m, se=se)
json.dump(out, open('v02_simcheck.json', 'w'), indent=1, default=float)
