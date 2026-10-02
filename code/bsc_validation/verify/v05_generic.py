import numpy as np
from vlat import *
import vevents as E
r = np.load('res_holdout.npz'); evs = ['T1','T2','T5','C1']
def pmf_rr(rho, n1, n2, K):
    # independent Bernoulli with p_near = rho * p_far, expected total K
    pf = K/(n1*rho+n2); pn = min(rho*pf, 1-1e-9)
    if rho*pf > 1:  # cap
        pn = 1-1e-9; pf = (K-n1*pn)/n2
    return cond_pmf(np.full(n1,pn), np.full(n2,pf), K)
def pooled(pms):
    p0 = [np.clip(r[f'{e}_M0'],1e-300,None) for e in evs]; pm = [np.clip(p,1e-300,None) for p in pms]
    g = np.meshgrid(*[np.arange(len(p)) for p in p0], indexing='ij')
    d = sum(np.log(pm[a][g[a]])-np.log(p0[a][g[a]]) for a in range(4)); w = np.prod([p0[a][g[a]] for a in range(4)], axis=0)
    dobs = sum(np.log(pm[a][E.OBS[e][0]])-np.log(p0[a][E.OBS[e][0]]) for a,e in enumerate(evs))
    return dobs, w[d >= dobs-1e-9].sum()/w.sum()
print('M2 (verifier re-implementation): pooled', pooled([r[f'{e}_M2'] for e in evs]))
for rho in (1.5, 2, 2.5, 3, 4, 5, 8):
    pms = []; ps = []
    for e in evs:
        k1,n1,k2,n2 = E.OBS[e]; pm = pmf_rr(rho, n1, n2, k1+k2); pms.append(pm); ps.append(two_sided(pm, k1))
    d, p = pooled(pms)
    print(f'constant RR={rho}: pooled delta {d:+.2f}, p under M0 {p:.5f}; A2 p-values', np.round(ps,3), 'failures', [e for e,q in zip(evs,ps) if q<0.05])
