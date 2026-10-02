import numpy as np, sys, json
from scipy import stats
from vlat import *
import vevents as E
rng = np.random.default_rng(4242)
t4 = np.load('res_T4.npz'); r1 = np.load('res_R1_primary.npz')
logD = t4['logD']; step = logD[1]-logD[0]
def draws(w, n):
    w = w/w.sum(); i = rng.choice(len(w), n, p=w); return 10**(logD[i] + rng.uniform(-step/2, step/2, n))
ND = 600
post = {('veh','M1'): draws(np.exp(t4['ll1']-t4['ll1'].max()), ND),
        ('veh','M2'): draws(np.exp(t4['ll2']-t4['ll2'].max()).sum(axis=1), ND),
        ('room','M1'): draws(np.exp(r1['ll1']-r1['ll1'].max()), ND),
        ('room','M2'): draws(np.exp(r1['ll2']-r1['ll2'].max()), ND)}
cls = dict(T1='veh', T2='veh', T5='veh', C1='room')
res = {}
for ev in ('T1','T2','T5','C1'):
    k1, n1, k2, n2 = E.OBS[ev]; K = k1+k2
    h = stats.hypergeom.pmf(np.arange(n1+1), n1+n2, n1, K)
    res[ev] = {'M0': h}
    for m in ('M1','M2'):
        pmf = np.zeros(n1+1); rr = []; epsl = []; quanta = []
        for i in range(ND):
            cfg = E.B[ev](rng); D = post[(cls[ev], m)][i]
            X = np.clip(expo(cfg['lat'], m, D, cfg['kappa'], cfg['segs'], cfg['src'], cfg['rec']), 1e-300, None)
            eps = prof_eps(X, K); p = -np.expm1(-eps*X)
            pmf += cond_pmf(p[cfg['near']], p[~cfg['near']], K)
            rr.append(p[cfg['near']].mean()/p[~cfg['near']].mean()); epsl.append(eps)
        pmf /= ND; res[ev][m] = pmf
        res[ev][m+'_rr'] = np.quantile(rr, [.05,.5,.95]); res[ev][m+'_eps'] = np.quantile(epsl, [.05,.5,.95])
    print(ev, 'obs', E.OBS[ev])
    for m in ('M0','M1','M2'):
        pm = res[ev][m]; ek = (np.arange(n1+1)*pm).sum()
        print(f'   {m}: E[k_near]={ek:.2f} modelRR(from Ek)={(ek/n1)/((K-ek)/n2):.2f} p2={two_sided(pm,k1):.4f} logS={np.log(pm[k1]):.3f} delta={np.log(pm[k1])-np.log(res[ev]["M0"][k1]):+.3f}',
              '' if m=='M0' else f'RRq={np.round(res[ev][m+"_rr"],2)} eps q={res[ev][m+"_eps"]}')
# pooled exact test
def pooled(m):
    evs = ['T1','T2','T5','C1']
    pm = [np.clip(res[e][m],1e-300,None) for e in evs]; p0 = [np.clip(res[e]['M0'],1e-300,None) for e in evs]
    g = np.meshgrid(*[np.arange(len(p)) for p in p0], indexing='ij')
    d = sum(np.log(pm[a][g[a]])-np.log(p0[a][g[a]]) for a in range(4)); w = np.prod([p0[a][g[a]] for a in range(4)], axis=0)
    dobs = sum(np.log(pm[a][E.OBS[e][0]])-np.log(p0[a][E.OBS[e][0]]) for a,e in enumerate(evs))
    return dobs, w[d >= dobs-1e-9].sum()/w.sum()
for m in ('M1','M2'): print('pooled', m, pooled(m))
np.savez('res_holdout.npz', **{f'{e}_{m}': res[e][m] for e in res for m in ('M0','M1','M2')})
