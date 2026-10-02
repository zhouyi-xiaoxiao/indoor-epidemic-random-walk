import numpy as np
from scipy import stats
from vlat import *
import vevents as E
rng = np.random.default_rng(5)
t4 = np.load('res_T4.npz'); r1 = np.load('res_R1_primary.npz'); logD = t4['logD']
# --- C1: per-D conditional likelihood, implied emission, share of posterior with absurd emission
k1,n1,k2,n2 = E.OBS['C1']; K = 12
cfgs = [E.C1(rng) for _ in range(60)]
wr = np.exp(r1['ll2']-r1['ll2'].max()); wr /= wr.sum()
rows = []
for lD in logD:
    pm = np.zeros(n1+1); Em = []
    for cfg in cfgs:
        X = np.clip(expo(cfg['lat'], 'M2', 10**lD, cfg['kappa'], cfg['segs'], cfg['src'], cfg['rec']), 1e-300, None)
        eps = prof_eps(X, K); p = -np.expm1(-eps*X); pm += cond_pmf(p[cfg['near']], p[~cfg['near']], K)/len(cfgs); Em.append(eps*3.0/0.5)
    ek = (np.arange(n1+1)*pm).sum(); rows.append((lD, np.median(Em), (ek/n1)/((K-ek)/n2), pm[k1], two_sided(pm,k1)))
rows = np.array(rows)
for r in rows[::4]: print('C1 D=%8.3g  median implied E=%10.3g q/h  model RR=%5.2f  pmf(obs)=%.3f  p=%.3f' % (10**r[0], r[1], r[2], r[3], r[4]))
for thr in (1e3, 1e4, 1e5):
    print(f'posterior mass (R1-calibrated D_air) with C1 median implied emission > {thr:g} q/h: {wr[rows[:,1]>thr].sum():.3f}')
ok = rows[:,1] <= 1e4; w2 = wr*ok; w2 /= w2.sum()
pmf_ok = None
print('restricted to implied E <= 1e4 q/h: D >=', 10**logD[ok].min())
# --- vehicle class: per-event conditional likelihood as function of D (kappa, seating marginalised) and common-D LR test
cl = {}
for ev in ('T1','T2','T5'):
    k1,n1,k2,n2 = E.OBS[ev]; K = k1+k2; cf = [E.B[ev](rng) for _ in range(120)]; L = []
    for lD in logD[::2]:
        pm = 0.0
        for cfg in cf:
            X = np.clip(expo(cfg['lat'], 'M2', 10**lD, cfg['kappa'], cfg['segs'], cfg['src'], cfg['rec']), 1e-300, None)
            eps = prof_eps(X, K); p = -np.expm1(-eps*X); pm += cond_pmf(p[cfg['near']], p[~cfg['near']], K)[k1]/len(cf)
        L.append(np.log(pm))
    cl[ev] = np.array(L); print(ev, 'D preferred alone', 10**logD[::2][cl[ev].argmax()], 'max logL', cl[ev].max())
ll = t4['ll2']; m = ll.max(); T4m = (m + np.log(np.exp(ll-m).mean(axis=1)))[::2]
print('T4 preferred', 10**logD[::2][T4m.argmax()], T4m.max())
tot = T4m + cl['T1'] + cl['T2'] + cl['T5']; sep = T4m.max()+cl['T1'].max()+cl['T2'].max()+cl['T5'].max()
print('sum of maxima', sep, 'common-D max', tot.max(), 'at D', 10**logD[::2][tot.argmax()], 'LR', 2*(sep-tot.max()), 'p(3 df)', stats.chi2.sf(2*(sep-tot.max()),3))
