import numpy as np
from vlat import *
import vevents as E
rng = np.random.default_rng(8)
t4 = np.load('res_T4.npz'); r1 = np.load('res_R1_primary.npz'); logD = t4['logD']
wr = np.exp(r1['ll2']-r1['ll2'].max()); wr /= wr.sum()
wv = np.exp(t4['ll2']-t4['ll2'].max()).sum(axis=1); wv /= wv.sum()
def curve(ev, n=80):
    k1,n1,k2,n2 = E.OBS[ev]; K = k1+k2; cf = [E.B[ev](rng) for _ in range(n)]; out = []
    for lD in logD:
        pm = np.zeros(n1+1)
        for cfg in cf:
            X = np.clip(expo(cfg['lat'], 'M2', 10**lD, cfg['kappa'], cfg['segs'], cfg['src'], cfg['rec']), 1e-300, None)
            eps = prof_eps(X, K); p = -np.expm1(-eps*X); pm += cond_pmf(p[cfg['near']], p[~cfg['near']], K)/n
        out.append(pm)
    return np.array(out)
c = curve('C1'); k1 = 8
tot = (wr[:,None]*c).sum(axis=0)
print('C1 mixture pmf(obs)=%.4f  two-sided p=%.3f' % (tot[k1], two_sided(tot,k1)))
for lo, hi in ((-2.01,-1.0),(-1.0,0.75),(0.75,4.01)):
    m = (logD>lo)&(logD<=hi)
    print('  D in (%.3g, %.3g]: posterior mass %.3f, contributes %.4f of pmf(obs) (%.0f%%); pmf(obs|range)=%.4f' % (10**lo, 10**hi, wr[m].sum(), (wr[m]*c[m,k1]).sum(), 100*(wr[m]*c[m,k1]).sum()/tot[k1], (wr[m]*c[m,k1]).sum()/wr[m].sum()))
im = wr.argmax(); print('  at posterior mode D=%.3g: pmf(obs)=%.4f, two-sided p=%.4f, E[k_near]=%.2f' % (10**logD[im], c[im,k1], two_sided(c[im],k1), (np.arange(11)*c[im]).sum()))
# with the vehicle-class (train) posterior instead (post hoc, exploratory)
tv = (wv[:,None]*c).sum(axis=0); print('C1 with the train-calibrated D_air (exploratory, post hoc): pmf(obs)=%.3f p=%.3f E[k_near]=%.2f' % (tv[k1], two_sided(tv,k1), (np.arange(11)*tv).sum()))
for ev in ('T5','T1','T2'):
    c5 = curve(ev, 60); k = E.OBS[ev][0]; im = wv.argmax()
    print(ev, 'at train posterior mode D=%.3g: pmf(obs)=%.4f p=%.4f' % (10**logD[im], c5[im,k], two_sided(c5[im],k)), '| at 5%%/95%% D: p=%.3f / %.3f' % (two_sided(c5[np.argmin(abs(logD-np.log10(9.7)))],k), two_sided(c5[np.argmin(abs(logD-np.log10(43)))],k)))
