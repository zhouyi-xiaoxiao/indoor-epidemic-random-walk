import numpy as np, json
from scipy import stats
from vlat import *
import vevents as E
t4 = np.load('res_T4.npz'); r1 = np.load('res_R1_primary.npz'); logD = t4['logD']; step = logD[1]-logD[0]
evs = ['T1','T2','T5','C1']; cls = dict(T1='veh', T2='veh', T5='veh', C1='room')
VOL = dict(T1=2.5/6*0.75*2.0, T2=0.5*(11.4/13)*(60.42/(11.4*2.5)), T5=0.9*1.1*2.2, C1=3.0)
def run(wveh, wroom, lg_veh=logD, lg_room=logD, ND=400, seed=1, kmul=None, label=''):
    rng = np.random.default_rng(seed)
    def dr(w, lg):
        w = w/w.sum(); i = rng.choice(len(w), ND, p=w); st = lg[1]-lg[0]; return 10**(lg[i]+rng.uniform(-st/2, st/2, ND))
    post = dict(veh=dr(wveh, lg_veh), room=dr(wroom, lg_room))
    pms, p0s, out = [], [], []
    for ev in evs:
        k1,n1,k2,n2 = E.OBS[ev]; K = k1+k2; pmf = np.zeros(n1+1); Eq = []
        for i in range(ND):
            cfg = E.B[ev](rng); kap = cfg['kappa']
            if kmul and ev in kmul: kap = (kap-0.93)*kmul[ev]+0.93
            X = np.clip(expo(cfg['lat'], 'M2', post[cls[ev]][i], kap, cfg['segs'], cfg['src'], cfg['rec']), 1e-300, None)
            eps = prof_eps(X, K); p = -np.expm1(-eps*X); pmf += cond_pmf(p[cfg['near']], p[~cfg['near']], K); Eq.append(eps*VOL[ev]/0.5)
        pmf /= ND; pms.append(pmf); p0s.append(stats.hypergeom.pmf(np.arange(n1+1), n1+n2, n1, K))
        ek = (np.arange(n1+1)*pmf).sum()
        out.append((ev, round((ek/n1)/((K-ek)/n2),2), round(two_sided(pmf,k1),4), np.quantile(Eq,[.05,.5,.95])))
    g = np.meshgrid(*[np.arange(len(p)) for p in p0s], indexing='ij')
    pm = [np.clip(p,1e-300,None) for p in pms]; p0 = [np.clip(p,1e-300,None) for p in p0s]
    d = sum(np.log(pm[a][g[a]])-np.log(p0[a][g[a]]) for a in range(4)); w = np.prod([p0[a][g[a]] for a in range(4)], axis=0)
    dobs = sum(np.log(pm[a][E.OBS[e][0]])-np.log(p0[a][E.OBS[e][0]]) for a,e in enumerate(evs))
    print(f'[{label}] pooled delta {dobs:+.2f} p {w[d>=dobs-1e-9].sum()/w.sum():.5f} | ' + ' | '.join(f'{e}: RR {r}, p {p}' for e,r,p,_ in out))
    return out
wv = np.exp(t4['ll2']-t4['ll2'].max()).sum(axis=1); wr = np.exp(r1['ll2']-r1['ll2'].max())
o = run(wv, wr, label='primary (verifier)')
print('  implied emission quanta/h under M2, 5/50/95%:'); [print('    ', e, np.array2string(q, precision=3)) for e,_,_,q in o]
# overdispersion-tempered train posterior (phi = deviance/df = 72.63/20)
phi = 72.63/20; wv_t = np.exp((t4['ll2']-t4['ll2'].max())/phi).sum(axis=1)
c = np.cumsum(wv_t/wv_t.sum()); print('tempered D_air 5/50/95', [round(10**np.interp(a,c,logD),1) for a in (.05,.5,.95)])
run(wv_t, wr, label='T4 posterior tempered by overdispersion 3.6')
# kappa prior sensitivity (ventilation ACH multiplied)
for km in (0.33, 3.0):
    run(wv, wr, kmul=dict(T1=km, T5=km, C1=km), label=f'ACH priors of T1,T5,C1 x{km}')
run(wv, wr, kmul=dict(T2=0.33), label='T2 ACH x0.33')
# S-adj / S-row
for name in ('S-adj','S-row'):
    l = np.load(f'res_T4_{name}_ll2.npy'); run(np.exp(l-l.max()).sum(axis=1), wr, lg_veh=logD[::2], label=name)
