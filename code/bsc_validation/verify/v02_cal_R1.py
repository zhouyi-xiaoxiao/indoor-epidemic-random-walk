import numpy as np, pandas as pd
from scipy import stats, optimize
from vlat import *
rng = np.random.default_rng(777)
# retyped by the re-check from Li et al. 2021 table as stored (table, patrons, overlap min); T04 empty, TA index table
tabs = [('TB',4,53),('TC',7,75),('T05',2,52),('T06',4,82),('T07',3,69),('T08',2,55),('T09',10,75),('T10',6,82),
        ('T11',7,70),('T12',2,64),('T13',6,50),('T14',3,61),('T15',8,82),('T16',5,48),('T17',5,23),('T18',5,77)]
chk = pd.read_csv('../../bsc_outbreaks/data/li2021_restaurant_tables.csv')
chk = chk[(chk.table!='TA')&(chk.patrons>0)]
assert [(a,int(b),int(c)) for a,b,c in zip(chk.table, chk.patrons, chk.overlap_with_table_A_min)] == tabs
print('patrons excluding A', sum(t[1] for t in tabs))
lat = Lat(17, 8, 1.0, 1.0)
XS = [1,4,7,10,12,15]; YS = [1,4,6]
names = [t[0] for t in tabs]; n = np.array([t[1] for t in tabs]); T = np.array([t[2] for t in tabs])/60
def layout():
    iA = rng.integers(1,5); side = rng.choice([-1,1])
    pos = {'TA':(iA,0),'TB':(iA+side,0),'TC':(iA-side,0),'T18':(iA,1)}
    free = [(i,j) for i in range(6) for j in range(3) if (i,j) not in pos.values()]
    rem = [t for t in names if t not in pos] + ['T04']
    perm = rng.permutation(len(free))
    for t,kk in zip(rem, perm): pos[t] = free[kk]
    s = np.array([lat.site(XS[pos[t][0]], YS[pos[t][1]]) for t in names])
    return dict(src=int(lat.site(XS[iA], YS[0])), sites=s, kappa=rng.uniform(0.56,0.77)+0.93)
lays = [layout() for _ in range(300)]
iB, iC = 0, 1; oth = np.arange(2,16)
def Xall(model, D):
    X = np.zeros((300,16))
    for a,ly in enumerate(lays):
        for t in np.unique(T):
            sel = T==t
            X[a,sel] = expo(lat, model, D, ly['kappa'], [t], ly['src'], ly['sites'][sel])
    return X
def loglik(X, eps, kb=(1,2,3), kc=(1,2)):
    p = np.clip(-np.expm1(-eps*X), 0, 1-1e-15)
    ll = (n[oth]*np.log1p(-p[:,oth])).sum(axis=1)
    PB = sum(stats.binom.pmf(q, n[iB], p[:,iB]) for q in kb); PC = sum(stats.binom.pmf(q, n[iC], p[:,iC]) for q in kc)
    ll = ll + np.log(np.clip(PB,1e-300,None)) + np.log(np.clip(PC,1e-300,None))
    m = ll.max(); return m + np.log(np.mean(np.exp(ll-m)))
def profile(X, **kw):
    f = lambda le: -loglik(X, np.exp(le), **kw)
    les = np.linspace(-10, 60, 141); v = [f(l) for l in les]; i = int(np.argmin(v))
    r = optimize.minimize_scalar(f, bounds=(les[max(i-1,0)], les[min(i+1,140)]), method='bounded', options={'xatol':1e-7})
    return -r.fun, np.exp(r.x)
logD = np.linspace(-2,4,61)
out = {}
for var, kw in (('primary',{}), ('upper',dict(kb=(3,),kc=(2,))), ('lower',dict(kb=(1,),kc=(1,)))):
    l0,e0 = profile(Xall('M0',1.0), **kw)
    r = {}
    for model in ('M1','M2'):
        ll = np.zeros(61); ee = np.zeros(61)
        for i,D in enumerate(10**logD):
            ll[i], ee[i] = profile(Xall(model, D), **kw)
        r[model] = (ll, ee)
        w = np.exp(ll-ll.max()); c = np.cumsum(w/w.sum())
        print(var, model, 'M0', round(l0,3), 'max', round(ll.max(),3), 'D', 10**logD[ll.argmax()], 'eps', ee[ll.argmax()],
              'q5/50/95', [round(10**np.interp(a,c,logD),3) for a in (.05,.5,.95)])
    out[var] = r
    np.savez(f'res_R1_{var}.npz', logD=logD, ll1=r['M1'][0], e1=r['M1'][1], ll2=r['M2'][0], e2=r['M2'][1], l0=l0, e0=e0)
