import numpy as np, pandas as pd
from vlat import *
rng = np.random.default_rng(31337)
df = pd.read_csv('../../bsc_outbreaks/data/park2020_callcentre_seats_digitized.csv')
nw = df[df.region=='north_wing'].reset_index(drop=True)
xy = nw[['x_pitch_units','y_pitch_units']].to_numpy(); case = nw.case.to_numpy().astype(int)
print('desks', len(nw), 'cases', case.sum(), ' south:', (df.region=='south_wing').sum(), df[df.region=='south_wing'].case.sum())
d = np.hypot(xy[:,None,0]-xy[None,:,0], xy[:,None,1]-xy[None,:,1]); A = ((d<=1.25)&(d>0)).astype(float)
def jz(c):
    return float(c@A@c)/2
# permutation moments by simulation (independent of analytic formula)
perm = np.array([jz(rng.permutation(case).astype(float)) for _ in range(20000)])
jobs = jz(case.astype(float)); m, s = perm.mean(), perm.std()
print('JJ obs', jobs, 'perm mean', m, 'sd', s, 'z', (jobs-m)/s, 'pairs', A.sum()/2)
def zstat(c, cache={}):
    k = int(c.sum())
    if k not in cache:
        lab = np.zeros(len(c)); lab[:k] = 1
        pp = np.array([jz(rng.permutation(lab)) for _ in range(3000)]); cache[k] = (pp.mean(), pp.std())
    return (jz(c.astype(float)) - cache[k][0])/cache[k][1]
zobs = (jobs-m)/s
def percolate(Xe, R):
    n = Xe.shape[0]; inf = np.zeros((R,n),bool); new = np.zeros((R,n),bool)
    new[np.arange(R), rng.integers(n,size=R)] = True; inf |= new
    while new.any():
        p = -np.expm1(-(new.astype(float)@Xe)); nxt = (rng.random((R,n))<p)&~inf; inf |= nxt; new = nxt
    return inf
def tune(X):
    lo = np.log(1/X.sum(axis=1).mean())-2; hi = lo+14
    for it in range(14):
        mid = (lo+hi)/2; fs = percolate(np.exp(mid)*X, 400).sum(axis=1); big = fs[fs>=40]
        mm = big.mean() if len(big)>=5 else 0
        if mm < 79: lo = mid
        else: hi = mid
    return np.exp((lo+hi)/2)
r1 = np.load('res_R1_primary.npz'); logD = r1['logD']; w = np.exp(r1['ll2']-r1['ll2'].max()); w/=w.sum()
x0, y0 = xy[:,0].min()-1, xy[:,1].min()-1
ix = np.rint(xy[:,0]-x0).astype(int); iy = np.rint(xy[:,1]-y0).astype(int)
print('lattice', ix.max()+2, iy.max()+2, 'distinct sites', len(set(zip(ix,iy))), 'of', len(ix))
res = []
ND = 40
allz = []
for i in range(ND):
    D = 10**(logD[rng.choice(61,p=w)] + rng.uniform(-0.05,0.05)); pitch = rng.uniform(1.0,1.6)
    kap = np.exp(rng.uniform(np.log(1),np.log(6)))+0.93
    lat = Lat(int(ix.max()+2), int(iy.max()+2), pitch, pitch)
    sites = lat.site(ix, iy); P = lat.V[sites]; X = np.clip((P*(1/(kap+D*lat.mu)))@P.T, 0, None); np.fill_diagonal(X, 0)
    eps = tune(X); inf = percolate(eps*X, 1500); fs = inf.sum(axis=1); keep = (fs>=69)&(fs<=89)
    z = np.array([zstat(c) for c in inf[keep]]); allz.append(z)
    res.append((D, kap, pitch, keep.sum(), np.median(z) if len(z) else np.nan, (z<=zobs).mean() if len(z) else np.nan))
    print(i, 'D=%.3g kap=%.2f pitch=%.2f keep=%d medz=%.2f P(z<=obs)=%.3f' % res[-1], flush=True)
zz = np.concatenate(allz); fr = np.array([r[5] for r in res if r[3]>0])
print('M2 pooled P(z<=obs)=', (zz<=zobs).mean(), 'q2.5/50/97.5', np.quantile(zz,[.025,.5,.975]), '| equal weight', fr.mean(), '| draws with P>0.025:', (fr>0.025).sum(), 'of', len(fr))
# M0
X = np.ones((137,137)); np.fill_diagonal(X,0); eps = tune(X); inf = percolate(eps*X, 20000); fs = inf.sum(axis=1); keep=(fs>=69)&(fs<=89)
z0 = np.array([zstat(c) for c in inf[keep][:4000]]); print('M0 P(z<=obs)', (z0<=zobs).mean(), np.quantile(z0,[.025,.5,.975]))
