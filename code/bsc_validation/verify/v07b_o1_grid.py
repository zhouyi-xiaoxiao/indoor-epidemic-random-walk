import numpy as np, pandas as pd, json
from vlat import *
rng = np.random.default_rng(2718)
df = pd.read_csv('../../bsc_outbreaks/data/park2020_callcentre_seats_digitized.csv')
nw = df[df.region=='north_wing'].reset_index(drop=True)
xy = nw[['x_pitch_units','y_pitch_units']].to_numpy(); case = nw.case.to_numpy().astype(float)
d = np.hypot(xy[:,None,0]-xy[None,:,0], xy[:,None,1]-xy[None,:,1]); A = ((d<=1.25)&(d>0)).astype(float)
def moments(k, cache={}):
    if k not in cache:
        lab = np.zeros(137); lab[:k] = 1
        pp = np.array([(lambda c: c@A@c/2)(rng.permutation(lab)) for _ in range(4000)]); cache[k] = (pp.mean(), pp.std())
    return cache[k]
def z(c):
    m, s = moments(int(c.sum())); return (c@A@c/2 - m)/s
zobs = z(case); print('zobs', zobs)
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
x0, y0 = xy[:,0].min()-1, xy[:,1].min()-1
ix = np.rint(xy[:,0]-x0).astype(int); iy = np.rint(xy[:,1]-y0).astype(int)
grid = np.round(np.arange(-2, 4.01, 0.25), 2); out = {}
for lD in grid:
    P = []; KF = []
    for rep in range(3):
        D = 10**lD; pitch = rng.uniform(1.0,1.6); kap = np.exp(rng.uniform(np.log(1),np.log(6)))+0.93
        lat = Lat(int(ix.max()+2), int(iy.max()+2), pitch, pitch); sites = lat.site(ix, iy); Pm = lat.V[sites]
        X = np.clip((Pm*(1/(kap+D*lat.mu)))@Pm.T, 0, None); np.fill_diagonal(X, 0)
        eps = tune(X); inf = percolate(eps*X, 1500); fs = inf.sum(axis=1); keep = (fs>=69)&(fs<=89)
        zz = np.array([z(c.astype(float)) for c in inf[keep]])
        P.append((zz<=zobs).mean() if len(zz) else 0.0); KF.append(keep.mean())
    out[float(lD)] = (float(np.mean(P)), float(np.mean(KF)))
    print(lD, out[float(lD)], flush=True)
json.dump(out, open('res_o1_grid.json','w'))
for name, f in (('re-check R1 fit', 'res_R1_primary.npz'),):
    r1 = np.load(f); logD = r1['logD']; w = np.exp(r1['ll2']-r1['ll2'].max()); w/=w.sum()
    Pg = np.interp(logD, grid, [out[float(g)][0] for g in grid]); Kg = np.interp(logD, grid, [out[float(g)][1] for g in grid])
    print(name, 'equal-weight P(z<=obs) =', (w*Pg).sum(), ' pooled-weight =', (w*Kg*Pg).sum()/(w*Kg).sum())
c = json.load(open('../data/01_cal_R1_primary.json')); ll = np.array(c['M2']['profile_logL']); w = np.exp(ll-ll.max()); w/=w.sum(); logD = np.array(c['logD'])
Pg = np.interp(logD, grid, [out[float(g)][0] for g in grid]); Kg = np.interp(logD, grid, [out[float(g)][1] for g in grid])
print('first-analysis R1 posterior: equal-weight P =', (w*Pg).sum(), ' pooled-weight =', (w*Kg*Pg).sum()/(w*Kg).sum())
# sensitivity to the upper limit of the flat prior on log10 D
for top in (2, 3, 4):
    m = logD <= top; ww = w[m]/w[m].sum()
    print('prior upper limit 10^%d: equal %.4f pooled %.4f' % (top, (ww*Pg[m]).sum(), (ww*Kg[m]*Pg[m]).sum()/(ww*Kg[m]).sum()))
# extend to 10^6 assuming flat profile beyond grid
ext = np.concatenate([w, np.full(20, w[-1])]); ext/=ext.sum(); Pe = np.concatenate([Pg, np.full(20, Pg[-1])]); Ke = np.concatenate([Kg, np.full(20, Kg[-1])])
print('prior upper limit 10^6: equal %.4f pooled %.4f' % ((ext*Pe).sum(), (ext*Ke*Pe).sum()/(ext*Ke).sum()))
