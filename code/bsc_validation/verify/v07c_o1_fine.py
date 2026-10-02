import numpy as np, pandas as pd
from vlat import *
exec(open('v07b_o1_grid.py').read().split("x0, y0 =")[0].replace("rng = np.random.default_rng(2718)", "rng = np.random.default_rng(11)"))
for f in (1, 2, 3):
    x0, y0 = xy[:,0].min()-1, xy[:,1].min()-1
    ix = np.rint((xy[:,0]-x0)*f).astype(int); iy = np.rint((xy[:,1]-y0)*f).astype(int)
    print('refinement', f, 'distinct sites', len(set(zip(ix,iy))), 'of 137; lattice', ix.max()+1+f, iy.max()+1+f)
    for D in (3.0, 30.0, 100.0):
        P = []; 
        for rep in range(3):
            pitch = rng.uniform(1.0,1.6)/f; kap = np.exp(rng.uniform(np.log(1),np.log(6)))+0.93
            lat = Lat(int(ix.max()+1+f), int(iy.max()+1+f), pitch, pitch); sites = lat.site(ix, iy); Pm = lat.V[sites]
            X = np.clip((Pm*(1/(kap+D*lat.mu)))@Pm.T, 0, None); np.fill_diagonal(X, 0)
            eps = tune(X); inf = percolate(eps*X, 1500); fs = inf.sum(axis=1); keep = (fs>=69)&(fs<=89)
            zz = np.array([z(c.astype(float)) for c in inf[keep]]); P.append(((zz<=zobs).mean(), np.median(zz)))
        print('   D=%g: P(z<=obs) per rep, median z:' % D, np.round(P,3).tolist(), flush=True)
