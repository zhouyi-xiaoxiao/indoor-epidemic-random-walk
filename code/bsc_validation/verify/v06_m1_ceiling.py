"""The re-check's own event simulator for the per-cell rule in the limit beta -> infinity:
all persons walk independently on the lattice (hop rate D/ax^2 across, D/ay^2 along, reflecting walls);
a susceptible is infected the first time it shares a site with the index. Seats reset per segment.
Also a finite-beta version (hazard beta * 1/n_total on the shared site; q=1) for a sanity check."""
import numpy as np, sys
import vevents as E
def sim(cfg, D, reps, rng, beta=np.inf):
    lat = cfg['lat']; nx, ny = lat.nx, lat.ny
    wx, wy = D/lat.ax**2, D/lat.ay**2; R = 2*wx+2*wy
    n = len(cfg['rec']); tot = np.zeros(reps, int); inf_all = np.zeros((reps, n), bool)
    for rep in range(reps):
        inf = np.zeros(n, bool)
        for T in cfg['segs']:
            pos = np.concatenate([[cfg['src']], cfg['rec']]).copy()
            ne = rng.poisson(R*T*(n+1))
            who = rng.integers(0, n+1, ne); u = rng.random(ne); times = np.sort(rng.random(ne))*T
            tprev = 0.0
            for e in range(ne):
                if np.isfinite(beta):
                    # exposure accumulates during [tprev, t) for those co-located with index
                    co = np.flatnonzero((pos[1:] == pos[0]) & ~inf)
                    if len(co):
                        ntot = (pos == pos[0]).sum()
                        pinf = -np.expm1(-beta/ntot*(times[e]-tprev))
                        inf[co[rng.random(len(co)) < pinf]] = True
                    tprev = times[e]
                w = who[e]; s = pos[w]; x, y = s % nx, s // nx; q = u[e]*R
                if q < wx: x2, y2 = x+1, y
                elif q < 2*wx: x2, y2 = x-1, y
                elif q < 2*wx+wy: x2, y2 = x, y+1
                else: x2, y2 = x, y-1
                if 0 <= x2 < nx and 0 <= y2 < ny:
                    pos[w] = y2*nx+x2
                    if not np.isfinite(beta):
                        if w == 0: inf |= (pos[1:] == pos[0])
                        elif pos[w] == pos[0]: inf[w-1] = True
            if np.isfinite(beta):
                co = np.flatnonzero((pos[1:] == pos[0]) & ~inf)
                if len(co):
                    ntot = (pos == pos[0]).sum(); pinf = -np.expm1(-beta/ntot*(T-tprev)); inf[co[rng.random(len(co)) < pinf]] = True
        tot[rep] = inf.sum(); inf_all[rep] = inf
    return tot, inf_all
rng = np.random.default_rng(99)
for ev, Ds in (('T1',[0.31,0.5,0.79]), ('T2',[0.31,0.5,0.79]), ('T5',[0.31,0.5,0.79]), ('C1',[0.2,1.0,14.0])):
    for D in Ds:
        tots = []
        for c in range(6):
            cfg = E.B[ev](rng); t,_ = sim(cfg, D, 60, rng); tots.append(t)
        t = np.concatenate(tots); K = cfg['K']
        print(f'{ev} D={D}: beta=inf mean total {t.mean():.2f} (sd {t.std():.2f}), max {t.max()}, P(total>={K}) = {(t>=K).mean():.4f}  [observed {K}/{len(cfg["rec"])}]', flush=True)
