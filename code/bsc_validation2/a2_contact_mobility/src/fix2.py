"""F2: home-anchored lattice walk with clustered homes (RWclus). Wraps a2lib.generate (frozen file untouched).
Parameters (M, p, b, sigma) fitted on TRAIN days only."""
import numpy as np
import a2lib as A
import fix1

BGRID = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
SGRID = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
_orig_generate = A.generate


def clustered_homes(rng, N, groups, Lx, Ly, sigma):
    if groups is None:
        return rng.integers(0, Lx, N), rng.integers(0, Ly, N)
    G = int(groups.max()) + 1
    gx = int(np.ceil(np.sqrt(G))); gy = int(np.ceil(G / gx))
    cx = (np.arange(G) % gx + 0.5) * Lx / gx; cy = (np.arange(G) // gx + 0.5) * Ly / gy
    hx = np.clip(np.floor(cx[groups] + sigma * rng.normal(size=N)), 0, Lx - 1).astype(np.int64)
    hy = np.clip(np.floor(cy[groups] + sigma * rng.normal(size=N)), 0, Ly - 1).astype(np.int64)
    return hx, hy


def generate(model, rng, pres, N, groups, cal, day_list, scale=1.0, extra=None):
    if model != 'RWclus':
        return _orig_generate(model, rng, pres, N, groups, cal, day_list, scale=scale, extra=extra)
    Lx, Ly = A.lattice_dims(extra['M'] / scale)
    home = clustered_homes(np.random.default_rng(int(rng.integers(1 << 30))), N, groups, Lx, Ly, extra['sigma'])
    return {d: A.rw_day(rng, pres[d], Lx, Ly, extra['p_fast'], home=home, bias=extra['bias'], start='home') for d in day_list}


A.generate = generate


def calfun(ds, pres, cal, train, seed):
    best = None; trace = []
    sgrid = SGRID if ds['groups'] is not None else [0.0]
    for b in BGRID:
        for s in sgrid:
            par, st, obs = fix1.fit_b(ds, pres, cal, train, seed, b, base=dict(sigma=s), model='RWclus')
            err = abs(np.log(st['S3_deg_mean'] / obs['S3_deg_mean']))
            if ds['groups'] is not None:
                err += abs(np.log(st['S6_within_frac'] / obs['S6_within_frac']))
            trace.append(dict(b=b, sigma=s, M=par['M'], p=par['p_fast'], deg=st['S3_deg_mean'], within=st['S6_within_frac'],
                              tot_ratio=st['S0_total'] / obs['S0_total'], dur=st['S1_mean_dur'], err=float(err)))
            if best is None or err < best[0]:
                best = (err, par)
    out = dict(best[1]); out['trace'] = trace
    return {'RWclus': out}
