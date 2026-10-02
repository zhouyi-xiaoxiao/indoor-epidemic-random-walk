"""F1: home-anchored lattice walk (RWhome). Parameters (M, p, b) fitted on TRAIN days only by simulation:
for each bias b on a fixed grid, (M, p) are iterated so that simulated total contact time and mean event
duration on the training days equal the observed ones; b is then chosen to match the training mean number of
distinct contacts per person-day (S3_deg_mean)."""
import numpy as np
import a2lib as A

BGRID = [0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7]
NITER = 4


def _q(p): return (1 - p) ** 2 + p * p / 4


def _sim(ds, pres, cal, train, par, seed, model='RWhome'):
    rng = np.random.default_rng([seed, 4242])
    gen = A.generate(model, rng, pres, ds['N'], ds['groups'], cal, train, extra=par)
    days = [gen.get(d) for d in range(len(ds['days']))]
    return A.stats(days, pres, ds['N'], ds['groups'], train)


def fit_b(ds, pres, cal, train, seed, b, base=None, model='RWhome'):
    obs = A.stats(ds['days'], pres, ds['N'], ds['groups'], train)
    par = dict(base or {}); par.update(M=cal['M'], p_fast=cal['p'], bias=b, home_seed=seed)
    for it in range(NITER):
        st = _sim(ds, pres, cal, train, par, seed + it, model)
        par['M'] = float(np.clip(par['M'] * st['S0_total'] / obs['S0_total'], 4, 1e7))
        brk = (1 - _q(par['p_fast'])) * st['S1_mean_dur'] / obs['S1_mean_dur']
        par['p_fast'] = A.p_from_q(1 - np.clip(brk, 1e-4, 0.8))
    st = _sim(ds, pres, cal, train, par, seed + 99, model)
    return par, st, obs


def calfun(ds, pres, cal, train, seed):
    best = None; trace = []
    for b in BGRID:
        par, st, obs = fit_b(ds, pres, cal, train, seed, b)
        err = abs(np.log(st['S3_deg_mean'] / obs['S3_deg_mean']))
        trace.append(dict(b=b, M=par['M'], p=par['p_fast'], deg=st['S3_deg_mean'], tot_ratio=st['S0_total'] / obs['S0_total'],
                          dur=st['S1_mean_dur'], err=float(err)))
        if best is None or err < best[0]:
            best = (err, par)
    out = dict(best[1]); out['trace'] = trace
    return {'RWhome': out}
