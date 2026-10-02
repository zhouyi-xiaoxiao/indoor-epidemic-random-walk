"""F3: clustered-home walk with slow home sites (RWstickZ: slow zone = all home sites, harmonic-mean interface rule;
RWstickO: a walker is slow only on its own home site). Wraps fix2 (frozen a2lib untouched).
Parameters (M, p, b, sigma, rho) fitted on TRAIN days only."""
import numpy as np
import a2lib as A
import fix1, fix2

RHO = [1.0, 0.3, 0.1, 0.03, 0.01]
_prev_generate = A.generate          # = fix2.generate


def rw_day_stick(rng, pr, Lx, Ly, p, rho, home, bias, variant):
    ids, a, b = pr['ids'], pr['a'], pr['b']; n = len(ids)
    T0 = int(a.min()); T1 = int(b.max()); T = T1 - T0 + 1
    hx = home[0][ids]; hy = home[1][ids]
    x = hx.copy(); y = hy.copy()
    slow = np.zeros(Lx * Ly, bool); slow[home[0] * Ly + home[1]] = True     # all individuals' homes
    hsite = hx * Ly + hy
    P = np.empty((T, n), dtype=np.int64)
    neg = -(np.arange(n) + 1)
    dxs = np.array([1, -1, 0, 0]); dys = np.array([0, 0, 1, -1])
    ps = rho * p
    for s in range(T0, T1 + 1):
        act = (a <= s) & (s <= b)
        mv = act & (a < s)
        d = rng.integers(0, 4, n)
        dx = dxs[d]; dy = dys[d]
        if bias > 0:
            hb = rng.random(n) < bias
            ex = hx - x; ey = hy - y
            usex = (ex != 0) & ((ey == 0) | (rng.random(n) < 0.5))
            bx = np.where(usex, np.sign(ex), 0); by = np.where(~usex, np.sign(ey), 0)
            dx = np.where(hb, bx, dx); dy = np.where(hb, by, dy)
        nx = x + dx; ny = y + dy
        ok = (nx >= 0) & (nx < Lx) & (ny >= 0) & (ny < Ly) & ((dx != 0) | (dy != 0))
        cur = x * Ly + y
        if variant == 'Z':
            nxt = np.clip(nx, 0, Lx - 1) * Ly + np.clip(ny, 0, Ly - 1)
            p0 = np.where(slow[cur], ps, p); p1 = np.where(slow[nxt], ps, p)
            acc = 2 * p0 * p1 / (p0 + p1)
        else:
            acc = np.where(cur == hsite, ps, p)
        go = mv & ok & (rng.random(n) < acc)
        x = np.where(go, nx, x); y = np.where(go, ny, y)
        P[s - T0] = np.where(act, x * Ly + y, neg)
    return A._contacts_from_positions(P, ids, T0)


def generate(model, rng, pres, N, groups, cal, day_list, scale=1.0, extra=None):
    if model not in ('RWstickZ', 'RWstickO'):
        return _prev_generate(model, rng, pres, N, groups, cal, day_list, scale=scale, extra=extra)
    Lx, Ly = A.lattice_dims(extra['M'] / scale)
    home = fix2.clustered_homes(np.random.default_rng(int(rng.integers(1 << 30))), N, groups, Lx, Ly, extra['sigma'])
    return {d: rw_day_stick(rng, pres[d], Lx, Ly, extra['p_fast'], extra['rho'], home, extra['bias'], model[-1]) for d in day_list}


A.generate = generate


def _err(st, obs, groups):
    e = abs(np.log(st['S3_deg_mean'] / obs['S3_deg_mean']))
    if groups is not None:
        e += abs(np.log(st['S6_within_frac'] / obs['S6_within_frac']))
    e += abs(np.log(max(st['S1_timefrac_ge15'], 0.002) / max(obs['S1_timefrac_ge15'], 0.002)))
    return float(e)


def _fit(ds, pres, cal, train, seed, model, b, s, rho, trace):
    par, st, obs = fix1.fit_b(ds, pres, cal, train, seed, b, base=dict(sigma=s, rho=rho), model=model)
    e = _err(st, obs, ds['groups'])
    trace.append(dict(b=b, sigma=s, rho=rho, M=par['M'], p=par['p_fast'], deg=st['S3_deg_mean'], within=st['S6_within_frac'],
                      long=st['S1_timefrac_ge15'], tot_ratio=st['S0_total'] / obs['S0_total'], dur=st['S1_mean_dur'], err=e))
    return e, par


def make_calfun(models, f2=None, refine=False):
    def calfun(ds, pres, cal, train, seed):
        out = {}
        f2par = (f2 or fix2.calfun)(ds, pres, cal, train, seed)['RWclus']
        out['RWclus'] = f2par
        b0, s0 = f2par['bias'], f2par['sigma']
        sgrid = fix2.SGRID if ds['groups'] is not None else [0.0]
        for model in models:
            trace = []; best = None
            for rho in RHO:                                    # stage B
                e, par = _fit(ds, pres, cal, train, seed, model, b0, s0, rho, trace)
                if best is None or e < best[0]: best = (e, par)
            rho = best[1]['rho']
            ib = fix2.BGRID.index(b0); isg = sgrid.index(s0)
            for b in fix2.BGRID[max(ib - 1, 0):ib + 2]:        # stage C
                for s in sgrid[max(isg - 1, 0):isg + 2]:
                    if b == b0 and s == s0: continue
                    e, par = _fit(ds, pres, cal, train, seed, model, b, s, rho, trace)
                    if e < best[0]: best = (e, par)
            if refine:                                         # stage D
                bb, ss, r0 = best[1]['bias'], best[1]['sigma'], best[1]['rho']
                for rho in (r0 * 3 ** 0.5, r0 / 3 ** 0.5):
                    if rho > 1: continue
                    e, par = _fit(ds, pres, cal, train, seed, model, bb, ss, rho, trace)
                    if e < best[0]: best = (e, par)
            o = dict(best[1]); o['fit_err'] = best[0]; o['trace'] = trace
            out[model] = o
        return out
    return calfun


def _f2_cached(tag):
    """Reuse an F2 fit already stored in out/<tag>_<ds>.json (identical deterministic procedure and seed)."""
    import json, os
    def f(ds, pres, cal, train, seed):
        for t in tag.split(','):
            p = f"{A.ROOT}/out/{t}_{ds['name']}.json"
            if os.path.exists(p):
                e = json.load(open(p))['extra']['RWclus']
                return {'RWclus': {k: v for k, v in e.items() if k != 'trace'}}
        return fix2.calfun(ds, pres, cal, train, seed)
    return f


calfun = make_calfun(['RWstickZ', 'RWstickO'], _f2_cached('s2b'))
calfun_refined = make_calfun(['RWstickZ', 'RWstickO'], _f2_cached('s2b'), refine=True)
calfun_Z = make_calfun(['RWstickZ'], _f2_cached('s2b,s3b'), refine=True)
calfun_O = make_calfun(['RWstickO'], _f2_cached('s2b,s3b'), refine=True)
