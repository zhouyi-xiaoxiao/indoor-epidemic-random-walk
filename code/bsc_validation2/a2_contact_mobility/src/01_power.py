"""Power analysis on SYNTHETIC data only (no real data read). For each synthetic 'truth' the frozen
pipeline is run and every candidate model is scored by the pre-registered rule."""
import sys, json, numpy as np, time
import a2lib as A, pipeline as PL

def synth_presence(rng, N, ndays):
    pres = []
    for d in range(ndays):
        ids = np.flatnonzero(rng.random(N) < 0.9)
        a = rng.integers(1440, 1620, len(ids)); b = rng.integers(2880, 3060, len(ids))
        pres.append(dict(ids=ids, a=a, b=b))
    return pres

def syn_structured(rng, pres, N, groups, c, ndays):
    """Out-of-family truth: fixed lognormal pair affinities x group factor, heavy-tailed (discrete Pareto) durations."""
    iu, ju = np.triu_indices(N, 1)
    aff = np.exp(rng.normal(0, 1.5, len(iu))) * np.where(groups[iu] == groups[ju], 8.0, 1.0)
    aff = aff / aff.mean()
    A_ = np.zeros((N, N)); A_[iu, ju] = aff
    out = []
    for d in range(ndays):
        pr = pres[d]; ids, a, b = pr['ids'], pr['a'], pr['b']; n = len(ids)
        i2, j2 = np.triu_indices(n, 1)
        lo = np.maximum(a[i2], a[j2]); hi = np.minimum(b[i2], b[j2]); L = np.clip(hi - lo + 1, 0, None)
        du_mean = 3.0
        ne = rng.poisson(c * A_[ids[i2], ids[j2]] * L / du_mean)
        rep = np.repeat(np.arange(len(i2)), ne)
        st = lo[rep] + np.floor(rng.random(len(rep)) * L[rep]).astype(np.int64)
        du = np.minimum(np.floor(rng.pareto(1.5, len(rep)) + 1).astype(np.int64), 400)   # mean ~3
        en = np.minimum(st + du, lo[rep] + L[rep])
        ln = en - st; tot = int(ln.sum())
        S = np.repeat(st, ln) + (np.arange(tot) - np.repeat(np.cumsum(ln) - ln, ln))
        I = np.repeat(ids[i2[rep]], ln); J = np.repeat(ids[j2[rep]], ln)
        k = np.unique(S * (N * N) + I * N + J)
        out.append(dict(s=k // (N * N), i=(k // N) % N, j=k % N))
    return out

def make_truth(kind, seed, N=120, ndays=4, c=6e-4):
    rng = np.random.default_rng([seed, 777])
    groups = np.arange(N) % 4
    pres = synth_presence(rng, N, ndays)
    cal = dict(c=c, c_in=c * 8 / (1 + 7 * 0.24), c_out=c / (1 + 7 * 0.24), mean_dur=2.5, M=1 / c, p=A.p_from_q(1 - 1 / 2.5))
    alld = list(range(ndays))
    if kind == 'SYN':
        days = syn_structured(rng, pres, N, groups, c, ndays)
    else:
        g = A.generate(kind, rng, pres, N, groups, cal, alld)
        days = [g[d] for d in alld]
    return dict(N=N, groups=groups, days=days)

if __name__ == '__main__':
    nrep = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    models = ['WM', 'BLOCK', 'RW0', 'RWhet']
    out = {}
    for kind in ['WM', 'BLOCK', 'RW0', 'SYN']:
        out[kind] = []
        for r in range(nrep):
            t0 = time.time()
            ds = make_truth(kind, r)
            res = PL.run_dataset(ds, models, G=4, G_epi=10, nrep_real=2000, nrep_model=200, seed=100 + r, verbose=False)
            out[kind].append(res)
            print(kind, r, f'{time.time()-t0:.0f}s', {m: {k: v for k, v in res['models'][m]['score'].items()} for m in models}, flush=True)
            PL.dump(out, A.ROOT + '/out/01_power.json')
