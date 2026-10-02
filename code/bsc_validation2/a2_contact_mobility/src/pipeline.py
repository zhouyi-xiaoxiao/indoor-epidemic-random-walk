"""Frozen analysis pipeline: calibrate on train days, predict test days, score, epidemic comparison."""
import json, time, numpy as np
import a2lib as A

STAT_KEYS = ['S0_total', 'S1_mean_dur', 'S1_frac_ge3', 'S1_timefrac_ge15', 'S2_burst', 'S2_median_gap', 'S2_recurrence',
             'S3_deg_mean', 'S3_deg_cv', 'S4_strength_cv', 'S5_persistence', 'S6_within_frac', 'S7_alpha', 'S7_n_range']


def _avg(list_of_stats):
    return {k: float(np.nanmean([s.get(k, np.nan) for s in list_of_stats])) for k in STAT_KEYS}


def run_dataset(ds, models, G=5, G_epi=10, nrep_real=2000, nrep_model=200, seed=1, extra=None, calfun=None, verbose=True, epi=True):
    """ds: dict(N, groups, days). Returns a JSON-able result dict. `extra[model]` = parameters of fix models
    (already fitted on train days by `calfun`)."""
    t0 = time.time()
    N, groups, days = ds['N'], ds['groups'], ds['days']
    pres = [A.presence(d) for d in days]
    train, test = A.split(len(days))
    cal = A.calibrate(days, pres, N, groups, train)
    obs = A.stats(days, pres, N, groups, test)
    obs_train = A.stats(days, pres, N, groups, train)
    res = dict(N=N, n_days=len(days), train=train, test=test, cal=cal, obs=obs, obs_train=obs_train,
               n_groups=(int(groups.max()) + 1 if groups is not None else 0), models={}, extra={})
    res['models']['REALtrain'] = dict(stats=obs_train, score=A.score(obs_train, obs))
    if calfun is not None:
        extra = calfun(ds, pres, cal, train, seed)
        res['extra'] = extra
    extra = extra or {}
    alld = list(range(len(days)))
    for mi, m in enumerate(models):
        rng = np.random.default_rng([seed, mi, 11])
        sts = []
        for g in range(G):
            gen = A.generate(m, rng, pres, N, groups, cal, alld, extra=extra.get(m))
            sts.append(A.stats([gen[d] for d in alld], pres, N, groups, test))
        st = _avg(sts)
        res['models'][m] = dict(stats=st, score=A.score(st, obs),
                                stats_sd={k: float(np.nanstd([s.get(k, np.nan) for s in sts])) for k in STAT_KEYS})
        if verbose: print(f"   {m}: stats done {time.time()-t0:.0f}s", flush=True)
    if epi:
        elig = np.unique(np.concatenate([pres[d]['ids'] for d in test])); ne = len(elig)
        tot_test = obs['S0_total']; ov_test = sum(A.pair_overlap(pres[d]) for d in test)
        scale = (tot_test / ov_test) / cal['c']
        res['epi_scale'] = scale; res['n_elig'] = ne
        real_by_day = {d: days[d] for d in test}
        res['epi'] = {}
        nets = {}
        for mi, m in enumerate(models):
            rng = np.random.default_rng([seed, mi, 22])
            # persistence into the first test day is irrelevant for looping the test days: generate test days only
            nets[m] = [A.generate(m, rng, pres, N, groups, cal, test, scale=scale, extra=extra.get(m)) for _ in range(G_epi)]
            res['models'][m]['epi_total_ratio'] = float(np.mean([sum(len(nn[d]['s']) for d in test) for nn in nets[m]]) / tot_test)
        for si, (R0, tau) in enumerate(A.EPI_SCEN):
            beta = A.epi_beta(tot_test, ne, len(test), R0, tau)
            key = f"R0={R0},tau={tau}d"
            fs, off = A.sir(real_by_day, test, N, elig, beta, tau * A.DAY, nrep_real, seed * 1000 + si)
            e_real = A.epi_summary(fs, off, ne)
            row = dict(beta=beta, REAL=e_real)
            for mi, m in enumerate(models):
                F, O = [], []
                for g, nn in enumerate(nets[m]):
                    f, o = A.sir(nn, test, N, elig, beta, tau * A.DAY, nrep_model, seed * 100000 + si * 1000 + mi * 100 + g)
                    F.append(f); O.append(o)
                e = A.epi_summary(np.concatenate(F), np.concatenate(O), ne)
                e['pass'] = A.epi_pass(e, e_real)
                row[m] = e
            res['epi'][key] = row
            if verbose: print(f"   epi {key} done {time.time()-t0:.0f}s", flush=True)
        for m in models:
            npass = sum(res['epi'][k][m]['pass'] for k in res['epi'])
            res['models'][m]['score']['E_epidemic'] = bool(npass >= 3)
            res['models'][m]['epi_pass_count'] = int(npass)
    res['seconds'] = time.time() - t0
    return res


S_FAMS = list(A.FAMILIES.keys())


def claim(res_by_dataset, model):
    """Pre-registered claim ladder for one model over a list of dataset results."""
    n = len(res_by_dataset); epass = 0; strong = 0; a2 = 0
    for r in res_by_dataset:
        sc = r['models'][model]['score']
        sp = [sc[f] for f in S_FAMS if sc[f] is not None]
        e = bool(sc.get('E_epidemic'))
        epass += e
        strong += (e and sum(sp) >= 6)
        a2 += (e and not r['models']['WM']['score'].get('E_epidemic'))
    need = int(np.ceil(2 * n / 3))
    if strong >= need: return 'A1'
    if a2 >= need: return 'A2'
    return 'A0/A3'


def dump(obj, path):
    def f(o):
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.bool_,)): return bool(o)
        raise TypeError(type(o))
    json.dump(obj, open(path, 'w'), indent=1, default=f)
