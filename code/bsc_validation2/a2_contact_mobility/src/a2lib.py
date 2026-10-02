"""Avenue 2 library: real contact data, contact-structure statistics, generative contact models
(well-mixed, block, lattice random walk variants), calibration, SIR on temporal networks, scoring.
Time unit everywhere: one 20 s interval; a day has DAY = 4320 intervals."""
import ctypes, io, os, zipfile
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, 'data', 'raw')
DAY = 4320

# ------------------------------------------------------------------ data
SPEC = {
    'InVS13':     dict(file='workplace_InVS_tij.dat.zip', member='tij_InVS.dat', meta='workplace_InVS_metadata.txt', off=0),
    'InVS15':     dict(file='workplace_InVS15_tij.dat.gz', meta='workplace_InVS15_metadata.txt', off=0),
    'LH10':       dict(file='hospital_lyon_contacts.dat.gz', meta=None, off=46800, inline_groups=True),   # t=0 is Monday 13:00
    'LyonSchool': dict(file='primaryschool.csv.gz', meta='primaryschool_metadata.txt', off=0),
    'SFHH':       dict(file='SFHH_tij.dat.gz', meta=None, off=0),
    'Thiers13':   dict(file='HighSchool2013_proximity_net.csv.gz', meta='HighSchool2013_metadata.txt', off=3600),  # unix, CET
}
DEV = ['InVS13', 'LyonSchool', 'LH10']
CONF = ['InVS15', 'Thiers13', 'SFHH']


def load(name):
    """Return dict(name, N, groups (int array or None), group_names, days=[dict(s,i,j)], day_index)."""
    sp = SPEC[name]
    path = os.path.join(RAW, sp['file'])
    if 'member' in sp:
        df = pd.read_csv(io.BytesIO(zipfile.ZipFile(path).read(sp['member'])), sep=r'\s+', header=None)
    else:
        df = pd.read_csv(path, sep=r'\s+', header=None)
    t = df[0].values.astype(np.int64) + sp['off']
    a = df[1].values.astype(np.int64); b = df[2].values.astype(np.int64)
    ids = np.unique(np.r_[a, b]); N = len(ids)
    ia = np.searchsorted(ids, a); ib = np.searchsorted(ids, b)
    groups = None; gnames = None
    lab = {}
    if sp.get('inline_groups'):
        for x, g in zip(np.r_[a, b], np.r_[df[3].values, df[4].values]):
            lab[int(x)] = g
    elif sp['meta']:
        m = pd.read_csv(os.path.join(RAW, sp['meta']), sep=r'\s+', header=None)
        lab = {int(x): g for x, g in zip(m[0].values, m[1].values)}
    if lab:
        gl = [str(lab.get(int(x), 'NA')) for x in ids]
        gnames = sorted(set(gl)); groups = np.array([gnames.index(g) for g in gl])
    day = t // 86400; s = ((t % 86400) // 20).astype(np.int64)
    lo = np.minimum(ia, ib); hi = np.maximum(ia, ib)
    keep = lo != hi
    day, s, lo, hi = day[keep], s[keep], lo[keep], hi[keep]
    udays, cnt = np.unique(day, return_counts=True)
    udays = udays[cnt >= 0.01 * cnt.sum()]            # drop near-empty days (weekends)
    days = []
    for d in udays:
        mk = day == d
        k = np.unique(s[mk] * (N * N) + lo[mk] * N + hi[mk])
        days.append(dict(s=k // (N * N), i=(k // N) % N, j=k % N))
    return dict(name=name, N=N, groups=groups, group_names=gnames, days=days, day_index=[int(d) for d in udays])


def split(ndays):
    """First ceil(D/2) days = train, last floor(D/2) days = test."""
    ntest = ndays // 2
    return list(range(ndays - ntest)), list(range(ndays - ntest, ndays))


def presence(day):
    """Presence windows from the individual's own first/last contact of the day."""
    ind = np.r_[day['i'], day['j']]; ss = np.r_[day['s'], day['s']]
    ids = np.unique(ind)
    k = np.searchsorted(ids, ind)
    a = np.full(len(ids), 10**9); b = np.full(len(ids), -1)
    np.minimum.at(a, k, ss); np.maximum.at(b, k, ss)
    return dict(ids=ids, a=a, b=b)


def pair_overlap(pr, groups=None):
    """Total pair-intervals of joint presence; optionally (within-group, between-group)."""
    a, b = pr['a'], pr['b']
    ov = np.minimum(b[:, None], b[None, :]) - np.maximum(a[:, None], a[None, :]) + 1
    ov = np.triu(np.clip(ov, 0, None), 1)
    if groups is None:
        return float(ov.sum())
    g = groups[pr['ids']]
    same = g[:, None] == g[None, :]
    return float(ov[same].sum()), float(ov[~same].sum())


# ------------------------------------------------------------------ statistics
def events(day, N):
    """Maximal runs of consecutive intervals for a pair -> (i, j, start, dur)."""
    key = day['i'] * N + day['j']
    o = np.lexsort((day['s'], key)); k = key[o]; ss = day['s'][o]
    if len(k) == 0:
        z = np.zeros(0, dtype=np.int64); return z, z, z, z
    new = np.r_[True, (k[1:] != k[:-1]) | (ss[1:] != ss[:-1] + 1)]
    idx = np.flatnonzero(new); dur = np.diff(np.r_[idx, len(k)])
    return k[idx] // N, k[idx] % N, ss[idx], dur


def _cv(x):
    x = np.asarray(x, float)
    return float(x.std() / x.mean()) if len(x) and x.mean() > 0 else np.nan


def stats(days, pres, N, groups, test, all_days_for_persistence=True, raw=False):
    """Contact-structure statistics on the test days. days/pres are lists over all retained days."""
    durs, gaps, degs, nev, npd = [], [], [], 0, 0
    strength = np.zeros(N); present_any = np.zeros(N, bool)
    tot = 0; win = 0.0; within = 0.0
    bins_n, bins_c = [], []
    for d in test:
        dy = days[d]; pr = pres[d]
        tot += len(dy['s'])
        ei, ej, es, ed = events(dy, N)
        durs.append(ed); nev += len(ed); npd += len(np.unique(ei * N + ej))
        np.add.at(strength, dy['i'], 1); np.add.at(strength, dy['j'], 1)
        present_any[pr['ids']] = True
        # inter-contact gaps per individual
        ind = np.r_[ei, ej]; st = np.r_[es, es]; en = np.r_[es + ed, es + ed]
        o = np.lexsort((st, ind)); ind, st, en = ind[o], st[o], en[o]
        bnd = np.flatnonzero(np.r_[True, ind[1:] != ind[:-1], True])
        for u in range(len(bnd) - 1):
            lo, hi = bnd[u], bnd[u + 1]
            if hi - lo < 2: continue
            cm = np.maximum.accumulate(en[lo:hi])
            g = st[lo + 1:hi] - cm[:-1]
            gaps.append(g[g > 0])
        # daily degree for person-days with window >= 90 intervals
        deg = np.zeros(N); pk = np.unique(dy['i'] * N + dy['j'])
        np.add.at(deg, pk // N, 1); np.add.at(deg, pk % N, 1)
        long_ = pr['ids'][(pr['b'] - pr['a'] + 1) >= 90]
        degs.append(deg[long_])
        if groups is not None:
            within += float((groups[dy['i']] == groups[dy['j']]).sum())
        # density scaling (descriptive): 20-min bins
        if len(dy['s']):
            s0 = int(pr['a'].min()); s1 = int(pr['b'].max())
            edges = np.arange(s0, s1 + 61, 60)
            c, _ = np.histogram(dy['s'], bins=edges)
            occ = np.zeros(s1 - s0 + 2)
            np.add.at(occ, pr['a'] - s0, 1); np.add.at(occ, pr['b'] - s0 + 1, -1)
            occ = np.cumsum(occ)[:-1]
            nb = np.array([occ[k * 60:(k + 1) * 60].mean() if len(occ[k * 60:(k + 1) * 60]) else 0 for k in range(len(c))])
            bins_n.append(nb); bins_c.append(c)
    dur = np.concatenate(durs) if durs else np.zeros(0)
    gap = np.concatenate(gaps) if gaps else np.zeros(0)
    deg = np.concatenate(degs) if degs else np.zeros(0)
    # persistence between consecutive retained days (d-1, d) with d a test day
    num = den = 0
    for d in test:
        if d == 0: continue
        A = set((days[d - 1]['i'] * N + days[d - 1]['j']).tolist())
        both = np.zeros(N, bool); both[np.intersect1d(pres[d - 1]['ids'], pres[d]['ids'])] = True
        pk = np.unique(days[d]['i'] * N + days[d]['j'])
        pk = pk[both[pk // N] & both[pk % N]]
        den += len(pk); num += sum(1 for x in pk.tolist() if x in A)
    out = dict(
        S0_total=float(tot),
        S1_mean_dur=float(dur.mean()) if len(dur) else np.nan,
        S1_frac_ge3=float((dur >= 3).mean()) if len(dur) else np.nan,
        S1_timefrac_ge15=float(dur[dur >= 15].sum() / dur.sum()) if len(dur) else np.nan,
        S2_burst=float((gap.std() - gap.mean()) / (gap.std() + gap.mean())) if len(gap) > 1 else np.nan,
        S2_median_gap=float(np.median(gap)) if len(gap) else np.nan,
        S2_recurrence=float(1 - npd / nev) if nev else np.nan,
        S3_deg_mean=float(deg.mean()) if len(deg) else np.nan,
        S3_deg_cv=_cv(deg),
        S4_strength_cv=_cv(strength[present_any]),
        S5_persistence=float(num / den) if den else np.nan,
        S6_within_frac=float(within / tot) if (groups is not None and tot) else np.nan,
    )
    if bins_n:
        n = np.concatenate(bins_n); c = np.concatenate(bins_c)
        ok = (n >= 10) & (c > 0)
        if ok.sum() >= 5 and np.ptp(np.log(n[ok])) > 0:
            out['S7_alpha'] = float(np.polyfit(np.log(n[ok]), np.log(c[ok]), 1)[0])
            out['S7_n_range'] = float(np.percentile(n[ok], 90) / np.percentile(n[ok], 10))
        else:
            out['S7_alpha'] = np.nan; out['S7_n_range'] = np.nan
    if raw:
        out['_dur'] = dur; out['_gap'] = gap; out['_deg'] = deg; out['_strength'] = strength[present_any]
    return out


# ------------------------------------------------------------------ scoring (pre-registered tolerances)
RATIO = ['S0_total', 'S1_mean_dur', 'S2_median_gap', 'S3_deg_mean', 'S3_deg_cv', 'S4_strength_cv']
FRAC = ['S1_frac_ge3', 'S1_timefrac_ge15', 'S2_recurrence', 'S5_persistence', 'S6_within_frac']
FAMILIES = {
    'S0_level': ['S0_total'],
    'S1_duration': ['S1_mean_dur', 'S1_timefrac_ge15'],
    'S2_intercontact': ['S2_burst', 'S2_recurrence'],
    'S3_degree': ['S3_deg_mean', 'S3_deg_cv'],
    'S4_heterogeneity': ['S4_strength_cv'],
    'S5_persistence': ['S5_persistence'],
    'S6_groups': ['S6_within_frac'],
}


def metric_pass(key, model, obs):
    if obs is None or model is None or np.isnan(obs) or np.isnan(model):
        return None
    if key == 'S2_burst':
        return bool(abs(model - obs) <= 0.10)
    r = model / obs if obs != 0 else np.inf
    ok = 0.8 <= r <= 1.25
    if key in FRAC:
        ok = ok or abs(model - obs) <= 0.03
    return bool(ok)


def score(model_stats, obs_stats):
    """Family passes: a family passes iff all its applicable metrics pass; None if not applicable."""
    res = {}
    for fam, keys in FAMILIES.items():
        ps = [metric_pass(k, model_stats.get(k, np.nan), obs_stats.get(k, np.nan)) for k in keys]
        ps = [p for p in ps if p is not None]
        res[fam] = (all(ps) if ps else None)
    return res


def epi_pass(m, o):
    """m, o: dicts with R_index, P_major, attack. Pass iff all three within tolerance."""
    r = m['R_index'] / o['R_index'] if o['R_index'] > 0 else np.inf
    return bool(0.9 <= r <= 1/0.9 and abs(m['P_major'] - o['P_major']) <= 0.05 and abs(m['attack'] - o['attack']) <= 0.05)


# ------------------------------------------------------------------ models
def _contacts_from_positions(P, ids, T0):
    """P[T, n] site index (negative = absent, unique). Returns s, i, j with global ids, i<j."""
    o = np.argsort(P, axis=1, kind='stable'); Ps = np.take_along_axis(P, o, axis=1)
    S, I, J = [], [], []
    k = 1
    while k < P.shape[1]:
        eq = (Ps[:, k:] == Ps[:, :-k]) & (Ps[:, k:] >= 0)
        if not eq.any(): break
        r, c = np.nonzero(eq)
        a = ids[o[r, c]]; b = ids[o[r, c + k]]
        S.append(r + T0); I.append(np.minimum(a, b)); J.append(np.maximum(a, b))
        k += 1
    if not S:
        z = np.zeros(0, dtype=np.int64); return dict(s=z, i=z.copy(), j=z.copy())
    return dict(s=np.concatenate(S).astype(np.int64), i=np.concatenate(I).astype(np.int64), j=np.concatenate(J).astype(np.int64))


def rw_day(rng, pr, Lx, Ly, p_fast, p_slow=None, f_slow=0.0, home=None, bias=0.0, start='uniform'):
    """Lattice random walk day. Each 20 s step a present walker attempts a move to one of the four
    von Neumann neighbours (uniformly, or towards its home site with probability `bias`); the move is
    accepted with probability H(p(r), p(r')) (harmonic mean, the interface rule of the model); attempts
    across a wall are rejected (reflecting boundary). Contact = same site at a sampling instant."""
    ids, a, b = pr['ids'], pr['a'], pr['b']; n = len(ids)
    T0 = int(a.min()); T1 = int(b.max()); T = T1 - T0 + 1
    xs = int(round(f_slow * Lx))                       # columns x < xs are the slow zone
    if home is not None and start == 'home':
        x = home[0][ids].copy(); y = home[1][ids].copy()
    else:
        x = rng.integers(0, Lx, n); y = rng.integers(0, Ly, n)
    if home is not None:
        hx = home[0][ids]; hy = home[1][ids]
    P = np.empty((T, n), dtype=np.int64)
    neg = -(np.arange(n) + 1)
    dxs = np.array([1, -1, 0, 0]); dys = np.array([0, 0, 1, -1])
    ps = p_fast if p_slow is None else p_slow
    for s in range(T0, T1 + 1):
        act = (a <= s) & (s <= b)
        mv = act & (a < s)
        d = rng.integers(0, 4, n)
        dx = dxs[d]; dy = dys[d]
        if home is not None and bias > 0:
            hb = rng.random(n) < bias
            ex = hx - x; ey = hy - y
            usex = (ex != 0) & ((ey == 0) | (rng.random(n) < 0.5))
            bx = np.where(usex, np.sign(ex), 0); by = np.where(~usex, np.sign(ey), 0)
            dx = np.where(hb, bx, dx); dy = np.where(hb, by, dy)
        nx = x + dx; ny = y + dy
        ok = (nx >= 0) & (nx < Lx) & (ny >= 0) & (ny < Ly) & ((dx != 0) | (dy != 0))
        if xs > 0:
            p0 = np.where(x < xs, ps, p_fast); p1 = np.where(nx < xs, ps, p_fast)
            acc = 2 * p0 * p1 / (p0 + p1)
        else:
            acc = p_fast
        go = mv & ok & (rng.random(n) < acc)
        x = np.where(go, nx, x); y = np.where(go, ny, y)
        P[s - T0] = np.where(act, x * Ly + y, neg)
    return _contacts_from_positions(P, ids, T0)


def pair_markov_day(rng, pr, groups, c_in, c_out, mean_dur):
    """Memoryless pair model (WM if c_in == c_out, BLOCK otherwise). For each jointly present pair,
    contact events start as a Poisson process such that the expected contact-time fraction is c;
    event durations are geometric with the given mean; overlapping events merge."""
    ids, a, b = pr['ids'], pr['a'], pr['b']; n = len(ids)
    iu, ju = np.triu_indices(n, 1)
    lo = np.maximum(a[iu], a[ju]); hi = np.minimum(b[iu], b[ju]); L = hi - lo + 1
    ok = L > 0; iu, ju, lo, L = iu[ok], ju[ok], lo[ok], L[ok]
    if groups is None:
        c = np.full(len(iu), c_in)
    else:
        g = groups[ids]; c = np.where(g[iu] == g[ju], c_in, c_out)
    ne = rng.poisson(c * L / mean_dur)
    rep = np.repeat(np.arange(len(iu)), ne)
    st = lo[rep] + np.floor(rng.random(len(rep)) * L[rep]).astype(np.int64)
    du = rng.geometric(1.0 / mean_dur, len(rep))
    en = np.minimum(st + du, lo[rep] + L[rep])
    tot = int((en - st).sum())
    S = np.repeat(st, en - st) + (np.arange(tot) - np.repeat(np.cumsum(en - st) - (en - st), en - st))
    I = np.repeat(ids[iu[rep]], en - st); J = np.repeat(ids[ju[rep]], en - st)
    N = int(ids.max()) + 1
    k = np.unique(S * (N * N) + np.minimum(I, J) * N + np.maximum(I, J))
    return dict(s=k // (N * N), i=(k // N) % N, j=k % N)


def p_from_q(q):
    """Persistence probability of a co-location, q = (1-p)^2 + p^2/4, inverted for p in (0, 0.8]."""
    q = float(np.clip(q, 0.2, 0.999999))
    return float((2 - np.sqrt(4 - 5 * (1 - q))) / 2.5)


def lattice_dims(M):
    M = max(int(round(M)), 4)
    Lx = max(int(round(np.sqrt(M))), 2); Ly = max(int(round(M / Lx)), 2)
    return Lx, Ly


def fit_two_geometric(dur):
    """MLE of a two-component geometric mixture by EM. Returns (w_slow, q_slow, q_fast), q_slow >= q_fast."""
    d = np.asarray(dur, float); vals, cnt = np.unique(d, return_counts=True)
    m = d.mean()
    q1, q2, w = 1 - 1 / (3 * m), max(1 - 1 / max(0.6 * m, 1.01), 0.01), 0.2
    for _ in range(2000):
        l1 = w * q1 ** (vals - 1) * (1 - q1); l2 = (1 - w) * q2 ** (vals - 1) * (1 - q2)
        r = l1 / (l1 + l2 + 1e-300)
        n1 = (r * cnt).sum(); n2 = ((1 - r) * cnt).sum()
        m1 = (r * cnt * vals).sum() / n1; m2 = ((1 - r) * cnt * vals).sum() / n2
        nq1, nq2, nw = 1 - 1 / max(m1, 1.0001), 1 - 1 / max(m2, 1.0001), n1 / (n1 + n2)
        if abs(nq1 - q1) + abs(nq2 - q2) + abs(nw - w) < 1e-10: break
        q1, q2, w = nq1, nq2, nw
    if q1 < q2: q1, q2, w = q2, q1, 1 - w
    return float(w), float(q1), float(q2)


def calibrate(days, pres, N, groups, train):
    """All base-model parameters from the training days only."""
    tot = sum(len(days[d]['s']) for d in train)
    ov = sum(pair_overlap(pres[d]) for d in train)
    dur = np.concatenate([events(days[d], N)[3] for d in train])
    cal = dict(c=tot / ov, mean_dur=float(dur.mean()), n_events=int(len(dur)))
    cal['M'] = 1.0 / cal['c']
    cal['p'] = p_from_q(1 - 1 / cal['mean_dur'])
    w, q1, q2 = fit_two_geometric(dur)
    cal['het'] = dict(w_slow=w, q_slow=q1, q_fast=q2, p_slow=p_from_q(q1), p_fast=p_from_q(q2))
    fs = (w / (1 - q1)) / (w / (1 - q1) + (1 - w) / (1 - q2))
    cal['het']['f_slow'] = float(fs)
    if groups is not None:
        win = sum(float((groups[days[d]['i']] == groups[days[d]['j']]).sum()) for d in train)
        oin = oout = 0.0
        for d in train:
            x, y = pair_overlap(pres[d], groups); oin += x; oout += y
        cal['c_in'] = win / oin if oin > 0 else cal['c']
        cal['c_out'] = (tot - win) / oout if oout > 0 else cal['c']
    else:
        cal['c_in'] = cal['c_out'] = cal['c']
    return cal


def generate(model, rng, pres, N, groups, cal, day_list, scale=1.0, extra=None):
    """Generate contacts for the listed days. `scale` multiplies the contact probability (level matching:
    c -> c*scale, i.e. M -> M/scale)."""
    out = {}
    if model in ('RW0', 'RWhet', 'RWhome'):
        Lx, Ly = lattice_dims(cal['M'] / scale if model != 'RWhome' else extra['M'] / scale)
    home = None
    if model == 'RWhome':
        home = make_homes(np.random.default_rng(extra.get('home_seed', 0) + int(rng.integers(1 << 30))), N, groups, Lx, Ly)
    for d in day_list:
        pr = pres[d]
        if model == 'WM':
            out[d] = pair_markov_day(rng, pr, None, cal['c'] * scale, cal['c'] * scale, cal['mean_dur'])
        elif model == 'BLOCK':
            out[d] = pair_markov_day(rng, pr, groups, cal['c_in'] * scale, cal['c_out'] * scale, cal['mean_dur'])
        elif model == 'RW0':
            out[d] = rw_day(rng, pr, Lx, Ly, cal['p'])
        elif model == 'RWhet':
            h = cal['het']
            out[d] = rw_day(rng, pr, Lx, Ly, h['p_fast'], h['p_slow'], h['f_slow'])
        elif model == 'RWhome':
            out[d] = rw_day(rng, pr, Lx, Ly, extra['p_fast'], extra.get('p_slow'), extra.get('f_slow', 0.0),
                            home=home, bias=extra['bias'], start='home')
        else:
            raise ValueError(model)
    return out


def make_homes(rng, N, groups, Lx, Ly):
    """Home sites: groups occupy contiguous vertical strips with width proportional to group size;
    each individual's home is a uniformly random site of its group's strip. No groups -> uniform."""
    hx = np.zeros(N, dtype=np.int64); hy = rng.integers(0, Ly, N)
    if groups is None:
        hx = rng.integers(0, Lx, N)
    else:
        G = int(groups.max()) + 1
        sizes = np.bincount(groups, minlength=G).astype(float)
        edges = np.round(np.cumsum(np.r_[0, sizes]) / sizes.sum() * Lx).astype(int)
        for g in range(G):
            lo, hi = edges[g], max(edges[g + 1], edges[g] + 1)
            m = groups == g
            hx[m] = rng.integers(lo, min(hi, Lx), m.sum()) if lo < Lx else Lx - 1
    return hx, hy


# ------------------------------------------------------------------ SIR
_lib = None


def _sir_lib():
    global _lib
    if _lib is None:
        so = os.path.join(HERE, 'libsir.dylib')
        if not os.path.exists(so):
            os.system(f'cc -O2 -shared -fPIC -o {so} {os.path.join(HERE, "sir.c")}')
        _lib = ctypes.CDLL(so)
        ip = ctypes.POINTER(ctypes.c_int)
        _lib.sir_loop.argtypes = [ctypes.c_int, ip, ip, ip, ctypes.c_int, ctypes.c_int, ip, ctypes.c_long,
                                  ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_uint64, ctypes.c_int, ip, ip]
    return _lib


def sir(contacts_by_day, day_list, N, elig, beta, tau, nrep, seed, maxloops=400):
    """SIR on the listed days looped periodically (period = len(day_list) days)."""
    S = np.concatenate([contacts_by_day[d]['s'] + k * DAY for k, d in enumerate(day_list)])
    I = np.concatenate([contacts_by_day[d]['i'] for d in day_list]); J = np.concatenate([contacts_by_day[d]['j'] for d in day_list])
    o = np.argsort(S, kind='stable')
    S = np.ascontiguousarray(S[o], dtype=np.int32); I = np.ascontiguousarray(I[o], dtype=np.int32); J = np.ascontiguousarray(J[o], dtype=np.int32)
    el = np.ascontiguousarray(elig, dtype=np.int32)
    fs = np.zeros(nrep, dtype=np.int32); off = np.zeros(nrep, dtype=np.int32)
    ip = ctypes.POINTER(ctypes.c_int)
    _sir_lib().sir_loop(len(S), S.ctypes.data_as(ip), I.ctypes.data_as(ip), J.ctypes.data_as(ip), int(N), len(el),
                        el.ctypes.data_as(ip), int(len(day_list) * DAY), float(beta), float(tau), int(nrep),
                        ctypes.c_uint64(seed), int(maxloops), fs.ctypes.data_as(ip), off.ctypes.data_as(ip))
    return fs, off


def epi_summary(fs, off, n_elig):
    ar = fs / n_elig
    major = ar >= 0.10
    return dict(R_index=float(off.mean()), P_major=float(major.mean()), attack=float(ar.mean()),
                attack_major=float(ar[major].mean()) if major.any() else float('nan'), nrep=int(len(fs)))


EPI_SCEN = [(1.5, 1.0), (3.0, 1.0), (1.5, 4.0), (3.0, 4.0)]   # (nominal R0, mean infectious period in days)


def epi_beta(total_test, n_elig, n_days, R0, tau_days):
    W = 2.0 * total_test / (n_elig * n_days * DAY)       # contact-intervals per person per interval of clock time
    return R0 / (W * tau_days * DAY)
