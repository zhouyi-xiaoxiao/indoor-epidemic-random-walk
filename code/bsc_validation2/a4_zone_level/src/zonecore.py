"""Zone-level (class / grade / school) chain-binomial model core for test T1 (Matsumoto influenza 2014/15).
hazard for a susceptible student of class c on day t:
    lam_c(t) = b . x_c(t) + a0 + a1 * Pcity_{-s}(t)
x: prevalence features built from the onset-to-onset kernel k (gamma mean 1.7 sd 1, daily bins, tau = 1..10).
Everything is at class-day level (all susceptibles of a class share the hazard)."""
import numpy as np, pandas as pd
from scipy import stats, optimize
T = 213                       # days; day 1 = 1 Oct 2014, day 213 = 1 May 2015
BREAK = np.arange(89, 101)    # onset days 28 Dec .. 8 Jan: no within-school transmission (Endo et al. 'duringbreak')
NEVER = 32767
LYON = np.array([0.7619836508613838, 0.1044788607478434, 0.13353748839077287])   # class, grade, school shares

PERPAIR = np.array([352.30174081237914, 46.42007434944238, 7.416012267087961])   # Lyon seconds per pair per day

def kernel(mean=1.7, sd=1.0, kmax=10):
    g = stats.gamma(a=mean**2/sd**2, scale=sd**2/mean)
    k = g.cdf(np.arange(1, kmax+1)) - g.cdf(np.arange(0, kmax)); return k/k.sum()
K = kernel()

def load_structure(path='data/derived/matsumoto_students.csv'):
    """structure only: school, grade, class ids of the respondents (NO outcomes)"""
    df = pd.read_csv(path, usecols=['schoolID', 'gradeID', 'classID'])
    return build_structure(df)

def build_structure(df):
    g = df.groupby(['schoolID', 'gradeID', 'classID']).size().reset_index(name='n')
    st = dict(school=g.schoolID.values, grade=g.gradeID.values, cls=g.classID.values, n=g.n.values.astype(float))
    key = g.schoolID.values*10+g.gradeID.values
    st['gkey'] = key
    st['n_grade'] = pd.Series(st['n']).groupby(key).transform('sum').values
    st['n_school'] = pd.Series(st['n']).groupby(st['school']).transform('sum').values
    st['N'] = st['n'].sum()
    return st

def load_outcomes(path='data/derived/matsumoto_students.csv'):
    df = pd.read_csv(path)
    st = build_structure(df)
    idx = {(s, g, c): i for i, (s, g, c) in enumerate(zip(st['school'], st['grade'], st['cls']))}
    ci = np.array([idx[(s, g, c)] for s, g, c in zip(df.schoolID, df.gradeID, df.classID)])
    on = df.onset.values.copy(); on[df.isinfected.values == 0] = NEVER
    Y = np.zeros((len(st['n']), T+1))
    ok = on <= T
    np.add.at(Y, (ci[ok], on[ok]), 1)
    return st, Y, dict(n_late=int(((on > T) & (on < NEVER)).sum()), n_inf=int((on < NEVER).sum()))

def group_sum(A, key):
    """sum rows of A within groups given by key, broadcast back"""
    u, inv = np.unique(key, return_inverse=True)
    G = np.zeros((len(u),) + A.shape[1:]); np.add.at(G, inv, A); return G[inv]

def features(st, Y):
    """Y[c,t] onsets (t = 0..T, column 0 unused). Returns dict of class-day arrays."""
    nC = Y.shape[0]
    Kc = np.zeros((nC, T+1))                       # infectious pressure sum_j k(t - t_j) by class
    for tau, k in enumerate(K, start=1):
        Kc[:, tau:] += k*Y[:, :T+1-tau]
    Kg = group_sum(Kc, st['gkey']); Ks = group_sum(Kc, st['school'])
    n, ng, ns = st['n'][:, None], st['n_grade'][:, None], st['n_school'][:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        Pc = np.where(n > 1, Kc/(n-1), 0.0)
        Pg = np.where(ng > n, (Kg-Kc)/(ng-n), 0.0)
        Ps = np.where(ns > ng, (Ks-Kg)/(ns-ng), 0.0)
        PH = np.where(ns > 1, Ks/(ns-1), 0.0)
    Ks_raw = Ks.copy()
    for P in (Pc, Pg, Ps, PH): P[:, BREAK] = 0.0
    # city level (4th mixing level, causal): kernel-weighted prevalence among respondents of all OTHER schools;
    # NOT switched off during the school break (community transmission continues)
    Kall = Kc.sum(0, keepdims=True)
    C = (Kall - Ks_raw)/(st['N'] - ns)
    cumY = np.cumsum(Y, 1)
    S = n - np.concatenate([np.zeros((nC, 1)), cumY[:, :-1]], 1)        # susceptible at start of day t
    Kraw = dict(Kc=Kc, Ks=Ks_raw, Kg=Kg)
    return dict(Pc=Pc, Pg=Pg, Ps=Ps, PH=PH, C=C, S=S, Y=Y, Kraw=Kraw)

def design(F, model, st=None, rho=None, w=LYON):
    one = np.ones_like(F['C'])
    if model == 'H': X = [F['PH']]
    elif model == 'X': X = [F['Pc']]
    elif model == 'Z': X = [w[0]*F['Pc']+w[1]*F['Pg']+w[2]*F['Ps']]
    elif model == 'F': X = [F['Pc'], F['Pg'], F['Ps']]
    elif model == 'Zdd':    # sensitivity: density-dependent transfer of Lyon per-pair contact times (class, grade, school)
        Kc, Kg, Ks = F['Kraw']['Kc'], F['Kraw']['Kg'], F['Kraw']['Ks']
        P = (PERPAIR[0]*Kc + PERPAIR[1]*(Kg-Kc) + PERPAIR[2]*(Ks-Kg))/PERPAIR[0]/22.2; P[:, BREAK] = 0.0
        X = [P]
    elif model == 'C':
        n, ns = st['n'][:, None], st['n_school'][:, None]
        Kc, Ks = F['Kraw']['Kc'], F['Kraw']['Ks']
        P = (rho*Kc + (Ks-Kc))/np.maximum(rho*(n-1)+(ns-n), 1e-12); P[:, BREAK] = 0.0
        X = [P]
    return np.stack(X+[one, F['C']], -1)[:, 1:, :]      # drop day 0

def _nll(theta, X, S, Y):
    lam = X@theta
    lam = np.maximum(lam, 1e-300)
    em = np.exp(-lam); p = -np.expm1(-lam)
    ll = -(lam*(S-Y)).sum() + (Y*np.log(p)).sum()
    r = -(S-Y) + Y*em/p
    return -ll, -(X*r[..., None]).sum((0, 1))

def loglik(theta, X, S, Y, per_class=False):
    lam = np.maximum(X@theta, 1e-300)
    ll = -(lam*(S-Y)) + Y*np.log(-np.expm1(-lam))
    return ll.sum(1) if per_class else ll.sum()

def fit(X, S, Y, starts=None):
    p = X.shape[-1]; best = None
    if starts is None:
        starts = [np.r_[np.full(p-2, 0.5/(p-2)), 1e-4, 0.5], np.r_[np.full(p-2, 0.1), 1e-3, 1.0], np.r_[np.full(p-2, 1.0/(p-2)), 1e-5, 0.1]]
    for x0 in starts:
        r = optimize.minimize(_nll, x0, args=(X, S, Y), jac=True, method='L-BFGS-B', bounds=[(0, None)]*(p-2)+[(1e-9, None), (0, None)],
                              options=dict(maxiter=2000, ftol=1e-13, gtol=1e-8))
        if best is None or r.fun < best.fun: best = r
    return best.x, -best.fun

def fit_model(F, st, model, mask=None, w=LYON):
    """mask: boolean over classes (calibration set). Returns dict(theta, ll, rho)"""
    S, Y = F['S'][:, 1:], F['Y'][:, 1:]
    m = slice(None) if mask is None else mask
    if model != 'C':
        X = design(F, model, st, w=w)
        th, ll = fit(X[m], S[m], Y[m]); return dict(theta=th, ll=ll, rho=None)
    def prof(lr):
        X = design(F, 'C', st, rho=np.exp(lr)); return -fit(X[m], S[m], Y[m])[1]
    grid = np.linspace(0, np.log(3000), 9); vals = [prof(g) for g in grid]; j = int(np.argmin(vals))
    lo, hi = grid[max(j-1, 0)], grid[min(j+1, len(grid)-1)]
    r = optimize.minimize_scalar(prof, bounds=(lo, hi), method='bounded', options=dict(xatol=0.03))
    lr = r.x if r.fun < vals[j] else grid[j]
    X = design(F, 'C', st, rho=np.exp(lr)); th, ll = fit(X[m], S[m], Y[m])
    return dict(theta=th, ll=ll, rho=float(np.exp(lr)))

def heldout_ll(F, st, model, res, mask, w=LYON):
    """per-class log predictive likelihood on classes in mask, parameters from res"""
    X = design(F, model, st, rho=res['rho'], w=w)
    return loglik(res['theta'], X[mask], F['S'][mask, 1:], F['Y'][mask, 1:], per_class=True)

def simulate(st, model, theta, rng, w=LYON, rho=None):
    """forward chain-binomial simulation of the whole city. theta = (b..., a0, a1) as in design();
    external hazard = a0 + a1 * (kernel-weighted prevalence in all other schools), endogenous."""
    nC = len(st['n']); n = st['n']; ng = st['n_grade']; ns = st['n_school']; N = st['N']
    Y = np.zeros((nC, T+1)); Kc = np.zeros((nC, T+1+len(K))); S = n.copy()
    gk_u, gk_inv = np.unique(st['gkey'], return_inverse=True); sk_u, sk_inv = np.unique(st['school'], return_inverse=True)
    isbreak = np.zeros(T+1, bool); isbreak[BREAK] = True
    b = np.atleast_1d(theta[:-2]); a0, a1 = theta[-2], theta[-1]
    for t in range(1, T+1):
        kc = Kc[:, t]; ks = np.bincount(sk_inv, kc)[sk_inv]
        lam = a0 + a1*(kc.sum()-ks)/(N-ns)
        if not isbreak[t]:
            kg = np.bincount(gk_inv, kc)[gk_inv]
            with np.errstate(divide='ignore', invalid='ignore'):
                if model == 'H': P = np.where(ns > 1, ks/(ns-1), 0)*b[0]
                elif model == 'X': P = np.where(n > 1, kc/(n-1), 0)*b[0]
                elif model == 'C': P = b[0]*(rho*kc+(ks-kc))/np.maximum(rho*(n-1)+(ns-n), 1e-12)
                else:
                    ww = b[0]*np.asarray(w) if model == 'Z' else b
                    P = ww[0]*np.where(n > 1, kc/(n-1), 0)+ww[1]*np.where(ng > n, (kg-kc)/(ng-n), 0)+ww[2]*np.where(ns > ng, (ks-kg)/(ns-ng), 0)
            lam = lam + P
        y = rng.binomial(S.astype(int), -np.expm1(-lam)); Y[:, t] = y; S = S - y
        Kc[:, t+1:t+1+len(K)] += y[:, None]*K[None, :]
    return Y

def dispersion_stat(st, Y, mask=None):
    """Pearson heterogeneity of class attack rates around the school attack rate, pooled; ~1 under homogeneity"""
    I = Y.sum(1); n = st['n']; sch = st['school']
    m = np.ones(len(n), bool) if mask is None else mask
    num = 0.0; df = 0
    for s in np.unique(sch[m]):
        j = m & (sch == s)
        if j.sum() < 2: continue
        p = I[j].sum()/n[j].sum()
        if p <= 0 or p >= 1: continue
        num += ((I[j]-n[j]*p)**2/(n[j]*p*(1-p))).sum(); df += j.sum()-1
    return num/df

def tv(a, b): return 0.5*np.abs(np.asarray(a)-np.asarray(b)).sum()
