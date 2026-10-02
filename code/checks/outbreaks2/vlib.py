"""Re-implementation with separately written code for the adversarial verification of the second outbreak test.

Nothing is imported from the first analysis's code except the transcribed seat-map strings
(src/a1/digitized.py), which were re-read against the published figures by the re-check.

Kernel: M2 exposure X = int_0^T (T-s) e^{-kappa s} G(r,s|r0) ds computed by UNIFORMIZATION
(a series with non-negative terms, component-wise accurate), different from both the first analysis's
cosine-mode sum and its image/Bessel sum.  For very large D (lambda*T > CUT) a numerical
eigen-decomposition is used (exposures are then nearly flat and well resolved).
"""
import os, sys, itertools
import numpy as np
from scipy.special import gammainc
from scipy.optimize import brentq
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bsc_validation2", "a1_more_outbreaks", "src"))
from a1 import digitized as dg          # strings only

LOGD = np.round(np.linspace(-2.0, 4.0, 61), 6)
LOGRHO = np.round(np.linspace(0.0, 2.5, 51), 6)
DEP = 0.93
X33 = [0, 1, 2, 4, 5, 6]
X242 = [0, 1, 3, 4, 5, 6, 8, 9]
X343 = [0, 1, 2, 4, 5, 6, 7, 9, 10, 11]
CUT = 2.0e5
SEEDV = 777001


# ----------------------------------------------------------------------------- kernels
def _lap1d(n):
    L = np.zeros((n, n))
    for i in range(n - 1):
        L[i, i + 1] = L[i + 1, i] = 1.0
        L[i, i] -= 1.0
        L[i + 1, i + 1] -= 1.0
    return L


def _w(L, T):
    x = L * T
    out = np.where(x < 1e-4, T * T * (0.5 - x / 6 + x * x / 24), (x - 1 + np.exp(-x)) / np.where(L > 0, L, 1) ** 2)
    return out


def kern_eigh(nx, ny, ax, ay, srcsets, T, kappas, D):
    """X[kappa, S, site] by numerical eigen-decomposition (site = y*nx + x)."""
    mx, Vx = np.linalg.eigh(-_lap1d(nx))
    my, Vy = np.linalg.eigh(-_lap1d(ny))
    mx = np.clip(mx, 0, None); my = np.clip(my, 0, None)
    out = np.zeros((len(kappas), len(srcsets), nx * ny))
    for ik, kap in enumerate(kappas):
        W = _w(kap + D * (my[:, None] / ay ** 2 + mx[None, :] / ax ** 2), T)      # (iy, ix)
        for s, srcs in enumerate(srcsets):
            acc = np.zeros((ny, nx))
            for (x0, y0) in srcs:
                acc += (Vy * Vy[y0][None, :]) @ W @ (Vx * Vx[x0][None, :]).T
            out[ik, s] = acc.reshape(-1)
    return out


def kern_table(nx, ny, ax, ay, srcsets, T, kappas, Dgrid, verbose=False):
    """X[iD, ik, S, site] for one segment of duration T.  srcsets: list of lists of (x0, y0)."""
    kappas = np.asarray(kappas, float)
    lam0 = 2.0 / ax ** 2 + 2.0 / ay ** 2
    px, py = (1 / ax ** 2) / lam0, (1 / ay ** 2) / lam0
    S = len(srcsets)
    out = np.zeros((len(Dgrid), len(kappas), S, nx * ny))
    uni = [i for i, D in enumerate(Dgrid) if D * lam0 * T <= CUT]
    for i, D in enumerate(Dgrid):
        if i not in uni:
            out[i] = kern_eigh(nx, ny, ax, ay, srcsets, T, kappas, D)
    if not uni:
        return out
    lam = np.array([Dgrid[i] * lam0 for i in uni])                      # (nD,)
    kmax = (lam * T + 14 * np.sqrt(lam * T) + 80).astype(int)
    mu = kappas[None, :] + lam[:, None]                                 # (nD, nk)
    u = np.zeros((S, ny, nx))
    for s, srcs in enumerate(srcsets):
        for (x0, y0) in srcs:
            u[s, y0, x0] += 1.0
    nbx = np.full(nx, 2.0); nbx[0] = nbx[-1] = 1.0
    nby = np.full(ny, 2.0); nby[0] = nby[-1] = 1.0
    stay = np.clip(1.0 - px * nbx[None, :] - py * nby[:, None], 0.0, None)
    acc = np.zeros((len(uni), len(kappas), S, ny * nx))
    CH = 400
    k0 = 0
    K = int(kmax.max())
    while k0 <= K:
        ks = np.arange(k0, min(k0 + CH, K + 1))
        U = np.zeros((len(ks), S, ny, nx))
        for j in range(len(ks)):
            U[j] = u
            v = stay[None] * u
            v[:, :, 1:] += px * u[:, :, :-1]
            v[:, :, :-1] += px * u[:, :, 1:]
            v[:, 1:, :] += py * u[:, :-1, :]
            v[:, :-1, :] += py * u[:, 1:, :]
            u = v
        act = np.flatnonzero(kmax >= k0)
        if len(act):
            m = mu[act][:, :, None]                                     # (a, nk, 1)
            l = lam[act][:, None, None]
            kk = ks[None, None, :]
            c = np.exp(kk * np.log(l / m)) / m * (T * gammainc(kk + 1, m * T) - (kk + 1) / m * gammainc(kk + 2, m * T))
            c = np.clip(c, 0.0, None)                                   # (a, nk, nks)
            acc[act] += np.einsum("ack,ksn->acsn", c, U.reshape(len(ks), S, ny * nx), optimize=True)
        k0 += CH
    for j, i in enumerate(uni):
        out[i] = acc[j]
    return out


# ----------------------------------------------------------------------------- statistics
def compositions(sizes, K):
    rng = [range(0, min(s, K) + 1) for s in sizes[:-1]]
    out = [h + (K - sum(h),) for h in itertools.product(*rng) if 0 <= K - sum(h) <= sizes[-1]]
    return np.array(out, int)


def pb(p, K):
    q = np.zeros(K + 1); q[0] = 1.0
    for pi in p:
        q[1:] = q[1:] * (1 - pi) + q[:-1] * pi
        q[0] *= (1 - pi)
    return q


def eps_profile(X, K):
    f = lambda le: np.sum(-np.expm1(-np.exp(le) * X)) - K
    lo = np.log(K / X.sum()) - 1.0
    hi = lo + 2.0
    while f(hi) < 0:
        hi += 2.0
        if hi > 600:
            return np.inf
    return np.exp(brentq(f, lo, hi, xtol=1e-12, rtol=1e-12))


def cond_pmf_p(p, bins, sizes, K, comps):
    pr = np.ones(len(comps))
    for b in range(len(sizes)):
        pr *= pb(p[bins == b], K)[comps[:, b]]
    s = pr.sum()
    return pr / s if s > 0 else np.full(len(comps), 1.0 / len(comps))


def cond_pmf(X, bins, sizes, K, comps):
    X = np.clip(np.asarray(X, float), 1e-300, None)
    eps = eps_profile(X, K)
    p = np.ones_like(X) if not np.isfinite(eps) else -np.expm1(-eps * X)
    return cond_pmf_p(np.clip(p, 0, 1 - 1e-15), bins, sizes, K, comps)


def cond_pmf_trueCRR(near, rho, bins, sizes, K, comps):
    """p_near = rho*q, p_far = q with expected total K (exact constant risk ratio)."""
    n1, n2 = near.sum(), (~near).sum()
    q = K / (rho * n1 + n2)
    q = min(q, (1 - 1e-12) / rho)
    return cond_pmf_p(np.where(near, rho * q, q), bins, sizes, K, comps)


def two_sided(pmf, i):
    return float(pmf[pmf <= pmf[i] * (1 + 1e-9) + 1e-15].sum())


def kappa_nodes(kind, lo=10.0, hi=30.0, n=7):
    q = (np.arange(n) + 0.5) / n
    if kind == "air":
        return np.exp(np.log(lo) + q * (np.log(hi) - np.log(lo))) + DEP
    if kind == "ward":
        return 7.79 * np.exp(0.2 * norm.ppf(q)) + DEP
    if kind == "plant":
        return 0.3 + 0.7 * q + DEP
    raise ValueError


# ----------------------------------------------------------------------------- geometry
def seatmap_cells(ev):
    """list of (x, y, symbol) in normalised symbols I / o / C / P(non-case in primary) / other."""
    if ev == "F1":
        tr = {"X": "o", "/": "o", "G": "C", "I": "I", ".": "e", "-": "-"}
        cells = [(X33[i], j, tr[ch]) for i, k in enumerate("ABCDEF") for j, ch in enumerate(dg.OLSEN[k])]
        return cells, 7, 22, 3.0
    if ev == "F3":
        cells = [(X343[i], j, ch) for i, s in enumerate(dg.BAKER) for j, ch in enumerate(s)]
        return cells, 12, 15, 13.0
    if ev == "F4":
        cells = [(X242[i], r - 1, ch) for r, s in dg.YOUNG.items() for i, ch in enumerate(s)]
        return cells, 10, 47, 9.5
    if ev == "F5":
        cells = [(X33[i], dg.TOYO_ROWS.index(r), ch) for r, s in dg.TOYO.items() for i, ch in enumerate(s)]
        return cells, 7, 30, 2.0
    if ev == "F6":
        cells = [(X242[i], j, ("C" if ch == "P" else ch)) for i, s in enumerate(dg.SPEAKE) for j, ch in enumerate(s)]
        return cells, 10, 17, 5.0
    if ev == "F8":
        cells = [(X33[i], dg.HOEHL_ROWS.index(r), ch) for r, s in dg.HOEHL.items() for i, ch in enumerate(s)]
        return cells, 7, 31, 4.67
    raise ValueError(ev)


def build(ev, n_cfg=None, seed=SEEDV, edges=(1, 2, 5), hh=False):
    """Return dict(nx, ny, ax, ay, T, nseg, srcsets, kappas, cfgs=[(iS, ik, rec_sites)], bins, sizes, obs, near, xy)"""
    rng = np.random.default_rng([seed, hash(ev) % 1000 if False else sum(map(ord, ev))])
    g = dict(ev=ev)
    if ev in ("F1", "F3", "F4", "F5", "F6", "F8"):
        cells, nx, ny, T = seatmap_cells(ev)
        if ev == "F5" and hh:
            drop = {(X33["ABCFGH".index(s[-1])], dg.TOYO_ROWS.index(int(s[:-1]))) for s in ["14F", "14G", "14H", "25G"]}
            cells = [c for c in cells if (c[0], c[1]) not in drop]
        src = [(x, y) for x, y, s in cells if s == "I"]
        sus = [(x, y, s) for x, y, s in cells if s in "oCP"]
        ry = np.array([y for x, y, s in sus]); sy = np.array([y for x, y in src])
        dist = np.abs(ry[:, None] - sy[None, :]).min(1)
        bins = np.searchsorted(np.asarray(edges), dist, side="left")
        ycase = np.array([s == "C" for x, y, s in sus])
        rec = np.array([y * nx + x for x, y, s in sus])
        kn = kappa_nodes("air")
        nb = len(edges) + 1
        g.update(nx=nx, ny=ny, ax=0.5, ay=0.8, T=T, nseg=1, srcsets=[src], kappas=kn,
                 cfgs=[(0, k, rec) for k in range(len(kn))], bins=bins,
                 obs=[int(ycase[bins == b].sum()) for b in range(nb)], dist=dist, ycase=ycase,
                 near_bins=[b for b in range(nb) if (edges[b] if b < len(edges) else 99) <= 2],
                 xy=np.array([(x * 0.5, y * 0.8) for x, y, s in sus]), srcxy=[np.array(src) * [0.5, 0.8]])
    elif ev == "F2":
        nx, ny = 12, 14
        kn = kappa_nodes("air")
        srcsets, cfgs = [], []
        for i in range(n_cfg or 200):
            r0 = int(rng.integers(2, 12)); x0 = int(rng.choice(X343))
            near = [(x, r) for r in range(r0 - 2, r0 + 3) for x in X343 if (x, r) != (x0, r0)]
            far = [(x, r) for r in range(ny) if abs(r - r0) > 2 for x in X343]
            a = [near[k] for k in rng.choice(len(near), 13, replace=False)]
            b = [far[k] for k in rng.choice(len(far), 55, replace=False)]
            if [(x0, r0)] not in srcsets:
                srcsets.append([(x0, r0)])
            cfgs.append((srcsets.index([(x0, r0)]), i % 7, np.array([r * nx + x for x, r in a + b])))
        g.update(nx=nx, ny=ny, ax=0.5, ay=0.8, T=8.75, nseg=1, srcsets=srcsets, kappas=kn, cfgs=cfgs,
                 bins=np.array([0] * 13 + [1] * 55), obs=[4, 2], near_bins=[0])
    elif ev == "F7":
        nx, ny = 12, 34
        Lx = dict(zip("ABCDEFGHJK", X343)); row = lambda r: r - 17
        src = [(Lx["G"], row(26)), (Lx["D"], row(26))]
        near = ["26A", "26C", "27D", "27K", "24C", "24D", "24E", "24F", "24G", "28A", "28D", "28G", "28K"]
        far = [(x, row(r)) for r in list(range(17, 23)) + list(range(31, 51)) for x in X343]
        kn = kappa_nodes("air")
        cfgs = []
        for i in range(n_cfg or 200):
            b = [far[k] for k in rng.choice(len(far), 71, replace=False)]
            rec = [row(int(s[:-1])) * nx + Lx[s[-1]] for s in near] + [y * nx + x for x, y in b]
            cfgs.append((0, i % 7, np.array(rec)))
        g.update(nx=nx, ny=ny, ax=0.5, ay=0.8, T=18.0, nseg=1, srcsets=[src], kappas=kn, cfgs=cfgs,
                 bins=np.array([0] * 4 + [1] * 9 + [2] * 71), obs=[1, 3, 0], near_bins=[0, 1])
    elif ev in ("W1", "W1x", "W2"):
        nx, ny = 21, 16
        bed = {}
        for i, b in enumerate(["13", "14", "15", "16", "16x"]): bed[b] = (2 * i, 1)
        for i, b in enumerate(["12", "11", "10", "9", "9x"]): bed[b] = (2 * i, 5)
        for i, b in enumerate(["17x", "17", "18", "19", "20"]): bed[b] = (12 + 2 * i, 1)
        for i, b in enumerate(["24x", "24", "23", "22", "21"]): bed[b] = (12 + 2 * i, 5)
        for i, b in enumerate(["5", "6", "7", "8"]): bed[b] = (2 * i, 11)
        for i, b in enumerate(["4", "3", "2", "1", "1x"]): bed[b] = (2 * i, 15)
        for i, b in enumerate(["34", "33"]): bed[b] = (2 * i, 8)
        for i, b in enumerate(["25x", "25", "26", "27", "28"]): bed[b] = (12 + 2 * i, 11)
        for i, b in enumerate(["32x", "32", "31", "30", "29"]): bed[b] = (12 + 2 * i, 15)
        same = ["9", "9x", "10", "12", "13", "14", "15", "16", "16x"]
        adj = ["17", "18", "19", "20", "21", "22", "23", "24", "17x", "24x"]
        distb = [b for b in bed if b not in same + adj + ["11"]]
        st = lambda b: bed[b][1] * nx + bed[b][0]
        kn = kappa_nodes("ward")
        if ev == "W1":
            cub = ["9", "9x", "13", "14", "15", "16", "16x"]; oth = adj + distb
            cfgs = []
            for i in range(n_cfg or 200):
                rec = [st("12")] * 3 + [st(cub[k]) for k in rng.integers(0, len(cub), 8)] + [st(oth[k]) for k in rng.integers(0, len(oth), 8)]
                cfgs.append((0, i % 7, np.array(rec)))
            g.update(T=40 / 60, nseg=1, cfgs=cfgs, bins=np.array([0] * 3 + [1] * 8 + [2] * 8), obs=[3, 4, 0], near_bins=[0, 1])
        elif ev == "W1x":      # exact bed assignments read by the re-check from Wong 2004 Fig. 4
            rec = [st("12")] * 3 + [st(b) for b in ["14", "15", "15", "16", "16x", "16x", "16x", "9"]] \
                + [st(b) for b in ["17x", "17x", "17x", "23", "23", "25x", "4", "30"]]
            g.update(T=40 / 60, nseg=1, cfgs=[(0, k, np.array(rec)) for k in range(7)],
                     bins=np.array([0] * 3 + [1] * 8 + [2] * 8), obs=[3, 4, 0], near_bins=[0, 1])
        else:
            rec = [st(same[i % len(same)]) for i in range(20)] + [st(adj[i % len(adj)]) for i in range(21)] + [st(distb[i % len(distb)]) for i in range(33)]
            g.update(T=100.0, nseg=1, cfgs=[(0, k, np.array(rec)) for k in range(7)],
                     bins=np.array([0] * 20 + [1] * 21 + [2] * 33), obs=[13, 11, 6], near_bins=[0])
        g.update(nx=nx, ny=ny, ax=1.0, ay=1.0, srcsets=[[bed["11"]]], kappas=kn)
    elif ev == "P1":
        nx, ny = 12, 32
        sizes = [9, 17, 22, 30]; ed = [0.0, 4.0, 8.0, 12.0, 29.0]
        xs, ys = np.meshgrid(np.arange(nx), np.arange(ny), indexing="xy"); xs, ys = xs.ravel(), ys.ravel()
        kn = kappa_nodes("plant")
        srcsets, cfgs, dd = [], [], []
        for i in range(n_cfg or 300):
            x0, y0 = int(rng.integers(0, 2)), int(rng.integers(0, 16))
            d = np.hypot(xs - x0, ys - y0)
            rec = []
            for b in range(4):
                cand = np.flatnonzero((d > ed[b]) & (d <= ed[b + 1]))
                rec += list(cand[rng.integers(0, len(cand), sizes[b])])
            if [(x0, y0)] not in srcsets:
                srcsets.append([(x0, y0)])
            cfgs.append((srcsets.index([(x0, y0)]), i % 7, np.array(rec)))
            dd.append(d[np.array(rec)])
        g.update(nx=nx, ny=ny, ax=1.0, ay=1.0, T=8.0, nseg=3, srcsets=srcsets, kappas=kn, cfgs=cfgs,
                 bins=np.repeat(np.arange(4), sizes), obs=[5, 12, 1, 2], near_bins=[0, 1], dists=dd)
    else:
        raise ValueError(ev)
    nb = int(g["bins"].max()) + 1
    g["sizes"] = [int((g["bins"] == b).sum()) for b in range(nb)]
    g["K"] = int(sum(g["obs"]))
    return g


def m2_table(g, Dgrid=None, cache=None):
    """pmf[iD, comp] of M2 averaged over configurations/kappa nodes; also returns comps."""
    Dgrid = 10.0 ** LOGD if Dgrid is None else Dgrid
    if cache and os.path.exists(cache):
        z = np.load(cache)
        return z["pmf"], z["comps"]
    X = kern_table(g["nx"], g["ny"], g["ax"], g["ay"], g["srcsets"], g["T"], g["kappas"], Dgrid) * g["nseg"]
    comps = compositions(g["sizes"], g["K"])
    pmf = np.zeros((len(Dgrid), len(comps)))
    for i in range(len(Dgrid)):
        for (iS, ik, rec) in g["cfgs"]:
            pmf[i] += cond_pmf(X[i, ik, iS, rec], g["bins"], g["sizes"], g["K"], comps)
    pmf /= len(g["cfgs"])
    if cache:
        np.savez_compressed(cache, pmf=pmf, comps=comps, Xmin=X.min(), Xmax=X.max())
    return pmf, comps


def crr_table(g, comps, near_bins=None, true=False):
    near = np.isin(g["bins"], g["near_bins"] if near_bins is None else near_bins)
    out = np.zeros((len(LOGRHO), len(comps)))
    for i, lr in enumerate(LOGRHO):
        if true:
            out[i] = cond_pmf_trueCRR(near, 10.0 ** lr, g["bins"], g["sizes"], g["K"], comps)
        else:
            out[i] = cond_pmf(np.where(near, 10.0 ** lr, 1.0), g["bins"], g["sizes"], g["K"], comps)
    return out


def m0_pmf(g, comps):
    return cond_pmf(np.ones(len(g["bins"])), g["bins"], g["sizes"], g["K"], comps)


def idx_of(comps, obs):
    return int(np.flatnonzero((comps == np.asarray(obs)).all(1))[0])
