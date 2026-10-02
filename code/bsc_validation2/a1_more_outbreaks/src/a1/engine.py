"""Conditional (given the event total) predictive distributions of bin-count vectors."""
from __future__ import annotations
import itertools
import numpy as np
from . import lattice as L

LOGD = np.round(np.linspace(-2.0, 4.0, 61), 6)
LOGRHO = np.round(np.linspace(0.0, 2.5, 51), 6)
WGRID = np.array([0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95, 1.0])
LOGD3 = LOGD[::2]


def compositions(sizes, K):
    """All vectors k with sum K and 0 <= k_b <= sizes[b]; array (n_comp, n_bins)."""
    rng = [range(0, min(s, K) + 1) for s in sizes[:-1]]
    out = []
    for head in itertools.product(*rng):
        last = K - sum(head)
        if 0 <= last <= sizes[-1]:
            out.append(head + (last,))
    return np.array(out, dtype=int)


def pb_trunc(p, K):
    """First K+1 terms of the Poisson-binomial pmf of independent Bernoulli(p_i)."""
    pmf = np.zeros(K + 1)
    pmf[0] = 1.0
    for pi in p:
        pmf[1:] = pmf[1:] * (1.0 - pi) + pmf[:-1] * pi
        pmf[0] *= (1.0 - pi)
    return pmf


def profile_eps(X, K):
    """eps with sum(1 - exp(-eps X)) = K (Newton from below; f is concave increasing)."""
    X = np.asarray(X, float)
    n = len(X)
    if K <= 0:
        return 0.0
    if K >= n:
        return np.inf
    eps = K / X.sum()
    for _ in range(200):
        e = np.exp(-eps * X)
        f = n - e.sum() - K
        if abs(f) < 1e-10 * max(K, 1):
            break
        fp = (X * e).sum()
        if fp <= 0 or not np.isfinite(fp):
            break
        step = -f / fp
        eps += step
        if step < 1e-14 * eps:
            break
    return eps


def cond_pmf(X, bins, sizes, K, comps):
    """P(bin vector | total K) over `comps`, for exposures X (eps profiled)."""
    X = np.clip(X, 1e-300, None)
    eps = profile_eps(X, K)
    if not np.isfinite(eps):
        p = np.ones_like(X)
    else:
        p = -np.expm1(-eps * X)
    p = np.clip(p, 0.0, 1.0 - 1e-15)
    pr = np.ones(len(comps))
    for b in range(len(sizes)):
        pb = pb_trunc(p[bins == b], K)
        pr *= pb[comps[:, b]]
    s = pr.sum()
    return pr / s if s > 0 else np.full(len(comps), 1.0 / len(comps))


def x_m2(d, D, model="M2"):
    lat = d["lat"]
    w = np.zeros(lat.M)
    for T in d["segs"]:
        w += L.weights(lat, model, D, d["kappa"], float(T))
    v = (lat.Phi[d["src"]] * w).sum(axis=0)
    X = lat.Phi[d["rec"]] @ v
    xbar = len(d["src"]) * w[0] / lat.M            # lattice average (k = 0 mode)
    return X, xbar


def table(meta, draws, model, comps):
    """Predictive pmf over comps at every grid value of the model's parameter(s)."""
    sizes, K = meta["sizes"], meta["K"]
    nd = len(draws)
    if model == "M0":
        return cond_pmf(np.ones(sum(sizes)), draws[0]["bins"], sizes, K, comps)[None, :]
    if model == "CRR":
        near = np.isin(draws[0]["bins"], meta["near_bins"])
        out = np.zeros((len(LOGRHO), len(comps)))
        for i, lr in enumerate(LOGRHO):
            out[i] = cond_pmf(np.where(near, 10.0 ** lr, 1.0), draws[0]["bins"], sizes, K, comps)
        return out
    if model in ("M2", "M1"):
        out = np.zeros((len(LOGD), len(comps)))
        for i, ld in enumerate(LOGD):
            for d in draws:
                X, _ = x_m2(d, 10.0 ** ld, model)
                out[i] += cond_pmf(X, d["bins"], sizes, K, comps)
        return out / nd
    if model == "M3":
        out = np.zeros((len(LOGD3), len(WGRID), len(comps)))
        for i, ld in enumerate(LOGD3):
            for d in draws:
                X, xbar = x_m2(d, 10.0 ** ld, "M2")
                for j, w in enumerate(WGRID):
                    out[i, j] += cond_pmf((1 - w) * X + w * xbar, d["bins"], sizes, K, comps)
        return out / nd
    raise ValueError(model)


def seat_logscore(meta, draws, model, theta):
    """log P(observed set of infected seats | K), averaged over kappa nodes (seat-map events)."""
    y, K = meta["y"], meta["K"]
    vals = []
    for d in draws:
        if model == "M0":
            X = np.ones(len(y))
        elif model == "CRR":
            X = np.where(np.isin(d["bins"], meta["near_bins"]), theta, 1.0)
        elif model == "M3":
            X0, xbar = x_m2(d, theta[0], "M2")
            X = (1 - theta[1]) * X0 + theta[1] * xbar
        else:
            X, _ = x_m2(d, theta, model)
        X = np.clip(X, 1e-300, None)
        eps = profile_eps(X, K)
        p = np.clip(-np.expm1(-eps * X), 1e-300, 1 - 1e-15)
        ll = np.log(p[y]).sum() + np.log1p(-p[~y]).sum() - np.log(pb_trunc(p, K)[K])
        vals.append(ll)
    vals = np.array(vals)
    return float(np.log(np.mean(np.exp(vals - vals.max()))) + vals.max())
