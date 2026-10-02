"""Statistics for the validation: conditional (given the event total) stratum
distributions, exact binomial limits, predictive p-values."""
from __future__ import annotations

import numpy as np
from scipy import optimize, stats


def pb_pmf(p: np.ndarray) -> np.ndarray:
    """Poisson-binomial pmf of sum of independent Bernoulli(p_i)."""
    pmf = np.zeros(len(p) + 1)
    pmf[0] = 1.0
    for i, pi in enumerate(p):
        pmf[1:i + 2] = pmf[1:i + 2] * (1 - pi) + pmf[0:i + 1] * pi
        pmf[0] *= (1 - pi)
    return pmf


def profile_eps(X: np.ndarray, K: float) -> float:
    """eps such that sum_j (1 - exp(-eps X_j)) = K."""
    X = np.asarray(X, float)
    n = len(X)
    if K <= 0:
        return 0.0
    if K >= n:
        return np.inf
    f = lambda le: np.sum(-np.expm1(-np.exp(le) * X)) - K
    lo = np.log(K / X.sum()) - 1.0
    hi = lo + 2.0
    while f(hi) < 0:
        hi += 2.0
        if hi > 200:
            return np.inf
    while f(lo) > 0:
        lo -= 2.0
    return float(np.exp(optimize.brentq(f, lo, hi, xtol=1e-10)))


def cond_pmf(p_near: np.ndarray, p_far: np.ndarray, K: int) -> np.ndarray:
    """P(k_near = k | k_near + k_far = K), k = 0..len(p_near)."""
    a = pb_pmf(np.asarray(p_near))
    b = pb_pmf(np.asarray(p_far))
    n1, n2 = len(p_near), len(p_far)
    out = np.zeros(n1 + 1)
    for k in range(n1 + 1):
        if 0 <= K - k <= n2:
            out[k] = a[k] * b[K - k]
    s = out.sum()
    return out / s if s > 0 else out


def hypergeom_pmf(n1: int, n2: int, K: int) -> np.ndarray:
    k = np.arange(n1 + 1)
    return stats.hypergeom.pmf(k, n1 + n2, n1, K)


def two_sided_p(pmf: np.ndarray, k_obs: int) -> float:
    """Sum of pmf over outcomes no more probable than the observed one."""
    tol = 1e-12
    return float(pmf[pmf <= pmf[k_obs] * (1 + 1e-9) + tol].sum())


def pmf_interval(pmf: np.ndarray, level: float = 0.90):
    """Central interval (equal tails) of a discrete pmf: (lo, hi) inclusive."""
    c = np.cumsum(pmf)
    a = (1 - level) / 2
    lo = int(np.searchsorted(c, a, side="left"))
    hi = int(np.searchsorted(c, 1 - a, side="left"))
    return lo, min(hi, len(pmf) - 1)


def rr(k1, n1, k2, n2, cc: float = 0.5):
    """Risk ratio with continuity correction when a cell is zero."""
    k1 = np.asarray(k1, float)
    k2 = np.asarray(k2, float)
    zero = (k1 == 0) | (k2 == 0)
    a = np.where(zero, k1 + cc, k1)
    b = np.where(zero, k2 + cc, k2)
    m1 = np.where(zero, n1 + cc, n1)
    m2 = np.where(zero, n2 + cc, n2)
    return (a / m1) / (b / m2)


def se_log_rr(k1, n1, k2, n2, cc: float = 0.5):
    if k1 == 0 or k2 == 0:
        k1, k2, n1, n2 = k1 + cc, k2 + cc, n1 + cc, n2 + cc
    return float(np.sqrt(1 / k1 - 1 / n1 + 1 / k2 - 1 / n2))


def clopper_pearson(k: int, n: int, level: float = 0.95):
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else stats.beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a, k + 1, n - k)
    return float(lo), float(hi)


def upper_bound_one_sided(k: int, n: int, level: float = 0.95):
    return 1.0 if k == n else float(stats.beta.ppf(level, k + 1, n - k))


def binom_loglik(k, n, p):
    p = np.clip(p, 1e-300, 1 - 1e-16)
    return float(np.sum(stats.binom.logpmf(k, n, p)))


def joincount_moments(A: np.ndarray, k: int):
    """Mean and variance of the number of case-case joins when k of n nodes
    are labelled at random (non-free sampling).  A: symmetric 0/1 adjacency."""
    n = A.shape[0]
    deg = A.sum(axis=1)
    P = A.sum() / 2.0                         # number of joins
    # pairs of joins sharing one node, and disjoint pairs
    S1 = float((deg * (deg - 1)).sum()) / 2.0  # unordered pairs of joins sharing a node
    S0 = P * (P - 1) / 2.0 - S1                # unordered disjoint pairs
    p2 = k * (k - 1) / (n * (n - 1))
    p3 = p2 * (k - 2) / (n - 2)
    p4 = p3 * (k - 3) / (n - 3)
    mean = P * p2
    ex2 = P * p2 + 2 * S1 * p3 + 2 * S0 * p4
    return mean, ex2 - mean ** 2


def joincount_z(A: np.ndarray, case: np.ndarray):
    case = case.astype(float)
    jj = float(case @ A @ case) / 2.0
    m, v = joincount_moments(A, int(case.sum()))
    return (jj - m) / np.sqrt(v), jj, m, np.sqrt(v)
