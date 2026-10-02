"""Which interval method reproduces the CIs printed in Hu et al. Table 1? (vectorised search)"""
import numpy as np, json
from scipy import stats
from pathlib import Path
Z = 1.959964
def wilson(k, n):
    p = k / n; d = 1 + Z * Z / n; c = (p + Z * Z / (2 * n)) / d
    h = Z * np.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d; return c - h, c + h
def cp(k, n):
    return stats.beta.ppf(.025, k, n - k + 1), stats.beta.ppf(.975, k + 1, n - k)
def wald(k, n):
    p = k / n; s = Z * np.sqrt(p * (1 - p) / n); return np.maximum(0, p - s), p + s
def logit(k, n):
    p = k / n; se = np.sqrt(1 / k + 1 / (n - k)); l = np.log(p / (1 - p))
    return 1 / (1 + np.exp(-(l - Z * se))), 1 / (1 + np.exp(-(l + Z * se)))
def logwald(k, n):
    p = k / n; se = np.sqrt(1 / k - 1 / n); return p * np.exp(-Z * se), p * np.exp(Z * se)
METH = dict(wilson=wilson, cp=cp, wald=wald, logit=logit, logwald=logwald)
cells = {"r0a(3.53)": (3.53, 2.89, 4.31), "r0b(1.65)": (1.65, 1.18, 2.31), "r0c(0.38)": (0.38, 0.18, 0.78), "r0d(0.38)": (0.38, 0.19, 0.79),
         "r0e(0.29)": (0.29, 0.10, 0.85), "r1c0": (0.21, 0.11, 0.38), "r1c4": (0.03, 0.00, 0.16), "r1c5": (0.05, 0.00, 0.30),
         "r3c5": (0.06, 0.00, 0.36), "row0": (1.53, 1.30, 1.80), "overall": (0.32, 0.28, 0.36), "col1": (0.68, 0.56, 0.81)}
res = {}
for name, (p, lo, hi) in cells.items():
    res[name] = {}
    for mname, m in METH.items():
        best = None
        for k in range(1, 320):
            n_lo = int(np.floor(k / ((p + 0.005) / 100))); n_hi = int(np.ceil(k / (max(p - 0.005, 1e-4) / 100)))
            n = np.arange(max(n_lo, k + 1), n_hi + 1, dtype=float)
            if len(n) == 0 or len(n) > 40000: continue
            ok = np.round(100 * k / n, 2) == p
            if not ok.any(): continue
            n = n[ok]; l, h = m(float(k), n)
            err = np.abs(100 * l - lo) + np.abs(100 * h - hi); j = int(err.argmin())
            if best is None or err[j] < best[0]: best = (float(err[j]), k, int(n[j]), float(100 * l[j]), float(100 * h[j]))
        res[name][mname] = best
        print(name, mname, best, flush=True)
json.dump(res, open(Path(__file__).resolve().parents[1] / "results" / "train_ci_method_exploration.json", "w"), indent=1)
