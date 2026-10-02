"""Self-tests on synthetic inputs (no outbreak outcome is read)."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from scipy import stats as st, linalg
from a1 import engine as G, lattice as L, stats as S
rng = np.random.default_rng(1)
out = {}
# 1. truncated Poisson-binomial vs full
p = rng.uniform(0, 0.4, 30)
out["pb_trunc_maxerr"] = float(np.abs(G.pb_trunc(p, 8) - S.pb_pmf(p)[:9]).max())
# 2. Newton profile vs brentq
X = rng.gamma(1.0, 1.0, 50)
out["profile_eps_relerr"] = float(abs(G.profile_eps(X, 7) / S.profile_eps(X, 7) - 1))
# 3. equal exposures -> multivariate hypergeometric
sizes = [5, 7, 20]; K = 6
comps = G.compositions(sizes, K)
bins = np.repeat(np.arange(3), sizes)
pm = G.cond_pmf(np.ones(32), bins, sizes, K, comps)
ref = np.array([st.multivariate_hypergeom.pmf(c, sizes, K) for c in comps])
out["hypergeom_maxerr"] = float(np.abs(pm - ref).max()); out["pmf_sum"] = float(pm.sum())
# 4. two bins vs round-1 cond_pmf
X = rng.gamma(1.0, 1.0, 12); b2 = np.array([0] * 5 + [1] * 7); c2 = G.compositions([5, 7], 4)
pm2 = G.cond_pmf(X, b2, [5, 7], 4, c2)
eps = S.profile_eps(X, 4); pp = -np.expm1(-eps * X)
ref2 = S.cond_pmf(pp[:5], pp[5:], 4)
out["round1_cond_maxerr"] = float(max(abs(pm2[i] - ref2[c[0]]) for i, c in enumerate(c2)))
# 5. rejection sampling check of the conditional vector pmf
X = rng.gamma(0.5, 1.0, 32); pm = G.cond_pmf(X, bins, sizes, K, comps)
eps = G.profile_eps(X, K); pp = -np.expm1(-eps * X)
sim = rng.random((400000, 32)) < pp
keep = sim[sim.sum(1) == K]
cnt = np.stack([keep[:, bins == b].sum(1) for b in range(3)], 1)
emp = np.array([(cnt == c).all(1).mean() for c in comps])
out["rejection_maxerr"] = float(np.abs(emp - pm).max()); out["rejection_n"] = int(len(keep))
# 6. M2 kernel vs quadrature of the matrix exponential; M2 -> M0 as D -> infinity; lattice mean
lat = L.Lattice(5, 6, 0.5, 0.8)
D, kap, T = 3.0, 4.0, 2.0
d = dict(lat=lat, src=np.array([7]), rec=np.arange(30), kappa=kap, segs=[T])
Xm, xbar = G.x_m2(d, D)
Q = L.generator(lat, D)
s = np.linspace(0, T, 4001)
vals = np.array([(T - si) * np.exp(-kap * si) * linalg.expm(Q * si)[:, 7] for si in s[::20]])
ref = np.trapezoid(vals, s[::20], axis=0)
out["m2_quadrature_relerr"] = float(np.abs(Xm / ref - 1).max())
out["m2_mean_vs_xbar"] = float(abs(Xm.mean() / xbar - 1))
Xinf, xb = G.x_m2(d, 1e9)
out["m2_Dinf_spread"] = float(Xinf.max() / Xinf.min() - 1)
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(ROOT, "out", "00_selftest.json"), "w"), indent=1)
