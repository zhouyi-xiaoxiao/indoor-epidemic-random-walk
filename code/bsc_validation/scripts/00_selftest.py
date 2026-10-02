"""Software verification on synthetic inputs only (no outbreak count is read).

A  closed-form kernels against brute-force quadrature of the matrix exponential
B  nesting: M2 -> M0 as D -> infinity; spatial mean of M2 equals M0
C  conditional stratum distribution against enumeration and the hypergeometric law
D  join-count moments against permutation
E  M1 closed form against the verified exact simulator (two people, then a full cabin)
"""
import itertools

import numpy as np
from scipy.linalg import expm

import _common as C
from valmod import lattice as L
from valmod import stats as S
from valmod.exact import simulate_m1

out = {"stamp": C.stamp()}
rng = np.random.default_rng(C.SEED)

# ---- A ---------------------------------------------------------------------
lat = L.Lattice(nx=4, ny=5, ax=0.5, ay=0.9)
D, kappa, T = 0.8, 2.5, 1.7
Q = L.generator(lat, D)
xs, ws = np.polynomial.legendre.leggauss(120)
t = 0.5 * T * (xs + 1)
w = 0.5 * T * ws
src = 6
X1 = np.zeros(lat.M)
X2 = np.zeros(lat.M)
for ti, wi in zip(t, w):
    X1 += wi * expm(Q * 2 * ti)[:, src]
    X2 += wi * (T - ti) * np.exp(-kappa * ti) * expm(Q * ti)[:, src]
rec = np.arange(lat.M)
e1 = np.abs(L.exposure(lat, "M1", D, kappa, [T], src, rec) - X1).max()
e2 = np.abs(L.exposure(lat, "M2", D, kappa, [T], src, rec) - X2).max()
out["A_max_abs_err_M1"] = float(e1)
out["A_max_abs_err_M2"] = float(e2)
assert e1 < 1e-9 and e2 < 1e-9, (e1, e2)

# ---- B ---------------------------------------------------------------------
x0 = L.exposure(lat, "M0", D, kappa, [T], src, rec)
xinf = L.exposure(lat, "M2", 1e9, kappa, [T], src, rec)
xm = L.exposure(lat, "M2", D, kappa, [T], src, rec)
out["B_M2_Dinf_minus_M0"] = float(np.abs(xinf - x0).max())
out["B_mean_M2_minus_M0"] = float(abs(xm.mean() - x0[0]))
assert out["B_M2_Dinf_minus_M0"] < 1e-8 and out["B_mean_M2_minus_M0"] < 1e-12

# ---- C ---------------------------------------------------------------------
p1 = rng.uniform(0.05, 0.9, 5)
p2 = rng.uniform(0.05, 0.9, 6)
K = 4
pmf = S.cond_pmf(p1, p2, K)
brute = np.zeros(6)
p = np.concatenate([p1, p2])
for bits in itertools.product([0, 1], repeat=11):
    b = np.array(bits)
    if b.sum() == K:
        brute[b[:5].sum()] += np.prod(np.where(b == 1, p, 1 - p))
brute /= brute.sum()
out["C_cond_pmf_err"] = float(np.abs(pmf - brute).max())
ph = S.cond_pmf(np.full(7, 0.3), np.full(9, 0.3), 5)
out["C_hypergeom_err"] = float(np.abs(ph - S.hypergeom_pmf(7, 9, 5)).max())
eps = S.profile_eps(np.array([1.0, 2.0, 0.5, 3.0]), 2.2)
out["C_profile_eps_resid"] = float(np.sum(1 - np.exp(-eps * np.array([1.0, 2.0, 0.5, 3.0]))) - 2.2)
assert out["C_cond_pmf_err"] < 1e-12 and out["C_hypergeom_err"] < 1e-12

# ---- D ---------------------------------------------------------------------
n = 30
A = (rng.random((n, n)) < 0.15).astype(float)
A = np.triu(A, 1)
A = A + A.T
k = 12
m, v = S.joincount_moments(A, k)
sims = []
for _ in range(40000):
    c = np.zeros(n)
    c[rng.choice(n, k, replace=False)] = 1
    sims.append(c @ A @ c / 2)
sims = np.array(sims)
out["D_mean_theory_sim"] = [float(m), float(sims.mean())]
out["D_var_theory_sim"] = [float(v), float(sims.var())]
assert abs(m - sims.mean()) < 0.05 and abs(v - sims.var()) / v < 0.05

# ---- E ---------------------------------------------------------------------
# two people on a 4 x 5 lattice: N_r = 2 whenever they meet, pair hazard beta/2
lat = L.Lattice(nx=4, ny=5, ax=0.5, ay=0.9)
src, recs = 6, np.array([7])
T, D = 2.0, 0.4
rows = []
for beta in (0.2, 2.0, 20.0):
    X = L.exposure(lat, "M1", D, 0.0, [T], src, recs)[0]
    p_cf = 1 - np.exp(-0.5 * beta * X)
    p_sim, _ = simulate_m1(lat, src, recs, beta, D, T, n_rep=60000, seed=11)
    se = np.sqrt(p_sim[0] * (1 - p_sim[0]) / 60000)
    rows.append(dict(beta=beta, closed_form=float(p_cf), exact=float(p_sim[0]), se=float(se)))
out["E_pair"] = rows
# small-hazard limit must agree; large hazard: closed form is an upper bound
assert abs(rows[0]["closed_form"] - rows[0]["exact"]) < 4 * rows[0]["se"] + 0.002
assert rows[2]["closed_form"] >= rows[2]["exact"] - 4 * rows[2]["se"]

# full cabin, 30 people on a 5 x 8 lattice, per-cell rule beta q n_I / n_total
lat = L.Lattice(nx=5, ny=8, ax=0.5, ay=0.8)
seats = rng.choice(lat.M, 30, replace=False)
src, recs = int(seats[0]), seats[1:]
T, D, beta = 1.5, 3.0, 1.0
X = L.exposure(lat, "M1", D, 0.0, [T], src, recs)
p_sim, inf = simulate_m1(lat, src, recs, beta, D, T, n_rep=20000, seed=12)
eps_fit = S.profile_eps(X, p_sim.sum())
p_cf = 1 - np.exp(-eps_fit * X)
out["E_cabin"] = dict(total_exact=float(p_sim.sum()), eps_over_beta=float(eps_fit / beta),
                      max_abs_diff=float(np.abs(p_cf - p_sim).max()),
                      corr=float(np.corrcoef(p_cf, p_sim)[0, 1]))

C.save_json("00_selftest.json", out)
for k_, v_ in out.items():
    print(k_, v_)
print("SELFTEST OK")
