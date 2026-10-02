"""Checks with separately written code of the A3 numbers (different numerical routes, no import of a3lib prediction code
except the event geometry tables).  Writes results/05_verify.json."""
import itertools, json, os, sys
import numpy as np, pandas as pd
from scipy import stats
from scipy.special import comb
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "bsc_validation2", "a3_tracer_physics")
sys.path.insert(0, os.path.join(ROOT, "src"))
import a3lib as A
from a3lib import E, L
out = {}
rng = np.random.default_rng(777)
pred = A.load_json("02_predictions.json")["events"]; ker = A.load_json("01_kernels.json")

# 1. steady M2 field: mode sum vs dense linear solve (kappa I - Q) c = delta
ev = A.t2_event(); lat = ev["lat"]
D, kap = 1258.9, 6.5
Q = L.generator(lat, D)
c_lin = np.linalg.solve(kap * np.eye(lat.M) - Q, np.eye(lat.M)[:, ev["src_test"]])
c_mode = A.steady_conc(lat, D, kap, ev["src_test"], np.arange(lat.M))
out["steady_field_max_rel_diff"] = float(np.max(np.abs(c_lin - c_mode) / c_lin))

# 2. coach D fit by a separate optimiser on the linear-solve field
df = ev["df"]; m = df[df["s3_measured"].notna() & (df.seat != "12D")]
y = np.log(m["s3_measured"].to_numpy()); rec = m["site"].to_numpy()
def sse(logD):
    c = np.linalg.solve(kap * np.eye(lat.M) - L.generator(lat, 10 ** logD), np.eye(lat.M)[:, ev["src_test"]])[rec]
    r = y - np.log(c); r -= r.mean(); return float(r @ r)
grid = np.arange(-1, 5.01, 0.05); v = np.array([sse(g) for g in grid])
out["coach_D_separate"] = float(10 ** grid[int(np.argmin(v))]); out["coach_D_a3"] = ker["coach"]["D_hat"]
out["coach_sse_at_min_and_wellmixed"] = [float(v.min()), float(v[-1])]

# 3. T5: exact conditional distribution by enumerating all C(20,12) case sets (K-A 5L and 5A-mirrored)
ev5, kern = A.t5_direct_kernels(False)
near = ev5["near"]
def enum_pmf(X, K=12):
    from scipy.optimize import brentq
    X = np.clip(X, 1e-300, None)
    eps = np.exp(brentq(lambda le: np.sum(-np.expm1(-np.exp(le) * X)) - K, -40, 40))
    p = -np.expm1(-eps * X); odds = p / (1 - p)
    pmf = np.zeros(near.sum() + 1); tot = 0.0
    idx = np.arange(len(X))
    for s in itertools.combinations(idx, K):
        w = np.prod(odds[list(s)]); pmf[int(near[list(s)].sum())] += w; tot += w
    return pmf / tot
pm = np.mean([enum_pmf(v[0]) for v in kern.values()], axis=0)
out["T5_P_pmf_enum_vs_a3_maxdiff"] = float(np.max(np.abs(pm - np.array(pred["T5"]["pmf"]["P"]))))
out["T5_P_prob_obs_enum"] = float(pm[11])

# 4. T2: Monte Carlo of the conditional distribution (rejection on the total) for the noise-free K-A kernel
X = ev["sus"]["s4_cfd"].to_numpy()
from scipy.optimize import brentq
eps = np.exp(brentq(lambda le: np.sum(-np.expm1(-np.exp(le) * X)) - 7, -40, 40))
p = -np.expm1(-eps * X)
draws = rng.random((400000, len(X))) < p
keep = draws.sum(axis=1) == 7
kn = draws[keep][:, ev["near"]].sum(axis=1)
mc = np.bincount(kn, minlength=20)[:len(pred["T2"]["pmf"]["KA_nonoise"])] / keep.sum()
out["T2_KA_pmf_mc_vs_a3_maxdiff"] = float(np.max(np.abs(mc - np.array(pred["T2"]["pmf"]["KA_nonoise"])))); out["T2_mc_n"] = int(keep.sum())
out["T2_M0_hypergeom_p_obs"] = float(stats.hypergeom.pmf(3, 45, 19, 7))

# 5. R1: closed form for K = 2 (two cases, both at near tables)
tb = A.r1_tables(); Xt = (tb["tracer"] * tb["T_h"]).to_numpy(); npat = tb["patrons"].to_numpy(); nearT = tb["neighbour_class"].eq("immediate").to_numpy()
Xp = np.repeat(Xt, npat); nP = np.repeat(nearT, npat)
eps = np.exp(brentq(lambda le: np.sum(-np.expm1(-np.exp(le) * Xp)) - 2, -40, 40)); o = (lambda q: q / (1 - q))(-np.expm1(-eps * Xp))
e2 = lambda a: (a.sum() ** 2 - (a ** 2).sum()) / 2
out["R1_K2_prob_both_near_closed_form"] = float(e2(o[nP]) / e2(o)); out["R1_K2_a3"] = pred["R1"]["pmf"]["P"][2]
out["R1_K2_M0"] = float(comb(16, 2) / comb(79, 2))

# 6. CRR pmf against direct binomial product
n1, n2, K, rho = 12, 8, 12, 3.0
q = K / (n1 * rho + n2); pn = min(1.0, rho * q)
if rho * q > 1: q = (K - n1) / n2
w = np.array([stats.binom.pmf(k, n1, pn) * stats.binom.pmf(K - k, n2, q) for k in range(n1 + 1)]); w /= w.sum()
out["CRR3_T5_maxdiff"] = float(np.max(np.abs(w - np.array(pred["T5"]["pmf"]["CRR3"]))))

# 7. T4 observed deviance for M0 by hand
d4 = E.t4_data(); mask = ~((d4["dr"] == 0) & (d4["dc"] == 1)).to_numpy(); k = d4["k"].to_numpy()[mask]; n = d4["n"].to_numpy()[mask]
ex = k.sum() * n / n.sum(); out["T4_M0_deviance_by_hand"] = float(2 * np.sum(np.where(k > 0, k * np.log(np.where(k > 0, k, 1) / ex), 0)))
out["T4_total_and_exposed"] = [int(k.sum()), int(n.sum())]

# 8. exact p-value of the pooled difference P - M0 (full precision)
ev_ = ["T2", "T1", "T5", "C1", "R1"]; obs = {"T2": 3, "T1": 14, "T5": 11, "C1": 8, "R1": 2}
pmP = {e: np.array(pred[e]["pmf"]["P"]) for e in ev_}; pm0 = {e: np.array(pred[e]["pmf"]["M0"]) for e in ev_}
d_obs = sum(np.log(pmP[e][obs[e]]) - np.log(pm0[e][obs[e]]) for e in ev_); tot = 0.0
for combo in itertools.product(*[np.flatnonzero(pm0[e] > 0) for e in ev_]):
    dd = sum(np.log(max(pmP[e][k], 1e-300)) - np.log(pm0[e][k]) for e, k in zip(ev_, combo))
    if dd >= d_obs - 1e-12: tot += np.prod([pm0[e][k] for e, k in zip(ev_, combo)])
out["pooled_P_vs_M0"] = dict(delta=float(d_obs), p_exact=float(tot))
# without the flight, and without the classroom
for drop in ("T5", "C1", "R1"):
    ee = [e for e in ev_ if e != drop]
    d_o = sum(np.log(pmP[e][obs[e]]) - np.log(pm0[e][obs[e]]) for e in ee); tot = 0.0
    for combo in itertools.product(*[np.flatnonzero(pm0[e] > 0) for e in ee]):
        dd = sum(np.log(max(pmP[e][k], 1e-300)) - np.log(pm0[e][k]) for e, k in zip(ee, combo))
        if dd >= d_o - 1e-12: tot += np.prod([pm0[e][k] for e, k in zip(ee, combo)])
    out["pooled_P_vs_M0_without_" + drop] = dict(delta=float(d_o), p_exact=float(tot))
A.save_json("05_verify.json", out)
print(json.dumps(out, indent=1))
