"""Checks with separately written code of the EXT-F numbers: exposure by the closed form B^-2 (exp(BT) - I - BT) with a dense matrix exponential
(not the mode sum), conditional distribution by brute-force enumeration / Monte Carlo.  Writes results/E3_ext_verify.json."""
import os, sys, itertools, json
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(os.path.dirname(os.path.dirname(HERE)), "bsc_validation2", "a3_tracer_physics")
sys.path.insert(0, os.path.join(ROOT, "src"))
import a3lib as A
from a3lib import L
sys.path.insert(0, os.path.normpath(os.path.join(ROOT, "..", "a1_more_outbreaks", "src")))
from a1 import events as E1
T = np.load(os.path.join(A.RES, "E1_ext_tables.npz")); out = {}
kap_nodes = 30.0 + 6.0 * (np.arange(7) + 0.5) / 7 + 0.93
def X_expm(lat, D, kap, Tt, src, rec):
    """X = int_0^T (T - s) exp((Q - kappa) s) ds v = B^-2 (exp(B T) - I - B T) v, dense linear algebra."""
    B = L.generator(lat, D) - kap * np.eye(lat.M); v = np.zeros(lat.M); np.add.at(v, np.atleast_1d(src), 1.0)
    y = expm(B * Tt) @ v - v - Tt * (B @ v)
    return np.linalg.solve(B, np.linalg.solve(B, y))[rec]
def probs(X, K):
    eps = np.exp(brentq(lambda le: np.sum(-np.expm1(-np.exp(le) * X)) - K, -60, 60)); return -np.expm1(-eps * X)
# F8: K = 2, exact enumeration over pairs
meta, draws = E1.build("F8"); d = draws[0]; bins = d["bins"]; sizes = meta["sizes"]; comps = T["F8|comps"]
acc = np.zeros(len(comps)); xdiff = 0
for kap in kap_nodes:
    for D in (707.945784384138, 112.2018454301963):
        X = X_expm(d["lat"], D, kap, meta["T"], d["src"], d["rec"])
        dd = dict(d); dd["kappa"] = kap
        from a1 import engine as G
        Xm, _ = G.x_m2(dd, D, "M2"); xdiff = max(xdiff, float(np.max(np.abs(X / Xm - 1))))
        p = probs(X, 2); o = p / (1 - p); pm = np.zeros(len(comps)); tot = 0
        for i, j in itertools.combinations(range(len(X)), 2):
            v = np.bincount([bins[i], bins[j]], minlength=len(sizes)); k = int(np.flatnonzero((comps == v).all(axis=1))[0])
            pm[k] += o[i] * o[j]; tot += o[i] * o[j]
        acc += pm / tot
acc /= 14
out["F8_exposure_expm_vs_modesum_max_rel"] = xdiff
out["F8_P_pmf_enum_vs_table_maxdiff"] = float(np.max(np.abs(acc - T["F8|P"])))
# F5: Monte Carlo (rejection on the total) for D = 708, middle node
meta, draws = E1.build("F5"); d = draws[0]; bins = d["bins"]; comps = T["F5|comps"]; K = meta["K"]
rng = np.random.default_rng(31)
X = X_expm(d["lat"], 707.945784384138, kap_nodes[3], meta["T"], d["src"], d["rec"]); p = probs(X, K)
cnt = np.zeros(len(comps)); n = 0
key = {tuple(c): i for i, c in enumerate(comps)}
for _ in range(40):
    dr = rng.random((50000, len(X))) < p; keep = dr[dr.sum(axis=1) == K]
    v = np.stack([keep[:, bins == b].sum(axis=1) for b in range(len(meta["sizes"]))], axis=1)
    for row in v: cnt[key[tuple(row)]] += 1
    n += len(v)
from a1 import engine as G
dd = dict(d); dd["kappa"] = kap_nodes[3]; ref = G.cond_pmf(G.x_m2(dd, 707.945784384138, "M2")[0], bins, meta["sizes"], K, comps)
out["F5_mc_n"] = int(n); out["F5_D708_node3_mc_vs_condpmf_maxdiff"] = float(np.max(np.abs(cnt / n - ref)))
obs = np.array(meta["obs"]); i = key[tuple(obs)]; out["F5_obs_prob_mc_vs_condpmf"] = [float(cnt[i] / n), float(ref[i])]
# M0 by hand: multivariate hypergeometric for F7
from scipy.special import comb
s7 = [4, 9, 71]; o7 = [1, 3, 0]; out["F7_M0_logp_by_hand"] = float(np.log(comb(4, 1) * comb(9, 3) * comb(71, 0) / comb(84, 4)))
ev = A.load_json("E2_ext_evaluation.json"); out["F7_M0_logp_stored"] = ev["events"]["F7"]["models"]["M0"]["logp"]
# CRR3, F8 by hand (K = 2): odds weights
w = np.where(np.isin(E1.build("F8")[1][0]["bins"], [0, 1]), 3.0, 1.0); eps = brentq(lambda e: np.sum(e * w) - 2, 0, 1); pp = eps * w
# note: a1's CRR uses 1 - exp(-eps X) with X = rho or 1; recompute that way
pp = probs(w, 2); o = pp / (1 - pp); e2 = lambda a: (a.sum() ** 2 - (a ** 2).sum()) / 2
b0 = E1.build("F8")[1][0]["bins"] == 0
out["F8_CRR3_logp_by_hand"] = float(np.log(e2(o[b0]) / e2(o))); out["F8_CRR3_logp_stored"] = ev["events"]["F8"]["models"]["CRR3"]["logp"]
json.dump(out, open(os.path.join(A.RES, "E3_ext_verify.json"), "w"), indent=1); print(json.dumps(out, indent=1))
