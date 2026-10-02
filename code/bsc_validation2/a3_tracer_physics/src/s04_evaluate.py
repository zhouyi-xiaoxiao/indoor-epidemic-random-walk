"""Step 4: comparison with the outbreak counts (first script that reads them).
Writes results/04_evaluation.json."""
import itertools, json, os
import numpy as np, pandas as pd
import a3lib as A
from a3lib import E, L, S

pred = A.load_json("02_predictions.json")["events"]
power = A.load_json("03_power.json")
ho = json.load(open(os.path.join(A.VAL1, "data", "heldout_outcomes.json")))
obs = {"T2": 3, "T1": ho["shape_primary"]["T1"]["near_k"], "T5": 11, "C1": ho["shape_primary"]["C1"]["near_k"], "R1": 2}
assert ho["shape_primary"]["T2"]["near_k"] == 3 and ho["shape_primary"]["T5"]["near_k"] == 11
EV = ["T2", "T1", "T5", "C1", "R1"]
out = {"stamp": A.stamp(), "observed_near": obs, "events": {}}

def table(e, k, d=None):
    d = pred[e] if d is None else d
    r = {}
    for m, pmf in d["pmf"].items():
        pmf = np.array(pmf)
        lo, hi = S.pmf_interval(pmf, 0.95)
        r[m] = dict(p_two_sided=A.two_sided_p(pmf, k), logscore=float(np.log(max(pmf[k], 1e-300))), prob_obs=float(pmf[k]),
                    mean_near=float((pmf * np.arange(len(pmf))).sum()), int95=[int(lo), int(hi)], rr_expected=d["rr_expected"].get(m))
    return dict(n_near=d["n_near"], n_far=d["n_far"], K=d["K"], k_obs=int(k),
                rr_obs=float(S.rr(k, d["n_near"], d["K"] - k, d["n_far"])), models=r)

for e in EV:
    out["events"][e] = table(e, obs[e])
out["events"]["T2_side"] = table("T2_side", ho["shape_secondary"]["T2_driver_side"]["near_k"])
out["events"]["T5_no2A"] = table("T5_no2A", 11)
out["events"]["R1_byK"] = {K: table("R1", int(K), pred["R1_byK"][K]) for K in pred["R1_byK"]}

# ---------------- pooled scores and exact tests
def joint(models_by_event):
    """exact distribution of the pooled log-score difference under each truth"""
    return None
def pooled(pkey, events=EV, competitor_keys=("M0", "CRR2", "CRR3"), kobs=None):
    kobs = obs if kobs is None else kobs
    pmP = {e: np.array(pred[e]["pmf"][pkey[e] if isinstance(pkey, dict) else pkey]) for e in events}
    res = {}
    for j in competitor_keys:
        pmJ = {e: np.array(pred[e]["pmf"][j]) for e in events}
        d_obs = sum(np.log(max(pmP[e][kobs[e]], 1e-300)) - np.log(max(pmJ[e][kobs[e]], 1e-300)) for e in events)
        supp = [np.flatnonzero(pmJ[e] > 0) for e in events]
        tot = 0.0
        for combo in itertools.product(*supp):
            wj = np.prod([pmJ[e][k] for e, k in zip(events, combo)])
            dd = sum(np.log(max(pmP[e][k], 1e-300)) - np.log(pmJ[e][k]) for e, k in zip(events, combo))
            if dd >= d_obs - 1e-12: tot += wj
        res[j] = dict(delta=float(d_obs), p_exact_under_competitor=float(tot),
                      per_event={e: float(np.log(max(pmP[e][kobs[e]], 1e-300)) - np.log(max(pmJ[e][kobs[e]], 1e-300))) for e in events})
    return res
out["pooled"] = {"P": pooled("P")}
for tag in ("P_lin", "P_gam"):
    out["pooled"][tag] = pooled(tag)
out["pooled"]["S-form (K-B for T2, T5)"] = pooled({"T2": "KB_measD", "T1": "P", "T5": "KB_measD", "C1": "P", "R1": "P"})
out["pooled"]["M2cal_mode (T2,T1,T5,C1 only)"] = pooled("M2cal_mode", events=["T2", "T1", "T5", "C1"])
out["pooled"]["P (T2,T1,T5,C1 only)"] = pooled("P", events=["T2", "T1", "T5", "C1"])
out["pooled"]["CRR2 vs M0"] = pooled("CRR2", competitor_keys=("M0",))
out["pooled"]["CRR3 vs M0"] = pooled("CRR3", competitor_keys=("M0",))
# CRR-best (information only)
best = None
for rho in np.arange(1.0, 8.01, 0.1):
    ll = sum(np.log(A.crr_pmf(pred[e]["n_near"], pred[e]["n_far"], pred[e]["K"], rho)[obs[e]]) for e in EV)
    if best is None or ll > best[1]: best = (float(rho), float(ll))
out["crr_best"] = dict(rho=best[0], loglik=best[1], loglik_P=float(sum(out["events"][e]["models"]["P"]["logscore"] for e in EV)),
                       loglik_M0=float(sum(out["events"][e]["models"]["M0"]["logscore"] for e in EV)))

# ---------------- seat-level conditional likelihoods
seat = {}
ev2 = A.t2_event(); z2 = np.load(A.RES + "/02_t2_seatprobs.npz")
m0_2 = -np.log(float(S.stats.binom.pmf(0, 1, 0)) * 1.0) if False else None
from scipy.special import comb, logsumexp
def mixll(P, case, K):
    return float(logsumexp([A.seat_loglik(p, case, K) for p in P]) - np.log(len(P)))
seat["T2"] = dict(M0=float(-np.log(comb(45, 7))), KA_nonoise=A.seat_loglik(z2["KA_nonoise"], ev2["case"], 7),
                  P=mixll(z2["KA_noise_draws"], ev2["case"], 7), KB=mixll(z2["KB_draws"], ev2["case"], 7),
                  KA_measured_points=A.seat_loglik(z2["KA_measured_points"], ev2["case"], 7))
z3 = np.load(A.RES + "/02_t3_seatprobs.npz"); c3 = np.isin(z3["seats"], ["1B", "5A"])
seat["T3"] = dict(M0=float(-np.log(comb(len(c3), 2))), KA_nonoise=A.seat_loglik(z3["KA_nonoise"], c3, 2), P=mixll(z3["KA_noise_draws"], c3, 2),
                  n=int(len(c3)), s4_cases=dict(zip(z3["seats"][c3].tolist(), np.array(pred["T3"]["s4"])[c3].tolist())),
                  rank_of_cases=[int((np.array(pred["T3"]["s4"]) > v).sum() + 1) for v in np.array(pred["T3"]["s4"])[c3]])
ev5 = A.t5_event(); z5 = np.load(A.RES + "/02_t5_seatprobs.npz")
seat["T5"] = dict(M0=float(-np.log(comb(20, 12))), KA_5L=A.seat_loglik(z5["KA_5L"], ev5["case"], 12),
                  KA_5A_mirrored=A.seat_loglik(z5["KA_5A-mirrored"], ev5["case"], 12),
                  P=mixll(np.stack([z5["KA_5L"], z5["KA_5A-mirrored"]]), ev5["case"], 12), KB=mixll(z5["KB_draws"], ev5["case"], 12))
for e in seat:
    for k in list(seat[e]):
        if k not in ("M0", "n", "s4_cases", "rank_of_cases"): seat[e]["delta_" + k] = seat[e][k] - seat[e]["M0"]
out["seat_level"] = seat

# ---------------- T4
d4 = E.t4_data(); mask = ~((d4["dr"] == 0) & (d4["dc"] == 1)).to_numpy()
k4 = d4["k"].to_numpy()[mask].astype(float); n4 = d4["n"].to_numpy()[mask].astype(float); K4 = int(k4.sum())
assert K4 == 138
t4 = pred["T4"]
def dev(k, pi):
    e = k.sum() * pi
    return float(2 * np.sum(np.where(k > 0, k * np.log(np.where(k > 0, k, 1) / e), 0.0)))
def ll_mix(pis):
    a = np.log(np.array(pis)) @ k4
    return float(logsumexp(a) - np.log(len(a)))
null_dev = np.load(A.RES + "/03_t4_null_deviance.npy")
rngb = np.random.default_rng(A.SEED + 4)
T4 = {"K": K4, "models": {}}
for name in ("P", "KB_end_release", "M2cal_mode"):
    pis = np.array(t4[name]["pi_draws"]); pibar = pis.mean(axis=0)
    idx = rngb.integers(len(pis), size=2000)
    nd = np.array([dev(rngb.multinomial(K4, pis[i]).astype(float), pibar) for i in idx])
    T4["models"][name] = dict(loglik=ll_mix(pis), deviance=dev(k4, pibar), p_boot=float((nd >= dev(k4, pibar)).mean()),
                              p_chi2_df21=float(S.stats.chi2.sf(dev(k4, pibar), len(k4) - 1)), expected=(K4 * pibar).tolist())
for name in ("M0", "CRR2", "CRR3"):
    pi = np.array(t4[name]["pi"])
    nd = np.array([dev(rngb.multinomial(K4, pi).astype(float), pi) for _ in range(2000)])
    T4["models"][name] = dict(loglik=float(k4 @ np.log(pi)), deviance=dev(k4, pi), p_boot=float((nd >= dev(k4, pi)).mean()),
                              p_chi2_df21=float(S.stats.chi2.sf(dev(k4, pi), len(k4) - 1)), expected=(K4 * pi).tolist())
for m in T4["models"]:
    T4["models"][m]["delta_vs_M0"] = T4["models"][m]["loglik"] - T4["models"]["M0"]["loglik"]
dr = d4["dr"].to_numpy()[mask]; dc = d4["dc"].to_numpy()[mask]
T4["cells"] = dict(dr=dr.tolist(), dc=dc.tolist(), k=k4.tolist(), n=n4.tolist())
T4["row_profile"] = {int(r): dict(k=float(k4[dr == r].sum()), n=float(n4[dr == r].sum()), rate_pct=float(100 * k4[dr == r].sum() / n4[dr == r].sum()),
                                  **{m: float(np.array(T4["models"][m]["expected"])[dr == r].sum()) for m in T4["models"]}) for r in range(4)}
out["events"]["T4"] = T4

# ---------------- O1
xy, case = E.o1_desks(); Aadj = E.o1_adjacency(xy, 1.25)
z_obs, jj, mu, sd = S.joincount_z(Aadj, case)
oz = np.load(A.RES + "/02_o1_z.npz")
O1 = dict(z_obs=float(z_obs), case_case_pairs=float(jj), null_mean=float(mu), null_sd=float(sd), models={})
for m in ("P", "M0", "KB_Dx2", "KB_Dhalf"):
    z = oz[m]; dr_ = oz[m + "_draw"]
    fr = np.array([(z[dr_ == j] <= z_obs).mean() for j in np.unique(dr_)])
    lo, hi = np.quantile(z, [0.025, 0.975])
    O1["models"][m] = dict(P_le_obs_pooled=float((z <= z_obs).mean()), P_le_obs_equal_weight=float(fr.mean()), int95=[float(lo), float(hi)],
                           median=float(np.median(z)), inside95=bool(lo <= z_obs <= hi), n=int(len(z)), n_draws=int(len(fr)),
                           draws_with_P_gt_0025=int((fr > 0.025).sum()))
out["events"]["O1"] = O1

# ---------------- decision
fails = []
for e in ["T2", "T1", "T5", "C1"]:
    if out["events"][e]["models"]["P"]["p_two_sided"] < 0.05: fails.append(e)
r1p = [out["events"]["R1_byK"][K]["models"]["P"]["p_two_sided"] for K in ("2", "3", "4", "5")]
if all(p < 0.05 for p in r1p): fails.append("R1")
if T4["models"]["P"]["p_boot"] < 0.05: fails.append("T4")
pO = O1["models"]["P"]["P_le_obs_pooled"]
if not (0.025 <= pO <= 0.975): fails.append("O1")
po = out["pooled"]["P"]
D2 = bool(po["M0"]["delta"] > 0 and po["M0"]["p_exact_under_competitor"] < 0.05 and T4["models"]["P"]["delta_vs_M0"] > 0)
D3 = bool(all(po[j]["delta"] > 0 and po[j]["p_exact_under_competitor"] < 0.05 for j in ("CRR2", "CRR3")))
if len(fails) >= 2: oc = "S5"
elif len(fails) == 1: oc = "S4"
elif D2 and D3: oc = "S1"
elif D2: oc = "S2"
else: oc = "S3"
out["decision"] = dict(D1_failures=fails, R1_p_by_K=dict(zip(("2", "3", "4", "5"), r1p)), D2=D2, D3=D3, outcome=oc)

# ---------------- same rule applied to the competitors (for comparison; adequacy only)
comp = {}
for m in ("M0", "CRR2", "CRR3"):
    f = [e for e in ["T2", "T1", "T5", "C1"] if out["events"][e]["models"][m]["p_two_sided"] < 0.05]
    if all(out["events"]["R1_byK"][K]["models"][m]["p_two_sided"] < 0.05 for K in ("2", "3", "4", "5")): f.append("R1")
    if T4["models"][m]["p_boot"] < 0.05: f.append("T4")
    if m == "M0" and not (0.025 <= O1["models"]["M0"]["P_le_obs_pooled"] <= 0.975): f.append("O1")
    comp[m] = f
out["competitor_adequacy_failures"] = comp
A.save_json("04_evaluation.json", out)

# ---------------- print
for e in EV + ["T2_side", "T5_no2A"]:
    t = out["events"][e]
    print(f"\n{e}: observed {t['k_obs']}/{t['n_near']} vs {t['K'] - t['k_obs']}/{t['n_far']}  RR_obs {t['rr_obs']:.2f}")
    for m, r in t["models"].items():
        print(f"   {m:26s} RR_exp {(r['rr_expected'] if r['rr_expected'] is not None else float('nan')):.2f}  E[near] {r['mean_near']:.2f}  95% {r['int95']}  p {r['p_two_sided']:.4f}  logscore {r['logscore']:.3f}")
print("\nR1 by K:", {K: {m: round(v["p_two_sided"], 4) for m, v in out["events"]["R1_byK"][K]["models"].items()} for K in out["events"]["R1_byK"]})
print("\npooled:")
for tag, r in out["pooled"].items():
    print("  ", tag, {j: (round(v["delta"], 3), round(v["p_exact_under_competitor"], 5)) for j, v in r.items()})
print("   per event P:", {j: {e: round(x, 2) for e, x in v["per_event"].items()} for j, v in out["pooled"]["P"].items()})
print("crr best", out["crr_best"])
print("\nseat level", json.dumps(seat, indent=1, default=float))
print("\nT4:", {m: {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != "expected"} for m, r in T4["models"].items()})
print("row profile", json.dumps(T4["row_profile"], indent=1))
print("\nO1:", json.dumps(O1, indent=1))
print("\nDECISION", out["decision"], "competitors:", comp)
