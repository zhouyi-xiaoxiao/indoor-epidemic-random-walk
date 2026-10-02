"""Step 3: power of the decision rule, before unblinding.  Reads predictions only (no outcome).
Writes results/03_power.json and PREREG_ADDENDUM_power.md."""
import itertools
import numpy as np
import a3lib as A

pred = A.load_json("02_predictions.json")["events"]
oz = np.load(A.RES + "/02_o1_z.npz")
rng = np.random.default_rng(A.SEED + 3)
EV = ["T2", "T1", "T5", "C1", "R1"]
TRUTHS = ["P", "M0", "CRR2", "CRR3"]
pm = {e: {m: np.array(pred[e]["pmf"][m]) for m in TRUTHS} for e in EV}
supp = {e: np.flatnonzero(sum(pm[e][m] for m in TRUTHS) > 0) for e in EV}
tsp = {e: {m: np.array([A.two_sided_p(pm[e][m], k) for k in range(len(pm[e][m]))]) for m in TRUTHS} for e in EV}
lg = {e: {m: np.log(np.clip(pm[e][m], 1e-300, None)) for m in TRUTHS} for e in EV}

# ---- per-event discrimination
disc = {}
for e in EV:
    disc[e] = {f"P(reject {j} | {t} true)": float(pm[e][t][tsp[e][j] < 0.05].sum()) for t in TRUTHS for j in TRUTHS}

# ---- joint enumeration over the five stratum events
grids = [supp[e] for e in EV]
combos = np.array(list(itertools.product(*grids)))              # (N, 5)
def col(fn):
    return np.stack([fn(e)[combos[:, i]] for i, e in enumerate(EV)], axis=1)
logw = {t: col(lambda e: lg[e][t]).sum(axis=1) for t in TRUTHS}
w = {t: np.exp(logw[t]) for t in TRUTHS}
for t in TRUTHS: assert abs(w[t].sum() - 1) < 1e-9, (t, w[t].sum())
nfailP = (col(lambda e: tsp[e]["P"]) < 0.05).sum(axis=1)
delta = {j: logw["P"] - logw[j] for j in ("M0", "CRR2", "CRR3")}
def crit(j):   # smallest c with P_j(delta >= c) < 0.05  -> reject j when delta > c95
    o = np.argsort(delta[j]); cw = np.cumsum(w[j][o])
    return float(delta[j][o][np.searchsorted(cw, 0.95)])
c95 = {j: crit(j) for j in delta}
passD2s = (delta["M0"] > 0) & (delta["M0"] > c95["M0"])
passD3 = (delta["CRR2"] > 0) & (delta["CRR2"] > c95["CRR2"]) & (delta["CRR3"] > 0) & (delta["CRR3"] > c95["CRR3"])

# ---- T4: adequacy (deviance, parametric bootstrap) and sign of the log-score difference
t4 = pred["T4"]; K4 = 138
piP = np.array(t4["P"]["pi_draws"]); pibar = piP.mean(axis=0)
pi_truth = {"P": piP, "M0": np.array(t4["M0"]["pi"])[None, :], "CRR2": np.array(t4["CRR2"]["pi"])[None, :], "CRR3": np.array(t4["CRR3"]["pi"])[None, :]}
def dev(k, pi):
    e = k.sum(axis=-1, keepdims=True) * pi
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(k > 0, k * np.log(k / e), 0.0)
    return 2 * t.sum(axis=-1)
def loglikP(k):       # log mean over draws of prod pi^k
    a = k @ np.log(piP).T                                         # (R, draws)
    mx = a.max(axis=1, keepdims=True)
    return (mx + np.log(np.exp(a - mx).mean(axis=1, keepdims=True)))[:, 0]
def sim(pis, R):
    idx = rng.integers(len(pis), size=R)
    return np.array([rng.multinomial(K4, pis[i]) for i in idx])
null_dev = dev(sim(piP, 4000), pibar)
t4res = {}
for t in TRUTHS:
    k = sim(pi_truth[t], 2000)
    pval = (null_dev[None, :] >= dev(k, pibar)[:, None]).mean(axis=1)
    okA = pval >= 0.05
    dl = loglikP(k) - k @ np.log(np.array(t4["M0"]["pi"]))
    t4res[t] = dict(P_adequate=float(okA.mean()), P_delta_pos=float((dl > 0).mean()),
                    joint={"adequate&pos": float((okA & (dl > 0)).mean()), "adequate&neg": float((okA & (dl <= 0)).mean()),
                           "fail&pos": float((~okA & (dl > 0)).mean()), "fail&neg": float((~okA & (dl <= 0)).mean())})
np.save(A.RES + "/03_t4_null_deviance.npy", null_dev)

# ---- O1
zP, zM0 = oz["P"], oz["M0"]
lo, hi = np.quantile(zP, [0.025, 0.975])
o1 = {"P": 0.95, "M0": float(((zM0 >= lo) & (zM0 <= hi)).mean())}
o1["CRR2"] = o1["CRR3"] = o1["M0"]

# ---- outcome probabilities
res = {}
for t in TRUTHS:
    pf = np.array([w[t][nfailP == i].sum() for i in range(6)])          # failures among the 5 stratum events
    out = dict(S1=0.0, S2=0.0, S3=0.0, S4=0.0, S5=0.0)
    for a4, lab in ((True, "adequate"), (False, "fail")):
        for pos in (True, False):
            p4 = t4res[t]["joint"][f"{lab}&{'pos' if pos else 'neg'}"]
            for aO in (True, False):
                pO = o1[t] if aO else 1 - o1[t]
                extra = (0 if a4 else 1) + (0 if aO else 1)
                for nf in range(6):
                    sel = nfailP == nf
                    tot = nf + extra
                    if tot >= 2:
                        out["S5"] += p4 * pO * w[t][sel].sum()
                    elif tot == 1:
                        out["S4"] += p4 * pO * w[t][sel].sum()
                    else:
                        d2 = passD2s[sel] & pos
                        out["S1"] += p4 * pO * w[t][sel][d2 & passD3[sel]].sum()
                        out["S2"] += p4 * pO * w[t][sel][d2 & ~passD3[sel]].sum()
                        out["S3"] += p4 * pO * w[t][sel][~d2].sum()
    res[t] = dict(outcomes=out, P_no_stratum_failure=float(pf[0]), P_D2_stratum=float(w[t][passD2s].sum()),
                  P_D3=float(w[t][passD3].sum()), T4=t4res[t], O1_adequate=o1[t])
out = dict(stamp=A.stamp(), c95=c95, per_event_discrimination=disc, truths=res, o1_interval_P=[float(lo), float(hi)],
           n_joint_outcomes=int(len(combos)))
A.save_json("03_power.json", out)

L_ = ["# Addendum: power of the decision rule (computed before unblinding)", "",
      f"Written {A.stamp()} by `src/s03_power.py` from `results/02_predictions.json` (kernels from physics data,",
      "strata sizes and event totals). No stratum outcome was read. Exact enumeration over the joint distribution",
      f"of the five stratum events ({len(combos)} outcomes); T4 by 2,000 simulated tables per truth; O1 from the simulated",
      "join-count distributions.", "",
      "## Probability of each outcome of the decision rule", "",
      "| Truth | S1 | S2 | S3 | S4 (one failure) | S5 (two or more) | no failure among 5 stratum events | D2 (stratum part) | D3 | T4 adequate | T4 score above M0 | O1 adequate |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for t in TRUTHS:
    r = res[t]; o = r["outcomes"]
    L_.append(f"| {t} | {o['S1']:.3f} | {o['S2']:.3f} | {o['S3']:.3f} | {o['S4']:.3f} | {o['S5']:.3f} | {r['P_no_stratum_failure']:.3f} | {r['P_D2_stratum']:.3f} | {r['P_D3']:.3f} | {r['T4']['P_adequate']:.3f} | {r['T4']['P_delta_pos']:.3f} | {r['O1_adequate']:.3f} |")
L_ += ["", "Reading: row P is the probability that the physics kernel earns each outcome if it is true; rows M0, CRR2,",
       "CRR3 are the probabilities that it earns them when the truth is well mixed or a constant near/far risk ratio",
       "of 2 or 3 (for these truths O1 is generated under M0 and T4 with near = same row).", "",
       f"Critical values of the pooled log-score difference (95th percentile under the competitor): {', '.join(f'{k}: {v:+.2f}' for k, v in c95.items())}.", "",
       "## Which events can discriminate", "",
       "Probability that the model in the column is rejected by the two-sided 5% test when the model in the row is true.", "",
       "| Event | Truth | reject P | reject M0 | reject CRR2 | reject CRR3 |", "|---|---|---|---|---|---|"]
for e in EV:
    for t in TRUTHS:
        d = disc[e]
        L_.append(f"| {e} | {t} | " + " | ".join(f"{d[f'P(reject {j} | {t} true)']:.3f}" for j in TRUTHS) + " |")
open(A.ROOT + "/PREREG_ADDENDUM_power.md", "w").write("\n".join(L_) + "\n")
print("\n".join(L_))
