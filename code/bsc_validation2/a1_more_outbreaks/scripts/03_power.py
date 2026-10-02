"""Power of the decision rule (pre-registration 0.3, 7).  Reads the calibration
posteriors and the hold-out predictive tables (which depend on hold-out geometry and
event totals only).  Does NOT read any hold-out stratum count."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from a1 import events as E, analysis as A
cal = json.load(open(os.path.join(ROOT, "out", "02_calibration.json")))
N = 200000
hold = [(e, "air") for e in E.HOLD["air"]] + [(e, "room") for e in E.HOLD["room"]]
pm = {m: [A.predictive(e, m, np.array(cal["classes"][c][m]["w"])) for e, c in hold] for m in A.MODELS}
lg = lambda a: np.log(np.clip(a, 1e-300, None))


def adequacy_all(model_pm, idx):
    """per event: two-sided p of sampled outcomes under model_pm."""
    ps = []
    for a, i in zip(model_pm, idx):
        order = np.argsort(a)
        cum = np.cumsum(a[order])
        rank = np.empty(len(a), int); rank[order] = np.arange(len(a))
        # p = sum of pmf <= pmf[i] (ties included)
        srt = a[order]
        hi = np.searchsorted(srt, a[i] * (1 + 1e-9) + 1e-15, side="right") - 1
        ps.append(cum[hi])
    return np.array(ps)


res = dict(n_mc=N, hold=[e for e, _ in hold], rules={})
for target in ("M2", "M3", "M1"):
    idx = {t: A.sample_idx(pm[t], N, seed=11 + k) for k, t in enumerate(["M0", "CRR", target])}
    stat = {}
    for t in idx:
        d0 = sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[target], pm["M0"], idx[t]))
        dc = sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[target], pm["CRR"], idx[t]))
        ad = adequacy_all(pm[target], idx[t])
        stat[t] = dict(d0=d0, dc=dc, nfail=(ad < 0.05).sum(axis=0))
    c0 = float(np.quantile(stat["M0"]["d0"], 0.95))           # B1 critical value
    cC = float(np.quantile(stat["CRR"]["dc"], 0.95))          # B2 critical value
    cM = float(np.quantile(stat[target]["dc"], 0.05))         # mirror critical value
    r = dict(crit_B1=c0, crit_B2=cC, crit_mirror=cM, truth={})
    for t in idx:
        s = stat[t]
        B1 = (s["d0"] > max(c0, 0)); B2 = (s["dc"] > max(cC, 0)); mir = (s["dc"] < min(cM, 0))
        B3 = s["nfail"] == 0; B3a = s["nfail"] <= 1
        r["truth"][t] = dict(P_B1=float(B1.mean()), P_B2=float(B2.mean()), P_mirror=float(mir.mean()),
                             P_B3_all=float(B3.mean()), P_B3_atmost1=float(B3a.mean()),
                             P_V1=float((B1 & B2 & B3).mean()), P_B1_B2=float((B1 & B2).mean()),
                             P_V2=float((B1 & B3a & ~B2 & ~mir & (s["nfail"] < 2)).mean()),
                             P_V3=float((B1 & mir & (s["nfail"] < 2)).mean()), P_V4=float((~B1 & (s["nfail"] < 2)).mean()),
                             P_V5=float((s["nfail"] >= 2).mean()),
                             mean_d0=float(s["d0"].mean()), mean_dc=float(s["dc"].mean()))
    # class-specific power of B1 & B2 (aircraft events only)
    k = len(E.HOLD["air"])
    sub = {}
    for t in idx:
        d0 = sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[target][:k], pm["M0"][:k], idx[t][:k]))
        dc = sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[target][:k], pm["CRR"][:k], idx[t][:k]))
        sub[t] = (d0, dc)
    a0 = float(np.quantile(sub["M0"][0], 0.95)); aC = float(np.quantile(sub["CRR"][1], 0.95))
    r["air_only"] = dict(crit_B1=a0, crit_B2=aC,
                         **{t: dict(P_B1=float((sub[t][0] > max(a0, 0)).mean()), P_B2=float((sub[t][1] > max(aC, 0)).mean()),
                                    P_B1_B2=float(((sub[t][0] > max(a0, 0)) & (sub[t][1] > max(aC, 0))).mean())) for t in idx})
    res["rules"][target] = r
    print(target, json.dumps(r, indent=1))
# expected bin counts under each model (shape of predictions, no outcome)
res["expected"] = {e: {m: A.expected_counts(A.tab(e), p).tolist() for m, p in ((m, pm[m][k]) for m in A.MODELS)}
                   for k, (e, _) in enumerate(hold)}
res["sizes"] = {e: A.tab(e)["sizes"].tolist() for e, _ in hold}
res["K"] = {e: int(A.tab(e)["K"]) for e, _ in hold}
json.dump(res, open(os.path.join(ROOT, "out", "03_power.json"), "w"), indent=1)
print(json.dumps(res["expected"], indent=0))
