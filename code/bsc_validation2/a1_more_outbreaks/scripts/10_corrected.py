"""Post-hoc re-run of the registered pipeline with the round-off-safe M2 kernel (POSTHOC_LOG P5).
Same events, strata, split, grids, statistics and decision rule; only the M2/M3 tables differ."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from scipy import stats as st
from a1 import events as E, analysis as A, engine as G, pipeline as P
obs = {e: E.build(e, {"n_cfg": 1})[0]["obs"] for e in E.EVENTS}
TAGS = {e: "nfix" for e in E.EVENTS}
ALL = {"air": E.CAL["air"] + E.HOLD["air"], "room": E.CAL["room"] + E.HOLD["room"]}
lg = P.lg
res = {}
def keep(r):
    return dict(pooled=r["pooled"], post=r["post"], events={e: dict(obs=v["obs"], sizes=v["sizes"], **{m: v[m] for m in A.MODELS}) for e, v in r["events"].items()})
def show(name, r, targets=("M2", "M3")):
    for t in targets:
        a = r["pooled"][t]["all"]; d = r["pooled"][t]["decision"]
        print("%-28s %s d0 %+.2f (p %.2g) dC %+.2f (pB2 %.3g, mirror %.2g) %s fails %s" % (name, t, a["delta0"], a["p_B1"], a["deltaC"], a["p_B2"], a["p_mirror"], d["outcome"], d["B3_failures"]))
        for sub in ("air", "room"):
            if sub in r["pooled"][t]:
                s = r["pooled"][t][sub]; print("      %s: d0 %+.2f (p %.2g) dC %+.2f (pB2 %.3g, mirror %.3g)" % (sub, s["delta0"], s["p_B1"], s["deltaC"], s["p_B2"], s["p_mirror"]))
    for e, v in r["events"].items():
        print("   ", e, v["obs"], {m: (round(v[m]["logp"], 2), float("%.2g" % v[m]["p_adeq"]), [round(x, 1) for x in v[m]["expected"]]) for m in ("M0", "CRR", "M2", "M3")})
    print("    post", {k: {kk: (float("%.3g" % vv) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in r["post"].items() if k.split("|")[1] in ("M2", "M3")})

r = P.run(E.CAL, E.HOLD, obs, tags=TAGS); res["primary_corrected"] = keep(r); show("primary (corrected kernel)", r)
r = P.run({"air": E.CAL["air"], "room": ["W1"]}, E.HOLD, obs, tags=TAGS); res["S-noW2_corrected"] = keep(r); show("S-noW2 (corrected)", r)
r = P.run(E.HOLD, E.CAL, obs, tags=TAGS); res["S-swap_corrected"] = keep(r); show("S-swap (corrected)", r)
# frozen transfer of round-1 posteriors
R1DIR = os.path.join(ROOT, "..", "..", "bsc_validation", "data")
t4 = json.load(open(os.path.join(R1DIR, "01_cal_T4_primary.json"))); r1 = json.load(open(os.path.join(R1DIR, "01_cal_R1_primary.json")))
def wts(ll):
    ll = np.array(ll); w = np.exp(ll - ll.max()); w = w.sum(axis=1) if w.ndim == 2 else w; return w / w.sum()
W = {("air", "M2"): wts(t4["M2"]["profile_logL"]), ("room", "M2"): wts(r1["M2"]["profile_logL"])}
r = P.run(E.CAL, ALL, obs, tags=TAGS, weights=W, targets=("M2",)); res["S-T4_corrected"] = keep(r); show("S-T4 all 11 (corrected)", r, ("M2",))
r = P.run(E.CAL, {"air": E.HOLD["air"]}, obs, tags=TAGS, weights=W, targets=("M2",)); res["S-T4_covid_flights_corrected"] = keep(r); show("S-T4 covid flights (corr.)", r, ("M2",))
# power of the rule under the corrected calibration (hold-out outcomes not used)
post = {(c, m): A.posterior(m, E.CAL[c], obs, TAGS)[0] for c in E.CAL for m in A.MODELS}
evs = [(e, c) for c in E.HOLD for e in E.HOLD[c]]
pm = {m: [A.predictive(e, m, post[c, m], "nfix") for e, c in evs] for m in A.MODELS}
n = 200000; power = {}
for target in ("M2", "M3"):
    crit0 = np.quantile(A.delta_reference(pm[target], pm["M0"], pm["M0"], n, 11), 0.95)
    critC = np.quantile(A.delta_reference(pm[target], pm["CRR"], pm["CRR"], n, 12), 0.95)
    critM = np.quantile(A.delta_reference(pm[target], pm["CRR"], pm[target], n, 13), 0.05)
    row = {}
    for truth in ("M0", "CRR", target):
        idx = A.sample_idx(pm[truth], n, 21)
        d0 = sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[target], pm["M0"], idx))
        dc = sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[target], pm["CRR"], idx))
        nf = sum((np.array([A.two_sided_p(p, j) for j in range(len(p))])[i] < 0.05).astype(int) for p, i in zip(pm[target], idx))
        B1 = (d0 > max(0, crit0)); B2 = (dc > max(0, critC)); mir = (dc < min(0, critM))
        V5 = nf >= 2; V4 = ~V5 & ~B1; V3 = ~V5 & B1 & mir; V1 = ~V5 & B1 & ~mir & B2 & (nf == 0); V2 = ~V5 & B1 & ~mir & ~V1
        row[truth] = {k: float(v.mean()) for k, v in dict(B1=B1, B2=B2, B1B2=B1 & B2, mirror=mir, V1=V1, V2=V2, V3=V3, V4=V4, V5=V5).items()}
    power[target] = row
    print("power", target, json.dumps({k: {a: round(b, 3) for a, b in v.items()} for k, v in row.items()}))
res["power_corrected"] = power
# leave-one-out and per-event D
loo = {}
for cls in ALL:
    for e in ALL[cls]:
        rr = P.run({cls: [x for x in ALL[cls] if x != e]}, {cls: [e]}, obs, tags=TAGS, n_mc=50000, targets=("M2",))
        loo[e] = dict(obs=obs[e], **{m: rr["events"][e][m] for m in A.MODELS}, post=rr["post"])
res["S-LOO_corrected"] = loo
for cls in list(ALL) + ["both"]:
    ev = ALL[cls] if cls != "both" else ALL["air"] + ALL["room"]
    pmm = {m: [] for m in A.MODELS}; io = []
    for e in ev:
        c = E.CLS[e]
        for m in A.MODELS:
            pmm[m].append(A.predictive(e, m, A.posterior(m, [x for x in ALL[c] if x != e], obs, TAGS)[0], "nfix"))
        io.append(A.comp_index(A.tab(e, "nfix"), obs[e]))
    o = {}
    for t in ("M2", "M3"):
        d0 = float(sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pmm[t], pmm["M0"], io))); dc = float(sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pmm[t], pmm["CRR"], io)))
        o[t] = dict(delta0=d0, deltaC=dc, p_B1=float((A.delta_reference(pmm[t], pmm["M0"], pmm["M0"], 200000, 7) >= d0).mean()),
                    p_B2=float((A.delta_reference(pmm[t], pmm["CRR"], pmm["CRR"], 200000, 8) >= dc).mean()),
                    p_mirror=float((A.delta_reference(pmm[t], pmm["CRR"], pmm[t], 200000, 9) <= dc).mean()))
    res["S-LOO_pooled_%s_corrected" % cls] = o; print("LOO pooled", cls, json.dumps({t: {k: round(v, 4) for k, v in x.items()} for t, x in o.items()}))
for e, v in loo.items():
    print("LOO", e, v["obs"], {m: (round(v[m]["logp"], 2), float("%.2g" % v[m]["p_adeq"])) for m in ("M0", "CRR", "M2", "M3")})
pref = {}
for cls in ALL:
    curves = {e: A.event_curve(e, "M2", obs[e], "nfix") for e in ALL[cls]}
    tot = sum(curves.values()); lr = 2 * (sum(c.max() for c in curves.values()) - tot.max())
    def ci(c):
        ok = np.flatnonzero(c >= c.max() - 1.92); return float(10 ** G.LOGD[ok[0]]), float(10 ** G.LOGD[ok[-1]])
    pref[cls] = dict(per_event={e: dict(D_hat=float(10 ** G.LOGD[int(np.argmax(c))]), D_lo=ci(c)[0], D_hi=ci(c)[1], gain_over_wellmixed=float(c.max() - c[-1]),
                                         logL_at_common=float(c[int(np.argmax(tot))]), logL_max=float(c.max())) for e, c in curves.items()},
                     D_common=float(10 ** G.LOGD[int(np.argmax(tot))]), D_common_ci=ci(tot), LR=float(lr), df=len(curves) - 1, p=float(st.chi2.sf(lr, len(curves) - 1)),
                     curves={e: c.tolist() for e, c in curves.items()})
    bp = {}
    for pth in sorted(set(E.PATH[e] for e in ALL[cls])):
        ev = [e for e in ALL[cls] if E.PATH[e] == pth]; t = sum(curves[e] for e in ev)
        bp[pth] = dict(events=ev, D_hat=float(10 ** G.LOGD[int(np.argmax(t))]), D_lo=ci(t)[0], D_hi=ci(t)[1], gain_over_wellmixed=float(t.max() - t[-1]))
    pref[cls]["by_pathogen"] = bp
    print(cls, "common D %.3g" % pref[cls]["D_common"], pref[cls]["D_common_ci"], "LR %.2f df %d p %.4g" % (lr, len(curves) - 1, pref[cls]["p"]))
    for e, x in pref[cls]["per_event"].items(): print("    ", e, {k: float("%.3g" % v) for k, v in x.items()})
    print("     by pathogen", {k: (float("%.3g" % v["D_hat"]), float("%.3g" % v["D_lo"]), float("%.3g" % v["D_hi"]), round(v["gain_over_wellmixed"], 2)) for k, v in bp.items()})
res["preferred_D_corrected"] = pref
json.dump(res, open(os.path.join(ROOT, "out", "10_corrected.json"), "w"), indent=1)
