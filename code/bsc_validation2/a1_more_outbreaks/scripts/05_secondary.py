"""Secondary analyses (pre-registration 8).  Reported, never used to upgrade the primary outcome."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from scipy import stats as st
from a1 import events as E, analysis as A, engine as G, pipeline as P
R1DIR = os.path.join(ROOT, "..", "..", "bsc_validation", "data")
lg = P.lg
obs = {e: E.build(e, {"n_cfg": 1})[0]["obs"] for e in E.EVENTS}
ALL = {"air": E.CAL["air"] + E.HOLD["air"], "room": E.CAL["room"] + E.HOLD["room"]}
res = {}
slim = lambda r: dict(pooled=r["pooled"], events={e: {m: dict(logp=v[m]["logp"], p_adeq=v[m]["p_adeq"]) for m in A.MODELS}
                                                   | dict(obs=v["obs"], sizes=v["sizes"]) for e, v in r["events"].items()}, post=r["post"])


def crr_pmf(ev, rho, opt=None):
    meta, draws = E.build(ev, dict(opt or {}, n_cfg=1))
    t = A.tab(ev)
    near = np.isin(draws[0]["bins"], meta["near_bins"])
    return G.cond_pmf(np.where(near, rho, 1.0), draws[0]["bins"], meta["sizes"], meta["K"], t["comps"])


# ---- CRR with the published constant 2.4, primary split
fixed = {e: crr_pmf(e, 2.4) for e in E.EVENTS}
res["CRR_2.4_primary_split"] = slim(P.run(E.CAL, E.HOLD, obs, crr_fixed=fixed, targets=("M2", "M3")))

# ---- S-T4: frozen transfer of the round-1 posteriors
t4 = json.load(open(os.path.join(R1DIR, "01_cal_T4_primary.json")))
r1 = json.load(open(os.path.join(R1DIR, "01_cal_R1_primary.json")))
assert np.allclose(t4["logD"], G.LOGD) and np.allclose(r1["logD"], G.LOGD)
def wts(ll):
    ll = np.array(ll)
    w = np.exp(ll - ll.max())
    w = w.sum(axis=1) if w.ndim == 2 else w
    return w / w.sum()
W = {("air", "M2"): wts(t4["M2"]["profile_logL"]), ("air", "M1"): wts(t4["M1"]["profile_logL"]),
     ("room", "M2"): wts(r1["M2"]["profile_logL"]), ("room", "M1"): wts(r1["M1"]["profile_logL"])}
res["S-T4"] = slim(P.run(E.CAL, ALL, obs, weights=W, crr_fixed=fixed, targets=("M2", "M1")))
res["S-T4"]["round1_posterior"] = {"%s|%s" % k: A.summarize(k[1], v) for k, v in W.items()}
res["S-T4_newflights_only"] = slim(P.run(E.CAL, {"air": ALL["air"]}, obs, weights=W, crr_fixed=fixed, targets=("M2", "M1")))
res["S-T4_covid_flights"] = slim(P.run(E.CAL, {"air": E.HOLD["air"]}, obs, weights=W, crr_fixed=fixed, targets=("M2", "M1")))

# ---- S-swap
res["S-swap"] = slim(P.run(E.HOLD, E.CAL, obs))

# ---- S-LOO, per-event preferred D, common-D test
loo = {}
for cls in ALL:
    for e in ALL[cls]:
        r = P.run({cls: [x for x in ALL[cls] if x != e]}, {cls: [e]}, obs, n_mc=50000)
        loo[e] = dict(obs=obs[e], sizes=r["events"][e]["sizes"],
                      **{m: r["events"][e][m] for m in A.MODELS},
                      pooled={t: r["pooled"][t]["all"] for t in ("M2", "M3", "M1")},
                      post={k: v for k, v in r["post"].items()})
res["S-LOO"] = loo
# pooled LOO statistic with reference distributions
for cls in list(ALL) + ["both"]:
    evs = ALL[cls] if cls != "both" else ALL["air"] + ALL["room"]
    pm = {m: [] for m in A.MODELS}
    io = []
    for e in evs:
        c = E.CLS[e]
        for m in A.MODELS:
            w = A.posterior(m, [x for x in ALL[c] if x != e], obs)[0]
            pm[m].append(A.predictive(e, m, w))
        io.append(A.comp_index(A.tab(e), obs[e]))
    out = {}
    for t in ("M2", "M3", "M1"):
        d0 = float(sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[t], pm["M0"], io)))
        dc = float(sum(lg(a[i]) - lg(b[i]) for a, b, i in zip(pm[t], pm["CRR"], io)))
        out[t] = dict(delta0=d0, deltaC=dc,
                      p_B1=float((A.delta_reference(pm[t], pm["M0"], pm["M0"], 200000, 7) >= d0).mean()),
                      p_B2=float((A.delta_reference(pm[t], pm["CRR"], pm["CRR"], 200000, 8) >= dc).mean()),
                      p_mirror=float((A.delta_reference(pm[t], pm["CRR"], pm[t], 200000, 9) <= dc).mean()))
    res["S-LOO_pooled_" + cls] = out

pref = {}
for cls in ALL:
    for model in ("M2", "M1"):
        curves = {e: A.event_curve(e, model, obs[e]) for e in ALL[cls]}
        tot = sum(curves.values())
        lr = 2 * (sum(c.max() for c in curves.values()) - tot.max())
        pref["%s|%s" % (cls, model)] = dict(
            per_event={e: dict(D_hat=float(10 ** G.LOGD[int(np.argmax(c))]), logL_max=float(c.max()),
                               logL_at_wellmixed=float(c[-1]), logL_at_common=float(c[int(np.argmax(tot))]),
                               D_lo=float(10 ** G.LOGD[np.flatnonzero(c >= c.max() - 1.92)[0]]),
                               D_hi=float(10 ** G.LOGD[np.flatnonzero(c >= c.max() - 1.92)[-1]])) for e, c in curves.items()},
            D_common=float(10 ** G.LOGD[int(np.argmax(tot))]), LR=float(lr), df=len(curves) - 1,
            p=float(st.chi2.sf(lr, len(curves) - 1)), curves={e: c.tolist() for e, c in curves.items()})
        # by pathogen
        bp = {}
        for pth in sorted(set(E.PATH[e] for e in ALL[cls])):
            ev = [e for e in ALL[cls] if E.PATH[e] == pth]
            t = sum(curves[e] for e in ev)
            ok = np.flatnonzero(t >= t.max() - 1.92)
            bp[pth] = dict(events=ev, D_hat=float(10 ** G.LOGD[int(np.argmax(t))]), D_lo=float(10 ** G.LOGD[ok[0]]),
                           D_hi=float(10 ** G.LOGD[ok[-1]]), gain_over_wellmixed=float(t.max() - t[-1]))
        pref["%s|%s" % (cls, model)]["by_pathogen"] = bp
res["preferred_D"] = pref

# ---- S-seat: seat-level conditional log score, hold-out seat-map events, primary posterior
cal = json.load(open(os.path.join(ROOT, "out", "02_calibration.json")))
seat = {}
for e in ("F5", "F6", "F8"):
    meta, draws = E.build(e)
    r = {"M0": G.seat_logscore(meta, draws, "M0", None)}
    for model, grid in (("CRR", 10 ** G.LOGRHO), ("M2", 10 ** G.LOGD), ("M1", 10 ** G.LOGD)):
        w = np.array(cal["classes"]["air"][model]["w"])
        v = np.array([G.seat_logscore(meta, draws, model, th) if w[i] > 1e-5 else -np.inf for i, th in enumerate(grid)])
        r[model] = float(np.log(np.sum(w * np.exp(v - v[np.isfinite(v)].max()))) + v[np.isfinite(v)].max())
    w = np.array(cal["classes"]["air"]["M3"]["w"])
    grid = [(10 ** d, ww) for d in G.LOGD3 for ww in G.WGRID]
    v = np.array([G.seat_logscore(meta, draws, "M3", th) if w[i] > 1e-5 else -np.inf for i, th in enumerate(grid)])
    r["M3"] = float(np.log(np.sum(w * np.exp(v - v[np.isfinite(v)].max()))) + v[np.isfinite(v)].max())
    seat[e] = r
seat["pooled"] = {m: float(sum(seat[e][m] for e in ("F5", "F6", "F8"))) for m in A.MODELS}
res["S-seat"] = seat
json.dump(res, open(os.path.join(ROOT, "out", "05_secondary.json"), "w"), indent=1)

def show(name):
    r = res[name]
    for t, v in r["pooled"].items():
        print(name, t, {k: (round(x, 4) if isinstance(x, float) else x) for k, x in v["all"].items() if k in ("delta0", "p_B1", "deltaC", "p_B2", "p_mirror")}, v["decision"])
    print("   ", {e: {m: round(v[m]["p_adeq"], 3) for m in ("M0", "CRR", "M2")} for e, v in r["events"].items()})
for n in ("CRR_2.4_primary_split", "S-T4", "S-T4_newflights_only", "S-T4_covid_flights", "S-swap"):
    show(n)
print("T4 posterior", res["S-T4"]["round1_posterior"])
for e, v in loo.items():
    print("LOO", e, v["obs"], {m: (round(v[m]["logp"], 2), round(v[m]["p_adeq"], 3)) for m in A.MODELS}, "d0 %.2f dC %.2f" % (v["pooled"]["M2"]["delta0"], v["pooled"]["M2"]["deltaC"]))
for k in res:
    if k.startswith("S-LOO_pooled"): print(k, json.dumps(res[k]))
for k, v in pref.items():
    print(k, "common D", v["D_common"], "LR", round(v["LR"], 2), "p", round(v["p"], 4), {e: (x["D_hat"], x["D_lo"], x["D_hi"]) for e, x in v["per_event"].items()})
    print("  by pathogen", v["by_pathogen"])
print("seat", json.dumps(seat))
