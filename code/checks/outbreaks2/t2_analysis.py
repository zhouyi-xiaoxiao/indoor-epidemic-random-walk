"""Recomputation, with separately written code, of the A1 primary analysis, power, and additional checks."""
import numpy as np, json, os, sys
import vlib as V

EV = ["F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "W1", "W2", "P1"]
CLS = {e: ("air" if e[0] == "F" else "room") for e in EV}
CAL = {"air": ["F1", "F2", "F3", "F4"], "room": ["W1", "W2"]}
HOLD = {"air": ["F5", "F6", "F7", "F8"], "room": ["P1"]}
lg = lambda a: np.log(np.clip(a, 1e-300, None))
NMC = 400000
out = {}

G, T = {}, {}
for e in EV + ["W1x"]:
    g = V.build(e)
    pmf, comps = V.m2_table(g, cache="cache/m2_%s.npz" % e)
    G[e] = g
    T[e] = dict(comps=comps, M2=pmf, M0=V.m0_pmf(g, comps)[None, :], CRR=V.crr_table(g, comps),
                CRRtrue=V.crr_table(g, comps, true=True), i=V.idx_of(comps, g["obs"]))


def post(model, events, tabs=T, mask=None):
    ll = sum(lg(tabs[e][model][:, tabs[e]["i"]]) for e in events)
    if mask is not None:
        ll = np.where(mask, ll, -np.inf)
    w = np.exp(ll - ll.max())
    return w / w.sum(), ll


def summ(w, x):
    c = np.cumsum(w)
    return dict(mode=float(10 ** x[np.argmax(w)]), q05=float(10 ** np.interp(.05, c, x)), q50=float(10 ** np.interp(.5, c, x)),
                q95=float(10 ** np.interp(.95, c, x)))


def refdist(a, b, truth, seed):
    rng = np.random.default_rng(seed)
    tot = np.zeros(NMC)
    for x, y, t in zip(a, b, truth):
        idx = rng.choice(len(t), size=NMC, p=t / t.sum())
        tot += lg(x[idx]) - lg(y[idx])
    return tot


def run(cal, hold, target="M2", comp="CRR", tabs=T, weights=None, seed=5, mask=None, label=""):
    w = {}
    for c in cal:
        for m in (target, comp):
            w[c, m] = weights[c, m] if weights and (c, m) in weights else post(m, cal[c], tabs, mask if m == target else None)[0]
    evs = [(e, c) for c in hold for e in hold[c]]
    pm = {m: [w[c, m] @ tabs[e][m] for e, c in evs] for m in (target, comp)}
    pm["M0"] = [tabs[e]["M0"][0] for e, c in evs]
    io = [tabs[e]["i"] for e, c in evs]
    r = dict(events={}, post={"%s|%s" % (c, m): summ(w[c, m], V.LOGD if m == "M2" else V.LOGRHO) for (c, m) in w if m in ("M2", "CRR", "CRRtrue")})
    for k, (e, c) in enumerate(evs):
        r["events"][e] = {m: dict(logp=float(lg(pm[m][k][io[k]])), p_adeq=V.two_sided(pm[m][k], io[k]),
                                  exp=np.round(pm[m][k] @ tabs[e]["comps"], 2).tolist()) for m in pm}
        r["events"][e]["obs"] = G[e]["obs"] if e in G else None
    for name, ks in (("all", range(len(evs))), ("air", [k for k, (e, c) in enumerate(evs) if c == "air"]),
                     ("room", [k for k, (e, c) in enumerate(evs) if c == "room"])):
        ks = list(ks)
        if not ks:
            continue
        a = [pm[target][k] for k in ks]; b0 = [pm["M0"][k] for k in ks]; bc = [pm[comp][k] for k in ks]
        d0 = float(sum(lg(a[j][io[k]]) - lg(b0[j][io[k]]) for j, k in enumerate(ks)))
        dc = float(sum(lg(a[j][io[k]]) - lg(bc[j][io[k]]) for j, k in enumerate(ks)))
        r0 = refdist(a, b0, b0, seed); rc = refdist(a, bc, bc, seed + 1); rm = refdist(a, bc, a, seed + 2)
        r[name] = dict(delta0=d0, p_B1=float((r0 >= d0 - 1e-12).mean()), deltaC=dc, p_B2=float((rc >= dc - 1e-12).mean()),
                       p_mirror=float((rm <= dc + 1e-12).mean()), E_dC_if_target=float(rm.mean()), E_dC_if_comp=float(rc.mean()))
    fails = [e for e, c in evs if r["events"][e][target]["p_adeq"] < 0.05]
    al = r["all"]
    B1 = al["delta0"] > 0 and al["p_B1"] < .05; B2 = al["deltaC"] > 0 and al["p_B2"] < .05; mir = al["deltaC"] < 0 and al["p_mirror"] < .05
    V_ = "V5" if len(fails) >= 2 else "V4" if not B1 else "V3" if mir else "V1" if (B2 and not fails) else "V2"
    r["decision"] = dict(B1=bool(B1), B2=bool(B2), mirror=bool(mir), B3_failures=fails, outcome=V_)
    r["_pm"] = pm; r["_evs"] = evs
    if label:
        print("%-34s d0 %+6.2f (p %.2g) dC %+6.2f (pB2 %.3g, mirror %.3g) | air d0 %+.2f dC %+.2f (pB2 %.3g mir %.3g) | %s fails %s"
              % (label, al["delta0"], al["p_B1"], al["deltaC"], al["p_B2"], al["p_mirror"],
                 r.get("air", {}).get("delta0", np.nan), r.get("air", {}).get("deltaC", np.nan), r.get("air", {}).get("p_B2", np.nan),
                 r.get("air", {}).get("p_mirror", np.nan), V_, fails), flush=True)
    return r


def strip(r):
    return {k: v for k, v in r.items() if not k.startswith("_")}


# ------------------------------------------------------------------ 1. primary
prim = run(CAL, HOLD, label="PRIMARY (re-check kernel)")
out["primary"] = strip(prim)
for e in prim["events"]:
    print("   ", e, G[e]["sizes"], G[e]["obs"], {m: (round(v["logp"], 2), round(v["p_adeq"], 4), v["exp"]) for m, v in prim["events"][e].items() if m != "obs"})
print("    post", prim["post"])
out["calib_maxloglik"] = {c: {m: float(post(m, CAL[c])[1].max()) for m in ("M0", "CRR", "M2")} for c in CAL}
print("    max loglik", out["calib_maxloglik"])

# ------------------------------------------------------------------ 2. power
def power(pm, evs, target="M2", comp="CRR", seed=50):
    res = {}
    samp = {}
    for kt, t in enumerate(["M0", comp, target]):
        rng = np.random.default_rng(seed + kt)
        samp[t] = [rng.choice(len(p), size=NMC, p=p / p.sum()) for p in pm[t]]
    def stats(idx, ks):
        d0 = sum(lg(pm[target][k][idx[k]]) - lg(pm["M0"][k][idx[k]]) for k in ks)
        dc = sum(lg(pm[target][k][idx[k]]) - lg(pm[comp][k][idx[k]]) for k in ks)
        return d0, dc
    def adeq(idx):
        nf = np.zeros(NMC, int)
        for k, p in enumerate(pm[target]):
            order = np.argsort(p); cum = np.cumsum(p[order])
            pv = cum[np.searchsorted(p[order], p[idx[k]] * (1 + 1e-9) + 1e-15, side="right") - 1]
            nf += pv < 0.05
        return nf
    for name, ks in (("all", list(range(len(evs)))), ("air", [k for k, (e, c) in enumerate(evs) if c == "air"])):
        st = {t: stats(samp[t], ks) for t in samp}
        c0 = np.quantile(st["M0"][0], .95); cC = np.quantile(st[comp][1], .95); cM = np.quantile(st[target][1], .05)
        res[name] = dict(crit_B1=float(c0), crit_B2=float(cC), crit_mirror=float(cM))
        for t in samp:
            d0, dc = st[t]
            B1 = d0 > max(c0, 0); B2 = dc > max(cC, 0); mir = dc < min(cM, 0)
            rr = dict(P_B1=float(B1.mean()), P_B2=float(B2.mean()), P_B1B2=float((B1 & B2).mean()), P_mirror=float(mir.mean()))
            if name == "all":
                nf = adeq(samp[t])
                rr.update(P_B3_all=float((nf == 0).mean()), P_V1=float((B1 & B2 & (nf == 0)).mean()), P_V5=float((nf >= 2).mean()),
                          P_V3=float((B1 & mir & (nf < 2)).mean()), P_V4=float((~B1 & (nf < 2)).mean()))
            res[name][t] = rr
    return res

pw = power(prim["_pm"], prim["_evs"])
out["power"] = pw
print("POWER all:", json.dumps(pw["all"]))
print("POWER air:", json.dumps(pw["air"]))

# ------------------------------------------------------------------ 3. alternative generic competitors
# (a) exact constant risk ratio instead of the first analysis's cloglog form
r = run(CAL, HOLD, comp="CRRtrue", label="vs exact CRR (p_near = rho q)")
out["vs_trueCRR"] = strip(r)
# (b) CRR with other near zones (aircraft: <=1 row, <=5 rows ; P1: <=4 m, <=12 m)
for nm, nb_air, nb_p1 in (("near<=1row / P1<=4m", [0], [0]), ("near<=5rows / P1<=12m", [0, 1, 2], [0, 1, 2])):
    T2 = {}
    for e in EV:
        T2[e] = dict(T[e])
        if e in ("F1", "F3", "F4", "F5", "F6", "F8"):
            T2[e]["CRR"] = V.crr_table(G[e], T[e]["comps"], near_bins=nb_air)
        if e == "F7":
            T2[e]["CRR"] = V.crr_table(G[e], T[e]["comps"], near_bins=[0] if nb_air == [0] else [0, 1])
        if e == "P1":
            T2[e]["CRR"] = V.crr_table(G[e], T[e]["comps"], near_bins=nb_p1)
    r = run(CAL, HOLD, tabs=T2, label="CRR " + nm)
    out["CRR_" + nm] = strip(r)
# (c) fixed published rho = 2.4
wfix = {(c, "CRR"): (np.arange(len(V.LOGRHO)) == np.argmin(np.abs(10 ** V.LOGRHO - 2.4))).astype(float) for c in CAL}
r = run(CAL, HOLD, weights=wfix, label="CRR fixed rho=2.4 (grid 2.40)")
out["CRR_fixed2.4"] = strip(r)

# (d) generic static distance kernel  X_j = sum_src exp(-r/ell)   (no lattice, no airborne physics)
LOGELL = np.round(np.linspace(-0.5, 2.0, 51), 6)
def exp_table(g):
    nx = g["nx"]; comps = compositions = None
    comps = T[g["ev"]]["comps"]
    tab = np.zeros((len(LOGELL), len(comps)))
    for (iS, ik, rec) in g["cfgs"]:
        src = np.array(g["srcsets"][iS], float) * [g["ax"], g["ay"]]
        xy = np.c_[(rec % nx) * g["ax"], (rec // nx) * g["ay"]]
        r_ = np.sqrt(((xy[:, None, :] - src[None, :, :]) ** 2).sum(-1))
        for i, le in enumerate(LOGELL):
            tab[i] += V.cond_pmf(np.exp(-r_ / 10 ** le).sum(1), g["bins"], g["sizes"], g["K"], comps)
    return tab / len(g["cfgs"])
def row_table(g):
    """exponential in ROW distance only (aircraft) / same Euclid for rooms"""
    nx = g["nx"]; comps = T[g["ev"]]["comps"]
    tab = np.zeros((len(LOGELL), len(comps)))
    for (iS, ik, rec) in g["cfgs"]:
        sy = np.array([y for x, y in g["srcsets"][iS]], float) * g["ay"]
        ry = (rec // nx) * g["ay"]
        d = np.abs(ry[:, None] - sy[None, :])
        for i, le in enumerate(LOGELL):
            tab[i] += V.cond_pmf(np.exp(-d / 10 ** le).sum(1), g["bins"], g["sizes"], g["K"], comps)
    return tab / len(g["cfgs"])
for e in EV:
    T[e]["EXP"] = exp_table(G[e])
    T[e]["ROW"] = row_table(G[e]) if e[0] == "F" else T[e]["EXP"]
for gen in ("EXP", "ROW"):
    r = run(CAL, HOLD, target=gen, label="generic %s kernel as TARGET vs CRR" % gen)
    out["generic_%s_target" % gen] = strip(r)
    we, _ = post(gen, CAL["air"]); print("    %s air posterior ell (m):" % gen, summ(we, LOGELL))
    r = run(CAL, HOLD, comp=gen, label="M2 vs generic %s kernel" % gen)
    out["M2_vs_%s" % gen] = strip(r)

# ------------------------------------------------------------------ 4. secondary: swap, LOO, per-event D, heterogeneity
r = run({"air": HOLD["air"], "room": HOLD["room"]}, {"air": CAL["air"], "room": CAL["room"]}, label="S-swap")
out["swap"] = strip(r)
for e in r["events"]:
    print("   ", e, {m: (round(v["logp"], 2), round(v["p_adeq"], 4)) for m, v in r["events"][e].items() if m != "obs"})
air = [e for e in EV if e[0] == "F"]
d0 = dc = 0; loo = {}
for e in air:
    rr = run({"air": [x for x in air if x != e]}, {"air": [e]}, seed=9)
    loo[e] = dict(d0=rr["all"]["delta0"], dC=rr["all"]["deltaC"], p_adeq=rr["events"][e]["M2"]["p_adeq"], p_adeq_CRR=rr["events"][e]["CRR"]["p_adeq"])
    d0 += rr["all"]["delta0"]; dc += rr["all"]["deltaC"]
out["LOO_air"] = dict(per_event=loo, delta0=d0, deltaC=dc)
print("LOO air: d0 %+.2f dC %+.2f" % (d0, dc), {e: (round(v["dC"], 2), round(v["p_adeq"], 3)) for e, v in loo.items()})
from scipy.stats import chi2
per = {}
for e in EV:
    ll = lg(T[e]["M2"][:, T[e]["i"]])
    ok = np.flatnonzero(ll >= ll.max() - 1.92)
    per[e] = dict(Dhat=float(10 ** V.LOGD[ll.argmax()]), lo=float(10 ** V.LOGD[ok.min()]), hi=float(10 ** V.LOGD[ok.max()]),
                  gain_over_M0=float(ll.max() - lg(T[e]["M0"][0, T[e]["i"]])),
                  crr_gain_over_M0=float(lg(T[e]["CRR"][:, T[e]["i"]]).max() - lg(T[e]["M0"][0, T[e]["i"]])))
out["per_event_D"] = per
print("per-event D:", {e: (round(v["Dhat"], 3), round(v["lo"], 3), round(v["hi"], 1), round(v["gain_over_M0"], 2), round(v["crr_gain_over_M0"], 2)) for e, v in per.items()})
het = {}
for name, evs in (("air", air), ("room", ["W1", "W2", "P1"]), ("air_covid", ["F5", "F6", "F7", "F8"])):
    lls = [lg(T[e]["M2"][:, T[e]["i"]]) for e in evs]
    LR = 2 * (sum(l.max() for l in lls) - sum(lls).max())
    het[name] = dict(LR=float(LR), df=len(evs) - 1, p=float(chi2.sf(LR, len(evs) - 1)), D_common=float(10 ** V.LOGD[sum(lls).argmax()]))
out["heterogeneity"] = het
print("heterogeneity:", het)

# ------------------------------------------------------------------ 5. S-T4 frozen transfer (round-1 train posterior)
try:
    c = json.load(open("../../bsc_validation/data/01_cal_T4_primary.json"))
    print("round-1 T4 keys:", list(c.keys())[:12])
    out["T4_keys"] = list(c.keys())[:20]
except Exception as ex:
    print("T4 load failed", ex)

# ------------------------------------------------------------------ 6. forking paths
# (a) restrict D prior range; (b) households collapsed; (c) two-bin strata; (d) exact W1 positions
for nm, mask in (("D prior 1..1000", (V.LOGD >= 0) & (V.LOGD <= 3)), ("D prior 10..10000", V.LOGD >= 1)):
    r = run(CAL, HOLD, mask=mask, label="M2 " + nm); out["prior_" + nm] = strip(r)
g = V.build("F5", hh=True); pmf, comps = V.m2_table(g, cache="cache/m2_F5hh.npz")
Th = dict(T); G["F5"] = g
Th["F5"] = dict(comps=comps, M2=pmf, M0=V.m0_pmf(g, comps)[None, :], CRR=V.crr_table(g, comps), i=V.idx_of(comps, g["obs"]))
r = run(CAL, HOLD, tabs=Th, label="S-hh (F5 households collapsed)"); out["S_hh"] = strip(r)
print("    F5hh", g["sizes"], g["obs"], {m: (round(v["logp"], 2), round(v["p_adeq"], 4)) for m, v in r["events"]["F5"].items() if m != "obs"})
G["F5"] = V.build("F5")
Tx = dict(T); Tx["W1"] = T["W1x"]
r = run(CAL, HOLD, tabs=Tx, label="W1 exact bed assignments"); out["W1_exact"] = strip(r)
r = run({"air": CAL["air"], "room": ["W1"]}, HOLD, label="S-noW2"); out["S_noW2"] = strip(r)
r = run({"air": CAL["air"], "room": ["W2"]}, HOLD, label="room calibrated on W2 only"); out["room_W2only"] = strip(r)

# ------------------------------------------------------------------ 7. descriptive near/far RR
def nf(e):
    g = G[e]; nb = np.isin(np.arange(len(g["sizes"])), g["near_bins"]); s = np.array(g["sizes"]); o = np.array(g["obs"])
    return int(o[nb].sum()), int(s[nb].sum()), int(o[~nb].sum()), int(s[~nb].sum())
def pool(evs):
    y, v = [], []
    for e in evs:
        a, n1, b, n2 = nf(e)
        if a == 0 or b == 0:
            a, b, n1, n2 = a + .5, b + .5, n1 + .5, n2 + .5
        y.append(np.log((a / n1) / (b / n2))); v.append(1 / a - 1 / n1 + 1 / b - 1 / n2)
    y, v = np.array(y), np.array(v); w = 1 / v
    fe = (w * y).sum() / w.sum(); Q = (w * (y - fe) ** 2).sum(); k = len(y)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    wr = 1 / (v + tau2); re = (wr * y).sum() / wr.sum(); se = 1 / np.sqrt(wr.sum())
    return dict(RR_re=float(np.exp(re)), lo=float(np.exp(re - 1.96 * se)), hi=float(np.exp(re + 1.96 * se)), RR_fe=float(np.exp(fe)),
                Q=float(Q), df=k - 1, pQ=float(chi2.sf(Q, k - 1)) if k > 1 else None)
out["nearfar"] = {e: nf(e) for e in EV}
out["RR"] = {"all": pool(EV), "air": pool(air), "covid": pool(["F5", "F6", "F7", "F8", "P1"]), "sars1": pool(["F1", "W1", "W2"]),
             "flu": pool(["F3", "F4"]), "tb": pool(["F2"]), "gradeA": pool([e for e in EV if e != "W2"])}
print("near/far:", out["nearfar"]); print("RR:", json.dumps(out["RR"]))
json.dump(out, open("t2_results.json", "w"), indent=1, default=float)
