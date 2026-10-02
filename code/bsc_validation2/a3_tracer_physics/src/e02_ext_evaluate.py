"""EXT-F step 2: compare with the observed bin vectors of F5-F8 (first script that reads them).
Writes results/E2_ext_evaluation.json."""
import os, sys, json
import numpy as np
import a3lib as A
sys.path.insert(0, os.path.normpath(os.path.join(A.ROOT, "..", "a1_more_outbreaks", "src")))
from a1 import events as E1
FL = ["F5", "F6", "F7", "F8"]; MODELS = ["P", "P_econ", "P_D708", "P_D112", "M0", "CRR2", "CRR3"]
pred = A.load_json("E1_ext_predictions.json"); T = np.load(os.path.join(A.RES, "E1_ext_tables.npz"))
lg = lambda a: np.log(np.clip(a, 1e-300, None))
def two_sided(p, i): return float(p[p <= p[i] * (1 + 1e-12)].sum())
out = {"stamp": A.stamp(), "events": {}}
idx = {}
for ev in FL:
    meta, _ = E1.build(ev, {"n_cfg": 1})
    obs = list(map(int, meta["obs"])); comps = T[f"{ev}|comps"]
    i = int(np.flatnonzero((comps == np.array(obs)).all(axis=1))[0]); idx[ev] = i
    nb = np.isin(np.arange(len(obs)), meta["near_bins"]); sz = np.array(meta["sizes"])
    kn, kf = int(np.array(obs)[nb].sum()), int(np.array(obs)[~nb].sum()); n1, n2 = int(sz[nb].sum()), int(sz[~nb].sum())
    r = dict(label=meta["label"], sizes=meta["sizes"], obs=obs, near=[kn, n1], far=[kf, n2],
             rr_obs=(kn / n1) / (kf / n2) if kf > 0 else None, rr_obs_cc=((kn + 0.5) / (n1 + 1)) / ((kf + 0.5) / (n2 + 1)), models={})
    for m in MODELS:
        p = T[f"{ev}|{m}"]
        r["models"][m] = dict(logp=float(lg(p[i])), p_two_sided=two_sided(p, i), rr_expected=pred[ev]["rr_expected"][m],
                              expected=(p[:, None] * comps).sum(axis=0).tolist())
    out["events"][ev] = r
rng = np.random.default_rng(A.SEED + 41); NMC = 200000
out["pooled"] = {}
for target in ("P", "P_econ", "P_D708", "P_D112"):
    res = {}
    for c in ("M0", "CRR2", "CRR3"):
        per = {ev: float(lg(T[f"{ev}|{target}"][idx[ev]]) - lg(T[f"{ev}|{c}"][idx[ev]])) for ev in FL}
        d = sum(per.values())
        ref = sum((lambda j, ev=ev: lg(T[f"{ev}|{target}"][j]) - lg(T[f"{ev}|{c}"][j]))(rng.choice(len(T[f"{ev}|{c}"]), NMC, p=T[f"{ev}|{c}"] / T[f"{ev}|{c}"].sum())) for ev in FL)
        mir = sum((lambda j, ev=ev: lg(T[f"{ev}|{target}"][j]) - lg(T[f"{ev}|{c}"][j]))(rng.choice(len(T[f"{ev}|{target}"]), NMC, p=T[f"{ev}|{target}"] / T[f"{ev}|{target}"].sum())) for ev in FL)
        res[c] = dict(delta=float(d), per_event=per, p_under_competitor=float((ref >= d - 1e-12).mean()),
                      p_low_under_target=float((mir <= d + 1e-12).mean()), E_delta_if_target=float(mir.mean()), E_delta_if_competitor=float(ref.mean()))
    fails = [ev for ev in FL if out["events"][ev]["models"][target]["p_two_sided"] < 0.05]
    B0 = res["M0"]["delta"] > 0 and res["M0"]["p_under_competitor"] < 0.05
    BC = all(res[c]["delta"] > 0 and res[c]["p_under_competitor"] < 0.05 for c in ("CRR2", "CRR3"))
    o = "E5" if len(fails) >= 2 else "E4" if len(fails) == 1 else "E3" if not B0 else "E1" if BC else "E2"
    res["decision"] = dict(failures=fails, beats_M0=bool(B0), beats_CRR=bool(BC), outcome=o)
    out["pooled"][target] = res
out["competitor_failures"] = {m: [ev for ev in FL if out["events"][ev]["models"][m]["p_two_sided"] < 0.05] for m in ("M0", "CRR2", "CRR3")}
# reference: the second outbreak test's interim hold-out file (calibrated M2 and calibrated CRR), copied not recomputed
try:
    h = json.load(open(os.path.join(A.ROOT, "..", "a1_more_outbreaks", "out", "04_holdout.json")))
    out["a1_reference"] = {ev: {m: dict(logp=h["events"][ev][m]["logp"], p=h["events"][ev][m]["p_adeq"]) for m in ("M0", "CRR", "M2")} for ev in FL}
    out["a1_reference"]["post"] = {k: v for k, v in h["post"].items() if k.startswith("air")}
except Exception as e:
    out["a1_reference"] = str(e)
A.save_json("E2_ext_evaluation.json", out)
for ev in FL:
    r = out["events"][ev]; print(ev, r["label"], r["obs"], "near", r["near"], "far", r["far"], {m: (round(v["logp"], 2), round(v["p_two_sided"], 4), round(v["rr_expected"], 2)) for m, v in r["models"].items()})
for t, r in out["pooled"].items():
    print(t, {c: (round(r[c]["delta"], 2), r[c]["p_under_competitor"], round(r[c]["p_low_under_target"], 4)) for c in ("M0", "CRR2", "CRR3")}, r["decision"])
print("competitor failures", out["competitor_failures"]); print(json.dumps(out["a1_reference"], indent=0)[:1500])
