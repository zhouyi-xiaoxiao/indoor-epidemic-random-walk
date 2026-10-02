"""Pre-registration 6.5: observed and predicted near/far risk ratios per event, pooled by inverse-variance
weights on the log scale, Cochran's Q, by pathogen and by setting.  Predicted values: primary posterior
for hold-out events, leave-one-out posterior for calibration events (flagged)."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from scipy import stats as st
from a1 import events as E, analysis as A
C = json.load(open(os.path.join(ROOT, "out", "10_corrected.json")))   # round-off-safe kernel (POSTHOC_LOG P5)
hold = C["primary_corrected"]["events"]
loo = C["S-LOO_corrected"]
rows = {}
for e in E.EVENTS:
    meta = E.build(e, {"n_cfg": 1})[0]
    a, n1, c, n0 = A.near_far(meta["sizes"], meta["obs"], meta["near_bins"])
    cc = 0.5 if min(a, c, n1 - a, n0 - c) == 0 else 0.0
    rr = ((a + cc) / (n1 + cc)) / ((c + cc) / (n0 + cc))
    se = np.sqrt(1 / (a + cc) - 1 / (n1 + cc) + 1 / (c + cc) - 1 / (n0 + cc))
    src = hold[e] if e in hold else loo[e]
    pred = {}
    for m in A.MODELS:
        ex = np.array(src[m]["expected"])
        nbm = np.isin(np.arange(len(ex)), meta["near_bins"]); pa, pc = ex[nbm].sum(), ex[~nbm].sum()
        pred[m] = float((pa / n1) / (pc / n0)) if pc > 0 else float("inf")
    rows[e] = dict(label=meta["label"], pathogen=meta["pathogen"], cls=meta["cls"], role="hold-out" if e in hold else "calibration (LOO prediction)",
                   near=[a, n1], far=[c, n0], cc=cc, RR=float(rr), logRR=float(np.log(rr)), se=float(se),
                   lo=float(np.exp(np.log(rr) - 1.96 * se)), hi=float(np.exp(np.log(rr) + 1.96 * se)),
                   fisher_p=float(st.fisher_exact([[a, n1 - a], [c, n0 - c]])[1]), pred_RR=pred)
def pool(evs):
    y = np.array([rows[e]["logRR"] for e in evs]); w = 1 / np.array([rows[e]["se"] for e in evs]) ** 2
    m = (w * y).sum() / w.sum(); s = 1 / np.sqrt(w.sum()); Q = float((w * (y - m) ** 2).sum()); df = len(evs) - 1
    # DerSimonian-Laird random effects
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())) if df > 0 else 0.0
    wr = 1 / (1 / w + tau2); mr = (wr * y).sum() / wr.sum(); sr = 1 / np.sqrt(wr.sum())
    # Mantel-Haenszel (no continuity correction needed)
    num = sum(rows[e]["near"][0] * rows[e]["far"][1] / (rows[e]["near"][1] + rows[e]["far"][1]) for e in evs)
    den = sum(rows[e]["far"][0] * rows[e]["near"][1] / (rows[e]["near"][1] + rows[e]["far"][1]) for e in evs)
    return dict(events=evs, n=len(evs), RR_fixed=float(np.exp(m)), lo=float(np.exp(m - 1.96 * s)), hi=float(np.exp(m + 1.96 * s)),
                Q=Q, df=df, p_Q=float(st.chi2.sf(Q, df)) if df > 0 else None, I2=float(max(0, (Q - df) / Q)) if df > 0 and Q > 0 else 0.0,
                tau2=float(tau2), RR_random=float(np.exp(mr)), lo_r=float(np.exp(mr - 1.96 * sr)), hi_r=float(np.exp(mr + 1.96 * sr)),
                RR_MH=float(num / den) if den > 0 else None)
groups = {"all 11 events": E.EVENTS, "aircraft (8)": E.EVENTS[:8], "rooms (3)": E.EVENTS[8:],
          "new grade-A only (10, without W2)": [e for e in E.EVENTS if e != "W2"],
          "calibration aircraft (pre-2020)": E.CAL["air"], "hold-out aircraft (SARS-CoV-2)": E.HOLD["air"]}
for p in sorted(set(E.PATH.values())):
    groups["pathogen: " + p] = [e for e in E.EVENTS if E.PATH[e] == p]
pooled = {k: pool(v) for k, v in groups.items()}
json.dump(dict(events=rows, pooled=pooled), open(os.path.join(ROOT, "out", "07_descriptive.json"), "w"), indent=1)
for e, r in rows.items():
    print("%-3s %-34s %-24s %3d/%-3d vs %3d/%-3d RR %5.2f (%.2f-%.2f) Fisher p %.4f | pred M0 %.2f CRR %.2f M2 %.2f M3 %.2f M1 %.2f | %s" % (
        e, r["label"], r["pathogen"], *r["near"], *r["far"], r["RR"], r["lo"], r["hi"], r["fisher_p"],
        *[r["pred_RR"][m] for m in ("M0", "CRR", "M2", "M3", "M1")], r["role"]))
for k, v in pooled.items():
    print("%-40s n=%d fixed RR %.2f (%.2f-%.2f) Q %.1f df %d p %s I2 %.2f | random %.2f (%.2f-%.2f) | MH %.2f" % (
        k, v["n"], v["RR_fixed"], v["lo"], v["hi"], v["Q"], v["df"], "%.3f" % v["p_Q"] if v["p_Q"] is not None else "-", v["I2"],
        v["RR_random"], v["lo_r"], v["hi_r"], v["RR_MH"]))
