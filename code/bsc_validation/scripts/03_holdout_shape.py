"""Hold-out shape tests: compare the stored predictive distributions with the
observed stratum counts.  This is the first script that reads
data/heldout_outcomes.json.

Usage: 03_holdout_shape.py [tag]    (tag selects data/02_predictive_<tag>.json)
"""
import sys

import numpy as np

import _common as C
from valmod import decision as DC

EVENTS = ["T1", "T2", "T5", "C1"]


def main(tag="primary"):
    pred = C.load_json("02_predictive.json" if tag == "primary" else f"02_predictive_{tag}.json")
    obs = C.load_json("heldout_outcomes.json")
    try:
        ex = C.load_json("02b_m1_exact.json")["events"] if tag == "primary" else {}
    except FileNotFoundError:
        ex = {}
    out = {"stamp": C.stamp(), "tag": tag, "events": {}, "pooled": {}, "secondary": {}}
    pm = {m: [] for m in ("M0", "M1", "M2", "M1cf")}
    kobs = []
    for ev in EVENTS:
        e = pred["events"][ev]
        o = obs["shape_primary"][ev]
        n1, n2, K = e["n_near"], e["n_far"], e["K"]
        assert (n1, n2) == (o["near_n"], o["far_n"]) and K == o["near_k"] + o["far_k"]
        k = o["near_k"]
        kobs.append(k)
        rec = {"label": o["label"], "obs": o, "n_near": n1, "n_far": n2, "K": K}
        pmfs = {"M0": e["M0"]["pmf"], "M2": e["M2"]["pmf"], "M1cf": e["M1"]["pmf"]}
        use_exact = bool(ex.get(ev, {}).get("replace_closed_form", False))
        if use_exact:
            pe = np.array(ex[ev]["pmf_sim"], float)
            # smooth empty cells of the Monte Carlo pmf with half a count
            pe = (pe * ex[ev]["n_cond"] + 0.5 / len(pe)) / (ex[ev]["n_cond"] + 0.5)
            pmfs["M1"] = (pe / pe.sum()).tolist()
        else:
            pmfs["M1"] = e["M1"]["pmf"]
        rec["M1_uses_exact_simulator"] = use_exact
        for m, p in pmfs.items():
            rec[m] = DC.event_row(p, n1, n2, K, k)
            pm[m].append(p)
        for m in ("M1", "M2", "M1cf"):
            rec[m]["delta_vs_M0"] = rec[m]["log_score"] - rec["M0"]["log_score"]
            rec[m]["label"] = DC.label(rec[m]["p_two_sided"], rec["M0"]["p_two_sided"])
        out["events"][ev] = rec
    for m in ("M1", "M2", "M1cf"):
        d, p, per = DC.pooled_delta_test(pm[m], pm["M0"], kobs)
        nfail = sum(out["events"][ev][m]["p_two_sided"] < 0.05 for ev in EVENTS)
        out["pooled"][m] = dict(delta=d, p_under_M0=p, per_event=dict(zip(EVENTS, per)),
                                A1_pass=bool(d > 0 and p < 0.05), A2_failures_single_events=int(nfail),
                                failed=[ev for ev in EVENTS if out["events"][ev][m]["p_two_sided"] < 0.05])
    out["pooled"]["M0"] = dict(failed=[ev for ev in EVENTS if out["events"][ev]["M0"]["p_two_sided"] < 0.05])
    for key, e in pred.get("secondary", {}).items():
        o = obs["shape_secondary"][key]
        n1, n2, K, k = e["n_near"], e["n_far"], e["K"], o["near_k"]
        assert (n1, n2) == (o["near_n"], o["far_n"]) and K == o["near_k"] + o["far_k"]
        out["secondary"][key] = {m: DC.event_row(e[m]["pmf"], n1, n2, K, k) for m in ("M0", "M1", "M2")}
        out["secondary"][key]["obs"] = o
    C.save_json("03_holdout_shape.json" if tag == "primary" else f"03_holdout_shape_{tag}.json", out)
    print(f"{'event':5s} {'obs':>12s} {'RRobs':>6s} | " + " | ".join(
        f"{m:>4s}: RRpred (90% PI)   p    logS" for m in ("M0", "M1", "M2")))
    for ev in EVENTS:
        r = out["events"][ev]
        o = r["obs"]
        s = f"{ev:5s} {o['near_k']:>2d}/{o['near_n']:<2d} {o['far_k']:>2d}/{o['far_n']:<2d} {r['M0']['rr_obs']:6.2f} | "
        for m in ("M0", "M1", "M2"):
            x = r[m]
            s += f"{x['rr_pred']:7.2f} ({x['rr_pred_90'][0]:.2f}-{x['rr_pred_90'][1]:.2f}) {x['p_two_sided']:.4f} {x['log_score']:7.2f} | "
        print(s + ("  [M1 exact]" if r["M1_uses_exact_simulator"] else ""))
    for m in ("M1", "M2", "M1cf"):
        print(m, out["pooled"][m])
    for key, r in out["secondary"].items():
        print(key, {m: (round(r[m]["rr_pred"], 2), round(r[m]["p_two_sided"], 4)) for m in ("M0", "M1", "M2")},
              "obs RR", round(r["M0"]["rr_obs"], 2))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "primary")
