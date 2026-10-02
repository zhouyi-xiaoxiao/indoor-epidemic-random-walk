"""One full run of the primary pipeline for a given calibration/hold-out assignment and
set of table tags (used by the primary analysis and by every sensitivity variant)."""
from __future__ import annotations
import numpy as np
from . import events as E, analysis as A

lg = lambda a: np.log(np.clip(a, 1e-300, None))


def run(cal, hold, obs, tags=None, weights=None, n_mc=200000, seed=101, targets=("M2", "M3", "M1"),
        crr_fixed=None):
    """cal/hold: dict cls -> event list.  obs: dict ev -> bin counts.  tags: dict ev -> table tag.
    weights: optional dict (cls, model) -> grid weights that replace the calibration posterior.
    crr_fixed: optional dict ev -> pmf replacing the CRR predictive (published constant)."""
    tags = tags or {}
    tg = lambda e: tags.get(e, "primary")
    post = {}
    for cls in cal:
        for m in A.MODELS:
            if weights and (cls, m) in weights:
                post[cls, m] = np.asarray(weights[cls, m])
            else:
                post[cls, m] = A.posterior(m, cal[cls], obs, tags)[0]
    evs = [(e, cls) for cls in hold for e in hold[cls]]
    pm = {m: [A.predictive(e, m, post[c, m], tg(e)) for e, c in evs] for m in A.MODELS}
    if crr_fixed:
        pm["CRR"] = [crr_fixed[e] for e, _ in evs]
    idx_obs = [A.comp_index(A.tab(e, tg(e)), obs[e]) for e, _ in evs]
    out = dict(events={}, pooled={}, post={"%s|%s" % k: A.summarize(k[1], v) for k, v in post.items() if k[1] != "M0"})
    for k, (e, c) in enumerate(evs):
        t = A.tab(e, tg(e))
        r = dict(cls=c, sizes=t["sizes"].tolist(), obs=list(map(int, obs[e])), K=int(t["K"]))
        for m in A.MODELS:
            p = pm[m][k]
            r[m] = dict(logp=float(lg(p[idx_obs[k]])), p_adeq=A.two_sided_p(p, idx_obs[k]),
                        expected=A.expected_counts(t, p).tolist())
        out["events"][e] = r
    sub = {"all": list(range(len(evs))), "air": [k for k, (e, c) in enumerate(evs) if c == "air"],
           "room": [k for k, (e, c) in enumerate(evs) if c == "room"]}
    for target in targets:
        res = {}
        for name, ks in sub.items():
            if not ks:
                continue
            a = [pm[target][k] for k in ks]; b0 = [pm["M0"][k] for k in ks]; bc = [pm["CRR"][k] for k in ks]
            io = [idx_obs[k] for k in ks]
            d0 = float(sum(lg(x[i]) - lg(y[i]) for x, y, i in zip(a, b0, io)))
            dc = float(sum(lg(x[i]) - lg(y[i]) for x, y, i in zip(a, bc, io)))
            ref0 = A.delta_reference(a, b0, b0, n_mc, seed)
            refc = A.delta_reference(a, bc, bc, n_mc, seed + 1)
            refm = A.delta_reference(a, bc, a, n_mc, seed + 2)
            refm0 = A.delta_reference(a, b0, a, n_mc, seed + 3)
            res[name] = dict(delta0=d0, p_B1=float((ref0 >= d0 - 1e-12).mean()),
                             deltaC=dc, p_B2=float((refc >= dc - 1e-12).mean()),
                             p_mirror=float((refm <= dc + 1e-12).mean()),
                             p_delta0_low_under_target=float((refm0 <= d0 + 1e-12).mean()),
                             E_delta0_if_target=float(refm0.mean()), E_deltaC_if_target=float(refm.mean()),
                             E_deltaC_if_CRR=float(refc.mean()), E_delta0_if_M0=float(ref0.mean()),
                             n_events=len(ks))
        fails = [e for (e, c) in evs if out["events"][e][target]["p_adeq"] < 0.05]
        al = res["all"]
        B1 = al["delta0"] > 0 and al["p_B1"] < 0.05
        B2 = al["deltaC"] > 0 and al["p_B2"] < 0.05
        mir = al["deltaC"] < 0 and al["p_mirror"] < 0.05
        if len(fails) >= 2: V = "V5"
        elif not B1: V = "V4"
        elif mir: V = "V3"
        elif B2 and not fails: V = "V1"
        else: V = "V2"
        res["decision"] = dict(B1=bool(B1), B2=bool(B2), mirror=bool(mir), B3_failures=fails, outcome=V)
        out["pooled"][target] = res
    return out
