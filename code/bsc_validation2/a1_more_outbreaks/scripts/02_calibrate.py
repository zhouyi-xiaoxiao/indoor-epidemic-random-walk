"""Calibration of each model's class parameter on the calibration events (pre-registration 5).
Reads the bin counts of calibration events only."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from a1 import events as E, analysis as A, engine as G
obs, sizes = {}, {}
for cls in ("air", "room"):
    for e in E.CAL[cls]:
        m, _ = E.build(e, {"n_cfg": 1})
        obs[e], sizes[e] = m["obs"], m["sizes"]
res = dict(obs=obs, sizes=sizes, classes={})
for cls in ("air", "room"):
    r = {}
    for model in A.MODELS:
        w, ll = A.posterior(model, E.CAL[cls], obs)
        r[model] = dict(w=w.tolist(), logL_grid=ll.tolist(), logL_max=float(ll.max()),
                        logL_marg=float(np.log(np.mean(np.exp(ll - ll.max()))) + ll.max()),
                        per_event={e: A.event_curve(e, model, obs[e]).tolist() for e in E.CAL[cls]},
                        summary=A.summarize(model, w))
        # in-sample adequacy and expected counts at the posterior
        r[model]["fit"] = {}
        for e in E.CAL[cls]:
            t = A.tab(e); pm = A.predictive(e, model, w); i = A.comp_index(t, obs[e])
            r[model]["fit"][e] = dict(logp=float(np.log(pm[i])), p_adeq=A.two_sided_p(pm, i),
                                      expected=A.expected_counts(t, pm).tolist())
        print(cls, model, "logLmax %.2f" % ll.max(), r[model]["summary"],
              {e: round(r[model]["fit"][e]["p_adeq"], 3) for e in E.CAL[cls]})
    res["classes"][cls] = r
json.dump(res, open(os.path.join(ROOT, "out", "02_calibration.json"), "w"), indent=1)
