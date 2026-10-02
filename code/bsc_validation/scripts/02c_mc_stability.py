"""Monte Carlo stability of the hold-out shape results: repeat the M2 (and the
closed-form M1) predictive computation with five other seeds for the posterior
draws and the seating configurations."""
import numpy as np

import _common as C
from valmod import decision as DC
from valmod import events as E
from valmod import predict as P
from valmod import stats as S

obs = C.load_json("heldout_outcomes.json")["shape_primary"]
EVENTS = ["T1", "T2", "T5", "C1"]
out = {"stamp": C.stamp(), "seeds": [], "M2": [], "M1cf": []}
for s in range(1, 6):
    row = {"M1": {}, "M2": {}}
    pm = {"M0": [], "M1": [], "M2": []}
    kobs = []
    for ev in EVENTS:
        cls = E.CLASS[ev]
        n1, n2 = E.STRATA_N[ev]
        K = E.TOTALS[ev]
        k = obs[ev]["near_k"]
        kobs.append(k)
        pm["M0"].append(S.hypergeom_pmf(n1, n2, K))
        for m in ("M1", "M2"):
            post = P.posterior_D(cls, m, n=400, seed=C.SEED + 1000 * s)
            r = P.shape_predictive(ev, m, post["D"], n_draw=300, seed=C.SEED + 2000 * s)
            pm[m].append(r["pmf"])
            row[m][ev] = S.two_sided_p(r["pmf"] / r["pmf"].sum(), k)
    for m, key in (("M2", "M2"), ("M1", "M1cf")):
        d, p, _ = DC.pooled_delta_test(pm[m], pm["M0"], kobs)
        out[key].append(dict(p=row[m], delta=d, p_pooled=p))
    out["seeds"].append(s)
    print(s, {ev: round(row["M2"][ev], 4) for ev in EVENTS}, round(out["M2"][-1]["delta"], 2),
          out["M2"][-1]["p_pooled"], flush=True)
    C.save_json("02c_mc_stability.json", out)
