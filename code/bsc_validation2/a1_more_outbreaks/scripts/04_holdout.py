"""Primary hold-out evaluation (pre-registration 6-7).  First script that reads hold-out stratum counts."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from a1 import events as E, pipeline as P
obs = {e: E.build(e, {"n_cfg": 1})[0]["obs"] for e in E.EVENTS}
res = P.run(E.CAL, E.HOLD, obs)
json.dump(res, open(os.path.join(ROOT, "out", "04_holdout.json"), "w"), indent=1)
for e, r in res["events"].items():
    print(e, r["sizes"], r["obs"], {m: (round(r[m]["logp"], 2), round(r[m]["p_adeq"], 4), [round(x, 1) for x in r[m]["expected"]]) for m in ("M0", "CRR", "M2", "M3", "M1")})
print(json.dumps(res["pooled"], indent=1))
