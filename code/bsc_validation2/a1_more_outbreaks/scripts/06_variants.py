"""Sensitivity variants of pre-registration section 8 (S-out, S-src, S-hh, S-pitch, S-kappa, S-noW2, S-int).
Each variant re-runs the whole primary pipeline (calibration + hold-out + decision) with the variant tables."""
import sys, os, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
from a1 import events as E, analysis as A, pipeline as P
base = {e: E.build(e, {"n_cfg": 1})[0]["obs"] for e in E.EVENTS}
AIR = E.CAL["air"] + E.HOLD["air"]
V = {  # name: (dict ev -> (opt, tag))
 "S-out F5 confirmed+probable": {"F5": ({"outcome": "wide"}, "wide")},
 "S-out F6 flight-associated only (8)": {"F6": ({"outcome": "wide"}, "wide")},
 "S-out F3 with suspected case": {"F3": ({"outcome": "wide"}, "wide")},
 "S-out F7 with passenger G": {"F7": ({"outcome": "wide"}, "wide")},
 "S-out all four": {e: ({"outcome": "wide"}, "wide") for e in ("F5", "F6", "F3", "F7")},
 "S-src F6 three symptomatic sources": {"F6": ({"src": "alt1"}, "src1")},
 "S-src F6 six infectious": {"F6": ({"src": "alt2"}, "src2")},
 "S-src F3 suspected cases as sources": {"F3": ({"src": "alt1"}, "src1")},
 "S-src F7 passenger A only": {"F7": ({"src": "alt1"}, "src1")},
 "S-hh F5 households collapsed": {"F5": ({"hh": True}, "hh")},
 "S-int F1 interviewed only": {"F1": ({"interviewed": True}, "int")},
 "S-pitch 0.75 m": {e: ({"ay": 0.75}, "ay075") for e in AIR},
 "S-pitch 0.86 m": {e: ({"ay": 0.86}, "ay086") for e in AIR},
 "S-kappa cabin ACH 5-15": {e: ({"cabin_ach": [5, 15]}, "ach5") for e in AIR},
}
res = {}
def one(name, cal, hold, obs, tags):
    r = P.run(cal, hold, obs, tags=tags)
    res[name] = dict(pooled=r["pooled"], post=r["post"],
                     events={e: dict(obs=v["obs"], sizes=v["sizes"], **{m: dict(logp=v[m]["logp"], p_adeq=v[m]["p_adeq"]) for m in A.MODELS}) for e, v in r["events"].items()})
    for t in ("M2", "M3", "M1"):
        a = r["pooled"][t]["all"]; ai = r["pooled"][t].get("air", {})
        print("%-38s %s d0 %+.2f (p %.4f) dC %+.2f (pB2 %.3f, mirror %.4f) | air d0 %+.2f dC %+.2f mirror %.3f | %s fails %s" % (
            name, t, a["delta0"], a["p_B1"], a["deltaC"], a["p_B2"], a["p_mirror"], ai.get("delta0", np.nan), ai.get("deltaC", np.nan),
            ai.get("p_mirror", np.nan), r["pooled"][t]["decision"]["outcome"], r["pooled"][t]["decision"]["B3_failures"]), flush=True)
    print("     post", {k: v["map"] for k, v in r["post"].items() if k.split("|")[1] in ("M2", "CRR")},
          {e: (v["obs"], round(v["M2"]["p_adeq"], 3)) for e, v in r["events"].items()})
one("primary (reference)", E.CAL, E.HOLD, base, {})
for name, spec in V.items():
    obs = dict(base); tags = {}
    for e, (opt, tag) in spec.items():
        obs[e] = E.build(e, dict(opt, n_cfg=1))[0]["obs"]; tags[e] = tag
    one(name, E.CAL, E.HOLD, obs, tags)
one("S-noW2 room class calibrated on W1 only", {"air": E.CAL["air"], "room": ["W1"]}, E.HOLD, base, {})
json.dump(res, open(os.path.join(ROOT, "out", "06_variants.json"), "w"), indent=1)
