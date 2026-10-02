"""Pre-registered sensitivity analyses (section 9).  Each variant recomputes the
hold-out conditional predictive distributions and the shape tests.  M1 uses the
closed form here (label M1cf); see POSTHOC_LOG P4.  Results never replace the
primary analysis."""
import numpy as np

import _common as C
from valmod import decision as DC
from valmod import events as E
from valmod import predict as P
from valmod import stats as S

EVENTS = ["T1", "T2", "T5", "C1"]
obs = C.load_json("heldout_outcomes.json")["shape_primary"]
N_DRAW = 300


def aniso_posterior(model, n=400, seed=0):
    cal = C.load_json("01_cal_T4_S-aniso.json")
    rng = np.random.default_rng(seed)
    ll = np.array(cal[model]["profile_logL"])
    th = np.array(cal["theta"])
    logD = np.array(cal["logD"])
    w = np.exp(ll - ll.max())
    w /= w.sum()
    idx = np.unravel_index(rng.choice(w.size, size=n, p=w.ravel()), w.shape)
    step = logD[1] - logD[0]
    return dict(D=10 ** (logD[idx[1]] + rng.uniform(-step / 2, step / 2, n)), theta=th[idx[0]],
                theta_hat=cal[model]["theta_hat"], D_hat=cal[model]["D_hat"], logL=cal[model]["logL"])


def run(tag, veh="primary", room="primary", kw=None, theta=False, nu=False):
    kw = kw or {}
    res = {"tag": tag, "events": {}}
    pm = {m: [] for m in ("M0", "M1", "M2")}
    kobs = []
    for ev in EVENTS:
        cls = E.CLASS[ev]
        n1, n2 = E.STRATA_N[ev]
        K = E.TOTALS[ev]
        k = obs[ev]["near_k"]
        kobs.append(k)
        row = {}
        pmfs = {"M0": S.hypergeom_pmf(n1, n2, K)}
        for m in ("M1", "M2"):
            th, nuv = None, None
            if cls == "vehicle" and theta:
                post = aniso_posterior(m, seed=C.SEED + 7)
                th = post["theta"]
                row[m + "_theta_hat"] = post["theta_hat"]
            else:
                post = P.posterior_D(cls, m, veh if cls == "vehicle" else room, n=400, seed=C.SEED + 7)
            if nu and cls == "vehicle":
                cal = P.load_cal("vehicle")
                c = cal["cells"]
                i01 = [i for i in range(len(c["dr"])) if c["dr"][i] == 0 and c["dc"][i] == 1][0]
                nuv = (c["k"][i01] / c["n"][i01]) / cal[m]["p"][i01]
                row[m + "_nu"] = nuv
            r = P.shape_predictive(ev, m, post["D"], n_draw=N_DRAW, seed=C.SEED + 13,
                                   builder_kwargs=kw.get(ev, {}), theta_draws=th, nu=nuv)
            pmfs[m] = r["pmf"]
        for m, p in pmfs.items():
            x = DC.event_row(p, n1, n2, K, k)
            row[m] = dict(rr_pred=x["rr_pred"], rr_pred_90=x["rr_pred_90"], p=x["p_two_sided"],
                          log_score=x["log_score"])
            pm[m].append(p)
        res["events"][ev] = row
    for m in ("M1", "M2"):
        d, p, per = DC.pooled_delta_test(pm[m], pm["M0"], kobs)
        res[m] = dict(delta=d, p_under_M0=p, A1=bool(d > 0 and p < 0.05),
                      failed=[ev for ev in EVENTS if res["events"][ev][m]["p"] < 0.05])
    print(tag, {m: (round(res[m]["delta"], 2), round(res[m]["p_under_M0"], 5), res[m]["failed"]) for m in ("M1", "M2")},
          {ev: {m: (round(res["events"][ev][m]["rr_pred"], 2), round(res["events"][ev][m]["p"], 4)) for m in ("M1", "M2")}
           for ev in EVENTS}, flush=True)
    return res


if __name__ == "__main__":
    out = {"stamp": C.stamp(), "variants": {}}
    jobs = [
        ("primary (closed form for M1)", {}),
        ("S-adj", dict(veh="S-adj")),
        ("S-row", dict(veh="S-row")),
        ("S-logit", dict(veh="S-logit")),
        ("S-pitch", dict(veh="S-pitch")),
        ("S-aniso", dict(theta=True)),
        ("S-R1-upper", dict(room="S-R1-upper")),
        ("S-R1-lower", dict(room="S-R1-lower")),
        ("S-5K", dict(kw={"T5": {"index_seat": (5, 3)}})),
        ("S-Luo", dict(kw={"T2": {"minutes": 150.0}})),
        ("S-corner", dict(kw={"C1": {"corner": True}})),
        ("S-nu", dict(nu=True)),
    ]
    for tag, args in jobs:
        out["variants"][tag] = run(tag, **args)
        C.save_json("09_sensitivity.json", out)
