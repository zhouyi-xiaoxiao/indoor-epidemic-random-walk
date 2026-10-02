"""Secondary analysis: leave-one-event-out over the stratified events of each
setting class (pre-registration 5.5).

For every event the conditional likelihood of its observed stratum split is
computed on a grid of D (seating, index seat and kappa marginalised).  For the
calibration events the profile likelihood of script 01 is used.  Leaving one
event out, D is weighted by the product of the other events' likelihoods (flat
prior on log10 D) and the left-out event is scored.

M1 is evaluated with the closed form (see POSTHOC_LOG P4) and is reported for
information only.
"""
import numpy as np

import _common as C
from valmod import events as E
from valmod import predict as P
from valmod import stats as S

LOGD = np.linspace(-2, 4, 31)
N_DRAW = 120
obs = C.load_json("heldout_outcomes.json")["shape_primary"]


def cal_curve(cls, model):
    cal = P.load_cal(cls)
    logD = np.array(cal["logD"])
    m = cal[model]
    ll = np.array(m["marg_logL"]) if (cls == "vehicle" and model == "M2") else np.array(m["profile_logL"])
    return np.interp(LOGD, logD, ll), cal["M0"]["logL"]


def main():
    out = {"stamp": C.stamp(), "logD": LOGD.tolist(), "classes": {}}
    groups = {"vehicle": ["T4", "T1", "T2", "T5"], "room": ["R1", "C1"]}
    for cls, evs in groups.items():
        res = {}
        for model in ("M1", "M2"):
            ll = {}
            ll0 = {}
            for ev in evs:
                if ev in ("T4", "R1"):
                    ll[ev], ll0[ev] = cal_curve(cls, model)
                else:
                    n1, n2 = E.STRATA_N[ev]
                    K = E.TOTALS[ev]
                    k = obs[ev]["near_k"]
                    v = []
                    for ld in LOGD:
                        r = P.shape_predictive(ev, model, None, n_draw=N_DRAW, seed=C.SEED + 71,
                                               fixed_D=10 ** ld)
                        v.append(np.log(max(r["pmf"][k], 1e-300)))
                    ll[ev] = np.array(v)
                    ll0[ev] = float(np.log(S.hypergeom_pmf(n1, n2, K)[k]))
                print(cls, model, ev, "max logL", round(float(ll[ev].max()), 2), "at D",
                      round(10 ** LOGD[int(np.argmax(ll[ev]))], 3), "M0", round(ll0[ev], 2), flush=True)
            rows = {}
            for ev in evs:
                others = [e for e in evs if e != ev]
                lw = sum(ll[e] for e in others)
                w = np.exp(lw - lw.max())
                w /= w.sum()
                score = float(np.log(np.sum(w * np.exp(ll[ev] - ll[ev].max()))) + ll[ev].max())
                rows[ev] = dict(score=score, score_M0=ll0[ev], delta=score - ll0[ev],
                                D_train_mode=float(10 ** LOGD[int(np.argmax(w))]),
                                D_own_mode=float(10 ** LOGD[int(np.argmax(ll[ev]))]),
                                own_max=float(ll[ev].max()))
            # joint fit of all events
            lw = sum(ll[e] for e in evs)
            res[model] = dict(loo=rows, curves={e: ll[e].tolist() for e in evs},
                              joint_D_mode=float(10 ** LOGD[int(np.argmax(lw))]),
                              joint_max=float(lw.max()), joint_M0=float(sum(ll0.values())),
                              sum_own_max=float(sum(ll[e].max() for e in evs)))
            out["classes"][cls] = res
            C.save_json("07_loo.json", out)
            print(cls, model, {e: (round(r["delta"], 2), r["D_train_mode"], r["D_own_mode"]) for e, r in rows.items()},
                  "joint D", res[model]["joint_D_mode"], flush=True)


if __name__ == "__main__":
    main()
