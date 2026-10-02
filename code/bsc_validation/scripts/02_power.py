"""Hold-out shape predictions (conditional on event totals) and the power of the
decision rule.  Reads calibration files, event geometry, stratum sizes and
totals.  Does NOT read data/heldout_outcomes.json.

Outputs
  data/02_predictive.json        predictive pmfs of k_near | K for M0, M1, M2
  data/02_power.json             exact power of A1 and A2 (four single-event shape tests)
  PREREG_ADDENDUM_A1_power.md    (written by hand from 02_power.json)
"""
import itertools
import sys

import numpy as np

import _common as C
from valmod import events as E
from valmod import level as LV
from valmod import predict as P
from valmod import stats as S

N_DRAW = 300
EVENTS = ["T1", "T2", "T5", "C1"]


def predictive_all(variant="primary", extra=None, tag="primary"):
    out = {"stamp": C.stamp(), "variant": variant, "n_draw": N_DRAW, "events": {}}
    post = {}
    for cls in ("vehicle", "room"):
        for m in ("M1", "M2"):
            pd_ = P.posterior_D(cls, m, variant if cls == "vehicle" else
                                (variant if variant.startswith("S-R1") else "primary"),
                                n=400, seed=C.SEED + 7)
            post[(cls, m)] = pd_
            out.setdefault("posterior", {})[f"{cls}_{m}"] = P.summarize_posterior(pd_["w_logD"], pd_["logD"])
    for ev in EVENTS:
        cls = E.CLASS[ev]
        n1, n2 = E.STRATA_N[ev]
        K = E.TOTALS[ev]
        rec = {"n_near": n1, "n_far": n2, "K": K, "M0": {"pmf": S.hypergeom_pmf(n1, n2, K).tolist()}}
        for m in ("M1", "M2"):
            kw = (extra or {}).get(ev, {})
            r = P.shape_predictive(ev, m, post[(cls, m)]["D"], n_draw=N_DRAW, seed=C.SEED + 13,
                                   builder_kwargs=kw)
            rec[m] = {"pmf": r["pmf"].tolist(),
                      "rr_expected_q": np.quantile(r["rr_expected"], [0.05, 0.5, 0.95]).tolist(),
                      "p_near_mean": float(r["p_near"].mean()), "p_far_mean": float(r["p_far"].mean())}
        out["events"][ev] = rec
        print(ev, {m: np.round(rec[m]["rr_expected_q"], 2).tolist() for m in ("M1", "M2")}, flush=True)
    # secondary splits
    sec = {}
    for ev, split, n1, n2 in (("T1", "rows6-10", 23, 44), ("T2", "driver_side", 24, 20)):
        cls = E.CLASS[ev]
        rec = {"n_near": n1, "n_far": n2, "K": E.TOTALS[ev],
               "M0": {"pmf": S.hypergeom_pmf(n1, n2, E.TOTALS[ev]).tolist()}}
        for m in ("M1", "M2"):
            r = P.shape_predictive(ev, m, post[(cls, m)]["D"], n_draw=N_DRAW, seed=C.SEED + 17,
                                   split=split)
            rec[m] = {"pmf": r["pmf"].tolist(),
                      "rr_expected_q": np.quantile(r["rr_expected"], [0.05, 0.5, 0.95]).tolist()}
        sec[f"{ev}_{split}"] = rec
    out["secondary"] = sec
    C.save_json(f"02_predictive_{tag}.json" if tag != "primary" else "02_predictive.json", out)
    return out


def power(pred):
    """Exact enumeration over the outcome space of the four events."""
    evs = EVENTS
    pm = {m: [np.clip(np.array(pred["events"][e][m]["pmf"]), 1e-300, None) for e in evs]
          for m in ("M0", "M1", "M2")}
    sizes = [len(x) for x in pm["M0"]]
    grids = np.meshgrid(*[np.arange(s) for s in sizes], indexing="ij")
    res = {}

    def joint(m):
        j = np.ones(sizes)
        for a, g in enumerate(grids):
            j = j * pm[m][a][g]
        return j / j.sum()

    J = {m: joint(m) for m in pm}
    for cand in ("M1", "M2"):
        delta = np.zeros(sizes)
        for a, g in enumerate(grids):
            delta += np.log(pm[cand][a][g]) - np.log(pm["M0"][a][g])
        # A1: one-sided p under M0 of Delta >= observed, and Delta > 0
        d = np.round(delta.ravel(), 9)
        w0 = J["M0"].ravel()
        order = np.argsort(-d, kind="stable")
        ds = d[order]
        cum = np.cumsum(w0[order])
        last = np.searchsorted(-ds, -d, side="right") - 1    # P_M0(Delta >= d), ties included
        pval = cum[last].reshape(sizes)
        a1 = (pval < 0.05) & (delta > 0)
        # A2: per-event two-sided predictive p >= 0.05
        fails = np.zeros(sizes, int)
        for a, g in enumerate(grids):
            pe = np.array([S.two_sided_p(pm[cand][a] / pm[cand][a].sum(), k) for k in range(sizes[a])])
            fails += (pe[g] < 0.05).astype(int)
        # same adequacy test applied to M0 itself (how often would the null be 'adequate')
        fails0 = np.zeros(sizes, int)
        for a, g in enumerate(grids):
            pe = np.array([S.two_sided_p(pm["M0"][a] / pm["M0"][a].sum(), k) for k in range(sizes[a])])
            fails0 += (pe[g] < 0.05).astype(int)
        r = {}
        for truth in ("M0", "M1", "M2"):
            w = J[truth]
            r[truth] = dict(
                P_A1=float(w[a1].sum()),
                P_A2_nofail=float(w[fails == 0].sum()),
                P_A2_le1=float(w[fails <= 1].sum()),
                P_W1_like=float(w[a1 & (fails == 0)].sum()),
                P_W2_like=float(w[a1 & (fails <= 1)].sum()),
                P_W3_like=float(w[(~a1) & (fails <= 1)].sum()),
                P_W4_like=float(w[fails >= 2].sum()),
                P_null_adequate=float(w[fails0 == 0].sum()),
                E_delta=float((w * delta).sum()),
            )
        res[cand] = r
        res[cand + "_per_event_power"] = {}
        for a, e in enumerate(evs):
            # probability that M0 is rejected by its own two-sided test when cand is true
            p0 = pm["M0"][a] / pm["M0"][a].sum()
            pc = pm[cand][a] / pm[cand][a].sum()
            rej0 = np.array([S.two_sided_p(p0, k) < 0.05 for k in range(sizes[a])])
            rejc = np.array([S.two_sided_p(pc, k) < 0.05 for k in range(sizes[a])])
            res[cand + "_per_event_power"][e] = dict(
                P_reject_M0_if_model_true=float(pc[rej0].sum()),
                P_reject_model_if_M0_true=float(p0[rejc].sum()),
                P_reject_model_if_model_true=float(pc[rejc].sum()))
    return res


def level_dispersion_power(n_sim=20000, sigma=0.94):
    """Power of the dispersion test of section 7.2 if the Wells-Riley scaling is
    true (between-index log-sd sigma)."""
    rng = np.random.default_rng(C.SEED + 23)
    evs = ["H1", "O2", "R2", "C5"]
    n = {"H1": 60, "O2": 12, "R2": 13, "C5": 217}
    perms = list(itertools.permutations(range(4)))
    hits = 0
    hits_const = 0
    # structural factors evaluated at prior medians of kappa (as the analysis will do)
    s_med = np.array([LV.wr_factor(LV.LEVEL_EVENTS[e]["V"],
                                   np.median(E.draw_kappa(e, np.random.default_rng(1), 2000)),
                                   LV.LEVEL_EVENTS[e]["segments"]) for e in evs])
    for _ in range(n_sim):
        s, th = [], []
        for e in evs:
            ev = LV.LEVEL_EVENTS[e]
            kap = E.draw_kappa(e, rng)
            s.append(LV.wr_factor(ev["V"], kap, ev["segments"]))
            th.append(sum(ev["segments"]))
        s, th = np.array(s), np.array(th)
        logE = rng.normal(0, sigma, 4)
        logE[3] = np.log(np.mean(np.exp(rng.normal(0, sigma, 8))))        # 8 instructors pooled
        # typical superspreading level so that attack rates are in the observed range
        dose = np.exp(logE + np.log(0.5 / np.median(s))) * s
        k = np.array([rng.binomial(n[e], 1 - np.exp(-d)) for e, d in zip(evs, dose)])
        k = np.clip(k, 1, np.array([n[e] for e in evs]) - 1)
        h = -np.log1p(-k / np.array([n[e] for e in evs]))
        sd_obs = np.std(np.log(h / s_med))
        sd_perm = np.array([np.std(np.log(h / s_med[list(p)])) for p in perms])
        p = np.mean(sd_perm <= sd_obs + 1e-12)
        hits += p < 0.05
        sd_const = np.std(np.log(h / th))
        hits_const += sd_obs < sd_const
    return dict(power_perm_test=hits / n_sim, P_sd_smaller_than_constant_hazard=hits_const / n_sim,
                min_attainable_p=1 / 24, sigma=sigma, n_sim=n_sim)


if __name__ == "__main__":
    pred = predictive_all()
    pw = power(pred)
    pw["level_dispersion"] = level_dispersion_power()
    pw["stamp"] = C.stamp()
    C.save_json("02_power.json", pw)
    import json
    print(json.dumps(pw, indent=1))
