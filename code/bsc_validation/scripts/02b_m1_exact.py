"""Pre-registered check of the M1 closed form against the verified exact simulator
(pre-registration 3.2).  Reads calibration files, geometry and event totals only.

For each hold-out shape event and for 24 posterior draws (D, seating
configuration, index seat):
  * closed form: eps profiled so that the expected total equals K;
  * exact simulator (per-cell rule beta q n_I/n_total, only the index infectious,
    no recovery): beta tuned by bisection so that the mean total equals K, then
    a production run; stratum attack rates and the distribution of k_near among
    replicates whose total equals K.

Output: data/02b_m1_exact.json
"""
import numpy as np

import _common as C
from valmod import events as E
from valmod import lattice as L
from valmod import predict as P
from valmod import stats as S
from valmod.exact import simulate_m1

N_DRAW = 24
N_TUNE = 1500
N_PROD = 8000


def sim_event(cfg, beta, D, n_rep, seed):
    inf = None
    for s, T in enumerate(cfg["segments"]):
        _, a = simulate_m1(cfg["lat"], cfg["src"], cfg["rec"], beta, D, T, n_rep=n_rep,
                           seed=seed + 1000 * s)
        inf = a if inf is None else (inf | a)
    return inf


def tune_beta(cfg, D, K, eps_cf, seed):
    lo, hi = np.log(eps_cf) - 1.0, np.log(eps_cf) + 9.0
    tot_hi = sim_event(cfg, np.exp(hi), D, N_TUNE, seed).sum(axis=1).mean()
    if tot_hi < K:
        return np.exp(hi), False
    for it in range(9):
        mid = 0.5 * (lo + hi)
        tot = sim_event(cfg, np.exp(mid), D, N_TUNE, seed + 7 * it).sum(axis=1).mean()
        if tot < K:
            lo = mid
        else:
            hi = mid
    return float(np.exp(0.5 * (lo + hi))), True


def main():
    out = {"stamp": C.stamp(), "n_draw": N_DRAW, "events": {}}
    for ev in ["T1", "T2", "T5", "C1"]:
        cls = E.CLASS[ev]
        post = P.posterior_D(cls, "M1", n=400, seed=C.SEED + 7)
        cal = P.load_cal(cls)
        D_mode = cal["M1"]["D_hat"]
        rng = np.random.default_rng(C.SEED + 31)
        n1, n2 = E.STRATA_N[ev]
        K = E.TOTALS[ev]
        pmf_sim = np.zeros(n1 + 1)
        pmf_cf = np.zeros(n1 + 1)
        rows = []
        n_cond = 0
        for i in range(N_DRAW):
            cfg = E.BUILDERS[ev](rng)
            D = D_mode if i < 4 else float(post["D"][i])
            X = np.clip(L.exposure(cfg["lat"], "M1", D, 0.0, cfg["segments"], cfg["src"], cfg["rec"]),
                        1e-300, None)
            eps = S.profile_eps(X, K)
            p_cf = -np.expm1(-eps * X)
            near = cfg["near"]
            beta, ok = tune_beta(cfg, D, K, eps, seed=C.SEED + 100 * i)
            inf = sim_event(cfg, beta, D, N_PROD, seed=C.SEED + 100 * i + 55)
            p_sim = inf.mean(axis=0)
            tot = inf.sum(axis=1)
            sel = tot == K
            kn = inf[sel][:, near].sum(axis=1)
            pm = np.bincount(kn, minlength=n1 + 1).astype(float)
            n_cond += int(sel.sum())
            pmf_sim += pm                    # pooled counts: draws weighted by P(total = K | draw)
            pmf_cf += S.cond_pmf(p_cf[near], p_cf[~near], K)
            rows.append(dict(D=D, beta=beta, reached_K=bool(ok), eps_cf=eps,
                             mean_total_sim=float(tot.mean()),
                             near_cf=float(p_cf[near].mean()), far_cf=float(p_cf[~near].mean()),
                             near_sim=float(p_sim[near].mean()), far_sim=float(p_sim[~near].mean()),
                             n_cond=int(sel.sum()), p_total_ge_K=float((tot >= K).mean()),
                             total_q=np.quantile(tot, [0.5, 0.95, 0.999]).tolist(), max_total=int(tot.max())))
            print(ev, i, {k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v)
                          for k, v in rows[-1].items() if k != "total_q"}, flush=True)
            out["events"][ev] = dict(rows=rows, pmf_sim=(pmf_sim / max(pmf_sim.sum(), 1)).tolist(),
                                     pmf_cf=(pmf_cf / (i + 1)).tolist(), n_cond=n_cond,
                                     n_near=n1, n_far=n2, K=K)
            C.save_json("02b_m1_exact.json", out)
        r = out["events"][ev]
        d_near = max(abs(x["near_cf"] - x["near_sim"]) for x in rows)
        d_far = max(abs(x["far_cf"] - x["far_sim"]) for x in rows)
        r["max_abs_diff_stratum_AR"] = max(d_near, d_far)
        r["replace_closed_form"] = bool(r["max_abs_diff_stratum_AR"] > 0.03)
        C.save_json("02b_m1_exact.json", out)
        print(ev, "max |diff|", r["max_abs_diff_stratum_AR"], "replace:", r["replace_closed_form"], flush=True)


if __name__ == "__main__":
    main()
