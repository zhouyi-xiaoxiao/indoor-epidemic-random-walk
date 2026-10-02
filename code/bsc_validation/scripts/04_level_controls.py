"""Level checks and negative controls (pre-registration sections 5.3, 5.4, 7 and
addendum A1.5).

1  implied emission rate per event (Wells-Riley scaling shared by M0 and M2),
   implied hazard per hour, implied M1 pair hazard; dispersion test
2  tied-index check T2 -> T3
3  unselected-cohort transfer T4 -> C3
4  negative controls S1, C3, T5 premium economy, T7
"""
import itertools

import numpy as np
from scipy import stats

import _common as C
from valmod import events as E
from valmod import lattice as L
from valmod import level as LV
from valmod import predict as P
from valmod import stats as S
from valmod.exact import simulate_m1

rng = np.random.default_rng(C.SEED + 41)
obs = C.load_json("heldout_outcomes.json")
out = {"stamp": C.stamp()}
SIGMA = 0.94
N = 20000

# ---------------------------------------------------------------------------
# 1 implied infectiousness and dispersion
# ---------------------------------------------------------------------------
lev = {}
for ev in ["H1", "O2", "R2", "C5", "T2", "T3"]:
    d = LV.LEVEL_EVENTS[ev]
    o = obs["level"][ev]
    kap = E.draw_kappa(ev, rng, 4000)
    s = LV.wr_factor(d["V"], kap, d["segments"])
    lo, hi = S.clopper_pearson(o["k"], o["n"])
    dose, dlo, dhi = [LV.implied_dose(x * o["n"], o["n"]) for x in (o["k"] / o["n"], max(lo, 1e-9), min(hi, 1 - 1e-9))]
    s_med = float(np.median(s))
    T = sum(d["segments"])
    lev[ev] = dict(label=d["label"], k=o["k"], n=o["n"], T_h=T, V=d["V"], kappa_median=float(np.median(kap)),
                   wr_factor_median=s_med, wr_factor_90=np.quantile(s, [0.05, 0.95]).tolist(),
                   dose=float(dose), dose_ci=[float(dlo), float(dhi)],
                   E_implied=float(dose / s_med),
                   E_ci_binomial=[float(dlo / s_med), float(dhi / s_med)],
                   E_range_all=[float(dlo / np.quantile(s, 0.95)), float(dhi / np.quantile(s, 0.05))],
                   hazard_per_h=float(dose / T),
                   m1_scaling=float(T / d["area"]), m1_implied=float(dose / (T / d["area"])))
    if "k_alt" in o:
        da = LV.implied_dose(o["k_alt"], o.get("n_alt", o["n"]))
        lev[ev]["E_implied_alt"] = float(da / s_med)
out["implied"] = lev

evs = ["H1", "O2", "R2", "C5"]
dose = np.array([lev[e]["dose"] for e in evs])
fac = {"wells_riley": np.array([lev[e]["wr_factor_median"] for e in evs]),
       "m1_T_over_area": np.array([lev[e]["m1_scaling"] for e in evs]),
       "constant_per_hour": np.array([lev[e]["T_h"] for e in evs])}
perms = list(itertools.permutations(range(4)))
disp = {}
for name, f in fac.items():
    sd = float(np.std(np.log(dose / f)))
    sdp = np.array([np.std(np.log(dose / f[list(p)])) for p in perms])
    disp[name] = dict(sd_log=sd, gsd=float(np.exp(sd)), perm_p=float(np.mean(sdp <= sd + 1e-12)))
out["dispersion"] = disp
# continuity with the original protocol: slope of log observed dose on log predicted (common E)
x = np.log(fac["wells_riley"])
y = np.log(dose)
slope, icpt, r, p, se = stats.linregress(x, y)
out["original_slope"] = dict(slope=float(slope), se=float(se), r=float(r), n=4,
                             note="log observed cumulative hazard on log Wells-Riley factor, four events")

# ---------------------------------------------------------------------------
# 2 tied-index check T2 -> T3
# ---------------------------------------------------------------------------
tied = {}
leps = np.linspace(-12, 12, 481)
for model in ("M0", "M1", "M2"):
    post = None if model == "M0" else P.posterior_D("vehicle", model, n=400, seed=C.SEED + 7)
    pm3 = np.zeros(18)
    pm3_alt = np.zeros(13)
    for i in range(200):
        c2 = E.build_T2(rng)
        c3 = E.build_T3(rng)
        D = 1.0 if model == "M0" else float(post["D"][i])
        X2 = np.clip(L.exposure(c2["lat"], model, D, c2["kappa"], c2["segments"], c2["src"], c2["rec"]), 1e-300, None)
        X3 = np.clip(L.exposure(c3["lat"], model, D, c3["kappa"], c3["segments"], c3["src"], c3["rec"]), 1e-300, None)
        if model == "M1":
            # same pair hazard, corrected for the occupancy term E[1/N] of the per-cell rule
            pass
        # likelihood of the T2 total (7 of 45 passengers; the driver adds one exposed, 0 cases)
        ll = np.array([np.log(max(S.pb_pmf(-np.expm1(-np.exp(le) * X2))[7], 1e-300)) for le in leps[::4]])
        w = np.exp(ll - ll.max())
        w /= w.sum()
        for le, wi in zip(leps[::4], w):
            if wi < 1e-6:
                continue
            p3 = -np.expm1(-np.exp(le) * X3)
            pm3 += wi * S.pb_pmf(p3)
            pm3_alt += wi * S.pb_pmf(p3[rng.choice(17, 12, replace=False)])
    pm3 /= pm3.sum()
    pm3_alt /= pm3_alt.sum()
    lo, hi = S.pmf_interval(pm3, 0.90)
    lo_a, hi_a = S.pmf_interval(pm3_alt, 0.90)
    k3 = obs["level"]["T3"]["k"]
    tied[model] = dict(pmf=pm3.tolist(), mean=float((np.arange(18) * pm3).sum()), pi90=[lo, hi],
                       inside=bool(lo <= k3 <= hi), p_two_sided=S.two_sided_p(pm3, k3),
                       alt_n12=dict(mean=float((np.arange(13) * pm3_alt).sum()), pi90=[lo_a, hi_a],
                                    inside=bool(lo_a <= k3 <= hi_a)))
    print("tied", model, tied[model]["mean"], tied[model]["pi90"], tied[model]["inside"], flush=True)
out["tied_T2_T3"] = tied

# ---------------------------------------------------------------------------
# 3-4 T4-anchored infectiousness; C3, S1, T7
# ---------------------------------------------------------------------------
cal = P.load_cal("vehicle")
post2 = P.posterior_D("vehicle", "M2", n=400, seed=C.SEED + 7)
E_mean_T4 = post2["eps_cal"] * cal["M2"]["cell_volume"] / E.BREATH          # quanta/h, mean over index cases
out["E_T4"] = dict(q=np.quantile(E_mean_T4, [0.05, 0.5, 0.95]).tolist(),
                   note="mean emission rate over train index cases (M2 calibration, kappa prior 5-20 ACH + 0.93)")
post1 = P.posterior_D("vehicle", "M1", n=400, seed=C.SEED + 7)
h_T4 = post1["eps_cal"]                                                    # pair hazard per hour co-located
out["h_T4_M1"] = dict(q=np.quantile(h_T4, [0.05, 0.5, 0.95]).tolist(), cell_area_m2=0.5)


def c3_sim(mask_factor, model="M2"):
    """Predictive distribution of the pooled number of cases among 728 contacts of 51 index cases."""
    NC = 4000
    ks = np.zeros(NC, int)
    n_i = rng.multinomial(728, np.ones(51) / 51, size=NC)
    for a in range(NC):
        V = np.exp(rng.uniform(np.log(180), np.log(350), 51))
        kap = E.draw_kappa("C3", rng, 51)
        T = rng.uniform(6.5, 13.0, 51)
        nd = np.maximum(np.rint(T / 6.5), 1)
        if model == "M2":
            Em = E_mean_T4[rng.integers(len(E_mean_T4))]
            Ei = np.exp(rng.normal(LV.lognormal_from_mean(Em, SIGMA), SIGMA, 51))
            s = E.BREATH / V * (T / kap - nd * (-np.expm1(-kap * T / nd)) / kap ** 2)
            dose = mask_factor * Ei * s
        else:
            hm = h_T4[rng.integers(len(h_T4))]
            hi_ = np.exp(rng.normal(LV.lognormal_from_mean(hm, SIGMA), SIGMA, 51))
            dose = mask_factor * hi_ * T / (V / 3.0)                        # T / M, M = area / (1 m)^2
        ks[a] = rng.binomial(n_i[a], -np.expm1(-dose)).sum()
    return ks


c3 = {}
k_obs = obs["level"]["C3"]["k"]
ub_c3 = S.upper_bound_one_sided(5, 728)
for model in ("M2", "M1"):
    for name, f in (("masks_0.35", 0.35), ("masks_0.2", 0.2), ("masks_0.6", 0.6), ("no_masks", 1.0)):
        ks = c3_sim(f, model)
        lo, hi = np.quantile(ks, [0.05, 0.95])
        c3[f"{model}_{name}"] = dict(mean=float(ks.mean()), median=float(np.median(ks)), pi90=[float(lo), float(hi)],
                                     inside=bool(lo <= k_obs <= hi), P_exceed_bound=float((ks / 728 > ub_c3).mean()),
                                     ar_mean_pct=float(100 * ks.mean() / 728))
        print("C3", model, name, c3[f"{model}_{name}"], flush=True)
out["C3"] = c3
out["C3"]["bound_pct"] = 100 * ub_c3

# S1 supermarket customers (maximal-exposure assumption)
ub_s1 = S.upper_bound_one_sided(0, 8224)
V, Tdw = 6000.0, 27 / 60
kap = E.draw_kappa("S1", rng, N)
Em = E_mean_T4[rng.integers(len(E_mean_T4), size=N)]
Ei = np.exp(rng.normal(LV.lognormal_from_mean(Em, SIGMA), SIGMA))
ar_single = -np.expm1(-Ei * E.BREATH * Tdw / (V * kap))              # one infectious employee, one draw of E
ar_pooled = Em * E.BREATH * Tdw / (V * kap)                          # averaged over many index-days
hm = h_T4[rng.integers(len(h_T4), size=N)]
ar_m1 = hm * Tdw / (1500.0 / 4.0)                                    # T / M, supermarket scene a = 2 m
out["S1"] = dict(bound_pct=100 * ub_s1,
                 M2_single_index=dict(median_pct=float(100 * np.median(ar_single)),
                                      q95_pct=float(100 * np.quantile(ar_single, 0.95)),
                                      P_exceed=float((ar_single > ub_s1).mean())),
                 M2_pooled=dict(median_pct=float(100 * np.median(ar_pooled)), P_exceed=float((ar_pooled > ub_s1).mean())),
                 M1_pooled=dict(median_pct=float(100 * np.median(ar_m1)), P_exceed=float((ar_m1 > ub_s1).mean())))
print("S1", out["S1"], flush=True)

# T7 metro ride
V7, T7 = 148.5, 20 / 60
kap = E.draw_kappa("T7", rng, N)
s7 = LV.wr_factor(V7, kap, [T7])
k7 = rng.binomial(309, -np.expm1(-Ei * s7))
k7_m1 = rng.binomial(309, -np.expm1(-np.exp(rng.normal(LV.lognormal_from_mean(hm, SIGMA), SIGMA)) * T7 / (45 * 6)))
out["T7"] = dict(M2=dict(mean_cases=float(k7.mean()), P_ge5=float((k7 >= 5).mean())),
                 M1=dict(mean_cases=float(k7_m1.mean()), P_ge5=float((k7_m1 >= 5).mean())))
print("T7", out["T7"], flush=True)

# T5 premium economy (0/35), M2 at the eps profiled from the business-class total
ub_pe = S.upper_bound_one_sided(0, 35)
kpe, ppe = [], []
for i in range(300):
    c5 = E.build_T5(rng)
    D = float(post2["D"][i])
    X = np.clip(L.exposure(c5["lat"], "M2", D, c5["kappa"], c5["segments"], c5["src"], c5["rec"]), 1e-300, None)
    eps = S.profile_eps(X, 12)
    Xpe = np.clip(L.exposure(c5["lat"], "M2", D, c5["kappa"], c5["segments"], c5["src"], c5["pe_sites"]), 0, None)
    p = float(np.mean(-np.expm1(-eps * Xpe)))
    ppe.append(p)
    kpe.append(rng.binomial(35, p))
kpe, ppe = np.array(kpe), np.array(ppe)
out["T5_PE"] = dict(bound_pct=100 * ub_pe, M2=dict(ar_median_pct=float(100 * np.median(ppe)),
                                                   ar_q95_pct=float(100 * np.quantile(ppe, 0.95)),
                                                   mean_cases=float(kpe.mean()),
                                                   P_exceed=float((kpe / 35 > ub_pe).mean()),
                                                   P_zero=float((kpe == 0).mean())))
# M1 by the exact simulator: business passengers plus 35 premium-economy passengers, beta at its cap
pe_rates = []
for i in range(6):
    c5 = E.build_T5(rng)
    pe_people = np.concatenate([c5["pe_sites"], c5["pe_sites"][:5]])
    rec = np.concatenate([c5["rec"], pe_people])
    p, _ = simulate_m1(c5["lat"], c5["src"], rec, 1e5, float(post1["D"][i]), 10.0, n_rep=4000, seed=C.SEED + 900 + i)
    pe_rates.append(float(p[20:].mean()))
out["T5_PE"]["M1_exact_beta_cap"] = dict(ar_mean_pct=float(100 * np.mean(pe_rates)),
                                         P_exceed=float(np.mean(np.array(pe_rates) > ub_pe)))
print("T5 PE", out["T5_PE"], flush=True)

C.save_json("04_level_controls.json", out)
print({k: (round(v["E_implied"], 1), [round(x, 1) for x in v["E_range_all"]]) for k, v in lev.items()})
print(disp)
