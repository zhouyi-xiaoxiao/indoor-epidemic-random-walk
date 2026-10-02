"""Descriptive (not part of the decision): physics D against the D each event prefers when fitted alone,
mixing lengths, implied emission, near/far tracer contrast.  Writes results/05_descriptive.json."""
import numpy as np, pandas as pd
import a3lib as A
from a3lib import E, L, S
ker = A.load_json("01_kernels.json"); ev = A.load_json("04_evaluation.json"); out = {}
pref = {"T5 flight": 0.1, "T4 trains": 25.1, "T1 bus": 63.1, "T2 coach": 6.31e3, "R1 restaurant": 0.398}   # first validation, Table V16
out["D_table"] = [
    dict(enclosure="coach (Ou 2022, ethane, 8 points)", D_phys=ker["coach"]["D_hat"], lo=ker["coach"]["boot"]["q05"], hi=ker["coach"]["boot"]["q95"], kappa_tracer=6.5,
         D_outbreak_alone=pref["T2 coach"], D_calibrated=25.1, event="T2 coach"),
    dict(enclosure="wide-body forward cabin (Kinahan 2021, 1 um particles, 2 releases)", D_phys=ker["cabin"]["zeros_kept"]["D_geomean"],
         lo=ker["cabin"]["zeros_kept"]["fits"]["5A"]["D_hat"], hi=ker["cabin"]["zeros_kept"]["fits"]["5L"]["D_hat"], kappa_tracer=35.0,
         D_outbreak_alone=pref["T5 flight"], D_calibrated=25.1, event="T5 flight"),
    dict(enclosure="rail carriage (Woodward 2022, salt aerosol, middle release)", D_phys=ker["train"]["B_middle"]["D_hat"], lo=ker["train"]["B_middle"]["boot"]["q05"],
         hi=ker["train"]["B_middle"]["boot"]["q95"], kappa_tracer=13.0, D_outbreak_alone=pref["T4 trains"], D_calibrated=25.1, event="T4 trains"),
    dict(enclosure="rail carriage, end release", D_phys=ker["train"]["A_end"]["D_hat"], lo=ker["train"]["A_end"]["boot"]["q05"], hi=ker["train"]["A_end"]["boot"]["q95"],
         kappa_tracer=13.0, D_outbreak_alone=pref["T4 trains"], D_calibrated=25.1, event="T4 trains"),
    dict(enclosure="classroom 429 m3 (Cheng 2011 relation, 2-12 ACH)", D_phys=float(np.sqrt(76.8 * 372.6)), lo=76.8 / 2, hi=372.6 * 2, kappa_tracer=None,
         D_outbreak_alone=None, D_calibrated=0.5, event="C1 classroom"),
    dict(enclosure="open-plan office wing (Cheng 2011 relation, 1-6 ACH)", D_phys=float(np.sqrt(91.1 * 376.5)), lo=91.1 / 2, hi=376.5 * 2, kappa_tracer=None,
         D_outbreak_alone=None, D_calibrated=0.5, event="O1 call centre"),
]
for r in out["D_table"]:
    if r["kappa_tracer"]:
        r["mixing_length_m"] = float(np.sqrt(r["D_phys"] / r["kappa_tracer"])); r["D_phys_m2_per_s"] = r["D_phys"] / 3600
# restaurant: tracer D by the same estimator is not possible without table coordinates; report the contrast instead
tb = A.r1_tables()
out["R1_tracer"] = dict(near_mean=float(tb[tb.neighbour_class == "immediate"].tracer.mean()), remote_mean=float(tb[tb.neighbour_class == "remote"].tracer.mean()),
                        remote_min=float(tb[tb.neighbour_class == "remote"].tracer.min()), remote_max=float(tb[tb.neighbour_class == "remote"].tracer.max()),
                        aerosol_near_mean=float(tb[tb.neighbour_class == "immediate"].norm_predicted_exposure.mean()),
                        aerosol_remote_mean=float(tb[tb.neighbour_class == "remote"].norm_predicted_exposure.mean()))
# dose ratio needed: what near/far exposure ratio would the observed counts need under the exponential law?
need = {}
for e, (k1, n1, k2, n2) in {"T2": (3, 19, 4, 26), "T1": (14, 33, 9, 34), "T5": (11, 12, 1, 8), "C1": (8, 10, 4, 14)}.items():
    need[e] = dict(dose_ratio_needed=float(np.log(1 - k1 / n1) / np.log(1 - k2 / n2)), risk_ratio_obs=float((k1 / n1) / (k2 / n2)),
                   rr_predicted_P=ev["events"][e]["models"]["P"]["rr_expected"])
out["dose_ratio_needed"] = need
# implied emission for K-B events (E = eps * cell volume / breathing rate), medians over draws
rng = np.random.default_rng(A.SEED + 5); Db = np.array(ker["coach"]["boot_draws"]); em = {}
v = []
for i in range(200):
    cfg = E.build_T1(rng); X = L.exposure(cfg["lat"], "M2", float(Db[i]), cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"])
    v.append(S.profile_eps(X, 23) * cfg["lat"].cell_volume / E.BREATH)
em["T1 bus"] = np.quantile(v, [0.05, 0.5, 0.95]).tolist()
t2 = A.t2_event(); kap = E.draw_kappa("T2", rng, 200)
v = [S.profile_eps(L.exposure(t2["lat"], "M2", float(Db[i]), float(kap[i]), [200 / 60], t2["src"], t2["sus"]["site"].to_numpy()), 7) * t2["lat"].cell_volume / E.BREATH for i in range(200)]
em["T2 coach"] = np.quantile(v, [0.05, 0.5, 0.95]).tolist()
t5 = A.t5_event(); D2 = [ker["cabin"]["zeros_kept"]["fits"][r]["D_hat"] for r in ("5L", "5A")]
v = [S.profile_eps(L.exposure(t5["lat"], "M2", D2[i % 2], rng.uniform(30, 36) + A.LOSS, [10.0], t5["src"], t5["rec"]), 12) * t5["lat"].cell_volume / E.BREATH for i in range(100)]
em["T5 flight (K-B)"] = np.quantile(v, [0.05, 0.5, 0.95]).tolist()
v = []
for i in range(200):
    cfg = E.build_C1(rng); lat = cfg["lat"]; ach = cfg["kappa"] - A.LOSS
    Dr = (0.52 * ach + 0.31) * (lat.M * lat.cell_volume) ** (2 / 3) * np.exp(0.42 * rng.standard_normal())
    v.append(S.profile_eps(L.exposure(lat, "M2", Dr, cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"]), 12) * lat.cell_volume / E.BREATH)
em["C1 classroom"] = np.quantile(v, [0.05, 0.5, 0.95]).tolist()
out["implied_emission_quanta_per_h"] = em
A.save_json("05_descriptive.json", out)
import json; print(json.dumps(out, indent=1))
