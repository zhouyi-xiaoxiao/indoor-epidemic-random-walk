"""Step 2: predictive distributions of every model for every event.  Reads kernels (step 1), strata sizes
and totals.  Does NOT read any stratum outcome (heldout_outcomes.json is never opened; the seat-level case
labels inside a3lib are not used here).  Writes results/02_predictions.json and results/02_o1_z.npz."""
import sys
import numpy as np, pandas as pd
import a3lib as A
from a3lib import E, L, S

t0 = A.time.time()
ker = A.load_json("01_kernels.json")
rng = np.random.default_rng(A.SEED + 2)
DOSES = ("exp", "lin", "gam")
out = {"stamp": A.stamp(), "events": {}}

def mix(pmfs):
    return np.mean(np.array(pmfs), axis=0)

def base_models(n1, n2, K):
    return {"M0": A.m0_pmf(n1, n2, K), "CRR2": A.crr_pmf(n1, n2, K, 2.0), "CRR3": A.crr_pmf(n1, n2, K, 3.0)}

def pack(n1, n2, K, pm, extra=None):
    d = dict(n_near=int(n1), n_far=int(n2), K=int(K), pmf={k: np.asarray(v).tolist() for k, v in pm.items()},
             rr_expected={k: A.expected_rr(np.asarray(v), n1, n2, K) for k, v in pm.items()})
    if extra: d.update(extra)
    return d

# ============================================================ T2 Hunan coach
ev = A.t2_event(); sus = ev["sus"]; lat = ev["lat"]; near = ev["near"]; K = ev["K"]
n1, n2 = int(near.sum()), int((~near).sum())
X0 = sus["s4_cfd"].to_numpy()
sig = ker["coach"]["sigma_log_meas_vs_cfd"]
noise = np.exp(sig * np.random.default_rng(A.SEED + 21).standard_normal((300, len(X0))))
pm = base_models(n1, n2, K)
p_seat = {}
for dose in DOSES:
    tag = "" if dose == "exp" else "_" + dose
    pm["P" + tag] = mix([A.near_pmf(X0 * z, near, K, dose)[0] for z in noise])
    pm["KA_nonoise" + tag] = A.near_pmf(X0, near, K, dose)[0]
p_seat["KA_nonoise"] = A.probs(X0, K)[0]
p_seat["KA_noise_draws"] = np.array([A.probs(X0 * z, K)[0] for z in noise])
# K-B with D from the measured tracer points
Db = np.array(ker["coach"]["boot_draws"])
kap = E.draw_kappa("T2", np.random.default_rng(A.SEED + 22), 200)
XB = [L.exposure(lat, "M2", float(Db[i]), float(kap[i]), [200 / 60], ev["src"], sus["site"].to_numpy()) for i in range(200)]
pm["KB_measD"] = mix([A.near_pmf(x, near, K)[0] for x in XB])
p_seat["KB_draws"] = np.array([A.probs(x, K)[0] for x in XB])
XBc = [L.exposure(lat, "M2", ker["coach"]["D_hat_cfd_S4A"], float(kap[i]), [200 / 60], ev["src"], sus["site"].to_numpy()) for i in range(200)]
pm["KB_cfdD"] = mix([A.near_pmf(x, near, K)[0] for x in XBc])
# measured points only: nearest measured point (physical distance on the lattice)
df = ev["df"]; m = df[df["s3_measured"].notna()]
mx, my = lat.xy_m(m["site"].to_numpy()); sx, sy = lat.xy_m(sus["site"].to_numpy())
nn = np.argmin(np.hypot(sx[:, None] - mx[None, :], sy[:, None] - my[None, :]), axis=1)
Xm = m["s3_measured"].to_numpy()[nn]
pm["KA_measured_points"] = A.near_pmf(Xm, near, K)[0]
p_seat["KA_measured_points"] = A.probs(Xm, K)[0]
# first validation's calibrated kernel for reference (posterior mode D = 25 m2/h)
pm["M2cal_mode"] = mix([A.near_pmf(L.exposure(lat, "M2", 25.1, float(k_), [200 / 60], ev["src"], sus["site"].to_numpy()), near, K)[0] for k_ in kap[:50]])
out["events"]["T2"] = pack(n1, n2, K, pm, dict(kernel="K-A Fig. S4A CFD + lognormal noise", sigma_log=sig))
# secondary split: driver side vs opposite side (13E excluded)
v = ev["side_valid"]; ds = ev["driver_side"][v]
n1s, n2s = int(ds.sum()), int((~ds).sum())
pms = base_models(n1s, n2s, K)
pms["P"] = mix([A.near_pmf((X0 * z)[v], ds, K)[0] for z in noise])
pms["KA_nonoise"] = A.near_pmf(X0[v], ds, K)[0]
pms["KB_measD"] = mix([A.near_pmf(x[v], ds, K)[0] for x in XB])
out["events"]["T2_side"] = pack(n1s, n2s, K, pms)
np.savez(A.RES + "/02_t2_seatprobs.npz", **p_seat)
print("T2 done", round(A.time.time() - t0, 1), flush=True)

# ============================================================ T3 minibus (seat level only)
b2 = pd.read_csv(A.DER + "/ou2022_B2_seats.csv")
s3 = b2[b2.s4_cfd.notna() & ~b2.seat.isin(["D", "C"])].reset_index(drop=True)
sig3 = ker["minibus"]["sigma_log_passenger_seats"]
nz3 = np.exp(sig3 * np.random.default_rng(A.SEED + 23).standard_normal((300, len(s3))))
np.savez(A.RES + "/02_t3_seatprobs.npz", seats=s3.seat.to_numpy().astype(str), KA_nonoise=A.probs(s3.s4_cfd.to_numpy(), 2)[0],
         KA_noise_draws=np.array([A.probs(s3.s4_cfd.to_numpy() * z, 2)[0] for z in nz3]))
out["events"]["T3"] = dict(n=int(len(s3)), K=2, seats=s3.seat.tolist(), s4=s3.s4_cfd.tolist())

# ============================================================ T1 Zhejiang bus (K-B, D from the coach)
r1 = np.random.default_rng(A.SEED + 24)
pmP = {d: [] for d in DOSES}; pmC = []; pmcal = []
for i in range(300):
    cfg = E.build_T1(r1)
    X = L.exposure(cfg["lat"], "M2", float(Db[i % 200]), cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"])
    for d in DOSES:
        pmP[d].append(A.near_pmf(X, cfg["near"], cfg["K"], d)[0])
    Xc = L.exposure(cfg["lat"], "M2", ker["coach"]["D_hat_cfd_S4A"], cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"])
    pmC.append(A.near_pmf(Xc, cfg["near"], cfg["K"])[0])
    Xk = L.exposure(cfg["lat"], "M2", 25.1, cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"])
    pmcal.append(A.near_pmf(Xk, cfg["near"], cfg["K"])[0])
pm = base_models(33, 34, 23)
pm.update({"P": mix(pmP["exp"]), "P_lin": mix(pmP["lin"]), "P_gam": mix(pmP["gam"]), "KB_cfdD": mix(pmC), "M2cal_mode": mix(pmcal)})
out["events"]["T1"] = pack(33, 34, 23, pm, dict(kernel="K-B, D from coach tracer (bootstrap)"))
print("T1 done", round(A.time.time() - t0, 1), flush=True)

# ============================================================ T5 flight VN54 (K-A, Kinahan forward cabin)
ev5, kern = A.t5_direct_kernels(False)
_, kern_zm = A.t5_direct_kernels(True)
near5 = ev5["near"]; K5 = ev5["K"]; n1, n2 = int(near5.sum()), int((~near5).sum())
pm = base_models(n1, n2, K5)
p5 = {}
for d in DOSES:
    tag = "" if d == "exp" else "_" + d
    pm["P" + tag] = mix([A.near_pmf(v[0], near5, K5, d)[0] for v in kern.values()])
for k_, v in kern.items():
    pm["KA_" + k_] = A.near_pmf(v[0], near5, K5)[0]
    p5["KA_" + k_] = A.probs(v[0], K5)[0]
pm["KA_zeros_missing"] = mix([A.near_pmf(v[0], near5, K5)[0] for v in kern_zm.values()])
# K-B
r5 = np.random.default_rng(A.SEED + 25)
D2 = [ker["cabin"]["zeros_kept"]["fits"][r]["D_hat"] for r in ("5L", "5A")]
kap5 = r5.uniform(30, 36, 100) + A.LOSS
kap5old = np.exp(r5.uniform(np.log(10), np.log(30), 100)) + A.LOSS
def kb5(Dv, kv):
    return A.near_pmf(L.exposure(ev5["lat"], "M2", Dv, kv, [10.0], ev5["src"], ev5["rec"]), near5, K5)[0]
pm["KB_measD"] = mix([kb5(D2[i % 2], kap5[i]) for i in range(100)])
pm["KB_measD_oldvent"] = mix([kb5(D2[i % 2], kap5old[i]) for i in range(100)])
pm["M2cal_mode"] = mix([kb5(25.1, kap5old[i]) for i in range(100)])
p5["KB_draws"] = np.array([A.probs(L.exposure(ev5["lat"], "M2", D2[i % 2], kap5[i], [10.0], ev5["src"], ev5["rec"]), K5)[0] for i in range(100)])
# S-2A: far stratum without the passenger lost to follow-up
keep = np.array([s != "2A" for s in ev5["seats"]])
pm2a = base_models(12, 7, K5)
pm2a["P"] = mix([A.near_pmf(v[0][keep], near5[keep], K5)[0] for v in kern.values()])
out["events"]["T5"] = pack(n1, n2, K5, pm, dict(kernel="K-A, 777 forward cabin releases 5L and 5A (mirrored)",
                                                 n_imputed={k_: v[1] for k_, v in kern.items()}))
out["events"]["T5_no2A"] = pack(12, 7, K5, pm2a)
# post-hoc S-5G: centre-seat release used unshifted (index treated as sitting in 5G)
mp = A.kinahan_maps(False)["5G"]
vals = []
for r, c in zip(ev5["row"], ev5["col"]):
    cc = "L" if c == "K" else c
    v_ = mp.loc[r, cc] if (r in mp.index and cc in mp.columns) else np.nan
    if not np.isfinite(v_):
        v_ = np.nanmean(mp.loc[r].to_numpy())
    vals.append(v_)
out["events"]["T5"]["pmf"]["KA_5G_posthoc"] = A.near_pmf(np.array(vals), near5, K5)[0].tolist()
np.savez(A.RES + "/02_t5_seatprobs.npz", seats=np.array(ev5["seats"]), **p5)
print("T5 done", round(A.time.time() - t0, 1), flush=True)

# ============================================================ C1 classroom (K-B, Cheng relation)
rc = np.random.default_rng(A.SEED + 26)
pmP = {k_: [] for k_ in ("exp", "lin", "gam", "x2", "d2")}; pmcal = []; Dc = []
for i in range(300):
    cfg = E.build_C1(rc)
    lat1 = cfg["lat"]; V = lat1.M * lat1.cell_volume
    ach = cfg["kappa"] - A.LOSS
    Dr = (0.52 * ach + 0.31) * V ** (2 / 3) * np.exp(0.42 * rc.standard_normal())
    Dc.append(Dr)
    X = L.exposure(lat1, "M2", Dr, cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"])
    for d in DOSES:
        pmP[d].append(A.near_pmf(X, cfg["near"], cfg["K"], d)[0])
    for tag, f in (("x2", 2.0), ("d2", 0.5)):
        pmP[tag].append(A.near_pmf(L.exposure(lat1, "M2", Dr * f, cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"]), cfg["near"], cfg["K"])[0])
    pmcal.append(A.near_pmf(L.exposure(lat1, "M2", 0.5, cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"]), cfg["near"], cfg["K"])[0])
pm = base_models(10, 14, 12)
pm.update({"P": mix(pmP["exp"]), "P_lin": mix(pmP["lin"]), "P_gam": mix(pmP["gam"]), "KB_Dx2": mix(pmP["x2"]), "KB_Dhalf": mix(pmP["d2"]), "M2cal_mode": mix(pmcal)})
out["events"]["C1"] = pack(10, 14, 12, pm, dict(kernel="K-B, D from Cheng relation", D_q=np.quantile(Dc, [0.05, 0.5, 0.95]).tolist()))
print("C1 done", round(A.time.time() - t0, 1), flush=True)

# ============================================================ R1 restaurant (K-A, Li Table 3)
tb = A.r1_tables()
nearT = tb["neighbour_class"].eq("immediate").to_numpy()
npat = tb["patrons"].to_numpy()
Xt = (tb["tracer"] * tb["T_h"]).to_numpy()
Xe = tb["norm_predicted_exposure"].to_numpy()
Xp = np.repeat(Xt, npat); Xpe = np.repeat(Xe, npat); nearP = np.repeat(nearT, npat)
n1, n2 = int(nearP.sum()), int((~nearP).sum())
assert (n1, n2) == (16, 63)
r1e = {}
for K in (2, 3, 4, 5):
    pm = base_models(n1, n2, K)
    for d in DOSES:
        pm["P" if d == "exp" else "P_" + d] = A.near_pmf(Xp, nearP, K, d)[0]
    pm["KA_cfd_aerosol_exposure"] = A.near_pmf(Xpe, nearP, K)[0]
    r1e[K] = pack(n1, n2, K, pm)
out["events"]["R1"] = r1e[2]; out["events"]["R1_byK"] = {str(k_): v for k_, v in r1e.items()}
out["events"]["R1"]["tracer_times_overlap"] = dict(zip(tb["table"], np.round(Xt, 3).tolist()))

# ============================================================ T4 trains (K-B, D from the rail carriage)
lat4 = E.t4_lattice(1.0); sites4, rows4, cols4 = E.t4_seats(lat4)
d4 = E.t4_data()                         # cell sizes n (and counts k, which are NOT used in this script)
dr_ = np.abs(rows4[:, None] - rows4[None, :]); dc_ = np.abs(cols4[:, None] - cols4[None, :])
cell_pairs = [np.flatnonzero(((dr_ == r) & (dc_ == c)).ravel()) for r, c in zip(d4["dr"], d4["dc"])]
mask = ~((d4["dr"] == 0) & (d4["dc"] == 1)).to_numpy()
n4 = d4["n"].to_numpy().astype(float)
Tn, Tw = E.t4_duration_nodes(16); kq, kw = E.kappa_quadrature("T4", 7)
allidx = np.concatenate(cell_pairs); lab = np.concatenate([np.full(len(ix), c) for c, ix in enumerate(cell_pairs)])
cnt = np.bincount(lab, minlength=len(cell_pairs)).astype(float)
K4 = None  # total is read in step 4; conditional cell probabilities do not depend on it in the rare-event regime
def t4_pi(D, kappa, total=138.0):
    Xall = np.stack([np.clip(L.exposure_matrix(lat4, "M2", D, kappa, [T], sites4).ravel()[allidx], 0, None) for T in Tn])
    def prob(eps):
        q = Tw @ (-np.expm1(-eps * Xall)); return np.bincount(lab, weights=q, minlength=len(cell_pairs)) / cnt
    f = lambda le: np.sum((n4 * prob(np.exp(le)))[mask]) - total
    lo, hi = -5.0, 5.0
    while f(lo) > 0: lo -= 3
    while f(hi) < 0: hi += 3
    eps = np.exp(A.optimize.brentq(f, lo, hi, xtol=1e-10))
    p = prob(eps); w = (n4 * p)[mask]
    return w / w.sum(), p[mask], eps
t4 = {"cells": dict(dr=d4["dr"][mask].tolist(), dc=d4["dc"][mask].tolist(), n=n4[mask].tolist())}
for name, panel in (("P", "B_middle"), ("KB_end_release", "A_end")):
    bd = np.array(ker["train"][panel]["boot_draws"]); Dq = np.quantile(bd, (np.arange(20) + 0.5) / 20)
    pis = np.array([[t4_pi(float(D_), float(k_))[0] for k_ in kq] for D_ in Dq])     # (20, 7, 22)
    t4[name] = dict(pi_draws=pis.reshape(-1, pis.shape[-1]).tolist(), D_draws=Dq.tolist())
    print("T4", name, round(A.time.time() - t0, 1), flush=True)
t4["M2cal_mode"] = dict(pi_draws=[t4_pi(25.1, float(k_))[0].tolist() for k_ in kq])
t4["M0"] = dict(pi=(n4[mask] / n4[mask].sum()).tolist())
for rho in (2.0, 3.0):
    w = n4[mask] * np.where(d4["dr"][mask].to_numpy() == 0, rho, 1.0)
    t4["CRR%d" % rho] = dict(pi=(w / w.sum()).tolist())
out["events"]["T4"] = t4
A.save_json("02_predictions.json", out)
print("T4 done", round(A.time.time() - t0, 1), flush=True)

# ============================================================ O1 call centre (K-B, Cheng relation)
xy, _case_not_used = E.o1_desks()
Aadj = E.o1_adjacency(xy, 1.25)
R_TUNE, R_PROD, N_PAR = 400, 1500, 60
def percolate(Xe, rng, R):
    n = Xe.shape[0]; inf = np.zeros((R, n), bool); new = np.zeros((R, n), bool)
    new[np.arange(R), rng.integers(n, size=R)] = True; inf |= new
    while new.any():
        p = -np.expm1(-(new.astype(float) @ Xe)); nxt = (rng.random((R, n)) < p) & ~inf; inf |= nxt; new = nxt
    return inf
def tune_eps(X, rng, target=79.0, thresh=40):
    c = np.log(1.0 / X.sum(axis=1).mean()); lo, hi = c - 2.0, c + 12.0
    for it in range(14):
        mid = 0.5 * (lo + hi); fs = percolate(np.exp(mid) * X, rng, R_TUNE).sum(axis=1); big = fs[fs >= thresh]
        if (big.mean() if len(big) >= 5 else 0.0) < target: lo = mid
        else: hi = mid
    return float(np.exp(0.5 * (lo + hi)))
def kern(cfg, D):
    lat = cfg["lat"]; w = L.steady_weights(lat, "M2", D, cfg["kappa"]); Ph = lat.Phi[cfg["sites"]]
    X = np.clip((Ph * w) @ Ph.T, 0, None); np.fill_diagonal(X, 0.0); return X
ro = np.random.default_rng(A.SEED + 27)
zs = {"P": [], "M0": [], "KB_Dx2": [], "KB_Dhalf": []}; rowsO = []
for i in range(N_PAR):
    cfg = E.build_O1(ro); lat = cfg["lat"]; V = lat.M * lat.cell_volume; ach = cfg["kappa"] - A.LOSS
    Dr = (0.52 * ach + 0.31) * V ** (2 / 3) * np.exp(0.42 * ro.standard_normal())
    for tag, Dv in (("P", Dr), ("KB_Dx2", 2 * Dr), ("KB_Dhalf", 0.5 * Dr), ("M0", 1e9)):
        if tag == "M0" and i >= 20: continue
        X = kern(cfg, Dv); eps = tune_eps(X, ro)
        inf = percolate(eps * X, ro, R_PROD); fs = inf.sum(axis=1); keep = (fs >= 69) & (fs <= 89)
        z = np.array([S.joincount_z(Aadj, c)[0] for c in inf[keep]])
        zs[tag].append(z)
        if tag == "P": rowsO.append(dict(D=float(Dr), ach=float(ach), pitch=float(cfg["pitch"]), V=float(V), n_keep=int(keep.sum()), z_mean=float(z.mean()) if len(z) else None))
    if i % 10 == 9: print("O1", i + 1, round(A.time.time() - t0, 1), flush=True)
np.savez(A.RES + "/02_o1_z.npz", **{k_: np.concatenate(v) for k_, v in zs.items()},
         **{k_ + "_draw": np.concatenate([np.full(len(z), j) for j, z in enumerate(v)]) for k_, v in zs.items()})
out["events"]["O1"] = dict(rows=rowsO, q={k_: np.quantile(np.concatenate(v), [0.025, 0.05, 0.5, 0.95, 0.975]).tolist() for k_, v in zs.items()},
                           n={k_: int(sum(len(z) for z in v)) for k_, v in zs.items()})
A.save_json("02_predictions.json", out)
print("all done", round(A.time.time() - t0, 1))
for e_, d_ in out["events"].items():
    if "rr_expected" in d_: print(e_, d_["n_near"], d_["n_far"], d_["K"], {k_: round(v, 2) for k_, v in d_["rr_expected"].items()})
print("O1 quantiles", out["events"]["O1"]["q"])
