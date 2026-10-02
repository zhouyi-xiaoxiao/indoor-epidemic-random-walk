"""Step 1 (physics only): estimate the kernels.  Reads tracer data, never an outbreak count.
Writes results/01_kernels.json."""
import numpy as np, pandas as pd
import a3lib as A
from a3lib import E, L, S

rng = np.random.default_rng(A.SEED + 1)
out = {"stamp": A.stamp(), "logD_grid": [float(A.LOGD[0]), float(A.LOGD[-1]), 0.05]}

def summ(Db):
    q = np.quantile(np.log10(Db), [0.05, 0.5, 0.95])
    return dict(q05=float(10 ** q[0]), q50=float(10 ** q[1]), q95=float(10 ** q[2]),
                frac_at_upper_edge=float(np.mean(Db >= A.DGRID[-1] * 0.999)), frac_at_lower_edge=float(np.mean(Db <= A.DGRID[0] * 1.001)))

# ---------------- P1 coach
ev = A.t2_event(); df = ev["df"]
m = df[df["s3_measured"].notna() & (df.seat != "12D")]          # measured points, reference seat excluded
kap_test = rng.uniform(5.79, 7.25, 200)
D_hat, sse = A.fit_D(ev["lat"], ev["src_test"], m["site"].to_numpy(), m["s3_measured"].to_numpy(), 6.5)
Db = A.boot_D(ev["lat"], ev["src_test"], m["site"].to_numpy(), m["s3_measured"].to_numpy(), kap_test, rng)
both = df[df["s3_measured"].notna() & df["s3_cfd"].notna() & (df.seat != "12D")]
lr = np.log(both["s3_measured"] / both["s3_cfd"])
c = df[df["s4_cfd"].notna() & (df.seat != "12C")]
D_cfd, sse_cfd = A.fit_D(ev["lat"], ev["src"], c["site"].to_numpy(), c["s4_cfd"].to_numpy(), 4.8)
sus = ev["sus"]
out["coach"] = dict(n_measured=int(len(m)), seats_measured=m.seat.tolist(), measured=m.s3_measured.tolist(),
                    D_hat=D_hat, sse=sse.tolist(), boot=summ(Db), boot_draws=Db.tolist(),
                    sigma_log_meas_vs_cfd=float(lr.std(ddof=1)), mean_log_meas_vs_cfd=float(lr.mean()),
                    D_hat_cfd_S4A=D_cfd, sse_cfd=sse_cfd.tolist(),
                    s4_mean_rear=float(sus.s4_cfd[ev["near"]].mean()), s4_mean_front=float(sus.s4_cfd[~ev["near"]].mean()),
                    measured_rear_mean=float(m[m.row >= 8].s3_measured.mean()), measured_front_mean=float(m[m.row < 8].s3_measured.mean()))
print("coach", {k: v for k, v in out["coach"].items() if k not in ("sse", "sse_cfd", "boot_draws")})

# ---------------- P1 minibus (K-A only; sigma for the noise)
b2 = pd.read_csv(A.DER + "/ou2022_B2_seats.csv")
bb = b2[b2.s3_measured.notna() & b2.s3_cfd.notna() & (b2.seat != "4B")]
lr2 = np.log(bb.s3_measured / bb.s3_cfd)
out["minibus"] = dict(sigma_log_meas_vs_cfd=float(lr2.std(ddof=1)), n=int(len(bb)),
                      sigma_log_passenger_seats=float(lr2[bb.seat != "D"].std(ddof=1)))
print("minibus", out["minibus"])

# ---------------- P2 aircraft cabin
out["cabin"] = {}
for zm in (False, True):
    ev5, kern = A.t5_direct_kernels(zero_missing=zm)
    maps = A.kinahan_maps(zm)
    fits = {}
    lat = ev5["lat"]
    colidx = {"A": 0, "D": 2, "G": 3, "L": 5}
    for rel in ("5L", "5A"):
        mp = maps[rel].stack().reset_index(); mp.columns = ["row", "col", "count"]
        mp = mp[(mp["count"] > 0) & ~((mp.row == int(rel[:-1])) & (mp.col == rel[-1]))]
        rec = lat.site(mp.col.map(colidx).to_numpy(), mp.row.to_numpy().astype(int) - 1 + E.T5_ROW0)
        src = int(lat.site(colidx[rel[-1]], int(rel[:-1]) - 1 + E.T5_ROW0))
        Dh, sse = A.fit_D(lat, src, rec, mp["count"].to_numpy(), 35.0)
        fits[rel] = dict(D_hat=Dh, n=int(len(mp)), sse_min=float(sse.min()), sse_wellmixed=float(sse[-1]), sse=sse.tolist())
    Dg = float(np.exp(np.mean([np.log(f["D_hat"]) for f in fits.values()])))
    near, far = ev5["near"], ~ev5["near"]
    out["cabin"]["zeros_missing" if zm else "zeros_kept"] = dict(
        fits=fits, D_geomean=Dg,
        direct={k: dict(values=dict(zip(ev5["seats"], np.round(v[0], 1).tolist())), n_imputed=v[1],
                        near_far_mean_ratio=float(v[0][near].mean() / v[0][far].mean())) for k, v in kern.items()})
    print("cabin", "zeros missing" if zm else "zeros kept", {r: (f["D_hat"], f["n"]) for r, f in fits.items()}, "geomean", Dg,
          {k: (v[1], round(float(v[0][near].mean() / v[0][far].mean()), 2)) for k, v in kern.items()})

# ---------------- P3 rail carriage
pts = pd.read_csv(A.DER + "/woodward2022_fig9_points.csv")
latw = L.Lattice(nx=6, ny=20, ax=2.69 / 6, ay=19.85 / 20, H=2.02)
out["train"] = {}
for panel, xsrc in (("B_middle", -0.135), ("A_end", -0.392)):
    d = pts[(pts.panel == panel) & (~pts.near_source)]
    ys = np.clip(np.floor((d.x_over_L.to_numpy() + 0.5) * 20).astype(int), 0, 19)
    xs = np.where(d.colour == "red", 1, 4)                     # y < 0 side / y >= 0 side
    rec = latw.site(xs, ys)
    src = int(latw.site(1, int(np.floor((xsrc + 0.5) * 20))))
    Dh, sse = A.fit_D(latw, src, rec, d.c_norm.to_numpy(), 13.0)
    Db = A.boot_D(latw, src, rec, d.c_norm.to_numpy(), rng.uniform(11, 15, 200), rng)
    ell = float(np.sqrt(Dh / 13.0))
    out["train"][panel] = dict(D_hat=Dh, n=int(len(d)), sse=sse.tolist(), boot=summ(Db), boot_draws=Db.tolist(), mixing_length_m=ell)
    print("train", panel, Dh, summ(Db), "ell", round(ell, 2), "n", len(d))
cur = pd.read_csv(A.DER + "/woodward2022_fig9_curve.csv")
for panel, xsrc in (("B_middle", -0.135), ("A_end", -0.392)):
    d = cur[cur.panel == panel]
    d = d[np.abs(d.x_over_L - xsrc) > 0.02].iloc[::4]
    ys = np.clip(np.floor((d.x_over_L.to_numpy() + 0.5) * 20).astype(int), 0, 19)
    rec = latw.site(np.full(len(ys), 2), ys)
    src = int(latw.site(1, int(np.floor((xsrc + 0.5) * 20))))
    Dh, sse = A.fit_D(latw, src, rec, np.clip(d.c_norm.to_numpy(), 0.05, None), 13.0)
    out["train"][panel]["D_hat_from_1D_curve"] = Dh
    print("train curve check", panel, Dh)

# ---------------- P5 room relation (Cheng et al. 2011)
out["room_relation"] = dict(slope=0.52, intercept=0.31, sigma_log=0.42,
                            note="K/L^2 = 0.52 ACH + 0.31 per hour, L = V^(1/3); R^2 = 0.92, n = 11")
for nm, V, ach in (("C1 classroom", 11 * 13 * 3.0, (2, 12)), ("O1 call centre (pitch 1.3 m)", None, (1, 6))):
    if V is None:
        cfg = E.build_O1(np.random.default_rng(0)); V = cfg["lat"].M * 1.3 ** 2 * 2.7
    lo, hi = [(0.52 * a + 0.31) * V ** (2 / 3) for a in ach]
    out["room_relation"][nm] = dict(V=float(V), D_at_low_ACH=float(lo), D_at_high_ACH=float(hi))
    print(nm, "V", round(V), "D range m2/h", round(lo, 1), round(hi, 1))
A.save_json("01_kernels.json", out)
