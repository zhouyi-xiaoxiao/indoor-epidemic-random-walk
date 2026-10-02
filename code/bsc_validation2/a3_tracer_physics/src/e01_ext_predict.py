"""EXT-F step 1 (PREREG_EXT_flights.md): predictive tables of P, P-econ, M0, CRR(2), CRR(3) for the four
SARS-CoV-2 flights of avenue 1, and the power of the rule.  Uses the totals K and bin sizes only; the observed
bin vectors (meta['obs'], meta['y']) are not read.  Writes results/E1_ext_predictions.json, E1_ext_tables.npz."""
import os, sys
import numpy as np, pandas as pd
import a3lib as A
from a3lib import L
sys.path.insert(0, os.path.normpath(os.path.join(A.ROOT, "..", "a1_more_outbreaks", "src")))
from a1 import events as E1, engine as G

FL = ["F5", "F6", "F7", "F8"]
ker = A.load_json("01_kernels.json")
out = {"stamp": A.stamp()}

# ---------------- physics D, economy sections (sensitivity P-econ)
df = pd.read_csv(os.path.join(A.DER, "kinahan_tidy.csv"))
g = df[(df.kind == "B") & (df["mask"] == "No") & df.row.notna()].copy()
g["airframe"] = g.airframe.astype(str)
LAY = {"343": (12, dict(zip("ABCDEFGJKL", [0, 1, 2, 4, 5, 6, 7, 9, 10, 11]))),
       "242": (10, dict(zip("ABDEFGKL", [0, 1, 3, 4, 5, 6, 8, 9]))),
       "232": (9, dict(zip("ABDEFKL", [0, 1, 3, 4, 5, 7, 8])))}
econ = []
for (af, sec), d in g.groupby(["airframe", "section"]):
    if af == "777" and sec == "FWD":
        continue
    cols = set(d.col.unique())
    lay = "343" if cols & {"C", "J"} else ("242" if "G" in cols else "232")
    nx, cmap = LAY[lay]
    r0, r1 = int(d.row.min()), int(d.row.max())
    lat = L.Lattice(nx=nx, ny=r1 - r0 + 21, ax=0.5, ay=0.8, H=2.0)
    kap = 35.0 if af == "777" else 32.0
    for rel, dd in d.groupby("release"):
        m = dd.groupby(["row", "col"])["count"].mean().reset_index()
        rr, rc = int(rel[:-1]), rel[-1]
        m = m[(m["count"] > 0) & ~((m.row == rr) & (m.col == rc)) & m.col.isin(list(cmap))]
        rec = lat.site(m.col.map(cmap).to_numpy(), m.row.to_numpy().astype(int) - r0 + 10)
        src = int(lat.site(cmap[rc], rr - r0 + 10))
        Dh, sse = A.fit_D(lat, src, rec, m["count"].to_numpy(), kap)
        econ.append(dict(airframe=af, section=sec, layout=lay, release=rel, n=int(len(m)), D_hat=Dh,
                         sse_min=float(sse.min()), sse_wellmixed=float(sse[-1])))
        print(af, sec, lay, rel, len(m), Dh, flush=True)
out["econ_fits"] = econ
D_econ = np.array([e["D_hat"] for e in econ])
out["D_econ_geomean"] = float(np.exp(np.log(D_econ).mean())); out["D_econ_range"] = [float(D_econ.min()), float(D_econ.max())]
D_fwd = np.array([ker["cabin"]["zeros_kept"]["fits"][r]["D_hat"] for r in ("5L", "5A")])
out["D_fwd"] = D_fwd.tolist()
kap_nodes = 30.0 + 6.0 * (np.arange(7) + 0.5) / 7 + E1.DEP

# ---------------- predictive tables
tabs = {}
for ev in FL:
    meta, draws = E1.build(ev)
    sizes, K = meta["sizes"], meta["K"]                    # K only; meta['obs'] is not used here
    comps = G.compositions(sizes, K)
    bins0 = draws[0]["bins"]; near = np.isin(bins0, meta["near_bins"])
    pm = {"M0": G.cond_pmf(np.ones(sum(sizes)), bins0, sizes, K, comps)}
    for rho in (2.0, 3.0):
        pm["CRR%d" % rho] = G.cond_pmf(np.where(near, rho, 1.0), bins0, sizes, K, comps)
    seatmap = len(draws) <= 7
    cfgs = [draws[0]] * 7 if seatmap else draws
    for name, Ds in (("P", D_fwd), ("P_econ", D_econ), ("P_D708", D_fwd[:1]), ("P_D112", D_fwd[1:])):
        acc = np.zeros(len(comps)); n = 0
        for i, d in enumerate(cfgs):
            dd = dict(d); dd["kappa"] = float(kap_nodes[i % 7])
            for D in Ds:
                X, _ = G.x_m2(dd, float(D), "M2")
                acc += G.cond_pmf(X, dd["bins"], sizes, K, comps); n += 1
        pm[name] = acc / n
    nb = np.isin(np.arange(len(sizes)), meta["near_bins"])
    n1, n2 = int(np.array(sizes)[nb].sum()), int(np.array(sizes)[~nb].sum())
    rr = {}
    for k_, p in pm.items():
        e1 = float((p * comps[:, nb].sum(axis=1)).sum()); rr[k_] = (e1 / n1) / ((K - e1) / n2)
    tabs[ev] = dict(comps=comps, **pm)
    out[ev] = dict(label=meta["label"], sizes=sizes, K=int(K), near_bins=list(map(int, meta["near_bins"])), n_near=n1, n_far=n2,
                   n_comp=int(len(comps)), rr_expected=rr, T_h=meta["T"], n_sources=int(len(draws[0]["src"])))
    print(ev, meta["label"], sizes, K, {k_: round(v, 2) for k_, v in rr.items()}, flush=True)
np.savez(os.path.join(A.RES, "E1_ext_tables.npz"), **{f"{ev}|{k_}": v for ev in FL for k_, v in tabs[ev].items()})

# ---------------- power of the rule (no observed vector)
def pvals(p):                       # two-sided p of every composition under table p
    o = np.argsort(p); c = np.cumsum(p[o]); pv = np.empty_like(p)
    # sum of probabilities <= p_i (ties included)
    pv[o] = c; 
    srt = p[o]; idx = np.searchsorted(srt, p * (1 + 1e-12), side="right") - 1
    return c[idx]
lg = lambda a: np.log(np.clip(a, 1e-300, None))
rng = np.random.default_rng(A.SEED + 40); NMC = 200000
draw = {t: {ev: rng.choice(len(tabs[ev]["comps"]), size=NMC, p=tabs[ev][t] / tabs[ev][t].sum()) for ev in FL} for t in ("P", "M0", "CRR2", "CRR3")}
pw = {}
for target in ("P", "P_econ"):
    pvt = {ev: pvals(tabs[ev][target]) for ev in FL}
    delta = lambda t, c: sum(lg(tabs[ev][target][draw[t][ev]]) - lg(tabs[ev][c][draw[t][ev]]) for ev in FL)
    c95 = {c: float(np.quantile(delta(c, c), 0.95)) for c in ("M0", "CRR2", "CRR3")}
    res = {}
    for t in ("P", "M0", "CRR2", "CRR3"):
        nf = sum((pvt[ev][draw[t][ev]] < 0.05).astype(int) for ev in FL)
        d0, d2, d3 = delta(t, "M0"), delta(t, "CRR2"), delta(t, "CRR3")
        B0 = (d0 > 0) & (d0 > c95["M0"]); BC = (d2 > 0) & (d2 > c95["CRR2"]) & (d3 > 0) & (d3 > c95["CRR3"])
        o = np.where(nf >= 2, 5, np.where(nf == 1, 4, np.where(~B0, 3, np.where(BC, 1, 2))))
        res[t] = dict(outcomes={"E%d" % k: float((o == k).mean()) for k in (1, 2, 3, 4, 5)}, beats_M0=float(B0.mean()), beats_CRR=float(BC.mean()),
                      per_event_reject={ev: float((pvt[ev][draw[t][ev]] < 0.05).mean()) for ev in FL})
    pw[target] = dict(c95=c95, truths=res)
    print(target, "c95", {k_: round(v, 2) for k_, v in c95.items()})
    for t in res: print("  truth", t, {k_: round(v, 3) for k_, v in res[t]["outcomes"].items()}, "beats M0", round(res[t]["beats_M0"], 3), "beats CRR", round(res[t]["beats_CRR"], 3))
out["power"] = pw
A.save_json("E1_ext_predictions.json", out)
