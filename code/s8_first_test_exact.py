"""s8_first_test_exact.py -- the first outbreak test re-run with a round-off-safe propagator.

Why.  The first test computes exposures with the cosine-mode sum of code/bsc_validation/src/valmod/
lattice.py.  Its absolute error reaches about 1e-10 of the largest exposure (kernel.rows.spectral_min_rel),
so receivers whose exposure is smaller than that (far desks of the classroom at mixing coefficients of
0.1 m^2/h and below) get round-off noise of either sign instead of values of order 1e-16 to 1e-33.  The second outbreak test found the same defect in its
copy of the code and replaced such values by a positive image sum
(code/bsc_validation2/a1_more_outbreaks/src/a1/engine_fix.py).  This script applies the same rule to the
first test: the spectral value is kept wherever it is resolved (exposure above 1e-8 of the largest exposure of
the lattice); elsewhere it is replaced by the image-sum value

    g(x, s | x0) = sum_m e^{-2ws} [ I_|x-x0-2nm|(2ws) + I_|x+x0+1-2nm|(2ws) ],   w = D / a^2   (per axis),
    M2: X = int_0^T (T-s) e^{-kappa s} gx gy ds,      M1: X = (1/2) int_0^{2T} gx gy ds,

a quadrature of a positive integrand.  The calibrations are not re-run: the error they could carry is below
1e-7 in the log-likelihood (code/checks/outbreaks2/t6_first_set_floor.json).

What is recomputed (same seeds, draws and builders as the stored analysis; nothing else changed):
  kernel    check of the image sum against a third method (uniformisation, from the check of the second
            test) on classroom layouts, and the number of non-positive spectral exposures;
  primary   predictive pmfs of the four held-out events for M1 (closed form) and M2 (seeds of 02_power.py),
            the event rows of 03_holdout_shape.py and the pooled test A1, with and without single events;
  power     the power table of 02_power.py with the corrected pmfs;
  breakdown share of the classroom's predictive probability by range of D (as validation_checks/v11);
  loo       the per-event likelihood curves of the room class (07_loo.py) and the common-coefficient test;
  mc        Monte Carlo stability over five other seeds (02c_mc_stability.py), M2;
  sens      the twelve rows of the sensitivity table (09_sensitivity.py);
  vent      ventilation priors of the bus, the flight and the classroom scaled by 1/3 and 3, and the emission
            rate implied by the classroom (as validation_checks/v09_sens.py, with the analysis's own
            builders);
  figures   Figs. "hold-out events of the first test" and "log-likelihood of the mixing coefficient for each
            event" redrawn from the corrected values (English and Chinese).
Output: ../data/s8_first_test_exact.json (checkpointed after every part; a finished part is not recomputed
        unless it is named on the command line or --force is given),
        ../figures/figV2_holdout_shape_{en,zh}.pdf, ../figures/figV8_mixing_heterogeneity_{en,zh}.pdf.
Nothing is written to the research directories.

    nice -n 19 python s8_first_test_exact.py [--force] [part ...]     (all parts: about 2 minutes)
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
from scipy.special import ive

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]                      # repository root
V1 = ROOT / "code/bsc_validation"
OUT = HERE.parent / "data" / "s8_first_test_exact.json"
FIG = HERE.parent / "figures"
sys.path.insert(0, str(V1 / "scripts"))
sys.path.insert(0, str(V1 / "src"))
import _common as C                         # noqa: E402  (reads only)
from valmod import decision as DC           # noqa: E402
from valmod import events as E              # noqa: E402
from valmod import lattice as L             # noqa: E402
from valmod import predict as P             # noqa: E402
from valmod import stats as S               # noqa: E402

EVENTS = ["T1", "T2", "T5", "C1"]
SPECTRAL = L.exposure
THRESH = 1e-8                               # as in the second test
NQ = 240
GL = np.polynomial.legendre.leggauss(NQ)
MIMG = np.arange(-3, 4)


# ------------------------------------------------------------------------------------- image-sum propagator
def g1d(n, x0, w, s):
    """(n, len(s)) 1-D reflecting-lattice propagator from site x0, as a sum of positive image terms."""
    x = np.arange(n)
    z = 2.0 * w * s
    o1 = np.abs(x[:, None] - x0 - 2 * n * MIMG[None, :])
    o2 = np.abs(x[:, None] + x0 + 1 - 2 * n * MIMG[None, :])
    return ive(o1[:, :, None], z[None, None, :]).sum(1) + ive(o2[:, :, None], z[None, None, :]).sum(1)


def x_images(lat, model, D, kappa, segments, src):
    X = np.zeros(lat.M)
    x0, y0 = int(src % lat.nx), int(src // lat.nx)
    for T in np.atleast_1d(segments):
        T = float(T)
        if model == "M2":
            lo, hi = np.log(T * 1e-9), np.log(T)
        else:                                               # M1: (1/2) int_0^{2T} P_s ds
            lo, hi = np.log(2 * T * 1e-9), np.log(2 * T)
        s = np.exp(0.5 * (hi - lo) * GL[0] + 0.5 * (hi + lo))
        wq = 0.5 * (hi - lo) * GL[1] * s
        wq = wq * ((T - s) * np.exp(-kappa * s) if model == "M2" else 0.5)
        gx = g1d(lat.nx, x0, D / lat.ax ** 2, s)
        gy = g1d(lat.ny, y0, D / lat.ay ** 2 * lat.theta_y, s)
        X += np.einsum("yq,xq,q->yx", gy, gx, wq).reshape(-1)
    return X


STATS = {"calls": 0, "calls_patched": 0, "values_replaced": 0}


def exposure_fixed(lat, model, D, kappa, segments, src, rec):
    rec = np.asarray(rec)
    if model == "M0":
        return SPECTRAL(lat, model, D, kappa, segments, src, rec)
    w = np.zeros(lat.M)
    for T in np.atleast_1d(segments):
        w += L.weights(lat, model, D, kappa, float(T))
    Xall = lat.Phi @ (w * lat.Phi[src])
    Xs = Xall[rec]
    bad = Xs < THRESH * Xall.max()
    STATS["calls"] += 1
    if bad.any():
        STATS["calls_patched"] += 1
        STATS["values_replaced"] += int(bad.sum())
        Xi = x_images(lat, model, D, kappa, segments, src)[rec]
        Xs = np.where(bad, Xi, Xs)
    return Xs


def use(kind):
    L.exposure = exposure_fixed if kind == "fixed" else SPECTRAL


# ---------------------------------------------------------------------------------------------- helpers
def load_out():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def save(res):
    OUT.parent.mkdir(exist_ok=True)
    tmp = OUT.with_suffix(".json.part")
    tmp.write_text(json.dumps(res, indent=1, default=float) + "\n")
    os.replace(tmp, OUT)


def r(x, n=6):
    return float(np.round(float(x), n))


OBS = C.load_json("heldout_outcomes.json")["shape_primary"]
KOBS = {ev: OBS[ev]["near_k"] for ev in EVENTS}
PRED = C.load_json("02_predictive.json")
EXACT_M1 = C.load_json("02b_m1_exact.json")["events"]


def m1_exact_pmf(ev):
    ex = EXACT_M1.get(ev, {})
    if not ex.get("replace_closed_form", False):
        return None
    pe = np.array(ex["pmf_sim"], float)
    pe = (pe * ex["n_cond"] + 0.5 / len(pe)) / (ex["n_cond"] + 0.5)
    return pe / pe.sum()


def pooled(pm_model, evs):
    pm0 = [PRED["events"][e]["M0"]["pmf"] for e in evs]
    d, p, per = DC.pooled_delta_test([pm_model[e] for e in evs], pm0, [KOBS[e] for e in evs])
    return dict(delta=r(d), p_under_M0=float(p), per_event={e: r(x) for e, x in zip(evs, per)})


# ------------------------------------------------------------------------------------------------- parts
def part_kernel(res):
    sys.path.insert(0, str((ROOT / "code/checks/outbreaks2/vlib.py").parent))   # the re-check's library
    import vlib as V                                                      # uniformisation (positive sums)
    rng = np.random.default_rng(3)
    rows = []
    for trial in range(3):
        cfg = E.build_C1(rng)
        lat, src, rec = cfg["lat"], cfg["src"], cfg["rec"]
        s0 = (int(src % lat.nx), int(src // lat.nx))
        for D in (0.01, 0.03, 0.1, 0.5, 1.0, 10.0):
            Xs = SPECTRAL(lat, "M2", D, cfg["kappa"], cfg["segments"], src, rec)
            Xf = exposure_fixed(lat, "M2", D, cfg["kappa"], cfg["segments"], src, rec)
            Xu = np.zeros(lat.M)
            for T in cfg["segments"]:
                Xu += V.kern_table(lat.nx, lat.ny, lat.ax, lat.ay, [[s0]], float(T), [cfg["kappa"]], [D])[0, 0, 0]
            Xu = Xu[rec]
            Xi = x_images(lat, "M2", D, cfg["kappa"], cfg["segments"], src)[rec]
            rows.append(dict(layout=trial, D=D, n_receivers=int(len(rec)),
                             spectral_nonpositive=int((Xs <= 0).sum()),
                             spectral_min_rel=float(Xs.min() / Xs.max()),
                             fixed_min_rel=float(Xf.min() / Xf.max()),
                             images_vs_uniformisation_max_rel_err=float(np.max(np.abs(Xi - Xu) / Xu)),
                             fixed_vs_uniformisation_max_rel_err=float(np.max(np.abs(Xf - Xu) / Xu))))
    res["kernel"] = dict(note="classroom C1 layouts, M2; relative errors are |a-b|/b per receiver", rows=rows)


def part_primary(res):
    out = {"events": {}, "pmf_change": {}}
    pm = {"M2": {}, "M1cf": {}, "M1": {}}
    post = {}
    for cls in ("vehicle", "room"):
        for m in ("M1", "M2"):
            post[(cls, m)] = P.posterior_D(cls, m, "primary", n=400, seed=C.SEED + 7)
    use("fixed")
    for ev in EVENTS:
        n1, n2 = E.STRATA_N[ev]
        K = E.TOTALS[ev]
        k = KOBS[ev]
        rec = {"n_near": n1, "n_far": n2, "K": K}
        for m in ("M1", "M2"):
            STATS.update(calls=0, calls_patched=0, values_replaced=0)
            rr = P.shape_predictive(ev, m, post[(E.CLASS[ev], m)]["D"], n_draw=300, seed=C.SEED + 13)
            p = np.asarray(rr["pmf"], float)
            stored = np.asarray(PRED["events"][ev][m]["pmf"], float)
            key = "M2" if m == "M2" else "M1cf"
            pm[key][ev] = p.tolist()
            out["pmf_change"][f"{ev}.{key}"] = dict(max_abs_change=float(np.max(np.abs(p - stored))),
                                                     draws_with_replaced_values=STATS["calls_patched"],
                                                     values_replaced=STATS["values_replaced"])
            rec[key] = DC.event_row(p, n1, n2, K, k)
            rec[key]["rr_expected_q"] = np.quantile(rr["rr_expected"], [0.05, 0.5, 0.95]).tolist()
            rec[key]["pmf"] = p.tolist()
            rec[key + "_stored"] = DC.event_row(stored, n1, n2, K, k)
        pe = m1_exact_pmf(ev)
        pm["M1"][ev] = (pe if pe is not None else np.asarray(pm["M1cf"][ev])).tolist()
        rec["M1"] = DC.event_row(pm["M1"][ev], n1, n2, K, k)
        rec["M1_uses_exact_simulator"] = pe is not None
        rec["M0"] = DC.event_row(PRED["events"][ev]["M0"]["pmf"], n1, n2, K, k)
        for m in ("M1", "M2", "M1cf"):
            rec[m]["delta_vs_M0"] = rec[m]["log_score"] - rec["M0"]["log_score"]
        out["events"][ev] = rec
        print(ev, {m: (round(rec[m]["rr_pred"], 2), round(rec[m]["p_two_sided"], 4), round(rec[m]["delta_vs_M0"], 2))
                   for m in ("M1", "M2")}, out["pmf_change"][f"{ev}.M2"], flush=True)
    use("spectral")
    out["pooled"] = {}
    for m in ("M1", "M2", "M1cf"):
        sets = {"all": EVENTS, "without_T5": ["T1", "T2", "C1"], "without_C1": ["T1", "T2", "T5"],
                "without_T5_C1": ["T1", "T2"]}
        out["pooled"][m] = {k: pooled(pm[m], v) for k, v in sets.items()}
        a = out["pooled"][m]["all"]
        out["pooled"][m]["A1_pass"] = bool(a["delta"] > 0 and a["p_under_M0"] < 0.05)
        out["pooled"][m]["failed"] = [ev for ev in EVENTS if out["events"][ev][m]["p_two_sided"] < 0.05]
    stored3 = C.load_json("03_holdout_shape.json")["pooled"]
    out["pooled_stored"] = {m: dict(delta=stored3[m]["delta"], p_under_M0=stored3[m]["p_under_M0"],
                                    failed=stored3[m]["failed"]) for m in ("M1", "M2", "M1cf")}
    res["primary"] = out


def part_power(res):
    spec = importlib.util.spec_from_file_location("power02", V1 / "scripts" / "02_power.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    pred = json.loads(json.dumps(PRED))
    for ev in EVENTS:
        pred["events"][ev]["M2"]["pmf"] = res["primary"]["events"][ev]["M2"]["pmf"]
        pred["events"][ev]["M1"]["pmf"] = res["primary"]["events"][ev]["M1cf"]["pmf"]
    pw = mod.power(pred)
    stored = C.load_json("02_power.json")
    keys = ("P_A1", "P_A2_nofail", "P_A2_le1", "P_W1_like")
    res["power"] = {c: {t: {k: pw[c][t][k] for k in keys} for t in ("M0", "M1", "M2")} for c in ("M1", "M2")}
    res["power_stored"] = {c: {t: {k: stored[c][t][k] for k in keys} for t in ("M0", "M1", "M2")} for c in ("M1", "M2")}


def part_breakdown(res):
    cal = P.load_cal("room")
    logD = np.array(cal["logD"])
    ll = np.array(cal["M2"]["profile_logL"])
    w = np.exp(ll - ll.max())
    w /= w.sum()
    k1 = KOBS["C1"]
    curves = {}
    for kind in ("spectral", "fixed"):
        use(kind)
        cur = []
        for lD in logD:
            rr = P.shape_predictive("C1", "M2", None, n_draw=80, seed=C.SEED + 91, fixed_D=10 ** lD)
            cur.append((np.asarray(rr["pmf"]) / np.sum(rr["pmf"])).tolist())
        curves[kind] = np.array(cur)
    use("spectral")
    out = {"logD": logD.tolist(), "posterior_w": w.tolist()}
    for kind, c in curves.items():
        tot = (w[:, None] * c).sum(0)
        d = {"pmf_obs": float(tot[k1]), "p_two_sided": S.two_sided_p(tot, k1), "ranges": []}
        for lo, hi in ((-2.01, -1.0), (-1.0, 0.75), (0.75, 4.01)):
            m = (logD > lo) & (logD <= hi)
            d["ranges"].append(dict(D_lo=10 ** lo, D_hi=10 ** hi, posterior_mass=float(w[m].sum()),
                                    share_of_pmf_obs=float((w[m] * c[m, k1]).sum() / tot[k1])))
        im = int(w.argmax())
        d["mode"] = dict(D=float(10 ** logD[im]), pmf_obs=float(c[im, k1]), p_two_sided=S.two_sided_p(c[im], k1),
                         k_mean=float((np.arange(c.shape[1]) * c[im]).sum()))
        d["p_by_D"] = [S.two_sided_p(c[i], k1) for i in range(len(logD))]
        out[kind] = d
    res["breakdown"] = out


def part_loo(res):
    lo = C.load_json("07_loo.json")
    LOGD = np.array(lo["logD"])
    k = KOBS["C1"]
    n1, n2 = E.STRATA_N["C1"]
    ll0 = float(np.log(S.hypergeom_pmf(n1, n2, E.TOTALS["C1"])[k]))
    use("fixed")
    v = []
    for ld in LOGD:
        rr = P.shape_predictive("C1", "M2", None, n_draw=120, seed=C.SEED + 71, fixed_D=10 ** ld)
        v.append(float(np.log(max(rr["pmf"][k], 1e-300))))
    use("spectral")
    room = lo["classes"]["room"]["M2"]
    curves = {"R1": np.array(room["curves"]["R1"]), "C1": np.array(v)}
    ll0s = {"R1": room["loo"]["R1"]["score_M0"], "C1": ll0}
    rows = {}
    for ev in ("R1", "C1"):
        other = curves["C1" if ev == "R1" else "R1"]
        wo = np.exp(other - other.max())
        wo /= wo.sum()
        sc = float(np.log(np.sum(wo * np.exp(curves[ev] - curves[ev].max()))) + curves[ev].max())
        rows[ev] = dict(score=sc, score_M0=ll0s[ev], delta=sc - ll0s[ev], D_own_mode=float(10 ** LOGD[curves[ev].argmax()]),
                        own_max=float(curves[ev].max()))
    joint = curves["R1"] + curves["C1"]
    sep = curves["R1"].max() + curves["C1"].max()
    from scipy import stats as st
    LR = 2 * (sep - joint.max())
    stored_LR = 2 * (room["sum_own_max"] - room["joint_max"])
    res["loo_room"] = dict(logD=LOGD.tolist(), C1_curve=v, C1_curve_stored=room["curves"]["C1"], loo=rows,
                           joint_D_mode=float(10 ** LOGD[joint.argmax()]), LR=float(LR), p=float(st.chi2.sf(LR, 1)),
                           LR_stored=float(stored_LR), p_stored=float(st.chi2.sf(stored_LR, 1)),
                           C1_own_max_D_stored=float(10 ** LOGD[int(np.argmax(room["curves"]["C1"]))]),
                           C1_at_stored_max_rel=float(v[int(np.argmax(room["curves"]["C1"]))] - max(v)))


def part_mc(res):
    use("fixed")
    out = []
    for s in range(1, 6):
        pm = {}
        p = {}
        for ev in EVENTS:
            post = P.posterior_D(E.CLASS[ev], "M2", n=400, seed=C.SEED + 1000 * s)
            rr = P.shape_predictive(ev, "M2", post["D"], n_draw=300, seed=C.SEED + 2000 * s)
            pm[ev] = (np.asarray(rr["pmf"]) / np.sum(rr["pmf"])).tolist()
            p[ev] = S.two_sided_p(np.asarray(pm[ev]), KOBS[ev])
        po = pooled(pm, EVENTS)
        out.append(dict(seed=s, p=p, delta=po["delta"], p_pooled=po["p_under_M0"]))
        print("mc", s, {e: round(x, 4) for e, x in p.items()}, po["delta"], flush=True)
    use("spectral")
    res["mc"] = out


def part_sens(res):
    spec = importlib.util.spec_from_file_location("sens09", V1 / "scripts" / "09_sensitivity.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    jobs = [("primary (closed form for M1)", {}), ("S-adj", dict(veh="S-adj")), ("S-row", dict(veh="S-row")),
            ("S-logit", dict(veh="S-logit")), ("S-pitch", dict(veh="S-pitch")), ("S-aniso", dict(theta=True)),
            ("S-R1-upper", dict(room="S-R1-upper")), ("S-R1-lower", dict(room="S-R1-lower")),
            ("S-5K", dict(kw={"T5": {"index_seat": (5, 3)}})), ("S-Luo", dict(kw={"T2": {"minutes": 150.0}})),
            ("S-corner", dict(kw={"C1": {"corner": True}})), ("S-nu", dict(nu=True))]
    done = res.get("sens", {})
    use("fixed")
    for tag, args in jobs:
        if tag in done:
            continue
        done[tag] = mod.run(tag, **args)
        res["sens"] = done
        save(res)
    use("spectral")
    res["sens"] = done


def part_vent(res):
    out = {}
    base = dict(E.BUILDERS)

    def scaled(km, evs):
        def wrap(b):
            def f(rng, **kw):
                cfg = b(rng, **kw)
                cfg["kappa"] = (cfg["kappa"] - 0.93) * km + 0.93
                return cfg
            return f
        return {e: (wrap(base[e]) if e in evs else base[e]) for e in base}

    post = {cls: P.posterior_D(cls, "M2", "primary", n=400, seed=C.SEED + 7) for cls in ("vehicle", "room")}
    use("fixed")
    for km in (None, 0.33, 3.0):
        E.BUILDERS.clear()
        E.BUILDERS.update(base if km is None else scaled(km, ("T1", "T5", "C1")))
        pm, row = {}, {}
        for ev in EVENTS:
            rr = P.shape_predictive(ev, "M2", post[E.CLASS[ev]]["D"], n_draw=300, seed=C.SEED + 13)
            pm[ev] = (np.asarray(rr["pmf"]) / np.sum(rr["pmf"])).tolist()
            row[ev] = dict(p=S.two_sided_p(np.asarray(pm[ev]), KOBS[ev]),
                           rr_pred=DC.event_row(pm[ev], *E.STRATA_N[ev], E.TOTALS[ev], KOBS[ev])["rr_pred"])
            if km is None and ev == "C1":
                cfg = base["C1"](np.random.default_rng(0))
                vol = cfg["lat"].cell_volume
                em = np.asarray(rr["eps"]) * vol / 0.5
                row[ev]["emission_q05_q50_q95"] = np.quantile(em, [0.05, 0.5, 0.95]).tolist()
                row[ev]["emission_share_above_1e4"] = float(np.mean(em > 1e4))
                row[ev]["cell_volume_m3"] = vol
        po = pooled(pm, EVENTS)
        out["primary" if km is None else f"ACH x{km}"] = dict(events=row, delta=po["delta"], p_pooled=po["p_under_M0"])
        print("vent", km, {e: round(v["p"], 4) for e, v in row.items()}, po["delta"], flush=True)
    E.BUILDERS.clear()
    E.BUILDERS.update(base)
    use("spectral")
    res["vent"] = out


def part_figures(res):
    spec = importlib.util.spec_from_file_location("figs20", V1 / "scripts" / "20_figures.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    pred = json.loads(json.dumps(PRED))
    sh = C.load_json("03_holdout_shape.json")
    for ev in EVENTS:
        pred["events"][ev]["M2"]["pmf"] = res["primary"]["events"][ev]["M2"]["pmf"]
        for m in ("M1", "M2"):
            sh["events"][ev][m]["p_two_sided"] = res["primary"]["events"][ev][m]["p_two_sided"]
    lo = C.load_json("07_loo.json")
    lo["classes"]["room"]["M2"]["curves"]["C1"] = res["loo_room"]["C1_curve"]
    repl = {"02_predictive.json": pred, "03_holdout_shape.json": sh, "07_loo.json": lo}
    orig = mod.C.load_json
    mod.C.load_json = lambda name: repl[name] if name in repl else orig(name)
    for fn, name in ((mod.fig_holdout_pmf, "figV2_holdout_shape"), (mod.fig_heterogeneity, "figV8_mixing_heterogeneity")):
        for lang in ("en", "zh"):
            mod.PL.setup(lang)
            fig = fn(lang)
            fig.savefig(FIG / f"{name}_{lang}.pdf", bbox_inches="tight")
            mod.plt.close(fig)
    mod.C.load_json = orig
    res["figures"] = [f"figures/{n}_{l}.pdf" for n in ("figV2_holdout_shape", "figV8_mixing_heterogeneity")
                      for l in ("en", "zh")]


PARTS = [("kernel", part_kernel), ("primary", part_primary), ("power", part_power), ("breakdown", part_breakdown),
         ("loo_room", part_loo), ("mc", part_mc), ("sens", part_sens), ("vent", part_vent), ("figures", part_figures)]

if __name__ == "__main__":
    force = "--force" in sys.argv
    want = [a for a in sys.argv[1:] if a != "--force"]
    res = load_out()
    if force:                                # recompute every part (or every named part) from scratch
        want = want or [name for name, _ in PARTS]
    for name, fn in PARTS:
        if want and name not in want:
            continue
        if name in res and name not in want and not (name == "sens" and len(res["sens"]) < 12):
            print("done:", name)
            continue
        if force:
            res.pop(name, None)              # parts that checkpoint internally (power, sens) start from scratch
        t0 = time.time()
        fn(res)
        res.setdefault("seconds", {})[name] = round(time.time() - t0, 1)
        res["stamp"] = time.strftime("%Y-%m-%d %H:%M:%S %Z")
        save(res)
        print(f"part {name}: {time.time() - t0:.1f} s", flush=True)
