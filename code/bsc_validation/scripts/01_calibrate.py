"""Calibration on the calibration records only (pre-registration 5.1).

Vehicle-cabin class: T4 high-speed trains (binomial on reconstructed cells).
Seated-room class:   R1 Guangzhou restaurant (interval-censored tables B and C).

For every model and every pre-registered likelihood variant the script stores
the profile log-likelihood on the log10 D grid (eps profiled out; kappa of the
calibration event integrated over its prior for T4 and drawn per layout for
R1), the posterior weights, the fitted cell/table probabilities at the mode and
AIC/BIC.  No hold-out record is read.

Usage: 01_calibrate.py [variant ...]   (default: all variants)
"""
import sys

import numpy as np
from scipy import optimize, stats

import _common as C
from valmod import events as E
from valmod import lattice as L

LOGD = np.linspace(-2, 4, 61)
DGRID = 10.0 ** LOGD


# =============================================================================
# T4
# =============================================================================
def t4_setup(row_pitch=1.0, theta_y=1.0):
    lat = E.t4_lattice(row_pitch, theta_y)
    sites, rows, cols = E.t4_seats(lat)
    df = E.t4_data()
    dr = np.abs(rows[:, None] - rows[None, :])
    dc = np.abs(cols[:, None] - cols[None, :])
    cell_pairs = []
    for r, c in zip(df["dr"], df["dc"]):
        cell_pairs.append(np.flatnonzero(((dr == r) & (dc == c)).ravel()))
    return lat, sites, df, cell_pairs


def t4_cellprob(lat, sites, cell_pairs, model, D, kappa, Tn, Tw):
    """Returns f(eps) -> model attack rate per cell (duration-averaged)."""
    allidx = np.concatenate(cell_pairs)
    lab = np.concatenate([np.full(len(ix), c) for c, ix in enumerate(cell_pairs)])
    cnt = np.bincount(lab, minlength=len(cell_pairs)).astype(float)
    Xall = np.stack([np.clip(L.exposure_matrix(lat, model, D, kappa, [T], sites).ravel()[allidx],
                             0.0, None) for T in Tn])                 # (nT, npairs)
    Tw = np.asarray(Tw)

    def prob(eps):
        q = Tw @ (-np.expm1(-eps * Xall))                             # duration-averaged
        return np.bincount(lab, weights=q, minlength=len(cell_pairs)) / cnt
    return prob


def t4_profile(prob, k, n, mask):
    def nll(le):
        p = np.clip(prob(np.exp(le)), 1e-300, 1 - 1e-12)
        return -np.sum(stats.binom.logpmf(k[mask], n[mask], p[mask]))
    # coarse bracket then Brent
    les = np.linspace(-12, 25, 38)
    v = np.array([nll(le) for le in les])
    i = int(np.argmin(v))
    lo, hi = les[max(i - 1, 0)], les[min(i + 1, len(les) - 1)]
    res = optimize.minimize_scalar(nll, bounds=(lo, hi), method="bounded",
                                   options={"xatol": 1e-6})
    return -res.fun, float(np.exp(res.x))


def t4_logit_profile(prob, df, mask):
    """S-logit: normal likelihood on the logit scale, cells with lo > 0."""
    lo, hi, pct = df["lo"].to_numpy() / 100, df["hi"].to_numpy() / 100, df["pct"].to_numpy() / 100
    ok = mask & (lo > 0)
    lg = lambda p: np.log(p / (1 - p))
    se = (lg(hi[ok]) - lg(lo[ok])) / 3.92
    y = lg(pct[ok])

    def nll(le):
        p = np.clip(prob(np.exp(le)), 1e-300, 1 - 1e-12)
        return 0.5 * np.sum(((lg(p[ok]) - y) / se) ** 2) + np.sum(np.log(se * np.sqrt(2 * np.pi)))
    les = np.linspace(-12, 25, 38)
    v = np.array([nll(le) for le in les])
    i = int(np.argmin(v))
    res = optimize.minimize_scalar(nll, bounds=(les[max(i - 1, 0)], les[min(i + 1, len(les) - 1)]),
                                   method="bounded", options={"xatol": 1e-6})
    return -res.fun, float(np.exp(res.x)), int(ok.sum())


def calibrate_T4(variant):
    row_pitch = 0.4 if variant == "S-pitch" else 1.0
    lat, sites, df, cell_pairs = t4_setup(row_pitch)
    k, n = df["k"].to_numpy(), df["n"].to_numpy()
    dr, dc = df["dr"].to_numpy(), df["dc"].to_numpy()
    if variant == "S-adj":
        mask = np.ones(len(df), bool)
    elif variant == "S-row":
        mask = dr > 0
    else:
        mask = ~((dr == 0) & (dc == 1))
    Tn, Tw = E.t4_duration_nodes(16)
    kq, kw = E.kappa_quadrature("T4", 7)
    use_logit = variant == "S-logit"
    out = dict(variant=variant, logD=LOGD.tolist(), kappa_nodes=kq.tolist(),
               cells=df.to_dict(orient="list"), mask=mask.tolist(), row_pitch=row_pitch,
               n_cells=int(mask.sum()))

    def fit(prob):
        if use_logit:
            ll, eps, nc = t4_logit_profile(prob, df, mask)
            out["n_cells"] = nc
            return ll, eps
        return t4_profile(prob, k, n, mask)

    # M0
    prob0 = t4_cellprob(lat, sites, cell_pairs, "M0", 1.0, kq[3], Tn, Tw)
    ll0, eps0 = fit(prob0)
    out["M0"] = dict(logL=ll0, eps=eps0, p=prob0(eps0).tolist(), n_par=1)
    # M1
    ll1, eps1 = np.zeros(len(DGRID)), np.zeros(len(DGRID))
    for i, D in enumerate(DGRID):
        prob = t4_cellprob(lat, sites, cell_pairs, "M1", D, 0.0, Tn, Tw)
        ll1[i], eps1[i] = fit(prob)
    ib = int(np.argmax(ll1))
    pb = t4_cellprob(lat, sites, cell_pairs, "M1", DGRID[ib], 0.0, Tn, Tw)(eps1[ib])
    out["M1"] = dict(profile_logL=ll1.tolist(), eps=eps1.tolist(), logL=float(ll1[ib]),
                     D_hat=float(DGRID[ib]), eps_hat=float(eps1[ib]), p=pb.tolist(), n_par=2)
    C.save_json(f"01_cal_T4_{variant}.json", out)
    # M2
    ll2 = np.zeros((len(DGRID), len(kq)))
    eps2 = np.zeros_like(ll2)
    for i, D in enumerate(DGRID):
        for j, kap in enumerate(kq):
            prob = t4_cellprob(lat, sites, cell_pairs, "M2", D, kap, Tn, Tw)
            ll2[i, j], eps2[i, j] = fit(prob)
    ib, jb = np.unravel_index(np.argmax(ll2), ll2.shape)
    pb = t4_cellprob(lat, sites, cell_pairs, "M2", DGRID[ib], kq[jb], Tn, Tw)(eps2[ib, jb])
    # marginal over kappa prior (log-mean-exp)
    mx = ll2.max()
    marg = mx + np.log(np.sum(np.exp(ll2 - mx) * kw[None, :], axis=1))
    out["M2"] = dict(profile_logL=ll2.tolist(), eps=eps2.tolist(), logL=float(ll2[ib, jb]),
                     D_hat=float(DGRID[ib]), kappa_hat=float(kq[jb]),
                     eps_hat=float(eps2[ib, jb]), p=pb.tolist(), n_par=2,
                     marg_logL=marg.tolist(), cell_volume=lat.cell_volume)
    # saturated log-likelihood for the deviance
    if not use_logit:
        ps = k[mask] / n[mask]
        out["logL_saturated"] = float(np.sum(stats.binom.logpmf(k[mask], n[mask], ps)))
    C.save_json(f"01_cal_T4_{variant}.json", out)
    print(f"T4 {variant}: M0 {ll0:.2f} | M1 {out['M1']['logL']:.2f} at D={out['M1']['D_hat']:.3g} | "
          f"M2 {out['M2']['logL']:.2f} at D={out['M2']['D_hat']:.3g}, kappa={out['M2']['kappa_hat']:.2f}",
          flush=True)
    return out


def calibrate_T4_aniso():
    """S-aniso: fitted factor theta <= 1 on the between-row hop rate."""
    thetas = 10.0 ** np.linspace(-3, 0, 13)
    Tn, Tw = E.t4_duration_nodes(16)
    kq, kw = E.kappa_quadrature("T4", 7)
    Dg = DGRID[::2]
    res = dict(variant="S-aniso", logD=np.log10(Dg).tolist(), theta=thetas.tolist(),
               kappa_nodes=kq.tolist())
    ll1 = np.zeros((len(thetas), len(Dg)))
    e1 = np.zeros_like(ll1)
    ll2 = np.zeros((len(thetas), len(Dg), len(kq)))
    e2 = np.zeros_like(ll2)
    for a, th in enumerate(thetas):
        lat, sites, df, cell_pairs = t4_setup(1.0, th)
        k, n = df["k"].to_numpy(), df["n"].to_numpy()
        mask = ~((df["dr"].to_numpy() == 0) & (df["dc"].to_numpy() == 1))
        for i, D in enumerate(Dg):
            ll1[a, i], e1[a, i] = t4_profile(
                t4_cellprob(lat, sites, cell_pairs, "M1", D, 0.0, Tn, Tw), k, n, mask)
            for j in (1, 3, 5):
                ll2[a, i, j], e2[a, i, j] = t4_profile(
                    t4_cellprob(lat, sites, cell_pairs, "M2", D, kq[j], Tn, Tw), k, n, mask)
        res["M1"] = dict(profile_logL=ll1.tolist(), eps=e1.tolist())
        res["M2"] = dict(profile_logL=ll2[:, :, [1, 3, 5]].tolist(), eps=e2[:, :, [1, 3, 5]].tolist(),
                         kappa_used=kq[[1, 3, 5]].tolist())
        C.save_json("01_cal_T4_S-aniso.json", res)
        print("aniso theta", th, ll1[a].max(), ll2[a][:, [1, 3, 5]].max(), flush=True)
    for m, ll in (("M1", ll1), ("M2", ll2[:, :, [1, 3, 5]])):
        idx = np.unravel_index(np.argmax(ll), ll.shape)
        res[m]["logL"] = float(ll[idx])
        res[m]["theta_hat"] = float(thetas[idx[0]])
        res[m]["D_hat"] = float(Dg[idx[1]])
        res[m]["n_par"] = 3
        if m == "M2":
            res[m]["kappa_hat"] = float(kq[[1, 3, 5]][idx[2]])
            res[m]["eps_hat"] = float(e2[:, :, [1, 3, 5]][idx])
        else:
            res[m]["eps_hat"] = float(e1[idx])
    C.save_json("01_cal_T4_S-aniso.json", res)
    return res


# =============================================================================
# R1
# =============================================================================
def r1_layouts(n_lay=300):
    rng = np.random.default_rng(C.SEED + 101)
    return [E.build_R1(rng) for _ in range(n_lay)]


def r1_exposures(lays, model, D):
    """X[layout, table]"""
    X = np.zeros((len(lays), len(lays[0]["sites"])))
    for a, ly in enumerate(lays):
        lat = ly["lat"]
        for T in np.unique(ly["T"]):
            sel = ly["T"] == T
            X[a, sel] = L.exposure(lat, model, D, ly["kappa"], [T], ly["src"], ly["sites"][sel])
    return np.clip(X, 0.0, None)


def r1_loglik(X, eps, lays, variant):
    """log of the layout-averaged likelihood."""
    ly = lays[0]
    n = ly["n"]
    names = ly["name"]
    iB, iC = names.index("TB"), names.index("TC")
    others = np.array([i for i in range(len(names)) if i not in (iB, iC)])
    p = -np.expm1(-eps * X)                                  # (layout, table)
    p = np.clip(p, 0.0, 1 - 1e-15)
    ll = np.sum(n[others][None, :] * np.log1p(-p[:, others]), axis=1)
    if variant == "S-R1-upper":
        kb, kc = [3], [2]
    elif variant == "S-R1-lower":
        kb, kc = [1], [1]
    else:
        kb, kc = [1, 2, 3], [1, 2]
    PB = sum(stats.binom.pmf(kk, n[iB], p[:, iB]) for kk in kb)
    PC = sum(stats.binom.pmf(kk, n[iC], p[:, iC]) for kk in kc)
    ll = ll + np.log(np.clip(PB, 1e-300, None)) + np.log(np.clip(PC, 1e-300, None))
    m = ll.max()
    return m + np.log(np.mean(np.exp(ll - m)))


def r1_profile(X, lays, variant):
    f = lambda le: -r1_loglik(X, np.exp(le), lays, variant)
    les = np.linspace(-10, 60, 71)
    v = np.array([f(le) for le in les])
    i = int(np.argmin(v))
    res = optimize.minimize_scalar(f, bounds=(les[max(i - 1, 0)], les[min(i + 1, len(les) - 1)]),
                                   method="bounded", options={"xatol": 1e-6})
    return -res.fun, float(np.exp(res.x))


def calibrate_R1(variant):
    lays = r1_layouts()
    ly = lays[0]
    out = dict(variant=variant, logD=LOGD.tolist(), tables=ly["name"], n=ly["n"].tolist(),
               T=ly["T"].tolist(), cls=ly["cls"], n_layouts=len(lays),
               cell_volume=ly["lat"].cell_volume)
    X0 = r1_exposures(lays, "M0", 1.0)
    ll0, eps0 = r1_profile(X0, lays, variant)
    out["M0"] = dict(logL=ll0, eps=eps0, n_par=1,
                     p=(-np.expm1(-eps0 * X0)).mean(axis=0).tolist())
    for model in ("M1", "M2"):
        ll, ep = np.zeros(len(DGRID)), np.zeros(len(DGRID))
        for i, D in enumerate(DGRID):
            X = r1_exposures(lays, model, D)
            ll[i], ep[i] = r1_profile(X, lays, variant)
        ib = int(np.argmax(ll))
        Xb = r1_exposures(lays, model, DGRID[ib])
        out[model] = dict(profile_logL=ll.tolist(), eps=ep.tolist(), logL=float(ll[ib]),
                          D_hat=float(DGRID[ib]), eps_hat=float(ep[ib]), n_par=2,
                          p=(-np.expm1(-ep[ib] * Xb)).mean(axis=0).tolist(),
                          x_rel_src=(Xb / L.exposure(ly["lat"], model, DGRID[ib], ly["kappa"],
                                                     [82 / 60], ly["src"], [ly["src"]])[0]
                                     ).mean(axis=0).tolist())
        print(f"R1 {variant} {model}: logL {ll[ib]:.3f} at D={DGRID[ib]:.3g} (M0 {ll0:.3f})", flush=True)
    C.save_json(f"01_cal_R1_{variant}.json", out)
    return out


if __name__ == "__main__":
    todo = sys.argv[1:] or ["T4:primary", "R1:primary", "T4:S-adj", "T4:S-row", "T4:S-logit",
                            "T4:S-pitch", "R1:S-R1-upper", "R1:S-R1-lower", "T4:S-aniso"]
    for job in todo:
        ev, var = job.split(":")
        print(C.stamp(), "start", job, flush=True)
        if ev == "T4" and var == "S-aniso":
            calibrate_T4_aniso()
        elif ev == "T4":
            calibrate_T4(var)
        else:
            calibrate_R1(var)
        print(C.stamp(), "done", job, flush=True)
