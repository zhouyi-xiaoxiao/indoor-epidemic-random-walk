"""O1 Seoul call centre, north wing: predictive distribution of the within-wing
join-count statistic under M0, M1 and M2 (epidemic percolation, pre-registration
section 4), then comparison with the digitized seat map.

p_ij = 1 - exp(-eps X_ij) for an infectious desk i and a susceptible desk j;
X from the room-class posterior of D, the O1 ventilation prior and the pitch
prior.  eps is set by bisection so that outbreaks reaching at least 40 cases
have a mean final size of 79; outbreaks with 69-89 cases are retained.
"""
import numpy as np

import _common as C
from valmod import events as E
from valmod import lattice as L
from valmod import predict as P
from valmod import stats as S

N_PAR = 40
R_TUNE = 400
R_PROD = 1500


def percolate(Xe, rng, R):
    """Xe = eps * X (n x n, zero diagonal).  Returns boolean (R, n) final infected."""
    n = Xe.shape[0]
    inf = np.zeros((R, n), bool)
    new = np.zeros((R, n), bool)
    idx = rng.integers(n, size=R)
    new[np.arange(R), idx] = True
    inf |= new
    while new.any():
        dose = new.astype(float) @ Xe
        p = -np.expm1(-dose)
        nxt = (rng.random((R, n)) < p) & ~inf
        inf |= nxt
        new = nxt
    return inf


def tune_eps(X, rng, target=79.0, thresh=40):
    lo, hi = np.log(1.0 / X.sum(axis=1).mean()) - 2.0, np.log(1.0 / X.sum(axis=1).mean()) + 12.0
    for it in range(14):
        mid = 0.5 * (lo + hi)
        fs = percolate(np.exp(mid) * X, rng, R_TUNE).sum(axis=1)
        big = fs[fs >= thresh]
        m = big.mean() if len(big) >= 5 else 0.0
        if m < target:
            lo = mid
        else:
            hi = mid
    return float(np.exp(0.5 * (lo + hi)))


def kernel(cfg, model, D):
    lat = cfg["lat"]
    w = L.steady_weights(lat, model, D, cfg["kappa"], T=8.0)
    Ph = lat.Phi[cfg["sites"]]
    X = (Ph * w) @ Ph.T
    X = np.clip(X, 0.0, None)
    np.fill_diagonal(X, 0.0)
    return X


def run_model(model, A, rng):
    post = P.posterior_D("room", model, n=400, seed=C.SEED + 7) if model != "M0" else None
    zs, sizes, rows, zlist = [], [], [], []
    n_par = 1 if model == "M0" else N_PAR
    for i in range(n_par):
        cfg = E.build_O1(rng)
        D = 1.0 if model == "M0" else float(post["D"][i])
        X = kernel(cfg, model, D)
        eps = tune_eps(X, rng)
        inf = percolate(eps * X, rng, R_PROD * (20 if model == "M0" else 1))
        fs = inf.sum(axis=1)
        keep = (fs >= 69) & (fs <= 89)
        z = [S.joincount_z(A, c)[0] for c in inf[keep]]
        zs.extend(z)
        zlist.append(np.array(z))
        sizes.extend(fs[keep].tolist())
        rows.append(dict(D=D, kappa=cfg["kappa"], pitch=cfg["pitch"], eps=eps, n_keep=int(keep.sum()),
                         z_mean=float(np.mean(z)) if len(z) else None,
                         frac_big=float((fs >= 40).mean())))
        print(model, i, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()},
              flush=True)
    return np.array(zs), rows, zlist


def main():
    xy, case = E.o1_desks()
    A = E.o1_adjacency(xy, 1.25)
    out = {"stamp": C.stamp(), "n_desks": int(len(case)), "n_pairs": int(A.sum() / 2)}
    rng = np.random.default_rng(C.SEED + 61)
    pred, per_draw = {}, {}
    for model in ("M0", "M1", "M2"):
        z, rows, zlist = run_model(model, A, rng)
        pred[model] = z
        per_draw[model] = zlist
        out[model] = dict(n=int(len(z)), q=np.quantile(z, [0.025, 0.05, 0.5, 0.95, 0.975]).tolist(),
                          mean=float(z.mean()), sd=float(z.std()), rows=rows,
                          z_sample=np.round(z[:: max(1, len(z) // 4000)], 3).tolist())
        C.save_json("06_callcentre_pred.json", out)
    # ---- comparison with the data (after the predictive distributions are stored)
    z_obs, jj, m, sd = S.joincount_z(A, case)
    out["observed"] = dict(z=float(z_obs), case_case_pairs=float(jj), null_mean=float(m), null_sd=float(sd),
                           cases=int(case.sum()))
    for model in ("M0", "M1", "M2"):
        z = pred[model]
        lo, hi = np.quantile(z, [0.025, 0.975])
        out[model]["P_z_le_obs"] = float((z <= z_obs).mean())
        out[model]["inside_95"] = bool(lo <= z_obs <= hi)
        # equal weight per posterior draw (instead of pooling the retained outbreaks)
        fr = np.array([(zz <= z_obs).mean() for zz in per_draw[model] if len(zz)])
        out[model]["P_z_le_obs_equal_weight"] = float(fr.mean())
        out[model]["inside_95_equal_weight"] = bool(0.025 <= fr.mean() <= 0.975)
        out[model]["n_draws_with_P_gt_0.025"] = int((fr > 0.025).sum())
        out[model]["n_draws"] = int(len(fr))
        for r, f in zip(out[model]["rows"], fr):
            r["P_z_le_obs"] = float(f)
    C.save_json("06_callcentre.json", out)
    print("observed", out["observed"])
    for model in ("M0", "M1", "M2"):
        print(model, np.round(out[model]["q"], 2), "P(z<=obs)", out[model]["P_z_le_obs"],
              "inside95", out[model]["inside_95"], "| equal weight:", out[model]["P_z_le_obs_equal_weight"],
              out[model]["n_draws_with_P_gt_0.025"], "of", out[model]["n_draws"])


if __name__ == "__main__":
    main()
