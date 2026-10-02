"""End-to-end cross-check of one hold-out prediction by a separate route.

T2 (Hunan coach), M2 with D = 25 m^2/h and one seating configuration:
  * exposure by brute-force quadrature of the matrix exponential of the lattice
    generator (no mode expansion);
  * conditional distribution of the rear-stratum count given 7 cases by Monte
    Carlo rejection sampling of independent Bernoulli outcomes (no convolution).
Same for T5 (flight) and C1 (classroom) at the posterior-mode D of their class.
"""
import numpy as np
from scipy.linalg import expm

import _common as C
from valmod import events as E
from valmod import lattice as L
from valmod import predict as P
from valmod import stats as S

out = {"stamp": C.stamp(), "events": {}}
xs, ws = np.polynomial.legendre.leggauss(80)
for ev, cls in (("T2", "vehicle"), ("T5", "vehicle"), ("C1", "room")):
    rng = np.random.default_rng(C.SEED + 5)
    cfg = E.BUILDERS[ev](rng)
    D = P.load_cal(cls)["M2"]["D_hat"]
    lat, kap = cfg["lat"], cfg["kappa"]
    Q = L.generator(lat, D)
    e0 = np.zeros(lat.M)
    e0[cfg["src"]] = 1.0
    Xb = np.zeros(lat.M)
    for T in cfg["segments"]:
        # split the time integral to resolve the fast initial transient
        edges = np.concatenate([[0.0], np.geomspace(min(1e-3, T / 10), T, 12)])
        for a, b in zip(edges[:-1], edges[1:]):
            t = 0.5 * (b - a) * (xs + 1) + a
            w = 0.5 * (b - a) * ws
            for ti, wi in zip(t, w):
                Xb += wi * (T - ti) * np.exp(-kap * ti) * (expm(Q * ti) @ e0)
    Xs = L.exposure(lat, "M2", D, kap, cfg["segments"], cfg["src"], cfg["rec"])
    err = float(np.abs(Xb[cfg["rec"]] - Xs).max() / Xs.max())
    K = cfg["K"]
    eps = S.profile_eps(Xs, K)
    p = -np.expm1(-eps * Xs)
    pmf = S.cond_pmf(p[cfg["near"]], p[~cfg["near"]], K)
    draws = rng.random((400000, len(p))) < p
    sel = draws.sum(axis=1) == K
    kn = draws[sel][:, cfg["near"]].sum(axis=1)
    mc = np.bincount(kn, minlength=len(pmf)) / sel.sum()
    out["events"][ev] = dict(D=D, kappa=kap, rel_err_exposure=err, n_accepted=int(sel.sum()),
                             max_abs_diff_pmf=float(np.abs(mc - pmf).max()),
                             mc_se_max=float(np.sqrt(0.25 / sel.sum())))
    print(ev, out["events"][ev], flush=True)
C.save_json("08_bruteforce_check.json", out)
