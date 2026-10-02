"""Posterior draws from the calibration files and hold-out shape predictions."""
from __future__ import annotations

import json
import os

import numpy as np

from . import events as E
from . import lattice as L
from . import stats as S

DATA = os.path.join(E.ROOT, "data")
CAL_FILE = {"vehicle": "01_cal_T4_{v}.json", "room": "01_cal_R1_{v}.json"}


def load_cal(cls: str, variant: str = "primary"):
    with open(os.path.join(DATA, CAL_FILE[cls].format(v=variant))) as f:
        return json.load(f)


def posterior_D(cls: str, model: str, variant: str = "primary", n: int = 400,
                seed: int = 0, cal=None):
    """Draws of D (m^2/h) from the profile-likelihood posterior (flat prior on
    log10 D); for the vehicle class and M2 also the paired kappa_T4 and eps."""
    cal = load_cal(cls, variant) if cal is None else cal
    rng = np.random.default_rng(seed)
    logD = np.array(cal["logD"])
    step = logD[1] - logD[0]
    m = cal[model]
    if cls == "vehicle" and model == "M2":
        ll = np.array(m["profile_logL"])                       # (D, kappa)
        w = np.exp(ll - ll.max())
        w = w / w.sum()
        flat = rng.choice(w.size, size=n, p=w.ravel())
        i, j = np.unravel_index(flat, w.shape)
        D = 10 ** (logD[i] + rng.uniform(-step / 2, step / 2, n))
        return dict(D=D, kappa_cal=np.array(cal["kappa_nodes"])[j],
                    eps_cal=np.array(m["eps"])[i, j], w_logD=w.sum(axis=1), logD=logD)
    ll = np.array(m["profile_logL"])
    w = np.exp(ll - ll.max())
    w = w / w.sum()
    i = rng.choice(len(w), size=n, p=w)
    D = 10 ** (logD[i] + rng.uniform(-step / 2, step / 2, n))
    return dict(D=D, eps_cal=np.array(m["eps"])[i], w_logD=w, logD=logD)


def summarize_posterior(w, logD):
    c = np.cumsum(w)
    q = lambda a: float(10 ** np.interp(a, c, logD))
    return dict(mode=float(10 ** logD[int(np.argmax(w))]), q05=q(0.05), q50=q(0.5), q95=q(0.95))


def shape_predictive(ev: str, model: str, Ddraws, n_draw: int = 300, seed: int = 0,
                     split: str = "primary", builder_kwargs=None, K=None,
                     fixed_D=None, theta_draws=None, nu=None):
    """Predictive distribution of the number of cases in the near stratum,
    conditional on the event total.

    Returns dict(pmf, n_near, n_far, K, rr_expected (per draw), p_near, p_far).
    """
    rng = np.random.default_rng(seed)
    builder = E.BUILDERS[ev]
    kw = builder_kwargs or {}
    pmf = None
    rr_e, pn, pf, eps_l = [], [], [], []
    for i in range(n_draw):
        cfg = builder(rng, **kw)
        Kc = cfg["K"] if K is None else K
        if split == "primary":
            near, valid = cfg["near"], np.ones(len(cfg["rec"]), bool)
        else:
            near = cfg["alt"][split]
            valid = cfg["alt"].get("side_valid", np.ones(len(cfg["rec"]), bool)) \
                if split == "driver_side" else np.ones(len(cfg["rec"]), bool)
        D = float(Ddraws[i % len(Ddraws)]) if fixed_D is None else fixed_D
        lat = cfg["lat"]
        if theta_draws is not None:          # S-aniso: seat-back factor on the along-cabin hop rate
            lat = L.Lattice(lat.nx, lat.ny, lat.ax, lat.ay, lat.H,
                            theta_y=float(theta_draws[i % len(theta_draws)]))
        X = L.exposure(lat, model, D, cfg["kappa"], cfg["segments"], cfg["src"], cfg["rec"])
        X = np.clip(X, 1e-300, None)
        if nu is not None:                   # S-nu: near-field factor on seats adjacent to the index
            sx, sy = cfg["src"] % lat.nx, cfg["src"] // lat.nx
            rx, ry = cfg["rec"] % lat.nx, cfg["rec"] // lat.nx
            X = np.where((ry == sy) & (np.abs(rx - sx) == 1), nu * X, X)
        eps = S.profile_eps(X, Kc)
        p = -np.expm1(-eps * X)
        a, b = p[near & valid], p[~near & valid]
        pm = S.cond_pmf(a, b, Kc)
        pmf = pm if pmf is None else pmf + pm
        rr_e.append(a.mean() / max(b.mean(), 1e-300))
        pn.append(a.mean())
        pf.append(b.mean())
        eps_l.append(eps)
    pmf = pmf / n_draw
    return dict(pmf=pmf, n_near=int(len(a)), n_far=int(len(b)), K=int(Kc),
                rr_expected=np.array(rr_e), p_near=np.array(pn), p_far=np.array(pf),
                eps=np.array(eps_l))


def rr_predictive_interval(pmf, n1, n2, K, level=0.90):
    """Predictive interval of the risk ratio implied by the conditional pmf."""
    k = np.arange(len(pmf))
    r = S.rr(k, n1, K - k, n2)
    order = np.argsort(r)
    c = np.cumsum(pmf[order])
    a = (1 - level) / 2
    lo = r[order][np.searchsorted(c, a)]
    hi = r[order][min(np.searchsorted(c, 1 - a), len(c) - 1)]
    med = r[order][min(np.searchsorted(c, 0.5), len(c) - 1)]
    return float(lo), float(med), float(hi)
