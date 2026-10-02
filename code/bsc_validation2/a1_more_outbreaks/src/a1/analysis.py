"""Calibration, prediction, decision statistics from the predictive tables."""
from __future__ import annotations
import json, os
import numpy as np
from . import events as E, engine as G

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
OUT = os.path.join(ROOT, "out")
MODELS = ["M0", "CRR", "M2", "M1", "M3"]
_cache = {}


def tab(ev, tag="primary"):
    key = (ev, tag)
    if key not in _cache:
        z = np.load(os.path.join(OUT, "tab_%s_%s.npz" % (ev, tag)))
        _cache[key] = {k: z[k] for k in z.files}
    return _cache[key]


def flat(t, model):
    a = t[model]
    return a.reshape(-1, a.shape[-1])            # (n_theta, n_comp)


def comp_index(t, obs):
    return int(np.flatnonzero((t["comps"] == np.asarray(obs)).all(axis=1))[0])


def grid_labels(model):
    if model == "M0":
        return [dict()]
    if model == "CRR":
        return [dict(rho=float(10 ** x)) for x in G.LOGRHO]
    if model in ("M2", "M1"):
        return [dict(D=float(10 ** x)) for x in G.LOGD]
    return [dict(D=float(10 ** x), w=float(w)) for x in G.LOGD3 for w in G.WGRID]


def event_curve(ev, model, obs, tag="primary"):
    """log P(obs | K, theta) on the model's grid."""
    t = tab(ev, tag)
    return np.log(np.clip(flat(t, model)[:, comp_index(t, obs)], 1e-300, None))


def posterior(model, cal_events, obs, tags=None):
    """Grid posterior (flat prior) from the calibration events. obs: dict ev -> bin counts."""
    tags = tags or {}
    ll = sum(event_curve(e, model, obs[e], tags.get(e, "primary")) for e in cal_events)
    w = np.exp(ll - ll.max())
    return w / w.sum(), ll


def predictive(ev, model, w, tag="primary"):
    return w @ flat(tab(ev, tag), model)


def two_sided_p(pmf, i):
    return float(pmf[pmf <= pmf[i] * (1 + 1e-9) + 1e-15].sum())


def summarize(model, w):
    lab = grid_labels(model)
    out = dict(map=lab[int(np.argmax(w))])
    if model in ("M2", "M1", "CRR"):
        x = G.LOGRHO if model == "CRR" else G.LOGD
        c = np.cumsum(w)
        q = lambda a: float(10 ** np.interp(a, c, x))
        out.update(q05=q(0.05), q50=q(0.5), q95=q(0.95))
    if model == "M3":
        W = w.reshape(len(G.LOGD3), len(G.WGRID))
        cw = np.cumsum(W.sum(0)); cd = np.cumsum(W.sum(1))
        out.update(w_q05=float(np.interp(0.05, cw, G.WGRID)), w_q50=float(np.interp(0.5, cw, G.WGRID)),
                   w_q95=float(np.interp(0.95, cw, G.WGRID)),
                   D_q05=float(10 ** np.interp(0.05, cd, G.LOGD3)), D_q50=float(10 ** np.interp(0.5, cd, G.LOGD3)),
                   D_q95=float(10 ** np.interp(0.95, cd, G.LOGD3)))
    return out


def near_far(sizes, counts, near_bins):
    nb = np.isin(np.arange(len(sizes)), near_bins)
    s, c = np.asarray(sizes), np.asarray(counts)
    return int(c[nb].sum()), int(s[nb].sum()), int(c[~nb].sum()), int(s[~nb].sum())


def expected_counts(t, pmf):
    return pmf @ t["comps"]


def delta_reference(pm_a, pm_b, pm_truth, n=200000, seed=0):
    """Samples of sum_e [log pm_a(k_e) - log pm_b(k_e)] with k_e ~ pm_truth_e (lists over events)."""
    rng = np.random.default_rng(seed)
    tot = np.zeros(n)
    for a, b, tr in zip(pm_a, pm_b, pm_truth):
        idx = rng.choice(len(tr), size=n, p=tr / tr.sum())
        tot += np.log(np.clip(a[idx], 1e-300, None)) - np.log(np.clip(b[idx], 1e-300, None))
    return tot


def sample_idx(pm_truth, n, seed):
    rng = np.random.default_rng(seed)
    return [rng.choice(len(tr), size=n, p=tr / tr.sum()) for tr in pm_truth]
