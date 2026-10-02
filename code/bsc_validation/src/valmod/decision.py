"""Decision-rule computations (pre-registration section 8)."""
from __future__ import annotations

import numpy as np

from . import stats as S


def pooled_delta_test(pm_model, pm_null, k_obs):
    """Exact one-sided test of the pooled log-score difference.

    pm_model, pm_null: lists of pmfs (one per event) of k_near | K.
    Returns Delta_obs, p = P_null(Delta >= Delta_obs), per-event deltas.
    """
    pm_model = [np.clip(np.asarray(p, float), 1e-300, None) for p in pm_model]
    pm_null = [np.clip(np.asarray(p, float), 1e-300, None) for p in pm_null]
    sizes = [len(p) for p in pm_null]
    grids = np.meshgrid(*[np.arange(s) for s in sizes], indexing="ij")
    delta = np.zeros(sizes)
    w0 = np.ones(sizes)
    for a, g in enumerate(grids):
        delta += np.log(pm_model[a][g]) - np.log(pm_null[a][g])
        w0 = w0 * pm_null[a][g]
    w0 = w0 / w0.sum()
    per = [float(np.log(pm_model[a][k]) - np.log(pm_null[a][k])) for a, k in enumerate(k_obs)]
    d_obs = float(sum(per))
    p = float(w0[np.round(delta, 9) >= round(d_obs, 9)].sum())
    return d_obs, p, per


def event_row(pmf, n1, n2, K, k_obs):
    pmf = np.asarray(pmf, float)
    pmf = pmf / pmf.sum()
    k = np.arange(len(pmf))
    ek = float((k * pmf).sum())
    rr_pred = (ek / n1) / max((K - ek) / n2, 1e-12)
    order = np.argsort(S.rr(k, n1, K - k, n2))
    r_sorted = S.rr(k, n1, K - k, n2)[order]
    c = np.cumsum(pmf[order])
    lo = float(r_sorted[np.searchsorted(c, 0.05)])
    hi = float(r_sorted[min(np.searchsorted(c, 0.95), len(c) - 1)])
    klo, khi = S.pmf_interval(pmf, 0.90)
    rr_obs = float(S.rr(k_obs, n1, K - k_obs, n2))
    se = S.se_log_rr(k_obs, n1, K - k_obs, n2)
    return dict(log_score=float(np.log(max(pmf[k_obs], 1e-300))), p_two_sided=S.two_sided_p(pmf, k_obs),
                k_pred_mean=ek, k_pred_90=[int(klo), int(khi)], rr_pred=float(rr_pred),
                rr_pred_90=[lo, hi], rr_obs=rr_obs, se_log_rr_obs=se,
                z_original=float((np.log(rr_pred) - np.log(rr_obs)) / se),
                pmf_at_obs=float(pmf[k_obs]))


def label(p_model, p_null):
    if p_model < 0.05:
        return "not reproduced"
    return "reproduced, discriminating" if p_null < 0.05 else "reproduced, non-discriminating"
