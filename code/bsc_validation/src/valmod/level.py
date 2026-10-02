"""Level (attack-rate) calculations: Wells-Riley structural factors and implied
emission rates.  For a group spread over a room the spatial mean of the M2
exposure equals the M0 value, so M0 and M2 share these formulas."""
from __future__ import annotations

import numpy as np

from . import events as E

# volume (m^3), floor area (m^2), exposure segments (h), number exposed
LEVEL_EVENTS = {
    "H1": dict(V=810.0, area=180.0, segments=[2.5], label="Skagit choir"),
    "O2": dict(V=700.0 * 2.7, area=700.0, segments=[11.0], label="Zurich office"),
    "R2": dict(V=96.6 * 3.2, area=96.6, segments=[21 / 60], label="Jeonju restaurant"),
    "C5": dict(V=60.0 * 3.0, area=60.0, segments=[50 / 60, 50 / 60], label="Cheonan fitness classes",
               n_index=8),
    "T2": dict(V=60.42, area=11.4 * 2.5, segments=[200 / 60], label="Hunan coach"),
    "T3": dict(V=21.69, area=5.5 * 2.5, segments=[1.0], label="Hunan minibus"),
}


def wr_factor(V, kappa, segments, breath=E.BREATH):
    """Dose per unit emission rate (h / quantum): b/V * sum_T [T/k - (1-e^{-kT})/k^2]."""
    kappa = np.asarray(kappa, float)
    s = 0.0
    for T in segments:
        s = s + (T / kappa - (-np.expm1(-kappa * T)) / kappa ** 2)
    return breath / V * s


def implied_dose(k, n):
    return -np.log1p(-k / n)


def lognormal_from_mean(mean, sigma):
    """mu of a lognormal with given arithmetic mean and log-sd."""
    return np.log(mean) - 0.5 * sigma ** 2
