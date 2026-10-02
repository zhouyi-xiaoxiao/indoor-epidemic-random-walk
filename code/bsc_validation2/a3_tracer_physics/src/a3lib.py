"""Shared code of the A3 analysis (physics-calibrated non-local kernel).

Imports the lattice propagator, the event builders and the statistics of the
first validation read-only (``../../bsc_validation/src/valmod``).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VAL1 = os.path.normpath(os.path.join(ROOT, "..", "..", "bsc_validation"))
OUTB = os.path.normpath(os.path.join(ROOT, "..", "..", "bsc_outbreaks"))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(ROOT, ".mplcache"))
sys.path.insert(0, os.path.join(VAL1, "src"))
from valmod import events as E          # noqa: E402
from valmod import lattice as L         # noqa: E402
from valmod import stats as S           # noqa: E402

SEED = 20261001
DER = os.path.join(ROOT, "data", "derived")
RES = os.path.join(ROOT, "results")
LOSS = E.LOSS
LOGD = np.round(np.arange(-1.0, 5.0001, 0.05), 2)       # pre-registered grid 0.1 .. 1e5 m^2/h
DGRID = 10.0 ** LOGD


def save_json(name, obj):
    p = os.path.join(RES, name)
    with open(p + ".part", "w") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    os.replace(p + ".part", p)
    return p


def load_json(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def stamp():
    return time.strftime("%Y-%m-%d %H:%M:%S %Z")


# ---------------------------------------------------------------------------
# steady-state M2 field and the tracer fit
# ---------------------------------------------------------------------------
def steady_conc(lat, D, kappa, src, rec):
    """Steady concentration (arbitrary units) at sites rec for a unit source at src."""
    w = L.steady_weights(lat, "M2", D, kappa)
    return np.clip(lat.Phi[np.asarray(rec)] @ (w * lat.Phi[src]), 1e-300, None)


def fit_D(lat, src, rec, c_obs, kappa, weights=None):
    """Least squares of log model concentration on log tracer value with a free
    additive constant, on the pre-registered D grid.  Returns (D_hat, sse array)."""
    y = np.log(np.asarray(c_obs, float))
    w = np.ones(len(y)) if weights is None else np.asarray(weights, float)
    sse = np.empty(len(DGRID))
    for i, D in enumerate(DGRID):
        m = np.log(steady_conc(lat, D, kappa, src, rec))
        r = y - m
        r = r - np.average(r, weights=w)
        sse[i] = np.sum(w * r * r)
    return float(DGRID[int(np.argmin(sse))]), sse


def boot_D(lat, src, rec, c_obs, kappa_draws, rng, n=200):
    rec = np.asarray(rec)
    c_obs = np.asarray(c_obs, float)
    out = np.empty(n)
    for b in range(n):
        ix = rng.integers(len(rec), size=len(rec))
        while len(np.unique(rec[ix])) < 2:
            ix = rng.integers(len(rec), size=len(rec))
        out[b] = fit_D(lat, src, rec[ix], c_obs[ix], float(kappa_draws[b % len(kappa_draws)]))[0]
    return out


# ---------------------------------------------------------------------------
# dose-response and conditional distributions
# ---------------------------------------------------------------------------
def probs(X, K, dose="exp"):
    """Infection probabilities with the level parameter profiled to the total K."""
    X = np.clip(np.asarray(X, float), 1e-300, None)
    n = len(X)
    if dose == "exp":
        eps = S.profile_eps(X, K)
        return -np.expm1(-eps * X), eps
    if dose == "lin":                      # low-dose limit, p proportional to X (capped at 1)
        f = lambda le: np.sum(np.minimum(1.0, np.exp(le) * X)) - K
    elif dose == "gam":                    # exponentially distributed susceptibility
        f = lambda le: np.sum(np.exp(le) * X / (1.0 + np.exp(le) * X)) - K
    else:
        raise ValueError(dose)
    lo = np.log(K / X.sum()) - 2.0
    hi = lo + 4.0
    while f(hi) < 0:
        hi += 2.0
    while f(lo) > 0:
        lo -= 2.0
    eps = float(np.exp(optimize.brentq(f, lo, hi, xtol=1e-12)))
    p = np.minimum(1.0, eps * X) if dose == "lin" else eps * X / (1.0 + eps * X)
    return p, eps


def near_pmf(X, near, K, dose="exp"):
    p, eps = probs(X, K, dose)
    near = np.asarray(near, bool)
    return S.cond_pmf(p[near], p[~near], K), p, eps


def crr_pmf(n1, n2, K, rho):
    """Constant near/far risk ratio rho, expected total K."""
    q = K / (n1 * rho + n2)
    if rho * q > 1.0:
        q = (K - n1) / n2
        pn = 1.0
    else:
        pn = rho * q
    return S.cond_pmf(np.full(n1, pn), np.full(n2, q), K)


def m0_pmf(n1, n2, K):
    return S.hypergeom_pmf(n1, n2, K)


def two_sided_p(pmf, k):
    return S.two_sided_p(np.asarray(pmf), int(k))


def expected_rr(pmf, n1, n2, K):
    k = np.arange(len(pmf))
    e1 = float((pmf * k).sum())
    return (e1 / n1) / max((K - e1) / n2, 1e-300)


def seat_loglik(p, case_mask, K):
    """log P(this set of K cases | exactly K cases) for independent Bernoulli(p)."""
    p = np.clip(np.asarray(p, float), 1e-300, 1 - 1e-15)
    odds = p / (1 - p)
    e = np.zeros(K + 1)                     # elementary symmetric polynomials of the odds
    e[0] = 1.0
    for o in odds:
        e[1:] = e[1:] + o * e[:-1]
    return float(np.sum(np.log(odds[np.asarray(case_mask, bool)])) - np.log(e[K]))


# ---------------------------------------------------------------------------
# event geometries with explicit seat maps
# ---------------------------------------------------------------------------
T2_COLS = {"A": 0, "B": 1, "E": 2, "C": 3, "D": 4}
T2_CASES = ["1D", "5C", "6A", "6D", "9C", "9D", "13D"]


def t2_event():
    """Hunan coach from Ou et al. Fig. S4A.  Lattice of the first validation."""
    df = pd.read_csv(os.path.join(DER, "ou2022_B1_seats.csv"))
    df = df[df["seat"].str.fullmatch(r"\d+[A-E]")].copy()
    lat = L.Lattice(nx=5, ny=13, ax=0.5, ay=11.4 / 13, H=60.42 / (11.4 * 2.5))
    pitch = (df.loc[df.seat == "1D", "x_fig"].iloc[0] - df.loc[df.seat == "13D", "x_fig"].iloc[0]) / 12.0
    x1 = df.loc[df.seat == "1D", "x_fig"].iloc[0]
    df["lrow"] = np.clip(np.rint((x1 - df["x_fig"]) / pitch), 0, 12).astype(int)
    df["lcol"] = df["col"].map(T2_COLS)
    df["site"] = lat.site(df["lcol"].to_numpy(), df["lrow"].to_numpy())
    assert df["site"].is_unique
    sus = df[(df.seat != "12D") & (df.seat != "8C")].reset_index(drop=True)
    assert len(sus) == 45 and sus["s4_cfd"].notna().all()
    near = (sus["row"] >= 8).to_numpy()
    side_valid = (sus["seat"] != "13E").to_numpy()
    driver_side = sus["col"].isin(["C", "D"]).to_numpy()
    case = sus["seat"].isin(T2_CASES).to_numpy()
    assert near.sum() == 19 and case.sum() == 7 and (driver_side & side_valid).sum() == 24
    return dict(lat=lat, df=df, sus=sus, near=near, case=case, driver_side=driver_side, side_valid=side_valid,
                src=int(df.loc[df.seat == "12D", "site"].iloc[0]), src_test=int(df.loc[df.seat == "12C", "site"].iloc[0]),
                K=7)


T5_OCC = {"K": [2, 3, 4, 5], "G": [1, 3, 4, 5, 6, 7], "D": [1, 3, 4, 5, 6, 7], "A": [1, 2, 3, 4, 5]}
T5_CASES = ["3K", "4K", "4G", "5G", "6G", "7G", "3D", "4D", "5D", "6D", "7D", "5A"]
T5_NEAR = ["3K", "4K", "3G", "4G", "5G", "6G", "7G", "3D", "4D", "5D", "6D", "7D"]
T5_COLIDX = {"A": 0, "D": 2, "G": 3, "K": 5}


def t5_event():
    """VN54 business class, seat map digitised from Khanh et al. Fig. 1."""
    lat = L.Lattice(nx=6, ny=46, ax=0.9, ay=1.1, H=2.2)
    seats = [f"{r}{c}" for c, rows in T5_OCC.items() for r in rows if f"{r}{c}" != "5K"]
    row = np.array([int(s[:-1]) for s in seats])
    col = np.array([s[-1] for s in seats])
    site = lat.site(np.array([T5_COLIDX[c] for c in col]), row - 1 + E.T5_ROW0)
    near = np.isin(seats, T5_NEAR)
    case = np.isin(seats, T5_CASES)
    assert len(seats) == 20 and near.sum() == 12 and case.sum() == 12 and (near & case).sum() == 11
    return dict(lat=lat, seats=seats, row=row, col=col, rec=site, near=near, case=case,
                src=int(lat.site(T5_COLIDX["K"], 5 - 1 + E.T5_ROW0)), K=12)


def kinahan_maps(zero_missing=False, section="FWD", airframe="777"):
    """Mean integrated count per sensor seat for each unmasked breathing release."""
    df = pd.read_csv(os.path.join(DER, "kinahan_tidy.csv"))
    g = df[(df.airframe.astype(str) == airframe) & (df.section == section) & (df.kind == "B") & (df["mask"] == "No")]
    g = g[g["col"].notna() & (g["col"] != "")].copy()
    if zero_missing:
        g.loc[g["count"] == 0, "count"] = np.nan
    out = {}
    for rel, d in g.groupby("release"):
        out[rel] = d.groupby(["row", "col"])["count"].mean().unstack()
    return out


def t5_direct_kernels(zero_missing=False):
    """K-A kernels for an index in 5K: release 5L as is, release 5A mirrored.
    Returns dict name -> (values per VN54 susceptible seat, number imputed)."""
    ev = t5_event()
    maps = kinahan_maps(zero_missing)
    mirror = {"A": "L", "D": "G", "G": "D", "L": "A"}
    out = {}
    for name, rel, mir in (("5L", "5L", False), ("5A-mirrored", "5A", True)):
        m = maps[rel]
        def get(r, c):
            cc = "L" if c == "K" else c
            if mir:
                cc = mirror[cc]
            try:
                return m.loc[r, cc]
            except KeyError:
                return np.nan
        vals, nimp = [], 0
        vmirror = {"A": "K", "K": "A", "D": "G", "G": "D"}          # VN54 column letters
        for r, c in zip(ev["row"], ev["col"]):
            v = get(r, c)
            if not np.isfinite(v):
                nimp += 1
                v = get(r, vmirror[c])                               # mirror seat, same row, same map
                if not np.isfinite(v):
                    rowvals = np.array([get(r, x) for x in ("A", "D", "G", "K")], float)
                    v = np.nanmean(rowvals) if np.isfinite(rowvals).any() else np.nan
                if not np.isfinite(v):
                    adj = np.array([get(rr, x) for rr in (r - 1, r + 1) for x in ("A", "D", "G", "K")], float)
                    v = np.nanmean(adj)
            vals.append(v)
        out[name] = (np.array(vals, float), nimp)
    return ev, out


# ---------------------------------------------------------------------------
# R1 restaurant
# ---------------------------------------------------------------------------
def r1_tables():
    df = pd.read_csv(os.path.join(OUTB, "data", "li2021_restaurant_tables.csv"))
    df = df[(df["table"] != "TA") & (df["patrons"] > 0)].reset_index(drop=True)
    df["tracer"] = df["norm_measured_tracer"].fillna(df["norm_predicted_tracer"])
    df["T_h"] = df["overlap_with_table_A_min"] / 60.0
    return df
