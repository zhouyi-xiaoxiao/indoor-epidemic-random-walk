"""Observation models of the outbreak records (pre-registration section 4).

Every builder returns, for one random draw, a dict with

    lat       Lattice
    src       site index of the index case
    rec       site indices of the susceptibles
    near      boolean mask over ``rec`` (primary split)
    alt       dict of alternative masks (secondary splits), may be empty
    segments  list of exposure durations in hours
    kappa     removal rate (1/h) drawn from the prior of the event
    K         total number of secondary cases in the event (conditioning variable)

No stratum-level outcome is stored here; hold-out stratum counts live in
``data/heldout_outcomes.json`` and are read by the evaluation scripts only.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from .lattice import Lattice

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
OUTB = os.path.normpath(os.path.join(ROOT, "..", "bsc_outbreaks"))

LOSS = 0.93            # deposition 0.3/h + inactivation 0.63/h
BREATH = 0.5           # m^3/h reference breathing rate
SEED = 20261001

# ventilation priors (air changes per hour); kappa = ACH + LOSS
VENT = {
    "T4": ("loguniform", 5.0, 20.0),
    "T1": ("loguniform", 0.5, 5.0),
    "T2": ("lognormal", 4.8, 0.2),
    "T3": ("lognormal", 9.6, 0.2),
    "T5": ("loguniform", 10.0, 30.0),
    "R1": ("uniform", 0.56, 0.77),
    "C1": ("loguniform", 2.0, 12.0),
    "O1": ("loguniform", 1.0, 6.0),
    "H1": ("uniform", 0.3, 1.0),
    "O2": ("lognormal", 1.0, 0.3),
    "R2": ("loguniform", 0.1, 0.5),
    "C5": ("loguniform", 0.5, 5.0),
    "C3": ("loguniform", 1.0, 6.0),
    "S1": ("loguniform", 1.0, 4.0),
    "T7": ("loguniform", 5.0, 20.0),
}


def draw_kappa(ev: str, rng, size=None):
    kind, a, b = VENT[ev]
    if kind == "loguniform":
        v = np.exp(rng.uniform(np.log(a), np.log(b), size))
    elif kind == "uniform":
        v = rng.uniform(a, b, size)
    else:
        v = a * np.exp(b * rng.standard_normal(size))
    return v + LOSS


def kappa_quadrature(ev: str, n: int = 7):
    """Nodes and weights of the ventilation prior (midpoint rule in the
    natural scale of the prior)."""
    kind, a, b = VENT[ev]
    u = (np.arange(n) + 0.5) / n
    if kind == "loguniform":
        v = np.exp(np.log(a) + u * (np.log(b) - np.log(a)))
    elif kind == "uniform":
        v = a + u * (b - a)
    else:
        from scipy.stats import norm
        v = a * np.exp(b * norm.ppf(u))
    return v + LOSS, np.full(n, 1.0 / n)


# ---------------------------------------------------------------------------
# T4 high-speed trains (calibration, vehicle class)
# ---------------------------------------------------------------------------
T4_ROWS = 17
T4_SEATCOLS = np.array([0, 1, 2, 4, 5])          # A B C | aisle(3) | D F
T4_T_MEAN, T4_T_SD = 2.1, 1.8


def t4_lattice(row_pitch: float = 1.0, theta_y: float = 1.0) -> Lattice:
    return Lattice(nx=6, ny=T4_ROWS, ax=0.5, ay=row_pitch, H=2.4, theta_y=theta_y)


def t4_seats(lat: Lattice):
    rows = np.repeat(np.arange(T4_ROWS), len(T4_SEATCOLS))
    cols = np.tile(T4_SEATCOLS, T4_ROWS)
    return lat.site(cols, rows), rows, cols


def t4_duration_nodes(n: int = 16):
    """Equal-probability nodes of the gamma co-travel-time distribution."""
    from scipy.stats import gamma
    shape = (T4_T_MEAN / T4_T_SD) ** 2
    scale = T4_T_SD ** 2 / T4_T_MEAN
    u = (np.arange(n) + 0.5) / n
    return gamma.ppf(u, shape, scale=scale), np.full(n, 1.0 / n)


def t4_data():
    """Reconstructed cell counts (verifier).  Row-0 cells are printed one column
    to the left on the publisher page: printed_col c means c + 1 seats apart."""
    df = pd.read_csv(os.path.join(OUTB, "verify", "results",
                                  "train_matrix_reconstruction.csv"))
    df = df.rename(columns={"row": "dr"})
    df["dc"] = np.where(df["dr"] == 0, df["printed_col"] + 1, df["printed_col"])
    pct = pd.read_csv(os.path.join(OUTB, "data", "hu2021_train_attack_matrix.csv"))
    pct = pct[pct["rows_apart"].astype(str).str.isdigit()
              & pct["cols_apart"].astype(str).str.isdigit()].copy()
    pct["dr"] = pct["rows_apart"].astype(int)
    pct["dc"] = pct["cols_apart"].astype(int)
    df = df.merge(pct[["dr", "dc", "attack_rate_pct", "ci95_low_pct", "ci95_high_pct"]],
                  on=["dr", "dc"], how="left")
    assert np.allclose(df["pct"], df["attack_rate_pct"]), "train table alignment"
    return df[["dr", "dc", "k", "n", "pct", "lo", "hi"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# T1 Zhejiang bus
# ---------------------------------------------------------------------------
def build_T1(rng):
    lat = Lattice(nx=6, ny=15, ax=2.5 / 6, ay=0.75, H=2.0)
    seatcols = np.array([0, 1, 2, 4, 5])
    rows = np.repeat(np.arange(15), 5)
    cols = np.tile(seatcols, 15)
    src_mask = (rows == 7) & (cols == 1)                 # row 8, middle of 3-seat side
    zone_a = (rows >= 5) & (rows <= 9)                   # rows 6-10
    zone_b = (rows == 4) | (rows == 10)                  # rows 5 and 11
    zone_c = ~(zone_a | zone_b)
    occ = np.ones(75, bool)
    occ[src_mask] = False
    cand = np.flatnonzero(zone_a & ~src_mask)            # 24 seats, 23 occupied
    occ[rng.choice(cand, 1, replace=False)] = False
    cand = np.flatnonzero(zone_c)                        # 40 seats, 34 occupied
    occ[rng.choice(cand, 6, replace=False)] = False
    rec = lat.site(cols[occ], rows[occ])
    near = (zone_a | zone_b)[occ]
    alt = {"rows6-10": zone_a[occ]}
    assert near.sum() == 33 and (~near).sum() == 34 and alt["rows6-10"].sum() == 23
    return dict(lat=lat, src=int(lat.site(1, 7)), rec=rec, near=near, alt=alt,
                segments=[50 / 60, 50 / 60], kappa=float(draw_kappa("T1", rng)), K=23,
                rows=rows[occ], cols=cols[occ])


# ---------------------------------------------------------------------------
# T2 Hunan coach
# ---------------------------------------------------------------------------
def build_T2(rng, minutes: float = 200.0):
    lat = Lattice(nx=5, ny=13, ax=0.5, ay=11.4 / 13, H=60.42 / (11.4 * 2.5))
    rows, cols = [], []
    for r in range(13):
        for c in (0, 1, 3, 4):
            rows.append(r)
            cols.append(c)
    rows.append(12)
    cols.append(2)                                        # seat 13E on the aisle site
    rows, cols = np.array(rows), np.array(cols)
    is_src = (rows == 11) & (cols == 4)                   # 12D
    is_8c = (rows == 7) & (cols == 3)
    occ = ~(is_src | is_8c)
    opp = cols <= 1
    front = rows <= 6
    occ[rng.choice(np.flatnonzero(opp & front), 2, replace=False)] = False
    occ[rng.choice(np.flatnonzero(opp & ~front), 4, replace=False)] = False
    rec = lat.site(cols[occ], rows[occ])
    near = ~front[occ]                                    # rear rows 8-13
    r_occ, c_occ = rows[occ], cols[occ]
    is13e = (r_occ == 12) & (c_occ == 2)
    alt = {"driver_side": (c_occ >= 3), "side_valid": ~is13e}
    assert near.sum() == 19 and (~near).sum() == 26
    assert (alt["driver_side"] & ~is13e).sum() == 24 and ((c_occ <= 1)).sum() == 20
    return dict(lat=lat, src=int(lat.site(4, 11)), rec=rec, near=near, alt=alt,
                segments=[minutes / 60], kappa=float(draw_kappa("T2", rng)), K=7,
                rows=r_occ, cols=c_occ)


# ---------------------------------------------------------------------------
# T3 minibus (tied to T2)
# ---------------------------------------------------------------------------
def build_T3(rng):
    lat = Lattice(nx=4, ny=6, ax=2.5 / 4, ay=5.5 / 6, H=21.69 / (5.5 * 2.5))
    rows = np.repeat(np.arange(6), 3)
    cols = np.tile(np.array([0, 1, 3]), 6)
    src_i = rng.choice(np.flatnonzero(rows == 3))
    occ = np.ones(18, bool)
    occ[src_i] = False
    rec = lat.site(cols[occ], rows[occ])
    return dict(lat=lat, src=int(lat.site(cols[src_i], rows[src_i])), rec=rec,
                near=np.ones(17, bool), alt={}, segments=[1.0],
                kappa=float(draw_kappa("T3", rng)), K=2)


# ---------------------------------------------------------------------------
# T5 flight VN54
# ---------------------------------------------------------------------------
T5_ROW0 = 2                                               # business rows start here
T5_SEATCOLS = np.array([0, 2, 3, 5])                      # A | aisle | D G | aisle | K
T5_PE_ROWS = np.arange(T5_ROW0 + 7 + 3, T5_ROW0 + 7 + 3 + 5)


def build_T5(rng, index_seat=None):
    lat = Lattice(nx=6, ny=46, ax=0.9, ay=1.1, H=2.2)
    rows = np.repeat(np.arange(7), 4) + T5_ROW0
    cols = np.tile(T5_SEATCOLS, 7)
    if index_seat is None:
        i0 = int(rng.integers(28))
    else:                                                 # (row 1..7, column index 0..3)
        i0 = (index_seat[0] - 1) * 4 + index_seat[1]
    others = np.delete(np.arange(28), i0)
    occ_i = rng.choice(others, 20, replace=False)
    x, y = (cols + 0.5) * lat.ax, (rows + 0.5) * lat.ay
    d = np.hypot(x[occ_i] - x[i0], y[occ_i] - y[i0]) + 1e-6 * rng.random(20)
    order = np.argsort(d)
    near = np.zeros(20, bool)
    near[order[:12]] = True
    rec = lat.site(cols[occ_i], rows[occ_i])
    pe = lat.site(np.tile(np.arange(6), 5), np.repeat(T5_PE_ROWS, 6))
    return dict(lat=lat, src=int(lat.site(cols[i0], rows[i0])), rec=rec, near=near,
                alt={}, segments=[10.0], kappa=float(draw_kappa("T5", rng)), K=12,
                pe_sites=pe, dist=d)


# ---------------------------------------------------------------------------
# R1 Guangzhou restaurant (calibration, seated-room class)
# ---------------------------------------------------------------------------
R1_XS = np.array([1, 4, 7, 10, 12, 15])
R1_YS = np.array([1, 4, 6])


def r1_tables():
    df = pd.read_csv(os.path.join(OUTB, "data", "li2021_restaurant_tables.csv"))
    df["T_h"] = df["overlap_with_table_A_min"] / 60.0
    return df


def build_R1(rng):
    """One random table layout.  Returns per-table arrays (excluding TA and the
    empty T04): site, n patrons, duration, name, class."""
    lat = Lattice(nx=17, ny=8, ax=1.0, ay=1.0, H=3.14)
    df = r1_tables()
    iA = int(rng.integers(1, 5))
    side = int(rng.choice([-1, 1]))
    pos = {"TA": (iA, 0), "TB": (iA + side, 0), "TC": (iA - side, 0), "T18": (iA, 1)}
    free = [(i, j) for i in range(6) for j in range(3) if (i, j) not in pos.values()]
    remote = [t for t in df["table"] if t not in pos]
    perm = rng.permutation(len(free))
    for t, k in zip(remote, perm):
        pos[t] = free[k]
    keep = df[(df["table"] != "TA") & (df["patrons"] > 0)].reset_index(drop=True)
    sites = np.array([lat.site(R1_XS[pos[t][0]], R1_YS[pos[t][1]]) for t in keep["table"]])
    return dict(lat=lat, src=int(lat.site(R1_XS[iA], R1_YS[0])), sites=sites,
                n=keep["patrons"].to_numpy(), T=keep["T_h"].to_numpy(),
                name=keep["table"].to_list(), cls=keep["neighbour_class"].to_list(),
                zone=keep["zone"].to_list(), kappa=float(draw_kappa("R1", rng)))


# ---------------------------------------------------------------------------
# C1 Marin classroom
# ---------------------------------------------------------------------------
def build_C1(rng, corner: bool = False):
    lat = Lattice(nx=11, ny=13, ax=1.0, ay=1.0, H=3.0)
    xs = np.array([1, 3, 5, 7, 9])
    ys = np.array([2, 4, 6, 8, 10])
    cols = np.tile(xs, 5)
    rows = np.repeat(ys, 5)
    front = rows <= 4
    occ = np.ones(25, bool)
    occ[rng.choice(np.flatnonzero(~front), 1)] = False
    rec = lat.site(cols[occ], rows[occ])
    src = lat.site(0, 0) if corner else lat.site(5, 0)
    return dict(lat=lat, src=int(src), rec=rec, near=front[occ], alt={},
                segments=[6.5, 6.5], kappa=float(draw_kappa("C1", rng)), K=12)


# ---------------------------------------------------------------------------
# O1 Seoul call centre, north wing
# ---------------------------------------------------------------------------
def o1_desks():
    df = pd.read_csv(os.path.join(OUTB, "data", "park2020_callcentre_seats_digitized.csv"))
    nw = df[df["region"] == "north_wing"].reset_index(drop=True)
    xy = nw[["x_pitch_units", "y_pitch_units"]].to_numpy()
    return xy, nw["case"].to_numpy().astype(int)


def o1_adjacency(xy, radius: float = 1.25):
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    A = ((d <= radius) & (d > 0)).astype(float)
    return A


def build_O1(rng):
    xy, _ = o1_desks()
    pitch = float(rng.uniform(1.0, 1.6))
    x0, y0 = xy[:, 0].min() - 1.0, xy[:, 1].min() - 1.0
    ix = np.rint(xy[:, 0] - x0).astype(int)
    iy = np.rint(xy[:, 1] - y0).astype(int)
    lat = Lattice(nx=int(ix.max() + 2), ny=int(iy.max() + 2), ax=pitch, ay=pitch, H=2.7)
    return dict(lat=lat, sites=lat.site(ix, iy), xy=xy, pitch=pitch,
                kappa=float(draw_kappa("O1", rng)))


BUILDERS = {"T1": build_T1, "T2": build_T2, "T5": build_T5, "C1": build_C1}
STRATA_N = {"T1": (33, 34), "T2": (19, 26), "T5": (12, 8), "C1": (10, 14)}
TOTALS = {"T1": 23, "T2": 7, "T5": 12, "C1": 12}
CLASS = {"T1": "vehicle", "T2": "vehicle", "T5": "vehicle", "T3": "vehicle",
         "C1": "room", "O1": "room", "R1": "room", "T4": "vehicle"}
