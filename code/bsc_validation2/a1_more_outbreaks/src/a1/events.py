"""Event geometry for avenue A1 (pre-registration sections 2-4).

`build(ev, opt)` returns
    meta : dict(cls, pathogen, n_bins, sizes (per bin), K, obs (bin counts), near_bins, label, T)
    draws: list of dict(lat, src (site array), rec (site array), bins (int array), kappa, segs)
Seat-map events have one geometry and 7 kappa nodes; zonal events have random
configurations (seed fixed), each paired with one kappa node in turn.
Outcome counts enter only through meta['obs'] / meta['y'], which the kernel
tables (01_tables.py) never read apart from the total K.
"""
from __future__ import annotations
import numpy as np
from . import digitized as dg
from .lattice import Lattice

SEED = 20261002
DEP = 0.93                      # deposition + inactivation, 1/h
AX, AY, HCAB = 0.5, 0.8, 2.0
EVENTS = ["F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "W1", "W2", "P1"]
CLS = {e: "air" for e in EVENTS[:8]} | {"W1": "room", "W2": "room", "P1": "room"}
CAL = {"air": ["F1", "F2", "F3", "F4"], "room": ["W1", "W2"]}
HOLD = {"air": ["F5", "F6", "F7", "F8"], "room": ["P1"]}
PATH = {"F1": "SARS-CoV-1", "F2": "M. tuberculosis", "F3": "influenza A(H1N1)pdm09",
        "F4": "influenza A(H1N1)pdm09", "F5": "SARS-CoV-2", "F6": "SARS-CoV-2", "F7": "SARS-CoV-2",
        "F8": "SARS-CoV-2", "W1": "SARS-CoV-1", "W2": "SARS-CoV-1", "P1": "SARS-CoV-2"}
LABEL = {"F1": "CA112 (Olsen 2003)", "F2": "Chicago-Honolulu (Kenyon 1996)", "F3": "LA-Auckland (Baker 2010)",
         "F4": "Cancun-Birmingham (Young 2014)", "F5": "Naha (Toyokawa 2022)", "F6": "Sydney-Perth (Speake 2020)",
         "F7": "EK448 (Swadi 2021)", "F8": "Tel Aviv-Frankfurt (Hoehl 2020)", "W1": "Ward 8A students (Wong 2004)",
         "W2": "Ward 8A inpatients (Yu 2005)", "P1": "Meat plant (Günther 2020)"}
DUR = {"F1": 3.0, "F2": 8.75, "F3": 13.0, "F4": 9.5, "F5": 2.0, "F6": 5.0, "F7": 18.0, "F8": 4.67}
X33 = [0, 1, 2, 4, 5, 6]                    # 3-3, nx = 7
X242 = [0, 1, 3, 4, 5, 6, 8, 9]             # 2-4-2, nx = 10
X343 = [0, 1, 2, 4, 5, 6, 7, 9, 10, 11]     # 3-4-3, nx = 12


def kappa_nodes(kind: str, n: int = 7, opt=None):
    q = (np.arange(n) + 0.5) / n
    opt = opt or {}
    if kind == "air":
        lo, hi = opt.get("cabin_ach", (10.0, 30.0))
        return np.exp(np.log(lo) + q * (np.log(hi) - np.log(lo))) + DEP
    if kind == "ward":
        from scipy.stats import norm
        return 7.79 * np.exp(0.2 * norm.ppf(q)) + DEP
    if kind == "plant":
        return 0.3 + q * 0.7 + DEP
    raise ValueError(kind)


def row_bins(dist, edges):
    """edges e.g. [1, 2, 5]: bins {<=1}, {2}, {3-5}, {>5}."""
    return np.searchsorted(np.asarray(edges), dist, side="left")


def _seatmap(cells, nx, ny, ay, opt, T, edges=(1, 2, 5)):
    """cells: list of (x, y, symbol) with symbols I (source), o/C/P (susceptible), others ignored."""
    lat = Lattice(nx, ny, AX, ay, HCAB)
    src = np.array([lat.site(x, y) for x, y, s in cells if s == "I"])
    sus = [(x, y, s) for x, y, s in cells if s in "oCP"]
    rec = np.array([lat.site(x, y) for x, y, s in sus])
    sy = np.array([y for x, y, s in cells if s == "I"])
    ry = np.array([y for x, y, s in sus])
    dist = np.abs(ry[:, None] - sy[None, :]).min(axis=1)
    bins = row_bins(dist, edges)
    wide = opt.get("outcome") == "wide"
    y = np.array([s == "C" or (wide and s == "P") for x, y_, s in sus])
    draws = [dict(lat=lat, src=src, rec=rec, bins=bins, kappa=k, segs=[T]) for k in kappa_nodes("air", opt=opt)]
    return draws, bins, y, dist


def _cells_lines(lines, xmap):
    """lines: list of strings (one per seat line, characters = rows)."""
    return [(xmap[i], j, ch) for i, s in enumerate(lines) for j, ch in enumerate(s)]


def _cells_rows(rows, order, xmap):
    """rows: dict row -> string (characters = seats)."""
    return [(xmap[i], order.index(r), ch) for r, s in rows.items() for i, ch in enumerate(s)]


def build(ev: str, opt=None):
    opt = dict(opt or {})
    rng = np.random.default_rng([SEED, EVENTS.index(ev)])
    ay = opt.get("ay", AY)
    meta = dict(ev=ev, cls=CLS[ev], pathogen=PATH[ev], label=LABEL[ev])
    y = None
    if ev == "F1":
        sym = {"X": "o", "/": "o", "G": "C", "I": "I", ".": "e", "-": "-"}
        if opt.get("interviewed"):
            sym["/"] = "x"
        lines = ["".join(sym[c] for c in dg.OLSEN[k]) for k in "ABCDEF"]
        draws, bins, y, dist = _seatmap(_cells_lines(lines, X33), 7, 22, ay, opt, DUR[ev])
    elif ev == "F3":
        lines = list(dg.BAKER)
        if opt.get("src") == "alt1":
            lines = [s.replace("j", "I") for s in lines]
        draws, bins, y, dist = _seatmap(_cells_lines(lines, X343), 12, 15, ay, opt, DUR[ev])
    elif ev == "F4":
        draws, bins, y, dist = _seatmap(_cells_rows(dg.YOUNG, list(range(1, 48)), X242), 10, 47, ay, opt, DUR[ev])
    elif ev == "F5":
        rows = dict(dg.TOYO)
        if opt.get("hh"):
            keep = {"14B", "25F"}
            for hh in dg.TOYO_HOUSEHOLDS:
                for seat in hh:
                    if seat not in keep:
                        r, c = int(seat[:-1]), "ABCFGH".index(seat[-1])
                        rows[r] = rows[r][:c] + "x" + rows[r][c + 1:]
        draws, bins, y, dist = _seatmap(_cells_rows(rows, dg.TOYO_ROWS, X33), 7, 30, ay, opt, DUR[ev])
    elif ev == "F6":
        lines = [list(s) for s in dg.SPEAKE]
        if opt.get("src") == "alt1":            # three symptomatic sources only
            for i, s in enumerate(lines):
                for j, ch in enumerate(s):
                    if ch == "I" and (i, j) not in dg.SPEAKE_SYMPT:
                        lines[i][j] = "j"
        if opt.get("src") == "alt2":            # all six infectious passengers
            for (i, j) in dg.SPEAKE_B1_INFECTIOUS:
                lines[i][j] = "I"
        lines = ["".join(s) for s in lines]
        # primary outcome: all 11 secondary cases; S-out: flight-associated only (8)
        lines = [s.replace("P", "o" if opt.get("outcome") == "wide" else "C") for s in lines]
        draws, bins, y, dist = _seatmap(_cells_lines(lines, X242), 10, 17, ay, opt, DUR[ev])
    elif ev == "F8":
        draws, bins, y, dist = _seatmap(_cells_rows(dg.HOEHL, dg.HOEHL_ROWS, X33), 7, 31, ay, opt, DUR[ev])
    elif ev == "F2":
        lat = Lattice(12, 14, AX, ay, HCAB)
        kn = kappa_nodes("air", opt=opt)
        draws = []
        for i in range(opt.get("n_cfg", 200)):
            r0 = rng.integers(2, 12)
            x0 = rng.choice(X343)
            near_seats = [(x, r) for r in range(r0 - 2, r0 + 3) for x in X343 if not (x == x0 and r == r0)]
            far_seats = [(x, r) for r in range(14) if abs(r - r0) > 2 for x in X343]
            a = [near_seats[k] for k in rng.choice(len(near_seats), 13, replace=False)]
            b = [far_seats[k] for k in rng.choice(len(far_seats), 55, replace=False)]
            rec = np.array([lat.site(x, r) for x, r in a + b])
            draws.append(dict(lat=lat, src=np.array([lat.site(x0, r0)]), rec=rec,
                              bins=np.array([0] * 13 + [1] * 55), kappa=kn[i % len(kn)], segs=[DUR[ev]]))
        bins = draws[0]["bins"]
        meta["obs"] = [4, 2]
    elif ev == "F7":
        lat = Lattice(12, 34, AX, ay, HCAB)
        L = dict(zip("ABCDEFGHJK", X343))
        row = lambda r: r - 17
        srcs = ["26G", "26D"] if opt.get("src") != "alt1" else ["26G"]
        near = ["26A", "26C", "27D", "27K"] + ["24C", "24D", "24E", "24F", "24G", "28A", "28D", "28G", "28K"]
        b_near = [0] * 4 + [1] * 9
        if opt.get("src") == "alt1":             # passenger B becomes a susceptible non-case at 26D
            near = ["26D"] + near
            b_near = [0] + b_near
        far_seats = [(x, row(r)) for r in list(range(17, 23)) + list(range(31, 51)) for x in X343]
        kn = kappa_nodes("air", opt=opt)
        draws = []
        for i in range(opt.get("n_cfg", 200)):
            b = [far_seats[k] for k in rng.choice(len(far_seats), 71, replace=False)]
            rec = np.array([lat.site(L[s[-1]], row(int(s[:-1]))) for s in near] + [lat.site(x, r) for x, r in b])
            draws.append(dict(lat=lat, src=np.array([lat.site(L[s[-1]], row(int(s[:-1]))) for s in srcs]),
                              rec=rec, bins=np.array(b_near + [2] * 71), kappa=kn[i % len(kn)], segs=[DUR[ev]]))
        bins = draws[0]["bins"]
        obs = [1, 3, 0]
        if opt.get("outcome") == "wide":
            obs = [1, 4, 0]
        meta["obs"] = obs
    elif ev in ("W1", "W2"):
        lat = Lattice(21, 16, 1.0, 1.0, 3.0)
        bed = {}
        for i, b in enumerate(["13", "14", "15", "16", "16x"]):
            bed[b] = (2 * i, 1)
        for i, b in enumerate(["12", "11", "10", "9", "9x"]):
            bed[b] = (2 * i, 5)
        for i, b in enumerate(["17x", "17", "18", "19", "20"]):
            bed[b] = (12 + 2 * i, 1)
        for i, b in enumerate(["24x", "24", "23", "22", "21"]):
            bed[b] = (12 + 2 * i, 5)
        for i, b in enumerate(["5", "6", "7", "8"]):
            bed[b] = (2 * i, 11)
        for i, b in enumerate(["4", "3", "2", "1", "1x"]):
            bed[b] = (2 * i, 15)
        for i, b in enumerate(["34", "33"]):
            bed[b] = (2 * i, 8)
        for i, b in enumerate(["25x", "25", "26", "27", "28"]):
            bed[b] = (12 + 2 * i, 11)
        for i, b in enumerate(["32x", "32", "31", "30", "29"]):
            bed[b] = (12 + 2 * i, 15)
        same = ["9", "9x", "10", "12", "13", "14", "15", "16", "16x"]
        adj = ["17", "18", "19", "20", "21", "22", "23", "24", "17x", "24x"]
        dist_beds = [b for b in bed if b not in same + adj + ["11"]]
        src = np.array([lat.site(*bed["11"])])
        kn = kappa_nodes("ward")
        st = lambda b: lat.site(*bed[b])
        if ev == "W1":
            cub_other = ["9", "9x", "13", "14", "15", "16", "16x"]
            others = adj + dist_beds
            draws = []
            for i in range(opt.get("n_cfg", 200)):
                rec = [st("12")] * 3 + [st(cub_other[k]) for k in rng.integers(0, len(cub_other), 8)] \
                    + [st(others[k]) for k in rng.integers(0, len(others), 8)]
                draws.append(dict(lat=lat, src=src, rec=np.array(rec), bins=np.array([0] * 3 + [1] * 8 + [2] * 8),
                                  kappa=kn[i % len(kn)], segs=[40.0 / 60.0]))
            meta["obs"] = [3, 4, 0]
            meta["near_bins"] = [0, 1]
        else:
            rec = [st(same[i % len(same)]) for i in range(20)] + [st(adj[i % len(adj)]) for i in range(21)] \
                + [st(dist_beds[i % len(dist_beds)]) for i in range(33)]
            draws = [dict(lat=lat, src=src, rec=np.array(rec), bins=np.array([0] * 20 + [1] * 21 + [2] * 33),
                          kappa=k, segs=[100.0]) for k in kn]
            meta["obs"] = [13, 11, 6]
            meta["near_bins"] = [0]
        bins = draws[0]["bins"]
    elif ev == "P1":
        lat = Lattice(12, 32, 1.0, 1.0, 6.1)
        sizes = [9, 17, 22, 30]
        edges = [0.0, 4.0, 8.0, 12.0, 29.0]
        xs, ys = np.meshgrid(np.arange(12), np.arange(32), indexing="xy")
        xs, ys = xs.ravel(), ys.ravel()
        kn = kappa_nodes("plant")
        draws = []
        for i in range(opt.get("n_cfg", 300)):
            x0, y0 = int(rng.integers(0, 2)), int(rng.integers(0, 16))
            d = np.hypot(xs - x0, ys - y0)
            rec = []
            for b in range(4):
                cand = np.flatnonzero((d > edges[b]) & (d <= edges[b + 1]))
                rec += list(lat.site(xs[cand], ys[cand])[rng.integers(0, len(cand), sizes[b])])
            draws.append(dict(lat=lat, src=np.array([lat.site(x0, y0)]), rec=np.array(rec),
                              bins=np.repeat(np.arange(4), sizes), kappa=kn[i % len(kn)], segs=[8.0, 8.0, 8.0]))
        bins = draws[0]["bins"]
        meta["obs"] = [5, 12, 1, 2]
        meta["near_bins"] = [0, 1]
    else:
        raise ValueError(ev)

    nb = int(bins.max()) + 1
    meta["n_bins"] = nb
    meta["sizes"] = [int((bins == b).sum()) for b in range(nb)]
    if y is not None:
        meta["y"] = y
        meta["dist"] = dist
        meta["obs"] = [int(y[bins == b].sum()) for b in range(nb)]
    if "near_bins" not in meta:                      # aircraft: rows <= 2
        meta["near_bins"] = [0] if ev == "F2" else [0, 1]
    meta["K"] = int(sum(meta["obs"]))
    meta["T"] = float(sum(draws[0]["segs"]))
    return meta, draws
