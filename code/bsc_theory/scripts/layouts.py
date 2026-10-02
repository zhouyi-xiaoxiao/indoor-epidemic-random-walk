"""layouts.py -- explicit office layouts used in the worked examples.

Grid: nx = 20, ny = 13 (a = 1.5 m), as in the office scene.
Zone codes: 'W' workstations, 'C' corridor, 'M' meeting room, 'K' kitchen, 'B' obstacle/wall.
Zone fractions (of all 260 cells): W 60% (156), C 15% (39), M 10% (26), K 5% (13), B 10% (26).
Zone parameters:
    D/D0: W 0.3, C 1.5, M 0.5, K 1.0       q: W 1.4, C 0.8, M 1.5, K 1.2
"""
from __future__ import annotations

import numpy as np

NX, NY = 20, 13
D_OFFICE = {"W": 0.3, "C": 1.5, "M": 0.5, "K": 1.0}
Q_OFFICE = {"W": 1.4, "C": 0.8, "M": 1.5, "K": 1.2}
D_TWO = {"W": 0.3, "C": 2.0}
Q_TWO = {"W": 1.5, "C": 0.5}


def _blank(fill="W"):
    return np.full((NY, NX), fill, dtype="<U1")


def office_O1():
    """Meeting room in a corner behind walls with one door, kitchen in the opposite corner,
    T-shaped corridor, cabinets scattered among desks.  Exact zone counts."""
    z = _blank("W")
    # meeting room (25 interior cells + door cell) and its walls
    z[8:13, 15:20] = "M"
    z[7, 14:20] = "B"           # wall below the meeting room
    z[8:13, 14] = "B"           # wall on its left ...
    z[10, 14] = "M"             # ... with a door cell
    # kitchen (13 cells), open to the corridor
    z[0:2, 15:20] = "K"
    z[2, 17:20] = "K"
    # corridor (39 cells)
    z[:, 13] = "C"              # vertical spine
    z[6, 0:13] = "C"            # horizontal spine, left
    z[6, 14:20] = "C"           # horizontal spine, right
    z[0:5, 14] = "C"
    z[2, 15:17] = "C"
    # cabinets / pillars among the desks (16 cells)
    for y in (1, 4, 8, 11):
        for x in (1, 4, 7, 10):
            z[y, x] = "B"
    return z


def office_O2():
    """Alternative arrangement: long meeting room along the top wall, kitchen top right,
    corridor = one horizontal spine + two vertical aisles.  Exact zone counts."""
    z = _blank("W")
    z[11:13, 0:13] = "M"        # meeting room 13 x 2
    z[10, 0:13] = "B"           # its wall ...
    z[10, 6] = "C"              # ... with a door (corridor cell)
    z[11:13, 13] = "B"
    z[11:13, 14:20] = "K"       # kitchen 12 cells
    z[10, 14] = "K"             # + 1
    z[9, :] = "C"               # horizontal spine
    z[0:9, 4] = "C"             # aisles
    z[0:9, 9] = "C"
    for y in (2, 5, 7):
        for x in (1, 7, 12, 16):
            z[y, x] = "B"
    return z


def two_zone(kind="side"):
    """Two-zone office (no obstacles): ~70% workspace, ~30% corridor."""
    z = _blank("W")
    if kind == "side":          # corridor band along one wall (4 rows = 80 cells, 30.8%)
        z[0:4, :] = "C"
    elif kind == "aisles":      # desk blocks 3 rows deep separated by aisles (rows 0,4,8,12 = 80 cells)
        z[[0, 4, 8, 12], :] = "C"
    elif kind == "grid":        # vertical aisles every 4 columns + one cross aisle (80 cells)
        z[:, [2, 6, 10, 14, 18]] = "C"
        z[6, :] = "C"
    elif kind == "scatter":     # corridor cells scattered at random (78 cells) - the well-interleaved extreme
        rng = np.random.default_rng(1)
        idx = rng.choice(NX * NY, 78, replace=False)
        z.ravel()[idx] = "C"
    else:
        raise ValueError(kind)
    return z


def counts(z):
    u, c = np.unique(z, return_counts=True)
    return dict(zip(u.tolist(), c.tolist()))


def fields(z, Dmap, qmap):
    """mask (accessible), D grid, q grid from a zone array."""
    mask = z != "B"
    D = np.zeros(z.shape)
    q = np.zeros(z.shape)
    for k, v in Dmap.items():
        D[z == k] = v
    for k, v in qmap.items():
        q[z == k] = v
    return mask, D, q


if __name__ == "__main__":
    for name, z in [("O1", office_O1()), ("O2", office_O2())] + [(k, two_zone(k)) for k in ("side", "aisles", "grid", "scatter")]:
        print(name, counts(z))
        for row in z[::-1]:
            print("   " + "".join(row))
