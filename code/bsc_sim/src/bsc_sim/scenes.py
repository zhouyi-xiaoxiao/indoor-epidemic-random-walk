"""Scene data: lattice geometry, zones, D(r), q(r), occupancy.

A :class:`Scene` holds everything the theory and the simulator need:

* ``access``  boolean (ny, nx) mask of walkable cells (Omega = Lambda minus B)
* ``site_xy`` (M, 2) integer coordinates of the M walkable cells
* ``D``       (M,) diffusion coefficient in lattice^2/day
* ``q``       (M,) dimensionless infection efficiency
* ``N``       number of people

Units: length in lattice cells (cell side ``a_m`` metres), time in days.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field

import numpy as np

SCENE_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scenes")
)
SCENE_NAMES = ["office", "supermarket", "classroom", "metro"]


@dataclass
class Scene:
    name: str
    nx: int
    ny: int
    a_m: float
    zone: np.ndarray            # (ny, nx) array of single-character zone codes
    access: np.ndarray          # (ny, nx) bool
    D_grid: np.ndarray          # (ny, nx) float, 0 on barriers
    q_grid: np.ndarray          # (ny, nx) float, 0 on barriers
    N: int
    meta: dict = field(default_factory=dict)

    # ---- derived -----------------------------------------------------
    def __post_init__(self):
        ys, xs = np.nonzero(self.access)
        self.site_xy = np.stack([xs, ys], axis=1)            # (M, 2)
        self.M = len(xs)
        self.site_index = -np.ones((self.ny, self.nx), dtype=np.int64)
        self.site_index[ys, xs] = np.arange(self.M)
        self.D = self.D_grid[ys, xs].astype(float)
        self.q = self.q_grid[ys, xs].astype(float)
        self.site_zone = self.zone[ys, xs]

    @property
    def area_m2(self) -> float:
        return self.nx * self.ny * self.a_m ** 2

    @property
    def occupancy(self) -> float:
        """Mean number of people per walkable cell."""
        return self.N / self.M

    def to_grid(self, v, fill=np.nan):
        """Scatter a per-site vector back onto the (ny, nx) grid."""
        g = np.full((self.ny, self.nx), fill, dtype=float)
        g[self.site_xy[:, 1], self.site_xy[:, 0]] = v
        return g

    def copy_with(self, **kw) -> "Scene":
        d = dict(name=self.name, nx=self.nx, ny=self.ny, a_m=self.a_m,
                 zone=self.zone.copy(), access=self.access.copy(),
                 D_grid=self.D_grid.copy(), q_grid=self.q_grid.copy(),
                 N=self.N, meta=dict(self.meta))
        d.update(kw)
        return Scene(**d)


def keep_largest_component(access: np.ndarray) -> np.ndarray:
    """Return a mask with only the largest 4-connected walkable component."""
    ny, nx = access.shape
    lab = -np.ones((ny, nx), dtype=int)
    sizes = []
    for y0 in range(ny):
        for x0 in range(nx):
            if access[y0, x0] and lab[y0, x0] < 0:
                k = len(sizes)
                stack = [(y0, x0)]
                lab[y0, x0] = k
                n = 0
                while stack:
                    y, x = stack.pop()
                    n += 1
                    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        yy, xx = y + dy, x + dx
                        if 0 <= yy < ny and 0 <= xx < nx and access[yy, xx] \
                                and lab[yy, xx] < 0:
                            lab[yy, xx] = k
                            stack.append((yy, xx))
                sizes.append(n)
    if not sizes:
        return access.copy()
    return lab == int(np.argmax(sizes))


def is_connected(access: np.ndarray) -> bool:
    return bool(keep_largest_component(access).sum() == access.sum())


def load_scene(name: str, N: int | None = None, D0: float = 1.0) -> Scene:
    """Load a Chapter-3 scene from ``scenes/<name>.json``.

    ``D0`` is the reference diffusion coefficient in lattice^2/day (thesis
    default 1.0).  For the metro, D depends on the crowd density
    rho = N / area (persons/m^2) through thesis eq. (3-22).
    """
    with open(os.path.join(SCENE_DIR, name + ".json"), encoding="utf-8") as f:
        sc = json.load(f)
    rows = sc["map"]
    ny, nx = len(rows), len(rows[0])
    zone = np.array([list(r) for r in rows])
    N = int(sc["N_default"] if N is None else N)
    area = sc["size_m"][0] * sc["size_m"][1]
    rho = N / area
    access = np.ones((ny, nx), dtype=bool)
    D = np.zeros((ny, nx))
    q = np.zeros((ny, nx))
    for z, zp in sc["zones"].items():
        m = zone == z
        if zp.get("barrier", False):
            access[m] = False
            continue
        if "D_coef" in zp:           # density-dependent mobility (metro)
            D[m] = D0 * zp["D_coef"] / rho ** zp["D_exp"]
        else:
            D[m] = D0 * zp["D"]
        q[m] = zp["q"]
    if not is_connected(access):
        raise ValueError(f"scene {name}: walkable region is not connected")
    meta = {k: v for k, v in sc.items() if k not in ("map",)}
    meta["D0"] = D0
    meta["rho_per_m2"] = rho
    return Scene(name=name, nx=nx, ny=ny, a_m=sc["a_m"], zone=zone,
                 access=access, D_grid=D, q_grid=q, N=N, meta=meta)


def uniform_room(nx: int, ny: int, N: int, D: float = 1.0, q: float = 1.0,
                 a_m: float = 1.0, barriers: np.ndarray | None = None,
                 name: str = "room") -> Scene:
    """Homogeneous rectangular room, optionally with a barrier mask."""
    access = np.ones((ny, nx), dtype=bool)
    if barriers is not None:
        access &= ~barriers
        access = keep_largest_component(access)
    zone = np.where(access, "U", "B")
    return Scene(name=name, nx=nx, ny=ny, a_m=a_m, zone=zone, access=access,
                 D_grid=np.where(access, D, 0.0), q_grid=np.where(access, q, 0.0),
                 N=N, meta={})


def nested_barriers(nx: int, ny: int, fractions, seed: int = 0):
    """Nested random barrier masks: the barriers at a lower density are a
    subset of those at a higher density (same random permutation of cells).

    After placing the barriers, walkable cells cut off from the largest
    component are also closed, so the realised barrier fraction can exceed the
    nominal one.  To keep the family nested, pockets are closed for the
    densest mask first and every sparser mask keeps the largest component of
    its own walkable set.  Returns ``{fraction: barrier_mask}``.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(nx * ny)
    out = {}
    for f in sorted(fractions):
        k = int(round(f * nx * ny))
        b = np.zeros(nx * ny, dtype=bool)
        b[perm[:k]] = True
        b = b.reshape(ny, nx)
        acc = keep_largest_component(~b)
        out[f] = ~acc
    return out
