"""Python front end of the exact stochastic simulator (C core, ctypes).

    res = simulate(scene, beta=0.5, gamma=0.14, n_rep=1000, seed=1)

The C source ``gillespie.c`` is compiled on first use with the system C
compiler into ``build/libgillespie.dylib``.  A slow pure-Python reference
implementation of the same algorithm is in ``refsim.py``; the two are
cross-validated in ``experiments/01_validate.py``.
"""
from __future__ import annotations

import ctypes as C
import os
import subprocess
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .scenes import Scene
from .theory import contact_kernel, hop_rates

_HERE = os.path.dirname(os.path.abspath(__file__))
_BUILD = os.path.normpath(os.path.join(_HERE, "..", "..", "build"))
_LIB = None


def _lib():
    global _LIB
    if _LIB is not None:
        return _LIB
    src = os.path.join(_HERE, "gillespie.c")
    out = os.path.join(_BUILD, "libgillespie.dylib")
    os.makedirs(_BUILD, exist_ok=True)
    if (not os.path.exists(out)) or os.path.getmtime(out) < os.path.getmtime(src):
        subprocess.run(["cc", "-O2", "-shared", "-fPIC", "-o", out, src, "-lm"],
                       check=True)
    lib = C.CDLL(out)
    ip = np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS")
    dp = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
    fp = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS")
    lib.run_batch.restype = C.c_int
    lib.run_batch.argtypes = [
        C.c_int, ip, ip,                       # M, nbr_ptr, nbr_idx
        C.c_int, dp,                           # n_prof, W
        ip, ip, dp,                            # kernel CSR
        dp, C.c_int, C.c_double,               # bq, mode, gamma
        C.c_int, dp, ip, dp,                   # schedule
        C.c_int, C.c_int, ip, C.c_int, ip,     # N, n_rep, init_pos, n_index, index_agent
        C.c_int, C.c_double, C.c_double, C.c_int, C.c_uint64,
        C.c_int, ip,                           # snapshots
        ip, ip,                                # out_S, out_I
        fp, fp, ip, ip, ip,                    # t_inf, t_rec, site_inf, gen, infector
        dp, dp,                                # index_exposure, t_ext
        np.ctypeslib.ndpointer(np.uint16, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(np.int64, flags="C_CONTIGUOUS"),
    ]
    _LIB = lib
    return lib


@dataclass
class SimResult:
    t: np.ndarray            # (n_out,)
    S: np.ndarray            # (n_rep, n_out)
    I: np.ndarray            # (n_rep, n_out)
    t_inf: np.ndarray        # (n_rep, N), -1 if never infected
    t_rec: np.ndarray
    site_inf: np.ndarray
    gen: np.ndarray
    infector: np.ndarray
    index_agent: np.ndarray  # (n_rep, n_index)
    index_exposure: np.ndarray
    t_ext: np.ndarray
    snap: np.ndarray         # (n_rep, n_snap, M) infectives per cell
    snap_times: np.ndarray
    n_events: np.ndarray
    init_pos: np.ndarray
    N: int

    # ---- per-replicate summaries ------------------------------------
    @property
    def final_size(self):
        """Number ever infected, including the index case(s)."""
        return (self.t_inf >= 0).sum(axis=1)

    @property
    def attack_rate(self):
        return self.final_size / self.N

    @property
    def peak_I(self):
        return self.I.max(axis=1)

    @property
    def peak_time(self):
        return self.t[self.I.argmax(axis=1)]

    @property
    def index_offspring(self):
        """Secondary cases caused directly by the first index case."""
        idx0 = self.index_agent[:, 0]
        return (self.infector == idx0[:, None]).sum(axis=1)

    def generation_sizes(self, gmax=6):
        out = np.zeros((self.gen.shape[0], gmax + 1), dtype=np.int64)
        for g in range(gmax + 1):
            out[:, g] = (self.gen == g).sum(axis=1)
        return out

    def major(self, threshold_frac=0.1):
        """Boolean mask of major outbreaks (attack rate >= threshold)."""
        return self.attack_rate >= threshold_frac


def _kernel_arrays(scene: Scene, P, beta, rho, q):
    P = sp.csr_matrix(P)
    P.sort_indices()
    rowq = np.repeat(beta * q / rho, np.diff(P.indptr))
    return (P.indptr.astype(np.int32), P.indices.astype(np.int32),
            np.ascontiguousarray(P.data * rowq, dtype=np.float64))


def draw_initial_positions(M: int, N: int, n_rep: int, rng, sites=None):
    """iid uniform positions on the walkable cells (the stationary law of the
    symmetric walk), optionally restricted to a subset ``sites``."""
    if sites is None:
        return rng.integers(0, M, size=(n_rep, N), dtype=np.int32)
    sites = np.asarray(sites, dtype=np.int32)
    return sites[rng.integers(0, len(sites), size=(n_rep, N))].astype(np.int32)


def simulate(scene: Scene, beta: float, gamma: float, n_rep: int = 100,
             seed: int = 0, t_max: float = 150.0, dt_out: float = 0.25,
             kernel_radius: float = 0.0, P=None, rho_ref: float | None = None,
             mode: str = "kernel", N: int | None = None, q=None, D=None,
             init_pos=None, index_site=None, n_index: int = 1, gmax: int = -1,
             schedule=None, snap_times=()) -> SimResult:
    """Run ``n_rep`` independent epidemics.

    Parameters
    ----------
    kernel_radius, P : contact kernel (radius in lattice units, or an explicit
        row-stochastic matrix).  Default: same-cell rule.
    rho_ref : None for the frequency-dependent calibration
        rho_bar = (N-1)/|Omega|; a number fixes the pair hazard
        (density-dependent transmission, R0 proportional to occupancy).
    mode : 'kernel' (pair hazard) or 'percell' (per-cell rule: same cell,
        beta q n_I / n_total).
    index_site : None -> the index case is person 0 at its random position
        (uniform over walkable cells); an int -> person 0 is moved to that
        cell; an array of M weights -> person 0's cell is drawn from it.
    gmax : highest infectious generation (-1: no cap).  gmax=0 lets only the
        index case transmit, which measures its secondary cases cleanly.
    schedule : None, or dict(period=.., segments=[(t_end, D_array_or_None, m)])
        piecewise-constant periodic schedule of mobility profile and
        transmission multiplier.
    snap_times : times at which the per-cell number of infectives is stored.
    """
    lib = _lib()
    M = scene.M
    N = scene.N if N is None else int(N)
    q = scene.q if q is None else np.asarray(q, float)
    rng = np.random.default_rng(seed)

    # hop-rate profiles
    ptr, idx, W0 = hop_rates(scene, D)
    profiles = [W0]
    if schedule is None:
        seg_end = np.array([1.0])
        seg_prof = np.array([0], dtype=np.int32)
        seg_m = np.array([1.0])
    else:
        seg_end, seg_prof, seg_m = [], [], []
        for t_end, Dseg, m in schedule["segments"]:
            if Dseg is None:
                k = 0
            else:
                profiles.append(hop_rates(scene, Dseg)[2])
                k = len(profiles) - 1
            seg_end.append(t_end)
            seg_prof.append(k)
            seg_m.append(m)
        assert abs(seg_end[-1] - schedule["period"]) < 1e-12
        seg_end = np.array(seg_end, float)
        seg_prof = np.array(seg_prof, dtype=np.int32)
        seg_m = np.array(seg_m, float)
    W = np.ascontiguousarray(np.concatenate(profiles))

    # contact kernel and pair hazards
    if P is None:
        P = contact_kernel(scene, kernel_radius)
    rho = (N - 1) / M if rho_ref is None else float(rho_ref)
    kptr, kidx, kw = _kernel_arrays(scene, P, beta, rho, q)
    bq = np.ascontiguousarray(beta * q, dtype=np.float64)

    # initial conditions
    if init_pos is None:
        init_pos = draw_initial_positions(M, N, n_rep, rng)
    init_pos = np.ascontiguousarray(init_pos, dtype=np.int32).copy()
    assert init_pos.shape == (n_rep, N)
    index_agent = np.tile(np.arange(n_index, dtype=np.int32), (n_rep, 1))
    if index_site is not None:
        if np.ndim(index_site) == 0:
            init_pos[:, 0] = int(index_site)
        else:
            w = np.asarray(index_site, float)
            init_pos[:, 0] = rng.choice(M, size=n_rep, p=w / w.sum())
    index_agent = np.ascontiguousarray(index_agent)

    n_out = int(round(t_max / dt_out)) + 1
    t = np.arange(n_out) * dt_out
    snap_idx = np.ascontiguousarray(
        np.round(np.asarray(snap_times, float) / dt_out).astype(np.int32))
    n_snap = len(snap_idx)
    out_S = np.zeros((n_rep, n_out), dtype=np.int32)
    out_I = np.zeros((n_rep, n_out), dtype=np.int32)
    t_inf = np.zeros((n_rep, N), dtype=np.float32)
    t_rec = np.zeros((n_rep, N), dtype=np.float32)
    site_inf = np.zeros((n_rep, N), dtype=np.int32)
    gen = np.zeros((n_rep, N), dtype=np.int32)
    infector = np.zeros((n_rep, N), dtype=np.int32)
    expo = np.zeros(n_rep)
    t_ext = np.zeros(n_rep)
    snap = np.zeros((n_rep, max(n_snap, 1), M), dtype=np.uint16)
    nev = np.zeros(n_rep, dtype=np.int64)

    seed64 = int(np.random.SeedSequence(seed).generate_state(1, np.uint64)[0])
    rc = lib.run_batch(
        M, ptr, idx, len(profiles), W, kptr, kidx, kw, bq,
        0 if mode == "kernel" else 1, float(gamma),
        len(seg_end), seg_end, seg_prof, seg_m,
        N, n_rep, init_pos, n_index, index_agent,
        int(gmax), float(t_max), float(dt_out), n_out, C.c_uint64(seed64),
        n_snap, snap_idx if n_snap else np.zeros(1, dtype=np.int32),
        out_S, out_I, t_inf, t_rec, site_inf, gen, infector, expo, t_ext,
        snap, nev)
    assert rc == 0
    return SimResult(t=t, S=out_S, I=out_I, t_inf=t_inf, t_rec=t_rec,
                     site_inf=site_inf, gen=gen, infector=infector,
                     index_agent=index_agent, index_exposure=expo,
                     t_ext=t_ext, snap=snap[:, :n_snap], snap_times=t[snap_idx]
                     if n_snap else np.array([]), n_events=nev,
                     init_pos=init_pos, N=N)


# ----------------------------------------------------------------------------
# Statistics helpers
# ----------------------------------------------------------------------------
def wilson(k: int, n: int, z: float = 1.96):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, c - h, c + h


def mean_ci(x, z: float = 1.96):
    """Mean and normal-approximation 95 % interval."""
    x = np.asarray(x, float)
    n = len(x)
    if n == 0:
        return (np.nan, np.nan, np.nan)
    m = x.mean()
    se = x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    return m, m - z * se, m + z * se


def summarize(res: SimResult, threshold_frac: float = 0.1) -> dict:
    """Outbreak statistics of a batch of replicates."""
    n = len(res.final_size)
    maj = res.major(threshold_frac)
    k = int(maj.sum())
    p, lo, hi = wilson(k, n)
    out = dict(n_rep=n, N=res.N, n_major=k, p_major=p, p_major_lo=lo,
               p_major_hi=hi, threshold_frac=threshold_frac)
    ar = res.attack_rate
    out["attack_all_mean"] = float(ar.mean())
    out["attack_all_sd"] = float(ar.std(ddof=1))
    off = res.index_offspring
    m, lo, hi = mean_ci(off)
    out.update(index_offspring_mean=m, index_offspring_lo=lo,
               index_offspring_hi=hi, index_offspring_var=float(off.var(ddof=1)))
    m, lo, hi = mean_ci(res.index_exposure)
    out.update(index_exposure_mean=m, index_exposure_lo=lo, index_exposure_hi=hi)
    if k > 0:
        for name, arr in (("attack", ar[maj]),
                          ("peak_prev", res.peak_I[maj] / res.N),
                          ("peak_time", res.peak_time[maj]),
                          ("duration", res.t_ext[maj])):
            arr = arr[np.isfinite(arr)]
            m, lo, hi = mean_ci(arr)
            out[name + "_major_mean"] = m
            out[name + "_major_lo"] = lo
            out[name + "_major_hi"] = hi
            out[name + "_major_sd"] = float(arr.std(ddof=1)) if len(arr) > 1 else np.nan
            out[name + "_major_q05"] = float(np.quantile(arr, 0.05)) if len(arr) else np.nan
            out[name + "_major_q95"] = float(np.quantile(arr, 0.95)) if len(arr) else np.nan
    out["n_unfinished"] = int(np.isnan(res.t_ext).sum())
    out["events_per_rep"] = float(res.n_events.mean())
    return out
