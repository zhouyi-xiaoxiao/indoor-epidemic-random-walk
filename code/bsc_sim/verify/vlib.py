"""Re-implementation with separately written code for the check of the simulations.
Nothing here imports the bsc_sim package of the first analysis.
"""
import ctypes as C
import json
import os
import subprocess

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
SCDIR = os.path.join(HERE, "..", "scenes")
BETA, GAMMA = 0.5, 0.14


class Room:
    """walkable cells indexed 0..M-1 ; D, q per cell ; N people."""

    def __init__(self, access, Dg, qg, N, a_m=1.0, zone=None):
        self.access = np.asarray(access, bool)
        self.ny, self.nx = self.access.shape
        self.idx = -np.ones((self.ny, self.nx), int)
        cells = [(r, c) for r in range(self.ny) for c in range(self.nx) if self.access[r, c]]
        for k, (r, c) in enumerate(cells):
            self.idx[r, c] = k
        self.cells = cells
        self.M = len(cells)
        self.D = np.array([Dg[r][c] for r, c in cells], float)
        self.q = np.array([qg[r][c] for r, c in cells], float)
        self.N = N
        self.a_m = a_m
        self.zone = None if zone is None else np.array([zone[r][c] for r, c in cells])


def load(name, D0=1.0, N=None):
    d = json.load(open(os.path.join(SCDIR, name + ".json"), encoding="utf-8"))
    rows = d["map"]
    N = d["N_default"] if N is None else N
    rho = N / (d["size_m"][0] * d["size_m"][1])
    ny, nx = len(rows), len(rows[0])
    acc = np.zeros((ny, nx), bool)
    Dg = np.zeros((ny, nx))
    qg = np.zeros((ny, nx))
    for r in range(ny):
        for c in range(nx):
            z = d["zones"][rows[r][c]]
            if z.get("barrier"):
                continue
            acc[r, c] = True
            if "D_coef" in z:
                Dg[r, c] = D0 * z["D_coef"] / rho ** z["D_exp"]
            else:
                Dg[r, c] = D0 * z["D"]
            qg[r, c] = z["q"]
    return Room(acc, Dg, qg, N, d["a_m"], zone=rows)


def uniform(nx, ny, N, D=1.0, q=1.0, barriers=None):
    acc = np.ones((ny, nx), bool)
    if barriers is not None:
        acc &= ~barriers
    return Room(acc, np.full((ny, nx), D), np.full((ny, nx), q), N)


def neighbours(room):
    """nbr[M,4] (-1 if none), w[M,4] harmonic-mean hop rates."""
    nbr = -np.ones((room.M, 4), np.int32)
    w = np.zeros((room.M, 4))
    for k, (r, c) in enumerate(room.cells):
        for j, (dr, dc) in enumerate(((0, 1), (0, -1), (1, 0), (-1, 0))):
            rr, cc = r + dr, c + dc
            if 0 <= rr < room.ny and 0 <= cc < room.nx and room.idx[rr, cc] >= 0:
                t = room.idx[rr, cc]
                nbr[k, j] = t
                a, b = room.D[k], room.D[t]
                w[k, j] = 2 * a * b / (a + b) if a + b > 0 else 0.0
    return nbr, w


def gen_matrix(room):
    """dense generator, Lm[y,x] = rate x->y, zero column sums."""
    nbr, w = neighbours(room)
    Lm = np.zeros((room.M, room.M))
    for x in range(room.M):
        for j in range(4):
            y = nbr[x, j]
            if y >= 0:
                Lm[y, x] += w[x, j]
                Lm[x, x] -= w[x, j]
    return Lm


def kernel(room, radius_cells):
    """row-stochastic P[x,y], uniform over walkable cells with centre distance <= radius."""
    M = room.M
    P = np.zeros((M, M))
    R = int(np.floor(radius_cells + 1e-9))
    for k, (r, c) in enumerate(room.cells):
        for dr in range(-R, R + 1):
            for dc in range(-R, R + 1):
                if dr * dr + dc * dc <= radius_cells ** 2 + 1e-9:
                    rr, cc = r + dr, c + dc
                    if 0 <= rr < room.ny and 0 <= cc < room.nx and room.idx[rr, cc] >= 0:
                        P[k, room.idx[rr, cc]] = 1.0
    return P / P.sum(axis=1, keepdims=True)


def ngm(room, P=None, beta=BETA, gamma=GAMMA, q=None, Lm=None):
    Lm = gen_matrix(room) if Lm is None else Lm
    q = room.q if q is None else q
    G = np.linalg.inv(gamma * np.eye(room.M) - Lm)
    K = (beta * q)[:, None] * G
    if P is not None:
        K = P.T @ K
    return K


def rho(K):
    return float(np.max(np.abs(np.linalg.eigvals(K))))


def pair_solve(room, P=None, beta=BETA, gamma=GAMMA, N=None, rho_ref=None, q=None, Lm=None):
    """p[x,y]: prob. an infective born at x infects a person then at y. Direct sparse LU
    for small M, GMRES/CG otherwise.  Separate from the code of the first analysis."""
    M = room.M
    N = room.N if N is None else N
    q = room.q if q is None else q
    Lm = gen_matrix(room) if Lm is None else Lm
    P = np.eye(M) if P is None else P
    rb = (N - 1) / M if rho_ref is None else rho_ref
    h = (beta * q / rb)[:, None] * P          # h[x,y]
    Ls = sp.csr_matrix(Lm)
    I = sp.identity(M, format="csr")
    # backward generator on pair space = L^T (+) L^T ; L symmetric here
    L2 = sp.kron(Ls.T, I) + sp.kron(I, Ls.T)
    A = (gamma * sp.identity(M * M) + sp.diags(h.ravel()) - L2).tocsc()
    b = h.ravel().copy()                       # p solves (gamma + h - L2) p = h
    if M * M <= 12000:
        p = spla.spsolve(A, b)
    else:
        d = A.diagonal()
        p, info = spla.cg(A.tocsr(), b, rtol=1e-10, maxiter=100000,
                          M=spla.LinearOperator(A.shape, matvec=lambda v: v / d))
        assert info == 0
    return p.reshape(M, M)


# ------------------------------------------------------------- C simulator
_lib = None


def lib():
    global _lib
    if _lib is None:
        src = os.path.join(HERE, "vsim.c")
        out = os.path.join(HERE, "libvsim.dylib")
        if not os.path.exists(out) or os.path.getmtime(out) < os.path.getmtime(src):
            subprocess.run(["cc", "-O2", "-shared", "-fPIC", "-o", out, src, "-lm"], check=True)
        _lib = C.CDLL(out)
        ip = np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS")
        dp = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
        _lib.vrun.restype = C.c_int
        _lib.vrun.argtypes = [C.c_int, ip, dp, ip, ip, dp, dp, C.c_int, C.c_double,
                              C.c_int, C.c_int, ip, ip, C.c_int, C.c_double, C.c_uint64,
                              C.c_int, dp, dp, dp,
                              ip, ip, ip, dp, dp, dp, ip, ip, dp, C.c_double, C.c_int,
                              C.c_int, ip]
    return _lib


def sim(room, n_rep, seed, P=None, mode=0, gmax=-1, beta=BETA, gamma=GAMMA, N=None,
        rho_ref=None, q=None, index_pos=None, init_pos=None, t_max=2000.0,
        sched=None, n_index=1, dt_curve=0.5, n_curve=0, move_all=False):
    """mode 0: pair hazard beta q(x) P(y|x)/rho_bar ; mode 1: per-cell rule.
    sched: list of (t_end, hop_scale, m) over one period (period = last t_end).
    index_pos: None (uniform) or int array (n_rep,) of cells for person 0.
    init_pos: None or int array (n_rep, N) of all starting cells."""
    M = room.M
    N = room.N if N is None else N
    q = room.q if q is None else q
    nbr, w = neighbours(room)
    P = np.eye(M) if P is None else P
    rb = (N - 1) / M if rho_ref is None else rho_ref
    Ps = sp.csr_matrix(P)
    Ps.sort_indices()
    kw = Ps.data * np.repeat(beta * q / rb, np.diff(Ps.indptr))
    if sched is None:
        sched = [(1.0, 1.0, 1.0)]
    se = np.array([s[0] for s in sched], float)
    sh = np.array([s[1] for s in sched], float)
    sm = np.array([s[2] for s in sched], float)
    if init_pos is None:
        ipos = -np.ones((n_rep, N), np.int32)
        if index_pos is not None:
            ipos[:, 0] = index_pos
    else:
        ipos = np.ascontiguousarray(init_pos, np.int32)
    fs = np.zeros(n_rep, np.int32)
    off = np.zeros(n_rep, np.int32)
    pk = np.zeros(n_rep, np.int32)
    pkt = np.zeros(n_rep)
    dur = np.zeros(n_rep)
    expo = np.zeros(n_rep)
    ncv = max(n_curve, 1)
    g2 = np.zeros(n_rep * 8, np.int32)
    cg = np.zeros(4 * M, np.int32)
    cv = np.zeros(n_rep * ncv)
    nunf = np.zeros(1, np.int32)
    move_all = 1 if (mode == 1 or move_all) else 0
    lib().vrun(M, np.ascontiguousarray(nbr.ravel()), np.ascontiguousarray(w.ravel()),
               Ps.indptr.astype(np.int32), Ps.indices.astype(np.int32),
               np.ascontiguousarray(kw), np.ascontiguousarray(beta * q), mode, gamma,
               N, n_rep, np.ascontiguousarray(ipos.ravel()), np.array([n_index], np.int32), gmax,
               t_max, C.c_uint64(seed), len(sched), se, sh, sm,
               fs, off, pk, pkt, dur, expo, g2, cg, cv, dt_curve, n_curve, move_all, nunf)
    return dict(final=fs, off=off, peak=pk, peak_t=pkt, dur=dur, expo=expo,
                gens=g2.reshape(n_rep, 8), cellgen=cg.reshape(4, M),
                curve=cv.reshape(n_rep, ncv), unfinished=int(nunf[0]), N=N)


def mci(x):
    x = np.asarray(x, float)
    m = x.mean()
    se = x.std(ddof=1) / np.sqrt(len(x))
    return m, se


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, c - h, c + h
