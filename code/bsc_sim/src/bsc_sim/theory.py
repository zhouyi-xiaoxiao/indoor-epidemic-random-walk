"""Deterministic theory for the lattice random-walk SIR model.

Model:

* People are independent continuous-time random walkers on the walkable cells
  Omega of a rectangular lattice.  The hop rate between neighbouring cells is
  the harmonic mean

      W(x -> y) = 2 D(x) D(y) / (D(x) + D(y))      [1/day, lattice units]

  so the generator L (L[y, x] = W(x -> y), columns sum to zero) is symmetric
  and the stationary occupation is uniform on Omega.

* An infective at cell x makes infectious contacts at rate beta*q(x).  A
  contact lands on cell y with probability P(y | x) (the contact kernel; the
  same-cell rule is P = identity).  At the
  individual level the pair hazard is

      h(x, y) = beta * q(x) * P(y | x) / rho_bar,   rho_bar = (N - 1) / |Omega|

  so that the expected contact rate of an infective in a fully susceptible
  room is exactly beta*q(x) (frequency-dependent calibration).

* Recovery at rate gamma.

Mean-field (many people per contact neighbourhood) this gives

      dI/dt = L I + F(S) I - gamma I,     F[y, x] = beta q(x) P(y|x) S_y / rho_bar

and the next-generation matrix

      K = P^T diag(beta q) (gamma I - L)^{-1},     R0 = spectral radius of K.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import solve_ivp

from .scenes import Scene

NBRS = ((1, 0), (-1, 0), (0, 1), (0, -1))


# ----------------------------------------------------------------------------
# Movement
# ----------------------------------------------------------------------------
def hop_rates(scene: Scene, D: np.ndarray | None = None):
    """Neighbour lists and harmonic-mean hop rates.

    Returns (ptr, idx, W): CSR arrays over walkable sites; for site s the
    neighbours are idx[ptr[s]:ptr[s+1]] with rates W[ptr[s]:ptr[s+1]].
    """
    D = scene.D if D is None else np.asarray(D, float)
    ptr = [0]
    idx, W = [], []
    for s in range(scene.M):
        x, y = scene.site_xy[s]
        for dx, dy in NBRS:
            xx, yy = x + dx, y + dy
            if 0 <= xx < scene.nx and 0 <= yy < scene.ny:
                t = scene.site_index[yy, xx]
                if t >= 0:
                    d1, d2 = D[s], D[t]
                    w = 2.0 * d1 * d2 / (d1 + d2) if (d1 + d2) > 0 else 0.0
                    idx.append(t)
                    W.append(w)
        ptr.append(len(idx))
    return (np.array(ptr, dtype=np.int32), np.array(idx, dtype=np.int32),
            np.array(W, dtype=float))


def generator(scene: Scene, D: np.ndarray | None = None) -> sp.csr_matrix:
    """Random-walk generator L (symmetric, zero column sums)."""
    ptr, idx, W = hop_rates(scene, D)
    M = scene.M
    rows = np.repeat(np.arange(M), np.diff(ptr))
    A = sp.csr_matrix((W, (idx, rows)), shape=(M, M))      # A[y, x] = W(x->y)
    out = np.asarray(A.sum(axis=0)).ravel()
    return (A - sp.diags(out)).tocsr()


# ----------------------------------------------------------------------------
# Contact kernel
# ----------------------------------------------------------------------------
def contact_kernel(scene: Scene, radius: float = 0.0) -> sp.csr_matrix:
    """Row-stochastic contact kernel P[x, y] = P(y | x).

    Uniform over the walkable cells whose centres lie within Euclidean
    distance ``radius`` (in lattice units) of x.  radius = 0 gives the
    same-cell rule (identity).  Line of sight is not checked.
    """
    M = scene.M
    if radius <= 0:
        return sp.identity(M, format="csr")
    r = int(np.floor(radius + 1e-9))
    offs = [(dx, dy) for dx in range(-r, r + 1) for dy in range(-r, r + 1)
            if dx * dx + dy * dy <= radius * radius + 1e-9]
    rows, cols = [], []
    for s in range(M):
        x, y = scene.site_xy[s]
        for dx, dy in offs:
            xx, yy = x + dx, y + dy
            if 0 <= xx < scene.nx and 0 <= yy < scene.ny:
                t = scene.site_index[yy, xx]
                if t >= 0:
                    rows.append(s)
                    cols.append(t)
    P = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(M, M))
    P = sp.diags(1.0 / np.asarray(P.sum(axis=1)).ravel()) @ P
    return P.tocsr()


def kernel_radius_cells(scene: Scene, ell_c_m: float) -> float:
    """Contact radius in lattice units for a physical contact range (metres)."""
    return ell_c_m / scene.a_m


# ----------------------------------------------------------------------------
# Next-generation matrix
# ----------------------------------------------------------------------------
def ngm(scene: Scene, beta: float, gamma: float, P=None, L=None, q=None):
    """Dense next-generation matrix K = P^T diag(beta q) (gamma I - L)^{-1}."""
    L = generator(scene) if L is None else L
    q = scene.q if q is None else q
    M = scene.M
    G = np.linalg.inv(gamma * np.eye(M) - L.toarray())
    K = (beta * q)[:, None] * G
    if P is not None:
        K = P.T @ K
    return np.asarray(K)


def R0_dense(scene: Scene, beta: float, gamma: float, P=None, L=None, q=None,
             return_vectors: bool = False):
    """R0 = rho(K) by a dense eigen-solve; optionally the Perron vectors.

    Returns R0, or (R0, v, u) where v is the right Perron vector (distribution
    of the birth cell of a typical case, normalised to sum 1) and u is the left
    Perron vector (reproductive value of a case born at each cell, normalised
    so that u . v = 1).
    """
    K = ngm(scene, beta, gamma, P=P, L=L, q=q)
    if not return_vectors:
        ev = np.linalg.eigvals(K)
        return float(np.max(ev.real))
    w, vl, vr = sla.eig(K, left=True, right=True)
    k = int(np.argmax(w.real))
    v = np.abs(vr[:, k].real)
    v /= v.sum()
    u = np.abs(vl[:, k].real)
    u /= (u @ v)
    return float(w[k].real), v, u


def R0_sparse(scene: Scene, beta: float, gamma: float, P=None, L=None, q=None,
              tol: float = 1e-10, maxit: int = 10000, return_vector=False):
    """R0 by power iteration with one sparse LU of (gamma I - L).

    This is the 'spectral evaluation' that is benchmarked against Monte Carlo.
    """
    L = generator(scene) if L is None else L
    q = scene.q if q is None else q
    M = scene.M
    lu = spla.splu((gamma * sp.identity(M, format="csc") - L.tocsc()))
    bq = beta * q
    PT = None if P is None else P.T.tocsr()
    v = np.ones(M) / M
    lam = 0.0
    for it in range(maxit):
        w = bq * lu.solve(v)
        if PT is not None:
            w = PT @ w
        new = w.sum()
        w /= new
        if abs(new - lam) < tol * new:
            lam = new
            v = w
            break
        lam, v = new, w
    if return_vector:
        return float(lam), v, it + 1
    return float(lam)


def R0_bounds(scene: Scene, beta: float, gamma: float, q=None):
    """Well-mixed value beta<q>/gamma and the frozen value beta max(q)/gamma."""
    q = scene.q if q is None else q
    return beta * q.mean() / gamma, beta * q.max() / gamma


def growth_rate(scene: Scene, beta: float, gamma: float, P=None, L=None, q=None):
    """Principal eigenvalue of the linearised operator L + P^T diag(beta q) - gamma."""
    L = generator(scene) if L is None else L
    q = scene.q if q is None else q
    A = L.toarray()
    F = np.diag(beta * q)
    if P is not None:
        F = P.T @ F
    ev = np.linalg.eigvals(A + F - gamma * np.eye(scene.M))
    return float(np.max(ev.real))


def offspring_mean_by_birth_cell(scene, beta, gamma, P=None, L=None, q=None):
    """Expected number of secondary cases of a case born at each cell
    (column sums of K).  Its uniform average is exactly beta<q>/gamma."""
    K = ngm(scene, beta, gamma, P=P, L=L, q=q)
    return K.sum(axis=0)


# ----------------------------------------------------------------------------
# Spatial branching process: probability of a major outbreak
# ----------------------------------------------------------------------------
def extinction_probability(scene: Scene, beta: float, gamma: float, P=None,
                           L=None, q=None, tol=1e-12, maxit=100000):
    """Extinction probability s(x) of the line started by one case born at x.

    Mean-field branching approximation: along its random walk X_t an infective
    makes contacts as a Poisson process of rate beta q(X_t); each contact starts
    an independent line at a cell drawn from P(.|X_t).  With an exponential
    infectious period the Feynman-Kac formula gives the fixed point

        s = gamma * (gamma I + diag(beta q (1 - P s)) - L)^{-1} 1 .

    Iterated from s = 0 it converges monotonically to the minimal solution.
    Returns s (per walkable cell); the major-outbreak probability of an index
    case at x is 1 - s(x).
    """
    L = generator(scene) if L is None else L
    q = scene.q if q is None else q
    M = scene.M
    Ld = L.toarray()
    s = np.zeros(M)
    one = np.ones(M)
    for it in range(maxit):
        Ps = s if P is None else P @ s
        A = gamma * np.eye(M) + np.diag(beta * q * (1.0 - Ps)) - Ld
        new = gamma * np.linalg.solve(A, one)
        if np.max(np.abs(new - s)) < tol:
            s = new
            break
        s = new
    return np.clip(s, 0.0, 1.0)


# ----------------------------------------------------------------------------
# Mean-field reaction-diffusion ODE on the lattice
# ----------------------------------------------------------------------------
def meanfield_ode(scene: Scene, beta: float, gamma: float, t_eval, I0_site=None,
                  I0=1.0, P=None, L=None, q=None, N=None, rho_ref=None,
                  m_of_t=None, rtol=1e-8, atol=1e-10):
    """Integrate the lattice reaction-diffusion SIR system.

    State: expected numbers S_x, I_x of susceptible / infective people per
    cell.  ``I0`` infectives are placed at ``I0_site`` (index, or an array of
    weights summing to one, or None for uniform); the remaining N - I0 people
    are susceptible and uniformly spread.

        dS/dt = L S - S * (P^T (beta q I)) / rho
        dI/dt = L I + S * (P^T (beta q I)) / rho - gamma I

    with rho = (N-1)/|Omega| (frequency-dependent calibration) unless
    ``rho_ref`` is given (density-dependent: fixed pair hazard).
    Returns dict with t, S, I, R totals and the final per-cell cumulative
    incidence.
    """
    L = generator(scene) if L is None else L
    q = scene.q if q is None else q
    N = scene.N if N is None else N
    M = scene.M
    rho = (N - 1) / M if rho_ref is None else rho_ref
    Lc = L.tocsr()
    PT = None if P is None else P.T.tocsr()
    if I0_site is None:
        w = np.ones(M) / M
    elif np.ndim(I0_site) == 0:
        w = np.zeros(M)
        w[int(I0_site)] = 1.0
    else:
        w = np.asarray(I0_site, float)
        w = w / w.sum()
    I_init = I0 * w
    S_init = (N - I0) * np.ones(M) / M
    y0 = np.concatenate([S_init, I_init, np.zeros(M)])
    bq = beta * q

    def rhs(t, y):
        S, I = y[:M], y[M:2 * M]
        e = bq * I
        if PT is not None:
            e = PT @ e
        m = 1.0 if m_of_t is None else m_of_t(t)
        inf = m * S * e / rho
        return np.concatenate([Lc @ S - inf, Lc @ I + inf - gamma * I, inf])

    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), y0, t_eval=t_eval,
                    method="LSODA", rtol=rtol, atol=atol)
    S = sol.y[:M].sum(axis=0)
    I = sol.y[M:2 * M].sum(axis=0)
    C = sol.y[2 * M:]
    return dict(t=sol.t, S=S, I=I, R=N - S - I, cum_inc_site=C[:, -1],
                I_site=sol.y[M:2 * M], S_site=sol.y[:M])


def wellmixed_final_size(R0: float) -> float:
    """Final-size fraction z solving z = 1 - exp(-R0 z) (0 if R0 <= 1)."""
    if R0 <= 1.0:
        return 0.0
    z = 1.0 - 1.0 / R0 if R0 < 1.5 else 0.9
    for _ in range(200):
        z = z - (z - 1.0 + np.exp(-R0 * z)) / (1.0 - R0 * np.exp(-R0 * z))
    return float(z)


# ----------------------------------------------------------------------------
# Pair-level theory for discrete individuals
# ----------------------------------------------------------------------------
def _pair_operator(scene, beta, gamma, N, L, q, P, rho_ref):
    """A = gamma + H - (L (+) L) on pair states (x, y) -> x*M + y, and the
    diagonal pair hazard h(x, y) = beta q(x) P(y|x) / rho_bar."""
    L = generator(scene) if L is None else L
    q = scene.q if q is None else q
    N = scene.N if N is None else N
    M = scene.M
    rho = (N - 1) / M if rho_ref is None else rho_ref
    Lc = sp.csr_matrix(L)
    I = sp.identity(M, format="csr")
    L2 = sp.kron(Lc, I, format="csr") + sp.kron(I, Lc, format="csr")
    Pk = sp.identity(M, format="csr") if P is None else sp.csr_matrix(P)
    h = (sp.diags(beta * q / rho) @ Pk).toarray().ravel()     # (M*M,)
    A = (gamma * sp.identity(M * M, format="csr") + sp.diags(h) - L2).tocsr()
    return A, h, N, M


def _cg(A, b, rtol=1e-10):
    # A is symmetric positive definite: conjugate gradients with a Jacobi
    # preconditioner (a sparse LU of this 4-D lattice operator would fill in).
    d = A.diagonal()
    Minv = spla.LinearOperator(A.shape, matvec=lambda v: v / d)
    u, info = spla.cg(A, b, rtol=rtol, maxiter=50000, M=Minv)
    if info != 0:
        raise RuntimeError(f"pair solve did not converge (info={info})")
    return u


def pair_infection_probability(scene: Scene, beta: float, gamma: float,
                               N: int | None = None, L=None, q=None, P=None,
                               rho_ref=None):
    """Exact probability p[x, y] that an infective born at cell x infects one
    given other person who is at cell y at that moment (discrete people).

    The pair (X_t, Y_t) is a Markov chain on Omega x Omega with generator
    L (+) L.  While the pair is in state (x, y) the infective transmits at the
    pair hazard h(x, y) = beta q(x) P(y|x) / rho_bar, and the infectious period
    is Exp(gamma).  Feynman-Kac gives

        1 - p = gamma * (gamma + H - L (+) L)^{-1} 1 .

    Infection is a partially absorbing 'defect' on the pair lattice.  The expected number
    of secondary cases of an index case born at x, with the other N-1 people
    uniformly spread and no competition from later cases, is

        R1(x) = (N - 1) * mean_y p[x, y].
    """
    A, h, N, M = _pair_operator(scene, beta, gamma, N, L, q, P, rho_ref)
    u = _cg(A, gamma * np.ones(M * M))
    return 1.0 - u.reshape(M, M)


def pair_ngm(scene: Scene, beta: float, gamma: float, N: int | None = None,
             L=None, q=None, P=None, rho_ref=None, progress=None):
    """Pair-level next-generation matrix for discrete people.

    K1[z, x] = expected number of people infected *at cell z* by a case born
    at cell x, when each of the other N-1 people starts at an independent
    uniform cell and can be infected only once.  Uses M symmetric solves on the
    pair lattice.  Column sums are R1(x); the spectral radius is the
    pair-level reproduction number (it still neglects that a new case starts
    next to its infector and siblings).
    """
    A, h, N, M = _pair_operator(scene, beta, gamma, N, L, q, P, rho_ref)
    H = h.reshape(M, M)
    K1 = np.zeros((M, M))
    for z in range(M):
        b = np.zeros((M, M))
        b[:, z] = H[:, z]
        if not b.any():
            continue
        w = _cg(A, b.ravel(), rtol=1e-9).reshape(M, M)
        K1[z, :] = (N - 1) / M * w.sum(axis=1)
        if progress is not None and z % 50 == 0:
            progress(z, M)
    return K1


def pair_R0_homogeneous(Lx: int, Ly: int, N: int, beta_q: float, gamma: float,
                        D: float, boundary: str = "reflecting"):
    """Closed form for a homogeneous room and the same-cell rule:

        R1 = (beta q / gamma) / (1 + (beta q / rho_bar) * g0),

    where g0 = int_0^inf e^{-gamma t} P(X_t = Y_t | X_0 = Y_0) dt is the
    lattice Green's function of the relative walk at the origin.  Because the
    walk is symmetric, P(X_t = Y_t | X_0 = Y_0 = z) = P_{2t}(z | z), so

        g0 = (1/|Omega|) * sum_k 1 / (gamma + 2 mu_k)

    with mu_k the eigenvalues of -L: for reflecting walls
    mu = 2D(1 - cos(n pi / Lx)) + 2D(1 - cos(m pi / Ly)), n = 0..Lx-1; for
    periodic walls the same with 2 pi n / Lx.  The formula is exact on the
    periodic lattice and uses the cell-averaged g0 for reflecting walls.
    """
    f = np.pi if boundary == "reflecting" else 2 * np.pi
    mx = 2 * D * (1 - np.cos(f * np.arange(Lx) / Lx))
    my = 2 * D * (1 - np.cos(f * np.arange(Ly) / Ly))
    mu = mx[:, None] + my[None, :]
    g0 = np.mean(1.0 / (gamma + 2 * mu))
    rho = (N - 1) / (Lx * Ly)
    return (beta_q / gamma) / (1.0 + beta_q / rho * g0), g0
