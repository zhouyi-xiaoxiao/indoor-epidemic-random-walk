"""ngm_lattice.py -- next-generation-matrix (NGM) tools for the lattice SIR
reaction-diffusion model.

Conventions
-----------
* Cells: the accessible cells of a 2-D boolean mask of shape (ny, nx), indexed
  0..n-1 in row-major order.  Obstacles / walls are simply not accessible.
* Movement generator M (n x n, sparse CSR):
      M[i, j] = rate at which ONE person at cell j jumps to cell i   (i != j)
      M[j, j] = -(total jump rate out of j)
  so that for occupation numbers N (column vector)  dN/dt = M N ; the columns of M
  sum to zero (people are conserved: reflecting walls), i.e. 1^T M = 0.
* Leaky room: exit (door) rates kappa_r >= 0 remove people:  M_kappa = M - diag(kappa).
* Linearised infected dynamics at the disease-free state:  dI/dt = (M + beta Q - gamma I) I,
  Q = diag(q).  Transmission F = beta Q, transition V = gamma I - M (+ diag kappa).
  NGM  K = F V^{-1} = beta Q (gamma I - M)^{-1}.

Two hopping conventions:
  'symmetric'  (harmonic mean): per-neighbour rate w_ij = w_ji = h(D_i, D_j)/a^2, h = harmonic
               mean (Fickian limit div(D grad rho)); stationary distribution uniform.
  'departure'  (departure-site rate, nearest neighbours): rate i->j = D_i/a^2 per accessible
               neighbour (Ito limit Laplacian(D rho)); reversible with pi_i proportional to 1/D_i.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components

# --------------------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------------------


class Lattice:
    """Accessible cells of a 2-D mask with nearest-neighbour (von Neumann) edges."""

    def __init__(self, mask):
        mask = np.asarray(mask, dtype=bool)
        self.mask = mask
        self.ny, self.nx = mask.shape
        self.index = -np.ones(mask.shape, dtype=int)
        self.cells = np.argwhere(mask)  # (y, x), row-major
        self.n = len(self.cells)
        self.index[mask] = np.arange(self.n)
        edges = []
        for dy, dx in [(0, 1), (1, 0)]:
            both = mask[: self.ny - dy, : self.nx - dx] & mask[dy:, dx:]
            ys, xs = np.nonzero(both)
            i = self.index[ys, xs]
            j = self.index[ys + dy, xs + dx]
            edges.append(np.stack([i, j], axis=1))
        self.edges = np.concatenate(edges) if edges else np.zeros((0, 2), int)

    # convenience -------------------------------------------------------------------
    def vec(self, grid):
        """2-D array (ny, nx) -> vector over accessible cells."""
        return np.asarray(grid, dtype=float)[self.mask]

    def grid(self, v, fill=np.nan):
        g = np.full(self.mask.shape, fill, dtype=float)
        g[self.mask] = v
        return g

    def adjacency(self, weights=None):
        w = np.ones(len(self.edges)) if weights is None else weights
        A = sp.coo_matrix((w, (self.edges[:, 0], self.edges[:, 1])), shape=(self.n, self.n))
        return (A + A.T).tocsr()

    def n_components(self, weights=None):
        A = self.adjacency(weights)
        if weights is not None:
            A.data = (A.data > 0).astype(float)
            A.eliminate_zeros()
        return connected_components(A, directed=False)[0]

    def edge_index(self, cell_a, cell_b):
        """index of the undirected edge between two (y,x) cells, or -1."""
        ia, ib = self.index[cell_a], self.index[cell_b]
        hit = np.nonzero(((self.edges[:, 0] == ia) & (self.edges[:, 1] == ib))
                         | ((self.edges[:, 0] == ib) & (self.edges[:, 1] == ia)))[0]
        return int(hit[0]) if len(hit) else -1


# --------------------------------------------------------------------------------------
# generators
# --------------------------------------------------------------------------------------


def generator_from_edge_rates(n, edges, w_fwd, w_bwd):
    """Generator with rate w_fwd for edges[:,0] -> edges[:,1] and w_bwd for the reverse."""
    src = np.concatenate([edges[:, 0], edges[:, 1]])
    dst = np.concatenate([edges[:, 1], edges[:, 0]])
    val = np.concatenate([w_fwd, w_bwd]).astype(float)
    keep = val > 0
    Moff = sp.coo_matrix((val[keep], (dst[keep], src[keep])), shape=(n, n)).tocsr()
    out = np.asarray(Moff.sum(axis=0)).ravel()  # column sums = total out-rate of each source
    return (Moff - sp.diags(out)).tocsr()


def interface_mean(Di, Dj, kind="harmonic"):
    if kind == "harmonic":
        return 2.0 * Di * Dj / (Di + Dj)
    if kind == "arithmetic":
        return 0.5 * (Di + Dj)
    if kind == "geometric":
        return np.sqrt(Di * Dj)
    raise ValueError(kind)


def generator_symmetric(lat: Lattice, D, a=1.0, kind="harmonic", edge_scale=None):
    """Symmetric per-neighbour rates h(D_i,D_j)/a^2.
    edge_scale (len = #edges) multiplies individual edge rates (0 = partition wall on that edge)."""
    D = np.asarray(D, float)
    i, j = lat.edges[:, 0], lat.edges[:, 1]
    w = interface_mean(D[i], D[j], kind) / a**2
    if edge_scale is not None:
        w = w * edge_scale
    return generator_from_edge_rates(lat.n, lat.edges, w, w)


def generator_departure(lat: Lattice, D, a=1.0, edge_scale=None):
    """Departure-site rule, uniform over accessible neighbours, per-neighbour rate D_departure/a^2."""
    D = np.asarray(D, float)
    i, j = lat.edges[:, 0], lat.edges[:, 1]
    wf, wb = D[i] / a**2, D[j] / a**2
    if edge_scale is not None:
        wf, wb = wf * edge_scale, wb * edge_scale
    return generator_from_edge_rates(lat.n, lat.edges, wf, wb)


def stationary(M):
    """Stationary distribution pi (M pi = 0, sum pi = 1) of an irreducible generator."""
    n = M.shape[0]
    A = M.toarray() if sp.issparse(M) else np.asarray(M)
    # replace one equation by normalisation
    B = A.copy()
    B[-1, :] = 1.0
    rhs = np.zeros(n)
    rhs[-1] = 1.0
    pi = np.linalg.solve(B, rhs)
    return pi


def stationary_sparse(M):
    n = M.shape[0]
    B = sp.lil_matrix(M)
    B[n - 1, :] = np.ones(n)
    rhs = np.zeros(n)
    rhs[-1] = 1.0
    return spla.spsolve(B.tocsc(), rhs)


# --------------------------------------------------------------------------------------
# dense spectral quantities (n up to a few thousand)
# --------------------------------------------------------------------------------------


def _dense(M):
    return M.toarray() if sp.issparse(M) else np.asarray(M, float)


def spectral_abscissa(J):
    """s(J) = max real part of the spectrum (dense)."""
    ev = np.linalg.eigvals(_dense(J))
    return float(np.max(ev.real))


def perron_metzler(J):
    """Perron root s(J) of an irreducible Metzler matrix with positive right/left eigenvectors."""
    J = _dense(J)
    ev, VL, VR = sla.eig(J, left=True, right=True)
    k = int(np.argmax(ev.real))
    s = ev[k]
    u = np.real(VR[:, k])
    l = np.real(VL[:, k])
    u = u / u.sum()
    l = l / l.sum()
    return float(s.real), u, l, float(abs(s.imag))


def ngm_dense(M, q, beta, gamma, kappa=None):
    """K = beta Q (gamma I - M + diag kappa)^{-1} (dense)."""
    A = _dense(M)
    n = A.shape[0]
    V = gamma * np.eye(n) - A
    if kappa is not None:
        V = V + np.diag(kappa)
    Vinv = np.linalg.inv(V)
    return beta * np.asarray(q, float)[:, None] * Vinv, Vinv


def R0_dense(M, q, beta, gamma, kappa=None, vectors=False):
    """Spectral radius of the NGM, with right (w) and left (z) Perron vectors."""
    K, Vinv = ngm_dense(M, q, beta, gamma, kappa)
    ev, VL, VR = sla.eig(K, left=True, right=True)
    k = int(np.argmax(np.abs(ev)))
    R0 = float(ev[k].real)
    if not vectors:
        return R0
    w = np.abs(np.real(VR[:, k]))
    z = np.abs(np.real(VL[:, k]))
    return R0, w / w.sum(), z / z.sum(), K


def growth_rate_dense(M, q, beta, gamma, kappa=None):
    """lambda_1 = s(M - diag kappa + beta Q - gamma I) and its Perron vectors (u right, l left)."""
    J = _dense(M) + np.diag(beta * np.asarray(q, float) - gamma)
    if kappa is not None:
        J = J - np.diag(kappa)
    return perron_metzler(J)


# --------------------------------------------------------------------------------------
# sparse spectral quantities
# --------------------------------------------------------------------------------------


def R0_sparse(M, q, beta, gamma, kappa=None, tol=1e-12, return_vec=False):
    """rho(beta Q V^{-1}) via ARPACK on the operator x -> beta q * V^{-1} x (sparse LU of V)."""
    n = M.shape[0]
    V = gamma * sp.identity(n, format="csc") - M.tocsc()
    if kappa is not None:
        V = V + sp.diags(kappa)
    lu = spla.splu(V.tocsc())
    q = np.asarray(q, float)
    op = spla.LinearOperator((n, n), matvec=lambda x: beta * q * lu.solve(np.asarray(x, float)), dtype=float)
    vals, vecs = spla.eigs(op, k=min(3, n - 2), which="LM", tol=tol, maxiter=20000)
    k = int(np.argmax(np.abs(vals)))
    R0 = float(vals[k].real)
    if return_vec:
        w = np.abs(np.real(vecs[:, k]))
        return R0, w / w.sum()
    return R0


def R0_reversible_sparse(M, q, beta, gamma, pi=None, kappa=None, tol=1e-12, return_vec=False):
    """For a reversible generator (detailed balance w.r.t. pi): symmetric generalised eigenproblem
       beta Q y = R (gamma I - S) y,  S = Pi^{-1/2} M Pi^{1/2} symmetric.  Uses eigsh (Lanczos)."""
    n = M.shape[0]
    if pi is None:
        pi = np.full(n, 1.0 / n)
    sq = np.sqrt(pi)
    S = sp.diags(1.0 / sq) @ M @ sp.diags(sq)
    S = 0.5 * (S + S.T)
    B = gamma * sp.identity(n) - S
    if kappa is not None:
        B = B + sp.diags(kappa)
    A = sp.diags(beta * np.asarray(q, float))
    vals, vecs = spla.eigsh(A.tocsc(), k=1, M=B.tocsc(), which="LA", tol=tol)
    R0 = float(vals[0])
    if return_vec:
        y = np.abs(vecs[:, 0])
        return R0, y  # y in symmetrised coordinates
    return R0


def growth_rate_sparse(M, q, beta, gamma, kappa=None, tol=1e-12, return_vec=False):
    """s(J), J = M - diag kappa + beta Q - gamma I, by shift-invert with sigma above s(J)."""
    n = M.shape[0]
    q = np.asarray(q, float)
    J = M.tocsc() + sp.diags(beta * q - gamma)
    if kappa is not None:
        J = J - sp.diags(kappa)
    sigma = beta * q.max() - gamma + 1e-3 * max(1.0, abs(beta * q.max()))
    vals, vecs = spla.eigs(J, k=1, sigma=sigma, which="LM", tol=tol, maxiter=20000)
    s = float(vals[0].real)
    if return_vec:
        u = np.abs(np.real(vecs[:, 0]))
        return s, u / u.sum()
    return s


# --------------------------------------------------------------------------------------
# derived quantities used in the theorems
# --------------------------------------------------------------------------------------


def bounds(q, beta, gamma, pi):
    q = np.asarray(q, float)
    return beta * float(pi @ q) / gamma, beta * float(q.max()) / gamma


def root_characterisation(M, q, beta, gamma, R, kappa=None):
    """s(M - diag kappa + beta Q / R - gamma I): zero exactly at R = R0 (Theorem 1)."""
    J = _dense(M) + np.diag(beta * np.asarray(q, float) / R - gamma)
    if kappa is not None:
        J = J - np.diag(kappa)
    return spectral_abscissa(J)


def mu_leak(M, kappa):
    """mu_kappa = -s(M - diag kappa) (dense): principal decay rate of the leaky room."""
    return -spectral_abscissa(_dense(M) - np.diag(kappa))


def mu_zone(M, Z):
    """mu_Z = -s(M[Z,Z]) (dense): exit rate of the principal quasi-stationary mode of zone Z
    (reflecting at walls, absorbing where Z borders other accessible cells)."""
    A = _dense(M)[np.ix_(Z, Z)]
    return -spectral_abscissa(A)


def large_D_coefficient(M0, q, pi=None):
    """C in  R0(D) = beta<q>_pi/gamma + beta C /(D <q>_pi) + O(D^-2)   (M = D M0).
    C = q^T x,  (-M0) x = Q pi - <q>_pi pi,  1^T x = 0."""
    A = _dense(M0)
    n = A.shape[0]
    if pi is None:
        pi = stationary(A)
    q = np.asarray(q, float)
    qbar = float(pi @ q)
    y = q * pi - qbar * pi
    # solve (-A) x = y with 1^T x = 0 via bordered system
    Bm = np.zeros((n + 1, n + 1))
    Bm[:n, :n] = -A
    Bm[:n, n] = pi  # Lagrange column (range of -A is {1^T y = 0}; pi spans its complement direction)
    Bm[n, :n] = 1.0
    rhs = np.concatenate([y, [0.0]])
    sol = np.linalg.solve(Bm, rhs)
    x = sol[:n]
    return float(q @ x), qbar


def ritz_R0(M_sym, q, beta, gamma, n_modes):
    """Galerkin/Ritz R0 in the span of the n_modes slowest eigenvectors phi_n of -M (symmetric M).
    Returns (R0_ritz, diag_terms) where diag_terms[n] = beta<phi_n|q|phi_n>/(gamma+mu_n)
    are lower bounds for R0 (largest diagonal modal term)."""
    A = _dense(M_sym)
    mu, phi = np.linalg.eigh(-A)  # ascending, mu_0 = 0
    q = np.asarray(q, float)
    P = phi[:, :n_modes]
    Lam = P.T @ (q[:, None] * P)  # Lambda_nm = <phi_n|q|phi_m>
    Bd = gamma + mu[:n_modes]
    Ms = beta * Lam / np.sqrt(np.outer(Bd, Bd))
    R_ritz = float(np.linalg.eigvalsh(Ms)[-1])
    diag_terms = beta * np.diag(Lam) / Bd
    return R_ritz, diag_terms, mu, phi


def elasticity(w, z):
    """e_r = z_r w_r / (z . w): elasticity of R0 w.r.t. q_r; sums to 1."""
    e = z * w
    return e / e.sum()
