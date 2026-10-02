"""Re-implementation with separately written code (re-check). No import from ../src/valmod.
Lattice walk propagator by numerical eigendecomposition of the generator."""
import numpy as np
from scipy import stats, optimize

class Lat:
    def __init__(self, nx, ny, ax, ay):
        self.nx, self.ny, self.ax, self.ay = nx, ny, ax, ay
        M = nx*ny
        Lp = np.zeros((M, M))   # unit-D generator (negative semidefinite)
        for y in range(ny):
            for x in range(nx):
                s = y*nx+x
                if x+1 < nx:
                    t = s+1; w = 1/ax**2
                    Lp[s,t]+=w; Lp[t,s]+=w; Lp[s,s]-=w; Lp[t,t]-=w
                if y+1 < ny:
                    t = s+nx; w = 1/ay**2
                    Lp[s,t]+=w; Lp[t,s]+=w; Lp[s,s]-=w; Lp[t,t]-=w
        ev, V = np.linalg.eigh(-Lp)
        ev[np.abs(ev) < 1e-10] = 0.0
        self.mu, self.V, self.M = ev, V, M
    def site(self, x, y):
        return np.asarray(y)*self.nx+np.asarray(x)

def w_m2(L, T):
    L = np.asarray(L, float); x = L*T
    out = np.where(x < 1e-4, T*T*(0.5 - x/6 + x*x/24), (x - 1 + np.exp(-x))/np.where(L>0, L, 1)**2)
    return out

def w_m1(lam, T):
    lam = np.asarray(lam, float); x = 2*lam*T
    return np.where(x < 1e-8, T*(1-x/2), -np.expm1(-x)/np.where(lam>0, 2*lam, 1))

def modew(lat, model, D, kappa, segs):
    w = np.zeros(lat.M)
    for T in np.atleast_1d(segs):
        if model == 'M0':
            w0 = np.zeros(lat.M); w0[lat.mu == 0] = w_m2(np.array([kappa]), T)[0]; w += w0
        elif model == 'M1':
            w += w_m1(D*lat.mu, T)
        else:
            w += w_m2(kappa + D*lat.mu, T)
    return w

def expo(lat, model, D, kappa, segs, src, rec):
    w = modew(lat, model, D, kappa, segs)
    return np.clip(lat.V[np.asarray(rec)] @ (w*lat.V[src]), 0, None)

def expo_mat(lat, model, D, kappa, segs, sites):
    w = modew(lat, model, D, kappa, segs)
    P = lat.V[np.asarray(sites)]
    return np.clip((P*w) @ P.T, 0, None)

def pb(p):
    pmf = np.array([1.0])
    for q in p:
        pmf = np.convolve(pmf, [1-q, q])
    return pmf

def cond_pmf(pn, pf, K):
    a, b = pb(pn), pb(pf)
    out = np.zeros(len(pn)+1)
    for k in range(len(pn)+1):
        if 0 <= K-k <= len(pf):
            out[k] = a[k]*b[K-k]
    return out/out.sum()

def prof_eps(X, K):
    f = lambda le: np.sum(-np.expm1(-np.exp(le)*X)) - K
    lo = np.log(K/X.sum()) - 2; hi = lo + 4
    while f(hi) < 0:
        hi += 4
        if hi > 600: return np.inf
    return np.exp(optimize.brentq(f, lo, hi, xtol=1e-12))

def two_sided(pmf, k):
    return float(pmf[pmf <= pmf[k]*(1+1e-9)].sum())
