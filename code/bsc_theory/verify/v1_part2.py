"""Independent adversarial re-verification of Theorems 1-5 (own implementation; does NOT import scripts/ngm_lattice.py)."""
import numpy as np, scipy.linalg as sla, json, time, sys
from scipy.optimize import minimize
rng = np.random.default_rng(424242)
t0 = time.time()
OUT = {}
def P(*a):
    print(*a, flush=True)

# ---------------- generators (convention: M[r,s] = rate s->r, columns sum to 0) ----------------
def lattice(ny, nx, pobs=0.0):
    while True:
        mask = rng.random((ny, nx)) >= pobs
        idx = -np.ones((ny, nx), int); idx[mask] = np.arange(mask.sum())
        n = int(mask.sum())
        if n < 4: continue
        E = []
        for y in range(ny):
            for x in range(nx):
                if not mask[y, x]: continue
                if x + 1 < nx and mask[y, x + 1]: E.append((idx[y, x], idx[y, x + 1]))
                if y + 1 < ny and mask[y + 1, x]: E.append((idx[y, x], idx[y + 1, x]))
        # connectivity
        adj = np.zeros((n, n), bool)
        for i, j in E: adj[i, j] = adj[j, i] = True
        seen = {0}; st = [0]
        while st:
            k = st.pop()
            for j in np.nonzero(adj[k])[0]:
                if j not in seen: seen.add(int(j)); st.append(int(j))
        if len(seen) == n: return n, np.array(E)

def from_rates(n, src, dst, w):
    M = np.zeros((n, n))
    np.add.at(M, (dst, src), w)
    M -= np.diag(M.sum(axis=0))
    return M

def gen(kind):
    if kind == 'dense':
        n = int(rng.integers(3, 12))
        W = rng.exponential(1.0, (n, n)) * (rng.random((n, n)) < 0.7)
        np.fill_diagonal(W, 0)
        # ensure irreducible: add a directed cycle
        for i in range(n): W[(i + 1) % n, i] += rng.exponential(0.5) + 1e-3
        M = W - np.diag(W.sum(axis=0)); return M, None
    n, E = lattice(int(rng.integers(2, 8)), int(rng.integers(2, 8)), rng.choice([0, 0.1, 0.25]))
    D = np.exp(rng.normal(0, 1.2, n))
    i, j = E[:, 0], E[:, 1]
    if kind == 'sym':
        w = 2 * D[i] * D[j] / (D[i] + D[j])
        M = from_rates(n, np.r_[i, j], np.r_[j, i], np.r_[w, w])
    elif kind == 'dep':
        M = from_rates(n, np.r_[i, j], np.r_[j, i], np.r_[D[i], D[j]])
    elif kind == 'nonrev':
        M = from_rates(n, np.r_[i, j], np.r_[j, i], np.exp(rng.normal(0, 1.5, 2 * len(E))))
    return M, (E, D)

def randq(n):
    c = rng.integers(0, 4)
    if c == 0: q = rng.exponential(1.0, n)
    elif c == 1: q = rng.choice([0.5, 0.8, 1.2, 1.5, 2.0], n)
    elif c == 2: q = rng.exponential(1.0, n) * (rng.random(n) < 0.5); q[rng.integers(n)] += 0.3
    else: q = 1 + 1e-3 * rng.normal(size=n)       # nearly constant (hard case for strictness)
    return np.abs(q)

def sab(A): return np.max(np.linalg.eigvals(A).real)
def statdist(M):
    w, v = np.linalg.eig(M); k = np.argmin(np.abs(w)); p = np.real(v[:, k]); return p / p.sum()
def R0f(M, q, b, g, kap=None):
    n = len(q); V = g * np.eye(n) - M
    if kap is not None: V = V + np.diag(kap)
    K = b * q[:, None] * np.linalg.inv(V)
    ev = np.linalg.eigvals(K); return float(np.max(np.abs(ev))), K

KINDS = ['sym', 'dep', 'nonrev', 'dense']
b3 = 0; v_mu = -1; v_mono = -1; v_lo = -1; v_up = -1; lim = 0; mupos = 1e9
for t in range(400):
    kind = KINDS[t % 4]; M, _ = gen(kind); n = M.shape[0]; b = np.exp(rng.normal()); g = np.exp(rng.normal(-1, 1)); pi = statdist(M)
    kap = rng.exponential(1.0, n) * (rng.random(n) < 0.3)
    if kap.sum() == 0: kap[rng.integers(n)] = 1.0
    mu = -sab(M - np.diag(kap)); mupos = min(mupos, mu)
    q0 = 1.1; R, _ = R0f(M, np.full(n, q0), b, g, kap); b3 = max(b3, abs(R - b * q0 / (g + mu)) / R)
    v_mu = max(v_mu, (mu - pi @ kap) / (pi @ kap))
    kap2 = kap.copy(); kap2[rng.integers(n)] += rng.exponential(); v_mono = max(v_mono, mu - (-sab(M - np.diag(kap2))))
    q = randq(n); R, _ = R0f(M, q, b, g, kap)
    w, vl, vr = sla.eig(M - np.diag(kap), left=True, right=True); k = np.argmax(w.real); u = np.abs(vr[:, k].real); l = np.abs(vl[:, k].real)
    nu = u * l / (u @ l)
    v_lo = max(v_lo, (b * (nu @ q) / (g + mu) - R) / R); v_up = max(v_up, (R - b * q.max() / (g + mu)) / R)
    Dl = np.nonzero(kap > 0)[0]; rest = np.setdiff1d(np.arange(n), Dl)
    if len(rest) > 0:
        mu_abs = -sab(M[np.ix_(rest, rest)]); mu_big = -sab(M - 1e6 * np.diag((kap > 0).astype(float)))
        lim = max(lim, abs(mu_big - mu_abs) / max(mu_abs, 1e-12)); viol_abs = max(globals().get('viol_abs', -1), (mu_big - mu_abs) / max(mu_abs, 1e-12))
P('mu(k=1e6) - mu_abs max rel (should be <= ~1e-9):', viol_abs); P(f'(b) uniform q leaky: rel err {b3:.2e}; mu<=<kappa>_pi violation {v_mu:.2e}; mu>0 min {mupos:.2e}; mono violation {v_mono:.2e}; het bounds violations lower {v_lo:.2e} upper {v_up:.2e}; absorbing limit err {lim:.2e}')
for L in (3, 5, 10, 13, 20, 40, 100):
    w = 1.0; A = -2 * w * np.eye(L) + w * (np.eye(L, k=1) + np.eye(L, k=-1))
    mu = -sab(A); cf = 2 * w * (1 - np.cos(np.pi / (L + 1)))
    P(f'   corridor L={L}: mu={mu:.6f} closed {cf:.6f} pi^2 D/L^2={np.pi**2/L**2:.6f}  ratio mu/(pi^2/L^2)={mu/(np.pi**2/L**2):.3f}')
OUT['T4'] = dict(b3=b3, v_mu=v_mu, v_lo=v_lo, v_up=v_up)

# ------------------------------------------------- Theorem 5
P('== Thm 5 ==')
e_sum = 0; e_fd = 0; e_sym = 0; e_w = 0; e_rev = 0; e_prof = 0
for t in range(240):
    kind = KINDS[t % 4]; M, _ = gen(kind); n = M.shape[0]; q = rng.exponential(1.0, n) + 0.05; b = np.exp(rng.normal()); g = np.exp(rng.normal(-1, 1))
    R, K = R0f(M, q, b, g)
    wv, vl, vr = sla.eig(K, left=True, right=True); k = np.argmax(np.abs(wv)); w = np.abs(vr[:, k].real); z = np.abs(vl[:, k].real)
    w /= w.sum(); e = z * w / (z @ w); e_sum = max(e_sum, abs(e.sum() - 1))
    r = int(rng.integers(n)); h = 1e-6; q2 = q.copy(); q2[r] *= (1 + h); q3 = q.copy(); q3[r] *= (1 - h)
    fd = (np.log(R0f(M, q2, b, g)[0]) - np.log(R0f(M, q3, b, g)[0])) / (np.log(1 + h) - np.log(1 - h))
    e_fd = max(e_fd, abs(fd - e[r]))
    # generalised eigenvector x: bQx = R0 (g - M) x
    V = g * np.eye(n) - M; x = np.linalg.solve(V, w); x = np.abs(x)
    wq = q * x; e_w = max(e_w, np.max(np.abs(wq / wq.sum() - w)))
    if kind in ('sym', 'dep'):
        pi = statdist(M); er = q * x**2 / pi; er /= er.sum(); e_rev = max(e_rev, np.max(np.abs(er - e)))
    # early-time profile
    J = M + np.diag(b * q - g); wj, vj = np.linalg.eig(J); kk = np.argmax(wj.real); u = np.abs(vj[:, kk].real); u /= u.sum()
    I0 = rng.random(n); gap = wj.real[kk] - np.sort(wj.real)[-2]
    if gap > 1e-3:
        T = 40 / gap; It = sla.expm((J - wj.real[kk] * np.eye(n)) * T) @ I0; e_prof = max(e_prof, np.max(np.abs(It / It.sum() - u)))
P(f'|sum e - 1| {e_sum:.2e}; FD elasticity err {e_fd:.2e}; w ~ Qx err {e_w:.2e}; reversible e ~ q x^2/pi err {e_rev:.2e}; profile -> u err {e_prof:.2e}')
OUT['T5'] = dict(e_sum=e_sum, e_fd=e_fd, e_rev=e_rev)
P(f'runtime {time.time()-t0:.1f}s')
json.dump(OUT, open('v1_part2.json', 'w'), indent=1, default=float)
