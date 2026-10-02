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
# ------------------------------------------------- Theorem 1
P('== Thm 1 ==')
mx_root = 0; sign_fail = 0; cnt = 0; mx_lap = 0
for t in range(600):
    kind = KINDS[t % 4]; M, _ = gen(kind); n = M.shape[0]; q = randq(n)
    b = np.exp(rng.normal(0, 1)); g = np.exp(rng.normal(-1, 1))
    R, K = R0f(M, q, b, g)
    mx_root = max(mx_root, abs(sab(M + np.diag(b * q / R - g))))
    lam = sab(M + np.diag(b * q - g))
    if abs(R - 1) > 1e-9 and np.sign(R - 1) != np.sign(lam): sign_fail += 1
    cnt += 1
    if t < 40:  # Laplace transform of propagator via eigen-decomposition free quadrature: integral e^{-gt} e^{Mt} dt
        from scipy.integrate import quad_vec
        Vinv = np.linalg.inv(g * np.eye(n) - M)
        I, err = quad_vec(lambda s: np.exp(-g * s) * sla.expm(M * s), 0, 60 / g, epsrel=1e-9)
        mx_lap = max(mx_lap, np.max(np.abs(I - Vinv) / Vinv))
P(f'instances {cnt}: max|s(M+bQ/R0-g)| = {mx_root:.2e}; sign failures {sign_fail}; Laplace-vs-resolvent max rel err {mx_lap:.2e}')
OUT['T1'] = dict(root=mx_root, sign_fail=sign_fail, laplace=mx_lap)

# ------------------------------------------------- Theorem 1' and Cor 1.3 (reversible)
P("== Thm 1' / Cor 1.3 ==")
mx_ray = 0; diag_excess = -1; nonmono = 0; mx_full = 0; asym = 0
for t in range(300):
    kind = ['sym', 'dep'][t % 2]; M, (E, D) = gen(kind); n = M.shape[0]; q = randq(n)
    b = np.exp(rng.normal(0, 1)); g = np.exp(rng.normal(-1, 1))
    pi = statdist(M)
    if kind == 'dep':
        assert np.allclose(pi, (1 / D) / (1 / D).sum(), rtol=1e-8), 'pi not ~1/D'
    else:
        assert np.allclose(pi, 1 / n)
    A = M * np.sqrt(pi)[None, :] / np.sqrt(pi)[:, None]
    asym = max(asym, np.max(np.abs(A - A.T)))
    A = (A + A.T) / 2
    R, K = R0f(M, q, b, g)
    Rray = sla.eigh(np.diag(b * q), g * np.eye(n) - A, eigvals_only=True)[-1]
    mx_ray = max(mx_ray, abs(Rray - R) / R)
    mu, phi = np.linalg.eigh(-A)
    Lam = phi.T @ (q[:, None] * phi)
    Rm = b * Lam / np.sqrt(np.outer(g + mu, g + mu))
    diag_excess = max(diag_excess, (np.max(np.diag(Rm)) - R) / R)
    prev = -1
    for N in range(1, n + 1):
        val = np.linalg.eigvalsh(Rm[:N, :N])[-1]
        if val < prev - 1e-11 * R: nonmono += 1
        prev = val
    mx_full = max(mx_full, abs(prev - R) / R)
    # first Galerkin = b<q>_pi/g
    assert abs(Rm[0, 0] - b * (pi @ q) / g) < 1e-9 * R
P(f'Rayleigh vs NGM rel err {mx_ray:.2e}; asym of Pi^-1/2 M Pi^1/2 {asym:.2e}; max rel (diag modal max - R0) {diag_excess:.2e}; Galerkin non-monotone {nonmono}; full basis err {mx_full:.2e}')
OUT['T1p'] = dict(ray=mx_ray, diag_excess=diag_excess, nonmono=nonmono)

# ------------------------------------------------- Theorem 2
P('== Thm 2 ==')
vl = vu = -1; gl = gu = 1e9; cq = 0; vl_lam = vu_lam = -1
for t in range(3000):
    kind = KINDS[t % 4]; M, _ = gen(kind); n = M.shape[0]; q = randq(n)
    b = np.exp(rng.normal(0, 1)); g = np.exp(rng.normal(-1, 1)); pi = statdist(M)
    R, _ = R0f(M, q, b, g); lo = b * (pi @ q) / g; hi = b * q.max() / g
    lam = sab(M + np.diag(b * q - g))
    vl = max(vl, (lo - R) / R); vu = max(vu, (R - hi) / R)
    vl_lam = max(vl_lam, (b * (pi @ q) - g) - lam); vu_lam = max(vu_lam, lam - (b * q.max() - g))
    if q.max() - q.min() > 1e-9:
        gl = min(gl, (R - lo) / R / ((q.max() - q.min()) / q.max()) ** 2); gu = min(gu, (hi - R) / R)
    if t % 10 == 0:
        R2, _ = R0f(M, np.full(n, 1.3), b, g); cq = max(cq, abs(R2 - b * 1.3 / g) / R2)
P(f'max rel violation lower {vl:.2e}, upper {vu:.2e}; lambda1 bounds violation {vl_lam:.2e}, {vu_lam:.2e}; min normalised lower gap {gl:.2e}; min upper gap {gu:.2e}; const q err {cq:.2e}')
OUT['T2'] = dict(vl=vl, vu=vu)

# adversarial: try to MINIMISE (R0 - lower)/lower over q and nonreversible rates for small n
def adv_T2():
    worst = 1e9
    for trial in range(60):
        n = int(rng.integers(3, 6))
        def f(z):
            W = np.exp(z[:n * n]).reshape(n, n); np.fill_diagonal(W, 0)
            M = W - np.diag(W.sum(axis=0)); q = np.exp(z[n * n:n * n + n]); g = 0.3; b = 1.0
            pi = statdist(M); R, _ = R0f(M, q, b, g)
            return R / (b * (pi @ q) / g) - 1.0
        z0 = rng.normal(0, 1, n * n + n)
        res = minimize(f, z0, method='Nelder-Mead', options=dict(maxiter=1500, xatol=1e-6, fatol=1e-14))
        worst = min(worst, res.fun)
    return worst
w2 = adv_T2(); P(f'adversarial min of R0/(b<q>/g)-1 over nonreversible M, q (Nelder-Mead, 60 starts): {w2:.3e}  (must be >= 0)')
OUT['T2_adv'] = w2

# ------------------------------------------------- Theorem 3
P('== Thm 3 ==')
mx_inc = -1; nonstrict = 0; liminf = 0; lim0 = 0; conv = -1
Ds = np.logspace(-4, 4, 49)
for t in range(400):
    kind = KINDS[t % 4]; M0, _ = gen(kind); n = M0.shape[0]; q = randq(n)
    b = np.exp(rng.normal(0, 1)); g = np.exp(rng.normal(-1, 1)); pi = statdist(M0)
    sc = np.max(np.abs(np.diag(M0)))
    Rs = np.array([R0f(D * M0 * g / sc, q, b, g)[0] for D in Ds])
    d = np.diff(Rs) / Rs[:-1]; mx_inc = max(mx_inc, d.max())
    if q.max() - q.min() > 1e-6 * q.max():
        nonstrict += int(np.sum(d[8:40] >= 0))
    lam = np.array([sab(D * M0 + np.diag(b * q - g)) for D in np.linspace(0.05, 5, 30)])
    conv = max(conv, np.max(-(lam[2:] - 2 * lam[1:-1] + lam[:-2])), np.max(np.diff(lam)))
    Rbig = R0f(1e7 * M0 * g / sc, q, b, g)[0]; Rsm = R0f(1e-8 * M0 * g / sc, q, b, g)[0]
    liminf = max(liminf, abs(Rbig - b * (pi @ q) / g) / Rbig); lim0 = max(lim0, abs(Rsm - b * q.max() / g) / Rsm)
P(f'max rel increase of R0 in D {mx_inc:.2e}; strict-decrease failures (mid-range) {nonstrict}; D->inf err {liminf:.2e}; D->0 err {lim0:.2e}; lambda1 convexity/monotone violation {conv:.2e}')
OUT['T3'] = dict(mx_inc=mx_inc, nonstrict=nonstrict)

def adv_T3():
    worst = -1e9
    for trial in range(80):
        n = int(rng.integers(3, 6))
        def f(z):
            W = np.exp(z[:n * n]).reshape(n, n); np.fill_diagonal(W, 0)
            M0 = W - np.diag(W.sum(axis=0)); q = np.exp(z[n * n:n * n + n]); g = 0.3; b = 1.0
            D1 = np.exp(z[-2]); D2 = D1 * (1 + np.exp(z[-1]))
            return -(R0f(D2 * M0, q, b, g)[0] - R0f(D1 * M0, q, b, g)[0])   # want to maximise increase
        z0 = rng.normal(0, 1.5, n * n + n + 2)
        res = minimize(f, z0, method='Nelder-Mead', options=dict(maxiter=2000, xatol=1e-6, fatol=1e-15))
        worst = max(worst, -res.fun)
    return worst
w3 = adv_T3(); P(f'adversarial max of R0(D2)-R0(D1), D2>D1, nonreversible (80 starts): {w3:.3e}  (must be <= 0)')
OUT['T3_adv'] = w3

# ------------------------------------------------- Prop 3.3
P('== Prop 3.3 ==')
ratL = []; ratS = []; Cdiff = 0; Cneg = 0
for t in range(300):
    kind = KINDS[t % 4]; M0, _ = gen(kind); n = M0.shape[0]; q = randq(n)
    if q.max() - q.min() < 1e-2: q = q + rng.random(n)
    b = 1.0; g = 0.5; pi = statdist(M0); qb = pi @ q
    G = np.linalg.pinv(-M0)  # NOT the group inverse in general; build group inverse explicitly:
    Pj = np.outer(pi, np.ones(n)); Gsharp = np.linalg.inv(-M0 + Pj) - Pj
    x = Gsharp @ (q * pi - qb * pi); C = q @ x
    if C < -1e-12: Cneg += 1
    if kind == 'sym':
        dq = q - q.mean(); Cdiff = max(Cdiff, abs(C - dq @ G @ dq / n) / abs(C))
    sc = np.max(np.abs(np.diag(M0)))
    def errL(D): return abs(R0f(D * M0, q, b, g)[0] - (b * qb / g + b * C / (qb * D)))
    e1, e2 = errL(300 * g / sc * 10), errL(3000 * g / sc * 10)
    if e2 > 1e-13: ratL.append(e1 / e2)
    Z = np.nonzero(q >= q.max() - 1e-12)[0]; muZ = -sab(M0[np.ix_(Z, Z)])
    def errS(D): return abs(R0f(D * M0, q, b, g)[0] - (b * q.max() / g) * (1 - D * muZ / g))
    e1, e2 = errS(1e-2 * g / sc), errS(1e-3 * g / sc)
    if e2 > 1e-14: ratS.append(e1 / e2)
P(f'large-D error ratio per decade: median {np.median(ratL):.1f} [{np.percentile(ratL,5):.1f},{np.percentile(ratL,95):.1f}] n={len(ratL)}; small-D: median {np.median(ratS):.1f} [{np.percentile(ratS,5):.1f},{np.percentile(ratS,95):.1f}] n={len(ratS)}; C<0 count {Cneg}; sym closed form diff {Cdiff:.2e}')
OUT['P33'] = dict(L=float(np.median(ratL)), S=float(np.median(ratS)), Cneg=Cneg)

# ------------------------------------------------- Prop 3.4 zone bound, Thm 3.5 edge monotone
P('== Prop 3.4 / Thm 3.5 ==')
zv = -1; nz = 0
for t in range(600):
    kind = KINDS[t % 4]; M, _ = gen(kind); n = M.shape[0]; q = randq(n); b = np.exp(rng.normal()); g = np.exp(rng.normal(-1, 1))
    R, _ = R0f(M, q, b, g)
    for _ in range(3):
        k = int(rng.integers(1, n)); Z = rng.choice(n, k, replace=False)
        muZ = -sab(M[np.ix_(Z, Z)]); bd = b * q[Z].min() / (g + muZ); zv = max(zv, (bd - R) / R); nz += 1
        assert muZ >= -1e-12
P(f'zone bound: {nz} zones, max rel violation {zv:.2e}')
ev = -1; ne = 0; dep_up = dep_dn = 0
for t in range(400):
    M, (E, D) = gen('sym'); n = M.shape[0]; q = randq(n); b = np.exp(rng.normal()); g = np.exp(rng.normal(-1, 1))
    R, _ = R0f(M, q, b, g); lam = sab(M + np.diag(b * q - g))
    for _ in range(3):
        e = E[rng.integers(len(E))]; dw = rng.exponential(1.0)
        M2 = M.copy(); i, j = e; M2[i, j] += dw; M2[j, i] += dw; M2[i, i] -= dw; M2[j, j] -= dw
        R2, _ = R0f(M2, q, b, g); ev = max(ev, (R2 - R) / R, sab(M2 + np.diag(b * q - g)) - lam); ne += 1
P(f'edge monotonicity (symmetric): {ne} trials, max increase {ev:.2e}')
# does single-edge monotonicity FAIL for reversible-but-not-symmetric / nonreversible? (the theorem claims only symmetric)
for t in range(300):
    n, E = lattice(int(rng.integers(2, 7)), int(rng.integers(2, 7))); D = np.exp(rng.normal(0, 1.2, n)); i, j = E[:, 0], E[:, 1]
    M = from_rates(n, np.r_[i, j], np.r_[j, i], np.r_[D[i], D[j]]); q = randq(n); b = 1.0; g = 0.3
    D2 = D.copy(); D2[rng.integers(n)] *= 2; M2 = from_rates(n, np.r_[i, j], np.r_[j, i], np.r_[D2[i], D2[j]])
    d = R0f(M2, q, b, g)[0] - R0f(M, q, b, g)[0]
    if d > 1e-12: dep_up += 1
    elif d < -1e-12: dep_dn += 1
P(f'departure rule, raise D at one cell: R0 up {dep_up}, down {dep_dn} of 300')
OUT['P34'] = zv; OUT['T35'] = ev

# ------------------------------------------------- Theorem 4
P('== Thm 4 ==')
a_err = 0; n_err = 0
for t in range(400):
    kind = KINDS[t % 4]; M, _ = gen(kind); n = M.shape[0]; b = np.exp(rng.normal()); g = np.exp(rng.normal(-1, 1)); q0 = rng.exponential() + 0.1
    R, K = R0f(M, np.full(n, q0), b, g); a_err = max(a_err, abs(R - b * q0 / g) / R); n_err = max(n_err, np.max(np.abs(K.sum(axis=0) - b * q0 / g)) / R)
P(f'(a) uniform q: rel err {a_err:.2e}; offspring number n(s) const err {n_err:.2e}')
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
        mu_abs = -sab(M[np.ix_(rest, rest)]); mu_big = -sab(M - 1e9 * np.diag((kap > 0).astype(float)))
        lim = max(lim, abs(mu_big - mu_abs) / max(mu_abs, 1e-12)); assert mu_big <= mu_abs * (1 + 1e-9) + 1e-12
P(f'(b) uniform q leaky: rel err {b3:.2e}; mu<=<kappa>_pi violation {v_mu:.2e}; mu>0 min {mupos:.2e}; mono violation {v_mono:.2e}; het bounds violations lower {v_lo:.2e} upper {v_up:.2e}; absorbing limit err {lim:.2e}')
for L in (3, 5, 10, 13, 20, 40, 100):
    w = 1.0; A = -2 * w * np.eye(L) + w * (np.eye(L, k=1) + np.eye(L, k=-1))
    mu = -sab(A); cf = 2 * w * (1 - np.cos(np.pi / (L + 1)))
    P(f'   corridor L={L}: mu={mu:.6f} closed {cf:.6f} pi^2 D/L^2={np.pi**2/L**2:.6f}  ratio mu/(pi^2/L^2)={mu/(np.pi**2/L**2):.3f}')
OUT['T4'] = dict(a=a_err, b3=b3, v_mu=v_mu, v_lo=v_lo, v_up=v_up)

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
json.dump(OUT, open(__file__.replace('.py', '.json'), 'w'), indent=1, default=float)
