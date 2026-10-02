"""Recomputation, with separately written code, of the worked examples. Own generator/eigen code; only the zone arrays (data) come from scripts/layouts.py."""
import sys, numpy as np, scipy.linalg as sla, scipy.sparse as sp, scipy.sparse.linalg as spla, json, time
sys.path.insert(0, '../scripts')
import layouts as lay
B_, G_ = 0.5, 0.14
t0 = time.time()
def P(*a): print(*a, flush=True)

def build(z, Dmap, qmap, rule='sym', D0=1.0, closed=None, return_all=False):
    ny, nx = z.shape; mask = z != 'B'; idx = -np.ones(z.shape, int); idx[mask] = np.arange(mask.sum()); n = int(mask.sum())
    D = np.zeros(n); q = np.zeros(n); zone = np.empty(n, dtype='<U1')
    for y in range(ny):
        for x in range(nx):
            if mask[y, x]:
                D[idx[y, x]] = Dmap[z[y, x]] * D0; q[idx[y, x]] = qmap[z[y, x]]; zone[idx[y, x]] = z[y, x]
    M = np.zeros((n, n)); edges = []
    for y in range(ny):
        for x in range(nx):
            if not mask[y, x]: continue
            for dy, dx in ((0, 1), (1, 0)):
                yy, xx = y + dy, x + dx
                if yy < ny and xx < nx and mask[yy, xx]:
                    i, j = idx[y, x], idx[yy, xx]; edges.append((i, j))
    for k, (i, j) in enumerate(edges):
        if closed is not None and closed[k]: continue
        if rule == 'sym':
            w = 2 * D[i] * D[j] / (D[i] + D[j]); M[j, i] += w; M[i, j] += w
        else:
            M[j, i] += D[i]; M[i, j] += D[j]     # M[dest, src] = rate src->dest = D_src
    M -= np.diag(M.sum(axis=0))
    if return_all: return M, q, zone, np.array(edges), idx, mask, D
    return M, q, zone

def R0(M, q, b=B_, g=G_, vec=False):
    n = len(q); K = b * q[:, None] * np.linalg.inv(g * np.eye(n) - M)
    if not vec: return float(np.max(np.abs(np.linalg.eigvals(K))))
    w, vl, vr = sla.eig(K, left=True, right=True); k = np.argmax(np.abs(w))
    return float(w[k].real), np.abs(vr[:, k].real), np.abs(vl[:, k].real)
def sab(A): return float(np.max(np.linalg.eigvals(A).real))
def stat(M):
    n = M.shape[0]; A = M.copy(); A[-1, :] = 1; r = np.zeros(n); r[-1] = 1; return np.linalg.solve(A, r)
def connected(M):
    n = M.shape[0]; A = (np.abs(M) > 0); seen = np.zeros(n, bool); seen[0] = True; st = [0]
    while st:
        k = st.pop()
        for j in np.nonzero(A[:, k] | A[k, :])[0]:
            if not seen[j]: seen[j] = True; st.append(j)
    return seen.all()

RES = {}
P('== E1 two-zone ==')
for kind in ('side', 'aisles', 'grid', 'scatter'):
    z = lay.two_zone(kind); u, c = np.unique(z, return_counts=True)
    for rule in ('sym', 'dep'):
        M, q, zone = build(z, lay.D_TWO, lay.Q_TWO, rule); pi = stat(M)
        r = R0(M, q); lo = B_ * (pi @ q) / G_; lam = sab(M + np.diag(B_ * q - G_))
        P(f'  {kind:8s} {rule}: counts {dict(zip(u, c))} <q>_pi={pi@q:.4f} lower {lo:.3f} R0 {r:.4f} = {r/(B_/G_):.4f} b/g  lam1 {lam:.4f}')
        RES[f'E1_{kind}_{rule}'] = r

P('== E2 office ==')
for name, z in (('O1', lay.office_O1()), ('O2', lay.office_O2())):
    u, c = np.unique(z, return_counts=True); P(f'  {name} counts {dict(zip(u, c))}')
    for rule in ('sym', 'dep'):
        M, q, zone = build(z, lay.D_OFFICE, lay.Q_OFFICE, rule); n = len(q); assert connected(M)
        pi = stat(M); r, w, zl = R0(M, q, vec=True); lam = sab(M + np.diag(B_ * q - G_))
        P(f'  {name} {rule}: n={n} <q>_pi={pi@q:.4f} bounds [{B_*(pi@q)/G_:.3f},{B_*q.max()/G_:.3f}] R0={r:.4f} lam1={lam:.4f} doubling {np.log(2)/lam:.2f} d')
        RES[f'E2_{name}_{rule}'] = r
        if rule == 'sym':
            mu, phi = np.linalg.eigh(-M); Lam = phi.T @ (q[:, None] * phi); Rm = B_ * Lam / np.sqrt(np.outer(G_ + mu, G_ + mu))
            dm = np.max(np.diag(Rm)); P(f'     modal diag max {dm:.4f} at n={np.argmax(np.diag(Rm))} mu1={mu[1]:.5f}; Galerkin N=1,2,3,5,10: ' + ', '.join(f'{np.linalg.eigvalsh(Rm[:N,:N])[-1]:.4f}' for N in (1, 2, 3, 5, 10)))
            for zc in 'MWK':
                Z = np.nonzero(zone == zc)[0]; muZ = -sab(M[np.ix_(Z, Z)]); P(f'     zone {zc}: mu_Z={muZ:.5f} (1/mu = {1/muZ:.1f} d) bound {B_*q[Z].min()/(G_+muZ):.4f}')
            dq = q - q.mean(); C = dq @ np.linalg.pinv(-M) @ dq / n
            Zm = np.nonzero(q >= q.max() - 1e-12)[0]; muZ0 = -sab(M[np.ix_(Zm, Zm)])
            P(f'     large-D: {B_*q.mean()/G_:.4f} + {B_*C/q.mean():.4f}/D0 ; small-D: {B_*q.max()/G_:.4f}(1 - {muZ0/G_:.4f} D0)')
            e = zl * w / (zl @ w); P('     elasticity share: ' + ', '.join(f'{zc}: {e[zone==zc].sum():.4f}' for zc in 'WCMK'))
            ev = np.sort(np.linalg.eigvalsh(M + np.diag(B_ * q - G_))); P(f'     spectral gap of J: {ev[-1]-ev[-2]:.4f}')
            M100, _, _ = build(z, lay.D_OFFICE, lay.Q_OFFICE, 'sym', D0=100.0); r100 = R0(M100, q); P(f'     D0=100: R0={r100:.4f} (+{100*(r100/(B_*q.mean()/G_)-1):.3f}%)')
            for D0 in (2.43e3, 4.87e4):
                Mx, _, _ = build(z, lay.D_OFFICE, lay.Q_OFFICE, 'sym', D0=D0); rx = R0(Mx, q); P(f'     D0={D0:g}: R0={rx:.5f} (+{100*(rx/(B_*q.mean()/G_)-1):.4f}%)')

P('== E3 partitions (O1, sym) ==')
z = lay.office_O1(); M, q, zone, edges, idx, mask, D = build(z, lay.D_OFFICE, lay.Q_OFFICE, 'sym', return_all=True); base = R0(M, q)
ww = np.array([(zone[i] == 'W') and (zone[j] == 'W') for i, j in edges]); P(f'  base {base:.4f}; desk-desk edges {ww.sum()}')
rng = np.random.default_rng(7)
def walls(p):
    closed = np.zeros(len(edges), bool); cand = rng.permutation(np.nonzero(ww)[0]); target = int(round(p * len(cand))); k = 0
    Mcur = M.copy()
    for e in cand:
        if k >= target: break
        closed[e] = True; Mt, _, _ = build(z, lay.D_OFFICE, lay.Q_OFFICE, 'sym', closed=closed)
        if connected(Mt): k += 1
        else: closed[e] = False
    Mt, _, _ = build(z, lay.D_OFFICE, lay.Q_OFFICE, 'sym', closed=closed); return R0(Mt, q), k
for p in (0.1, 0.3, 0.5, 1.0):
    vals = []; ks = []
    for _ in range(4 if p < 1 else 2):
        v, k = walls(p); vals.append(v); ks.append(k)
    P(f'  walls p={p}: closed {np.mean(ks):.0f}; R0 mean {np.mean(vals):.4f} min {np.min(vals):.4f} -> {100*(np.mean(vals)/base-1):+.2f}%')
    RES[f'E3_walls_{p}'] = float(np.mean(vals))
vals = []; Wc = np.argwhere(z == 'W')
while len(vals) < 12:
    z2 = z.copy(); pick = Wc[rng.choice(len(Wc), 39, replace=False)]; z2[pick[:, 0], pick[:, 1]] = 'B'
    M2, q2, _ = build(z2, lay.D_OFFICE, lay.Q_OFFICE, 'sym')
    if connected(M2): vals.append(R0(M2, q2))
P(f'  39 desks->obstacles: R0 {np.mean(vals):.4f} +- {np.std(vals):.4f} [{np.min(vals):.4f},{np.max(vals):.4f}] -> {100*(np.mean(vals)/base-1):+.2f}%')
q2 = q.copy(); q2[zone == 'W'] *= 0.5; P(f'  q_W->0.5q_W at D0=1: {R0(M, q2):.4f}')
M100, _, _ = build(z, lay.D_OFFICE, lay.Q_OFFICE, 'sym', D0=100.0); P(f'  same at D0=100: {R0(M100, q2):.4f} from {R0(M100, q):.4f}')

P('== E4 corridor D ==')
for name, z, Dm, qm in (('O1', lay.office_O1(), lay.D_OFFICE, lay.Q_OFFICE), ('aisles', lay.two_zone('aisles'), lay.D_TWO, lay.Q_TWO), ('side', lay.two_zone('side'), lay.D_TWO, lay.Q_TWO)):
    for rule in ('sym', 'dep'):
        out = []
        for c in (1, 0.1, 0.01, 10):
            Dm2 = dict(Dm); Dm2['C'] = Dm['C'] * c; M, q, _ = build(z, Dm2, qm, rule); out.append(R0(M, q))
        P(f'  {name} {rule}: x1 {out[0]:.4f} x0.1 {out[1]:.4f} x0.01 {out[2]:.4f} x10 {out[3]:.4f}')
        RES[f'E4_{name}_{rule}'] = out

P('== E5 targeting (O1 sym) ==')
z = lay.office_O1()
for D0 in (1.0, 100.0):
    M, q, zone = build(z, lay.D_OFFICE, lay.Q_OFFICE, 'sym', D0=D0); n = len(q)
    r, w, zl = R0(M, q, vec=True); e = zl * w; order_e = np.argsort(-e)
    rr = np.random.default_rng(11)
    for f in (0.05, 0.1, 0.2, 0.5):
        k = int(round(f * n)); q2 = q.copy(); q2[order_e[:k]] *= 0.5; one = R0(M, q2)
        rv = []
        for _ in range(15):
            q3 = q.copy(); q3[rr.choice(n, k, replace=False)] *= 0.5; rv.append(R0(M, q3))
        # q-ranked with random tie-breaks
        oq = np.lexsort((rr.random(n), -q)); q4 = q.copy(); q4[oq[:k]] *= 0.5
        P(f'  D0={D0:g} f={f}: one-shot elasticity {one:.4f}; q-ranked {R0(M, q4):.4f}; random {np.mean(rv):.4f} +- {np.std(rv):.4f}')
    # adaptive greedy 12 cells per step
    q2 = q.copy(); treated = np.zeros(n, bool); kd = 0; step = int(round(0.05 * n)); out = {}
    while kd < n * 0.55:
        r2, w2, z2 = R0(M, q2, vec=True); ec = z2 * w2; ec[treated] = -1; pick = np.argsort(-ec)[:step]; q2[pick] *= 0.5; treated[pick] = True; kd += len(pick); out[kd] = R0(M, q2)
    P(f'  D0={D0:g} adaptive greedy: ' + ', '.join(f'{k}/{n}={100*k/n:.1f}%: {v:.4f}' for k, v in out.items()))
    # FAIR comparison: adaptive greedy with exactly k cells
    for f in (0.1, 0.2, 0.5):
        k = int(round(f * n)); q2 = q.copy(); treated = np.zeros(n, bool); kd = 0
        while kd < k:
            r2, w2, z2 = R0(M, q2, vec=True); ec = z2 * w2; ec[treated] = -1; pick = np.argsort(-ec)[:min(step, k - kd)]; q2[pick] *= 0.5; treated[pick] = True; kd += len(pick)
        P(f'     adaptive greedy with exactly {k} cells ({100*k/n:.1f}%): {R0(M, q2):.4f}')

P('== E6 room size (D=0.5, uniform q) ==')
def square(L, D=0.5):
    n = L * L; I = sp.identity(L); T = sp.diags([np.ones(L - 1), np.ones(L - 1)], [1, -1]); A = sp.kron(I, T) + sp.kron(T, I)
    deg = np.asarray(A.sum(axis=0)).ravel(); return (D * (A - sp.diags(deg))).tocsc()
def mu_of(Mk, n):
    if n <= 400: return -sab(Mk.toarray())
    val = spla.eigsh(-Mk, k=1, sigma=-1e-9, which='LM', return_eigenvectors=False); return float(val[0])
for L in (3, 13, 20, 40):
    M = square(L); n = L * L
    kap = np.zeros(n); kap[L // 2] = 0.5                # door cell in middle of first row (wall)
    mu1 = mu_of((M - sp.diags(kap)).tocsc(), n)
    kap2 = np.zeros(n); kap2[:L] = 0.5                  # one side fully open
    mu2 = mu_of((M - sp.diags(kap2)).tocsc(), n)
    Md = M.toarray() if n <= 1700 else None
    refl = R0(Md, np.ones(n)) / (B_ / G_) if Md is not None else float('nan')
    P(f'  L={L}: reflecting {refl:.6f}; one door mu={mu1:.5f} ratio {G_/(G_+mu1):.4f}; open side mu={mu2:.5f} ratio {G_/(G_+mu2):.4f}')
for L in (8, 16, 32, 64, 128):
    M = square(L, 1.0); n = L * L; door = L // 2
    keep = np.setdiff1d(np.arange(n), [door]); Mabs = M[keep][:, keep]
    mu_abs = float(spla.eigsh(-Mabs, k=1, sigma=-1e-12, which='LM', return_eigenvectors=False)[0])
    c = L / np.exp(np.pi / (L * L * mu_abs))
    kap = np.zeros(n); kap[door] = 1.0; mu_k = float(spla.eigsh(-(M - sp.diags(kap)).tocsc(), k=1, sigma=-1e-12, which='LM', return_eigenvectors=False)[0])
    ser = 1 / (1 / mu_abs + n / 1.0)
    P(f'  one door L={L}: L^2 mu_abs={L*L*mu_abs:.4f}, c={c:.3f}; mu_kappa={mu_k:.4e}; series {ser:.4e} ({100*(ser/mu_k-1):+.2f}%)')
for nm, T, qq in (('office', 480, 1.2), ('supermarket', 27, 1.2), ('classroom', 60, 1.5), ('metro', 20, 2.5)):
    P(f'  dwell {nm}: {B_*qq/(G_+1440/T):.4f}; closed {B_*qq/G_:.3f}')
from scipy.optimize import brentq
P(f'  SIR final size at R0=6: {brentq(lambda a: 1-a-np.exp(-6*a), 0.5, 1):.4f}')
P(f'runtime {time.time()-t0:.1f}s')
json.dump(RES, open('v2_examples.json', 'w'), indent=1, default=float)
