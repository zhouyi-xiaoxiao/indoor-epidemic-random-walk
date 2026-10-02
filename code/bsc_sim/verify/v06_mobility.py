"""Office: pair-level NGM K1 by a separately written (forward) formulation, rho(K1), counted cases, gen ratio, hotspots."""
import json, os, time, sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, vlib as V
B, G = V.BETA, V.GAMMA

def pair_ngm_forward(room, P=None, N=None, rho_ref=None, q=None):
    M = room.M; N = room.N if N is None else N; q = room.q if q is None else q
    Lm = sp.csr_matrix(V.gen_matrix(room))
    Pm = np.eye(M) if P is None else P
    rb = (N - 1) / M if rho_ref is None else rho_ref
    h = (B * q / rb)[:, None] * Pm
    I = sp.identity(M, format='csr')
    A = (G * sp.identity(M * M) + sp.diags(h.ravel()) - sp.kron(Lm, I) - sp.kron(I, Lm)).tocsr()
    d = A.diagonal()
    pre = spla.LinearOperator(A.shape, matvec=lambda v: v / d)
    K1 = np.zeros((M, M))
    for x in range(M):
        src = np.zeros((M, M)); src[x, :] = 1.0          # index at x, other anywhere (weight 1 each)
        phi, info = spla.cg(A, src.ravel(), rtol=1e-9, maxiter=200000, M=pre)
        assert info == 0
        phi = phi.reshape(M, M)                           # occupation density of surviving pair
        K1[:, x] = (N - 1) / M * (phi * h).sum(axis=0)   # infections landing at cell z = other's cell
    return K1

def perron(K):
    w, Vv = np.linalg.eig(K); k = int(np.argmax(w.real)); v = np.abs(Vv[:, k].real)
    return float(w[k].real), v / v.sum()

if __name__ == '__main__':
    out = json.load(open('v06_mobility.json')) if os.path.exists('v06_mobility.json') else {}
    for D0 in (1.0, 10.0, 0.03, 100.0):
        key = f'D0={D0}'
        if key in out: continue
        t0 = time.time()
        r = V.load('office', D0)
        K1 = pair_ngm_forward(r)
        np.save(f'v06_K1_office_D{D0}.npy', K1)
        R1, v1 = perron(K1)
        R1u = K1.sum(0).mean()
        K = V.ngm(r)
        tk = time.time() - t0
        n = 40000 if D0 <= 10 else 8000
        rng = np.random.default_rng(int(D0 * 100) + 5)
        ip = rng.choice(r.M, size=n, p=v1)
        s = V.sim(r, n, 5550001 + int(D0 * 100), gmax=0, index_pos=ip)
        m, se = V.mci(s['off'])
        su = V.sim(r, n, 5560001 + int(D0 * 100), gmax=0)
        mu, seu = V.mci(su['off'])
        s2 = V.sim(r, n, 5570001 + int(D0 * 100), gmax=1, index_pos=ip)
        g = s2['gens']; ratio = g[:, 2].sum() / g[:, 1].sum()
        bs = [g[ii, 2].sum() / g[ii, 1].sum() for ii in (rng.integers(0, n, n) for _ in range(200))]
        out[key] = dict(R0=V.rho(K), R1=R1, R1_uniform=R1u, sim_perron=m, sim_perron_se=se, sim_uniform=mu, sim_uniform_se=seu,
                        gen2_gen1=float(ratio), gen2_gen1_sd=float(np.std(bs)), n=n, k1_secs=tk, secs=time.time() - t0)
        print(key, {a: round(b, 4) for a, b in out[key].items()}, flush=True)
        json.dump(out, open('v06_mobility.json', 'w'), indent=1, default=float)
