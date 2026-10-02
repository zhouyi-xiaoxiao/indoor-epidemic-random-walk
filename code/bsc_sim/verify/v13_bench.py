"""Benchmark, with separately written code, of 'spectral R0' vs 'best-case Monte Carlo for the same number', MC in C."""
import ctypes as C, time, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, vlib as V
B, G = V.BETA, V.GAMMA
lib = C.CDLL('./libvwalk.dylib')
ip = np.ctypeslib.ndpointer(np.int32, flags='C_CONTIGUOUS'); dp = np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS')
lib.walk.argtypes = [C.c_int, ip, dp, dp, dp, C.c_double, C.c_int, ip, C.c_uint64, dp, dp]; lib.walk.restype = None
out = {}
for D0 in (1.0, 10.0):
    r = V.load('office', D0); M = r.M
    nbr, w = V.neighbours(r); wt = w.sum(1)
    Ls = sp.csc_matrix(V.gen_matrix(r))
    def spectral():
        lu = spla.splu((G * sp.identity(M, format='csc') - Ls).tocsc())
        v = np.ones(M) / M; lam = 0
        for it in range(100000):
            x = B * r.q * lu.solve(v); new = x.sum(); x /= new
            if abs(new - lam) < 1e-10 * new: break
            lam, v = new, x
        return new, it
    ts = []
    for _ in range(9):
        t0 = time.process_time(); val, its = spectral(); ts.append(time.process_time() - t0)
    t_spec = float(np.median(ts))
    Rd = V.rho(V.ngm(r))
    ev = np.sort(np.abs(np.linalg.eigvals(V.ngm(r))))[::-1]; ratio = ev[1] / ev[0]
    burn = int(np.ceil(np.log(1e-3) / np.log(ratio)))
    rng = np.random.default_rng(3); n_w, Gn = 20000, 30
    nbf = np.ascontiguousarray(nbr.ravel()); wf = np.ascontiguousarray(w.ravel()); bq = np.ascontiguousarray(B * r.q)
    t0 = time.process_time(); dist = np.ones(M) / M; ks = []
    for g in range(burn + Gn):
        starts = rng.choice(M, size=n_w, p=dist).astype(np.int32)
        ws = np.zeros(n_w); cell = np.zeros(M)
        lib.walk(M, nbf, wf, wt, bq, G, n_w, starts, C.c_uint64(1000 + g), ws, cell)
        if g >= burn: ks.append(ws.mean())
        dist = cell / cell.sum()
    t_mc = time.process_time() - t0
    ks = np.array(ks); est = ks.mean(); rel = ks.std(ddof=1) / np.sqrt(len(ks)) / est
    # also the walker-level relative SE (CV/sqrt(n))
    rel_w = ws.std(ddof=1) / ws.mean() / np.sqrt(n_w * Gn)
    t1 = t_mc * (rel / 0.01) ** 2
    out[D0] = dict(R0=val, its=its, dense=Rd, t_spec=t_spec, ratio=float(ratio), burn=burn, mc=float(est), rel=float(rel), rel_walker=float(rel_w), t_mc=t_mc,
                   t_mc_1pct=float(t1), speedup_1pct=float(t1 / t_spec), speedup_01pct=float(t1 * 100 / t_spec),
                   t_mc_1pct_noburn=float(t_mc * Gn / (burn + Gn) * (rel / 0.01) ** 2))
    print(D0, out[D0], flush=True)
json.dump(out, open('v13_bench.json', 'w'), indent=1)
