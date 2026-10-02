"""Occupancy + interventions (office, D0=1) with the second simulator and the pair solver of the re-check."""
import json, os, time, copy, numpy as np, vlib as V
B, G = V.BETA, V.GAMMA
out = json.load(open('v10.json')) if os.path.exists('v10.json') else {}
def run(key, r, N, q=None, beta=B, rho_ref=None, n=4000, P=None, scale=1.0):
    if key in out: return
    t0 = time.time()
    q = r.q if q is None else q
    R0 = V.rho(V.ngm(r, P, beta=beta, q=q)) * scale
    p = V.pair_solve(r, P, beta=beta, N=N, rho_ref=rho_ref, q=q)
    R1 = (N - 1) * p.mean()
    s0 = V.sim(r, n, 1234567 + len(out), P=P, gmax=0, beta=beta, N=N, rho_ref=rho_ref, q=q)
    m, se = V.mci(s0['off'])
    s = V.sim(r, n, 2234567 + len(out), P=P, beta=beta, N=N, rho_ref=rho_ref, q=q, t_max=800.0, dt_curve=0.5, n_curve=1601)
    ar = s['final'] / N; maj = ar >= 0.1; k = int(maj.sum())
    pm, lo, hi = V.wilson(k, n)
    row = dict(N=N, cells=r.M, R0=R0, R1=R1, counted=m, counted_se=se, p_major=pm, p_lo=lo, p_hi=hi, attack_all=float(ar.mean()))
    if k > 1:
        cv = s['curve'][maj]
        row.update(attack_major=float(ar[maj].mean()), peak_major=float((cv.max(1) / N).mean()), peak_day_major=float((cv.argmax(1) * 0.5).mean()),
                   peak_day_se=float((cv.argmax(1) * 0.5).std(ddof=1) / np.sqrt(k)))
    row['secs'] = time.time() - t0
    out[key] = row
    print(key, {a: round(float(b), 4) for a, b in row.items()}, flush=True)
    json.dump(out, open('v10.json', 'w'), indent=1, default=float)
base = V.load('office', 1.0)
rho100 = 99 / base.M
for N in (50, 100, 400):
    r = V.load('office', 1.0, N=N)
    run(f'occ|FD|N={N}', r, N)
    if N != 100:
        run(f'occ|DD|N={N}', r, N, rho_ref=rho100, scale=(N - 1) / 99)
run('vent q*0.5', base, 100, q=0.5 * base.q)
run('masks beta*0.3', base, 100, beta=0.3 * B)
# partitions: 39 random W cells closed, floor kept connected (layout drawn for the re-check)
def connected(acc):
    ys, xs = np.nonzero(acc); seen = np.zeros(acc.shape, bool); st = [(ys[0], xs[0])]; seen[ys[0], xs[0]] = True; c = 0
    while st:
        y, x = st.pop(); c += 1
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            yy, xx = y + dy, x + dx
            if 0 <= yy < acc.shape[0] and 0 <= xx < acc.shape[1] and acc[yy, xx] and not seen[yy, xx]:
                seen[yy, xx] = True; st.append((yy, xx))
    return c == acc.sum()
for lay in range(3):
    rng = np.random.default_rng(77 + lay)
    acc = base.access.copy()
    cand = [c for k, c in enumerate(base.cells) if base.zone[k] == 'W']
    added = 0
    for j in rng.permutation(len(cand)):
        if added == 39: break
        y, x = cand[j]; acc[y, x] = False
        if connected(acc): added += 1
        else: acc[y, x] = True
    Dg = np.zeros(acc.shape); qg = np.zeros(acc.shape)
    for k, (y, x) in enumerate(base.cells): Dg[y, x] = base.D[k]; qg[y, x] = base.q[k]
    part = V.Room(acc, Dg, qg, 100, base.a_m)
    run(f'partitions lay{lay}|FD', part, 100)
    run(f'combination lay{lay}|FD', part, 50, q=0.5 * part.q)
    run(f'combination lay{lay}|DD', part, 50, q=0.5 * part.q, rho_ref=rho100, scale=(49 / part.M) / rho100)
