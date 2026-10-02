"""Barrier experiment, layouts and simulator of the re-check. Identical starts across arms + a stationary-start control."""
import json, time, numpy as np, vlib as V
B, G = V.BETA, V.GAMMA
NX, NY, N, Q, D = 20, 13, 100, 1.2, 1.0
FR = [0.0, 0.10, 0.30]
def largest(acc):
    ny, nx = acc.shape; lab = -np.ones(acc.shape, int); best = None; k = 0
    for y0 in range(ny):
        for x0 in range(nx):
            if acc[y0, x0] and lab[y0, x0] < 0:
                st = [(y0, x0)]; lab[y0, x0] = k; cells = []
                while st:
                    y, x = st.pop(); cells.append((y, x))
                    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        yy, xx = y + dy, x + dx
                        if 0 <= yy < ny and 0 <= xx < nx and acc[yy, xx] and lab[yy, xx] < 0:
                            lab[yy, xx] = k; st.append((yy, xx))
                if best is None or len(cells) > len(best): best = cells
                k += 1
    m = np.zeros(acc.shape, bool)
    for y, x in best: m[y, x] = True
    return m
out = {}
agg = {f: dict(ar=[], maj=[], pk=[], pkt=[], ar_stat=[], R1=[], Rsim=[], cells=[], R0=[]) for f in FR}
t0 = time.time()
lay = 0; seed = 100
while lay < 8:
    seed += 1
    rng = np.random.default_rng(seed)
    perm = rng.permutation([i for i in range(NX * NY) if i != 6 * NX + 10])
    accs = {}
    for f in FR:
        b = np.zeros(NX * NY, bool); b[perm[:int(round(f * NX * NY))]] = True
        accs[f] = largest(~b.reshape(NY, NX))
    if not accs[0.30][6, 10]: continue
    if any((accs[0.30] & ~accs[f]).any() for f in FR): continue
    ys, xs = np.nonzero(accs[0.30])
    nrep = 1000
    pick = rng.integers(0, len(xs), size=(nrep, N))
    sy, sx = ys[pick], xs[pick]; sy[:, 0] = 6; sx[:, 0] = 10
    for f in FR:
        r = V.uniform(NX, NY, N, D=D, q=Q, barriers=~accs[f])
        init = r.idx[sy, sx]
        assert (init >= 0).all()
        s = V.sim(r, nrep, 9000 + 17 * lay + int(100 * f), init_pos=init, t_max=600.0, dt_curve=0.5, n_curve=1201)
        ar = s['final'] / N; maj = ar >= 0.1
        cv = s['curve'][maj]
        a = agg[f]
        a['ar'].append(ar); a['maj'].append(maj); a['pk'].append(cv.max(1) / N); a['pkt'].append(cv.argmax(1) * 0.5)
        # stationary start control: everyone uniform on own walkable set, index at centre
        c = int(r.idx[6, 10])
        s2 = V.sim(r, nrep, 9500 + 17 * lay + int(100 * f), index_pos=np.full(nrep, c), t_max=600.0)
        a['ar_stat'].append(s2['final'] / N)
        p = V.pair_solve(r)
        a['R1'].append((N - 1) * p.mean()); a['cells'].append(r.M); a['R0'].append(V.rho(V.ngm(r)))
        s3 = V.sim(r, 2000, 9700 + 17 * lay + int(100 * f), gmax=0)
        a['Rsim'].append(s3['off'].mean())
    lay += 1
    print('layout', lay, 'done', time.time() - t0, flush=True)
for f in FR:
    a = agg[f]
    AR = np.array(a['ar']); MJ = np.concatenate(a['maj'])
    out[f'{f}'] = dict(cells=float(np.mean(a['cells'])), R0=float(np.mean(a['R0'])), R1=float(np.mean(a['R1'])), Rsim=float(np.mean(a['Rsim'])),
                       attack=float(AR.mean()), attack_se_layout=float(AR.mean(1).std(ddof=1) / np.sqrt(len(AR))),
                       p_major=float(MJ.mean()), attack_stationary=float(np.mean(a['ar_stat'])),
                       peak_major=float(np.concatenate(a['pk']).mean()), peak_day_major=float(np.concatenate(a['pkt']).mean()),
                       peak_day_se=float(np.concatenate(a['pkt']).std(ddof=1) / np.sqrt(MJ.sum())))
    print(f, out[f'{f}'])
d = np.array(agg[0.10]['ar']) - np.array(agg[0.0]['ar'])
out['paired_10_minus_0'] = dict(mean=float(d.mean()), se_runs=float(d.std(ddof=1) / np.sqrt(d.size)), se_layouts=float(d.mean(1).std(ddof=1) / np.sqrt(d.shape[0])))
print(out['paired_10_minus_0'])
json.dump(out, open('v08_obstacles.json', 'w'), indent=1)
