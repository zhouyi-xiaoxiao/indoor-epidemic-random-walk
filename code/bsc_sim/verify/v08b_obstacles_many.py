"""40 barrier layouts x 1000 reps, identical starts, arms 0/0.10/0.30: layout-level uncertainty of the effects."""
import json, time, numpy as np, vlib as V
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

res = {f: dict(ar=[], pkt=[], pk=[], pm=[]) for f in FR}
dpk = []; dar = []; dpkt_pair = []
lay = 0; seed = 5000; t0 = time.time()
while lay < 40:
    seed += 1
    rng = np.random.default_rng(seed)
    perm = rng.permutation([i for i in range(NX * NY) if i != 6 * NX + 10])
    accs = {}
    for f in FR:
        b = np.zeros(NX * NY, bool); b[perm[:int(round(f * NX * NY))]] = True
        accs[f] = largest(~b.reshape(NY, NX))
    if not accs[0.30][6, 10]: continue
    if any((accs[0.30] & ~accs[f]).any() for f in FR): continue
    ys, xs = np.nonzero(accs[0.30]); nrep = 1000
    pick = rng.integers(0, len(xs), size=(nrep, N))
    sy, sx = ys[pick], xs[pick]; sy[:, 0] = 6; sx[:, 0] = 10
    tmp = {}
    for f in FR:
        r = V.uniform(NX, NY, N, D=D, q=Q, barriers=~accs[f])
        s = V.sim(r, nrep, 31000 + 13 * lay + int(100 * f), init_pos=r.idx[sy, sx], t_max=600.0, dt_curve=0.5, n_curve=1201)
        ar = s['final'] / N; maj = ar >= 0.1; cv = s['curve']
        tmp[f] = (ar, maj, cv.argmax(1) * 0.5, cv.max(1) / N)
        res[f]['ar'].append(ar.mean()); res[f]['pkt'].append(tmp[f][2][maj].mean()); res[f]['pk'].append(tmp[f][3][maj].mean()); res[f]['pm'].append(maj.mean())
    lay += 1
out = {}
for f in FR:
    out[str(f)] = {k: (float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(len(v)))) for k, v in res[f].items()}
    print(f, out[str(f)])
for f in (0.10, 0.30):
    for k in ('ar', 'pkt', 'pk', 'pm'):
        d = np.array(res[f][k]) - np.array(res[0.0][k])
        out[f'diff_{f}_{k}'] = (float(d.mean()), float(d.std(ddof=1) / np.sqrt(len(d))))
        print('diff', f, k, out[f'diff_{f}_{k}'])
print(time.time() - t0)
json.dump(out, open('v08b_obstacles_many.json', 'w'), indent=1)
