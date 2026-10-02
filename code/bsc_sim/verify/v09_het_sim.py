import json, os, time, numpy as np, vlib as V
from v06_mobility import pair_ngm_forward, perron
from v09_het_det import twozone
out = json.load(open('v09_het.json')) if os.path.exists('v09_het.json') else {}
ct = 1 - 0.5 / 1.2
jobs = [(1.0, 100, 0.0, 'zoneD'), (1.0, 100, ct, 'zoneD'), (1.0, 100, ct, 'uniformD'), (1.0, 100, 0.0, 'uniformD'),
        (1.0, 100, -2.0, 'zoneD'), (10.0, 100, 0.0, 'zoneD'), (10.0, 100, ct, 'zoneD'),
        (1.0, 1000, 0.0, 'zoneD'), (1.0, 1000, ct, 'zoneD'), (1.0, 1000, 1.0, 'zoneD')]
for D0, N, c, var in jobs:
    key = f'{var}|D0={D0}|N={N}|c={c:.3f}'
    if key in out: continue
    t0 = time.time()
    r = twozone(c, D0, N, var)
    K1 = pair_ngm_forward(r)
    R1, v1 = perron(K1)
    n = 20000 if N == 100 else 3000
    if D0 > 1: n = 6000
    rng = np.random.default_rng(42)
    ip = rng.choice(r.M, size=n, p=v1)
    s = V.sim(r, n, 660001 + len(out), gmax=0, index_pos=ip)
    m, se = V.mci(s['off'])
    su = V.sim(r, n, 670001 + len(out), gmax=0)
    mu, seu = V.mci(su['off'])
    row = dict(R0=V.rho(V.ngm(r)), R1=R1, R1_uniform=float(K1.sum(0).mean()), sim_perron=m, se=se, sim_uniform=mu, se_u=seu, n=n)
    if N == 100:
        nf = 4000 if D0 == 1 else 1500
        f = V.sim(r, nf, 680001 + len(out), t_max=800.0)
        ar = f['final'] / N
        p, lo, hi = V.wilson(int((ar >= 0.1).sum()), nf)
        row.update(p_major=p, p_lo=lo, p_hi=hi, attack_all=float(ar.mean()), attack_se=float(ar.std(ddof=1) / np.sqrt(nf)))
    row['secs'] = time.time() - t0
    out[key] = row
    print(key, {a: round(float(b), 4) for a, b in row.items()}, flush=True)
    json.dump(out, open('v09_het.json', 'w'), indent=1, default=float)
