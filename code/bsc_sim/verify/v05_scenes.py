"""Scene epidemics with the second simulator: P(major), attack, peak, peak day, duration, counted R."""
import json, time, sys, numpy as np, vlib as V
B, G = V.BETA, V.GAMMA
out = json.load(open('v05_scenes.json')) if __import__('os').path.exists('v05_scenes.json') else {}
jobs = [('office', 1.0, '1m', 4000), ('supermarket', 1.0, '1m', 4000), ('classroom', 1.0, '1m', 4000),
        ('metro', 1.0, '1m', 3000), ('metro', 1.0, 'same', 3000), ('classroom', 1.0, 'same', 4000),
        ('office', 10.0, '1m', 2000), ('supermarket', 10.0, '1m', 1500), ('classroom', 10.0, '1m', 2000)]
for name, D0, kn, n in jobs:
    key = f'{name}|{D0}|{kn}'
    if key in out: continue
    t0 = time.time()
    r = V.load(name, D0)
    rad = 1.0 / r.a_m
    P = V.kernel(r, rad) if (kn == '1m' and rad >= 1) else None
    s = V.sim(r, n, 990001 + len(out), P=P, t_max=400.0, dt_curve=0.5, n_curve=801)
    ar = s['final'] / r.N
    maj = ar >= 0.10
    k = int(maj.sum())
    p, lo, hi = V.wilson(k, n)
    row = dict(n=n, p_major=p, p_lo=lo, p_hi=hi, attack_all=float(ar.mean()), unfinished=s['unfinished'])
    if k > 1:
        cv = s['curve'][maj]
        gp = cv.max(axis=1) / r.N
        gt = cv.argmax(axis=1) * 0.5
        for nm, arr in (('attack_major', ar[maj]), ('peak_exact_major', s['peak'][maj] / r.N),
                        ('peak_grid_major', gp), ('peak_day_grid_major', gt), ('peak_day_exact_major', s['peak_t'][maj]),
                        ('duration_major', s['dur'][maj])):
            m, se = V.mci(arr)
            row[nm] = m; row[nm + '_se'] = se
    g = s['gens']
    row['gen2_over_gen1'] = float(g[:, 2].sum() / max(g[:, 1].sum(), 1))
    s0 = V.sim(r, n, 880001 + len(out), P=P, gmax=0, t_max=400.0)
    m, se = V.mci(s0['off'])
    row['R_alone'] = m; row['R_alone_se'] = se
    row['secs'] = time.time() - t0
    out[key] = row
    print(key, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in row.items()}, flush=True)
    json.dump(out, open('v05_scenes.json', 'w'), indent=1, default=float)
