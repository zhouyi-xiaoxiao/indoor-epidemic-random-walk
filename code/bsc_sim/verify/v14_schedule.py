import sys, json, numpy as np, vlib as V
out = {}
r = V.load('office', 1.0)
def summ(s, N):
    ar = s['final'] / N; maj = ar >= 0.1; cv = s['curve'][maj]
    return dict(p=float(maj.mean()), attack_major=float(ar[maj].mean()), peak=float((cv.max(1) / N).mean()), day=float((cv.argmax(1) * 0.5).mean()),
                day_se=float((cv.argmax(1) * 0.5).std(ddof=1) / np.sqrt(maj.sum())), peak_se=float((cv.max(1) / N).std(ddof=1) / np.sqrt(maj.sum())))
a = V.sim(r, 6000, 31, sched=[(0.5, 1.8, 1.3), (1.0, 0.2, 0.7)], t_max=600.0, dt_curve=0.5, n_curve=1201)
b = V.sim(r, 6000, 32, sched=[(1.0, 1.0, 1.0)], t_max=600.0, dt_curve=0.5, n_curve=1201)
out['check_varying'] = summ(a, 100); out['check_mean'] = summ(b, 100)
print('second simulator, office: varying', out['check_varying']); print('second simulator, office: mean   ', out['check_mean'])
# their metro arms with a new seed
sys.path.insert(0, '../src'); sys.path.insert(0, '../experiments')
from bsc_sim import theory as th
from bsc_sim.scenes import load_scene
from bsc_sim.sim import simulate, summarize
import importlib
m9 = importlib.import_module('09_metro_time')
hours = np.arange(24) + 0.5; rho_h = m9.rho_of_hour(hours); m_h = np.where(rho_h >= 3.0, 1.0, 0.8)
for D0 in (1.0,):
    sc = load_scene('metro', D0=D0); P = th.contact_kernel(sc, 1.0 / sc.a_m)
    zones = sc.meta['zones']
    coef = np.array([zones[z]['D_coef'] for z in sc.site_zone]); expo = np.array([zones[z]['D_exp'] for z in sc.site_zone])
    D_h = [D0 * coef / rr ** expo for rr in rho_h]; D_mean = np.mean(D_h, axis=0)
    for arm, cfg in (('varying', dict(D=D_h[0], schedule=dict(period=1.0, segments=[((k + 1) / 24.0, D_h[k], float(m_h[k])) for k in range(24)]))),
                     ('mean', dict(D=D_mean, schedule=dict(period=1.0, segments=[(1.0, D_mean, float(m_h.mean()))])))):
        res = simulate(sc, 0.5, 0.14, n_rep=2000, seed=4321, P=P, D=cfg['D'], schedule=cfg['schedule'], t_max=300.0, dt_out=0.25)
        s = summarize(res)
        out[f'their_metro_{arm}'] = {k: s[k] for k in ('p_major', 'attack_major_mean', 'peak_prev_major_mean', 'peak_time_major_mean')}
        print('THEIR sim metro D0=1, new seed', arm, out[f'their_metro_{arm}'])
json.dump(out, open('v14_schedule.json', 'w'), indent=1, default=float)
