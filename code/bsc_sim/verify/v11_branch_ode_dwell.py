import json, time, numpy as np, vlib as V
from scipy.integrate import solve_ivp
B, G = V.BETA, V.GAMMA
out = {}
def branching(r, P=None, gamma=G):
    Lm = V.gen_matrix(r); M = r.M; s = np.zeros(M)
    for it in range(100000):
        Ps = s if P is None else P @ s
        new = gamma * np.linalg.solve(gamma * np.eye(M) + np.diag(B * r.q * (1 - Ps)) - Lm.T, np.ones(M))
        if np.abs(new - s).max() < 1e-13: break
        s = new
    return 1 - s.mean()
def ode(r, P=None, T=400.0):
    Lm = V.gen_matrix(r); M = r.M; N = r.N; rb = (N - 1) / M
    y0 = np.concatenate([np.full(M, (N - 1) / M), np.full(M, 1 / M)])
    PT = None if P is None else P.T
    def f(t, y):
        S, I = y[:M], y[M:]
        e = B * r.q * I
        if PT is not None: e = PT @ e
        inf = S * e / rb
        return np.concatenate([Lm @ S - inf, Lm @ I + inf - G * I])
    t = np.arange(0, T + 0.5, 0.5)
    sol = solve_ivp(f, (0, T), y0, t_eval=t, method='LSODA', rtol=1e-8, atol=1e-10)
    I = sol.y[M:].sum(0); S = sol.y[:M].sum(0)
    return 1 - S[-1] / N, I.max() / N, t[I.argmax()]
for name in ('office', 'supermarket', 'classroom', 'metro'):
    for D0 in (1.0, 10.0):
        r = V.load(name, D0); rad = 1 / r.a_m
        P = V.kernel(r, rad) if rad >= 1 else None
        pb = branching(r, P)
        o = ode(r, P)
        out[f'{name}|{D0}'] = dict(branch=pb, ode_attack=o[0], ode_peak=o[1], ode_day=o[2])
        print(name, D0, 'branching P(major) %.3f ; ODE attack %.3f peak %.3f day %.1f' % (pb, *o), flush=True)
# ---- dwell time
r = V.load('office', 1.0); w = 1 / 3
p = V.pair_solve(r, gamma=G / w); R1dc = 99 * p.mean()
R0dc = V.rho(V.ngm(r, gamma=G / w))
s = V.sim(r, 40000, 90210, gmax=0, sched=[(w, 1.0, 1.0), (1.0, 0.0, 0.0)], t_max=600.0)
m, se = V.mci(s['off'])
print('office 8h/day D0=1: R0 duty %.3f pair R1 duty %.3f counted with explicit schedule %.4f +- %.4f' % (R0dc, R1dc, m, se))
out['office_duty_D1'] = dict(R0=R0dc, R1=R1dc, sim=m, se=se)
f = V.sim(r, 4000, 90211, sched=[(w, 1.0, 1.0), (1.0, 0.0, 0.0)], t_max=3000.0)
ar = f['final'] / 100
print('   full epidemics with schedule: P(attack>=10%%) = %.3f, unfinished %d' % ((ar >= 0.1).mean(), f['unfinished']))
out['office_duty_D1_pmajor'] = float((ar >= 0.1).mean())
# per visit
for name, hours in (('office', 8.0), ('supermarket', 27 / 60), ('classroom', 1.0), ('metro', 20 / 60)):
    r = V.load(name, 1.0); rad = 1 / r.a_m
    P = V.kernel(r, rad) if rad >= 1 else None
    T = hours / 24
    s = V.sim(r, 400000, 5150 + int(hours * 100), P=P, gmax=0, t_max=T)
    m, se = V.mci(s['off'])
    mf = B * r.q.mean() * (1 - np.exp(-G * T)) / G
    print(name, 'per visit (%.2f h): mean-field %.4f counted %.4f +- %.4f' % (hours, mf, m, se), flush=True)
    out[f'visit_{name}'] = dict(mf=mf, sim=m, se=se)
json.dump(out, open('v11.json', 'w'), indent=1, default=float)
