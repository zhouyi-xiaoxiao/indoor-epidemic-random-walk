import json, time, numpy as np, vlib as V
out = {}
t0 = time.time()
r = V.load('office', 1.0)
s = V.sim(r, 20000, 4242, mode=1)
m, se = V.mci(s['off']); 
out['office'] = dict(off=m, se=se, mean_final=float(s['final'].mean()), max_final=int(s['final'].max()), p10=float((s['final'] >= 10).mean()))
print('office per-cell rule', out['office'], time.time() - t0)
r = V.uniform(20, 20, 100)
s = V.sim(r, 20000, 4243, mode=1)
m, se = V.mci(s['off'])
out['uniform20'] = dict(off=m, se=se, mean_final=float(s['final'].mean()), max_final=int(s['final'].max()))
print('uniform 20x20 per-cell rule', out['uniform20'], time.time() - t0)
for name in ('supermarket', 'classroom', 'metro'):
    r = V.load(name, 1.0)
    s = V.sim(r, 4000, 4244, mode=1)
    m, se = V.mci(s['off'])
    out[name] = dict(off=m, se=se, attack=float(s['final'].mean() / r.N), p10=float((s['final'] >= 0.1 * r.N).mean()))
    print(name, out[name], time.time() - t0)
json.dump(out, open('v03_percell_rule.json', 'w'), indent=1)
