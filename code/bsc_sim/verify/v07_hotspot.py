import json, numpy as np, vlib as V, time
B, G = V.BETA, V.GAMMA
r = V.load('office', 1.0)
K = V.ngm(r); K1 = np.load('v06_K1_office_D1.0.npy')
u = np.ones(r.M) / r.M
out = {}
t0 = time.time()
s = V.sim(r, 40000, 7770001, gmax=2)
cg = s['cellgen']
for g in (1, 2, 3):
    obs = cg[g - 1] / cg[g - 1].sum()
    mf = np.linalg.matrix_power(K, g) @ u; mf /= mf.sum()
    pr = np.linalg.matrix_power(K1, g) @ u; pr /= pr.sum()
    out[f'gen{g}'] = dict(n=int(cg[g - 1].sum()), r_mf=float(np.corrcoef(obs, mf)[0, 1]), r_pair=float(np.corrcoef(obs, pr)[0, 1]),
                          per_run=float(cg[g - 1].sum() / 40000))
    print('gmax=2 gen', g, out[f'gen{g}'])
# first generation alone
s = V.sim(r, 100000, 7770002, gmax=0)
obs = s['cellgen'][0] / s['cellgen'][0].sum()
mf = K @ u; mf /= mf.sum(); pr = K1 @ u; pr /= pr.sum()
# ceiling: correlation expected if pr were exact, by multinomial resampling
rng = np.random.default_rng(1); n1 = int(s['cellgen'][0].sum())
ceil = np.mean([np.corrcoef(rng.multinomial(n1, pr) / n1, pr)[0, 1] for _ in range(200)])
out['gen1_alone'] = dict(n=n1, r_mf=float(np.corrcoef(obs, mf)[0, 1]), r_pair=float(np.corrcoef(obs, pr)[0, 1]), ceiling=float(ceil))
print('gen1 alone', out['gen1_alone'])
for z in 'WCMK':
    m = r.zone == z
    print(z, 'obs %.4f pair %.4f mf %.4f  area %.4f' % (obs[m].sum(), pr[m].sum(), mf[m].sum(), m.mean()))
    out['zone_' + z] = dict(obs=float(obs[m].sum()), pair=float(pr[m].sum()), mf=float(mf[m].sum()), area=float(m.mean()))
json.dump(out, open('v07_hotspot.json', 'w'), indent=1)
print(time.time() - t0)
