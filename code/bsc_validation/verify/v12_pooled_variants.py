import numpy as np, json
import vevents as E
d = json.load(open('../data/02_predictive.json'))['events']
def pooled(evs, pmM=None):
    pm = [np.clip(np.array(d[e]['M2']['pmf'] if (pmM is None or e not in pmM) else pmM[e]),1e-300,None) for e in evs]; p0 = [np.clip(np.array(d[e]['M0']['pmf']),1e-300,None) for e in evs]
    g = np.meshgrid(*[np.arange(len(p)) for p in p0], indexing='ij')
    dd = sum(np.log(pm[a][g[a]])-np.log(p0[a][g[a]]) for a in range(len(evs))); w = np.prod([p0[a][g[a]] for a in range(len(evs))], axis=0)
    dobs = sum(np.log(pm[a][E.OBS[e][0]])-np.log(p0[a][E.OBS[e][0]]) for a,e in enumerate(evs))
    return round(float(dobs),2), float(w[dd>=dobs-1e-9].sum()/w.sum())
print('first-analysis pmfs, all four:', pooled(['T1','T2','T5','C1']))
print('without C1:', pooled(['T1','T2','T5']))
print('without T5:', pooled(['T1','T2','C1']))
print('without T5 and C1:', pooled(['T1','T2']))
print('T5 only:', pooled(['T5']), ' C1 only:', pooled(['C1']))
for e in ('T1','T2','T5','C1'):
    pm = np.array(d[e]['M2']['pmf']); k = E.OBS[e][0]; print(e, 'first-analysis pmf(obs)', pm[k], 'M0', np.array(d[e]['M0']['pmf'])[k])
