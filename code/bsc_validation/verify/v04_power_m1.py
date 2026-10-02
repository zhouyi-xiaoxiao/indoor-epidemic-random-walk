import numpy as np
from vlat import *
import vevents as E
r = np.load('res_holdout.npz')
evs = ['T1','T2','T5','C1']
pm = {m: [np.clip(r[f'{e}_{m}'],1e-300,None) for e in evs] for m in ('M0','M1','M2')}
sizes = [len(p) for p in pm['M0']]
g = np.meshgrid(*[np.arange(s) for s in sizes], indexing='ij')
J = {m: np.prod([pm[m][a][g[a]] for a in range(4)], axis=0) for m in pm}
for m in J: J[m] /= J[m].sum()
for cand in ('M1','M2'):
    d = sum(np.log(pm[cand][a][g[a]])-np.log(pm['M0'][a][g[a]]) for a in range(4))
    flat = d.ravel(); o = np.argsort(-flat); cum = np.cumsum(J['M0'].ravel()[o]); pv = np.empty_like(flat)
    # ties-inclusive p
    ds = flat[o]; last = np.searchsorted(-ds, -flat, side='right')-1; pv = cum[last].reshape(sizes)
    a1 = (pv < 0.05) & (d > 0)
    fails = np.zeros(sizes, int)
    for a in range(4):
        pe = np.array([two_sided(pm[cand][a]/pm[cand][a].sum(), k) for k in range(sizes[a])])
        fails += (pe[g[a]] < 0.05)
    for truth in ('M0','M1','M2'):
        w = J[truth]
        print(cand, 'truth', truth, 'P(A1)=%.4f P(noA2fail)=%.4f P(<=1 fail)=%.4f P(A1&nofail)=%.4f' % (w[a1].sum(), w[fails==0].sum(), w[fails<=1].sum(), w[a1&(fails==0)].sum()))
# what is the weakest monotone alternative that passes A1 given these observations? (any model with pmf = M2?) -> instead compute
# pooled p using only T5 (Fisher-type) to show where the significance comes from
import scipy.stats as st
for e in evs:
    k1,n1,k2,n2 = E.OBS[e]; print(e, 'one-sided Fisher (near>far) p =', st.fisher_exact([[k1,n1-k1],[k2,n2-k2]], alternative='greater')[1])
