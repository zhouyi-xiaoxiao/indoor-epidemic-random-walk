import numpy as np, json, pandas as pd
from scipy import stats, optimize
from vlat import *
# ---- retyped by the verifier from Hu et al. printed percentages; counts from reconstruction file
df = pd.read_csv('../../bsc_outbreaks/verify/results/train_matrix_reconstruction.csv')
df['dr'] = df['row']; df['dc'] = np.where(df['row']==0, df['printed_col']+1, df['printed_col'])
print('sum k, n', df.k.sum(), df.n.sum(), ' max |k/n - pct|', np.abs(100*df.k/df.n - df.pct).max())
lat = Lat(6, 17, 0.5, 1.0)
cols = np.array([0,1,2,4,5]); R = np.repeat(np.arange(17), 5); Cc = np.tile(cols, 17)
sites = lat.site(Cc, R)
DR = np.abs(R[:,None]-R[None,:]); DC = np.abs(Cc[:,None]-Cc[None,:])
cellmask = [((DR==r)&(DC==c)) for r,c in zip(df.dr, df.dc)]
print('pairs per cell', [int(m.sum()) for m in cellmask])
from scipy.stats import gamma
sh = (2.1/1.8)**2; sc = 1.8**2/2.1
u = (np.arange(16)+0.5)/16; Tn = gamma.ppf(u, sh, scale=sc)
kq = np.exp(np.log(5)+ (np.arange(7)+0.5)/7*(np.log(20)-np.log(5))) + 0.93
k = df.k.to_numpy(); n = df.n.to_numpy()
use = ~((df.dr==0)&(df.dc==1)).to_numpy()
def cellX(model, D, kap):
    out = np.zeros((16, len(df)))  # mean of (1-exp(-eps X)) must be done pairwise -> keep pair values
    Xs = [expo_mat(lat, model, D, kap, [T], sites) for T in Tn]
    return Xs
def prob(Xs, eps):
    p = np.zeros(len(df))
    for X in Xs:
        q = -np.expm1(-eps*X)
        p += np.array([q[m].mean() for m in cellmask])/16
    return p
def profile(Xs, use):
    def nll(le):
        p = np.clip(prob(Xs, np.exp(le)), 1e-300, 1-1e-12)
        return -stats.binom.logpmf(k[use], n[use], p[use]).sum()
    les = np.linspace(-12, 25, 38); v = [nll(l) for l in les]; i = int(np.argmin(v))
    r = optimize.minimize_scalar(nll, bounds=(les[max(i-1,0)], les[min(i+1,37)]), method='bounded', options={'xatol':1e-7})
    return -r.fun, np.exp(r.x)
res = {}
X0 = cellX('M0', 1.0, kq[3]); ll0, e0 = profile(X0, use); print('M0', ll0, e0, prob(X0,e0)[:3])
logD = np.linspace(-2, 4, 61)
ll1 = []; e1=[]
for D in 10**logD:
    l,e = profile(cellX('M1', D, 0.0), use); ll1.append(l); e1.append(e)
ll1 = np.array(ll1); print('M1 max', ll1.max(), 10**logD[ll1.argmax()])
ll2 = np.zeros((61,7)); e2 = np.zeros((61,7))
for i,D in enumerate(10**logD):
    for j,kap in enumerate(kq):
        ll2[i,j], e2[i,j] = profile(cellX('M2', D, kap), use)
i,j = np.unravel_index(ll2.argmax(), ll2.shape); print('M2 max', ll2.max(), 10**logD[i], kq[j])
sat = stats.binom.logpmf(k[use], n[use], k[use]/n[use]).sum()
print('deviance M0,M1,M2', 2*(sat-ll0), 2*(sat-ll1.max()), 2*(sat-ll2.max()), 'df', use.sum()-1, use.sum()-2)
for d,dfree in ((2*(sat-ll0), use.sum()-1), (2*(sat-ll1.max()), use.sum()-2), (2*(sat-ll2.max()), use.sum()-2)):
    print('  p', stats.chi2.sf(d, dfree))
# posterior summaries
def q(w, a):
    c = np.cumsum(w/w.sum()); return 10**np.interp(a, c, logD)
w1 = np.exp(ll1-ll1.max()); w2 = np.exp(ll2-ll2.max()).sum(axis=1)
print('M1 D 5/50/95', q(w1,.05), q(w1,.5), q(w1,.95)); print('M2 D 5/50/95', q(w2,.05), q(w2,.5), q(w2,.95))
# fitted cell probs at mode
p1 = prob(cellX('M1', 10**logD[ll1.argmax()], 0.0), e1[int(ll1.argmax())])
p2 = prob(cellX('M2', 10**logD[i], kq[j]), e2[i,j])
tab = df[['dr','dc','k','n']].copy(); tab['obs%']=100*k/n; tab['M1%']=100*p1; tab['M2%']=100*p2
print(tab.round(3).to_string())
# variants
for name, um in (('S-adj', np.ones(len(df),bool)), ('S-row', (df.dr>0).to_numpy())):
    l0,_ = profile(X0, um)
    l1v = np.array([profile(cellX('M1', D, 0.0), um)[0] for D in 10**logD[::2]])
    l2v = np.array([[profile(cellX('M2', D, kap), um) for kap in kq] for D in 10**logD[::2]])
    np.save(f'res_T4_{name}_ll2.npy', l2v[:,:,0]); 
    print(name, 'M0', l0, 'M1', l1v.max(), 10**logD[::2][l1v.argmax()], 'M2', l2v[:,:,0].max(), 10**logD[::2][np.unravel_index(l2v[:,:,0].argmax(), l2v[:,:,0].shape)[0]])
np.savez('res_T4.npz', logD=logD, ll1=ll1, e1=np.array(e1), ll2=ll2, e2=e2, kq=kq, ll0=ll0, e0=e0)
