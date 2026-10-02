import numpy as np
rng = np.random.default_rng(3)
t4 = np.load('res_T4.npz'); ll = t4['ll2']; w = np.exp(ll-ll.max()); w/=w.sum(); e2 = t4['e2']
flat = rng.choice(w.size, 4000, p=w.ravel()); Emean = e2.ravel()[flat]*(0.5*1.0*2.4)/0.5
print('T4 mean emission q/h 5/50/95:', np.quantile(Emean,[.05,.5,.95]))
def sim(mask, n=20000):
    out = np.zeros(n, int)
    ncon = np.array([15]*14 + [14]*37); assert ncon.sum()==728
    for r in range(n):
        Em = Emean[rng.integers(len(Emean))]; E = np.exp(np.log(Em)-0.5*0.94**2 + 0.94*rng.standard_normal(51))
        V = np.exp(rng.uniform(np.log(180),np.log(350),51)); T = rng.uniform(6.5,13,51); k = np.exp(rng.uniform(0,np.log(6),51))+0.93
        dose = mask*E*0.5/V*(T/k-(1-np.exp(-k*T))/k**2)
        out[r] = rng.binomial(ncon, 1-np.exp(-dose)).sum()
    return out
for m in (0.35, 0.2, 0.6, 1.0):
    x = sim(m); print(f'mask {m}: mean {x.mean():.2f}, P(X>=5)={np.mean(x>=5):.3f}, P(X<=5)={np.mean(x<=5):.3f}, 5%/95% quantiles {np.quantile(x,[.05,.95])}')
