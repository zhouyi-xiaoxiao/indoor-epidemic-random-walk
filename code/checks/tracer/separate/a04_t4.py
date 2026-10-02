import numpy as np,json
from advlib import *
from scipy.stats import gamma,chi2
dr=np.array([0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3]); dc=np.array([2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5])
k=np.array([33,7,7,3,10,13,5,3,1,1,11,9,8,6,3,3,2,2,4,3,3,1.]); n=np.array([1998,1848,1824,1033,4839,5420,3653,3524,3537,1884,4372,5179,3427,3637,3263,1760,4048,4286,3112,2936,2936,1570.])
assert k.sum()==138
lat=Lat(6,17,0.5,1.0); rows=np.repeat(np.arange(17),5); cols=np.tile([0,1,2,4,5],17); sites=lat.site(cols,rows)
DR=abs(rows[:,None]-rows[None,:]); DC=abs(cols[:,None]-cols[None,:])
sh=(2.1/1.8)**2; sc=1.8**2/2.1; Tn=gamma.ppf((np.arange(16)+.5)/16,sh,scale=sc)
def pi(D,kap):
    Xc=np.zeros(22)
    for T in Tn:
        X=lat.expo_mat(D,kap,T,sites)
        Xc+=np.array([X[(DR==a)&(DC==b)].mean() for a,b in zip(dr,dc)])/16
    w=n*Xc; return w/w.sum()       # rare-event (linear) limit
def dev(k,p): e=k.sum()*p; return 2*np.sum(np.where(k>0,k*np.log(np.where(k>0,k,1)/e),0))
kq=np.exp(np.log(5)+(np.arange(7)+.5)/7*np.log(4))+.93
for D in(891.,355.,25.1,2.0,0.5):
    p=np.mean([pi(D,kk) for kk in kq],axis=0)
    print('D',D,'dev',round(dev(k,p),1),'chi2 p',chi2.sf(dev(k,p),21),'loglik-M0',round(k @ np.log(p)-k @ np.log(n/n.sum()),2),'E same row',round(138*p[dr==0].sum(),1))
for rho in(1,2,3,5,8):
    w=n*np.where(dr==0,rho,1); p=w/w.sum(); print('CRR',rho,'dev',round(dev(k,p),1),'E row0',round(138*p[dr==0].sum(),1))
# saturated-by-row model (4 free row levels): is ANY row-only model adequate?
w=n*np.array([k[dr==r].sum()/n[dr==r].sum() for r in dr]); p=w/w.sum(); print('row-saturated dev',round(dev(k,p),1),'df',18,chi2.sf(dev(k,p),18))
