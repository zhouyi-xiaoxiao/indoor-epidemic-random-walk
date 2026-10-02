import numpy as np, vlib as V
from scipy.stats import chi2 as c2
B, G = V.BETA, V.GAMMA
exec(open('v02_simcheck.py').read().split("N = 40")[0].split("# (a)")[1].split("\n",1)[1])
N=40; r=V.uniform(1,1,N); pex=exact_fs(N,B/(N-1),G)
for seed in (1,2,3,4,5):
    n=200000
    s=V.sim(r,n,seed*1000003)
    emp=np.bincount(s['final'],minlength=N+1); E=pex*n
    chi=0;dof=-1;ce=co=0
    for k in range(1,N+1):
        ce+=E[k];co+=emp[k]
        if ce>=5: chi+=(co-ce)**2/ce;dof+=1;ce=co=0
    m,se=V.mci(s['final'])
    print(seed,'mean %.3f+-%.3f chi2 %.1f/%d p=%.3f'%(m,se,chi,dof,c2.sf(chi,dof)))
