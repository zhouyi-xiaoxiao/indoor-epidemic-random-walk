import numpy as np, vcore as v, json
from scipy import stats
from v03_floor_power import fitb
LY=np.array([0.7619836508613838,0.1044788607478434,0.13353748839077287])
rng=np.random.default_rng(5); M=400000
eps=np.exp(rng.uniform(np.log(3),np.log(30),M))/540; kap=np.exp(rng.uniform(0,np.log(5),M))
L11=-np.log(1-94/216); c=L11/((1-eps)*94/216+kap*eps*94/1143); mu=927*(1-np.exp(-c*kap*eps*94/1143)); cnt=rng.poisson(mu)
print('Seoul mu q',np.quantile(mu,[.025,.5,.975]).round(2),'count q',np.quantile(cnt,[.025,.5,.975]),'P(<=3)',(cnt<=3).mean(),'wm',927*97/1143,stats.binom.cdf(3,927,97/1143))
print('generic pass prob',np.log(8.77/0.62)/np.log(40000))
# prior sensitivity: how wide can the eps prior be moved and still contain 3
for lo,hi in ((1,10),(10,60),(30,120)):
    e=np.exp(rng.uniform(np.log(lo),np.log(hi),M))/540; cc=L11/((1-e)*94/216+kap*e*94/1143); k=rng.poisson(927*(1-np.exp(-cc*kap*e*94/1143))); print('minutes',lo,hi,np.quantile(k,[.025,.5,.975]))
cl,Y=v.load()
def run(tag,**kw):
    D=v.Data(cl,Y,**kw); S=D.S[:,1:]; Yd=Y[:,1:]
    thZ,lZ=fitb(D.feats('Z',w=LY),S,Yd,1e-9); thF,lF=fitb(D.feats('F'),S,Yd,1e-9); w=thF[:3]/thF[:3].sum()
    lH=fitb(D.feats('H'),S,Yd,1e-9)[1]; lX=fitb(D.feats('X'),S,Yd,1e-9)[1]
    print(tag,'w',w.round(3),'TV',round(.5*abs(w-LY).sum(),3),'LR',round(2*(lF-lZ),1),'Z-H',round(lZ-lH,1),'Z-X',round(lZ-lX,1))
run('base'); run('nobreak',brk=[]); run('kern2.5',K=v.kern(2.5,1.5)); run('kern1.2/0.8',K=v.kern(1.2,0.8)); run('break 88-102',brk=range(88,103))
# drop multi-class restriction: only grades with >=2 classes
