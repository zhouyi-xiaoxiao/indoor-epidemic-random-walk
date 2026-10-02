import json,numpy as np,pandas as pd,itertools
from advlib import *
rng=np.random.default_rng(777)
PM={}; OBS={'T2':3,'T1':14,'T5':11,'C1':8,'R1':2}
def base(n1,n2,K): return {'M0':m0(n1,n2,K),'CRR2':crr(n1,n2,K,2),'CRR3':crr(n1,n2,K,3)}
# ---------- T2: Fig S4A values typed again for the re-check from the rendered figure
S4={'13D':.74,'11D':.70,'10D':.58,'9D':.58,'8D':.44,'7D':.45,'6D':.43,'5D':.43,'4D':.42,'3D':.41,'2D':.41,'1D':.41,
 '13C':.66,'12C':1.00,'11C':.55,'10C':.52,'9C':.51,'7C':.46,'6C':.43,'5C':.43,'4C':.42,'3C':.41,'2C':.41,'1C':.41,'13E':.56,
 '13B':.47,'12B':.53,'11B':.46,'10B':.44,'7B':.42,'6B':.42,'5B':.42,'4B':.41,'3B':.41,'1B':.41,
 '13A':.47,'12A':.48,'11A':.45,'10A':.44,'7A':.42,'6A':.42,'5A':.42,'4A':.41,'3A':.41,'1A':.41}
assert len(S4)==45
seats=list(S4); X=np.array([S4[s] for s in seats]); near=np.array([int(s[:-1])>=8 for s in seats]); assert near.sum()==19
cases=['1D','5C','6A','6D','9C','9D','13D']; assert sum(near[seats.index(c)] for c in cases)==3
# sigma from S3A measured/CFD typed again for the re-check
ms={'1D':(.70,.56),'5C':(.76,.58),'6B':(.82,.58),'9C':(.92,.75),'11B':(.99,.62),'11C':(1.03,.68),'13A':(.72,.64),'13D':(1.00,2.00)}
sig=np.std([np.log(a/b) for a,b in ms.values()],ddof=1); print('sigma_log',sig)
pm=base(19,26,7); pm['KA']=nearpmf(X,near,7)
pm['P']=np.mean([nearpmf(X*np.exp(sig*rng.standard_normal(45)),near,7) for _ in range(300)],axis=0)
# K-B with D=1259 on coach lattice, kappa 4.8+0.93, T=200/60
coach=Lat(5,13,0.5,11.4/13); cc={'A':0,'B':1,'E':2,'C':3,'D':4}
rec=[coach.site(cc[s[-1]],int(s[:-1])-1) for s in seats]; src=coach.site(4,11)
for D in(1259,25.1,6.3): pm['KB_D%g'%D]=nearpmf(coach.expo(D,5.73,200/60,src,rec),near,7)
PM['T2']=pm
# ---------- T5
km=json.load(open('kinahan_fwd_means.json'))
OCC={'K':[2,3,4],'G':[1,3,4,5,6,7],'D':[1,3,4,5,6,7],'A':[1,2,3,4,5]}
s5=[f'{r}{c}' for c,rs in OCC.items() for r in rs]; NEAR=['3K','4K','3G','4G','5G','6G','7G','3D','4D','5D','6D','7D']
near5=np.array([s in NEAR for s in s5]); assert len(s5)==20 and near5.sum()==12
def kern(rel,mirror,zero_missing=False):
    mp={'A':'L','D':'G','G':'D','L':'A'}
    def get(r,c):
        c='L' if c=='K' else c
        if mirror: c=mp[c]
        v=km[rel].get(f'{r}{c}'); 
        if v is None or (zero_missing and v==0): return np.nan
        return v
    vm={'A':'K','K':'A','D':'G','G':'D'}; out=[]
    for s in s5:
        r,c=int(s[:-1]),s[-1]; v=get(r,c)
        if not np.isfinite(v): v=get(r,vm[c])
        if not np.isfinite(v): v=np.nanmean([get(r,x) for x in 'ADGK'])
        out.append(v)
    return np.array(out)
k5L,k5A=kern('5L',False),kern('5A',True)
print('5A-mirrored',dict(zip(s5,np.round(k5A,1))))
pm=base(12,8,12)
def npz(X,near,K,dose='exp'): return nearpmf(np.clip(X,1e-300,None),near,K,dose)
pm['KA_5L']=npz(k5L,near5,12); pm['KA_5A']=npz(k5A,near5,12); pm['P']=(pm['KA_5L']+pm['KA_5A'])/2
pm['P_zeromiss']=(npz(kern('5L',False,True),near5,12)+npz(kern('5A',True,True),near5,12))/2
# pooled-kernel alternative: geometric / arithmetic mean of normalised kernels instead of mixture of pmfs
pm['KA_avgkernel']=npz(k5L/k5L.mean()+k5A/k5A.mean(),near5,12)
cab=Lat(6,46,0.9,1.1); ci={'A':0,'D':2,'G':3,'K':5}
rec5=[cab.site(ci[s[-1]],int(s[:-1])+1) for s in s5]; src5=cab.site(5,6)
pm['KB']=np.mean([npz(cab.expo(D,k+0.93,10.0,src5,rec5),near5,12) for D in(707.9,112.2) for k in np.linspace(30,36,7)],axis=0)
pm['KB_D0.1']=npz(cab.expo(0.1,20.93,10.0,src5,rec5),near5,12)
PM['T5']=pm
# ---------- R1 (typed from Li Table 1 and 3 fetched again for the re-check): table: patrons, overlap, measured, cfd, aerosol
R={'TB':(4,53,.87,1.04,.76),'TC':(7,75,.98,.93,.89),'T05':(2,52,None,.62,.07),'T06':(4,82,None,.47,.13),'T07':(3,69,None,.42,.04),'T08':(2,55,None,.42,.06),
 'T09':(10,75,None,.32,.04),'T10':(6,82,.55,.52,.08),'T11':(7,70,None,.57,.12),'T12':(2,64,None,.50,.09),'T13':(6,50,None,.55,.05),'T14':(3,61,None,.63,.11),
 'T15':(8,82,.58,.54,.23),'T16':(5,48,.70,.56,.06),'T17':(5,23,.86,.75,.47),'T18':(5,77,.73,.85,.40)}
Xr=np.concatenate([np.full(n,(m if m is not None else c)*t/60) for n,t,m,c,a in R.values()])
nr=np.concatenate([np.full(n,k in('TB','TC','T18')) for k,(n,*_) in R.items()]); assert nr.sum()==16 and (~nr).sum()==63
R1K={}
for K in(2,3,4,5):
    d=base(16,63,K); d['P']=nearpmf(Xr,nr,K); R1K[K]=d
PM['R1']=R1K[2]
# ---------- T1 (own geometry draw) and C1
def t1(D):
    lat=Lat(6,15,2.5/6,0.75); sc=[0,1,2,4,5]; out=[]
    rows=np.repeat(np.arange(15),5); cols=np.tile(sc,15)
    for i in range(150):
        occ=np.ones(75,bool); srcm=(rows==7)&(cols==1); occ[srcm]=False
        za=(rows>=5)&(rows<=9); zb=(rows==4)|(rows==10); zc=~(za|zb)
        occ[rng.choice(np.flatnonzero(za&~srcm),1)]=False; occ[rng.choice(np.flatnonzero(zc),6,replace=False)]=False
        kap=np.exp(rng.uniform(np.log(.5),np.log(5)))+.93
        Dd=D if np.isscalar(D) else D[i%len(D)]
        out.append(nearpmf(lat.expo(Dd,kap,[50/60,50/60],lat.site(1,7),lat.site(cols[occ],rows[occ])),(za|zb)[occ],23))
    return np.mean(out,axis=0)
kj=json.load(open('../../../bsc_validation2/a3_tracer_physics/results/01_kernels.json')); Db=np.array(kj['coach']['boot_draws'])
pm=base(33,34,23); pm['P']=t1(Db); pm['P_D1259']=t1(1259.); PM['T1']=pm
def c1(f=1.0):
    lat=Lat(11,13,1,1); xs=np.array([1,3,5,7,9]); ys=np.array([2,4,6,8,10]); cols=np.tile(xs,5); rows=np.repeat(ys,5); front=rows<=4; out=[];Ds=[]
    for i in range(300):
        occ=np.ones(25,bool); occ[rng.choice(np.flatnonzero(~front),1)]=False
        ach=np.exp(rng.uniform(np.log(2),np.log(12))); D=(0.52*ach+0.31)*(143*3.0)**(2/3)*np.exp(.42*rng.standard_normal())*f; Ds.append(D)
        out.append(nearpmf(lat.expo(D,ach+.93,[6.5,6.5],lat.site(5,0),lat.site(cols[occ],rows[occ])),front[occ],12))
    return np.mean(out,axis=0),np.median(Ds)
pm=base(10,14,12); pm['P'],Dm=c1(); print('C1 median D',Dm); PM['C1']=pm
# ---------- report
N={'T2':(19,26,7),'T1':(33,34,23),'T5':(12,8,12),'C1':(10,14,12),'R1':(16,63,2)}
for e in PM:
    print(e,'obs',OBS[e])
    for m,p in PM[e].items(): print('   %-14s p=%.4f  logscore=%.3f  RR=%.2f'%(m,p2(p,OBS[e]),np.log(max(p[OBS[e]],1e-300)),rr(p,*N[e])))
print('R1 by K (obs k=K):',{K:{m:round(p2(p,K),4) for m,p in d.items()} for K,d in R1K.items()})
EV=['T2','T1','T5','C1','R1']
def pooled(pk,comp,evs=EV):
    d=sum(np.log(PM[e][pk][OBS[e]])-np.log(PM[e][comp][OBS[e]]) for e in evs); tot=0
    for combo in itertools.product(*[np.flatnonzero(PM[e][comp]>0) for e in evs]):
        w=np.prod([PM[e][comp][k] for e,k in zip(evs,combo)])
        dd=sum(np.log(max(PM[e][pk][k],1e-300))-np.log(PM[e][comp][k]) for e,k in zip(evs,combo))
        if dd>=d-1e-12: tot+=w
    return round(d,3),tot
for c in('M0','CRR2','CRR3'): print('pooled P vs',c,pooled('P',c),' without T5:',pooled('P',c,['T2','T1','C1','R1']))
print('pooled CRR2 vs M0',pooled('CRR2','M0'),'CRR3 vs M0',pooled('CRR3','M0'))
json.dump({e:{m:p.tolist() for m,p in d.items()} for e,d in PM.items()},open('pmfs.json','w'))
