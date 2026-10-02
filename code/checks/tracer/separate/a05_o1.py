import numpy as np,pandas as pd
from advlib import *
df=pd.read_csv('../../../bsc_outbreaks/data/park2020_callcentre_seats_digitized.csv'); nw=df[df.region=='north_wing']
xy=nw[['x_pitch_units','y_pitch_units']].values; case=nw.case.values.astype(int); n=len(xy); print(n,case.sum())
d=np.hypot(xy[:,None,0]-xy[None,:,0],xy[:,None,1]-xy[None,:,1]); A=((d<=1.25)&(d>0)).astype(float)
def z(c,perm_rng=None):
    return c@A@c/2
rng=np.random.default_rng(5)
jj=z(case); perm=np.array([z(rng.permutation(case)) for _ in range(20000)]); mu,sd=perm.mean(),perm.std()
print('obs joins',jj,'perm mean',mu,'sd',sd,'z',(jj-mu)/sd)
ix=np.rint(xy[:,0]-xy[:,0].min()+1).astype(int); iy=np.rint(xy[:,1]-xy[:,1].min()+1).astype(int)
def simz(D,ach,pitch,R=1500):
    lat=Lat(ix.max()+2,iy.max()+2,pitch,pitch); s=lat.site(ix,iy); w=1/(ach+.93+D*lat.lam); P=lat.V[s]; X=np.clip((P*w)@P.T,0,None); np.fill_diagonal(X,0)
    def perc(eps,R):
        inf=np.zeros((R,n),bool); new=np.zeros((R,n),bool); new[np.arange(R),rng.integers(n,size=R)]=True; inf|=new
        while new.any():
            p=1-np.exp(-eps*(new.astype(float)@X)); nxt=(rng.random((R,n))<p)&~inf; inf|=nxt; new=nxt
        return inf
    lo,hi=np.log(1/X.sum(1).mean())-2,np.log(1/X.sum(1).mean())+12
    for _ in range(14):
        mid=(lo+hi)/2; fs=perc(np.exp(mid),400).sum(1); big=fs[fs>=40]
        if (big.mean() if len(big)>=5 else 0)<79: lo=mid
        else: hi=mid
    inf=perc(np.exp((lo+hi)/2),R); fs=inf.sum(1); keep=inf[(fs>=69)&(fs<=89)]
    zs=[]
    for c in keep:
        c=c.astype(float); k=int(c.sum()); 
        # moments by formula-free approach: use permutation mean/sd scaled approx for k -> use exact moments
        deg=A.sum(1); Pn=A.sum()/2; S1=(deg*(deg-1)).sum()/2; S0=Pn*(Pn-1)/2-S1
        p2_=k*(k-1)/(n*(n-1)); p3=p2_*(k-2)/(n-2); p4=p3*(k-3)/(n-3); m=Pn*p2_; v=Pn*p2_+2*S1*p3+2*S0*p4-m*m
        zs.append((c@A@c/2-m)/np.sqrt(v))
    return np.array(zs)
zo=(jj-mu)/sd
for D,ach,pitch in((185,2.45,1.3),(91,1.0,1.3),(376,6.0,1.3),(185,2.45,1.0),(1e9,2.45,1.3),(0.5,2.45,1.3)):
    zs=simz(D,ach,pitch); print('D',D,'ach',ach,'pitch',pitch,'n',len(zs),'median z',round(np.median(zs),2),'q2.5,97.5',np.round(np.quantile(zs,[.025,.975]),2),'P(z<=obs)',round((zs<=zo).mean(),3))
