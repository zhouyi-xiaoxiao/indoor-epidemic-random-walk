"""Re-implementation with separately written code (re-check). Class-day chain binomial. Written from PREREGISTRATION.md sec 2.3 only."""
import numpy as np, pandas as pd
from scipy import stats, optimize
T=213
def kern(mean=1.7,sd=1.0,m=10):
    sh=(mean/sd)**2; sc=sd**2/mean
    c=stats.gamma.cdf(np.arange(0,m+1),sh,scale=sc); k=np.diff(c); return k/k.sum()
def load(path='v_students.csv'):
    df=pd.read_csv(path)
    df['cid']=df.groupby(['schoolID','gradeID','classID']).ngroup()
    cl=df.groupby('cid').agg(school=('schoolID','first'),grade=('gradeID','first'),n=('onset','size')).reset_index()
    Y=np.zeros((len(cl),T+1))
    inf=df[df.isinfected==1]
    for c,t in zip(inf.cid.values,inf.onset.values): Y[c,t]+=1
    return cl,Y
class Data:
    def __init__(s,cl,Y,K=None,brk=range(89,101)):
        s.cl=cl; s.Y=Y; K=kern() if K is None else K; s.K=K
        nC=len(cl); s.n=cl.n.values.astype(float); s.sch=cl.school.values; s.gk=(cl.school*100+cl.grade).values
        A=np.zeros((nC,T+1))
        for t in range(1,T+1):
            for tau in range(1,len(K)+1):
                if t-tau>=1: A[:,t]+=K[tau-1]*Y[:,t-tau]
        s.A=A
        def gs(v,key):
            return pd.DataFrame(v).groupby(key).transform('sum').values
        s.Ag=gs(A,s.gk); s.As=gs(A,s.sch)
        s.ng=pd.Series(s.n).groupby(s.gk).transform('sum').values; s.ns=pd.Series(s.n).groupby(s.sch).transform('sum').values
        s.N=s.n.sum()
        n,ng,ns=s.n[:,None],s.ng[:,None],s.ns[:,None]
        def sd(a,b): 
            out=np.zeros_like(a); np.divide(a,b,out=out,where=np.broadcast_to(b>0,a.shape)); return out
        s.Pc=sd(A,n-1); s.Pg=sd(s.Ag-A,ng-n); s.Ps=sd(s.As-s.Ag,ns-ng); s.PH=sd(s.As,ns-1)
        s.city=(A.sum(0,keepdims=True)-s.As)/(s.N-ns)
        s.inbreak=np.zeros(T+1,bool); s.inbreak[list(brk)]=True
        cum=np.cumsum(Y,1); s.S=n-(cum-Y)   # susceptibles at start of day
    def x(s,model,w=None,rho=None):
        if model=='H': X=s.PH.copy()
        elif model=='X': X=s.Pc.copy()
        elif model=='Z': X=w[0]*s.Pc+w[1]*s.Pg+w[2]*s.Ps
        elif model=='C':
            n,ns=s.n[:,None],s.ns[:,None]; X=(rho*s.A+(s.As-s.A))/(rho*(n-1)+ns-n)
        X[:,s.inbreak]=0; return X
    def feats(s,model,w=None,rho=None):
        if model=='F':
            L=[s.Pc.copy(),s.Pg.copy(),s.Ps.copy()]
            for a in L: a[:,s.inbreak]=0
        else: L=[s.x(model,w,rho)]
        L+= [np.ones_like(s.city*s.Pc), np.broadcast_to(s.city,s.Pc.shape)]
        return np.stack(L,-1)[:,1:]
def ll_cls(th,X,S,Y):
    lam=np.maximum(X@th,1e-300)
    return (-(S-Y)*lam+Y*np.log1p(-np.exp(-lam))).sum(1)
def fit(X,S,Y):
    p=X.shape[-1]; best=None
    def f(u):
        th=np.exp(u); lam=np.maximum(X@th,1e-300); e=np.exp(-lam); pr=1-e
        l=(-(S-Y)*lam+Y*np.log(pr)).sum(); r=-(S-Y)+Y*e/pr
        return -l, -((X*r[...,None]).sum((0,1)))*th
    for a0 in (-9,-12):
      for b0 in (0.3,1.0):
        u0=np.log(np.r_[np.full(p-2,b0/(p-2)),np.exp(a0),0.3])
        r=optimize.minimize(f,u0,jac=True,method='BFGS',options=dict(gtol=1e-6,maxiter=3000))
        if best is None or r.fun<best.fun: best=r
    return np.exp(best.x),-best.fun
def fitC(D,S,Y,m):
    def prof(lr): return -fit(D.feats('C',rho=np.exp(lr))[m],S[m],Y[m])[1]
    r=optimize.minimize_scalar(prof,bounds=(0,np.log(3000)),method='bounded',options=dict(xatol=0.01))
    rho=float(np.exp(r.x)); th,l=fit(D.feats('C',rho=rho)[m],S[m],Y[m]); return th,l,rho
def disp(cl,Y):
    I=Y.sum(1); n=cl.n.values; num=0; df=0
    for s in np.unique(cl.school.values):
        j=cl.school.values==s
        if j.sum()<2: continue
        p=I[j].sum()/n[j].sum()
        if p<=0 or p>=1: continue
        num+=((I[j]-n[j]*p)**2/(n[j]*p*(1-p))).sum(); df+=j.sum()-1
    return num/df
def simulate(cl,model,th,rng,w=None,rho=None,K=None,brk=range(89,101)):
    K=kern() if K is None else K
    n=cl.n.values.astype(float); nC=len(n); sch=pd.factorize(cl.school)[0]; gk=pd.factorize(cl.school*100+cl.grade)[0]
    ng=np.bincount(gk,n)[gk]; ns=np.bincount(sch,n)[sch]; N=n.sum()
    Y=np.zeros((nC,T+1)); S=n.copy().astype(int); b=th[:-2]; a0,a1=th[-2:]; brk=set(brk)
    for t in range(1,T+1):
        a=np.zeros(nC)
        for tau in range(1,len(K)+1):
            if t-tau>=1: a+=K[tau-1]*Y[:,t-tau]
        as_=np.bincount(sch,a)[sch]; ag=np.bincount(gk,a)[gk]
        lam=a0+a1*(a.sum()-as_)/(N-ns)
        if t not in brk:
            dv=lambda u,v: np.where(v>0,u/np.where(v>0,v,1),0)
            if model=='H': lam=lam+b[0]*dv(as_,ns-1)
            elif model=='X': lam=lam+b[0]*dv(a,n-1)
            elif model=='C': lam=lam+b[0]*(rho*a+as_-a)/(rho*(n-1)+ns-n)
            else:
                ww=b[0]*np.asarray(w) if model=='Z' else b
                lam=lam+ww[0]*dv(a,n-1)+ww[1]*dv(ag-a,ng-n)+ww[2]*dv(as_-ag,ns-ng)
        y=rng.binomial(S,1-np.exp(-lam)); Y[:,t]=y; S=S-y
    return Y
