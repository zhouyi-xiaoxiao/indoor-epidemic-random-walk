"""Re-implementation with separately written code (re-check). No import from valmod / a3lib."""
import numpy as np, itertools
from scipy import optimize, stats
LOGD=np.round(np.arange(-1,5.0001,0.05),2); DGRID=10**LOGD
class Lat:
    def __init__(s,nx,ny,ax,ay):
        s.nx,s.ny,s.ax,s.ay=nx,ny,ax,ay; M=nx*ny
        Lp=np.zeros((M,M))
        for y in range(ny):
            for x in range(nx):
                i=y*nx+x
                for dx,dy,w in((1,0,1/ax**2),(-1,0,1/ax**2),(0,1,1/ay**2),(0,-1,1/ay**2)):
                    xx,yy=x+dx,y+dy
                    if 0<=xx<nx and 0<=yy<ny:
                        j=yy*nx+xx; Lp[j,i]+=w; Lp[i,i]-=w
        s.lam,s.V=np.linalg.eigh(-Lp)      # -Lap = V diag(lam) V^T, lam>=0
        s.lam=np.clip(s.lam,0,None)
    def site(s,x,y): return np.asarray(y)*s.nx+np.asarray(x)
    def steady(s,D,kap,src,rec):
        w=1/(kap+D*s.lam); return s.V[np.asarray(rec)]@(w*s.V[src])
    def expo(s,D,kap,T,src,rec):
        """int_0^T (T-u) exp(-kap u) G(u) du, summed over segments T"""
        L=kap+D*s.lam; w=np.zeros_like(L)
        for t in np.atleast_1d(T): w+=(L*t-1+np.exp(-L*t))/L**2
        return s.V[np.asarray(rec)]@(w*s.V[src])
    def expo_mat(s,D,kap,T,sites):
        L=kap+D*s.lam; w=np.zeros_like(L)
        for t in np.atleast_1d(T): w+=(L*t-1+np.exp(-L*t))/L**2
        P=s.V[np.asarray(sites)]; return (P*w)@P.T
def fitD(lat,src,rec,c,kap):
    y=np.log(c); sse=[]
    for D in DGRID:
        m=np.log(np.clip(lat.steady(D,kap,src,rec),1e-300,None)); r=y-m; r-=r.mean(); sse.append((r*r).sum())
    sse=np.array(sse); return DGRID[sse.argmin()],sse
def pb(p):
    f=np.array([1.0])
    for q in p: f=np.convolve(f,[1-q,q])
    return f
def cond(pn,pf,K):
    a,b=pb(pn),pb(pf); out=np.array([a[k]*b[K-k] if 0<=K-k<len(b) else 0. for k in range(len(a))]); return out/out.sum()
def prof(X,K,dose='exp'):
    X=np.asarray(X,float)
    if dose=='exp': g=lambda e:1-np.exp(-e*X)
    elif dose=='lin': g=lambda e:np.minimum(1,e*X)
    else: g=lambda e:e*X/(1+e*X)
    Xp=X[X>0]
    le=optimize.brentq(lambda l:g(np.exp(l)).sum()-K,np.log(1e-12/Xp.max()),np.log(1e6/Xp.min()),xtol=1e-13)
    return g(np.exp(le)),np.exp(le)
def nearpmf(X,near,K,dose='exp'):
    p,_=prof(X,K,dose); return cond(p[near],p[~near],K)
def crr(n1,n2,K,rho):
    q=K/(n1*rho+n2); pn=rho*q
    if pn>1: pn=1; q=(K-n1)/n2
    return cond(np.full(n1,pn),np.full(n2,q),K)
def m0(n1,n2,K): return stats.hypergeom.pmf(np.arange(n1+1),n1+n2,n1,K)
def p2(pmf,k): return float(pmf[pmf<=pmf[k]*(1+1e-9)+1e-12].sum())
def rr(pmf,n1,n2,K):
    e=(pmf*np.arange(len(pmf))).sum(); return (e/n1)/((K-e)/n2)
