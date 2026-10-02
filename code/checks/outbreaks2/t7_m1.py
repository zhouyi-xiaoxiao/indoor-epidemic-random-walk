import json, numpy as np, os
from scipy.special import gammainc
src=open('t2_analysis.py').read().split("# ------------------------------------------------------------------ 1. primary")[0]
exec(src)
def kern_m1(g, Dgrid):
    nx,ny,ax,ay,T=g['nx'],g['ny'],g['ax'],g['ay'],g['T']
    lam0=2/ax**2+2/ay**2; px,py=(1/ax**2)/lam0,(1/ay**2)/lam0
    S=len(g['srcsets']); out=np.zeros((len(Dgrid),S,nx*ny))
    mx,Vx=np.linalg.eigh(-V._lap1d(nx)); my,Vy=np.linalg.eigh(-V._lap1d(ny)); mx=np.clip(mx,0,None); my=np.clip(my,0,None)
    uni=[i for i,D in enumerate(Dgrid) if 2*D*lam0*T<=V.CUT]
    for i,D in enumerate(Dgrid):
        if i in uni: continue
        lamk=D*(my[:,None]/ay**2+mx[None,:]/ax**2); x=2*lamk*T
        W=np.where(x<1e-8,T,-np.expm1(-x)/np.where(lamk>0,2*lamk,1))
        for s,srcs in enumerate(g['srcsets']):
            acc=np.zeros((ny,nx))
            for (x0,y0) in srcs: acc+=(Vy*Vy[y0][None,:])@W@(Vx*Vx[x0][None,:]).T
            out[i,s]=acc.reshape(-1)
    if uni:
        lam=np.array([2*Dgrid[i]*lam0 for i in uni]); kmax=(lam*T+14*np.sqrt(lam*T)+80).astype(int)
        u=np.zeros((S,ny,nx))
        for s,srcs in enumerate(g['srcsets']):
            for (x0,y0) in srcs: u[s,y0,x0]+=1
        nbx=np.full(nx,2.);nbx[0]=nbx[-1]=1;nby=np.full(ny,2.);nby[0]=nby[-1]=1
        stay=np.clip(1-px*nbx[None,:]-py*nby[:,None],0,None)
        acc=np.zeros((len(uni),S,ny*nx)); CH=400; k0=0; K=int(kmax.max())
        while k0<=K:
            ks=np.arange(k0,min(k0+CH,K+1)); U=np.zeros((len(ks),S,ny,nx))
            for j in range(len(ks)):
                U[j]=u; v=stay[None]*u
                v[:,:,1:]+=px*u[:,:,:-1]; v[:,:,:-1]+=px*u[:,:,1:]; v[:,1:,:]+=py*u[:,:-1,:]; v[:,:-1,:]+=py*u[:,1:,:]; u=v
            act=np.flatnonzero(kmax>=k0)
            if len(act):
                c=gammainc(ks[None,:]+1,lam[act][:,None]*T)/lam[act][:,None]
                acc[act]+=np.einsum('ak,ksn->asn',c,U.reshape(len(ks),S,ny*nx),optimize=True)
            k0+=CH
        for j,i in enumerate(uni): out[i]=acc[j]
    return out*g['nseg']
Dg=10.0**V.LOGD
for e in EV:
    fn='cache/m1_%s.npz'%e
    if os.path.exists(fn): T[e]['M1']=np.load(fn)['pmf']; continue
    g=G[e]; X=kern_m1(g,Dg); comps=T[e]['comps']; pmf=np.zeros((len(Dg),len(comps)))
    for i in range(len(Dg)):
        for (iS,ik,rec) in g['cfgs']: pmf[i]+=V.cond_pmf(X[i,iS,rec],g['bins'],g['sizes'],g['K'],comps)
    pmf/=len(g['cfgs']); np.savez_compressed(fn,pmf=pmf); T[e]['M1']=pmf; print(e,'done',flush=True)
_summ=summ
r=run(CAL,HOLD,target='M1',label='M1cf (positive kernel)')
for e in r['events']: print('  ',e,{m:(round(v['logp'],2),round(v['p_adeq'],4),v['exp']) for m,v in r['events'][e].items() if m!='obs'})
for c in CAL:
    w,ll=post('M1',CAL[c]); print(c,'M1 posterior',summ(w,V.LOGD),'maxll',ll.max())
json.dump(strip(r),open('t7_m1.json','w'),indent=1,default=float)
