import numpy as np, pandas as pd, json, time, vcore as v, sys
from scipy import optimize
LY=np.array([0.7619836508613838,0.1044788607478434,0.13353748839077287])
FLOOR=1e-9
def fitb(X,S,Y,floor):
    p=X.shape[-1]; best=None
    def f(th):
        lam=np.maximum(X@th,1e-300); e=np.exp(-lam); pr=-np.expm1(-lam)
        l=(-(S-Y)*lam+Y*np.log(pr)).sum(); r=-(S-Y)+Y*e/pr
        return -l,-(X*r[...,None]).sum((0,1))
    for x0 in ([*([0.5/(p-2)]*(p-2)),1e-4,0.5],[*([1.0/(p-2)]*(p-2)),1e-5,0.1],[*([0.1]*(p-2)),1e-3,1.0]):
        r=optimize.minimize(f,np.array(x0),jac=True,method='L-BFGS-B',bounds=[(0,None)]*(p-2)+[(floor,None),(0,None)],options=dict(maxiter=3000,ftol=1e-14,gtol=1e-9))
        if best is None or r.fun<best.fun: best=r
    return best.x,-best.fun
def pipeline(cl,Y,rng,floor=1e-9,nb=1000,nf=60,extra=()):
    D=v.Data(cl,Y); S=D.S[:,1:]; Yd=Y[:,1:]; odd=cl.school.values%2==1; schools=np.unique(cl.school.values)
    specs={'Z':('Z',LY),'H':('H',None),'X':('X',None)}; specs.update({k:('Z',np.array(w)) for k,w in extra})
    ho={k:np.zeros(len(cl)) for k in specs}
    for cal in (odd,~odd):
        for k,(m,w) in specs.items():
            X=D.feats(m,w=w); th,_=fitb(X[cal],S[cal],Yd[cal],floor); ho[k][~cal]=v.ll_cls(th,X[~cal],S[~cal],Yd[~cal])
    ps={k:np.array([ho[k][cl.school.values==s].sum() for s in schools]) for k in ho}
    idx=rng.integers(0,len(schools),(nb,len(schools))); o={}
    for k in ps:
        if k=='Z': continue
        d=ps['Z']-ps[k]; bs=d[idx].sum(1); o['D_'+k]=[float(d.sum()),float(np.quantile(bs,.025)),float(np.quantile(bs,.975))]
    o['tier1']=bool(o['D_H'][1]>0 and o['D_X'][1]>0)
    thZ,lZ=fitb(D.feats('Z',w=LY),S,Yd,floor); thF,lF=fitb(D.feats('F'),S,Yd,floor)
    w=thF[:3]/thF[:3].sum(); o['w']=w.tolist(); o['TV']=float(.5*np.abs(w-LY).sum()); o['LR']=float(2*(lF-lZ)); o['tier2a']=bool(o['LR']<=5.991); o['tier2b']=bool(o['TV']<=.15)
    if nf:
        ds=[v.disp(cl,v.simulate(cl,'Z',thZ,rng,w=LY)) for _ in range(nf)]; ob=v.disp(cl,Y); o['disp']=[float(ob),float(np.quantile(ds,.025)),float(np.quantile(ds,.975))]
        o['tier3']=bool(o['disp'][1]<=ob<=o['disp'][2]); o['overall']=bool(o['tier1'] and o['tier2b'] and o['tier3'])
    return o
if __name__=='__main__':
    cl,Y=v.load(); t0=time.time()
    if sys.argv[1]=='floor':
        res={}
        for fl in (1e-9,1e-7,1e-6,1e-5):
            res[str(fl)]=pipeline(cl,Y,np.random.default_rng(1),floor=fl,nb=10000,nf=0,extra=[('N702',[.7,.1,.2]),('N811',[.8,.1,.1]),('N2/3',[2/3,1/6,1/6])]); print(fl,json.dumps(res[str(fl)]),flush=True)
        json.dump(res,open('v_floor.json','w'),indent=1)
    else:
        A=dict(beta=0.8,a0=2e-5,a1=0.42); out=[]
        jobs=[('Z',None)]*10+[('H',None)]*10+[('X',None)]*10+[('C',37.)]*16+[('C',10.)]*8+[('C',100.)]*8+[('F',[.66,.10,.24])]*10
        for i,(tr,par) in enumerate(jobs):
            rng=np.random.default_rng([4242,i])
            if tr=='C': Ys=v.simulate(cl,'C',np.array([A['beta'],A['a0'],A['a1']]),rng,rho=par)
            elif tr=='F': Ys=v.simulate(cl,'F',np.r_[A['beta']*np.array(par),A['a0'],A['a1']],rng)
            else: Ys=v.simulate(cl,tr,np.array([A['beta'],A['a0'],A['a1']]),rng,w=LY)
            o=pipeline(cl,Ys,rng); o['truth']=tr; o['par']=par; o['AR']=float(Ys.sum()/cl.n.sum()); out.append(o)
            json.dump(out,open('v_power.json','w')); print(i,tr,par,o['tier1'],o['tier2a'],o['tier2b'],o['tier3'],o['overall'],round(o['TV'],3),round(o['LR'],1),round(time.time()-t0),flush=True)
