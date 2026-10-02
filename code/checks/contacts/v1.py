import sys,json,numpy as np,pandas as pd,time
import vlib as V
name=sys.argv[1]; scen=[(1.5,1.0),(3.0,4.0)] if len(sys.argv)<3 else eval(sys.argv[2])
D,N,grp,nd=V.load(name)
ntest=nd//2; train=list(range(nd-ntest)); test=list(range(nd-ntest,nd))
pres={d:V.presence(D[D.day==d]) for d in range(nd)}
tot_tr=sum((D.day==d).sum() for d in train); ov_tr=sum(V.overlap(pres[d]) for d in train)
c=tot_tr/ov_tr; dur=np.concatenate([V.events(D[D.day==d]).dur.values for d in train]); m=dur.mean()
p=V.p_from_q(1-1/m); M=1/c
obs=V.stats(D,pres,N,grp,test)
out=dict(name=name,N=N,nd=nd,c=c,M=M,mean_dur=m,p=p,obs=obs)
print(json.dumps(out,default=float),flush=True)
# RW0 structure: 3 realisations of test days (+ previous day for persistence)
def dims(M):
    M=max(int(round(M)),4); Lx=max(int(round(np.sqrt(M))),2); Ly=max(int(round(M/Lx)),2); return Lx,Ly
rng=np.random.default_rng(4242)
need=sorted(set(test)|{test[0]-1})
sts=[]
for g in range(3):
    Lx,Ly=dims(M); fr=[]
    for d in need:
        x=V.rw0_day(rng,pres[d],Lx,Ly,p); x['day']=d; fr.append(x)
    G=pd.concat(fr); sts.append(V.stats(G,pres,N,grp,test))
rw={k:float(np.mean([s[k] for s in sts])) for k in sts[0]}
out['RW0']=rw; print('RW0',json.dumps(rw),flush=True)
# epidemic
elig=np.unique(np.concatenate([pres[d].id.values for d in test])); ne=len(elig)
tot=obs['S0_total']; ov=sum(V.overlap(pres[d]) for d in test); scale=(tot/ov)/c
P=len(test)*V.DAY
def flat(G):
    S=np.concatenate([G[G.day==d].s.values+k*V.DAY for k,d in enumerate(test)]); I=np.concatenate([G[G.day==d].i.values for d in test]);J=np.concatenate([G[G.day==d].j.values for d in test]); return S,I,J
out['epi']={}
for R0,tau in scen:
    W=2*tot/(ne*len(test)*V.DAY); beta=R0/(W*tau*V.DAY)
    t0=time.time(); fs,off=V.sir_vec(*flat(D),N,elig,P,beta,tau*V.DAY,2000,rng)
    real=dict(R=off.mean(),Pm=(fs/ne>=0.1).mean(),att=(fs/ne).mean())
    F=[];O=[];tr=[]
    for g in range(4):
        Lx,Ly=dims(M/scale); fr=[]
        for d in test:
            x=V.rw0_day(rng,pres[d],Lx,Ly,p); x['day']=d; fr.append(x)
        G=pd.concat(fr); tr.append(len(G)/tot)
        f,o=V.sir_vec(*flat(G),N,elig,P,beta,tau*V.DAY,500,rng); F.append(f);O.append(o)
    F=np.concatenate(F);O=np.concatenate(O)
    mod=dict(R=O.mean(),Pm=(F/ne>=0.1).mean(),att=(F/ne).mean(),totratio=np.mean(tr))
    out['epi'][f'{R0},{tau}']=dict(beta=beta,real=real,RW0=mod)
    print(R0,tau,'beta',beta,'real',real,'RW0',mod,f'{time.time()-t0:.0f}s',flush=True)
json.dump(out,open(f'v1_{name}.json','w'),indent=1,default=float)
