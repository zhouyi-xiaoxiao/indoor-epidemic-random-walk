"""The re-check's extra generic competitor: STATIC = pair-specific contact rate learned on train days (no space, no walk),
plus E noise ceiling (real train days vs real test days at equal nominal R0)."""
import sys,json,numpy as np
sys.path.insert(0,'../../bsc_validation2/a2_contact_mobility/src'); import a2lib as A
def static_day(rng,pr,C,mean_dur,N):
    ids,a,b=pr['ids'],pr['a'],pr['b'];n=len(ids)
    iu,ju=np.triu_indices(n,1)
    lo=np.maximum(a[iu],a[ju]);hi=np.minimum(b[iu],b[ju]);L=hi-lo+1
    ok=L>0;iu,ju,lo,L=iu[ok],ju[ok],lo[ok],L[ok]
    c=C[ids[iu],ids[ju]]
    ne=rng.poisson(c*L/mean_dur); rep=np.repeat(np.arange(len(iu)),ne)
    st=lo[rep]+np.floor(rng.random(len(rep))*L[rep]).astype(np.int64)
    du=rng.geometric(1/mean_dur,len(rep)); en=np.minimum(st+du,lo[rep]+L[rep])
    S=np.concatenate([np.arange(x,y) for x,y in zip(st,en)]) if len(st) else np.zeros(0,np.int64)
    I=np.repeat(ids[iu[rep]],en-st);J=np.repeat(ids[ju[rep]],en-st)
    k=np.unique(S*(N*N)+I*N+J)
    return dict(s=k//(N*N),i=(k//N)%N,j=k%N)
out={}
for name in sys.argv[1].split(','):
    ds=A.load(name);N,groups,days=ds['N'],ds['groups'],ds['days']
    pres=[A.presence(d) for d in days];train,test=A.split(len(days))
    cal=A.calibrate(days,pres,N,groups,train);obs=A.stats(days,pres,N,groups,test)
    W=np.zeros((N,N));O=np.zeros((N,N))
    for d in train:
        np.add.at(W,(days[d]['i'],days[d]['j']),1)
        pr=pres[d];a,b=pr['a'],pr['b'];ov=np.clip(np.minimum(b[:,None],b[None,:])-np.maximum(a[:,None],a[None,:])+1,0,None)
        O[np.ix_(pr['ids'],pr['ids'])]+=ov
    C=np.where(O>0,W/np.maximum(O,1),cal['c']); C=np.triu(C,1); C=C+C.T
    exp_tot=lambda dl:sum(float(np.triu(C[np.ix_(pres[d]['ids'],pres[d]['ids'])]*np.clip(np.minimum(pres[d]['b'][:,None],pres[d]['b'][None,:])-np.maximum(pres[d]['a'][:,None],pres[d]['a'][None,:])+1,0,None),1).sum()) for d in dl)
    rng=np.random.default_rng([20261001,99])
    alld=list(range(len(days)));sts=[]
    for g in range(5):
        gen=[static_day(rng,pres[d],C,cal['mean_dur'],N) for d in alld]
        sts.append(A.stats(gen,pres,N,groups,test))
    st={k:float(np.nanmean([s.get(k,np.nan) for s in sts])) for k in sts[0]}
    sc=A.score(st,obs)
    elig=np.unique(np.concatenate([pres[d]['ids'] for d in test]));ne=len(elig)
    tot=obs['S0_total'];scale=tot/exp_tot(test)
    nets=[{d:static_day(rng,pres[d],C*scale,cal['mean_dur'],N) for d in test} for _ in range(10)]
    totratio=np.mean([sum(len(nn[d]['s']) for d in test) for nn in nets])/tot
    # ceiling: real train days looped, own beta
    eligT=np.unique(np.concatenate([pres[d]['ids'] for d in train]));totT=sum(len(days[d]['s']) for d in train)
    res=dict(score=sc,stats=st,scale=scale,totratio=totratio,epi={});npass=0;npassT=0
    for si,(R0,tau) in enumerate(A.EPI_SCEN):
        beta=A.epi_beta(tot,ne,len(test),R0,tau)
        fs,off=A.sir({d:days[d] for d in test},test,N,elig,beta,tau*A.DAY,4000,20261001*1000+si)
        real=A.epi_summary(fs,off,ne)
        F=[];Of=[]
        for g,nn in enumerate(nets):
            f,o=A.sir(nn,test,N,elig,beta,tau*A.DAY,400,777000+si*100+g);F.append(f);Of.append(o)
        e=A.epi_summary(np.concatenate(F),np.concatenate(Of),ne);e['pass']=A.epi_pass(e,real);npass+=e['pass']
        betaT=A.epi_beta(totT,len(eligT),len(train),R0,tau)
        fs,off=A.sir({d:days[d] for d in train},train,N,eligT,betaT,tau*A.DAY,4000,555000+si)
        eT=A.epi_summary(fs,off,len(eligT));eT['pass']=A.epi_pass(eT,real);npassT+=eT['pass']
        # MC noise: second independent real run
        fs,off=A.sir({d:days[d] for d in test},test,N,elig,beta,tau*A.DAY,4000,31337+si)
        r2=A.epi_summary(fs,off,ne)
        res['epi'][f'{R0},{tau}']=dict(real=real,real2=r2,STATIC=e,TRAINREPLAY=eT)
        print(name,R0,tau,'real R %.3f P %.3f A %.3f | rerun %.3f %.3f %.3f | STATIC %.3f %.3f %.3f %s | TRAINdays %.3f %.3f %.3f %s'%(real['R_index'],real['P_major'],real['attack'],r2['R_index'],r2['P_major'],r2['attack'],e['R_index'],e['P_major'],e['attack'],e['pass'],eT['R_index'],eT['P_major'],eT['attack'],eT['pass']),flush=True)
    res['E_pass_count']=npass;res['train_replay_pass_count']=npassT
    print(name,'STATIC S-fam',{k:v for k,v in sc.items()},'E',npass,'/4 totratio %.3f'%totratio,'| train-replay E',npassT,'/4')
    print('   stats',{k:round(v,3) for k,v in st.items() if not k.startswith('S7')})
    out[name]=res
    json.dump(out,open(f'v4_static_{sys.argv[1].replace(",","_")}.json','w'),indent=1,default=float)
