"""Re-implementation with separately written code (re-check). Does not import a2lib."""
import numpy as np, pandas as pd, zipfile, io, os
RAW=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','bsc_validation2','a2_contact_mobility','data','raw')
SPEC={'InVS13':('workplace_InVS_tij.dat.zip','workplace_InVS_metadata.txt',0),
'InVS15':('workplace_InVS15_tij.dat.gz','workplace_InVS15_metadata.txt',0),
'LH10':('hospital_lyon_contacts.dat.gz',None,46800),
'LyonSchool':('primaryschool.csv.gz','primaryschool_metadata.txt',0),
'SFHH':('SFHH_tij.dat.gz',None,0),
'Thiers13':('HighSchool2013_proximity_net.csv.gz','HighSchool2013_metadata.txt',3600)}
DAY=4320
def load(name):
    f,meta,off=SPEC[name]
    p=os.path.join(RAW,f)
    if f.endswith('.zip'):
        z=zipfile.ZipFile(p); df=pd.read_csv(io.BytesIO(z.read(z.namelist()[0])),sep=r'\s+',header=None)
    else: df=pd.read_csv(p,sep=r'\s+',header=None)
    lab={}
    if name=='LH10':
        for a,g in zip(df[1],df[3]): lab[a]=g
        for a,g in zip(df[2],df[4]): lab[a]=g
    elif meta:
        m=pd.read_csv(os.path.join(RAW,meta),sep=r'\s+',header=None); lab=dict(zip(m[0],m[1]))
    ids=sorted(set(df[1])|set(df[2])); ix={v:k for k,v in enumerate(ids)}
    d=pd.DataFrame({'t':df[0]+off,'a':df[1].map(ix),'b':df[2].map(ix)})
    d=d[d.a!=d.b]
    d['i']=d[['a','b']].min(axis=1); d['j']=d[['a','b']].max(axis=1)
    d['day']=d.t//86400; d['s']=(d.t%86400)//20
    cnt=d.groupby('day').size(); keep=cnt[cnt>=0.01*cnt.sum()].index
    d=d[d.day.isin(keep)][['day','s','i','j']].drop_duplicates()
    dmap={v:k for k,v in enumerate(sorted(keep))}; d['day']=d.day.map(dmap)
    grp=None
    if lab:
        names=sorted(set(str(lab.get(x,'NA')) for x in ids)); grp=np.array([names.index(str(lab.get(x,'NA'))) for x in ids])
    return d.reset_index(drop=True),len(ids),grp,len(keep)
def presence(dd):
    """dd: one day's frame -> frame id,a,b"""
    x=pd.concat([dd[['i','s']].rename(columns={'i':'id'}),dd[['j','s']].rename(columns={'j':'id'})])
    g=x.groupby('id').s.agg(['min','max']).reset_index(); g.columns=['id','a','b']; return g
def overlap(pr):
    a=pr.a.values;b=pr.b.values
    ov=np.minimum.outer(b,b)-np.maximum.outer(a,a)+1
    return np.triu(np.clip(ov,0,None),1).sum()
def events(dd):
    x=dd.sort_values(['i','j','s']); new=(x.i.diff()!=0)|(x.j.diff()!=0)|(x.s.diff()!=1)
    x=x.assign(ev=new.cumsum())
    e=x.groupby('ev').agg(i=('i','first'),j=('j','first'),st=('s','min'),dur=('s','size')).reset_index(drop=True)
    return e
def stats(D,pres,N,grp,test):
    """D: frame all days; pres: dict day->presence"""
    durs=[];gaps=[];degs=[];nev=0;npair=0;tot=0;within=0
    strength=np.zeros(N);pany=np.zeros(N,bool)
    for d in test:
        dd=D[D.day==d]; tot+=len(dd); e=events(dd); durs.append(e.dur.values); nev+=len(e); npair+=len(dd[['i','j']].drop_duplicates())
        np.add.at(strength,dd.i.values,1);np.add.at(strength,dd.j.values,1); pany[pres[d].id.values]=True
        # gaps = holes in union of each person's contact intervals
        pe=pd.concat([e[['i','st','dur']].rename(columns={'i':'id'}),e[['j','st','dur']].rename(columns={'j':'id'})])
        for _,g in pe.groupby('id'):
            g=g.sort_values('st'); end=-1
            for st,du in zip(g.st.values,g.dur.values):
                if end>=0 and st>end: gaps.append(st-end)
                end=max(end,st+du)
        pr=pres[d]; lg=set(pr.id[(pr.b-pr.a+1)>=90])
        prs=dd[['i','j']].drop_duplicates(); deg=pd.concat([prs.i,prs.j]).value_counts()
        degs+= [deg.get(k,0) for k in lg]
        if grp is not None: within+=(grp[dd.i.values]==grp[dd.j.values]).sum()
    num=den=0
    for d in test:
        if d==0: continue
        A=set(map(tuple,D[D.day==d-1][['i','j']].drop_duplicates().values))
        both=set(pres[d-1].id)&set(pres[d].id)
        for i,j in D[D.day==d][['i','j']].drop_duplicates().values:
            if i in both and j in both: den+=1; num+=((i,j) in A)
    dur=np.concatenate(durs);gap=np.array(gaps,float);deg=np.array(degs,float);st=strength[pany]
    return dict(S0_total=tot,S1_mean_dur=dur.mean(),S1_timefrac_ge15=dur[dur>=15].sum()/dur.sum(),
        S2_burst=(gap.std()-gap.mean())/(gap.std()+gap.mean()),S2_recurrence=1-npair/nev,
        S3_deg_mean=deg.mean(),S3_deg_cv=deg.std()/deg.mean(),S4_strength_cv=st.std()/st.mean(),
        S5_persistence=num/den if den else np.nan,S6_within_frac=within/tot if grp is not None else np.nan)
def p_from_q(q):
    # solve (1-p)^2+p^2/4=q -> 1.25p^2-2p+1-q=0
    return (2-np.sqrt(4-5*(1-q)))/2.5
def rw0_day(rng,pr,Lx,Ly,p):
    ids=pr.id.values;a=pr.a.values;b=pr.b.values;n=len(ids)
    x=rng.integers(0,Lx,n);y=rng.integers(0,Ly,n); rows=[]
    for s in range(a.min(),b.max()+1):
        act=(a<=s)&(s<=b); mv=(a<s)&(s<=b)
        d=rng.integers(0,4,n); dx=np.array([1,-1,0,0])[d]; dy=np.array([0,0,1,-1])[d]
        nx=x+dx;ny=y+dy; ok=mv&(rng.random(n)<p)&(nx>=0)&(nx<Lx)&(ny>=0)&(ny<Ly)
        x=np.where(ok,nx,x);y=np.where(ok,ny,y)
        site=x*Ly+y; ia=np.flatnonzero(act); o=ia[np.argsort(site[ia],kind='stable')]; ss=site[o]
        # groups of equal site
        brk=np.flatnonzero(np.r_[True,ss[1:]!=ss[:-1],True])
        for u in np.flatnonzero(np.diff(brk)>1):
            mem=np.sort(ids[o[brk[u]:brk[u+1]]])
            for q1 in range(len(mem)):
                for q2 in range(q1+1,len(mem)): rows.append((s,mem[q1],mem[q2]))
    return pd.DataFrame(rows,columns=['s','i','j'])
def sir_vec(S,I,J,N,elig,P,beta,tau,nrep,rng,maxloops=400):
    idx=rng.choice(elig,nrep); t0=rng.random(nrep)*P; R=np.arange(nrep)
    tinf=np.full((N,nrep),np.inf);trec=np.full((N,nrep),np.inf)
    tinf[idx,R]=t0;trec[idx,R]=t0+rng.exponential(tau,nrep)
    off=np.zeros(nrep,int);maxrec=trec[idx,R].max()
    o=np.argsort(S,kind='stable');S=S[o];I=I[o];J=J[o]
    for L in range(maxloops):
        base=L*P
        if base>maxrec:break
        for k in range(len(S)):
            T=base+S[k]
            if T>maxrec:break
            a=I[k];b=J[k]
            ia=(tinf[a]<T)&(trec[a]>T); ib=(tinf[b]<T)&(trec[b]>T)
            for (u,v,iu,iv) in ((a,b,ia,ib),(b,a,ib,ia)):
                c=iu&~iv&np.isinf(tinf[v])
                if c.any():
                    w=np.flatnonzero(c); w=w[rng.random(len(w))<beta]
                    if len(w):
                        tinf[v,w]=T; r=T+rng.exponential(tau,len(w)); trec[v,w]=r; maxrec=max(maxrec,r.max())
                        off[w]+= (idx[w]==u)
    fs=np.isfinite(tinf).sum(0)
    return fs,off
