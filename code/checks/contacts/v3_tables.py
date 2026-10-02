import json,numpy as np,sys
sys.path.insert(0,'../../bsc_validation2/a2_contact_mobility/src'); import a2lib as A
DS=['InVS13','LyonSchool','LH10','InVS15','Thiers13','SFHH']
tags={'base':('s1','s3a'),'RWclus':('s2b','s3b'),'RWstickZ':('s2d','s3c'),'RWstickO':('s2d','s3d')}
def get(m,n):
    t=tags.get(m,tags['base'])[0 if n in A.DEV else 1]; return json.load(open(f'../../bsc_validation2/a2_contact_mobility/out/{t}_{n}.json'))
FAM=list(A.FAMILIES)
allr={}
for m in ['WM','BLOCK','RW0','RWhet','RWclus','RWstickZ','RWstickO']:
    sf=[];ep=[];rr=[];da=[];dp=[];hi=0;cells=0;lr=[];adA=[];adP=[];tr=[]
    fam={f:0 for f in FAM}
    for n in DS:
        r=get(m,n); sc=r['models'][m]['score']
        # recompute score independently from stats
        st=r['models'][m]['stats']; ob=r['obs']
        def mp(k):
            a,b=st.get(k),ob.get(k)
            if a is None or b is None or a!=a or b!=b: return None
            if k=='S2_burst': return abs(a-b)<=0.1
            ok=0.8<=a/b<=1.25
            if k in A.FRAC: ok=ok or abs(a-b)<=0.03
            return ok
        myS=0
        for f,ks in A.FAMILIES.items():
            ps=[mp(k) for k in ks]; ps=[p for p in ps if p is not None]
            v=all(ps) if ps else None
            assert v==sc[f],(m,n,f)
            if v: myS+=1; fam[f]+=1
        sf.append(myS); npass=0; rat=[];a_=[];p_=[]
        for k,row in r['epi'].items():
            e,o=row[m],row['REAL']; q=e['R_index']/o['R_index']
            ok=0.9<=q<=1/0.9 and abs(e['P_major']-o['P_major'])<=0.05 and abs(e['attack']-o['attack'])<=0.05
            assert ok==e['pass']; npass+=ok; rat.append(q); a_.append(e['attack']-o['attack']); p_.append(e['P_major']-o['P_major'])
            hi+=q>1; cells+=1; lr.append(abs(np.log(q))); adA.append(abs(a_[-1])); adP.append(abs(p_[-1]))
        ep.append(npass); rr.append(np.mean(rat)); da.append(np.mean(a_)); dp.append(np.mean(p_)); tr.append(r['models'][m]['epi_total_ratio'])
        allr[(m,n)]=(min(rat),max(rat))
    print(m,'S',sf,'E',ep,'Rratio',np.round(rr,2),'dAtt',np.round(da,3),'dP',np.round(dp,3),f'hi {hi}/{cells}','|logR| %.3f |dA| %.3f |dP| %.3f'%(np.mean(lr),np.mean(adA),np.mean(adP)),'totratio',np.round(tr,2))
    print('    fam',fam, 'cell range R ratio',round(min(v[0] for k,v in allr.items() if k[0]==m),2),round(max(v[1] for k,v in allr.items() if k[0]==m),2))
# REALtrain
print('REALtrain',{n:[f for f in FAM if get('base',n)['models']['REALtrain']['score'][f] is False] for n in DS})
# real identical across tags
for n in A.CONF:
    rs=[json.load(open(f'../../bsc_validation2/a2_contact_mobility/out/s3{t}_{n}.json')) for t in 'abcd']
    print(n,'real identical',all(json.dumps({k:v['REAL'] for k,v in r['epi'].items()})==json.dumps({k:v['REAL'] for k,v in rs[0]['epi'].items()}) for r in rs))
# S7
for n in DS:
    r=get('base',n); print(n,'S7 obs',round(r['obs']['S7_alpha'],2),'range',round(r['obs']['S7_n_range'],2),'RW0',round(r['models']['RW0']['stats']['S7_alpha'],2),'WM',round(r['models']['WM']['stats']['S7_alpha'],2))
# fits
for n in DS:
    for m in ['RWclus','RWstickZ','RWstickO']:
        e=get(m,n)['extra'][m]; print(n,m,{k:(round(v,4) if isinstance(v,float) else v) for k,v in e.items() if k!='trace'})
