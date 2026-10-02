"""Re-run the first analysis's frozen pipeline to check stored outputs are reproducible from frozen code."""
import sys,json,importlib,numpy as np
sys.path.insert(0,'../../bsc_validation2/a2_contact_mobility/src')
import a2lib as A, pipeline as PL
tag,n,models=sys.argv[1],sys.argv[2],sys.argv[3].split(',')
calfun=importlib.import_module(sys.argv[4]).calfun if len(sys.argv)>4 else None
ds=A.load(n)
res=PL.run_dataset(ds,models,G=5,G_epi=10,nrep_real=4000,nrep_model=400,seed=20261001,calfun=calfun,verbose=False)
PL.dump(res,f'rerun_{tag}_{n}.json')
old=json.load(open(f'../../bsc_validation2/a2_contact_mobility/out/{tag}_{n}.json')); new=json.load(open(f'rerun_{tag}_{n}.json'))
def cmp(a,b,path=''):
    d=0
    if isinstance(a,dict):
        for k in a:
            if k in('seconds',): continue
            d+=cmp(a[k],b.get(k),path+'/'+str(k))
    elif isinstance(a,list):
        for i,(x,y) in enumerate(zip(a,b)): d+=cmp(x,y,path+f'[{i}]')
    else:
        same = (a==b) or (isinstance(a,float) and isinstance(b,float) and (abs(a-b)<=1e-9*max(1,abs(a)) or (a!=a and b!=b)))
        if not same:
            d=1
            if cmp.n<8: print('DIFF',path,a,b); cmp.n+=1
    return d
cmp.n=0
print(tag,n,'n_diffs',cmp(old,new))
