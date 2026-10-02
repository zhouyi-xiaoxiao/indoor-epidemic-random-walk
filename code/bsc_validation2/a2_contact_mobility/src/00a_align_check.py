"""Format check only: are f2f and co-presence clocks/IDs aligned? No model, no contact statistics."""
import numpy as np, pandas as pd, gzip, zipfile, io, json
R='data/raw/'; C='data/copres/co-presence/'
def f2f(name):
    if name=='InVS13':
        z=zipfile.ZipFile(R+'workplace_InVS_tij.dat.zip'); return pd.read_csv(io.BytesIO(z.read('tij_InVS.dat')),sep=r'\s+',header=None,usecols=[0,1,2]).values
    fn={'InVS15':'workplace_InVS15_tij.dat.gz','SFHH':'SFHH_tij.dat.gz','LH10':'hospital_lyon_contacts.dat.gz','LyonSchool':'primaryschool.csv.gz','Thiers13':'HighSchool2013_proximity_net.csv.gz'}[name]
    return pd.read_csv(R+fn,sep=r'\s+',header=None,usecols=[0,1,2]).values
out={}
for name in ['InVS13','InVS15','LH10','LyonSchool','SFHH','Thiers13']:
    a=f2f(name); c=pd.read_csv(C+f'tij_pres_{name}.dat',sep=r'\s+',header=None).values
    ia=set(np.unique(a[:,1:])); ic=set(np.unique(c[:,1:]))
    key=lambda t,i,j:(t.astype(np.int64)<<24)+(np.minimum(i,j).astype(np.int64)<<12)+np.maximum(i,j)
    kc=np.unique(key(c[:,0],c[:,1],c[:,2]))
    t0a,t0c=a[:,0].min(),c[:,0].min()
    res={}
    for d in sorted(set([0,-20,20,int(t0c-t0a),int(t0c-t0a)-20,int(t0c-t0a)+20, int(c[:,0].max()-a[:,0].max())])):
        ka=key(a[:,0]+d,a[:,1],a[:,2]); res[d]=float(np.isin(ka,kc).mean())
    out[name]=dict(n_f2f=len(a),n_cop=len(c),ids_f2f=len(ia),ids_cop=len(ic),ids_common=len(ia&ic),t_f2f=[int(a[:,0].min()),int(a[:,0].max())],t_cop=[int(c[:,0].min()),int(c[:,0].max())],frac_f2f_in_cop_by_offset=res)
    print(name,out[name],flush=True)
json.dump(out,open('out/00a_align_check.json','w'),indent=1)
