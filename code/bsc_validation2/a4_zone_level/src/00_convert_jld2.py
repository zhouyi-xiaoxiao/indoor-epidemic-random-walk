"""Convert Endo et al. anonymizedstudents.jld2 (HDF5/JLD2) to CSV without Julia.
Reads each referenced 89-byte compound record with h5py low-level API.
Prints STRUCTURE ONLY (no class-level outcome tabulation) so that the pre-registration stays model-blind."""
import sys, struct, hashlib
sys.path.insert(0, 'vendor')
import h5py, numpy as np, pandas as pd
P='data/raw/matsumoto/data_anonymizedstudents.jld2'
f=h5py.File(P,'r'); d=f['anonymizedstudents']; refs=d[:]
rows=[]
for r in refs:
    o=f[r]; raw=np.empty((),dtype='V89'); o.id.read(h5py.h5s.ALL,h5py.h5s.ALL,raw,mtype=o.id.get_type())
    b=raw.tobytes()
    isinf=b[0]; onset,school,grade,cls,sex=struct.unpack('<5q',b[1:41])
    classsize,nclasses,hhlr,sw=struct.unpack('<4d',b[57:89])
    rows.append((isinf,onset,school,grade,cls,sex,classsize,nclasses,hhlr,sw))
df=pd.DataFrame(rows,columns='isinfected onset schoolID gradeID classID sex classsize nclasses hh_loglr sampleweight'.split())
df.to_csv('data/derived/matsumoto_students.csv',index=False)
print('n',len(df),'sha256 raw',hashlib.sha256(open(P,'rb').read()).hexdigest())
print('schools',df.schoolID.nunique(),'grades',sorted(df.gradeID.unique()),'classIDs',sorted(df.classID.unique()))
g=df.groupby(['schoolID','gradeID','classID']).size()
print('n classes',len(g),'class size (respondents) quantiles',g.quantile([0,.1,.5,.9,1]).tolist())
print('reported classsize quantiles',df.classsize.quantile([0,.1,.5,.9,1]).tolist())
print('nclasses values',sorted(df.nclasses.unique()))
print('students per school', df.groupby('schoolID').size().tolist())
print('sex',df.sex.value_counts().to_dict(),'sampleweight range',df.sampleweight.min(),df.sampleweight.max())
print('onset coding: min/max over all',df.onset.min(),df.onset.max(),' n infected total',int(df.isinfected.sum()))
