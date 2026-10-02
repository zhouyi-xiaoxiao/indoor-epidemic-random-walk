import h5py, numpy as np, pandas as pd
f=h5py.File('raw/m.jld2','r'); d=f['anonymizedstudents']; refs=d[:]
o=f[refs[0]]; tp=o.id.get_type()
print('compound size',tp.get_size(),'members',[(tp.get_member_name(i),tp.get_member_offset(i),tp.get_member_type(i).get_size()) for i in range(tp.get_nmembers())])
names=[tp.get_member_name(i).decode() for i in range(tp.get_nmembers())]; offs=[tp.get_member_offset(i) for i in range(tp.get_nmembers())]
rows=[]
for r in refs:
    o=f[r]; raw=np.empty((),dtype='V%d'%tp.get_size()); o.id.read(h5py.h5s.ALL,h5py.h5s.ALL,raw,mtype=o.id.get_type()); b=raw.tobytes()
    rec={}
    for n,of in zip(names,offs):
        if n=='isinfected': rec[n]=b[of]
        elif n in('onset','schoolID','gradeID','classID','sex'): rec[n]=int.from_bytes(b[of:of+8],'little',signed=True)
        elif n in('classsize','nclasses','householdlikelihoodratio','sampleweight'): rec[n]=np.frombuffer(b[of:of+8],'<f8')[0]
    rows.append(rec)
df=pd.DataFrame(rows); df.to_csv('v_students.csv',index=False)
a=pd.read_csv('../../bsc_validation2/a4_zone_level/data/derived/matsumoto_students.csv')
print(len(df),len(a),'identical core cols:',all((df[c].values==a[c].values).all() for c in ['isinfected','onset','schoolID','gradeID','classID']))
print('infected',df.isinfected.sum(),'schools',df.schoolID.nunique(),'classes',df.groupby(['schoolID','gradeID','classID']).ngroups)
g=df.groupby(['schoolID','gradeID','classID']).agg(n=('onset','size'),cs=('classsize','first'),ncl=('nclasses','first'),csn=('classsize','nunique')).reset_index()
nc=g.groupby(['schoolID','gradeID']).classID.transform('size')
print('nclasses field == classes seen:',(nc==g.ncl).mean(),' classsize>=respondents:',(g.cs>=g.n).mean(),'resp/classsize overall',g.n.sum()/g.cs.sum(), 'sum classsize',g.cs.sum())
print('single-class grades: classes',(nc==1).sum(),'students',g.n[nc==1].sum())
inf=df[df.isinfected==1]; print('onset range infected',inf.onset.min(),inf.onset.max(),'uninfected onset values',df[df.isinfected==0].onset.unique()[:5])
print('onset in break 89-100:',((inf.onset>=89)&(inf.onset<=100)).sum())
