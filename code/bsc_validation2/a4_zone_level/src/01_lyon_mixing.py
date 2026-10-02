"""Mixing weights from the SocioPatterns Lyon primary-school contact data (independent of any outbreak).
Students only (teachers dropped). Contact time = 20 s per record. Output: data/derived/lyon_mixing.json"""
import pandas as pd, numpy as np, json
c=pd.read_csv('data/raw/sociopatterns/primaryschool.csv.gz',sep='\t',header=None,names=['t','i','j','ci','cj'])
m=pd.read_csv('data/raw/sociopatterns/primaryschool_metadata.txt',sep='\t',header=None,names=['id','cls','sex'])
stu=m[m.cls!='Teachers']; ncls=stu.groupby('cls').size()
c=c[(c.ci!='Teachers')&(c.cj!='Teachers')].copy()
c['day']=(c.t>100000).astype(int)+1
def stratum(a,b): return 'class' if a==b else ('grade' if a[0]==b[0] else 'school')
c['str']=[stratum(a,b) for a,b in zip(c.ci,c.cj)]
N=len(stu); out={'n_students':int(N),'class_sizes':ncls.to_dict(),'t_range':[int(c.t.min()),int(c.t.max())]}
def weights(df,ndays):
    # person-seconds: each record contributes 20 s to each of the two students
    tot=df.groupby('str').size()*20*2/N/ndays   # mean contact seconds per student per day, by stratum
    return tot
w=weights(c,2); out['sec_per_student_day']=w.to_dict(); out['shares']=(w/w.sum()).to_dict()
for d in (1,2):
    wd=weights(c[c.day==d],1); out[f'shares_day{d}']=(wd/wd.sum()).to_dict()
# per-pair mean seconds/day
npairs={'class':sum(n*(n-1)/2 for n in ncls),'grade':sum(ncls[g+'A']*ncls[g+'B'] for g in '12345')}
npairs['school']=N*(N-1)/2-npairs['class']-npairs['grade']
out['per_pair_sec_day']={k:float(c[c.str==k].shape[0]*20/2/npairs[k]) for k in npairs}
# bootstrap over students is awkward; give per-class variation of the within-class share instead
sh=[]
for k in ncls.index:
    a=c[(c.ci==k)|(c.cj==k)]
    own=2*(a.str=='class').sum(); g=(a.str=='grade').sum(); s=(a.str=='school').sum()
    sh.append([own/(own+g+s),g/(own+g+s),s/(own+g+s)])
sh=np.array(sh); out['per_class_shares_min']=sh.min(0).tolist(); out['per_class_shares_max']=sh.max(0).tolist()
# number-of-contact-events alternative is not computed (duration is the quantity of the model)
json.dump(out,open('data/derived/lyon_mixing.json','w'),indent=1); print(json.dumps(out,indent=1))
