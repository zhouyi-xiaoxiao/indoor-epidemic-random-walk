import pandas as pd, numpy as np, json
c=pd.read_csv('raw/primaryschool.csv.gz',sep='\t',header=None,names=['t','i','j','ci','cj'])
m=pd.read_csv('raw/meta.txt',sep='\t',header=None,names=['id','cls','sex'])
print(m.cls.value_counts().to_dict()); print('t range',c.t.min(),c.t.max(), 'gap check', np.sort(c.t.unique())[np.argmax(np.diff(np.sort(c.t.unique())))])
c=c[(c.ci!='Teachers')&(c.cj!='Teachers')]
lev=np.where(c.ci==c.cj,'class',np.where(c.ci.str[0]==c.cj.str[0],'grade','school')); c=c.assign(lev=lev,day=np.where(c.t>90000,2,1))
d=c.groupby('lev').size(); print('duration shares',(d/d.sum()).round(4).to_dict(), 'sec/student/day',(d*40/232/2).round(0).to_dict())
for dd in (1,2):
    x=c[c.day==dd].groupby('lev').size(); print('day',dd,(x/x.sum()).round(4).to_dict())
# alternative contact metrics (forking-path sensitivity): distinct pairs per day; contact events
pairs=c.drop_duplicates(['day','i','j']).groupby('lev').size(); print('distinct-pair-day shares',(pairs/pairs.sum()).round(4).to_dict())
# pairs with >= 1 min cumulative / day
pd_=c.groupby(['day','i','j','lev']).size().reset_index(name='k'); 
for thr in (3,15):
    y=pd_[pd_.k>=thr].groupby('lev').size(); print('pairs with >=%ds/day shares'%(thr*20),(y/y.sum()).round(4).to_dict())
# excluding lunch/break? time-of-day: class hours only is not defined in data; skip
