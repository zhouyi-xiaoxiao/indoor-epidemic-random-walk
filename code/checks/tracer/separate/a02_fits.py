import json,numpy as np,pandas as pd
from advlib import *
res={}
# ---- cabin
km=json.load(open('kinahan_fwd_means.json'))   # values are half (double block) -> scale irrelevant for log fit; zeros stay zeros
cab=Lat(6,46,0.9,1.1); ci={'A':0,'D':2,'G':3,'L':5}
def cabfit(rel,drop=(),kap=35.0,zero_missing=False):
    seats=[s for s,v in km[rel].items() if v is not None and v>0 and s!=rel and s not in drop]
    rec=[cab.site(ci[s[-1]],int(s[:-1])-1+2) for s in seats]; c=[km[rel][s] for s in seats]
    src=cab.site(ci[rel[-1]],int(rel[:-1])-1+2)
    D,sse=fitD(cab,src,rec,c,kap); return D,len(seats),sse.min(),sse[-1]
for rel in('5L','5A','5G'):
    print('cabin',rel,cabfit(rel),'| without extra sensors 5A/5G/5L/7A:',cabfit(rel,drop=('5A','5G','5L','7A'))[:2],'| kap 30:',cabfit(rel,kap=30)[0])
# ---- coach
b1=pd.read_csv('../../../bsc_validation2/a3_tracer_physics/data/derived/ou2022_B1_seats.csv')
coach=Lat(5,13,0.5,11.4/13); cc={'A':0,'B':1,'E':2,'C':3,'D':4}
meas={'1D':.70,'5C':.76,'6B':.82,'9C':.92,'11B':.99,'11C':1.03,'13A':.72,'13D':1.00}   # typed again for the re-check from rendered Fig S3A
rec=[coach.site(cc[s[-1]],int(s[:-1])-1) for s in meas]
D,sse=fitD(coach,coach.site(cc['C'],11),rec,list(meas.values()),6.5)
print('coach D',D,'sse min',sse.min(),'sse wellmixed',sse[-1],'sse at D=25',sse[np.argmin(abs(DGRID-25.1))],'D range within +10% of min sse',DGRID[sse<=1.1*sse.min()][[0,-1]])
# incl 12D reference seat
m2=dict(meas); m2['12D']=1.0
rec2=[coach.site(cc[s[-1]],int(s[:-1])-1) for s in m2]
print('coach incl 12D',fitD(coach,coach.site(3,11),rec2,list(m2.values()),6.5)[0])
# ---- train
pts=pd.read_csv('../../../bsc_validation2/a3_tracer_physics/data/derived/woodward2022_fig9_points.csv'); print(pts.head(3).to_string()); print(pts.groupby(['panel','near_source']).size())
tr=Lat(6,20,2.69/6,19.85/20)
for panel,xs in(('B_middle',-0.135),('A_end',-0.392)):
    d=pts[(pts.panel==panel)&(~pts.near_source)]
    ys=np.clip(np.floor((d.x_over_L.values+.5)*20).astype(int),0,19); xsn=np.where(d.colour=='red',1,4)
    src=tr.site(1,int(np.floor((xs+.5)*20)))
    for kap in(11,13,15):
        D,sse=fitD(tr,src,tr.site(xsn,ys),d.c_norm.values,kap); print('train',panel,'kap',kap,'D',D,'ell',np.sqrt(D/kap))
