/* SIR on a periodic (looped) temporal contact list. Exact per-contact Bernoulli transmission,
   exponential infectious period. Time unit = 20 s interval. */
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
static inline uint64_t nxt(uint64_t *s){ uint64_t z=(*s+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
static inline double u01(uint64_t *s){ return ((nxt(s)>>11)+0.5)*(1.0/9007199254740992.0);}
void sir_loop(int nc,const int *s,const int *ci,const int *cj,int N,int nel,const int *elig,
              long P,double beta,double tau,int nrep,uint64_t seed,int maxloops,int *finalsize,int *idxoff){
  double *tinf=malloc(sizeof(double)*N),*trec=malloc(sizeof(double)*N);
  uint64_t st=seed;
  for(int r=0;r<nrep;r++){
    for(int k=0;k<N;k++){tinf[k]=INFINITY;trec[k]=INFINITY;}
    int idx=elig[(int)(u01(&st)*nel)]; double t0=u01(&st)*P;
    tinf[idx]=t0; trec[idx]=t0-tau*log(u01(&st)); double maxrec=trec[idx];
    int ninf=1,off=0,done=0;
    int lo=0,hi=nc; while(lo<hi){int m=(lo+hi)/2; if(s[m]>t0) hi=m; else lo=m+1;}
    int k0=lo;
    for(int L=0;L<maxloops && !done;L++){
      double base=(double)L*P;
      for(int k=(L==0?k0:0);k<nc;k++){
        double T=base+s[k];
        if(T>maxrec){done=1;break;}
        int a=ci[k],b=cj[k];
        int ia=(tinf[a]<T && trec[a]>T), ib=(tinf[b]<T && trec[b]>T);
        if(ia==ib) continue;
        if(ib){int t=a;a=b;b=t;}
        if(tinf[b]!=INFINITY) continue;
        if(u01(&st)<beta){ tinf[b]=T; trec[b]=T-tau*log(u01(&st)); if(trec[b]>maxrec)maxrec=trec[b]; ninf++; if(a==idx)off++; }
      }
      if(base+P>maxrec) done=1;
    }
    finalsize[r]=ninf; idxoff[r]=off;
  }
  free(tinf);free(trec);
}
