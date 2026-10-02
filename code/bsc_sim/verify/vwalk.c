#include <math.h>
#include <stdint.h>
static uint64_t st;
static inline uint64_t nx(void){uint64_t z=(st+=0x9E3779B97F4A7C15ULL);z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;return z^(z>>31);}
static inline double U(void){return (nx()>>11)*(1.0/9007199254740992.0);}
/* n walkers killed at rate gamma; accumulates bq*dt ; start cells given */
void walk(int M,const int*nbr,const double*w,const double*wtot,const double*bq,double gamma,int n,const int*start,uint64_t seed,double*wsum,double*cell){
    st=seed;
    for(int i=0;i<n;i++){
        int x=start[i]; double life=-log(1.0-U())/gamma, acc=0.0;
        while(1){
            double hold = wtot[x]>0 ? -log(1.0-U())/wtot[x] : 1e300;
            if(hold>=life){acc+=bq[x]*life; cell[x]+=bq[x]*life; break;}
            acc+=bq[x]*hold; cell[x]+=bq[x]*hold; life-=hold;
            double r=U()*wtot[x],a=0.0; int y=x;
            for(int j=0;j<4;j++){ if(nbr[4*x+j]<0) continue; y=nbr[4*x+j]; a+=w[4*x+j]; if(a>r) break; }
            x=y;
        }
        wsum[i]=acc;
    }
}
