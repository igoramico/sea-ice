int i,L;
int time,NTIME;
int nnn;

int printfreq,volfreq,dumpfreq;

int restore;
int start_time;

double dt;

double *h;
double *Jdet,*Jnoise;
double *rhs,*rhsold;
double *Br;

double surftens;
double Temp;
double cutoff;

#ifdef DISJOINING
double Ham1;
#endif

#ifdef DISJOINING_CONJOINING
double Ham1,Ham2;
#endif 

#ifdef INIT_SINE
int KMIN,KMAX;
double h0,dh0;
#endif

#ifdef INIT_RIM
double hmin,hmax;
double rim_delta0,rim_delta1,rim_lambda,rim_x0;
#endif

char basename[500],OutDir[500];
