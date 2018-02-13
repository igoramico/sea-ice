int i,ip,im;
int L;
int iter,NTIME;

int printfreq,volfreq,dumpfreq;

int restore;
int start_time;

double dt;

double *h;
double *hrhs,*hrhsold;
double *sigma,*psi;
double *fR;

double S,sigma0;

double hmin;

int kk,KMIN,KMAX;
double h0,eps,pert;
double phi,epsf;

#ifdef MELT_PONDS
double *w;
double *wrhs,*wrhsold;
double *flux;
#endif

char basename[500],OutDir[500];
char icename[500],pondname[500];
