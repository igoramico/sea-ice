int i,ip,im;
int j,jp,jm;
int idx;
int LX,LY,TOTSIZE;
int iter,NTIME;

int printfreq,volfreq,dumpfreq;

int restore;
int start_time;

double dt;

double *h;
double *hrhs,*hrhsold;
double *psi;
double *sigmax,*sigmay;
double *fR;

double S,sigma0,psigma;

int RC;

double hmin;

int kk,KMIN,KMAX;
double h0,eps;

#ifdef MELT_PONDS
double *w;
double *wrhs,*wrhsold;
double *fluxx,*fluxy;
double *s;
double *tsx,*tsy;

double alpha1,alphad,kappa;
double ts0x,ts0y;

double mluethje,mpluethje,wmaxfr;
#endif

char basename[500],OutDir[500];
char icename[500],pondname[500];
char fname[500];

FILE *finp;
