int i,ip,im;
int j,jp,jm;
int idx;
int LX,LY,TOTSIZE;
int iter,NTIME;

int printfreq,volfreq,dumpfreq;

int restore;
int start_time;

double dt,dx;

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
double *Flux;
double *s;
double *tsx,*tsy;

double s0;
double wmin,wmin_melt,wminss;

int freq_super_seeper;

double alpha1,alphad,kappa;
double ts0x,ts0y,tsmod;

double mluethje,mpluethje,wmaxfr;
double ws_luethje;

double WWSmax;

#endif

#ifdef TEMPERATURE_FIELD

double *t;

#endif

char basename[500],OutDir[500];
char icename[500],pondname[500];
char fname[500];

FILE *finp;
