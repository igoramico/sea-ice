#include <time.h>
#include "variables.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Der (double *h, int pos, int order) {

  double y;

  switch(order) {
   case 1: 
     y = (-h[(pos+2)] + 8.*h[(pos+1)] - 8.*h[(pos-1)] + h[(pos-2)])/12.;
     break;
   case 2:
     y = (-h[(pos+2)] + 16.*h[(pos+1)] - 30.*h[pos] + 16.*h[(pos-1)] - h[(pos-2)])/12.;
     break;
   case 3:
     y = (-h[(pos+3)] + 8.*h[(pos+2)] - 13.*h[(pos+1)] + 13.*h[(pos-1)] - 8.*h[(pos-2)] + h[(pos-3)])/8.;
     break;
   case 4:
     y = (-h[(pos+3)] + 12.*h[(pos+2)] - 39.*h[(pos+1)] + 56.*h[pos] - 39.*h[(pos-1)] + 12.*h[(pos-2)] - h[(pos-3)])/6.;
     break;
  }

  return y;

}

float ran1(long idum) {

    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    if (idum <= 0 || !iy) { 
	if (-(idum) < 1) idum=1; 
	else idum = -(idum);
      for (j=NTAB+7;j>=0;j--) { 
	  k=(idum)/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	if (j < NTAB) iv[j] = idum;
      }
      iy=iv[0];
    }
    k=(idum)/IQ; 
    idum=IA*(idum-k*IQ)-IR*k; 
    if (idum < 0) idum += IM; 
    j=iy/NDIV; 
    iy=iv[j]; 		      
    iv[j] = idum; 		      							       
    if ((temp=AM*iy) > RNMX) return RNMX; 	
    else return temp;
}

double gasdev(int idum) {  

  static int iset=0; 
  static float gset; 
  float fac,rsq,v1,v2; 

  if (idum < 0) iset=0; 
  if (iset == 0) {  
    do { 
      v1=2.0*ran1(idum)-1.0; 
      v2=2.0*ran1(idum)-1.0; 
      rsq=v1*v1+v2*v2; 
    } 
    while (rsq >= 1.0 || rsq == 0.0);  
    fac=sqrt(-2.0*log(rsq)/rsq); 
    gset=v1*fac; 
    iset=1;  
    return (double) v2*fac; 
  } else {  
    iset=0; 
    return (double) gset; 
  } 
}

void dump_config (int t, char name[], double *h) {

  FILE *fout;
  char filename[500];

  sprintf(filename,"%s/%s_config_t%d",OutDir,name,t);
  fout = fopen(filename,"wb");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  fwrite(h,sizeof(double),LP6,fout);
  fclose(fout);

}

void restore_config (int t, char name[], double *h) {

  FILE *fout;
  char filename[500];
  int status;

  sprintf(filename,"%s/%s_config_t%d",OutDir,name,t);
  fout = fopen(filename,"rb");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  status = fread(h,sizeof(double),LP6,fout);
  if(!status) { 
    fprintf(stderr,"Error! File %s could not be read! \n",filename);
   fprintf(stderr,"Restore of dumped configuration aborted!\n");
   exit(2);
  }

  fclose(fout);

}

void generate_stochastic (int idum, double *Br) {

  int i;

  for(i=0;i<LP6;i++) {
    Br[i] = gasdev(idum);
  }

}

void initialisation (double *g, double *phi) {

  int i,j,h;

  for(i=1;i<=L;i++) {
   for(j=1;j<=L;j++) {

     index = j + i*LP2;
     for(h=0;h<NHEIGHTS;h++) {
g[index].v[h] = 
     }

   }
  }

#ifdef INIT_SINE

  int i,kk;
  double nm;

  srand48(time(NULL));

  nm = (double)(KMAX-KMIN+1);

  for(i=3;i<=(L+2);i++) {
   h[i] = h0;
   for(kk=KMIN;kk<=KMAX;kk++) {
#ifdef INIT_SINE_RANDOM_PHASES
     h[i] += (dh0/nm)*sin(2.*M_PI*kk*((double)(i-3))/((double)(L))+M_PI*drand48());
#else
     h[i] += (dh0/nm)*sin(2.*M_PI*kk*((double)(i-3))/((double)(L)));
#endif
   }
  }

#endif

#ifdef INIT_RIM

  int i;
  double x,f1,f2;

  for(i=3;i<=(L+2);i++) {
  x = ((double)(i-3));
  f1 = hmin + 0.5*(hmax*(1.+exp(-(x-rim_x0)/rim_lambda)) - hmin)*(1.+tanh((x-rim_x0)/rim_delta1));
  f2 = hmin + 0.5*(hmax*(1.+exp((x+rim_x0)/rim_lambda)) - hmin)*(1.-tanh((x+rim_x0)/rim_delta1));
  h[i] = f2 + 0.5*(f1 - f2)*tanh(x/rim_delta0);
}

#endif

}

void bc (double *h) {

#ifdef BC_PERIODIC 
  h[0] = h[L];
  h[1] = h[(L+1)];
  h[2] = h[(L+2)];

  h[(L+3)] = h[3];
  h[(L+4)] = h[4];
  h[(L+5)] = h[5];
#endif

}

double Pi (double x, int order) {

  double y;

#ifdef DISJOINING 
  switch(order) {
   case 0: 
     y = -Ham1/(x*x*x);
     break;
   case 1:
     y = 3.*Ham1/(x*x*x*x);
     break;
   case 2:
     y = -12.*Ham1/(x*x*x*x*x);
     break;
  }
#endif

#ifdef DISJOINING_CONJOINING
  switch(order) {
   case 0: 
     y = -Ham1/(x*x*x) + Ham2/(x*x);
     break;
   case 1:
     y = 3.*Ham1/(x*x*x*x) - 2.*Ham2/(x*x*x);
     break;
   case 2:
     y = -12.*Ham1/(x*x*x*x*x) + 6.*Ham2/(x*x*x*x);
     break;
  }
#endif

  return y;

}

void deterministic_flux (int nnn, double *h, double *Jdet, int tstep) {

  int i;
  double gp;

  gravnow = grav*cos(OmegaGrav*tstep);

  for(i=3;i<=(L+2);i++) {

    gp = Pi(h[i],1)*Der(h,i,1) - surftens*Der(h,i,3);

    Jdet[i] = (1./(3.*mu))*h[i]*h[i]*h[i]*(gp + gravnow*Der(h,i,1));

  }

}

void fluctuating_flux (int nnn, double *h, double *Br, double *Jnoise) {

  int i,j;
  double invdt,fact;
  double Noise;
  double gp;
  double Temph;

  invdt = 1./(sqrt(dt));
  fact = sqrt(2./((double)L));

  fact *= invdt;

  for(i=3;i<=(L+2);i++) {

    /* Wiener term and derivative */
    Noise = 0.0;
    for(j=1;j<=L;j++) {
      Noise += sin(j*M_PI*((double)(i-3))/((double)L))*Br[j];
    }
    Noise *= fact;

    gp = Pi(h[i],1)*Der(h,i,1) - surftens*Der(h,i,3);

    Temph = sqrt(Temp)*(pow(h[i],((4.*nnn+3.)/2.)));

    Jnoise[i] = (pow(gp,(2.*nnn)))*Temph*Noise;

  }

}

void compute_rhs (double *Jdet, double *Jnoise, double *rhs) {

  int i;

  for(i=3;i<=(L+2);i++) {

    rhs[i] = Der(Jdet,i,1) + Der(Jnoise,i,1);

  }

}

void time_marching (double *h, double *rhs, double *rhsold) {

#ifdef ADAMS_BASHFORTH
  for(i=3;i<=(L+2);i++) {

    h[i] = h[i] + 1.5*dt*rhs[i] - 0.5*dt*rhsold[i];

  }
#endif

}

void total_volume (int t, double *h, FILE *fvol) {

  int i;
  double vol;

  vol = 0.0;
  for(i=3;i<=(L+2);i++) {
    vol += h[i];
  }
  fprintf(fvol," %d %g \n",t,vol);

}

void stabiliser (double *h) {

  int i;

  for(i=3;i<=(L+2);i++) {

    if(h[i]<cutoff) { h[i] = cutoff; }

  }

}

void print_f (int t, double *h) {

  FILE *fout;
  char filename[500];

  sprintf(filename,"%s/%s_t%d.dat",OutDir,basename,t);
  fout = fopen(filename,"w");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  
  for(i=3;i<=(L+2);i++) {

    fprintf(fout," %g %g \n",((double)(i-2)),h[i]);

  }

  fclose(fout);

}

void copy_array (double *v1, double *v2, int N) {

  int i;

  for(i=0;i<N;i++) {
    v1[i] = v2[i];
  }

}


