#include <time.h>
#include "variables.h"

double gasdev() {

  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;

  if (iset == 0) {
    do {
      v1=2.0*drand48()-1.0;
      v2=2.0*drand48()-1.0;
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

void dump_config (int t, char basename[], double *h) {

  FILE *fout;
  char filename[500];

  sprintf(filename,"%s/%s_config_t%d",OutDir,basename,t);
  fout = fopen(filename,"wb");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  fwrite(h,sizeof(double),L,fout);
  fclose(fout);
  
}

void restore_config (int t, char basename[], double *h) {

  FILE *fout;
  char filename[500];
  int status;

  sprintf(filename,"%s/%s_config_t%d",OutDir,basename,t);
  fout = fopen(filename,"rb");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  status = fread(h,sizeof(double),L,fout);
  if(!status) { 
    fprintf(stderr,"Error! File %s could not be read! \n",filename);
   fprintf(stderr,"Restore of dumped configuration aborted!\n");
   exit(2);
  }

  fclose(fout);

}

#ifdef MELT_PONDS
void initialisation (double *h, double *w) {


}
#else
void initialisation (double *h) {
  
#ifdef INIT_ICE_RANDOM
    for(i=0;i<L;i++) {
      /* to be continued */
    }
#endif

#ifdef INIT_ICE_FOURIER
    for(i=0;i<L;i++) {
      pert = 0.0;
      for(kk=KMIN;kk<=KMAX;kk++) {
	r = drand48();
	phi = epsf*(2.0*r-1.0);
	pert += sin(2.0*kk*pi/L+phi);
      }
      h[i] = h0 + eps*pert;
    }
#endif

#ifdef INIT_ICE_FROM_INPUT
    fin = fopen("init_ice_topography.inp","r");
    if(fin==NULL) {
      fprintf(stderr,"Error! Initial topography file not found!! \n");
      exit(2);
    } else {
      for(i=0;i<L;i++) {
	fscanf(fin,"%lf",&h[i]);
      }
      fclose(fin);
    }
#endif

}

#endif


/*
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
*/

#ifdef MELT_PONDS
void compute_fR (double *h, double *w, double *fR) {

  for(i=0;i<L;i++) {
    if(h[j]>hmin) {
     fR[i] = S/h[i]; 
    } else {
      fR[i]=0.0;
    }
  }

}
#else
void compute_fR (double *h, double *fR) {

  for(i=0;i<L;i++) {
    if(h[i]>hmin) {
     fR[i] = S/h[i]; 
    } else {
      fR[i]=0.0;
    }
  }

}
#endif

#ifdef MELT_PONDS
void compute_sigma (double *h, double *w, double *sigma) {

  for(i=0;i<L;i++) {
  sigma[i] = sigma0*gasdev();
  }
  
}
#else
void compute_sigma (double *h, double *sigma) {


  for(i=0;i<L;i++) {
  sigma[i] = sigma0*gasdev();
  }
  
}
#endif


void compute_hrhs (double *fR, double *psi, double *hrhs) {


  for(i=0;i<L;i++) {

    hrhs[i] = -fR[i] + psi[i];

  }

}

#ifdef MELT_PONDS
void compute_wrhs (double *fR, double *flux, double *s, double *wrhs) {


  for(i=0;i<L;i++) {

    wrhs[i] = fR[i] - flux[i] + s[i];

  }

}
#endif

void compute_psi (double *sigma, double *psi) {

  for(i=0;i<L;i++) {

    ip = (i + 1 + L)%L;
    im = (i - 1 + L)%L;

    psi[i] = 0.5 * (sigma[ip] - sigma[im]); 

  }
  
}

#ifdef MELT_PONDS
void compute_flux (double *h, double *w, double *flux) {

}  
#endif

void time_marching (double *h, double *rhs, double *rhsold) {

#ifdef ADAMS_BASHFORTH
  for(i=0;i<L;i++) {

    h[i] = h[i] + 1.5*dt*hrhs[i] - 0.5*dt*hrhsold[i];

  }
#endif

}

void total_volume (int t, double *h, FILE *fvol) {

  double vol;

  vol = 0.0;
  for(i=0;i<L;i++) {
    vol += h[i];
  }
  fprintf(fvol," %d %g \n",t,vol);

}

/*
void stabiliser (double *h) {

  int i;

  for(i=3;i<=(L+2);i++) {

    if(h[i]<cutoff) { h[i] = cutoff; }

  }

}
*/
 
void print_f (int t, double *h, char basename[]) {

  FILE *fout;
  char filename[500];

  sprintf(filename,"%s/%s_t%d.dat",OutDir,basename,t);
  fout = fopen(filename,"w");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  
  for(i=0;i<L;i++) {

    fprintf(fout," %g %g \n",((double)i),h[i]);

  }

  fclose(fout);

}

void copy_array (double *v1, double *v2, int N) {

  for(i=0;i<N;i++) {
    v1[i] = v2[i];
  }

}


