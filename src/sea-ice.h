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

  double r,epsh0;
  
#ifdef INIT_ICE_RANDOM

  epsh0 = eps*h0;

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
      idx = j + i*LY;
      r = (2.0*drand48() - 1.0);
      h[idx] = h0 + epsh0*r;
     }
    }
#endif

#ifdef INIT_ICE_FOURIER

    double phi;
    double pert,epsf;
    int NK;
    double ff=7./8.;
    
    epsf  = eps*M_PI;
    epsh0 = eps*h0;

    NK = KMAX-KMIN+1;
    
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       idx = j + i*LY;
      pert = 0.0;
      for(kk=KMIN;kk<=KMAX;kk++) {
	r = drand48();
	phi = pow(-1.0,kk)*ff*(kk-KMIN)*M_PI/NK;
	pert += sin((2.0*(double)kk*M_PI*(double)i/((double)LX))+phi)*sin((2.0*(double)kk*M_PI*(double)j/((double)LY))+phi);
      }
      h[idx] = h0 + epsh0*pert/sqrt(NK);
     }
    }
#endif

#ifdef INIT_ICE_FROM_INPUT

    fin = fopen("init_ice_topography.inp","r");
    if(fin==NULL) {
      fprintf(stderr,"Error! Initial topography file not found!! \n");
      exit(2);
    } else {
      for(i=0;i<LX;i++) {
       for(j=0;j<LY;j++) {
	 idx = j + i*LY;
	fscanf(fin,"%lf",&h[idx]);
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

  /*** to be implemented ***/

}
#else
void compute_fR (double *h, double *fR) {

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
    if(h[idx]>hmin) {
     fR[idx] = S/h[idx]; 
    } else {
      fR[idx]=0.0;
    }
  }

}
#endif

#ifdef MELT_PONDS
void compute_sigma (double *h, double *w, double *sigma) {

  /*** to be implemented ***/
}
#else
void compute_sigma (double *h, double *sigma) {

#ifdef MECHANICAL_CORRELATED
  int ii,jj;
  int irc,jrc,idxrc;
  double xh,yh,xt,yt;
  double dist;

  for(i=0;i<LX;i++) {
    xh = (double)i;
   for(j=0;j<LY;j++) {
     yh = (double)j;
     idx = j + i*LY;
     sigma[i]=0.0;
    for(ii=-RC;ii<=RC;ii++) {
      xt = (double)ii;
    for(jj=-RC;jj<=RC;jj++) {
      yt = (double)jj;
      dist = sqrt((xh-xt)*(xh-xt)+(yh-yt)*(yh-yt));
      if(dist<=((double)RC)&&dist!=0.0) {
	irc = (i + ii + LX)%LX;
	jrc = (j + jj + LY)%LY;
	idxrc = jrc + irc*LY;
        sigma[idx] += sigma0*drand48()*h[idx]*h[idxrc]/pow(dist,psigma);
      }
     }
    }
   }
  }
#else

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
    idx = j + i*LY;
    sigma[idx] = sigma0*gasdev();
   }
  }
#endif
  
}
#endif


void compute_hrhs (double *fR, double *psi, double *hrhs) {


  for(i=0;i<LX;i++) {
  for(j=0;j<LY;j++) {
    idx = j +i*LY;
    hrhs[idx] = -fR[idx] + psi[idx];

  }
  }

}

#ifdef MELT_PONDS
void compute_wrhs (double *fR, double *flux, double *s, double *wrhs) {


  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
    wrhs[idx] = fR[idx] - flux[idx] + s[idx];

  }
  }

}
#endif

void compute_psi (double *sigma, double *psi) {

  int idxp,idxm;
  double dsxdx,dsydy;

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
    if(h[idx]>hmin) {     
     ip = (i + 1 + LX)%LX;
     im = (i - 1 + LX)%LX;
     idxp = j + ip*LY;
     idxm = j + im*LY;
     dsxdx = 0.5 * (sigma[idxp].x - sigma[idxm].x);

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxp = jp + i*LY;
     idxm = jm + i*LY;
     dsydy = 0.5 * (sigma[idxp].y - sigma[idxm].y);


     psi[idx] = dsxdx + dsydy;
    } else {
      psi[idx] = 0.0;
    }
    
   }
  }
  
}

#ifdef MELT_PONDS
void compute_flux (double *h, double *w, double *flux) {
  /*** to be implemented ***/
}  
#endif

void time_marching (double *h, double *rhs, double *rhsold) {

#ifdef ADAMS_BASHFORTH
  for(i=0;i<L;i++) {

    h[i] = h[i] + 1.5*dt*hrhs[i] - 0.5*dt*hrhsold[i];

  }
#endif

#ifdef CRANK_NICHOLSON
  /*** to be implemented ***/
#endif 

}

void total_volume (int t, double *h, FILE *fvol) {

  double vol;

  vol = 0.0;
  for(i=0;i<LX;i++) {
   for(j=0;j<LX;j++) {
     idx = j + i*LY;
    vol += h[idx];
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
  
  for(i=0;i<LX;i++) {
  for(j=0;j<LY;j++) {
    idx = j + i*LY;
    fprintf(fout," %g %g \n",((double)i),((double)j),h[idx]);

  }
  fprintf(fout,"\n");
  }
  fclose(fout);

}

void copy_array (double *v1, double *v2, int N) {

  for(i=0;i<N;i++) {
    v1[i] = v2[i];
  }

}


