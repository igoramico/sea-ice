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
  fwrite(h,sizeof(double),(LX*LY),fout);
  fclose(fout);
  
}

void restore_config (int t, char basename[], double *h) {

  FILE *fout;
  char filename[500];
  int status;

  sprintf(filename,"%s/%s_config_t%d",OutDir,basename,t);
  fout = fopen(filename,"rb");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  status = fread(h,sizeof(double),(LX*LY),fout);
  if(!status) { 
    fprintf(stderr,"Error! File %s could not be read! \n",filename);
   fprintf(stderr,"Restore of dumped configuration aborted!\n");
   exit(2);
  }

  fclose(fout);

}

#ifdef MELT_PONDS
void initialisation (double *h, double *w) {

  double r,epsh0;

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
      idx = j + i*LY;
      h[idx] = h0;
      w[idx] = 0.0;
     }
    }

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
    }
#endif

}
#else
void initialisation (double *h) {

  double r,epsh0;

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
      idx = j + i*LY;
      h[idx] = h0;
     }
    }
    
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

#ifdef FR_LUETHJE_ET_AL
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     if(h[idx]>hmin) {
      if(w[idx]<wmaxfr) {
	fR[idx] = (1.0 + mpluethje/mluethje*w[idx]/wmaxfr)*mluethje; 
       } else {
	fR[idx] = (1.0 + mpluethje/mluethje)*mluethje; 
      }
     } else {
      fR[idx]=0.0;
     }
   }
  }
#endif

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
}  
#endif

/* 
#ifdef MELT_PONDS
void compute_sigma (double *h, double *w, double *sigmax, double *sigmay) {

  to be implemented
}
#else
*/

void compute_sigma (double *h, double *sigmax, double *sigmay) {

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
     sigmax[idx]=sigmay[idx]=0.0;
    for(ii=-RC;ii<=RC;ii++) {
      xt = (double)ii;
    for(jj=-RC;jj<=RC;jj++) {
      yt = (double)jj;
      dist = sqrt((xh-xt)*(xh-xt)+(yh-yt)*(yh-yt));
      if(dist<=((double)RC)&&dist!=0.0) {
	irc = (i + ii + LX)%LX;
	jrc = (j + jj + LY)%LY;
	idxrc = jrc + irc*LY;
        sigmax[idx] += sigma0*drand48()*h[idx]*h[idxrc]*(xh-xt)/pow(dist,psigma);
        sigmay[idx] += sigma0*drand48()*h[idx]*h[idxrc]*(yh-yt)/pow(dist,psigma);
      }
     }
    }
   }
  }
#else

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
    idx = j + i*LY;
    sigmax[idx] = sigma0*gasdev();
    sigmay[idx] = sigma0*gasdev();
   }
  }
#endif
  
}
//#endif


void compute_hrhs (double *fR, double *psi, double *hrhs) {


  for(i=0;i<LX;i++) {
  for(j=0;j<LY;j++) {
    idx = j +i*LY;
    hrhs[idx] = -fR[idx] + psi[idx];

  }
  }

}

#ifdef MELT_PONDS

void compute_flux (double *h, double *w, double *fluxx, double *fluxy) {

  int idxp,idxm;
  double dhx,dwx;
  double dhy,dwy;

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {

     idx = j + i*LY;

     ip = (i + 1 + LX)%LX;
     im = (i - 1 + LX)%LX;
     idxp = j + ip*LY;
     idxm = j + im*LY;

     dhx = 0.5*(h[idxp] - h[idxm]);
     dwx = 0.5*(w[idxp] - w[idxm]);

     fluxx[idx] = -alpha1*w[idx]*w[idx]*(dhx + dwx);

#ifdef WIND_SHEAR 
     fluxx[idx] += tsx[idx]*w[idx];
#endif 

#ifdef LATERAL_DRAINAGE 
     fluxx[idx] += -alphad*(dhx + dwx);
#endif 

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxp = jp + i*LY;
     idxm = jm + i*LY;

     dhy = 0.5*(h[idxp] - h[idxm]);
     dwy = 0.5*(w[idxp] - w[idxm]);

     fluxy[idx] = -alpha1*w[idx]*w[idx]*(dhy + dwy);  

#ifdef WIND_SHEAR 
     fluxy[idx] += tsy[idx]*w[idx];
#endif 

#ifdef LATERAL_DRAINAGE 
     fluxy[idx] += -alphad*(dhy + dwy);
#endif 

       }
  }

}  

#ifdef SEEPAGE
void compute_seepage (double *h, double *w) {

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {

     idx = j + i*LY;

     if(h[idx]>hmin) {
      s[idx] = kappa*w[idx]/h[idx];
     } else {
       s[idx] = 0.0;
     }

   }
  }

}
#endif 


void compute_wrhs (double *fR, double *fluxx, double *fluxy, double *s, double *wrhs) {

  int idxp,idxm;
  double dfluxx,dfluxy;
  double div;

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {

     idx = j + i*LY;

     ip = (i + 1 + LX)%LX;
     im = (i - 1 + LX)%LX;
     idxp = j + ip*LY;
     idxm = j + im*LY;

     dfluxx = 0.5*(fluxx[idxp] - fluxx[idxm]);

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxp = jp + i*LY;
     idxm = jm + i*LY;

     dfluxy = 0.5*(fluxy[idxp] - fluxy[idxm]);

     div = dfluxx + dfluxy;

     wrhs[idx] = fR[idx] - div + s[idx];

  }
  }

}
#endif /* ifdef melt ponds */

void compute_psi (double *sigmax, double *sigmay, double *psi) {

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
     dsxdx = 0.5 * (sigmax[idxp] - sigmax[idxm]);

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxp = jp + i*LY;
     idxm = jm + i*LY;
     dsydy = 0.5 * (sigmay[idxp] - sigmay[idxm]);


     psi[idx] = dsxdx + dsydy;
    } else {
      psi[idx] = 0.0;
    }
    
   }
  }
  
}

#ifdef CRANK_NICOLSON
  /*** to be implemented ***/
#elif defined RUNGE_KUTTA

void time_marching (double *h, double *k1, double *k2, double *k3, double *k4) {

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     h[idx] = h[idx] + (1./6.*k1[idx] + 1./3.*k2[idx] + 1./3.*k3[idx] + 1./6.*k4[idx]);

  }
 }

}

#else
void time_marching (double *h, double *rhs, double *rhsold) {

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     h[idx] = h[idx] + 1.5*dt*hrhs[idx] - 0.5*dt*hrhsold[idx];

  }
 }
}
#endif 



void total_volume (int t, double *h, FILE *fvol) {

  double vol;

  vol = 0.0;
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
    vol += h[idx];
  }
  }
   fprintf(fvol," %d %g \n",t,vol);
   fflush(fvol);
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
    fprintf(fout," %g %g %g \n",((double)i),((double)j),h[idx]);

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


