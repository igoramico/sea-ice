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

  epsh0 = eps*h0;
  
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
      idx = j + i*LY;
      h[idx] = h0;
      w[idx] = 0.0;
     }
    }

#ifdef INIT_ICE_RANDOM

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
      idx = j + i*LY;
      r = (2.0*drand48() - 1.0);
      h[idx] = h0 + epsh0*r;
     }
    }

#endif

#ifdef INIT_ICE_SINGLEMODE

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       idx = j + i*LY;
       h[idx] = h0 + epsh0*sin((2.0*(double)KMIN*M_PI*(double)i/((double)LX)))*sin((2.0*(double)KMIN*M_PI*(double)j/((double)LY)));
     }
    }
    
#endif
    
#ifdef INIT_ICE_FOURIER

    double phi;
    double pert,epsf;
    int NK;
    double ff=7./8.;
    
    epsf  = eps*M_PI;

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


#ifdef INIT_ICE_FROM_FILE

    sprintf(fname,"init_ice_topography.inp");
    finp = fopen(fname,"r");
    if(finp==NULL) {
      fprintf(stderr,"Error! Initial topography file %s not found!! \n",fname);
      exit(2);
    } else {
      for(i=0;i<LX;i++) {
       for(j=0;j<LY;j++) {
	 idx = j + i*LY;
	fscanf(finp,"%lf",&h[idx]);
       }
      }
      fclose(finp);
    }


#endif

#ifdef INIT_MELT_PONDS_FROM_FILE

    sprintf(fname,"init_melt_water_distribution.inp");
    finp = fopen(fname,"r");
    if(finp==NULL) {
      fprintf(stderr,"Error! Initial melt water distribution file '%s' not found!! \n",fname);
      exit(2);
    } else {
      for(i=0;i<LX;i++) {
       for(j=0;j<LY;j++) {
	 idx = j + i*LY;
	fscanf(finp,"%lf",&w[idx]);
      }
      }
     fclose(finp); 
    }

#endif
    
}

#else

void initialisation (double *h) {

  double r,epsh0;

  epsh0 = eps*h0;
  
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
      idx = j + i*LY;
      h[idx] = h0;
     }
    }
    
#ifdef INIT_ICE_RANDOM

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
      idx = j + i*LY;
      r = (2.0*drand48() - 1.0);
      h[idx] = h0 + epsh0*r;
     }
    }
#endif

#ifdef INIT_ICE_SINGLEMODE

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       idx = j + i*LY;
       h[idx] = h0 + epsh0*sin((2.0*(double)KMIN*M_PI*(double)i/((double)LX)))*sin((2.0*(double)KMIN*M_PI*(double)j/((double)LY)));
     }
    }
    
#endif
    
#ifdef INIT_ICE_FOURIER

    double phi;
    double pert,epsf;
    int NK;
    double ff=7./8.;
    
    epsf  = eps*M_PI;

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

    finp = fopen("init_ice_topography.inp","r");
    if(finp==NULL) {
      fprintf(stderr,"Error! Initial topography file not found!! \n");
      exit(2);
    } else {
      for(i=0;i<LX;i++) {
       for(j=0;j<LY;j++) {
	 idx = j + i*LY;
	fscanf(finp,"%lf",&h[idx]);
      }
    }
      fclose(finp);
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

#ifdef TEMPERATURE_FIELD

#endif

#ifdef MELT_PONDS
void compute_fR (double *h, double *w, double *fR) {

  double meltrate,meltenh;
  
#ifdef MELTING_LATERAL  
  int idxpx,idxmx;
  int idxpy,idxmy;
  double dwx,dwy;
  double normal;
  double fVert,fLongX,fLongY;
#endif
  
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     fR[idx] = 0.0;
   }
  }
     
#ifdef MELTING_LUETHJE
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     if((h[idx]>hmin)&&(mluethje!=0.0)) {
       if(w[idx]>wmin_melt) {

#ifdef MELTING_WIND_SHEAR
	 
	 tsmod = sqrt(ts0x*ts0x + ts0y*ts0y);
	 mpluethje *= (1.0/(1.0 + ws_luethje*pow(tsmod,(3./4.))*pow(w[idx],(7./4.)))
		       + ws_luethje*tsmod*(w[idx]/WWSmax)*(w[idx]/WWSmax));

#endif	 
	 
#ifdef MELTING_LUETHJE_TWO_SEVENTH

         meltenh = (1.0 + mpluethje/mluethje*pow((wmaxfr/w[idx]),(1./7.)));

#else

	 if(w[idx]<wmaxfr) {
           meltenh = (1.0 + mpluethje/mluethje*w[idx]/wmaxfr);
         } else {
	   meltenh = (1.0 + mpluethje/mluethje);
         }

         meltrate = meltenh*mluethje;

#endif /* ifdef melting_luethje_two_seventh */
	 	 
#ifdef MELTING_LATERAL

     ip = (i + 1 + LX)%LX;
     im = (i - 1 + LX)%LX;
     idxpx = j + ip*LY;
     idxmx = j + im*LY;

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxpy = jp + i*LY;
     idxmy = jm + i*LY;
    
     dwx = 0.5*(w[idxpx] - w[idxmx])/dx;
     dwy = 0.5*(w[idxpy] - w[idxmy])/dx;
    
     normal = sqrt(1.0 + dwx*dwx + dwy*dwy);
     fVert  = 1.0/normal;
     fLongX = fabs(dwx)/normal;
     fLongY = fabs(dwy)/normal;

     fR[idx] += fVert*meltrate;

     if(h[idxpx] > h[idx]) {
       fR[idxpx] += fLongX*meltrate;
     }
     if(h[idxmx] > h[idx]) {
       fR[idxmx] += fLongX*meltrate;
     }

     if(h[idxpy] > h[idx]) {
       fR[idxpy] += fLongY*meltrate;
     }
     if(h[idxmy] > h[idx]) {
       fR[idxmy] += fLongY*meltrate;
     }
     
#else

     fR[idx] = meltrate;
     
#endif /* ifdef melting_lateral */     

       } else {
	fR[idx] = mluethje;
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
     fR[idx] = 0.0;
   }
  }
  
#ifdef MELTING_WETTLAUFER  
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
  
}


#endif /* endif MELT_PONDS */

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

#ifdef VOLUMETRIC_FLUXES
void compute_flux (double *h, double *w, double *Flux) {

  int idxpx,idxmx;
  int idxpy,idxmy;
  double wm;
  double Jpx,Jmx,Jpy,Jmy;
  
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {

     idx = j + i*LY;

     ip = (i + 1 + LX)%LX;
     im = (i - 1 + LX)%LX;
     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;

     idxpx = j + ip*LY;
     idxmx = j + im*LY;

     idxpy = jp + i*LY;
     idxmy = jm + i*LY;

     wm  = 0.5 * (w[idxpx] + w[idx]);
     Jpx = alpha1 * wm * wm * wm * (w[idxpx] + h[idxpx] - w[idx] - h[idx]); 
#ifdef WIND_SHEAR 
     Jpx -= tsx[idx]*wm; 
#endif 

     wm  = 0.5 * (w[idxpy] + w[idx]);
     Jpy = alpha1 * wm * wm * wm * (w[idxpy] + h[idxpy] - w[idx] - h[idx]); 
#ifdef WIND_SHEAR 
     Jpy -= tsy[idx]*wm; 
#endif 

     wm  = 0.5 * (w[idxmy] + w[idx]);
     Jmy = alpha1 * wm * wm * wm * (w[idxmy] + h[idxmy] - w[idx] - h[idx]); 
#ifdef WIND_SHEAR 
     Jmy += tsy[idx]*wm; 
#endif 

     wm  = 0.5 * (w[idxmx] + w[idx]);
     Jmx = alpha1 * wm * wm * wm * (w[idxmx] + h[idxmx] - w[idx] - h[idx]); 
#ifdef WIND_SHEAR 
     Jmx += tsx[idx]*wm; 
#endif 
     
     Flux[idx] = (Jpx + Jpy + Jmx + Jmy)/dx;
     
   }
  }

}  
#else
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

     dhx = 0.5*(h[idxp] - h[idxm])/dx;
     dwx = 0.5*(w[idxp] - w[idxm])/dx;

     fluxx[idx] = -alpha1*w[idx]*w[idx]*(dhx + dwx);

#ifdef WIND_SHEAR 
     fluxx[idx] += tsx[idx]*w[idx];
#endif 

#ifdef LATERAL_DRAINAGE
     fluxx[idx] += -alphad * w[idx] * (dhx + dwx);
#endif 

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxp = jp + i*LY;
     idxm = jm + i*LY;

     dhy = 0.5*(h[idxp] - h[idxm])/dx;
     dwy = 0.5*(w[idxp] - w[idxm])/dx;

     fluxy[idx] = -alpha1*w[idx]*w[idx]*(dhy + dwy);  

#ifdef WIND_SHEAR 
     fluxy[idx] += tsy[idx]*w[idx];
#endif 

#ifdef LATERAL_DRAINAGE 
     fluxy[idx] += -alphad * w[idx] * (dhy + dwy);
#endif 
 
       }
  }

}  
#endif /* if VOLUMETRIC_FLUXES */


void compute_seepage (double *h, double *w, double *s) {

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {

     idx = j + i*LY;

     s[idx] = 0.0;

#ifdef SEEPAGE_CONSTANT
     if(w[idx]>0.0) {
       s[idx] = s0;
     } else {
       s[idx] = 0.0;
     }
#endif

#ifdef SEEPAGE_DARCY
     if(h[idx]>hmin) {
      s[idx] = kappa*w[idx]/h[idx];
     } else {
       s[idx] = 0.0;
     }
#endif
     
   }
  }

}

#ifdef MELT_PONDS_SUPERSEEPER

void melt_ponds_superseeper(double *w, double wminss) {

    
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     if(w[idx]<wminss) {
       w[idx] = 0.0;
     }
   }
  }

}
#endif

#ifdef VOLUMETRIC_FLUXES
void compute_wrhs (double *fR, double *Flux, double *s, double *wrhs) {

  int idxp,idxm;
  double dfluxx,dfluxy;
  double div;

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {

     idx = j + i*LY;

     wrhs[idx] = fR[idx] + Flux[idx] + s[idx];

   }
  }

}
#else
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

     dfluxx = 0.5*(fluxx[idxp] - fluxx[idxm])/dx;

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxp = jp + i*LY;
     idxm = jm + i*LY;

     dfluxy = 0.5*(fluxy[idxp] - fluxy[idxm])/dx;

     div = dfluxx + dfluxy;
     wrhs[idx] = fR[idx] - div + s[idx];
   }
  }

}
#endif

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
     dsxdx = 0.5 * (sigmax[idxp] - sigmax[idxm])/dx;

     jp = (j + 1 + LY)%LY;
     jm = (j - 1 + LY)%LY;
     idxp = jp + i*LY;
     idxm = jm + i*LY;
     dsydy = 0.5 * (sigmay[idxp] - sigmay[idxm])/dx;


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
     h[idx] = h[idx] + 1.5*dt*rhs[idx] - 0.5*dt*rhsold[idx];
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


void stabiliser (double *h, double cutoff) {

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
    if(h[idx]<cutoff) { h[idx] = cutoff; }

   }
  }

}


#ifdef MELT_PONDS
void print_f (int t, double *h, double *w, char icebasename[], char pondbasename[]) {

  FILE *fout;
  char filename[500];

  sprintf(filename,"%s/%s_%s_t%d.dat",OutDir,icebasename,pondbasename,t);
  fout = fopen(filename,"w");
  if(fout==NULL) { fprintf(stderr,"Error! File %s could not be opened! \n",filename); exit(2); }
  
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     fprintf(fout," %g %g %g %g \n",((double)i),((double)j),h[idx],w[idx]);
   }
   fprintf(fout,"\n");
  }

  fclose(fout);

}
#else
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
#endif

void copy_array (double *v1, double *v2, int N) {

  for(i=0;i<N;i++) {
    v1[i] = v2[i];
  }

}




