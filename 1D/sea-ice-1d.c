#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define L 1000

void compute_fR (double *h, double *w, double *fR);

int main (int argc, char **argv) {

  int    i,ip,im;
  double *h,*w;
  double *flux,*fR,*s,*psi;
  double *hrhs,*wrhs;
  double *hrhsold,*wrhsold;

  int kk,KMIN,KMAX;
  double h0,eps,pert;
  double phi,epsf;

  for(time=0;time<=NTIME;time++) {

    /*** initialisation ***/
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

#ifdef INIT_PONDS_NONE
    for(i=0;i<L;i++) {

      w[i] = 0.0;

    }
#endif 

#ifdef INIT_PONDS_RANDOM
    while(init_ponds_volume>0.0) {
      r = drand48();
      i = (int)floor(r*L);
      r = drand48();
      /* to be continued */
    }
#endif

   /*** evaluate rhs's ***/

   for(i=0;i<L;i++) {

     ip = (i+1+L)%L;
     im = (i-1+L)%L;

     hrhs[i] = -fR[i] + psi[i];
     flux[i] = 0.5*(w[i]*(u[ip]-u[im])+u[i]*(w[ip]-w[im]));
     wrhs[i] = fR[i] - s[i] - flux[i];

   }

   /*** time marching ***/

   for(i=0;i<L;i++) {

     h[i] = h[i] + 1.5*hrhs[i] - 0.5*hrhsold[i];
     w[i] = w[i] + 1.5*wrhs[i] - 0.5*wrhsold[i];

     hrhsold[i] = hrhs[i];
     wrhsold[i] = wrhs[i];

   }

   if((time%printfreq)==0) {
     sprintf(fname,"%s_t%d.dat",basename,time);
     fout = fopen(fname,"w");
     for(i=0;i<L;i++) {
       fprintf(fout," %d %g %g \n",i,h[i],w[i]);
     }
     fclose(fout);
   }

  } /* time loop */

  return 0;

}

void compute_fR (double *h, double *w, double *fR) {

  int j;

#ifdef NO_MELT_PONDS
  for(j=0;j<L;j++) {
    fR[j] = fRin; 
  }
#endif

#ifdef EISENMANN_WETTLAUFER

#endif 

}
