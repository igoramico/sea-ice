#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "sea-ice-1d.h"
#include "variables.h"

#define L 1000

void compute_fR (double *h, double *w, double *fR);

void compute_sigma (double *h, double *w, double *fR);

void compute_psi (double *h, double *w, double *fR);


int main (int argc, char **argv) {

  int i,nparam;
  char finput[500];
  char fvolname[500];
  FILE *fin,*fvol;
  char str[500],value[500];
  param_type *Parameters;  
  int status;
  int seed;

  int    i,ip,im;
  double *h,*w;
  double *flux,*fR,*s;
  double *sigma,*psi;
  double *hrhs,*wrhs;
  double *hrhsold,*wrhsold;

  int kk,KMIN,KMAX;
  double h0,eps,pert;
  double phi,epsf;

  if(argc!=2) {
    fprintf(stderr,"Error! Usage is: %s <input-file> \n",argv[0]);
    exit(2);
  }

  strcpy(finput,argv[1]);  
  fin = fopen(finput,"r");
  if(fin==NULL) { 
   fprintf(stderr,"Error! File %s cannot be opened...! \n",finput); 
   exit(2);
  }

  i=0;
  while((fscanf(fin,"%s %s",str,value)==2)) {
    ++i;
  }

  fclose(fin);

  nparam = i;

  Parameters = (param_type *)malloc(sizeof(param_type)*nparam);

  i=0;
  fin = fopen(finput,"r");
  while((fscanf(fin,"%s %s",str,value)==2)) {
    strcpy(Parameters[i].name,str);
    strcpy(Parameters[i].value,value);
    ++i;
  }
  fclose(fin);

  for(i=0;i<nparam;i++) {

    if((strcmp("system_size",Parameters[i].name))==0) {
      L = atoi(Parameters[i].value);
    }
    
    if((strcmp("iterations",Parameters[i].name))==0) {
      NTSTEPS = atoi(Parameters[i].value);
    }

    if((strcmp("time_step",Parameters[i].name))==0) {
      dt = atof(Parameters[i].value);
    }
        
    if((strcmp("dump_frequency",Parameters[i].name))==0) {
      printfreq = atoi(Parameters[i].value);
    }

    if((strcmp("Stefan_number",Parameters[i].name))==0) {
      S = atof(Parameters[i].value);
    }

    if((strcmp("mechanical_noise_amplitude",Parameters[i].name))==0) {
      sigma0 = atof(Parameters[i].value);
    }

    if((strcmp("base_file_name_ice_field",Parameters[i].name))==0) {
      strcpy(icename,Parameters[i].value);
    }

    if((strcmp("base_file_name_pond_field",Parameters[i].name))==0) {
      strcpy(pondname,Parameters[i].value);
    }
    
    if((strcmp("output_directory",Parameters[i].name))==0) {
      strcpy(OutDir,Parameters[i].value);
    }

    if((strcmp("restore",Parameters[i].name))==0) {
      if((strcmp("no",Parameters[i].value))==0) {
	restore = 0;
      } else {
	restore = 1;
      }
    }

    if((strcmp("seed_random_generator",Parameters[i].name))==0) {
      seed = atoi(Parameters[i].value);
    }

  }

  status = mkdir(OutDir, S_IRWXU | S_IRWXG );

  h       = (double *)malloc(sizeof(double)*L);
  hrhs    = (double *)malloc(sizeof(double)*L);
  hrhsold = (double *)malloc(sizeof(double)*L);

#ifdef MELT_PONDS  
  w       = (double *)malloc(sizeof(double)*L);
  wrhs    = (double *)malloc(sizeof(double)*L);
  wrhsold = (double *)malloc(sizeof(double)*L);
#endif
  
  if(!restore) {
#ifdef MELT_PONDS
    initialisation(h,w);
#else
    initialisation(h);
#endif    
    fprintf(stdout,"Configuration successfully initialised! \n");
  } else {
    /*
    restore_config(start_time,"h",h);
    restore_config(start_time,"rhs",rhs);
    restore_config(start_time,"rhsold",rhsold);
    fprintf(stdout,"Configuration successfully restored from a dump at t=%d! \n",start_time);
    */
  }

  fvol = fopen("total_ice_volume.dat","w");
  
  for(time=0;time<=NTIME;time++) {

   /*** evaluate rhs's ***/

 #ifdef MELT_PONDS   
    compute_fR(h,w,fR);
    compute_sigma(h,w,sigma);
 #else
    compute_fR(h,fR);
    compute_sigma(h,sigma);
 #endif
    compute_psi(sigma,psi);
    compute_hrhs(fR,psi,hrhs);
#ifdef MELT_PONDS
    compute_wrhs(fR,flux,s,wrhs);
    time_marching(w,wrhs,wrhsold);
#endif    
    time_marching(h,hrhs,hrhsold);
    copy_array(hrhsold,hrhs,L);
#ifdef MELT_PONDS
    copy_array(wrhsold,wrhs,L);
#endif
    
   if((time%printfreq)==0) {
     print_f(time,h,icename);
#ifdef MELT_PONDS
     print_f(time,w,pondname);
#endif
   }

   total_volume(time,h,fvol);
   
  } /* time loop */

  return 0;

}

