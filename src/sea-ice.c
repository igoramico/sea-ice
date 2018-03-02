#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include "variables.h"
#include "sea-ice.h"

typedef struct {

  char name[500];
  char value[500];
  
} param_type;


int main (int argc, char **argv) {

  int i,nparam;
  char finput[500];
  char fvolname[500];
  FILE *fin,*fvol,*fwvol;
  char str[500],value[500];
  param_type *Parameters;  
  int status;
  int seed;

  srand48(time(NULL));
  
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

    if((strcmp("system_size_x",Parameters[i].name))==0) {
      LX = atoi(Parameters[i].value);
    }

    if((strcmp("system_size_y",Parameters[i].name))==0) {
      LY = atoi(Parameters[i].value);
    }

    if((strcmp("minimum_allowed_height",Parameters[i].name))==0) {
      hmin = atof(Parameters[i].value);
    }

    if((strcmp("initial_height",Parameters[i].name))==0) {
      h0 = atof(Parameters[i].value);
    }
    
    if((strcmp("initial_perturbation_amplitude",Parameters[i].name))==0) {
      eps = atof(Parameters[i].value);
    }

    if((strcmp("iterations",Parameters[i].name))==0) {
      NTIME = atoi(Parameters[i].value);
    }

    if((strcmp("time_step",Parameters[i].name))==0) {
      dt = atof(Parameters[i].value);
    }

    if((strcmp("k_min_mode",Parameters[i].name))==0) {
      KMIN = atoi(Parameters[i].value);
    }

    if((strcmp("k_max_mode",Parameters[i].name))==0) {
      KMAX = atoi(Parameters[i].value);
    }

    if((strcmp("dump_frequency",Parameters[i].name))==0) {
      printfreq = atoi(Parameters[i].value);
    }

    if((strcmp("volume_check_frequency",Parameters[i].name))==0) {
      volfreq = atoi(Parameters[i].value);
    }
    
    if((strcmp("Stefan_number",Parameters[i].name))==0) {
      S = atof(Parameters[i].value);
    }

    if((strcmp("mechanical_noise_amplitude",Parameters[i].name))==0) {
      sigma0 = atof(Parameters[i].value);
    }

    if((strcmp("mechanical_noise_kernel_power",Parameters[i].name))==0) {
      psigma = atof(Parameters[i].value);
    }

    if((strcmp("mechanical_noise_kernel_range",Parameters[i].name))==0) {
      RC = atoi(Parameters[i].value);
    }

#ifdef MELT_PONDS
    
    if((strcmp("melting_luethje_m",Parameters[i].name))==0) {
      mluethje = atoi(Parameters[i].value);
    }

    if((strcmp("melting_luethje_mp",Parameters[i].name))==0) {
      mpluethje = atoi(Parameters[i].value);
    }

    if((strcmp("melting_luethje_wmax",Parameters[i].name))==0) {
      wmaxfr = atoi(Parameters[i].value);
    }
    
    if((strcmp("meltwater_flux_alpha1",Parameters[i].name))==0) {
      alpha1 = atoi(Parameters[i].value);
    }

    if((strcmp("meltwater_flux_wind_x",Parameters[i].name))==0) {
      ts0x = atoi(Parameters[i].value);
    }

    if((strcmp("meltwater_flux_wind_y",Parameters[i].name))==0) {
      ts0y = atoi(Parameters[i].value);
    }

    if((strcmp("vertical_seepage_rate_kappa",Parameters[i].name))==0) {
      kappa = atoi(Parameters[i].value);
    }

    if((strcmp("lateral_drainage_rate_alphad",Parameters[i].name))==0) {
      alphad = atoi(Parameters[i].value);
    }

    if((strcmp("base_file_name_pond_field",Parameters[i].name))==0) {
      strcpy(pondname,Parameters[i].value);
    }
    
#endif
    
    if((strcmp("base_file_name_ice_field",Parameters[i].name))==0) {
      strcpy(icename,Parameters[i].value);
    }
    
    if((strcmp("output_directory",Parameters[i].name))==0) {
      strcpy(OutDir,Parameters[i].value);
    }

    if((strcmp("restore",Parameters[i].name))==0) {
      if((strcmp("no",Parameters[i].value))==0) {
	restore = 0;
      } else {
	restore = 1;
        if((strcmp("start_time",Parameters[i].name))==0) {
          start_time = atoi(Parameters[i].value);
        }
      }
    }

    if((strcmp("seed_random_generator",Parameters[i].name))==0) {
      seed = atoi(Parameters[i].value);
    }

#ifdef MELT_PONDS

    if((strcmp("melt_ponds_alpha1_flux",Parameters[i].name))==0) {
      alpha1 = atof(Parameters[i].value);
    }

    if((strcmp("melt_ponds_alphad_lateral_drainage",Parameters[i].name))==0) {
      alphad = atof(Parameters[i].value);
    }

    if((strcmp("melt_ponds_kappa_seepage",Parameters[i].name))==0) {
      kappa = atof(Parameters[i].value);
    }

#endif 

  }

  TOTSIZE = LX*LY;
  
  status = mkdir(OutDir, S_IRWXU | S_IRWXG );

  h       = (double *)malloc(sizeof(double)*TOTSIZE);
  hrhs    = (double *)malloc(sizeof(double)*TOTSIZE);
  hrhsold = (double *)malloc(sizeof(double)*TOTSIZE);
  
  psi     = (double *)malloc(sizeof(double)*TOTSIZE);
  sigmax  = (double *)malloc(sizeof(double)*TOTSIZE);
  sigmay  = (double *)malloc(sizeof(double)*TOTSIZE);
  fR      = (double *)malloc(sizeof(double)*TOTSIZE);

#ifdef MELT_PONDS  
  w       = (double *)malloc(sizeof(double)*TOTSIZE);
  wrhs    = (double *)malloc(sizeof(double)*TOTSIZE);
  wrhsold = (double *)malloc(sizeof(double)*TOTSIZE);

  fluxx   = (double *)malloc(sizeof(double)*TOTSIZE);
  fluxy   = (double *)malloc(sizeof(double)*TOTSIZE);

  s       = (double *)malloc(sizeof(double)*TOTSIZE);
#endif
  
  if(!restore) {
    start_time = 0;
#ifdef MELT_PONDS
    initialisation(h,w);
#else
    initialisation(h);
#endif    
    fprintf(stdout,"Configuration successfully initialised! \n");
  } else {
    restore_config(start_time,"h",h);
    restore_config(start_time,"hrhs",hrhs);
#ifndef RUNGE_KUTTA
    restore_config(start_time,"hrhsold",hrhsold);
#endif

#ifdef MELT_PONDS 
    restore_config(start_time,"w",w);
    restore_config(start_time,"wrhs",wrhs);
#ifndef RUNGE_KUTTA
    restore_config(start_time,"wrhsold",wrhsold);
#endif

#endif

    fprintf(stdout,"Configuration successfully restored from a dump at t=%d! \n",start_time);
  }

  fvol = fopen("total_ice_volume.dat","a");

#ifdef MELT_PONDS
  fwvol = fopen("total_water_volume.dat","a");
#endif
  
  for(iter=start_time;iter<=(start_time+NTIME);iter++) {

   if((iter%printfreq)==0) {
     print_f(iter,h,icename);
#ifdef MELT_PONDS
     print_f(iter,w,pondname);
#endif
   }

   if((iter%volfreq)==0) {
    total_volume(iter,h,fvol);
#ifdef MELT_PONDS
    total_volume(iter,w,fwvol);
#endif    
   }

    
    fprintf(stdout,"Starting calculations at iteration #%d...",iter);
    
   /*** evaluate rhs's ***/

#ifdef RUNGE_KUTTA

    /* Runge-Kutta time marching to be implemented... */

#else /* else time marching */

 #ifdef MELT_PONDS   
    compute_fR(h,w,fR);
    //    compute_sigma(h,w,sigmax,sigmay);
 #else
    compute_fR(h,fR);
    //    compute_sigma(h,sigmax,sigmay);
 #endif
    compute_sigma(h,sigmax,sigmay);
    compute_psi(sigmax,sigmay,psi);
    compute_hrhs(fR,psi,hrhs);

 #ifdef MELT_PONDS
    compute_flux(h,w,fluxx,fluxy);
    compute_wrhs(fR,fluxx,fluxy,s,wrhs);
    if(iter==0) {
      copy_array(wrhsold,wrhs,TOTSIZE);
    }    
    time_marching(w,wrhs,wrhsold);
 #endif
    if(iter==0) {
      copy_array(hrhsold,hrhs,TOTSIZE);
    }    
    time_marching(h,hrhs,hrhsold);
    copy_array(hrhsold,hrhs,TOTSIZE);
 #ifdef MELT_PONDS
    copy_array(wrhsold,wrhs,TOTSIZE);
 #endif

#endif /* endif time marching */
    
   fprintf(stdout,"Done!\n");

  } /* time loop */

  fclose(fvol);
#ifdef MELT_PONDS  
  fclose(fwvol);
#endif
  
  return 0;

}

