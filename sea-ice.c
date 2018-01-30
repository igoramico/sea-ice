#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "sea-ice.h" 
#include "variables.h"

typedef struct {

  char name[500];
  char value[500];

} param_type;

int main (int argc, char **argv) {

  int i,nparam;
  char finput[500];
  char fvolname[500];
  FILE *fin,*fvol;
  char str[500],value[500];
  param_type *Parameters;  
  int status;
  int seed;

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

    if((strcmp("temperature",Parameters[i].name))==0) {
      Temp = atof(Parameters[i].value);
    }

    if((strcmp("iterations",Parameters[i].name))==0) {
      NTSTEPS = atoi(Parameters[i].value);
    }

    if((strcmp("dump_frequency",Parameters[i].name))==0) {
      printfreq = atoi(Parameters[i].value);
    }

    if((strcmp("dump_config_frequency",Parameters[i].name))==0) {
      dumpfreq = atoi(Parameters[i].value);
    }

    if((strcmp("normalization_check_frequency",Parameters[i].name))==0) {
      volfreq = atoi(Parameters[i].value);
    }

    if((strcmp("system_size",Parameters[i].name))==0) {
      L = atoi(Parameters[i].value);
      LP6 = L + 6;
    }

    if((strcmp("time_step",Parameters[i].name))==0) {
      dt = atof(Parameters[i].value);
    }

    if((strcmp("base_file_name",Parameters[i].name))==0) {
      strcpy(basename,Parameters[i].value);
    }

    if((strcmp("restore",Parameters[i].name))==0) {
      if((strcmp("no",Parameters[i].value))==0) {
	restore = 0;
      } else {
	restore = 1;
      }
    }

    if((strcmp("output_directory",Parameters[i].name))==0) {
      strcpy(OutDir,Parameters[i].value);
    }

    if((strcmp("seed_random_generator",Parameters[i].name))==0) {
      seed = atoi(Parameters[i].value);
    }

    if((strcmp("starting_time",Parameters[i].name))==0) {
      start_time = atoi(Parameters[i].value);
    }

  }

  status = mkdir(OutDir, S_IRWXU | S_IRWXG );

  phi       = (double *)malloc(sizeof(double)*LP2*LP2);
  phirhs    = (double *)malloc(sizeof(double)*LP2*LP2);
  phirhsold = (double *)malloc(sizeof(double)*LP2*LP2);

  g       = (double *)malloc(sizeof(double)*LP2*LP2);
  grhs    = (double *)malloc(sizeof(double)*LP2*LP2);
  grhsold = (double *)malloc(sizeof(double)*LP2*LP2);

  if(!restore) {
    initialisation(g,phi);
    fprintf(stdout,"Configuration successfully initialised! \n");
  } else {
    /*
    restore_config(start_time,"h",h);
    restore_config(start_time,"rhs",rhs);
    restore_config(start_time,"rhsold",rhsold);
    fprintf(stdout,"Configuration successfully restored from a dump at t=%d! \n",start_time);
    */
  }

  print_phi(start_time,phi);
  print_g(start_time,g);

  sprintf(fvolname,"%s/gnorm_ti%d-tf%d.dat",OutDir,start_time,NTSTEPS);
  fvol = fopen(fvolname,"a");
  total_volume(start_time,g,fvol);
  fflush(fvol);

  for(iter=(start_time+1);iter<=NTSTEPS;iter++) {

   fprintf(stdout,"Starting %d time iteration... \n",iter);

   bc(h);
   deterministic_flux(nnn,h,Jdet);
   generate_stochastic(seed,Br);
   fluctuating_flux(nnn,h,Br,Jnoise);
   bc(Jdet);
   bc(Jnoise);
   compute_grhs(g,phi,grhs);
   compute_phirhs(g,phi,phirhs);

   if((!restore)&&(iter!=(start_time+1))) {
     copy_array(grhsold,grhs,LP2,LP2);
     copy_array(phirhsold,phirhs,LP2,LP2);
   }

   /*
   if((!((iter-1)%dumpfreq))&&((iter-1)!=start_time)) {
     dump_config((iter-1),"h",h);
     dump_config((iter-1),"rhs",rhs);
     dump_config((iter-1),"rhsold",rhsold);
   }
   */

   time_marching(g,phi,grhs,grhsold,phirhs,phirhsold);

   if(iter!=NTSTEPS) {
     copy_array(grhsold,grhs,LP2,LP2);
     copy_array(phirhsold,phirhs,LP2,LP2);
   }

   if(!(iter%printfreq)) {
     print_phi(iter,phi);
     print_g(iter,g);
   }

   if(!(iter%volfreq)) {
     total_volume(iter,g,fvol);
     fflush(fvol);
   }

   fprintf(stdout,"Completed. \n");

  }

  fclose(fvol);
  /*
  dump_config(NTSTEPS,"h",h);
  dump_config(NTSTEPS,"rhs",rhs);
  dump_config(NTSTEPS,"rhsold",rhsold);
  */

  fprintf(stdout,"Performed %d time iterations. Programme ended succesfully! \n",NTSTEPS);
  fprintf(stdout,"Good bye! \n");

  return 0;

}
