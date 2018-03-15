#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char **argv) {

  int i,j,n;
  int LX,LY;
  int BIN;
  double h;
  double max,min,delta;
  double *pdf;
  FILE *fin;
  
  if(argc!=6) {
    fprintf(stderr,"Error! Usage is: %s <start time> <end time> <time step> <LX> <LY> <BIN> \n",argv[0]);
    exit(2);
  }

  LX  = atoi(argv[1]);
  LY  = atoi(argv[2]);
  BIN = atoi(argv[5]);


  for(it=STIME;it<=ETIME;it+=TSTEP) {

#ifdef MELT_PONDS    
    sprintf(fname,"%s/%s_%s_t%d.dat",path,icename,pondname,it); 
#else
    sprintf(fname,"%s/%s_t%d.dat",path,icename,it); 
#endif

    fin = fopen(fname,"r");
    
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       idx = j + i*LY;

       fscanf(fin," %d %d %lf %lf \n",&idum,&idum,h[idx],w[idx]);
       if(h[idx]>hmax) {
	 hmax = h[idx];
       }
       if(h[idx]<hmin) {
	 hmin = h[idx];
       }
       
   delta = (max - min)/((double)BIN);

  fin = fopen("init_ice_topography.inp","r");
  
  pdf = (double *)malloc(sizeof(double)*BIN);
  for(n=0;n<BIN;n++) { pdf[n] = 0.0; }
  
  for(i=0;i<LX;i++) {
    for(j=0;j<LY;j++) {
      fscanf(fin,"%lf",&h);
      n = (int)floor((h-min)/delta);
      pdf[n] += 1.0;
    }
  }

  for(n=0;n<BIN;n++) {
    fprintf(stdout," %g %g \n",(min+(n+0.5)*delta),pdf[n]/(LX*LY));
  }
  
  free(pdf);
  
  return 0;

}
