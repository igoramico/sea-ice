#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

int main (int argc, char **argv) {

  int i,j,idx;
  int LX,LY;
  double *h,*w;
  FILE *fin;

#ifdef INIT_MELT_WATER_ZERO
  
  if(argc!=3) {
    fprintf(stderr,"Error! Usage is: %s <size x> <size y> \n",argv[0]);
    exit(2);
  }

  LX = atoi(argv[1]);
  LY = atoi(argv[2]);

  w = (double *)malloc(sizeof(double)*LX*LY);

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     w[idx] = 0.0;
   }
  }

  fin = fopen("init_melt_water_distribution.inp","w");
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     fprintf(fin,"%g\n",w[idx]);
   }
  }
  fclose(fin);

  free(w);

#endif
  
  
#ifdef INIT_MELT_WATER_PATCH  

  double x,x1,x2,dx;
  double y,y1,y2,dy;
  double w0;
  
  if(argc!=10) {
   fprintf(stderr,"Error! Usage is: %s <size x> <size y> <[x1,x2]> <[y1,y2]> <dx> <dy> <water level> \n",argv[0]);
    exit(2);
  }

  LX = atoi(argv[1]);
  LY = atoi(argv[2]);
  x1 = atof(argv[3]);
  x2 = atof(argv[4]);
  y1 = atof(argv[5]);
  y2 = atof(argv[6]);
  dx = atof(argv[7]);
  dy = atof(argv[8]);
  w0 = atof(argv[9]);

  w = (double *)malloc(sizeof(double)*LX*LY);
  
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     x = (double)i;
     y = (double)j;
     if((x>(x1-5.0*dx)&&x<(x2+5.0*dx))&&(y>(y1-5.0*dy)&&y<(y2+5.0*dy))) {
       w[idx] = w0*0.25*(1.0+tanh((x-x1)/dx)*tanh((x2-x)/dx))*(1.0+tanh((y-y1)/dy)*tanh((y2-y)/dy));
     } else {
       w[idx] = 0.0;
     }
   }
  }

  fin = fopen("init_melt_water_distribution.inp","w");
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     fprintf(fin,"%g\n",w[idx]);
   }
  }
  fclose(fin);

  free(w);
  
#endif

#ifdef INIT_ICE_GAUSSIAN

  double hmean,hvar;
  double hmin,hmax;
  double r;
  int ok;
  
  if(argc!=7) {
   fprintf(stderr,"Error! Usage is: %s <size x> <size y> <mean> <standard deviation> <maximum> <minimum> \n",argv[0]);
    exit(2);
  }

  LX    = atoi(argv[1]);
  LY    = atoi(argv[2]);
  hmean = atof(argv[3]);
  hvar  = atof(argv[4]);
  hmax  = atof(argv[5]);
  hmin  = atof(argv[6]);

  h = (double *)malloc(sizeof(double)*LX*LY);

  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     ok = 0;
     while(ok==0) {
      r = hvar*gasdev() + hmean;
      if(r<=hmax&&r>=hmin) {
	ok=1;
      }
     }
     h[idx] = r;
   }
  }

  fin = fopen("init_ice_topography.inp","w");
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i*LY;
     fprintf(fin,"%g\n",h[idx]);
   }
  }
  fclose(fin);

  free(h);

#endif
  
  return 0;

}
