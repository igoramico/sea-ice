#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char **argv) {

  int i,j,idx;
  int LX,LY;
  double x,x1,x2,dx;
  double y,y1,y2,dy;
  double *w;
  double w0;
  FILE *fin;

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
  
  return 0;

}
