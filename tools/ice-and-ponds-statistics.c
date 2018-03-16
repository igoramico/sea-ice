#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main (int argc, char **argv) {

  int i,j,idx;
  int LX,LY;
  int n,BIN;
  int it,STIME,ETIME,TSTEP;
  int idum;
  double *h;
  double hmax,hmin,hdelta;
  double *pdfh;
  FILE *fin,*fout;
  char path[128];
  char fname[128];
  char icename[128];
#ifdef MET_PONDS
  double *w;
  char pondmame[128];
  double wmax,wmin,wdelta;
  double *pdfw;
#endif
  
#ifdef MELT_PONDS
  
  if(argc!=10) {
    fprintf(stderr,"Error! Usage is: %s <path> <ice basename> <pond basename> <start time> <end time> <time step> <LX> <LY> <BIN> \n",argv[0]);
    exit(2);
  }

  strcpy(path,argv[1]);
  strcpy(icename,argv[2]);
  strcpy(pondname,argv[3]);
  STIME = atoi(argv[4]);
  ETIME = atoi(argv[5]);
  TSTEP = atoi(argv[6]);
  LX    = atoi(argv[7]);
  LY    = atoi(argv[8]);
  BIN   = atoi(argv[9]);

#else

  if(argc!=9) {
    fprintf(stderr,"Error! Usage is: %s <path> <ice basename> <start time> <end time> <time step> <LX> <LY> <BIN> \n",argv[0]);
    exit(2);
  }

  strcpy(path,argv[1]);
  strcpy(icename,argv[2]);
  STIME = atoi(argv[3]);
  ETIME = atoi(argv[4]);
  TSTEP = atoi(argv[5]);
  LX    = atoi(argv[6]);
  LY    = atoi(argv[7]);
  BIN   = atoi(argv[8]);

#endif  
  //  totpondarea = 0.0;
  //  avepondarea = 0.0;

  pdfh = (double *)malloc(sizeof(double)*BIN);
#ifdef MELT_PONDS
  pdfw = (double *)malloc(sizeof(double)*BIN);
#endif
  
  for(it=STIME;it<=ETIME;it+=TSTEP) {

#ifdef MELT_PONDS
    for(n=0;n<BIN;n++) { pdfw[n]=pdfh[n]=0.0; }
    sprintf(fname,"%s/%s_%s_t%d.dat",path,icename,pondname,it); 
#else
    for(n=0;n<BIN;n++) { pdfh[n]=0.0; }
    sprintf(fname,"%s/%s_t%d.dat",path,icename,it);     
#endif

    fin = fopen(fname,"r");
    
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       idx = j + i*LY;

#ifdef MELT_PONDS
       fscanf(fin," %d %d %lf %lf \n",&idum,&idum,h[idx],w[idx]);
#else
       fscanf(fin," %d %d %lf \n",&idum,&idum,h[idx]);
#endif
       
       if(h[idx]>hmax) {
	 hmax = h[idx];
       }
       if(h[idx]<hmin) {
	 hmin = h[idx];
       }

#ifdef MELT_PONDS       
       if(w[idx]>wmax) {
	 wmax = w[idx];
       }
       if(w[idx]<wmin) {
	 wmin = w[idx];
       }
#endif
       
     }
    }

    fclose(fin);
    
    hmax *= 1.1;
    hmin *= 1.1;
    
    hdelta = (hmax - hmin)/((double)BIN);

#ifdef MELT_PONDS    
    wmax *= 1.1;
    wmin *= 1.1;
    
    wdelta = (wmax - wmin)/((double)BIN);
#endif
    
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       idx = j + i*LY;

       n = (int)floor((h[idx]-hmin)/hdelta);
       pdfh[n] += 1.0;

#ifdef MELT_PONDS       
       n = (int)floor((w[idx]-wmin)/wdelta);
       pdfw[n] += 1.0;
#endif
      
     }
    }


    sprintf(fname,"pdf_%s_t%d.dat",icename,it); 
    fout = fopen(fname,"w");
    for(n=0;n<BIN;n++) {
      fprintf(stdout," %g %g \n",(hmin+(n+0.5)*hdelta),pdfh[n]/(LX*LY));
    }    
    fclose(fout);
    
#ifdef MELT_PONDS    

    sprintf(fname,"pdf_%s_t%d.dat",pondname,it); 
    fout = fopen(fname,"w");
    for(n=0;n<BIN;n++) {
      fprintf(stdout," %g %g \n",(wmin+(n+0.5)*wdelta),pdfw[n]/(LX*LY));
    }    
    fclose(fout);
    
#endif
    
    
  }
  
  free(pdfh);
#ifdef MELT_PONDS  
  free(pdfw);
#endif
  
  return 0;

}
