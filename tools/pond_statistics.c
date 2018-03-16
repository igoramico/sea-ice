#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "hk.h"

int *labels;
int  n_labels = 0;     /* length of the labels array */

/*  uf_find returns the canonical label for the equivalence class containing x */

int uf_find(int x) {
  int y = x;
  while (labels[y] != y)
    y = labels[y];
  
  while (labels[x] != x) {
    int z = labels[x];
    labels[x] = y;
    x = z;
  }
  return y;
}

/*  uf_union joins two equivalence classes and returns the canonical label of the resulting class. */

int uf_union(int x, int y) {
  return labels[uf_find(x)] = uf_find(y);
}

/*  uf_make_set creates a new equivalence class and returns its label */

int uf_make_set(void) {
  labels[0] ++;
  assert(labels[0] < n_labels);
  labels[labels[0]] = labels[0];
  return labels[0];
}

/*  uf_intitialize sets up the data structures needed by the union-find implementation. */

void uf_initialize(int max_labels) {
  n_labels = max_labels;
  labels = calloc(sizeof(int), n_labels);
  labels[0] = 0;
}

/*  uf_done frees the memory used by the union-find data structures */

void uf_done(void) {
  n_labels = 0;
  free(labels);
  labels = 0;
}

/* End Union-Find implementation */

#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

/* print_matrix prints out a matrix that is set up in the "pointer to pointers" scheme
   (aka, an array of arrays); this is incompatible with C's usual representation of 2D
   arrays, but allows for 2D arrays with dimensions determined at run-time */

void print_matrix(int **matrix, int m, int n) {
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++)
      printf("%3d ",matrix[i][j]);
    printf("\n");
  }
}


/* Label the clusters in "matrix".  Return the total number of clusters found. */

int hoshen_kopelman(int **matrix, int m, int n) {
  
  uf_initialize(m * n / 2);
  
  /* scan the matrix */
  
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {                        // if occupied ...

	int up = (i==0 ? 0 : matrix[i-1][j]);    //  look up  
	int left = (j==0 ? 0 : matrix[i][j-1]);  //  look left
	
	switch (!!up + !!left) {
	  
	case 0:
	  matrix[i][j] = uf_make_set();      // a new cluster
	  break;
	  
	case 1:                              // part of an existing cluster
	  matrix[i][j] = max(up,left);       // whichever is nonzero is labelled
	  break;
	  
	case 2:                              // this site binds two clusters
	  matrix[i][j] = uf_union(up, left);
	  break;
	}
	
      }
  
  /* apply the relabeling to the matrix */

  /* This is a little bit sneaky.. we create a mapping from the canonical labels
     determined by union/find into a new set of canonical labels, which are 
     guaranteed to be sequential. */
  
  int *new_labels = calloc(sizeof(int), n_labels); // allocate array, initialized to zero
  
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {
	int x = uf_find(matrix[i][j]);
	if (new_labels[x] == 0) {
	  new_labels[0]++;
	  new_labels[x] = new_labels[0];
	}
	matrix[i][j] = new_labels[x];
      }
 
  int total_clusters = new_labels[0];

  free(new_labels);
  uf_done();

  return total_clusters;
}

/* This procedure checks to see that any occupied neighbors of an occupied site
   have the same label. */

void check_labelling(int **matrix, int m, int n) {
  int N,S,E,W;
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (matrix[i][j]) {
	N = ( i==0 ? 0 : matrix[i-1][j] );
	S = ( i==m-1 ? 0 : matrix[i+1][j] );
	E = ( j==n-1 ? 0 : matrix[i][j+1] );
	W = ( j==0 ? 0 : matrix[i][j-1] );
	
	assert( N==0 || matrix[i][j]==N );
	assert( S==0 || matrix[i][j]==S );
	assert( E==0 || matrix[i][j]==E );
	assert( W==0 || matrix[i][j]==W );
      }
}


int main(int argc, char **argv) {

  int i,j,idx,n;
  int ip,im,jp,jm;
  int labh;
  int labE,labN,labW,labS;
  int idum;
  int LX,LY;
  int time;
  int clusters;
  int **matrix;
  double *w;
  double threshold;
  //  int *area,*perimeter;
  double *area;
  double *xcm,*ycm;
  double *r2mean;
  double *xmean,*ymean;
  char path[128],fname[128];
  FILE *fin,*fout,*fcheck;
  double xm,ym,r2m,Rg;
  double ddum;
  
  if(argc!=6) {
    fprintf(stderr,"Error! Usage is: %s <input path> <size x> <size y> <time> <threshold> \n",argv[0]);
    exit(2);
  }

  strcpy(path,argv[1]);
  LX        = atoi(argv[2]);
  LY        = atoi(argv[3]);
  time      = atoi(argv[4]);
  threshold = atof(argv[5]);
  
  matrix = (int **)calloc(LX, sizeof(int*));
  w = (double *)malloc(sizeof(double)*LX*LY);

  /* initialize the matrix */
  
  sprintf(fname,"%s/h_w_t%d.dat",path,time);
  fin = fopen(fname,"r");
  if(fin==NULL) {
    fprintf(stderr,"Error! File %s not found!!\n",fname);
    exit(2);
  }
  
  for(i=0;i<LX;i++) {
   for(j=0;j<LY;j++) {
     idx = j + i * LY;
     fscanf(fin," %d %d %lf %lf \n",&idum,&idum,&ddum,&w[idx]); 
   }
  }
  fclose(fin);
  
     for (i=0; i<LX; i++) {
      matrix[i] = (int *)calloc(LY, sizeof(int));
      for (j=0; j<LY; j++) {
        idx = j + i * LY;
	if(w[idx]>threshold) {
 	 matrix[i][j] = 1;
	} else {
	 matrix[i][j] = 0;
	}	
      }
     }

     free(w);
     
    /* Process the matrix */
    
    clusters = hoshen_kopelman(matrix,LX,LY);

    //    perimeter = (int *)malloc(sizeof(int)*(clusters+1));
    //    area      = (int *)malloc(sizeof(int)*(clusters+1));
    area      = (double *)malloc(sizeof(double)*(clusters+1));

    xcm      = (double *)malloc(sizeof(double)*(clusters+1));
    ycm      = (double *)malloc(sizeof(double)*(clusters+1));
  
    xmean      = (double *)malloc(sizeof(double)*(clusters+1));
    ymean      = (double *)malloc(sizeof(double)*(clusters+1));

    r2mean      = (double *)malloc(sizeof(double)*(clusters+1));

    /*
    for(n=0;n<=clusters;n++) {
      area[n] = 0.0;
      perimeter[n] = 0;
    }
    */
    
    check_labelling(matrix,LX,LY);

    fcheck = fopen("check.dat","w");
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       fprintf(fcheck," %d %d %d \n",i,j,matrix[i][j]);
     }
     fprintf(fcheck,"\n");
    }
    fclose(fcheck);
    
    printf("HK reports %d clusters found\n",clusters);

    /* compute clusters area and perimeter */

    /*
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       
       ip = (i + 1 + LX) % LX; 
       im = (i - 1 + LX) % LX; 
       jp = (j + 1 + LY) % LY; 
       jm = (j - 1 + LY) % LY; 

       labh = matrix[i][j];
       labE = matrix[ip][j];
       labN = matrix[i][jp];
       labW = matrix[im][j];
       labS = matrix[i][jm];

       if(labh!=0) {
        if((labE!=labh)||(labN!=labh)||(labW!=labh)||(labS!=labh)) {
 	  ++perimeter[labh];
        }
        ++area[labh];
       }
       
     }
    }
    */

    for(n=0;n<=clusters;n++) {
      area[n] = 0.0;
      xcm[n] = 0.0;
      ycm[n] = 0.0;
      xmean[n] = 0.0;
      ymean[n] = 0.0;
      r2mean[n] = 0.0;
    }
    
    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       n = matrix[i][j];
       xcm[n] += (double)i;
       ycm[n] += (double)j;
       area[n] += 1.0;
     }
    }

    for(i=0;i<LX;i++) {
     for(j=0;j<LY;j++) {
       n = matrix[i][j];
       xmean[n] += ((double)i - xcm[n]);   
       ymean[n] += ((double)j - ycm[n]);   
       r2mean[n] += ((double)i - xcm[n])*((double)i - xcm[n]) + ((double)j - ycm[n])*((double)j - ycm[n]);
     }
    }

    sprintf(fname,"clusters_t%d.dat",time);
    fout = fopen(fname,"w");
    for(n=0;n<=clusters;n++) {

      r2m = r2mean[n]/area[n];
      xm = xmean[n]/area[n];
      ym = ymean[n]/area[n];
      
      Rg = sqrt(r2m - (xm*xm + ym*ym));
      fprintf(fout,"%d %g %g \n",n,area[n],Rg);

    }
    fclose(fout);

    for(n=0;n<clusters;n++) {
      if(area[n]>amax) {
	amax = area[n];
      }
    }

    delta = amax/((double)BIN);

    for(n=0;n<clusters;n++) {
      nc = (int)floor(area[n]/delta);
      pdf[nc] += 1.0;
      atot += area[n];
    }

    
    /* print out cluster informations */

    /*
    sprintf(fname,"clusters_t%d.dat",time);
    fout = fopen(fname,"w");
    for(n=1; n<=clusters; n++) {
      fprintf(fout," %d %d %d \n",n,area[n],perimeter[n]);
    }
    fclose(fout);
    */
    
    for (i=0; i<LX; i++)
      free(matrix[i]);
    free(matrix);
  
    free(area);
    //    free(perimeter);
    free(xmean);
    free(ymean);
    free(r2mean);
    free(xcm);
    free(ycm);
    
  return 0;
  
}


