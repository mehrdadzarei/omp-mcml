/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Some routines modified from Numerical Recipes in C,
 *	including error report, array or matrix declaration
 *	and releasing.
 ****/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/***********************************************************
 *	Report error message to stderr, then exit the program
 *	with signal 1.
 ****/
void nrerror(char error_text[])
     
{
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

/***********************************************************
 *	Allocate an array with index from nl to nh inclusive.
 *
 *	Original matrix and vector from Numerical Recipes in C
 *	don't initialize the elements to zero. This will
 *	be accomplished by the following functions. 
 ****/
double *AllocVector(short nl, short nh)
{
  double *v;
  short i;
  
  v=(double *)calloc((unsigned) (nh-nl+1), sizeof(double));
  if (!v) nrerror("allocation failure in vector()");
  
  v -= nl;
  return v;
}

/***********************************************************
 *	Allocate a matrix with row index from nrl to nrh 
 *	inclusive, and column index from ncl to nch
 *	inclusive.
 ****/
double **AllocMatrix(short nrl,short nrh,
					 short ncl,short nch)
{
  short i,j;
  double **m;
  
  m=(double **) calloc((unsigned) (nrh-nrl+1), sizeof(double*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
    if (!m[i]) nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  
  return m;
}

double ***AllocCubic(short nrl,short nrh, short ncl,short nch, short nzl,short nzh){
  short i,j,k;
  double ***m;
  
  m=(double ***) calloc((unsigned) (nrh-nrl+1), sizeof(double**));
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double **) calloc((unsigned) (nch-ncl+1), sizeof(double*));
    m[i] -= ncl;
  }
  
  for(i=nrl;i<=nrh;i++) 
	for(j=ncl;j<=nch;j++) {
        m[i][j]=(double *) calloc((unsigned) (nzh-nzl+1), sizeof(double));
        m[i][j] -= nzl;
	}
  
  return m;
}

/***********************************************************
 *	Release the memory.
 ****/
void FreeVector(double *v,short nl,short nh)
{
  free((char*) (v+nl));
}

/***********************************************************
 *	Release the memory.
 ****/
void FreeMatrix(double **m,short nrl,short nrh,
			    short ncl,short nch)
{
  short i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void FreeCubic(double ***m,short nxl,short nxh,short nyl,short nyh, short ncl,short nch){
  short i,j;
  
  for(i=nxh;i>=nxl;i--){
    for(j=nyh;j>=nyl;j--) 
	  free((char*) (m[i][j]+ncl));
	free((char*) (m[i]+nyl));
  }
  free((char*) (m+nxl));
}
