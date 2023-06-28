#include "header.h"

#undef __FUNC__
#define __FUNC__ "mtranspose"

double *mtranspose(int n,double *A)
{
  int i,j;
  double *C=0;
  C=createMemMore(double,sqr(n));
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
	    C[i+n*j]=A[j+n*i];
	  }
  }
  return(C);
}

#undef __FUNC__
#define __FUNC__ "maddition"

void *maddition(int n,double *A,double *B)
{
  int i,j;
  double *C=0;
  C=createMemMore(double,sqr(n));

  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
  	  C[i+n*j]=A[i+n*j]+B[i+n*j];
  	}
  }
  return(C);
}

#undef __FUNC__
#define __FUNC__ "msubtraction"

void *msubtraction(int n,double *A,double *B)
{
  int i,j;
  double *C=0;
  C=createMemMore(double,sqr(n));
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
  	  C[i+n*j]=A[i+n*j]-B[i+n*j];
  	}
  }
  return(C);
}

#undef __FUNC__
#define __FUNC__ "mmultiplication"

double *mmultiplication(int n,double *A,double b)
{
  int i,j;
  double *C=0; //double C[sqr(n)];
  C=createMemMore(double,sqr(n));
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
      C[i+n*j]=A[i+n*j]*b;
	  }
  }
  return(C);
}

#undef __FUNC__
#define __FUNC__ "mproduct"

void *mproduct(int n,double *A,double *B)
{
  int i,j,k;
  double *C=0;
  C=createMemMore(double,n*n);
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
	    C[i+n*j]=0.0;
	    for (k=0;k<n;k++)
	      {
	        C[i+n*j]=C[i+n*j]+A[i+n*k]*B[k+n*j];
	      }
	  }
  }
  return(C);
}
