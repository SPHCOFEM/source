#undef __FUNC__
#define __FUNC__ "matrix_transpose"

void matrix_transpose(int n,
  double *A,double *B)     
{
  int i,j;
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
	    B[i+n*j]=A[j+n*i];
	  }
  }
}

#undef __FUNC__
#define __FUNC__ "matrices_addition"

void matrices_addition(int n,
  double *A,double *B,double *C)     
{
  int i,j;
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
  	  C[i+n*j]=A[i+n*j]+B[i+n*j];
  	}
  }
}

#undef __FUNC__
#define __FUNC__ "matrices_multiplication"

void matrix_multiplication(int n,
  double *A,double c,double *C)     
{
  int i,j;
  
  for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
	  {
	    A[i+n*j]=c*A[i+n*j];
	  }
  }
}

#undef __FUNC__
#define __FUNC__ "matrices_product"

void matrices_product(int n,
  double *A,double *B,double *C)
{
  int i,j,k;
  
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
}
