#include "header.h"
#include "matrices_operations.h"

#undef __FUNC__
#define __FUNC__ "fem_triangle"

void fem_triangle(double *Me, double *Ke,

  double x1,double x2,double x3,double y1,double y2,double y3,

  double rho,double t,double E,double nu)
{
  double m,k;
  double *Si,*Sit,*I;
  double *aux;

  Si=createMemMore(double,36);
  Sit=createMemMore(double,36);
  I=createMemMore(double,36);
  aux=createMemMore(double,36);

  m=rho*t*fabs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  k=t*E/(1.0-sqr(nu))*fabs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);

  memset(Si,0,36*sizeof(double));
  memset(I,0,36*sizeof(double));

  Si[0+6*0] = (x2*y3-x3*y2)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[0+6*1] = -(x1*y3-x3*y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[0+6*2] = (x1*y2-x2*y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[1+6*0] = (-y3+y2)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[1+6*1] = -(-y3+y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[1+6*2] = (-y2+y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[2+6*0] = -(-x3+x2)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[2+6*1] = (-x3+x1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[2+6*2] = -(-x2+x1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[3+6*3] = (x2*y3-x3*y2)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[3+6*4] = -(x1*y3-x3*y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[3+6*5] = (x1*y2-x2*y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[4+6*3] = (-y3+y2)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[4+6*4] = -(-y3+y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[4+6*5] = (-y2+y1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[5+6*3] = -(-x3+x2)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[5+6*4] = (-x3+x1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  Si[5+6*5] = -(-x2+x1)/(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);

  I[0+6*0] = 1.0/2.0;
  I[0+6*1] = x1/6+x3/6+x2/6;
  I[0+6*2] = x1/2+y3/6-y1/3+y2/6;
  I[1+6*0] = x1/6+x3/6+x2/6;
  I[1+6*1] = x1*x1/12+x2*x1/12+x3*x1/12+x2*x2/12+x3*x2/12+x3*x3/12;
  I[1+6*2] = x3*y3/12-x3*y1/8+x1*y3/24-x1*y1/12+x3*x1/6+x1*x1/6+x2*y2/12-x2*y1/8+x2*y3/24+x3*y2/24+x1*y2/24+x2*x1/6;
  I[2+6*0] = x1/2+y3/6-y1/3+y2/6;
  I[2+6*1] = x3*y3/12-x3*y1/8+x1*y3/24-x1*y1/12+x3*x1/6+x1*x1/6+x2*y2/12-x2*y1/8+x2*y3/24+x3*y2/24+x1*y2/24+x2*x1/6;
  I[2+6*2] = y1*y1/4-2.0/3.0*x1*y1-y1*y3/4-y2*y1/4+x1*x1/2+x1*y3/3+y3*y3/12+x1*y2/3+y2*y3/12+y2*y2/12;
  I[3+6*3] = 1.0/2.0;
  I[3+6*4] = x1/6+x3/6+x2/6;
  I[3+6*5] = x1/2+y3/6-y1/3+y2/6;
  I[4+6*3] = x1/6+x3/6+x2/6;
  I[4+6*4] = x1*x1/12+x2*x1/12+x3*x1/12+x2*x2/12+x3*x2/12+x3*x3/12;
  I[4+6*5] = x3*y3/12-x3*y1/8+x1*y3/24-x1*y1/12+x3*x1/6+x1*x1/6+x2*y2/12-x2*y1/8+x2*y3/24+x3*y2/24+x1*y2/24+x2*x1/6;
  I[5+6*3] = x1/2+y3/6-y1/3+y2/6;
  I[5+6*4] = x3*y3/12-x3*y1/8+x1*y3/24-x1*y1/12+x3*x1/6+x1*x1/6+x2*y2/12-x2*y1/8+x2*y3/24+x3*y2/24+x1*y2/24+x2*x1/6;
  I[5+6*5] = y1*y1/4-2.0/3.0*x1*y1-y1*y3/4-y2*y1/4+x1*x1/2+x1*y3/3+y3*y3/12+x1*y2/3+y2*y3/12+y2*y2/12;

  matrix_transpose(6,
		Si,
		Sit);

  matrix_multiplication(6,
		I,m,
		I);
  
  matrices_product(6,
		Sit,I,
		aux);
  
  matrices_product(6,
		aux,Si,
		Me);

  memset(I,0,36*sizeof(double));

  I[1+6*1] = 1.0/2.0;
  I[2+6*2] = 1.0/2.0;
  I[2+6*4] = 1.0/2.0;
  I[4+6*2] = 1.0/2.0;
  I[4+6*4] = 1.0/2.0;
  I[5+6*5] = 1.0/2.0;
	  
  matrix_multiplication(6,
		I,k,
		I);
  
  matrices_product(6,
		Sit,I,
		aux);
  
  matrices_product(6,
		aux,Si,
		Ke);
  
  freeMem(Si);
  freeMem(Sit);
  freeMem(I);
  freeMem(aux);
}
