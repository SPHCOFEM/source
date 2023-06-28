#include "header.h"
#include "matrices_operations.h"

#undef __FUNC__
#define __FUNC__ "fem_tet"

void fem_tetrahedron(double *Me, double *Ke,
		     
  double x1,double x2,double x3,double x4,
  double y1,double y2,double y3,double y4,
  double z1,double z2,double z3,double z4,

  double rho,double E,double nu)
{
  double m,k;
  double *Si,*Sit,*I;
  double *aux;

  Si=createMemMore(double,144);
  Sit=createMemMore(double,144);
  I=createMemMore(double,144);
  aux=createMemMore(double,144);

  memset(Si,0,144*sizeof(double));
  memset(I,0,144*sizeof(double));

  m=rho*fabs(x2*y1*z3-x2*y3*z1-x1*y2*z3+x1*y3*z2+x4*y2*z3-x3*y1*z2+x3*y2*z1-x1*y3*z4+x1*y4*z3+x3*y1*z4-x3*y4*z1-x4*y1*z3+x4*y3*z1+x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2-x4*y3*z2+x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1);
  k=E/(1.0-2.0*nu)/(1.0+nu)*fabs(x2*y1*z3-x2*y3*z1-x1*y2*z3+x1*y3*z2+x4*y2*z3-x3*y1*z2+x3*y2*z1-x1*y3*z4+x1*y4*z3+x3*y1*z4-x3*y4*z1-x4*y1*z3+x4*y3*z1+x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2-x4*y3*z2+x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1);

  Si[0+12*0] = -(x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[0+12*1] = (x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[0+12*2] = -(x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[0+12*3] = (x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[1+12*0] = (y3*z4-y4*z3-y2*z4+y4*z2+y2*z3-y3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[1+12*1] = (-y3*z4+y4*z3+y1*z4-y4*z1-y1*z3+y3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[1+12*2] = -(-y2*z4+y4*z2+y1*z4-y4*z1-y1*z2+y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[1+12*3] = (y1*z3+y2*z1+y3*z2-y1*z2-y3*z1-y2*z3)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[2+12*0] = -(x3*z4-x4*z3-x2*z4+x4*z2+x2*z3-x3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[2+12*1] = -(-x3*z4+x4*z3+x1*z4-x4*z1-x1*z3+x3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[2+12*2] = (-x2*z4+x4*z2+x1*z4-x4*z1-x1*z2+x2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[2+12*3] = -(-x2*z3+x3*z2+x1*z3-x3*z1-x1*z2+x2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[3+12*0] = (x3*y4-x4*y3-x2*y4+x4*y2+x2*y3-x3*y2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[3+12*1] = -(x3*y4-x4*y3-x1*y4+x4*y1+x1*y3-x3*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[3+12*2] = (x2*y4-x4*y2-x1*y4+x4*y1+x1*y2-x2*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[3+12*3] = -(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[4+12*4] = -(x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[4+12*5] = (x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[4+12*6] = -(x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[4+12*7] = (x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[5+12*4] = (y3*z4-y4*z3-y2*z4+y4*z2+y2*z3-y3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[5+12*5] = (-y3*z4+y4*z3+y1*z4-y4*z1-y1*z3+y3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[5+12*6] = -(-y2*z4+y4*z2+y1*z4-y4*z1-y1*z2+y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[5+12*7] = (y1*z3+y2*z1+y3*z2-y1*z2-y3*z1-y2*z3)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[6+12*4] = -(x3*z4-x4*z3-x2*z4+x4*z2+x2*z3-x3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[6+12*5] = -(-x3*z4+x4*z3+x1*z4-x4*z1-x1*z3+x3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[6+12*6] = (-x2*z4+x4*z2+x1*z4-x4*z1-x1*z2+x2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[6+12*7] = -(-x2*z3+x3*z2+x1*z3-x3*z1-x1*z2+x2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[7+12*4] = (x3*y4-x4*y3-x2*y4+x4*y2+x2*y3-x3*y2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[7+12*5] = -(x3*y4-x4*y3-x1*y4+x4*y1+x1*y3-x3*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[7+12*6] = (x2*y4-x4*y2-x1*y4+x4*y1+x1*y2-x2*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[7+12*7] = -(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[8+12*8] = -(x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[8+12*9] = (x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[8+12*10] = -(x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[8+12*11] = (x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[9+12*8] = (y3*z4-y4*z3-y2*z4+y4*z2+y2*z3-y3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[9+12*9] = (-y3*z4+y4*z3+y1*z4-y4*z1-y1*z3+y3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[9+12*10] = -(-y2*z4+y4*z2+y1*z4-y4*z1-y1*z2+y2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[9+12*11] = (y1*z3+y2*z1+y3*z2-y1*z2-y3*z1-y2*z3)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[10+12*8] = -(x3*z4-x4*z3-x2*z4+x4*z2+x2*z3-x3*z2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[10+12*9] = -(-x3*z4+x4*z3+x1*z4-x4*z1-x1*z3+x3*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[10+12*10] = (-x2*z4+x4*z2+x1*z4-x4*z1-x1*z2+x2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[10+12*11] = -(-x2*z3+x3*z2+x1*z3-x3*z1-x1*z2+x2*z1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[11+12*8] = (x3*y4-x4*y3-x2*y4+x4*y2+x2*y3-x3*y2)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[11+12*9] = -(x3*y4-x4*y3-x1*y4+x4*y1+x1*y3-x3*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[11+12*10] = (x2*y4-x4*y2-x1*y4+x4*y1+x1*y2-x2*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);
  Si[11+12*11] = -(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)/(-x2*y3*z4+x2*y4*z3+x3*y2*z4-x3*y4*z2-x4*y2*z3+x4*y3*z2+x1*y3*z4-x1*y4*z3-x3*y1*z4+x3*y4*z1+x4*y1*z3-x4*y3*z1-x1*y2*z4+x1*y4*z2+x2*y1*z4-x2*y4*z1-x4*y1*z2+x4*y2*z1+x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1);

  I[0+12*0] = 1.0/6.0;
  I[0+12*1] = x1/24+x4/24+x3/24+x2/24;
  I[0+12*2] = x1/6-y1/8+y2/24+y3/24+y4/24;
  I[0+12*3] = z1/24+z4/24+z3/24+z2/24;
  I[1+12*0] = x1/24+x4/24+x3/24+x2/24;
  I[1+12*1] = x2*x1/60+x4*x3/60+x4*x2/60+x2*x3/60+x3*x3/60+x2*x2/60+x1*x1/60+x4*x1/60+x4*x4/60+x3*x1/60;
  I[1+12*2] = x2*y3/120+x2*x1/24+x4*y2/120+x1*y2/120+x3*y2/120+x2*y2/60+x2*y4/120-x2*y1/30+x1*y4/120+x4*y4/60-x4*y1/30-x1*y1/40+x1*x1/24+x4*x1/24+x4*y3/120+x1*y3/120+x3*y3/60+x3*x1/24+x3*y4/120-x3*y1/30;
  I[1+12*3] = x4*z4/60+x4*z1/120+x1*z4/120+x1*z1/60+x3*z2/120+x4*z2/120+x2*z4/120+x2*z1/120+x1*z2/120+x2*z3/120+x2*z2/60+x3*z4/120+x3*z3/60+x1*z3/120+x3*z1/120+x4*z3/120;
  I[2+12*0] = x1/6-y1/8+y2/24+y3/24+y4/24;
  I[2+12*1] = x2*y3/120+x2*x1/24+x4*y2/120+x1*y2/120+x3*y2/120+x2*y2/60+x2*y4/120-x2*y1/30+x1*y4/120+x4*y4/60-x4*y1/30-x1*y1/40+x1*x1/24+x4*x1/24+x4*y3/120+x1*y3/120+x3*y3/60+x3*x1/24+x3*y4/120-x3*y1/30;
  I[2+12*2] = x1*y2/12+x1*y4/12-x1*y1/4+x1*x1/6+y3*y4/60+y4*y4/60+y1*y1/10+y2*y2/60+y3*y3/60-y1*y3/15-y1*y4/15+x1*y3/12+y4*y2/60+y2*y3/60-y1*y2/15;
  I[2+12*3] = x1*z4/24+x1*z1/24+y3*z1/120+y3*z3/60+y4*z4/60+y4*z1/120-y1*z4/30-y1*z1/40+y3*z4/120+y4*z3/120-y1*z3/30+y2*z1/120+y3*z2/120+y2*z3/120+y2*z4/120+y2*z2/60+x1*z2/24-y1*z2/30+y4*z2/120+x1*z3/24;
  I[3+12*0] = z1/24+z4/24+z3/24+z2/24;
  I[3+12*1] = x4*z4/60+x4*z1/120+x1*z4/120+x1*z1/60+x3*z2/120+x4*z2/120+x2*z4/120+x2*z1/120+x1*z2/120+x2*z3/120+x2*z2/60+x3*z4/120+x3*z3/60+x1*z3/120+x3*z1/120+x4*z3/120;
  I[3+12*2] = x1*z4/24+x1*z1/24+y3*z1/120+y3*z3/60+y4*z4/60+y4*z1/120-y1*z4/30-y1*z1/40+y3*z4/120+y4*z3/120-y1*z3/30+y2*z1/120+y3*z2/120+y2*z3/120+y2*z4/120+y2*z2/60+x1*z2/24-y1*z2/30+y4*z2/120+x1*z3/24;
  I[3+12*3] = z2*z2/60+z3*z3/60+z4*z4/60+z1*z1/60+z1*z3/60+z1*z4/60+z4*z3/60+z1*z2/60+z4*z2/60+z2*z3/60;
  I[4+12*4] = 1.0/6.0;
  I[4+12*5] = x1/24+x4/24+x3/24+x2/24;
  I[4+12*6] = x1/6-y1/8+y2/24+y3/24+y4/24;
  I[4+12*7] = z1/24+z4/24+z3/24+z2/24;
  I[5+12*4] = x1/24+x4/24+x3/24+x2/24;
  I[5+12*5] = x2*x1/60+x4*x3/60+x4*x2/60+x2*x3/60+x3*x3/60+x2*x2/60+x1*x1/60+x4*x1/60+x4*x4/60+x3*x1/60;
  I[5+12*6] = x2*y3/120+x2*x1/24+x4*y2/120+x1*y2/120+x3*y2/120+x2*y2/60+x2*y4/120-x2*y1/30+x1*y4/120+x4*y4/60-x4*y1/30-x1*y1/40+x1*x1/24+x4*x1/24+x4*y3/120+x1*y3/120+x3*y3/60+x3*x1/24+x3*y4/120-x3*y1/30;
  I[5+12*7] = x4*z4/60+x4*z1/120+x1*z4/120+x1*z1/60+x3*z2/120+x4*z2/120+x2*z4/120+x2*z1/120+x1*z2/120+x2*z3/120+x2*z2/60+x3*z4/120+x3*z3/60+x1*z3/120+x3*z1/120+x4*z3/120;
  I[6+12*4] = x1/6-y1/8+y2/24+y3/24+y4/24;
  I[6+12*5] = x2*y3/120+x2*x1/24+x4*y2/120+x1*y2/120+x3*y2/120+x2*y2/60+x2*y4/120-x2*y1/30+x1*y4/120+x4*y4/60-x4*y1/30-x1*y1/40+x1*x1/24+x4*x1/24+x4*y3/120+x1*y3/120+x3*y3/60+x3*x1/24+x3*y4/120-x3*y1/30;
  I[6+12*6] = x1*y2/12+x1*y4/12-x1*y1/4+x1*x1/6+y3*y4/60+y4*y4/60+y1*y1/10+y2*y2/60+y3*y3/60-y1*y3/15-y1*y4/15+x1*y3/12+y4*y2/60+y2*y3/60-y1*y2/15;
  I[6+12*7] = x1*z4/24+x1*z1/24+y3*z1/120+y3*z3/60+y4*z4/60+y4*z1/120-y1*z4/30-y1*z1/40+y3*z4/120+y4*z3/120-y1*z3/30+y2*z1/120+y3*z2/120+y2*z3/120+y2*z4/120+y2*z2/60+x1*z2/24-y1*z2/30+y4*z2/120+x1*z3/24;
  I[7+12*4] = z1/24+z4/24+z3/24+z2/24;
  I[7+12*5] = x4*z4/60+x4*z1/120+x1*z4/120+x1*z1/60+x3*z2/120+x4*z2/120+x2*z4/120+x2*z1/120+x1*z2/120+x2*z3/120+x2*z2/60+x3*z4/120+x3*z3/60+x1*z3/120+x3*z1/120+x4*z3/120;
  I[7+12*6] = x1*z4/24+x1*z1/24+y3*z1/120+y3*z3/60+y4*z4/60+y4*z1/120-y1*z4/30-y1*z1/40+y3*z4/120+y4*z3/120-y1*z3/30+y2*z1/120+y3*z2/120+y2*z3/120+y2*z4/120+y2*z2/60+x1*z2/24-y1*z2/30+y4*z2/120+x1*z3/24;
  I[7+12*7] = z2*z2/60+z3*z3/60+z4*z4/60+z1*z1/60+z1*z3/60+z1*z4/60+z4*z3/60+z1*z2/60+z4*z2/60+z2*z3/60;
  I[8+12*8] = 1.0/6.0;
  I[8+12*9] = x1/24+x4/24+x3/24+x2/24;
  I[8+12*10] = x1/6-y1/8+y2/24+y3/24+y4/24;
  I[8+12*11] = z1/24+z4/24+z3/24+z2/24;
  I[9+12*8] = x1/24+x4/24+x3/24+x2/24;
  I[9+12*9] = x2*x1/60+x4*x3/60+x4*x2/60+x2*x3/60+x3*x3/60+x2*x2/60+x1*x1/60+x4*x1/60+x4*x4/60+x3*x1/60;
  I[9+12*10] = x2*y3/120+x2*x1/24+x4*y2/120+x1*y2/120+x3*y2/120+x2*y2/60+x2*y4/120-x2*y1/30+x1*y4/120+x4*y4/60-x4*y1/30-x1*y1/40+x1*x1/24+x4*x1/24+x4*y3/120+x1*y3/120+x3*y3/60+x3*x1/24+x3*y4/120-x3*y1/30;
  I[9+12*11] = x4*z4/60+x4*z1/120+x1*z4/120+x1*z1/60+x3*z2/120+x4*z2/120+x2*z4/120+x2*z1/120+x1*z2/120+x2*z3/120+x2*z2/60+x3*z4/120+x3*z3/60+x1*z3/120+x3*z1/120+x4*z3/120;
  I[10+12*8] = x1/6-y1/8+y2/24+y3/24+y4/24;
  I[10+12*9] = x2*y3/120+x2*x1/24+x4*y2/120+x1*y2/120+x3*y2/120+x2*y2/60+x2*y4/120-x2*y1/30+x1*y4/120+x4*y4/60-x4*y1/30-x1*y1/40+x1*x1/24+x4*x1/24+x4*y3/120+x1*y3/120+x3*y3/60+x3*x1/24+x3*y4/120-x3*y1/30;
  I[10+12*10] = x1*y2/12+x1*y4/12-x1*y1/4+x1*x1/6+y3*y4/60+y4*y4/60+y1*y1/10+y2*y2/60+y3*y3/60-y1*y3/15-y1*y4/15+x1*y3/12+y4*y2/60+y2*y3/60-y1*y2/15;
  I[10+12*11] = x1*z4/24+x1*z1/24+y3*z1/120+y3*z3/60+y4*z4/60+y4*z1/120-y1*z4/30-y1*z1/40+y3*z4/120+y4*z3/120-y1*z3/30+y2*z1/120+y3*z2/120+y2*z3/120+y2*z4/120+y2*z2/60+x1*z2/24-y1*z2/30+y4*z2/120+x1*z3/24;
  I[11+12*8] = z1/24+z4/24+z3/24+z2/24;
  I[11+12*9] = x4*z4/60+x4*z1/120+x1*z4/120+x1*z1/60+x3*z2/120+x4*z2/120+x2*z4/120+x2*z1/120+x1*z2/120+x2*z3/120+x2*z2/60+x3*z4/120+x3*z3/60+x1*z3/120+x3*z1/120+x4*z3/120;
  I[11+12*10] = x1*z4/24+x1*z1/24+y3*z1/120+y3*z3/60+y4*z4/60+y4*z1/120-y1*z4/30-y1*z1/40+y3*z4/120+y4*z3/120-y1*z3/30+y2*z1/120+y3*z2/120+y2*z3/120+y2*z4/120+y2*z2/60+x1*z2/24-y1*z2/30+y4*z2/120+x1*z3/24;
  I[11+12*11] = z2*z2/60+z3*z3/60+z4*z4/60+z1*z1/60+z1*z3/60+z1*z4/60+z4*z3/60+z1*z2/60+z4*z2/60+z2*z3/60;
  
  matrix_transpose(12,
		Si,
		Sit);

  matrix_multiplication(12,
		I,m,
		I);
  
  matrices_product(12,
		Sit,I,
		aux);
  
  matrices_product(12,
		aux,Si,
		Me);

  memset(I,0,144*sizeof(double));

  I[1+12*1] = 1.0/6.0-nu/6;
  I[1+12*6] = nu/6;
  I[1+12*11] = nu/6;
  I[2+12*2] = 1.0/12.0-nu/6;
  I[2+12*5] = 1.0/12.0-nu/6;
  I[3+12*3] = 1.0/12.0-nu/6;
  I[3+12*9] = 1.0/12.0-nu/6;
  I[5+12*2] = 1.0/12.0-nu/6;
  I[5+12*5] = 1.0/12.0-nu/6;
  I[6+12*1] = nu/6;
  I[6+12*6] = 1.0/6.0-nu/6;
  I[6+12*11] = nu/6;
  I[7+12*7] = 1.0/12.0-nu/6;
  I[7+12*10] = 1.0/12.0-nu/6;
  I[9+12*3] = 1.0/12.0-nu/6;
  I[9+12*9] = 1.0/12.0-nu/6;
  I[10+12*7] = 1.0/12.0-nu/6;
  I[10+12*10] = 1.0/12.0-nu/6;
  I[11+12*1] = nu/6;
  I[11+12*6] = nu/6;
  I[11+12*11] = 1.0/6.0-nu/6;
  
  matrix_multiplication(12,
		I,k,
		I);

  matrices_product(12,
		Sit,I,
		aux);
  
  matrices_product(12,
		aux,Si,
		Ke);

  freeMem(Si);
  freeMem(Sit);
  freeMem(I);
  freeMem(aux);
}
