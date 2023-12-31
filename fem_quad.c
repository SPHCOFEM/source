#include "header.h"
#include "matrices_operations.h"

#undef __FUNC__
#define __FUNC__ "fem_quad"

void fem_quad(double *Me, double *Ke,

  double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4,

  double rho,double t,double E,double nu)
{
  double ma,mb,ka,kb;
  double *Si,*Sit,*Ia,*Ib,*I;
  double *aux;

  Si=createMemMore(double,64);
  Sit=createMemMore(double,64);
  Ia=createMemMore(double,64);
  Ib=createMemMore(double,64);
  I=createMemMore(double,64);
  aux=createMemMore(double,64);
  
  ma=rho*t*fabs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  mb=rho*t*fabs(x3*y4-x4*y3-x1*y4+x4*y1+x1*y3-x3*y1);
  ka=t*E/(1.0-nu*nu)*fabs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
  kb=t*E/(1.0-nu*nu)*fabs(x3*y4-x4*y3-x1*y4+x4*y1+x1*y3-x3*y1);
  
  memset(Si,0,64*sizeof(double));
  memset(Ia,0,64*sizeof(double));
  memset(Ib,0,64*sizeof(double));

  Si[0+8*0] = -(x2*x4*y3*y4-x2*x3*y3*y4-x3*y2*x4*y4+x3*x2*y2*y4+x4*x3*y2*y3-x4*x2*y2*y3)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[0+8*1] = (x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[0+8*2] = -(x1*y2*x4*y4-x1*x2*y2*y4-x2*x4*y4*y1+x2*x1*y4*y1+x4*x2*y2*y1-x4*x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[0+8*3] = (x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[1+8*0] = -(-x4*y3*y4+x3*y3*y4+y2*x4*y4-x2*y2*y4-x3*y2*y3+x2*y2*y3)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[1+8*1] = (-x4*y3*y4+x3*y3*y4+x4*y4*y1-x1*y4*y1-x3*y3*y1+x1*y3*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[1+8*2] = -(-y2*x4*y4+x2*y2*y4+x4*y4*y1-x1*y4*y1-x2*y2*y1+x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[1+8*3] = (-x3*y2*y3+x2*y2*y3+x3*y3*y1-x1*y3*y1-x2*y2*y1+x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[2+8*0] = (-x4*x3*y4+x4*x3*y3+x2*x4*y4-x2*y2*x4-x3*x2*y3+x3*x2*y2)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[2+8*1] = -(-x4*x3*y4+x4*x3*y3+x4*x1*y4-x4*x1*y1-x3*x1*y3+x3*x1*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[2+8*2] = (-x2*x4*y4+x2*y2*x4+x4*x1*y4-x4*x1*y1-x1*x2*y2+x1*x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[2+8*3] = -(-x3*x2*y3+x3*x2*y2+x3*x1*y3-x3*x1*y1-x1*x2*y2+x1*x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[3+8*0] = (x3*y4-x4*y3-x2*y4+y2*x4+x2*y3-x3*y2)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[3+8*1] = -(x3*y4-x4*y3-x1*y4+x4*y1+x1*y3-x3*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[3+8*2] = (x2*y4-y2*x4-x1*y4+x4*y1+x1*y2-x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[3+8*3] = -(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[4+8*4] = -(x2*x4*y3*y4-x2*x3*y3*y4-x3*y2*x4*y4+x3*x2*y2*y4+x4*x3*y2*y3-x4*x2*y2*y3)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);          Si[4+8*5] = (x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[4+8*6] = -(x1*y2*x4*y4-x1*x2*y2*y4-x2*x4*y4*y1+x2*x1*y4*y1+x4*x2*y2*y1-x4*x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[4+8*7] = (x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[5+8*4] = -(-x4*y3*y4+x3*y3*y4+y2*x4*y4-x2*y2*y4-x3*y2*y3+x2*y2*y3)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[5+8*5] = (-x4*y3*y4+x3*y3*y4+x4*y4*y1-x1*y4*y1-x3*y3*y1+x1*y3*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[5+8*6] = -(-y2*x4*y4+x2*y2*y4+x4*y4*y1-x1*y4*y1-x2*y2*y1+x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[5+8*7] = (-x3*y2*y3+x2*y2*y3+x3*y3*y1-x1*y3*y1-x2*y2*y1+x1*y2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[6+8*4] = (-x4*x3*y4+x4*x3*y3+x2*x4*y4-x2*y2*x4-x3*x2*y3+x3*x2*y2)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);          Si[6+8*5] = -(-x4*x3*y4+x4*x3*y3+x4*x1*y4-x4*x1*y1-x3*x1*y3+x3*x1*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[6+8*6] = (-x2*x4*y4+x2*y2*x4+x4*x1*y4-x4*x1*y1-x1*x2*y2+x1*x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[6+8*7] = -(-x3*x2*y3+x3*x2*y2+x3*x1*y3-x3*x1*y1-x1*x2*y2+x1*x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[7+8*4] = (x3*y4-x4*y3-x2*y4+y2*x4+x2*y3-x3*y2)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[7+8*5] = -(x3*y4-x4*y3-x1*y4+x4*y1+x1*y3-x3*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[7+8*6] = (x2*y4-y2*x4-x1*y4+x4*y1+x1*y2-x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);
  Si[7+8*7] = -(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)/(-x2*x4*y3*y4+x2*x3*y3*y4+x3*y2*x4*y4-x3*x2*y2*y4-x4*x3*y2*y3+x4*x2*y2*y3+x4*x1*y3*y4-x3*x1*y3*y4-x3*x4*y4*y1+x3*x1*y4*y1+x3*x4*y3*y1-x4*x1*y3*y1-x1*y2*x4*y4+x1*x2*y2*y4+x2*x4*y4*y1-x2*x1*y4*y1-x4*x2*y2*y1+x4*x1*y2*y1+x3*x1*y2*y3-x2*x1*y2*y3-x2*x3*y3*y1+x2*x1*y3*y1+x2*x3*y2*y1-x3*x1*y2*y1);

  Ia[0+8*0] = 1.0/2.0;
  Ia[0+8*1] = x1/6+x3/6+x2/6;
  Ia[0+8*2] = x1/2+y3/6-y1/3+y2/6;
  Ia[0+8*3] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[1+8*0] = x1/6+x3/6+x2/6;
  Ia[1+8*1] = x1*x1/12+x1*x2/12+x3*x1/12+x2*x2/12+x3*x2/12+x3*x3/12;
  Ia[1+8*2] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[1+8*3] = x3*x2*y2/30+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x2*x2*y2/20+x3*x2*x1/12+x3*x2*y3/30-x3*x2*y1/15+x2*x2*y3/60-x2*x2*y1/15+x2*x2*x1/12+x3*x3*y2/60+x1*x1*y2/60+x1*x2*y3/60+x1*x2*y2/30-x1*x2*y1/20+x2*x1*x1/12+x3*x1*y2/60;
  Ia[2+8*0] = x1/2+y3/6-y1/3+y2/6;
  Ia[2+8*1] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[2+8*2] = y1*y1/4-2.0/3.0*x1*y1-y2*y1/4-y1*y3/4+x1*y2/3+y2*y2/12+y2*y3/12+x1*y3/3+x1*x1/2+y3*y3/12;
  Ia[2+8*3] = -2.0/15.0*x2*y2*y1-x1*y3*y1/20+x2*y1*y1/10+x1*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x2*y3*y3/60+x1*y2*y2/60+x3*y2*y2/60+x2*y2*y2/20+x1*x1*y2/12+x2*y2*y3/30+x1*y2*y3/60+x1*x2*y3/12+x1*x2*y2/6-x3*y2*y1/15-x2*y3*y1/15+x3*y2*y3/30-x1*x2*y1/4+x2*x1*x1/6+x3*x1*y2/12-x1*y2*y1/20;
  Ia[3+8*0] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[3+8*1] = x3*x2*y2/30+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x2*x2*y2/20+x3*x2*x1/12+x3*x2*y3/30-x3*x2*y1/15+x2*x2*y3/60-x2*x2*y1/15+x2*x2*x1/12+x3*x3*y2/60+x1*x1*y2/60+x1*x2*y3/60+x1*x2*y2/30-x1*x2*y1/20+x2*x1*x1/12+x3*x1*y2/60;
  Ia[3+8*2] = -2.0/15.0*x2*y2*y1-x1*y3*y1/20+x2*y1*y1/10+x1*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x2*y3*y3/60+x1*y2*y2/60+x3*y2*y2/60+x2*y2*y2/20+x1*x1*y2/12+x2*y2*y3/30+x1*y2*y3/60+x1*x2*y3/12+x1*x2*y2/6-x3*y2*y1/15-x2*y3*y1/15+x3*y2*y3/30-x1*x2*y1/4+x2*x1*x1/6+x3*x1*y2/12-x1*y2*y1/20;
  Ia[3+8*3] = (x2*x2*y2*y3/60-x2*x2*y2*y1/12-x3*x3*y3*y1/12+x2*x1*y2*y2/60+x2*x3*x1*x1/12+x3*x1*x1*y3/15+x2*x2*x1*y3/30+x3*x1*y1*y1/30-x1*x1*y3*y1/60+x3*x3*y3*y3/30+x3*x3*y1*y1/18+x1*x1*y3*y3/180+x1*x1*y1*y1/60+x1*x1*x1*y3/30-x1*x1*x1*y1/15+x3*x1*y3*y3/60-2.0/45.0*x3*x1*y3*y1+x3*x3*x1*y3/10-2.0/15.0*x3*x3*x1*y1-x3*x1*x1*y1/10+pow(x1,4.0)/12+x3*x3*x1*x1/12+x3*x1*x1*x1/12+x2*x3*y1*y1/18-2.0/15.0*x2*x2*x1*y1+x2*x1*y1*y1/30+x3*x3*y2*y3/60+x2*x2*y1*y1/18+x1*x1*y2*y2/180-x3*x1*y2*y1/45)+x2*x3*x1*y2/15-x3*x3*y2*y1/36-x2*x1*y3*y1/45+x2*x2*y3*y3/180+x1*x1*y2*y3/180-x2*x2*y3*y1/36+x3*x1*y2*y2/180+x3*x1*x1*y2/30+x3*x3*x1*y2/30-x2*x1*x1*y1/10+x2*x1*x1*y3/30+x2*x2*x1*y2/10-x1*x1*y2*y1/60+x2*x3*y2*y2/60+x2*x1*x1*y2/15+x2*x3*y3*y3/60+x2*x1*y3*y3/180-x2*x3*y2*y1/18+x2*x3*x1*y3/15-2.0/15.0*x2*x3*x1*y1-2.0/45.0*x2*x1*y2*y1+x2*x1*y2*y3/90+x2*x3*y2*y3/45+x3*x3*y2*y2/180+x3*x1*y2*y3/90+x1*x1*x1*y2/30-x2*x3*y3*y1/18+x2*x1*x1*x1/12+x2*x2*x1*x1/12+x2*x2*y2*y2/30;
  Ia[4+8*4] = 1.0/2.0;
  Ia[4+8*5] = x1/6+x3/6+x2/6;
  Ia[4+8*6] = x1/2+y3/6-y1/3+y2/6;
  Ia[4+8*7] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[5+8*4] = x1/6+x3/6+x2/6;
  Ia[5+8*5] = x1*x1/12+x1*x2/12+x3*x1/12+x2*x2/12+x3*x2/12+x3*x3/12;
  Ia[5+8*6] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[5+8*7] = x3*x2*y2/30+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x2*x2*y2/20+x3*x2*x1/12+x3*x2*y3/30-x3*x2*y1/15+x2*x2*y3/60-x2*x2*y1/15+x2*x2*x1/12+x3*x3*y2/60+x1*x1*y2/60+x1*x2*y3/60+x1*x2*y2/30-x1*x2*y1/20+x2*x1*x1/12+x3*x1*y2/60;
  Ia[6+8*4] = x1/2+y3/6-y1/3+y2/6;
  Ia[6+8*5] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[6+8*6] = y1*y1/4-2.0/3.0*x1*y1-y2*y1/4-y1*y3/4+x1*y2/3+y2*y2/12+y2*y3/12+x1*y3/3+x1*x1/2+y3*y3/12;
  Ia[6+8*7] = -2.0/15.0*x2*y2*y1-x1*y3*y1/20+x2*y1*y1/10+x1*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x2*y3*y3/60+x1*y2*y2/60+x3*y2*y2/60+x2*y2*y2/20+x1*x1*y2/12+x2*y2*y3/30+x1*y2*y3/60+x1*x2*y3/12+x1*x2*y2/6-x3*y2*y1/15-x2*y3*y1/15+x3*y2*y3/30-x1*x2*y1/4+x2*x1*x1/6+x3*x1*y2/12-x1*y2*y1/20;
  Ia[7+8*4] = x1*x1/6+x1*y2/24+x1*y3/24-x1*y1/12+x3*x1/6+x3*y3/12-x3*y1/8-x2*y1/8+x2*y3/24+x3*y2/24+x2*y2/12+x1*x2/6;
  Ia[7+8*5] = x3*x2*y2/30+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x2*x2*y2/20+x3*x2*x1/12+x3*x2*y3/30-x3*x2*y1/15+x2*x2*y3/60-x2*x2*y1/15+x2*x2*x1/12+x3*x3*y2/60+x1*x1*y2/60+x1*x2*y3/60+x1*x2*y2/30-x1*x2*y1/20+x2*x1*x1/12+x3*x1*y2/60;
  Ia[7+8*6] = -2.0/15.0*x2*y2*y1-x1*y3*y1/20+x2*y1*y1/10+x1*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x2*y3*y3/60+x1*y2*y2/60+x3*y2*y2/60+x2*y2*y2/20+x1*x1*y2/12+x2*y2*y3/30+x1*y2*y3/60+x1*x2*y3/12+x1*x2*y2/6-x3*y2*y1/15-x2*y3*y1/15+x3*y2*y3/30-x1*x2*y1/4+x2*x1*x1/6+x3*x1*y2/12-x1*y2*y1/20;
  Ia[7+8*7] = (x2*x2*y2*y3/60-x2*x2*y2*y1/12-x3*x3*y3*y1/12+x2*x1*y2*y2/60+x2*x3*x1*x1/12+x3*x1*x1*y3/15+x2*x2*x1*y3/30+x3*x1*y1*y1/30-x1*x1*y3*y1/60+x3*x3*y3*y3/30+x3*x3*y1*y1/18+x1*x1*y3*y3/180+x1*x1*y1*y1/60+x1*x1*x1*y3/30-x1*x1*x1*y1/15+x3*x1*y3*y3/60-2.0/45.0*x3*x1*y3*y1+x3*x3*x1*y3/10-2.0/15.0*x3*x3*x1*y1-x3*x1*x1*y1/10+pow(x1,4.0)/12+x3*x3*x1*x1/12+x3*x1*x1*x1/12+x2*x3*y1*y1/18-2.0/15.0*x2*x2*x1*y1+x2*x1*y1*y1/30+x3*x3*y2*y3/60+x2*x2*y1*y1/18+x1*x1*y2*y2/180-x3*x1*y2*y1/45)+x2*x3*x1*y2/15-x3*x3*y2*y1/36-x2*x1*y3*y1/45+x2*x2*y3*y3/180+x1*x1*y2*y3/180-x2*x2*y3*y1/36+x3*x1*y2*y2/180+x3*x1*x1*y2/30+x3*x3*x1*y2/30-x2*x1*x1*y1/10+x2*x1*x1*y3/30+x2*x2*x1*y2/10-x1*x1*y2*y1/60+x2*x3*y2*y2/60+x2*x1*x1*y2/15+x2*x3*y3*y3/60+x2*x1*y3*y3/180-x2*x3*y2*y1/18+x2*x3*x1*y3/15-2.0/15.0*x2*x3*x1*y1-2.0/45.0*x2*x1*y2*y1+x2*x1*y2*y3/90+x2*x3*y2*y3/45+x3*x3*y2*y2/180+x3*x1*y2*y3/90+x1*x1*x1*y2/30-x2*x3*y3*y1/18+x2*x1*x1*x1/12+x2*x2*x1*x1/12+x2*x2*y2*y2/30;

  Ib[0+8*0] = 1.0/2.0;
  Ib[0+8*1] = x1/6+x4/6+x3/6;
  Ib[0+8*2] = x1/2+y4/6-y1/3+y3/6;
  Ib[0+8*3] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[1+8*0] = x1/6+x4/6+x3/6;
  Ib[1+8*1] = x1*x1/12+x3*x1/12+x4*x1/12+x3*x3/12+x4*x3/12+x4*x4/12;
  Ib[1+8*2] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[1+8*3] = x4*x1*y4/30-x4*x1*y1/20+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20+x4*x1*x1/12-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x4*x4*y4/20-x4*x4*y1/15+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x4*x4*x1/12+x1*x1*y4/60+x4*x3*x1/12+x4*x3*y3/30-x4*x3*y1/15+x4*x3*y4/30+x3*x1*y4/60+x4*x1*y3/60+x3*x3*y4/60+x4*x4*y3/60;
  Ib[2+8*0] = x1/2+y4/6-y1/3+y3/6;
  Ib[2+8*1] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[2+8*2] = y1*y1/4-2.0/3.0*x1*y1-y1*y3/4-y1*y4/4+x1*y4/3+x1*x1/2+y3*y3/12+y3*y4/12+y4*y4/12+x1*y3/3;
  Ib[2+8*3] = -x1*y3*y1/20-x1*y4*y1/20-2.0/15.0*x4*y4*y1+x4*x1*y4/6-x4*x1*y1/4+x4*y4*y4/20+x4*y1*y1/10+x1*x1*x1/6+x4*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x1*y4*y4/60+x1*x1*y4/12+x3*x1*y4/12-x3*y4*y1/15+x1*y3*y4/60+x3*y3*y4/30+x4*x1*y3/12+x4*y3*y4/30-x4*y3*y1/15+x3*y4*y4/60+x4*y3*y3/60;
  Ib[3+8*0] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[3+8*1] = x4*x1*y4/30-x4*x1*y1/20+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20+x4*x1*x1/12-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x4*x4*y4/20-x4*x4*y1/15+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x4*x4*x1/12+x1*x1*y4/60+x4*x3*x1/12+x4*x3*y3/30-x4*x3*y1/15+x4*x3*y4/30+x3*x1*y4/60+x4*x1*y3/60+x3*x3*y4/60+x4*x4*y3/60;
  Ib[3+8*2] = -x1*y3*y1/20-x1*y4*y1/20-2.0/15.0*x4*y4*y1+x4*x1*y4/6-x4*x1*y1/4+x4*y4*y4/20+x4*y1*y1/10+x1*x1*x1/6+x4*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x1*y4*y4/60+x1*x1*y4/12+x3*x1*y4/12-x3*y4*y1/15+x1*y3*y4/60+x3*y3*y4/30+x4*x1*y3/12+x4*y3*y4/30-x4*y3*y1/15+x3*y4*y4/60+x4*y3*y3/60;
  Ib[3+8*3] = x4*x1*x1*x1/12+x4*x4*y4*y4/30+x4*x4*y1*y1/18+x1*x1*y4*y4/180+(-x3*x3*y3*y1/12-x3*x3*y4*y1/36+x3*x1*x1*y3/15+x3*x1*y1*y1/30-x1*x1*y3*y1/60+x4*x4*x1*y4/10-2.0/15.0*x4*x4*x1*y1+x3*x3*y3*y3/30+x3*x3*y1*y1/18+x1*x1*y3*y3/180+x1*x1*y1*y1/60+x1*x1*x1*y3/30-x1*x1*x1*y1/15+x3*x1*y3*y3/60-2.0/45.0*x3*x1*y3*y1+x3*x3*x1*y3/10-2.0/15.0*x3*x3*x1*y1-x3*x1*x1*y1/10+pow(x1,4.0)/12+x3*x3*x1*x1/12+x3*x1*x1*x1/12+x4*x1*x1*y4/15-x4*x1*x1*y1/10+x1*x1*x1*y4/30-x4*x4*y4*y1/12-2.0/45.0*x4*x1*y4*y1+x4*x1*y4*y4/60+x4*x1*y1*y1/30-x1*x1*y4*y1/60+x4*x4*x1*x1/12)-x3*x4*y3*y1/18+x3*x1*y3*y4/90-2.0/15.0*x3*x4*x1*y1-x3*x4*y4*y1/18+x3*x1*y4*y4/180+x1*x1*y3*y4/180+x4*x4*y3*y4/60-x4*x4*y3*y1/36+x3*x4*x1*x1/12+x4*x1*y3*y3/180+x3*x4*y3*y3/60+x3*x3*x1*y4/30+x3*x1*x1*y4/30-x4*x1*y3*y1/45+x4*x1*y3*y4/90+x3*x3*y4*y4/180+x3*x4*y3*y4/45-x3*x1*y4*y1/45+x3*x4*x1*y3/15+x4*x4*y3*y3/180+x4*x1*x1*y3/30+x3*x3*y3*y4/60+x4*x4*x1*y3/30+x3*x4*y4*y4/60+x3*x4*y1*y1/18+x3*x4*x1*y4/15;
  Ib[4+8*4] = 1.0/2.0;
  Ib[4+8*5] = x1/6+x4/6+x3/6;
  Ib[4+8*6] = x1/2+y4/6-y1/3+y3/6;
  Ib[4+8*7] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[5+8*4] = x1/6+x4/6+x3/6;
  Ib[5+8*5] = x1*x1/12+x3*x1/12+x4*x1/12+x3*x3/12+x4*x3/12+x4*x4/12;
  Ib[5+8*6] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[5+8*7] = x4*x1*y4/30-x4*x1*y1/20+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20+x4*x1*x1/12-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x4*x4*y4/20-x4*x4*y1/15+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x4*x4*x1/12+x1*x1*y4/60+x4*x3*x1/12+x4*x3*y3/30-x4*x3*y1/15+x4*x3*y4/30+x3*x1*y4/60+x4*x1*y3/60+x3*x3*y4/60+x4*x4*y3/60;
  Ib[6+8*4] = x1/2+y4/6-y1/3+y3/6;
  Ib[6+8*5] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[6+8*6] = y1*y1/4-2.0/3.0*x1*y1-y1*y3/4-y1*y4/4+x1*y4/3+x1*x1/2+y3*y3/12+y3*y4/12+y4*y4/12+x1*y3/3;
  Ib[6+8*7] = -x1*y3*y1/20-x1*y4*y1/20-2.0/15.0*x4*y4*y1+x4*x1*y4/6-x4*x1*y1/4+x4*y4*y4/20+x4*y1*y1/10+x1*x1*x1/6+x4*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x1*y4*y4/60+x1*x1*y4/12+x3*x1*y4/12-x3*y4*y1/15+x1*y3*y4/60+x3*y3*y4/30+x4*x1*y3/12+x4*y3*y4/30-x4*y3*y1/15+x3*y4*y4/60+x4*y3*y3/60;
  Ib[7+8*4] = x1*x1/6+x1*y3/24-x1*y1/12+x1*y4/24+x3*x1/6+x4*x1/6+x3*y3/12-x3*y1/8+x4*y4/12-x4*y1/8+x3*y4/24+x4*y3/24;
  Ib[7+8*5] = x4*x1*y4/30-x4*x1*y1/20+x1*x1*x1/12-x3*x3*y1/15+x3*x3*y3/20+x4*x1*x1/12-x1*x1*y1/30+x1*x1*y3/60+x3*x3*x1/12+x4*x4*y4/20-x4*x4*y1/15+x3*x1*x1/12+x3*x1*y3/30-x3*x1*y1/20+x4*x4*x1/12+x1*x1*y4/60+x4*x3*x1/12+x4*x3*y3/30-x4*x3*y1/15+x4*x3*y4/30+x3*x1*y4/60+x4*x1*y3/60+x3*x3*y4/60+x4*x4*y3/60;
  Ib[7+8*6] = -x1*y3*y1/20-x1*y4*y1/20-2.0/15.0*x4*y4*y1+x4*x1*y4/6-x4*x1*y1/4+x4*y4*y4/20+x4*y1*y1/10+x1*x1*x1/6+x4*x1*x1/6-x1*x1*y1/6+x1*y1*y1/20+x1*x1*y3/12+x1*y3*y3/60+x3*y3*y3/20+x3*y1*y1/10+x3*x1*x1/6-2.0/15.0*x3*y3*y1+x3*x1*y3/6-x3*x1*y1/4+x1*y4*y4/60+x1*x1*y4/12+x3*x1*y4/12-x3*y4*y1/15+x1*y3*y4/60+x3*y3*y4/30+x4*x1*y3/12+x4*y3*y4/30-x4*y3*y1/15+x3*y4*y4/60+x4*y3*y3/60;
  Ib[7+8*7] = x4*x1*x1*x1/12+x4*x4*y4*y4/30+x4*x4*y1*y1/18+x1*x1*y4*y4/180+(-x3*x3*y3*y1/12-x3*x3*y4*y1/36+x3*x1*x1*y3/15+x3*x1*y1*y1/30-x1*x1*y3*y1/60+x4*x4*x1*y4/10-2.0/15.0*x4*x4*x1*y1+x3*x3*y3*y3/30+x3*x3*y1*y1/18+x1*x1*y3*y3/180+x1*x1*y1*y1/60+x1*x1*x1*y3/30-x1*x1*x1*y1/15+x3*x1*y3*y3/60-2.0/45.0*x3*x1*y3*y1+x3*x3*x1*y3/10-2.0/15.0*x3*x3*x1*y1-x3*x1*x1*y1/10+pow(x1,4.0)/12+x3*x3*x1*x1/12+x3*x1*x1*x1/12+x4*x1*x1*y4/15-x4*x1*x1*y1/10+x1*x1*x1*y4/30-x4*x4*y4*y1/12-2.0/45.0*x4*x1*y4*y1+x4*x1*y4*y4/60+x4*x1*y1*y1/30-x1*x1*y4*y1/60+x4*x4*x1*x1/12)-x3*x4*y3*y1/18+x3*x1*y3*y4/90-2.0/15.0*x3*x4*x1*y1-x3*x4*y4*y1/18+x3*x1*y4*y4/180+x1*x1*y3*y4/180+x4*x4*y3*y4/60-x4*x4*y3*y1/36+x3*x4*x1*x1/12+x4*x1*y3*y3/180+x3*x4*y3*y3/60+x3*x3*x1*y4/30+x3*x1*x1*y4/30-x4*x1*y3*y1/45+x4*x1*y3*y4/90+x3*x3*y4*y4/180+x3*x4*y3*y4/45-x3*x1*y4*y1/45+x3*x4*x1*y3/15+x4*x4*y3*y3/180+x4*x1*x1*y3/30+x3*x3*y3*y4/60+x4*x4*x1*y3/30+x3*x4*y4*y4/60+x3*x4*y1*y1/18+x3*x4*x1*y4/15;

  matrix_transpose(8,
		Si,
		Sit);
  
  matrix_multiplication(8,
		Ia,ma,
		Ia);
  
  matrix_multiplication(8,
		Ib,mb,
		Ib);

  matrices_addition(8,
		Ia,Ib,
		I);
  
  matrices_product(8,
		Sit,I,
		aux);
  
  matrices_product(8,
		aux,Si,
		Me);

  memset(Ia,0,64*sizeof(double));
  memset(Ib,0,64*sizeof(double));

  Ia[1+8*1] = 1.0/2.0;
  Ia[1+8*3] = x1/2+y3/6-y1/3+y2/6;
  Ia[1+8*6] = nu/2;
  Ia[1+8*7] = nu*x1/6+nu*x3/6+nu*x2/6;
  Ia[2+8*2] = 1.0/4.0-nu/4;
  Ia[2+8*3] = x1/12-nu*x2/12-nu*x1/12+x3/12+x2/12-nu*x3/12;
  Ia[2+8*5] = 1.0/4.0-nu/4;
  Ia[2+8*7] = x1/4-y2*nu/12-nu*x1/4-y1/6+y1*nu/6+y3/12+y2/12-y3*nu/12;
  Ia[3+8*1] = x1/2+y3/6-y1/3+y2/6;
  Ia[3+8*2] = x1/12-nu*x2/12-nu*x1/12+x3/12+x2/12-nu*x3/12;
  Ia[3+8*3] = -nu*x1*x1/24+13.0/24.0*x1*x1+x1*y3/3-2.0/3.0*x1*y1+x1*y2/3+x3*x1/24-x1*nu*x2/24+x1*x2/24-nu*x3*x1/24+y2*y2/12+x2*x2/24-y1*y3/4-x2*x2*nu/24+y1*y1/4-y2*y1/4+y3*y3/12+y2*y3/12+x3*x3/24-x3*nu*x2/24+x3*x2/24-nu*x3*x3/24;
  Ia[3+8*5] = x1/12-nu*x2/12-nu*x1/12+x3/12+x2/12-nu*x3/12;
  Ia[3+8*6] = nu*x1/2+y3*nu/6-y1*nu/3+y2*nu/6;
  Ia[3+8*7] = nu*x1*y2/48+nu*x1*x1/12+nu*x2*y2/24+x1*x1/12+x1*y2/48+x1*y3/48-x1*y1/24+x1*nu*x2/12+x3*x1/12+x3*y3/24-x3*y1/16+nu*x3*x1/12+nu*x3*y3/24-nu*x3*y1/16-nu*x1*y1/24+nu*x1*y3/48-x2*y1/16+x2*y3/48+x3*y2/48+x2*y2/24+x1*x2/12+nu*x2*y3/48+nu*x3*y2/48-nu*x2*y1/16;
  Ia[5+8*2] = 1.0/4.0-nu/4;
  Ia[5+8*3] = x1/12-nu*x2/12-nu*x1/12+x3/12+x2/12-nu*x3/12;
  Ia[5+8*5] = 1.0/4.0-nu/4;
  Ia[5+8*7] = x1/4-y2*nu/12-nu*x1/4-y1/6+y1*nu/6+y3/12+y2/12-y3*nu/12;
  Ia[6+8*1] = nu/2;
  Ia[6+8*3] = nu*x1/2+y3*nu/6-y1*nu/3+y2*nu/6;
  Ia[6+8*6] = 1.0/2.0;
  Ia[6+8*7] = x1/6+x3/6+x2/6;
  Ia[7+8*1] = nu*x1/6+nu*x3/6+nu*x2/6;
  Ia[7+8*2] = x1/4-y2*nu/12-nu*x1/4-y1/6+y1*nu/6+y3/12+y2/12-y3*nu/12;
  Ia[7+8*3] = nu*x1*y2/48+nu*x1*x1/12+nu*x2*y2/24+x1*x1/12+x1*y2/48+x1*y3/48-x1*y1/24+x1*nu*x2/12+x3*x1/12+x3*y3/24-x3*y1/16+nu*x3*x1/12+nu*x3*y3/24-nu*x3*y1/16-nu*x1*y1/24+nu*x1*y3/48-x2*y1/16+x2*y3/48+x3*y2/48+x2*y2/24+x1*x2/12+nu*x2*y3/48+nu*x3*y2/48-nu*x2*y1/16;
  Ia[7+8*5] = x1/4-y2*nu/12-nu*x1/4-y1/6+y1*nu/6+y3/12+y2/12-y3*nu/12;
  Ia[7+8*6] = x1/6+x3/6+x2/6;
  Ia[7+8*7] = -nu*x1*y2/6-nu*x1*x1/4+x3*x3/12+x1*x1/3+y3*y3/24+y1*y1/8+y2*y3/24+x1*y2/6+x1*y3/6-x1*y1/3-y2*y1/8-y1*y3/8+x2*x2/12+y2*y2/24+x3*x1/12-y2*y3*nu/24-y3*y3*nu/24-y2*y2*nu/24+y3*y1*nu/8+y2*y1*nu/8-y1*y1*nu/8+nu*x1*y1/3-nu*x1*y3/6+x3*x2/12+x1*x2/12;

  Ib[1+8*1] = 1.0/2.0;
  Ib[1+8*3] = x1/2+y4/6-y1/3+y3/6;
  Ib[1+8*6] = nu/2;
  Ib[1+8*7] = nu*x1/6+nu*x4/6+nu*x3/6;
  Ib[2+8*2] = 1.0/4.0-nu/4;
  Ib[2+8*3] = x1/12-nu*x3/12-nu*x1/12+x4/12+x3/12-nu*x4/12;
  Ib[2+8*5] = 1.0/4.0-nu/4;
  Ib[2+8*7] = x1/4-y3*nu/12-nu*x1/4-y1/6+y1*nu/6+y4/12+y3/12-y4*nu/12;
  Ib[3+8*1] = x1/2+y4/6-y1/3+y3/6;
  Ib[3+8*2] = x1/12-nu*x3/12-nu*x1/12+x4/12+x3/12-nu*x4/12;
  Ib[3+8*3] = -nu*x1*x1/24+13.0/24.0*x1*x1-2.0/3.0*x1*y1+x1*y3/3+x4*x1/24+x3*x1/24-nu*x3*x1/24+x1*y4/3-nu*x4*x1/24-y1*y4/4-y1*y3/4+y1*y1/4+y4*y4/12+y3*y3/12+x3*x3/24-nu*x3*x3/24+y3*y4/12+x4*x4/24+x4*x3/24-x4*nu*x3/24-nu*x4*x4/24;
  Ib[3+8*5] = x1/12-nu*x3/12-nu*x1/12+x4/12+x3/12-nu*x4/12;
  Ib[3+8*6] = nu*x1/2+y4*nu/6-y1*nu/3+y3*nu/6;
  Ib[3+8*7] = x4*y4*nu/24-x4*y1*nu/16+nu*x1*x1/12+x1*x1/12+nu*x1*y4/48+nu*x4*x1/12+x1*y3/48-x1*y1/24+x1*y4/48+x3*x1/12+x4*x1/12+x3*y3/24-x3*y1/16+nu*x3*x1/12+nu*x3*y3/24-nu*x3*y1/16-nu*x1*y1/24+nu*x1*y3/48+x4*y4/24-x4*y1/16+x3*y4/48+x4*y3/48+nu*x4*y3/48+nu*x3*y4/48;
  Ib[5+8*2] = 1.0/4.0-nu/4;
  Ib[5+8*3] = x1/12-nu*x3/12-nu*x1/12+x4/12+x3/12-nu*x4/12;
  Ib[5+8*5] = 1.0/4.0-nu/4;
  Ib[5+8*7] = x1/4-y3*nu/12-nu*x1/4-y1/6+y1*nu/6+y4/12+y3/12-y4*nu/12;
  Ib[6+8*1] = nu/2;
  Ib[6+8*3] = nu*x1/2+y4*nu/6-y1*nu/3+y3*nu/6;
  Ib[6+8*6] = 1.0/2.0;
  Ib[6+8*7] = x1/6+x4/6+x3/6;
  Ib[7+8*1] = nu*x1/6+nu*x4/6+nu*x3/6;
  Ib[7+8*2] = x1/4-y3*nu/12-nu*x1/4-y1/6+y1*nu/6+y4/12+y3/12-y4*nu/12;
  Ib[7+8*3] = x4*y4*nu/24-x4*y1*nu/16+nu*x1*x1/12+x1*x1/12+nu*x1*y4/48+nu*x4*x1/12+x1*y3/48-x1*y1/24+x1*y4/48+x3*x1/12+x4*x1/12+x3*y3/24-x3*y1/16+nu*x3*x1/12+nu*x3*y3/24-nu*x3*y1/16-nu*x1*y1/24+nu*x1*y3/48+x4*y4/24-x4*y1/16+x3*y4/48+x4*y3/48+nu*x4*y3/48+nu*x3*y4/48;
  Ib[7+8*5] = x1/4-y3*nu/12-nu*x1/4-y1/6+y1*nu/6+y4/12+y3/12-y4*nu/12;
  Ib[7+8*6] = x1/6+x4/6+x3/6;
  Ib[7+8*7] = -y4*y4*nu/24-y3*y4*nu/24-nu*x1*x1/4+x3*x3/12+x1*x1/3+y3*y3/24+y1*y1/8-nu*x1*y4/6+x1*y3/6-x1*y1/3-y1*y3/8+y4*y1*nu/8-y1*y4/8+y3*y4/24+x1*y4/6+x3*x1/12-y3*y3*nu/24+y3*y1*nu/8-y1*y1*nu/8+x4*x4/12+x4*x1/12+y4*y4/24+nu*x1*y1/3-nu*x1*y3/6+x4*x3/12;
  
  matrix_multiplication(8,
		Ia,ka,
		Ia);
  
  matrix_multiplication(8,
		Ib,kb,
		Ib);
  
  matrices_addition(8,
		Ia,Ib,
		I);
  
  matrices_product(8,
		Sit,I,
		aux);
  
  matrices_product(8,
		aux,Si,
		Ke);

  freeMem(Si);
  freeMem(Sit);
  freeMem(Ia);
  freeMem(Ib);
  freeMem(I);
  freeMem(aux);
}
