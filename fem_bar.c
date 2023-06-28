#include "header.h"

#undef __FUNC__
#define __FUNC__ "fem_bar"

void fem_bar(double *Me,double *Ke,

  double x1,double x2,double y1,double y2,

  double A,double rho,double E)
{
  double c=0.0,s=0.0,t1,t2,l;

  t1=x2-x1;
  t2=y2-y1;

  l=norm(t1,t2,0.0);
  
  if (t1==0)
  {
    c=0.0;
    if (t2>0) s=1.0;
    else s=-1.0;
  }
  if (t2==0)
  {
    s=0.0;
    if (t1>0.0) c=1.0;
    else c=-1.0;
  }
  if ((t1!=0.0)&&(t2!=0.0))
  {
    s=fabs(t2/l);
    c=fabs(t1/l);
    if (t1<0.0)
	  {
	    c=-c;
	    if (t2<0.0) s=-s;
	  }
    if (t1>0.0)
	  {
	    c=-c;
	    if (t2<0.0) s=-s;
	  }
  }

  Me[0+4*0] = rho*A*l*(35.0*c*c+39.0*s*s)/105;
  Me[0+4*1] = rho*A*l*(35.0*c*c+27.0*s*s)/210;
  Me[0+4*2] = -4.0/105.0*c*rho*A*l*s;
  Me[0+4*3] = 4.0/105.0*c*rho*A*l*s;
  Me[1+4*0] = rho*A*l*(35.0*c*c+27.0*s*s)/210;
  Me[1+4*1] = rho*A*l*(35.0*c*c+39.0*s*s)/105;
  Me[1+4*2] = 4.0/105.0*c*rho*A*l*s;
  Me[1+4*3] = -4.0/105.0*c*rho*A*l*s;
  Me[2+4*0] = -4.0/105.0*c*rho*A*l*s;
  Me[2+4*1] = 4.0/105.0*c*rho*A*l*s;
  Me[2+4*2] = rho*A*l*(35.0*s*s+39.0*c*c)/105;
  Me[2+4*3] = rho*A*l*(35.0*s*s+27.0*c*c)/210;
  Me[3+4*0] = 4.0/105.0*c*rho*A*l*s;
  Me[3+4*1] = -4.0/105.0*c*rho*A*l*s;
  Me[3+4*2] = rho*A*l*(35.0*s*s+27.0*c*c)/210;
  Me[3+4*3] = rho*A*l*(35.0*s*s+39.0*c*c)/105;

  Ke[0+4*0] = c*c*E*A/l;
  Ke[0+4*1] = -c*c*E*A/l;
  Ke[0+4*2] = c*E*A/l*s;
  Ke[0+4*3] = -c*E*A/l*s;
  Ke[1+4*0] = -c*c*E*A/l;
  Ke[1+4*1] = c*c*E*A/l;
  Ke[1+4*2] = -c*E*A/l*s;
  Ke[1+4*3] = c*E*A/l*s;
  Ke[2+4*0] = c*E*A/l*s;
  Ke[2+4*1] = -c*E*A/l*s;
  Ke[2+4*2] = s*s*E*A/l;
  Ke[2+4*3] = -s*s*E*A/l;
  Ke[3+4*0] = -c*E*A/l*s;
  Ke[3+4*1] = c*E*A/l*s;
  Ke[3+4*2] = -s*s*E*A/l;
  Ke[3+4*3] = s*s*E*A/l;
}
