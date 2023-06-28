#include "header.h"

#undef __FUNC__
#define __FUNC__ "kernel"

double kernel(double s,double h,double sigma,int dim)
{
	double w=0.0;
  	int kernel=1;

  	switch (kernel)
  	{
		case 1:
		{
			if ((s>=0.0)&&(s<1.0)) w=1.0-1.5*sqr(s)+0.75*pow(s,(double)3);
			if ((s>=1.0)&&(s<2.0)) w=0.25*pow(2.0-s,3.0);
			w=w*sigma/pow(h,(double)dim);
			break;
		}
		case 2: /* cubic B-spline with a fixed number of neighbours (~ 13) */
		{
			if ((s>=0.0)&&(s<1.0)) w=2.0/3.0-sqr(s)+1.0/2.0*pow(s,(double)3);
			if ((s>=1.0)&&(s<2.0)) w=1.0/6.0*pow(2.0-s,3.0);
			w=15.0*w/7.0/PI/sqr(h);
			break;
		}
		case 3: /* Wendland C4 function with varying number of neighbours */
		{
			if ((s>=0.0)&&(s<1.0)) w=pow(1.0-s,6.0)*(1.0+6.0*s+35.0/3.0*sqr(s));
			w=9.0*w/PI/sqr(h);
			break;
		}
	}
  	return(w);
}

#undef __FUNC__
#define __FUNC__ "gradient_of_kernel"

double gradient_of_kernel(double s,double h,double sigma,int dim)
{
	double dw=0.0;
  	int kernel=1;

  	switch (kernel)
	{
		case 1:
	  	{
			if ((s>=0.0)&&(s<1.0)) dw=-3.0*s+2.25*sqr(s);
			if ((s>=1.0)&&(s<2.0)) dw=-0.75*sqr(2.0-s);
				dw=dw*sigma/pow(h,(double)dim);
			break;
		}
		case 2: /* cubic B-spline with a fixed number of neighbours (~ 13) */
		{
			if ((s>=0.0)&&(s<1.0)) dw=-2.0*s+3.0/2.0*sqr(s);
			if ((s>=1.0)&&(s<2.0)) dw=-1.0/2.0*sqr(2.0-s);
			dw=15.0*dw/7.0/PI/sqr(h);
			break;
		}
		case 3: /* Wendland C4 function with varying number of neighbours */
	  	{
			if ((s>=0.0)&&(s<1.0)) dw=-6.0*pow(1.0-s,5.0)*(1.0+6.0*s+35.0/3.0*sqr(s))+pow(1.0-s,6.0)*(6.0+70.0/3.0*s);
			dw=9.0*dw/PI/sqr(h);
			break;
	  	}
	}
  	return(dw);
}

#undef __FUNC__
#define __FUNC__ "gradient_of_kernel"

double second_derivative_of_kernel(double s,double h,double sigma,int dim)
{
	double d2w=0.0;
  	int kernel=1;

  	switch (kernel)
	{
		case 1:
	  	{
			if ((s>=0.0)&&(s<1.0)) d2w=-3.0+4.5*s;
			if ((s>=1.0)&&(s<2.0)) d2w=1.5*(2.0-s);
				d2w=d2w*sigma/pow(h,(double)dim);
			break;
	  	}
		case 2: /* cubic B-spline with a fixed number of neighbours (~ 13) */
	  	{
			if ((s>=0.0)&&(s<1.0)) d2w=-2.0+3.0*s;
			if ((s>=1.0)&&(s<2.0)) d2w=2.0-s;
			d2w=15.0*d2w/7.0/PI/sqr(h);
			break;
	  	}
		case 3: /* Wendland C4 function with varying number of neighbours */
	  	{
			if ((s>=0.0)&&(s<1.0)) d2w=6.0/5.0*pow(1.0-s,4.0)*(1.0+6.0*s+35.0/3.0*sqr(s))-6.0*pow(1.0-s,5.0)*(6.0+70.0/3.0*s)-6.0*pow(1.0-s,5.0)*(6.0+70.0/3.0*s)+pow(1.0-s,6.0)*70.0/3.0;
			d2w=9.0*d2w/PI/sqr(h);
			break;
	  	}
	}
  	return(d2w);
}
