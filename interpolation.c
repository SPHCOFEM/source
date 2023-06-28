#include "header.h"

#undef __FUNC__
#define __FUNC__ "interpolation"

double interpolation(double xi,int fun,int n_u,
	int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u)

{
	int i,fi=0; // if fi is not set, warning that may be uninitialized
	double yi=0.0;

	/* find function */
	for (i=0;i<n_u;i++)
	{
		if (num_u[i]==fun)
		{
			fi=i;
			break;
		}
	}

	if (xi<=fxi_u[fun_loc[fi]]) // xi below interval range -> f(xi) = f(x1)
	{
		yi=fyi_u[fun_loc[fi]];
	}
	else if (xi>fxi_u[fun_loc[fi]+fun_u[fi]-1]) // xi over interval range -> f(xi) = f(x2)
	{
		yi=fyi_u[fun_loc[fi]+fun_u[fi]-1];
	}
	else
	{
		for (i=1;i<fun_u[fi];i++)
		{
			if ((xi>fxi_u[i+fun_loc[fi]-1])&&(xi<=fxi_u[i+fun_loc[fi]])) // xi in range [x1, x2]
			{
				yi=fyi_u[i+fun_loc[fi]-1]+
					(xi-fxi_u[i+fun_loc[fi]-1])/
						(fxi_u[i+fun_loc[fi]]-fxi_u[i+fun_loc[fi]-1])*
							(fyi_u[i+fun_loc[fi]]-fyi_u[i+fun_loc[fi]-1]); // linear interpolation
				break;
			}
		}
	}

  	return(yi);
}
