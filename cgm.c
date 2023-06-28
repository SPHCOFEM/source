#include "header.h"

#undef __FUNC__
#define __FUNC__ "cgm"

void cgm(int n_n,

  int n_b,int *num_b,
  double *x_b_0,
	double *x_b,double *v_b,double *a_b,

	double *M,double *B,double *K,double *f,

	int n_o_0_b,int *ind_o_0_b,int *o_0_b,

	int n_r_0_b,int *ind_r_0_b,

	int dim)
{
  int i,j,k;

  double *r=0,*v=0,*Ar=0,*Av=0,*g=0;
  double nom,den,al,be;
  double nrm;

  /* allocation of auxiliary variables */
  if (n_n>0)
  {
    r=createMemMore(double,n_n);
    v=createMemMore(double,n_n);
    Ar=createMemMore(double,n_n);
    Av=createMemMore(double,n_n);
    g=createMemMore(double,n_n);
  }
  
  /* right-hand side construction */
  for (i=0;i<n_n;i++)
  {
    g[i]=f[i];
    for (j=0;j<n_n;j++)
	  {
	    g[i]=g[i]-B[i+j*n_n]*v_b[j]-K[i+j*n_n]*(x_b[j]-x_b_0[j]);
	  }
  }

  /* known values - boundary conditions */
  for (i=0;i<n_o_0_b;i++)
  {
    for (k=0;k<n_n;k++)
	  {
	    if (o_0_b[i])
	    {
	      M[ind_o_0_b[i]+k*n_n]=0.0;
	      M[k+ind_o_0_b[i]*n_n]=0.0;
	    }
	  
	    if ((dim>1)&&(o_0_b[i+n_o_0_b]))
	    {
	      M[ind_o_0_b[i]+n_b+k*n_n]=0.0;
	      M[k+(ind_o_0_b[i]+n_b)*n_n]=0.0;
	    }

	    if ((dim>2)&&(o_0_b[i+2*n_o_0_b]))
	    {
	      M[ind_o_0_b[i]+2*n_b+k*n_n]=0.0;
	      M[k+(ind_o_0_b[i]+2*n_b)*n_n]=0.0;
	    }
	  }

    if (o_0_b[i]) g[ind_o_0_b[i]]=0.0;
    if ((dim>1)&&(o_0_b[i+n_o_0_b])) g[ind_o_0_b[i]+n_b]=0.0;
    if ((dim>2)&&(o_0_b[i+2*n_o_0_b])) g[ind_o_0_b[i]+2*n_b]=0.0;

    if (o_0_b[i]) M[ind_o_0_b[i]+ind_o_0_b[i]*n_n]=1.0;
    if ((dim>1)&&(o_0_b[i+n_o_0_b])) M[ind_o_0_b[i]+n_b+(ind_o_0_b[i]+n_b)*n_n]=1.0;
    if ((dim>2)&&(o_0_b[i+2*n_o_0_b])) M[ind_o_0_b[i]+2*n_b+(ind_o_0_b[i]+2*n_b)*n_n]=1.0;
  }

  /* known values - rigid bodies */
  for (i=0;i<n_r_0_b;i++)
  {
    for (k=0;k<n_n;k++)
	  {
      M[ind_r_0_b[i]+k*n_n]=0.0;
      M[k+ind_r_0_b[i]*n_n]=0.0;
	  
	    if (dim>1)
	    {
	      M[ind_r_0_b[i]+n_b+k*n_n]=0.0;
	      M[k+(ind_r_0_b[i]+n_b)*n_n]=0.0;
	    }

	    if (dim>2)
	    {
	      M[ind_r_0_b[i]+2*n_b+k*n_n]=0.0;
	      M[k+(ind_r_0_b[i]+2*n_b)*n_n]=0.0;
	    }
	  }

    g[ind_r_0_b[i]]=0.0;
    if (dim>1) g[ind_r_0_b[i]+n_b]=0.0;
    if (dim>2) g[ind_r_0_b[i]+2*n_b]=0.0;

    M[ind_r_0_b[i]+ind_r_0_b[i]*n_n]=1.0;
    if (dim>1) M[ind_r_0_b[i]+n_b+(ind_r_0_b[i]+n_b)*n_n]=1.0;
    if (dim>2) M[ind_r_0_b[i]+2*n_b+(ind_r_0_b[i]+2*n_b)*n_n]=1.0;
  }

  /* printing matrices for checking
  for (i=0;i<n_n;i++)
    {
      for (j=0;j<n_n;j++)
	{
	  printf("%f ",M[i+n_n*j]);
	}
      printf("\n");
    }
  printf("\n");
  for (i=0;i<n_n;i++)
    {
      for (j=0;j<n_n;j++)
	{
	  printf("%f ",B[i+n_n*j]);
	}
      printf("\n");
    }
  printf("\n");
  for (i=0;i<n_n;i++)
    {
      for (j=0;j<n_n;j++)
	{
	  printf("%f ",K[i+n_n*j]);
	}
      printf("\n");
    }
  printf("\n");
  for (i=0;i<n_n;i++)
    {
      printf("%f\n",f[i]);
    }
  exit(0);
  */

  /* cgm */
  nrm=0.0;
  for (i=0;i<n_n;i++)
  {
    r[i]=-g[i];
    for (j=0;j<n_n;j++)
	  {
	    r[i]=r[i]+M[i+j*n_n]*a_b[j];
	  }
    v[i]=-r[i];
    nrm=nrm+sqr(r[i]);
  }
  nrm=sqrt(nrm);

  while (nrm>=EPS)
  {
    for (i=0;i<n_n;i++)
	  {
	    Av[i]=0.0;
	    for (j=0;j<n_n;j++)
	    {
	      Av[i]=Av[i]+M[i+j*n_n]*v[j];
	    }
	  }

    nom=0.0;
    den=0.0;
    for (i=0;i<n_n;i++)
	  {
	    nom=nom+r[i]*v[i];
	    den=den+Av[i]*v[i];
	  }
    al=-nom/den;

    nrm=0.0;
    for (i=0;i<n_n;i++)
	  {
	    a_b[i]=a_b[i]+al*v[i];
	    nrm=nrm+sqr(r[i]+al*Av[i]);
	  }
    nrm=sqrt(nrm);

    if (nrm<EPS)
	  {
	    break;
	  }

    for (i=0;i<n_n;i++)
	  {
	    r[i]=r[i]+al*Av[i];
	  }
      
    for (i=0;i<n_n;i++)
	  {
      Ar[i]=0.0;
      Av[i]=0.0;
      for (j=0;j<n_n;j++)
	    {
	      Ar[i]=Ar[i]+M[i+j*n_n]*r[j];
	      Av[i]=Av[i]+M[i+j*n_n]*v[j];
	    }
  	}

    nom=0.0;
    den=0.0;
    for (i=0;i<n_n;i++)
	  {
	    nom=nom+Ar[i]*v[i];
	    den=den+Av[i]*v[i];
	  }
    be=nom/den;

    for (i=0;i<n_n;i++)
	  {
	    v[i]=-r[i]+be*v[i];
	  }
  }
  
  /* free - fem - auxiliary variables */
  if (n_n>0)
  {
    freeMem(r);
    freeMem(v);
    freeMem(Ar);
    freeMem(Av);
    freeMem(g);
  }
}
