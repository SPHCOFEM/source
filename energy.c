#include "header.h"

#undef __FUNC__
#define __FUNC__ "energy"

void energy(int n_s,
	
  double *x_s,double *v_s,
	double *m_s,double *u_s,

	int n_b,
	double *x_b,double *x_b_0,double *v_b,

	int n_n,
	double *M,double *B,double *K,double *f,
	    
	double *kine_s,double *inne_s,double *pote_s,
	double *kine_b,double *defo_b,double *disi_b,double *pote_b,
	double *tote,
	    
	double ax,double ay,double az)
{
  int i,j;
  
  double *kine_v_b=0,*disi_v_b=0,*defo_v_b=0;

  if (n_n>0)
  {
    kine_v_b=createMemMore(double,n_n);
    disi_v_b=createMemMore(double,n_n);
    defo_v_b=createMemMore(double,n_n);
  }
  
  *kine_s=0.0;
  *inne_s=0.0;
  *pote_s=0.0;

  *kine_b=0.0;
  *defo_b=0.0;
  *disi_b=0.0;
  *pote_b=0.0;

  for (i=0;i<n_s;i++)
  {
    *kine_s=*kine_s+0.5*m_s[i]*(sqr(v_s[i])+sqr(v_s[i+n_s])+sqr(v_s[i+2*n_s]));
    *inne_s=*inne_s+m_s[i]*u_s[i];
    *pote_s=*pote_s+m_s[i]*(ax*x_s[i]+ay*x_s[i+n_s]+az*x_s[i+2*n_s]);
  }

  for (i=0;i<n_n;i++)
  {
    kine_v_b[i]=0.0;
    defo_v_b[i]=0.0;
    disi_v_b[i]=0.0;

    for (j=0;j<n_n;j++)
  	{
	    kine_v_b[i]=kine_v_b[i]+0.5*v_b[j]*M[i+j*n_n];
	    defo_v_b[i]=defo_v_b[i]+0.5*(x_b[j]-x_b_0[j])*K[i+j*n_n];
  	  disi_v_b[i]=disi_v_b[i]+0.5*v_b[j]*B[i+j*n_n];
	  }
    *kine_b=*kine_b+v_b[i]*kine_v_b[i];
    *defo_b=*defo_b+(x_b[i]-x_b_0[i])*defo_v_b[i];
    *disi_b=*disi_b+v_b[i]*disi_v_b[i];
    *pote_b=*pote_b+x_b[i]*f[i];
  }

  if (n_n>0)
  {
    freeMem(kine_v_b);
    freeMem(disi_v_b);
    freeMem(defo_v_b);
  }

  *tote=*kine_s+*inne_s+*pote_s+*kine_b+*defo_b+*disi_b+*pote_b;
}
