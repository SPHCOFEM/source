#include "header.h"
#include "contact_forces.h"

#undef __FUNC__
#define __FUNC__ "contacts"

void contacts(int *type_m,int *domain_m,
	
	int n_c,double *f_c,
	int *seg_c,int *slave_c,int *master_c,int *sw_c,int *mat_c,
	double *m_slave,double *m_master,double *ct_c,double *klin_c,double *knon_c,double *kf_c,double *kd_c,

	int n_s,double *h_s,
	double *x_s,double *v_s,double *f_s,

	int n_b,
	double *x_b,double *v_b,double *f_b,

	int n_e,int *nod_e,

	//double t, // t because of printing force
	double dt,

	int dim)
{
	int k,count_slave=0,count_master=0;

  	/* setting contact forces to zeros */
  	memset(f_c,0,3*n_c*sizeof(double));
	
  	/* domain versus boundary contact */
  	for (k=0;k<n_c;k++)
    {
		/* domain node versus boundary element */
		if (((type_m[mat_c[k]]<4)||(type_m[mat_c[k]]>6))&& // domain
			((type_m[mat_c[k+n_c]]>3)&&(type_m[mat_c[k+n_c]]<7))&& // boundary
			 (domain_m[mat_c[k]]==domain_m[mat_c[k+n_c]])) // both materials in same domain
		{
			contact_forces(n_c,f_c,0, // 0 = SPH versus FEM
				count_slave,count_master,sw_c[k],seg_c,slave_c,master_c,
				ct_c[k],klin_c[k],knon_c[k],kf_c[k],kd_c[k],m_slave,m_master,
			 	n_s,
			 	x_s,v_s,f_s,
				h_s,
			 	n_b,
			 	x_b,v_b,f_b,
			 	n_e,nod_e,
			 	//t, // t because of printing force
				dt,
			 	dim,
			 	k);
		}
    	/* boundary node versus boundary element */
		else if (((type_m[mat_c[k]]>3)&&(type_m[mat_c[k]]<7))&& // boundary
				 ((type_m[mat_c[k+n_c]]>3)&&(type_m[mat_c[k+n_c]]<7))&& // boundary
				  (domain_m[mat_c[k]]==domain_m[mat_c[k+n_c]])) // both materials in same domain
		{
			contact_forces(n_c,f_c,1, // 1 = FEM versus FEM
				count_slave,count_master,sw_c[k],seg_c,slave_c,master_c,
				ct_c[k],klin_c[k],knon_c[k],kf_c[k],kd_c[k],m_slave,m_master,
				n_b,
				x_b,v_b,f_b,
				h_s,
				n_b,
				x_b,v_b,f_b,
				n_e,nod_e,
				//t, // t because of printing force
				dt,
				dim,
				k);
		}
		else
		{
			/* domain node versus domain node -> acceleration.c */
		}
    	count_slave=count_slave+seg_c[k];
    	count_master=count_master+seg_c[k+n_c];
    }
}
