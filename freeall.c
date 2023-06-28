#include "header.h"

#undef __FUNC__
#define __FUNC__ "freeall"

void freeall(int n_m,int *num_m,int *type_m,int *domain_m,
  double *rho_m,double *mu_m,double *T_m,double *gamma_m,double *kappa_m,
	     
  int n_u,int *num_u,int *fun_u,
  double *fmx_u,double *fmy_u,double *fdx_u,double *fdy_u,double *fxi_u,double *fyi_u,

  int n_c,int *num_c,int *seg_c,int *slave_c,int *master_c,int *sw_c,int *mat_c,
  double *m_slave,double *m_master,double *ct_c,double *klin_c,double *knon_c,double *kf_c,double *kd_c,double *f_c,
	     
  int n_s,int ix,int *num_s,int *mat_s,
  double *x_s_0,double *c_s_0,double *h_s_0,
  double *x_s,double *v_s,double *dv_s,
  //double *x_s_a,double *v_s_a,
  double *x_s_o,double *v_s_o,
  
	//double *gradvx_s,double *gradvy_s,double *gradvz_s,

  double *V_s,double *m_s,
  double *rho_s,double *u_s,double *p_s,double *c_s,double *h_s,double *mu_s,
  //double *rho_s_a,double *u_s_a,//double *p_s_a,double *c_s_a,double *h_s_a,
  double *rho_s_o,double *u_s_o,//double *p_s_o,double *c_s_o,double *h_s_o,

  double *drhodt_s,double *dudt_s,
  double *drhodt_s_o,double *dudt_s_o,

  double *a_s,
  double *a_s_o,
  double *f_s,

  double *e_s,double *S_s,
  //double *e_s_a,double *S_s_a,
  double *e_s_o,double *S_s_o,
  
  double *dedt_s,double *dSdt_s,
  double *dedt_s_o,double *dSdt_s_o,double *O_s,
  //double *O_s_o,

  int n_b,int *num_b,int *mat_b,
  double *x_b_0,
  double *x_b,double *v_b,
  //double *x_b_a,double *v_b_a,
  double *x_b_o,double *v_b_o,

  double *a_b,
  double *a_b_o,
  double *f_b,

  int nnopt,int *neig_s,int *nn_s,

  int n_a_0_s,int *ind_a_0_s,int *type_a_0_s,int *frame_a_0_s,
  double *a_0_s,

  int n_f_0_s,int *ind_f_0_s,int *type_f_0_s,int *frame_f_0_s,
  double *f_0_s,double *M_0_s,

  int n_o_0_s,int *ind_o_0_s,int *type_o_0_s,int *o_0_s,

  int n_a_0_b,int *ind_a_0_b,int *type_a_0_b,int *frame_a_0_b,
  double *a_0_b,

  int n_f_0_b,int *ind_f_0_b,int *type_f_0_b,int *frame_f_0_b,
  double *f_0_b,double *M_0_b,

  int n_o_0_b,int *ind_o_0_b,int *type_o_0_b,int *o_0_b,

  int n_d_0_b,int *ind_d_0_b,
  double *d_0_b,

  int *constrained_s,int *constrained_b,
  int *constrained_r,int *constrained_r_frame,
  //int n_unc_s,int n_unc_b,
  //int *unconstrained_s,int *unconstrained_b,

  int n_r,
  int *num_r,int *type_r,int *mat_r,int *COG_r,int *N1_r,int *N2_r,int *N3_r,
  int *COG_r_type,int *N1_r_type,int *N2_r_type,int *N3_r_type,
	double *m_r,double *I1_r,double *I2_r,double *I3_r,double *D1_r,double *D2_r,double *D3_r,
  double *x_r_0,double *psi_r_0,
  double *x_r,double *v_r,double *psi_r,double *o_r,
  //double *x_r_a,double *v_r_a,double *psi_r_a,double *o_r_a,
  double *x_r_o,double *v_r_o,double *psi_r_o,double *o_r_o,
  //double *x1_r,double *x2_r,double *x3_r,
  //double *x1_r_o,double *x2_r_o,double *x3_r_o,
  double *a_r,double *alpha_r,
  double *a_r_o,double *alpha_r_o,
  double *f_r,double *M_r,
  int *type_a_0_r,int *frame_a_0_r,
  double *a_0_r,
  int *type_f_0_r,int *frame_f_0_r,
  double *f_0_r,double *M_0_r,
	double *u1_r,double *u2_r,double *u3_r,
	int *constrained_s_rb,int *constrained_b_rb,
  int *rb_s,int *rb_b,
	int *ind_r_0_s,int *ind_r_0_b,
  double *l_s,double *l_b,

  int n_n,
  double *M,double *B,double *K,double *f,

  int n_e,
  int *nod_e,int *num_e,int *mat_e,    
  double *V_e)
{
  if (n_u>0)
  {
    freeMem(num_u);
    freeMem(fun_u);
    freeMem(fmx_u);
    freeMem(fmy_u);
    freeMem(fdx_u);
    freeMem(fdy_u);
    freeMem(fxi_u);
    freeMem(fyi_u);
  }

  if (n_m>0)
  {
    freeMem(num_m);
    freeMem(type_m);
    freeMem(domain_m);
    freeMem(rho_m);
    freeMem(mu_m);
    freeMem(T_m);
    freeMem(gamma_m);
    freeMem(kappa_m);
  }

  if (n_s>0)
  {
    freeMem(mat_s);
    freeMem(num_s);
    freeMem(x_s);
    freeMem(x_s_0);
    freeMem(x_s_o);
    //freeMem(x_s_a);
    freeMem(v_s);
    freeMem(v_s_o);
    //freeMem(v_s_a);
    if (ix) freeMem(dv_s);
    //freeMem(gradvx_s);
    //freeMem(gradvy_s);
    //freeMem(gradvz_s);
    freeMem(a_s);
    freeMem(a_s_o);
    freeMem(f_s);
    freeMem(rho_s);
    freeMem(rho_s_o);
    //freeMem(rho_s_a);
    freeMem(drhodt_s);
    freeMem(drhodt_s_o);
    freeMem(u_s);
    freeMem(u_s_o);
    //freeMem(u_s_a);
    freeMem(dudt_s);
    freeMem(dudt_s_o);
    freeMem(mu_s);
    freeMem(V_s);
    freeMem(m_s);
    freeMem(p_s);
    //freeMem(p_s_o);
    //freeMem(p_s_a);
    freeMem(c_s);
    freeMem(c_s_0);
    //freeMem(c_s_o);
    //freeMem(c_s_a);
    freeMem(h_s);
    freeMem(h_s_0);
    //freeMem(h_s_o);
    //freeMem(h_s_a);
    freeMem(S_s);
    freeMem(S_s_o);
    //freeMem(S_s_a);
    freeMem(dSdt_s);
    freeMem(e_s);
    freeMem(e_s_o);
    //freeMem(e_s_a);
    freeMem(O_s);
    //freeMem(O_s_o);
    freeMem(dedt_s);
    freeMem(dedt_s_o);
    if (nnopt)
	  {
	    freeMem(neig_s);
	    freeMem(nn_s);
	  }
    freeMem(constrained_s);
  }

  if (n_f_0_s>0)
  {
    freeMem(ind_f_0_s);
    freeMem(type_f_0_s);
    freeMem(frame_f_0_s);
    freeMem(f_0_s);
    freeMem(M_0_s);
  }

  if (n_a_0_s>0)
  {
    freeMem(ind_a_0_s);
    freeMem(type_a_0_s);
    freeMem(frame_a_0_s);
    freeMem(a_0_s);
  }
  
  if (n_o_0_s>0)
  {
    freeMem(ind_o_0_s);
    freeMem(type_o_0_s);
    freeMem(o_0_s);
  }

  //if (n_unc_s>0) freeMem(unconstrained_s);

  if (n_b>0)
  {
    freeMem(num_b);
    freeMem(mat_b);
    freeMem(x_b);
    freeMem(x_b_0);
    freeMem(x_b_o);
    //freeMem(x_b_a);
    freeMem(v_b);
    freeMem(v_b_o);
    //freeMem(v_b_a);
    freeMem(a_b);
    freeMem(a_b_o);
    freeMem(f_b);
    freeMem(constrained_b);
  }

  if (n_f_0_b>0)
  {
    freeMem(ind_f_0_b);
    freeMem(type_f_0_b);
    freeMem(frame_f_0_b);
    freeMem(f_0_b);
    freeMem(M_0_b);
  }

  if (n_a_0_b>0)
  {
    freeMem(ind_a_0_b);
    freeMem(type_a_0_b);
    freeMem(frame_a_0_b);
    freeMem(a_0_b);
  }
  
  if (n_o_0_b>0)
  {
    freeMem(ind_o_0_b);
    freeMem(type_o_0_b);
    freeMem(o_0_b);
  }

  if (n_o_0_b>0)
  {
    freeMem(ind_d_0_b);
    freeMem(d_0_b);
  }

  //if (n_unc_b>0) freeMem(unconstrained_b);

  if (n_n>0)
  {
    freeMem(M);
    freeMem(B);
    freeMem(K);
    freeMem(f);
  }

  if (n_e>0)
  {
    freeMem(nod_e);
    freeMem(mat_e);
    freeMem(num_e);
    freeMem(V_e);
  }

  if (n_c>0)
  {
    freeMem(seg_c);
    freeMem(slave_c);
    freeMem(master_c);
    freeMem(mat_c);
    freeMem(num_c);
    freeMem(sw_c);
    freeMem(m_slave);
    freeMem(m_master);
    freeMem(ct_c);
    freeMem(klin_c);
    freeMem(knon_c);
    freeMem(kf_c);
    freeMem(kd_c);
    freeMem(f_c);
  }

  if (n_r>0)
  {
    freeMem(m_r);
    freeMem(I1_r);
    freeMem(I2_r);
    freeMem(I3_r);
    freeMem(D1_r);
    freeMem(D2_r);
    freeMem(D3_r);

    freeMem(num_r);
    freeMem(type_r);
    freeMem(mat_r);
    freeMem(COG_r);
    freeMem(COG_r_type);
    freeMem(N1_r);
    freeMem(N1_r_type);
    freeMem(N2_r);
    freeMem(N2_r_type);
    freeMem(N3_r);
    freeMem(N3_r_type);

    freeMem(x_r);
    freeMem(x_r_0);
    freeMem(x_r_o);
    //freeMem(x_r_a);
    freeMem(psi_r);
    freeMem(psi_r_0);
    //freeMem(psi_r_a);
    freeMem(psi_r_o);
    freeMem(v_r);
    freeMem(v_r_o);
    //freeMem(v_r_a);
    freeMem(o_r);
    freeMem(o_r_o);
    //freeMem(o_r_a);
    freeMem(a_r);
    freeMem(a_r_o);
    freeMem(alpha_r);
    freeMem(alpha_r_o);
    freeMem(f_r);
    freeMem(M_r);
    freeMem(a_0_r);
    freeMem(type_a_0_r);
    freeMem(frame_a_0_r);
    freeMem(f_0_r);
    freeMem(M_0_r);
    freeMem(type_f_0_r);
    freeMem(frame_f_0_r);

    //freeMem(x1_r);
    //freeMem(x2_r);
    //freeMem(x3_r);
    freeMem(u1_r);
    freeMem(u2_r);
    freeMem(u3_r);

    freeMem(x_r_0);
    freeMem(psi_r_0);

	  freeMem(constrained_s_rb);
    freeMem(constrained_b_rb);
    freeMem(rb_s);
    freeMem(rb_b);
	  freeMem(l_s);
    freeMem(l_b);
	  freeMem(ind_r_0_s);
    freeMem(ind_r_0_b);

    freeMem(constrained_r);
    freeMem(constrained_r_frame);
  }
}
