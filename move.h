void move(int *num_m,int *type_m,
	double *rho_m,double *T_m,double *gamma_m,double *kappa_m,double *mu_m,
	double *coef1_m,double *coef2_m,double *coef3_m,double *coef4_m,

	int n_u,int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u,

	int n_s,int *num_s,int *mat_s,
	double *m_s,double *x_s_0,double *c_s_0,
	double *p_s,double *c_s,double *h_s,
	double *dv_s,//double *O_s,//double *gradvx_s,double *gradvy_s,double *gradvz_s,

	double *x_s,double *v_s,double *a_s,
	double *rho_s,double *u_s,
	double *drhodt_s,double *dudt_s,
	double *e_s,double *S_s,
	double *dedt_s,double *dSdt_s,
	
	double *x_s_o,double *v_s_o,double *a_s_o,
	double *rho_s_o,double *u_s_o,
	double *drhodt_s_o,double *dudt_s_o,
	double *e_s_o,double *S_s_o,
	double *dedt_s_o,double *dSdt_s_o,

	int n_b,int *num_b,double *x_b_0,
	double *x_b,double *v_b,double *a_b,
	double *x_b_o,double *v_b_o,double *a_b_o,

	int n_o_0_s,int *o_0_s,int *ind_o_0_s,int *type_o_0_s,
	int n_o_0_b,int *o_0_b,int *ind_o_0_b,int *type_o_0_b,

  	int n_r,int *num_r,int *type_r,
	int *COG_r,int *N1_r,int *N2_r,int *N3_r,
	//int *COG_r_type,
	int *N1_r_type,int *N2_r_type,int *N3_r_type,
  	double *x_r_0,double *psi_r_0,
  	double *x_r,double *v_r,double *psi_r,double *o_r,double *a_r,double *alpha_r,
  	double *x_r_o,double *v_r_o,double *psi_r_o,double *o_r_o,double *a_r_o,double *alpha_r_o,
	//double *x1_r,double *x2_r,double *x3_r,
	//double *x1_r_o,double *x2_r_o,double *x3_r_o,
	double *u1_r,double *u2_r,double *u3_r,
	int *constrained_s_rb,int *constrained_b_rb,
	int *rb_s,int *rb_b,
	int *ind_r_0_s,int *ind_r_0_b,
	//int n_r_0_s,int n_r_0_b,
	double *l_s,double *l_b,

	int dim,

	int ix,
	double xeps,

	int *constrained_s,int *constrained_b,int *constrained_r,
	int *constrained_r_frame,

	int integration,int step,

	double t,double dt,

	double v_max,

	double sigma,

	double h0);
