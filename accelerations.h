void accelerations(int *type_m,int *domain_m,
	double *rho_m,double *T_m,double *mu_m,double *kappa_m,

	int n_u,int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u,

	int n_c,int *mat_c,int *sw_c,int *seg_c,int *slave_c,int *master_c,
	double *ct_c,double *klin_c,double *knon_c,double *kf_c,double *kd_c,

	int n_s,int *mat_s,
	double *m_s,double *mu_s,double *x_s_0,double *h_s_0,
	double *p_s,double *c_s,double *h_s,
	double *dv_s,double *O_s,//double *gradvx_s,double *gradvy_s,double *gradvz_s,

	double *x_s,double *v_s,
	double *a_s,double *f_s,
	double *rho_s,double *u_s,
	double *drhodt_s,double *dudt_s,
	double *S_s,
	double *dedt_s,double *dSdt_s,

	int nnopt,int *nn_s,int *neig_s,

	int n_a_0_s,int *ind_a_0_s,int *type_a_0_s,
	double *a_0_s,

	int n_f_0_s,int *ind_f_0_s,int *type_f_0_s,
	double *f_0_s,

	int dim,

	int is,

	int ix,

	double t,double dt,

	double *v_max,

	double sigma,

	double alpha,double beta,double eta,
		   
	double zeta,double nas,double theta,

	double c0,double c1,

	int count_cycle, int cycle_contact,
	
	double ax,double ay,double az);
