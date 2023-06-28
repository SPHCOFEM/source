double time_step(int *type_m,
	double *rho_m,double *T_m,double *mu_m,double *kappa_m,

	int n_s,int *mat_s,
	double *x_s,double *v_s,
	//double *x_s_a,double *v_s_a,
	double *x_s_o,double *v_s_o,

	double *rho_s,double *u_s,double *p_s,double *c_s,double *h_s,
	//double *rho_s_a,double *u_s_a,//double *p_s_a,double *c_s_a,double *h_s_a,
	double *rho_s_o,double *u_s_o,//double *p_s_o,double *c_s_o,double *h_s_o,

	double *a_s,double *drhodt_s,double *dudt_s,
	double *a_s_o,double *drhodt_s_o,double *dudt_s_o,
	
	double *e_s,double *S_s,double *dedt_s,double *dSdt_s,
	//double *e_s_a,double *S_s_a,
	double *e_s_o,double *S_s_o,double *dedt_s_o,double *dSdt_s_o,

	int n_b,
	double *x_b,double *v_b,
	//double *x_b_a,double *v_b_a,
	double *x_b_o,double *v_b_o,

	double *a_b,
	double *a_b_o,

	int n_e,
	int *mat_e,int *nod_e,
	double *V_e,

	int n_r,
	double *x_r,double *psi_r,double *v_r,double *o_r,
	//double *x_r_a,double *psi_r_a,double *v_r_a,double *o_r_a,
	double *x_r_o,double *psi_r_o,double *v_r_o,double *o_r_o,
	double *a_r,double *alpha_r,
	double *a_r_o,double *alpha_r_o,

	int dim,

	double alpha,double beta,

	double cour,double kstab,

	int integration,

	double *v_max,

	double dr,

	double dt_max,

	double dt);
