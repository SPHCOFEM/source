void large_deformations(int *type_m,
	
	double *rho_m,double *T_m,double *mu_m,double *kappa_m,double *gamma_m,

	int n_u,int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u,

	int n_b,
	double *x_b,double *x_b_0,double *v_b,double *f_b,

	int n_e,int *nod_e,int *mat_e,

	int n_n,
	double *M,double *B,double *K,double *f,


	int n_a_0_b,int *ind_a_0_b,int *type_a_0_b,
	double *a_0_b,

	int n_f_0_b,int *ind_f_0_b,int *type_f_0_b,
	double *f_0_b,

	int n_d_0_b,int *ind_d_0_b,
	double *d_0_b,

	double t, double dt,

	int dim,
	 
	double c0,double c1,

	double ax,double ay,double az);

void mass_stiffness_damping(int *type_m,
	
	double *rho_m,double *T_m,double *mu_m,double *kappa_m,

	int n_b,
	double *x_b,double *x_b_0,double *v_b,

	int n_e,int *nod_e,int *mat_e,

	int n_n,
	double *M,double *B,double *K,

	double c0, double c1,

	int dim);

void mass_stiffness(int *type_m,
	
	double *rho_m,double *T_m,double *mu_m,double *kappa_m,

	int n_b,
	double *x_b,double *x_b_0,double *v_b,

	int n_e,int *nod_e,int *mat_e,

	int n_n,
	double *M,double *K,

	int dim);

void damping_rhs(int *type_m,
	double *gamma_m,

	int n_u,int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u,

	int n_b,
	double *v_b,double *f_b,

	int n_e,int *nod_e,int *mat_e,

	int n_n,
	double *M,double *B,double *K,double *f,

	int n_a_0_b,int *ind_a_0_b,int *type_a_0_b,
	double *a_0_b,

	int n_f_0_b,int *ind_f_0_b,int *type_f_0_b,
	double *f_0_b,

	int n_d_0_b,int *ind_d_0_b,
	double *d_0_b,

	double t, double dt,

	int dim,

	double c0,double c1,

	double ax,double ay,double az);

void rhs(int *type_m,

	int n_u,int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u,

	int n_b,
	double *f_b,

	int n_e,int *nod_e,int *mat_e,

	int n_n,
	double *M,double *f,

	int n_a_0_b,int *ind_a_0_b,int *type_a_0_b,
	double *a_0_b,

	int n_f_0_b,int *ind_f_0_b,int *type_f_0_b,
	double *f_0_b,

	double t, double dt,

	int dim,

	double ax,double ay,double az);
