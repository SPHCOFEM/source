void store(int n_s,
	double *m_s,
	double *x_s,double *v_s,double *dv_s,double *a_s,double *f_s,
	double *rho_s,double *drhodt_s,double *u_s,double *dudt_s,double *p_s,double *c_s,double *h_s,
	double *S_s,double *dSdt_s,double *e_s,double *dedt_s,double *O_s,

	int n_b,
	double *x_b,double *v_b,double *a_b,double *f_b,

    int n_r,
	double *x_r,double *v_r,double *psi_r,double *o_r,
	double *a_r,double *alpha_r,double *f_r,double *M_r,

	int n_c,
	double *f_c,

	double t,double dt,

	double kine_s,double inne_s,double pote_s,
	double kine_b,double pote_b,double defo_b,double disi_b,
	double tote,
	   
	int dim,

	int ix,

	char *output,

	int save_dt,
	   
	int save_kine_s,int save_inne_s,int save_pote_s,
	int save_kine_b,int save_disi_b,int save_defo_b,int save_pote_b,
	int save_tote,

	int save_v_s,int save_dv_s,int save_a_s,int save_f_s,

	int save_rho_s,int save_drhodt_s,int save_u_s,int save_dudt_s,int save_p_s,int save_c_s,int save_h_s,

	int save_S_s,int save_dSdt_s,int save_e_s,int save_dedt_s,int save_O_s,

	int save_v_b,int save_a_b,int save_f_b,
	
	int save_f_c,
	
    int save_x_r,int save_psi_r,int save_v_r,int save_o_r,
	int save_a_r,int save_alpha_r,int save_f_r,int save_M_r);
	