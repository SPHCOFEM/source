void rigid_bodies(int n_r,int dim,
	int *N1_r,int *N2_r,int *N3_r,
	int *N1_r_type,int *N2_r_type,int *N3_r_type,
	double *x_r,double *u1_r,double *u2_r,double *u3_r,
	double *m_r,double *I1_r,double *I2_r,double *I3_r,double *D1_r,double *D2_r,double *D3_r,
	int *rb_s,int *rb_b,
	int *ind_r_0_s,int *ind_r_0_b,
	double *a_r,double *alpha_r,
	double *f_r,double *M_r,
	int *type_a_0_r,int *frame_a_0_r,
	double *a_0_r,
	int *type_f_0_r,int *frame_f_0_r,
	double *f_0_r,double *M_0_r,	
	int n_s,int n_b,
	double *x_s,double *x_b,
	double *f_s,double *f_b,
	int n_u,int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u,
	double t,double dt,
	double ax,double ay,double az);