void contact_forces(int n_c,double *f_c,int f_t,

	int count_slave,int count_master,int sw,int *seg_c,int *slave_c,int *master_c,
	double ct,double klin,double knon,double kf,double kd,double *m_slave,double *m_master,

	int n_s,
	double *x_s,double *v_s,double *f_s,
	double *h_s,

	int n_b,
	double *x_b,double *v_b,double *f_b,

	int n_e,int *nod_e,

	//double t,
	double dt,

	int dim,

	int k);
