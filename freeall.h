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
    double *V_e);
