#include "header.h"

#undef __FUNC__
#define __FUNC__ "store"

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
    int save_a_r,int save_alpha_r,int save_f_r,int save_M_r)
{
  FILE *sph;

  double *ve=0,*vde=0,*vO=0,*vS=0,*vdS=0;
  int i;

  sph=fopen(output,"a");
  
  fwrite(&(t),sizeof(double),1,sph);
  if (save_dt) fwrite(&(dt),sizeof(double),1,sph);

  if (n_s>0)
  {
    if (save_kine_s) fwrite(&(kine_s),sizeof(double),1,sph);
    if (save_inne_s) fwrite(&(inne_s),sizeof(double),1,sph);
    if (save_pote_s) fwrite(&(pote_s),sizeof(double),1,sph);
  }

  if (n_b>0)
  {
    if (save_kine_b) fwrite(&(kine_b),sizeof(double),1,sph);
    if (save_disi_b) fwrite(&(disi_b),sizeof(double),1,sph);
    if (save_defo_b) fwrite(&(defo_b),sizeof(double),1,sph);
    if (save_pote_b) fwrite(&(pote_b),sizeof(double),1,sph);
  }

  if (save_tote) fwrite(&(tote),sizeof(double),1,sph);

  if (n_s>0)
  {
    fwrite(x_s,sizeof(double),dim*n_s,sph);
    if (save_v_s) fwrite(v_s,sizeof(double),dim*n_s,sph);
    if ((save_dv_s)&&(ix)) fwrite(dv_s,sizeof(double),dim*n_s,sph);
    if (save_a_s) fwrite(a_s,sizeof(double),dim*n_s,sph);
    if (save_f_s) fwrite(f_s,sizeof(double),dim*n_s,sph);
    if (save_rho_s) fwrite(rho_s,sizeof(double),n_s,sph);
    if (save_drhodt_s) fwrite(drhodt_s,sizeof(double),n_s,sph);
    if (save_u_s) fwrite(u_s,sizeof(double),n_s,sph);
    if (save_dudt_s) fwrite(dudt_s,sizeof(double),n_s,sph);
    if (save_p_s) fwrite(p_s,sizeof(double),n_s,sph);
    if (save_c_s) fwrite(c_s,sizeof(double),n_s,sph);
    if (save_h_s) fwrite(h_s,sizeof(double),n_s,sph);

    /* full 3x3 matrix (order of appearance needed)
    if (save_e_s) fwrite(e_s,sizeof(double),6*n_s,sph);
    if (save_dedt_s) fwrite(dedt_s,sizeof(double),6*n_s,sph);
    if (save_O_s) fwrite(O_s,sizeof(double),6*n_s,sph);
    if (save_S_s) fwrite(S_s,sizeof(double),6*n_s,sph);
    if (save_dSdt_s) fwrite(dSdt_s,sizeof(double),6*n_s,sph); */
    switch (dim)
    {
      case 1: // e_s, dedt_s scalars 1x1 [0], O_s, S_s, dSdt_s are zero (not relevant) []
      {
        if (save_e_s) fwrite(e_s,sizeof(double),n_s,sph);
        if (save_dedt_s) fwrite(dedt_s,sizeof(double),n_s,sph);

        break;
      }
      case 2: // e_s, dedt_s, S_s, dSdt_s symmetric 2x2 [0, 1, 3], O_s antisymmetric 2x2 [1] (scalar)
      {
        if (save_e_s) ve=createMemMore(double,3*n_s);
        if (save_dedt_s) vde=createMemMore(double,3*n_s);
        if (save_O_s) vO=createMemMore(double,n_s);
        if (save_S_s) vS=createMemMore(double,3*n_s);
        if (save_dSdt_s) vdS=createMemMore(double,3*n_s);

        for (i=0;i<n_s;i++)
        {
          if (save_e_s) 
          {
            ve[i]=e_s[i];
            ve[i+n_s]=e_s[i+n_s];
            ve[i+2*n_s]=e_s[i+3*n_s];
          }

          if (save_dedt_s)
          {
            vde[i]=dedt_s[i];
            vde[i+n_s]=dedt_s[i+n_s];
            vde[i+2*n_s]=dedt_s[i+3*n_s];
          }

          if (save_O_s) vO[i]=O_s[i+n_s];

          if (save_S_s)
          {
            vS[i]=S_s[i];
            vS[i+n_s]=S_s[i+n_s];
            vS[i+2*n_s]=S_s[i+3*n_s];
          }

          if (save_dSdt_s)
          {
            vdS[i]=dSdt_s[i];
            vdS[i+n_s]=dSdt_s[i+n_s];
            vdS[i+2*n_s]=dSdt_s[i+3*n_s];
          }
        }

        if (save_e_s)
        {
          fwrite(ve,sizeof(double),3*n_s,sph);
          freeMem(ve);
        }
        
        if (save_dedt_s)
        {
          fwrite(vde,sizeof(double),3*n_s,sph);
          freeMem(vde);
        }

        if (save_O_s)
        {
          fwrite(vO,sizeof(double),n_s,sph);
          freeMem(vO);
        }
        
        if (save_S_s)
        {
          fwrite(vS,sizeof(double),3*n_s,sph);
          freeMem(vS);
        } 

        if (save_dSdt_s)
        {
          fwrite(dSdt_s,sizeof(double),3*n_s,sph);
          freeMem(vdS);
        }
        break;
      }
      case 3: // e_s, dedt_s, S_s, dSdt_s symmetric 3x3 [0, 1, 2, 3, 4, 5], O_s antisymmetric 3x3 [1, 2, 4]
      {
        if (save_O_s) vO=createMemMore(double,3*n_s);

        for (i=0;i<n_s;i++)
        {
          if (save_O_s)
          {
            vO[i]=O_s[i+n_s];
            vO[i+n_s]=O_s[i+2*n_s];
            vO[i+2*n_s]=O_s[i+4*n_s];
          }
        }
        
        if (save_e_s) fwrite(e_s,sizeof(double),6*n_s,sph);
        if (save_dedt_s) fwrite(dedt_s,sizeof(double),6*n_s,sph);

        if (save_O_s)
        {
          fwrite(vO,sizeof(double),3*n_s,sph);
          freeMem(vO);
        }

        if (save_S_s) fwrite(S_s,sizeof(double),6*n_s,sph);
        if (save_dSdt_s) fwrite(dSdt_s,sizeof(double),6*n_s,sph);

        break;
      }
    }
  }

  if (n_b>0)
  {
    fwrite(x_b,sizeof(double),dim*n_b,sph);
    if (save_v_b) fwrite(v_b,sizeof(double),dim*n_b,sph);
    if (save_a_b) fwrite(a_b,sizeof(double),dim*n_b,sph);
    if (save_f_b) fwrite(f_b,sizeof(double),dim*n_b,sph);
  }
  
  if (n_r>0)
  {
    switch (dim)
    {
      case 1: // x_r, v_r, a_r scalars 1x1 [0], psi_r, o_r, alpha_r are zero (not relevant) []
      {
        if (save_x_r) fwrite(x_r,sizeof(double),dim*n_r,sph);
        if (save_v_r) fwrite(v_r,sizeof(double),dim*n_r,sph);
        if (save_a_r) fwrite(a_r,sizeof(double),dim*n_r,sph);
        if (save_f_r) fwrite(f_r,sizeof(double),dim*n_r,sph);

        break;
      }
      case 2: // x_r, v_r, a_r vectors 2x1 [0, 1], psi_r, o_r, alpha_r scalars 1x1 [0]
      {
        if (save_x_r) fwrite(x_r,sizeof(double),dim*n_r,sph);
        if (save_psi_r) fwrite(psi_r,sizeof(double),n_r,sph);
        if (save_v_r) fwrite(v_r,sizeof(double),dim*n_r,sph);
        if (save_o_r) fwrite(o_r,sizeof(double),n_r,sph);
        if (save_a_r) fwrite(a_r,sizeof(double),dim*n_r,sph);
        if (save_alpha_r) fwrite(alpha_r,sizeof(double),n_r,sph);
        if (save_f_r) fwrite(f_r,sizeof(double),dim*n_r,sph);
        if (save_M_r) fwrite(M_r,sizeof(double),n_r,sph);

        break;
      }
      case 3: // x_r, v_r, a_r vectors 3x1 [0, 1, 2], psi_r, o_r, alpha_r vectors 3x1 [0, 1, 2]
      {
        if (save_x_r) fwrite(x_r,sizeof(double),dim*n_r,sph);
        if (save_psi_r) fwrite(psi_r,sizeof(double),dim*n_r,sph);
        if (save_v_r) fwrite(v_r,sizeof(double),dim*n_r,sph);
        if (save_o_r) fwrite(o_r,sizeof(double),dim*n_r,sph);
        if (save_a_r) fwrite(a_r,sizeof(double),dim*n_r,sph);
        if (save_alpha_r) fwrite(alpha_r,sizeof(double),dim*n_r,sph);
        if (save_f_r) fwrite(f_r,sizeof(double),dim*n_r,sph);
        if (save_M_r) fwrite(M_r,sizeof(double),dim*n_r,sph);

        break;
      }
    }
  }
  
  if ((n_c>0)&&(save_f_c)) fwrite(f_c,sizeof(double),dim*n_c,sph);
  
  fclose(sph);

  fprintf(stdout,
	  "done"
	  "\n"
	  "\n"
	  );
}
