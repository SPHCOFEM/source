#include "header.h"
#include "kernel.h"
#include "interpolation.h"

#undef __FUNC__
#define __FUNC__ "accelerations"

void accelerations(int *type_m,int *domain_m,
	double *rho_m,double *T_m,double *mu_m,double *kappa_m,

	int n_u,int *num_u,int *fun_u,int *fun_loc,double *fxi_u,double *fyi_u,

	int n_c,int *mat_c,int *sw_c,int *seg_c,int *slave_c,int *master_c,
	double *ct_c,double *klin_c,double *knon_c,double *kf_c,double *kd_c,

	int n_s,int *mat_s,
	double *m_s,double *mu_s,double *x_s_0,double *h_s_0,
	double *p_s,double *c_s,double *h_s,
	double *dv_s,double *O_s,
	//double *gradvx_s,double *gradvy_s,double *gradvz_s,

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
	
	double ax,double ay,double az)
{
  	int i,j,k,count_neig=0;
  	int nn1,nn2;
	int contact,icontact,contact_num;
	double de,ct,klin,knon,kf,p,vr,fnx,fny,fnz,fn,ftx=0.0,fty=0.0,ftz=0.0,ft;
	double tijx,tijy,tijz,sijx,sijy,sijz,sij;
  	double hij,hij0,muij;
	double cij,rhoij;
  	double rijx,rijy,rijz,rij,rij0;
  	double v,vijx,vijy,vijz,vij;
  	double w,w0,dw,d2w,gradwx,gradwy,gradwz,grad2wx,grad2wy,grad2wz;
  	double afx,afy,afz,vfx,vfy,vfz,pfx,pfy,pfz,rfx,rfy,rfz;
	double Pijx,Pijy,Pijz,Fijx,Fijy,Fijz;
  	double vijrij,av,as;
  	double *gradt,*grad2t,gradtn;
	double *gradvx,*gradvy,*gradvz;
	double a1,a2,a3,f1,f2,f3;

	/* set particle forces to zero 
	-> set in sphcofem.c */
	//memset(f_s,0,3*n_s*sizeof(double));

	/* setting velocity gradients to zero */
	gradvx=createMemMore(double,3*n_s);
	gradvy=createMemMore(double,3*n_s);
	gradvz=createMemMore(double,3*n_s);
	memset(gradvx,0,3*n_s*sizeof(double));
	memset(gradvy,0,3*n_s*sizeof(double));
	memset(gradvz,0,3*n_s*sizeof(double));

	/* setting tension gradients to zero */
	gradt=createMemMore(double,3*n_s);
	grad2t=createMemMore(double,3*n_s);
	memset(gradt,0,3*n_s*sizeof(double));
	memset(grad2t,0,3*n_s*sizeof(double));

	/* setting deformation rate and 
	   stress tensor rate to zero */
	memset(dedt_s,0,6*n_s*sizeof(double));
	memset(dSdt_s,0,6*n_s*sizeof(double));
	
	/* setting domain derivatives to zero */
	memset(drhodt_s,0,n_s*sizeof(double));
	memset(dudt_s,0,n_s*sizeof(double));

	/* xsph */
	if (ix) memset(dv_s,0,3*n_s*sizeof(double));

	/* set v_max to very small number */
	*v_max=1e-9;

  	/* acceleration of domain particles */
  	for (i=0;i<n_s;i++)
  	{
		/* maximum velocity for corrector movement */
		v=norm(v_s[i],v_s[i+n_s],v_s[i+2*n_s]);
		if (v>*v_max) *v_max=v;

		/* optimized nearest neighbour search */
		if (nnopt) /* optimized search */
		{
			nn1=0;
			nn2=nn_s[i];
		}
    	else /* full search */
  		{
			nn1=i+1;
			nn2=n_s;
		}

    	/* loop over nearest neighbours */
    	for (k=nn1;k<nn2;k++)
		{
			/* optimized nearest neighbour search */
			if (nnopt) j=neig_s[k+count_neig]; // optimized search
			else j=k; // full search

			/* particles in same domain */
			if (domain_m[mat_s[i]]==domain_m[mat_s[j]])
			{
				/* partciles mutual position vector */
				rijx=x_s[i]-x_s[j];
				rijy=x_s[i+n_s]-x_s[j+n_s];
				rijz=x_s[i+2*n_s]-x_s[j+2*n_s];

				/* partciles mutual velocity vector */
				vijx=v_s[i]-v_s[j];
				vijy=v_s[i+n_s]-v_s[j+n_s];
				vijz=v_s[i+2*n_s]-v_s[j+2*n_s];

				rij=norm(rijx,rijy,rijz);
				vij=norm(vijx,vijy,vijz);
				hij=(h_s[i]+h_s[j])/2.0;

				cij=(c_s[i]+c_s[j])/2.0;
				rhoij=(rho_s[i]+rho_s[j])/2.0;

				/* checking contact between 2 domains*/
				icontact=0;
				for (contact=0;contact<n_c;contact++)
				{
					if (((mat_s[i]==mat_c[contact])&&(mat_s[j]==mat_c[contact+n_c]))||
						((mat_s[i]==mat_c[contact+n_c])&&(mat_s[j]==mat_c[contact])))
					{
						icontact=1;
						contact_num=contact; // particles in contact
						break;
					}
				}

				if (!icontact) // no contact between 2 domains defined
				{
					if (rij<2.0*hij) // particles within smoothing length
					{
						w=kernel(rij/hij,hij,sigma,dim);
						dw=gradient_of_kernel(rij/hij,hij,sigma,dim)/hij;
						d2w=second_derivative_of_kernel(rij/hij,hij,sigma,dim)/sqr(hij);
				
						rij0=norm(x_s_0[i]-x_s_0[j],x_s_0[i+n_s]-x_s_0[j+n_s],x_s_0[i+2*n_s]-x_s_0[j+2*n_s]);
						hij0=(h_s_0[i]+h_s_0[j])/2.0;
						w0=kernel(rij0/hij0,hij0,sigma,dim);
						//it might happen that rij is outside hij => w0=0 => NaN in as
						//w0=0.001;

						/* gradient of kernel */
						gradwx=rijx*dw/rij;
						gradwy=rijy*dw/rij;
						gradwz=rijz*dw/rij;

						/* gradient of kernel second derivative */
						grad2wx=rijx*d2w/rij;
						grad2wy=rijy*d2w/rij;
						grad2wz=rijz*d2w/rij;
					
						/* compression (bulk stress) */
						afx=(p_s[i]/sqr(rho_s[i])+p_s[j]/sqr(rho_s[j]))*gradwx;
						afy=(p_s[i]/sqr(rho_s[i])+p_s[j]/sqr(rho_s[j]))*gradwy;
						afz=(p_s[i]/sqr(rho_s[i])+p_s[j]/sqr(rho_s[j]))*gradwz;
				
						/* deviatoric (shear) stress for solids */
						if ((type_m[mat_s[i]]==0)||
							(type_m[mat_s[i]]==7)|| // (type_m[mat_s[i]]==8)||(type_m[mat_s[i]]==9) // solid
							(type_m[mat_s[i]]==27))
						{
							vfx=(S_s[i]/sqr(rho_s[i])+S_s[j]/sqr(rho_s[j]))*gradwx+
								(S_s[i+n_s]/sqr(rho_s[i])+S_s[j+n_s]/sqr(rho_s[j]))*gradwy+
								(S_s[i+2*n_s]/sqr(rho_s[i])+S_s[j+2*n_s]/sqr(rho_s[j]))*gradwz;

							vfy=(S_s[i+n_s]/sqr(rho_s[i])+S_s[j+n_s]/sqr(rho_s[j]))*gradwx+
								(S_s[i+3*n_s]/sqr(rho_s[i])+S_s[j+3*n_s]/sqr(rho_s[j]))*gradwy+
								(S_s[i+4*n_s]/sqr(rho_s[i])+S_s[j+4*n_s]/sqr(rho_s[j]))*gradwz;

							vfz=(S_s[i+2*n_s]/sqr(rho_s[i])+S_s[j+2*n_s]/sqr(rho_s[j]))*gradwx+
								(S_s[i+4*n_s]/sqr(rho_s[i])+S_s[j+4*n_s]/sqr(rho_s[j]))*gradwy+
								(S_s[i+5*n_s]/sqr(rho_s[i])+S_s[j+5*n_s]/sqr(rho_s[j]))*gradwz;
						}
						else /* deviatoric (shear) stress for non-solids if no viscosity model chosen */
						{
							vfx=0.0;
							vfy=0.0;
							vfz=0.0;
						}

						/* viscosity for fluids (Macia, 2011) */
						if ((type_m[mat_s[i]]==1)||(type_m[mat_s[i]]==2)||(type_m[mat_s[i]]==3)|| // fluid
							(type_m[mat_s[i]]==12)||(type_m[mat_s[i]]==13)|| // tension only
							(type_m[mat_s[i]]==23)) // Mie-Grueneisen liquid
						{
							switch (is) // viscosity model switch for materials 1, 2, 3, 12, 13 and 23
							{
								case 1: // second order viscous term
								{
									//vfx=m_s[j]*(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j]*(sqr(rij)+sqr(eta*hij)))*(rijx*gradwx+rijy*gradwy+rijz*gradwz)*vijx;
									//vfy=m_s[j]*(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j]*(sqr(rij)+sqr(eta*hij)))*(rijx*gradwx+rijy*gradwy+rijz*gradwz)*vijy;
									//vfz=m_s[j]*(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j]*(sqr(rij)+sqr(eta*hij)))*(rijx*gradwx+rijy*gradwy+rijz*gradwz)*vijz;
									vfx=(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j]*(sqr(rij)+sqr(eta*hij)))*(rijx*gradwx+rijy*gradwy+rijz*gradwz)*vijx;
									vfy=(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j]*(sqr(rij)+sqr(eta*hij)))*(rijx*gradwx+rijy*gradwy+rijz*gradwz)*vijy;
									vfz=(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j]*(sqr(rij)+sqr(eta*hij)))*(rijx*gradwx+rijy*gradwy+rijz*gradwz)*vijz;
								}
								case 2: // Monaghan, Cleary, Gingold (2006)
								{
									vfx=4.9633/(rho_s[i]*rho_s[j])*4.0*(mu_s[i]*mu_s[j])/(mu_s[i]+mu_s[j])*(rijx*vijx+rijy*vijy+rijz*vijz)/(sqr(rij)+sqr(eta*hij))*gradwx;
									vfy=4.9633/(rho_s[i]*rho_s[j])*4.0*(mu_s[i]*mu_s[j])/(mu_s[i]+mu_s[j])*(rijx*vijx+rijy*vijy+rijz*vijz)/(sqr(rij)+sqr(eta*hij))*gradwy;
									vfz=4.9633/(rho_s[i]*rho_s[j])*4.0*(mu_s[i]*mu_s[j])/(mu_s[i]+mu_s[j])*(rijx*vijx+rijy*vijy+rijz*vijz)/(sqr(rij)+sqr(eta*hij))*gradwz;
								}
								case 3: // Morris et al. (1997)
								{
									vfx=-(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j])*vijx*(rijx*gradwx+rijy*gradwy+rijz*gradwz)/(sqr(rij)+sqr(eta*hij));
									vfy=-(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j])*vijy*(rijx*gradwx+rijy*gradwy+rijz*gradwz)/(sqr(rij)+sqr(eta*hij));
									vfz=-(mu_s[i]+mu_s[j])/(rho_s[i]*rho_s[j])*vijz*(rijx*gradwx+rijy*gradwy+rijz*gradwz)/(sqr(rij)+sqr(eta*hij));
								}
								case 4: // Takeda et al. (1994)
								{
									vfx=(mu_s[i]+mu_s[j])/2.0/rho_s[i]/rho_s[j]*(7.0/3.0*(-2.0*w/sqr(hij))*vijx+((rijx*vijx+rijy*vijy+rijz*vijz)/3.0*rijx+vijx*sqr(rij))*4.0*w/sqr(sqr(hij)));
									vfy=(mu_s[i]+mu_s[j])/2.0/rho_s[i]/rho_s[j]*(7.0/3.0*(-2.0*w/sqr(hij))*vijy+((rijx*vijx+rijy*vijy+rijz*vijz)/3.0*rijy+vijy*sqr(rij))*4.0*w/sqr(sqr(hij)));
									vfz=(mu_s[i]+mu_s[j])/2.0/rho_s[i]/rho_s[j]*(7.0/3.0*(-2.2*w/sqr(hij))*vijz+((rijx*vijx+rijy*vijy+rijz*vijz)/3.0*rijz+vijz*sqr(rij))*4.0*w/sqr(sqr(hij)));
								}
								case 5: // Onderik et al. (2007)
								{
									vfx=(mu_s[i]+mu_s[j])/2.0/rho_s[i]/rho_s[j]*vijx*grad2wx;
									vfy=(mu_s[i]+mu_s[j])/2.0/rho_s[i]/rho_s[j]*vijy*grad2wy;
									vfz=(mu_s[i]+mu_s[j])/2.0/rho_s[i]/rho_s[j]*vijz*grad2wz;
								}
								case 6: // Monaghan and Gingold (1983)
								{
									vfx=(2.0*dim+2.0)*(mu_s[i]+mu_s[j])/rho_s[i]*(vijx*rijx+vijy*rijy+vijz*rijz)/(sqr(rij)+sqr(eta*hij))*gradwx*vijx;
									vfy=(2.0*dim+2.0)*(mu_s[i]+mu_s[j])/rho_s[i]*(vijx*rijx+vijy*rijy+vijz*rijz)/(sqr(rij)+sqr(eta*hij))*gradwy*vijy;
									vfz=(2.0*dim+2.0)*(mu_s[i]+mu_s[j])/rho_s[i]*(vijx*rijx+vijy*rijy+vijz*rijz)/(sqr(rij)+sqr(eta*hij))*gradwz*vijz;
								}
								//default: // artificial viscosity only
							}
						}
					
						// always calculate artificial viscosity
						vijrij=(vijx*rijx+vijy*rijy+vijz*rijz);

						/* artificial viscosity for gas and fluid */
						if (vijrij<0.0) /* particles going against */
						{
							/*
							Hyncik (2002): alpha ~ 1.2 and beta ~ 1.5
							Dalrymple (2005): alpha in <0.01, 0.1> and beta = 0
							*/
							muij=hij*vijrij/(sqr(rij)+sqr(eta*hij));
							av=(beta*sqr(muij)-alpha*muij*cij)/rhoij;
							/*
							Fuller (2010) : av=-alpha*hij*vijrij*(cij-2*hij*vijrij/sqr(rij))/rhoij/sqr(rij); => myalpha = alpha, mybeta = 2 * alpha
							Lobovsky (2009) : av = -ksi1*cij*muij/2/rho+ksi2*sqr(muij)/2/rho => myalpha = ksi1 / 2, mybeta = ksi2 / 2
							*/
						}
						else  /* particles going apart */
						{
							av=0.0;
						}

						/* artificial viscosity term */
						pfx=av*gradwx;
						pfy=av*gradwy;
						pfz=av*gradwz;

						/* artificial stress (optimum nas = 4, zeta = 0.3) */
						if ((type_m[mat_s[i]]==0)||
							(type_m[mat_s[i]]==7)|| // (type_m[mat_s[i]]==8)||(type_m[mat_s[i]]==9)||
							(type_m[mat_s[i]]==27)) /* solid */
						{
							/*
							Gray (2001, p. 6649): if (-p + S) < 0 => system is in compression => as = 0
												  if (-p + S) > 0 => system is in tension => as < -3 / (3 * n + 2) * |p - S| / rho ^ 2

							as is zero for (-p + S) < 0 (Gray, 2001, p. 6649)
							as is negative for (-p + S) > 0 (Gray, 2001, p. 6649)
							
							zeta = 0.2 (Monaghan, 2000)
							(-p + S) < 0 for liquids => zeta = 0.01, otherwise particles have slight tendency to form linear structures (Monaghan, 2000)

							zeta in [0.3, 0.8] (Maurel, 2006)
							nas in {3, 4} (Maurel, 2006)

							av acts againts compression => +
							as acts againts tension => -
							*/
							as=-zeta*((-p_s[i]+S_s[i])/sqr(rho_s[i])+(-p_s[j]+S_s[j])/sqr(rho_s[j]))*pow(w/w0,nas); /* Lobovsky (2007, 2009) */
							if (as>0.0) as=theta*as/zeta; // Monaghan, 2000
						}
						else as=0.0; /* fluid - Lobovsky et al. (2007) */

						/* artificial pressure term */
						rfx=as*gradwx;
						rfy=as*gradwy;
						rfz=as*gradwz;
				
						/* tension only for materials type 12 and 13 */
						if ((type_m[mat_s[i]]==12)||(type_m[mat_s[i]]==13))
						{
							gradt[i]=gradt[i]+m_s[j]/rho_s[j]*c_s[j]*gradwx;
							gradt[i+n_s]=gradt[i+n_s]+m_s[j]/rho_s[j]*c_s[j]*gradwy;
							gradt[i+2*n_s]=gradt[i+2*n_s]+m_s[j]/rho_s[j]*c_s[j]*gradwz;
						
							grad2t[i]=grad2t[i]+m_s[j]/rho_s[j]*c_s[j]*grad2wx;
							grad2t[i+n_s]=grad2t[i+n_s]+m_s[j]/rho_s[j]*c_s[j]*grad2wy;
							grad2t[i+2*n_s]=grad2t[i+2*n_s]+m_s[j]/rho_s[j]*c_s[j]*grad2wz;
						}

						/* right-hand side
						f = - pressure term + stress term =
						= - (pressure + artificial viscosity) + (stress + artificial stress) = 
						= - (pressure + artificial viscosity) + (- pressure + shear stress + artificial stress)
						= - (af + pf) + (vf + rf) = 
						= - (p + pi) + (-p + S + r) =
						= - Pij + Fij
						*/

						/* pressure term */
						Pijx=-m_s[i]*m_s[j]*(afx+pfx);
						Pijy=-m_s[i]*m_s[j]*(afy+pfy);
						Pijz=-m_s[i]*m_s[j]*(afz+pfz);

						/* viscosity term */
						Fijx=m_s[i]*m_s[j]*(vfx+rfx);
						Fijy=m_s[i]*m_s[j]*(vfy+rfy);
						Fijz=m_s[i]*m_s[j]*(vfz+rfz);

						/* inter-particle forces i->j */
						f_s[i]=f_s[i]+Pijx+Fijx;
						f_s[i+n_s]=f_s[i+n_s]+Pijy+Fijy;
						f_s[i+2*n_s]=f_s[i+2*n_s]+Pijz+Fijz;

						/* continuity equation and internal energy i->j */
						drhodt_s[i]=drhodt_s[i]+m_s[j]*(vijx*gradwx+vijy*gradwy+vijz*gradwz);
						dudt_s[i]=dudt_s[i]+m_s[j]*p_s[i]/sqr(rho_s[i])*(vijx*gradwx+vijy*gradwy+vijz*gradwz);

						/* damping in energy term i->j
						dudt_s[i]=dudt_s[i]+m_s[j]*(av+p_s[i]/sqr(rho_s[i]))*(vijx*gradwx+vijy*gradwy+vijz*gradwz);
						*/

						/* velocity change i->j */
						gradvx[i]=gradvx[i]-m_s[j]/rho_s[j]*vijx*gradwx;
						gradvx[i+n_s]=gradvx[i+n_s]-m_s[j]/rho_s[j]*vijx*gradwy;
						gradvx[i+2*n_s]=gradvx[i+2*n_s]-m_s[j]/rho_s[j]*vijx*gradwz;
						gradvy[i]=gradvy[i]-m_s[j]/rho_s[j]*vijy*gradwx;
						gradvy[i+n_s]=gradvy[i+n_s]-m_s[j]/rho_s[j]*vijy*gradwy;
						gradvy[i+2*n_s]=gradvy[i+2*n_s]-m_s[j]/rho_s[j]*vijy*gradwz;
						gradvz[i]=gradvz[i]-m_s[j]/rho_s[j]*vijz*gradwx;
						gradvz[i+n_s]=gradvz[i+n_s]-m_s[j]/rho_s[j]*vijz*gradwy;
						gradvz[i+2*n_s]=gradvz[i+2*n_s]-m_s[j]/rho_s[j]*vijz*gradwz;

						/* inter-particle forces j->i */
						f_s[j]=f_s[j]-Pijx-Fijx;
						f_s[j+n_s]=f_s[j+n_s]-Pijy-Fijy;
						f_s[j+2*n_s]=f_s[j+2*n_s]-Pijz-Fijz;

						/* continuity equation and internal energy j->i */
						drhodt_s[j]=drhodt_s[j]+m_s[i]*(vijx*gradwx+vijy*gradwy+vijz*gradwz);
						dudt_s[j]=dudt_s[j]+m_s[i]*p_s[j]/sqr(rho_s[j])*(vijx*gradwx+vijy*gradwy+vijz*gradwz);

						/* damping in energy term j->i
						dudt_s[j]=dudt_s[j]+m_s[i]*(av+p_s[j]/sqr(rho_s[j]))*(vijx*gradwx+vijy*gradwy+vijz*gradwz);
						*/

						/* velocity change j->i */
						gradvx[j]=gradvx[j]-m_s[i]/rho_s[i]*vijx*gradwx;
						gradvx[j+n_s]=gradvx[j+n_s]-m_s[i]/rho_s[i]*vijx*gradwy;
						gradvx[j+2*n_s]=gradvx[j+2*n_s]-m_s[i]/rho_s[i]*vijx*gradwz;
						gradvy[j]=gradvy[j]-m_s[i]/rho_s[i]*vijy*gradwx;
						gradvy[j+n_s]=gradvy[j+n_s]-m_s[i]/rho_s[i]*vijy*gradwy;
						gradvy[j+2*n_s]=gradvy[j+2*n_s]-m_s[i]/rho_s[i]*vijy*gradwz;
						gradvz[j]=gradvz[j]-m_s[i]/rho_s[i]*vijz*gradwx;
						gradvz[j+n_s]=gradvz[j+n_s]-m_s[i]/rho_s[i]*vijz*gradwy;
						gradvz[j+2*n_s]=gradvz[j+2*n_s]-m_s[i]/rho_s[i]*vijz*gradwz;

						/* xsph */
						if (ix)
						{
							/* Gray (2001, p. 6644)
							dv_s[i]=dv_s[i]-m_s[j]*vijx/rhoij*w;
							dv_s[i+n_s]=dv_s[i+n_s]-m_s[j]*vijy/rhoij*w;
							dv_s[i+2*n_s]=dv_s[i+2*n_s]-m_s[j]*vijz/rhoij*w;

							dv_s[j]=dv_s[j]+m_s[i]*vijx/rhoij*w;
							dv_s[j+n_s]=dv_s[j+n_s]+m_s[i]*vijy/rhoij*w;
							dv_s[j+2*n_s]=dv_s[j+2*n_s]+m_s[i]*vijz/rhoij*w;
							*/
							/* opposite sign */
							dv_s[i]=dv_s[i]+m_s[j]*vijx/rhoij*w;
							dv_s[i+n_s]=dv_s[i+n_s]+m_s[j]*vijy/rhoij*w;
							dv_s[i+2*n_s]=dv_s[i+2*n_s]+m_s[j]*vijz/rhoij*w;

							dv_s[j]=dv_s[j]-m_s[i]*vijx/rhoij*w;
							dv_s[j+n_s]=dv_s[j+n_s]-m_s[i]*vijy/rhoij*w;
							dv_s[j+2*n_s]=dv_s[j+2*n_s]-m_s[i]*vijz/rhoij*w;
						}
					}
				}
				else if ((sw_c[contact_num]>0)&&(!(count_cycle%cycle_contact))) // contact between 2 domains defined by contact frequence
				{
					/* contact thickness
					-> negative contact thickness and SPH to FEM -> smoothing length
					-> otherwise contact thickness defined */
					if (ct_c[contact_num]<0) ct=-2.0*ct_c[contact_num]*hij;
					else ct=ct_c[contact_num];

					/* penetration in direction (rijx,rijy,rijz) */
					p=ct-rij;

					if (p>0) // penetration
					{
						/* linear and nonlinear stiffness */
						klin=klin_c[contact_num];
						knon=knon_c[contact_num];

						/* normal contact force */
						if (klin>EPS) // penalty coefficient given
						{
							fn=klin*p/ct+knon*p*p*p/(ct*ct*ct);
						}
						else // penalty coefficient calculated
						{
							fn=m_s[i]*m_s[j]/(m_s[i]+m_s[j])/sqr(dt)*p;
		
							if (knon>EPS) // non-linear penalty active
							{
								fn=fn+m_s[i]*m_s[j]/sqr(dt)*pow(p,knon);
							}
						}

						/* normal force */
						fnx=fn*rijx/rij;
						fny=fn*rijy/rij;
						fnz=fn*rijz/rij;

						/* friction force */
						if (vijx*rijx+vijy*rijy+vijz*rijz>0) // particles mutual velocity not parallel to particles connection line
						{
							/* friction coefficient  */
							kf=kf_c[contact_num];

							/* tangent friction force */
							ft=kf*fn;

							switch (dim)
							{
								case 1: /* 1D */
								{
									/* no tangential motion */
									break;
								}
								case 2: /* 2D */
								{
									/*  fn acting vector is oposite to rij = (rijx,rijy)
									    ft perpendicular to clock wise fn || (-rijy,rijx)
										projection of (vijx,vijy) to (-rijy,rijx) = (-vijx*rijy+vijy*rijx)/vij/rij */

									/* direction against projection of (vijx,vijy) to (-rijy,rijx) = -(-vijx*rijy+vijy*rijx)/vij/rij */
									if (vij>0)
									{
										vr=(vijx*rijy-vijy*rijx)/vij/rij; // relative velocity ubit vector

										ftx=-vr*ft*rijy/rij;
										fty=vr*ft*rijx/rij;
									}
									break;
								}
								case 3: /* 3D */
								{
									/* rij = (rijx,rijy,rijz) is tangent plane normal vector
									-> vij = (vijx,vijy,vijz) projected to plane defined by normal vector rij
									-> vector product tij = vij x rij lies in plane defined by rij
									-> vector product sij = tij x rij is opposite to projection of vij */

									tijx=vijy*rijz-vijz*rijy;
									tijy=vijz*rijx-vijx*rijz;
									tijz=vijx*rijy-vijy*rijz;

									sijx=tijy*rijz-tijz*rijy;
									sijy=tijz*rijx-tijx*rijz;
									sijz=tijx*rijy-tijy*rijz;

									sij=norm(sijx,sijy,sijz);

									if (sij>0)
									{
										ftx=ft*sijx/sij;
										fty=ft*sijy/sij;
										ftz=ft*sijz/sij;
									}
									break;
								}
							}
						}
						
						/* contact force */
						Fijx=fnx+ftx;
						Fijy=fny+fty;
						Fijz=fnz+ftz;

						/* inter-particle contact forces i->j */
						f_s[i]=f_s[i]+Fijx;
						f_s[i+n_s]=f_s[i+n_s]+Fijy;
						f_s[i+2*n_s]=f_s[i+2*n_s]+Fijz;

						/* inter-particle contact forces j->i */
						f_s[j]=f_s[j]-Fijx;
						f_s[j+n_s]=f_s[j+n_s]-Fijy;
						f_s[j+2*n_s]=f_s[j+2*n_s]-Fijz;
					}
				}
				//else // no interation between particles
				//{
				//	// particles not influencing each other neither by SPH or contact
				//}
			}
		}

		/* deformation rate tensor */
		dedt_s[i]=gradvx[i]; // 1D, 2D, 3D
		if (dim>1) // 2D, 3D
		{
			dedt_s[i+n_s]=0.5*(gradvx[i+n_s]+gradvy[i]);
			dedt_s[i+3*n_s]=gradvy[i+n_s];
		}
		if (dim>2) // 3D
		{
			dedt_s[i+2*n_s]=0.5*(gradvx[i+2*n_s]+gradvz[i]);
			dedt_s[i+4*n_s]=0.5*(gradvy[i+2*n_s]+gradvz[i+n_s]);
			dedt_s[i+5*n_s]=gradvz[i+2*n_s];
		}

		/* deformation rate tensor trace */
		de=(dedt_s[i]+dedt_s[i+3*n_s]+dedt_s[i+5*n_s])/(double)dim;

		/* rotation tensor 
		-> O_s[i]=0.0;O_s[i+5*n_s]=0.0;O_s[i+3*n_s]=0.0; */
		if (dim>1) O_s[i+n_s]=0.5*(gradvx[i+n_s]-gradvy[i]);
		if (dim>2)
		{
			O_s[i+2*n_s]=0.5*(gradvx[i+2*n_s]-gradvz[i]);
			O_s[i+4*n_s]=0.5*(gradvy[i+2*n_s]-gradvz[i+n_s]);
		}

		/* deviatoric stress rate tensor for elastic materials */
		if ((type_m[mat_s[i]]==0)|| // rigid body
			(type_m[mat_s[i]]==7)|| // (type_m[mat_s[i]]==8)||(type_m[mat_s[i]]==9)|| // solid
			(type_m[mat_s[i]]==27)) // Mie-Grueneisen solid
		{
			dSdt_s[i]=-2.0*mu_m[mat_s[i]]*(de-dedt_s[i])+2.0*O_s[i+n_s]*S_s[i+n_s]+2.0*O_s[i+2*n_s]*S_s[i+2*n_s];
			dSdt_s[i+n_s]=2.0*dedt_s[i+n_s]*mu_m[mat_s[i]]+O_s[i+n_s]*S_s[i+3*n_s]-O_s[i+n_s]*S_s[i]+O_s[i+2*n_s]*S_s[i+4*n_s]+O_s[i+4*n_s]*S_s[i+2*n_s];
			dSdt_s[i+2*n_s]=2.0*dedt_s[i+2*n_s]*mu_m[mat_s[i]]+O_s[i+n_s]*S_s[i+4*n_s]-O_s[i+2*n_s]*S_s[i]-O_s[i+4*n_s]*S_s[i+n_s]+O_s[i+2*n_s]*S_s[i+5*n_s];
			dSdt_s[i+3*n_s]=-2*mu_m[mat_s[i]]*(de-dedt_s[i+3*n_s])+2*O_s[i+4*n_s]*S_s[i+4*n_s]-2*O_s[i+n_s]*S_s[i+n_s];
			dSdt_s[i+4*n_s]=2.0*dedt_s[i+4*n_s]*mu_m[mat_s[i]]+O_s[i+4*n_s]*S_s[i+5*n_s]-O_s[i+2*n_s]*S_s[i+n_s]-O_s[i+4*n_s]*S_s[i+3*n_s]-O_s[i+n_s]*S_s[i+2*n_s];
			dSdt_s[i+5*n_s]=-2.0*mu_m[mat_s[i]]*(de-dedt_s[i+5*n_s])-2.0*O_s[i+2*n_s]*S_s[i+2*n_s]-2.0*O_s[i+4*n_s]*S_s[i+4*n_s];
		}

      	/* tension only material */
		if ((type_m[mat_s[i]]==12)||(type_m[mat_s[i]]==13))
		{
	  		gradtn=norm(gradt[i],gradt[i+n_s],gradt[i+2*n_s]);
	  		if (gradtn>EPS)
	    	{
	      		f_s[i]=f_s[i]-kappa_m[mat_s[i]]*grad2t[i]*gradt[i]/gradtn;
	      		f_s[i+n_s]=f_s[i+n_s]-kappa_m[mat_s[i]]*grad2t[i+n_s]*gradt[i+n_s]/gradtn;
	      		f_s[i+2*n_s]=f_s[i+2*n_s]-kappa_m[mat_s[i]]*grad2t[i+2*n_s]*gradt[i+2*n_s]/gradtn;
	    	}
		}
		//else // other possiblke materials
		//{
		//	// other possiblke materials
		//}

		/* particle accelerations */
		a_s[i]=f_s[i]/m_s[i]+ax;
		a_s[i+n_s]=f_s[i+n_s]/m_s[i]+ay;
		a_s[i+2*n_s]=f_s[i+2*n_s]/m_s[i]+az;
		
    	/* optimized nearest neighbour search */
		if (nnopt) count_neig=count_neig+nn_s[i];
    }
    
  	/* free gradients */
	freeMem(gradt);
	freeMem(grad2t);
	freeMem(gradvx);
	freeMem(gradvy);
	freeMem(gradvz);
	
  	/* acceleration field */
  	for (i=0;i<n_a_0_s;i++)
    {
		if (type_a_0_s[i]==0) // constant
		{
			a1=a_0_s[i];
			a2=a_0_s[i+n_a_0_s];
			a3=a_0_s[i+2*n_a_0_s];
		}
		else // function
		{
			a1=interpolation(t-dt,(int)a_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			a2=interpolation(t-dt,(int)a_0_s[i+n_a_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			a3=interpolation(t-dt,(int)a_0_s[i+2*n_a_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
		}
      	a_s[ind_a_0_s[i]]=a_s[ind_a_0_s[i]]+a1;
      	a_s[ind_a_0_s[i]+n_s]=a_s[ind_a_0_s[i]+n_s]+a2;
      	a_s[ind_a_0_s[i]+2*n_s]=a_s[ind_a_0_s[i]+2*n_s]+a3;
    }

  	/* concentrated forces */
  	for (i=0;i<n_f_0_s;i++)
    {
		if (type_f_0_s[i]==0) // constant
		{
			f1=f_0_s[i];
			f2=f_0_s[i+n_f_0_s];
			f3=f_0_s[i+2*n_f_0_s];
		}
		else // function
		{
			f1=interpolation(t-dt,(int)f_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f2=interpolation(t-dt,(int)f_0_s[i+n_f_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f3=interpolation(t-dt,(int)f_0_s[i+2*n_f_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
		}
      	a_s[ind_f_0_s[i]]=a_s[ind_f_0_s[i]]+f1/m_s[ind_a_0_s[i]];
      	a_s[ind_f_0_s[i]+n_s]=a_s[ind_f_0_s[i]+n_s]+f2/m_s[ind_a_0_s[i]];
      	a_s[ind_f_0_s[i]+2*n_s]=a_s[ind_f_0_s[i]+2*n_s]+f3/m_s[ind_a_0_s[i]];
    }
}
