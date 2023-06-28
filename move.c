#include "header.h"
#include "kernel.h"
#include "interpolation.h"
#include "user_defined_material.h"
#include "spatial_rotation.h"

#undef __FUNC__
#define __FUNC__ "move"

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

	double h0)
{
  	int i,j,k;
	double Eta,En,s,G0;
	double C0,C1,C2,C3,C4,C5;

	double u0,u1,u2;
	double v1,v2;
	double a1;

	int count_rb_s=0,count_rb_b=0;
	double psi,theta,phi,R[9],RO[9],ROO[9];//*R,*RO,*ROO;
	double utx=1.0,uty=0.0,utz=0.0,utn=1.0;
	double vtx=0.0,vty=1.0,vtz=0.0,vtn=1.0;
	double wtx=0.0,wty=0.0,wtz=1.0,wtn=1.0;
	double urx=1.0,ury=0.0,urz=0.0,urn=1.0;
	double vrx=0.0,vry=1.0,vrz=0.0,vrn=1.0;
	double wrx=0.0,wry=0.0,wrz=1.0,wrn=1.0;

	//double dx,dy,dz;
	//double C1,C2,C3,S1,S2,S3;

  	/* test on function output 
	if (step==1)
	{
		FILE *out;

		if (t==0) out=fopen("fun_out.txt","w");
		else out=fopen("fun_out.txt","a");

		int f=3;
		fprintf(out,"%d,%f,%f\n", // out,"f%d(%f) = %f\n",
			f,t,interpolation(t,f,n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u));

		fclose(out);
	}
	*/

  	/* boundary conditions on SPH */
	for (i=0;i<n_o_0_s;i++)
	{
		/* particle index */
		j=ind_o_0_s[i];

		/* boundary condition type */
		switch (type_o_0_s[i])
		{
			case -1: // keeps initial conditions
			{
				if (o_0_s[i]) // x-displacement
				{
					x_s[j]=x_s[j]+v_s[j]*dt; // +0.5*m_s[j]*ax*sqr(dt);
					//v_s[j]=v_s[j]+m_s[j]*ax*dt;
					a_s[j]=0.0; // m_s[j]*ax;
				}

				if (o_0_s[i+n_o_0_s]) // y-displacement
				{
					x_s[j+n_s]=x_s[j+n_s]+v_s[j+n_s]*dt; // +0.5*m_s[j]*ay*sqr(dt);
					//v_s[j+n_s]=v_s[j+n_s]+m_s[j]*ay*dt;
					a_s[j+n_s]=0.0; // m_s[j]*ay;
				}

				if (o_0_s[i+2*n_o_0_s]) // z-displacement
				{
					x_s[j+2*n_s]=x_s[j+2*n_s]+v_s[j+2*n_s]*dt; // +0.5*m_s[j]*az*sqr(dt);
					//v_s[j+2*n_s]=v_s[j+2*n_s]+m_s[j]*az*dt;
					a_s[j+2*n_s]=0.0; // m_s[j]*az;
				}

				if (constrained_s_rb[j]) // particle is affiliated to rigid body
				{
					/* constrained_s_rb[j] happens only for COG particle
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_s_rb[j]-1; // center of gravity index

					/* centre of gravity particle constrained by boundary conditions */
					if (constrained_r[k]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_s[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_s[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_s[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
						
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_s[COG_r[k]+n_s]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_s[COG_r[k]+n_s]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_s[COG_r[k]+n_s]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_s[COG_r[k]+2*n_s]; // else x_r[k+2*n_r] calculated by rigid body dynamics
						v_r[k+2*n_r]=v_s[COG_r[k]+2*n_s]; // else v_r[k+2*n_r] calculated by rigid body dynamics
						a_r[k+2*n_r]=a_s[COG_r[k]+2*n_s]; // else a_r[k+2*n_r] calculated by rigid body dynamics
					}

					if (o_0_s[i+3*n_o_0_s]) // x-angle
					{
						psi_r[k]=psi_r[k]+o_r[k]*dt; // +0.5*a1*sqr(dt);
						//o_r[k]=o_r[k]+a1*dt;
						alpha_r[k]=0.0;
					}

					if (o_0_s[i+4*n_o_0_s]) // y-angle
					{
						psi_r[k+n_r]=psi_r[k+n_r]+o_r[k+n_r]*dt; // +0.5*a1*sqr(dt);
						//o_r[k+n_r]=o_r[k+n_r]+a1*dt;
						alpha_r[k+n_r]=0.0;
					}

					if (o_0_s[i+5*n_o_0_s]) // z-angle
					{
						psi_r[k+2*n_r]=psi_r[k+2*n_r]+o_r[k+2*n_r]*dt; // +0.5*a1*sqr(dt);
						//o_r[k+2*n_r]=o_r[k+2*n_r]+a1*dt;
						alpha_r[k+2*n_r]=0.0;
					}
				}
				break;
			}
			case 1: // displacement boundary condition
			{
				if (o_0_s[i]>0) // x-displacement o_0_s[i](t) (-1 means fixed)
				{
					u0=interpolation(t-dt,o_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
					u1=interpolation(t,o_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
					u2=interpolation(t+dt,o_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

					//x_s[j]=x_s_0[j]+u1; // total Lagrangean
					x_s[j]=x_s[j]+u1-u0; // updated Lagrangean
					v_s[j]=(u2-u1)/dt;
					a_s[j]=(u2-2.0*u1+u0)/sqr(dt); // ((u2-u1)/dt-(u1-u0)/dt)/dt;
				}

				if (o_0_s[i+n_o_0_s]>0) // y-displacement o_0_s[i+n_o_0_s](t) (-1 means fixed)
				{
					u0=interpolation(t-dt,o_0_s[i+n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
					u1=interpolation(t,o_0_s[i+n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
					u2=interpolation(t+dt,o_0_s[i+n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

					//x_s[j+n_s]=x_s_0[j+n_s]+u1; // total Lagrangean
					x_s[j+n_s]=x_s[j+n_s]+u1-u0; // updated Lagrangean
					v_s[j+n_s]=(u2-u1)/dt;
					a_s[j+n_s]=(u2-2.0*u1+u0)/sqr(dt); // ((u2-u1)/dt-(u1-u0)/dt)/dt;
				}

				if (o_0_s[i+2*n_o_0_s]>0) // z-displacement o_0_s[i+2*n_o_0_s](t) (-1 means fixed)
				{
					u0=interpolation(t-dt,o_0_s[i+2*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
					u1=interpolation(t,o_0_s[i+2*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
					u2=interpolation(t+dt,o_0_s[i+2*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

					//x_s[j+2*n_s]=x_s_0[j+2*n_s]+u1; // total Lagrangean
					x_s[j+2*n_s]=x_s[j+2*n_s]+u1-u0; // updated Lagrangean
					v_s[j+2*n_s]=(u2-u1)/dt;
					a_s[j+2*n_s]=(u2-2.0*u1+u0)/sqr(dt); // ((u2-u1)/dt-(u1-u0)/dt)/dt;
				}

				if (constrained_s_rb[j]) // particle is affiliated to rigid body
				{
					/* constrained_s_rb[j] happens only for COG particle
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_s_rb[j]-1; // center of gravity index

					/* centre of gravity particle constrained by boundary conditions */
					if (constrained_r[k]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_s[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_s[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_s[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_s[COG_r[k]+n_s]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_s[COG_r[k]+n_s]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_s[COG_r[k]+n_s]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_s[COG_r[k]+2*n_s]; // else x_r[k+2*n_r] calculated by rigid body dynamics
						v_r[k+2*n_r]=v_s[COG_r[k]+2*n_s]; // else v_r[k+2*n_r] calculated by rigid body dynamics
						a_r[k+2*n_r]=a_s[COG_r[k]+2*n_s]; // else a_r[k+2*n_r] calculated by rigid body dynamics
					}

					if (o_0_s[i+3*n_o_0_s]>0) // x-angle (-1 means fixed)
					{
						u0=interpolation(t-dt,o_0_s[i+3*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
						u1=interpolation(t,o_0_s[i+3*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
						u2=interpolation(t+dt,o_0_s[i+3*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

						//psi_r[k]=psi_r_0[k]+u1; // total Lagrangean
						psi_r[k]=psi_r[k]+u1-u0; // updated Lagrangean
						o_r[k]=(u2-u1)/dt;
						alpha_r[k]=(u2-2.0*u1+u0)/sqr(dt);
					}

					if (o_0_s[i+4*n_o_0_s]>0) // y-angle (-1 means fixed)
					{
						u0=interpolation(t-dt,o_0_s[i+4*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
						u1=interpolation(t,o_0_s[i+4*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
						u2=interpolation(t+dt,o_0_s[i+4*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

						//psi_r[k+n_r]=psi_r_0[k+n_r]+u1; // total Lagrangean
						psi_r[k+n_r]=psi_r[k+n_r]+u1-u0; // updated Lagrangean
						o_r[k+n_r]=(u2-u1)/dt;
						alpha_r[k+n_r]=(u2-2.0*u1+u0)/sqr(dt);
					}

					if (o_0_s[i+5*n_o_0_s]>0) // z-angle (-1 means fixed)
					{
						u0=interpolation(t-dt,o_0_s[i+5*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
						u1=interpolation(t,o_0_s[i+5*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
						u2=interpolation(t+dt,o_0_s[i+5*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

						//psi_r[k+2*n_r]=psi_r_0[k+2*n_r]+u1; // total Lagrangean
						psi_r[k+2*n_r]=psi_r[k+2*n_r]+u1-u0; // updated Lagrangean
						o_r[k+2*n_r]=(u2-u1)/dt;
						alpha_r[k+2*n_r]=(u2-2.0*u1+u0)/sqr(dt);
					}
				}
				break;
			}
			case 2: // velocity boundary condition
			{
				if (o_0_s[i]>0) // x-velocity o_0_s[i](t) (-1 means fixed)
				{
					v1=interpolation(t,o_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
					v2=interpolation(t+dt,o_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

					x_s[j]=x_s[j]+v1*dt;
					v_s[j]=v1;
					a_s[j]=(v2-v1)/dt;
				}

				if (o_0_s[i+n_o_0_s]>0) // y-velocity o_0_s[i+n_o_0_s](t) (-1 means fixed)
				{
					v1=interpolation(t,o_0_s[i+n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
					v2=interpolation(t+dt,o_0_s[i+n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

					x_s[j+n_s]=x_s[j+n_s]+v1*dt;
					v_s[j+n_s]=v1;
					a_s[j+n_s]=(v2-v1)/dt;
				}

				if (o_0_s[i+2*n_o_0_s]>0) // z-velocity o_0_s[i+2*n_o_0_s](t) (-1 means fixed)
				{
					v1=interpolation(t,o_0_s[i+2*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
					v2=interpolation(t+dt,o_0_s[i+2*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

					x_s[j+2*n_s]=x_s[j+2*n_s]+v1*dt;
					v_s[j+2*n_s]=v1;
					a_s[j+2*n_s]=(v2-v1)/dt;
				}

				if (constrained_s_rb[j]) // particle is affiliated to rigid body
				{
					/* constrained_s_rb[j] happens only for COG particle
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_s_rb[j]-1; // center of gravity index

					/* centre of gravity particle constrained by boundary conditions */
					if (constrained_r[k]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_s[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_s[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_s[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_s[COG_r[k]+n_s]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_s[COG_r[k]+n_s]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_s[COG_r[k]+n_s]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_s[COG_r[k]+2*n_s]; // else x_r[k+2*n_r] calculated by rigid body dynamics
						v_r[k+2*n_r]=v_s[COG_r[k]+2*n_s]; // else v_r[k+2*n_r] calculated by rigid body dynamics
						a_r[k+2*n_r]=a_s[COG_r[k]+2*n_s]; // else a_r[k+2*n_r] calculated by rigid body dynamics
					}

					if (o_0_s[i+3*n_o_0_s]>0) // x-rotational velocity (-1 means fixed)
					{
						v1=interpolation(t,o_0_s[i+3*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
						v2=interpolation(t+dt,o_0_s[i+3*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

						psi_r[k]=psi_r[k]+v1*dt;
						o_r[k]=v1;
						alpha_r[k]=(v2-v1)/dt;
					}

					if (o_0_s[i+4*n_o_0_s]>0) // y-rotational velocity (-1 means fixed)
					{
						v1=interpolation(t,o_0_s[i+4*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
						v2=interpolation(t+dt,o_0_s[i+4*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

						psi_r[k+n_r]=psi_r[k+n_r]+v1*dt;
						o_r[k+n_r]=v1;
						alpha_r[k+n_r]=(v2-v1)/dt;
					}

					if (o_0_s[i+5*n_o_0_s]>0) // z-rotational velocity (-1 means fixed)
					{
						v1=interpolation(t,o_0_s[i+5*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
						v2=interpolation(t+dt,o_0_s[i+5*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

						psi_r[k+2*n_r]=psi_r[k+2*n_r]+v1*dt;
						o_r[k+2*n_r]=v1;
						alpha_r[k+2*n_r]=(v2-v1)/dt;
					}
				}
				break;
			}
			case 3: // acceleration boundary condition
			{
				if (o_0_s[i]>0) // x-acceleration o_0_s[i](t) (-1 means fixed)
				{
					a1=interpolation(t,o_0_s[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

					x_s[j]=x_s[j]+v_s[j]*dt+0.5*a1*sqr(dt);
					v_s[j]=v_s[j]+a1*dt;
					a_s[j]=a1;
				}

				if (o_0_s[i+n_o_0_s]>0) // y-acceleration o_0_s[i+n_o_0_s](t) (-1 means fixed)
				{
					a1=interpolation(t,o_0_s[i+n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

					x_s[j+n_s]=x_s[j+n_s]+v_s[j+n_s]*dt+0.5*a1*sqr(dt);
					v_s[j+n_s]=v_s[j+n_s]+a1*dt;
					a_s[j+n_s]=a1;
				}

				if (o_0_s[i+2*n_o_0_s]>0) // z-acceleration o_0_s[i+2*n_o_0_s](t) (-1 means fixed)
				{
					a1=interpolation(t,o_0_s[i+2*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

					x_s[j+2*n_s]=x_s[j+2*n_s]+v_s[j+2*n_s]*dt+0.5*a1*sqr(dt);
					v_s[j+2*n_s]=v_s[j+2*n_s]+a1*dt;
					a_s[j+2*n_s]=a1;
				}

				if (constrained_s_rb[j]) // particle is affiliated to rigid body
				{
					/* constrained_s_rb[j] happens only for COG particle
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_s_rb[j]-1; // center of gravity index

					/* centre of gravity particle constrained by boundary conditions */
					if (constrained_r[k]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_s[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_s[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_s[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_s[COG_r[k]+n_s]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_s[COG_r[k]+n_s]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_s[COG_r[k]+n_s]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_s[COG_r[k]+2*n_s]; // else x_r[k+2*n_r] calculated by rigid body dynamics
						v_r[k+2*n_r]=v_s[COG_r[k]+2*n_s]; // else v_r[k+2*n_r] calculated by rigid body dynamics
						a_r[k+2*n_r]=a_s[COG_r[k]+2*n_s]; // else a_r[k+2*n_r] calculated by rigid body dynamics
					}

					if (o_0_s[i+3*n_o_0_s]>0) // x-rotational acceleration (-1 means fixed)
					{
						a1=interpolation(t,o_0_s[i+3*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

						psi_r[k]=psi_r[k]+o_r[k]*dt+0.5*a1*sqr(dt);
						o_r[k]=o_r[k]+a1*dt;
						alpha_r[k]=a1;
					}

					if (o_0_s[i+4*n_o_0_s]>0) // y-rotational acceleration (-1 means fixed)
					{
						a1=interpolation(t,o_0_s[i+4*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

						psi_r[k+n_r]=psi_r[k+n_r]+o_r[k+n_r]*dt+0.5*a1*sqr(dt);
						o_r[k+n_r]=o_r[k+n_r]+a1*dt;
						alpha_r[k+n_r]=a1;
					}

					if (o_0_s[i+5*n_o_0_s]>0) // z-rotational acceleration (-1 means fixed)
					{
						a1=interpolation(t,o_0_s[i+5*n_o_0_s],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

						psi_r[k+2*n_r]=psi_r[k+2*n_r]+o_r[k+2*n_r]*dt+0.5*a1*sqr(dt);
						o_r[k+2*n_r]=o_r[k+2*n_r]+a1*dt;
						alpha_r[k+2*n_r]=a1;
					}
				}
				break;
			}
			/*
			default: // case 0: fixed boundary condition -> not necessary, nothing moves (default)
			{
				if (o_0_s[i]) // fixed x-displacement
				{
					x_s[j]=x_s_0[j];
					v_s[j]=0.0;
					a_s[j]=0.0;
				}

				if (o_0_s[i+n_o_0_s]) // fixed y-displacement
				{
					x_s[j+n_s]=x_s_0[j+n_s];
					v_s[j+n_s]=0.0;
					a_s[j+n_s]=0.0;
				}

				if (o_0_s[i+2*n_o_0_s]) // fixed z-displacement
				{
					x_s[j+2*n_s]=x_s_0[j+2*n_s];
					v_s[j+2*n_s]=0.0;
					a_s[j+2*n_s]=0.0;
				}

				if (constrained_s_rb[j]) // particle is affiliated to rigid body
				{
					// constrained_s_rb[j] happens only for COG particle
					// because constraints are allowed
					// and checked only on COG
					k=constrained_s_rb[j]-1; // center of gravity index

					// centre of gravity particle constrained by boundary conditions
					if (constrained_r[k]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_s[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_s[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_s[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_s[COG_r[k]+n_s]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_s[COG_r[k]+n_s]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_s[COG_r[k]+n_s]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_s[COG_r[k]+2*n_s]; // else x_r[k+2*n_r] calculated by rigid body dynamics
						v_r[k+2*n_r]=v_s[COG_r[k]+2*n_s]; // else v_r[k+2*n_r] calculated by rigid body dynamics
						a_r[k+2*n_r]=a_s[COG_r[k]+2*n_s]; // else a_r[k+2*n_r] calculated by rigid body dynamics
					}

					if (o_0_s[i+3*n_o_0_s]) // fixed x-angle
					{
						psi_r[k]=psi_r_0[k];
						o_r[k]=0.0;
						alpha_r[k]=0.0;
					}

					if (o_0_s[i+4*n_o_0_s]) // fixed y-angle
					{
						psi_r[k+n_r]=psi_r_0[k+n_r];
						o_r[k+n_r]=0.0;
						alpha_r[k+n_r]=0.0;
					}

					if (o_0_s[i+5*n_o_0_s]) // fixed z-angle
					{
						psi_r[k+2*n_r]=psi_r_0[k+2*n_r];
						o_r[k+2*n_r]=0.0;
						alpha_r[k+2*n_r]=0.0;
					}
				}
				break;
			}
			*/
		}
	}

  	/* boundary conditions on FEM */
	for (i=0;i<n_o_0_b;i++)
	{
		/* node index */
		j=ind_o_0_b[i];

		/* boundary condition type */
		switch (type_o_0_b[i])
		{
			case -1: // keeps initial conditions
			{
				if (o_0_b[i]) // x-displacement
				{
					x_b[j]=x_b[j]+v_b[j]*dt; // +0.5*m_b[j]*ax*sqr(dt);
					//v_b[j]=v_b[j]+m_b[j]*ax*dt;
					a_b[j]=0.0; // m_b[j]*ax;
				}

				if (o_0_b[i+n_o_0_b]) // y-displacement
				{
					x_b[j+n_b]=x_b[j+n_b]+v_b[j+n_b]*dt; // +0.5*m_b[j]*ay*sqr(dt);
					//v_b[j+n_b]=v_b[j+n_b]+m_b[j]*ay*dt;
					a_b[j+n_b]=0.0; // m_b[j]*ay;
				}

				if (o_0_b[i+2*n_o_0_b]) // z-displacement
				{
					x_b[j+2*n_b]=x_b[j+2*n_b]+v_b[j+2*n_b]*dt; // +0.5*m_b[j]*az*sqr(dt);
					//v_b[j+2*n_b]=v_b[j+2*n_b]+m_b[j]*az*dt;
					a_b[j+2*n_b]=0.0; // m_b[j]*az;
				}

				if (constrained_b_rb[j]) // node is affiliated to rigid body
				{
					/* constrained_b_rb[j] happens only for COG node
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_b_rb[j]-1; // center of gravity index

					/* centre of gravity node not constrained by boundary conditions */
					if (constrained_r[i]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_b[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_b[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_b[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_b[COG_r[k]+n_b]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_b[COG_r[k]+n_b]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_b[COG_r[k]+n_b]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_b[COG_r[k]+2*n_b]; // else x_r[k+2*n_r] calculated by rigid body dynamics		
						v_r[k+2*n_r]=v_b[COG_r[k]+2*n_b]; // else v_r[k+2*n_r] calculated by rigid body dynamics		
						a_r[k+2*n_r]=a_b[COG_r[k]+2*n_b]; // else a_r[k+2*n_r] calculated by rigid body dynamics		
					}

					if (o_0_b[i+3*n_o_0_b]) // x-angle
					{
						psi_r[k]=psi_r[k]+o_r[k]*dt; // +0.5*a1*sqr(dt);
						//o_r[k]=o_r[k]+a1*dt;
						alpha_r[k]=0.0;
					}

					if (o_0_b[i+4*n_o_0_b]) // y-angle
					{
						psi_r[k+n_r]=psi_r[k+n_r]+o_r[k+n_r]*dt; // +0.5*a1*sqr(dt);
						//o_r[k+n_r]=o_r[k+n_r]+a1*dt;
						alpha_r[k+n_r]=0.0;
					}

					if (o_0_b[i+5*n_o_0_b]) // z-angle
					{
						psi_r[k+2*n_r]=psi_r[k+2*n_r]+o_r[k+2*n_r]*dt; // +0.5*a1*sqr(dt);
						//o_r[k+2*n_r]=o_r[k+2*n_r]+a1*dt;
						alpha_r[k+2*n_r]=0.0;
					}
				}
				break;
			}
			case 1: // displacement boundary condition
			{
				if (o_0_b[i]>0) // x-displacement o_0_b[i](t) (-1 means fixed)
				{
					u0=interpolation(t-dt,o_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
					u1=interpolation(t,o_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
					u2=interpolation(t+dt,o_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

					//x_b[j]=x_b_0[j]+u1; // total Lagrangean
					x_b[j]=x_b[j]+u1-u0; // updated Lagrangean
					v_b[j]=(u2-u1)/dt;
					a_b[j]=(u2-2.0*u1+u0)/sqr(dt); // ((u2-u1)/dt-(u1-u0)/dt)/dt;
				}

				if (o_0_b[i+n_o_0_b]>0) // y-displacement o_0_b[i+n_o_0_b](t) (-1 means fixed)
				{
					u0=interpolation(t-dt,o_0_b[i+n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
					u1=interpolation(t,o_0_b[i+n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
					u2=interpolation(t+dt,o_0_b[i+n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

					//x_b[j+n_b]=x_b_0[j+n_b]+u1; // total Lagrangean
					x_b[j+n_b]=x_b[j+n_b]+u1-u0; // updated Lagrangean
					v_b[j+n_b]=(u2-u1)/dt;
					a_b[j+n_b]=(u2-2.0*u1+u0)/sqr(dt); // ((u2-u1)/dt-(u1-u0)/dt)/dt;
				}

				if (o_0_b[i+2*n_o_0_b]>0) // z-displacement o_0_b[i+2*n_o_0_b](t) (-1 means fixed)
				{
					u0=interpolation(t-dt,o_0_b[i+2*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
					u1=interpolation(t,o_0_b[i+2*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
					u2=interpolation(t+dt,o_0_b[i+2*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

					//x_b[j+2*n_b]=x_b_0[j+2*n_b]+u1; // total Lagrangean
					x_b[j+2*n_b]=x_b[j+2*n_b]+u1-u0; // updated Lagrangean
					v_b[j+2*n_b]=(u2-u1)/dt;
					a_b[j+2*n_b]=(u2-2.0*u1+u0)/sqr(dt); // ((u2-u1)/dt-(u1-u0)/dt)/dt;
				}

				if (constrained_b_rb[j]) // node is affiliated to rigid body
				{
					/* constrained_b_rb[j] happens only for COG node
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_b_rb[j]-1; // center of gravity index]

					/* centre of gravity node not constrained by boundary conditions */
					if (constrained_r[i]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_b[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_b[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_b[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_b[COG_r[k]+n_b]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_b[COG_r[k]+n_b]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_b[COG_r[k]+n_b]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_b[COG_r[k]+2*n_b]; // else x_r[k+2*n_r] calculated by rigid body dynamics		
						v_r[k+2*n_r]=v_b[COG_r[k]+2*n_b]; // else v_r[k+2*n_r] calculated by rigid body dynamics		
						a_r[k+2*n_r]=a_b[COG_r[k]+2*n_b]; // else a_r[k+2*n_r] calculated by rigid body dynamics		
					}

					if (o_0_b[i+3*n_o_0_b]>0) // x-angle (-1 means fixed)
					{
						u0=interpolation(t-dt,o_0_b[i+3*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
						u1=interpolation(t,o_0_b[i+3*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
						u2=interpolation(t+dt,o_0_b[i+3*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

						//psi_r[k]=psi_r_0[k]+u1; // total Lagrangean
						psi_r[k]=psi_r[k]+u1-u0; // updated Lagrangean
						o_r[k]=(u2-u1)/dt;
						alpha_r[k]=(u2-2.0*u1+u0)/sqr(dt);
					}

					if (o_0_b[i+4*n_o_0_b]>0) // y-angle (-1 means fixed)
					{
						u0=interpolation(t-dt,o_0_b[i+4*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
						u1=interpolation(t,o_0_b[i+4*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
						u2=interpolation(t+dt,o_0_b[i+4*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

						//psi_r[k+n_r]=psi_r_0[k+n_r]+u1; // total Lagrangean
						psi_r[k+n_r]=psi_r[k+n_r]+u1-u0; // updated Lagrangean
						o_r[k+n_r]=(u2-u1)/dt;
						alpha_r[k+n_r]=(u2-2.0*u1+u0)/sqr(dt);
					}

					if (o_0_b[i+5*n_o_0_b]>0) // z-angle (-1 means fixed)
					{
						u0=interpolation(t-dt,o_0_b[i+5*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t - dt)
						u1=interpolation(t,o_0_b[i+5*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t)
						u2=interpolation(t+dt,o_0_b[i+5*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // u(t + dt)

						//psi_r[k+2*n_r]=psi_r_0[k+2*n_r]+u1; // total Lagrangean
						psi_r[k+2*n_r]=psi_r[k+2*n_r]+u1-u0; // updated Lagrangean`
						o_r[k+2*n_r]=(u2-u1)/dt;
						alpha_r[k+2*n_r]=(u2-2.0*u1+u0)/sqr(dt);
					}
				}
				break;
			}
			case 2: // velocity boundary condition
			{
				if (o_0_b[i]>0) // x-velocity o_0_b[i](t) (-1 means fixed)
				{
					v1=interpolation(t,o_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
					v2=interpolation(t+dt,o_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

					x_b[j]=x_b[j]+v1*dt;
					v_b[j]=v1;
					a_b[j]=(v2-v1)/dt;
				}

				if (o_0_b[i+n_o_0_b]>0) // y-velocity o_0_b[i+n_o_0_b](t) (-1 means fixed)
				{
					v1=interpolation(t,o_0_b[i+n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
					v2=interpolation(t+dt,o_0_b[i+n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

					x_b[j+n_b]=x_b[j+n_b]+v1*dt;
					v_b[j+n_b]=v1;
					a_b[j+n_b]=(v2-v1)/dt;
				}

				if (o_0_b[i+2*n_o_0_b]>0) // z-velocity o_0_b[i+2*n_o_0_b](t) (-1 means fixed)
				{
					v1=interpolation(t,o_0_b[i+2*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
					v2=interpolation(t+dt,o_0_b[i+2*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

					x_b[j+2*n_b]=x_b[j+2*n_b]+v1*dt;
					v_b[j+2*n_b]=v1;
					a_b[j+2*n_b]=(v2-v1)/dt;
				}

				if (constrained_b_rb[j]) // node is affiliated to rigid body
				{
					/* constrained_b_rb[j] happens only for COG node
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_b_rb[j]-1; // center of gravity index

					/* centre of gravity node not constrained by boundary conditions */
					if (constrained_r[i]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_b[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_b[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_b[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_b[COG_r[k]+n_b]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_b[COG_r[k]+n_b]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_b[COG_r[k]+n_b]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_b[COG_r[k]+2*n_b]; // else x_r[k+2*n_r] calculated by rigid body dynamics		
						v_r[k+2*n_r]=v_b[COG_r[k]+2*n_b]; // else v_r[k+2*n_r] calculated by rigid body dynamics		
						a_r[k+2*n_r]=a_b[COG_r[k]+2*n_b]; // else a_r[k+2*n_r] calculated by rigid body dynamics		
					}

					if (o_0_b[i+3*n_o_0_b]>0) // x-rotational velocity (-1 means fixed)
					{
						v1=interpolation(t,o_0_b[i+3*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
						v2=interpolation(t+dt,o_0_b[i+3*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

						psi_r[k]=psi_r[k]+v1*dt;
						o_r[k]=v1;
						alpha_r[k]=(v2-v1)/dt;
					}

					if (o_0_b[i+4*n_o_0_b]>0) // y-rotational velocity (-1 means fixed)
					{
						v1=interpolation(t,o_0_b[i+4*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
						v2=interpolation(t+dt,o_0_b[i+4*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

						psi_r[k+n_r]=psi_r[k+n_r]+v1*dt;
						o_r[k+n_r]=v1;
						alpha_r[k+n_r]=(v2-v1)/dt;
					}

					if (o_0_b[i+5*n_o_0_b]>0) // z-rotational velocity (-1 means fixed)
					{
						v1=interpolation(t,o_0_b[i+5*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t)
						v2=interpolation(t+dt,o_0_b[i+5*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // v(t + dt)

						psi_r[k+2*n_r]=psi_r[k+2*n_r]+v1*dt;
						o_r[k+2*n_r]=v1;
						alpha_r[k+2*n_r]=(v2-v1)/dt;
					}
				}
				break;
			}
			case 3: // acceleration boundary condition
			{
				if (o_0_b[i]>0) // x-acceleration o_0_b[i](t) (-1 means fixed)
				{
					a1=interpolation(t,o_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

					x_b[j]=x_b[j]+v_b[j]*dt+0.5*a1*sqr(dt);
					v_b[j]=v_b[j]+a1*dt;
					a_b[j]=a1;
				}

				if (o_0_b[i+n_o_0_b]>0) // y-acceleration o_0_b[i+n_o_0_b](t) (-1 means fixed)
				{
					a1=interpolation(t,o_0_b[i+n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

					x_b[j+n_b]=x_b[j+n_b]+v_b[j+n_b]*dt+0.5*a1*sqr(dt);
					v_b[j+n_b]=v_b[j+n_b]+a1*dt;
					a_b[j+n_b]=a1;
				}

				if (o_0_b[i+2*n_o_0_b]>0) // z-acceleration o_0_b[i+2*n_o_0_b](t) (-1 means fixed)
				{
					a1=interpolation(t,o_0_b[i+2*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

					x_b[j+2*n_b]=x_b[j+2*n_b]+v_b[j+2*n_b]*dt+0.5*a1*sqr(dt);
					v_b[j+2*n_b]=v_b[j+2*n_b]+a1*dt;
					a_b[j+2*n_b]=a1;
				}

				if (constrained_b_rb[j]) // node is affiliated to rigid body
				{
					/* constrained_b_rb[j] happens only for COG node
					   because constraints are allowed
					   and checked only on COG */
					k=constrained_b_rb[j]-1; // center of gravity index

					/* centre of gravity node not constrained by boundary conditions */
					if (constrained_r[i]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_b[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_b[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_b[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_b[COG_r[k]+n_b]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_b[COG_r[k]+n_b]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_b[COG_r[k]+n_b]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_b[COG_r[k]+2*n_b]; // else x_r[k+2*n_r] calculated by rigid body dynamics		
						v_r[k+2*n_r]=v_b[COG_r[k]+2*n_b]; // else v_r[k+2*n_r] calculated by rigid body dynamics		
						a_r[k+2*n_r]=a_b[COG_r[k]+2*n_b]; // else a_r[k+2*n_r] calculated by rigid body dynamics		
					}

					if (o_0_b[i+3*n_o_0_b]>0) // x-rotational acceleration (-1 means fixed)
					{
						a1=interpolation(t,o_0_b[i+3*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

						psi_r[k]=psi_r[k]+o_r[k]*dt+0.5*a1*sqr(dt);
						o_r[k]=o_r[k]+a1*dt;
						alpha_r[k]=a1;
					}

					if (o_0_b[i+4*n_o_0_b]>0) // y-rotational acceleration (-1 means fixed)
					{
						a1=interpolation(t,o_0_b[i+4*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

						psi_r[k+n_r]=psi_r[k+n_r]+o_r[k+n_r]*dt+0.5*a1*sqr(dt);
						o_r[k+n_r]=o_r[k+n_r]+a1*dt;
						alpha_r[k+n_r]=a1;
					}

					if (o_0_b[i+5*n_o_0_b]>0) // fixed z-rotational acceleration (-1 means fixed)
					{
						a1=interpolation(t,o_0_b[i+5*n_o_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u); // a(t)

						psi_r[k+2*n_r]=psi_r[k+2*n_r]+o_r[k+2*n_r]*dt+0.5*a1*sqr(dt);
						o_r[k+2*n_r]=o_r[k+2*n_r]+a1*dt;
						alpha_r[k+2*n_r]=a1;
					}
				}
				break;
			}
			/*
			default: //case 0: fixed boundary condition -> not necessary, nothing moves (default)
			{
				if (o_0_b[i]) // fixed x-displacement
				{
					x_b[j]=x_b_0[j];
					v_b[j]=0.0;
					a_b[j]=0.0;
				}

				if (o_0_b[i+n_o_0_b]) // fixed y-displacement
				{
					x_b[j+n_b]=x_b_0[j+n_b];
					v_b[j+n_b]=0.0;
					a_b[j+n_b]=0.0;
				}

				if (o_0_b[i+2*n_o_0_b]) // fixed z-displacement
				{
					x_b[j+2*n_b]=x_b_0[j+2*n_b];
					v_b[j+2*n_b]=0.0;
					a_b[j+2*n_b]=0.0;
				}

				if (constrained_b_rb[j]) // node is affiliated to rigid body
				{
					// constrained_b_rb[j] happens only for COG node
					// because constraints are allowed
					// and checked only on COG
					k=constrained_b_rb[j]-1; // center of gravity index

					// centre of gravity node not constrained by boundary conditions
					if (constrained_r[i]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k]=x_b[COG_r[k]]; // else x_r[k] calculated by rigid body dynamics
						v_r[k]=v_b[COG_r[k]]; // else v_r[k] calculated by rigid body dynamics
						a_r[k]=a_b[COG_r[k]]; // else a_r[k] calculated by rigid body dynamics
					}
					if (constrained_r[i+n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+n_r]=x_b[COG_r[k]+n_b]; // else x_r[k+n_r] calculated by rigid body dynamics
						v_r[k+n_r]=v_b[COG_r[k]+n_b]; // else v_r[k+n_r] calculated by rigid body dynamics
						a_r[k+n_r]=a_b[COG_r[k]+n_b]; // else a_r[k+n_r] calculated by rigid body dynamics
					}
					if (constrained_r[i+2*n_r]) // this conditon probably redundant as taken from constrained_s
					{
						x_r[k+2*n_r]=x_b[COG_r[k]+2*n_b]; // else x_r[k+2*n_r] calculated by rigid body dynamics		
						v_r[k+2*n_r]=v_b[COG_r[k]+2*n_b]; // else v_r[k+2*n_r] calculated by rigid body dynamics		
						a_r[k+2*n_r]=a_b[COG_r[k]+2*n_b]; // else a_r[k+2*n_r] calculated by rigid body dynamics		
					}

					if (o_0_b[i+3*n_o_0_b]) // fixed x-angle
					{
						psi_r[k]=psi_r_0[k];
						o_r[k]=0.0;
						alpha_r[k]=0.0;
					}

					if (o_0_b[i+4*n_o_0_b]) // fixed y-angle
					{
						psi_r[k+n_r]=psi_r_0[k+n_r];
						o_r[k+n_r]=0.0;
						alpha_r[k+n_r]=0.0;
					}

					if (o_0_b[i+5*n_o_0_b]) // fixed z-angle
					{
						psi_r[k+2*n_r]=psi_r_0[k+2*n_r];
						o_r[k+2*n_r]=0.0;
						alpha_r[k+2*n_r]=0.0;
					}
				}
				break;
			}
			*/
		}
	}
	
	/* move particles not constrained in rigid body 
	-> constrained_s_rb = particle affiliation to rigid body (rigid body number)
	-> type_r[constrained_s_rb] = rigid body type (0 / 1) */
	for (i=0;i<n_s;i++)
    {
		/* particle rigid body number 
		-> k = 0, 1, 2, ... if particle in rigid body
		-> k = -1 otherwise */
		k=constrained_s_rb[i]-1;

		/* corrector step in (2) predictor-corrector averages 
		accelerations from predictor and corrector */
		if ((step==1)&&(integration==2))
		{
			/* particle not constrained by rigid body */
			if (!constrained_s_rb[i])
			{
				/* particle not constrained by boundary conditions */
				if (!constrained_s[i]) a_s[i]=(a_s[i]+a_s_o[i])/2.0;
				if (!constrained_s[i+n_s]) a_s[i+n_s]=(a_s[i+n_s]+a_s_o[i+n_s])/2.0;
				if (!constrained_s[i+2*n_s]) a_s[i+2*n_s]=(a_s[i+2*n_s]+a_s_o[i+2*n_s])/2.0;
			}

			/* density updated only for rigid body type 0 for smoothing
			   -> type_r[k] does not exist when no rigid body present -> second in condition 
			      (constrained_s_rb[i] = rigid body number from 1)
			   -> density not stable for broken dam
			if ((!constrained_s_rb[i])||(!type_r[k]))
			{
				drhodt_s[i]=(drhodt_s[i]+drhodt_s_o[i])/2.0;
				dudt_s[i]=(dudt_s[i]+dudt_s_o[i])/2.0;
			}
			*/
		}

		/* no predictor step for (0) Euler,
		predictor step for (1) central acceleraion, (2) predictor-corrector, (3) predictor-corrector leapfrog
		corrector step for (3) predictor-corrector leapfrog does not update x_s (Gray, Monaghan, Swift, 2006, p. 6653)
		-> corrector step for  (1) central acceleration, (2) predictor-corrector with averaged V_s
			-> for step = 0 except integration = 0 (condition in sphcofem.c)
			-> for step = 1 for integrations = 0, 1 and 2 */
		if ((step!=1)||(integration!=3))
		{
			/* particle not constrained by rigid body */
			if (!constrained_s_rb[i])
			{
				/* particle not constrained by boundary conditions */
				if (!constrained_s[i]) // x-translation
				{
					x_s[i]=x_s[i]+v_s[i]*dt+0.5*a_s[i]*sqr(dt);
					v_s[i]=v_s[i]+a_s[i]*dt;
				}
				if (!constrained_s[i+n_s]) // y-translation
				{
					x_s[i+n_s]=x_s[i+n_s]+v_s[i+n_s]*dt+0.5*a_s[i+n_s]*sqr(dt);
					v_s[i+n_s]=v_s[i+n_s]+a_s[i+n_s]*dt;
				}
				if (!constrained_s[i+2*n_s]) // z-translation
				{
					x_s[i+2*n_s]=x_s[i+2*n_s]+v_s[i+2*n_s]*dt+0.5*a_s[i+2*n_s]*sqr(dt);
					v_s[i+2*n_s]=v_s[i+2*n_s]+a_s[i+2*n_s]*dt;
				}
			}
			
			/* density updated only for rigid body type 0 for smoothing
			   -> type_r[k] does not exist when no rigid body present -> second in condition 
			      (constrained_s_rb[i] = rigid body number from 1) */
			if ((!constrained_s_rb[i])||(!type_r[k])) rho_s[i]=rho_s[i]+drhodt_s[i]*dt;
		}
		/* corrector step in (3) predictor-corrector updates only particles velocities
		   and densities from predictor by dV-dV_o (Gray, Monaghan, Swift, 2006, p. 6653) */
		else // if ((step==1)&&(integration==3))
		{
			/* x_s[i] holds from predictor
				x_s[i+n_s] holds from predictor
				x_s[i+2*n_s] holds from predictor */
			/* particle not constrained by rigid body */
			if (!constrained_s_rb[i])
			{
				/* particle not constrained by boundary condition */
				if (!constrained_s[i]) v_s[i]=v_s[i]+(a_s[i]-a_s_o[i])*dt;
				if (!constrained_s[i+n_s]) v_s[i+n_s]=v_s[i+n_s]+(a_s[i+n_s]-a_s_o[i+n_s])*dt;
				if (!constrained_s[i+2*n_s]) v_s[i+2*n_s]=v_s[i+2*n_s]+(a_s[i+2*n_s]-a_s_o[i+2*n_s])*dt;
			}

			/* density updated only for rigid body type 0 for smoothing 
			   -> type_r[k] does not exist when no rigid body present -> second in condition
			      (constrained_s_rb[i] = rigid body number from 1) */
			if ((!constrained_s_rb[i])||(!type_r[k])) rho_s[i]=rho_s[i]+(drhodt_s[i]-drhodt_s_o[i])*dt;
		}

		/* xsph */
		if (((ix==1)&&(step==0))||((ix==2)&&(step==1))||(ix==3))
		{
			/* particle not constrained by rigid body */
			if (!constrained_s_rb[i])
			{
				/* particle not constrained by boundary conditions */
				if (!constrained_s[i]) v_s[i]=v_s[i]+xeps*dv_s[i]; // x-translation
				if (!constrained_s[i+n_s]) v_s[i+n_s]=v_s[i+n_s]+xeps*dv_s[i+n_s]; // y-translation
				if (!constrained_s[i+2*n_s]) v_s[i+2*n_s]=v_s[i+2*n_s]+xeps*dv_s[i+2*n_s]; // z-translation
			}
		}

		/* test on large particle translation */
		if (norm(x_s[i],x_s[i+n_s],x_s[i+2*n_s])-norm(x_s_0[i],x_s_0[i+n_s],x_s_0[i+2*n_s])>MAX)
		{
			fprintf(stdout,
				"ERROR: large translation of domain node number %d!"
				"\n"
				"\n"
				"error termination"
				"\n",
				num_s[i]);
			exit(0);
		}
	
		/* test on negative density */
		if (rho_s[i]<=0.0)
		{
			fprintf(stdout,
				"ERROR: domain node number %d has zero or negative density!"
				"\n"
				"\n"
				"error termination"
				"\n",
				num_s[i]);
			exit(0);
		}

		/* rigid body type 0
		   0 -> state variables updated (for smoothing)
		   1 -> state variables not updated 
		-> type_r[k] does not exist when no rigid body present -> second in condition
		   (constrained_s_rb[i] = rigid body number from 1) */
		if ((!constrained_s_rb[i])||(!type_r[k]))
		{
			/* smoothing length */
			h_s[i]=h0*pow(m_s[i]/rho_s[i],1.0/(double)dim);

			/* deformation tensor 
			e_s[i]=e_s_o[i]+dedt_s[i]*dt;
			e_s[i+n_s]=e_s_o[i+n_s]+dedt_s[i+n_s]*dt;
			e_s[i+2*n_s]=e_s_o[i+2*n_s]+dedt_s[i+2*n_s]*dt;
			e_s[i+3*n_s]=e_s_o[i+3*n_s]+dedt_s[i+3*n_s]*dt;
			e_s[i+4*n_s]=e_s_o[i+4*n_s]+dedt_s[i+4*n_s]*dt;
			e_s[i+5*n_s]=e_s_o[i+5*n_s]+dedt_s[i+5*n_s]*dt;
			*/
			e_s[i]=e_s[i]+dedt_s[i]*dt;
			e_s[i+n_s]=e_s[i+n_s]+dedt_s[i+n_s]*dt;
			e_s[i+2*n_s]=e_s[i+2*n_s]+dedt_s[i+2*n_s]*dt;
			e_s[i+3*n_s]=e_s[i+3*n_s]+dedt_s[i+3*n_s]*dt;
			e_s[i+4*n_s]=e_s[i+4*n_s]+dedt_s[i+4*n_s]*dt;
			e_s[i+5*n_s]=e_s[i+5*n_s]+dedt_s[i+5*n_s]*dt;

			/* deviatoric stress tensor 
			S_s[i]=S_s_o[i]+dSdt_s[i]*dt;
			S_s[i+n_s]=S_s_o[i+n_s]+dSdt_s[i+n_s]*dt;
			S_s[i+2*n_s]=S_s_o[i+2*n_s]+dSdt_s[i+2*n_s]*dt;
			S_s[i+3*n_s]=S_s_o[i+3*n_s]+dSdt_s[i+3*n_s]*dt;
			S_s[i+4*n_s]=S_s_o[i+4*n_s]+dSdt_s[i+4*n_s]*dt;
			S_s[i+5*n_s]=S_s_o[i+5*n_s]+dSdt_s[i+5*n_s]*dt;
			*/
			S_s[i]=S_s[i]+dSdt_s[i]*dt;
			S_s[i+n_s]=S_s[i+n_s]+dSdt_s[i+n_s]*dt;
			S_s[i+2*n_s]=S_s[i+2*n_s]+dSdt_s[i+2*n_s]*dt;
			S_s[i+3*n_s]=S_s[i+3*n_s]+dSdt_s[i+3*n_s]*dt;
			S_s[i+4*n_s]=S_s[i+4*n_s]+dSdt_s[i+4*n_s]*dt;
			S_s[i+5*n_s]=S_s[i+5*n_s]+dSdt_s[i+5*n_s]*dt;

			/* internal energy */
			u_s[i]=u_s[i]+(dudt_s[i]+S_s[i]*dedt_s[i]+2*S_s[i+n_s]*dedt_s[i+2*n_s]+2*S_s[i+2*n_s]*dedt_s[i+2*n_s]+
				S_s[i+3*n_s]*dedt_s[i+3*n_s]+2*S_s[i+4*n_s]*dedt_s[i+4*n_s]+2*S_s[i+5*n_s]*dedt_s[i+5*n_s])*dt;

			/* equation of state for different types of materials
			EOS depends on v_max for materials 3, 13, 23 and 27
			-> additional loop is necessary for calculating v_max
			
			v_max is calculated (only for particles)
			-> at time_step.c for predictor movement
			-> at acceleration.c for corrector movement */

			/* state quantities calculated regardless presence
			in rigid bodies because of neighbour smoothing*/
			switch (type_m[mat_s[i]])
			{
				case 0: /* rigid body */
				{
					/* same as material type 7 for rigid body type 0 */
					if (!type_r[k]) // rigid body type 0 as material type 7
						p_s[i]=kappa_m[mat_s[i]]+rho_m[mat_s[i]]*sqr(c_s_0[i])/gamma_m[mat_s[i]]*
							(pow(rho_s[i]/rho_m[mat_s[i]],gamma_m[mat_s[i]])-1.0); // Monaghan (2000)
					else p_s[i]=0.0; // rigid body type 1 keeps initial conditions
					break;
				}
				case 1: /* gas */
				{
					c_s[i]=sqrt((kappa_m[mat_s[i]]-1.0)*u_s[i]);
					p_s[i]=rho_s[i]*(kappa_m[mat_s[i]]-1.0)*u_s[i];
					break;
				}
				case 2: /* liquid */
				{
					/* initial sound speed calculated in the beginning */
					//c_s[i]=gamma_m[mat_s[i]]/rho_m[mat_s[i]];
					//c_s[i]=sqrt(T_m[mat_s[i]]*gamma_m[mat_s[i]]/rho_m[mat_s[i]]);
					//c_s[i]=sqrt(T_m[mat_s[i]]/rho_m[mat_s[i]]);
					
					/* gamma = 7 (Dalrymple, 2005) - reference value in most cases
						gamma = 2 (Ellero, 2007) */
					p_s[i]=kappa_m[mat_s[i]]+T_m[mat_s[i]]*(pow((rho_s[i]/rho_m[mat_s[i]]),gamma_m[mat_s[i]])-1.0);
					break;
				}
				case 3: /* SPH */
				{
					c_s[i]=10.0*v_max;
					p_s[i]=kappa_m[mat_s[i]]+100.0*rho_m[mat_s[i]]*sqr(v_max)/gamma_m[mat_s[i]]*
						(pow((rho_s[i]/rho_m[mat_s[i]]),gamma_m[mat_s[i]])-1.0);
					/* Lobovsky
					p_s[i]=100.0*rho_m[mat_s[i]]/gamma_m[mat_s[i]]*
						(pow((rho_s[i]/rho_m[mat_s[i]]),gamma_m[mat_s[i]])-1.0); */
					break;
				}
				case 7: /* SPS (Gray, 2001) */
				{
					/* SPS = SPH with shear (linear elastodynamics) 
					-> initial sound speed calculated in the beginning */
					//p_s[i]=kappa_m[mat_s[i]]+sqr(c_s_0[i])*(rho_s[i]-rho_m[mat_s[i]]);
					p_s[i]=kappa_m[mat_s[i]]+rho_m[mat_s[i]]*sqr(c_s_0[i])/gamma_m[mat_s[i]]*
						(pow(rho_s[i]/rho_m[mat_s[i]],gamma_m[mat_s[i]])-1.0); // Monaghan (2000)
					break;
				}
				case 8: /* Neo-Hookean */
				{
					// not implemented yet
					break;
				}
				case 9: /* used-defined */
				{
					user_defined_material(dim,num_m,type_m,
						rho_m,T_m,gamma_m,kappa_m,mu_m,
						p_s,
						e_s,S_s,//O_s,
						dedt_s,dSdt_s,
						coef1_m,coef2_m,coef3_m,coef4_m);
					break;
				}
				case 12: /* liquid with tension */
				{
					/* same as material type 2, but coefficient kappa is used for tension modulus */
					c_s[i]=sqrt(T_m[mat_s[i]]/rho_m[mat_s[i]]);
					p_s[i]=T_m[mat_s[i]]*(pow((rho_s[i]/rho_m[mat_s[i]]),gamma_m[mat_s[i]])-1.0);
					break;
				}
				case 13: /* SPH with tension */
				{
					/* same as material type 3, but coefficient kappa is used for tension modulus */
					c_s[i]=10.0*v_max;
					p_s[i]=100.0*rho_m[mat_s[i]]*sqr(v_max)/gamma_m[mat_s[i]]*(pow((rho_s[i]/rho_m[mat_s[i]]),gamma_m[mat_s[i]])-1.0);
					break;
				}
				case 23: // Mie-Grueneisen EOS for liquid
				case 27: // Mie-Grueneisen EOS for solid
				{
					c_s[i]=10.0*v_max; // sound speed
					s=coef1_m[mat_s[i]]; // linear coefficient of Hugoniot (1.87)
					G0=coef2_m[mat_s[i]]; // Grueneisen constant (0.17)

					// EOS coefficients
					C0=kappa_m[mat_s[i]]; // initial pressure
					C1=rho_m[mat_s[i]]*sqr(c_s_0[i])-G0/2.0*kappa_m[mat_s[i]];
					C2=rho_m[mat_s[i]]*sqr(c_s_0[i])*(2.0*s-1.0)-G0/2.0*rho_m[mat_s[i]]*sqr(c_s_0[i]);
					C3=rho_m[mat_s[i]]*sqr(c_s_0[i])*(2.0*s-1.0)*(3.0*s-1.0)-G0/2.0*rho_m[mat_s[i]]*sqr(c_s_0[i])*(2.0*s-1.0);
					C4=G0;C5=G0;En=0.0; // energy per unit of initial volume
					
					//Eta=1.0-rho_s[i]/rho_m[mat_s[i]]; // (Awoukeng-Goumtcha, 2014)
					Eta=rho_s[i]/rho_m[mat_s[i]]-1.0; // (Taddei, 2015; Frissane, 2019)

					/* (Taddei, 2015; Frissane, 2019) */
					if (Eta>=0) p_s[i]=C0+C1*Eta+C2*sqr(Eta)+C3*pow(Eta,3.0)+(C4+C5*Eta)*En;
					else p_s[i]=C0+C1*Eta;
					break;
				}
			}
		}
	}
  
  	/* move boundary */
	for (i=0;i<n_b;i++)
    {
		/* node not constrained in rigid body
		-> constrained_b_rb = node affiliation to rigid body (rigid body number)
		-> type_r[constrained_s_rb] = rigid body type (0 / 1) */
		if (!constrained_b_rb[i])
		{
			/* corrector step in (2) predictor-corrector averages 
			   accelerations from predictor and corrector */
			if ((step==1)&&(integration==2))
			{
				/* node not constrained by boundary conditions */
				if (!constrained_b[i]) a_b[i]=(a_b[i]+a_b_o[i])/2.0;
				if (!constrained_b[i+n_b]) a_b[i+n_b]=(a_b[i+n_b]+a_b_o[i+n_b])/2.0;
				if (!constrained_b[i+2*n_b]) a_b[i+2*n_b]=(a_b[i+2*n_b]+a_b_o[i+2*n_b])/2.0;
			}

			/* no predictor step for (0) Euler,
			   predictor step for (1) central acceleraion, (2) predictor-corrector, (3) predictor-corrector leapfrog
			   corrector step for (3) predictor-corrector leapfrog does not update x_s (Gray, Monaghan, Swift, 2006, p. 6653)
			   -> corrector step for  (1) central acceleration, (2) predictor-corrector with averaged V_b
			      -> for step = 0 except integration = 0 (condition in sphcofem.c)
				  -> for step = 1 for integrations = 0, 1 and 2 */
			if ((step!=1)||(integration!=3))
			{
				/* node not constrained by boundary conditions */
				if (!constrained_b[i]) // x-translation
				{
					x_b[i]=x_b[i]+v_b[i]*dt+0.5*a_b[i]*sqr(dt);
					v_b[i]=v_b[i]+a_b[i]*dt;
				}
				if (!constrained_b[i+n_b]) // y-translation
				{
					x_b[i+n_b]=x_b[i+n_b]+v_b[i+n_b]*dt+0.5*a_b[i+n_b]*sqr(dt);
					v_b[i+n_b]=v_b[i+n_b]+a_b[i+n_b]*dt;
				}
				if (!constrained_b[i+2*n_b]) // z-translation
				{
					x_b[i+2*n_b]=x_b[i+2*n_b]+v_b[i+2*n_b]*dt+0.5*a_b[i+2*n_b]*sqr(dt);
					v_b[i+2*n_b]=v_b[i+2*n_b]+a_b[i+2*n_b]*dt;
				}
			}
			/* corrector step in (3) predictor-corrector updates only nodal velocities
			   from predictor (Gray, Monaghan, Swift, 2006, p. 6653) */
			else // if ((step==1)&&(integration==3))
			{
				/* x_b[i] holds from predictor
				   x_b[i+n_b] holds from predictor
				   x_b[i+2*n_b] holds from predictor */
				/* node not constrained by boundary conditions */
				if (!constrained_b[i]) v_b[i]=v_b[i]+(a_b[i]-a_b_o[i])*dt;
				if (!constrained_b[i+n_b]) v_b[i+n_b]=v_b[i+n_b]+(a_b[i+n_b]-a_b_o[i+n_b])*dt;
				if (!constrained_b[i+2*n_b]) v_b[i+2*n_b]=v_b[i+2*n_b]+(a_b[i+2*n_b]-a_b_o[i+2*n_b])*dt;
			}

			/* test on large node translation */
			if (norm(x_b[i],x_b[i+n_b],x_b[i+2*n_b])-norm(x_b_0[i],x_b_0[i+n_b],x_b_0[i+2*n_b])>MAX)
			{
				fprintf(stdout,
					"ERROR: large translation of boundary node number %d!"
					"\n"
					"\n"
					"error termination"
					"\n",
					num_b[i]);
				exit(0);
			}
		}
	}

	/* move rigid bodies including particles and nodes */
	for (i=0;i<n_r;i++)
    {
		/* corrector step in (2) predictor-corrector averages 
		   accelerations from predictor and corrector */
      	if ((step==1)&&(integration==2))
	  	{
			if (!constrained_r[i]) a_r[i]=(a_r[i]+a_r_o[i])/2.0;
            if (!constrained_r[i+n_r]) a_r[i+n_r]=(a_r[i+n_r]+a_r_o[i+n_r])/2.0;
            if (!constrained_r[i+2*n_r]) a_r[i+2*n_r]=(a_r[i+2*n_r]+a_r_o[i+2*n_r])/2.0;
			if (!constrained_r[i+3*n_r]) alpha_r[i]=(alpha_r[i]+alpha_r_o[i])/2.0;
            if (!constrained_r[i+4*n_r]) alpha_r[i+n_r]=(alpha_r[i+n_r]+alpha_r_o[i+n_r])/2.0;
            if (!constrained_r[i+5*n_r]) alpha_r[i+2*n_r]=(alpha_r[i+2*n_r]+alpha_r_o[i+2*n_r])/2.0;
        }

		/* no predictor step for (0) Euler,
		   predictor step for (1) central acceleraion, (2) predictor-corrector, (3) predictor-corrector leapfrog
		   corrector step for (3) predictor-corrector leapfrog does not update x_s (Gray, Monaghan, Swift, 2006, p. 6653)
		   -> corrector step for  (1) central acceleration, (2) predictor-corrector with averaged V_r
			  -> for step = 0 except integration = 0 (condition in sphcofem.c)
			  -> for step = 1 for integrations = 0, 1 and 2 */
      	if ((step!=1)||(integration!=3))
	  	{
			/* centre of gravity not constrained by boundary conditions */
          	if (!constrained_r[i]) // x-translation
			{
				x_r[i]=x_r[i]+v_r[i]*dt+0.5*a_r[i]*sqr(dt);
				v_r[i]=v_r[i]+a_r[i]*dt;
			}
          	if (!constrained_r[i+n_r]) // y-translation
			{
				x_r[i+n_r]=x_r[i+n_r]+v_r[i+n_r]*dt+0.5*a_r[i+n_r]*sqr(dt);
				v_r[i+n_r]=v_r[i+n_r]+a_r[i+n_r]*dt;
			}
          	if (!constrained_r[i+2*n_r]) // z-translation
			{
				x_r[i+2*n_r]=x_r[i+2*n_r]+v_r[i+2*n_r]*dt+0.5*a_r[i+2*n_r]*sqr(dt);
				v_r[i+2*n_r]=v_r[i+2*n_r]+a_r[i+2*n_r]*dt;
			}
          	if (!constrained_r[i+3*n_r]) // x-rotation
			{
				psi_r[i]=psi_r[i]+o_r[i]*dt+0.5*alpha_r[i]*sqr(dt);
				o_r[i]=o_r[i]+alpha_r[i]*dt;
			}
          	if (!constrained_r[i+4*n_r]) // y-rotation
			{
				psi_r[i+n_r]=psi_r[i+n_r]+o_r[i+n_r]*dt+0.5*alpha_r[i+n_r]*sqr(dt);
				o_r[i+n_r]=o_r[i+n_r]+alpha_r[i+n_r]*dt;
			}
          	if (!constrained_r[i+5*n_r]) // z-rotation
			{
				psi_r[i+2*n_r]=psi_r[i+2*n_r]+o_r[i+2*n_r]*dt+0.5*alpha_r[i+2*n_r]*sqr(dt);
				o_r[i+2*n_r]=o_r[i+2*n_r]+alpha_r[i+2*n_r]*dt;
			}
        }
      	/* corrector step in (3) predictor-corrector updates only centre of gravity velocities 
		   from predictor (Gray, Monaghan, Swift, 2006, p. 6653) */
		else // if ((step==1)&&(integration==3))
		{
			/* x_r[i] holds from predictor
			   x_r[i+n_r] holds from predictor
			   x_r[i+2*n_r] holds from predictor 
			   psi_r[i] holds from predictor
			   psi_r[i+n_r] holds from predictor
			   psi_r[i+2*n_r] holds from predictor */
			/* centre of gravity not constrained by boundary conditions */
			if (!constrained_r[i]) v_r[i]=v_r[i]+(a_r[i]-a_r_o[i])*dt;
			if (!constrained_r[i+n_r]) v_r[i+n_r]=v_r[i+n_r]+(a_r[i+n_r]-a_r_o[i+n_r])*dt;
			if (!constrained_r[i+2*n_r]) v_r[i+2*n_r]=v_r[i+2*n_r]+(a_r[i+2*n_r]-a_r_o[i+2*n_r])*dt;
			if (!constrained_r[i+3*n_r]) o_r[i]=o_r[i]+(alpha_r[i]-alpha_r_o[i])*dt;
			if (!constrained_r[i+4*n_r]) o_r[i+n_r]=o_r[i+n_r]+(alpha_r[i+n_r]-alpha_r_o[i+n_r])*dt;
			if (!constrained_r[i+5*n_r]) o_r[i+2*n_r]=o_r[i+2*n_r]+(alpha_r[i+2*n_r]-alpha_r_o[i+2*n_r])*dt;
		}

		/* test on large centre of gravity translation */
		if (norm(x_r[i],x_r[i+n_r],x_r[i+2*n_r])-norm(x_r_0[i],x_r_0[i+n_r],x_r_0[i+2*n_r])>MAX)
		{
			fprintf(stdout,
				"ERROR: large translation of rigid body %d!"
				"\n"
				"\n"
				"error termination"
				"\n",
				num_r[i]);
			exit(0);
		}

		/* test on large centre of gravity rotation */
		if (norm(psi_r[i],psi_r[i+n_r],psi_r[i+2*n_r])-norm(psi_r_0[i],psi_r_0[i+n_r],psi_r_0[i+2*n_r])>MAX)
		{
			fprintf(stdout,
				"ERROR: large rotation of rigid body %d!"
				"\n"
				"\n"
				"error termination"
				"\n",
				num_r[i]);
			exit(0);
		}

		/* local coordinate vectors */
		switch (constrained_r_frame[i])
		{
			case 1: // translation in local coordinate system
			{
				if (dim>1)
				{
					/* translation */
					utx=u1_r[i];uty=u1_r[i+n_r];utz=u1_r[i+2*n_r];
					utn=norm(utx,uty,utz);
					utx=utx/utn;uty=uty/utn;utz=utz/utn;

					/* rotation */
					urx=1.0;ury=0.0;urz=0.0;
				}
				if (dim>2)
				{
					/* translation */
					vtx=u2_r[i];vty=u2_r[i+n_r];vtz=u2_r[i+2*n_r];
					vtn=norm(vtx,vty,vtz);
					vtx=vtx/vtn;vty=vty/vtn;vtz=vtz/vtn;

					wtx=u3_r[i];wty=u3_r[i+n_r];wtz=u3_r[i+2*n_r];
					wtn=norm(wtx,wty,wtz);
					wtx=wtx/wtn;wty=wty/wtn;wtz=wtz/wtn;

					/* rotation */
					vrx=0.0;vry=1.0;vrz=0.0;
					wrx=0.0;wry=0.0;wrz=1.0;
				}
				break;
			}
			case 2: // rotation in local coordinate system
			{
				if (dim>1)
				{
					/* translation */
					utx=1.0;uty=0.0;utz=0.0;

					/* rotation */
					urx=u1_r[i];ury=u1_r[i+n_r];urz=u1_r[i+2*n_r];
					urn=norm(urx,ury,urz);
					urx=urx/urn;ury=ury/urn;urz=urz/urn;
				}
				if (dim>2)
				{
					/* translation */
					vtx=0.0;vty=1.0;vtz=0.0;
					wtx=0.0;wty=0.0;wtz=1.0;

					/* rotation */
					vrx=u2_r[i];vry=u2_r[i+n_r];vrz=u2_r[i+2*n_r];
					vrn=norm(vrx,vry,vrz);
					vrx=vrx/vrn;vry=vry/vrn;vrz=vrz/vrn;

					wrx=u3_r[i];wry=u3_r[i+n_r];wrz=u3_r[i+2*n_r];
					wrn=norm(wrx,wry,wrz);
					wrx=wrx/wrn;wry=wry/wrn;wrz=wrz/wrn;
				}
				break;
			}
			case 3: // translation and rotation in local coordinate system
			{
				if (dim>1)
				{
					/* translation */
					utx=u1_r[i];uty=u1_r[i+n_r];utz=u1_r[i+2*n_r];
					utn=norm(utx,uty,utz);
					utx=utx/utn;uty=uty/utn;utz=utz/utn;

					/* rotation */
					urx=u1_r[i];ury=u1_r[i+n_r];urz=u1_r[i+2*n_r];
					urn=norm(urx,ury,urz);
					urx=urx/urn;ury=ury/urn;urz=urz/urn;
				}
				if (dim>2)
				{
					/* translation */
					vtx=u2_r[i];vty=u2_r[i+n_r];vtz=u2_r[i+2*n_r];
					vtn=norm(vtx,vty,vtz);
					vtx=vtx/vtn;vty=vty/vtn;vtz=vtz/vtn;

					wtx=u3_r[i];wty=u3_r[i+n_r];wtz=u3_r[i+2*n_r];
					wtn=norm(wtx,wty,wtz);
					wtx=wtx/wtn;wty=wty/wtn;wtz=wtz/wtn;

					/* rotation */
					vrx=u2_r[i];vry=u2_r[i+n_r];vrz=u2_r[i+2*n_r];
					vrn=norm(vrx,vry,vrz);
					vrx=vrx/vrn;vry=vry/vrn;vrz=vrz/vrn;

					wrx=u3_r[i];wry=u3_r[i+n_r];wrz=u3_r[i+2*n_r];
					wrn=norm(wrx,wry,wrz);
					wrx=wrx/wrn;wry=wry/wrn;wrz=wrz/wrn;
				}
				break;
			}
			default: // translation and rotation in global coordinate system
			{
				if (dim>1)
				{
					/* translation */
					utx=1.0;uty=0.0;utz=0.0;

					/* rotation */
					urx=1.0;ury=0.0;urz=0.0;
				}
				if (dim>2)
				{
					/* translation */
					vtx=0.0;vty=1.0;vtz=0.0;
					wtx=0.0;wty=0.0;wtz=1.0;

					/* rotation */
					vrx=0.0;vry=1.0;vrz=0.0;
					wrx=0.0;wry=0.0;wrz=1.0;
				}
				break;
			}
		}

		/* to implement rotation (sines and cosines)
		C1=cos(psi_r[i]);C2=cos(psi_r[i+n_r]);C3=cos(psi_r[i+2*n_r]);
		S1=sin(psi_r[i]);S2=sin(psi_r[i+n_r]);S3=sin(psi_r[i+2*n_r]);
		.
		.
		.
		.
		.
		.
		.
		.
		.
		.
		*/

		/*
		Rotace
		-> 1. okolo globalni osy "od ucha k uchu"
		-> 2.  kolem "brada - sesule"

		% matice rotace kolem globalni z
		R  = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
		for i=1:r
			rotX(i,2:end) = (R*X(i,2:end)')'; % (R*X)
		end

		% pak vytahnu z X2 uhel n(smerovy vektor)
		V = rotX(2,2:end);
		n = (V/norm(V))';
		
		% pak rotX = A2*X2
		A2 = cs*(n*n') + cos(fi)*eye(3) + sin(fi)*[0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];    
		for i=1:r
			rotX(i,2:end) = (A2*rotX(i,2:end)')';
		end
		*/

		/* rotation angle */
		psi=psi_r[i]-psi_r_0[i]; // total Lagrangean
		//psi=psi_r[i]-psi_r_o[i]; // updated Lagrangean
		theta=psi_r[i+n_r]-psi_r_0[i+n_r];
		phi=psi_r[i+2*n_r]-psi_r_0[i+2*n_r];

		fprintf(stdout,"%f: %f, %f, %f\n",t,psi,theta,phi);
		
		/* rotation matrices 
		R=createMemMore(double,9);
		RO=createMemMore(double,9);
		ROO=createMemMore(double,9);
		*/
		spatial_rotation(dim,R,RO,ROO,
			psi,theta,phi, // rotations
			o_r[i],o_r[i+n_r],o_r[i+2*n_r], // anglular velocities
			alpha_r[i],alpha_r[i+n_r],alpha_r[i+2*n_r], // anglular accelerations
			urx,ury,urz, // first rotation axis
			vrx,vry,vrz, // second rotation axis
			wrx,wry,wrz); // third rotation axis

		if (i==1)
			fprintf(stdout,"%d, %f, %f, %f | %18.18f , %18.18f , %18.18f | = %f\n",
				step,t,psi,psi_r[i]-psi_r_o[i],urx,ury,urz,norm(urx,ury,urz));

		/*
		if (i==1) fprintf(stdout,
			"----------------------------------------------------------------------\n"
		   	"| %18.18f , %18.18f , %18.18f |\n"
		   	"| %18.18f , %18.18f , %18.18f |\n"
		   	"| %18.18f , %18.18f , %18.18f |\n"
			"----------------------------------------------------------------------\n",
		   	R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8]);
		*/

		/*
		constrained_s_rb ... particles affiliation to rigid body
		constrained_b_rb ... nodes affiliation to rigid body

		l_s ... particles coordinates in rigid body local coordinate system
		l_b ... nodes coordinates in rigid body local coordinate system

		rb_s ... number of particles of rigid body -> particles in ind_r_0_s[j+count_rb_s]
		rb_b ... number of nodes of rigid body -> nodes in ind_r_0_b[j+count_rb_b]
		*/

		/* update all particles and nodes */
		for (j=0;j<rb_s[i];j++) // particles
		{
			k=ind_r_0_s[j+count_rb_s]; // particle index

			/* particle position including COG
			-> COG can be unconstrained particle
			-> x_s in local coordinate system */

			/* total Lagrangean 
			dx=x_s_0[k]-x_r_0[i];
			dy=x_s_0[k+n_s]-x_r_0[i+n_r];
			dz=x_s_0[k+2*n_s]-x_r_0[i+2*n_r];
			*/

			/* position */
			x_s[k]=x_r[i]+R[0]*l_s[3*count_rb_s+j]+
						  R[1]*l_s[3*count_rb_s+j+rb_s[i]]+
						  R[2]*l_s[3*count_rb_s+j+2*rb_s[i]];
			x_s[k+n_s]=x_r[i+n_r]+R[3]*l_s[3*count_rb_s+j]+
								  R[4]*l_s[3*count_rb_s+j+rb_s[i]]+
								  R[5]*l_s[3*count_rb_s+j+2*rb_s[i]];
			x_s[k+2*n_s]=x_r[i+2*n_r]+R[6]*l_s[3*count_rb_s+j]+
									  R[7]*l_s[3*count_rb_s+j+rb_s[i]]+
									  R[8]*l_s[3*count_rb_s+j+2*rb_s[i]];

			/* velocity */
			v_s[k]=v_r[i]+RO[0]*l_s[3*count_rb_s+j]+
						  RO[1]*l_s[3*count_rb_s+j+rb_s[i]]+
						  RO[2]*l_s[3*count_rb_s+j+2*rb_s[i]];
			v_s[k+n_s]=v_r[i+n_r]+RO[3]*l_s[3*count_rb_s+j]+
								  RO[4]*l_s[3*count_rb_s+j+rb_s[i]]+
								  RO[5]*l_s[3*count_rb_s+j+2*rb_s[i]];
			v_s[k+2*n_s]=v_r[i+2*n_r]+RO[6]*l_s[3*count_rb_s+j]+
									  RO[7]*l_s[3*count_rb_s+j+rb_s[i]]+
									  RO[8]*l_s[3*count_rb_s+j+2*rb_s[i]];

			/* acceleration */
			a_s[k]=a_r[i]+ROO[0]*l_s[3*count_rb_s+j]+
						  ROO[1]*l_s[3*count_rb_s+j+rb_s[i]]+
						  ROO[2]*l_s[3*count_rb_s+j+2*rb_s[i]];
			a_s[k+n_s]=a_r[i+n_r]+ROO[3]*l_s[3*count_rb_s+j]+
								  ROO[4]*l_s[3*count_rb_s+j+rb_s[i]]+
								  ROO[5]*l_s[3*count_rb_s+j+2*rb_s[i]];
			a_s[k+2*n_s]=a_r[i+2*n_r]+ROO[6]*l_s[3*count_rb_s+j]+
									  ROO[7]*l_s[3*count_rb_s+j+rb_s[i]]+
									  ROO[8]*l_s[3*count_rb_s+j+2*rb_s[i]];

			/* updated Lagrangean 
			l_s[3*count_rb_s+j]=x_s[k]-x_r[i];
			l_s[3*count_rb_s+j+rb_s[i]]=x_s[k+n_s]-x_r[i+n_r];
			l_s[3*count_rb_s+j+2*rb_s[i]]=x_s[k+2*n_s]-x_r[i+2*n_r];
			*/
		}
		count_rb_s=count_rb_s+rb_s[i]; // update vector position
		
		for (j=0;j<rb_b[i];j++) // nodes
		{
			k=ind_r_0_b[j+count_rb_b]; // nodal index

			/* nodal position including COG
			-> COG can be unconstrained node
			-> x_b in local coordinate system */

			/* total Lagrangean 
			dx=x_b_0[k]-x_r_0[i];
			dy=x_b_0[k+n_b]-x_r_0[i+n_r];
			dz=x_b_0[k+2*n_b]-x_r_0[i+2*n_r];
			*/

			/* updated Lagrangean 
			dx=x_b[k]-x_r[i];
			dy=x_b[k+n_b]-x_r[i+n_r];
			dz=x_b[k+2*n_b]-x_r[i+2*n_r];
			*/
			
			/* position */
			x_b[k]=x_r[i]+R[0]*l_b[3*count_rb_b+j]+
						  R[1]*l_b[3*count_rb_b+j+rb_b[i]]+
						  R[2]*l_b[3*count_rb_b+j+2*rb_b[i]];
			x_b[k+n_b]=x_r[i+n_r]+R[3]*l_b[3*count_rb_b+j]+
						  		  R[4]*l_b[3*count_rb_b+j+rb_b[i]]+
								  R[5]*l_b[3*count_rb_b+j+2*rb_b[i]];
			x_b[k+2*n_b]=x_r[i+2*n_r]+R[6]*l_b[3*count_rb_b+j]+
									  R[7]*l_b[3*count_rb_b+j+rb_b[i]]+
									  R[8]*l_b[3*count_rb_b+j+2*rb_b[i]];

			/* velocity */
			v_b[k]=x_r[i]+RO[0]*l_b[3*count_rb_b+j]+
						  RO[1]*l_b[3*count_rb_b+j+rb_b[i]]+
						  RO[2]*l_b[3*count_rb_b+j+2*rb_b[i]];
			v_b[k+n_b]=x_r[i+n_r]+RO[3]*l_b[3*count_rb_b+j]+
								  RO[4]*l_b[3*count_rb_b+j+rb_b[i]]+
								  RO[5]*l_b[3*count_rb_b+j+2*rb_b[i]];
			v_b[k+2*n_b]=x_r[i+2*n_r]+RO[6]*l_b[3*count_rb_b+j]+
									  RO[7]*l_b[3*count_rb_b+j+rb_b[i]]+
									  RO[8]*l_b[3*count_rb_b+j+2*rb_b[i]];

			/* acceleration */
			a_b[k]=x_r[i]+ROO[0]*l_b[3*count_rb_b+j]+
						  ROO[1]*l_b[3*count_rb_b+j+rb_b[i]]+
						  ROO[2]*l_b[3*count_rb_b+j+2*rb_b[i]];
			a_b[k+n_b]=x_r[i+n_r]+ROO[3]*l_b[3*count_rb_b+j]+
								  ROO[4]*l_b[3*count_rb_b+j+rb_b[i]]+
								  ROO[5]*l_b[3*count_rb_b+j+2*rb_b[i]];
			a_b[k+2*n_b]=x_r[i+2*n_r]+ROO[6]*l_b[3*count_rb_b+j]+
									  ROO[7]*l_b[3*count_rb_b+j+rb_b[i]]+
									  ROO[8]*l_b[3*count_rb_b+j+2*rb_b[i]];

			/* updated Lagrangean 
			l_b[3*count_rb_b+j]=x_b[k]-x_r[i];
			l_b[3*count_rb_b+j+rb_b[i]]=x_b[k+n_b]-x_r[i+n_r];
			l_b[3*count_rb_b+j+2*rb_b[i]]=x_b[k+2*n_b]-x_r[i+2*n_r];
			*/

			/* R(psi).x = w.(w.x) + cos(psi).(w cross x) cross w + sin(psi).(w cross x) 
			   R(psi).x = w. scal + cos(psi).     cross1 cross w + sin(psi).     cross1
			   R(psi).x = w. scal + cos(psi).             cross2 + sin(psi).     cross1
			
			// x = (xx, xy, xz)
			double xx=l_b[3*count_rb_b+j];
			double xy=l_b[3*count_rb_b+j+n_r_0_b];
			double xz=l_b[3*count_rb_b+j+2*n_r_0_b];

			// scal = w.x
			double scal=urx*xx+ury*xy+urz*xz;

			// cross1 = w cross x
			double cross1x=ury*xz-urz*xy;
			double cross1y=urz*xx-urx*xz;
			double cross1z=urx*xy-ury*xx;

			// cross2 = cross1 cross w
			double cross2x=cross1y*urz-cross1z*ury;
			double cross2y=cross1z*urx-cross1x*urz;
			double cross2z=cross1x*ury-cross1y*urx;
			
			x_b[k]=x_r[i]+urx*scal+cos(psi)*cross2x+sin(psi)*cross1x;
			x_b[k+n_b]=x_r[i+n_r]+ury*scal+cos(psi)*cross2y+sin(psi)*cross1y;
			x_b[k+2*n_b]=x_r[i+2*n_r]+urz*scal+cos(psi)*cross2z+sin(psi)*cross1z;
			*/
		}
		count_rb_b=count_rb_b+rb_b[i]; // update vector position
	}
}
