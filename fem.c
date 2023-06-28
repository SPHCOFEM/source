#include "header.h"
#include "fem_bar.h"
#include "fem_beam.h"
#include "fem_triangle.h"
#include "fem_quad.h"
#include "fem_membrane.h"
#include "fem_shell.h"
#include "fem_tetrahedron.h"
#include "interpolation.h"

#undef __FUNC__
#define __FUNC__ "fem"

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

	double ax,double ay,double az)
{
	int i,j,k,h=0;
  
  	double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;

  	double *Me=0,*Ke=0;

  	int *ind=0;

  	double vavx,vavy,vavz;

	double a1,a2,a3,f1,f2,f3;

	/* average velocity for element damping */
	vavx=0.0;
	vavy=0.0;
	vavz=0.0;

	/* allocation - fem - set matrices to zero */
	memset(M,0,n_n*n_n*sizeof(double));
	memset(B,0,n_n*n_n*sizeof(double));
	memset(K,0,n_n*n_n*sizeof(double));
	memset(f,0,n_n*sizeof(double));

	for (i=0;i<n_e;i++)
    {
    	switch (dim)
		{
			case 1: /* 1d */
			{
				/* allocate element matrices */
				Me=createMemMore(double,1);
				Ke=createMemMore(double,1);

				Me[0]=rho_m[mat_e[i]];
				Ke[0]=T_m[mat_e[i]];

				/* element matrix size */
				h=1;

				/* allocate index vector */
				ind=createMemMore(int,h);
				ind[0]=nod_e[i];

				/* average velocity */
				vavx=v_b[nod_e[i]];

				break;
			}
			case 2: /* 2d */
			{	    
				switch (type_m[mat_e[i]])
				{
					case 4: /* no bending */
					{
						x1=x_b[nod_e[i]];
						y1=x_b[nod_e[i]+n_b];

						x2=x_b[nod_e[i+n_e]];
						y2=x_b[nod_e[i+n_e]+n_b];

						/*
						x10=x_b_0[nod_e[i]];
						y10=x_b_0[nod_e[i]+n_b];

						x20=x_b_0[nod_e[i+n_e]];
						y20=x_b_0[nod_e[i+n_e]+n_b];
						*/
				
						/* allocate element matrices */
						Me=createMemMore(double,16);
						Ke=createMemMore(double,16);

						/* small edeformation - original configuration
						eps_x=(norm(x2-x1,y2-y1,0.0)-norm(x20-x10,y20-y10,0.0))/norm(x2-x1,y2-y1,0.0);
						eps_y=-mu_m[mat_e[i]]*eps_x;
						A=sqr(sqrt(kappa_m[mat_e[i]]-mu_m[mat_e[i]]*eps_x));

						x1=x10;
						x2=x20;

						y1=y10;
						y2=y20;

						printf("l0 = %f, l = %f\n",norm(x20-x10,y20-y10,0.0),norm(x2-x1,y2-y1,0.0));
				
						printf("A0 = %f, A = %f\n",kappa_m[mat_e[i]],
							sqr(sqrt(kappa_m[mat_e[i]]-mu_m[mat_e[i]]*
							(norm(x2-x1,y2-y1,0.0)-norm(x20-x10,y20-y10,0.0))/norm(x2-x1,y2-y1,0.0))));
				
						fem_bar(Me,Ke,
							x1,x2,y1,y2,
							sqr(sqrt(kappa_m[mat_e[i]]-mu_m[mat_e[i]]*
							(norm(x2-x1,y2-y1,0.0)-norm(x20-x10,y20-y10,0.0))/norm(x2-x1,y2-y1,0.0))),
							rho_m[mat_e[i]]*norm(x20-x10,y20-y10,0.0)*kappa_m[mat_e[i]]/
							(norm(x2-x1,y2-y1,0.0)-norm(x20-x10,y20-y10,0.0)),
							T_m[mat_e[i]]);
						*/  

						fem_bar(Me,Ke,
							x1,x2,y1,y2,
							kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);
				
						/* element matrix size */
						h=4;
				
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;
				
						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]])/2.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b])/2.0;

						break;
					}
					case 5: /* bending active */
					{
						x1=x_b[nod_e[i]];
						y1=x_b[nod_e[i]+n_b];

						x2=x_b[nod_e[i+n_e]];
						y2=x_b[nod_e[i+n_e]+n_b];

						/* allocate element matrices */
						Me=createMemMore(double,16);
						Ke=createMemMore(double,16);

						fem_beam(Me,Ke,x1,x2,y1,y2,
							kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);

						/* element matrix size */
						h=4;
				
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;

						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]])/2.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b])/2.0;

						break;
					}
					case 6: /* full 2d */
					{
						if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
						{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
					
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,36);
							Ke=createMemMore(double,36);

							fem_triangle(Me,Ke,
								x1,x2,x3,y1,y2,y3,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=6;
				
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i]+n_b;
							ind[4]=nod_e[i+n_e]+n_b;
							ind[5]=nod_e[i+2*n_e]+n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]])/3.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b])/3.0;
						}
						else /* rectangle */
						{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
							x4=x_b[nod_e[i+3*n_e]];
					
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];
							y4=x_b[nod_e[i+3*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,64);
							Ke=createMemMore(double,64);
					
							fem_quad(Me,Ke,
								x1,x2,x3,x4,y1,y2,y3,y4,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=8;
				
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i+3*n_e];
							ind[4]=nod_e[i]+n_b;
							ind[5]=nod_e[i+n_e]+n_b;
							ind[6]=nod_e[i+2*n_e]+n_b;
							ind[7]=nod_e[i+3*n_e]+n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]]+v_b[nod_e[i+3*n_e]])/4.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b]+v_b[nod_e[i+3*n_e]+n_b])/4.0;
						}
						break;
					}
				}
				break;
			}
			case 3: /* 3d */
			{
				switch (type_m[mat_e[i]])
				{
					case 4: /* no bending */
					{
						if (nod_e[i+n_e]==nod_e[i+2*n_e]) /* bar */
						{
							x1=x_b[nod_e[i]];
							y1=x_b[nod_e[i]+n_b];

							x2=x_b[nod_e[i+n_e]];
							y2=x_b[nod_e[i+n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,16);
							Ke=createMemMore(double,16);

							fem_bar(Me,Ke,
								x1,x2,y1,y2,
								kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);
					
							/* element matrix size */
							h=4;
					
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i]+n_b;
							ind[3]=nod_e[i+n_e]+n_b;
					
							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]])/2.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b])/2.0;

							/* not implemented yet */
							/*
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
							.
							*/
						}
						else if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
						{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];

							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];		  

							z1=x_b[nod_e[i]+2*n_b];
							z2=x_b[nod_e[i+n_e]+2*n_b];
							z3=x_b[nod_e[i+2*n_e]+2*n_b];		  

							/* allocate element matrices */
							Me=createMemMore(double,81);
							Ke=createMemMore(double,81);

							fem_membrane(Me,Ke,
								x1,x2,x3,y1,y2,y3,z1,z2,z3,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=9;
				
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i]+n_b;
							ind[4]=nod_e[i+n_e]+n_b;
							ind[5]=nod_e[i+2*n_e]+n_b;
							ind[6]=nod_e[i]+2*n_b;
							ind[7]=nod_e[i+n_e]+2*n_b;
							ind[8]=nod_e[i+2*n_e]+2*n_b;
					
							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]])/3.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b])/3.0;
							vavz=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+2*n_b]+v_b[nod_e[i+2*n_e]+2*n_b])/3.0;

							/* not implemented yet */
							/*
							.
							.
							.
							.
							.
							. only as contact surface
							.
							.
							.
							.
							.
							*/
						}
						else /* rectangle */
						{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
							x4=x_b[nod_e[i+3*n_e]];
					
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];
							y4=x_b[nod_e[i+3*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,64);
							Ke=createMemMore(double,64);
					
							fem_quad(Me,Ke,
								x1,x2,x3,x4,y1,y2,y3,y4,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=8;
				
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i+3*n_e];
							ind[4]=nod_e[i]+n_b;
							ind[5]=nod_e[i+n_e]+n_b;
							ind[6]=nod_e[i+2*n_e]+n_b;
							ind[7]=nod_e[i+3*n_e]+n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]]+v_b[nod_e[i+3*n_e]])/4.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b]+v_b[nod_e[i+3*n_e]+n_b])/4.0;

							/* not implemented yet */
							/*
							.
							.
							.
							.
							.
							. only as contact surface (using triangle defined by nodes N1, N2, N3)
							.
							.
							.
							.
							.
							*/
						}
						break;
					}
					case 5: /* bending active */
					{
						if (nod_e[i+n_e]==nod_e[i+2*n_e]) /* bar */
						{
							x1=x_b[nod_e[i]];
							y1=x_b[nod_e[i]+n_b];

							x2=x_b[nod_e[i+n_e]];
							y2=x_b[nod_e[i+n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,16);
							Ke=createMemMore(double,16);

							fem_beam(Me,Ke,x1,x2,y1,y2,
								kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);

							/* element matrix size */
							h=4;
					
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i]+n_b;
							ind[3]=nod_e[i+n_e]+n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]])/2.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b])/2.0;

							/* not implemented yet */
							/*
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
							.
							*/
						}
						else if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
						{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];

							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];		  

							z1=x_b[nod_e[i]+2*n_b];
							z2=x_b[nod_e[i+n_e]+2*n_b];
							z3=x_b[nod_e[i+2*n_e]+2*n_b];		  

							/* allocate element matrices */
							Me=createMemMore(double,81);
							Ke=createMemMore(double,81);

							fem_shell(Me,Ke,
								x1,x2,x3,y1,y2,y3,z1,z2,z3,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=9;
							
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i]+n_b;
							ind[4]=nod_e[i+n_e]+n_b;
							ind[5]=nod_e[i+2*n_e]+n_b;
							ind[6]=nod_e[i]+2*n_b;
							ind[7]=nod_e[i+n_e]+2*n_b;
							ind[8]=nod_e[i+2*n_e]+2*n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]])/3.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b])/3.0;
							vavz=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+2*n_b]+v_b[nod_e[i+2*n_e]+2*n_b])/3.0;

							/* not implemented yet */
							/*
							.
							.
							.
							.
							.
							. only as contact surface
							.
							.
							.
							.
							.
							*/
						}
						else /* rectangle */
						{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
							x4=x_b[nod_e[i+3*n_e]];
					
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];
							y4=x_b[nod_e[i+3*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,64);
							Ke=createMemMore(double,64);
					
							fem_quad(Me,Ke,
								x1,x2,x3,x4,y1,y2,y3,y4,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=8;
				
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i+3*n_e];
							ind[4]=nod_e[i]+n_b;
							ind[5]=nod_e[i+n_e]+n_b;
							ind[6]=nod_e[i+2*n_e]+n_b;
							ind[7]=nod_e[i+3*n_e]+n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]]+v_b[nod_e[i+3*n_e]])/4.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b]+v_b[nod_e[i+3*n_e]+n_b])/4.0;

							/* not implemented yet */
							/*
							.
							.
							.
							.
							.
							. only as contact surface (using triangle defined by nodes N1, N2, N3)
							.
							.
							.
							.
							.
							*/
						}
						break;
					}
					case 6: /* tetrahedron */
					{
						x1=x_b[nod_e[i]];
						x2=x_b[nod_e[i+n_e]];
						x3=x_b[nod_e[i+2*n_e]];
						x4=x_b[nod_e[i+3*n_e]];

						y1=x_b[nod_e[i]+n_b];
						y2=x_b[nod_e[i+n_e]+n_b];
						y3=x_b[nod_e[i+2*n_e]+n_b];		  
						y4=x_b[nod_e[i+3*n_e]+n_b];		  

						z1=x_b[nod_e[i]+2*n_b];
						z2=x_b[nod_e[i+n_e]+2*n_b];
						z3=x_b[nod_e[i+2*n_e]+2*n_b];		  
						z4=x_b[nod_e[i+3*n_e]+2*n_b];		  

						/* allocate element matrices */
						Me=createMemMore(double,144);
						Ke=createMemMore(double,144);

						fem_tetrahedron(Me,Ke,
							x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,
							rho_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

						/* elemenent matrix size */
						h=12;
							
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i+3*n_e];
						ind[4]=nod_e[i]+n_b;
						ind[5]=nod_e[i+n_e]+n_b;
						ind[6]=nod_e[i+2*n_e]+n_b;
						ind[7]=nod_e[i+3*n_e]+n_b;
						ind[8]=nod_e[i]+2*n_b;
						ind[9]=nod_e[i+n_e]+2*n_b;
						ind[10]=nod_e[i+2*n_e]+2*n_b;
						ind[11]=nod_e[i+3*n_e]+2*n_b;

						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]]+v_b[nod_e[i+3*n_e]])/4.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b]+v_b[nod_e[i+3*n_e]+n_b])/4.0;
						vavz=(v_b[nod_e[i]+2*n_b]+v_b[nod_e[i+n_e]+2*n_b]+v_b[nod_e[i+2*n_e]+2*n_b]+v_b[nod_e[i+3*n_e]+2*n_b])/4.0;

						break;
					}
		      	}
	    		break;
	  		}
		}
      
    	for (j=0;j<h;j++)
		{
			for (k=0;k<h;k++)
	    	{
	    		M[ind[j]+n_n*ind[k]]=M[ind[j]+n_n*ind[k]]+Me[j+h*k];
	      		K[ind[j]+n_n*ind[k]]=K[ind[j]+n_n*ind[k]]+Ke[j+h*k];
	    	}
	  		B[ind[j]+n_n*ind[j]]=B[ind[j]+n_n*ind[j]]+gamma_m[mat_e[i]]*((j<n_n/dim)*vavx+
				((j>=n_n/dim)&&(j<2*n_n/dim))*vavy+
				((j>=2*n_n/dim)&&(j<3*n_n/dim))*vavz-
				v_b[ind[j]]);
	  		f[ind[j]]=f[ind[j]]+f_b[ind[j]];
		}

	    /* free index vector */
    	freeMem(ind);

    	/* free element matrices */
    	freeMem(Me);
    	freeMem(Ke);
	}

  	/* acceleration field */
  	for (i=0;i<n_a_0_b;i++)
    {
    	for (j=0;j<n_a_0_b;j++)
		{
			if (type_a_0_b[j]==0) // constant
			{
				a1=a_0_b[j];
				a2=a_0_b[j+n_a_0_b];
				a3=a_0_b[j+2*n_a_0_b];
			}
			else // function
			{
				a1=interpolation(t-dt,(int)a_0_b[j],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
				a2=interpolation(t-dt,(int)a_0_b[j+n_a_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
				a3=interpolation(t-dt,(int)a_0_b[j+2*n_a_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			}
			f[ind_a_0_b[i]]=f[ind_a_0_b[i]]+M[ind_a_0_b[i]+ind_a_0_b[j]*n_n]*((ind_a_0_b[j]<n_n/dim)*a1+
				((ind_a_0_b[j]>=n_n/dim)&&(ind_a_0_b[j]<2*n_n/dim))*a2+
				((ind_a_0_b[j]>=2*n_n/dim)&&(ind_a_0_b[j]<3*n_n/dim))*a3);
		}
    }

  	/* concentrated forces on boundary */
  	for (i=0;i<n_f_0_b;i++)
    {
		if (type_f_0_b[i]==0) // constant
		{
			f1=f_0_b[i];
			f2=f_0_b[i+n_f_0_b];
			f3=f_0_b[i+2*n_f_0_b];
		}
		else // function
		{
			f1=interpolation(t-dt,(int)f_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f2=interpolation(t-dt,(int)f_0_b[i+n_f_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f3=interpolation(t-dt,(int)f_0_b[i+2*n_f_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
		}
		f[ind_a_0_b[i]]=f[ind_a_0_b[i]]+(ind_a_0_b[i]<n_n/dim)*f1+
			((ind_a_0_b[i]>=n_n/dim)&&(ind_a_0_b[i]<2*n_n/dim))*f2+
			((ind_a_0_b[i]>=2*n_n/dim)&&(ind_a_0_b[i]<3*n_n/dim))*f3;
    }

  	/* material damping and acceleration field */  
  	for (i=0;i<n_n;i++)
  	{
		for (j=0;j<n_n;j++)
		{
	  		B[i+j*n_n]=B[i+j*n_n]+c0*M[i+j*n_n]+c1*K[i+j*n_n];

	  		f[i]=f[i]+M[i+j*n_n]*((j<n_n/dim)*ax+
				((j>=n_n/dim)&&(j<2*n_n/dim))*ay+
				((j>=2*n_n/dim)&&(j<3*n_n/dim))*az);
		}
	}

  	/* nodal damping */
  	for (i=0;i<n_d_0_b;i++)
    {
    	B[ind_d_0_b[i]+ind_d_0_b[i]*n_n]=B[ind_d_0_b[i]+ind_d_0_b[i]*n_n]+d_0_b[i];
      	if (dim>1) B[ind_d_0_b[i]+n_b+(ind_d_0_b[i]+n_b)*n_n]=B[ind_d_0_b[i]+n_b+(ind_d_0_b[i]+n_b)*n_n]+d_0_b[i+n_d_0_b];
      	if (dim>2) B[ind_d_0_b[i]+2*n_b+(ind_d_0_b[i]+2*n_b)*n_n]=B[ind_d_0_b[i]+2*n_b+(ind_d_0_b[i]+2*n_b)*n_n]+d_0_b[i+2*n_d_0_b];
    }
}

#undef __FUNC__
#define __FUNC__ "mass_stiffness_damping"

void mass_stiffness_damping(int *type_m,

    double *rho_m,double *T_m,double *mu_m,double *kappa_m,

    int n_b,
    double *x_b,double *x_b_0,double *v_b,

    int n_e,int *nod_e,int *mat_e,

    int n_n,
    double *M,double *B,double *K,

    double c0, double c1,

    int dim)
{
	int i,j=0,k,h=0;
  
	double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;

	double *Me=0,*Ke=0;

	int *ind=0;

	/* allocation - fem - set matrices to zero */
	memset(M,0,n_n*n_n*sizeof(double));
	memset(B,0,n_n*n_n*sizeof(double));
	memset(K,0,n_n*n_n*sizeof(double));

	for (i=0;i<n_e;i++)
    {
		switch (dim)
		{
			case 1: /* 1d */
	  		{
				/* allocate element matrices */
				Me=createMemMore(double,1);
				Ke=createMemMore(double,1);

				Me[0]=rho_m[mat_e[i]];
				Ke[0]=T_m[mat_e[i]];

				/* element matrix size */
				h=1;

				/* allocate index vector */
				ind=createMemMore(int,h);
				ind[0]=nod_e[i];

				break;
			}
			case 2: /* 2d */
	  		{	    
	    		switch (type_m[mat_e[i]])
				{
					case 4: /* no bending */
					{
						x1=x_b[nod_e[i]];
						y1=x_b[nod_e[i]+n_b];

						x2=x_b[nod_e[i+n_e]];
						y2=x_b[nod_e[i+n_e]+n_b];

						/* allocate element matrices */
						Me=createMemMore(double,16);
						Ke=createMemMore(double,16);

						fem_bar(Me,Ke,
					  		x1,x2,y1,y2,
			  				kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);
		  
						/* element matrix size */
						h=4;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;
						
						break;
					}
	      			case 5: /* bending active */
					{
						x1=x_b[nod_e[i]];
						y1=x_b[nod_e[i]+n_b];

						x2=x_b[nod_e[i+n_e]];
						y2=x_b[nod_e[i+n_e]+n_b];

						/* allocate element matrices */
						Me=createMemMore(double,16);
						Ke=createMemMore(double,16);

						fem_beam(Me,Ke,x1,x2,y1,y2,
							kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);

						/* element matrix size */
						h=4;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;
						break;
					}
	    		  	case 6: /* full 2d */
					{
		  				if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
		    			{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
							
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,36);
							Ke=createMemMore(double,36);

							fem_triangle(Me,Ke,
								x1,x2,x3,y1,y2,y3,
				   				rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=6;
						
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i]+n_b;
							ind[4]=nod_e[i+n_e]+n_b;
							ind[5]=nod_e[i+2*n_e]+n_b;
						}
		  				else /* rectangle */
		    			{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
							x4=x_b[nod_e[i+3*n_e]];
							
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];
							y4=x_b[nod_e[i+3*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,64);
							Ke=createMemMore(double,64);
							
							fem_quad(Me,Ke,
								x1,x2,x3,x4,y1,y2,y3,y4,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=8;
						
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i+3*n_e];
							ind[4]=nod_e[i]+n_b;
							ind[5]=nod_e[i+n_e]+n_b;
							ind[6]=nod_e[i+2*n_e]+n_b;
							ind[7]=nod_e[i+3*n_e]+n_b;
						}
		  				break;
					}
	      		}
	    		break;
	  		}
			case 3: /* 3d */
	  		{
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 4: /* no bending */
					{
						x1=x_b[nod_e[i]];
						x2=x_b[nod_e[i+n_e]];
						x3=x_b[nod_e[i+2*n_e]];

						y1=x_b[nod_e[i]+n_b];
						y2=x_b[nod_e[i+n_e]+n_b];
						y3=x_b[nod_e[i+2*n_e]+n_b];		  

						z1=x_b[nod_e[i]+2*n_b];
						z2=x_b[nod_e[i+n_e]+2*n_b];
						z3=x_b[nod_e[i+2*n_e]+2*n_b];		  

						/* allocate element matrices */
						Me=createMemMore(double,81);
						Ke=createMemMore(double,81);

						fem_membrane(Me,Ke,
							x1,x2,x3,y1,y2,y3,z1,z2,z3,
							rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

						/* element matrix size */
						h=9;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;

			  			break;
					}
	      			case 5: /* bending active */
					{
						x1=x_b[nod_e[i]];
						x2=x_b[nod_e[i+n_e]];
						x3=x_b[nod_e[i+2*n_e]];

						y1=x_b[nod_e[i]+n_b];
						y2=x_b[nod_e[i+n_e]+n_b];
						y3=x_b[nod_e[i+2*n_e]+n_b];		  

						z1=x_b[nod_e[i]+2*n_b];
						z2=x_b[nod_e[i+n_e]+2*n_b];
						z3=x_b[nod_e[i+2*n_e]+2*n_b];		  

						/* allocate element matrices */
						Me=createMemMore(double,81);
						Ke=createMemMore(double,81);

						fem_shell(Me,Ke,
			    			x1,x2,x3,y1,y2,y3,z1,z2,z3,
			    			rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

						/* element matrix size */
						h=9;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;

		  				break;
					}
	      			case 6: /* tetrahedron */
					{
						x1=x_b[nod_e[i]];
						x2=x_b[nod_e[i+n_e]];
						x3=x_b[nod_e[i+2*n_e]];
						x4=x_b[nod_e[i+3*n_e]];

						y1=x_b[nod_e[i]+n_b];
						y2=x_b[nod_e[i+n_e]+n_b];
						y3=x_b[nod_e[i+2*n_e]+n_b];		  
						y4=x_b[nod_e[i+3*n_e]+n_b];		  

						z1=x_b[nod_e[i]+2*n_b];
						z2=x_b[nod_e[i+n_e]+2*n_b];
						z3=x_b[nod_e[i+2*n_e]+2*n_b];		  
						z4=x_b[nod_e[i+3*n_e]+2*n_b];		  

						/* allocate element matrices */
						Me=createMemMore(double,144);
						Ke=createMemMore(double,144);

						fem_tetrahedron(Me,Ke,
							x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,
						  	rho_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

						/* elemenent matrix size */
						h=12;
							
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i+3*n_e];
						ind[4]=nod_e[i]+n_b;
						ind[5]=nod_e[i+n_e]+n_b;
						ind[6]=nod_e[i+2*n_e]+n_b;
						ind[7]=nod_e[i+3*n_e]+n_b;
						ind[8]=nod_e[i]+2*n_b;
						ind[9]=nod_e[i+n_e]+2*n_b;
						ind[10]=nod_e[i+2*n_e]+2*n_b;
						ind[11]=nod_e[i+3*n_e]+2*n_b;

			  			break;
					}
	      		}
	    		break;
	  		}
		}
      
      	for (j=0;j<h;j++)
		{
	  		for (k=0;k<h;k++)
	    	{
				M[ind[j]+n_n*ind[k]]=M[ind[j]+n_n*ind[k]]+Me[j+h*k];
				K[ind[j]+n_n*ind[k]]=K[ind[j]+n_n*ind[k]]+Ke[j+h*k];
	    	}
		}

      	/* free index vector */
      	freeMem(ind);

      	/* free element matrices */
      	freeMem(Me);
      	freeMem(Ke);
    }

  	for (i=0;i<n_n;i++)
    {
    	B[i+n_n*j]=c0*M[i+n_n*j]+c1*K[i+n_n*j];
    }
}

#undef __FUNC__
#define __FUNC__ "mass_stiffness"

void mass_stiffness(int *type_m,

    double *rho_m,double *T_m,double *mu_m,double *kappa_m,

    int n_b,
    double *x_b,double *x_b_0,double *v_b,

    int n_e,int *nod_e,int *mat_e,

    int n_n,
    double *M,double *K,

    int dim)
{
	int i,j,k,h=0;
  
  	double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;

  	double *Me=0,*Ke=0;

  	int *ind=0;

  	/* allocation - fem - set matrices to zero */
  	memset(M,0,n_n*n_n*sizeof(double));
  	memset(K,0,n_n*n_n*sizeof(double));

  	for (i=0;i<n_e;i++)
    {
      	switch (dim)
		{
			case 1: /* 1d */
			{
				/* allocate element matrices */
				Me=createMemMore(double,1);
				Ke=createMemMore(double,1);

				Me[0]=rho_m[mat_e[i]];
				Ke[0]=T_m[mat_e[i]];

				/* element matrix size */
				h=1;

				/* allocate index vector */
				ind=createMemMore(int,h);
				ind[0]=nod_e[i];

				break;
			}
			case 2: /* 2d */
		  	{	    
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 4: /* no bending */
					{
						x1=x_b[nod_e[i]];
						y1=x_b[nod_e[i]+n_b];

						x2=x_b[nod_e[i+n_e]];
						y2=x_b[nod_e[i+n_e]+n_b];

						/* allocate element matrices */
						Me=createMemMore(double,16);
						Ke=createMemMore(double,16);

						fem_bar(Me,Ke,
							x1,x2,y1,y2,
							kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);
						
						/* element matrix size */
						h=4;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;
						
		  				break;
					}
	      			case 5: /* bending active */
					{
						x1=x_b[nod_e[i]];
						y1=x_b[nod_e[i]+n_b];

						x2=x_b[nod_e[i+n_e]];
						y2=x_b[nod_e[i+n_e]+n_b];

						/* allocate element matrices */
						Me=createMemMore(double,16);
						Ke=createMemMore(double,16);

						fem_beam(Me,Ke,x1,x2,y1,y2,
							kappa_m[mat_e[i]],rho_m[mat_e[i]],T_m[mat_e[i]]);

						/* element matrix size */
						h=4;
						
		  				/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;
						break;
					}
	      			case 6: /* full 2d */
					{
		  				if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
		    			{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
							
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,36);
							Ke=createMemMore(double,36);

							fem_triangle(Me,Ke,
								x1,x2,x3,y1,y2,y3,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=6;
						
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i]+n_b;
							ind[4]=nod_e[i+n_e]+n_b;
							ind[5]=nod_e[i+2*n_e]+n_b;
						}
		  				else /* rectangle */
		    			{
							x1=x_b[nod_e[i]];
							x2=x_b[nod_e[i+n_e]];
							x3=x_b[nod_e[i+2*n_e]];
							x4=x_b[nod_e[i+3*n_e]];
							
							y1=x_b[nod_e[i]+n_b];
							y2=x_b[nod_e[i+n_e]+n_b];
							y3=x_b[nod_e[i+2*n_e]+n_b];
							y4=x_b[nod_e[i+3*n_e]+n_b];

							/* allocate element matrices */
							Me=createMemMore(double,64);
							Ke=createMemMore(double,64);
							
							fem_quad(Me,Ke,
								x1,x2,x3,x4,y1,y2,y3,y4,
								rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

							/* element matrix size */
							h=8;
						
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i+3*n_e];
							ind[4]=nod_e[i]+n_b;
							ind[5]=nod_e[i+n_e]+n_b;
							ind[6]=nod_e[i+2*n_e]+n_b;
							ind[7]=nod_e[i+3*n_e]+n_b;
		    			}
		  				break;
					}
	      		}
	    		break;
	  		}
			case 3: /* 3d */
	  		{
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 4: /* no bending */
					{
						x1=x_b[nod_e[i]];
						x2=x_b[nod_e[i+n_e]];
						x3=x_b[nod_e[i+2*n_e]];

						y1=x_b[nod_e[i]+n_b];
						y2=x_b[nod_e[i+n_e]+n_b];
						y3=x_b[nod_e[i+2*n_e]+n_b];		  

						z1=x_b[nod_e[i]+2*n_b];
						z2=x_b[nod_e[i+n_e]+2*n_b];
						z3=x_b[nod_e[i+2*n_e]+2*n_b];		  

						/* allocate element matrices */
						Me=createMemMore(double,81);
						Ke=createMemMore(double,81);

            	      	fem_membrane(Me,Ke,
			       			x1,x2,x3,y1,y2,y3,z1,z2,z3,
			       			rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

						/* element matrix size */
						h=9;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;

		  				break;
					}
	      			case 5: /* bending active */
					{
						x1=x_b[nod_e[i]];
						x2=x_b[nod_e[i+n_e]];
						x3=x_b[nod_e[i+2*n_e]];

						y1=x_b[nod_e[i]+n_b];
						y2=x_b[nod_e[i+n_e]+n_b];
						y3=x_b[nod_e[i+2*n_e]+n_b];		  

						z1=x_b[nod_e[i]+2*n_b];
						z2=x_b[nod_e[i+n_e]+2*n_b];
						z3=x_b[nod_e[i+2*n_e]+2*n_b];		  

						/* allocate element matrices */
						Me=createMemMore(double,81);
						Ke=createMemMore(double,81);

						fem_shell(Me,Ke,
							x1,x2,x3,y1,y2,y3,z1,z2,z3,
							rho_m[mat_e[i]],kappa_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

						/* element matrix size */
		  				h=9;
		  
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;

						break;
					}
	      			case 6: /* tetrahedron */
					{
						x1=x_b[nod_e[i]];
						x2=x_b[nod_e[i+n_e]];
						x3=x_b[nod_e[i+2*n_e]];
						x4=x_b[nod_e[i+3*n_e]];

						y1=x_b[nod_e[i]+n_b];
						y2=x_b[nod_e[i+n_e]+n_b];
						y3=x_b[nod_e[i+2*n_e]+n_b];		  
						y4=x_b[nod_e[i+3*n_e]+n_b];		  

						z1=x_b[nod_e[i]+2*n_b];
						z2=x_b[nod_e[i+n_e]+2*n_b];
						z3=x_b[nod_e[i+2*n_e]+2*n_b];		  
						z4=x_b[nod_e[i+3*n_e]+2*n_b];		  

						/* allocate element matrices */
						Me=createMemMore(double,144);
						Ke=createMemMore(double,144);

                  		fem_tetrahedron(Me,Ke,
				  			x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,
				  			rho_m[mat_e[i]],T_m[mat_e[i]],mu_m[mat_e[i]]);

						/* elemenent matrix size */
						h=12;
							
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i+3*n_e];
						ind[4]=nod_e[i]+n_b;
						ind[5]=nod_e[i+n_e]+n_b;
						ind[6]=nod_e[i+2*n_e]+n_b;
						ind[7]=nod_e[i+3*n_e]+n_b;
						ind[8]=nod_e[i]+2*n_b;
						ind[9]=nod_e[i+n_e]+2*n_b;
						ind[10]=nod_e[i+2*n_e]+2*n_b;
						ind[11]=nod_e[i+3*n_e]+2*n_b;

		  				break;
					}
	      		}
	    		break;
	  		}
		}
      
      	for (j=0;j<h;j++)
		{
	  		for (k=0;k<h;k++)
	    	{
	      		M[ind[j]+n_n*ind[k]]=M[ind[j]+n_n*ind[k]]+Me[j+h*k];
	      		K[ind[j]+n_n*ind[k]]=K[ind[j]+n_n*ind[k]]+Ke[j+h*k];
	    	}
		}

		/* free index vector */
		freeMem(ind);

		/* free element matrices */
		freeMem(Me);
		freeMem(Ke);
    }
}

#undef __FUNC__
#define __FUNC__ "damping_rhs"

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

	double ax,double ay,double az)
{
	int i,j,h=0;
	
	int *ind=0;

	double vavx,vavy,vavz;

	double a1,a2,a3,f1,f2,f3;

	/* average velocity for element damping */
	vavx=0.0;
	vavy=0.0;
	vavz=0.0;

	/* allocation - fem - set matrices to zero */
	memset(B,0,n_n*n_n*sizeof(double));
	memset(f,0,n_n*sizeof(double));

	for (i=0;i<n_e;i++)
    {
      	switch (dim)
		{
			case 1: /* 1d */
	  		{
				/* allocate index vector */
				ind=createMemMore(int,h);
				ind[0]=nod_e[i];

				/* average velocity */
				vavx=v_b[nod_e[i]];

				break;
	  		}
			case 2: /* 2d */
	  		{	    
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 4: /* no bending */
					{
						/* element matrix size */
						h=4;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;
						
						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]])/2.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b])/2.0;

						break;
					}
	      			case 5: /* bending active */
					{
						/* element matrix size */
						h=4;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;

						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]])/2.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b])/2.0;

						break;
					}
	      			case 6: /* full 2d */
					{
		  				if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
		    			{
							/* element matrix size */
							h=6;
						
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i]+n_b;
							ind[4]=nod_e[i+n_e]+n_b;
							ind[5]=nod_e[i+2*n_e]+n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]])/3.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b])/3.0;
		    			}
		  				else /* rectangle */
		    			{
		      				/* element matrix size */
		      				h=8;
		  
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i+3*n_e];
							ind[4]=nod_e[i]+n_b;
							ind[5]=nod_e[i+n_e]+n_b;
							ind[6]=nod_e[i+2*n_e]+n_b;
							ind[7]=nod_e[i+3*n_e]+n_b;

							/* average velocity */
							vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]]+v_b[nod_e[i+3*n_e]])/4.0;
							vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b]+v_b[nod_e[i+3*n_e]+n_b])/4.0;
		    			}
		  				break;
					}
	      		}
	    		break;
	  		}
			case 3: /* 3d */
	  		{
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 4: /* no bending */
					{
						/* element matrix size */
						h=9;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;
							
						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]])/3.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b])/3.0;
						vavz=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+2*n_b]+v_b[nod_e[i+2*n_e]+2*n_b])/3.0;

				  		break;
					}
	      			case 5: /* bending active */
					{
		  				/* element matrix size */
						h=9;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;

						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]])/3.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b])/3.0;
						vavz=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+2*n_b]+v_b[nod_e[i+2*n_e]+2*n_b])/3.0;

						break;
					}
	      			case 6: /* tetrahedron */
					{
						/* elemenent matrix size */
						h=12;
							
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i+3*n_e];
						ind[4]=nod_e[i]+n_b;
						ind[5]=nod_e[i+n_e]+n_b;
						ind[6]=nod_e[i+2*n_e]+n_b;
						ind[7]=nod_e[i+3*n_e]+n_b;
						ind[8]=nod_e[i]+2*n_b;
						ind[9]=nod_e[i+n_e]+2*n_b;
						ind[10]=nod_e[i+2*n_e]+2*n_b;
						ind[11]=nod_e[i+3*n_e]+2*n_b;

						/* average velocity */
						vavx=(v_b[nod_e[i]]+v_b[nod_e[i+n_e]]+v_b[nod_e[i+2*n_e]]+v_b[nod_e[i+3*n_e]])/4.0;
						vavy=(v_b[nod_e[i]+n_b]+v_b[nod_e[i+n_e]+n_b]+v_b[nod_e[i+2*n_e]+n_b]+v_b[nod_e[i+3*n_e]+n_b])/4.0;
						vavz=(v_b[nod_e[i]+2*n_b]+v_b[nod_e[i+n_e]+2*n_b]+v_b[nod_e[i+2*n_e]+2*n_b]+v_b[nod_e[i+3*n_e]+2*n_b])/4.0;

						break;
					}
	      		}
	    		break;
	  		}
		}
      
      	for (j=0;j<h;j++)
		{
	  		B[ind[j]+n_n*ind[j]]=B[ind[j]+n_n*ind[j]]+gamma_m[mat_e[i]]*((j<n_n/dim)*vavx+
			    ((j>=n_n/dim)&&(j<2*n_n/dim))*vavy+
			    ((j>=2*n_n/dim)&&(j<3*n_n/dim))*vavz-
			    v_b[ind[j]]);
	  		f[ind[j]]=f[ind[j]]+f_b[ind[j]];
		}

      	/* free index vector */
      	freeMem(ind);
    }

	/* acceleration field */
	for (i=0;i<n_a_0_b;i++)
    {
    	for (j=0;j<n_a_0_b;j++)
		{
			if (type_a_0_b[j]==0) // constant
			{
				a1=a_0_b[j];
				a2=a_0_b[j+n_a_0_b];
				a3=a_0_b[j+2*n_a_0_b];
			}
			else // function
			{
				a1=interpolation(t-dt,(int)a_0_b[j],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
				a2=interpolation(t-dt,(int)a_0_b[j+n_a_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
				a3=interpolation(t-dt,(int)a_0_b[j+2*n_a_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			}
			f[ind_a_0_b[i]]=f[ind_a_0_b[i]]+M[ind_a_0_b[i]+ind_a_0_b[j]*n_n]*((ind_a_0_b[j]<n_n/dim)*a1+
				((ind_a_0_b[j]>=n_n/dim)&&(ind_a_0_b[j]<2*n_n/dim))*a2+
				((ind_a_0_b[j]>=2*n_n/dim)&&(ind_a_0_b[j]<3*n_n/dim))*a3);
		}
    }

  	/* concentrated forces on boundary */
  	for (i=0;i<n_f_0_b;i++)
    {
		if (type_f_0_b[i]==0) // constant
		{
			f1=f_0_b[i];
			f2=f_0_b[i+n_f_0_b];
			f3=f_0_b[i+2*n_f_0_b];
		}
		else // function
		{
			f1=interpolation(t-dt,(int)f_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f2=interpolation(t-dt,(int)f_0_b[i+n_f_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f3=interpolation(t-dt,(int)f_0_b[i+2*n_f_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
		}
		f[ind_a_0_b[i]]=f[ind_a_0_b[i]]+(ind_a_0_b[i]<n_n/dim)*f1+
			((ind_a_0_b[i]>=n_n/dim)&&(ind_a_0_b[i]<2*n_n/dim))*f2+
			((ind_a_0_b[i]>=2*n_n/dim)&&(ind_a_0_b[i]<3*n_n/dim))*f3;
    }

	/* material damping and acceleration field */  
	for (i=0;i<n_n;i++)
    {
      	for (j=0;j<n_n;j++)
		{
	  		B[i+j*n_n]=B[i+j*n_n]+c0*M[i+j*n_n]+c1*K[i+j*n_n];

	  		f[i]=f[i]+M[i+j*n_n]*((j<n_n/dim)*ax+
				((j>=n_n/dim)&&(j<2*n_n/dim))*ay+
				((j>=2*n_n/dim)&&(j<3*n_n/dim))*az);
		}
    }

  	/* nodal damping */
  	for (i=0;i<n_d_0_b;i++)
    {
    	B[ind_d_0_b[i]+ind_d_0_b[i]*n_n]=B[ind_d_0_b[i]+ind_d_0_b[i]*n_n]+d_0_b[i];
      	if (dim>1) B[ind_d_0_b[i]+n_b+(ind_d_0_b[i]+n_b)*n_n]=B[ind_d_0_b[i]+n_b+(ind_d_0_b[i]+n_b)*n_n]+d_0_b[i+n_d_0_b];
      	if (dim>2) B[ind_d_0_b[i]+2*n_b+(ind_d_0_b[i]+2*n_b)*n_n]=B[ind_d_0_b[i]+2*n_b+(ind_d_0_b[i]+2*n_b)*n_n]+d_0_b[i+2*n_d_0_b];
    }
}

#undef __FUNC__
#define __FUNC__ "rhs"

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

	double ax,double ay,double az)
{
	int i,j,h=0;
	
	int *ind=0;

	double a1,a2,a3,f1,f2,f3;

	/* allocation - fem - set matrices to zero */
	memset(f,0,n_n*sizeof(double));

	for (i=0;i<n_e;i++)
    {
      	switch (dim)
		{
			case 1: /* 1d */
	  		{
				/* allocate index vector */
				ind=createMemMore(int,h);
				ind[0]=nod_e[i];

				break;
	  		}
			case 2: /* 2d */
	  		{	    
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 4: /* no bending */
					{
						/* element matrix size */
						h=4;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;

						break;
					}
	      			case 5: /* bending active */
					{
		  				/* element matrix size */
						h=4;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i]+n_b;
						ind[3]=nod_e[i+n_e]+n_b;

						break;
					}
	      			case 6: /* full 2d */
					{
		  				if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
		    			{
							/* element matrix size */
							h=6;
						
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i]+n_b;
							ind[4]=nod_e[i+n_e]+n_b;
							ind[5]=nod_e[i+2*n_e]+n_b;
		    			}
		  				else /* rectangle */
		    			{
							/* element matrix size */
							h=8;
						
							/* allocate index vector */
							ind=createMemMore(int,h);

							ind[0]=nod_e[i];
							ind[1]=nod_e[i+n_e];
							ind[2]=nod_e[i+2*n_e];
							ind[3]=nod_e[i+3*n_e];
							ind[4]=nod_e[i]+n_b;
							ind[5]=nod_e[i+n_e]+n_b;
							ind[6]=nod_e[i+2*n_e]+n_b;
							ind[7]=nod_e[i+3*n_e]+n_b;
		    			}
		  				break;
					}
	      		}
	    		break;
	  		}
			case 3: /* 3d */
	  		{
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 4: /* no bending */
					{
						/* element matrix size */
						h=9;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;

						break;
					}
	      			case 5: /* bending active */
					{
						/* element matrix size */
						h=9;
						
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i]+n_b;
						ind[4]=nod_e[i+n_e]+n_b;
						ind[5]=nod_e[i+2*n_e]+n_b;
						ind[6]=nod_e[i]+2*n_b;
						ind[7]=nod_e[i+n_e]+2*n_b;
						ind[8]=nod_e[i+2*n_e]+2*n_b;

						break;
					}
	      			case 6: /* tetrahedron */
					{
						/* elemenent matrix size */
						h=12;
							
						/* allocate index vector */
						ind=createMemMore(int,h);

						ind[0]=nod_e[i];
						ind[1]=nod_e[i+n_e];
						ind[2]=nod_e[i+2*n_e];
						ind[3]=nod_e[i+3*n_e];
						ind[4]=nod_e[i]+n_b;
						ind[5]=nod_e[i+n_e]+n_b;
						ind[6]=nod_e[i+2*n_e]+n_b;
						ind[7]=nod_e[i+3*n_e]+n_b;
						ind[8]=nod_e[i]+2*n_b;
						ind[9]=nod_e[i+n_e]+2*n_b;
						ind[10]=nod_e[i+2*n_e]+2*n_b;
						ind[11]=nod_e[i+3*n_e]+2*n_b;

						break;
					}
	      		}
	    		break;
	  		}
		}
      
      	for (j=0;j<h;j++)
		{
	  		f[ind[j]]=f[ind[j]]+f_b[ind[j]];
		}

		/* free index vector */
		freeMem(ind);
    }

  	/* concentrated forces on boundary */
  	for (i=0;i<n_f_0_b;i++)
    {
		if (type_f_0_b[i]==0) // constant
		{
			f1=f_0_b[i];
			f2=f_0_b[i+n_f_0_b];
			f3=f_0_b[i+2*n_f_0_b];
		}
		else // function
		{
			f1=interpolation(t-dt,(int)f_0_b[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f2=interpolation(t-dt,(int)f_0_b[i+n_f_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f3=interpolation(t-dt,(int)f_0_b[i+2*n_f_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
		}
		f[ind_a_0_b[i]]=f[ind_a_0_b[i]]+(ind_a_0_b[i]<n_n/dim)*f1+
			((ind_a_0_b[i]>=n_n/dim)&&(ind_a_0_b[i]<2*n_n/dim))*f2+
			((ind_a_0_b[i]>=2*n_n/dim)&&(ind_a_0_b[i]<3*n_n/dim))*f3;
    }

	/* acceleration field */
	for (i=0;i<n_a_0_b;i++)
    {
      	for (j=0;j<n_a_0_b;j++)
		{
			if (type_a_0_b[j]==0) // constant
			{
				a1=a_0_b[j];
				a2=a_0_b[j+n_a_0_b];
				a3=a_0_b[j+2*n_a_0_b];
			}
			else // function
			{
				a1=interpolation(t-dt,(int)a_0_b[j],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
				a2=interpolation(t-dt,(int)a_0_b[j+n_a_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
				a3=interpolation(t-dt,(int)a_0_b[j+2*n_a_0_b],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			}
			f[ind_a_0_b[i]]=f[ind_a_0_b[i]]+M[ind_a_0_b[i]+ind_a_0_b[j]*n_n]*((ind_a_0_b[j]<n_n/dim)*a1+
				((ind_a_0_b[j]>=n_n/dim)&&(ind_a_0_b[j]<2*n_n/dim))*a2+
				((ind_a_0_b[j]>=2*n_n/dim)&&(ind_a_0_b[j]<3*n_n/dim))*a3);
		}
    }
}
