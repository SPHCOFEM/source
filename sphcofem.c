/* Copyright 2000, 2002, 2015, 2020 - 2023 University of West Bohemia
BSD 3-Clause "New" or "Revised" License
*/

/* Errors
-> XSPH term has opposite sign
-> predictor-corrector (2) ustable for broken dam 
-> tetrahedron elasticity does not work
-> predictor-corrector leapfrog (3) does not update positions
*/

/* ToDo
-> rotation velocity (matrix multiplication) -> updated Lagrangean
-> hyperelastic material
-> contact damping
   -> update manual
	  -> klin
	  -> knon
	  -> kf
	  -> kd
-> test triangles and quads in 2D
-> rigid body element volume in 3D (mark difference between rectangle and tetrahedron)
-> rigid body dynamics (angular motion and forces and moments from particles and nodes to centre of gravity)
-> total Lagrangean formulation (initial kernel w and rho at t = 0)
-> model and results save and view in Python

-> t_s[n_s] as a particle type (0 = domain, 1 = interface, 2 = wall)
-> optionally change between small and large deformations
-> better nearest neighbour search (NNS)
-> more kernel types as input option

Leapfrog method: coefficient leapfrog = 1 (move velocity) / 2 (move position)
-> a[i]
-> v[i+1/2]=v[i-1/2]+a[i]*dt -> v[i+1]=v[i]+a[i+1/2]*dt
-> x[i+1]=x[i]+v[i+1/2]*dt   -> x[i+3/2]=x[i+1/2]+v[i+1]*dt
*/

/* Remarks
c_s[n_s] = 10.0 for gas shock tube convergence
Data (coordinates, velocity, accelerations, forces, moments, boundary conditions) are decalred as full 3D
-> spatial arrays of size 3 times number of values
-> calculations are restricted to particular dimension
e_s, dedt_s, O_s, S_s, dSdt_s declared as [6 x n] matrix [xx, xy, xz, yy, yz, zz]
([0, 1, 2, 3, 4, 5]) in all dimensions
-> for saving (store.c):
	1D: e_s, dedt_s scalars 1x1 [0], O_s, S_s, dSdt_s are zero (not relevant) []
	2D: e_s, dedt_s, S_s, dSdt_s symmetric 2x2 [0, 1, 3], O_s antisymmetric 2x2 [1] (scalar)
	3D: e_s, dedt_s, S_s, dSdt_s symmetric 3x3 [0, 1, 2, 3, 4, 5], O_s antisymmetric 3x3 [1, 2, 4]

For zero dt_save each time step saved
Maximum time reached at line 7439 must be plus epsilon
Positive dt_max sets initial time step in each loop, otherwise previous time step is taken as initial
Higher inital time step causes instability (optimum initial time step equals 0.1)
Energy (doubles) updated inside functions by addressing using &

Particle acceleration a_s[3*n_s] = f_s[3*n_s] / m_s[n_s]
Perticle force f_s[3*n_s] = m_s[n_s] * a_s[3*n_s]
Nodal acceleration a_b[n_b]
Nodal force f_b[n_b]
Rigid body acceleration a_r[3*n_r]
Rigid body force f_r[3*n_r]
Rigid body rotational acceleration alpha_r[3*n_r]
Rigid body moment M_r[3*n_r]

Functions, materials, domains, contacts, particles, nodes and elements cannot have negative numbers
Rigid body membranes or shells in 3D can be only triangles (otherwise 4 nodes define tetrahedron)
BOUNC type -1 keeps initial conditions
Element normal must be outside domain
Orthogonaity check must be 1e6*EPS
*/

/* Material models
1 : gas
2 : liquid
3 : SPH
4 : FEM bar (2D, 3D), membrane (3D)
5 : FEM beam (2D, 3D), shell (3D)
6 : FEM tria/quad (2D), tetrahedron (3D)
7 : SPS (SPH + shear stress)
8 : SPL (linear solid)
9 : SPY (hyperelasic solid)
12 : liquid (2) with tension
13 : SPH (3) with tension
23 : Mue-Grueneisen EOS (SPH)
27 : Mue-Grueneisen EOS with shear (SPS)
*/

/* Viscosity models
0 : artificial viscosity only (default)
1 : second order viscous term
2 : Monaghan-Cleary-Gingold (2006)
3 : Morris et al. (1997)
4 : Takeda et al. (1994)
5 : Onderik et al. (2007)
6 : Monaghan and Gingold (1983)
*/

/* Integration
0 : Euler method
1 : central acceleration (default)
2 : predictor-corrector
3 : predictor-corrector leapfrog
*/

/* XSPH
0 : off (default)
1 : predictor motion updated
2 : corrector motion updated
3 : both predictor and corrector motions updated
*/

#include "header.h"
#include "move.h"
#include "accelerations.h"
#include "rigid_bodies.h"
#include "contacts.h"
#include "nns.h"
#include "fem.h"
#include "cgm.h"
#include "store.h"
#include "freeall.h"
#include "volume.h"
#include "energy.h"
#include "time_step.h"
#include "interpolation.h"

#undef __FUNC__
#define __FUNC__ "main"

int main(int argc,char *argv[])
{
	/* file */
	FILE *sph,*msg,*sig;

	/* strings */
  	int ch; // n_max=0;
  	char keyword[10],output[80],textout[80],msgout[80],logout[80];
	char dt_max_type[20],dt_save_type[20];
	char data_check_type[20],data_print_type[20],mem_check_type[20];
	char nns_type[20],opt_type[20],int_scheme[40],xsph[20],is_type[40];
	char contact_type[40],ct_type[20];
	char rb_type[40];

	/* switches */
	char swon[10],swof[10];
	strcpy(swon,"1 (on)");
	strcpy(swof,"0 (off)");

	/* default values */
	strcpy(output,OUTFILE);
	strcpy(textout,TXTFILE);
	strcpy(msgout,MSGFILE);
	strcpy(logout,LOGNAME); 
	strcpy(dt_save_type,"0.0 (each cycle)");
	strcpy(dt_max_type,"0.0 (not bounded)");
	strcpy(data_check_type,swon);
	strcpy(data_print_type,"1 (standard output)");
	strcpy(mem_check_type,swon);
	strcpy(nns_type,"0 (each cycle)");
	strcpy(opt_type,"0.0 (ignored)");
	strcpy(int_scheme,"1 (central acceleration)");
	strcpy(xsph,swof);
	strcpy(is_type,"0 (artificial viscosity only)");

  	/* model */
  	int dim=3,integration=1,mem_check=1,is=0,ix=0,error=0; // default values
	int data_check=1,data_print=1,cycle_print=1,cycle_contact=1;  // default values
  	double cour=0.9,kstab=0.9,ax=0.0,ay=0.0,az=0.0; // default values

  	/* energy */
  	double kine_s,inne_s,pote_s;
  	double kine_b,defo_b,disi_b,pote_b;
  	double tote;
 
  	/* materials */
  	int *num_m=0,*domain_m=0,*type_m=0;
  	double *rho_m=0,*mu_m=0,*T_m=0,*gamma_m=0,*kappa_m=0;
  	double *coef1_m=0,*coef2_m=0,*coef3_m=0,*coef4_m=0;

	/* functions */
	int *num_u=0,*fun_u=0;
	int count_line_fun=0,count_fun_pairs=0,*fun_loc=0;
	double *fmx_u=0,*fmy_u=0,*fdx_u=0,*fdy_u=0,*fxi_u=0,*fyi_u=0;

  	/* contacts */
  	int *seg_c=0,*slave_c=0,*master_c=0;
  	int *mat_c=0,*num_c=0,*sw_c=0;
  	double *m_slave=0,*m_master=0;
  	double *ct_c=0,*klin_c=0,*knon_c=0,*kf_c=0,*kd_c=0;
  	double *f_c=0; // if *f_c is not set, warning that may be uninitialized

	/* constraints */
	int *constrained_s=0,*constrained_b=0,*constrained_r=0; // constraint attributes (0 or 1)
	int *constrained_s_frame=0,*constrained_b_frame=0;
	int *constrained_r_frame=0; // rigid body coordinate systems frames (0 or 1)
	//int i_unc_s=0,i_unc_b=0,n_unc_s=0,n_unc_b=0; // number of unconstrained particles and nodes
	//int *unconstrained_s=0,*unconstrained_b=0; // list of unconstrained particles and nodes

	/* rigid bodies */
	double *m_r=0,*I1_r=0,*I2_r=0,*I3_r=0,*D1_r=0,*D2_r=0,*D3_r=0;
	int *num_r=0,*type_r=0,*mat_r=0,*COG_r=0,*N1_r=0,*N2_r=0,*N3_r=0;
	int loc=0; // local coordinate frame
	int *COG_r_type=0,*N1_r_type=0,*N2_r_type=0,*N3_r_type=0;
  	double *x_r=0,*psi_r=0,*v_r=0,*o_r=0;
  	double *x_r_0=0,*psi_r_0=0;
	double *x_r_o=0,*psi_r_o=0,*v_r_o=0,*o_r_o=0;
	//double *x_r_a=0,*psi_r_a=0,*v_r_a=0,*o_r_a=0;
	double *a_r=0,*alpha_r=0;
	double *a_r_o=0,*alpha_r_o=0;
	double *f_r=0,*M_r=0;
	int *type_a_0_r=0,*frame_a_0_r=0;
	double *a_0_r=0;
	int *type_f_0_r=0,*frame_f_0_r=0;
	double *f_0_r=0,*M_0_r=0;	
	//double *x1_r=0,*x2_r=0,*x3_r=0;
	//double *x1_r_o=0,*x2_r_o=0,*x3_r_o=0;
	double *u1_r=0,*u2_r=0,*u3_r=0; // rigid bodies local axes system
	int *constrained_s_rb=0,*constrained_b_rb=0;
	int i_r_0_s=0,n_r_0_s=0,i_r_0_b=0,n_r_0_b=0;
	int *rb_s=0,*rb_b=0;
	double *l_s=0,*l_b=0;
	int *ind_r_0_s=0,*ind_r_0_b=0;
	int count_rb_s=0,count_rb_b=0;
	double ux=1.0,uy=0.0,uz=0.0,un=1.0;
	double vx=0.0,vy=1.0,vz=0.0,vn=1.0;
	double wx=0.0,wy=0.0,wz=0.0,wn=1.0;
	double rx=0.0,ry=0.0,rz=0.0,ox=0.0,oy=0.0,oz=0.0;

  	/* SPH */
  	int *mat_s=0,*num_s=0;
	//int *frame_v_0_s=0;
  	int *ind_f_0_s=0,*type_f_0_s=0,*frame_f_0_s=0;
	int *ind_a_0_s=0,*type_a_0_s=0,*frame_a_0_s=0;
	int *ind_o_0_s=0,*type_o_0_s=0,*o_0_s=0;
  	double *f_0_s=0,*M_0_s=0,*a_0_s=0,*dv_s=0;
	double *x_s_0=0,*V_s=0,*m_s=0,*mu_s=0,*c_s_0=0,*h_s_0=0;
  	double *x_s=0,*v_s=0;
	double *x_s_o=0,*v_s_o=0;
	//double *x_s_a=0,*v_s_a=0;
  	double *rho_s=0,*u_s=0,*p_s=0,*c_s=0,*h_s=0;
	double *rho_s_o=0,*u_s_o=0; //*p_s_o=0,*c_s_o=0,*h_s_o=0;
	//double *rho_s_a=0,*u_s_a=0; //*p_s_a=0,*c_s_a=0,*h_s_a=0;
	double *a_s=0;
	double *a_s_o=0;
	double *f_s=0;
	double *drhodt_s=0,*dudt_s=0;
  	double *drhodt_s_o=0,*dudt_s_o=0;
  	double *e_s=0,*S_s=0,*O_s=0;
	double *e_s_o=0,*S_s_o=0; //double *O_s_o=0;
	//double *e_s_a=0,*S_s_a=0;
	double *dedt_s=0,*dSdt_s=0;
	double *dedt_s_o=0,*dSdt_s_o=0;
	//double *gradvx_s=0,*gradvy_s=0,*gradvz_s=0;
  	double alpha=1.2,beta=1.5,eta=0.1,zeta=0.3,nas=4.0,theta=0.01; // default values
  	double v_max=0.0;
  	double sigma=1.0/PI,h0=1.4,xeps=0.5; // constants and default values

  	/* nearest neigbour optimization */
  	int nnopt=0; // default value
  	int *nn_s=0,*neig_s=0;
  	double opt=1.0; // default value

  	/* FEM */
  	int *num_b=0,*mat_b=0;
	//int *frame_v_0_b=0;
  	int *ind_f_0_b=0,*type_f_0_b=0,*frame_f_0_b=0;
	int *ind_a_0_b=0,*type_a_0_b=0,*frame_a_0_b=0,*ind_d_0_b=0,*ind_o_0_b=0,*type_o_0_b=0,*o_0_b=0;
  	double *f_0_b=0,*M_0_b=0,*a_0_b=0,*d_0_b=0;
	double *x_b_0=0,*m_b=0;
  	double *x_b=0,*v_b=0,*a_b=0,*f_b=0;
	double *x_b_o=0,*v_b_o=0,*a_b_o=0;
	//double *x_b_a=0,*v_b_a=0;
  	double *M=0,*B=0,*K=0,*f=0;
  	int *nod_e=0,*mat_e=0,*num_e=0;
  	double *V_e=0;
  	double c0=0.0,c1=0.0,dr=0.0;

  	/* initial and boundary conditions */
  	int *inpre=0;
  	double *p_0=0;
	int *invel=0,*invel_frame=0;
  	double *v_0=0,*omega_0=0;
	int *acfld=0,*acfld_type=0,*acfld_frame=0;
	double *a_0=0;
	int *force=0,*force_type=0,*force_frame=0;
	double *f_0=0,*M_0=0;
	int *ndamp=0;
	double *d_0=0;
	int *bounc=0,*bounc_frame=0,*bounc_type=0,*o_0=0;

  	/* save */
  	int save_dt=0;
  	int save_kine_s=0,save_inne_s=0,save_pote_s=0;
	int save_kine_b=0,save_disi_b=0,save_defo_b=0,save_pote_b=0;
  	int save_tote=0;
  	int save_v_s=0,save_dv_s=0,save_a_s=0,save_f_s=0;
  	int save_rho_s=0,save_drhodt_s=0,save_u_s=0,save_dudt_s=0,save_p_s=0,save_c_s=0,save_h_s=0;
  	int save_S_s=0,save_dSdt_s=0,save_e_s=0,save_dedt_s=0,save_O_s=0;
  	int save_v_b=0,save_a_b=0,save_f_b=0;
  	int save_f_c=0;
	int save_x_r=0,save_psi_r=0,save_v_r=0,save_o_r=0;
	int save_a_r=0,save_alpha_r=0,save_f_r=0,save_M_r=0;
  
	/* save message */
  	char save_dt_type[10];
  	char save_kine_s_type[10],save_inne_s_type[10],save_pote_s_type[10];
	char save_kine_b_type[10],save_disi_b_type[10],save_defo_b_type[10],save_pote_b_type[10];
  	char save_tote_type[10];
  	char save_v_s_type[10],save_dv_s_type[10],save_a_s_type[10],save_f_s_type[10];
  	char save_rho_s_type[10],save_drhodt_s_type[10],save_u_s_type[10],save_dudt_s_type[10];
	char save_p_s_type[10],save_c_s_type[10],save_h_s_type[10];
  	char save_S_s_type[10],save_dSdt_s_type[10],save_e_s_type[10],save_dedt_s_type[10],save_O_s_type[10];
  	char save_v_b_type[10],save_a_b_type[10],save_f_b_type[10];
  	char save_f_c_type[10];
  	char save_x_r_type[10],save_psi_r_type[10],save_v_r_type[10],save_o_r_type[10];
	char save_a_r_type[10],save_alpha_r_type[10],save_f_r_type[10],save_M_r_type[10];

	/* default values */
  	strcpy(save_dt_type,swof);
  
  	strcpy(save_kine_s_type,swof);
	strcpy(save_inne_s_type,swof);
	strcpy(save_pote_s_type,swof);;

	strcpy(save_kine_b_type,swof);
	strcpy(save_disi_b_type,swof);
	strcpy(save_defo_b_type,swof);
	strcpy(save_pote_b_type,swof);

  	strcpy(save_tote_type,swof);

  	strcpy(save_v_s_type,swof);
	strcpy(save_dv_s_type,swof);
	strcpy(save_a_s_type,swof);
	strcpy(save_f_s_type,swof);

  	strcpy(save_rho_s_type,swof);
	strcpy(save_drhodt_s_type,swof);
	strcpy(save_u_s_type,swof);
	strcpy(save_dudt_s_type,swof);

	strcpy(save_p_s_type,swof);
	strcpy(save_c_s_type,swof);
	strcpy(save_h_s_type,swof);

  	strcpy(save_S_s_type,swof);
	strcpy(save_dSdt_s_type,swof);
	strcpy(save_e_s_type,swof);
	strcpy(save_dedt_s_type,swof);
	strcpy(save_O_s_type,swof);

  	strcpy(save_v_b_type,swof);
	strcpy(save_a_b_type,swof);
	strcpy(save_f_b_type,swof);

  	strcpy(save_x_r_type,swof);
  	strcpy(save_psi_r_type,swof);
  	strcpy(save_v_r_type,swof);
  	strcpy(save_o_r_type,swof);
  	strcpy(save_a_r_type,swof);
  	strcpy(save_alpha_r_type,swof);
  	strcpy(save_f_r_type,swof);
  	strcpy(save_M_r_type,swof);

  	strcpy(save_f_c_type,swof);

  	/* counters */  
  	int i,j,k,l;
  	int count_line=1,count_cycle=0,count_save=0;
	int count_slave=0,count_master=0,count_char;//,count_neig=0;
  	int n_u=0,n_m=0,n_s=0,n_v=0,n_p=0,n_f=0,n_a=0,n_o=0,n_r=0,n_b=0,n_d=0,n_e=0,n_c=0,n_n=0;
  	int i_u=0,i_m=0,i_r=0,i_s=0,i_v=0,i_p=0,i_f=0,i_a=0,i_o=0,i_b=0,i_d=0,i_e=0,i_c=0;
  	int n_f_0_s=0,n_a_0_s=0,n_o_0_s=0;
  	int n_f_0_b=0,n_a_0_b=0,n_d_0_b=0,n_o_0_b=0;
  	int i_f_0_s=0,i_a_0_s=0,i_o_0_s=0;
  	int i_f_0_b=0,i_a_0_b=0,i_d_0_b=0,i_o_0_b=0;

	/* check switches */
  	int iwarning=0,ioutput=0,idim=0,itmax,itime=0,igacc=0,ioptim=0,isph=0,ifem=0,isave=0;
	int imate=0,inode=0,ibfun=0,icontact=0;

  	/* time */
  	int time_h,time_m;
  	double t=0.0,dt,t_max,dt_save,t_save=0.0,dt_init=0.1,dt_max=0.0,time_s; // default values

	/* stopwatch */
  	clock_t t_start=clock();
  	//long clk_tck=sysconf(_SC_CLK_TCK);

  	/* integration parameters */
  	double prei=0.5,cori=1.0; // default values
		  
	/* auxiliary variables */
  	double aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8,aux9,aux10;
  	double aux11,aux12,aux13,aux14,aux15,aux16,aux17,aux18,aux19,aux20;
  	double aux21,aux22,aux23,aux24,aux25,aux26,aux27,aux28,aux29,aux30;
	double aux31,aux32,aux33,aux34,aux35,aux36,aux37;
  
  	/* Mie-Gruenesen EOS */
  	double Eta,s=1.87,G0=0.17,En; // default values
  	double C0,C1,C2,C3,C4,C5;

  	io_alloc(&g_ioc);
  	io_setProgName(&g_ioc,argv[0]);
  	g_ioc.stdoutTermOutput=1;
  	g_ioc.stderrTermOutput=1;
  	if (atexit( &dummy_mem_stats ))
    {
		perror("error:");
		errput("cannot registrate atexit() function (%d)\n!", errno);
    }
	
	strcpy(g_ioc.stdoutLogFileName,logout);
  	strcpy(g_ioc.stderrLogFileName,logout);

  	/* error file */
  	msg=fopen(msgout,"w");
	
  	fprintf(msg,
		"---------------"
		"\n"
		"SPHCOFEM solver"
		"\n"
		"---------------"
		"\n"
		"analyzing input file ..."
		"\n"
		"\n");

  	/* temporary input file */
  	sph=fopen(TMP,"w");
  	while ((ch=getc(stdin))!=EOF) putc(ch,sph);
	/* copy without commnets - cannot allocate line number for errors 
  	while ((ch=getc(stdin))!=EOF)
	{
		if (ch=='$') while ((ch=getc(stdin))!='\n');
		else putc(ch,sph);
	}
	*/
  	fclose(sph);

   	/* detection of numbers and correct lines */
  	sph=fopen(TMP,"r");

  	while ((ch=getc(sph))!=EOF)
    {
		count_line++;
		
		if (ch!='$')
		{
			/* reading keyword of fixed length
			str[0]=ch;
	  		str[1]=getc(sph);
	  		str[2]=getc(sph);
	  		str[3]=getc(sph);
	  		str[4]=getc(sph);
	  		str[5]='\0';
			*/

			/* reading chars till space or end of line or end of file */
			count_char=0;

			while ((ch!=' ')&&(ch!='\n')&&(ch!=EOF))
			{
				keyword[count_char]=ch;
				ch=getc(sph);
				count_char++;
			}

			/* input file keyword */
			keyword[count_char]='\0';

	  		if (!strcmp(keyword,"NAME"))
			{
				ioutput=1;
				fscanf(sph,
					"%s\n",
					output);
			}
			else if (!strcmp(keyword,"DIM"))
			{
				idim=1;
				fscanf(sph,
					"%lf",
		     		&(aux1));
				dim=(int)aux1;
		    }
			else if (!strcmp(keyword,"TMAX"))
			{
				itmax=1;
				fscanf(sph,
					"%lf",
		     		&(t_max));
		    }
			else if (!strcmp(keyword,"TIME"))
			{
				itime=1;
				fscanf(sph,
					"%lf %lf %lf %lf %lf",
					&(dt_save),
					&(dt_init),
					&(dt_max),
					&(cour),
					&(kstab));
			}
			else if (!strcmp(keyword,"GACC"))
			{
				igacc=1;
				fscanf(sph,
					"%lf %lf %lf",
		     		&(ax),
		     		&(ay),
		     		&(az));
            }
			else if (!strcmp(keyword,"OPTIM"))
			{
				ioptim=1;
				fscanf(sph,
					"%lf %lf %lf %lf %lf %lf %lf %lf",
					&(aux1),
					&(aux2),
					&(aux3),
					&(aux4),
					&(opt),
					&(aux5),
					&(aux6),
					&(aux7));
				data_check=(int)aux1;
				data_print=(int)aux2;
				cycle_print=(int)aux3;
				nnopt=(int)aux4;
				cycle_contact=(int)aux5;
				integration=(int)aux6;
				mem_check=(int)aux7;
			}
			else if (!strcmp(keyword,"SPH"))
			{
				isph=1;
				fscanf(sph,
					"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&(aux1),
					&(alpha),
					&(beta),
					&(eta),
					&(zeta),
					&(nas),
					&(theta),
					&(aux2),
					&(xeps)
					);
				is=(int)aux1;
				ix=(int)aux2;
			}
			else if (!strcmp(keyword,"FEM"))
			{
				ifem=1;
				fscanf(sph,
					"%lf %lf",
					&(c0),
					&(c1));
			}
			else if (!strcmp(keyword,"SAVE"))
			{
				isave=1;
				fscanf(sph,
					"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
					"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
					"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
					"%lf %lf %lf %lf %lf %lf %lf",
					&(aux1),
					&(aux2),
					&(aux3),
					&(aux4),
					&(aux5),
					&(aux6),
					&(aux7),
					&(aux8),
					&(aux9),
					&(aux10),
					&(aux11),
					&(aux12),
					&(aux13),
					&(aux14),
					&(aux15),
					&(aux16),
					&(aux17),
					&(aux18),
					&(aux19),
					&(aux20),
					&(aux21),
					&(aux22),
					&(aux23),
					&(aux24),
					&(aux25),
					&(aux26),
					&(aux27),
					&(aux28),
					&(aux29),
					&(aux30),
					&(aux31),
					&(aux32),
					&(aux33),
					&(aux34),
					&(aux35),
					&(aux36),
					&(aux37));
				save_dt=(int)aux1;
			
				save_kine_s=(int)aux2;
				save_inne_s=(int)aux3;
				save_pote_s=(int)aux4;
				save_kine_b=(int)aux5;
				save_disi_b=(int)aux6;
				save_defo_b=(int)aux7;
				save_pote_b=(int)aux8;

				save_tote=(int)aux9;

				save_v_s=(int)aux10;
				save_dv_s=(int)aux11;
				save_a_s=(int)aux12;
				save_f_s=(int)aux13;

				save_rho_s=(int)aux14;
				save_drhodt_s=(int)aux15;
				save_u_s=(int)aux16;
				save_dudt_s=(int)aux17;
				save_p_s=(int)aux18;
				save_c_s=(int)aux19;
				save_h_s=(int)aux20;

				save_S_s=(int)aux21;
				save_dSdt_s=(int)aux22;
				save_e_s=(int)aux23;
				save_dedt_s=(int)aux24;
				save_O_s=(int)aux25;

				save_v_b=(int)aux26;
				save_a_b=(int)aux27;
				save_f_b=(int)aux28;

				save_x_r=(int)aux29;
				save_psi_r=(int)aux30;
				save_v_r=(int)aux31;
				save_o_r=(int)aux32;
				save_a_r=(int)aux33;
				save_alpha_r=(int)aux34;
				save_f_r=(int)aux35;
				save_M_r=(int)aux36;

				save_f_c=(int)aux37;
			}
			else if (!strcmp(keyword,"FUNCT"))
			{
				n_u++;
				fscanf(sph,
					"%lf %lf", // %lf %lf %lf %lf",
					&(aux1),
					&(aux2));
				count_line_fun=(int)aux2;
				count_fun_pairs=count_fun_pairs+count_line_fun;

				while ((ch=getc(sph))!='\n'); // reading until end of line

				for (i=0;i<count_line_fun;i++) // read all function pairs
				{
					count_line++;
					if (i<count_line_fun-1)
					{
						while ((ch=getc(sph))!='\n'); // reading until end of line
					}
				}
			}
			else if (!strcmp(keyword,"MATER"))
			{
				n_m++;
				fscanf(sph,
		     		"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		     		&(aux1),
		     		&(aux2),
		     		&(aux3),
		     		&(aux4),
		     		&(aux5),
		     		&(aux6),
		     		&(aux7),
		     		&(aux8),
		     		&(aux9),
		     		&(aux10),
		     		&(aux11));
				if ((int)aux3==0) n_r++; // rigid body
			}
			else if (!strcmp(keyword,"SNODE")) n_s++;
			else if (!strcmp(keyword,"BNODE")) n_b++;
			else if (!strcmp(keyword,"BELEM")) n_e++;
			else if (!strcmp(keyword,"INVEL")) n_v++;
			else if (!strcmp(keyword,"INPRE")) n_p++;
			else if (!strcmp(keyword,"FORCE")) n_f++;
			else if (!strcmp(keyword,"ACFLD")) n_a++;
			else if (!strcmp(keyword,"BOUNC")) n_o++;
			else if (!strcmp(keyword,"NDAMP")) n_d++;
			else if (!strcmp(keyword,"CONTACT")) n_c++;
			else if (!strcmp(keyword,"END")) break;
			else
			{
				fprintf(msg,
					"ERROR: unknown keyword %s on input line %d!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					keyword,count_line);
				fclose(sph);
				fclose(msg);
				remove(TMP);
				return(0);	  
			}
		}
	    while ((ch=getc(sph))!='\n'); // reading until end of line
	}
	fclose(sph);

  	/* test on output file */
  	if (!ioutput)
    {
		fprintf(msg,
			"WARNING: output file is not defined, default value %s is used!"
			"\n"
			"\n",
			output);
		iwarning=1;
    }
  	else 
	{
		strcpy(textout,output);
		strcpy(msgout,output);
		strcpy(logout,output);
		strcpy(&output[strlen(output)],".out\0");
		strcpy(&textout[strlen(textout)],".txt\0");
		strcpy(&msgout[strlen(msgout)],".msg\0");
		strcpy(&logout[strlen(logout)],".log\0");
	}

  	/* test on dimension */
  	if (!idim)
    {
    	fprintf(msg,
	    	"WARNING: dimension is not defined, default value %d is used!"
			"\n"
	      	"\n",
			dim);
		iwarning=1;
    }
	else
	{
		/* test on dimension */
		if (dim==0)
		{
			fprintf(msg,
				"WARNING: dim = %d",
				dim);
			dim=3;
			fprintf(msg,
				" -> default value %d is used!"
				"\n"
				"\n",
				dim);
			iwarning=1;
		}
		else if ((dim!=1)&&(dim!=2)&&(dim!=3))
		{
			fprintf(msg,
				"ERROR: dimension must be 1 or 2 or 3!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
	}

	switch (dim)
	{
		case 1:
		{
			sigma=2.0/3.0;
			h0=1.4;
			break;
		}
		case 2:
		{
			sigma=10.0/7.0/PI;
			h0=1.4*sqrt(PI);
			break;
		}
		case 3:
		{
			sigma=1.0/PI;
			h0=1.4;
			break;
		}
	}

  	/* test on termination time */
  	if (!itmax)
    {
    	fprintf(msg,
	    	"ERROR: termination time must be defined!"
			"\n"
			"\n"
			"error termination"
			"\n"
	      	"\n");
		fclose(msg);
      	remove(TMP);
      	return(0);
    }
	else
	{
		/* test on t_max */
		if (t_max<0.0)
		{
			fprintf(msg,
				"ERROR: termination time must be 0 or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
	}

  	/* test on time controls */
  	if (!itime)
    {
      	fprintf(msg,
	    	"WARNING: controls are not defined, default values are used!"
	      	"\n"
			"\n");
		iwarning=1;
    }
	else
	{
		/* test on dt_save */
		if (dt_save<0.0)
		{
			fprintf(stdout,
				"ERROR: dt_save must be 0 or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		else
		{
			if (dt_save==0.0)
			{
				strcpy(dt_save_type,"0.0 (each state)");
			}
			else
			{
				sprintf(dt_save_type,"%f",dt_save);
			}
		}

		/* test on dt_init */
		if (dt_init==0.0)
		{
			fprintf(msg,
				"WARNING: dt_init = %f",
				dt_init);
			dt_init=0.1;
			fprintf(msg,
				" -> default value %f is used!"
				"\n"
				"\n",
				dt_init);
			iwarning=1;
		}
		else if (dt_init<0.0)
		{
			fprintf(msg,
				"ERROR: dt_init must be 0 or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}

		/* test on dt_max */
		if (dt_max<0.0)
			{
			fprintf(msg,
				"ERROR: dt_max be zero or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		else
		{
			if (dt_max==0.0)
			{
				strcpy(dt_max_type,"0.0 (not bounded)");
			}
			else
			{
				sprintf(dt_max_type,"%f",dt_max);
			}
		}

		/* test on cour */
		if ((n_s>0)&&(cour==0.0))
		{
			fprintf(stdout,
				"WARNING: cour = %f",
				cour);
			cour=0.9;
			fprintf(stdout,
				" -> default value %f is used!"
				"\n"
				"\n",
				cour);
			iwarning=1;
		}

		/* test on kstab */
		if ((n_b>0)&&(kstab==0.0))
		{
			fprintf(stdout,
				"WARNING: kstab = %f",
				kstab);
			kstab=0.9;
			fprintf(stdout,
				" -> default value %f is used!"
				"\n"
				"\n",
				kstab);
			iwarning=1;
		}
	}

  	/* test on global acceleration field */
  	if (!igacc)
    {
      	fprintf(msg,
	    	"WARNING: global acceleration field is not defined, default values are used!"
			"\n"
	      	"\n");
		iwarning=1;
	}

  	/* test on optimization */
  	if (!ioptim)
    {
      	fprintf(msg,
	    	"WARNING: optimization controls are not defined, default values are used!"
			"\n"
	      	"\n");
		iwarning=1;
	}
	else
	{
		/* test on data_check */
		if ((data_check!=0)&&(data_check!=1))
		{
			fprintf(msg,
				"ERROR: data_check must be 0 or 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		else
		{
			if (data_check==1)
			{
				strcpy(data_check_type,swon);
			}
			else
			{
				strcpy(data_check_type,swof);
			}
		}
	
		/* test on data_print */
		if ((data_print!=0)&&(data_print!=1)&&(data_print!=2))
		{
			fprintf(msg,
				"ERROR: data_print must be 0 or 1 or 2!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		else
		{
			switch (data_print)
			{
				case 0:
				{
					strcpy(data_print_type,"0 (no output)");
					freopen("/dev/null", "w", stdout);
					break;
				}
				case 1:
				{
					strcpy(data_print_type,"1 (standard output)");
					break;
				}
				case 2:
				{
					strcpy(data_print_type,"2 (file output)");
					freopen(textout, "w", stdout);
					break;
				}
			}
		}

		/* test on cycle_print */
		if ((cycle_print<0)||(cycle_print!=(int)cycle_print))
		{
			fprintf(msg,
				"ERROR: cycle_print must be 0 or an integer value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}

		/* test on nearest neighbour search (0 = each cycle, 1 = once at beginning, N = each nopt cycles) */
		if ((nnopt!=0)&&(nnopt!=(int)nnopt))
		{
			fprintf(msg,
				"ERROR: nearest neighnour search method must be 0 or an integer value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}

		/* test on opt */
		if (nnopt>0) // clustering nns
		{
			if (opt<0.0)
			{
				fprintf(msg,
					"ERROR: coeffitient opt must be 0.0 or a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n");
				fclose(msg);
				remove(TMP);
				return(0);
			}
			else if (opt==0.0)
			{
				fprintf(msg,
					"WARNING: opt = %f",
					opt);
				opt=1.0;
				fprintf(msg,
					" -> default value %f is used!"
					"\n"
					"\n",
					opt);
				iwarning=1;
			}
		}

		switch (nnopt)
		{
			case 0:
			{
				strcpy(nns_type,"0 (no nns clustering)");
				strcpy(opt_type,"0.0 (ignored)");
				break;
			}
			case 1:
			{
				strcpy(nns_type,"1 (nns clustering at beginning)");
				sprintf(opt_type,"%f",opt);
				break;
			}
			default:
			{
				sprintf(nns_type,"%d (nns clustering each %d cycles)",nnopt,nnopt);
				sprintf(opt_type,"%f",opt);
			}
		}
	
		/* test on cycle_contact */
		if (cycle_contact!=(int)cycle_contact)
		{
			fprintf(msg,
				"ERROR: cycle_contact must be an integer value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		
		/* test on integration scheme (0 = Euler, 1 = central acceleration, 
		2 = predictor-corrector, 3 = predictor-corrector leapfrog) */
		if ((integration!=0)&&(integration!=1)&&(integration!=2)&&(integration!=3))
			{
			fprintf(msg,
				"ERROR: integration method must be 0 or 1 or 2 or 3, defaule value 1 is used!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		else
		{
			/* integration types and constants */
			switch (integration)
			{
				case 0:
				{
					//prei=0.0; // no predictor step
					strcpy(int_scheme,"0 (Euler method)");
					cori=1.0;
					break;
				}
				case 1:
				{
					strcpy(int_scheme,"1 (central acceleration)");
					prei=0.5;
					cori=1.0;
					break;
				}
				case 2:
				{
					strcpy(int_scheme,"2 (predictor-corrector)");
					prei=1.0;
					cori=1.0;
					break;
				}
				case 3:
				{
					strcpy(int_scheme,"3 (predictor-corrector leapfrog)");
					prei=1.0;
					cori=0.5;
					break;
				}
			}
		}

		/* test on mem_check */
		if ((mem_check!=0)&&(mem_check!=1))
		{
			fprintf(msg,
				"ERROR: mem_check must be 0 or 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		else
		{
			if (mem_check==1)
			{
				strcpy(mem_check_type,swon);
			}
			else
			{
				strcpy(mem_check_type,swof);
			}
		}
	}

  	/* test on SPH */
  	if ((n_s>0)&&(!isph))
    {
      	fprintf(msg,
	    	"WARNING: SPH controls are not defined, default values are used!"
			"\n"
	      	"\n");
		iwarning=1;
    }
	else
	{
		/* test on XSPH */
		if ((n_s>0)&&(ix!=0)&&(ix!=1)&&(ix!=2)&&(ix!=3))
		{
			fprintf(msg,
				"ERROR: XSPH option must be 0 or 1 or 2 or 3!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
		else
		{
			switch (ix)
			{
				case 0:
				{
					strcpy(xsph,swof);
					break;
				}
				case 1:
				{
					strcpy(xsph,"1 (active for predictor)");
					break;
				}
				case 2:
				{
					strcpy(xsph,"2 (active for corrector)");
					break;
				}
				case 3:
				{
					strcpy(xsph,"3 (active for predictor and corrector)");
					break;
				}
			}
		}

		/* test on xeps */
		if ((n_s>1)&&((ix==1)||(ix==2)||(ix==3))&&((xeps<0.0)||(xeps>1.0)))
		{
			fprintf(msg,
				"ERROR: coeffitient xeps must be between 0 and 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			fclose(msg);
	      	remove(TMP);
			return(0);
		}
	} 

  	/* test on FEM */
  	if ((n_e>0)&&(!ifem))
    {
      	fprintf(msg,
	    	"WARNING: finite elements controls are not defined, default values are used!"
			"\n"
	      	"\n");
		iwarning=1;
    }
      
 	/* test on save */
  	if (!isave)
    {
      	fprintf(msg,
	    	"WARNING: save controls are not defined, default values are used!"
			"\n"
	      	"\n");
		iwarning=1;
    }
	else
	{
		if (save_dt) strcpy(save_dt_type,swon);

		if (n_s>0)
		{
			if (save_v_s) strcpy(save_v_s_type,swon);
			if ((save_dv_s)&&(ix)) strcpy(save_dv_s_type,swon);
			if (save_a_s) strcpy(save_a_s_type,swon);
			if (save_f_s) strcpy(save_f_s_type,swon);

			if (save_rho_s) strcpy(save_rho_s_type,swon);
			if (save_drhodt_s) strcpy(save_drhodt_s_type,swon);
			if (save_u_s) strcpy(save_u_s_type,swon);
			if (save_dudt_s) strcpy(save_dudt_s_type,swon);
			if (save_p_s) strcpy(save_p_s_type,swon);
			if (save_c_s) strcpy(save_c_s_type,swon);
			if (save_h_s) strcpy(save_h_s_type,swon);

			if (save_e_s) strcpy(save_e_s_type,swon);
			if (save_dedt_s) strcpy(save_dedt_s_type,swon);
			if (save_O_s) strcpy(save_O_s_type,swon);
			if (save_S_s) strcpy(save_S_s_type,swon);
			if (save_dSdt_s) strcpy(save_dSdt_s_type,swon);

			if (save_kine_s) strcpy(save_kine_s_type,swon);
			if (save_inne_s) strcpy(save_inne_s_type,swon);
			if (save_pote_s) strcpy(save_pote_s_type,swon);
		}

		if (n_b>0)
		{
			if (save_v_b) strcpy(save_v_b_type,swon);
			if (save_a_b) strcpy(save_a_b_type,swon);
			if (save_f_b) strcpy(save_f_b_type,swon);

			if (save_kine_b) strcpy(save_kine_b_type,swon);
			if (save_disi_b) strcpy(save_disi_b_type,swon);
			if (save_defo_b) strcpy(save_defo_b_type,swon);
			if (save_pote_b) strcpy(save_pote_b_type,swon);
		}

		if (save_tote) strcpy(save_tote_type,swon);
		
		if ((n_c>0)&&(save_f_c)) strcpy(save_f_c_type,swon);

		if (n_r>0)
		{
			if (save_x_r) strcpy(save_x_r_type,swon);
			if (save_psi_r) strcpy(save_psi_r_type,swon);
			if (save_v_r) strcpy(save_v_r_type,swon);
			if (save_o_r) strcpy(save_o_r_type,swon);
			if (save_a_r) strcpy(save_a_r_type,swon);
			if (save_alpha_r) strcpy(save_alpha_r_type,swon);
			if (save_f_r) strcpy(save_f_r_type,swon);
			if (save_M_r) strcpy(save_M_r_type,swon);
		}
	}

  	/* test on materials */
  	if (n_m==0)
    {
    	fprintf(msg,
			"ERROR: there are no materials in the input file!"
	    	"\n"
			"\n"
			"error termination"
			"\n"
	    	"\n");
		fclose(msg);
      	remove(TMP);
      	return(0);
    }

  	/* test on nodes */
  	if ((n_s==0)&&(n_b==0))
    {
    	fprintf(msg,
	    	"ERROR: there are no particles and nodes in the input file!"
	      	"\n"
			"\n"
			"error termination"
			"\n"
	      	"\n");
		fclose(msg);
      	remove(TMP);
      	return(0);
    }

	/* renaming log file */
	strcpy(g_ioc.stdoutLogFileName,logout);
  	strcpy(g_ioc.stderrLogFileName,logout);
 
	/* checking memory consistency */
  	if (mem_check)
    {
		/* log file */
      	g_ioc.stdoutLog=1;
      	g_ioc.stderrLog=1;
    }
  	else
    {
    	/* no log file */
      	g_ioc.stdoutLog=0;
      	g_ioc.stderrLog=0;
    }
  
  	/* closing message file */
	fprintf(msg,
		" -> no errors found"
		"\n");

	if (!iwarning)
	{
		fprintf(msg,
			" -> no warnings found"
			"\n"
			"\n"
			"closing error file"
			"\n");
		fclose(msg);
	}
	else
	{
		fprintf(msg,
			"\n"
			" -> warnings found"
			"\n"
			"\n"
			"closing error file"
			"\n");
		fclose(msg);
	}

	/* renaming message file */
	rename(MSGFILE,msgout);

	/* freeing strings (except output) 
	free(keyword);
	free(textout);
	
	free(dt_max_type);
	free(dt_save_type);

	free(data_check_type);
	free(data_print_type);
	free(mem_check_type);

	free(nns_type);
	free(opt_type);
	free(int_scheme);
	free(xsph);
	free(is_type);

	free(contact_type);
	free(ct_type);

	free(save_dt_type);
	free(save_kine_s_type);
	free(save_inne_s_type);
	free(save_pote_s_type);
	free(save_kine_b_type);
	free(save_disi_b_type);
	free(save_defo_b_type);
	free(save_pote_b_type);
	free(save_tote_type);
	free(save_v_s_type);
	free(save_dv_s_type);
	free(save_a_s_type);
	free(save_f_s_type);
	free(save_rho_s_type);
	free(save_drhodt_s_type);
	free(save_u_s_type);
	free(save_dudt_s_type);
	free(save_p_s_type);
	free(save_c_s_type);
	free(save_h_s_type);
	free(save_S_s_type);
	free(save_dSdt_s_type);
	free(save_e_s_type);
	free(save_dedt_s_type);
	free(save_O_s_type);
	free(save_v_b_type);
	free(save_a_b_type);
	free(save_f_b_type);
	free(save_x_r_type);
	free(save_psi_r_type);
	free(save_v_r_type);
	free(save_o_r_type);
	free(save_a_r_type);
	free(save_alpha_r_type);
	free(save_f_r_type);
	free(save_M_r_type);
	free(save_f_c_type);
	*/

  	/* reading input */
  	fprintf(stdout,
		"---------------"
		"\n"
		"SPHCOFEM solver"
		"\n"
		"---------------"
		"\n"
		"allocating memory ..."
		"\n");

  	/* allocation - functions */
	if (n_u>0)
	{
		num_u=createMemMore(int,n_u);
		fun_u=createMemMore(int,n_u);
		fmx_u=createMemMore(double,n_u);
		fmy_u=createMemMore(double,n_u);
		fdx_u=createMemMore(double,n_u);
		fdy_u=createMemMore(double,n_u);
		fun_loc=createMemMore(int,n_u);
		fxi_u=createMemMore(double,count_fun_pairs);
		fyi_u=createMemMore(double,count_fun_pairs);
		count_fun_pairs=0; // reset all functions data pairs counter
	}

  	/* allocation - materials */
	if (n_m>0)
	{
		rho_m=createMemMore(double,n_m);
		mu_m=createMemMore(double,n_m);
		T_m=createMemMore(double,n_m);
		gamma_m=createMemMore(double,n_m);
		kappa_m=createMemMore(double,n_m);
		num_m=createMemMore(int,n_m);
		type_m=createMemMore(int,n_m);
		domain_m=createMemMore(int,n_m);
		coef1_m=createMemMore(double,n_m);
		coef2_m=createMemMore(double,n_m);
		coef3_m=createMemMore(double,n_m);
		coef4_m=createMemMore(double,n_m);
	}

  	/* allocation - SPH nodes */
  	if (n_s>0)
    {
		x_s=createMemMore(double,3*n_s);
      	x_s_0=createMemMore(double,3*n_s);
      	x_s_o=createMemMore(double,3*n_s);
      	//x_s_a=createMemMore(double,3*n_s);
      	v_s=createMemMore(double,3*n_s);
      	v_s_o=createMemMore(double,3*n_s);
      	//v_s_a=createMemMore(double,3*n_s);
      	a_s=createMemMore(double,3*n_s);
      	a_s_o=createMemMore(double,3*n_s);
      	f_s=createMemMore(double,3*n_s);
      	rho_s=createMemMore(double,n_s);
      	rho_s_o=createMemMore(double,n_s);
      	//rho_s_a=createMemMore(double,n_s);
      	drhodt_s=createMemMore(double,n_s);
      	drhodt_s_o=createMemMore(double,n_s);
      	u_s=createMemMore(double,n_s);
      	u_s_o=createMemMore(double,n_s);
      	//u_s_a=createMemMore(double,n_s);
      	dudt_s=createMemMore(double,n_s);
      	dudt_s_o=createMemMore(double,n_s);
      	mu_s=createMemMore(double,n_s);
      	V_s=createMemMore(double,n_s);
      	m_s=createMemMore(double,n_s);
      	p_s=createMemMore(double,n_s);
      	//p_s_o=createMemMore(double,n_s);
      	//p_s_a=createMemMore(double,n_s);
      	c_s=createMemMore(double,n_s);
      	c_s_0=createMemMore(double,n_s);
      	//c_s_o=createMemMore(double,n_s);
      	//c_s_a=createMemMore(double,n_s);
      	h_s=createMemMore(double,n_s);
      	h_s_0=createMemMore(double,n_s);
      	//h_s_o=createMemMore(double,n_s);
      	//h_s_a=createMemMore(double,n_s);
      	e_s=createMemMore(double,6*n_s);
      	e_s_o=createMemMore(double,6*n_s);
      	//e_s_a=createMemMore(double,6*n_s);
      	S_s=createMemMore(double,6*n_s);
      	S_s_o=createMemMore(double,6*n_s);
      	//S_s_a=createMemMore(double,6*n_s);
      	dedt_s=createMemMore(double,6*n_s);
      	dedt_s_o=createMemMore(double,6*n_s);
      	dSdt_s=createMemMore(double,6*n_s);
      	dSdt_s_o=createMemMore(double,6*n_s);
      	O_s=createMemMore(double,6*n_s);
      	//O_s_o=createMemMore(double,6*n_s);
      	//gradvx_s=createMemMore(double,3*n_s);
      	//gradvy_s=createMemMore(double,3*n_s);
      	//gradvz_s=createMemMore(double,3*n_s);
      	mat_s=createMemMore(int,n_s);
      	num_s=createMemMore(int,n_s);
		//frame_v_0_s=createMemMore(int,n_s);
		constrained_s=createMemMore(int,6*n_s);
		constrained_s_frame=createMemMore(int,n_s);
		constrained_s_rb=createMemMore(int,n_s);

      	/* setting domain vectors to zeros */
      	memset(v_s,0,3*n_s*sizeof(double));
      	memset(a_s,0,3*n_s*sizeof(double));
      	memset(f_s,0,3*n_s*sizeof(double));
      	memset(e_s,0,6*n_s*sizeof(double));
      	memset(dedt_s,0,6*n_s*sizeof(double));
      	memset(S_s,0,6*n_s*sizeof(double));
      	memset(dSdt_s,0,6*n_s*sizeof(double));
      	memset(O_s,0,6*n_s*sizeof(double));
      	//memset(gradvx_s,0,3*n_s*sizeof(double));
      	//memset(gradvy_s,0,3*n_s*sizeof(double));
      	//memset(gradvz_s,0,3*n_s*sizeof(double));
      	memset(drhodt_s,0,n_s*sizeof(double));
      	memset(dudt_s,0,n_s*sizeof(double));
		memset(constrained_s,0,6*n_s*sizeof(int));
		memset(constrained_s_frame,0,n_s*sizeof(int));
		memset(constrained_s_rb,0,n_s*sizeof(int));

	  	/* velocity corrections in case of XSPH */
	  	if (ix)
	    {
			dv_s=createMemMore(double,3*n_s);
      		memset(dv_s,0,3*n_s*sizeof(double));
    	}
    }
  
  	/* allocation - boundary nodes */
  	if (n_b>0)
    {
      	x_b=createMemMore(double,3*n_b);
      	x_b_0=createMemMore(double,3*n_b);
      	x_b_o=createMemMore(double,3*n_b);
      	//x_b_a=createMemMore(double,3*n_b);
      	v_b=createMemMore(double,3*n_b);
      	v_b_o=createMemMore(double,3*n_b);
		//v_b_a=createMemMore(double,3*n_b);
      	a_b=createMemMore(double,3*n_b);
      	a_b_o=createMemMore(double,3*n_b);
      	f_b=createMemMore(double,3*n_b);
		m_b=createMemMore(double,n_b);
      	num_b=createMemMore(int,n_b);
      	mat_b=createMemMore(int,n_b);
		//frame_v_0_b=createMemMore(int,n_b);
		constrained_b=createMemMore(int,6*n_b);
		constrained_b_frame=createMemMore(int,n_b);
		constrained_b_rb=createMemMore(int,n_b);

      	/* setting boundary vectors to zeros */
      	memset(v_b,0,3*n_b*sizeof(double));
      	memset(a_b,0,3*n_b*sizeof(double));
      	memset(f_b,0,3*n_b*sizeof(double));
      	memset(constrained_b,0,6*n_b*sizeof(int));
      	memset(constrained_b_frame,0,n_b*sizeof(int));
		memset(constrained_b_rb,0,n_b*sizeof(int));
    }
  
  	/* allocation - boundary elements */
  	if (n_e>0)
    {
      	V_e=createMemMore(double,n_e);
      	mat_e=createMemMore(int,n_e);
      	num_e=createMemMore(int,n_e);
      	nod_e=createMemMore(int,4*n_e);
    }
  
  	/* allocation - contacts */
  	if (n_c>0)
    {
      	seg_c=createMemMore(int,2*n_c);
      	ct_c=createMemMore(double,n_c);
      	klin_c=createMemMore(double,n_c);
      	knon_c=createMemMore(double,n_c);
      	kf_c=createMemMore(double,n_c);
      	kd_c=createMemMore(double,n_c);
      	num_c=createMemMore(int,n_c);
      	sw_c=createMemMore(int,n_c);
      	mat_c=createMemMore(int,2*n_c);
      	f_c=createMemMore(double,3*n_c);

      	/* setting contact vectors to zero */
      	memset(seg_c,0,2*n_c*sizeof(int));
      	memset(f_c,0,3*n_c*sizeof(double));
    }

  	/* allocation - pressure i1nitial conditions */
  	if (n_p>0)
    {
		inpre=createMemMore(int,n_p);
      	p_0=createMemMore(double,n_p);
    }

  	/* allocation - velocity initial conditions */
  	if (n_v>0)
    {
		invel=createMemMore(int,n_v);
      	invel_frame=createMemMore(int,n_v);
      	v_0=createMemMore(double,3*n_v);
      	omega_0=createMemMore(double,3*n_v);
    }

  	/* allocation - accelerations */
  	if (n_a>0)
    {
      	acfld=createMemMore(int,n_a);
      	acfld_type=createMemMore(int,n_a);
      	acfld_frame=createMemMore(int,n_a);
      	a_0=createMemMore(double,3*n_a);
    }
  
  	/* allocation - forces */
  	if (n_f>0)
    {
      	force=createMemMore(int,n_f);
      	force_type=createMemMore(int,n_f);
      	force_frame=createMemMore(int,n_f);
      	f_0=createMemMore(double,3*n_f);
      	M_0=createMemMore(double,3*n_f);
    }
  
  	/* allocation - nodal dampings */
  	if (n_d>0)
    {
      	ndamp=createMemMore(int,n_d);
      	d_0=createMemMore(double,3*n_d);
    }
  
  	/* allocation - boundary conditions */
  	if (n_o>0)
    {
      	bounc=createMemMore(int,n_o);
      	bounc_type=createMemMore(int,n_o);
		bounc_frame=createMemMore(int,n_o);
      	o_0=createMemMore(int,6*n_o);
    }
  
  	/* allocation - rigid bodies */
	if (n_r>0)
	{
		num_r=createMemMore(int,n_r);
		type_r=createMemMore(int,n_r);
		mat_r=createMemMore(int,n_r);
		m_r=createMemMore(double,n_r);
		I1_r=createMemMore(double,n_r);
		I2_r=createMemMore(double,n_r);
		I3_r=createMemMore(double,n_r);
		D1_r=createMemMore(double,n_r);
		D2_r=createMemMore(double,n_r);
		D3_r=createMemMore(double,n_r);
		COG_r=createMemMore(int,n_r);
		COG_r_type=createMemMore(int,n_r);
		N1_r=createMemMore(int,n_r);
		N1_r_type=createMemMore(int,n_r);
		N2_r=createMemMore(int,n_r);
		N2_r_type=createMemMore(int,n_r);
		N3_r=createMemMore(int,n_r);
		N3_r_type=createMemMore(int,n_r);
		x_r=createMemMore(double,3*n_r);
		x_r_0=createMemMore(double,3*n_r);
		x_r_o=createMemMore(double,3*n_r);
		//x_r_a=createMemMore(double,3*n_r);
		psi_r=createMemMore(double,3*n_r);
		psi_r_0=createMemMore(double,3*n_r);
		psi_r_o=createMemMore(double,3*n_r);
		//psi_r_a=createMemMore(double,3*n_r);
		v_r=createMemMore(double,3*n_r);
		v_r_o=createMemMore(double,3*n_r);
		//v_r_a=createMemMore(double,3*n_r);
		o_r=createMemMore(double,3*n_r);
		o_r_o=createMemMore(double,3*n_r);
		//o_r_a=createMemMore(double,3*n_r);
		a_r=createMemMore(double,3*n_r);
		a_r_o=createMemMore(double,3*n_r);
		alpha_r=createMemMore(double,3*n_r);
		alpha_r_o=createMemMore(double,3*n_r);
		f_r=createMemMore(double,3*n_r);
		M_r=createMemMore(double,3*n_r);
		type_a_0_r=createMemMore(int,n_r);
		frame_a_0_r=createMemMore(int,n_r);
		a_0_r=createMemMore(double,3*n_r);
		type_f_0_r=createMemMore(int,n_r);
		frame_f_0_r=createMemMore(int,n_r);
		f_0_r=createMemMore(double,3*n_r);
		M_0_r=createMemMore(double,3*n_r);
		//x1_r=createMemMore(double,3*n_r);
		//x1_r_o=createMemMore(double,3*n_r);
		//x2_r=createMemMore(double,3*n_r);
		//x2_r_o=createMemMore(double,3*n_r);
		//x3_r=createMemMore(double,3*n_r);
		//x3_r_o=createMemMore(double,3*n_r);
		u1_r=createMemMore(double,3*n_r);
		u2_r=createMemMore(double,3*n_r);
		u3_r=createMemMore(double,3*n_r);
		rb_s=createMemMore(int,n_r);
		rb_b=createMemMore(int,n_r);
		constrained_r=createMemMore(int,6*n_r);
		constrained_r_frame=createMemMore(int,n_r);

		/* setting number of rigid bodies vectors to zero */
		memset(o_r,0,3*n_r*sizeof(double));
		memset(psi_r,0,3*n_r*sizeof(double));
		memset(a_r,0,3*n_r*sizeof(double));
		memset(alpha_r,0,3*n_r*sizeof(double));
		memset(f_r,0,3*n_r*sizeof(double));
		memset(M_r,0,3*n_r*sizeof(double));
		memset(type_a_0_r,0,n_r*sizeof(int));
		memset(frame_a_0_r,0,n_r*sizeof(int));
		memset(a_0_r,0,3*n_r*sizeof(double));
		memset(type_f_0_r,0,n_r*sizeof(int));
		memset(frame_f_0_r,0,n_r*sizeof(int));
		memset(f_0_r,0,3*n_r*sizeof(double));
		memset(M_0_r,0,3*n_r*sizeof(double));
		memset(rb_s,0,n_r*sizeof(int));
		memset(rb_b,0,n_r*sizeof(int));
	}

  	/* check memory integrity */
  	if (mem_check)
    {
		checkMemoryIntegrity();
    }

  	fprintf(stdout,
		"\n"
		"domain and boundary reading ..."
		"\n");

  	sph=fopen(TMP,"r");
  	while ((ch=getc(sph))!=EOF)
    {
		if (ch!='$')
		{
			/* reading chars till space or end of line or end of file */
			count_char=0;

			while ((ch!=' ')&&(ch!='\n')&&(ch!=EOF))
			{
				keyword[count_char]=ch;
				ch=getc(sph);
				count_char++;
			}

			keyword[count_char]='\0';

			if (!strcmp(keyword,"FUNCT"))
			{
				fscanf(sph,
					"%lf %lf %lf %lf %lf %lf",
					&(aux1),
					&(aux2),
					&(fmx_u[i_u]),
					&(fmy_u[i_u]),
					&(fdx_u[i_u]),
					&(fdy_u[i_u]));
				num_u[i_u]=(int)aux1;
				fun_u[i_u]=(int)aux2;
				fun_loc[i_u]=count_fun_pairs;

				while ((ch=getc(sph))!='\n'); // reading until end of line

				for (i=0;i<fun_u[i_u];i++) // read all function data pairs
				{
					/* function data pair */
					fscanf(sph,"%lf %lf",
						&(fxi_u[i+fun_loc[i_u]]),
						&(fyi_u[i+fun_loc[i_u]]));

					if (i<fun_u[i_u]-1)
					{
						while ((ch=getc(sph))!='\n'); // reading until end of line
					}
				}
				count_fun_pairs=count_fun_pairs+fun_u[i_u];
				i_u++;
			}
	  		else if (!strcmp(keyword,"MATER"))
	    	{
	      		fscanf(sph,
		     		"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		     		&(aux1),
		     		&(aux2),
		     		&(aux3),
		     		&(rho_m[i_m]),
		     		&(mu_m[i_m]),
		     		&(T_m[i_m]),
		     		&(kappa_m[i_m]),
		     		&(gamma_m[i_m]),
		     		&(coef1_m[i_m]),
		     		&(coef2_m[i_m]),
		     		&(coef3_m[i_m]),
		     		&(coef4_m[i_m]));
				num_m[i_m]=(int)aux1;
				domain_m[i_m]=(int)aux2;
				type_m[i_m]=(int)aux3;

				if (type_m[i_m]==0) // rigid body
				{
					num_r[i_r]=i_r+1; // rigid body numbering from 1
					mat_r[i_r]=num_m[i_m]; // mat_r[i_r]=i_m; means direct renumbering
					if (rho_m[i_m]>0) // rigid body type 0
					{
						type_r[i_r]=0;
						// rho = rho_m[i_m];
						// G = mu_m[i_m];
						// K = T_m[i_m];
						// kappa_m[i_m];
						// gamma_m[i_m];
						m_r[i_r]=0.0; // mass calculated automatically
						I1_r[i_r]=0.0; // first principal moment of inertial calculated automatically
						I2_r[i_r]=0.0; // second principal moment of inertial calculated automatically
						I3_r[i_r]=0.0; // third principal moment of inertial calculated automatically
						D1_r[i_r]=0.0; // first deviation moment of inertial calculated automatically
						D2_r[i_r]=0.0; // second deviation moment of inertial calculated automatically
						D3_r[i_r]=0.0; // third deviation moment of inertial calculated automatically
					}
					else // rho_m[i_m] < 0 -> rigid body type 1
					{
						type_r[i_r]=1;
						//rho_m[i_m]=-rho_m[i_m]; // density for smoothing -> density as switch for rigid body type
						m_r[i_r]=mu_m[i_m]; // mass on next position
						I1_r[i_r]=T_m[i_m];
						I2_r[i_r]=kappa_m[i_m];
						I3_r[i_r]=gamma_m[i_m];
						D1_r[i_r]=0.0;
						D2_r[i_r]=0.0;
						D3_r[i_r]=0.0;
					}
					COG_r[i_r]=(int)coef1_m[i_m];
					N1_r[i_r]=(int)coef2_m[i_m];
					N2_r[i_r]=(int)coef3_m[i_m];
					N3_r[i_r]=(int)coef4_m[i_m];
					i_r++;
				}
				i_m++;
	    	}      
		  	else if (!strcmp(keyword,"CONTACT"))
		    {
	    		fscanf(sph,
					"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			     	&(aux1),
			     	&(aux2),
			     	&(aux3),
		    	 	&(aux4),
		     		&(ct_c[i_c]),
			     	&(klin_c[i_c]),
			     	&(knon_c[i_c]),
			     	&(kf_c[i_c]),
			     	&(kd_c[i_c]));
	      		num_c[i_c]=(int)aux1;
		      	mat_c[i_c]=(int)aux2;
		      	mat_c[i_c+n_c]=(int)aux3;
	    	  	sw_c[i_c]=(int)aux4;
	      		i_c++;
	    	}
		 	else if (!strcmp(keyword,"SNODE"))
		    {
				fscanf(sph,
					"%lf %lf %lf %lf %lf %lf",
			     	&(aux1),
			     	&(aux2),
		    	 	&(V_s[i_s]),
		     		&(x_s[i_s]),&(x_s[i_s+n_s]),&(x_s[i_s+2*n_s]));
		      	x_s_0[i_s]=x_s[i_s];
		      	x_s_0[i_s+n_s]=x_s[i_s+n_s];
	    	  	x_s_0[i_s+2*n_s]=x_s[i_s+2*n_s];
	      		num_s[i_s]=(int)aux1;
		      	mat_s[i_s]=(int)aux2;
		      	i_s++;
	    	}
	  		else if (!strcmp(keyword,"BNODE"))
		    {
				fscanf(sph,
			    	"%lf %lf %lf %lf",
		    	 	&(aux1),
		     		&(x_b[i_b]),&(x_b[i_b+n_b]),&(x_b[i_b+2*n_b]));
		      	num_b[i_b]=(int)aux1;
		      	x_b_0[i_b]=x_b[i_b];
	    	  	x_b_0[i_b+n_b]=x_b[i_b+n_b];
	      		x_b_0[i_b+2*n_b]=x_b[i_b+2*n_b];
		      	i_b++;
		    }
		  	else if (!strcmp(keyword,"BELEM"))
	    	{
	      		fscanf(sph,
			     	"%lf %lf %lf %lf %lf %lf",
				    &(aux1),
				    &(aux2),
		    	 	&(aux3),
		     		&(aux4),
			     	&(aux5),
			     	&(aux6));
	    	  	num_e[i_e]=(int)aux1;
	      		mat_e[i_e]=(int)aux2;
		      	nod_e[i_e]=(int)aux3;
		      	nod_e[i_e+n_e]=(int)aux4;
	    	  	nod_e[i_e+2*n_e]=(int)aux5;
	      		nod_e[i_e+3*n_e]=(int)aux6;
		      	i_e++;
	    	}
		  	else if (!strcmp(keyword,"INPRE"))
		    {
				fscanf(sph,"%lf %lf",
					&(aux1),
			     	&(p_0[i_p]));
	    	  	inpre[i_p]=(int)aux1;
	      		i_p++;
		    }
		  	else if (!strcmp(keyword,"INVEL"))
		    {
				fscanf(sph,"%lf %lf %lf %lf %lf %lf %lf %lf",
					&(aux1),
			     	&(v_0[i_v]),&(v_0[i_v+n_v]),&(v_0[i_v+2*n_v]),
			     	&(omega_0[i_v]),&(omega_0[i_v+n_v]),&(omega_0[i_v+2*n_v]),&(aux2));
	    	  	invel[i_v]=(int)aux1;
				invel_frame[i_v]=(int)aux2;
	      		if (sqrt(sqr(v_0[i_v])+sqr(v_0[i_v+n_v])+sqr(v_0[i_v+2*n_v]))>v_max)
				{
					v_max=sqrt(sqr(v_0[i_v])+sqr(v_0[i_v+n_v])+sqr(v_0[i_v+2*n_v]));
				}
				i_v++;
	    	}
		  	else if (!strcmp(keyword,"ACFLD"))
		    {
				fscanf(sph,"%lf %lf %lf %lf %lf %lf",
		    	 	&(aux1),
		    	 	&(aux2),
		     		&(a_0[i_a]),&(a_0[i_a+n_a]),&(a_0[i_a+2*n_a]),
					&(aux3));
		      	acfld[i_a]=(int)aux1;
		      	acfld_type[i_a]=(int)aux2;
		      	acfld_frame[i_a]=(int)aux3;
		      	i_a++;
	    	}
		  	else if (!strcmp(keyword,"FORCE"))
		    {
				fscanf(sph,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&(aux1),
					&(aux2),
		     		&(f_0[i_f]),&(f_0[i_f+n_f]),&(f_0[i_f+2*n_f]),
		     		&(M_0[i_f]),&(M_0[i_f+n_f]),&(M_0[i_f+2*n_f]),&(aux3));
		      	force[i_f]=(int)aux1;
		      	force_type[i_f]=(int)aux2;
				force_frame[i_f]=(int)aux3;
		      	i_f++;
	    	}
		  	else if (!strcmp(keyword,"NDAMP"))
		    {
				fscanf(sph,"%lf %lf %lf %lf",
					&(aux1),
			     	&(d_0[i_d]),
		    	 	&(d_0[i_d+n_d]),
		     		&(d_0[i_d+2*n_d]));
		      	ndamp[i_d]=(int)aux1;
		      	i_d++;
	    	}
		  	else if (!strcmp(keyword,"BOUNC"))
		    {
				fscanf(sph,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&(aux1),
			     	&(aux2),
			     	&(aux3),
			     	&(aux4),
			     	&(aux5),
			     	&(aux6),
			     	&(aux7),
			     	&(aux8),
		    	 	&(aux9));
		      	bounc[i_o]=(int)aux1;
				bounc_type[i_o]=(int)aux2;
		      	o_0[i_o]=(int)aux3;
	    	  	o_0[i_o+n_o]=(int)aux4;
	      		o_0[i_o+2*n_o]=(int)aux5;
	      		o_0[i_o+3*n_o]=(int)aux6;
	      		o_0[i_o+4*n_o]=(int)aux7;
	      		o_0[i_o+5*n_o]=(int)aux8;
				bounc_frame[i_o]=(int)aux9;
		      	i_o++;
		    }
	  		else if (!strcmp(keyword,"END")) break;
		}
		while ((ch=getc(sph))!='\n'); // reading until end of line
	}
  	fclose(sph);
  	remove(TMP);

	/* dimension reduction 
	-> free dimension input set to zero */
	/* particles position and rotation */
	
	/* particles coordinates */
	for (i=0;i<n_s;i++)
	{
		if (dim<3) // 2D
		{
			// rotation psi(0) = 0
			x_s[i+2*n_s]=0.0; // coordinate z
		}
		if (dim<2) // 2D -> 1D
		{
			// rotations psi(0) = 0
			x_s[i+n_s]=0.0; // coordinate y
		}
	}

	/* nodal coordinates */
	for (i=0;i<n_b;i++)
	{
		if (dim<3) // 2D
		{
			x_b[i+2*n_b]=0.0; // coordinate z
		}
		if (dim<2) // 2D -> 1D
		{
			x_b[i+n_b]=0.0; // coordinate y
		}
	}

	/* initial velocities and rotational velocities */
	for (i=0;i<n_v;i++)
	{
		if (dim<3) // 2D
		{
			v_0[i+2*n_v]=0.0; // translational velocity z
			omega_0[i+2*n_v]=0.0; // rotational velocity z
			omega_0[i+n_v]=0.0; // rotational velocity y
		}
		if (dim<2) // 2D -> 1D 
		{
			v_0[i+n_v]=0.0; // translational velocity y
			omega_0[i]=0.0; // rotational velocity x -> no rotation
		}
	}

	/* boundary acceleration fields */
	for (i=0;i<n_a;i++)
	{
		if (dim<3) // 2D
		{
			a_0[i+2*n_a]=0; // acceleration z
		}
		if (dim<2) // 2D -> 1D
		{
			a_0[i+n_a]=0; // acceleration y
		}
	}

	/* boundary forces and moments */
	for (i=0;i<n_f;i++)
	{
		if (dim<3) // 2D
		{
			f_0[i+2*n_f]=0; // force z
			M_0[i+2*n_f]=0; // moment z
			M_0[i+n_f]=0; // moment y
		}
		if (dim<2) // 2D -> 1D
		{
			f_0[i+n_f]=0; // force y
			M_0[i]=0; // moment x -> no moment in 1D
		}
	}

	/* nodal dampings */
	for (i=0;i<n_d;i++)
	{
		if (dim<3) // 2D
		{
			d_0[i+2*n_d]=0.0; // damping z
		}
		if (dim<2) // 2D -> 1D
		{
			d_0[i+n_d]=0.0; // damping y
		}
	}

	/* particles boundary conditions */
	for (i=0;i<n_o;i++)
	{
		if (dim<3) // 2D 
		{
			o_0[i+2*n_o]=0; // translation z
			o_0[i+5*n_o]=0; // rotation z
			o_0[i+4*n_o]=0; // rotation y
		}
		if (dim<2) // 2D -> 1D 
		{
			o_0[i+n_o]=0; // translation y
			o_0[i+3*n_o]=0; // rotation x -> no rotation
		}
	}

 	/* test on particles */
	for (i=0;i<n_s;i++)
    {
		/* test on zero or negative particle numbering */
		if (num_s[i]<=0)
		{
			fprintf(stdout,
				"ERROR: particle %d must be a positive integer!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_s[i]);
	  		return(0);
		}

		/* maximum node number */
		//if (num_s[i]>n_max) n_max=num_s[i];

		/* test on zero or negative volume */
		if (V_s[i]<=0.0)
		{
			fprintf(stdout,
				"ERROR: particle %d has zero or negative volume!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_s[i]);
	  		return(0);
		}

      	/* test on coincident particles */
      	for (j=i+1;j<n_s;j++)
		{
			if ((x_s[i]==x_s[j])&&(x_s[i+n_s]==x_s[j+n_s])&&(x_s[i+2*n_s]==x_s[j+2*n_s]))
			{
				fprintf(stdout,
		      		"ERROR: particle %d has coincident coordinates with partcile %d!"
					"\n"
		      		"\n",
		      		num_s[i],num_s[j]);
	      		return(0);
	    	}
		}
      
      	/* test on double snode numbering */
      	for (j=i+1;j<n_s;j++)
		{
			if (num_s[i]==num_s[j])
			{
				fprintf(stdout,
		      		"ERROR: particle %d occurs more than once!"
					"\n"
		      		"\n",
		      		num_s[i]);
	      		return(0);
	    	}
		}
      
    	/* test on double snode/bnode numbering */
    	for (j=0;j<n_b;j++)
		{
			if (num_s[i]==num_b[j])
	    	{
				fprintf(stdout,
		    		"ERROR: particle %d is the same as node %d!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		num_s[i],num_s[i]);
	      		return(0);
	    	}
		}
    }

  	/* test on nodes */
	for (i=0;i<n_b;i++)
    {
		/* test on zero or negative nodal numbering */
		if (num_b[i]<=0)
		{
			fprintf(stdout,
				"ERROR: node %d must be a positive integer!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_b[i]);
	  		return(0);
		}

		/* maximum node number */
		//if (num_b[i]>n_max) n_max=num_b[i];

		/* test on double bnode numbering */
      	for (j=i+1;j<n_b;j++)
		{
	  		if (num_b[i]==num_b[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: node %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		num_b[i]);
	      		return(0);
		    }
		}	      
    }

  	/* test on double element numbering */
  	for (i=0;i<n_e;i++)
    {
		/* test on zero or negative element numbering */
		if (num_e[i]<=0)
		{
			fprintf(stdout,
				"ERROR: element number %d must be a positive integer!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_e[i]);
	  		return(0);
		}

		for (j=i+1;j<n_e;j++)
		{
	  		if (num_e[i]==num_e[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: element %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		num_e[i]);
	      		return(0);
	    	}
		}
    }
  
  	/* rejecting free nodes and their properties */
  	if (!n_b) inode=1;
  	while (!inode)
    {
		/* checking for free node */
		for (i=0;i<n_b;i++)
		{
	  		inode=0;
	  		for (j=0;j<n_e;j++)
	    	{
	      		if ((num_b[i]==nod_e[j])||(num_b[i]==nod_e[j+n_e])||
					(num_b[i]==nod_e[j+2*n_e])||(num_b[i]==nod_e[j+3*n_e]))
				{
		  			inode=1;
		  			break;
				}
	    	}

	  		if (!inode)
	    	{
				/* free node can be COG or N1 or N2 or N3 */
				for (j=0;j<n_r;j++)
				{
					if ((COG_r[j]==num_b[i])||(N1_r[j]==num_b[i])||(N2_r[j]==num_b[i])||(N3_r[j]==num_b[i]))
					{
						constrained_b_rb[i]=num_r[j]; // num_r = rigid body number from 1)
						mat_b[i]=-1; // free node does not have material number
						inode=1;
						break;
					}
				}

				if (!inode)
				{
					fprintf(stdout,
						"node %d is free -> rejected!\n",
						num_b[i]);
			
					/* rejecting free initial velocities */
					for (j=0;j<n_v;j++)
					{
						if (num_b[i]==invel[j])
						{
							for (k=j;k<n_v-1;k++)
							{
								invel[j]=invel[j+1];
								v_0[j]=v_0[j+1];
								v_0[j+n_v-1]=v_0[j+1+n_v];
								v_0[j+2*n_v-2]=v_0[j+1+2*n_v];
								omega_0[j]=omega_0[j+1];
								omega_0[j+n_v-1]=omega_0[j+1+n_v];
								omega_0[j+2*n_v-2]=omega_0[j+1+2*n_v];
								invel_frame[j]=invel_frame[j+1];
							}
							n_v--;
							fprintf(stdout,
								" -> initial velocity to free node %d rejected!\n",
								num_b[i]);
							break;
						}
					}

					/* rejecting free initial acceleration fields */
					for (j=0;j<n_a;j++)
					{
						if (num_b[i]==acfld[j])
						{
							for (k=j;k<n_a-1;k++)
							{
								acfld[j]=acfld[j+1];
								acfld_type[j]=acfld_type[j+1];
								acfld_frame[j]=acfld_frame[j+1];
								a_0[j]=a_0[j+1];
								a_0[j+n_a-1]=a_0[j+1+n_a];
								a_0[j+2*n_a-2]=a_0[j+1+2*n_a];
							}
							n_a--;
							fprintf(stdout,
								" -> acceleration field to free node %d rejected!\n",
								num_b[i]);
							break;
						}
					}

					/* rejecting free nodal forces */
					for (j=0;j<n_f;j++)
					{
						if (num_b[i]==force[j])
						{
							for (k=j;k<n_f-1;k++)
							{
								force[j]=force[j+1];
								force_type[j]=force_type[j+1];
								f_0[j]=f_0[j+1];
								f_0[j+n_f-1]=f_0[j+1+n_f];
								f_0[j+2*n_f-2]=f_0[j+1+2*n_f];
								M_0[j]=M_0[j+1];
								M_0[j+n_f-1]=M_0[j+1+n_f];
								M_0[j+2*n_f-2]=M_0[j+1+2*n_f];
								force_frame[j]=force_frame[j+1];
							}
							n_f--;
							fprintf(stdout,
								" -> nodal force to free node %d rejected!\n",
								num_b[i]);
							break;
						}
					}

					/* rejecting free damping */
					for (j=0;j<n_d;j++)
					{
						if (num_b[i]==ndamp[j])
						{
							for (k=j;k<n_d-1;k++)
							{
								ndamp[j]=ndamp[j+1];
								d_0[j]=d_0[j+1];
							}
							n_d--;
							fprintf(stdout,
								" -> nodal damping to free node %d rejected!\n",
								num_b[i]);
							break;
						}
					}

					/* rejecting free boundary conditions */
					for (j=0;j<n_o;j++)
					{
						if (num_b[i]==bounc[j])
						{
							for (k=j;k<n_o-1;k++)
							{
								bounc[j]=bounc[j+1];
								bounc_type[j]=bounc_type[j+1];
								o_0[j]=o_0[j+1];
								o_0[j+n_o-1]=o_0[j+1+n_o];
								o_0[j+2*n_o-2]=o_0[j+1+2*n_o];
								o_0[j+3*n_o-2]=o_0[j+1+3*n_o];
								o_0[j+4*n_o-2]=o_0[j+1+4*n_o];
								o_0[j+5*n_o-2]=o_0[j+1+5*n_o];
								bounc_frame[j]=bounc_frame[j+1];
							}
							n_o--;
							fprintf(stdout,
								" -> boundary condition on free node %d rejected!\n",
								num_b[i]);
							break;
						}
					}

					/* rejecting free node */
					for (j=i;j<n_b-1;j++)
					{
						num_b[j]=num_b[j+1];
						x_b[j]=x_b[j+1];
						x_b_0[j]=x_b_0[j+1];
						x_b[j+n_b-1]=x_b[j+1+n_b];
						x_b_0[j+n_b-1]=x_b_0[j+1+n_b];
						x_b[j+2*n_b-2]=x_b[j+1+2*n_b];
						x_b_0[j+2*n_b-2]=x_b_0[j+1+2*n_b];
					}
					n_b--;
					break;
				}
		    }
		}
    }

  	/* test on nodes - second test after rejecting free nodes */
  	if ((n_s==0)&&(n_b==0))
    {
      	fprintf(stdout,
	      	"ERROR: there are no partciles and nodes in the input file!"
	      	"\n"
			"\n"
			"error termination"
			"\n"
	      	"\n");
      	return(0);
    }

  	/* dimension of boundary equation of motion matrices */
	if (n_e>0)
	{
		n_n=dim*n_b;

		M=createMemMore(double,n_n*n_n);
		B=createMemMore(double,n_n*n_n);
		K=createMemMore(double,n_n*n_n);
		f=createMemMore(double,n_n);

		/* setting matrices to zeros */
      	memset(M,0,n_n*n_n*sizeof(double));
      	memset(B,0,n_n*n_n*sizeof(double));
      	memset(K,0,n_n*n_n*sizeof(double));
      	memset(f,0,n_n*sizeof(double));
    }
  
  	/* printing all info */
  	fprintf(stdout,
	  	"\n"
	  	"file name"
	  	"\n"
	  	" -> %s"
	  	"\n"
	  	" -> dimension = %d"
	  	"\n"
	  	" -> termination time = %f"
	  	"\n"
	  	"\n",
	  	output,dim,t_max);

  	fprintf(stdout,
	  	"time controls"
	  	"\n"
	  	" -> saving time interval = %s"
	  	"\n"
	  	" -> initial time step = %f"
	  	"\n"
	  	" -> maximum time step = %s"
	  	"\n"
	  	" -> Courant's number = %f"
	  	"\n"
	  	" -> stabile time step factor = %f"
	  	"\n"
	  	"\n",
	  	dt_save_type,dt_init,dt_max_type,cour,kstab);

  	fprintf(stdout,
	  	"global acceleration field"
	  	"\n"
	  	" -> global acceleration field ax = %f"
	  	"\n"
	  	" -> global acceleration field ay = %f"
	  	"\n"
	  	" -> global acceleration field az = %f"
	  	"\n"
	  	"\n",
	  	ax,ay,az);
  
  	fprintf(stdout,
	  	"optimization control"
	  	"\n"
	  	" -> data check = %s"
	  	"\n"
	  	" -> data print = %s"
	  	"\n"
	  	" -> cycle print = %d"
	  	"\n"
	  	" -> nearest neighbour search = %s"
	  	"\n"
	  	" -> nns radius multiplier = %s"
	  	"\n"
	  	" -> contact frequency = %d"
	  	"\n"
	  	" -> integration scheme = %s"
	  	"\n"
	  	" -> memory check active = %s"
	  	"\n"
	  	"\n",	  
	  	data_check_type,data_print_type,cycle_print,
		nns_type,opt_type,
		cycle_contact,
		int_scheme,mem_check_type);

  	fprintf(stdout,
	  	"SPH"
	  	"\n"
	  	" -> viscosity type = %s"
	  	"\n"
	  	"\n"
	  	" -> artificial viscosity coeffitient alpha = %f"
	  	"\n"
	  	" -> artificial viscosity coeffitient beta = %f"
	  	"\n"
	  	" -> artificial viscosity coeffitient eta = %f"
	  	"\n"
	  	"\n"
	  	" -> artificial stress coeffitient zeta = %f"
	  	"\n"
	  	" -> artificial stress coeffitient n = %f"
	  	"\n"
	  	" -> artificial stress coeffitient theta = %f"
	  	"\n"
	  	"\n"
	  	" -> XSPH = %s"
	  	"\n"
	  	" -> XSPH epsilon = %f"
	  	"\n"
	  	"\n",
	  	is_type,alpha,beta,eta,zeta,nas,theta,xsph,xeps);
  
  	fprintf(stdout,
	  	"FEM"
	  	"\n"
	  	" -> damping coeffitient c0 = %f"
	  	"\n"
	  	" -> damping coeffitient c1 = %f"
	  	"\n"	  
	  	"\n",
	  	c0,c1);
  
  	fprintf(stdout,
	  	"save"
	  	"\n"
	  	" -> time step = %s"
	  	"\n"
	  	" -> particles kinetic energy = %s"
	  	"\n"
	  	" -> particles internal energy = %s"
	  	"\n"
	  	" -> particles potential energy = %s"
	  	"\n"
	  	" -> nodes kinetic energy = %s"
	  	"\n"
	  	" -> nodes disipation energy = %s"
	  	"\n"
	  	" -> nodes deformation energy = %s"
	  	"\n"
	  	" -> nodes potential energy = %s"
	  	"\n"
	  	" -> total energy = %s"
	  	"\n"
	  	" -> particles velocities = %s"
	  	"\n"
	  	" -> particles velocity corrections = %s"
	  	"\n"
	  	" -> particles accelerations = %s"
	  	"\n"
	  	" -> particles forces = %s"
	  	"\n"
	  	" -> particles densities = %s"
	  	"\n"
	  	" -> particles changes of densities = %s"
	  	"\n"
	  	" -> particles internal energies = %s"
	  	"\n"
	  	" -> particles changes of internal energies = %s"
	  	"\n"
	  	" -> particles pressures = %s"
	  	"\n"
	  	" -> particles sound speeds = %s"
	  	"\n"
	  	" -> particles smoothing lenghts = %s"
	  	"\n"
	  	" -> particles deviatoric stress tensors = %s"
	  	"\n"
	  	" -> particles change of deviatoric stress tensors = %s"
	  	"\n"
	  	" -> particles deformation tensors = %s"
	  	"\n"
	  	" -> particles rate of deformation tensors = %s"
	  	"\n"
	  	" -> particles rotation tensors = %s"
	  	"\n"
	  	" -> nodal velocities = %s"
	  	"\n"
	  	" -> nodal accelerations = %s"
	  	"\n"
	  	" -> nodal forces = %s"
	  	"\n"
	  	" -> rigid bodies rotation = %s"
	  	"\n"
	  	" -> rigid bodies translational velocity = %s"
	  	"\n"
	  	" -> rigid bodies rotational velocity = %s"
	  	"\n"
	  	" -> rigid bodies translational acceleration = %s"
	  	"\n"
	  	" -> rigid bodies rotational acceleration = %s"
	  	"\n"
	  	" -> rigid bodies forces = %s"
	  	"\n"
	  	" -> rigid bodies moments = %s"
	  	"\n"
	  	" -> contact forces = %s"
	  	"\n"
	  	" -> rigid bodies centre of gravity position= %s"
	  	"\n"
	  	"\n",
	  	save_dt_type,
	  	save_kine_s_type,save_inne_s_type,save_pote_s_type,
		save_kine_b_type,save_disi_b_type,save_defo_b_type,save_pote_b_type,
	  	save_tote_type,
	  	save_v_s_type,save_dv_s_type,save_a_s_type,save_f_s_type,
	  	save_rho_s_type,save_drhodt_s_type,save_u_s_type,save_dudt_s_type,
		save_p_s_type,save_c_s_type,save_h_s_type,
	  	save_S_s_type,save_dSdt_s_type,save_e_s_type,save_dedt_s_type,save_O_s_type,
	  	save_v_b_type,save_a_b_type,save_f_b_type,
		save_x_r_type,save_psi_r_type,save_v_r_type,save_o_r_type,save_a_r_type,
		save_alpha_r_type,save_f_r_type,save_M_r_type,
	  	save_f_c_type);

	/* test on is */
	if ((n_s>1)&&(is!=0)&&(is!=1)&&(is!=2)&&(is!=3)&&(is!=4)&&(is!=5)&&(is!=6))
	{
		fprintf(stdout,
			"ERROR: unknown viscosity type!"
			"\n"
			"\n"
			"error termination"
			"\n"
			"\n");
		return(0);
	}
	else
	{
		switch (is)
		{
			case 0:
			{
				strcpy(is_type,"0 (artificial viscosity only)");
				break;
			}
			case 1:
			{
				strcpy(is_type,"1 (second order viscous term)");
				break;
			}
			case 2:
			{
				strcpy(is_type,"2 (Monaghan-Cleary-Gingold, 2006)");
				break;
			}
			case 3:
			{
				strcpy(is_type,"3 (Morris et al., 1997)");
				break;
			}
			case 4:
			{
				strcpy(is_type,"4 (Takeda et al., 1994)");
				break;
			}
			case 5:
			{
				strcpy(is_type,"5 (Onderik et al., 2007)");
				break;
			}
			case 6:
			{
				strcpy(is_type,"6 (Monaghan and Gingold, 1983)");
				break;
			}
		}
	}

	/* test on data */
	if (data_check)
	{
		/* test on cour */
		if ((n_s>0)&&(cour<0.0))
		{
			fprintf(stdout,
				"ERROR: Courant number must be a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			return(0);
		}

		/* test on kstab */
		if ((n_b>0)&&(kstab<0.0))
		{
			fprintf(stdout,
				"ERROR: Stable time step coefficient must be a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			return(0);
		}

		/* test on alpha */
		if ((n_s>1)&&(is==0)&&(alpha<0.0))
		{
			fprintf(stdout,
				"ERROR: coeffitient alpha must be 0 or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			return(0);
		}

		/* test on beta */
		if ((n_s>1)&&(is==0)&&(beta<0.0))
		{
			fprintf(stdout,
				"ERROR: coeffitient beta must be 0 or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			return(0);
		}

		/* test on eta */
		if ((n_s>1)&&(is==0)&&(eta<=0.0))
		{
			fprintf(stdout,
				"ERROR: coeffitient eta must be a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			return(0);
		}

		/* test on c0 */
		if ((n_e>0)&&(c0<0.0))
		{
			fprintf(stdout,
				"ERROR: coeffitient c0 must be 0 or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			return(0);
		}

		/* test on c1 */
		if ((n_e>0)&&(c1<0.0))
		{
			fprintf(stdout,
				"ERROR: coeffitient c1 must be 0 or a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n");
			return(0);
		}
	}
  
  	/* functions */
  	for (i=0;i<n_u;i++)
    {
		fprintf(stdout,
			"function %d"
			"\n"
			" -> number of data pairs = %d"
			"\n"
			" -> function abscissa multiplier = %f"
			"\n"
			" -> function ordinate multiplier = %f"
			"\n"
			" -> function abscissa shift = %f"
			"\n"
			" -> function ordinate shift  = %f"
			"\n",
			num_u[i],fun_u[i],fmx_u[i],fmy_u[i],fdx_u[i],fdy_u[i]);

      	/* print data pairs */
		for (j=0;j<fun_u[i];j++)
		{
			fprintf(stdout,
				"     -> function data pair %d (%f,%f)",
				j,fxi_u[j+fun_loc[i]],fyi_u[j+fun_loc[i]]);

				fxi_u[j+fun_loc[i]]=fmx_u[i]*fxi_u[j+fun_loc[i]]+fdx_u[i];
				fyi_u[j+fun_loc[i]]=fmy_u[i]*fyi_u[j+fun_loc[i]]+fdy_u[i];

			fprintf(stdout,
				" -> (%f,%f)"
				"\n",
				fxi_u[j+fun_loc[i]],fyi_u[j+fun_loc[i]]);
		}
		fprintf(stdout,"\n");

		/* test on zero or negative function numbering */
		if (num_u[i]<=0)
		{
			fprintf(stdout,
				"ERROR: function number %d must be a positive integer!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_u[i]);
	  		return(0);
		}

		/* test on double numbering */      
      	for (j=i+1;j<n_u;j++)
	  	{
			if (num_u[i]==num_u[j])
			{
	      		fprintf(stdout,
		      		"ERROR: function %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		num_u[i]);
	      		return(0);
		  	}
		}

	  	/* test on at least 2 data pairs */
		if (fun_u[i]<2)
		{
			fprintf(stdout,
				"ERROR: function %d must have at least 2 data pairs!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				num_u[i]);
			return(0);
		}
			
	  	/* test on increasing abscissa */
		for (j=1;j<fun_u[i];j++)
		{
			if (fxi_u[j+fun_loc[i]]<=fxi_u[j+fun_loc[i]-1])
			{
				fprintf(stdout,
					"ERROR: function %d must have increasing abscissa!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_u[i]);
				return(0);
			}
		}
    }

  	/* materials */
  	for (i=0;i<n_m;i++)
    {
		/* test on zero or negative material numbering */
		if (num_m[i]<=0)
		{
			fprintf(stdout,
				"ERROR: material number %d must be a positive integer!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_m[i]);
	  		return(0);
		}

		/* test on zero or negative domain numbering */
		if (domain_m[i]<=0)
		{
			fprintf(stdout,
				"ERROR: material number %d must be a positive integer!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		domain_m[i]);
	  		return(0);
		}

      	/* test on type */
      	if ((type_m[i]!=0)&&
			(type_m[i]!=1)&&(type_m[i]!=2)&&(type_m[i]!=3)&&
			(type_m[i]!=4)&&(type_m[i]!=5)&&(type_m[i]!=6)&&
			(type_m[i]!=7)&& (type_m[i]!=8)&&(type_m[i]!=9)&&
			(type_m[i]!=12)&&(type_m[i]!=13)&&
			(type_m[i]!=23)&&(type_m[i]!=27))
	  	{
		  	fprintf(stdout,
		  		"ERROR: unknown material type %d of material %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		type_m[i],num_m[i]);
	  		return(0);
		}
      
		if ((type_m[i]==0)&&(rho_m[i]>0))
	  	{
			fprintf(stdout,
			  	"material number %d"
			  	"\n"
			  	" -> domain = %d"
			  	"\n"
			  	" -> type: rigid body 0"
		  		"\n"
			  	" -> density = %f"
			  	"\n"
		  		" -> shear modulus = %f"
		  		"\n"
		  		" -> bulk modulus = %f"
		  		"\n"
			  	" -> first free coefficient = %f"
		  		"\n"
			  	" -> second free coefficient = %f"
			  	"\n"
			  	" -> centre of gravity node = %d"
			  	"\n"
			  	" -> first principal axis node = %d"
			  	"\n"
			  	" -> second principal axis node = %d"
			  	"\n"
			  	" -> third principal axis node = %d"
			  	"\n"
			  	"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				(int)coef1_m[i],(int)coef2_m[i],(int)coef3_m[i],(int)coef4_m[i]);
		}
		if ((type_m[i]==0)&&(rho_m[i]<0))
	  	{
			fprintf(stdout,
			  	"material number %d"
			  	"\n"
			  	" -> domain = %d"
			  	"\n"
			  	" -> type: rigid body 1"
		  		"\n"
			  	" -> density = %f"
			  	"\n"
		  		" -> mass = %f"
		  		"\n"
		  		" -> first principal moment of inertia = %f"
		  		"\n"
			  	" -> second principal moment of inertia = %f"
		  		"\n"
			  	" -> third principal moment of inertia = %f"
			  	"\n"
			  	" -> centre of gravity node = %d"
			  	"\n"
			  	" -> first principal axis node = %d"
			  	"\n"
			  	" -> second principal axis node = %d"
			  	"\n"
			  	" -> third principal axis node = %d"
			  	"\n"
			  	"\n",
		  		num_m[i],domain_m[i],-rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				(int)coef1_m[i],(int)coef2_m[i],(int)coef3_m[i],(int)coef4_m[i]);
		}
		if (type_m[i]==1)
	  	{
			fprintf(stdout,
			  	"material number %d"
			  	"\n"
			  	" -> domain = %d"
			  	"\n"
			  	" -> type: gas"
		  		"\n"
			  	" -> density = %f"
			  	"\n"
			  	" -> viscosity = %f"
		  		"\n"
			  	" -> absolute temperature = %f"
			  	"\n"
			  	" -> adiabatic coeffitient = %f"
			  	"\n"
			  	" -> constant volume heat coefficient = %f"
			  	"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
			  	"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==2)
	  	{
			fprintf(stdout,
			  	"material number %d"
			  	"\n"
			  	" -> domain = %d"
			  	"\n"
			  	" -> type: liquid"
		  		"\n"
			  	" -> density = %f"
			  	"\n"
			  	" -> viscosity = %f"
		  		"\n"
			  	" -> bulk modulus = %f"
			  	"\n"
			  	" -> initial pressure = %f"
		  		"\n"
			  	" -> power in equation of state = %f"
			  	"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
			  	"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==3)
	  	{
			fprintf(stdout,
			  	"material number %d"
			  	"\n"
			  	" -> domain = %d"
			  	"\n"
			  	" -> type: SPH"
			  	"\n"
			  	" -> density = %f"
			  	"\n"
		  		" -> viscosity = %f"
			  	"\n"
			  	" -> first free coefficient = %f"
			  	"\n"
		  		" -> initial pressure = %f"
			  	"\n"
			  	" -> power in equation of state = %f"
			  	"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
			  	"\n",
			  	num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if ((type_m[i]>3)&&(type_m[i]<7))
	  	{
			fprintf(stdout,
			  	"material number %d"
			  	"\n"
			  	" -> domain = %d"
			  	"\n"
			  	" -> type: FEM"
			  	"\n"
			  	" -> density = %f"
			  	"\n"
			  	" -> Poisson's ratio = %f"
			  	"\n"
			  	" -> Young's modulus = %f"
			  	"\n"
			  	" -> thickness / area = %f"
			  	"\n"
			  	" -> damping = %f"
			  	"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
			  	"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==7)
	  	{
			fprintf(stdout,
			  	"material number %d"
		  		"\n"
			  	" -> domain = %d"
			  	"\n"
		  		" -> type: SPS"
		  		"\n"
		  		" -> density = %f"
		  		"\n"
		  		" -> shear modulus = %f"
		  		"\n"
		  		" -> bulk modulus = %f"
		  		"\n"
		  		" -> initial pressure = %f"
		  		"\n"
		  		" -> power in equation of state = %f"
		  		"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
		  		"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==8)
	  	{
			fprintf(stdout,
		  		"material number %d"
		  		"\n"
			  	" -> domain = %d"
			  	"\n"
		  		" -> type: Neo-hookean solid"
		  		"\n"
				" -> density = %f"
		  		"\n"
		  		" -> shear modulus = %f"
		  		"\n"
		  		" -> bulk modulus = %f"
		  		"\n"
		  		" -> initial pressure = %f"
		  		"\n"
		  		" -> power in equation of state = %f"
		  		"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
		  		"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==9)
	  	{
			fprintf(stdout,
		  		"material number %d"
		  		"\n"
			  	" -> domain = %d"
			  	"\n"
		  		" -> type: user-defined solid"
		  		"\n"
		  		" -> density = %f"
		  		"\n"
		  		" -> first free coefficient = %f"
		  		"\n"
		  		" -> second free coefficient = %f"
		  		"\n"
		  		" -> third free coefficient = %f"
		  		"\n"
		  		" -> fourth free coefficient = %f"
				"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
		  		"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==12)
	  	{
			fprintf(stdout,
		  		"material number %d"
		  		"\n"
			  	" -> domain = %d"
			  	"\n"
		  		" -> type: liquid with tension"
		  		"\n"
		  		" -> density = %f"
		  		"\n"
		  		" -> viscosity = %f"
		  		"\n"
		  		" -> bulk modulus = %f"
		  		"\n"
		  		" -> tension modulus = %f"
		  		"\n"
		  		" -> power in equation of state = %f"
		  		"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
		  		"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==13)
	  	{
			fprintf(stdout,
		  		"material number %d"
		  		"\n"
			  	" -> domain = %d"
			  	"\n"
		  		" -> type: SPH with tension"
		  		"\n"
		  		" -> density = %f"
		  		"\n"
		  		" -> viscosity = %f"
		  		"\n"
		  		" -> first free coefficient = %f"
		  		"\n"
		  		" -> tension modulus = %f"
		  		"\n"
		  		" -> power in equation of state = %f"
		  		"\n"
			  	" -> first auxiliary coefficient = %f"
			  	"\n"
			  	" -> second auxiliary coefficient = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
		  		"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==23)
	  	{
			fprintf(stdout,
		  		"material number %d"
		  		"\n"
			  	" -> domain = %d"
			  	"\n"
		  		" -> type: Mie-Grueneisen EOS for SPH"
		  		"\n"
		  		" -> density = %f"
		  		"\n"
		  		" -> viscosity = %f"
		  		"\n"
		  		" -> bulk modulus = %f"
		  		"\n"
		  		" -> initial pressure = %f"
		  		"\n"
		  		" -> first free coefficient = %f"
		  		"\n"
			  	" -> linear coefficient of Hugoniot = %f"
			  	"\n"
			  	" -> Grueneisen constant = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
		  		"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}
      	if (type_m[i]==27)
	  	{
			fprintf(stdout,
		  		"material number %d"
		  		"\n"
			  	" -> domain = %d"
			  	"\n"
		  		" -> type: Mie-Grueneisen EOS for SPS"
		  		"\n"
		  		" -> density = %f"
		  		"\n"
		  		" -> shear modulus = %f"
		  		"\n"
		  		" -> bulk modulus = %f"
		  		"\n"
		  		" -> initial pressure = %f"
		  		"\n"
		  		" -> first free coefficient = %f"
		  		"\n"
			  	" -> linear coefficient of Hugoniot = %f"
			  	"\n"
			  	" -> Grueneisen constant = %f"
			  	"\n"
			  	" -> third auxiliary coefficient = %f"
			  	"\n"
			  	" -> fourth auxiliary coefficient = %f"
			  	"\n"
		  		"\n",
		  		num_m[i],domain_m[i],rho_m[i],
				mu_m[i],T_m[i],kappa_m[i],gamma_m[i],
				coef1_m[i],coef2_m[i],coef3_m[i],coef4_m[i]);
		}

      	/* test on double numbering */      
      	for (j=i+1;j<n_m;j++)
	  	{
			if (num_m[i]==num_m[j])
			{
	      		fprintf(stdout,
		      		"ERROR: material number %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		num_m[i]);
	      		return(0);
		  	}
		}

		/* test on particles or elements in rigid body */
		if (type_m[i]==0)
		{
			for (j=0;j<n_s;j++)
			{
				if (mat_s[j]==num_m[i])
				{
					imate=1;
					break;
				}
			}
			if (!imate)
			{
				for (j=0;j<n_e;j++)
				{
					if (mat_e[j]==num_m[i])
					{
						imate=1;
						break;
					}
				}

			}
			if (!imate)
			{
				fprintf(stdout,
					"ERROR: material number %d has no domain particles or nodes!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			imate=0;
		}

		/* test on particles in domain material */
		if (((type_m[i]<4)||(type_m[i]>6))&&(type_m[i]!=0))
		{
			for (j=0;j<n_s;j++)
			{
				if (mat_s[j]==num_m[i])
				{
					imate=1;
					break;
				}
			}
			if (!imate)
			{
				fprintf(stdout,
					"ERROR: material number %d has no domain particles!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			imate=0;
		}

		/* test on elements in finite element material */
		if ((type_m[i]>3)&&(type_m[i]<7))
		{
			for (j=0;j<n_e;j++)
			{
				if (mat_e[j]==num_m[i])
				{
					imate=1;
					break;
				}
			}
			if (!imate)
			{
				fprintf(stdout,
					"ERROR: material number %d has no finite elements!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			imate=0;
		}

	  	/* test on data */
	  	if (data_check)
		{
			/* test on m */
			if (!type_m[i]) // rigid body
			{
				if (rho_m[i]==0.0) // rigid body positive density or negative mass
				{
					fprintf(stdout,
						"ERROR: density of rigid material %d must be a negative or a positive value!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						num_m[i]);
					return(0);
				}				
			}
			
			if ((type_m[i]>0)&&(rho_m[i]<=0.0)) // non-rigid material
			{
				fprintf(stdout,
					"ERROR: density of material %d must be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
		
			/* test on mu */
			if ((mu_m[i]<0.0)&&
				((type_m[i]<4)||(type_m[i]==12)||(type_m[i]==13)||(type_m[i]==23)||(type_m[i]==27)))
			{
				fprintf(stdout,
					"ERROR: viscosity of material %d must be 0 or positive!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((mu_m[i]<0.0)&&(mu_m[i]>0.5)&&(type_m[i]>3)&&(type_m[i]<7))
			{
				fprintf(stdout,
					"ERROR: Poisson's ratio of material %d must be between 0 and 0.5!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((mu_m[i]<0.0)&&(type_m[i]>3)&&(type_m[i]<7))
			{
				fprintf(stdout,
					"ERROR: Shear modulus of material %d be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
		
			/* test on T */
			if ((T_m[i]<=0.0)&&(type_m[i]==1))
			{
				fprintf(stdout,
					"ERROR: absolute temperature of material %d be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((T_m[i]<=0.0)&&(type_m[i]>3)&&(type_m[i]<7))
			{
				fprintf(stdout,
					"ERROR: Young's modulus of material %d be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((T_m[i]<=0.0)&&(type_m[i]>3)&&(type_m[i]<7))
			{
				fprintf(stdout,
					"ERROR: bulk modulus of material %d be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
		
			/* test on kappa */
			if ((type_m[i]==1)&&(kappa_m[i]<=1.0))
			{
				fprintf(stdout,
					"ERROR: specific heat ratio of material %d must be greater than 1!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
				num_m[i]);
				return(0);
			}
			if ((type_m[i]==4)&&(kappa_m[i]<=0.0))
			{
				fprintf(stdout,
					"ERROR: cross_sectional area of material %d be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((type_m[i]==5)&&(kappa_m[i]<=0.0))
			{
				fprintf(stdout,
					"ERROR: thickness of material %d be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}

			/* test on gamma */
			if ((type_m[i]==1)&&(gamma_m[i]<0.0))
			{
				fprintf(stdout,
					"ERROR: constant volume specific heat ratio of material %d must be zero or positive!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if (((type_m[i]==0)||
				(type_m[i]==2)||(type_m[i]==3)||
				(type_m[i]==7)|| // (type_m[i]==8)||(type_m[i]==9)||
				(type_m[i]==12)||(type_m[i]==13))
					&&(gamma_m[i]<1.0))
			{
				fprintf(stdout,
					"ERROR: power of material %d must be 1 or greater than 1!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((type_m[i]>3)&&(type_m[i]<7)&&(gamma_m[i]<0.0))
			{
				fprintf(stdout,
					"ERROR: damping of material %d must be greater than 0!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}

			/* test on K */
			if ((((type_m[i]==0)&&(rho_m[i]>0))||
				(type_m[i]==7))
				//((type_m[i]==7)||(type_m[i]==8)||(type_m[i]==9)))
				&&(T_m[i]<=0.0))
			{
				fprintf(stdout,
					"ERROR: bulk modulus of material %d must be greater than 0!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
				num_m[i]);
				return(0);
			}

			/* test on G */
			if (((type_m[i]==7)|| // (type_m[i]==8)|| // (type_m[i]==9)|| (type_m[i]==23)||
				(type_m[i]==27))&&
					(mu_m[i]<=0.0))
			{
				fprintf(stdout,
					"ERROR: shear modulus of material %d must be greater than 0!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
				num_m[i]);
				return(0);
			}

			/* test on G */
			if (((type_m[i]==23)||
				(type_m[i]==27))&&
					(mu_m[i]<=0.0))
			{
				fprintf(stdout,
					"ERROR: ciscosity of material %d must be greater than 0!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
				num_m[i]);
				return(0);
			}

			/* test on nu */
			double nu=(3.0 * T_m[i] - 2 * mu_m[i]) / 2 / (3 * T_m[i] + mu_m[i]);
			if ((type_m[i]==7) // ||(type_m[i]==8)) // ||(type_m[i]==9))
				&&((nu<0.0)||(nu>0.5)))
			{
				fprintf(stdout,
					"ERROR: Poisson's ratio s of material %d must be between 0 and 0.5!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
				num_m[i]);
				return(0);
			}
		}
    }

  	/* contacts */
  	for (i=0;i<n_c;i++)
    {
		/* test on zero or negative particle numbering */
		if (num_c[i]<=0)
		{
			fprintf(stdout,
				"ERROR: contact %d must be a positive integer!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}

		/* contact type */
		switch (sw_c[i])
		{
			case 0: 
			{
				strcpy(contact_type,"0 (null contact)");
				break;
			}
			case 1:
			{
				strcpy(contact_type,"1 (soft contact)");
				break;
			}
			case 2:
			{
				strcpy(contact_type,"2 (sliding without separation)");
				break;
			}
			case 3:
			{
				strcpy(contact_type,"3 (tied contact)");
				break;
			}
		}

		/* contact thickness */
		if (ct_c[i]<0)
		{
			sprintf(ct_type,"%f (times smoothing length)",-ct_c[i]);
		}
		else
		{
			sprintf(ct_type,"%f",ct_c[i]);
		}

      	fprintf(stdout,
	      	"contact number %d"
	      	"\n"
	      	" -> slave material = %d"
	      	"\n"
	      	" -> master material = %d"
	      	"\n"
	      	" -> contact type = %s"
	      	"\n"
	      	" -> ct = %s"
	      	"\n"
	      	" -> klin = %f"
	      	"\n"
	      	" -> knon = %f"
	      	"\n"
	      	" -> kf = %f"
	      	"\n"
	      	" -> kd = %f"
	      	"\n"
	      	"\n",
	      	num_c[i],mat_c[i],mat_c[i+n_c],contact_type,ct_type,klin_c[i],knon_c[i],kf_c[i],kd_c[i]);

      	/* test on both same materials */
      	if (mat_c[i]==mat_c[i+n_c])
		{
	  		fprintf(stdout,
		  		"ERROR: self-contact in contact %d is not supported!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
		  	return(0);
		}
	  
      	/* test on type */
      	if ((sw_c[i]!=0)&&(sw_c[i]!=1)&&(sw_c[i]!=2)&&(sw_c[i]!=3))
		{
	  		fprintf(stdout,
		  		"ERROR: unknown contact type in contact %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}

      	/* test on ct - for particles it can be negative*/
      	if ((mat_c[i]>3)&&(mat_c[i]<7)&&(ct_c[i]<0.0))
		{
	  		fprintf(stdout,
		  		"ERROR: contact thickness in contact %d must be a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}
	
      	/* test on klin */
      	if (klin_c[i]<0.0)
	  	{
		  	fprintf(stdout,
		  		"ERROR: linear stiffness coeffitient in contact %d must be a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}
	
      	/* test on knon */
      	if (knon_c[i]<0.0)
		{
	  		fprintf(stdout,
		  		"ERROR: nonlinear stiffnes coeffitient in contact %d must be a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
				num_c[i]);
	  	return(0);
		}
	
      	/* test on kf */
      	if (kf_c[i]<0.0)
		{
	  		fprintf(stdout,
		  		"ERROR: friciton coeffitient in contact %d must be a positive value!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}
      
      	/* test on kd */
      	if ((kd_c[i]<0.0)||(kd_c[i]>1.0))
		{
	  		fprintf(stdout,
		  		"ERROR: damping coeffitient in contact %d must be between 0 and 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}
      
      	/* test on double numbering */      
      	for (j=i+1;j<n_c;j++)
		{
	  		if (num_c[i]==num_c[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: contact %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		num_c[i]);
	      		return(0);
	    	}
		}
    }

  	/* renumbering of contact materials */
  	for (i=0;i<n_c;i++)
    {
      	for (j=0;j<n_m;j++)
		{
	  		if (mat_c[i]==num_m[j])
	    	{
	      		icontact=1;
	      		mat_c[i]=j;
	      		break;
	    	}
		}

      	if (!icontact)
		{
	  		fprintf(stdout,
		  		"ERROR: material number %d in contact %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		mat_c[i],num_c[i]);
	  		return(0);
		}

      	icontact=0;
      	for (j=0;j<n_m;j++)
		{
	  		if (mat_c[i+n_c]==num_m[j])
	    	{
	      		icontact=1;
	      		mat_c[i+n_c]=j;
	      		break;
	    	}
		}
      
      	if (!icontact)
		{
	  		fprintf(stdout,
		  		"ERROR: material number %d in contact %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		mat_c[i+n_c],num_c[i]);
	  		return(0);
		}
  
		/* test on particles as slave and master - after renumbering */
      	if ((sw_c[i]==0)&&(type_m[mat_c[i]]>3)&&(type_m[mat_c[i]]<7))
		{
	  		fprintf(stdout,
		  		"ERROR: contact %d must have particles as slave!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}
      	if ((sw_c[i]==0)&&(type_m[mat_c[i+n_c]]>3)&&(type_m[mat_c[i+n_c]]<7))
		{
	  		fprintf(stdout,
		  		"ERROR: contact %d must have partciles as master!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}
    }

  	/* renumbering of domain materials and initialization */
  	for (i=0;i<n_s;i++)
    {
      	for (j=0;j<n_m;j++)
		{
	  		if (mat_s[i]==num_m[j])
	    	{
	      		imate=1;
	      		mat_s[i]=j;
				rho_s[i]=rho_m[j]; // used as witch between rigid body types 0 and 1 below
				if ((type_m[j]==0)&&(rho_m[j]<0)) rho_m[j]=-rho_m[j];
				m_s[i]=rho_m[j]*V_s[i];
	      		u_s[i]=gamma_m[j]*T_m[j];
	      		mu_s[i]=mu_m[j];
	      		dudt_s[i]=0.0;
	      		drhodt_s[i]=0.0;
	      		h_s[i]=h0*pow(m_s[i]/rho_s[i],1.0/(double)dim);
	      		h_s_0[i]=h_s[i];

              	/* equation of state */
				switch (type_m[j])
              	{
              		case 0: /* rigid body */
              	  	{
						if (rho_s[i]>0) // rigid body type 0 -> c_s[i] and p_s[i] for smoothing (as for material 7)
						{
							c_s[i]=sqrt(T_m[j]/rho_m[j]);
							p_s[i]=kappa_m[j]+rho_m[j]*sqr(c_s[i])/gamma_m[j]*
								(pow(rho_s[i]/rho_m[j],gamma_m[j])-1.0);
						}
						else // rigid body type 1 -> c_s[i] and p_s[i] irrelevant
						{
							rho_s[i]=-rho_s[i]; // rigid body type 1 has negative density
							c_s[i]=10.0;
							p_s[i]=0.0;
						}
	              		break;
        	      	}
              		case 1: /* gas */
              	  	{
						//c_s[i]=10.0; for shock tube convergence, otherwise sqrt((kappa_m[j]-1.0)*u_s[i]);
              	    	c_s[i]=10.0;
              	    	p_s[i]=rho_s[i]*(kappa_m[j]-1.0)*u_s[i];
	              		break;
        	      	}
              		case 2: /* liquid */
              	  	{
	            		//c_s[i]=gamma_m[j]/rho_m[j];
	            		c_s[i]=sqrt(T_m[j]*gamma_m[j]/rho_m[j]);
	            		//p_s[i]=kappa_m[j]+gamma_m[j]*(rho_s[i]/rho_m[j]-1.0);
	            		p_s[i]=kappa_m[j]+T_m[j]*(pow((rho_s[i]/rho_m[j]),gamma_m[j])-1.0);
	              	    break;
              	  	}
              		case 3: /* SPH */
              	  	{
              	    	if (v_max==0) c_s[i]=10.0;
						else c_s[i]=10.0*v_max;
						p_s[i]=kappa_m[j]+100.0*rho_m[j]*sqr(v_max)/gamma_m[j]*
              	      		(pow((rho_s[i]/rho_m[j]),gamma_m[j])-1.0);
                    	break;
              	  	}
              		case 7: /* SPS = SPH with shear */
              		case 8: /* Neo-hookean */
              		case 9: /* user-defined */
              	  	{
              	    	c_s[i]=sqrt(T_m[j]/rho_m[j]);
						p_s[i]=kappa_m[j]+rho_m[j]*sqr(c_s[i])/gamma_m[j]*
							(pow(rho_s[i]/rho_m[j],gamma_m[j])-1.0);
                    	break;
              	  	}
              		case 12: /* SPH with tension */
              	  	{
						c_s[i]=sqrt(T_m[j]*gamma_m[j]/rho_m[j]);
						p_s[i]=T_m[j]*(pow(rho_s[i]/rho_m[j],gamma_m[j])-1.0);
              	    	break;
              	  	}
              		case 13: /* SPS with tension */
              	  	{
						c_s[i]=sqrt(T_m[j]*gamma_m[j]/rho_m[j]);
						p_s[i]=100.0*rho_m[j]*sqr(v_max)/gamma_m[j]*
							(pow(rho_s[i]/rho_m[j],gamma_m[j])-1.0);
                    	break;
              	  	}
              		case 23: /* Mie-Grueneisen EOS */
					/* {
              	    	if (v_max==0) c_s[i]=10.0;
						else c_s[i]=10.0*v_max;
					} */
              		case 27: /* SPS with Mie-Grueneisen EOS */
              	  	{
						// initial sound speed = sqrt(bulk modulus over initial density)
						//c_s[i]=sqrt(T_m[j]*gamma_m[j]/rho_m[j]);
						c_s[i]=sqrt(T_m[j]/rho_m[j]);

						s=coef1_m[mat_s[j]]; // linear coefficient of Hugoniot
						G0=coef2_m[mat_s[j]]; // Grueneisen constant

	            		C0=kappa_m[mat_s[j]]; // initial pressure
	            		C1=rho_m[mat_s[j]]*sqr(c_s[i])-G0/2.0*kappa_m[mat_s[j]];
	            		C2=rho_m[mat_s[j]]*sqr(c_s[i])*(2.0*s-1.0)-
							G0/2.0*rho_m[mat_s[j]]*sqr(c_s[i]);
	            		C3=rho_m[mat_s[j]]*sqr(c_s[i])*(2.0*s-1.0)*(3.0*s-1.0)-
							G0/2.0*rho_m[mat_s[j]]*sqr(c_s[i])*(2.0*s-1.0);
	            		C4=G0;C5=G0;En=0.0; // energy per unit of initial volume
	            		Eta=rho_s[i]/rho_m[mat_s[j]]-1.0;

						/* initial pressure (Taddei, 2015; Frissane, 2019) */
				    	if (Eta>=0) p_s[i]=C0+C1*Eta+C2*sqr(Eta)+C3*pow(Eta,3.0)+(C4+C5*Eta)*En;
						else p_s[i]=C0+C1*Eta;
                    	break;
              	  	}
              	}
				/* initial sound speed */
				c_s_0[i]=c_s[i];
	      		break;
		    }
		}
		if (!imate)
		{
		  	fprintf(stdout,
			  	"ERROR: material number %d in domain particle %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
			  	"\n",
			  	mat_s[i],num_s[i]);
		  	return(0);
		}
    	imate=0;

		if ((type_m[mat_s[i]]>3)&&(type_m[mat_s[i]]<7)) // FEM material without rigid body
		{
		  	fprintf(stdout,
			  	"ERROR: material number %d in domain particle %d is a FEM material!"
				"\n"
				"\n"
				"error termination"
				"\n"
			  	"\n",
			  	num_m[mat_s[i]],num_s[i]);
		  	return(0);
		}
    }

  	/* renumbering of element materials */
  	for (i=0;i<n_e;i++)
    {
      	for (j=0;j<n_m;j++)
		{
	  		if (mat_e[i]==num_m[j])
	    	{
	      		imate=1;
	      		mat_e[i]=j;
				break;
	    	}
		}
      	if (!imate)
		{
	  		fprintf(stdout,
		  		"ERROR: material number %d in element %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_m[mat_e[i]],num_b[i]);
	  		return(0);
		}
      	imate=0;

		if ((type_m[mat_e[i]]!=0)&&((type_m[mat_e[i]]<4)||(type_m[mat_e[i]]>6))) // SPH material without rigid body
		{
		  	fprintf(stdout,
			  	"ERROR: material number %d in element %d is a SPH material!"
				"\n"
				"\n"
				"error termination"
				"\n"
			  	"\n",
			  	num_m[mat_e[i]],num_e[i]);
		  	return(0);
		}
    }
  
  	/* renumbering of element nodes */
  	for (i=0;i<n_e;i++)
    {
      	for (j=0;j<n_b;j++)
		{
	  		if (nod_e[i]==num_b[j])
	    	{
	    	  	inode=1;
	      		nod_e[i]=j;
	      		break;
	    	}
		}
      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: node %d of element %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		nod_e[i],num_e[i]);
	  		return(0);
		}
      	inode=0;

      	for (j=0;j<n_b;j++)
		{
	  		if (nod_e[i+n_e]==num_b[j])
	    	{
	      		inode=1;
	      		nod_e[i+n_e]=j;
	      		break;
	    	}
		}
      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: node %d of element %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		nod_e[i+2*n_e],num_e[i]);
	  		return(0);
		}
      	inode=0;

      	for (j=0;j<n_b;j++)
		{
	  		if (nod_e[i+2*n_e]==num_b[j])
	    	{
	      		inode=1;
	      		nod_e[i+2*n_e]=j;
	      		break;
	    	}
		}
      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: node %d of element %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		nod_e[i+2*n_e],num_e[i]);
	  		return(0);
		}
      	inode=0;

      	for (j=0;j<n_b;j++)
		{
	  		if (nod_e[i+3*n_e]==num_b[j])
	    	{
	      		inode=1;
	      		nod_e[i+3*n_e]=j;
	      		break;
	    	}
		}
      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: node %d of element %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		nod_e[i+3*n_e],num_e[i]);
	  		return(0);
		}
    }

  	/* node numbers according to materials */
  	for (i=0;i<n_e;i++)
    {
      	switch(dim)
		{
			case 1:
	  		{
	    		nod_e[i+2*n_e]=nod_e[i+n_e];
	    		nod_e[i+3*n_e]=nod_e[i+n_e];
	    		break;
	  		}
			case 2:
	  		{
	    		if ((type_m[mat_e[i]]==4)||(type_m[mat_e[i]]==5))
	      		{
					nod_e[i+2*n_e]=nod_e[i+n_e];
					nod_e[i+3*n_e]=nod_e[i+n_e];
	      		}
	    		break;
	  		}
			case 3:
	  		{
	    		if ((type_m[mat_e[i]]==4)||(type_m[mat_e[i]]==5))
	      		{
					nod_e[i+3*n_e]=nod_e[i+2*n_e];
	      		}
	    		break;
	  		}
		}
    }

	/* element areas and test on element areas */
	if (n_e>0)
	{
	  	volume(type_m, // materials
		 	n_b,x_b, // initial values
		 	n_e,num_e,mat_e,nod_e, // elements
		 	V_e,
	 		dim);
	}

  	/* creating of contact slaves and masters */
  	for (i=0;i<n_c;i++)
    {
		/* slave segment */
      	if (type_m[mat_c[i]]<4||type_m[mat_c[i]]>6)
		{
	  		for (j=0;j<n_s;j++)
	    	{
	      		if (mat_s[j]==mat_c[i])
				{
		  			seg_c[i]++;
				}
	    	}
		}
    	else
		{
	  		for (j=0;j<n_b;j++)
	    	{
	      		for (k=0;k<n_e;k++)
				{
		  			if (((nod_e[k]==j)||(nod_e[k+n_e]==j)||(nod_e[k+2*n_e]==j)||
						(nod_e[k+3*n_e]==j))&&(mat_e[k]==mat_c[i]))
		    		{
		      			seg_c[i]++;
		      			break;
		    		}
				}
	    	}
		}

      	if (seg_c[i]==0)
		{
	  		fprintf(stdout,
		  		"ERROR: contact %d has no slave nodes!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}

      	count_slave=count_slave+seg_c[i];
	
      	/* master segment */
      	if (type_m[mat_c[i+n_c]]<4||type_m[mat_c[i+n_c]]>6)
		{
	  		for (j=0;j<n_s;j++)
	    	{
	      		if (mat_s[j]==mat_c[i+n_c])
				{
		  			seg_c[i+n_c]++;
				}
	    	}
		}
		else
		{
			for (j=0;j<n_e;j++)
			{
				if (mat_e[j]==mat_c[i+n_c])
					{
						seg_c[i+n_c]++;
					}
			}
		}
	  
      	if (seg_c[i+n_c]==0)
		{
	  		fprintf(stdout,
		  		"ERROR: material %d has no master elements!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		num_c[i]);
	  		return(0);
		}

      	count_master=count_master+seg_c[i+n_c];
    }

  	if (count_slave>0)
    {
      	m_slave=createMemMore(double,count_slave);
      	memset(m_slave,0,count_slave*sizeof(double));
      	slave_c=createMemMore(int,count_slave);
    }

  	if (count_master>0)
	{
		m_master=createMemMore(double,count_master);
		memset(m_master,0,count_master*sizeof(double));
		master_c=createMemMore(int,count_master);
    }

  	count_slave=0;
  	count_master=0;

  	/* set segment counter to zero */
  	for (i=0;i<n_c;i++)
    {
		/* slave segment */
      	if (type_m[mat_c[i]]<4||type_m[mat_c[i]]>6)
		{
	  		for (j=0;j<n_s;j++)
	    	{
	      		if (mat_s[j]==mat_c[i])
				{
		  			slave_c[count_slave]=j;
		  			m_slave[count_slave]=m_s[j];
		  			count_slave++;
				}
	    	}
		}
      	else
		{
	  		for (j=0;j<n_b;j++)
	    	{
	      		for (k=0;k<n_e;k++)
				{
		  			if (((nod_e[k]==j)||(nod_e[k+n_e]==j)||(nod_e[k+2*n_e]==j)||
						(nod_e[k+3*n_e]==j))&&(mat_e[k]==mat_c[i]))
		    		{
		      			slave_c[count_slave]=j;

		      			switch (dim)
						{
							case 1:
			  				{
			    				m_slave[count_slave]=m_slave[count_slave]+rho_m[mat_e[k]]*V_e[k];
			    				break;
			  				}
							case 2:
			  				{
			    				if (nod_e[k+n_e]==nod_e[k+2*n_e]) m_slave[count_slave]=m_slave[count_slave]+
									rho_m[mat_e[k]]*kappa_m[mat_e[k]]*V_e[k]/2.0;
			    				else
			      				{
									if (nod_e[k+2*n_e]==nod_e[k+3*n_e]) m_slave[count_slave]=
										m_slave[count_slave]+rho_m[mat_e[k]]*kappa_m[mat_e[k]]*V_e[k]/3.0;
									else m_slave[count_slave]=m_slave[count_slave]+
										rho_m[mat_e[k]]*kappa_m[mat_e[k]]*V_e[k]/4.0;
			      				}
			    				break;
			  				}
							case 3:
			  				{
			    				if (nod_e[k+2*n_e]==nod_e[k+3*n_e]) m_slave[count_slave]=
									m_slave[count_slave]+rho_m[mat_e[k]]*kappa_m[mat_e[k]]*V_e[k]/3.0;
			    				else m_slave[count_slave]=m_slave[count_slave]+rho_m[mat_e[k]]*V_e[k]/4.0;

			    				break;
			  				}
						}

		      			count_slave++;
		      			break;
		    		}
				}
	    	}
		}
      
      	/* master segment */
      	if (type_m[mat_c[i+n_c]]<4||type_m[mat_c[i+n_c]]>6)
		{
	  		for (j=0;j<n_s;j++)
	    	{
	      		if (mat_s[j]==mat_c[i+n_c])
				{
		  			master_c[count_master]=j;
		  			m_master[count_master]=m_s[j];
		  			count_master++;
				}
	    	}
		}
      	else
		{
			for (j=0;j<n_e;j++)
			{
				if (mat_e[j]==mat_c[i+n_c])
				{
					master_c[count_master]=j;
					m_master[count_master]=rho_m[mat_e[j]]*V_e[j];
					count_master++;
				}
			}
		}
    }

  	/* test on initial pressure */
  	for (i=0;i<n_p;i++)
  	{
		inode=0;
		for (j=0;j<n_s;j++)
		{
	  	if (inpre[i]==num_s[j])
	    	{
	      		inode=1;
	      		switch (type_m[mat_s[i]])
				{
					case 1:
					{
						u_s[j]=p_0[i]/(rho_m[mat_s[i]]*(kappa_m[mat_s[i]]-1.0));

						if (u_s[j]<=0)
						{
							fprintf(stdout,
								"ERROR: internal energy of particle %d be a positive value!"
								"\n"
								"\n"
								"error termination"
								"\n"
								"\n",
								num_s[i]);
							return(0);
						} 
						break;
					}
					case 2:
					case 3:
					{
						//rho_s[j]=rho_m[mat_s[i]]*pow((p_0[i]-kappa_m[mat_s[i]])/
						//	(100.0*sqr(v_max)*rho_m[mat_s[i]])*gamma_m[mat_s[i]]+1.0,1.0/gamma_m[mat_s[i]]);
						rho_s[j]=rho_m[mat_s[i]]*pow((p_0[i]-kappa_m[mat_s[i]])/
							(100.0*rho_m[mat_s[i]])*gamma_m[mat_s[i]]+1.0,1.0/gamma_m[mat_s[i]]);
						
						if (rho_s[j]<=0)
						{
							fprintf(stdout,
								"ERROR: density of particle %d be a positive value!"
								"\n"
								"\n"
								"error termination"
								"\n"
								"\n",
								num_s[i]);
							return(0);
						}
						break;
					}
					case 7:
					case 8:
					case 9:
					{
						rho_s[j]=rho_m[mat_s[i]]*pow((p_0[i]-kappa_m[mat_s[i]])/
							T_m[mat_s[i]]*gamma_m[mat_s[i]]+1.0,1.0/gamma_m[mat_s[i]]);

						if (rho_s[j]<=0)
						{
							fprintf(stdout,
								"ERROR: density of particle %d be a positive value!"
								"\n"
								"\n"
								"error termination"
								"\n"
								"\n",
								num_s[i]);
							return(0);
						}
						break;
					}
					case 12:
					case 13:
					{
						//rho_s[j]=rho_m[mat_s[i]]*pow((p_0[i]-kappa_m[mat_s[i]])/
						//	(100.0*sqr(v_max)*rho_m[mat_s[i]])*gamma_m[mat_s[i]]+1.0,1.0/gamma_m[mat_s[i]]);
						rho_s[j]=rho_m[mat_s[i]]*pow((p_0[i]-kappa_m[mat_s[i]])/
							(100.0*rho_m[mat_s[i]])*gamma_m[mat_s[i]]+1.0,1.0/gamma_m[mat_s[i]]);

						if (rho_s[j]<=0)
						{
							fprintf(stdout,
								"ERROR: density of particle %d be a positive value!"
								"\n"
								"\n"
								"error termination"
								"\n"
								"\n",
								num_s[i]);
							return(0);
						}
						break;
					}
					case 23:
					case 27:
					{
						//rho_s[j]=rho_m[mat_s[i]]*pow((p_0[i]-kappa_m[mat_s[i]])/
						//	(100.0*sqr(v_max)*rho_m[mat_s[i]])*gamma_m[mat_s[i]]+1.0,1.0/gamma_m[mat_s[i]]);
						rho_s[j]=rho_m[mat_s[i]]*pow((p_0[i]-kappa_m[mat_s[i]])/
							T_m[mat_s[i]]*gamma_m[mat_s[i]]+1.0,1.0/gamma_m[mat_s[i]]);

						if (rho_s[j]<=0)
						{
							fprintf(stdout,
								"ERROR: density of particle %d be a positive value!"
								"\n"
								"\n"
								"error termination"
								"\n"
								"\n",
								num_s[i]);
							return(0);
						}
						break;
					}
					break;
				}
			}
		}

    	if (!inode)
		{
			fprintf(stdout,
				"ERROR: particle %d in initial pressure condition is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				inpre[i]);
			return(0);
		}
      
    	/* test on double initial pressures */
    	for (j=i+1;j<n_p;j++)
		{
	  		if (inpre[i]==inpre[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: initial pressure number %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		inpre[i]);
	      		return(0);
	    	}
		}     
    }

  	/* test on initial velocities */
  	for (i=0;i<n_v;i++)
    {
		inode=0;
		for (j=0;j<n_s;j++)
		{
	  		if (invel[i]==num_s[j])
	    	{
	      		inode=1;
				//frame_v_0_s[j]=invel_frame[i];
	      		v_s[j]=v_0[i];
	      		v_s[j+n_s]=v_0[i+n_v];
	      		v_s[j+2*n_s]=v_0[i+2*n_v];
	      		break;
		    }
		}

      	for (j=0;j<n_b;j++)
		{
			if (invel[i]==num_b[j])
	    	{
	      		inode=1;
				//frame_v_0_b[j]=invel_frame[i];
	      		v_b[j]=v_0[i];
	      		v_b[j+n_b]=v_0[i+n_v];
	      		v_b[j+2*n_b]=v_0[i+2*n_v];
	      		break;
	    	}
		}

      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: node %d in initial velocity condition is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		invel[i]);
	  		return(0);
		}

      	/* test on double invel */
      	for (j=i+1;j<n_v;j++)
		{
	  		if (invel[i]==invel[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: initial velocity number %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		invel[i]);
	      		return(0);
	    	}
		}     
    
		/* test on invel frame */
		if ((invel_frame[i]!=0)&&(invel_frame[i]!=1))
		{
			fprintf(stdout,
				"ERROR: initial velocity %d frame must be 0 or 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				invel[i]);
			return(0);
		}
	}

  	/* test on acceleration fields */
  	for (i=0;i<n_a;i++)
    {
      	inode=0;
      	for (j=0;j<n_s;j++)
		{
	  		if (acfld[i]==num_s[j])
	    	{
	      		n_a_0_s++;
	      		inode=1;
	    	}
		}

      	for (j=0;j<n_b;j++)
		{
	  		if (acfld[i]==num_b[j])
	    	{
	      		n_a_0_b++;
	      		inode=1;
	    	}
		}

      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: accel node %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		acfld[i]);
	  		return(0);
		}
      
      	/* test on double acceleration fields */
      	for (j=i+1;j<n_a;j++)
		{
	  		if (acfld[i]==acfld[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: accel %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		acfld[i]);
	      		return(0);
	    	}
		}     

		/* test on functions */
		if ((acfld_type[i]!=0)&&(acfld_type[i]!=1))
		{			
			fprintf(stdout,
				"ERROR: frame %d type must be 0 or 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				acfld[i]);
			return(0);
		}
		else
		{
			if (acfld_type[i]==1) // acceleration field by fuctions
			{
				if ((int)a_0[i]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)a_0[i]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in x acceleration field number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)a_0[i],acfld[i]);
						return(0);
					}
				}

				if ((int)a_0[i+n_a]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)a_0[i+n_a]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in y acceleration field number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)a_0[i+n_a],acfld[i]);
						return(0);
					}
				}

				if ((int)a_0[i+2*n_a]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)a_0[i+2*n_a]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in y acceleration field number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)a_0[i+2*n_a],acfld[i]);
						return(0);
					}
				}
			} // else acceleration field by real numbers
		}
    }

  	/* allocation - concentrated accelerations */
  	if (n_a_0_s>0)
    {
      	ind_a_0_s=createMemMore(int,n_a_0_s);
      	type_a_0_s=createMemMore(int,n_a_0_s);
      	frame_a_0_s=createMemMore(int,n_a_0_s);
      	a_0_s=createMemMore(double,3*n_a_0_s);
    }
  
  	if (n_a_0_b>0)
    {
      	ind_a_0_b=createMemMore(int,n_a_0_b);
      	type_a_0_b=createMemMore(int,n_a_0_b);
      	frame_a_0_b=createMemMore(int,n_a_0_b);
      	a_0_b=createMemMore(double,3*n_a_0_b);
    }
  
  	/* concentrated accelerations */
  	for (i=0;i<n_a;i++)
    {
      	for (j=0;j<n_s;j++)
		{
	  		if (acfld[i]==num_s[j])
	    	{
	      		ind_a_0_s[i_a_0_s]=j;
	      		type_a_0_s[i_a_0_s]=acfld_type[i];
				frame_a_0_s[i_a_0_s]=acfld_frame[i];
	      		a_0_s[i_a_0_s]=a_0[i];
	      		a_0_s[i_a_0_s+n_a_0_s]=a_0[i+n_a];
	      		a_0_s[i_a_0_s+2*n_a_0_s]=a_0[i+2*n_a];
	      		i_a_0_s++;
	    	}
		}

      	for (j=0;j<n_b;j++)
		{
	  		if (acfld[i]==num_b[j])
	    	{
	      		ind_a_0_b[i_a_0_b]=j;
	      		type_a_0_b[i_a_0_b]=acfld_type[i];
	      		frame_a_0_b[i_a_0_b]=acfld_frame[i];
	      		a_0_b[i_a_0_b]=a_0[i];
	      		a_0_b[i_a_0_b+n_a_0_b]=a_0[i+n_a];
	      		a_0_b[i_a_0_b+2*n_a_0_b]=a_0[i+2*n_a];
	      		i_a_0_b++;    
	    	}
	  
		}
    }
  
  	/* test on force */
  	for (i=0;i<n_f;i++)
    {
      	inode=0;
      	for (j=0;j<n_s;j++)
		{
	  		if (force[i]==num_s[j])
	    	{
	      		n_f_0_s++;
	      		inode=1;
			}
		}

      	for (j=0;j<n_b;j++)
		{
	  		if (force[i]==num_b[j])
	    	{
	      		n_f_0_b++;
	      		inode=1;
	    	}
		}

      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: force node %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		force[i]);
	  		return(0);
		}
      
      	/* test on double force */
      	for (j=i+1;j<n_f;j++)
		{
			if (force[i]==force[j])
			{
				fprintf(stdout,
		      		"ERROR: force %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		force[i]);
	      		return(0);
	    	}
		}     

		/* test on functions */
		if ((force_type[i]!=0)&&(force_type[i]!=1))
		{			
			fprintf(stdout,
				"ERROR: frame %d type must be 0 or 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				force[i]);
			return(0);
		}
		else
		{
			if (force_type[i]==1) // force and moments by fuctions
			{
				if ((int)f_0[i]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)f_0[i]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in x force number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)f_0[i],force[i]);
						return(0);
					}
				}

				if ((int)f_0[i+n_f]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)f_0[i+n_f]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in y force number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)f_0[i+n_f],force[i]);
						return(0);
					}
				}

				if ((int)f_0[i+2*n_f]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)f_0[i+2*n_f]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in y force number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)f_0[i+2*n_f],force[i]);
						return(0);
					}
				}

				if ((int)M_0[i]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)M_0[i]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in x moment number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)M_0[i],force[i]);
						return(0);
					}
				}

				if ((int)M_0[i+n_f]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)M_0[i+n_f]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in y moment number %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)M_0[i+n_f],force[i]);
						return(0);
					}
				}

				if ((int)M_0[i+2*n_f]>0) // 0 means free
				{
					ibfun=0;
					for (k=0;k<n_u;k++)
					{
						if ((int)M_0[i+2*n_f]==num_u[j])
						{
							ibfun=1;
							break;
						}
					}
					if (!ibfun)
					{
						fprintf(stdout,
							"ERROR: function %d in y moment nummber %d does not exist!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							(int)M_0[i+2*n_f],force[i]);
						return(0);
					}
				}
			} // else force and moments by real numbers
		}
		
		/* test on moment frame */
		if ((force_frame[i]!=0)&&(force_frame[i]!=1))
		{
			fprintf(stdout,
				"ERROR: initial velocity %d frame must be 0 or 1!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				force[i]);
			return(0);
		}
    }

  	/* allocation - concentrated forces */
  	if (n_f_0_s>0)
    {
      	ind_f_0_s=createMemMore(int,n_f_0_s);
      	type_f_0_s=createMemMore(int,n_f_0_s);
      	f_0_s=createMemMore(double,3*n_f_0_s);
      	M_0_s=createMemMore(double,3*n_f_0_s);
    }

  	if (n_f_0_b>0)
    {
      	ind_f_0_b=createMemMore(int,n_f_0_b);
      	type_f_0_b=createMemMore(int,n_f_0_b);
      	f_0_b=createMemMore(double,3*n_f_0_b);
      	M_0_b=createMemMore(double,3*n_f_0_b);
    }

  	/* concentrated force */
  	for (i=0;i<n_f;i++)
    {
      	for (j=0;j<n_s;j++)
		{
	  		if (force[i]==num_s[j])
	    	{
	      		ind_f_0_s[i_f_0_s]=j;
	      		type_f_0_s[i_f_0_s]=force_type[i];
				frame_f_0_s[i_f_0_s]=force_frame[i];
	      		f_0_s[i_f_0_s]=f_0[i];
	      		f_0_s[i_f_0_s+n_f_0_s]=f_0[i+n_f];
	      		f_0_s[i_f_0_s+2*n_f_0_s]=f_0[i+2*n_f];
	      		M_0_s[i_f_0_s]=M_0[i];
	      		M_0_s[i_f_0_s+n_f_0_s]=M_0[i+n_f];
	      		M_0_s[i_f_0_s+2*n_f_0_s]=M_0[i+2*n_f];
	      		i_f_0_s++;
	    	}
		}

      	for (j=0;j<n_b;j++)
		{
	  		if (force[i]==num_b[j])
	    	{
	      		ind_f_0_b[i_f_0_b]=j;
	      		type_f_0_b[i_f_0_b]=force_type[i];
				frame_f_0_b[i_f_0_b]=force_frame[i];
	      		f_0_b[i_f_0_b]=f_0[i];
	      		f_0_b[i_f_0_b+n_f_0_b]=f_0[i+n_f];
	      		f_0_b[i_f_0_b+2*n_f_0_b]=f_0[i+2*n_f];
	      		M_0_b[i_f_0_b]=M_0[i];
	      		M_0_b[i_f_0_b+n_f_0_b]=M_0[i+n_f];
	      		M_0_b[i_f_0_b+2*n_f_0_b]=M_0[i+2*n_f];
	      		i_f_0_b++;    
	    	}
	  
		}
    }

  	/* test on damping */
  	for (i=0;i<n_d;i++)
    {
      	inode=0;
      	for (j=0;j<n_b;j++)
		{
	  		if (ndamp[i]==num_b[j])
	    	{
	      		n_d_0_b++;
	      		inode=1;
	    	}
		}

      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: node %d in damping condition is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		ndamp[i]);
	  		return(0);
		}
      
      	/* test on double damping */
      	for (j=i+1;j<n_d;j++)
		{
	  		if (ndamp[i]==ndamp[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: damping %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		ndamp[i]);
	      		return(0);
	    	}
		}     
    }
  
  	/* allocation - dampings */
  	if (n_d_0_b>0)
    {
      	ind_d_0_b=createMemMore(int,n_d_0_b);
      	d_0_b=createMemMore(double,3*n_d_0_b);
    }

  	/* dampings */
  	for (i=0;i<n_d;i++)
    {
      	for (j=0;j<n_b;j++)
		{
	  		if (ndamp[i]==num_b[j])
	    	{
	      		ind_d_0_b[i_d_0_b]=j;
	      		d_0_b[i_d_0_b]=d_0[i];
	      		d_0_b[i_d_0_b+n_d_0_b]=d_0[i+n_d];
	      		d_0_b[i_d_0_b+2*n_d_0_b]=d_0[i+2*n_d];
	      		i_d_0_b++;
	    	}
		}
    }

  	/* test on bounc */
  	for (i=0;i<n_o;i++)
    {
      	inode=0;
      	for (j=0;j<n_s;j++)
		{
	  		if (bounc[i]==num_s[j])
	    	{
				n_o_0_s++;
	      		inode=1;
				break;
	    	}
		}

      	if (!inode)
		{
			for (j=0;j<n_b;j++)
			{
				if (bounc[i]==num_b[j])
				{
					n_o_0_b++;
					inode=1;
					break;
				}
			}
		}

      	if (!inode)
		{
	  		fprintf(stdout,
		  		"ERROR: node %d in boundary condition is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
		  		"\n",
		  		bounc[i]);
	  		return(0);
		}

      	/* test on double bounc */
      	for (j=i+1;j<n_o;j++)
		{
	  		if (bounc[i]==bounc[j])
	    	{
	      		fprintf(stdout,
		      		"ERROR: boundary condition %d occurs more than once!"
					"\n"
					"\n"
					"error termination"
					"\n"
		      		"\n",
		      		bounc[i]);
	      		return(0);
	    	}
		
		}

		/* test on bounc type */
		if ((bounc_type[i]!=-1)&&(bounc_type[i]!=0)&&(bounc_type[i]!=1)&&(bounc_type[i]!=2)&&(bounc_type[i]!=3))
		{
			fprintf(stdout,
				"ERROR: boundary condition %d type must be 0 or 1 or 2 or 3!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				bounc[i]);
			return(0);
		}

      	/* test on bounc */
		if (bounc_type[i]<=0) // static of keeping initial conditions
		{
			/* x boundary condition */
			if ((o_0[i]!=0)&&(o_0[i]!=1)) // 0 means free, -1 means fixed
			{
				fprintf(stdout,
					"ERROR: x-constraint must be 0 or 1 in static boundary condition %d!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					bounc[i]);
				return(0);
			}

			/* y boundary condition */
			if ((o_0[i+n_o]!=0)&&(o_0[i+n_o]!=1)) // 0 means free, -1 means fixed
			{
					fprintf(stdout,
						"ERROR: y-constraint must be 0 or 1 in static boundary condition %d!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						bounc[i]);
					return(0);
			}

			/* z boundary condition */
			if ((o_0[i+2*n_o]!=0)&&(o_0[i+2*n_o]!=1)) // 0 means free, -1 means fixed
			{
					fprintf(stdout,
						"ERROR: z-constraint must be 0 or 1 in static boundary condition %d!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						bounc[i]);
					return(0);
			}

			/* angular x boundary condition */
			if ((o_0[i+3*n_o]!=0)&&(o_0[i+3*n_o]!=1)) // 0 means free, -1 means fixed
			{
					fprintf(stdout,
						"ERROR: psi-constraint must be 0 or 1 in static boundary condition %d!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						bounc[i]);
					return(0);
			}

			/* angular y boundary condition */
			if ((o_0[i+4*n_o]!=0)&&(o_0[i+4*n_o]!=1)) // 0 means free, -1 means fixed
			{
					fprintf(stdout,
						"ERROR: theta-constraint must be 0 or 1 in static boundary condition %d!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						bounc[i]);
					return(0);
			}

			/* angular z boundary condition */
			if ((o_0[i+5*n_o]!=0)&&(o_0[i+5*n_o]!=1)) // 0 means free, -1 means fixed
			{
					fprintf(stdout,
						"ERROR: phi-constraint must be 0 or 1 in static boundary condition %d!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						bounc[i]);
					return(0);
			}
		}
		else  // test on function in bounc
		{
			/* x boundary condition */
			if (o_0[i]>0) // 0 means free, -1 means fixed
			{
				ibfun=0;
				for (j=0;j<n_u;j++)
				{
					if (o_0[i]==num_u[j])
					{
						ibfun=1;
						break;
					}
				}

				if (!ibfun)
				{
					fprintf(stdout,
						"ERROR: function %d in x boundary condition %d does not exist!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						o_0[i],bounc[i]);
					return(0);
				}
			}

			/* y boundary condition */
			if (o_0[i+n_o]>0) // 0 means free, -1 means fixed
			{
				ibfun=0;
				for (j=0;j<n_u;j++)
				{
					if (o_0[i+n_o]==num_u[j])
					{
						ibfun=1;
						break;
					}
				}

				if (!ibfun)
				{
					fprintf(stdout,
						"ERROR: function %d in y boundary condition %d does not exist!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						o_0[i+n_o],bounc[i]);
					return(0);
				}
			}

			/* z boundary condition */
			if (o_0[i+2*n_o]>0) // 0 means free, -1 means fixed
			{
				ibfun=0;
				for (j=0;j<n_u;j++)
				{
					if (o_0[i+2*n_o]==num_u[j])
					{
						ibfun=1;
						break;
					}
				}

				if (!ibfun)
				{
					fprintf(stdout,
						"ERROR: function %d in z boundary condition %d does not exist!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						o_0[i+2*n_o],bounc[i]);
					return(0);
				}
			}

			/* angular x boundary condition */
			if (o_0[i+3*n_o]>0) // 0 means free, -1 means fixed
			{
				ibfun=0;
				for (j=0;j<n_u;j++)
				{
					if (o_0[i+3*n_o]==num_u[j])
					{
						ibfun=1;
						break;
					}
				}

				if (!ibfun)
				{
					fprintf(stdout,
						"ERROR: function %d in angular x boundary condition %d does not exist!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						o_0[i+3*n_o],bounc[i]);
					return(0);
				}
			}

			/* angular y boundary condition */
			if (o_0[i+4*n_o]>0) // 0 means free, -1 means fixed
			{
				ibfun=0;
				for (j=0;j<n_u;j++)
				{
					if (o_0[i+4*n_o]==num_u[j])
					{
						ibfun=1;
						break;
					}
				}

				if (!ibfun)
				{
					fprintf(stdout,
						"ERROR: function %d in angular y boundary condition %d does not exist!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						o_0[i+4*n_o],bounc[i]);
					return(0);
				}
			}

			/* angular z boundary condition */
			if (o_0[i+5*n_o]>0) // 0 means free, -1 means fixed
			{
				ibfun=0;
				for (j=0;j<n_u;j++)
				{
					if (o_0[i+5*n_o]==num_u[j])
					{
						ibfun=1;
						break;
					}
				}

				if (!ibfun)
				{
					fprintf(stdout,
						"ERROR: function %d in angular z boundary condition %d does not exist!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						o_0[i+5*n_o],bounc[i]);
					return(0);
				}
			}
		}

		/* test on boundary conditions frames */
		if ((bounc_frame[i]!=0)&&(bounc_frame[i]!=1)&&
			(bounc_frame[i]!=2)&&(bounc_frame[i]!=3))
		{
			fprintf(stdout,
				"ERROR: boundary condition %d frame must be 0 or 1 or 2 or 3!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				bounc[i]);
			return(0);
		}
    }

  	/* allocation - boundary conditions */
  	if (n_o_0_s>0)
    {
      	ind_o_0_s=createMemMore(int,n_o_0_s);
      	type_o_0_s=createMemMore(int,n_o_0_s);
      	o_0_s=createMemMore(int,6*n_o_0_s);
    }
  
  	if (n_o_0_b>0)
    {
      	ind_o_0_b=createMemMore(int,n_o_0_b);
      	type_o_0_b=createMemMore(int,n_o_0_b);
      	o_0_b=createMemMore(int,6*n_o_0_b);
    }
  
  	/* distributing boundary conditions between particles and nodes */
	for (i=0;i<n_o;i++)
    {
      	for (j=0;j<n_s;j++)
		{
	  		if (bounc[i]==num_s[j])
	    	{
	      		ind_o_0_s[i_o_0_s]=j;
	      		type_o_0_s[i_o_0_s]=bounc_type[i];

	      		o_0_s[i_o_0_s]=o_0[i];
	      		o_0_s[i_o_0_s+n_o_0_s]=o_0[i+n_o];
	      		o_0_s[i_o_0_s+2*n_o_0_s]=o_0[i+2*n_o];
	      		o_0_s[i_o_0_s+3*n_o_0_s]=o_0[i+3*n_o];
	      		o_0_s[i_o_0_s+4*n_o_0_s]=o_0[i+4*n_o];
	      		o_0_s[i_o_0_s+5*n_o_0_s]=o_0[i+5*n_o];

	      		/* 
				-> 0 = unconstrained
				-> other = constrained 
				*/
				constrained_s[j]=o_0_s[i_o_0_s]; // constraint 0/1 in motion x
	      		constrained_s[j+n_s]=o_0_s[i_o_0_s+n_o_0_s]; // constraint 0/1 in motion y
	      		constrained_s[j+2*n_s]=o_0_s[i_o_0_s+2*n_o_0_s]; // constraint 0/1 in motion z
	      		constrained_s[j+3*n_s]=o_0_s[i_o_0_s+3*n_o_0_s]; // constraint 0/1 in rotation x
	      		constrained_s[j+4*n_s]=o_0_s[i_o_0_s+4*n_o_0_s]; // constraint 0/1 in rotation y
	      		constrained_s[j+5*n_s]=o_0_s[i_o_0_s+5*n_o_0_s]; // constraint 0/1 in rotation z

				constrained_s_frame[j]=bounc_frame[i]; // coordinate frame

	      		i_o_0_s++;
				break;
	    	}
		}

      	for (j=0;j<n_b;j++)
		{
	  		if (bounc[i]==num_b[j])
	    	{
	      		ind_o_0_b[i_o_0_b]=j;
	      		type_o_0_b[i_o_0_b]=bounc_type[i];

	      		o_0_b[i_o_0_b]=o_0[i];
	      		o_0_b[i_o_0_b+n_o_0_b]=o_0[i+n_o];
	      		o_0_b[i_o_0_b+2*n_o_0_b]=o_0[i+2*n_o];
	      		o_0_b[i_o_0_b+3*n_o_0_b]=o_0[i+3*n_o];
	      		o_0_b[i_o_0_b+4*n_o_0_b]=o_0[i+4*n_o];
	      		o_0_b[i_o_0_b+5*n_o_0_b]=o_0[i+5*n_o];

	      		/* 
				-> 0 = unconstrained
				-> other = constrained 
				*/
	      		constrained_b[j]=o_0_b[i_o_0_b]; // constraint 0/1 in motion x
	      		constrained_b[j+n_b]=o_0_b[i_o_0_b+n_o_0_b]; // constraint 0/1 in motion y
	      		constrained_b[j+2*n_b]=o_0_b[i_o_0_b+2*n_o_0_b]; // constraint 0/1 in motion z
	      		constrained_b[j+3*n_b]=o_0_b[i_o_0_b+3*n_o_0_b]; // constraint 0/1 in rotation x
	      		constrained_b[j+4*n_b]=o_0_b[i_o_0_b+4*n_o_0_b]; // constraint 0/1 in rotation y
	      		constrained_b[j+5*n_b]=o_0_b[i_o_0_b+5*n_o_0_b]; // constraint 0/1 in rotation z

				constrained_b_frame[j]=bounc_frame[i]; // coordinate frame

	      		i_o_0_b++;
				break;
	    	}
	  
		}
    }

  	/* creating unconstrained list of particles 
	n_unc_s=n_s-n_o_0_s-n_r_0_s;
	if (n_unc_s>0)
	{
		unconstrained_s=createMemMore(int,n_unc_s);
		for (i=0;i<n_s;i++)
		{
			if (!constrained_s[i])
			{
				unconstrained_s[i_unc_s]=i;
				i_unc_s++;
			}
		}
	}
	*/

  	/* creating unconstrained list of nodes 
	n_unc_b=n_b-n_o_0_b-n_r_0_b;
	if (n_unc_b>0)
	{
		unconstrained_b=createMemMore(int,n_unc_b);
		for (i=0;i<n_b;i++)
		{
			if (!constrained_b[i])
			{
				unconstrained_b[i_unc_b]=i;
				i_unc_b++;
			}
		}
	}
	*/

	/* renumbering rigid body materials 
	-> not necessary - done during reading
	for (i=0;i<n_r;i++)
	{
		for (j=0;j<n_m;j++)
		{
			if (mat_r[i]==num_m[j]) mat_r[i]=j;
		}
	}
	*/

	/* rigid bodies */
  	for (i=0;i<n_r;i++)
    {
		/* rigid body type */
		if (type_r[i]==0.0) // mass and inertia calculated, nodes rewritten
		{
			//type_r[i]=0;
			strcpy(rb_type,"0 (calculated from mass distribution)");
			// mass m_r[i] is 0 for automatic calculation
		}
		else // mass, inertia and nodes given
		{
			//type_r[i]=1;
			strcpy(rb_type,"1 (given by user)");
			// mass m_r[i] is set
		}

		fprintf(stdout,
			"rigid body number %d"
			"\n"
			" -> type = %s"
			"\n"
			" -> material = %d"
			"\n"
			" -> centre of gravity node = %d"
			"\n"
			" -> first principal axis node = %d"
			"\n"
			" -> second principal axis node = %d"
			"\n"
			" -> third principal axis node = %d"
			"\n"
			"\n",
			num_r[i],rb_type,mat_r[i],
			COG_r[i],N1_r[i],N2_r[i],N3_r[i]);

		/* test on rigid body material */
		for (j=0;j<n_m;j++)
		{
			if (mat_r[i]==num_m[j])
			{
				imate=1;
				break;
			}
		}
		if (!imate)
		{
			fprintf(stdout,
				"ERROR: material %d in rigid body %d is not defined!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				num_m[j],num_r[i]);
			return(0);
		}
		else if (type_m[j]!=0) // test on rigid body material type
		{
			fprintf(stdout,
				"ERROR: material type %d in material %d of rigid body %d must be zero!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				type_m[j],num_m[j],num_r[i]);
			return(0);
		}
		imate=0;
	
		/* user defined mass and inertia */
		if ((data_check)&&(type_r[i]))
		{
			if ((dim>1)&&(I1_r[i]<=0.0)) // first principal moment of inertia
			{
				fprintf(stdout,
					"ERROR: first principal moment of inertia of rigid material %d must be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((dim>2)&&(I2_r[i]<=0.0)) // second principal moment of inertia
			{
				fprintf(stdout,
					"ERROR: second principal moment of inertia of rigid material %d must be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
			if ((dim>2)&&(I3_r[i]<=0.0)) // third principal moment of inertia
			{
				fprintf(stdout,
					"ERROR: third principal moment of inertia of rigid material %d must be a positive value!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_m[i]);
				return(0);
			}
		}

		/* particles and nodes affiliation to rigid body 
		-> constrained_s_rb = particle affiliation to rigid body (0 / rigid body number)
		-> constrained_b_rb = node affiliation to rigid body (0 / rigid body number) */
		
		/* counting particles per rigid body */
		for (j=0;j<n_s;j++)
		{
			/* particles affiliation to rigid body */
			if (num_m[mat_s[j]]==mat_r[i])
			{
				constrained_s_rb[j]=num_r[i];

				/* mass and inertia calculated */
				if (!type_r[i])
				{
					m_r[i]=m_r[i]+m_s[j];
					
					/* static moments for centre of gravity */
					if (dim>1) x_r[i]=x_r[i]+m_s[j]*x_s[j];
					if (dim>2)
					{
						x_r[i+n_r]=x_r[i+n_r]+m_s[j]*x_s[j+n_s];
						x_r[i+2*n_r]=x_r[i+2*n_r]+m_s[j]*x_s[j+2*n_s];
					}
				}
				/* increase number of particles in i-th rigid body
				and number of particles in all rigid bodies n_r_0_b */
				rb_s[i]++; // number of particles in rigid body i
				n_r_0_s++; // number of particles in all rigid bodies
			}
		}

		/* elements affiliation to rigid body */
		memset(m_b,0,n_b*sizeof(double));
		for (j=0;j<n_e;j++)
		{
			/* nodes affiliation to rigid body */
			if (num_m[mat_e[j]]==mat_r[i])
			{
				switch (dim)
				{
					case 1:
					{
						constrained_b_rb[nod_e[j]]=num_r[i];
						constrained_b_rb[nod_e[j+n_e]]=num_r[i];

						mat_b[nod_e[j]]=mat_e[j];
						mat_b[nod_e[j+n_e]]=mat_e[j];

						if (!type_r[i]) // mass and inertia calculated (negative density)
						{
							m_b[nod_e[j]]=m_b[nod_e[j]]+rho_m[mat_e[j]]*V_e[j]/2.0;
							m_b[nod_e[j+n_e]]=m_b[nod_e[j+n_e]]+rho_m[mat_e[j]]*V_e[j]/2.0;							
						}

						break;
					}
					case 2:
					{
						if ((nod_e[j+2*n_e]==nod_e[j+3*n_e])&&
							(nod_e[j+3*n_e]==nod_e[j+4*n_e])) // bar or beam
						{
							constrained_b_rb[nod_e[j]]=num_r[i];
							constrained_b_rb[nod_e[j+n_e]]=num_r[i];

							mat_b[nod_e[j]]=mat_e[j];
							mat_b[nod_e[j+n_e]]=mat_e[j];

							if (!type_r[i]) // mass and inertia calculated (negative density)
							{
								m_b[nod_e[j]]=m_b[nod_e[j]]+rho_m[mat_e[j]]*V_e[j]/2.0;
								m_b[nod_e[j+n_e]]=m_b[nod_e[j+n_e]]+rho_m[mat_e[j]]*V_e[j]/2.0;
							}
						}
						else if (nod_e[j+3*n_e]==nod_e[j+4*n_e]) // triangle
						{
							constrained_b_rb[nod_e[j]]=num_r[i];
							constrained_b_rb[nod_e[j+n_e]]=num_r[i];
							constrained_b_rb[nod_e[j+2*n_e]]=num_r[i];

							mat_b[nod_e[j]]=mat_e[j];
							mat_b[nod_e[j+n_e]]=mat_e[j];
							mat_b[nod_e[j+2*n_e]]=mat_e[j];

							if (!type_r[i]) // mass and inertia calculated (negative density)
							{
								m_b[nod_e[j]]=m_b[nod_e[j]]+rho_m[mat_e[j]]*V_e[j]/3.0;
								m_b[nod_e[j+n_e]]=m_b[nod_e[j+n_e]]+rho_m[mat_e[j]]*V_e[j]/3.0;
								m_b[nod_e[j+2*n_e]]=m_b[nod_e[j+2*n_e]]+rho_m[mat_e[j]]*V_e[j]/3.0;
							}
						}
						else // quad
						{
							constrained_b_rb[nod_e[j]]=num_r[i];
							constrained_b_rb[nod_e[j+n_e]]=num_r[i];
							constrained_b_rb[nod_e[j+2*n_e]]=num_r[i];
							constrained_b_rb[nod_e[j+3*n_e]]=num_r[i];

							mat_b[nod_e[j]]=mat_e[j];
							mat_b[nod_e[j+n_e]]=mat_e[j];
							mat_b[nod_e[j+2*n_e]]=mat_e[j];
							mat_b[nod_e[j+3*n_e]]=mat_e[j];

							if (!type_r[i]) // mass and inertia calculated (negative density)
							{
								m_b[nod_e[j]]=m_b[nod_e[j]]+rho_m[mat_e[j]]*V_e[j]/4.0;
								m_b[nod_e[j+n_e]]=m_b[nod_e[j+n_e]]+rho_m[mat_e[j]]*V_e[j]/4.0;
								m_b[nod_e[j+2*n_e]]=m_b[nod_e[j+2*n_e]]+rho_m[mat_e[j]]*V_e[j]/4.0;
								m_b[nod_e[j+3*n_e]]=m_b[nod_e[j+3*n_e]]+rho_m[mat_e[j]]*V_e[j]/4.0;
							}
						}
						break;
					}
					case 3:
					{
						if ((nod_e[j+2*n_e]==nod_e[j+3*n_e])&&
							(nod_e[j+3*n_e]==nod_e[j+4*n_e])) // bar or beam
						{
							constrained_b_rb[nod_e[j]]=num_r[i];
							constrained_b_rb[nod_e[j+n_e]]=num_r[i];

							mat_b[nod_e[j]]=mat_e[j];
							mat_b[nod_e[j+n_e]]=mat_e[j];

							if (!type_r[i]) // mass and inertia calculated (negative density)
							{
								m_b[nod_e[j]]=m_b[nod_e[j]]+rho_m[mat_e[j]]*V_e[j]/2.0;
								m_b[nod_e[j+n_e]]=m_b[nod_e[j+n_e]]+rho_m[mat_e[j]]*V_e[j]/2.0;
							}
						}
						else if (nod_e[j+3*n_e]==nod_e[j+4*n_e]) // triangle
						{
							constrained_b_rb[nod_e[j]]=num_r[i];
							constrained_b_rb[nod_e[j+n_e]]=num_r[i];
							constrained_b_rb[nod_e[j+2*n_e]]=num_r[i];

							mat_b[nod_e[j]]=mat_e[j];
							mat_b[nod_e[j+n_e]]=mat_e[j];
							mat_b[nod_e[j+2*n_e]]=mat_e[j];

							if (!type_r[i]) // mass and inertia calculated (negative density)
							{
								m_b[nod_e[j]]=m_b[nod_e[j]]+rho_m[mat_e[j]]*V_e[j]/3.0;
								m_b[nod_e[j+n_e]]=m_b[nod_e[j+n_e]]+rho_m[mat_e[j]]*V_e[j]/3.0;
								m_b[nod_e[j+2*n_e]]=m_b[nod_e[j+2*n_e]]+rho_m[mat_e[j]]*V_e[j]/3.0;
							}
						}
						else // quad or tetrahedron
						{
							constrained_b_rb[nod_e[j]]=num_r[i];
							constrained_b_rb[nod_e[j+n_e]]=num_r[i];
							constrained_b_rb[nod_e[j+2*n_e]]=num_r[i];
							constrained_b_rb[nod_e[j+3*n_e]]=num_r[i];

							mat_b[nod_e[j]]=mat_e[j];
							mat_b[nod_e[j+n_e]]=mat_e[j];
							mat_b[nod_e[j+2*n_e]]=mat_e[j];
							mat_b[nod_e[j+3*n_e]]=mat_e[j];

							if (!type_r[i]) // mass and inertia calculated (negative density)
							{
								m_b[nod_e[j]]=m_b[nod_e[j]]+rho_m[mat_e[j]]*V_e[j]/4.0;
								m_b[nod_e[j+n_e]]=m_b[nod_e[j+n_e]]+rho_m[mat_e[j]]*V_e[j]/4.0;
								m_b[nod_e[j+2*n_e]]=m_b[nod_e[j+2*n_e]]+rho_m[mat_e[j]]*V_e[j]/4.0;
								m_b[nod_e[j+3*n_e]]=m_b[nod_e[j+3*n_e]]+rho_m[mat_e[j]]*V_e[j]/4.0;
							}
						}
						break;
					}
				}
			}
		}

		/* counting nodes per rigid body */
		for (j=0;j<n_b;j++)
		{
			if (constrained_b_rb[j]==num_r[i])
			{
				/* mass and inertia calculated */
				if (!type_r[i])
				{
					m_r[i]=m_r[i]+m_b[j];

					/* static moments for centre of gravity */
					if (dim>1) x_r[i]=x_r[i]+m_b[j]*x_b[j];
					if (dim>2)
					{
						x_r[i+n_r]=x_r[i+n_r]+m_b[j]*x_b[j+n_b];
						x_r[i+2*n_r]=x_r[i+2*n_r]+m_b[j]*x_b[j+2*n_b];
					}
				}

				/* increase number of nodes in i-th rigid body
				and number of nodes in all rigid bodies n_r_0_b */
				rb_b[i]++; // number of nodes in i-th rigid body
				n_r_0_b++; // number of nodes in all rigid bodies
			}
		}

		/* centre of gravity */
		if (dim>1) x_r[i]=x_r[i]/m_r[i];
		if (dim>2)
		{
			x_r[i+n_r]=x_r[i+n_r]/m_r[i];
			x_r[i+2*n_r]=x_r[i+2*n_r]/m_r[i];
		}
	}

	/* counting particles and nodes per rigid body */
	if (n_r_0_s>0) ind_r_0_s=createMemMore(int,n_r_0_s); // n_r_0_s -> number of particles in rigid body
	if (n_r_0_b>0) ind_r_0_b=createMemMore(int,n_r_0_b); // n_r_0_b -> number of nodes in rigid body

  	/* distributing nodes per rigid bodies */
  	for (i=0;i<n_r;i++)
    {
		/*
		fprintf(stdout,
			" -> %d particles\n",
			rb_s[i]);
		*/

		for (j=0;j<n_s;j++) // particles in all rigid bodies
		{
			if (constrained_s_rb[j]==num_r[i]) // particle in rigid body
			{
				//fprintf(stdout,"     -> particle %d\n",num_s[j]);
				ind_r_0_s[i_r_0_s]=j;
				i_r_0_s++;

				/* moments of inertia and deviation moments from particles */
				if (!type_r[i]) // only for rigid body type 0
				{
					if (dim>1) I1_r[i]=I1_r[i]+m_s[j]*(sqr(x_s[j]-x_r[i])+sqr(x_s[j+n_s]-x_r[i+n_r]));
					if (dim>2)
					{
						I2_r[i]=I2_r[i]+m_s[j]*(sqr(x_s[j]-x_r[i])+sqr(x_s[j+2*n_s]-x_r[i+2*n_r]));
						I3_r[i]=I3_r[i]+m_s[j]*(sqr(x_s[j]-x_r[i])+sqr(x_s[j+n_s]-x_r[i+n_r]));

						D1_r[i]=D1_r[i]+m_s[j]*(x_s[j]-x_r[i])*(x_s[j+n_s]-x_r[i+n_r]);
						D2_r[i]=D2_r[i]+m_s[j]*(x_s[j]-x_r[i])*(x_s[j+2*n_s]-x_r[i+2*n_r]);
						D3_r[i]=D3_r[i]+m_s[j]*(x_s[j+n_s]-x_r[i+n_r])*(x_s[j+2*n_s]-x_r[i+2*n_r]);
					}

					/* new inertia axes coincident with global coordinate system */
					if (N1_r[i]==num_s[j]) // first principal axis particle
					{
						x_s[j]=x_r[i]+1.0;
						x_s[j+n_s]=x_r[i+n_r];
						x_s[j+2*n_s]=x_r[i+2*n_r];
					}
					if (N2_r[i]==num_s[j]) // second principal axis particle
					{
						x_s[j]=x_r[i];
						x_s[j+n_s]=x_r[i+n_r]+1.0;
						x_s[j+2*n_s]=x_r[i+2*n_r];
					}
					if (N3_r[i]==num_s[j]) // third principal axis particle
					{
						x_s[j]=x_r[i];
						x_s[j+n_s]=x_r[i+n_r];
						x_s[j+2*n_s]=x_r[i+2*n_r]+1.0;
					}
				}
			}
		}

		/*
		fprintf(stdout,
			" -> %d nodes\n",
			rb_b[i]);
		*/
			
		for (j=0;j<n_b;j++) // nodes in all rigid bodies
		{
			if (constrained_b_rb[j]==num_r[i]) // node in rigid body
			{
				//fprintf(stdout,"     -> node %d\n",num_b[j]);
				ind_r_0_b[i_r_0_b]=j;
				i_r_0_b++;
				
				/* moments of inertia and deviation moments from nodes */
				if (!type_r[i]) // only for rigid body type 0
				{
					if (dim>1) I1_r[i]=I1_r[i]+m_b[j]*(sqr(x_b[j])+sqr(x_b[j+n_b]));
					if (dim>2)
					{
						I2_r[i]=I2_r[i]+m_b[j]*(sqr(x_b[j]-x_r[i])+sqr(x_b[j+2*n_b]-x_r[i+2*n_r]));
						I3_r[i]=I3_r[i]+m_b[j]*(sqr(x_b[j]-x_r[i])+sqr(x_b[j+n_b]-x_r[i+n_r]));

						D1_r[i]=D1_r[i]+m_b[j]*(x_b[j]-x_r[i])*(x_b[j+n_b]-x_r[i+n_r]);
						D2_r[i]=D2_r[i]+m_b[j]*(x_b[j]-x_r[i])*(x_b[j+2*n_b]-x_r[i+2*n_r]);
						D3_r[i]=D3_r[i]+m_b[j]*(x_b[j+n_b]-x_r[i+n_r])*(x_b[j+2*n_b]-x_r[i+2*n_r]);
					}

					/* new inertia axes coincident with global coordinate system */
					if (N1_r[i]==num_b[j]) // first principal axis node
					{
						x_b[j]=x_r[i]+1.0;
						x_b[j+n_b]=x_r[i+n_r];
						x_b[j+2*n_b]=x_r[i+2*n_r];
					}
					if (N2_r[i]==num_b[j]) // second principal axis node
					{
						x_b[j]=x_r[i];
						x_b[j+n_b]=x_r[i+n_r]+1.0;
						x_b[j+2*n_b]=x_r[i+2*n_r];
					}
					if (N3_r[i]==num_b[j]) // third principal axis node
					{
						x_b[j]=x_r[i];
						x_b[j+n_b]=x_r[i+n_r];
						x_b[j+2*n_b]=x_r[i+2*n_r]+1.0;
					}
				}
			}
		}

		/* number of particles, nodes, mass and inertia */
		fprintf(stdout,
			" -> %d particles"
			"\n"
			" -> %d nodes"
			"\n"
			" -> mass = %f"
			"\n"
			" -> first principal moment of inertia = %f"
			"\n"
			" -> second principal moment of inertia = %f"
			"\n"
			" -> third principal moment of inertia = %f"
			"\n"
			" -> first deviation moment of inertia = %f"
			"\n"
			" -> second deviation moment of inertia = %f"
			"\n"
			" -> third deviation moment of inertia = %f"
			"\n"
			"\n",
			rb_s[i],rb_b[i],m_r[i],I1_r[i],I2_r[i],I3_r[i],D1_r[i],D2_r[i],D3_r[i]);
	}
	
	/* print list of particles and nodes for each rigid body 
	//count_rb_s=0; // counter of number of particles in rigid body
	//count_rb_b=0; // counter of number of nodes in rigid body
	for (i=0;i<n_r;i++)
	{
		fprintf(stdout,
			"rigid body number %d"
			"\n",
			num_r[i]);

		fprintf(stdout," -> %d particles\n",rb_s[i]);
		for (j=0;j<rb_s[i];j++) // particles
		{
			fprintf(stdout,"     -> particle %d\n",num_s[ind_r_0_s[j+count_rb_s]]);
		}
		count_rb_s=count_rb_s+rb_s[i];
		
		fprintf(stdout," -> %d nodes\n",rb_b[i]);
		for (j=0;j<rb_b[i];j++) // nodes
		{
			fprintf(stdout,"     -> node %d\n",num_b[ind_r_0_b[j+count_rb_b]]);
		}
		count_rb_b=count_rb_b+rb_b[i];
	}
	*/

	/* test on coincidence with initial and boundary conditions 
	-> INPRE p_s (including COG)
	-> INVEL v_s, v_b
	-> ACFLD a_0_s, a_0_b
	-> FORCE f_0_s, f_0_b
	-> NDAMP d_0_b
	-> BOUNC o_0_s, o_0_b */
	/* set in declarations */
	count_rb_s=0; // counter of number of particles in rigid body
	count_rb_b=0; // counter of number of nodes in rigid body
	for (i=0;i<n_r;i++)
	{
		/* particles in rigid body */
		for (j=0;j<rb_s[i];j++)
		{
			 /* rigid body particle index */
			 k=ind_r_0_s[j+count_rb_s];

			/* initial pressure coincident with rigid body particle 
			if (p_s[k]>0.0)
			{
				fprintf(stdout,
					"\n"
					"ERROR: initial pressure on particle %d is coincident with rigid body %d particle %d!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					inpre[k],num_r[i],num_s[k]);
				return(0);
			} */
			for (l=0;l<n_p;l++)
			{
				if (inpre[l]==k)
				{
					fprintf(stdout,
						"\n"
						"ERROR: initial pressure on particle %d is coincident with rigid body %d particle %d!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						inpre[k],num_r[i],num_s[k]);
					return(0);
				}
			}

			/* only centre of gravity allowed */
			if (k!=COG_r[i]) // no center of gravity
			{
				/* initial velocity coincident with rigid body particle */
				for (l=0;l<n_v;l++)
				{
					if ((invel[l]==k)&& // initial velocity on rigid body particle
						((norm(v_0[l],v_0[l+n_v],v_0[l+2*n_v])>0.0)|| // translational velocity
						(norm(omega_0[l],omega_0[l+n_v],omega_0[l+2*n_v])>0.0))) // rotational velocity
					{
						fprintf(stdout,
							"\n"
							"ERROR: initial velocity on particle %d is coincident with rigid body %d particle %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							invel[l],num_r[i],num_s[ind_r_0_s[j+count_rb_s]]);
						return(0);
					}
				}

				/* particles acceleration field */
				for (k=0;k<n_a_0_s;k++)
				{
					/* acceleration field coincident with rigid body particle */
					if ((acfld[k]==k)&& // acceleration field on rigid body particle
						(norm(a_0_s[k],a_0_s[k+n_a_0_s],a_0_s[k+2*n_a_0_s])>0.0))
					{
						fprintf(stdout,
							"\n"
							"ERROR: acceleration field on particle %d is coincident with rigid body %d particle %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							acfld[k],num_r[i],num_s[ind_r_0_s[j+count_rb_s]]);
						return(0);
					}
				}

				/* particles force */
				for (k=0;k<n_f_0_s;k++)
				{
					/* force coincident with rigid body particle */
					if ((force[k]==k)&& // force and moment on rigid body particle
						((norm(f_0_s[k],f_0_s[k+n_f_0_s],f_0_s[k+2*n_f_0_s])>0.0)|| // force
						(norm(M_0_s[k],M_0_s[k+n_f_0_s],M_0_s[k+2*n_f_0_s])>0.0))) // moment
					{
						fprintf(stdout,
							"\n"
							"ERROR: force on particle %d is coincident with rigid body %d particle %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							force[k],num_r[i],num_s[ind_r_0_s[j+count_rb_s]]);
						return(0);
					}
				}

				/* particles damping
				-> not for particles */

				/* particles boundary conditions */
				for (k=0;k<n_o_0_s;k++)
				{
					/* boundary condition coincident with rigid body particle */
					if ((bounc[k]==k)&& // boundary condition on rigid body particle
						(num_r[i]==constrained_s_rb[ind_o_0_s[k]])&& // particle affiliated to rigid body
						(ind_o_0_s[k]==ind_r_0_s[j+count_rb_s])) // boundary condition on particle affiliated to rigid body
					{
						fprintf(stdout,
							"\n"
							"ERROR: boundary condition on particle %d is coincident with rigid body %d particle %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							bounc[k],num_r[i],num_s[ind_r_0_s[j+count_rb_s]]);
						return(0);
					}
				}
			}
		}
		count_rb_s=count_rb_s+rb_s[i];

		/* nodes in rigid body */
		for (j=0;j<rb_b[i];j++)
		{
			 /* rigid body node index */
			 k=ind_r_0_b[j+count_rb_b];

			/* nodes initial pressure 
			-> not for nodes */

			/* only centre of gravity allowed*/
			if (k!=COG_r[i]) // no center of gravity
			{
				/* initial velocity coincident with rigid body node */
				for (k=0;k<n_v;k++)
				{
					if ((invel[k]==k)&& // initial velocity on rigid body node
						((norm(v_0[k],v_0[k+n_v],v_0[k+2*n_v])>0.0)|| // translational velocity
						(norm(omega_0[k],omega_0[k+n_v],omega_0[k+2*n_v])>0.0))) // rotational velocity
					{
						fprintf(stdout,
							"\n"
							"ERROR: initial velocity on node %d is coincident with rigid body %d node %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							invel[k],num_r[i],num_b[ind_r_0_b[j+count_rb_b]]);
						return(0);
					}
				}

				/* nodal acceleration field */
				for (k=0;k<n_a_0_b;k++)
				{
					/* acceleration field coincident with rigid body node */
					if ((acfld[k]==k)&& // acceleration field on rigid body node
						(norm(a_0_b[k],a_0_b[k+n_a_0_b],a_0_b[k+2*n_a_0_b])>0.0))
					{
						fprintf(stdout,
							"\n"
							"ERROR: acceleration field on node %d is coincident with rigid body %d node %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							acfld[k],num_r[i],num_b[ind_r_0_b[j+count_rb_b]]);
						return(0);
					}
				}

				/* nodal force */
				for (k=0;k<n_f_0_b;k++)
				{
					/* force coincident with rigid body node */
					if ((force[k]==k)&& // force and moment on rigid body node
						((norm(f_0_b[k],f_0_b[k+n_f_0_b],f_0_b[k+2*n_f_0_b])>0.0)|| // force
						(norm(M_0_b[k],M_0_b[k+n_f_0_b],M_0_b[k+2*n_f_0_b])>0.0))) // moment
					{
						fprintf(stdout,
							"\n"
							"ERROR: force on node %d is coincident with rigid body %d node %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							force[k],num_r[i],num_b[ind_r_0_b[j+count_rb_b]]);
						return(0);
					}
				}

				/* nodal damping */
				for (k=0;k<n_d_0_b;k++)
				{
					/* damping coincident with rigid body node */
					if ((ndamp[k]==k)&& // damping on rigid body node
						(norm(d_0_b[k],d_0_b[k+n_d_0_b],d_0_b[k+2*n_d_0_b])>0.0))
					{
						fprintf(stdout,
							"\n"
							"ERROR: damping on node %d is coincident with rigid body %d node %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							ndamp[k],num_r[i],num_b[ind_r_0_b[j+count_rb_b]]);
						return(0);
					}
				}

				/* nodal boundary conditions */
				for (k=0;k<n_o_0_b;k++)
				{
					/* boundary condition coincident with rigid body node */
					if ((bounc[k]!=COG_r[i])&& // boundary condition on rigid body node
						(num_r[i]==constrained_b_rb[ind_o_0_b[k]])&& // node affiliated to rigid body
						(ind_o_0_b[k]==ind_r_0_b[j+count_rb_b])) // boundary condition on node affiliated to rigid body
					{
						fprintf(stdout,
							"\n"
							"ERROR: boundary condition on node %d is coincident with rigid body %d node %d!"
							"\n"
							"\n"
							"error termination"
							"\n"
							"\n",
							bounc[k],num_r[j],num_b[ind_r_0_b[j+count_rb_b]]);
						return(0);
					}
				}
			}
		}
		count_rb_b=count_rb_b+rb_b[i];
	}
		
  	/* rigid bodies
	-> COG, N1, N2, N3 must be defined in advance otherwise boundary conditions cannot be applied
	-> COG, N1, N2, N3 must be particles or element nodes as free nodes are removed 
	l_s ... particles coordinates in rigid body local coordinate system
	l_b ... nodes coordinates in rigid body local coordinate system */
	if (n_r_0_s>0) l_s=createMemMore(double,3*n_r_0_s);
	if (n_r_0_b>0) l_b=createMemMore(double,3*n_r_0_b);
	
	/* counters of nodes and particles position */
	count_rb_s=0; // counter of number of particles in rigid body
	count_rb_b=0; // counter of number of nodes in rigid body
	for (i=0;i<n_r;i++)
	{
		/* test of COG existence among particles */
		inode=0;
		for (j=0;j<n_s;j++)
		{
			if (COG_r[i]==num_s[j])
			{
				COG_r[i]=j; // renumbering
				COG_r_type[i]=1; // particle
					
				constrained_r[i]=constrained_s[j];
				constrained_r[i+n_r]=constrained_s[j+n_s];
				constrained_r[i+2*n_r]=constrained_s[j+2*n_s];
				constrained_r[i+3*n_r]=constrained_s[j+3*n_s];
				constrained_r[i+4*n_r]=constrained_s[j+4*n_s];
				constrained_r[i+5*n_r]=constrained_s[j+5*n_s];
				constrained_r_frame[i]=constrained_s_frame[j];

				if (type_r[i]) // rigid body type 1 (for 0 it is calculated)
				{
					x_r[i]=x_s[j];
					x_r[i+n_r]=x_s[j+n_s];
					x_r[i+2*n_r]=x_s[j+2*n_s];
					x_r_0[i]=x_r[i];
					x_r_0[i+n_r]=x_r[i+n_r];
					x_r_0[i+2*n_r]=x_r[i+2*n_r];
				}

				psi_r[i]=0.0;
				psi_r[i+n_r]=0.0;
				psi_r[i+2*n_r]=0.0;
				psi_r_0[i]=psi_r[i];
				psi_r_0[i+n_r]=psi_r[i+n_r];
				psi_r_0[i+2*n_r]=psi_r[i+2*n_r];

				/* translational velocities 
				v_r[i]=v_s[j];
				v_r[i+n_r]=v_s[j+n_s];
				v_r[i+2*n_r]=v_s[j+2*n_s];
				*/
				
				/* velocities */
				for (k=0;k<n_v;k++)
				{
					if (invel[k]==num_s[j])
					{
						v_r[i]=v_0[k];
						v_r[i+n_r]=v_0[k+n_v];
						v_r[i+2*n_r]=v_0[k+2*n_v];
						o_r[i]=omega_0[k];
						o_r[i+n_r]=omega_0[k+n_v];
						o_r[i+2*n_r]=omega_0[k+2*n_v];

						if (invel_frame[k]) loc=invel_frame[k];
						else loc=0;
						
						break;
					}
				}

				/* acceleration field */
				for (k=0;k<n_a_0_s;k++)
				{
					if (ind_a_0_s[k]==num_s[j])
					{
						type_a_0_r[i]=type_a_0_s[k];
						frame_a_0_r[i]=frame_a_0_s[k];
						a_0_r[i]=a_0_s[k];
						a_0_r[i+n_r]=a_0_s[k+n_a_0_s];
						a_0_r[i+2*n_r]=a_0_s[k+2*n_a_0_s];
						break;
					}
				}

				/* concentrated forces and moments */
				for (k=0;k<n_f_0_s;k++)
				{
					if (ind_f_0_s[k]==num_s[j])
					{
						type_f_0_r[i]=type_f_0_s[k];
						frame_f_0_r[i]=frame_f_0_s[k];
						f_0_r[i]=f_0_s[k];
						f_0_r[i+n_r]=f_0_s[k+n_f_0_s];
						f_0_r[i+2*n_r]=f_0_s[k+2*n_f_0_s];
						M_0_r[i]=M_0_s[k];
						M_0_r[i+n_r]=M_0_s[k+n_f_0_s];
						M_0_r[i+2*n_r]=M_0_s[k+2*n_f_0_s];
						break;
					}
				}

				inode=1;
				break;
			}
		}
		if (!inode)
		{
			/* test of COG existence among nodes */
			for (j=0;j<n_b;j++)
			{
				if (COG_r[i]==num_b[j])
				{
					COG_r[i]=j; // renumbering
					COG_r_type[i]=0; // node
					
					constrained_r[i]=constrained_b[j];
					constrained_r[i+n_r]=constrained_b[j+n_b];
					constrained_r[i+2*n_r]=constrained_b[j+2*n_b];
					constrained_r[i+3*n_r]=constrained_b[j+3*n_b];
					constrained_r[i+4*n_r]=constrained_b[j+4*n_b];
					constrained_r[i+5*n_r]=constrained_b[j+5*n_b];
					constrained_r_frame[i]=constrained_b_frame[j];

					if (type_r[i]) // rigid body type 1 (for 0 it is calculated)
					{
						x_r[i]=x_b[j];
						x_r[i+n_r]=x_b[j+n_b];
						x_r[i+2*n_r]=x_b[j+2*n_b];
						x_r_0[i]=x_r[i];
						x_r_0[i+n_r]=x_r[i+n_r];
						x_r_0[i+2*n_r]=x_r[i+2*n_r];
					}

					psi_r[i]=0.0;
					psi_r[i+n_r]=0.0;
					psi_r[i+2*n_r]=0.0;
					psi_r_0[i]=psi_r[i];
					psi_r_0[i+n_r]=psi_r[i+n_r];
					psi_r_0[i+2*n_r]=psi_r[i+2*n_r];

					/* translational velocities 
					v_r[i]=v_b[j];
					v_r[i+n_r]=v_b[j+n_b];
					v_r[i+2*n_r]=v_b[j+2*n_b];
					*/

					/* velocities */
					for (k=0;k<n_v;k++)
					{
						if (invel[k]==num_b[j])
						{
							v_r[i]=v_0[k];
							v_r[i+n_r]=v_0[k+n_v];
							v_r[i+2*n_r]=v_0[k+2*n_v];
							o_r[i]=omega_0[k];
							o_r[i+n_r]=omega_0[k+n_v];
							o_r[i+2*n_r]=omega_0[k+2*n_v];

							if (invel_frame[k]) loc=1;
							else loc=0;

							break;
						}
					}
					
					/* acceleration field */
					for (k=0;k<n_a_0_b;k++)
					{
						if (ind_a_0_b[k]==num_b[j])
						{
							type_a_0_r[i]=type_a_0_b[k];
							frame_a_0_r[i]=frame_a_0_b[k];
							a_0_r[i]=a_0_b[k];
							a_0_r[i+n_r]=a_0_b[k+n_a_0_b];
							a_0_r[i+2*n_r]=a_0_b[k+2*n_a_0_b];
							break;
						}
					}

					/* concentrated forces and moments */
					for (k=0;k<n_f_0_b;k++)
					{
						if (ind_f_0_s[k]==num_b[j])
						{
							type_f_0_r[i]=type_f_0_b[k];
							frame_f_0_r[i]=frame_f_0_b[k];
							f_0_r[i]=f_0_b[k];
							f_0_r[i+n_r]=f_0_b[k+n_f_0_b];
							f_0_r[i+2*n_r]=f_0_b[k+2*n_f_0_b];
							M_0_r[i]=M_0_b[k];
							M_0_r[i+n_r]=M_0_b[k+n_f_0_b];
							M_0_r[i+2*n_r]=M_0_b[k+2*n_f_0_b];
							break;
						}
					}

					inode=1;
					break;
				}
			}
			if (!inode)
			{
				fprintf(stdout,
					"ERROR: COG node %d in rigid body %d is not defined!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					COG_r[i],num_r[i]);
				return(0);
			}
		}

		/* test of N1 existence among particles */
		if (dim>1)
		{
			inode=0;
			for (j=0;j<n_s;j++)
			{
				if (N1_r[i]==num_s[j])
				{
					N1_r[i]=j; // renumbering
					N1_r_type[i]=1; // particle

					//x1_r[i]=x_s[j];
					//x1_r[i+n_r]=x_s[j+n_s];
					//x1_r[i+2*n_r]=x_s[j+2*n_s];
					u1_r[i]=x_s[j]-x_r[i];
					u1_r[i+n_r]=x_s[j+n_s]-x_r[i+n_r];
					u1_r[i+2*n_r]=x_s[j+2*n_s]-x_r[i+2*n_r];

					inode=1;
					break;
				}
			}
			if (!inode)
			{
				/* test of N1 existence among nodes */
				for (j=0;j<n_b;j++)
				{
					if (N1_r[i]==num_b[j])
					{
						N1_r[i]=j; // renumbering
						N1_r_type[i]=0; // node

						//x1_r[i]=x_b[j];
						//x1_r[i+n_r]=x_b[j+n_b];
						//x1_r[i+2*n_r]=x_b[j+2*n_b];
						u1_r[i]=x_b[j]-x_r[i];
						u1_r[i+n_r]=x_b[j+n_b]-x_r[i+n_r];
						u1_r[i+2*n_r]=x_b[j+2*n_b]-x_r[i+2*n_r];

						inode=1;
						break;
					}
				}
				if (!inode)
				{
					fprintf(stdout,
						"ERROR: first principal axis node %d in rigid body %d is not defined!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						N1_r[i],num_r[i]);
					return(0);
				}
			}
		}

		/* test of N2 existence among particles */
		if (dim>2)
		{
			inode=0;
			for (j=0;j<n_s;j++)
			{
				if (N2_r[i]==num_s[j])
				{
					N2_r[i]=j; // renumbering
					N2_r_type[i]=1; // particle

					//x2_r[i]=x_s[j];
					//x2_r[i+n_r]=x_s[j+n_s];
					//x2_r[i+2*n_r]=x_s[j+2*n_s];
					u2_r[i]=x_s[j]-x_r[i];
					u2_r[i+n_r]=x_s[j+n_s]-x_r[i+n_r];
					u2_r[i+2*n_r]=x_s[j+2*n_s]-x_r[i+2*n_r];

					inode=1;
					break;
				}
			}
			if (!inode)
			{
				/* test of N2 existence among nodes */
				for (j=0;j<n_b;j++)
				{
					if (N2_r[i]==num_b[j])
					{
						N2_r[i]=j; // renumbering
						N2_r_type[i]=0; // node

						//x2_r[i]=x_b[j];
						//x2_r[i+n_r]=x_b[j+n_b];
						//x2_r[i+2*n_r]=x_b[j+2*n_b];
						u2_r[i]=x_b[j]-x_r[i];
						u2_r[i+n_r]=x_b[j+n_b]-x_r[i+n_r];
						u2_r[i+2*n_r]=x_b[j+2*n_b]-x_r[i+2*n_r];

						inode=1;
						break;
					}
				}
				if (!inode)
				{
					fprintf(stdout,
						"ERROR: second principal axis node %d in rigid body %d is not defined!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						N2_r[i],num_r[i]);
					return(0);
				}
			}
		}

		/* test of N3 existence among particles */
		if (dim>2)
		{
			inode=0;
			for (j=0;j<n_s;j++)
			{
				if (N3_r[i]==num_s[j])
				{
					N3_r[i]=j; // renumbering
					N3_r_type[i]=1; // particle

					//x3_r[i]=x_s[j];
					//x3_r[i+n_r]=x_s[j+n_s];
					//x3_r[i+2*n_r]=x_s[j+2*n_s];
					u3_r[i]=x_s[j]-x_r[i];
					u3_r[i+n_r]=x_s[j+n_s]-x_r[i+n_r];
					u3_r[i+2*n_r]=x_s[j+2*n_s]-x_r[i+2*n_r];

					inode=1;
					break;
				}
			}
			if (!inode)
			{
				/* test of N3 existence among nodes */
				for (j=0;j<n_b;j++)
				{
					if (N3_r[i]==num_b[j])
					{
						N3_r[i]=j; // renumbering
						N3_r_type[i]=0; // node

						//x3_r[i]=x_b[j];
						//x3_r[i+n_r]=x_b[j+n_b];
						//x3_r[i+2*n_r]=x_b[j+2*n_b];
						u3_r[i]=x_b[j]-x_r[i];
						u3_r[i+n_r]=x_b[j+n_b]-x_r[i+n_r];
						u3_r[i+2*n_r]=x_b[j+2*n_b]-x_r[i+2*n_r];

						inode=1;
						break;
					}
				}
				if (!inode)
				{
					fprintf(stdout,
						"ERROR: third principal axis node %d in rigid body %d is not defined!"
						"\n"
						"\n"
						"error termination"
						"\n"
						"\n",
						N3_r[i],num_r[i]);
					return(0);
				}
			}
		}

		/* rigid body local coordinates */ 
		for (j=0;j<rb_s[i];j++) // particles
		{
			k=ind_r_0_s[j+count_rb_s]; // particle index

			l_s[3*count_rb_s+j]=x_s[k]-x_r[i];
			l_s[3*count_rb_s+j+rb_s[i]]=x_s[k+n_s]-x_r[i+n_r];
			l_s[3*count_rb_s+j+2*rb_s[i]]=x_s[k+2*n_s]-x_r[i+2*n_r];
		}
		count_rb_s=count_rb_s+rb_s[i];

		for (j=0;j<rb_b[i];j++) // nodes
		{
			k=ind_r_0_b[j+count_rb_b]; // nodal index

			l_b[3*count_rb_b+j]=x_b[k]-x_r[i];
			l_b[3*count_rb_b+j+rb_b[i]]=x_b[k+n_b]-x_r[i+n_r];
			l_b[3*count_rb_b+j+2*rb_b[i]]=x_b[k+2*n_b]-x_r[i+2*n_r];
		}
		count_rb_b=count_rb_b+rb_b[i];
	}

	/* after renumbering */
	for (i=0;i<n_r;i++)
	{
		/* test on coincidence between centre of gravity and principal axes particles or nodes */
	   	if (COG_r[i]==N1_r[i])
		{
			fprintf(stdout,
				"ERROR: centre of gravity %d of rigid body %d is coincident with first principal axis node %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				COG_r[i],num_r[i],N1_r[i]);
			return(0);
		}
		
		if (COG_r[i]==N2_r[i])
		{
			fprintf(stdout,
				"ERROR: centre of gravity %d of rigid body %d is coincident with second principal axis node %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				COG_r[i],num_r[i],N2_r[i]);
			return(0);
		}
		
		if (COG_r[i]==N3_r[i])
		{
			fprintf(stdout,
				"ERROR: centre of gravity %d of rigid body %d is coincident with third principal axis node %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				COG_r[i],num_r[i],N3_r[i]);
			return(0);
		}
		
		if (N1_r[i]==N2_r[i])
		{
			fprintf(stdout,
				"ERROR: first principal axis node %d of rigid body %d is coincident with second principal axis node %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				N1_r[i],num_r[i],N2_r[i]);
			return(0);
		}
		
		if (N1_r[i]==N3_r[i])
		{
			fprintf(stdout,
				"ERROR: first principal axis node %d of rigid body %d is coincident with third principal axis node %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				N1_r[i],num_r[i],N3_r[i]);
			return(0);
		}
		
		if (N2_r[i]==N3_r[i])
		{
			fprintf(stdout,
				"ERROR: second principal axis node %d of rigid body %d is coincident with third principal axis node %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				N2_r[i],num_r[i],N3_r[i]);
			return(0);
		}

		/* test on COG, N1, N2 and N3 belong to single rigid body */
		if (COG_r_type[i]) // particle
		{
			j=num_s[COG_r[i]];
			k=mat_s[COG_r[i]];
		}
		else // node
		{
			j=num_b[COG_r[i]];
			k=mat_b[COG_r[i]];
		}

		if ((num_m[k]!=mat_r[i])&&(k!=-1)) // -1 for free nodes without material affiliation
		{
			fprintf(stdout,
				"ERROR: centre of gravity node %d of rigid body %d has not rigid body material %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				j,num_r[i],mat_r[i]);
			return(0);
		}
		
		if (N1_r_type[i]) // particle
		{
			j=num_s[N1_r[i]];
			k=mat_s[N1_r[i]];
		}
		else // node
		{
			j=num_b[N1_r[i]];
			k=mat_b[N1_r[i]];
		}

		if ((num_m[k]!=mat_r[i])&&(k!=-1)) // -1 for free nodes without material affiliation
		{
			fprintf(stdout,
				"ERROR: first principal axis node %d of rigid body %d has not rigid body material %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				j,num_r[i],mat_r[i]);
			return(0);
		}

		if (N2_r_type[i]) // particle
		{
			j=num_s[N2_r[i]];
			k=mat_s[N2_r[i]];
		}
		else // node
		{
			j=num_b[N2_r[i]];
			k=mat_b[N2_r[i]];
		}
		
		if ((num_m[k]!=mat_r[i])&&(k!=-1)) // -1 for free nodes without material affiliation
		{
			fprintf(stdout,
				"ERROR: second principal axis node %d of rigid body %d has not rigid body material %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				j,num_r[i],mat_r[i]);
			return(0);
		}

		if (N3_r_type[i]) // particle
		{
			j=num_s[N3_r[i]];
			k=mat_s[N3_r[i]];
		}
		else // node
		{
			j=num_b[N3_r[i]];
			k=mat_b[N3_r[i]];
		}
		
		if ((num_m[k]!=mat_r[i])&&(k!=-1)) // -1 for free nodes without material affiliation
		{
			fprintf(stdout,
				"ERROR: third principal axis node %d of rigid body %d has not rigid body material %d!"
				"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
				j,num_r[i],mat_r[i]);
			return(0);
		}	
	}

	/* check on principal inertia axes */
	for (i=0;i<n_r;i++)
	{
		if (dim>1)
		{
			ux=u1_r[i];uy=u1_r[i+n_r];uz=u1_r[i+2*n_r];
			un=norm(ux,uy,uz);
			ux=ux/un;uy=uy/un;uz=uz/un;
		}

		if (dim>2)
		{
			vx=u2_r[i];vy=u2_r[i+n_r];vz=u2_r[i+2*n_r];
			vn=norm(vx,vy,vz);
			vx=vx/vn;vy=vy/vn;vz=vz/vn;

			wx=u3_r[i];wy=u3_r[i+n_r];wz=u3_r[i+2*n_r];
			wn=norm(wx,wy,wz);
			wx=wx/wn;wy=wy/wn;wz=wz/wn;

			/* test on orthogonality of first and second principal axes 
			if (ux*vx+uy*vy+uz*vz>1e6*EPS)
			{
				fprintf(stdout,
					"ERROR: first and second principal inertia axes of rigid body %d are not orthogonal!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_r[i]);
				return(0);
			}*/

			/* test on orthogonality of first and third principal axes 
			if (ux*wx+uy*wy+uz*wz>1e6*EPS)
			{
				fprintf(stdout,
					"ERROR: first and third principal inertia axes of rigid body %d are not orthogonal!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_r[i]);
				return(0);
			}*/

			/* test on orthogonality of second and third principal axes 
			if (vx*wx+vy*wy+vz*wz>1e6*EPS)
			{
				fprintf(stdout,
					"ERROR: second and third principal inertia axes of rigid body %d are not orthogonal!"
					"\n"
					"\n"
					"error termination"
					"\n"
					"\n",
					num_r[i]);
				return(0);
			}*/
		}

		/* initial velocity frame */
		switch(loc)
		{
			case 0: // translational and rotatinal velocities in global coordinate system
			{
				/* translational velocity in global coordinate system 
				   -> calculating displacements in global coordinate system */

				/* rotational velocity in global coordinate system
				   -> from global to local coordinate system */
				ox=o_r[i];
				oy=o_r[i+n_r];
				oz=o_r[i+2*n_r];

				/* scalar product to local coordinate system vectors */
				o_r[i]=ox*ux+oy*uy+oz*uz;
				o_r[i+n_r]=ox*vx+oy*vy+oz*vz;
				o_r[i+2*n_r]=ox*wx+oy*wy+oz*wz;

				break;
			}
			case 1: // translational velocity in local axis system, rotational velocity in global coordinate system
			{
				/* translational velocity 
				   -> from local to global coordinate system */
				rx=v_r[i];
				ry=v_r[i+n_r];
				rz=v_r[i+2*n_r];

				/* scalar product to global coordinate system vectors */
				v_r[i]=rx*ux+ry*vx+rz*wx;
				v_r[i+n_r]=rx*uy+ry*vy+rz*wy;
				v_r[i+2*n_r]=rx*uz+ry*vz+rz*wz;

				/* rotational velocity in global coordinate system
				   -> from global to local coordinate system */
				ox=o_r[i];
				oy=o_r[i+n_r];
				oz=o_r[i+2*n_r];

				/* scalar product to local coordinate system vectors */
				o_r[i]=ox*ux+oy*uy+oz*uz;
				o_r[i+n_r]=ox*vx+oy*vy+oz*vz;
				o_r[i+2*n_r]=ox*wx+oy*wy+oz*wz;
				   
				break;
			}
			case 2: // translational velocity in global axis system, rotational velocity in local coordinate system
			{
				/* translational velocity in global coordinate system 
				   -> calculating displacements in global coordinate system */

				/* rotational velocity in local coordinate system
				   -> calculating rotations in local coordinate system */

				break;
			}
			case 3: // translational and rotatiobal velocities in local coordinate system
			{
				/* translational velocity 
				   -> from local to global coordinate system */
				rx=v_r[i];
				ry=v_r[i+n_r];
				rz=v_r[i+2*n_r];

				/* scalar product to global coordinate system vectors */
				v_r[i]=rx*ux+ry*vx+rz*wx;
				v_r[i+n_r]=rx*uy+ry*vy+rz*wy;
				v_r[i+2*n_r]=rx*uz+ry*vz+rz*wz;

				/* rotational velocity 
				   -> calculating rotations in local coordinate system */

				break;
			}
		}
	}

	/* freeing variables distributed to rigid bodies */
	if (n_s>0)
	{
		freeMem(constrained_s_frame);
		//freeMem(frame_v_0_s);
	}
	if (n_b>0)
	{
		freeMem(constrained_b_frame);
		//freeMem(frame_v_0_b);
	}
  	if (n_p>0)
    {
		freeMem(inpre);
		freeMem(p_0);
    }
  
  	if (n_v>0)
    {
		freeMem(invel);
		freeMem(invel_frame);
      	freeMem(v_0);
      	freeMem(omega_0);
    }
  	if (n_a>0)
    {
      	freeMem(acfld);
      	freeMem(acfld_type);
		freeMem(acfld_frame);
      	freeMem(a_0);
    }
  	if (n_f>0)
    {
      	freeMem(force);
      	freeMem(force_type);
		freeMem(force_frame);
      	freeMem(f_0);
		freeMem(M_0);
    }
  
  	if (n_d>0)
    {
      	freeMem(ndamp);
      	freeMem(d_0);
    }

  	if (n_o>0)
    {
      	freeMem(bounc);
      	freeMem(bounc_type);
      	freeMem(bounc_frame);
      	freeMem(o_0);
    }

  	/* optimized nearest neighbours search */
	if ((n_s>0)&&(nnopt>0)) // number of neighbours for each particle
	{
		nn_s=createMemMore(int,n_s);
	}

  	sph=fopen(output,"wb");

  	fwrite(&(dim),sizeof(int),1,sph);
  	fwrite(&(ix),sizeof(int),1,sph);

  	fwrite(&(n_m),sizeof(int),1,sph);
  	fwrite(&(n_s),sizeof(int),1,sph);
  	fwrite(&(n_b),sizeof(int),1,sph);
  	fwrite(&(n_e),sizeof(int),1,sph);
  	fwrite(&(n_r),sizeof(int),1,sph);
  	fwrite(&(n_c),sizeof(int),1,sph);

  	fwrite(num_m,sizeof(int),n_m,sph);
  	if (n_s>0) fwrite(num_s,sizeof(int),n_s,sph);
  	if (n_b>0) fwrite(num_b,sizeof(int),n_b,sph);
  	if (n_e>0) fwrite(num_e,sizeof(int),n_e,sph);
  
  	if (n_s>0) fwrite(mat_s,sizeof(int),n_s,sph);
  	if (n_e>0)
	{
		fwrite(mat_e,sizeof(int),n_e,sph);
		fwrite(nod_e,sizeof(int),4*n_e,sph);
	}

  	fwrite(&(save_dt),sizeof(int),1,sph);

  	fwrite(&(save_kine_s),sizeof(int),1,sph);
  	fwrite(&(save_inne_s),sizeof(int),1,sph);
  	fwrite(&(save_pote_s),sizeof(int),1,sph);
  	fwrite(&(save_kine_b),sizeof(int),1,sph);
  	fwrite(&(save_disi_b),sizeof(int),1,sph);
  	fwrite(&(save_defo_b),sizeof(int),1,sph);
  	fwrite(&(save_pote_b),sizeof(int),1,sph);

  	fwrite(&(save_tote),sizeof(int),1,sph);

  	fwrite(&(save_v_s),sizeof(int),1,sph);
  	fwrite(&(save_dv_s),sizeof(int),1,sph);
  	fwrite(&(save_a_s),sizeof(int),1,sph);
  	fwrite(&(save_f_s),sizeof(int),1,sph);

  	fwrite(&(save_rho_s),sizeof(int),1,sph);
  	fwrite(&(save_drhodt_s),sizeof(int),1,sph);
  	fwrite(&(save_u_s),sizeof(int),1,sph);
  	fwrite(&(save_dudt_s),sizeof(int),1,sph);
  	fwrite(&(save_p_s),sizeof(int),1,sph);
  	fwrite(&(save_c_s),sizeof(int),1,sph);
  	fwrite(&(save_h_s),sizeof(int),1,sph);

  	fwrite(&(save_S_s),sizeof(int),1,sph);
  	fwrite(&(save_dSdt_s),sizeof(int),1,sph);
  	fwrite(&(save_e_s),sizeof(int),1,sph);
  	fwrite(&(save_dedt_s),sizeof(int),1,sph);
  	fwrite(&(save_O_s),sizeof(int),1,sph);

  	fwrite(&(save_v_b),sizeof(int),1,sph);
  	fwrite(&(save_a_b),sizeof(int),1,sph);
  	fwrite(&(save_f_b),sizeof(int),1,sph);
  
  	fwrite(&(save_x_r),sizeof(int),1,sph);
  	fwrite(&(save_psi_r),sizeof(int),1,sph);
  	fwrite(&(save_v_r),sizeof(int),1,sph);
  	fwrite(&(save_o_r),sizeof(int),1,sph);
  	fwrite(&(save_a_r),sizeof(int),1,sph);
  	fwrite(&(save_alpha_r),sizeof(int),1,sph);
  	fwrite(&(save_f_r),sizeof(int),1,sph);
  	fwrite(&(save_M_r),sizeof(int),1,sph);

  	fwrite(&(save_f_c),sizeof(int),1,sph);

  	fclose(sph);
  
  	fprintf(stdout,
	  	"solution ..."
	  	"\n"
	  	"\n"
	  	" -> %d functions"
	  	"\n"
	  	" -> %d materials"
	  	"\n"
	  	" -> %d contacts"
	  	"\n"
	  	" -> %d particles"
	  	"\n"
	  	" -> %d nodes"
	  	"\n"
	  	" -> %d elements"
	  	"\n"
	  	" -> %d initial nodal pressures"
	  	"\n"
	  	" -> %d initial nodal velocities"
	  	"\n"
	  	" -> %d nodal acceleration fields"
	  	"\n"
	  	" -> %d nodal forces"
	  	"\n"
	  	" -> %d nodal dampings"
	  	"\n"
	  	" -> %d nodal boundary conditions"
	  	"\n"
	  	" -> %d rigid bodies"
	  	"\n"
	  	"\n",
	  	n_u,n_m,n_c,n_s,n_b,n_e,n_p,n_v,n_a,n_f,n_d,n_o,n_r);

	/* for FEM only */
	if (n_e>0)
	{
	  	/* mass and stiffness matrices in case of small deformations */
  		mass_stiffness(type_m, // materisla
			rho_m,T_m,mu_m,kappa_m,
			n_b, // nodes
			x_b,x_b_0,v_b, // initial values
			n_e,nod_e,mat_e, // elements
			n_n, // boundary equation of motion matrices
			M,K,dim);
	}
	
	/* initial time step */
	dt=dt_init;

  	while (t<t_max)
    {
      	/* LHy for debugging saving time step
      	fprintf(stdout,
			"---------------------------------------"
			"\n"
			"t = %f ..."
			"\n"
			"dt = %f ..."
			"\n"
			"t + dt = %f ..."
			"\n"
			"dt_save = %f ..."
			"\n"
			"Next t_save = %f ... \n",
			t,dt,t+dt,dt_save,t_save);
		*/

      	/* print cycle */
      	if ((cycle_print>0)&&(!(count_cycle%cycle_print)))
		{
	  		if (count_cycle==0)
			{
				fprintf(stdout,
			  		"cycle = %d, time = %f, dt = not set, ready = %6.2f%%"
			  		"\n",
		  			count_cycle,t,100*t/t_max);
			}
			else
			{
				fprintf(stdout,
			  		"cycle = %d, time = %f, dt = %10.9f, ready = %6.2f%%"
			  		"\n",
		  			count_cycle,t,dt,100*t/t_max);
			}
		}

      	/* save configuration */
      	if ((dt_save==0.0)||(dt_save<dt)||((t+dt>t_save)&&(t<=t_save)))
		{
			/*
			fprintf(stdout,
		  		"\n"
		  		" t = %f, t_save = %f, t + dt = %f"
				"\n",
		  		t,t_save,t+dt);
			*/
	    
	  		/* print */
			fprintf(stdout,
		  		"\n"
		  		" -> saving %d. state at time = %f ... ",
		  		count_save,t);
	    
		  	/* energy */
		  	energy(n_s, // particles
			x_s,v_s,
			m_s,u_s,
			n_b, // nodes
			x_b,x_b_0,v_b,
			n_n, // boundary equation of motion materices
			M,B,K,f,
			&kine_s,&inne_s,&pote_s,
			&kine_b,&defo_b,&disi_b,&pote_b,
			&tote,
			ax,ay,az);

		  	/* store to ouptut file */
		  	store(n_s, // particles
				m_s,
				x_s,v_s,dv_s,a_s,f_s,
				rho_s,drhodt_s,u_s,dudt_s,p_s,c_s,h_s,
				S_s,dSdt_s,e_s,dedt_s,O_s,
				n_b, // nodes
				x_b,v_b,a_b,f_b,
				n_r, // rigid bodies
				x_r,v_r,psi_r,o_r,a_r,alpha_r,f_r,M_r,
				n_c, // contacts
				f_c,
				t,dt,
				kine_s,inne_s,pote_s, // energy
				kine_b,pote_b,defo_b,disi_b,
				tote,
				dim,
				ix,
				output, // output file name
				save_dt, // saving switches
				save_kine_s,save_inne_s,save_pote_s,
				save_kine_b,save_disi_b,save_defo_b,save_pote_b,
				save_tote,
				save_v_s,save_dv_s,save_a_s,save_f_s,
				save_rho_s,save_drhodt_s,save_u_s,save_dudt_s,
				save_p_s,save_c_s,save_h_s,
				save_S_s,save_dSdt_s,save_e_s,save_dedt_s,save_O_s,
				save_v_b,save_a_b,save_f_b,
				save_f_c,
				save_x_r,save_psi_r,save_v_r,save_o_r,
				save_a_r,save_alpha_r,save_f_r,save_M_r);

	  		/* update t_save after the last time configuration in case dt_save<dt
	        -> not relevant for dt_save=0.0 (each time step saved) */
	  		if ((dt_save>0.0)&&(dt_save<dt))
	    	{
	      		while (t_save<t)
	        	{
	          		t_save=t_save+dt_save;
				}
	    	}

	  		/* jump over current time */  
	  		t_save=t_save+dt_save;

			/* wait for user to hit enter or another key */
			//system("pause");

			count_save++;
		}

		/* maximum time reached */
		if (t_save-EPS>t_max) // without substracting 1e-9 stops one step earlier
		{
			//fprintf(stdout,"%f,%f\n",t_save,t_max);
			break;
		}		

      	/* calculate time step and store previous values V_o:=V just in one loop over particles */
      	dt=time_step(type_m,
			rho_m,T_m,mu_m,kappa_m,
			n_s,mat_s,
			x_s,v_s, // actual value
			//x_s_a,v_s_a, // value to be updated by predictor
			x_s_o,v_s_o, // value from beginning of time step
			rho_s,u_s,p_s,c_s,h_s, // actual value
			//rho_s_a,u_s_a,//p_s_a,c_s_a,h_s_a,
			rho_s_o,u_s_o,//p_s_o,c_s_o,h_s_o, // value from beginning of time step
			a_s,drhodt_s,dudt_s, // actual value
			a_s_o,drhodt_s_o,dudt_s_o, // value from beginning of time step
			e_s,S_s,dedt_s,dSdt_s, // actual value
			//e_s_a,S_s_a, // value to be updated by predictor
			e_s_o,S_s_o,dedt_s_o,dSdt_s_o, // value from beginning of time step
			n_b,
			x_b,v_b, // actual value
			//x_b_a,v_b_a, // value to be updated by predictor
			x_b_o,v_b_o, // value from beginning of time step
			a_b, // actual value
			a_b_o, // value from beginning of time step
			n_e,
			mat_e,nod_e,
			V_e,
			n_r,
			x_r,psi_r,v_r,o_r, // actual value
			//x_r_a,psi_r_a,v_r_a,o_r_a, // value to be updated by predictor
			x_r_o,psi_r_o,v_r_o,o_r_o, // value from beginning of time step
			a_r,a_r, // actual value
			alpha_r_o,alpha_r_o, // value from beginning of time step
			dim,
			alpha,beta,
			cour,kstab,
			integration,
			&v_max,
			dr,
			dt_max,
			dt);

      	/* termination because of small time step */
      	if (dt<EPS)
		{
			error=1;
			
			fprintf(stdout,
		  		"\n"
		  		"ERROR: time step %f is too small!"
		  		"\n"
				"\n"
		  		"saving error state ..."
		  		"\n"
				"\n"
				"error termination"
				"\n"
				"\n",
		  		dt);

	  		/* store to error file */
	  		store(n_s, // particles
				m_s,
				x_s,v_s,dv_s,a_s,f_s,
				rho_s,drhodt_s,u_s,dudt_s,p_s,c_s,h_s,
				S_s,dSdt_s,e_s,dedt_s,O_s,
				n_b, // nodes
				x_b,v_b,a_b,f_b,
				n_r, // rigid bodies
				x_r,v_r,psi_r,o_r,a_r,alpha_r,f_r,M_r,
				n_c, // contacts
				f_c,
				t,dt,
				kine_s,inne_s,pote_s, // energy
				kine_b,pote_b,defo_b,disi_b,
				tote,
				dim,
				ix,
				output, // output file name
				save_dt, // saving switches
				save_kine_s,save_inne_s,save_pote_s,
				save_kine_b,save_disi_b,save_defo_b,save_pote_b,
				save_tote,
				save_v_s,save_dv_s,save_a_s,save_f_s,
				save_rho_s,save_drhodt_s,save_u_s,save_dudt_s,
				save_p_s,save_c_s,save_h_s,
				save_S_s,save_dSdt_s,save_e_s,save_dedt_s,save_O_s,
				save_v_b,save_a_b,save_f_b,
				save_f_c,
				save_x_r,save_psi_r,save_v_r,save_o_r,
				save_a_r,save_alpha_r,save_f_r,save_M_r);

		  break;
		}

      	/* termination on signal */
		sig=fopen(SIGFILE,"r");
      	if (sig)
		{
			fclose(sig);
			remove(SIGFILE);
			error=2;

			fprintf(stdout,
		  		"WARNING: termination on signal!"
		  		"\n"
				"\n"
		  		"saving signal state"
		  		"\n"
				"\n");

	  		/* store to error file */
	  		store(n_s, // particles
				m_s,
				x_s,v_s,dv_s,a_s,f_s,
				rho_s,drhodt_s,u_s,dudt_s,p_s,c_s,h_s,
				S_s,dSdt_s,e_s,dedt_s,O_s,
				n_b, // nodes
				x_b,v_b,a_b,f_b,
				n_r, // rigid bodies
				x_r,v_r,psi_r,o_r,a_r,alpha_r,f_r,M_r,
				n_c, // contacts
				f_c,
				t,dt,
				kine_s,inne_s,pote_s, // energy
				kine_b,pote_b,defo_b,disi_b,
				tote,
				dim,
				ix,
				output, // output file name
				save_dt, // saving switches
				save_kine_s,save_inne_s,save_pote_s,
				save_kine_b,save_disi_b,save_defo_b,save_pote_b,
				save_tote,
				save_v_s,save_dv_s,save_a_s,save_f_s,
				save_rho_s,save_drhodt_s,save_u_s,save_dudt_s,
				save_p_s,save_c_s,save_h_s,
				save_S_s,save_dSdt_s,save_e_s,save_dedt_s,save_O_s,
				save_v_b,save_a_b,save_f_b,
				save_f_c,
				save_x_r,save_psi_r,save_v_r,save_o_r,
				save_a_r,save_alpha_r,save_f_r,save_M_r);

		  break;
		}

  		/* nearest neighbours search
		-> each nnopt cycles */
  		if ((n_s>0)&&(((nnopt==1)&&(count_cycle==0))||
			((nnopt>1)&&(!(count_cycle%nnopt)))))
		{
			freeMem(neig_s);
			nns(n_s,opt,nn_s,x_s,h_s,&neig_s);
		}

		/* print records
    	int count_neig=0;
    	for (i=0;i<n_s;i++)
		{
			fprintf(stdout,"Particle %d has %d neighbours : ",i,nn_s[i]);
			for (j=0;j<nn_s[i];j++)
			{
				fprintf(stdout,"%d ",neig_s[count_neig+j]);
			}
			fprintf(stdout,"\n");
			count_neig=count_neig+nn_s[i];
		}
		return(0);
		*/

      	/* no predictor step for (0) Euler,
		predictor step (predicted values) for (1) central acceleration,
		(2) predictor-corrector and (3) predictor-corrector leapfrog 
		-> V is actual value of variable (x, v, rho, u, e, S, drho, du, de, dS) to be predicted and corrected
		-> V_o is value of variable stored at beginning of time step interval
			-> for variables updated during time integration (x, v, rho, u, e, S, drho, du, de, dS)
			-> not necessary for p, c, h (directly calculated state variables) */
      	if (integration>0) // central acceleraion (1), predictor-corrector (2), predictor-corrector leapfrog (3)
		{
	  		move(num_m,type_m, // materials
	       		rho_m,T_m,gamma_m,kappa_m,mu_m,
				coef1_m,coef2_m,coef3_m,coef4_m,
				n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u, // functions
		       	n_s,num_s,mat_s, // particles
				m_s,x_s_0,c_s_0, // initial values
				p_s,c_s,h_s, // state values
				dv_s,//O_s,//gradvx_s,gradvy_s,gradvz_s,
				/* --- */
	       		//x_s,v_s, // predicted values x_s and v_s
	       		x_s_o,v_s_o,
				a_s, // by a_s (same as a_s_o)
				//rho_s,u_s, // predicted values rho_s and u_s
				rho_s_o,u_s_o,
				drhodt_s,dudt_s, // by drhodt_s and dudt_s (same as drhodt_s_o and dudt_s_o)
	       		//e_s,S_s, // predited values e_s and S_s
	       		e_s_o,S_s_o,
				dedt_s,dSdt_s, // by dedt_s and dSdt_s (same as dedt_s_o and dSdt_s_o)
				/* --- */
	       		x_s_o,v_s_o, // beginning of time step interval
				a_s_o, // beginning of time step interval
				rho_s_o,u_s_o, // beginning of time step interval
	       		drhodt_s_o,dudt_s_o, // not relevant in predictor
	       		e_s_o,S_s_o, //O_s_o, // beginning of time step interval
				dedt_s_o,dSdt_s_o, // not relevant in predictor
	       		n_b,num_b,x_b_0, // nodes
				/* --- */
	       		//x_b,v_b, // predicted values x_b and v_b
	       		x_b_o,v_b_o,
				a_b, // by a_b (same as a_b_o)
	       		x_b_o,v_b_o, // beginning of time step interval
				a_b_o, // beginning of time step interval
				/* --- */
	       		n_o_0_s,o_0_s,ind_o_0_s,type_o_0_s, // boundary conditions on SPH
				n_o_0_b,o_0_b,ind_o_0_b,type_o_0_b, // boundary conditions on FEM
				n_r,num_r,type_r,COG_r,N1_r,N2_r,N3_r, // rigid bodies
				//COG_r_type,
				N1_r_type,N2_r_type,N3_r_type,
				x_r_0,psi_r_0, // initial values
				/* --- */
				//x_r,v_r,psi_r,o_r, // predicted values x_r, psi_r, v_r and o_r
				x_r_o,v_r_o,psi_r_o,o_r_o,
				a_r,alpha_r, // by a_r and alpha_r (same as a_r_o and alpha_r_o)
				/* --- */
				x_r_o,v_r_o,psi_r_o,o_r_o, // beginning of time step interval
				a_r_o,alpha_r_o, // beginning of time step interval
				//x1_r,x2_r,x3_r,
				//x1_r_o,x2_r_o,x3_r_o,
				u1_r,u2_r,u3_r, // principal moments of inertial vectors
				constrained_s_rb,constrained_b_rb, // particles and nodes constrained in rigid bodies
				rb_s,rb_b, // particles and nodes in rigid bodies
				ind_r_0_s,ind_r_0_b, // boundary conditions on rigid bodies
				//n_r_0_s,n_r_0_b,
				l_s,l_b, // particles and nodes in rigid bodies local coordinate systems
	       		dim, // dimension
	       		ix, // XSPH switch
	       		xeps, // XSPH constant
				constrained_s,constrained_b,constrained_r, // constraint attributes
				constrained_r_frame, // coordinate systems
	       		integration,0, // 0 = predictor (for integration = 1, 2, 3)
	       		t,prei*dt, // time and predictor time step
				v_max, // maximum velocity
	       		sigma, // parameters
	       		h0);
		}
      
      	/* full forces f_s, rhodt_s, dudt_s and f_b
		for (0) Euler, (1) central acceleration
		or (2) predictor-corrector or corrector forces
		for (3) predictor-corrector leapfrog 
		-> or contact between SPH and SPH */

		/* set boundary forces to zero 
		-> before contact */
		if (n_b>0) memset(f_b,0,3*n_b*sizeof(double));
			
		/* set particle forcse to zero 
		-> before contact */
		if (n_s>0) memset(f_s,0,3*n_s*sizeof(double));

      	/* contact forces
			-> between SPH and FEM
			-> between FEM and FEM
		-> each cycle_contact cycles */
  		if (!(count_cycle%cycle_contact))
		{
			contacts(type_m,domain_m, // materials
		       	n_c,f_c, // contacts
		       	seg_c,slave_c,master_c,sw_c,mat_c,
	    	   	m_slave,m_master,ct_c,klin_c,knon_c,kf_c,kd_c,
		       	n_s, // particles
				h_s, // smoothing length
				/**/
		       	//x_s,v_s,f_s, // predicted values -> particles forces
		       	x_s_o,v_s_o,f_s,
				/**/
	    	   	n_b, // nodes
				/**/
	       		//x_b,v_b,f_b, // predicted values -> nodal forces
	       		x_b_o,v_b_o,f_b,
				/**/
	       		n_e,nod_e, // elements
				dt, // time step
		       	dim); // dimension
		}

		/* for SPH only */
		if (n_s>0)
		{
			accelerations(type_m,domain_m, // materials
				rho_m,T_m,mu_m,kappa_m,
				n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u, // functions
				n_c,mat_c,sw_c,seg_c,slave_c,master_c, // contacts
				ct_c,klin_c,knon_c,kf_c,kd_c,
				n_s,mat_s, // particles
				m_s,mu_s,x_s_0,h_s_0, // initial values
				p_s,c_s,h_s, // state values
				dv_s,O_s,
				//gradvx_s,gradvy_s,gradvz_s,
				/**/
				//x_s,v_s, // predicted values -> particles forces and accelerations
				x_s_o,v_s_o,
				a_s,f_s, // to calculate particles forces and accelerations
				//rho_s,u_s, // predicted values -> particles density and internal energy
				rho_s_o,u_s_o, // to calculate particles density and internal energy
				drhodt_s,dudt_s, // calculated particles drho, du
				//S_s, // predicted particles shear stress to calculate pressure
				S_s_o,
				dedt_s,dSdt_s, // calculated particles de, dS
				/**/
				nnopt,nn_s,neig_s, // nearest neightbour search
				n_a_0_s,ind_a_0_s,type_a_0_s, // boundary accelerations
				a_0_s,
				n_f_0_s,ind_f_0_s,type_f_0_s, // boundary forces
				f_0_s,
				dim, // dimension
				is, // viscosity switch
				ix, // XSPH switch
				t,dt, // actual time and time step
				&v_max, // maximum velocity to be updated
				sigma, // parameters
				alpha,beta,eta, // artificial viscosity coefficients
				zeta,nas,theta, // artificial stress coefficients
				c0,c1, // finite element method coefficients
				count_cycle,cycle_contact, // counters
				ax,ay,az); // global acceleration field
		}

	  	/* for FEM only */
		if (n_e>0)
		{
	      	/* small deformations
			-> damping matrix and right-hand side */
    	  	damping_rhs(type_m, // materials
			  	gamma_m,
				n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u, // functions
			  	n_b, // nodes
				/* --- */
			  	//v_b, // predicted nodal velocities
			  	v_b_o,
				f_b, // to calculate nodal forces
				/* --- */
			  	n_e,nod_e,mat_e, // elements
		  		n_n,
			  	M,B,K,f, // boundary equation of motion materices
			  	n_a_0_b,ind_a_0_b,type_a_0_b, // boundary accelerations
			  	a_0_b,
				n_f_0_b,ind_f_0_b,type_f_0_b, // boundary forces
				f_0_b,
		  		n_d_0_b,ind_d_0_b, // boundary conditions
			  	d_0_b,
				t,dt, // time and time step
			  	dim, // dimension
			  	c0,c1, // finite element method coefficients
		  		ax,ay,az); // global acceleration field
			
  			/* conjugate gradient method for solving FEM equation of motion 
			-> M, B, K, f from mass_stiffness.c and damping_rhs.c
			-> initial guess a_b from previous step */
    	  	cgm(n_n,
	  			n_b,num_b, // nodes
				x_b_0, // initial values
				/* --- */
	  			//x_b,v_b, // predicted nodal positions and velocities
	  			x_b_o,v_b_o,
				/* --- */
				a_b, // to calculate nodal accelerations
		  		M,B,K,f, // boundary equation of motion materices
		  		n_o_0_b,ind_o_0_b,o_0_b, // boundary conditions
				n_r_0_b,ind_r_0_b,
	  			dim); // dimension
		}

		/* for rigid bodies only 
		-> rigid body motion */
		if (n_r>0) rigid_bodies(n_r,dim,
						N1_r,N2_r,N3_r,
						N1_r_type,N2_r_type,N3_r_type,
						x_r,u1_r,u2_r,u3_r,
						m_r,I1_r,I2_r,I3_r,D1_r,D2_r,D3_r,
						rb_s,rb_b,
						ind_r_0_s,ind_r_0_b,
						a_r,alpha_r,
						f_r,M_r,
						type_a_0_r,frame_a_0_r,
						a_0_r,
						type_f_0_r,frame_f_0_r,
						f_0_r,M_0_r,
						n_s,n_b,
						x_s,x_b,
						f_s,f_b,
						n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u,
						t,dt,
						ax,ay,az);

      	/* full step (corrected values) for (0) Euler
	  	or corrector step for (1) central acceleration or (2) predictor-corrector
		-> V is actual value of variable (x, v, rho, u, e, S, drho, du, de, dS) to be predicted and corrected
		-> V_o is value of variable stored at beginning of time step interval
			-> for variables updated during time integration (x, v, rho, u, e, S, drho, du, de, dS)
			-> not necessary for p, c, h (directly calculated state variables)
        for corrector step for (3) predictor-corrector leapfrog with only v_s, rho_s and u_s updated
            -> update V:=V-V_o
		  	-> return V:=V+V_o as V_o is corrected by (dV-dV_o) */
		move(num_m,type_m, // materials
			rho_m,T_m,gamma_m,kappa_m,mu_m,
			coef1_m,coef2_m,coef3_m,coef4_m,
			n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u, // functions
			n_s,num_s,mat_s, // particles
			m_s,x_s_0,c_s_0, // initial values
			p_s,c_s,h_s, // state values
			dv_s,//O_s,
			//gradvx_s,gradvy_s,gradvz_s,
			x_s,v_s, // corrected values x_s and v_s
			a_s, // by updated a_s
			rho_s,u_s, // corrected values rho_s and u_s
			drhodt_s,dudt_s, // by updated drhodt_s and dudt_s
			e_s,S_s, // corrected values e_s and S_s
			dedt_s,dSdt_s, // by updated dedt_s and dSdt_s
			x_s_o,v_s_o, // beginning of time step interval
			a_s_o, // beginning of time step interval
			rho_s_o,u_s_o, // beginning of time step interval
			drhodt_s_o,dudt_s_o, // for updating values for (3) predictor-corrector leapfrog
			e_s_o,S_s_o, // beginning of time step interval
			dedt_s_o,dSdt_s_o, // for updating values for (3) predictor-corrector leapfrog
			n_b,num_b,x_b_0, // nodes
			x_b,v_b, // corrected values x_v and v_b
			a_b, // by updated a_b
			x_b_o,v_b_o,a_b_o, // beginning of time step interval
			n_o_0_s,o_0_s,ind_o_0_s,type_o_0_s, // boundary conditions on SPH
			n_o_0_b,o_0_b,ind_o_0_b,type_o_0_b, // boundary conditions on FEM
			n_r,num_r,type_r,COG_r,N1_r,N2_r,N3_r, // rigid bodies
			//COG_r_type,
			N1_r_type,N2_r_type,N3_r_type,
			x_r_0,psi_r_0, // initial values
			x_r,v_r,psi_r,o_r, // corrected values x_r, v_r, psi_r and o_r
			a_r,alpha_r, // by updated a_r and alpha_r,
			x_r_o,v_r_o,psi_r_o,o_r_o, // beginning of time step interval
			a_r_o,alpha_r_o, // for updating values for (3) predictor-corrector leapfrog
			//x1_r,x2_r,x3_r,
			//x1_r_o,x2_r_o,x3_r_o,
			u1_r,u2_r,u3_r, // principal moments of inertial vectors
			constrained_s_rb,constrained_b_rb,  // particles and nodes constrained in rigid bodies
			rb_s,rb_b, // particles and nodes in rigid bodies
			ind_r_0_s,ind_r_0_b, // boundary conditions on rigid bodies
			//n_r_0_s,n_r_0_b,
			l_s,l_b, // particles and nodes in rigid bodies local coordinate systems
			dim, // dimensions
			ix, // XSPH switch
			xeps, // XSPH constant
			constrained_s,constrained_b,constrained_r, // constraint attributes
			constrained_r_frame, // coordinate systems
			integration,1, // 1 = corrector (for integration = 1, 2, 3)
			t,cori*dt, // time and corrector time step
			v_max, // maximum velocity
			sigma, // parameters
			h0); // smoothing lenght coefficient

	  	/* for FEM only */
		if (n_e>0)
		{
	      	/* element areas and test on element areas */
	      	volume(type_m, // materials
		    	n_b,x_b, // corrected values
	     		n_e,num_e,mat_e,nod_e, // elements
		     	V_e, // element length or area or volume
		     	dim); // dimension
		}

      	/* update time step */
		t=t+dt;

		/* update cycle count */
      	count_cycle++;
    }

  	/* freeing variables */
  	freeall(n_m,num_m,type_m,domain_m, // materials
		rho_m,mu_m,T_m,gamma_m,kappa_m,
		n_u,num_u,fun_u, // functions
		fmx_u,fmy_u,fdx_u,fdy_u,fxi_u,fyi_u,
		n_c,num_c,seg_c,slave_c,master_c,sw_c,mat_c, // contacts
		m_slave,m_master,ct_c,klin_c,knon_c,kf_c,kd_c,f_c,
		n_s,ix,num_s,mat_s, // particles
		x_s_0,c_s_0,h_s_0,
		x_s,v_s,dv_s,
		//x_s_a,v_s_a,
		x_s_o,v_s_o,
		//gradvx_s,gradvy_s,gradvz_s,
		V_s,m_s,
		rho_s,u_s,p_s,c_s,h_s,mu_s,
		//rho_s_a,u_s_a,//p_s_a,c_s_a,h_s_a,
		rho_s_o,u_s_o,//p_s_o,c_s_o,h_s_o,
		drhodt_s,dudt_s,
		drhodt_s_o,dudt_s_o,
		a_s,
		a_s_o,
		f_s,
		e_s,S_s,
		//e_s_a,S_s_a,
		e_s_o,S_s_o,
		dedt_s,dSdt_s,
		dedt_s_o,dSdt_s_o,O_s,
		//O_s_o,
		n_b,num_b,mat_b, // nodes
		x_b_0,
		x_b,v_b,
		//x_b_a,v_b_a,
		x_b_o,v_b_o,
		a_b,
		a_b_o,
		f_b,
		nnopt,neig_s,nn_s, // nearest neighbours
		n_a_0_s,ind_a_0_s,type_a_0_s,frame_a_0_s, // particles boundary accelerations
		a_0_s,
		n_f_0_s,ind_f_0_s,type_f_0_s,frame_f_0_s, // particles boundary forces
		f_0_s,M_0_s,
		n_o_0_s,ind_o_0_s,type_o_0_s, // particles boundary conditions
		o_0_s,
		n_a_0_b,ind_a_0_b,type_a_0_b,frame_a_0_b, // nodal boundary accelerations
		a_0_b,
		n_f_0_b,ind_f_0_b,type_f_0_b,frame_f_0_b, // nodal boundary forces
		f_0_b,M_0_b,
		n_o_0_b,ind_o_0_b,type_o_0_b, // nodal boundary conditions
		o_0_b,
		n_d_0_b,ind_d_0_b,
		d_0_b,
		constrained_s,constrained_b, // constraint attributes for particles and nodes
		constrained_r,constrained_r_frame, // constraint attributes for rigid bodies
		//n_unc_s,n_unc_b,
		//unconstrained_s,unconstrained_b,
		n_r, // rigid bodies
		num_r,type_r,mat_r,COG_r,N1_r,N2_r,N3_r,
		COG_r_type,N1_r_type,N2_r_type,N3_r_type,
		m_r,I1_r,I2_r,I3_r,D1_r,D2_r,D3_r,
		x_r_0,psi_r_0,
		x_r,v_r,psi_r,o_r,
		//x_r_a,v_r_a,psi_r_a,o_r_a,
		x_r_o,v_r_o,psi_r_o,o_r_o,
		//x1_r,x2_r,x3_r,
		//x1_r_o,x2_r_o,x3_r_o,
		a_r,alpha_r,
		a_r_o,alpha_r_o,
		f_r,M_r,
		type_a_0_r,frame_a_0_r,a_0_r,
		type_f_0_r,frame_f_0_r,f_0_r,M_0_r,
		u1_r,u2_r,u3_r,
		constrained_s_rb,constrained_b_rb,
		rb_s,rb_b,
		ind_r_0_s,ind_r_0_b,
		l_s,l_b,
		n_n, // boundary equation of motion materices
		M,B,K,f,
		n_e, // elements
		nod_e,num_e,mat_e,    
		V_e);

  	//time_s=(clock()-t_start)/(double)clk_tck;
  	time_s=(clock()-t_start)/(double)CLOCKS_PER_SEC;
  	time_h=time_s/3600;
  	time_s=time_s-time_h*3600.0;
  	time_m=time_s/60;
  	time_s=time_s-time_m*60.0;

	/* line space to output */
	if (!error) fprintf(stdout,"\n");

	/* final output */
  	fprintf(stdout,
	  	" -> %d states",
	  	count_save);
  
  	if (error==1)
    {
      	fprintf(stdout,
	      	" + 1 error state saved"
	      	"\n"
	      	"\n");
      
      	fprintf(stdout,
	      	"\n"
	      	"error termination"
	      	"\n"
	      	"\n");
    }
	else if (error==2)
	{
      	fprintf(stdout,
	      	" + 1 signal state saved"
	      	"\n"
	      	"\n");
      
      	fprintf(stdout,
	      	"termination on signal"
	      	"\n"
	      	"\n");
	}
  	else
    {
      	fprintf(stdout,
	      	" saved"
	      	"\n"
	      	"\n");
      
      	fprintf(stdout,
	      	"normal termination"
	      	"\n"
	      	"\n");
    }
  
  	/* time output */
	fprintf(stdout,
	  	"elapsed time = %d hours %d minutes %4.2f seconds"
	  	"\n"
	  	"\n",
	  	time_h,time_m,time_s);
 
  	/* termination */
	return(0);
}
