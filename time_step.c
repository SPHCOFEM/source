#include "header.h"

#undef __FUNC__
#define __FUNC__ "time_step"

double time_step(int *type_m,
	double *rho_m,double *T_m,double *mu_m,double *kappa_m,

	int n_s,int *mat_s,
	double *x_s,double *v_s,
	//double *x_s_a,double *v_s_a,
	double *x_s_o,double *v_s_o,

	double *rho_s,double *u_s,double *p_s,double *c_s,double *h_s,
	//double *rho_s_a,double *u_s_a,//double *p_s_a,double *c_s_a,double *h_s_a,
	double *rho_s_o,double *u_s_o,//double *p_s_o,double *c_s_o,double *h_s_o,

	double *a_s,double *drhodt_s,double *dudt_s,
	double *a_s_o,double *drhodt_s_o,double *dudt_s_o,
	
	double *e_s,double *S_s,double *dedt_s,double *dSdt_s,
	//double *e_s_a,double *S_s_a,
	double *e_s_o,double *S_s_o,double *dedt_s_o,double *dSdt_s_o,

	int n_b,
	double *x_b,double *v_b,
	//double *x_b_a,double *v_b_a,
	double *x_b_o,double *v_b_o,

	double *a_b,
	double *a_b_o,

	int n_e,
	int *mat_e,int *nod_e,
	double *V_e,

	int n_r,
	double *x_r,double *psi_r,double *v_r,double *o_r,
	//double *x_r_a,double *psi_r_a,double *v_r_a,double *o_r_a,
	double *x_r_o,double *psi_r_o,double *v_r_o,double *o_r_o,
	double *a_r,double *alpha_r,
	double *a_r_o,double *alpha_r_o,

	int dim,

	double alpha,double beta,

	double cour,double kstab,

	int integration,

	double *v_max,

	double dr,

	double dt_max,

	double dt)
{
	int i;
  
  	double v,dt_min=1e-9;
  	double aux1,aux2,aux3,aux4;
  	double nx,ny,nz,ux,uy,uz,vx,vy,vz;

	/* reset v_max to very small number */
	if (integration!=0)
	{
		*v_max=1e-9;
	}

	/* time step depending on previous one always decreases
	-> time step not depending on previous step can increase again
	-> but higher inital time step causes instability (optimum dt=0.1) */
	if (dt_max>0)
	{
		dt=dt_max;
	}
  
  	for (i=0;i<n_s;i++)
    {
		if (type_m[mat_s[i]]>0)
		{
			/* minimum time step */
			dt_min=cour*h_s[i]/(c_s[i]+0.6*(alpha*c_s[i]+beta*mu_m[mat_s[i]]));
			if (dt_min<dt) dt=dt_min;

			/* maximum velocity for predictor movement */
			if (integration!=0) // not necessary for Euler integration scheme
			{
				v=norm(v_s[i],v_s[i+n_s],v_s[i+2*n_s]);
				if (v>*v_max) *v_max=v;
			}
		}

		/* beginning of time step */
		x_s_o[i]=x_s[i];
		x_s_o[i+n_s]=x_s[i+n_s];
		x_s_o[i+2*n_s]=x_s[i+2*n_s];

		v_s_o[i]=v_s[i];
		v_s_o[i+n_s]=v_s[i+n_s];
		v_s_o[i+2*n_s]=v_s[i+2*n_s];

		a_s_o[i]=a_s[i];
		a_s_o[i+n_s]=a_s[i+n_s];
		a_s_o[i+2*n_s]=a_s[i+2*n_s];

		drhodt_s_o[i]=drhodt_s[i];
		dudt_s_o[i]=dudt_s[i];

		rho_s_o[i]=rho_s[i];
		u_s_o[i]=u_s[i];
		
		/* directly calculated
		p_s_o[i]=p_s[i];
		c_s_o[i]=c_s[i];
		h_s_o[i]=h_s[i];

		O_s_o[i]=O_s[i];
		O_s_o[i+n_s]=O_s[i+n_s];
		O_s_o[i+2*n_s]=O_s[i+2*n_s];
		O_s_o[i+3*n_s]=O_s[i+3*n_s];
		O_s_o[i+4*n_s]=O_s[i+4*n_s];
		O_s_o[i+5*n_s]=O_s[i+5*n_s];
		*/

		e_s_o[i]=e_s[i];
		e_s_o[i+n_s]=e_s[i+n_s];
		e_s_o[i+2*n_s]=e_s[i+2*n_s];
		e_s_o[i+3*n_s]=e_s[i+3*n_s];
		e_s_o[i+4*n_s]=e_s[i+4*n_s];
		e_s_o[i+5*n_s]=e_s[i+5*n_s];

		S_s_o[i]=S_s[i];
		S_s_o[i+n_s]=S_s[i+n_s];
		S_s_o[i+2*n_s]=S_s[i+2*n_s];
		S_s_o[i+3*n_s]=S_s[i+3*n_s];
		S_s_o[i+4*n_s]=S_s[i+4*n_s];
		S_s_o[i+5*n_s]=S_s[i+5*n_s];

		dedt_s_o[i]=dedt_s[i];
		dedt_s_o[i+n_s]=dedt_s[i+n_s];
		dedt_s_o[i+2*n_s]=dedt_s[i+2*n_s];
		dedt_s_o[i+3*n_s]=dedt_s[i+3*n_s];
		dedt_s_o[i+4*n_s]=dedt_s[i+4*n_s];
		dedt_s_o[i+5*n_s]=dedt_s[i+5*n_s];

		dSdt_s_o[i]=dSdt_s[i];
		dSdt_s_o[i+n_s]=dSdt_s[i+n_s];
		dSdt_s_o[i+2*n_s]=dSdt_s[i+2*n_s];
		dSdt_s_o[i+3*n_s]=dSdt_s[i+3*n_s];
		dSdt_s_o[i+4*n_s]=dSdt_s[i+4*n_s];
		dSdt_s_o[i+5*n_s]=dSdt_s[i+5*n_s];

		/* accelerations states 
		x_s_a[i]=x_s[i];
		x_s_a[i+n_s]=x_s[i+n_s];
		x_s_a[i+2*n_s]=x_s[i+2*n_s];

		v_s_a[i]=v_s[i];
		v_s_a[i+n_s]=v_s[i+n_s];
		v_s_a[i+2*n_s]=v_s[i+2*n_s];

		rho_s_a[i]=rho_s[i];
		u_s_a[i]=u_s[i];

		e_s_a[i]=e_s[i];
		e_s_a[i+n_s]=e_s[i+n_s];
		e_s_a[i+2*n_s]=e_s[i+2*n_s];
		e_s_a[i+3*n_s]=e_s[i+3*n_s];
		e_s_a[i+4*n_s]=e_s[i+4*n_s];
		e_s_a[i+5*n_s]=e_s[i+5*n_s];

		S_s_a[i]=S_s[i];
		S_s_a[i+n_s]=S_s[i+n_s];
		S_s_a[i+2*n_s]=S_s[i+2*n_s];
		S_s_a[i+3*n_s]=S_s[i+3*n_s];
		S_s_a[i+4*n_s]=S_s[i+4*n_s];
		S_s_a[i+5*n_s]=S_s[i+5*n_s];
		*/
	}

    for (i=0;i<n_e;i++)
	{
		if (type_m[mat_e[i]]>0)
		{
			switch (dim)
			{
				case 1:
				{
					break;
				}
				case 2:
				{
					switch (type_m[mat_e[i]])
					{
						case 4: /* bar */
						{
							dt_min=kstab*V_e[i]/sqrt(T_m[mat_e[i]]/rho_m[mat_e[i]])/(sqrt(1.0+sqr(dr))-dr);
							break;
						}
						case 5: /* beam */
						{
							aux1=kappa_m[mat_e[i]];
							aux2=aux1*sqr(kappa_m[mat_e[i]])/12.0;
							aux3=2.0*sqrt((aux2/aux1)*(12.0*aux2+4.0*aux1*sqr(V_e[i])))/(12.0*aux2+aux1*sqr(V_e[i]));
							aux4=V_e[i]/(aux3*sqrt(3.0));
							if (aux4<1) aux4=1.0;
							dt_min=kstab*V_e[i]/sqrt(T_m[mat_e[i]]/rho_m[mat_e[i]])/(sqrt(1.0+sqr(dr))-dr)*aux4;
							break;
						}
						case 6: /* triangle or rectangle */
						{
							aux1=norm(x_b[nod_e[i]]-x_b[nod_e[i+n_e]],
								x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b],
								x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b]);
							aux2=norm(x_b[nod_e[i+n_e]]-x_b[nod_e[i+2*n_e]],
								x_b[nod_e[i+n_e]+n_b]-x_b[nod_e[i+2*n_e]+n_b],
								x_b[nod_e[i+n_e]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b]);
							if (aux2>aux1) aux1=aux2;
							aux2=norm(x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]],
							x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b],
							x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b]);
							if (aux2>aux1) aux1=aux2;			
							if (nod_e[i+2*n_e]!=nod_e[i+3*n_e]) /* recatngle */
							{
								aux2=norm(x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]],
									x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b],
									x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b]);
									if (aux2>aux1) aux1=aux2;
							}
							dt_min=kstab*(V_e[i]/aux1)/(sqrt(T_m[mat_e[i]]/rho_m[mat_e[i]]));
							break;
						}
					}
					break;
				}
				case 3:
				{
					switch (type_m[mat_e[i]])
					{
						case 4: /* membrane */
						{
							aux1=norm(x_b[nod_e[i]]-x_b[nod_e[i+n_e]],
								x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b],
								x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b]);
							aux2=norm(x_b[nod_e[i+n_e]]-x_b[nod_e[i+2*n_e]],
								x_b[nod_e[i+n_e]+n_b]-x_b[nod_e[i+2*n_e]+n_b],
								x_b[nod_e[i+n_e]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b]);
							if (aux2>aux1) aux1=aux2;
							aux2=norm(x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]],
								x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b],
								x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b]);
							if (aux2>aux1) aux1=aux2;
							dt_min=kstab*(V_e[i]/aux1)/(sqrt(T_m[mat_e[i]]/rho_m[mat_e[i]]));
							break;
						}
						case 5: /* shell */
						{
							aux1=norm(x_b[nod_e[i]]-x_b[nod_e[i+n_e]],
								x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b],
								x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b]);
							aux2=norm(x_b[nod_e[i+n_e]]-x_b[nod_e[i+2*n_e]],
								x_b[nod_e[i+n_e]+n_b]-x_b[nod_e[i+2*n_e]+n_b],
								x_b[nod_e[i+n_e]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b]);
							if (aux2>aux1) aux1=aux2;
							aux2=norm(x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]],
								x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b],
								x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b]);
							if (aux2>aux1) aux1=aux2;
							aux3=V_e[i]/aux1/(kappa_m[mat_e[i]]*sqrt(3.0));
							if (aux3>1.0) aux3=1.0;  
							dt_min=kstab*(V_e[i]/aux1)/(sqrt(T_m[mat_e[i]]/rho_m[mat_e[i]]))*aux3;
							break;
						}
						case 6: /* tetrahedron */
						{
							ux=x_b[nod_e[i]]-x_b[nod_e[i+n_e]];
							uy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b];
							uz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b];

							vx=x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]];
							vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b];
							vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b];

							nx=uy*vz-vy*uz;
							ny=-ux*vz+vx*uz;
							nz=ux*vy-vx*uy;

							aux1=norm(nx,ny,nz)/2.0;

							ux=x_b[nod_e[i]]-x_b[nod_e[i+n_e]];
							uy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b];
							uz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b];

							vx=x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]];
							vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
							vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];

							nx=uy*vz-vy*uz;
							ny=-ux*vz+vx*uz;
							nz=ux*vy-vx*uy;

							aux2=norm(nx,ny,nz)/2.0;

							if (aux1<aux2) aux1=aux2;

							ux=x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]];
							uy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b];
							uz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b];

							vx=x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]];
							vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
							vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];

							nx=uy*vz-vy*uz;
							ny=-ux*vz+vx*uz;
							nz=ux*vy-vx*uy;

							aux2=norm(nx,ny,nz)/2.0;

							if (aux1<aux2) aux1=aux2;
							
							ux=x_b[nod_e[i+n_e]]-x_b[nod_e[i+2*n_e]];
							uy=x_b[nod_e[i+n_e]+n_b]-x_b[nod_e[i+2*n_e]+n_b];
							uz=x_b[nod_e[i+n_e]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b];

							vx=x_b[nod_e[i+n_e]]-x_b[nod_e[i+3*n_e]];
							vy=x_b[nod_e[i+n_e]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
							vz=x_b[nod_e[i+n_e]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];

							nx=uy*vz-vy*uz;
							ny=-ux*vz+vx*uz;
							nz=ux*vy-vx*uy;

							aux2=norm(nx,ny,nz)/2.0;

							if (aux1<aux2) aux1=aux2;

							aux2=T_m[mat_e[i]]*mu_m[mat_e[i]]/((1-2*mu_m[mat_e[i]])*(1+mu_m[mat_e[i]]));
							aux3=T_m[mat_e[i]]/(2*(1+mu_m[mat_e[i]]));
							aux4=0.0;		      
							dt_min=kstab*(V_e[i]/aux1)/(aux4+sqrt(sqr(sqrt((aux2+aux3)/rho_m[mat_e[i]]))+sqr(aux4)));

							break;
						}
					}
					break;
				}
			}
			if (dt_min<dt) dt=dt_min;
		}
	}

	/* minimum time step 
	if (dt<dt_min) dt=dt_min; // dt_init
	*/

	/* forced time step 
	if ((dt_max<0)&&(dt<-dt_max))
	{
		dt=-dt_max;
	}
	*/
      
    for (i=0;i<n_b;i++)
	{
      	/* beginning of time step */
		x_b_o[i]=x_b[i];
        x_b_o[i+n_b]=x_b[i+n_b];
        x_b_o[i+2*n_b]=x_b[i+2*n_b];

        v_b_o[i]=v_b[i];
        v_b_o[i+n_b]=v_b[i+n_b];
        v_b_o[i+2*n_b]=v_b[i+2*n_b];

      	a_b_o[i]=a_b[i];
      	a_b_o[i+n_b]=a_b[i+n_b];
      	a_b_o[i+2*n_b]=a_b[i+2*n_b];

		/* accelerations states 
		x_b_a[i]=x_b[i];
		x_b_a[i+n_b]=x_b[i+n_b];
		x_b_a[i+2*n_b]=x_b[i+2*n_b];

		v_b_a[i]=v_b[i];
		v_b_a[i+n_b]=v_b[i+n_b];
		v_b_a[i+2*n_b]=v_b[i+2*n_b];
		*/
	}


    for (i=0;i<n_r;i++)
	{
      	/* beginning of time step */
		x_r_o[i]=x_r[i];
        x_r_o[i+n_r]=x_r[i+n_r];
        x_r_o[i+2*n_r]=x_r[i+2*n_r];

		psi_r_o[i]=psi_r[i];
        psi_r_o[i+n_r]=psi_r[i+n_r];
        psi_r_o[i+2*n_r]=psi_r[i+2*n_r];

		v_r_o[i]=v_r[i];
        v_r_o[i+n_r]=v_r[i+n_r];
        v_r_o[i+2*n_r]=v_r[i+2*n_r];

		o_r_o[i]=o_r[i];
        o_r_o[i+n_r]=o_r[i+n_r];
        o_r_o[i+2*n_r]=o_r[i+2*n_r];

		a_r_o[i]=a_r[i];
        a_r_o[i+n_r]=a_r[i+n_r];
        a_r_o[i+2*n_r]=a_r[i+2*n_r];

		alpha_r_o[i]=alpha_r[i];
        alpha_r_o[i+n_r]=alpha_r[i+n_r];
        alpha_r_o[i+2*n_r]=alpha_r[i+2*n_r];

		/* not implemented 
		x1_r_o[i]=x1_r[i];
        x1_r_o[i+n_r]=x1_r[i+n_r];
        x1_r_o[i+2*n_r]=x1_r[i+2*n_r];

		x2_r_o[i]=x2_r[i];
        x2_r_o[i+n_r]=x2_r[i+n_r];
        x2_r_o[i+2*n_r]=x2_r[i+2*n_r];

		x3_r_o[i]=x3_r[i];
        x3_r_o[i+n_r]=x3_r[i+n_r];
        x3_r_o[i+2*n_r]=x3_r[i+2*n_r];
		*/

		/* accelerations states 
		x_r_a[i]=x_r[i];
        x_r_a[i+n_r]=x_r[i+n_r];
        x_r_a[i+2*n_r]=x_r[i+2*n_r];

		psi_r_a[i]=psi_r[i];
        psi_r_a[i+n_r]=psi_r[i+n_r];
        psi_r_a[i+2*n_r]=psi_r[i+2*n_r];

		v_r_a[i]=v_r[i];
        v_r_a[i+n_r]=v_r[i+n_r];
        v_r_a[i+2*n_r]=v_r[i+2*n_r];

		o_r_a[i]=o_r[i];
        o_r_a[i+n_r]=o_r[i+n_r];
        o_r_a[i+2*n_r]=o_r[i+2*n_r];
		*/
	}
    return(dt);
}
