#include "header.h"
#include "interpolation.h"

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
	double ax,double ay,double az)
{
	int i,j,k;
	int count_rb_s=0,count_rb_b=0;
	double afx=0.0,afy=0.0,afz=0.0;
	double ux=1.0,uy=0.0,uz=0.0,un=1.0,vx=0.0,vy=1.0,vz=0.0,vn=1.0,wx=0.0,wy=0.0,wz=1.0,wn=1.0;
	double rx=0.0,ry=0.0,rz=0.0,ox=0.0,oy=0.0,oz=0.0,fx=0.0,fy=0.0,fz=0.0,Mx=0.0,My=0.0,Mz=0.0;
	double Di[6]; // inversed inertia matrix

	/* set rigid bodies forces to zero */
	memset(a_r,0,3*n_r*sizeof(double));
	memset(alpha_r,0,3*n_r*sizeof(double));
	memset(f_r,0,3*n_r*sizeof(double));
	memset(M_r,0,3*n_r*sizeof(double));

	for (i=0;i<n_r;i++)
	{
		/* local axis system (u, v, w)
		-> calculated before x_r moves
		-> x_s and x_b in rigid bodies cannot be constrained 
			-> only COG can be constrained but by its x_s or x_b
		u = N1 - COG;
		v = N2 - COG;
		w = N3 - COG; */

		/* first principal axis
		-> N1_r_type = 1 / 0 (particle / node) */
		if (dim>1)
		{
			if (N1_r_type[i]) // particle
			{
				u1_r[i]=x_s[N1_r[i]]-x_r[i];
				u1_r[i+n_r]=x_s[N1_r[i]+n_s]-x_r[i+n_r];
				u1_r[i+2*n_r]=x_s[N1_r[i]+2*n_s]-x_r[i+2*n_r];
			}
			else // node
			{
				u1_r[i]=x_b[N1_r[i]]-x_r[i];
				u1_r[i+n_r]=x_b[N1_r[i]+n_b]-x_r[i+n_r];
				u1_r[i+2*n_r]=x_b[N1_r[i]+2*n_b]-x_r[i+2*n_r];
			}
			ux=u1_r[i];uy=u1_r[i+n_r];uz=u1_r[i+2*n_r];
			un=norm(ux,uy,uz);
			ux=ux/un;uy=uy/un;uz=uz/un;
		}

		if (dim>2)
		{
			/* second principal axis
			-> N2_r_type = 1 / 0 (particle / node) */
			if (N2_r_type[i]) // particle
			{
				u2_r[i]=x_s[N2_r[i]]-x_r[i];
				u2_r[i+n_r]=x_s[N2_r[i]+n_s]-x_r[i+n_r];
				u2_r[i+2*n_r]=x_s[N2_r[i]+2*n_s]-x_r[i+2*n_r];
			}
			else // node
			{
				u2_r[i]=x_b[N2_r[i]]-x_r[i];
				u2_r[i+n_r]=x_b[N2_r[i]+n_b]-x_r[i+n_r];
				u2_r[i+2*n_r]=x_b[N2_r[i]+2*n_b]-x_r[i+2*n_r];
			}
			vx=u2_r[i];vy=u2_r[i+n_r];vz=u2_r[i+2*n_r];
			vn=norm(vx,vy,vz);
			vx=vx/vn;vy=vy/vn;vz=vz/vn;

			/* third principal axis
			-> N3_r_type = 1 / 0 (particle / node) */
			if (N3_r_type[i]) // particle
			{
				u3_r[i]=x_s[N3_r[i]]-x_r[i];
				u3_r[i+n_r]=x_s[N3_r[i]+n_s]-x_r[i+n_r];
				u3_r[i+2*n_r]=x_s[N3_r[i]+2*n_s]-x_r[i+2*n_r];
			}
			else // node
			{
				u3_r[i]=x_b[N3_r[i]]-x_r[i];
				u3_r[i+n_r]=x_b[N3_r[i]+n_b]-x_r[i+n_r];
				u3_r[i+2*n_r]=x_b[N3_r[i]+2*n_b]-x_r[i+2*n_r];
			}
			wx=u3_r[i];wy=u3_r[i+n_r];wz=u3_r[i+2*n_r];
			wn=norm(wx,wy,wz);
			wx=wx/wn;wy=wy/wn;wz=wz/wn;
		}

		/* acceleration field */
		if (type_a_0_r[i]==0) // constant
		{
			afx=a_0_r[i];
			afy=a_0_r[i+n_r];
			afz=a_0_r[i+2*n_r];
		}
		else // function
		{
			afx=interpolation(t-dt,(int)a_0_r[i+n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			afy=interpolation(t-dt,(int)a_0_r[i+n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			afz=interpolation(t-dt,(int)a_0_r[i+2*n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
		}

		/* acceleration field velocity frame */
		switch(frame_a_0_r[i])
		{
			case 0: // translational and rotatinal velocities in global coordinate system
			{
				break;
			}
			case 1: // translational velocity in local axis system, rotational velocity in global coordinate system
			{
				/* acceleration field 
				-> from local to global coordinate system */
				rx=afx;
				ry=afy;
				rz=afz;
				
				afx=rx*ux+ry*vx+rz*wx;
				afy=rx*uy+ry*vy+rz*wz;
				afz=rx*uz+ry*vz+rz*wz;

				break;
			}
			case 2: // translational velocity in global axis system, rotational velocity in local coordinate system
			{
				/* rotational acceleration field does not exist */
				break;
			}
			case 3: // translational and rotatiobal velocities in local coordinate system
			{
				/* translational velocity 
				-> from local to global coordinate system */

				/* rotational acceleration field does not exist */
				break;
			}
		}

		/* concentrated forces and moments */
		if (type_f_0_r[i]==0) // constant
		{
			f_r[i]=f_0_r[i];
			f_r[i+n_r]=f_0_r[i+n_r];
			f_r[i+2*n_r]=f_0_r[i+2*n_r];
			M_r[i]=M_0_r[i];
			M_r[i+n_r]=M_0_r[i+n_r];
			M_r[i+2*n_r]=M_0_r[i+2*n_r];
		}
		else // function
		{
			f_r[i]=interpolation(t-dt,(int)f_0_r[i+n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f_r[i+n_r]=interpolation(t-dt,(int)f_0_r[i+n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			f_r[i+2*n_r]=interpolation(t-dt,(int)f_0_r[i+2*n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			M_r[i]=interpolation(t-dt,(int)M_0_r[i],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			M_r[i+n_r]=interpolation(t-dt,(int)M_0_r[i+n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
			M_r[i+2*n_r]=interpolation(t-dt,(int)M_0_r[i+2*n_r],n_u,num_u,fun_u,fun_loc,fxi_u,fyi_u);
		}
				
		/* forces and moments frame */
		switch(frame_f_0_r[i])
		{
			case 0: // forces and moments in global coordinate system
			{
				/* forces in global coordinate system 
				   -> calculating displacements in global coordinate system */

				/* moments
				   -> from global to local coordinate system 
				      -> calculating rotations in local coordinate system */
				Mx=M_r[i];
				My=M_r[i+n_r];
				Mz=M_r[i+2*n_r];

				M_r[i]=Mx*ux+My*uy+Mz*uz;
				M_r[i+n_r]=Mx*vx+My*vy+Mz*vz;
				M_r[i+2*n_r]=Mx*wx+My*wy+Mz*wz;

				break;
			}
			case 1: // forces in local axis system, moments in global coordinate system
			{
				/* forces 
				   -> from local to global coordinate system 
				    -> calculating displacements in global coordinate system */
				fx=f_r[i];
				fy=f_r[i+n_r];
				fz=f_r[i+2*n_r];
				
				f_r[i]=fx*ux+fy*vx+fz*wx;
				f_r[i+n_r]=fx*uy+fy*vy+fz*wy;
				f_r[i+2*n_r]=fx*uz+fy*vz+fz*wz;

				/* moments 
				   -> from global to local coordinate system 
				      -> calculating rotations in local coordinate system */
				Mx=M_r[i];
				My=M_r[i+n_r];
				Mz=M_r[i+2*n_r];

				M_r[i]=Mx*ux+My*uy+Mz*uz;
				M_r[i+n_r]=Mx*vx+My*vy+Mz*vz;
				M_r[i+2*n_r]=Mx*wx+My*wy+Mz*wz;

				break;
			}
			case 2: // forces in global axis system, moments in local coordinate system
			{
				/* forces in global coordinate system 
				   -> calculating displacements in global coordinate system */

				/* moments in local coordinate system 
				   -> calculating rotations in local coordinate system */

				break;
			}
			case 3: // forces and moments in local coordinate system
			{
				/* forces 
				-> from local to global coordinate system 
				   -> calculating displacements in global coordinate system */
				fx=f_r[i];
				fy=f_r[i+n_r];
				fz=f_r[i+2*n_r];
				
				f_r[i]=fx*ux+fy*vx+fz*wx;
				f_r[i+n_r]=fx*uy+fy*vy+fz*wy;
				f_r[i+2*n_r]=fx*uz+fy*vz+fz*wz;

				/* moments in local coordinate system 
				   -> calculating rotations in local coordinate system */

				break;
			}
		}

		/* contact forces and moments */
		for (j=0;j<rb_s[i];j++) // particles
		{
			k=ind_r_0_s[j+count_rb_s]; // particle index

			/* move force to COG and add moments
			   f_0_s inputs to accelerations.c (only to COG)
			   type_f_0_s (function or double) -> f_0_s only on COG by input definition
			-> f_r = f_r + f_s
			-> M_r = M_r + f_s cross COG
			   -> project to (u1_r, u2_r, u3_r) */
			
			/* external (contact) force projected to local coordinate frame 
			fx=f_s[k]*ux+f_s[k+n_s]*uy+f_s[k+2*n_s]*uz;
			fy=f_s[k]*vx+f_s[k+n_s]*vy+f_s[k+2*n_s]*vz;
			fz=f_s[k]*wx+f_s[k+n_s]*wy+f_s[k+2*n_s]*wz;

			f_r[i]=f_r[i]+fx;
			f_r[i+n_r]=f_r[i+n_r]+fy;
			f_r[i+2*n_r]=f_r[i+2*n_r]+fz;
			*/

			/* external (contact) force in global coordinate frame */
			f_r[i]=f_r[i]+f_s[k];
			f_r[i+n_r]=f_r[i+n_r]+f_s[k+n_s];
			f_r[i+2*n_r]=f_r[i+2*n_r]+f_s[k+2*n_s];

			/* added moment projected to local coordinate frame 
			(f_s[k],f_s[k+n_s],f_s[k+2*n_s])
			(x_s[k],x_s[k+n_k],x_s[k+2*n_k]) */

			/* cross product between force and particle vector in local coordinate frame */
			Mx=f_s[k+n_s]*(x_s[k+2*n_s]-x_r[i+2*n_r])-f_s[k+2*n_s]*(x_s[k+n_s]-x_r[i+n_r]);
			My=f_s[k+2*n_s]*(x_s[i]-x_r[i])-f_s[k]*(x_s[k+2*n_s]-x_r[i+2*n_r]);
			Mz=f_s[k]*(x_s[k+n_s]-x_r[i+n_r])-f_s[k+n_s]*(x_s[k]-x_r[i]);

			/* scalar product to local coordinate frame */
			ox=Mx*ux+My*uy+Mz*uz;
			oy=Mx*vx+My*vy+Mz*vz;
			oz=Mx*wx+My*wy+Mz*wz;

			M_r[i]=M_r[i]+ox;
			M_r[i+n_r]=M_r[i+n_r]+oy;
			M_r[i+2*n_r]=M_r[i+2*n_r]+oz;
		}
		count_rb_s=count_rb_s+rb_s[i]; // update vector position
		
		for (j=0;j<rb_b[i];j++) // nodes
		{
			k=ind_r_0_b[j+count_rb_b]; // nodal index

			/* move force to COG and add moments
			   f_0_b inputs to damping_rhs.c (only to COG)
			   type_f_0_b (function or double) -> f_0_b only on COG by input definition
			-> f_r = f_r + f_b
			-> M_r = M_r + f_b cross COG
			   -> project to (u1_r, u2_r, u3_r) */

			/* external (contact) force projected to local coordinate frame 
			fx=(f_b[k]*ux+f_b[k+n_b]*uy+f_b[k+2*n_b]*uz)/un;
			fy=(f_b[k]*vx+f_b[k+n_b]*vy+f_b[k+2*n_b]*vz)/vn;
			fz=(f_b[k]*wx+f_b[k+n_b]*wy+f_b[k+2*n_b]*wz)/wn;

			f_r[i]=f_r[i]+fx;
			f_r[i+n_r]=f_r[i+n_r]+fy;
			f_r[i+2*n_r]=f_r[i+2*n_r]+fz;
			*/

			/* external (contact) force in global coordinate frame */
			f_r[i]=f_r[i]+f_b[k];
			f_r[i+n_r]=f_r[i+n_r]+f_b[k+n_b];
			f_r[i+2*n_r]=f_r[i+2*n_r]+f_b[k+2*n_b];

			/* added moment projected to local coordinate frame
			(f_b[k],f_b[k+n_b],f_b[k+2*n_b])
			(x_b[k],x_b[k+n_k],x_b[k+2*n_k]) */

			/* cross product between force and nodal vector in local coordinate frame */
			Mx=f_b[k+n_b]*(x_b[k+2*n_b]-x_r[i+2*n_r])-f_b[k+2*n_b]*(x_b[k+n_b]-x_r[i+n_r]);
			My=f_b[k+2*n_b]*(x_b[i]-x_r[i])-f_b[k]*(x_b[k+2*n_b]-x_r[i+2*n_r]);
			Mz=f_b[k]*(x_b[k+n_b]-x_r[i+n_r])-f_b[k+n_b]*(x_b[k]-x_r[i]);

			/* scalar product to local coordinate frame */
			ox=Mx*ux+My*uy+Mz*uz;
			oy=Mx*vx+My*vy+Mz*vz;
			oz=Mx*wx+My*wy+Mz*wz;

			M_r[i]=M_r[i]+ox;
			M_r[i+n_r]=M_r[i+n_r]+oy;
			M_r[i+2*n_r]=M_r[i+2*n_r]+oz;
		}
		count_rb_b=count_rb_b+rb_b[i]; // update vector position

		/* rigid body dynamics 
		D = [  I1, -D1, -D2]
			[ -D1,  I2, -D3]
			[ -D2, -D3,  I3]
		
		pow(D,-1) = [(D3^2 - I2*I3)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3),
						-(D2*D3 + D1*I3)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3),
							-(D1*D3 + D2*I2)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3)]
					[-(D2*D3 + D1*I3)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3),
						(D2^2 - I1*I3)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3),
							-(D1*D2 + D3*I1)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3)]
					[-(D1*D3 + D2*I2)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3),
						-(D1*D2 + D3*I1)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3),
							(D1^2 - I1*I2)/(I3*D1^2 + 2*D1*D2*D3 + I2*D2^2 + I1*D3^2 - I1*I2*I3)] */

		/* translational acceleration */
		a_r[i]=f_r[i]/m_r[i]+ax+afx;
		a_r[i+n_r]=f_r[i+n_r]/m_r[i]+ay+afy;
		a_r[i+2*n_r]=f_r[i+2*n_r]/m_r[i]+az+afz;

		/* inversed inertia matrix */
		Di[0]=(sqr(D3_r[i])-I2_r[i]*I3_r[i])/(I3_r[i]*sqr(D1_r[i])+2*D1_r[i]*D2_r[i]*D3_r[i]+I2_r[i]*sqr(D2_r[i])+I1_r[i]*sqr(D3_r[i])-I1_r[i]*I2_r[i]*I3_r[i]);
		Di[1]=-(D2_r[i]*D3_r[i]+D1_r[i]*I3_r[i])/(I3_r[i]*sqr(D1_r[i])+2*D1_r[i]*D2_r[i]*D3_r[i]+I2_r[i]*sqr(D2_r[i])+I1_r[i]*sqr(D3_r[i])-I1_r[i]*I2_r[i]*I3_r[i]);
		Di[2]=-(D1_r[i]*D3_r[i]+D2_r[i]*I2_r[i])/(I3_r[i]*sqr(D1_r[i])+2*D1_r[i]*D2_r[i]*D3_r[i]+I2_r[i]*sqr(D2_r[i])+I1_r[i]*sqr(D3_r[i])-I1_r[i]*I2_r[i]*I3_r[i]);
		Di[3]=(sqr(D2_r[i])-I1_r[i]*I3_r[i])/(I3_r[i]*sqr(D1_r[i])+2*D1_r[i]*D2_r[i]*D3_r[i]+I2_r[i]*sqr(D2_r[i])+I1_r[i]*sqr(D3_r[i])-I1_r[i]*I2_r[i]*I3_r[i]);
		Di[4]=-(D1_r[i]*D2_r[i]+D3_r[i]*I1_r[i])/(I3_r[i]*sqr(D1_r[i])+2*D1_r[i]*D2_r[i]*D3_r[i]+I2_r[i]*sqr(D2_r[i])+I1_r[i]*sqr(D3_r[i])-I1_r[i]*I2_r[i]*I3_r[i]);
		Di[5]=(sqr(D1_r[i])-I1_r[i]*I2_r[i])/(I3_r[i]*sqr(D1_r[i])+2*D1_r[i]*D2_r[i]*D3_r[i]+I2_r[i]*sqr(D2_r[i])+I1_r[i]*sqr(D3_r[i])-I1_r[i]*I2_r[i]*I3_r[i]);

		/* rotational acceleration */
		alpha_r[i]=Di[0]*M_r[i]+Di[1]*M_r[i+n_r]+Di[2]*M_r[i+2*n_r];
		alpha_r[i+n_r]=Di[1]*M_r[i]+Di[3]*M_r[i+n_r]+Di[4]*M_r[i+2*n_r];
		alpha_r[i+2*n_r]=Di[2]*M_r[i]+Di[4]*M_r[i+n_r]+Di[5]*M_r[i+2*n_r];
	}
}
