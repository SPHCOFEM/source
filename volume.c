#include "header.h"

#undef __FUNC__
#define __FUNC__ "volume"

void volume(int *type_m,

    int n_b,
    double *x_b,

    int n_e,int *num_e,int *mat_e,int *nod_e,
    double *V_e,

    int dim)
{
	int i;
  
  	double ux,uy,uz,vx,vy,vz,nx,ny,nz;

  	for (i=0;i<n_e;i++)
    {
    	switch (dim)
		{
			case 1: // 1D
	  		{
	    		V_e[i]=fabs(x_b[nod_e[i+n_e]]-x_b[nod_e[i]]);

			    break;
	  		}
			case 2: // 2D
	  		{
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 0: /* rigid body */
					{
						if (nod_e[i+n_e]==nod_e[i+2*n_e]) /* bar or beam */
						{
							V_e[i]=norm(x_b[nod_e[i]]-x_b[nod_e[i+n_e]],
			      				x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b],
			      				x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b]);
						}
						else if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
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

							V_e[i]=norm(nx,ny,nz)/2.0;
						}
						else /* rectangle */
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

							V_e[i]=norm(nx,ny,nz)/2.0; /* triangle */

							if (nod_e[i+2*n_e]!=nod_e[i+3*n_e]) 
							{
								ux=x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]];
								uy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b];
								uz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b];

								vx=x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]];
								vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
								vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];

								nx=uy*vz-vy*uz;
								ny=-ux*vz+vx*uz;
								nz=ux*vy-vx*uy;
								
								V_e[i]=V_e[i]+norm(nx,ny,nz)/2.0; /* rectangle */
							}
						}
						break;
					}
					case 4: /* bar */
					/* If we do not use the break statement, 
					all statements after the matching label
					are also executed */
	      			case 5: /* beam */
					{
		  				V_e[i]=norm(x_b[nod_e[i]]-x_b[nod_e[i+n_e]],
			      			x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b],
			      			x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b]);

		  				break;
					}
	      			case 6: /* triangle or rectangle */
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

						V_e[i]=norm(nx,ny,nz)/2.0; /* triangle */

		  				if (nod_e[i+2*n_e]!=nod_e[i+3*n_e]) /* rectangle */
		    			{
							ux=x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]];
							uy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b];
							uz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b];

							vx=x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]];
							vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
							vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];

							nx=uy*vz-vy*uz;
							ny=-ux*vz+vx*uz;
							nz=ux*vy-vx*uy;
							
							V_e[i]=V_e[i]+norm(nx,ny,nz)/2.0;
		    			}
		  				break;
					}
	      		}
	    		break;
	  		}
			case 3: // 3D
	  		{
	    		switch (type_m[mat_e[i]])
	      		{
	      			case 0: /* rigid body */
					{
						if (nod_e[i+n_e]==nod_e[i+2*n_e]) /* bar or beam */
						{
							V_e[i]=norm(x_b[nod_e[i]]-x_b[nod_e[i+n_e]],
			      				x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b],
			      				x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b]);
						}
						else if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle */
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

							V_e[i]=norm(nx,ny,nz)/2.0;
						}
						else /* tetrahedron (no rectangle) */
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

							vx=x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]];
							vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
							vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];
							
							V_e[i]=fabs(nx*vx+ny*vy+nz*vz)/6.0;
						}
						break;
					}
	      			case 4: /* bar or shell */
					/* If we do not use the break statement, 
					all statements after the matching label
					are also executed */
	      			case 5: /* beam or membrane */
					{
						if (nod_e[i+n_e]==nod_e[i+2*n_e]) /* bar or beam */
						{
							V_e[i]=norm(x_b[nod_e[i]]-x_b[nod_e[i+n_e]],
			      				x_b[nod_e[i]+n_b]-x_b[nod_e[i+n_e]+n_b],
			      				x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+n_e]+2*n_b]);
						}
						else if (nod_e[i+2*n_e]==nod_e[i+3*n_e]) /* triangle or rectangle */
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

							V_e[i]=norm(nx,ny,nz)/2.0;  /* triangle */

							if (nod_e[i+2*n_e]!=nod_e[i+3*n_e]) /* rectangle */
							{
								ux=x_b[nod_e[i]]-x_b[nod_e[i+2*n_e]];
								uy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+2*n_e]+n_b];
								uz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+2*n_e]+2*n_b];

								vx=x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]];
								vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
								vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];

								nx=uy*vz-vy*uz;
								ny=-ux*vz+vx*uz;
								nz=ux*vy-vx*uy;
								
								V_e[i]=V_e[i]+norm(nx,ny,nz)/2.0;
							}
						}
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

						vx=x_b[nod_e[i]]-x_b[nod_e[i+3*n_e]];
						vy=x_b[nod_e[i]+n_b]-x_b[nod_e[i+3*n_e]+n_b];
						vz=x_b[nod_e[i]+2*n_b]-x_b[nod_e[i+3*n_e]+2*n_b];
						
						V_e[i]=fabs(nx*vx+ny*vy+nz*vz)/6.0;

						break;
					}
	      			break;
	      		}
	    		break;
	  		}
		}
      	
		/* test on zero or negative volume */
		if (V_e[i]<=0.0)
		{
	  		fprintf(stdout,
				"ERROR: boundary element number %d has zero or negative volume %f!"
			  	"\n",
		  		num_e[i],V_e[i]);
	  		exit(0);
		}
    }
}
