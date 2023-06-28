#include "header.h"

#undef __FUNC__
#define __FUNC__ "contact_forces"

void contact_forces(int n_c,double *f_c,int f_t, // (0 = SPH versus FEM, 1 = FEM versus FEM)
	
	int count_slave,int count_master,int sw,int *seg_c,int *slave_c,int *master_c,
	double ct,double klin,double knon,double kf,double kd,double *m_slave,double *m_master,

	int n_s,
	double *x_s,double *v_s,double *f_s,
	double *h_s,

	int n_b,
	double *x_b,double *v_b,double *f_b,

	int n_e,int *nod_e,

	//double t, // t because of printing force
	double dt,

	int dim,

	int k)
{
	int i,j,ctriangle,ntriangle=1;
	int N1,N2,N3;
  	double s,p,fn,ft;
  	double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  	double n1,n2,n3,n11,n12,n13,t1,t2,t3;
  	double u1,u2,u3,v1,v2,v3;
  	double d1,d2,d3,k1,k2,k3;
  	double rt,rn,vr,at_s,an_s;
  	double vx1,vx2,vx3,vy1,vy2,vy3,vz1,vz2,vz3,vt1,vt2,vt3;
	
  	for (i=0;i<seg_c[k];i++)
    {
		/* contact thickness */
		if (ct<0)
		{
			if (!f_t) // negative contact thickness and SPH to FEM -> smoothing length
			{
				ct=-ct*h_s[slave_c[i+count_slave]];
			}
			else // negative contact thickness and FEM to FEM -> ct
			{
				// not yet implemented
			}
		}

		/* particles in same domain */
		//if (domain_m[mat_s[i]]==domain_m[mat_s[j]])
		//{
		for (j=0;j<seg_c[k+n_c];j++)
		{
			/* reset forces for case of no penetration */
			fn=0.0;
			ft=0.0;

			switch (dim)
			{
				case 1: /* 1D */
				{
					x1=x_b[nod_e[master_c[j+count_master]]];
					x2=x_b[nod_e[master_c[j+count_master]+n_e]];

					/* normal vector */
					n1=(x2-x1)/fabs(x2-x1);

					/* penetration */
					p=n1*(x_s[slave_c[i+count_slave]]-x1)+ct;

					if ((p>0.0)&&(p<=2.0*ct))
					{
						switch (sw)
						{
							case 1: /* soft contact */
							{
								if (klin>EPS)
								{
									fn=klin*p/ct+knon*p*p*p/(ct*ct*ct);
								}
								else
								{
									fn=m_slave[i+count_slave]*m_master[j+count_master]
										/(m_slave[i+count_slave]+m_master[j+count_master])
										/sqr(dt)*p;
				
									if (knon>EPS)
									{
										fn=fn+m_slave[i+count_slave]*m_master[j+count_master]/sqr(dt)*pow(p,knon);
									}
								}
								f_s[slave_c[i+count_slave]]=f_s[slave_c[i+count_slave]]-n1*fn;

								break;
							}
							case 2: /* sliding without separation - irrelevant in 1D */
							case 3: /* tied contact */
								{
									an_s=(v_b[nod_e[master_c[j+count_master]]]-v_s[slave_c[i+count_slave]])/dt;

									if (klin>EPS)
									{
										fn=-klin*an_s-knon*an_s*an_s*an_s;
									}
									else
									{
										fn=-m_slave[i+count_slave]*(an_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave]);

										if (knon>EPS)
										{
											fn=fn-m_slave[i+count_slave]*pow((an_s-f_s[slave_c[i+count_slave]]
												/m_slave[i+count_slave]),knon);
										}
									}
									f_s[slave_c[i+count_slave]]=f_s[slave_c[i+count_slave]]-fn;

									break;
								}
						}
			
						/* force on boundary element */
						f_b[nod_e[master_c[j+count_master]]]=f_b[nod_e[master_c[j+count_master]]]+n1*fn;
			
						/* contribution to the contact force */
						f_c[k]=f_c[k]+fn*n1;
					}
					break;
				}
				case 2: /* 2D */
				{
					/* nodes N1 and N2 */
					N1=nod_e[master_c[j+count_master]];
					N2=nod_e[master_c[j+count_master]+n_e];

					x1=x_b[N1];
					y1=x_b[N1+n_b];

					x2=x_b[N2];
					y2=x_b[N2+n_b];

					t1=x2-x1;
					t2=y2-y1;

					/* node located in element space */
					if ((t1*x_s[slave_c[i+count_slave]]+t2*x_s[slave_c[i+count_slave]+n_s]-t1*x1-t2*y1>=0.0)&&
						(t1*x_s[slave_c[i+count_slave]]+t2*x_s[slave_c[i+count_slave]+n_s]-t1*x2-t2*y2<=0.0))
					{
						n1=y1-y2;
						n2=x2-x1;

						/* normal vector */
						rn=norm(n1,n2,0.0);
						if (rn<EPS) rn=1.0; // zero normal vector

						/* penetration */
						p=(n1*x_s[slave_c[i+count_slave]]+n2*x_s[slave_c[i+count_slave]+n_s]-n1*x1-n2*y1)/rn+ct;

						if ((p>0.0)&&(p<=2.0*ct))
						{
							/* element velocity */
							vx1=v_b[N1];
							vy1=v_b[N1+n_b];

							vx2=v_b[N2];
							vy2=v_b[N2+n_b];

							vt1=(vx1+vx2)/2.0;
							vt2=(vy1+vy2)/2.0;

							/* relative motion */
							vr=v_s[slave_c[i+count_slave]]*vt1+v_s[slave_c[i+count_slave]+n_s]*vt2;
							if (vr > 0)
							{
								vr=1.0;
							}
							else
							{
								/* vr=1.0; rev 8. 1. 2015 - corrected direction */
								vr=-1.0;
							}
			
							/* projection of domain node */
							s=(n1*x1+n2*y1-n1*x_s[slave_c[i+count_slave]]-n2*x_s[slave_c[i+count_slave]+n_s])/(sqr(n1)+sqr(n2));

							x3=x_s[slave_c[i+count_slave]]+n1*s;
							y3=x_s[slave_c[i+count_slave]+n_s]+n2*s;

							/* unit tangent vector - velocity projection */
							rt=norm(t1,t2,0.0);
							if (rt<EPS) // zero tangent vector
							{
								vr=0.0; // multiplier -> ft = 0
								rt=1.0; // because of division by rt
							}

							k1=norm(x3-x2,y3-y2,0.0)/rt;
							k2=norm(x3-x1,y3-y1,0.0)/rt;

							switch (sw)
							{
								case 1: /* soft contact */
								{
									if (klin>EPS)
									{
										fn=klin*p/ct+knon*p*p*p/(ct*ct*ct);
									}
									else
									{
										/*fn=m_slave[i+count_slave]/sqr(dt)*(p-an_s*dt/2.0);*/
										fn=m_slave[i+count_slave]*m_master[j+count_master]
											/(m_slave[i+count_slave]+m_master[j+count_master])
											/sqr(dt)*p;
				
										if (knon>EPS)
										{
											fn=fn+m_slave[i+count_slave]*m_master[j+count_master]/sqr(dt)*pow(p,knon);
										}
									}
					
									/* friction force */
									//if (rt>0) // for rt=0.0 -> rt=1.0
									ft=vr*kf*fn;

									f_s[slave_c[i+count_slave]]=f_s[slave_c[i+count_slave]]-(fn*n1/rn+ft*t1/rt);
									f_s[slave_c[i+count_slave]+n_s]=f_s[slave_c[i+count_slave]+n_s]-(fn*n2/rn+ft*t2/rt);

									break;
								}
								case 2: /* sliding without separation */
								{
									v1=k1*v_b[N1]+k2*v_b[N2];
									v2=k1*v_b[N1+n_b]+k2*v_b[N2+n_b];

									an_s=(n1*(v1-v_s[slave_c[i+count_slave]])+
										n2*(v2-v_s[slave_c[i+count_slave]+n_s]))/rn/dt;
				
									at_s=(t1*f_s[slave_c[i+count_slave]]+
										t2*f_s[slave_c[i+count_slave]+n_s])/rt/m_slave[i+count_slave];

									if (klin>EPS)
									{
										fn=-klin*an_s-knon*an_s*an_s*an_s;
									}
									else
									{
										fn=-m_slave[i+count_slave]*(an_s-(n1*f_s[slave_c[i+count_slave]]+
											n2*f_s[slave_c[i+count_slave]+n_s])/m_slave[i+count_slave]/rn);

										if (knon>EPS)
										{
											fn=fn-m_slave[i+count_slave]*pow((an_s-(n1*f_s[slave_c[i+count_slave]]+
												n2*f_s[slave_c[i+count_slave]+n_s])/m_slave[i+count_slave]/rn),knon);
										}
									}

									/* friction force */
									ft=kf*fn;
			
									f_s[slave_c[i+count_slave]]=m_slave[i+count_slave]*(an_s*n1/rn+(at_s-ft/m_slave[i+count_slave])*t1/rt);
									f_s[slave_c[i+count_slave]+n_s]=m_slave[i+count_slave]*(an_s*n2/rn+(at_s-ft/m_slave[i+count_slave])*t2/rt);

									break;
								}
								case 3: /* tied contact */
								{
									v1=k1*v_b[N1]+k2*v_b[N2];
									v2=k1*v_b[N1+n_b]+k2*v_b[N2+n_b];

									at_s=(v1-v_s[slave_c[i+count_slave]])/dt;
									an_s=(v2-v_s[slave_c[i+count_slave]+n_s])/dt;

									if (klin>EPS)
									{
										fn=-klin*an_s-knon*an_s*an_s*an_s;
										ft=-klin*at_s-knon*at_s*at_s*at_s;
									}
									else
									{
										fn=-m_slave[i+count_slave]*(n1*(at_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave])+
											n2*(an_s-f_s[slave_c[i+count_slave]+n_s]/m_slave[i+count_slave]))/rn;

										ft=-m_slave[i+count_slave]*(t1*(at_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave])+
											t2*(an_s-f_s[slave_c[i+count_slave]+n_s]/m_slave[i+count_slave]))/rt;

										if (knon>EPS)
										{
											fn=fn-pow((n1*(at_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave])+
												n2*(an_s-f_s[slave_c[i+count_slave]+n_s]/m_slave[i+count_slave]))/rn,knon);

											ft=ft-pow((t1*(at_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave])+
												t2*(an_s-f_s[slave_c[i+count_slave]+n_s]/m_slave[i+count_slave]))/rt,knon);
										}
									}
									f_s[slave_c[i+count_slave]]=m_slave[i+count_slave]*at_s;
									f_s[slave_c[i+count_slave]+n_s]=m_slave[i+count_slave]*an_s;

									break;
								}
							}

							/* forces on boundary element */
							f_b[N1]=f_b[N1]+k1*(fn*n1/rn+ft*t1/rt);
							f_b[N1+n_b]=f_b[N1+n_b]+k1*(fn*n2/rn+ft*t2/rt);
			
							f_b[N2]=f_b[N2]+k2*(fn*n1/rn+ft*t1/rt);
							f_b[N2+n_b]=f_b[N2+n_b]+k2*(fn*n2/rn+ft*t2/rt);
			
							/* contribution to the contact force */
							f_c[k]=f_c[k]+fn*n1/rn+ft*t1/rt;
							f_c[k+n_c]=f_c[k+n_c]+fn*n2/rn+ft*t2/rt;
						}
					}
					break;
				}
				case 3: /* 3D */
				{
					/* ignore bars and beams */
					if (nod_e[master_c[j+count_master]+n_e]!=nod_e[master_c[j+count_master]+2*n_e])
					{
						/* rectangle */
						if (nod_e[master_c[j+count_master]+2*n_e]!=nod_e[master_c[j+count_master]+3*n_e])
						{
							ntriangle=2; /* 2 triangles below */
						}

						/* loop over triangles */
						for (ctriangle=0;ctriangle<ntriangle;ctriangle++)
						{
							if (ctriangle==0) // first triangle
							{
								/* nodes N1, N2, N3 */
								N1=nod_e[master_c[j+count_master]];
								N2=nod_e[master_c[j+count_master]+n_e];
								N3=nod_e[master_c[j+count_master]+2*n_e];
							}
							else // second triangle
							{
								/* nodes N1, N2, N3 */
								N1=nod_e[master_c[j+count_master]];
								N2=nod_e[master_c[j+count_master]+2*n_e];
								N3=nod_e[master_c[j+count_master]+3*n_e];
							}

							x1=x_b[N1];
							y1=x_b[N1+n_b];
							z1=x_b[N1+2*n_b];

							x2=x_b[N2];
							y2=x_b[N2+n_b];
							z2=x_b[N2+2*n_b];

							x3=x_b[N3];
							y3=x_b[N3+n_b];
							z3=x_b[N3+2*n_b];

							u1=x2-x1;
							u2=y2-y1;
							u3=z2-z1;

							v1=x3-x1;
							v2=y3-y1;
							v3=z3-z1;

							/* element normal vector */
							n1=u2*v3-u3*v2;
							n2=u3*v1-u1*v3;
							n3=u1*v2-u2*v1;

							/* element corridor side 1 */
							u1=x2-x3;
							u2=y2-y3;
							u3=z2-z3;

							n11=u2*n3-u3*n2;
							n12=u3*n1-u1*n3;
							n13=u1*n2-u2*n1;

							d1=(n11*x_s[slave_c[i+count_slave]]+n12*x_s[slave_c[i+count_slave]+n_s]+n13*x_s[slave_c[i+count_slave]+2*n_s]
								-n11*x2-n12*y2-n13*z2)/
								norm(n11,n12,n13);

							/* element corridor side 2 */
							u1=x3-x1;
							u2=y3-y1;
							u3=z3-z1;

							n11=u2*n3-u3*n2;
							n12=u3*n1-u1*n3;
							n13=u1*n2-u2*n1;
					
							d2=(n11*x_s[slave_c[i+count_slave]]+n12*x_s[slave_c[i+count_slave]+n_s]+n13*x_s[slave_c[i+count_slave]+2*n_s]
								-n11*x1-n12*y1-n13*z1)/
								norm(n11,n12,n13);

							/* element corridor side 3 */
							u1=x1-x2;
							u2=y1-y2;
							u3=z1-z2;

							n11=u2*n3-u3*n2;
							n12=u3*n1-u1*n3;
							n13=u1*n2-u2*n1;

							d3=(n11*x_s[slave_c[i+count_slave]]+n12*x_s[slave_c[i+count_slave]+n_s]+n13*x_s[slave_c[i+count_slave]+2*n_s]
								-n11*x1-n12*y1-n13*z1)/
								norm(n11,n12,n13);

							/* node located in element space defined by 3 element corridor sides */
							if ((d1>0.0)&&(d2>0.0)&&(d3>0.0))
							{
								/* normal vector*/
								rn=norm(n1,n2,n3);
								if (rn<EPS) rn=1.0; // zero normal vector

								/* penetration */
								p=ct+(n1*x_s[slave_c[i+count_slave]]+n2*x_s[slave_c[i+count_slave]+n_s]+n3*x_s[slave_c[i+count_slave]+2*n_s]
								-n1*x1-n2*y1-n3*z1)/rn;

								if ((p>0.0)&&(p<=2*ct))
								{
									/* projection of domain node */
									s=(n1*x1+n2*y1+n3*z1-n1*x_s[slave_c[i+count_slave]]-n2*x_s[slave_c[i+count_slave]+n_s]
										-n3*x_s[slave_c[i+count_slave]+2*n_s])/sqr(rn);

									x1=x_s[slave_c[i+count_slave]]+n1*s;
									x2=x_s[slave_c[i+count_slave]+n_s]+n2*s;
									x3=x_s[slave_c[i+count_slave]+2*n_s]+n3*s;

									u1=x1-x_b[N2];
									u2=x2-x_b[N2+n_b];
									u3=x3-x_b[N2+2*n_b];
									
									v1=x1-x_b[N3];
									v2=x2-x_b[N3+n_b];
									v3=x3-x_b[N3+2*n_b];
									
									n11=u2*v3-u3*v2;
									n12=-u1*v3+u3*v1;
									n13=u1*v2-u2*v1;

									k1=norm(n11,n12,n13)/rn;

									u1=x1-x_b[N1];
									u2=x2-x_b[N1+n_b];
									u3=x3-x_b[N1+2*n_b];

									v1=x1-x_b[N3];
									v2=x2-x_b[N3+n_b];
									v3=x3-x_b[N3+2*n_b];

									n11=u2*v3-v2*u3;
									n12=-u1*v3+v1*u3;
									n13=u1*v2-v1*u2;

									k2=norm(n11,n12,n13)/rn;

									u1=x1-x_b[N1];
									u2=x2-x_b[N1+n_b];
									u3=x3-x_b[N1+2*n_b];

									v1=x1-x_b[N2];
									v2=x2-x_b[N2+n_b];
									v3=x3-x_b[N2+2*n_b];

									n11=u2*v3-v2*u3;
									n12=-u1*v3+v1*u3;
									n13=u1*v2-v1*u2;

									k3=norm(n11,n12,n13)/rn;
									
									/* element velocity */
									vx1=v_b[N1];
									vy1=v_b[N1+n_b];
									vz1=v_b[N1+2*n_b];

									vx2=v_b[N2];
									vy2=v_b[N2+n_b];
									vz2=v_b[N2+2*n_b];

									vx3=v_b[N3];
									vy3=v_b[N3+n_b];
									vz3=v_b[N3+2*n_b];

									vt1=(vx1+vx2+vx3)/3.0;
									vt2=(vy1+vy2+vy3)/3.0;
									vt3=(vz1+vz2+vz3)/3.0;

									v1=v_s[slave_c[i+count_slave]]-vt1;
									v2=v_s[slave_c[i+count_slave]+n_s]-vt2;
									v3=v_s[slave_c[i+count_slave]+2*n_s]-vt3;

									u1=n2*v3-n3*v2;
									u2=n3*v1-n1*v3;
									u3=n1*v2-n2*v1;

									t1=u2*n3-n2*u3;
									t2=u3*n1-u1*n3;
									t3=u1*n2-u2*n1;

									/* unit tangent vector - velocity projection */
									rt=norm(t1,t2,t3);
									if (rt<EPS) // zero tangent vector
									{
										vr=0.0; // multiplier -> ft = 0
										rt=1.0; // because of division by rt
									}
									else vr=1.0; // analogy to 2D (relative motion)

									switch (sw)
									{
										case 1: /* soft contact */
										{
											if (klin>EPS)
											{
												fn=klin*p/ct+knon*p*p*p/(ct*ct*ct);
											}
											else
											{
												fn=m_slave[i+count_slave]*m_master[j+count_master]/(m_slave[i+count_slave]+m_master[j+count_master])
												/sqr(dt)*p;
							
												if (knon>EPS)
												{
													fn=fn+m_slave[i+count_slave]*m_master[j+count_master]/sqr(dt)*pow(p,knon);
												}
											}

											/* friction force */
											// if (rt>0) // for rt=0.0 -> rt=1.0
											ft=vr*kf*fn;

											/*
											f_s[slave_c[i+count_slave]]=f_s[slave_c[i+count_slave]]-(fn*n1+ft*t1/rt);
											f_s[slave_c[i+count_slave]+n_s]=f_s[slave_c[i+count_slave]+n_s]-(fn*n2+ft*t2/rt);
											f_s[slave_c[i+count_slave]+2*n_s]=f_s[slave_c[i+count_slave]+2*n_s]-(fn*n3+ft*t3/rt);
											*/

											f_s[slave_c[i+count_slave]]=f_s[slave_c[i+count_slave]]-(fn*n1/rn+ft*t1/rt);
											f_s[slave_c[i+count_slave]+n_s]=f_s[slave_c[i+count_slave]+n_s]-(fn*n2/rn+ft*t2/rt);
											f_s[slave_c[i+count_slave]+2*n_s]=f_s[slave_c[i+count_slave]+2*n_s]-(fn*n3/rn+ft*t3/rt);

											break;
										}
										case 2: /* sliding without separation */
										{
											v1=k1*v_b[N1]+
											k2*v_b[N2]+
											k3*v_b[N3];

											v2=k1*v_b[N1+n_b]+
											k2*v_b[N2+n_b]+
											k3*v_b[N3+n_b];

											v3=k1*v_b[N1+2*n_b]+
											k2*v_b[N2+2*n_b]+
											k3*v_b[N3+2*n_b];

											an_s=(n1*(v1-v_s[slave_c[i+count_slave]])+
												n2*(v2-v_s[slave_c[i+count_slave]+n_s])+
												n3*(v3-v_s[slave_c[i+count_slave]+2*n_s]))/rn/dt;
						
											at_s=(t1*f_s[slave_c[i+count_slave]]+
												t2*f_s[slave_c[i+count_slave]+n_s]+
												t3*f_s[slave_c[i+count_slave]+2*n_s])/m_slave[i+count_slave]/rt;

											if (klin>EPS)
											{
												fn=-klin*an_s-knon*an_s*an_s*an_s;
												ft=-klin*at_s-knon*at_s*at_s*at_s;
											}
											else
											{
												fn=-m_slave[i+count_slave]*(an_s-(n1*f_s[slave_c[i+count_slave]]+
													n2*f_s[slave_c[i+count_slave]+n_s]+
													n3*f_s[slave_c[i+count_slave]+2*n_s])/m_slave[i+count_slave]/rn);

												if (knon>EPS)
												{
													fn=fn-pow((n1*(at_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave])+
														n2*(an_s-f_s[slave_c[i+count_slave]+n_s]/m_slave[i+count_slave])+
														n3*(an_s-f_s[slave_c[i+count_slave]+2*n_s]/m_slave[i+count_slave]))/rn,knon);
												}
											}

											/* friction force */
											ft=kf*fn;

											f_s[slave_c[i+count_slave]]=m_slave[i+count_slave]*(n1*an_s/rn+t1*(at_s-ft)/rt);
											f_s[slave_c[i+count_slave]+n_s]=m_slave[i+count_slave]*(n2*an_s/rn+t2*(at_s-ft)/rt);
											f_s[slave_c[i+count_slave]+2*n_s]=m_slave[i+count_slave]*(n3*an_s/rn+t3*(at_s-ft)/rt);

											break;
										}
										case 3: /* tied contact */
										{
											v1=k1*v_b[N1]+
											k2*v_b[N2]+
											k3*v_b[N3];

											v2=k1*v_b[N1+n_b]+
											k2*v_b[N2+n_b]+
											k3*v_b[N3+n_b];

											v3=k1*v_b[N1+2*n_b]+
											k2*v_b[N2+2*n_b]+
											k3*v_b[N3+2*n_b];

											an_s=(n1*(v1-v_s[slave_c[i+count_slave]])+
												n2*(v2-v_s[slave_c[i+count_slave]+n_s])+
												n3*(v3-v_s[slave_c[i+count_slave]+2*n_s]))/rn/dt;
						
											at_s=(t1*(v1-v_s[slave_c[i+count_slave]])+
												t2*(v2-v_s[slave_c[i+count_slave]+n_s])+
												t3*(v3-v_s[slave_c[i+count_slave]+2*n_s]))/rt/dt;

											if (klin>EPS)
											{
												fn=-klin*an_s-knon*an_s*an_s*an_s;
												ft=-klin*at_s-knon*at_s*at_s*at_s;
											}
											else
											{
												fn=-m_slave[i+count_slave]*(an_s-(n1*f_s[slave_c[i+count_slave]]+
													n2*f_s[slave_c[i+count_slave]+n_s]+
													n3*f_s[slave_c[i+count_slave]+2*n_s])/m_slave[i+count_slave]/rn);
						
												ft=-m_slave[i+count_slave]*(at_s-(t1*f_s[slave_c[i+count_slave]]+
													t2*f_s[slave_c[i+count_slave]+n_s]+
													t3*f_s[slave_c[i+count_slave]+2*n_s])/m_slave[i+count_slave]/rt);

												if (knon>EPS)
												{
													fn=fn-m_slave[i+count_slave]*pow((n1*(at_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave])+
														n2*(an_s-f_s[slave_c[i+count_slave]+n_s]/m_slave[i+count_slave])+
														n3*(an_s-f_s[slave_c[i+count_slave]+2*n_s]/m_slave[i+count_slave]))
														/rn,knon);
							
													ft=ft-m_slave[i+count_slave]*pow((t1*(at_s-f_s[slave_c[i+count_slave]]/m_slave[i+count_slave])+
														t2*(an_s-f_s[slave_c[i+count_slave]+n_s]/m_slave[i+count_slave])+
														t3*(an_s-f_s[slave_c[i+count_slave]+2*n_s]/m_slave[i+count_slave]))
														/rt,knon);
												}
											}

											f_s[slave_c[i+count_slave]]=m_slave[i+count_slave]*(n1*an_s/rn+t1*at_s/rt);
											f_s[slave_c[i+count_slave]+n_s]=m_slave[i+count_slave]*(n2*an_s/rn+t2*at_s/rt);
											f_s[slave_c[i+count_slave]+2*n_s]=m_slave[i+count_slave]*(n3*an_s/rn+t3*at_s/rt);

											break;
										}
									}

									/* forces on boundary element */
									f_b[N1]=
										f_b[N1]+k1*(fn*n1/rn+ft*t1/rt);
									f_b[N1+n_b]=
										f_b[N1+n_b]+k1*(fn*n2/rn+ft*t2/rt);
									f_b[N1+2*n_b]=
										f_b[N1+2*n_b]+k1*(fn*n3/rn+ft*t3/rt);

									f_b[N2]=
										f_b[N2]+k2*(fn*n1/rn+ft*t1/rt);
									f_b[N2+n_b]=
										f_b[N2+n_b]+k2*(fn*n2/rn+ft*t2/rt);
									f_b[N2+2*n_b]=
										f_b[N2+2*n_b]+k2*(fn*n3/rn+ft*t3/rt);

									f_b[N3]=
										f_b[N3]+k3*(fn*n1/rn+ft*t1/rt);
									f_b[N3+n_b]=
										f_b[N3+n_b]+k3*(fn*n2/rn+ft*t2/rt);
									f_b[N3+2*n_b]=
										f_b[N3+2*n_b]+k3*(fn*n3/rn+ft*t3/rt);

									/*  contact force in case of ballistics (LHy for debugging)
									f_cx=f_cx+k1*(fn*n1/rn+ft*t1/rt)+k2*(fn*n1/rn+ft*t1/rt)+k3*(fn*n1/rn+ft*t1/rt);
									f_cy=f_cy+k1*(fn*n2/rn+ft*t2/rt)+k2*(fn*n2/rn+ft*t2/rt)+k3*(fn*n2/rn+ft*t2/rt);
									f_cz=f_cz+k1*(fn*n3/rn+ft*t3/rt)+k2*(fn*n3/rn+ft*t3/rt)+k3*(fn*n3/rn+ft*t3/rt);
									*/
					
									/* contribution to the contact force */
									f_c[k]=f_c[k]+fn*n1/rn+ft*t1/rt;
									f_c[k+n_c]=f_c[k+n_c]+fn*n2/rn+ft*t2/rt;
									f_c[k+2*n_c]=f_c[k+2*n_c]+fn*n3/rn+ft*t3/rt;
								}
							}
						}
					}    
					break;
				}
			}
		}
		//}
    }
}
