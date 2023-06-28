#include "header.h"

#undef __FUNC__
#define __FUNC__ "move"

void user_defined_material(int dim,int *num_m,int *type_m,
	double *rho_m,double *T_m,double *gamma_m,double *kappa_m,double *mu_m,
    double *p_s,
    double *e_s,double *S_s,//double *O_s,
    double *dedt_s,double *dSdt_s,
	double *coef1_m,double *coef2_m,double *coef3_m,double *coef4_m)
{
    /* tensors
    11 -> xx -> [i]
    12 -> xy -> [i+n_s]
    13 -> xz -> [i+2*n_s]
    22 -> yy -> [i+3*n_s]
    23 -> yz -> [i+4*n_s]
    33 -> zz -> [i+5*n_s]
    */

    /* rotation tensor
    Omega = | 0   Oxy Oxz |
            | Oyx   0 Oyz |
            | Ozx Ozy   0 |

    Oxx = 0.0;
    Oxy = O_s[i+n_s];
    Oxz = O_s[i+2*n_s];

    Oyx = - Oxy;
    Syy = 0.0;
    Oyz = O_s[i+4*n_s];

    Ozx = - Oxz;
    Ozy = - Oyz;
    Ozz = 0.0;
    */

    /* stress tensor
    Sigma = | Sxx Sxy Sxz |
            | Syx Syy Syz |
            | Szx Szy Szz |

    Sxx = - p_s[i] + S_s[i];
    Sxy = S_s[i+n_s];
    Sxz = S_s[i+2*n_s];

    Syx = Sxy;
    Syy = - p_s[i] + S_s[i+3*n_s];
    Syz = S_s[i+4*n_s];

    Szx = Sxz;
    Szy = Syz;
    Szz = - p_s[i] + S_s[i+5*n_s];
    */

    /* option to calculate stress tensor directly
    K = E / 3 / (1 - 2 * nu) // bulk modulus T_m
    G = E / 2 / (1 + nu)  // shear modulus mu_m

    Sigma_ij = S_ii * delta_ij + D_ij
        e_ij = e_ii * delta_ij + d_ij
        S_ii = K * e_ii
        D_ij = G * d_ij
        p = - S_ii / 3

    Sigma_ij = -p * delta_ij + S_ij
    */
   
    switch (dim)
    {
        case 1:
        {
            /* stress tensor direct calculation
            p_s[i]=-T_m[mat_s[i]]*e_s[i];
            S_s[i]=0.0;
            */
            break;
        }
        case 2:
        {
            /* stress tensor direct calculation
            e=(e_s[i]+e_s[i+3*n_s])/2.0;
            p_s[i]=-T_m[mat_s[i]]*e;
            S_s[i]=mu_m[mat_s[i]]*(e_s[i]-e);
            S_s[i+n_s]=mu_m[mat_s[i]]*(e_s[i+n_s]-e);
            S_s[i+3*n_s]=mu_m[mat_s[i]]*(e_s[i+3*n_s]-e);
            */
            break;
        }
        case 3:
        {
            /* stress tensor direct calculation
            e=(e_s[i]+e_s[i+3*n_s]+e_s[i+5*n_s])/3.0;
            p_s[i]=-T_m[mat_s[i]]*e;
            S_s[i]=mu_m[mat_s[i]]*(e_s[i]-e);
            S_s[i+n_s]=mu_m[mat_s[i]]*(e_s[i+n_s]-e);
            S_s[i+2*n_s]=mu_m[mat_s[i]]*(e_s[i+2*n_s]-e);
            S_s[i+3*n_s]=mu_m[mat_s[i]]*(e_s[i+3*n_s]-e);
            S_s[i+4*n_s]=mu_m[mat_s[i]]*(e_s[i+4*n_s]-e);
            S_s[i+5*n_s]=mu_m[mat_s[i]]*(e_s[i+5*n_s]-e);
            */
            break;
        }
    }

    /* deformation tensor trace 
    de=(dedt_s[i]+dedt_s[i+3*n_s]+dedt_s[i+5*n_s])/(double)dim; */

    /* deviatoric stress rate tensor based on deformation rate (Fuller, 2010)
    dSdt_ij = 2 * mu * (dedt_ij - delta_ij * dedt_ii / 3) + O_ik * S_kj - S_ik * O_kj */

    /* From symbolic Matlab:
    dS_11=2*O_12*S_12+2*O_13*S_13-2*mu*(e_11/3-de_11+e_22/3+e_33/3);
    dS_12=O_12*S_22-O_12*S_11+O_13*S_23+O_23*S_13+2*de_12*mu;
    dS_13=O_12*S_23-O_13*S_11-O_23*S_12+O_13*S_33+2*de_13*mu;
    dS_22=2*O_23*S_23-2*O_12*S_12-2*mu*(e_11/3-de_22+e_22/3+e_33/3);
    dS_23=O_23*S_33-O_13*S_12-O_23*S_22-O_12*S_13+2*de_23*mu;
    dS_33=-2*O_13*S_13-2*O_23*S_23-2*mu*(e_11/3-de_33+e_22/3+e_33/3);

    E=createMemMore(double,9);
    E[0]=dedt_s[i]-de/3.0;E[1]=dedt_s[i+n_s];E[2]=dedt_s[i+2*n_s];
    E[3]=de[1];E[4]=dedt_s[i+3*n_s]-de/3.0;E[5]=dedt_s[i+4*n_s];
    E[6]=de[2];E[7]=de[5];E[8]=dedt_s[i+5*n_s]-de/3.0;

    O=createMemMore(double,9);
    O[0]=O_s[i];O[1]=O_s[i+n_s];O[2]=O_s[i+2*n_s];
    O[3]=-O[1];O[4]=O_s[i+3*n_s];O[5]=O_s[i+4*n_s];
    O[6]=-O[2];O[7]=-O[5];O[8]=O_s[i+5*n_s];

    S=createMemMore(double,9);
    S[0]=S_s[i];S[1]=S_s[i+n_s];S[2]=S_s[i+2*n_s];
    S[3]=S[1];S[4]=S_s[i+3*n_s];S[5]=S_s[i+4*n_s];
    S[6]=S[2];S[7]=S[5];S[8]=S_s[i+5*n_s];

    dS=msubtraction(maddition(9,mmultiplication(9,2.0*mu,E),mproduct(9,O,S)),mproduct(9,S,O));
    */
}

