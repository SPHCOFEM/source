#include "header.h"
#include "matrices.h"

void spatial_rotation(int dim,
     double *R,double *RO,double *ROO,
	double psi,double theta,double phi, // rotation angles
	double ox,double oy,double oz, // anglular velocities
	double ax,double ay,double az, // anglular accelerations
	double urx,double ury,double urz, // first rotation axis
	double vrx,double vry,double vrz, // second rotation axis
	double wrx,double wry,double wrz) // third rotation axis
{
     memset(R,0,9*sizeof(double)); // rotation matrix
     memset(RO,0,9*sizeof(double)); // velocity matrix
     memset(ROO,0,9*sizeof(double)); // acceleration matrix
    
     switch (dim)
     {
          case 1: // 1D
          {
               // no rotation in 1D
               break;
          }
          case 2: // 2D
          {
               R[0]=cos(psi);
               R[1]=-sin(psi);
               R[3]=sin(psi);
               R[4]=cos(psi);

               RO[0]=R[1]*oz;
               RO[1]=-R[0]*oz;
               RO[3]=R[4]*oz;
               RO[4]=-R[3]*oz;

               ROO[0]=-R[0]*sqr(oz)+R[1]*az;
               ROO[1]=-R[1]*sqr(oz)-R[0]*az;
               ROO[3]=-R[3]*sqr(oz)+R[4]*az;
               ROO[4]=-R[4]*sqr(oz)-R[3]*az;

               break;
          }
          case 3: // 3D
          {
               /* rigid body motion
               position X = Xc + R * X0 
               -> Rba = Rab ^ (-1) = Rab'

               https://en.wikipedia.org/wiki/Rotation_matrix 
               R(psi).x = u.(u.x) + cos(psi).(u cross x) cross u + sin(psi).(u cross x)
               -> for big arrays declare in main and use pointers
               -> for small arrays declare double R[9]
               */
               double R1[9]; // first (psi) rotation
               R1[0]=cos(psi)+sqr(urx)*(1.0-cos(psi));
               R1[1]=urx*ury*(1.0-cos(psi))-urz*sin(psi);
               R1[2]=urx*urz*(1.0-cos(psi))+ury*sin(psi);
               R1[3]=ury*urx*(1.0-cos(psi))+urz*sin(psi);
               R1[4]=cos(psi)+sqr(ury)*(1.0-cos(psi));
               R1[5]=ury*urz*(1.0-cos(psi))-urx*sin(psi);
               R1[6]=urz*urx*(1.0-cos(psi))-ury*sin(psi);
               R1[7]=urz*ury*(1.0-cos(psi))+urx*sin(psi);
               R1[8]=cos(psi)+sqr(urz)*(1.0-cos(psi));

               double R2[9]; // second (theta) rotation
               R2[0]=cos(theta)+sqr(vrx)*(1.0-cos(theta));
               R2[1]=vrx*vry*(1.0-cos(theta))-vrz*sin(theta);
               R2[2]=vrx*vrz*(1.0-cos(theta))+vry*sin(theta);
               R2[3]=vry*vrx*(1.0-cos(theta))+vrz*sin(theta);
               R2[4]=cos(theta)+sqr(vry)*(1.0-cos(theta));
               R2[5]=vry*vrz*(1.0-cos(theta))-vrx*sin(theta);
               R2[6]=vrz*vrx*(1.0-cos(theta))-vry*sin(theta);
               R2[7]=vrz*vry*(1.0-cos(theta))+vrx*sin(theta);
               R2[8]=cos(theta)+sqr(vrz)*(1.0-cos(theta));

               double R3[9]; // third (phi) rotation
               R3[0]=cos(phi)+sqr(wrx)*(1.0-cos(phi));
               R3[1]=wrx*wry*(1.0-cos(phi))-wrz*sin(phi);
               R3[2]=wrx*wrz*(1.0-cos(phi))+wry*sin(phi);
               R3[3]=wry*wrx*(1.0-cos(phi))+wrz*sin(phi);
               R3[4]=cos(phi)+sqr(wry)*(1.0-cos(phi));
               R3[5]=wry*wrz*(1.0-cos(phi))-wrx*sin(phi);
               R3[6]=wrz*wrx*(1.0-cos(phi))-wry*sin(phi);
               R3[7]=wrz*wry*(1.0-cos(phi))+wrx*sin(phi);
               R3[8]=cos(phi)+sqr(wrz)*(1.0-cos(phi));

               /* R(psi, theta, phi) = R(phi) * R(theta) * R(psi) */
               R[0]=R1[0]*(R2[0]*R3[0]+R2[3]*R3[1]+R2[6]*R3[2])+
                    R1[3]*(R2[1]*R3[0]+R2[4]*R3[1]+R2[7]*R3[2])+
                    R1[6]*(R2[2]*R3[0]+R2[5]*R3[1]+R2[8]*R3[2]);
               R[1]=R1[1]*(R2[0]*R3[0]+R2[3]*R3[1]+R2[6]*R3[2])+
                    R1[4]*(R2[1]*R3[0]+R2[4]*R3[1]+R2[7]*R3[2])+
                    R1[7]*(R2[2]*R3[0]+R2[5]*R3[1]+R2[8]*R3[2]);
               R[2]=R1[2]*(R2[0]*R3[0]+R2[3]*R3[1]+R2[6]*R3[2])+
                    R1[5]*(R2[1]*R3[0]+R2[4]*R3[1]+R2[7]*R3[2])+
                    R1[8]*(R2[2]*R3[0]+R2[5]*R3[1]+R2[8]*R3[2]);
               R[3]=R1[0]*(R2[0]*R3[3]+R2[3]*R3[4]+R2[6]*R3[5])+
                    R1[3]*(R2[1]*R3[3]+R2[4]*R3[4]+R2[7]*R3[5])+
                    R1[6]*(R2[2]*R3[3]+R2[5]*R3[4]+R2[8]*R3[5]);
               R[4]=R1[1]*(R2[0]*R3[3]+R2[3]*R3[4]+R2[6]*R3[5])+
                    R1[4]*(R2[1]*R3[3]+R2[4]*R3[4]+R2[7]*R3[5])+
                    R1[7]*(R2[2]*R3[3]+R2[5]*R3[4]+R2[8]*R3[5]);
               R[5]=R1[2]*(R2[0]*R3[3]+R2[3]*R3[4]+R2[6]*R3[5])+
                    R1[5]*(R2[1]*R3[3]+R2[4]*R3[4]+R2[7]*R3[5])+
                    R1[8]*(R2[2]*R3[3]+R2[5]*R3[4]+R2[8]*R3[5]);
               R[6]=R1[0]*(R2[0]*R3[6]+R2[3]*R3[7]+R2[6]*R3[8])+
                    R1[3]*(R2[1]*R3[6]+R2[4]*R3[7]+R2[7]*R3[8])+
                    R1[6]*(R2[2]*R3[6]+R2[5]*R3[7]+R2[8]*R3[8]);
               R[7]=R1[1]*(R2[0]*R3[6]+R2[3]*R3[7]+R2[6]*R3[8])+
                    R1[4]*(R2[1]*R3[6]+R2[4]*R3[7]+R2[7]*R3[8])+
                    R1[7]*(R2[2]*R3[6]+R2[5]*R3[7]+R2[8]*R3[8]);
               R[8]=R1[2]*(R2[0]*R3[6]+R2[3]*R3[7]+R2[6]*R3[8])+
                    R1[5]*(R2[1]*R3[6]+R2[4]*R3[7]+R2[7]*R3[8])+
                    R1[8]*(R2[2]*R3[6]+R2[5]*R3[7]+R2[8]*R3[8]);
               //R=mproduct(dim,mproduct(dim,R1,R2),R3);

               /* velocity V = Vc + dR * X0 = Vc + R * O * X0 = Vc + RO * X0
               -> O = angular velocity matrix = R' * dR = [  0, -oz,  oy]
                                                            [ oz,   0, -ox]
                                                            [-oy,  ox,   0]
               RO = R * O */
               RO[0]=R[1]*oz-R[2]*oy;
               RO[1]=R[2]*ox-R[0]*oz;
               RO[2]=R[0]*oy-R[1]*ox;
               RO[3]=R[4]*oz-R[5]*oy;
               RO[4]=R[5]*ox-R[3]*oz;
               RO[5]=R[3]*oy-R[4]*ox;
               RO[6]=R[7]*oz-R[8]*oy;
               RO[7]=R[8]*ox-R[6]*oz;
               RO[8]=R[6]*oy-R[7]*ox;

               /* acceleration A = Ac + R * (O ^ 2 + dO) * X0 = Ac + ROO * X0
               -> dO = angular acceleration matrix = [  0, -az,  ay]
                                                       [ az,   0, -ax]
                                                       [-ay,  ax,   0]
               ROO = R * (O ^ 2 + dO) */
               ROO[0]=R[1]*(az+ox*oy)-R[0]*(sqr(oy)+sqr(oz))-R[2]*(ay-ox*oz);
               ROO[1]=R[2]*(ax+oy*oz)-R[0]*(az-ox*oy)-R[1]*(sqr(ox)+sqr(oz));
               ROO[2]=R[0]*(ay+ox*oz)-R[2]*(sqr(ox)+sqr(oy))-R[1]*(ax-oy*oz);
               ROO[3]=R[4]*(az+ox*oy)-R[3]*(sqr(oy)+sqr(oz))-R[5]*(ay-ox*oz);
               ROO[4]=R[5]*(ax+oy*oz)-R[3]*(az-ox*oy)-R[4]*(sqr(ox)+sqr(oz));
               ROO[5]=R[3]*(ay+ox*oz)-R[5]*(sqr(ox)+sqr(oy))-R[4]*(ax-oy*oz);
               ROO[6]=R[7]*(az+ox*oy)-R[6]*(sqr(oy)+sqr(oz))-R[8]*(ay-ox*oz);
               ROO[7]=R[8]*(ax+oy*oz)-R[6]*(az-ox*oy)-R[7]*(sqr(ox)+sqr(oz));
               ROO[8]=R[6]*(ay+ox*oz)-R[8]*(sqr(ox)+sqr(oy))-R[7]*(ax-oy*oz);

               break;
          }
     }
}
