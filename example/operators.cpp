#include "operators.h"
#include<cmath>
void P2M_0(double x, double y, double z, double q, double * M) {
M[0] += q;

}
void M2M_0(double x, double y, double z, double * M, double * Ms) {
#pragma omp atomic
Ms[0] += M[0];

}

void M2L_0(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[1];
D[0] = (1 / (R));
#pragma omp atomic
L[0] += D[0]*M[0];

}

void L2L_0(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += L[0];

}

void L2P_0(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += L[0];
#pragma omp atomic
F[1] += 0;
#pragma omp atomic
F[2] += 0;
#pragma omp atomic
F[3] += 0;

}

void M2P_0(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*M[0]/(R*R*R);
#pragma omp atomic
F[0] += M[0]/R;
#pragma omp atomic
F[1] += Ftmp0*x;
#pragma omp atomic
F[2] += Ftmp0*y;
#pragma omp atomic
F[3] += Ftmp0*z;

}

void P2P(double x, double y, double z, double * S, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*S[0]/(R*R*R);
#pragma omp atomic
F[0] += S[0]/R;
#pragma omp atomic
F[1] += Ftmp0*x;
#pragma omp atomic
F[2] += Ftmp0*y;
#pragma omp atomic
F[3] += Ftmp0*z;

}

void P2M_1(double x, double y, double z, double q, double * M) {
M[0] += q;
M[1] += -q*x;
M[2] += -q*y;
M[3] += -q*z;

}
void M2M_1(double x, double y, double z, double * M, double * Ms) {
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += x*M[0] + M[1];
#pragma omp atomic
Ms[2] += y*M[0] + M[2];
#pragma omp atomic
Ms[3] += z*M[0] + M[3];

}

void M2L_1(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[4];
double Dtmp0 = 1.0/(R*R*R);
D[0] = (1 / (R));
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3];
#pragma omp atomic
L[1] += D[1]*M[0];
#pragma omp atomic
L[2] += D[2]*M[0];
#pragma omp atomic
L[3] += D[3]*M[0];

}

void L2L_1(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += x*L[1] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += L[1];
#pragma omp atomic
Ls[2] += L[2];
#pragma omp atomic
Ls[3] += L[3];

}

void L2P_1(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += x*L[1] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
F[1] += -L[1];
#pragma omp atomic
F[2] += -L[2];
#pragma omp atomic
F[3] += -L[3];

}

void M2P_1(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*M[1];
double Ftmp1 = (1 / (R*R));
double Ftmp2 = Ftmp1*x;
double Ftmp3 = 1.0*M[2];
double Ftmp4 = Ftmp1*y;
double Ftmp5 = 1.0*M[3];
double Ftmp6 = (1 / (R*R*R));
double Ftmp7 = 1.0*M[0];
double Ftmp8 = 3.0*Ftmp2;
double Ftmp9 = Ftmp8*y;
double Ftmp10 = z*M[3];
double Ftmp11 = 3.0*Ftmp1;
double Ftmp12 = 3.0*Ftmp4;
#pragma omp atomic
F[0] += (-Ftmp0*Ftmp2 - Ftmp1*Ftmp5*z - Ftmp3*Ftmp4 + M[0])/R;
#pragma omp atomic
F[1] += Ftmp6*(Ftmp0 - Ftmp10*Ftmp8 - Ftmp11*(x*x)*M[1] + Ftmp7*x - Ftmp9*M[2]);
#pragma omp atomic
F[2] += Ftmp6*(-Ftmp10*Ftmp12 - Ftmp11*(y*y)*M[2] + Ftmp3 + Ftmp7*y - Ftmp9*M[1]);
#pragma omp atomic
F[3] += Ftmp6*(-Ftmp11*(z*z)*M[3] - Ftmp12*z*M[2] + Ftmp5 + Ftmp7*z - Ftmp8*z*M[1]);

}

void P2M_2(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = (1.0/2.0)*q;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -q*z;
M[4] += Mtmp2*(x*x);
M[5] += Mtmp0*y;
M[6] += Mtmp0*z;
M[7] += Mtmp2*(y*y);
M[8] += Mtmp1*z;
M[9] += Mtmp2*(z*z);

}
void M2M_2(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = (1.0/2.0)*M[0];
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += z*M[0] + M[3];
#pragma omp atomic
Ms[4] += Mstmp2*(x*x) + x*M[1] + M[4];
#pragma omp atomic
Ms[5] += Mstmp0*y + x*M[2] + y*M[1] + M[5];
#pragma omp atomic
Ms[6] += Mstmp0*z + x*M[3] + z*M[1] + M[6];
#pragma omp atomic
Ms[7] += Mstmp2*(y*y) + y*M[2] + M[7];
#pragma omp atomic
Ms[8] += Mstmp1*z + y*M[3] + z*M[2] + M[8];
#pragma omp atomic
Ms[9] += Mstmp2*(z*z) + z*M[3] + M[9];

}

void M2L_2(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[10];
double Dtmp0 = (1 / (R*R*R));
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = 3.0/(R*R);
double Dtmp3 = (1 / (R*R*R*R*R));
double Dtmp4 = 3.0*Dtmp3*x;
D[0] = (1 / (R));
D[1] = -Dtmp1*x;
D[2] = -Dtmp1*y;
D[3] = -Dtmp1*z;
D[4] = Dtmp0*(Dtmp2*(x*x) - 1.0);
D[5] = Dtmp4*y;
D[6] = Dtmp4*z;
D[7] = Dtmp0*(Dtmp2*(y*y) - 1.0);
D[8] = 3.0*Dtmp3*y*z;
D[9] = -D[4] - D[7];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3];
#pragma omp atomic
L[4] += D[4]*M[0];
#pragma omp atomic
L[5] += D[5]*M[0];
#pragma omp atomic
L[6] += D[6]*M[0];
#pragma omp atomic
L[7] += D[7]*M[0];
#pragma omp atomic
L[8] += D[8]*M[0];
#pragma omp atomic
L[9] += D[9]*M[0];

}

void L2L_2(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp2*y + (1.0/2.0)*(x*x)*L[4] + x*L[1] + (1.0/2.0)*(y*y)*L[7] + y*L[2] + (1.0/2.0)*(z*z)*L[9] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp2 + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += L[4];
#pragma omp atomic
Ls[5] += L[5];
#pragma omp atomic
Ls[6] += L[6];
#pragma omp atomic
Ls[7] += L[7];
#pragma omp atomic
Ls[8] += L[8];
#pragma omp atomic
Ls[9] += L[9];

}

void L2P_2(double x, double y, double z, double * L, double * F) {
double Ftmp0 = y*L[5];
double Ftmp1 = z*L[6];
double Ftmp2 = z*L[8];
#pragma omp atomic
F[0] += Ftmp0*x + Ftmp1*x + Ftmp2*y + (1.0/2.0)*(x*x)*L[4] + x*L[1] + (1.0/2.0)*(y*y)*L[7] + y*L[2] + (1.0/2.0)*(z*z)*L[9] + z*L[3] + L[0];
#pragma omp atomic
F[1] += -Ftmp0 - Ftmp1 - x*L[4] - L[1];
#pragma omp atomic
F[2] += -Ftmp2 - x*L[5] - y*L[7] - L[2];
#pragma omp atomic
F[3] += -x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_2(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = (1 / (R*R));
double Ftmp1 = 1.0*M[1];
double Ftmp2 = 1.0*M[2];
double Ftmp3 = 1.0*M[3];
double Ftmp4 = (1 / (R*R*R*R));
double Ftmp5 = 3.0*Ftmp4;
double Ftmp6 = x*y;
double Ftmp7 = Ftmp6*M[5];
double Ftmp8 = x*M[6];
double Ftmp9 = Ftmp5*z;
double Ftmp10 = y*M[8];
double Ftmp11 = (x*x);
double Ftmp12 = 3.0*Ftmp0;
double Ftmp13 = Ftmp11*Ftmp12;
double Ftmp14 = Ftmp0*M[4];
double Ftmp15 = (y*y);
double Ftmp16 = Ftmp12*Ftmp15;
double Ftmp17 = Ftmp0*M[7];
double Ftmp18 = (z*z);
double Ftmp19 = Ftmp12*Ftmp18;
double Ftmp20 = Ftmp0*M[9];
double Ftmp21 = (1 / (R*R*R));
double Ftmp22 = 1.0*M[0];
double Ftmp23 = Ftmp12*M[5];
double Ftmp24 = Ftmp12*z;
double Ftmp25 = Ftmp12*Ftmp6;
double Ftmp26 = Ftmp24*M[3];
double Ftmp27 = 15.0*Ftmp4;
double Ftmp28 = Ftmp27*z;
double Ftmp29 = Ftmp11*Ftmp27;
double Ftmp30 = 15.0*Ftmp0;
double Ftmp31 = -Ftmp15*Ftmp30;
double Ftmp32 = Ftmp17*(Ftmp31 + 3.0);
double Ftmp33 = -Ftmp18*Ftmp30;
double Ftmp34 = Ftmp20*(Ftmp33 + 3.0);
double Ftmp35 = -Ftmp11*Ftmp30;
double Ftmp36 = Ftmp15*Ftmp27;
double Ftmp37 = Ftmp14*(Ftmp35 + 3.0);
double Ftmp38 = Ftmp18*Ftmp27;
#pragma omp atomic
F[0] += (-Ftmp0*Ftmp1*x - Ftmp0*Ftmp2*y - Ftmp0*Ftmp3*z + Ftmp10*Ftmp9 - Ftmp14*(-Ftmp13 + 1.0) - Ftmp17*(-Ftmp16 + 1.0) - Ftmp20*(-Ftmp19 + 1.0) + Ftmp5*Ftmp7 + Ftmp8*Ftmp9 + M[0])/R;
#pragma omp atomic
F[1] += Ftmp21*(Ftmp1 + Ftmp10*Ftmp28*x - Ftmp13*M[1] - Ftmp14*x*(Ftmp35 + 9.0) + Ftmp22*x - Ftmp23*y - Ftmp24*M[6] - Ftmp25*M[2] - Ftmp26*x + Ftmp29*y*M[5] + Ftmp29*z*M[6] - Ftmp32*x - Ftmp34*x);
#pragma omp atomic
F[2] += Ftmp21*(-Ftmp16*M[2] - Ftmp17*y*(Ftmp31 + 9.0) + Ftmp2 + Ftmp22*y - Ftmp23*x - Ftmp24*M[8] - Ftmp25*M[1] - Ftmp26*y + Ftmp28*Ftmp8*y - Ftmp34*y + Ftmp36*x*M[5] + Ftmp36*z*M[8] - Ftmp37*y);
#pragma omp atomic
F[3] += Ftmp21*(-Ftmp10*Ftmp12 + Ftmp10*Ftmp38 - Ftmp12*Ftmp8 - Ftmp19*M[3] - Ftmp20*z*(Ftmp33 + 9.0) + Ftmp22*z - Ftmp24*x*M[1] - Ftmp24*y*M[2] + Ftmp28*Ftmp7 + Ftmp3 - Ftmp32*z - Ftmp37*z + Ftmp38*Ftmp8);

}

void P2M_3(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = (x*x);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = (y*y);
double Mtmp7 = (z*z);
double Mtmp8 = (1.0/6.0)*q;
double Mtmp9 = (1.0/2.0)*Mtmp3;
double Mtmp10 = (1.0/2.0)*Mtmp0;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -Mtmp2;
M[4] += Mtmp3*Mtmp4;
M[5] += Mtmp5;
M[6] += Mtmp0*z;
M[7] += Mtmp4*Mtmp6;
M[8] += Mtmp1*z;
M[9] += Mtmp4*Mtmp7;
M[10] += -Mtmp8*(x*x*x);
M[11] += -Mtmp1*Mtmp9;
M[12] += -Mtmp2*Mtmp9;
M[13] += -Mtmp10*Mtmp6;
M[14] += -Mtmp5*z;
M[15] += -Mtmp10*Mtmp7;
M[16] += -Mtmp8*(y*y*y);
M[17] += -1.0/2.0*Mtmp2*Mtmp6;
M[18] += -1.0/2.0*Mtmp1*Mtmp7;
M[19] += -Mtmp8*(z*z*z);

}
void M2M_3(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = z*M[0];
double Mstmp3 = x*M[1];
double Mstmp4 = (x*x);
double Mstmp5 = (1.0/2.0)*M[0];
double Mstmp6 = x*M[2];
double Mstmp7 = y*M[1];
double Mstmp8 = Mstmp0*y;
double Mstmp9 = x*M[3];
double Mstmp10 = y*M[2];
double Mstmp11 = (y*y);
double Mstmp12 = y*M[3];
double Mstmp13 = (z*z);
double Mstmp14 = (1.0/2.0)*Mstmp4;
double Mstmp15 = (1.0/6.0)*M[0];
double Mstmp16 = (1.0/2.0)*M[1];
double Mstmp17 = (1.0/2.0)*Mstmp11;
double Mstmp18 = (1.0/2.0)*Mstmp13;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4*Mstmp5 + M[4];
#pragma omp atomic
Ms[5] += Mstmp6 + Mstmp7 + Mstmp8 + M[5];
#pragma omp atomic
Ms[6] += Mstmp0*z + Mstmp9 + z*M[1] + M[6];
#pragma omp atomic
Ms[7] += Mstmp10 + Mstmp11*Mstmp5 + M[7];
#pragma omp atomic
Ms[8] += Mstmp1*z + Mstmp12 + z*M[2] + M[8];
#pragma omp atomic
Ms[9] += Mstmp13*Mstmp5 + z*M[3] + M[9];
#pragma omp atomic
Ms[10] += Mstmp14*M[1] + Mstmp15*(x*x*x) + x*M[4] + M[10];
#pragma omp atomic
Ms[11] += Mstmp1*Mstmp14 + Mstmp14*M[2] + Mstmp3*y + x*M[5] + y*M[4] + M[11];
#pragma omp atomic
Ms[12] += Mstmp14*Mstmp2 + Mstmp14*M[3] + Mstmp3*z + x*M[6] + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += Mstmp0*Mstmp17 + Mstmp11*Mstmp16 + Mstmp6*y + x*M[7] + y*M[5] + M[13];
#pragma omp atomic
Ms[14] += Mstmp6*z + Mstmp7*z + Mstmp8*z + Mstmp9*y + x*M[8] + y*M[6] + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += Mstmp0*Mstmp18 + Mstmp13*Mstmp16 + Mstmp9*z + x*M[9] + z*M[6] + M[15];
#pragma omp atomic
Ms[16] += Mstmp15*(y*y*y) + Mstmp17*M[2] + y*M[7] + M[16];
#pragma omp atomic
Ms[17] += Mstmp10*z + Mstmp17*Mstmp2 + Mstmp17*M[3] + y*M[8] + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += Mstmp1*Mstmp18 + Mstmp12*z + Mstmp18*M[2] + y*M[9] + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += Mstmp15*(z*z*z) + Mstmp18*M[3] + z*M[9] + M[19];

}

void M2L_3(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[20];
double Dtmp0 = (1 / (R*R*R));
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = (x*x);
double Dtmp3 = (1 / (R*R));
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = (1 / (R*R*R*R*R));
double Dtmp6 = Dtmp5*x;
double Dtmp7 = 3.0*Dtmp6;
double Dtmp8 = (y*y);
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 15.0*Dtmp3;
double Dtmp11 = -Dtmp10*Dtmp2;
double Dtmp12 = Dtmp5*(Dtmp11 + 3.0);
double Dtmp13 = -Dtmp10*Dtmp8;
double Dtmp14 = Dtmp13 + 3.0;
D[0] = (1 / (R));
D[1] = -Dtmp1*x;
D[2] = -Dtmp1*y;
D[3] = -Dtmp1*z;
D[4] = Dtmp0*(Dtmp2*Dtmp4 - 1.0);
D[5] = Dtmp7*y;
D[6] = Dtmp7*z;
D[7] = Dtmp0*(Dtmp4*Dtmp8 - 1.0);
D[8] = 3.0*Dtmp9*z;
D[9] = -D[4] - D[7];
D[10] = Dtmp6*(Dtmp11 + 9.0);
D[11] = Dtmp12*y;
D[12] = Dtmp12*z;
D[13] = 1.0*Dtmp14*Dtmp6;
D[14] = -15.0*x*y*z/(R*R*R*R*R*R*R);
D[15] = -D[10] - D[13];
D[16] = Dtmp9*(Dtmp13 + 9.0);
D[17] = Dtmp14*Dtmp5*z;
D[18] = -D[11] - D[16];
D[19] = -D[12] - D[17];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[15]*M[9];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3] + D[11]*M[4] + D[13]*M[5] + D[14]*M[6] + D[16]*M[7] + D[17]*M[8] + D[18]*M[9];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3] + D[12]*M[4] + D[14]*M[5] + D[15]*M[6] + D[17]*M[7] + D[18]*M[8] + D[19]*M[9];
#pragma omp atomic
L[4] += D[4]*M[0] + D[10]*M[1] + D[11]*M[2] + D[12]*M[3];
#pragma omp atomic
L[5] += D[5]*M[0] + D[11]*M[1] + D[13]*M[2] + D[14]*M[3];
#pragma omp atomic
L[6] += D[6]*M[0] + D[12]*M[1] + D[14]*M[2] + D[15]*M[3];
#pragma omp atomic
L[7] += D[7]*M[0] + D[13]*M[1] + D[16]*M[2] + D[17]*M[3];
#pragma omp atomic
L[8] += D[8]*M[0] + D[14]*M[1] + D[17]*M[2] + D[18]*M[3];
#pragma omp atomic
L[9] += D[9]*M[0] + D[15]*M[1] + D[18]*M[2] + D[19]*M[3];
#pragma omp atomic
L[10] += D[10]*M[0];
#pragma omp atomic
L[11] += D[11]*M[0];
#pragma omp atomic
L[12] += D[12]*M[0];
#pragma omp atomic
L[13] += D[13]*M[0];
#pragma omp atomic
L[14] += D[14]*M[0];
#pragma omp atomic
L[15] += D[15]*M[0];
#pragma omp atomic
L[16] += D[16]*M[0];
#pragma omp atomic
L[17] += D[17]*M[0];
#pragma omp atomic
L[18] += D[18]*M[0];
#pragma omp atomic
L[19] += D[19]*M[0];

}

void L2L_3(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (1.0/2.0)*(x*x);
double Lstmp6 = (1.0/2.0)*(y*y);
double Lstmp7 = (1.0/2.0)*(z*z);
double Lstmp8 = x*L[13];
double Lstmp9 = x*L[15];
double Lstmp10 = y*L[11];
double Lstmp11 = z*L[12];
double Lstmp12 = y*L[18];
double Lstmp13 = z*L[17];
double Lstmp14 = y*L[13];
double Lstmp15 = y*L[14];
double Lstmp16 = z*L[15];
double Lstmp17 = z*L[18];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp5 + Lstmp11*Lstmp5 + Lstmp12*Lstmp7 + Lstmp13*Lstmp6 + Lstmp2*y + Lstmp4*x + Lstmp5*L[4] + Lstmp6*Lstmp8 + Lstmp6*L[7] + Lstmp7*Lstmp9 + Lstmp7*L[9] + (1.0/6.0)*(x*x*x)*L[10] + x*L[1] + (1.0/6.0)*(y*y*y)*L[16] + y*L[2] + (1.0/6.0)*(z*z*z)*L[19] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*x + Lstmp11*x + Lstmp4 + Lstmp5*L[10] + Lstmp6*L[13] + Lstmp7*L[15] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp13*y + Lstmp14*x + Lstmp2 + Lstmp3*x + Lstmp5*L[11] + Lstmp6*L[16] + Lstmp7*L[18] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp15*x + Lstmp16*x + Lstmp17*y + Lstmp5*L[12] + Lstmp6*L[17] + Lstmp7*L[19] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp10 + Lstmp11 + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp14 + Lstmp3 + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp15 + Lstmp16 + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp13 + Lstmp8 + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp17 + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp12 + Lstmp9 + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += L[10];
#pragma omp atomic
Ls[11] += L[11];
#pragma omp atomic
Ls[12] += L[12];
#pragma omp atomic
Ls[13] += L[13];
#pragma omp atomic
Ls[14] += L[14];
#pragma omp atomic
Ls[15] += L[15];
#pragma omp atomic
Ls[16] += L[16];
#pragma omp atomic
Ls[17] += L[17];
#pragma omp atomic
Ls[18] += L[18];
#pragma omp atomic
Ls[19] += L[19];

}

void L2P_3(double x, double y, double z, double * L, double * F) {
double Ftmp0 = y*L[5];
double Ftmp1 = z*L[6];
double Ftmp2 = z*L[8];
double Ftmp3 = x*y;
double Ftmp4 = Ftmp3*L[14];
double Ftmp5 = (1.0/2.0)*(x*x);
double Ftmp6 = (1.0/2.0)*(y*y);
double Ftmp7 = (1.0/2.0)*(z*z);
double Ftmp8 = Ftmp6*L[13];
double Ftmp9 = Ftmp7*L[15];
double Ftmp10 = Ftmp5*L[11];
double Ftmp11 = Ftmp5*L[12];
double Ftmp12 = Ftmp7*L[18];
double Ftmp13 = Ftmp6*L[17];
double Ftmp14 = x*z;
double Ftmp15 = y*z;
#pragma omp atomic
F[0] += Ftmp0*x + Ftmp1*x + Ftmp10*y + Ftmp11*z + Ftmp12*y + Ftmp13*z + Ftmp2*y + Ftmp4*z + Ftmp5*L[4] + Ftmp6*L[7] + Ftmp7*L[9] + Ftmp8*x + Ftmp9*x + (1.0/6.0)*(x*x*x)*L[10] + x*L[1] + (1.0/6.0)*(y*y*y)*L[16] + y*L[2] + (1.0/6.0)*(z*z*z)*L[19] + z*L[3] + L[0];
#pragma omp atomic
F[1] += -Ftmp0 - Ftmp1 - Ftmp14*L[12] - Ftmp15*L[14] - Ftmp3*L[11] - Ftmp5*L[10] - Ftmp8 - Ftmp9 - x*L[4] - L[1];
#pragma omp atomic
F[2] += -Ftmp10 - Ftmp12 - Ftmp14*L[14] - Ftmp15*L[17] - Ftmp2 - Ftmp3*L[13] - Ftmp6*L[16] - x*L[5] - y*L[7] - L[2];
#pragma omp atomic
F[3] += -Ftmp11 - Ftmp13 - Ftmp14*L[15] - Ftmp15*L[18] - Ftmp4 - Ftmp7*L[19] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_3(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = (1 / (R*R));
double Ftmp1 = 1.0*M[1];
double Ftmp2 = 1.0*M[2];
double Ftmp3 = 1.0*M[3];
double Ftmp4 = (1 / (R*R*R*R));
double Ftmp5 = Ftmp4*x;
double Ftmp6 = 3.0*Ftmp5;
double Ftmp7 = y*M[5];
double Ftmp8 = z*M[6];
double Ftmp9 = Ftmp4*y;
double Ftmp10 = Ftmp9*z;
double Ftmp11 = 15.0*M[14];
double Ftmp12 = (1 / (R*R*R*R*R*R));
double Ftmp13 = Ftmp12*y*z;
double Ftmp14 = (x*x);
double Ftmp15 = 3.0*Ftmp0;
double Ftmp16 = Ftmp14*Ftmp15;
double Ftmp17 = Ftmp0*M[4];
double Ftmp18 = (y*y);
double Ftmp19 = Ftmp15*Ftmp18;
double Ftmp20 = Ftmp0*M[7];
double Ftmp21 = (z*z);
double Ftmp22 = Ftmp15*Ftmp21;
double Ftmp23 = Ftmp0*M[9];
double Ftmp24 = 15.0*Ftmp0;
double Ftmp25 = -Ftmp14*Ftmp24;
double Ftmp26 = Ftmp25 + 9.0;
double Ftmp27 = Ftmp26*M[10];
double Ftmp28 = Ftmp25 + 3.0;
double Ftmp29 = Ftmp9*M[11];
double Ftmp30 = -Ftmp18*Ftmp24;
double Ftmp31 = Ftmp30 + 9.0;
double Ftmp32 = Ftmp31*M[16];
double Ftmp33 = Ftmp4*z;
double Ftmp34 = Ftmp30 + 3.0;
double Ftmp35 = -Ftmp21*Ftmp24;
double Ftmp36 = Ftmp35 + 9.0;
double Ftmp37 = Ftmp36*M[19];
double Ftmp38 = 1.0*Ftmp5;
double Ftmp39 = Ftmp38*M[13];
double Ftmp40 = Ftmp35 + 3.0;
double Ftmp41 = Ftmp38*M[15];
double Ftmp42 = 1.0*M[18];
double Ftmp43 = Ftmp42*Ftmp9;
double Ftmp44 = (1 / (R*R*R));
double Ftmp45 = 1.0*M[0];
double Ftmp46 = Ftmp15*M[5];
double Ftmp47 = Ftmp15*z;
double Ftmp48 = Ftmp15*x;
double Ftmp49 = Ftmp48*y;
double Ftmp50 = Ftmp47*M[3];
double Ftmp51 = y*M[8];
double Ftmp52 = Ftmp5*z;
double Ftmp53 = 15.0*Ftmp52;
double Ftmp54 = 15.0*Ftmp14;
double Ftmp55 = 105.0*M[14];
double Ftmp56 = Ftmp0*Ftmp34;
double Ftmp57 = Ftmp0*Ftmp40;
double Ftmp58 = Ftmp20*Ftmp34;
double Ftmp59 = Ftmp23*Ftmp40;
double Ftmp60 = 105.0*Ftmp0;
double Ftmp61 = -Ftmp14*Ftmp60;
double Ftmp62 = Ftmp61 + 45.0;
double Ftmp63 = Ftmp5*y;
double Ftmp64 = Ftmp62*Ftmp63;
double Ftmp65 = -Ftmp18*Ftmp60;
double Ftmp66 = Ftmp65 + 45.0;
double Ftmp67 = Ftmp66*M[16];
double Ftmp68 = Ftmp65 + 15.0;
double Ftmp69 = Ftmp68*M[17];
double Ftmp70 = Ftmp52*Ftmp62;
double Ftmp71 = -Ftmp21*Ftmp60;
double Ftmp72 = Ftmp71 + 45.0;
double Ftmp73 = Ftmp72*M[19];
double Ftmp74 = Ftmp14*Ftmp4;
double Ftmp75 = Ftmp71 + 15.0;
double Ftmp76 = Ftmp75*y;
double Ftmp77 = 1.0*Ftmp74;
double Ftmp78 = 15.0*Ftmp18;
double Ftmp79 = Ftmp12*Ftmp55*x;
double Ftmp80 = Ftmp0*Ftmp28;
double Ftmp81 = Ftmp17*Ftmp28;
double Ftmp82 = Ftmp61 + 15.0;
double Ftmp83 = Ftmp82*M[12];
double Ftmp84 = Ftmp18*Ftmp4;
double Ftmp85 = 15.0*Ftmp21;
double Ftmp86 = Ftmp21*Ftmp4;
double Ftmp87 = Ftmp72*z;
#pragma omp atomic
F[0] += (-Ftmp0*Ftmp1*x - Ftmp0*Ftmp2*y - Ftmp0*Ftmp3*z + 3.0*Ftmp10*M[8] - Ftmp11*Ftmp13*x - Ftmp17*(-Ftmp16 + 1.0) - Ftmp20*(-Ftmp19 + 1.0) - Ftmp23*(-Ftmp22 + 1.0) + Ftmp27*Ftmp5 + Ftmp28*Ftmp29 + Ftmp28*Ftmp33*M[12] + Ftmp32*Ftmp9 + Ftmp33*Ftmp34*M[17] + Ftmp33*Ftmp37 + Ftmp34*Ftmp39 + Ftmp40*Ftmp41 + Ftmp40*Ftmp43 + Ftmp6*Ftmp7 + Ftmp6*Ftmp8 + M[0])/R;
#pragma omp atomic
F[1] += Ftmp44*(-Ftmp0*Ftmp27 + Ftmp1 + Ftmp10*Ftmp11 - Ftmp13*Ftmp14*Ftmp55 - Ftmp16*M[1] - Ftmp17*Ftmp26*x + Ftmp33*Ftmp54*M[6] + Ftmp38*Ftmp76*M[18] + Ftmp45*x - Ftmp46*y - Ftmp47*M[6] - Ftmp49*M[2] - Ftmp50*x + Ftmp51*Ftmp53 + Ftmp52*Ftmp69 + Ftmp52*Ftmp73 + Ftmp54*Ftmp9*M[5] - Ftmp56*M[13] - Ftmp57*M[15] - Ftmp58*x - Ftmp59*x + Ftmp63*Ftmp67 + Ftmp64*M[11] + Ftmp68*Ftmp77*M[13] + Ftmp70*M[12] + Ftmp74*(Ftmp61 + 75.0)*M[10] + Ftmp75*Ftmp77*M[15]);
#pragma omp atomic
F[2] += Ftmp44*(-Ftmp0*Ftmp32 + Ftmp10*Ftmp66*M[17] + Ftmp10*Ftmp73 + Ftmp10*Ftmp83 + Ftmp11*Ftmp52 - Ftmp18*Ftmp79*z - Ftmp19*M[2] + Ftmp2 - Ftmp20*Ftmp31*y + Ftmp33*Ftmp78*M[8] + Ftmp39*Ftmp66*y + Ftmp41*Ftmp76 + Ftmp42*Ftmp75*Ftmp84 + Ftmp45*y - Ftmp46*x - Ftmp47*M[8] - Ftmp49*M[1] + Ftmp5*Ftmp78*M[5] - Ftmp50*y - Ftmp57*M[18] - Ftmp59*y + 15.0*Ftmp63*Ftmp8 + Ftmp64*M[10] - Ftmp80*M[11] - Ftmp81*y + Ftmp82*Ftmp84*M[11] + Ftmp84*(Ftmp65 + 75.0)*M[16]);
#pragma omp atomic
F[3] += Ftmp44*(-Ftmp0*Ftmp37 + Ftmp10*Ftmp67 + Ftmp11*Ftmp63 - Ftmp15*Ftmp51 - Ftmp21*Ftmp79*y - Ftmp22*M[3] - Ftmp23*Ftmp36*z + Ftmp29*Ftmp82*z + Ftmp3 + Ftmp39*Ftmp68*z + Ftmp41*Ftmp87 + Ftmp43*Ftmp87 + Ftmp45*z - Ftmp47*x*M[1] - Ftmp47*y*M[2] - Ftmp48*M[6] + Ftmp5*Ftmp85*M[6] + Ftmp53*Ftmp7 - Ftmp56*M[17] - Ftmp58*z + Ftmp69*Ftmp86 + Ftmp70*M[10] - Ftmp80*M[12] - Ftmp81*z + Ftmp83*Ftmp86 + Ftmp85*Ftmp9*M[8] + Ftmp86*(Ftmp71 + 75.0)*M[19]);

}

void P2M_4(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = (x*x);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = (y*y);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = (z*z);
double Mtmp10 = (x*x*x);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = (y*y*y);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = (z*z*z);
double Mtmp18 = (1.0/24.0)*q;
double Mtmp19 = (1.0/6.0)*Mtmp10;
double Mtmp20 = (1.0/4.0)*Mtmp3*q;
double Mtmp21 = (1.0/6.0)*Mtmp0;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -Mtmp2;
M[4] += Mtmp3*Mtmp4;
M[5] += Mtmp5;
M[6] += Mtmp6;
M[7] += Mtmp4*Mtmp7;
M[8] += Mtmp8;
M[9] += Mtmp4*Mtmp9;
M[10] += -Mtmp10*Mtmp11;
M[11] += -Mtmp1*Mtmp12;
M[12] += -Mtmp12*Mtmp2;
M[13] += -Mtmp13*Mtmp7;
M[14] += -Mtmp5*z;
M[15] += -Mtmp13*Mtmp9;
M[16] += -Mtmp11*Mtmp14;
M[17] += -Mtmp15*Mtmp2;
M[18] += -Mtmp1*Mtmp16;
M[19] += -Mtmp11*Mtmp17;
M[20] += Mtmp18*(x*x*x*x);
M[21] += Mtmp1*Mtmp19;
M[22] += Mtmp19*Mtmp2;
M[23] += Mtmp20*Mtmp7;
M[24] += Mtmp12*Mtmp8;
M[25] += Mtmp20*Mtmp9;
M[26] += Mtmp14*Mtmp21;
M[27] += Mtmp15*Mtmp6;
M[28] += Mtmp16*Mtmp5;
M[29] += Mtmp17*Mtmp21;
M[30] += Mtmp18*(y*y*y*y);
M[31] += (1.0/6.0)*Mtmp14*Mtmp2;
M[32] += (1.0/4.0)*Mtmp7*Mtmp9*q;
M[33] += (1.0/6.0)*Mtmp1*Mtmp17;
M[34] += Mtmp18*(z*z*z*z);

}
void M2M_4(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = z*M[0];
double Mstmp3 = x*M[1];
double Mstmp4 = (x*x);
double Mstmp5 = (1.0/2.0)*M[0];
double Mstmp6 = x*M[2];
double Mstmp7 = y*M[1];
double Mstmp8 = Mstmp0*y;
double Mstmp9 = x*M[3];
double Mstmp10 = z*M[1];
double Mstmp11 = Mstmp0*z;
double Mstmp12 = y*M[2];
double Mstmp13 = (y*y);
double Mstmp14 = y*M[3];
double Mstmp15 = z*M[2];
double Mstmp16 = Mstmp1*z;
double Mstmp17 = z*M[3];
double Mstmp18 = (z*z);
double Mstmp19 = x*M[4];
double Mstmp20 = (1.0/2.0)*Mstmp4;
double Mstmp21 = (x*x*x);
double Mstmp22 = (1.0/6.0)*M[0];
double Mstmp23 = x*M[5];
double Mstmp24 = y*M[4];
double Mstmp25 = Mstmp3*y;
double Mstmp26 = x*M[6];
double Mstmp27 = x*M[7];
double Mstmp28 = y*M[5];
double Mstmp29 = Mstmp6*y;
double Mstmp30 = (1.0/2.0)*M[1];
double Mstmp31 = (1.0/2.0)*Mstmp13;
double Mstmp32 = x*M[8];
double Mstmp33 = y*M[6];
double Mstmp34 = Mstmp9*y;
double Mstmp35 = x*M[9];
double Mstmp36 = (1.0/2.0)*Mstmp18;
double Mstmp37 = y*M[7];
double Mstmp38 = (y*y*y);
double Mstmp39 = y*M[8];
double Mstmp40 = y*M[9];
double Mstmp41 = (z*z*z);
double Mstmp42 = (1.0/6.0)*Mstmp21;
double Mstmp43 = (1.0/24.0)*M[0];
double Mstmp44 = (1.0/4.0)*Mstmp4*M[0];
double Mstmp45 = (1.0/6.0)*M[1];
double Mstmp46 = (1.0/6.0)*Mstmp38;
double Mstmp47 = (1.0/6.0)*Mstmp41;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4*Mstmp5 + M[4];
#pragma omp atomic
Ms[5] += Mstmp6 + Mstmp7 + Mstmp8 + M[5];
#pragma omp atomic
Ms[6] += Mstmp10 + Mstmp11 + Mstmp9 + M[6];
#pragma omp atomic
Ms[7] += Mstmp12 + Mstmp13*Mstmp5 + M[7];
#pragma omp atomic
Ms[8] += Mstmp14 + Mstmp15 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp17 + Mstmp18*Mstmp5 + M[9];
#pragma omp atomic
Ms[10] += Mstmp19 + Mstmp20*M[1] + Mstmp21*Mstmp22 + M[10];
#pragma omp atomic
Ms[11] += Mstmp1*Mstmp20 + Mstmp20*M[2] + Mstmp23 + Mstmp24 + Mstmp25 + M[11];
#pragma omp atomic
Ms[12] += Mstmp2*Mstmp20 + Mstmp20*M[3] + Mstmp26 + Mstmp3*z + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += Mstmp0*Mstmp31 + Mstmp13*Mstmp30 + Mstmp27 + Mstmp28 + Mstmp29 + M[13];
#pragma omp atomic
Ms[14] += Mstmp32 + Mstmp33 + Mstmp34 + Mstmp6*z + Mstmp7*z + Mstmp8*z + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += Mstmp0*Mstmp36 + Mstmp18*Mstmp30 + Mstmp35 + Mstmp9*z + z*M[6] + M[15];
#pragma omp atomic
Ms[16] += Mstmp22*Mstmp38 + Mstmp31*M[2] + Mstmp37 + M[16];
#pragma omp atomic
Ms[17] += Mstmp12*z + Mstmp2*Mstmp31 + Mstmp31*M[3] + Mstmp39 + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += Mstmp1*Mstmp36 + Mstmp14*z + Mstmp36*M[2] + Mstmp40 + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += Mstmp22*Mstmp41 + Mstmp36*M[3] + z*M[9] + M[19];
#pragma omp atomic
Ms[20] += Mstmp20*M[4] + Mstmp42*M[1] + Mstmp43*(x*x*x*x) + x*M[10] + M[20];
#pragma omp atomic
Ms[21] += Mstmp1*Mstmp42 + Mstmp19*y + Mstmp20*Mstmp7 + Mstmp20*M[5] + Mstmp42*M[2] + x*M[11] + y*M[10] + M[21];
#pragma omp atomic
Ms[22] += Mstmp10*Mstmp20 + Mstmp19*z + Mstmp2*Mstmp42 + Mstmp20*M[6] + Mstmp42*M[3] + x*M[12] + z*M[10] + M[22];
#pragma omp atomic
Ms[23] += Mstmp12*Mstmp20 + Mstmp13*Mstmp44 + Mstmp20*M[7] + Mstmp23*y + Mstmp3*Mstmp31 + Mstmp31*M[4] + x*M[13] + y*M[11] + M[23];
#pragma omp atomic
Ms[24] += Mstmp14*Mstmp20 + Mstmp15*Mstmp20 + Mstmp16*Mstmp20 + Mstmp20*M[8] + Mstmp23*z + Mstmp24*z + Mstmp25*z + Mstmp26*y + x*M[14] + y*M[12] + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += Mstmp17*Mstmp20 + Mstmp18*Mstmp44 + Mstmp20*M[9] + Mstmp26*z + Mstmp3*Mstmp36 + Mstmp36*M[4] + x*M[15] + z*M[12] + M[25];
#pragma omp atomic
Ms[26] += Mstmp0*Mstmp46 + Mstmp27*y + Mstmp31*Mstmp6 + Mstmp31*M[5] + Mstmp38*Mstmp45 + x*M[16] + y*M[13] + M[26];
#pragma omp atomic
Ms[27] += Mstmp10*Mstmp31 + Mstmp11*Mstmp31 + Mstmp27*z + Mstmp28*z + Mstmp29*z + Mstmp31*Mstmp9 + Mstmp31*M[6] + Mstmp32*y + x*M[17] + y*M[14] + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += Mstmp32*z + Mstmp33*z + Mstmp34*z + Mstmp35*y + Mstmp36*Mstmp6 + Mstmp36*Mstmp7 + Mstmp36*Mstmp8 + Mstmp36*M[5] + x*M[18] + y*M[15] + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += Mstmp0*Mstmp47 + Mstmp35*z + Mstmp36*Mstmp9 + Mstmp36*M[6] + Mstmp41*Mstmp45 + x*M[19] + z*M[15] + M[29];
#pragma omp atomic
Ms[30] += Mstmp31*M[7] + Mstmp43*(y*y*y*y) + Mstmp46*M[2] + y*M[16] + M[30];
#pragma omp atomic
Ms[31] += Mstmp15*Mstmp31 + Mstmp2*Mstmp46 + Mstmp31*M[8] + Mstmp37*z + Mstmp46*M[3] + y*M[17] + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += Mstmp12*Mstmp36 + (1.0/4.0)*Mstmp13*Mstmp18*M[0] + Mstmp17*Mstmp31 + Mstmp31*M[9] + Mstmp36*M[7] + Mstmp39*z + y*M[18] + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += Mstmp1*Mstmp47 + Mstmp14*Mstmp36 + Mstmp36*M[8] + Mstmp40*z + Mstmp47*M[2] + y*M[19] + z*M[18] + M[33];
#pragma omp atomic
Ms[34] += Mstmp36*M[9] + Mstmp43*(z*z*z*z) + Mstmp47*M[3] + z*M[19] + M[34];

}

void M2L_4(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[35];
double Dtmp0 = (1 / (R*R*R));
double Dtmp1 = 1.0*Dtmp0;
double Dtmp2 = (x*x);
double Dtmp3 = (1 / (R*R));
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = (1 / (R*R*R*R*R));
double Dtmp6 = Dtmp5*x;
double Dtmp7 = 3.0*Dtmp6;
double Dtmp8 = (y*y);
double Dtmp9 = Dtmp5*y;
double Dtmp10 = 15.0*Dtmp3;
double Dtmp11 = -Dtmp10*Dtmp2;
double Dtmp12 = Dtmp11 + 3.0;
double Dtmp13 = Dtmp12*Dtmp5;
double Dtmp14 = -Dtmp10*Dtmp8;
double Dtmp15 = Dtmp14 + 3.0;
double Dtmp16 = (1 / (R*R*R*R*R*R*R));
double Dtmp17 = Dtmp16*y;
double Dtmp18 = Dtmp17*z;
double Dtmp19 = 105.0/(R*R*R*R);
double Dtmp20 = Dtmp2*Dtmp3;
double Dtmp21 = -105.0*Dtmp20;
double Dtmp22 = x*(Dtmp21 + 45.0);
double Dtmp23 = Dtmp16*z;
double Dtmp24 = Dtmp3*Dtmp8;
double Dtmp25 = -105.0*Dtmp24;
double Dtmp26 = Dtmp25 + 45.0;
double Dtmp27 = 1.0*x;
D[0] = (1 / (R));
D[1] = -Dtmp1*x;
D[2] = -Dtmp1*y;
D[3] = -Dtmp1*z;
D[4] = Dtmp0*(Dtmp2*Dtmp4 - 1.0);
D[5] = Dtmp7*y;
D[6] = Dtmp7*z;
D[7] = Dtmp0*(Dtmp4*Dtmp8 - 1.0);
D[8] = 3.0*Dtmp9*z;
D[9] = -D[4] - D[7];
D[10] = Dtmp6*(Dtmp11 + 9.0);
D[11] = Dtmp13*y;
D[12] = Dtmp13*z;
D[13] = 1.0*Dtmp15*Dtmp6;
D[14] = -15.0*Dtmp18*x;
D[15] = -D[10] - D[13];
D[16] = Dtmp9*(Dtmp14 + 9.0);
D[17] = Dtmp15*Dtmp5*z;
D[18] = -D[11] - D[16];
D[19] = -D[12] - D[17];
D[20] = Dtmp5*(Dtmp19*(x*x*x*x) - 90.0*Dtmp20 + 9.0);
D[21] = -Dtmp17*Dtmp22;
D[22] = -Dtmp22*Dtmp23;
D[23] = Dtmp5*(Dtmp12 + Dtmp14 + Dtmp19*Dtmp2*Dtmp8);
D[24] = -Dtmp18*(Dtmp21 + 15.0);
D[25] = -D[20] - D[23];
D[26] = -Dtmp17*Dtmp26*Dtmp27;
D[27] = -Dtmp23*Dtmp27*(Dtmp25 + 15.0);
D[28] = -D[21] - D[26];
D[29] = -D[22] - D[27];
D[30] = Dtmp5*(Dtmp19*(y*y*y*y) - 90.0*Dtmp24 + 9.0);
D[31] = -Dtmp18*Dtmp26;
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -D[25] - D[32];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[15]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[29]*M[19];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3] + D[11]*M[4] + D[13]*M[5] + D[14]*M[6] + D[16]*M[7] + D[17]*M[8] + D[18]*M[9] + D[21]*M[10] + D[23]*M[11] + D[24]*M[12] + D[26]*M[13] + D[27]*M[14] + D[28]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[33]*M[19];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3] + D[12]*M[4] + D[14]*M[5] + D[15]*M[6] + D[17]*M[7] + D[18]*M[8] + D[19]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[29]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19];
#pragma omp atomic
L[4] += D[4]*M[0] + D[10]*M[1] + D[11]*M[2] + D[12]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[25]*M[9];
#pragma omp atomic
L[5] += D[5]*M[0] + D[11]*M[1] + D[13]*M[2] + D[14]*M[3] + D[21]*M[4] + D[23]*M[5] + D[24]*M[6] + D[26]*M[7] + D[27]*M[8] + D[28]*M[9];
#pragma omp atomic
L[6] += D[6]*M[0] + D[12]*M[1] + D[14]*M[2] + D[15]*M[3] + D[22]*M[4] + D[24]*M[5] + D[25]*M[6] + D[27]*M[7] + D[28]*M[8] + D[29]*M[9];
#pragma omp atomic
L[7] += D[7]*M[0] + D[13]*M[1] + D[16]*M[2] + D[17]*M[3] + D[23]*M[4] + D[26]*M[5] + D[27]*M[6] + D[30]*M[7] + D[31]*M[8] + D[32]*M[9];
#pragma omp atomic
L[8] += D[8]*M[0] + D[14]*M[1] + D[17]*M[2] + D[18]*M[3] + D[24]*M[4] + D[27]*M[5] + D[28]*M[6] + D[31]*M[7] + D[32]*M[8] + D[33]*M[9];
#pragma omp atomic
L[9] += D[9]*M[0] + D[15]*M[1] + D[18]*M[2] + D[19]*M[3] + D[25]*M[4] + D[28]*M[5] + D[29]*M[6] + D[32]*M[7] + D[33]*M[8] + D[34]*M[9];
#pragma omp atomic
L[10] += D[10]*M[0] + D[20]*M[1] + D[21]*M[2] + D[22]*M[3];
#pragma omp atomic
L[11] += D[11]*M[0] + D[21]*M[1] + D[23]*M[2] + D[24]*M[3];
#pragma omp atomic
L[12] += D[12]*M[0] + D[22]*M[1] + D[24]*M[2] + D[25]*M[3];
#pragma omp atomic
L[13] += D[13]*M[0] + D[23]*M[1] + D[26]*M[2] + D[27]*M[3];
#pragma omp atomic
L[14] += D[14]*M[0] + D[24]*M[1] + D[27]*M[2] + D[28]*M[3];
#pragma omp atomic
L[15] += D[15]*M[0] + D[25]*M[1] + D[28]*M[2] + D[29]*M[3];
#pragma omp atomic
L[16] += D[16]*M[0] + D[26]*M[1] + D[30]*M[2] + D[31]*M[3];
#pragma omp atomic
L[17] += D[17]*M[0] + D[27]*M[1] + D[31]*M[2] + D[32]*M[3];
#pragma omp atomic
L[18] += D[18]*M[0] + D[28]*M[1] + D[32]*M[2] + D[33]*M[3];
#pragma omp atomic
L[19] += D[19]*M[0] + D[29]*M[1] + D[33]*M[2] + D[34]*M[3];
#pragma omp atomic
L[20] += D[20]*M[0];
#pragma omp atomic
L[21] += D[21]*M[0];
#pragma omp atomic
L[22] += D[22]*M[0];
#pragma omp atomic
L[23] += D[23]*M[0];
#pragma omp atomic
L[24] += D[24]*M[0];
#pragma omp atomic
L[25] += D[25]*M[0];
#pragma omp atomic
L[26] += D[26]*M[0];
#pragma omp atomic
L[27] += D[27]*M[0];
#pragma omp atomic
L[28] += D[28]*M[0];
#pragma omp atomic
L[29] += D[29]*M[0];
#pragma omp atomic
L[30] += D[30]*M[0];
#pragma omp atomic
L[31] += D[31]*M[0];
#pragma omp atomic
L[32] += D[32]*M[0];
#pragma omp atomic
L[33] += D[33]*M[0];
#pragma omp atomic
L[34] += D[34]*M[0];

}

void L2L_4(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (x*x);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = (1.0/6.0)*(x*x*x);
double Lstmp8 = (y*y);
double Lstmp9 = (1.0/2.0)*Lstmp8;
double Lstmp10 = (1.0/6.0)*(y*y*y);
double Lstmp11 = (z*z);
double Lstmp12 = (1.0/2.0)*Lstmp11;
double Lstmp13 = (1.0/6.0)*(z*z*z);
double Lstmp14 = x*L[13];
double Lstmp15 = x*L[26];
double Lstmp16 = x*L[15];
double Lstmp17 = x*L[29];
double Lstmp18 = y*L[11];
double Lstmp19 = z*L[12];
double Lstmp20 = y*L[21];
double Lstmp21 = z*L[22];
double Lstmp22 = y*L[18];
double Lstmp23 = y*L[33];
double Lstmp24 = z*L[17];
double Lstmp25 = z*L[31];
double Lstmp26 = y*L[28];
double Lstmp27 = Lstmp26*x;
double Lstmp28 = z*L[27];
double Lstmp29 = Lstmp28*x;
double Lstmp30 = z*L[24];
double Lstmp31 = Lstmp30*y;
double Lstmp32 = (1.0/4.0)*Lstmp5;
double Lstmp33 = x*L[23];
double Lstmp34 = x*L[25];
double Lstmp35 = y*L[13];
double Lstmp36 = Lstmp28*y;
double Lstmp37 = x*L[28];
double Lstmp38 = y*L[23];
double Lstmp39 = y*L[32];
double Lstmp40 = y*L[14];
double Lstmp41 = z*L[15];
double Lstmp42 = z*L[18];
double Lstmp43 = z*L[28];
double Lstmp44 = Lstmp43*y;
double Lstmp45 = x*L[27];
double Lstmp46 = y*L[24];
double Lstmp47 = z*L[25];
double Lstmp48 = z*L[32];
double Lstmp49 = y*L[26];
double Lstmp50 = y*L[27];
double Lstmp51 = z*L[29];
double Lstmp52 = z*L[33];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp15 + Lstmp10*Lstmp25 + Lstmp10*L[16] + Lstmp11*Lstmp32*L[25] + (1.0/4.0)*Lstmp11*Lstmp8*L[32] + Lstmp12*Lstmp16 + Lstmp12*Lstmp22 + Lstmp12*Lstmp27 + Lstmp12*L[9] + Lstmp13*Lstmp17 + Lstmp13*Lstmp23 + Lstmp13*L[19] + Lstmp14*Lstmp9 + Lstmp18*Lstmp6 + Lstmp19*Lstmp6 + Lstmp2*y + Lstmp20*Lstmp7 + Lstmp21*Lstmp7 + Lstmp24*Lstmp9 + Lstmp29*Lstmp9 + Lstmp31*Lstmp6 + Lstmp32*Lstmp8*L[23] + Lstmp4*x + Lstmp6*L[4] + Lstmp7*L[10] + Lstmp9*L[7] + (1.0/24.0)*(x*x*x*x)*L[20] + x*L[1] + (1.0/24.0)*(y*y*y*y)*L[30] + y*L[2] + (1.0/24.0)*(z*z*z*z)*L[34] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*L[26] + Lstmp12*Lstmp26 + Lstmp12*Lstmp34 + Lstmp12*L[15] + Lstmp13*L[29] + Lstmp18*x + Lstmp19*x + Lstmp20*Lstmp6 + Lstmp21*Lstmp6 + Lstmp28*Lstmp9 + Lstmp31*x + Lstmp33*Lstmp9 + Lstmp4 + Lstmp6*L[10] + Lstmp7*L[20] + Lstmp9*L[13] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*L[30] + Lstmp12*Lstmp37 + Lstmp12*Lstmp39 + Lstmp12*L[18] + Lstmp13*L[33] + Lstmp15*Lstmp9 + Lstmp2 + Lstmp24*y + Lstmp25*Lstmp9 + Lstmp3*x + Lstmp30*Lstmp6 + Lstmp35*x + Lstmp36*x + Lstmp38*Lstmp6 + Lstmp6*L[11] + Lstmp7*L[21] + Lstmp9*L[16] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*L[31] + Lstmp12*Lstmp17 + Lstmp12*Lstmp23 + Lstmp12*L[19] + Lstmp13*L[34] + Lstmp40*x + Lstmp41*x + Lstmp42*y + Lstmp44*x + Lstmp45*Lstmp9 + Lstmp46*Lstmp6 + Lstmp47*Lstmp6 + Lstmp48*Lstmp9 + Lstmp6*L[12] + Lstmp7*L[22] + Lstmp9*L[17] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp12*L[25] + Lstmp18 + Lstmp19 + Lstmp20*x + Lstmp21*x + Lstmp31 + Lstmp6*L[20] + Lstmp9*L[23] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp12*L[28] + Lstmp3 + Lstmp30*x + Lstmp35 + Lstmp36 + Lstmp38*x + Lstmp6*L[21] + Lstmp9*L[26] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp12*L[29] + Lstmp40 + Lstmp41 + Lstmp44 + Lstmp46*x + Lstmp47*x + Lstmp6*L[22] + Lstmp9*L[27] + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp12*L[32] + Lstmp14 + Lstmp24 + Lstmp25*y + Lstmp29 + Lstmp49*x + Lstmp6*L[23] + Lstmp9*L[30] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp12*L[33] + Lstmp42 + Lstmp43*x + Lstmp48*y + Lstmp50*x + Lstmp6*L[24] + Lstmp9*L[31] + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp12*L[34] + Lstmp16 + Lstmp22 + Lstmp27 + Lstmp51*x + Lstmp52*y + Lstmp6*L[25] + Lstmp9*L[32] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp20 + Lstmp21 + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp30 + Lstmp38 + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp46 + Lstmp47 + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp28 + Lstmp33 + Lstmp49 + L[13];
#pragma omp atomic
Ls[14] += Lstmp43 + Lstmp50 + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp26 + Lstmp34 + Lstmp51 + L[15];
#pragma omp atomic
Ls[16] += Lstmp15 + Lstmp25 + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp45 + Lstmp48 + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp37 + Lstmp39 + Lstmp52 + L[18];
#pragma omp atomic
Ls[19] += Lstmp17 + Lstmp23 + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += L[20];
#pragma omp atomic
Ls[21] += L[21];
#pragma omp atomic
Ls[22] += L[22];
#pragma omp atomic
Ls[23] += L[23];
#pragma omp atomic
Ls[24] += L[24];
#pragma omp atomic
Ls[25] += L[25];
#pragma omp atomic
Ls[26] += L[26];
#pragma omp atomic
Ls[27] += L[27];
#pragma omp atomic
Ls[28] += L[28];
#pragma omp atomic
Ls[29] += L[29];
#pragma omp atomic
Ls[30] += L[30];
#pragma omp atomic
Ls[31] += L[31];
#pragma omp atomic
Ls[32] += L[32];
#pragma omp atomic
Ls[33] += L[33];
#pragma omp atomic
Ls[34] += L[34];

}

void L2P_4(double x, double y, double z, double * L, double * F) {
double Ftmp0 = y*L[5];
double Ftmp1 = z*L[6];
double Ftmp2 = z*L[8];
double Ftmp3 = x*y;
double Ftmp4 = Ftmp3*L[14];
double Ftmp5 = (x*x);
double Ftmp6 = (1.0/2.0)*Ftmp5;
double Ftmp7 = (1.0/6.0)*(x*x*x);
double Ftmp8 = (y*y);
double Ftmp9 = (1.0/2.0)*Ftmp8;
double Ftmp10 = (1.0/6.0)*(y*y*y);
double Ftmp11 = (z*z);
double Ftmp12 = (1.0/2.0)*Ftmp11;
double Ftmp13 = (1.0/6.0)*(z*z*z);
double Ftmp14 = Ftmp9*L[13];
double Ftmp15 = Ftmp10*L[26];
double Ftmp16 = Ftmp12*L[15];
double Ftmp17 = Ftmp13*L[29];
double Ftmp18 = Ftmp6*L[11];
double Ftmp19 = Ftmp6*L[12];
double Ftmp20 = Ftmp7*L[21];
double Ftmp21 = Ftmp7*L[22];
double Ftmp22 = Ftmp12*L[18];
double Ftmp23 = Ftmp13*L[33];
double Ftmp24 = Ftmp9*L[17];
double Ftmp25 = Ftmp10*L[31];
double Ftmp26 = x*z;
double Ftmp27 = y*z;
double Ftmp28 = (1.0/4.0)*Ftmp5;
double Ftmp29 = Ftmp3*z;
double Ftmp30 = Ftmp9*x;
double Ftmp31 = Ftmp12*x;
double Ftmp32 = Ftmp6*y;
double Ftmp33 = Ftmp6*z;
double Ftmp34 = Ftmp12*y;
double Ftmp35 = Ftmp9*z;
#pragma omp atomic
F[0] += Ftmp0*x + Ftmp1*x + Ftmp10*L[16] + Ftmp11*Ftmp28*L[25] + (1.0/4.0)*Ftmp11*Ftmp8*L[32] + Ftmp12*Ftmp3*L[28] + Ftmp12*L[9] + Ftmp13*L[19] + Ftmp14*x + Ftmp15*x + Ftmp16*x + Ftmp17*x + Ftmp18*y + Ftmp19*z + Ftmp2*y + Ftmp20*y + Ftmp21*z + Ftmp22*y + Ftmp23*y + Ftmp24*z + Ftmp25*z + Ftmp26*Ftmp9*L[27] + Ftmp27*Ftmp6*L[24] + Ftmp28*Ftmp8*L[23] + Ftmp4*z + Ftmp6*L[4] + Ftmp7*L[10] + Ftmp9*L[7] + (1.0/24.0)*(x*x*x*x)*L[20] + x*L[1] + (1.0/24.0)*(y*y*y*y)*L[30] + y*L[2] + (1.0/24.0)*(z*z*z*z)*L[34] + z*L[3] + L[0];
#pragma omp atomic
F[1] += -Ftmp0 - Ftmp1 - Ftmp14 - Ftmp15 - Ftmp16 - Ftmp17 - Ftmp26*L[12] - Ftmp27*L[14] - Ftmp29*L[24] - Ftmp3*L[11] - Ftmp30*L[23] - Ftmp31*L[25] - Ftmp32*L[21] - Ftmp33*L[22] - Ftmp34*L[28] - Ftmp35*L[27] - Ftmp6*L[10] - Ftmp7*L[20] - x*L[4] - L[1];
#pragma omp atomic
F[2] += -Ftmp10*L[30] - Ftmp18 - Ftmp2 - Ftmp20 - Ftmp22 - Ftmp23 - Ftmp26*L[14] - Ftmp27*L[17] - Ftmp29*L[27] - Ftmp3*L[13] - Ftmp30*L[26] - Ftmp31*L[28] - Ftmp32*L[23] - Ftmp33*L[24] - Ftmp34*L[32] - Ftmp35*L[31] - Ftmp9*L[16] - x*L[5] - y*L[7] - L[2];
#pragma omp atomic
F[3] += -Ftmp12*L[19] - Ftmp13*L[34] - Ftmp19 - Ftmp21 - Ftmp24 - Ftmp25 - Ftmp26*L[15] - Ftmp27*L[18] - Ftmp29*L[28] - Ftmp30*L[27] - Ftmp31*L[29] - Ftmp32*L[24] - Ftmp33*L[25] - Ftmp34*L[33] - Ftmp35*L[32] - Ftmp4 - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void M2P_4(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = (1 / (R*R));
double Ftmp1 = 1.0*M[1];
double Ftmp2 = 1.0*M[2];
double Ftmp3 = 1.0*M[3];
double Ftmp4 = (1 / (R*R*R*R));
double Ftmp5 = Ftmp4*x;
double Ftmp6 = 3.0*Ftmp5;
double Ftmp7 = y*M[5];
double Ftmp8 = z*M[6];
double Ftmp9 = Ftmp4*y;
double Ftmp10 = 3.0*Ftmp9;
double Ftmp11 = (1 / (R*R*R*R*R*R));
double Ftmp12 = Ftmp11*x;
double Ftmp13 = Ftmp12*y;
double Ftmp14 = z*M[14];
double Ftmp15 = (x*x);
double Ftmp16 = 3.0*Ftmp0;
double Ftmp17 = Ftmp15*Ftmp16;
double Ftmp18 = Ftmp0*M[4];
double Ftmp19 = (y*y);
double Ftmp20 = Ftmp16*Ftmp19;
double Ftmp21 = Ftmp0*M[7];
double Ftmp22 = (z*z);
double Ftmp23 = Ftmp16*Ftmp22;
double Ftmp24 = Ftmp0*M[9];
double Ftmp25 = 15.0*Ftmp0;
double Ftmp26 = -Ftmp15*Ftmp25;
double Ftmp27 = Ftmp26 + 9.0;
double Ftmp28 = Ftmp27*M[10];
double Ftmp29 = Ftmp26 + 3.0;
double Ftmp30 = Ftmp9*M[11];
double Ftmp31 = -Ftmp19*Ftmp25;
double Ftmp32 = Ftmp31 + 9.0;
double Ftmp33 = Ftmp32*M[16];
double Ftmp34 = Ftmp4*z;
double Ftmp35 = Ftmp31 + 3.0;
double Ftmp36 = -Ftmp22*Ftmp25;
double Ftmp37 = Ftmp36 + 9.0;
double Ftmp38 = Ftmp37*M[19];
double Ftmp39 = 1.0*Ftmp5;
double Ftmp40 = Ftmp39*M[13];
double Ftmp41 = Ftmp36 + 3.0;
double Ftmp42 = Ftmp39*M[15];
double Ftmp43 = 1.0*Ftmp9;
double Ftmp44 = Ftmp43*M[18];
double Ftmp45 = Ftmp0*Ftmp22;
double Ftmp46 = (-35.0*Ftmp45 + 5.0)*M[28];
double Ftmp47 = 3.0*Ftmp13;
double Ftmp48 = 105.0*Ftmp0;
double Ftmp49 = -Ftmp15*Ftmp48;
double Ftmp50 = Ftmp49 + 45.0;
double Ftmp51 = Ftmp50*M[21];
double Ftmp52 = -Ftmp19*Ftmp48;
double Ftmp53 = Ftmp52 + 45.0;
double Ftmp54 = Ftmp53*M[26];
double Ftmp55 = 1.0*Ftmp13;
double Ftmp56 = Ftmp52 + 15.0;
double Ftmp57 = Ftmp56*M[27];
double Ftmp58 = Ftmp12*z;
double Ftmp59 = 1.0*Ftmp58;
double Ftmp60 = Ftmp50*M[22];
double Ftmp61 = -Ftmp22*Ftmp48;
double Ftmp62 = Ftmp61 + 45.0;
double Ftmp63 = Ftmp62*M[29];
double Ftmp64 = Ftmp49 + 15.0;
double Ftmp65 = Ftmp64*M[24];
double Ftmp66 = Ftmp11*y;
double Ftmp67 = Ftmp66*z;
double Ftmp68 = Ftmp53*M[31];
double Ftmp69 = Ftmp62*M[33];
double Ftmp70 = (x*x*x*x);
double Ftmp71 = 105.0*Ftmp4;
double Ftmp72 = 90.0*Ftmp0;
double Ftmp73 = Ftmp4*M[20];
double Ftmp74 = (y*y*y*y);
double Ftmp75 = Ftmp4*M[30];
double Ftmp76 = (z*z*z*z);
double Ftmp77 = Ftmp4*M[34];
double Ftmp78 = Ftmp15*Ftmp71;
double Ftmp79 = Ftmp4*M[23];
double Ftmp80 = Ftmp4*M[25];
double Ftmp81 = Ftmp19*Ftmp22;
double Ftmp82 = Ftmp4*M[32];
double Ftmp83 = (1 / (R*R*R));
double Ftmp84 = 1.0*M[0];
double Ftmp85 = Ftmp16*M[5];
double Ftmp86 = Ftmp16*z;
double Ftmp87 = 15.0*Ftmp9;
double Ftmp88 = Ftmp16*x;
double Ftmp89 = Ftmp88*y;
double Ftmp90 = Ftmp86*M[3];
double Ftmp91 = y*M[8];
double Ftmp92 = Ftmp5*z;
double Ftmp93 = 15.0*Ftmp92;
double Ftmp94 = 15.0*Ftmp34;
double Ftmp95 = 105.0*M[14];
double Ftmp96 = Ftmp0*Ftmp35;
double Ftmp97 = Ftmp0*Ftmp41;
double Ftmp98 = 1.0*Ftmp34;
double Ftmp99 = Ftmp62*Ftmp98;
double Ftmp100 = Ftmp21*Ftmp35;
double Ftmp101 = Ftmp24*Ftmp41;
double Ftmp102 = Ftmp5*y;
double Ftmp103 = Ftmp102*Ftmp50;
double Ftmp104 = Ftmp53*M[16];
double Ftmp105 = Ftmp56*M[17];
double Ftmp106 = Ftmp50*Ftmp92;
double Ftmp107 = Ftmp62*M[19];
double Ftmp108 = Ftmp15*Ftmp4;
double Ftmp109 = Ftmp61 + 15.0;
double Ftmp110 = Ftmp109*M[18];
double Ftmp111 = 1.0*Ftmp108;
double Ftmp112 = Ftmp0*Ftmp15;
double Ftmp113 = -945.0*Ftmp112;
double Ftmp114 = Ftmp113 + 315.0;
double Ftmp115 = Ftmp13*z;
double Ftmp116 = Ftmp114*Ftmp115;
double Ftmp117 = Ftmp0*Ftmp19;
double Ftmp118 = -945.0*Ftmp117;
double Ftmp119 = Ftmp118 + 315.0;
double Ftmp120 = Ftmp119*M[31];
double Ftmp121 = -945.0*Ftmp45;
double Ftmp122 = Ftmp121 + 315.0;
double Ftmp123 = Ftmp55*z;
double Ftmp124 = -315.0*Ftmp45;
double Ftmp125 = 3.0*(Ftmp124 + 35.0)*M[28];
double Ftmp126 = 1.0*Ftmp15;
double Ftmp127 = Ftmp119*M[26];
double Ftmp128 = Ftmp15*(Ftmp113 + 525.0);
double Ftmp129 = Ftmp11*z;
double Ftmp130 = Ftmp126*Ftmp129;
double Ftmp131 = (Ftmp118 + 105.0)*M[27];
double Ftmp132 = Ftmp122*M[29];
double Ftmp133 = 945.0*Ftmp4;
double Ftmp134 = Ftmp133*Ftmp74;
double Ftmp135 = Ftmp75*(-630.0*Ftmp117 + Ftmp134 + 45.0);
double Ftmp136 = Ftmp133*Ftmp76;
double Ftmp137 = Ftmp77*(Ftmp136 - 630.0*Ftmp45 + 45.0);
double Ftmp138 = Ftmp133*Ftmp70;
double Ftmp139 = Ftmp133*Ftmp81;
double Ftmp140 = -315.0*Ftmp117;
double Ftmp141 = Ftmp133*Ftmp15;
double Ftmp142 = Ftmp141*Ftmp19;
double Ftmp143 = Ftmp141*Ftmp22;
double Ftmp144 = 15.0*Ftmp5;
double Ftmp145 = 15.0*Ftmp102;
double Ftmp146 = Ftmp0*Ftmp29;
double Ftmp147 = Ftmp18*Ftmp29;
double Ftmp148 = Ftmp64*M[12];
double Ftmp149 = Ftmp9*z;
double Ftmp150 = Ftmp19*Ftmp4;
double Ftmp151 = Ftmp12*Ftmp19;
double Ftmp152 = Ftmp114*Ftmp12;
double Ftmp153 = Ftmp118 + 525.0;
double Ftmp154 = Ftmp129*Ftmp19;
double Ftmp155 = (Ftmp113 + 105.0)*M[24];
double Ftmp156 = 1.0*M[33];
double Ftmp157 = Ftmp73*(-630.0*Ftmp112 + Ftmp138 + 45.0);
double Ftmp158 = -315.0*Ftmp112;
double Ftmp159 = Ftmp22*Ftmp4;
double Ftmp160 = Ftmp62*z;
double Ftmp161 = 1.0*Ftmp12*Ftmp22;
double Ftmp162 = Ftmp121 + 525.0;
double Ftmp163 = Ftmp22*Ftmp66;
#pragma omp atomic
F[0] += (-Ftmp0*Ftmp1*x - Ftmp0*Ftmp2*y - Ftmp0*Ftmp3*z + Ftmp10*z*M[8] - 15.0*Ftmp13*Ftmp14 - Ftmp13*Ftmp51 - Ftmp18*(-Ftmp17 + 1.0) - Ftmp21*(-Ftmp20 + 1.0) - Ftmp24*(-Ftmp23 + 1.0) + Ftmp28*Ftmp5 + Ftmp29*Ftmp30 + Ftmp29*Ftmp34*M[12] + Ftmp33*Ftmp9 + Ftmp34*Ftmp35*M[17] + Ftmp34*Ftmp38 + Ftmp35*Ftmp40 + Ftmp41*Ftmp42 + Ftmp41*Ftmp44 - Ftmp46*Ftmp47 - Ftmp54*Ftmp55 - Ftmp57*Ftmp59 - Ftmp58*Ftmp60 - Ftmp59*Ftmp63 + Ftmp6*Ftmp7 + Ftmp6*Ftmp8 - Ftmp65*Ftmp67 - Ftmp67*Ftmp68 - 1.0*Ftmp67*Ftmp69 + Ftmp73*(-Ftmp15*Ftmp72 + Ftmp70*Ftmp71 + 9.0) + Ftmp75*(-Ftmp19*Ftmp72 + Ftmp71*Ftmp74 + 9.0) + Ftmp77*(-Ftmp22*Ftmp72 + Ftmp71*Ftmp76 + 9.0) + Ftmp79*(Ftmp19*Ftmp78 + Ftmp29 + Ftmp31) + Ftmp80*(Ftmp22*Ftmp78 + Ftmp29 + Ftmp36) + Ftmp82*(Ftmp35 + Ftmp36 + Ftmp71*Ftmp81) + M[0])/R;
#pragma omp atomic
F[1] += Ftmp83*(-Ftmp0*Ftmp28 + Ftmp1 + Ftmp10*Ftmp46 - Ftmp100*x - Ftmp101*x + Ftmp102*Ftmp104 + Ftmp103*M[11] + Ftmp105*Ftmp92 + Ftmp106*M[12] + Ftmp107*Ftmp92 + Ftmp108*(Ftmp49 + 75.0)*M[10] + Ftmp109*Ftmp111*M[15] + Ftmp110*Ftmp39*y + Ftmp111*Ftmp56*M[13] - Ftmp115*Ftmp120 - Ftmp116*M[24] - Ftmp122*Ftmp123*M[33] - Ftmp125*Ftmp15*Ftmp66 - Ftmp126*Ftmp127*Ftmp66 - Ftmp128*Ftmp129*M[22] - Ftmp128*Ftmp66*M[21] - Ftmp130*Ftmp131 - Ftmp130*Ftmp132 + Ftmp135*x + Ftmp137*x + Ftmp14*Ftmp87 - Ftmp15*Ftmp67*Ftmp95 + Ftmp15*Ftmp87*M[5] + Ftmp15*Ftmp94*M[6] - Ftmp17*M[1] - Ftmp18*Ftmp27*x + Ftmp34*Ftmp60 + Ftmp43*Ftmp54 + Ftmp51*Ftmp9 + Ftmp57*Ftmp98 + Ftmp73*x*(-1050.0*Ftmp112 + Ftmp138 + 225.0) + Ftmp79*x*(Ftmp140 + Ftmp142 + Ftmp50) + Ftmp80*x*(Ftmp124 + Ftmp143 + Ftmp50) + Ftmp82*x*(Ftmp139 + Ftmp56 + Ftmp61) + Ftmp84*x - Ftmp85*y - Ftmp86*M[6] - Ftmp89*M[2] - Ftmp90*x + Ftmp91*Ftmp93 - Ftmp96*M[13] - Ftmp97*M[15] + Ftmp99*M[29]);
#pragma omp atomic
F[2] += Ftmp83*(-Ftmp0*Ftmp33 - Ftmp101*y + Ftmp103*M[10] + Ftmp107*Ftmp149 + Ftmp109*Ftmp42*y + 1.0*Ftmp110*Ftmp150 - Ftmp116*M[22] - Ftmp119*Ftmp123*M[27] - Ftmp122*Ftmp154*Ftmp156 - Ftmp123*Ftmp132 - Ftmp125*Ftmp151 + Ftmp137*y + Ftmp14*Ftmp144 + Ftmp144*Ftmp19*M[5] + Ftmp145*Ftmp8 - Ftmp146*M[11] - Ftmp147*y + Ftmp148*Ftmp149 + Ftmp149*Ftmp53*M[17] + Ftmp150*Ftmp64*M[11] + Ftmp150*(Ftmp52 + 75.0)*M[16] - 1.0*Ftmp151*Ftmp153*M[26] - Ftmp152*Ftmp19*M[21] - Ftmp153*Ftmp154*M[31] - Ftmp154*Ftmp155 + Ftmp157*y - Ftmp19*Ftmp58*Ftmp95 + Ftmp19*Ftmp94*M[8] + Ftmp2 - Ftmp20*M[2] - Ftmp21*Ftmp32*y + Ftmp34*Ftmp65 + Ftmp34*Ftmp68 + Ftmp39*Ftmp54 + Ftmp40*Ftmp53*y + Ftmp46*Ftmp6 + Ftmp5*Ftmp51 + Ftmp75*y*(-1050.0*Ftmp117 + Ftmp134 + 225.0) + Ftmp79*y*(Ftmp142 + Ftmp158 + Ftmp53) + Ftmp80*y*(Ftmp143 + Ftmp61 + Ftmp64) + Ftmp82*y*(Ftmp124 + Ftmp139 + Ftmp53) + Ftmp84*y - Ftmp85*x - Ftmp86*M[8] - Ftmp89*M[1] - Ftmp90*y - Ftmp97*M[18] + Ftmp99*M[33]);
#pragma omp atomic
F[3] += Ftmp83*(-Ftmp0*Ftmp38 - Ftmp100*z + Ftmp104*Ftmp149 + Ftmp105*Ftmp159 + Ftmp106*M[10] - Ftmp116*M[21] - Ftmp120*Ftmp163 - Ftmp123*Ftmp127 - Ftmp13*Ftmp22*Ftmp95 - Ftmp131*Ftmp161 + Ftmp135*z + Ftmp144*Ftmp22*M[6] + Ftmp145*M[14] - Ftmp146*M[12] - Ftmp147*z + Ftmp148*Ftmp159 - Ftmp152*Ftmp22*M[22] - Ftmp155*Ftmp163 - Ftmp156*Ftmp162*Ftmp163 + Ftmp157*z + Ftmp159*(Ftmp61 + 75.0)*M[19] - Ftmp16*Ftmp91 + Ftmp160*Ftmp42 + Ftmp160*Ftmp44 - Ftmp161*Ftmp162*M[29] + Ftmp22*Ftmp87*M[8] - Ftmp23*M[3] - Ftmp24*Ftmp37*z + Ftmp3 + Ftmp30*Ftmp64*z + Ftmp39*Ftmp57 + Ftmp39*Ftmp63 + Ftmp40*Ftmp56*z + Ftmp43*Ftmp69 - Ftmp47*z*(Ftmp124 + 105.0)*M[28] + Ftmp5*Ftmp60 + Ftmp65*Ftmp9 + Ftmp68*Ftmp9 + Ftmp7*Ftmp93 + Ftmp77*z*(Ftmp136 - 1050.0*Ftmp45 + 225.0) + Ftmp79*z*(Ftmp142 + Ftmp52 + Ftmp64) + Ftmp80*z*(Ftmp143 + Ftmp158 + Ftmp62) + Ftmp82*z*(Ftmp139 + Ftmp140 + Ftmp62) + Ftmp84*z - Ftmp86*x*M[1] - Ftmp86*y*M[2] - Ftmp88*M[6] - Ftmp96*M[17]);

}

void P2M(double x, double y, double z, double q, double * M, int order) {
switch (order) {
  case 0:
    P2M_0(x, y, z, q, M);
    break;
  case 1:
    P2M_1(x, y, z, q, M);
    break;
  case 2:
    P2M_2(x, y, z, q, M);
    break;
  case 3:
    P2M_3(x, y, z, q, M);
    break;
  case 4:
    P2M_4(x, y, z, q, M);
    break;
  }
}
void M2M(double x, double y, double z, double * M, double * Ms, int order) {
switch (order) {
  case 0:
    M2M_0(x, y, z, M, Ms);
    break;
  case 1:
    M2M_1(x, y, z, M, Ms);
    break;
  case 2:
    M2M_2(x, y, z, M, Ms);
    break;
  case 3:
    M2M_3(x, y, z, M, Ms);
    break;
  case 4:
    M2M_4(x, y, z, M, Ms);
    break;
  }
}
void M2L(double x, double y, double z, double * M, double * L, int order) {
switch (order) {
  case 0:
    M2L_0(x, y, z, M, L);
    break;
  case 1:
    M2L_1(x, y, z, M, L);
    break;
  case 2:
    M2L_2(x, y, z, M, L);
    break;
  case 3:
    M2L_3(x, y, z, M, L);
    break;
  case 4:
    M2L_4(x, y, z, M, L);
    break;
  }
}
void L2L(double x, double y, double z, double * L, double * Ls, int order) {
switch (order) {
  case 0:
    L2L_0(x, y, z, L, Ls);
    break;
  case 1:
    L2L_1(x, y, z, L, Ls);
    break;
  case 2:
    L2L_2(x, y, z, L, Ls);
    break;
  case 3:
    L2L_3(x, y, z, L, Ls);
    break;
  case 4:
    L2L_4(x, y, z, L, Ls);
    break;
  }
}
void L2P(double x, double y, double z, double * L, double * F, int order) {
switch (order) {
  case 0:
    L2P_0(x, y, z, L, F);
    break;
  case 1:
    L2P_1(x, y, z, L, F);
    break;
  case 2:
    L2P_2(x, y, z, L, F);
    break;
  case 3:
    L2P_3(x, y, z, L, F);
    break;
  case 4:
    L2P_4(x, y, z, L, F);
    break;
  }
}
void M2P(double x, double y, double z, double * M, double * F, int order) {
switch (order) {
  case 0:
    M2P_0(x, y, z, M, F);
    break;
  case 1:
    M2P_1(x, y, z, M, F);
    break;
  case 2:
    M2P_2(x, y, z, M, F);
    break;
  case 3:
    M2P_3(x, y, z, M, F);
    break;
  case 4:
    M2P_4(x, y, z, M, F);
    break;
  }
}
