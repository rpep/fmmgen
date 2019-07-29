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
#pragma omp atomic
L[0] += M[0]/R;

}

void L2L_0(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += L[0];

}

void L2P_0(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += L[0];

}

void M2P_0(double x, double y, double z, double * M, double * F) {
#pragma omp atomic
F[0] += pow((x*x) + (y*y) + (z*z), -0.5)*M[0];

}

void P2P(double x, double y, double z, double * S, double * F) {
double R = sqrt(x*x + y*y + z*z);
#pragma omp atomic
F[0] += S[0]/R;

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
double x0;
double x1;
double x2;
double x3;
x0 = 1.0/(R*R*R);
x1 = x*x0;
x2 = x0*y;
x3 = x0*z;
#pragma omp atomic
L[0] += -x1*M[1] - x2*M[2] - x3*M[3] + M[0]/R;
#pragma omp atomic
L[1] += -x1*M[0];
#pragma omp atomic
L[2] += -x2*M[0];
#pragma omp atomic
L[3] += -x3*M[0];

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

}

void M2P_1(double x, double y, double z, double * M, double * F) {
double x0;
double x1;
x0 = (x*x) + (y*y) + (z*z);
x1 = 1.0*pow(x0, -1.5);
#pragma omp atomic
F[0] += -x*x1*M[1] + pow(x0, -0.5)*M[0] - x1*y*M[2] - x1*z*M[3];

}

void P2M_2(double x, double y, double z, double q, double * M) {
double x0;
double x1;
x0 = q*x;
x1 = q*y;
M[0] += q;
M[1] += -x0;
M[2] += -x1;
M[3] += -q*z;
M[4] += q*(x*x);
M[5] += x0*y;
M[6] += x0*z;
M[7] += q*(y*y);
M[8] += x1*z;
M[9] += q*(z*z);
}
void M2M_2(double x, double y, double z, double * M, double * Ms) {
double x0;
double x1;
x0 = x*M[0];
x1 = y*M[0];
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += x0 + M[1];
#pragma omp atomic
Ms[2] += x1 + M[2];
#pragma omp atomic
Ms[3] += z*M[0] + M[3];
#pragma omp atomic
Ms[4] += (x*x)*M[0] + x*M[1] + M[4];
#pragma omp atomic
Ms[5] += x*M[2] + x0*y + y*M[1] + M[5];
#pragma omp atomic
Ms[6] += x*M[3] + x0*z + z*M[1] + M[6];
#pragma omp atomic
Ms[7] += (y*y)*M[0] + y*M[2] + M[7];
#pragma omp atomic
Ms[8] += x1*z + y*M[3] + z*M[2] + M[8];
#pragma omp atomic
Ms[9] += (z*z)*M[0] + z*M[3] + M[9];

}

void M2L_2(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*x1;
x3 = x1*y;
x4 = x1*z;
x5 = 3.0/(R*R*R*R*R);
x6 = x*x5;
x7 = x6*y;
x8 = x6*z;
x9 = x5*y*z;
x10 = -x1;
x11 = (x*x)*x5;
x12 = x10 + x11;
x13 = x5*(y*y);
x14 = x10 + x13;
x15 = 2.0*x0 - x11 - x13;
#pragma omp atomic
L[0] += (1.0/2.0)*x12*M[4] + (1.0/2.0)*x14*M[7] + (1.0/2.0)*x15*M[9] - x2*M[1] - x3*M[2] - x4*M[3] + x7*M[5] + x8*M[6] + x9*M[8] + M[0]/R;
#pragma omp atomic
L[1] += x12*M[1] - x2*M[0] + x7*M[2] + x8*M[3];
#pragma omp atomic
L[2] += x14*M[2] - x3*M[0] + x7*M[1] + x9*M[3];
#pragma omp atomic
L[3] += x15*M[3] - x4*M[0] + x8*M[1] + x9*M[2];
#pragma omp atomic
L[4] += x12*M[0];
#pragma omp atomic
L[5] += x7*M[0];
#pragma omp atomic
L[6] += x8*M[0];
#pragma omp atomic
L[7] += x14*M[0];
#pragma omp atomic
L[8] += x9*M[0];
#pragma omp atomic
L[9] += x15*M[0];

}

void L2L_2(double x, double y, double z, double * L, double * Ls) {
double x0;
double x1;
double x2;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
#pragma omp atomic
Ls[0] += (1.0/2.0)*(x*x)*L[4] + x*x0 + x*x1 + x*L[1] + x2*y + (1.0/2.0)*(y*y)*L[7] + y*L[2] + (1.0/2.0)*(z*z)*L[9] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*L[4] + x0 + x1 + L[1];
#pragma omp atomic
Ls[2] += x*L[5] + x2 + y*L[7] + L[2];
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
#pragma omp atomic
F[0] += (1.0/2.0)*(x*x)*L[4] + x*y*L[5] + x*z*L[6] + x*L[1] + (1.0/2.0)*(y*y)*L[7] + y*z*L[8] + y*L[2] + (1.0/2.0)*(z*z)*L[9] + z*L[3] + L[0];

}

void M2P_2(double x, double y, double z, double * M, double * F) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = 3.0*pow(x2, -2.5);
x6 = x*x5;
x7 = -x4;
x8 = x0*x5;
x9 = x1*x5;
#pragma omp atomic
F[0] += -x*x4*M[1] + pow(x2, -0.5)*M[0] - x4*y*M[2] - x4*z*M[3] + x5*y*z*M[8] + x6*y*M[5] + x6*z*M[6] + (1.0/2.0)*(x7 + x8)*M[4] + (1.0/2.0)*(x7 + x9)*M[7] + (1.0/2.0)*(2.0*x3 - x8 - x9)*M[9];

}

void P2M_3(double x, double y, double z, double q, double * M) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = x0*y;
x5 = (y*y);
x6 = (z*z);
M[0] += q;
M[1] += -x0;
M[2] += -x1;
M[3] += -x2;
M[4] += q*x3;
M[5] += x4;
M[6] += x0*z;
M[7] += q*x5;
M[8] += x1*z;
M[9] += q*x6;
M[10] += -q*(x*x*x);
M[11] += -x1*x3;
M[12] += -x2*x3;
M[13] += -x0*x5;
M[14] += -x4*z;
M[15] += -x0*x6;
M[16] += -q*(y*y*y);
M[17] += -x2*x5;
M[18] += -x1*x6;
M[19] += -q*(z*z*z);
}
void M2M_3(double x, double y, double z, double * M, double * Ms) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
x0 = x*M[0];
x1 = y*M[0];
x2 = z*M[0];
x3 = x*M[1];
x4 = (x*x);
x5 = x*M[2];
x6 = y*M[1];
x7 = x0*y;
x8 = x*M[3];
x9 = y*M[2];
x10 = (y*y);
x11 = y*M[3];
x12 = (z*z);
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += x0 + M[1];
#pragma omp atomic
Ms[2] += x1 + M[2];
#pragma omp atomic
Ms[3] += x2 + M[3];
#pragma omp atomic
Ms[4] += x3 + x4*M[0] + M[4];
#pragma omp atomic
Ms[5] += x5 + x6 + x7 + M[5];
#pragma omp atomic
Ms[6] += x0*z + x8 + z*M[1] + M[6];
#pragma omp atomic
Ms[7] += x10*M[0] + x9 + M[7];
#pragma omp atomic
Ms[8] += x1*z + x11 + z*M[2] + M[8];
#pragma omp atomic
Ms[9] += x12*M[0] + z*M[3] + M[9];
#pragma omp atomic
Ms[10] += (x*x*x)*M[0] + x*M[4] + x4*M[1] + M[10];
#pragma omp atomic
Ms[11] += x*M[5] + x1*x4 + x3*y + x4*M[2] + y*M[4] + M[11];
#pragma omp atomic
Ms[12] += x*M[6] + x2*x4 + x3*z + x4*M[3] + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += x*M[7] + x0*x10 + x10*M[1] + x5*y + y*M[5] + M[13];
#pragma omp atomic
Ms[14] += x*M[8] + x5*z + x6*z + x7*z + x8*y + y*M[6] + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += x*M[9] + x0*x12 + x12*M[1] + x8*z + z*M[6] + M[15];
#pragma omp atomic
Ms[16] += x10*M[2] + (y*y*y)*M[0] + y*M[7] + M[16];
#pragma omp atomic
Ms[17] += x10*x2 + x10*M[3] + x9*z + y*M[8] + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += x1*x12 + x11*z + x12*M[2] + y*M[9] + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += x12*M[3] + (z*z*z)*M[0] + z*M[9] + M[19];

}

void M2L_3(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*x1;
x3 = x1*y;
x4 = x1*z;
x5 = (1 / (R*R*R*R*R));
x6 = 3.0*x5;
x7 = x*x6;
x8 = x7*y;
x9 = x7*z;
x10 = y*z;
x11 = x10*x6;
x12 = 15.0/pow(R, 7);
x13 = x*x10*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x6;
x17 = x14 + x16;
x18 = (1.0/2.0)*M[4];
x19 = (y*y);
x20 = x19*x6;
x21 = x14 + x20;
x22 = (1.0/2.0)*M[7];
x23 = x12*x19;
x24 = x*(x23 - x6);
x25 = 0.5*x24;
x26 = 9.0*x5;
x27 = -x12*x15;
x28 = x*(x26 + x27);
x29 = x27 + x6;
x30 = x29*y;
x31 = -x23;
x32 = y*(x26 + x31);
x33 = x29*z;
x34 = z*(x31 + x6);
x35 = 2.0*x0 - x16 - x20;
x36 = (1.0/2.0)*M[9];
x37 = 1.0*x24;
x38 = -x28 + x37;
x39 = (1.0/2.0)*x38;
x40 = -x30 - x32;
x41 = (1.0/2.0)*x40;
x42 = -x33 - x34;
#pragma omp atomic
L[0] += x11*M[8] - x13*M[14] + x17*x18 - x2*M[1] + x21*x22 - x25*M[13] + (1.0/6.0)*x28*M[10] - x3*M[2] + (1.0/2.0)*x30*M[11] + (1.0/6.0)*x32*M[16] + (1.0/2.0)*x33*M[12] + (1.0/2.0)*x34*M[17] + x35*x36 + x39*M[15] - x4*M[3] + x41*M[18] + (1.0/6.0)*x42*M[19] + x8*M[5] + x9*M[6] + M[0]/R;
#pragma omp atomic
L[1] += -x13*M[8] + x17*M[1] + x18*x28 - x2*M[0] - x25*M[7] + x30*M[5] + x33*M[6] + x39*M[9] + x8*M[2] + x9*M[3];
#pragma omp atomic
L[2] += x11*M[3] - x13*M[6] + x18*x30 + x21*M[2] + x22*x32 - x3*M[0] + x34*M[8] - x37*M[5] + x41*M[9] + x8*M[1];
#pragma omp atomic
L[3] += x11*M[2] - x13*M[5] + x18*x33 + x22*x34 + x35*M[3] + x36*x42 + x38*M[6] - x4*M[0] + x40*M[8] + x9*M[1];
#pragma omp atomic
L[4] += x17*M[0] + x28*M[1] + x30*M[2] + x33*M[3];
#pragma omp atomic
L[5] += -x13*M[3] + x30*M[1] - x37*M[2] + x8*M[0];
#pragma omp atomic
L[6] += -x13*M[2] + x33*M[1] + x38*M[3] + x9*M[0];
#pragma omp atomic
L[7] += x21*M[0] + x32*M[2] + x34*M[3] - x37*M[1];
#pragma omp atomic
L[8] += x11*M[0] - x13*M[1] + x34*M[2] + x40*M[3];
#pragma omp atomic
L[9] += x35*M[0] + x38*M[1] + x40*M[2] + x42*M[3];
#pragma omp atomic
L[10] += x28*M[0];
#pragma omp atomic
L[11] += x30*M[0];
#pragma omp atomic
L[12] += x33*M[0];
#pragma omp atomic
L[13] += -x37*M[0];
#pragma omp atomic
L[14] += -x13*M[0];
#pragma omp atomic
L[15] += x38*M[0];
#pragma omp atomic
L[16] += x32*M[0];
#pragma omp atomic
L[17] += x34*M[0];
#pragma omp atomic
L[18] += x40*M[0];
#pragma omp atomic
L[19] += x42*M[0];

}

void L2L_3(double x, double y, double z, double * L, double * Ls) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
x3 = z*L[14];
x4 = x3*y;
x5 = (1.0/2.0)*(x*x);
x6 = (1.0/2.0)*(y*y);
x7 = (1.0/2.0)*(z*z);
x8 = x*L[13];
x9 = x*L[15];
x10 = y*L[11];
x11 = z*L[12];
x12 = y*L[18];
x13 = z*L[17];
x14 = y*L[13];
x15 = y*L[14];
x16 = z*L[15];
x17 = z*L[18];
#pragma omp atomic
Ls[0] += (1.0/6.0)*(x*x*x)*L[10] + x*x0 + x*x1 + x*x4 + x*L[1] + x10*x5 + x11*x5 + x12*x7 + x13*x6 + x2*y + x5*L[4] + x6*x8 + x6*L[7] + x7*x9 + x7*L[9] + (1.0/6.0)*(y*y*y)*L[16] + y*L[2] + (1.0/6.0)*(z*z*z)*L[19] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*x10 + x*x11 + x*L[4] + x0 + x1 + x4 + x5*L[10] + x6*L[13] + x7*L[15] + L[1];
#pragma omp atomic
Ls[2] += x*x14 + x*x3 + x*L[5] + x13*y + x2 + x5*L[11] + x6*L[16] + x7*L[18] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*x15 + x*x16 + x*L[6] + x17*y + x5*L[12] + x6*L[17] + x7*L[19] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += x*L[10] + x10 + x11 + L[4];
#pragma omp atomic
Ls[5] += x*L[11] + x14 + x3 + L[5];
#pragma omp atomic
Ls[6] += x*L[12] + x15 + x16 + L[6];
#pragma omp atomic
Ls[7] += x13 + x8 + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += x*L[14] + x17 + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += x12 + x9 + z*L[19] + L[9];
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
double x0;
double x1;
double x2;
double x3;
x0 = x*y;
x1 = (1.0/2.0)*(x*x);
x2 = (1.0/2.0)*(y*y);
x3 = (1.0/2.0)*(z*z);
#pragma omp atomic
F[0] += (1.0/6.0)*(x*x*x)*L[10] + x*x2*L[13] + x*x3*L[15] + x*z*L[6] + x*L[1] + x0*z*L[14] + x0*L[5] + x1*y*L[11] + x1*z*L[12] + x1*L[4] + x2*z*L[17] + x2*L[7] + x3*y*L[18] + x3*L[9] + (1.0/6.0)*(y*y*y)*L[16] + y*z*L[8] + y*L[2] + (1.0/6.0)*(z*z*z)*L[19] + z*L[3] + L[0];

}

void M2P_3(double x, double y, double z, double * M, double * F) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*x6;
x8 = y*z;
x9 = 15.0*pow(x2, -3.5);
x10 = -x4;
x11 = x0*x6;
x12 = x1*x6;
x13 = 9.0*x5;
x14 = -x0*x9;
x15 = x*(x13 + x14);
x16 = x1*x9;
x17 = x*(x16 - x6);
x18 = x14 + x6;
x19 = x18*y;
x20 = -x16;
x21 = y*(x13 + x20);
x22 = x18*z;
x23 = z*(x20 + x6);
#pragma omp atomic
F[0] += -x*x4*M[1] - x*x8*x9*M[14] + (1.0/6.0)*x15*M[10] - 0.5*x17*M[13] + (1.0/2.0)*x19*M[11] + pow(x2, -0.5)*M[0] + (1.0/6.0)*x21*M[16] + (1.0/2.0)*x22*M[12] + (1.0/2.0)*x23*M[17] - x4*y*M[2] - x4*z*M[3] + x6*x8*M[8] + x7*y*M[5] + x7*z*M[6] + (1.0/2.0)*(x10 + x11)*M[4] + (1.0/2.0)*(x10 + x12)*M[7] + (1.0/2.0)*(-x15 + 1.0*x17)*M[15] + (1.0/2.0)*(-x19 - x21)*M[18] + (1.0/6.0)*(-x22 - x23)*M[19] + (1.0/2.0)*(-x11 - x12 + 2.0*x3)*M[9];

}

void P2M_4(double x, double y, double z, double q, double * M) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = q*x3;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = q*x7;
x9 = x1*z;
x10 = (z*z);
x11 = (x*x*x);
x12 = (y*y*y);
x13 = (z*z*z);
M[0] += q;
M[1] += -x0;
M[2] += -x1;
M[3] += -x2;
M[4] += x4;
M[5] += x5;
M[6] += x6;
M[7] += x8;
M[8] += x9;
M[9] += q*x10;
M[10] += -q*x11;
M[11] += -x1*x3;
M[12] += -x2*x3;
M[13] += -x0*x7;
M[14] += -x5*z;
M[15] += -x0*x10;
M[16] += -q*x12;
M[17] += -x2*x7;
M[18] += -x1*x10;
M[19] += -q*x13;
M[20] += q*(x*x*x*x);
M[21] += x1*x11;
M[22] += x11*x2;
M[23] += x4*x7;
M[24] += x3*x9;
M[25] += x10*x4;
M[26] += x0*x12;
M[27] += x6*x7;
M[28] += x10*x5;
M[29] += x0*x13;
M[30] += q*(y*y*y*y);
M[31] += x12*x2;
M[32] += x10*x8;
M[33] += x1*x13;
M[34] += q*(z*z*z*z);
}
void M2M_4(double x, double y, double z, double * M, double * Ms) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
x0 = x*M[0];
x1 = y*M[0];
x2 = z*M[0];
x3 = x*M[1];
x4 = (x*x);
x5 = x4*M[0];
x6 = x*M[2];
x7 = y*M[1];
x8 = x0*y;
x9 = x*M[3];
x10 = z*M[1];
x11 = x0*z;
x12 = y*M[2];
x13 = (y*y);
x14 = x13*M[0];
x15 = y*M[3];
x16 = z*M[2];
x17 = x1*z;
x18 = z*M[3];
x19 = (z*z);
x20 = x*M[4];
x21 = (x*x*x);
x22 = x*M[5];
x23 = y*M[4];
x24 = x3*y;
x25 = x*M[6];
x26 = x*M[7];
x27 = y*M[5];
x28 = x6*y;
x29 = x*M[8];
x30 = y*M[6];
x31 = x9*y;
x32 = x*M[9];
x33 = y*M[7];
x34 = (y*y*y);
x35 = y*M[8];
x36 = y*M[9];
x37 = (z*z*z);
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += x0 + M[1];
#pragma omp atomic
Ms[2] += x1 + M[2];
#pragma omp atomic
Ms[3] += x2 + M[3];
#pragma omp atomic
Ms[4] += x3 + x5 + M[4];
#pragma omp atomic
Ms[5] += x6 + x7 + x8 + M[5];
#pragma omp atomic
Ms[6] += x10 + x11 + x9 + M[6];
#pragma omp atomic
Ms[7] += x12 + x14 + M[7];
#pragma omp atomic
Ms[8] += x15 + x16 + x17 + M[8];
#pragma omp atomic
Ms[9] += x18 + x19*M[0] + M[9];
#pragma omp atomic
Ms[10] += x20 + x21*M[0] + x4*M[1] + M[10];
#pragma omp atomic
Ms[11] += x1*x4 + x22 + x23 + x24 + x4*M[2] + M[11];
#pragma omp atomic
Ms[12] += x2*x4 + x25 + x3*z + x4*M[3] + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += x0*x13 + x13*M[1] + x26 + x27 + x28 + M[13];
#pragma omp atomic
Ms[14] += x29 + x30 + x31 + x6*z + x7*z + x8*z + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += x0*x19 + x19*M[1] + x32 + x9*z + z*M[6] + M[15];
#pragma omp atomic
Ms[16] += x13*M[2] + x33 + x34*M[0] + M[16];
#pragma omp atomic
Ms[17] += x12*z + x13*x2 + x13*M[3] + x35 + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += x1*x19 + x15*z + x19*M[2] + x36 + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += x19*M[3] + x37*M[0] + z*M[9] + M[19];
#pragma omp atomic
Ms[20] += (x*x*x*x)*M[0] + x*M[10] + x21*M[1] + x4*M[4] + M[20];
#pragma omp atomic
Ms[21] += x*M[11] + x1*x21 + x20*y + x21*M[2] + x4*x7 + x4*M[5] + y*M[10] + M[21];
#pragma omp atomic
Ms[22] += x*M[12] + x10*x4 + x2*x21 + x20*z + x21*M[3] + x4*M[6] + z*M[10] + M[22];
#pragma omp atomic
Ms[23] += x*M[13] + x12*x4 + x13*x3 + x13*x5 + x13*M[4] + x22*y + x4*M[7] + y*M[11] + M[23];
#pragma omp atomic
Ms[24] += x*M[14] + x15*x4 + x16*x4 + x17*x4 + x22*z + x23*z + x24*z + x25*y + x4*M[8] + y*M[12] + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += x*M[15] + x18*x4 + x19*x3 + x19*x5 + x19*M[4] + x25*z + x4*M[9] + z*M[12] + M[25];
#pragma omp atomic
Ms[26] += x*M[16] + x0*x34 + x13*x6 + x13*M[5] + x26*y + x34*M[1] + y*M[13] + M[26];
#pragma omp atomic
Ms[27] += x*M[17] + x10*x13 + x11*x13 + x13*x9 + x13*M[6] + x26*z + x27*z + x28*z + x29*y + y*M[14] + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += x*M[18] + x19*x6 + x19*x7 + x19*x8 + x19*M[5] + x29*z + x30*z + x31*z + x32*y + y*M[15] + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += x*M[19] + x0*x37 + x19*x9 + x19*M[6] + x32*z + x37*M[1] + z*M[15] + M[29];
#pragma omp atomic
Ms[30] += x13*M[7] + x34*M[2] + (y*y*y*y)*M[0] + y*M[16] + M[30];
#pragma omp atomic
Ms[31] += x13*x16 + x13*M[8] + x2*x34 + x33*z + x34*M[3] + y*M[17] + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += x12*x19 + x13*x18 + x13*M[9] + x14*x19 + x19*M[7] + x35*z + y*M[18] + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += x1*x37 + x15*x19 + x19*M[8] + x36*z + x37*M[2] + y*M[19] + z*M[18] + M[33];
#pragma omp atomic
Ms[34] += x19*M[9] + x37*M[3] + (z*z*z*z)*M[0] + z*M[19] + M[34];

}

void M2L_4(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
double x49;
double x50;
double x51;
double x52;
double x53;
double x54;
double x55;
double x56;
double x57;
double x58;
double x59;
double x60;
double x61;
double x62;
double x63;
double x64;
double x65;
double x66;
double x67;
double x68;
double x69;
double x70;
double x71;
double x72;
double x73;
double x74;
double x75;
double x76;
double x77;
double x78;
double x79;
double x80;
double x81;
double x82;
double x83;
double x84;
double x85;
double x86;
double x87;
double x88;
double x89;
double x90;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*x1;
x3 = x1*y;
x4 = x1*z;
x5 = (1 / (R*R*R*R*R));
x6 = 3.0*x5;
x7 = x*x6;
x8 = x7*y;
x9 = x7*z;
x10 = y*z;
x11 = x10*x6;
x12 = pow(R, -7);
x13 = 15.0*x12;
x14 = x*x10*x13;
x15 = -x1;
x16 = (x*x);
x17 = x16*x6;
x18 = x15 + x17;
x19 = (1.0/2.0)*M[4];
x20 = (y*y);
x21 = x20*x6;
x22 = x15 + x21;
x23 = (1.0/2.0)*M[7];
x24 = x13*x20;
x25 = x24 - x6;
x26 = 0.5*x;
x27 = x25*x26;
x28 = 9.0*x5;
x29 = x13*x16;
x30 = -x29;
x31 = x*(x28 + x30);
x32 = (1.0/6.0)*M[10];
x33 = x30 + x6;
x34 = x33*y;
x35 = (1.0/2.0)*M[11];
x36 = -x24;
x37 = y*(x28 + x36);
x38 = (1.0/6.0)*M[16];
x39 = x33*z;
x40 = (1.0/2.0)*M[12];
x41 = z*(x36 + x6);
x42 = (1.0/2.0)*M[17];
x43 = -45.0*x12;
x44 = pow(R, -9);
x45 = 105.0*x44;
x46 = x16*x45;
x47 = x*(x43 + x46);
x48 = x47*y;
x49 = x20*x45;
x50 = x43 + x49;
x51 = x50*y;
x52 = 0.16666666666666666*x*x51;
x53 = x47*z;
x54 = -x13;
x55 = z*(x49 + x54);
x56 = x26*x55;
x57 = x10*x50;
x58 = x10*(x46 + x54);
x59 = (x*x*x*x)*x45;
x60 = 90.0*x12;
x61 = -x16*x60 + x28 + x59;
x62 = x45*(y*y*y*y);
x63 = -x20*x60 + x28 + x62;
x64 = 2.0*x0 - x17 - x21;
x65 = (1.0/2.0)*M[9];
x66 = x20*x46;
x67 = x33 + x36 + x66;
x68 = 1.0*x;
x69 = x25*x68;
x70 = -x31 + x69;
x71 = (1.0/2.0)*x70;
x72 = -x34 - x37;
x73 = (1.0/2.0)*M[18];
x74 = -x39 - x41;
x75 = (1.0/6.0)*M[19];
x76 = x51*x68;
x77 = -x48 - x76;
x78 = (1.0/2.0)*x77;
x79 = x55*x68;
x80 = -x53 - x79;
x81 = -x57 - x58;
x82 = 105.0*x12;
x83 = -12.0*x5 - x66;
x84 = x20*x82 + x29 - x62 + x83;
x85 = x16*x82 + x24 - x59 + x83;
x86 = 120.0*x12;
x87 = 210.0*x16*x20*x44 - x16*x86 - x20*x86 + 24.0*x5 + x59 + x62;
x88 = (1.0/2.0)*x67;
x89 = (1.0/2.0)*M[15];
x90 = x26*x51;
#pragma omp atomic
L[0] += x11*M[8] - x14*M[14] + x18*x19 - x2*M[1] + x22*x23 - x27*M[13] - x3*M[2] + x31*x32 + x34*x35 + x37*x38 + x39*x40 - x4*M[3] + x41*x42 + (1.0/6.0)*x48*M[21] + x52*M[26] + (1.0/6.0)*x53*M[22] + x56*M[27] + (1.0/6.0)*x57*M[31] + (1.0/2.0)*x58*M[24] + (1.0/24.0)*x61*M[20] + (1.0/24.0)*x63*M[30] + x64*x65 + (1.0/4.0)*x67*M[23] + x71*M[15] + x72*x73 + x74*x75 + x78*M[28] + x8*M[5] + (1.0/6.0)*x80*M[29] + (1.0/6.0)*x81*M[33] + (1.0/4.0)*x84*M[32] + (1.0/4.0)*x85*M[25] + (1.0/24.0)*x87*M[34] + x9*M[6] + M[0]/R;
#pragma omp atomic
L[1] += -x14*M[8] + x18*M[1] + x19*x31 - x2*M[0] - x27*M[7] + x32*x61 + x34*M[5] + x35*x48 + x39*M[6] + x40*x53 + x52*M[16] + x56*M[17] + x58*M[14] + x71*M[9] + x73*x77 + x75*x80 + x8*M[2] + x85*x89 + x88*M[13] + x9*M[3];
#pragma omp atomic
L[2] += x11*M[3] - x14*M[6] + x19*x34 + x22*M[2] + x23*x37 - x3*M[0] + x32*x48 + x38*x63 + x40*x58 + x41*M[8] + x42*x57 + x65*x72 - x69*M[5] + x73*x84 + x75*x81 + x78*M[15] + x79*M[14] + x8*M[1] + x88*M[11] + x90*M[13];
#pragma omp atomic
L[3] += x11*M[2] - x14*M[5] + x19*x39 + x23*x41 + x32*x53 + x35*x58 + x38*x57 - x4*M[0] + x40*x85 + x42*x84 + x56*M[13] + x64*M[3] + x65*x74 + x70*M[6] + x72*M[8] + x73*x81 + x75*x87 + x77*M[14] + x80*x89 + x9*M[1];
#pragma omp atomic
L[4] += x18*M[0] + x19*x61 + x23*x67 + x31*M[1] + x34*M[2] + x39*M[3] + x48*M[5] + x53*M[6] + x58*M[8] + x65*x85;
#pragma omp atomic
L[5] += -x14*M[3] + x19*x48 + x34*M[1] + x58*M[6] + x67*M[5] - x69*M[2] + x78*M[9] + x79*M[8] + x8*M[0] + x90*M[7];
#pragma omp atomic
L[6] += -x14*M[2] + x19*x53 + x39*M[1] + x56*M[7] + x58*M[5] + x65*x80 + x70*M[3] + x77*M[8] + x85*M[6] + x9*M[0];
#pragma omp atomic
L[7] += x19*x67 + x22*M[0] + x23*x63 + x37*M[2] + x41*M[3] + x57*M[8] + x65*x84 - x69*M[1] + x76*M[5] + x79*M[6];
#pragma omp atomic
L[8] += x11*M[0] - x14*M[1] + x19*x58 + x23*x57 + x41*M[2] + x65*x81 + x72*M[3] + x77*M[6] + x79*M[5] + x84*M[8];
#pragma omp atomic
L[9] += x19*x85 + x23*x84 + x64*M[0] + x65*x87 + x70*M[1] + x72*M[2] + x74*M[3] + x77*M[5] + x80*M[6] + x81*M[8];
#pragma omp atomic
L[10] += x31*M[0] + x48*M[2] + x53*M[3] + x61*M[1];
#pragma omp atomic
L[11] += x34*M[0] + x48*M[1] + x58*M[3] + x67*M[2];
#pragma omp atomic
L[12] += x39*M[0] + x53*M[1] + x58*M[2] + x85*M[3];
#pragma omp atomic
L[13] += x67*M[1] - x69*M[0] + x76*M[2] + x79*M[3];
#pragma omp atomic
L[14] += -x14*M[0] + x58*M[1] + x77*M[3] + x79*M[2];
#pragma omp atomic
L[15] += x70*M[0] + x77*M[2] + x80*M[3] + x85*M[1];
#pragma omp atomic
L[16] += x37*M[0] + x57*M[3] + x63*M[2] + x76*M[1];
#pragma omp atomic
L[17] += x41*M[0] + x57*M[2] + x79*M[1] + x84*M[3];
#pragma omp atomic
L[18] += x72*M[0] + x77*M[1] + x81*M[3] + x84*M[2];
#pragma omp atomic
L[19] += x74*M[0] + x80*M[1] + x81*M[2] + x87*M[3];
#pragma omp atomic
L[20] += x61*M[0];
#pragma omp atomic
L[21] += x48*M[0];
#pragma omp atomic
L[22] += x53*M[0];
#pragma omp atomic
L[23] += x67*M[0];
#pragma omp atomic
L[24] += x58*M[0];
#pragma omp atomic
L[25] += x85*M[0];
#pragma omp atomic
L[26] += x76*M[0];
#pragma omp atomic
L[27] += x79*M[0];
#pragma omp atomic
L[28] += x77*M[0];
#pragma omp atomic
L[29] += x80*M[0];
#pragma omp atomic
L[30] += x63*M[0];
#pragma omp atomic
L[31] += x57*M[0];
#pragma omp atomic
L[32] += x84*M[0];
#pragma omp atomic
L[33] += x81*M[0];
#pragma omp atomic
L[34] += x87*M[0];

}

void L2L_4(double x, double y, double z, double * L, double * Ls) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
double x49;
double x50;
double x51;
double x52;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
x3 = z*L[14];
x4 = x3*y;
x5 = (x*x);
x6 = (1.0/2.0)*x5;
x7 = (1.0/6.0)*(x*x*x);
x8 = (y*y);
x9 = (1.0/2.0)*x8;
x10 = (1.0/6.0)*(y*y*y);
x11 = (z*z);
x12 = (1.0/2.0)*x11;
x13 = (1.0/6.0)*(z*z*z);
x14 = x*L[13];
x15 = x*L[26];
x16 = x*L[15];
x17 = x*L[29];
x18 = y*L[11];
x19 = z*L[12];
x20 = y*L[21];
x21 = z*L[22];
x22 = y*L[18];
x23 = y*L[33];
x24 = z*L[17];
x25 = z*L[31];
x26 = y*L[28];
x27 = x*x26;
x28 = z*L[27];
x29 = x*x28;
x30 = z*L[24];
x31 = x30*y;
x32 = (1.0/4.0)*x5;
x33 = x*L[23];
x34 = x*L[25];
x35 = y*L[13];
x36 = x28*y;
x37 = x*L[28];
x38 = y*L[23];
x39 = y*L[32];
x40 = y*L[14];
x41 = z*L[15];
x42 = z*L[18];
x43 = z*L[28];
x44 = x43*y;
x45 = x*L[27];
x46 = y*L[24];
x47 = z*L[25];
x48 = z*L[32];
x49 = y*L[26];
x50 = y*L[27];
x51 = z*L[29];
x52 = z*L[33];
#pragma omp atomic
Ls[0] += (1.0/24.0)*(x*x*x*x)*L[20] + x*x0 + x*x1 + x*x4 + x*L[1] + x10*x15 + x10*x25 + x10*L[16] + x11*x32*L[25] + (1.0/4.0)*x11*x8*L[32] + x12*x16 + x12*x22 + x12*x27 + x12*L[9] + x13*x17 + x13*x23 + x13*L[19] + x14*x9 + x18*x6 + x19*x6 + x2*y + x20*x7 + x21*x7 + x24*x9 + x29*x9 + x31*x6 + x32*x8*L[23] + x6*L[4] + x7*L[10] + x9*L[7] + (1.0/24.0)*(y*y*y*y)*L[30] + y*L[2] + (1.0/24.0)*(z*z*z*z)*L[34] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*x18 + x*x19 + x*x31 + x*L[4] + x0 + x1 + x10*L[26] + x12*x26 + x12*x34 + x12*L[15] + x13*L[29] + x20*x6 + x21*x6 + x28*x9 + x33*x9 + x4 + x6*L[10] + x7*L[20] + x9*L[13] + L[1];
#pragma omp atomic
Ls[2] += x*x3 + x*x35 + x*x36 + x*L[5] + x10*L[30] + x12*x37 + x12*x39 + x12*L[18] + x13*L[33] + x15*x9 + x2 + x24*y + x25*x9 + x30*x6 + x38*x6 + x6*L[11] + x7*L[21] + x9*L[16] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*x40 + x*x41 + x*x44 + x*L[6] + x10*L[31] + x12*x17 + x12*x23 + x12*L[19] + x13*L[34] + x42*y + x45*x9 + x46*x6 + x47*x6 + x48*x9 + x6*L[12] + x7*L[22] + x9*L[17] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += x*x20 + x*x21 + x*L[10] + x12*L[25] + x18 + x19 + x31 + x6*L[20] + x9*L[23] + L[4];
#pragma omp atomic
Ls[5] += x*x30 + x*x38 + x*L[11] + x12*L[28] + x3 + x35 + x36 + x6*L[21] + x9*L[26] + L[5];
#pragma omp atomic
Ls[6] += x*x46 + x*x47 + x*L[12] + x12*L[29] + x40 + x41 + x44 + x6*L[22] + x9*L[27] + L[6];
#pragma omp atomic
Ls[7] += x*x49 + x12*L[32] + x14 + x24 + x25*y + x29 + x6*L[23] + x9*L[30] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += x*x43 + x*x50 + x*L[14] + x12*L[33] + x42 + x48*y + x6*L[24] + x9*L[31] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += x*x51 + x12*L[34] + x16 + x22 + x27 + x52*y + x6*L[25] + x9*L[32] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += x*L[20] + x20 + x21 + L[10];
#pragma omp atomic
Ls[11] += x*L[21] + x30 + x38 + L[11];
#pragma omp atomic
Ls[12] += x*L[22] + x46 + x47 + L[12];
#pragma omp atomic
Ls[13] += x28 + x33 + x49 + L[13];
#pragma omp atomic
Ls[14] += x*L[24] + x43 + x50 + L[14];
#pragma omp atomic
Ls[15] += x26 + x34 + x51 + L[15];
#pragma omp atomic
Ls[16] += x15 + x25 + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += x45 + x48 + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += x37 + x39 + x52 + L[18];
#pragma omp atomic
Ls[19] += x17 + x23 + z*L[34] + L[19];
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
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
x0 = x*y;
x1 = x*z;
x2 = y*z;
x3 = (x*x);
x4 = (1.0/2.0)*x3;
x5 = (1.0/6.0)*(x*x*x);
x6 = (y*y);
x7 = (1.0/2.0)*x6;
x8 = (1.0/6.0)*(y*y*y);
x9 = (z*z);
x10 = (1.0/2.0)*x9;
x11 = (1.0/6.0)*(z*z*z);
x12 = (1.0/4.0)*x3;
#pragma omp atomic
F[0] += (1.0/24.0)*(x*x*x*x)*L[20] + x*x10*L[15] + x*x11*L[29] + x*x7*L[13] + x*x8*L[26] + x*L[1] + x0*x10*L[28] + x0*z*L[14] + x0*L[5] + x1*x7*L[27] + x1*L[6] + x10*y*L[18] + x10*L[9] + x11*y*L[33] + x11*L[19] + x12*x6*L[23] + x12*x9*L[25] + x2*x4*L[24] + x2*L[8] + x4*y*L[11] + x4*z*L[12] + x4*L[4] + x5*y*L[21] + x5*z*L[22] + x5*L[10] + (1.0/4.0)*x6*x9*L[32] + x7*z*L[17] + x7*L[7] + x8*z*L[31] + x8*L[16] + (1.0/24.0)*(y*y*y*y)*L[30] + y*L[2] + (1.0/24.0)*(z*z*z*z)*L[34] + z*L[3] + L[0];

}

void M2P_4(double x, double y, double z, double * M, double * F) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*x6;
x8 = y*z;
x9 = pow(x2, -3.5);
x10 = 15.0*x9;
x11 = -x4;
x12 = x0*x6;
x13 = x1*x6;
x14 = 9.0*x5;
x15 = x0*x10;
x16 = -x15;
x17 = x*(x14 + x16);
x18 = x1*x10;
x19 = x18 - x6;
x20 = 0.5*x;
x21 = x16 + x6;
x22 = x21*y;
x23 = -x18;
x24 = y*(x14 + x23);
x25 = x21*z;
x26 = z*(x23 + x6);
x27 = -45.0*x9;
x28 = pow(x2, -4.5);
x29 = 105.0*x28;
x30 = x0*x29;
x31 = x*(x27 + x30);
x32 = x31*y;
x33 = x1*x29;
x34 = x27 + x33;
x35 = x34*y;
x36 = x31*z;
x37 = -x10;
x38 = z*(x33 + x37);
x39 = x8*(x30 + x37);
x40 = x34*x8;
x41 = 90.0*x9;
x42 = (x*x*x*x)*x29;
x43 = x29*(y*y*y*y);
x44 = x1*x30;
x45 = 1.0*x;
x46 = 105.0*x9;
x47 = -x44 - 12.0*x5;
x48 = 120.0*x9;
#pragma omp atomic
F[0] += -x*x10*x8*M[14] + 0.16666666666666666*x*x35*M[26] - x*x4*M[1] + (1.0/6.0)*x17*M[10] - x19*x20*M[13] + pow(x2, -0.5)*M[0] + x20*x38*M[27] + (1.0/2.0)*x22*M[11] + (1.0/6.0)*x24*M[16] + (1.0/2.0)*x25*M[12] + (1.0/2.0)*x26*M[17] + (1.0/6.0)*x32*M[21] + (1.0/6.0)*x36*M[22] + (1.0/2.0)*x39*M[24] - x4*y*M[2] - x4*z*M[3] + (1.0/6.0)*x40*M[31] + x6*x8*M[8] + x7*y*M[5] + x7*z*M[6] + (1.0/2.0)*(x11 + x12)*M[4] + (1.0/2.0)*(x11 + x13)*M[7] + (1.0/2.0)*(-x17 + x19*x45)*M[15] + (1.0/2.0)*(-x22 - x24)*M[18] + (1.0/6.0)*(-x25 - x26)*M[19] + (1.0/2.0)*(-x32 - x35*x45)*M[28] + (1.0/6.0)*(-x36 - x38*x45)*M[29] + (1.0/6.0)*(-x39 - x40)*M[33] + (1.0/2.0)*(-x12 - x13 + 2.0*x3)*M[9] + (1.0/4.0)*(x21 + x23 + x44)*M[23] + (1.0/24.0)*(-x0*x41 + x14 + x42)*M[20] + (1.0/24.0)*(-x1*x41 + x14 + x43)*M[30] + (1.0/4.0)*(x0*x46 + x18 - x42 + x47)*M[25] + (1.0/4.0)*(x1*x46 + x15 - x43 + x47)*M[32] + (1.0/24.0)*(210.0*x0*x1*x28 - x0*x48 - x1*x48 + x42 + x43 + 24.0*x5)*M[34];

}

void P2M_5(double x, double y, double z, double q, double * M) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = q*x3;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = q*x7;
x9 = x1*z;
x10 = (z*z);
x11 = q*x10;
x12 = (x*x*x);
x13 = x1*x3;
x14 = x2*x3;
x15 = x0*x7;
x16 = (y*y*y);
x17 = (z*z*z);
x18 = (x*x*x*x);
x19 = (y*y*y*y);
x20 = (z*z*z*z);
M[0] += q;
M[1] += -x0;
M[2] += -x1;
M[3] += -x2;
M[4] += x4;
M[5] += x5;
M[6] += x6;
M[7] += x8;
M[8] += x9;
M[9] += x11;
M[10] += -q*x12;
M[11] += -x13;
M[12] += -x14;
M[13] += -x15;
M[14] += -x5*z;
M[15] += -x0*x10;
M[16] += -q*x16;
M[17] += -x2*x7;
M[18] += -x1*x10;
M[19] += -q*x17;
M[20] += q*x18;
M[21] += x1*x12;
M[22] += x12*x2;
M[23] += x4*x7;
M[24] += x3*x9;
M[25] += x10*x4;
M[26] += x0*x16;
M[27] += x6*x7;
M[28] += x10*x5;
M[29] += x0*x17;
M[30] += q*x19;
M[31] += x16*x2;
M[32] += x10*x8;
M[33] += x1*x17;
M[34] += q*x20;
M[35] += -q*(x*x*x*x*x);
M[36] += -x1*x18;
M[37] += -x18*x2;
M[38] += -x12*x8;
M[39] += -x12*x9;
M[40] += -x11*x12;
M[41] += -x16*x4;
M[42] += -x14*x7;
M[43] += -x10*x13;
M[44] += -x17*x4;
M[45] += -x0*x19;
M[46] += -x16*x6;
M[47] += -x10*x15;
M[48] += -x17*x5;
M[49] += -x0*x20;
M[50] += -q*(y*y*y*y*y);
M[51] += -x19*x2;
M[52] += -x11*x16;
M[53] += -x17*x8;
M[54] += -x1*x20;
M[55] += -q*(z*z*z*z*z);
}
void M2M_5(double x, double y, double z, double * M, double * Ms) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
double x49;
double x50;
double x51;
double x52;
double x53;
double x54;
double x55;
double x56;
double x57;
double x58;
double x59;
double x60;
double x61;
double x62;
double x63;
double x64;
double x65;
double x66;
double x67;
double x68;
double x69;
double x70;
double x71;
double x72;
double x73;
double x74;
double x75;
double x76;
double x77;
double x78;
double x79;
double x80;
double x81;
double x82;
double x83;
double x84;
double x85;
double x86;
double x87;
double x88;
x0 = x*M[0];
x1 = y*M[0];
x2 = z*M[0];
x3 = x*M[1];
x4 = (x*x);
x5 = x4*M[0];
x6 = x*M[2];
x7 = y*M[1];
x8 = x0*y;
x9 = x*M[3];
x10 = z*M[1];
x11 = x0*z;
x12 = y*M[2];
x13 = (y*y);
x14 = x13*M[0];
x15 = y*M[3];
x16 = z*M[2];
x17 = x1*z;
x18 = z*M[3];
x19 = (z*z);
x20 = x19*M[0];
x21 = x*M[4];
x22 = x4*M[1];
x23 = (x*x*x);
x24 = x*M[5];
x25 = y*M[4];
x26 = x3*y;
x27 = x4*M[2];
x28 = x1*x4;
x29 = x*M[6];
x30 = z*M[4];
x31 = x3*z;
x32 = x4*M[3];
x33 = x2*x4;
x34 = x*M[7];
x35 = y*M[5];
x36 = x6*y;
x37 = x13*M[1];
x38 = x0*x13;
x39 = x*M[8];
x40 = y*M[6];
x41 = z*M[5];
x42 = x9*y;
x43 = x6*z;
x44 = x7*z;
x45 = x*M[9];
x46 = z*M[6];
x47 = x9*z;
x48 = y*M[7];
x49 = x13*M[2];
x50 = (y*y*y);
x51 = y*M[8];
x52 = z*M[7];
x53 = x12*z;
x54 = x13*M[3];
x55 = y*M[9];
x56 = z*M[8];
x57 = x15*z;
x58 = z*M[9];
x59 = (z*z*z);
x60 = x*M[10];
x61 = (x*x*x*x);
x62 = x*M[11];
x63 = y*M[10];
x64 = x21*y;
x65 = x*M[12];
x66 = x*M[13];
x67 = y*M[11];
x68 = x24*y;
x69 = x*M[14];
x70 = y*M[12];
x71 = x29*y;
x72 = x*M[15];
x73 = x*M[16];
x74 = y*M[13];
x75 = x34*y;
x76 = x*M[17];
x77 = y*M[14];
x78 = x39*y;
x79 = x*M[18];
x80 = y*M[15];
x81 = x45*y;
x82 = x*M[19];
x83 = y*M[16];
x84 = (y*y*y*y);
x85 = y*M[17];
x86 = y*M[18];
x87 = y*M[19];
x88 = (z*z*z*z);
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += x0 + M[1];
#pragma omp atomic
Ms[2] += x1 + M[2];
#pragma omp atomic
Ms[3] += x2 + M[3];
#pragma omp atomic
Ms[4] += x3 + x5 + M[4];
#pragma omp atomic
Ms[5] += x6 + x7 + x8 + M[5];
#pragma omp atomic
Ms[6] += x10 + x11 + x9 + M[6];
#pragma omp atomic
Ms[7] += x12 + x14 + M[7];
#pragma omp atomic
Ms[8] += x15 + x16 + x17 + M[8];
#pragma omp atomic
Ms[9] += x18 + x20 + M[9];
#pragma omp atomic
Ms[10] += x21 + x22 + x23*M[0] + M[10];
#pragma omp atomic
Ms[11] += x24 + x25 + x26 + x27 + x28 + M[11];
#pragma omp atomic
Ms[12] += x29 + x30 + x31 + x32 + x33 + M[12];
#pragma omp atomic
Ms[13] += x34 + x35 + x36 + x37 + x38 + M[13];
#pragma omp atomic
Ms[14] += x39 + x40 + x41 + x42 + x43 + x44 + x8*z + M[14];
#pragma omp atomic
Ms[15] += x0*x19 + x19*M[1] + x45 + x46 + x47 + M[15];
#pragma omp atomic
Ms[16] += x48 + x49 + x50*M[0] + M[16];
#pragma omp atomic
Ms[17] += x13*x2 + x51 + x52 + x53 + x54 + M[17];
#pragma omp atomic
Ms[18] += x1*x19 + x19*M[2] + x55 + x56 + x57 + M[18];
#pragma omp atomic
Ms[19] += x19*M[3] + x58 + x59*M[0] + M[19];
#pragma omp atomic
Ms[20] += x23*M[1] + x4*M[4] + x60 + x61*M[0] + M[20];
#pragma omp atomic
Ms[21] += x1*x23 + x23*M[2] + x4*x7 + x4*M[5] + x62 + x63 + x64 + M[21];
#pragma omp atomic
Ms[22] += x10*x4 + x2*x23 + x21*z + x23*M[3] + x4*M[6] + x65 + z*M[10] + M[22];
#pragma omp atomic
Ms[23] += x12*x4 + x13*x3 + x13*x5 + x13*M[4] + x4*M[7] + x66 + x67 + x68 + M[23];
#pragma omp atomic
Ms[24] += x15*x4 + x16*x4 + x17*x4 + x24*z + x25*z + x26*z + x4*M[8] + x69 + x70 + x71 + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += x18*x4 + x19*x3 + x19*x5 + x19*M[4] + x29*z + x4*M[9] + x72 + z*M[12] + M[25];
#pragma omp atomic
Ms[26] += x0*x50 + x13*x6 + x13*M[5] + x50*M[1] + x73 + x74 + x75 + M[26];
#pragma omp atomic
Ms[27] += x10*x13 + x11*x13 + x13*x9 + x13*M[6] + x34*z + x35*z + x36*z + x76 + x77 + x78 + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += x19*x6 + x19*x7 + x19*x8 + x19*M[5] + x39*z + x40*z + x42*z + x79 + x80 + x81 + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += x0*x59 + x19*x9 + x19*M[6] + x45*z + x59*M[1] + x82 + z*M[15] + M[29];
#pragma omp atomic
Ms[30] += x13*M[7] + x50*M[2] + x83 + x84*M[0] + M[30];
#pragma omp atomic
Ms[31] += x13*x16 + x13*M[8] + x2*x50 + x48*z + x50*M[3] + x85 + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += x12*x19 + x13*x18 + x13*M[9] + x14*x19 + x19*M[7] + x51*z + x86 + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += x1*x59 + x15*x19 + x19*M[8] + x55*z + x59*M[2] + x87 + z*M[18] + M[33];
#pragma omp atomic
Ms[34] += x19*M[9] + x59*M[3] + x88*M[0] + z*M[19] + M[34];
#pragma omp atomic
Ms[35] += (x*x*x*x*x)*M[0] + x*M[20] + x23*M[4] + x4*M[10] + x61*M[1] + M[35];
#pragma omp atomic
Ms[36] += x*M[21] + x1*x61 + x23*x7 + x23*M[5] + x25*x4 + x4*M[11] + x60*y + x61*M[2] + y*M[20] + M[36];
#pragma omp atomic
Ms[37] += x*M[22] + x10*x23 + x2*x61 + x23*M[6] + x30*x4 + x4*M[12] + x60*z + x61*M[3] + z*M[20] + M[37];
#pragma omp atomic
Ms[38] += x*M[23] + x12*x23 + x13*x21 + x13*x22 + x13*M[10] + x14*x23 + x23*M[7] + x35*x4 + x4*M[13] + x62*y + y*M[21] + M[38];
#pragma omp atomic
Ms[39] += x*M[24] + x15*x23 + x16*x23 + x17*x23 + x23*M[8] + x4*x40 + x4*x41 + x4*x44 + x4*M[14] + x62*z + x63*z + x64*z + x65*y + y*M[22] + z*M[21] + M[39];
#pragma omp atomic
Ms[40] += x*M[25] + x18*x23 + x19*x21 + x19*x22 + x19*M[10] + x20*x23 + x23*M[9] + x4*x46 + x4*M[15] + x65*z + z*M[22] + M[40];
#pragma omp atomic
Ms[41] += x*M[26] + x13*x24 + x13*x27 + x13*M[11] + x3*x50 + x4*x48 + x4*M[16] + x5*x50 + x50*M[4] + x66*y + y*M[23] + M[41];
#pragma omp atomic
Ms[42] += x*M[27] + x13*x29 + x13*x30 + x13*x31 + x13*x32 + x13*x33 + x13*M[12] + x4*x51 + x4*x52 + x4*x53 + x4*M[17] + x66*z + x67*z + x68*z + x69*y + y*M[24] + z*M[23] + M[42];
#pragma omp atomic
Ms[43] += x*M[28] + x19*x24 + x19*x25 + x19*x26 + x19*x27 + x19*x28 + x19*M[11] + x4*x55 + x4*x56 + x4*x57 + x4*M[18] + x69*z + x70*z + x71*z + x72*y + y*M[25] + z*M[24] + M[43];
#pragma omp atomic
Ms[44] += x*M[29] + x19*x29 + x19*x32 + x19*M[12] + x3*x59 + x4*x58 + x4*M[19] + x5*x59 + x59*M[4] + x72*z + z*M[25] + M[44];
#pragma omp atomic
Ms[45] += x*M[30] + x0*x84 + x13*x34 + x13*M[13] + x50*x6 + x50*M[5] + x73*y + x84*M[1] + y*M[26] + M[45];
#pragma omp atomic
Ms[46] += x*M[31] + x10*x50 + x11*x50 + x13*x39 + x13*x41 + x13*x43 + x13*M[14] + x50*x9 + x50*M[6] + x73*z + x74*z + x75*z + x76*y + y*M[27] + z*M[26] + M[46];
#pragma omp atomic
Ms[47] += x*M[32] + x13*x45 + x13*x46 + x13*x47 + x13*M[15] + x19*x34 + x19*x35 + x19*x36 + x19*x37 + x19*x38 + x19*M[13] + x76*z + x77*z + x78*z + x79*y + y*M[28] + z*M[27] + M[47];
#pragma omp atomic
Ms[48] += x*M[33] + x19*x39 + x19*x40 + x19*x42 + x19*M[14] + x59*x6 + x59*x7 + x59*x8 + x59*M[5] + x79*z + x80*z + x81*z + x82*y + y*M[29] + z*M[28] + M[48];
#pragma omp atomic
Ms[49] += x*M[34] + x0*x88 + x19*x45 + x19*M[15] + x59*x9 + x59*M[6] + x82*z + x88*M[1] + z*M[29] + M[49];
#pragma omp atomic
Ms[50] += x13*M[16] + x50*M[7] + x84*M[2] + (y*y*y*y*y)*M[0] + y*M[30] + M[50];
#pragma omp atomic
Ms[51] += x13*x52 + x13*M[17] + x16*x50 + x2*x84 + x50*M[8] + x83*z + x84*M[3] + y*M[31] + z*M[30] + M[51];
#pragma omp atomic
Ms[52] += x13*x56 + x13*M[18] + x18*x50 + x19*x48 + x19*x49 + x19*M[16] + x20*x50 + x50*M[9] + x85*z + y*M[32] + z*M[31] + M[52];
#pragma omp atomic
Ms[53] += x12*x59 + x13*x58 + x13*M[19] + x14*x59 + x19*x51 + x19*x54 + x19*M[17] + x59*M[7] + x86*z + y*M[33] + z*M[32] + M[53];
#pragma omp atomic
Ms[54] += x1*x88 + x15*x59 + x19*x55 + x19*M[18] + x59*M[8] + x87*z + x88*M[2] + y*M[34] + z*M[33] + M[54];
#pragma omp atomic
Ms[55] += x19*M[19] + x59*M[9] + x88*M[3] + (z*z*z*z*z)*M[0] + z*M[34] + M[55];

}

void M2L_5(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
double x49;
double x50;
double x51;
double x52;
double x53;
double x54;
double x55;
double x56;
double x57;
double x58;
double x59;
double x60;
double x61;
double x62;
double x63;
double x64;
double x65;
double x66;
double x67;
double x68;
double x69;
double x70;
double x71;
double x72;
double x73;
double x74;
double x75;
double x76;
double x77;
double x78;
double x79;
double x80;
double x81;
double x82;
double x83;
double x84;
double x85;
double x86;
double x87;
double x88;
double x89;
double x90;
double x91;
double x92;
double x93;
double x94;
double x95;
double x96;
double x97;
double x98;
double x99;
double x100;
double x101;
double x102;
double x103;
double x104;
double x105;
double x106;
double x107;
double x108;
double x109;
double x110;
double x111;
double x112;
double x113;
double x114;
double x115;
double x116;
double x117;
double x118;
double x119;
double x120;
double x121;
double x122;
double x123;
double x124;
double x125;
double x126;
double x127;
double x128;
double x129;
double x130;
double x131;
double x132;
double x133;
double x134;
double x135;
double x136;
double x137;
double x138;
double x139;
double x140;
double x141;
double x142;
double x143;
double x144;
double x145;
double x146;
double x147;
double x148;
double x149;
double x150;
double x151;
double x152;
double x153;
double x154;
double x155;
double x156;
double x157;
double x158;
double x159;
double x160;
double x161;
double x162;
double x163;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*x1;
x3 = x1*y;
x4 = x1*z;
x5 = (1 / (R*R*R*R*R));
x6 = 3.0*x5;
x7 = x*x6;
x8 = x7*y;
x9 = x7*z;
x10 = y*z;
x11 = x10*x6;
x12 = pow(R, -7);
x13 = 15.0*x12;
x14 = x*x10;
x15 = x13*x14;
x16 = -x1;
x17 = (x*x);
x18 = x17*x6;
x19 = x16 + x18;
x20 = (1.0/2.0)*M[4];
x21 = (y*y);
x22 = x21*x6;
x23 = x16 + x22;
x24 = (1.0/2.0)*M[7];
x25 = x13*x21;
x26 = x25 - x6;
x27 = 0.5*x;
x28 = x26*x27;
x29 = 9.0*x5;
x30 = x13*x17;
x31 = -x30;
x32 = x*(x29 + x31);
x33 = (1.0/6.0)*M[10];
x34 = x31 + x6;
x35 = x34*y;
x36 = (1.0/2.0)*M[11];
x37 = -x25;
x38 = y*(x29 + x37);
x39 = (1.0/6.0)*M[16];
x40 = x34*z;
x41 = (1.0/2.0)*M[12];
x42 = z*(x37 + x6);
x43 = (1.0/2.0)*M[17];
x44 = 45.0*x12;
x45 = -x44;
x46 = pow(R, -9);
x47 = 105.0*x46;
x48 = x17*x47;
x49 = x45 + x48;
x50 = x*x49;
x51 = x50*y;
x52 = (1.0/6.0)*M[21];
x53 = 0.16666666666666666*M[26];
x54 = x21*x47;
x55 = x45 + x54;
x56 = x55*y;
x57 = x*x56;
x58 = x50*z;
x59 = (1.0/6.0)*M[22];
x60 = -x13;
x61 = x54 + x60;
x62 = x61*z;
x63 = x27*x62;
x64 = x10*x55;
x65 = (1.0/6.0)*M[31];
x66 = x10*(x48 + x60);
x67 = (1.0/2.0)*M[24];
x68 = 315.0*x46;
x69 = 945.0/pow(R, 11);
x70 = x21*x69 - x68;
x71 = x14*x70;
x72 = 0.16666666666666666*x71;
x73 = x17*x69;
x74 = x14*(x68 - x73);
x75 = (x*x*x*x);
x76 = x47*x75;
x77 = 90.0*x12;
x78 = -x17*x77 + x29 + x76;
x79 = (1.0/24.0)*M[20];
x80 = (y*y*y*y);
x81 = x47*x80;
x82 = -x21*x77 + x29 + x81;
x83 = (1.0/24.0)*M[30];
x84 = 2.0*x0 - x18 - x22;
x85 = (1.0/2.0)*M[9];
x86 = -225.0*x12;
x87 = -x69*x75;
x88 = 1050.0*x46;
x89 = x*(x17*x88 + x86 + x87);
x90 = x69*x80;
x91 = 630.0*x46;
x92 = x21*x91;
x93 = x44 + x90 - x92;
x94 = x*x93;
x95 = 0.041666666666666664*x94;
x96 = -x90;
x97 = y*(x21*x88 + x86 + x96);
x98 = x17*x91 + x45 + x87;
x99 = x98*y;
x100 = x98*z;
x101 = z*(x45 + x92 + x96);
x102 = x21*x48;
x103 = x102 + x34 + x37;
x104 = (1.0/4.0)*M[23];
x105 = -x21*x73;
x106 = x*(x105 + x21*x68 + x49);
x107 = y*(x105 + x17*x68 + x55);
x108 = z*(x105 + x48 + x61);
x109 = 1.0*x;
x110 = x109*x26;
x111 = x110 - x32;
x112 = (1.0/2.0)*x111;
x113 = -x35 - x38;
x114 = (1.0/2.0)*M[18];
x115 = -x40 - x42;
x116 = (1.0/6.0)*M[19];
x117 = x109*x56;
x118 = -x117 - x51;
x119 = (1.0/2.0)*M[28];
x120 = x109*x62;
x121 = -x120 - x58;
x122 = (1.0/6.0)*M[29];
x123 = -x64 - x66;
x124 = (1.0/6.0)*M[33];
x125 = x10*x109*x70;
x126 = x125 - x74;
x127 = 105.0*x12;
x128 = -x102 - 12.0*x5;
x129 = x127*x21 + x128 + x30 - x81;
x130 = (1.0/4.0)*M[32];
x131 = x127*x17 + x128 + x25 - x76;
x132 = (1.0/4.0)*M[25];
x133 = 120.0*x12;
x134 = -x133*x17 - x133*x21 + 210.0*x17*x21*x46 + 24.0*x5 + x76 + x81;
x135 = (1.0/24.0)*M[34];
x136 = -x106;
x137 = x136 - x89;
x138 = x109*x93;
x139 = x136 + x138;
x140 = (1.0/4.0)*x139;
x141 = -x107;
x142 = x141 - x97;
x143 = x141 - x99;
x144 = (1.0/4.0)*x143;
x145 = -x108;
x146 = -x100 + x145;
x147 = -x101 + x145;
x148 = 2*x107 + x97 + x99;
x149 = (1.0/24.0)*x148;
x150 = x100 + x101 + 2*x108;
x151 = 2*x106 - x138 + x89;
x152 = (1.0/24.0)*x151;
x153 = 0.16666666666666666*M[16];
x154 = (1.0/2.0)*x103;
x155 = (1.0/2.0)*x108;
x156 = (1.0/2.0)*M[15];
x157 = x27*x56;
x158 = 0.5*x71;
x159 = 0.16666666666666666*x94;
x160 = (1.0/2.0)*x139;
x161 = (1.0/2.0)*x143;
x162 = (1.0/2.0)*M[13];
x163 = x27*x93;
#pragma omp atomic
L[0] += (1.0/24.0)*x100*M[37] + (1.0/24.0)*x101*M[51] + x103*x104 + (1.0/12.0)*x106*M[38] + (1.0/12.0)*x107*M[41] + (1.0/4.0)*x108*M[42] + x11*M[8] + x112*M[15] + x113*x114 + x115*x116 + x118*x119 + x121*x122 + x123*x124 + (1.0/6.0)*x126*M[48] + x129*x130 + x131*x132 + x134*x135 + (1.0/12.0)*x137*M[40] + x140*M[47] + (1.0/12.0)*x142*M[52] + x144*M[43] + (1.0/12.0)*x146*M[44] + (1.0/12.0)*x147*M[53] + x149*M[54] - x15*M[14] + (1.0/120.0)*x150*M[55] + x152*M[49] + x19*x20 - x2*M[1] + x23*x24 - x28*M[13] - x3*M[2] + x32*x33 + x35*x36 + x38*x39 - x4*M[3] + x40*x41 + x42*x43 + x51*x52 + x53*x57 + x58*x59 + x63*M[27] + x64*x65 + x66*x67 - x72*M[46] + (1.0/6.0)*x74*M[39] + x78*x79 + x8*M[5] + x82*x83 + x84*x85 + (1.0/120.0)*x89*M[35] + x9*M[6] - x95*M[45] + (1.0/120.0)*x97*M[50] + (1.0/24.0)*x99*M[36] + M[0]/R;
#pragma omp atomic
L[1] += x100*x59 + x104*x106 + (1.0/6.0)*x107*M[26] + x112*M[9] + x114*x118 + x116*x121 + x119*x143 + x122*x146 + x124*x126 + x131*x156 + x132*x137 + x140*M[32] - x15*M[8] + x152*M[34] + x153*x57 + x154*M[13] + x155*M[27] + x19*M[1] - x2*M[0] + x20*x32 - x28*M[7] + x33*x78 + x35*M[5] + x36*x51 + x40*M[6] + x41*x58 + x52*x99 + x63*M[17] + x66*M[14] + x67*x74 - x72*M[31] + x79*x89 + x8*M[2] + x9*M[3] - x95*M[30];
#pragma omp atomic
L[2] += x101*x65 + x104*x107 + x106*x52 + x11*M[3] - x110*M[5] + x113*x85 + x114*x129 + x116*x123 + x118*x156 + x119*x139 + x120*M[14] + x122*x126 + x124*x147 + x130*x142 + x144*M[25] + x149*M[34] - x15*M[6] + x154*M[11] + x155*M[24] + x157*M[13] - x158*M[27] - x159*M[26] + x20*x35 + x23*M[2] + x24*x38 - x3*M[0] + x33*x51 + x39*x82 + x41*x66 + x42*M[8] + x43*x64 + x59*x74 + x79*x99 + x8*M[1] + x83*x97;
#pragma omp atomic
L[3] += x100*x79 + x101*x83 + x104*x108 + x11*M[2] + x111*M[6] + x113*M[8] + x114*x123 + x115*x85 + x116*x134 + x118*M[14] + x119*x126 + x121*x156 + x122*x151 + x124*x148 + x129*x43 + x130*x147 + x131*x41 + x132*x146 + x135*x150 + x137*x59 + x142*x65 - x15*M[5] + x160*M[27] + x161*M[24] + x20*x40 + x24*x42 + x33*x58 + x36*x66 + x39*x64 - x4*M[0] + x52*x74 - x53*x71 + x63*M[13] + x84*M[3] + x9*M[1];
#pragma omp atomic
L[4] += x100*x41 + x103*x24 + x106*x162 + x107*x39 + x108*x43 + x114*x143 + x116*x146 + x131*x85 + x137*x156 + x19*M[0] + x20*x78 + x32*M[1] + x33*x89 + x35*M[2] + x36*x99 + x40*M[3] + x51*M[5] + x58*M[6] + x66*M[8] + x74*M[14];
#pragma omp atomic
L[5] += x103*M[5] + x106*x36 + x107*x162 + x108*M[14] - x110*M[2] + x114*x139 + x116*x126 + x118*x85 + x120*M[8] + x143*x156 - x15*M[3] + x157*M[7] - x158*M[17] - x159*M[16] + x20*x51 + x33*x99 + x35*M[1] + x41*x74 + x66*M[6] + x8*M[0];
#pragma omp atomic
L[6] += x100*x33 + x111*M[3] + x114*x126 + x116*x151 + x118*M[8] + x121*x85 + x131*M[6] + x137*x41 + x143*M[14] + x146*x156 - x15*M[2] - x153*x71 + x155*M[13] + x160*M[17] + x20*x58 + x36*x74 + x40*M[1] + x63*M[7] + x66*M[5] + x9*M[0];
#pragma omp atomic
L[7] += x101*x43 + x103*x20 + x106*x33 + x107*x36 + x108*x41 - x110*M[1] + x114*x142 + x116*x147 + x117*M[5] + x120*M[6] - x125*M[14] + x129*x85 + x139*x156 - x163*M[13] + x23*M[0] + x24*x82 + x38*M[2] + x39*x97 + x42*M[3] + x64*M[8];
#pragma omp atomic
L[8] += x101*x39 + x108*x36 + x11*M[0] + x113*M[3] + x114*x147 + x116*x148 + x118*M[6] + x120*M[5] + x123*x85 + x126*x156 + x129*M[8] + x139*M[14] + x142*x43 - x15*M[1] - x158*M[13] + x161*M[12] + x20*x66 + x24*x64 + x33*x74 + x42*M[2];
#pragma omp atomic
L[9] += x111*M[1] + x113*M[2] + x114*x148 + x115*M[3] + x116*x150 + x118*M[5] + x121*M[6] + x123*M[8] + x126*M[14] + x129*x24 + x131*x20 + x134*x85 + x137*x33 + x142*x39 + x146*x41 + x147*x43 + x151*x156 + x160*M[13] + x161*M[11] + x84*M[0];
#pragma omp atomic
L[10] += x100*M[6] + x106*x24 + x137*x85 + x20*x89 + x32*M[0] + x51*M[2] + x58*M[3] + x74*M[8] + x78*M[1] + x99*M[5];
#pragma omp atomic
L[11] += x103*M[2] + x106*M[5] + x107*x24 + x108*M[8] + x143*x85 + x20*x99 + x35*M[0] + x51*M[1] + x66*M[3] + x74*M[6];
#pragma omp atomic
L[12] += x100*x20 + x108*x24 + x131*M[3] + x137*M[6] + x143*M[8] + x146*x85 + x40*M[0] + x58*M[1] + x66*M[2] + x74*M[5];
#pragma omp atomic
L[13] += x103*M[1] + x106*x20 + x107*M[5] + x108*M[6] - x110*M[0] + x117*M[2] + x120*M[3] - x125*M[8] + x139*x85 - x163*M[7];
#pragma omp atomic
L[14] += x108*M[5] + x118*M[3] + x120*M[2] + x126*x85 + x139*M[8] + x143*M[6] - x15*M[0] - x158*M[7] + x20*x74 + x66*M[1];
#pragma omp atomic
L[15] += x111*M[0] + x118*M[2] + x121*M[3] + x126*M[8] + x131*M[1] + x137*x20 + x139*x24 + x143*M[5] + x146*M[6] + x151*x85;
#pragma omp atomic
L[16] += x101*M[8] + x107*x20 + x117*M[1] - x125*M[6] - x138*M[5] + x142*x85 + x24*x97 + x38*M[0] + x64*M[3] + x82*M[2];
#pragma omp atomic
L[17] += x101*x24 + x108*x20 + x120*M[1] - x125*M[5] + x129*M[3] + x139*M[6] + x142*M[8] + x147*x85 + x42*M[0] + x64*M[2];
#pragma omp atomic
L[18] += x113*M[0] + x118*M[1] + x123*M[3] + x126*M[6] + x129*M[2] + x139*M[5] + x142*x24 + x143*x20 + x147*M[8] + x148*x85;
#pragma omp atomic
L[19] += x115*M[0] + x121*M[1] + x123*M[2] + x126*M[5] + x134*M[3] + x146*x20 + x147*x24 + x148*M[8] + x150*x85 + x151*M[6];
#pragma omp atomic
L[20] += x100*M[3] + x78*M[0] + x89*M[1] + x99*M[2];
#pragma omp atomic
L[21] += x106*M[2] + x51*M[0] + x74*M[3] + x99*M[1];
#pragma omp atomic
L[22] += x100*M[1] + x137*M[3] + x58*M[0] + x74*M[2];
#pragma omp atomic
L[23] += x103*M[0] + x106*M[1] + x107*M[2] + x108*M[3];
#pragma omp atomic
L[24] += x108*M[2] + x143*M[3] + x66*M[0] + x74*M[1];
#pragma omp atomic
L[25] += x131*M[0] + x137*M[1] + x143*M[2] + x146*M[3];
#pragma omp atomic
L[26] += x107*M[1] + x117*M[0] - x125*M[3] - x138*M[2];
#pragma omp atomic
L[27] += x108*M[1] + x120*M[0] - x125*M[2] + x139*M[3];
#pragma omp atomic
L[28] += x118*M[0] + x126*M[3] + x139*M[2] + x143*M[1];
#pragma omp atomic
L[29] += x121*M[0] + x126*M[2] + x146*M[1] + x151*M[3];
#pragma omp atomic
L[30] += x101*M[3] - x138*M[1] + x82*M[0] + x97*M[2];
#pragma omp atomic
L[31] += x101*M[2] - x125*M[1] + x142*M[3] + x64*M[0];
#pragma omp atomic
L[32] += x129*M[0] + x139*M[1] + x142*M[2] + x147*M[3];
#pragma omp atomic
L[33] += x123*M[0] + x126*M[1] + x147*M[2] + x148*M[3];
#pragma omp atomic
L[34] += x134*M[0] + x148*M[2] + x150*M[3] + x151*M[1];
#pragma omp atomic
L[35] += x89*M[0];
#pragma omp atomic
L[36] += x99*M[0];
#pragma omp atomic
L[37] += x100*M[0];
#pragma omp atomic
L[38] += x106*M[0];
#pragma omp atomic
L[39] += x74*M[0];
#pragma omp atomic
L[40] += x137*M[0];
#pragma omp atomic
L[41] += x107*M[0];
#pragma omp atomic
L[42] += x108*M[0];
#pragma omp atomic
L[43] += x143*M[0];
#pragma omp atomic
L[44] += x146*M[0];
#pragma omp atomic
L[45] += -x138*M[0];
#pragma omp atomic
L[46] += -x125*M[0];
#pragma omp atomic
L[47] += x139*M[0];
#pragma omp atomic
L[48] += x126*M[0];
#pragma omp atomic
L[49] += x151*M[0];
#pragma omp atomic
L[50] += x97*M[0];
#pragma omp atomic
L[51] += x101*M[0];
#pragma omp atomic
L[52] += x142*M[0];
#pragma omp atomic
L[53] += x147*M[0];
#pragma omp atomic
L[54] += x148*M[0];
#pragma omp atomic
L[55] += x150*M[0];

}

void L2L_5(double x, double y, double z, double * L, double * Ls) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
double x49;
double x50;
double x51;
double x52;
double x53;
double x54;
double x55;
double x56;
double x57;
double x58;
double x59;
double x60;
double x61;
double x62;
double x63;
double x64;
double x65;
double x66;
double x67;
double x68;
double x69;
double x70;
double x71;
double x72;
double x73;
double x74;
double x75;
double x76;
double x77;
double x78;
double x79;
double x80;
double x81;
double x82;
double x83;
double x84;
double x85;
double x86;
double x87;
double x88;
double x89;
double x90;
double x91;
double x92;
double x93;
double x94;
double x95;
double x96;
double x97;
double x98;
double x99;
double x100;
double x101;
double x102;
double x103;
double x104;
double x105;
double x106;
double x107;
double x108;
double x109;
double x110;
double x111;
double x112;
double x113;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
x3 = z*L[14];
x4 = x3*y;
x5 = (x*x);
x6 = (1.0/2.0)*x5;
x7 = (x*x*x);
x8 = (1.0/6.0)*x7;
x9 = (1.0/24.0)*(x*x*x*x);
x10 = (y*y);
x11 = (1.0/2.0)*x10;
x12 = (y*y*y);
x13 = (1.0/6.0)*x12;
x14 = (1.0/24.0)*(y*y*y*y);
x15 = (z*z);
x16 = (1.0/2.0)*x15;
x17 = (z*z*z);
x18 = (1.0/6.0)*x17;
x19 = (1.0/24.0)*(z*z*z*z);
x20 = x*L[13];
x21 = x*L[26];
x22 = x*L[45];
x23 = x*L[15];
x24 = x*L[29];
x25 = x*L[49];
x26 = y*L[11];
x27 = z*L[12];
x28 = y*L[21];
x29 = z*L[22];
x30 = y*L[36];
x31 = z*L[37];
x32 = y*L[18];
x33 = y*L[33];
x34 = y*L[54];
x35 = z*L[17];
x36 = z*L[31];
x37 = z*L[51];
x38 = y*L[28];
x39 = x*x38;
x40 = y*L[48];
x41 = x*x40;
x42 = z*L[27];
x43 = x*x42;
x44 = z*L[46];
x45 = x*x44;
x46 = z*L[24];
x47 = x46*y;
x48 = z*L[39];
x49 = x48*y;
x50 = (1.0/4.0)*x5;
x51 = x10*x50;
x52 = (1.0/12.0)*x5;
x53 = x15*x50;
x54 = (1.0/12.0)*x7;
x55 = (1.0/4.0)*x10*x15;
x56 = x*L[47];
x57 = y*L[43];
x58 = z*L[42];
x59 = x*L[23];
x60 = x*L[41];
x61 = x*L[25];
x62 = x*L[44];
x63 = x*x57;
x64 = x*x58;
x65 = y*L[13];
x66 = x42*y;
x67 = x*L[28];
x68 = x*L[48];
x69 = y*L[23];
x70 = y*L[38];
x71 = y*L[32];
x72 = y*L[53];
x73 = y*L[47];
x74 = x*x73;
x75 = x58*y;
x76 = y*L[14];
x77 = z*L[15];
x78 = z*L[18];
x79 = z*L[28];
x80 = x79*y;
x81 = x*L[27];
x82 = x*L[46];
x83 = y*L[24];
x84 = z*L[25];
x85 = y*L[39];
x86 = z*L[40];
x87 = z*L[32];
x88 = z*L[52];
x89 = z*L[47];
x90 = x*x89;
x91 = z*L[43];
x92 = x91*y;
x93 = x*L[38];
x94 = x*L[40];
x95 = x*L[43];
x96 = x*L[42];
x97 = y*L[26];
x98 = x44*y;
x99 = y*L[41];
x100 = y*L[52];
x101 = y*L[27];
x102 = x89*y;
x103 = y*L[42];
x104 = z*L[29];
x105 = z*L[33];
x106 = z*L[48];
x107 = x106*y;
x108 = z*L[44];
x109 = z*L[53];
x110 = y*L[45];
x111 = y*L[46];
x112 = z*L[49];
x113 = z*L[54];
#pragma omp atomic
Ls[0] += (1.0/120.0)*(x*x*x*x*x)*L[35] + x*x0 + x*x1 + x*x4 + x*L[1] + (1.0/12.0)*x10*x17*L[53] + x10*x54*L[38] + x11*x20 + x11*x35 + x11*x43 + x11*L[7] + (1.0/12.0)*x12*x15*L[52] + x12*x52*L[41] + x13*x21 + x13*x36 + x13*x45 + x13*L[16] + x14*x22 + x14*x37 + x14*L[30] + x15*x54*L[40] + x16*x23 + x16*x32 + x16*x39 + x16*L[9] + x17*x52*L[44] + x18*x24 + x18*x33 + x18*x41 + x18*L[19] + x19*x25 + x19*x34 + x19*L[34] + x2*y + x26*x6 + x27*x6 + x28*x8 + x29*x8 + x30*x9 + x31*x9 + x47*x6 + x49*x8 + x51*x58 + x51*L[23] + x53*x57 + x53*L[25] + x55*x56 + x55*L[32] + x6*L[4] + x8*L[10] + x9*L[20] + (1.0/120.0)*(y*y*y*y*y)*L[50] + y*L[2] + (1.0/120.0)*(z*z*z*z*z)*L[55] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*x26 + x*x27 + x*x47 + x*L[4] + x0 + x1 + x11*x42 + x11*x59 + x11*x64 + x11*L[13] + x13*x44 + x13*x60 + x13*L[26] + x14*L[45] + x16*x38 + x16*x61 + x16*x63 + x16*L[15] + x18*x40 + x18*x62 + x18*L[29] + x19*L[49] + x28*x6 + x29*x6 + x30*x8 + x31*x8 + x4 + x49*x6 + x51*L[38] + x53*L[40] + x55*L[47] + x6*L[10] + x8*L[20] + x9*L[35] + L[1];
#pragma omp atomic
Ls[2] += x*x3 + x*x65 + x*x66 + x*L[5] + x11*x21 + x11*x36 + x11*x45 + x11*L[16] + x13*x22 + x13*x37 + x13*L[30] + x14*L[50] + x16*x67 + x16*x71 + x16*x74 + x16*L[18] + x18*x68 + x18*x72 + x18*L[33] + x19*L[54] + x2 + x35*y + x46*x6 + x48*x8 + x51*L[41] + x53*L[43] + x55*L[52] + x6*x69 + x6*x75 + x6*L[11] + x70*x8 + x8*L[21] + x9*L[36] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*x76 + x*x77 + x*x80 + x*L[6] + x11*x81 + x11*x87 + x11*x90 + x11*L[17] + x13*x82 + x13*x88 + x13*L[31] + x14*L[51] + x16*x24 + x16*x33 + x16*x41 + x16*L[19] + x18*x25 + x18*x34 + x18*L[34] + x19*L[55] + x51*L[42] + x53*L[44] + x55*L[53] + x6*x83 + x6*x84 + x6*x92 + x6*L[12] + x78*y + x8*x85 + x8*x86 + x8*L[22] + x9*L[37] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += x*x28 + x*x29 + x*x49 + x*L[10] + x11*x58 + x11*x93 + x11*L[23] + x13*L[41] + x16*x57 + x16*x94 + x16*L[25] + x18*L[44] + x26 + x27 + x30*x6 + x31*x6 + x47 + x6*L[20] + x8*L[35] + L[4];
#pragma omp atomic
Ls[5] += x*x46 + x*x69 + x*x75 + x*L[11] + x11*x44 + x11*x60 + x11*L[26] + x13*L[45] + x16*x73 + x16*x95 + x16*L[28] + x18*L[48] + x3 + x48*x6 + x6*x70 + x6*L[21] + x65 + x66 + x8*L[36] + L[5];
#pragma omp atomic
Ls[6] += x*x83 + x*x84 + x*x92 + x*L[12] + x11*x89 + x11*x96 + x11*L[27] + x13*L[46] + x16*x40 + x16*x62 + x16*L[29] + x18*L[49] + x6*x85 + x6*x86 + x6*L[22] + x76 + x77 + x8*L[37] + x80 + L[6];
#pragma omp atomic
Ls[7] += x*x97 + x*x98 + x100*x16 + x11*x22 + x11*x37 + x11*L[30] + x13*L[50] + x16*x56 + x16*L[32] + x18*L[53] + x20 + x35 + x36*y + x43 + x58*x6 + x6*x99 + x6*L[23] + x8*L[38] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += x*x101 + x*x102 + x*x79 + x*L[14] + x103*x6 + x11*x82 + x11*x88 + x11*L[31] + x13*L[51] + x16*x68 + x16*x72 + x16*L[33] + x18*L[54] + x6*x91 + x6*L[24] + x78 + x8*L[39] + x87*y + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += x*x104 + x*x107 + x105*y + x108*x6 + x109*x11 + x11*x56 + x11*L[32] + x13*L[52] + x16*x25 + x16*x34 + x16*L[34] + x18*L[55] + x23 + x32 + x39 + x57*x6 + x6*L[25] + x8*L[40] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += x*x30 + x*x31 + x*L[20] + x11*L[38] + x16*L[40] + x28 + x29 + x49 + x6*L[35] + L[10];
#pragma omp atomic
Ls[11] += x*x48 + x*x70 + x*L[21] + x11*L[41] + x16*L[43] + x46 + x6*L[36] + x69 + x75 + L[11];
#pragma omp atomic
Ls[12] += x*x85 + x*x86 + x*L[22] + x11*L[42] + x16*L[44] + x6*L[37] + x83 + x84 + x92 + L[12];
#pragma omp atomic
Ls[13] += x*x99 + x11*L[45] + x16*L[47] + x42 + x59 + x6*L[38] + x64 + x97 + x98 + L[13];
#pragma omp atomic
Ls[14] += x*x103 + x*x91 + x*L[24] + x101 + x102 + x11*L[46] + x16*L[48] + x6*L[39] + x79 + L[14];
#pragma omp atomic
Ls[15] += x*x108 + x104 + x107 + x11*L[47] + x16*L[49] + x38 + x6*L[40] + x61 + x63 + L[15];
#pragma omp atomic
Ls[16] += x*x110 + x11*L[50] + x16*L[52] + x21 + x36 + x37*y + x45 + x6*L[41] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += x*x111 + x11*L[51] + x16*L[53] + x6*L[42] + x81 + x87 + x88*y + x90 + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += x*x106 + x105 + x109*y + x11*L[52] + x16*L[54] + x6*L[43] + x67 + x71 + x74 + L[18];
#pragma omp atomic
Ls[19] += x*x112 + x11*L[53] + x113*y + x16*L[55] + x24 + x33 + x41 + x6*L[44] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += x*L[35] + x30 + x31 + L[20];
#pragma omp atomic
Ls[21] += x*L[36] + x48 + x70 + L[21];
#pragma omp atomic
Ls[22] += x*L[37] + x85 + x86 + L[22];
#pragma omp atomic
Ls[23] += x58 + x93 + x99 + L[23];
#pragma omp atomic
Ls[24] += x*L[39] + x103 + x91 + L[24];
#pragma omp atomic
Ls[25] += x108 + x57 + x94 + L[25];
#pragma omp atomic
Ls[26] += x110 + x44 + x60 + L[26];
#pragma omp atomic
Ls[27] += x111 + x89 + x96 + L[27];
#pragma omp atomic
Ls[28] += x106 + x73 + x95 + L[28];
#pragma omp atomic
Ls[29] += x112 + x40 + x62 + L[29];
#pragma omp atomic
Ls[30] += x22 + x37 + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += x82 + x88 + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += x100 + x109 + x56 + L[32];
#pragma omp atomic
Ls[33] += x113 + x68 + x72 + L[33];
#pragma omp atomic
Ls[34] += x25 + x34 + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += L[35];
#pragma omp atomic
Ls[36] += L[36];
#pragma omp atomic
Ls[37] += L[37];
#pragma omp atomic
Ls[38] += L[38];
#pragma omp atomic
Ls[39] += L[39];
#pragma omp atomic
Ls[40] += L[40];
#pragma omp atomic
Ls[41] += L[41];
#pragma omp atomic
Ls[42] += L[42];
#pragma omp atomic
Ls[43] += L[43];
#pragma omp atomic
Ls[44] += L[44];
#pragma omp atomic
Ls[45] += L[45];
#pragma omp atomic
Ls[46] += L[46];
#pragma omp atomic
Ls[47] += L[47];
#pragma omp atomic
Ls[48] += L[48];
#pragma omp atomic
Ls[49] += L[49];
#pragma omp atomic
Ls[50] += L[50];
#pragma omp atomic
Ls[51] += L[51];
#pragma omp atomic
Ls[52] += L[52];
#pragma omp atomic
Ls[53] += L[53];
#pragma omp atomic
Ls[54] += L[54];
#pragma omp atomic
Ls[55] += L[55];

}

void L2P_5(double x, double y, double z, double * L, double * F) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
x0 = x*y;
x1 = x*z;
x2 = y*z;
x3 = (x*x);
x4 = (1.0/2.0)*x3;
x5 = (x*x*x);
x6 = (1.0/6.0)*x5;
x7 = (1.0/24.0)*(x*x*x*x);
x8 = (y*y);
x9 = (1.0/2.0)*x8;
x10 = (y*y*y);
x11 = (1.0/6.0)*x10;
x12 = (1.0/24.0)*(y*y*y*y);
x13 = (z*z);
x14 = (1.0/2.0)*x13;
x15 = (z*z*z);
x16 = (1.0/6.0)*x15;
x17 = (1.0/24.0)*(z*z*z*z);
x18 = (1.0/4.0)*x3;
x19 = x18*x8;
x20 = (1.0/12.0)*x3;
x21 = x13*x18;
x22 = (1.0/12.0)*x5;
x23 = (1.0/4.0)*x13*x8;
#pragma omp atomic
F[0] += (1.0/120.0)*(x*x*x*x*x)*L[35] + x*x11*L[26] + x*x12*L[45] + x*x14*L[15] + x*x16*L[29] + x*x17*L[49] + x*x23*L[47] + x*x9*L[13] + x*L[1] + x0*x14*L[28] + x0*x16*L[48] + x0*z*L[14] + x0*L[5] + x1*x11*L[46] + x1*x9*L[27] + x1*L[6] + (1.0/12.0)*x10*x13*L[52] + x10*x20*L[41] + x11*z*L[31] + x11*L[16] + x12*z*L[51] + x12*L[30] + x13*x22*L[40] + x14*y*L[18] + x14*L[9] + x15*x20*L[44] + (1.0/12.0)*x15*x8*L[53] + x16*y*L[33] + x16*L[19] + x17*y*L[54] + x17*L[34] + x19*z*L[42] + x19*L[23] + x2*x4*L[24] + x2*x6*L[39] + x2*L[8] + x21*y*L[43] + x21*L[25] + x22*x8*L[38] + x23*L[32] + x4*y*L[11] + x4*z*L[12] + x4*L[4] + x6*y*L[21] + x6*z*L[22] + x6*L[10] + x7*y*L[36] + x7*z*L[37] + x7*L[20] + x9*z*L[17] + x9*L[7] + (1.0/120.0)*(y*y*y*y*y)*L[50] + y*L[2] + (1.0/120.0)*(z*z*z*z*z)*L[55] + z*L[3] + L[0];

}

void M2P_5(double x, double y, double z, double * M, double * F) {
double x0;
double x1;
double x2;
double x3;
double x4;
double x5;
double x6;
double x7;
double x8;
double x9;
double x10;
double x11;
double x12;
double x13;
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
double x22;
double x23;
double x24;
double x25;
double x26;
double x27;
double x28;
double x29;
double x30;
double x31;
double x32;
double x33;
double x34;
double x35;
double x36;
double x37;
double x38;
double x39;
double x40;
double x41;
double x42;
double x43;
double x44;
double x45;
double x46;
double x47;
double x48;
double x49;
double x50;
double x51;
double x52;
double x53;
double x54;
double x55;
double x56;
double x57;
double x58;
double x59;
double x60;
double x61;
double x62;
double x63;
double x64;
double x65;
double x66;
double x67;
double x68;
double x69;
double x70;
double x71;
double x72;
double x73;
double x74;
double x75;
double x76;
double x77;
double x78;
double x79;
double x80;
double x81;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*x6;
x8 = y*z;
x9 = pow(x2, -3.5);
x10 = 15.0*x9;
x11 = x*x8;
x12 = -x4;
x13 = x0*x6;
x14 = x1*x6;
x15 = 9.0*x5;
x16 = x0*x10;
x17 = -x16;
x18 = x*(x15 + x17);
x19 = x1*x10;
x20 = x19 - x6;
x21 = 0.5*x;
x22 = x17 + x6;
x23 = x22*y;
x24 = -x19;
x25 = y*(x15 + x24);
x26 = x22*z;
x27 = z*(x24 + x6);
x28 = 45.0*x9;
x29 = -x28;
x30 = pow(x2, -4.5);
x31 = 105.0*x30;
x32 = x0*x31;
x33 = x29 + x32;
x34 = x*x33;
x35 = x34*y;
x36 = x1*x31;
x37 = x29 + x36;
x38 = x37*y;
x39 = x34*z;
x40 = -x10;
x41 = x36 + x40;
x42 = x41*z;
x43 = x8*(x32 + x40);
x44 = x37*x8;
x45 = 315.0*x30;
x46 = 945.0*pow(x2, -5.5);
x47 = x0*x46;
x48 = x11*(x45 - x47);
x49 = x1*x46 - x45;
x50 = 90.0*x9;
x51 = (x*x*x*x);
x52 = x31*x51;
x53 = (y*y*y*y);
x54 = x31*x53;
x55 = -225.0*x9;
x56 = 1050.0*x30;
x57 = -x46*x51;
x58 = x*(x0*x56 + x55 + x57);
x59 = 630.0*x30;
x60 = x1*x59;
x61 = x46*x53;
x62 = x28 - x60 + x61;
x63 = x0*x59 + x29 + x57;
x64 = x63*y;
x65 = -x61;
x66 = y*(x1*x56 + x55 + x65);
x67 = x63*z;
x68 = z*(x29 + x60 + x65);
x69 = x1*x32;
x70 = -x1*x47;
x71 = x*(x1*x45 + x33 + x70);
x72 = y*(x0*x45 + x37 + x70);
x73 = z*(x32 + x41 + x70);
x74 = 1.0*x;
x75 = 105.0*x9;
x76 = -12.0*x5 - x69;
x77 = 120.0*x9;
x78 = -x71;
x79 = x62*x74;
x80 = -x72;
x81 = -x73;
#pragma omp atomic
F[0] += 0.16666666666666666*x*x38*M[26] - x*x4*M[1] - 0.041666666666666664*x*x62*M[45] - x10*x11*M[14] - 0.16666666666666666*x11*x49*M[46] + (1.0/6.0)*x18*M[10] + pow(x2, -0.5)*M[0] - x20*x21*M[13] + x21*x42*M[27] + (1.0/2.0)*x23*M[11] + (1.0/6.0)*x25*M[16] + (1.0/2.0)*x26*M[12] + (1.0/2.0)*x27*M[17] + (1.0/6.0)*x35*M[21] + (1.0/6.0)*x39*M[22] - x4*y*M[2] - x4*z*M[3] + (1.0/2.0)*x43*M[24] + (1.0/6.0)*x44*M[31] + (1.0/6.0)*x48*M[39] + (1.0/120.0)*x58*M[35] + x6*x8*M[8] + (1.0/24.0)*x64*M[36] + (1.0/120.0)*x66*M[50] + (1.0/24.0)*x67*M[37] + (1.0/24.0)*x68*M[51] + x7*y*M[5] + x7*z*M[6] + (1.0/12.0)*x71*M[38] + (1.0/12.0)*x72*M[41] + (1.0/4.0)*x73*M[42] + (1.0/2.0)*(x12 + x13)*M[4] + (1.0/2.0)*(x12 + x14)*M[7] + (1.0/2.0)*(-x18 + x20*x74)*M[15] + (1.0/2.0)*(-x23 - x25)*M[18] + (1.0/6.0)*(-x26 - x27)*M[19] + (1.0/2.0)*(-x35 - x38*x74)*M[28] + (1.0/6.0)*(-x39 - x42*x74)*M[29] + (1.0/6.0)*(-x43 - x44)*M[33] + (1.0/6.0)*(-x48 + x49*x74*x8)*M[48] + (1.0/12.0)*(-x58 + x78)*M[40] + (1.0/4.0)*(-x64 + x80)*M[43] + (1.0/12.0)*(-x66 + x80)*M[52] + (1.0/12.0)*(-x67 + x81)*M[44] + (1.0/12.0)*(-x68 + x81)*M[53] + (1.0/4.0)*(x78 + x79)*M[47] + (1.0/2.0)*(-x13 - x14 + 2.0*x3)*M[9] + (1.0/4.0)*(x22 + x24 + x69)*M[23] + (1.0/24.0)*(x58 + 2*x71 - x79)*M[49] + (1.0/24.0)*(x64 + x66 + 2*x72)*M[54] + (1.0/120.0)*(x67 + x68 + 2*x73)*M[55] + (1.0/24.0)*(-x0*x50 + x15 + x52)*M[20] + (1.0/24.0)*(-x1*x50 + x15 + x54)*M[30] + (1.0/4.0)*(x0*x75 + x19 - x52 + x76)*M[25] + (1.0/4.0)*(x1*x75 + x16 - x54 + x76)*M[32] + (1.0/24.0)*(210.0*x0*x1*x30 - x0*x77 - x1*x77 + 24.0*x5 + x52 + x54)*M[34];

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
  case 5:
    P2M_5(x, y, z, q, M);
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
  case 5:
    M2M_5(x, y, z, M, Ms);
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
  case 5:
    M2L_5(x, y, z, M, L);
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
  case 5:
    L2L_5(x, y, z, L, Ls);
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
  case 5:
    L2P_5(x, y, z, L, F);
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
  case 5:
    M2P_5(x, y, z, M, F);
    break;
  }
}
