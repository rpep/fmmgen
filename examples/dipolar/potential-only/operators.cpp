#include "operators.h"
#include "math.h"
void P2M_1(double x, double y, double z, double q, double * M) {
M[0] += -q*x;
M[1] += -q*y;
M[2] += -q*z;
}
void M2M_1(double x, double y, double z, double * M, double * Ms) {
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];

}

void M2L_1(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double x0;
x0 = 1.0/(R*R*R);
#pragma omp atomic
L[0] += -x*x0*M[0] - x0*y*M[1] - x0*z*M[2];

}

void L2L_1(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += L[0];

}

void L2P_1(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += L[0];

}

void M2P_1(double x, double y, double z, double * M, double * F) {
double x0;
x0 = 1.0*pow((x*x) + (y*y) + (z*z), -1.5);
#pragma omp atomic
F[0] += -x*x0*M[0] - x0*y*M[1] - x0*z*M[2];

}

void P2P(double x, double y, double z, double * S, double * F) {
double R = sqrt(x*x + y*y + z*z);
double x0;
x0 = 1.0/(R*R*R);
#pragma omp atomic
F[0] += -x*x0*S[0] - x0*y*S[1] - x0*z*S[2];

}

void P2M_2(double x, double y, double z, double q, double * M) {
double x0;
double x1;
double x2;
x0 = q*x;
x1 = q*y;
x2 = (1.0/2.0)*q;
M[0] += -x0;
M[1] += -x1;
M[2] += -q*z;
M[3] += (x*x)*x2;
M[4] += x0*y;
M[5] += x0*z;
M[6] += x2*(y*y);
M[7] += x1*z;
M[8] += x2*(z*z);
}
void M2M_2(double x, double y, double z, double * M, double * Ms) {
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x*M[0] + M[3];
#pragma omp atomic
Ms[4] += x*M[1] + y*M[0] + M[4];
#pragma omp atomic
Ms[5] += x*M[2] + z*M[0] + M[5];
#pragma omp atomic
Ms[6] += y*M[1] + M[6];
#pragma omp atomic
Ms[7] += y*M[2] + z*M[1] + M[7];
#pragma omp atomic
Ms[8] += z*M[2] + M[8];

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
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*M[0];
x3 = y*M[1];
x4 = z*M[2];
x5 = 3.0/pow(R, 5);
x6 = x*x5;
x7 = x5*y;
x8 = -x1;
x9 = (x*x)*x5;
x10 = x8 + x9;
x11 = x5*(y*y);
x12 = x11 + x8;
x13 = 2.0*x0 - x11 - x9;
x14 = x5*z;
#pragma omp atomic
L[0] += -x1*x2 - x1*x3 - x1*x4 + x10*M[3] + x12*M[6] + x13*M[8] + x6*y*M[4] + x6*z*M[5] + x7*z*M[7];
#pragma omp atomic
L[1] += x10*M[0] + x3*x6 + x4*x6;
#pragma omp atomic
L[2] += x12*M[1] + x2*x7 + x4*x7;
#pragma omp atomic
L[3] += x13*M[2] + x14*x2 + x14*x3;

}

void L2L_2(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += x*L[1] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += L[1];
#pragma omp atomic
Ls[2] += L[2];
#pragma omp atomic
Ls[3] += L[3];

}

void L2P_2(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += x*L[1] + y*L[2] + z*L[3] + L[0];

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
F[0] += -x*x4*M[0] - x4*y*M[1] - x4*z*M[2] + x5*y*z*M[7] + x6*y*M[4] + x6*z*M[5] + (x7 + x8)*M[3] + (x7 + x9)*M[6] + (2.0*x3 - x8 - x9)*M[8];

}

void P2M_3(double x, double y, double z, double q, double * M) {
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
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = (y*y);
x7 = (z*z);
x8 = (1.0/6.0)*q;
x9 = (1.0/2.0)*x3;
x10 = (1.0/2.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x0*z;
M[6] += x4*x6;
M[7] += x1*z;
M[8] += x4*x7;
M[9] += -(x*x*x)*x8;
M[10] += -x1*x9;
M[11] += -x2*x9;
M[12] += -x10*x6;
M[13] += -x5*z;
M[14] += -x10*x7;
M[15] += -x8*(y*y*y);
M[16] += -1.0/2.0*x2*x6;
M[17] += -1.0/2.0*x1*x7;
M[18] += -x8*(z*z*z);
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
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = y*M[1];
x5 = y*M[2];
x6 = (1.0/2.0)*(x*x);
x7 = (y*y);
x8 = (1.0/2.0)*M[0];
x9 = (z*z);
x10 = (1.0/2.0)*x7;
x11 = (1.0/2.0)*x9;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + z*M[0] + M[5];
#pragma omp atomic
Ms[6] += x4 + M[6];
#pragma omp atomic
Ms[7] += x5 + z*M[1] + M[7];
#pragma omp atomic
Ms[8] += z*M[2] + M[8];
#pragma omp atomic
Ms[9] += x*M[3] + x6*M[0] + M[9];
#pragma omp atomic
Ms[10] += x*M[4] + x0*y + x6*M[1] + y*M[3] + M[10];
#pragma omp atomic
Ms[11] += x*M[5] + x0*z + x6*M[2] + z*M[3] + M[11];
#pragma omp atomic
Ms[12] += x*M[6] + x1*y + x7*x8 + y*M[4] + M[12];
#pragma omp atomic
Ms[13] += x*M[7] + x1*z + x2*z + x3*y + y*M[5] + z*M[4] + M[13];
#pragma omp atomic
Ms[14] += x*M[8] + x3*z + x8*x9 + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += x10*M[1] + y*M[6] + M[15];
#pragma omp atomic
Ms[16] += x10*M[2] + x4*z + y*M[7] + z*M[6] + M[16];
#pragma omp atomic
Ms[17] += x11*M[1] + x5*z + y*M[8] + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += x11*M[2] + z*M[8] + M[18];

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
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*M[0];
x3 = y*M[1];
x4 = z*M[2];
x5 = pow(R, -5);
x6 = 3.0*x5;
x7 = x*x6;
x8 = y*M[4];
x9 = z*M[5];
x10 = x6*y;
x11 = z*M[7];
x12 = 15.0/pow(R, 7);
x13 = x*x12*y;
x14 = -x1;
x15 = (x*x);
x16 = x15*x6;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x6;
x20 = x14 + x19;
x21 = 9.0*x5;
x22 = -x12*x15;
x23 = x*(x21 + x22);
x24 = x22 + x6;
x25 = x24*y;
x26 = x12*x18;
x27 = -x26;
x28 = y*(x21 + x27);
x29 = x24*z;
x30 = z*(x27 + x6);
x31 = 1.0*x*(x26 - x6);
x32 = 2.0*x0 - x16 - x19;
x33 = -x23 + x31;
x34 = -x25 - x28;
x35 = -x29 - x30;
x36 = x6*z;
x37 = x12*z;
x38 = x*x37;
#pragma omp atomic
L[0] += -x1*x2 - x1*x3 - x1*x4 + x10*x11 - x13*z*M[13] + x17*M[3] + x20*M[6] + x23*M[9] + x25*M[10] + x28*M[15] + x29*M[11] + x30*M[16] - x31*M[12] + x32*M[8] + x33*M[14] + x34*M[17] + x35*M[18] + x7*x8 + x7*x9;
#pragma omp atomic
L[1] += -x11*x13 + x17*M[0] + x23*M[3] + x25*M[4] + x29*M[5] + x3*x7 - x31*M[6] + x33*M[8] + x4*x7;
#pragma omp atomic
L[2] += x10*x2 + x10*x4 - x13*x9 + x20*M[1] + x25*M[3] + x28*M[6] + x30*M[7] - x31*M[4] + x34*M[8];
#pragma omp atomic
L[3] += x2*x36 + x29*M[3] + x3*x36 + x30*M[6] + x32*M[2] + x33*M[5] + x34*M[7] + x35*M[8] - x38*x8;
#pragma omp atomic
L[4] += x23*M[0] + x25*M[1] + x29*M[2];
#pragma omp atomic
L[5] += -x13*x4 + x25*M[0] - x31*M[1];
#pragma omp atomic
L[6] += x29*M[0] - x3*x38 + x33*M[2];
#pragma omp atomic
L[7] += x28*M[1] + x30*M[2] - x31*M[0];
#pragma omp atomic
L[8] += -x2*x37*y + x30*M[1] + x34*M[2];
#pragma omp atomic
L[9] += x33*M[0] + x34*M[1] + x35*M[2];

}

void L2L_3(double x, double y, double z, double * L, double * Ls) {
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

void L2P_3(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += (1.0/2.0)*(x*x)*L[4] + x*y*L[5] + x*z*L[6] + x*L[1] + (1.0/2.0)*(y*y)*L[7] + y*z*L[8] + y*L[2] + (1.0/2.0)*(z*z)*L[9] + z*L[3] + L[0];

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
x16 = x14 + x6;
x17 = x16*y;
x18 = x1*x9;
x19 = -x18;
x20 = y*(x13 + x19);
x21 = x16*z;
x22 = z*(x19 + x6);
x23 = 1.0*x*(x18 - x6);
#pragma omp atomic
F[0] += -x*x4*M[0] - x*x8*x9*M[13] + x15*M[9] + x17*M[10] + x20*M[15] + x21*M[11] + x22*M[16] - x23*M[12] - x4*y*M[1] - x4*z*M[2] + x6*x8*M[7] + x7*y*M[4] + x7*z*M[5] + (x10 + x11)*M[3] + (x10 + x12)*M[6] + (-x15 + x23)*M[14] + (-x17 - x20)*M[17] + (-x21 - x22)*M[18] + (-x11 - x12 + 2.0*x3)*M[8];

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
double x14;
double x15;
double x16;
double x17;
double x18;
double x19;
double x20;
double x21;
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = x1*z;
x9 = (z*z);
x10 = (x*x*x);
x11 = (1.0/6.0)*q;
x12 = (1.0/2.0)*x3;
x13 = (1.0/2.0)*x0;
x14 = (y*y*y);
x15 = (1.0/2.0)*x7;
x16 = (1.0/2.0)*x9;
x17 = (z*z*z);
x18 = (1.0/24.0)*q;
x19 = (1.0/6.0)*x10;
x20 = (1.0/4.0)*q*x3;
x21 = (1.0/6.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x6;
M[6] += x4*x7;
M[7] += x8;
M[8] += x4*x9;
M[9] += -x10*x11;
M[10] += -x1*x12;
M[11] += -x12*x2;
M[12] += -x13*x7;
M[13] += -x5*z;
M[14] += -x13*x9;
M[15] += -x11*x14;
M[16] += -x15*x2;
M[17] += -x1*x16;
M[18] += -x11*x17;
M[19] += (x*x*x*x)*x18;
M[20] += x1*x19;
M[21] += x19*x2;
M[22] += x20*x7;
M[23] += x12*x8;
M[24] += x20*x9;
M[25] += x14*x21;
M[26] += x15*x6;
M[27] += x16*x5;
M[28] += x17*x21;
M[29] += x18*(y*y*y*y);
M[30] += (1.0/6.0)*x14*x2;
M[31] += (1.0/4.0)*q*x7*x9;
M[32] += (1.0/6.0)*x1*x17;
M[33] += x18*(z*z*z*z);
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
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = z*M[0];
x5 = y*M[1];
x6 = y*M[2];
x7 = z*M[1];
x8 = z*M[2];
x9 = x*M[3];
x10 = (1.0/2.0)*(x*x);
x11 = x*M[4];
x12 = y*M[3];
x13 = x0*y;
x14 = x*M[5];
x15 = x*M[6];
x16 = y*M[4];
x17 = x1*y;
x18 = (y*y);
x19 = (1.0/2.0)*M[0];
x20 = x*M[7];
x21 = y*M[5];
x22 = x3*y;
x23 = x*M[8];
x24 = (z*z);
x25 = y*M[6];
x26 = (1.0/2.0)*x18;
x27 = y*M[7];
x28 = y*M[8];
x29 = (1.0/2.0)*x24;
x30 = (1.0/6.0)*(x*x*x);
x31 = (y*y*y);
x32 = (1.0/6.0)*M[0];
x33 = (z*z*z);
x34 = (1.0/6.0)*x31;
x35 = (1.0/6.0)*x33;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + x4 + M[5];
#pragma omp atomic
Ms[6] += x5 + M[6];
#pragma omp atomic
Ms[7] += x6 + x7 + M[7];
#pragma omp atomic
Ms[8] += x8 + M[8];
#pragma omp atomic
Ms[9] += x10*M[0] + x9 + M[9];
#pragma omp atomic
Ms[10] += x10*M[1] + x11 + x12 + x13 + M[10];
#pragma omp atomic
Ms[11] += x0*z + x10*M[2] + x14 + z*M[3] + M[11];
#pragma omp atomic
Ms[12] += x15 + x16 + x17 + x18*x19 + M[12];
#pragma omp atomic
Ms[13] += x1*z + x2*z + x20 + x21 + x22 + z*M[4] + M[13];
#pragma omp atomic
Ms[14] += x19*x24 + x23 + x3*z + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += x25 + x26*M[1] + M[15];
#pragma omp atomic
Ms[16] += x26*M[2] + x27 + x5*z + z*M[6] + M[16];
#pragma omp atomic
Ms[17] += x28 + x29*M[1] + x6*z + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += x29*M[2] + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += x*M[9] + x10*M[3] + x30*M[0] + M[19];
#pragma omp atomic
Ms[20] += x*M[10] + x10*x2 + x10*M[4] + x30*M[1] + x9*y + y*M[9] + M[20];
#pragma omp atomic
Ms[21] += x*M[11] + x10*x4 + x10*M[5] + x30*M[2] + x9*z + z*M[9] + M[21];
#pragma omp atomic
Ms[22] += x*M[12] + x0*x26 + x10*x5 + x10*M[6] + x11*y + x26*M[3] + y*M[10] + M[22];
#pragma omp atomic
Ms[23] += x*M[13] + x10*x6 + x10*x7 + x10*M[7] + x11*z + x12*z + x13*z + x14*y + y*M[11] + z*M[10] + M[23];
#pragma omp atomic
Ms[24] += x*M[14] + x0*x29 + x10*x8 + x10*M[8] + x14*z + x29*M[3] + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += x*M[15] + x1*x26 + x15*y + x26*M[4] + x31*x32 + y*M[12] + M[25];
#pragma omp atomic
Ms[26] += x*M[16] + x15*z + x16*z + x17*z + x20*y + x26*x3 + x26*x4 + x26*M[5] + y*M[13] + z*M[12] + M[26];
#pragma omp atomic
Ms[27] += x*M[17] + x1*x29 + x2*x29 + x20*z + x21*z + x22*z + x23*y + x29*M[4] + y*M[14] + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += x*M[18] + x23*z + x29*x3 + x29*M[5] + x32*x33 + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += x26*M[6] + x34*M[1] + y*M[15] + M[29];
#pragma omp atomic
Ms[30] += x25*z + x26*x7 + x26*M[7] + x34*M[2] + y*M[16] + z*M[15] + M[30];
#pragma omp atomic
Ms[31] += x26*x8 + x26*M[8] + x27*z + x29*x5 + x29*M[6] + y*M[17] + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += x28*z + x29*x6 + x29*M[7] + x35*M[1] + y*M[18] + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += x29*M[8] + x35*M[2] + z*M[18] + M[33];

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
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*M[0];
x3 = y*M[1];
x4 = z*M[2];
x5 = pow(R, -5);
x6 = 3.0*x5;
x7 = x*x6;
x8 = y*z;
x9 = x6*x8;
x10 = pow(R, -7);
x11 = 15.0*x10;
x12 = x11*x8;
x13 = x*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x6;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x6;
x20 = x14 + x19;
x21 = 9.0*x5;
x22 = x11*x15;
x23 = -x22;
x24 = x*(x21 + x23);
x25 = x23 + x6;
x26 = x25*y;
x27 = x11*x18;
x28 = -x27;
x29 = y*(x21 + x28);
x30 = x25*z;
x31 = z*(x28 + x6);
x32 = 1.0*x;
x33 = x32*(x27 - x6);
x34 = -45.0*x10;
x35 = pow(R, -9);
x36 = 105.0*x35;
x37 = x15*x36;
x38 = x*(x34 + x37);
x39 = x38*y;
x40 = x38*z;
x41 = x18*x36;
x42 = x34 + x41;
x43 = x42*x8;
x44 = -x11;
x45 = x8*(x37 + x44);
x46 = x32*x42*y;
x47 = x32*z*(x41 + x44);
x48 = (x*x*x*x)*x36;
x49 = 90.0*x10;
x50 = -x15*x49 + x21 + x48;
x51 = x36*(y*y*y*y);
x52 = -x18*x49 + x21 + x51;
x53 = 2.0*x0 - x16 - x19;
x54 = x18*x37;
x55 = x25 + x28 + x54;
x56 = -x24 + x33;
x57 = -x26 - x29;
x58 = -x30 - x31;
x59 = -x39 - x46;
x60 = -x40 - x47;
x61 = -x43 - x45;
x62 = 105.0*x10;
x63 = -12.0*x5 - x54;
x64 = x18*x62 + x22 - x51 + x63;
x65 = x15*x62 + x27 - x48 + x63;
x66 = 120.0*x10;
x67 = 210.0*x15*x18*x35 - x15*x66 - x18*x66 + x48 + 24.0*x5 + x51;
x68 = x2*x6;
#pragma omp atomic
L[0] += -x1*x2 - x1*x3 - x1*x4 - x13*M[13] + x17*M[3] + x20*M[6] + x24*M[9] + x26*M[10] + x29*M[15] + x30*M[11] + x31*M[16] - x33*M[12] + x39*M[20] + x40*M[21] + x43*M[30] + x45*M[23] + x46*M[25] + x47*M[26] + x50*M[19] + x52*M[29] + x53*M[8] + x55*M[22] + x56*M[14] + x57*M[17] + x58*M[18] + x59*M[27] + x60*M[28] + x61*M[32] + x64*M[31] + x65*M[24] + x67*M[33] + x7*y*M[4] + x7*z*M[5] + x9*M[7];
#pragma omp atomic
L[1] += -x13*M[7] + x17*M[0] + x24*M[3] + x26*M[4] + x3*x7 + x30*M[5] - x33*M[6] + x39*M[10] + x4*x7 + x40*M[11] + x45*M[13] + x46*M[15] + x47*M[16] + x50*M[9] + x55*M[12] + x56*M[8] + x59*M[17] + x60*M[18] + x65*M[14];
#pragma omp atomic
L[2] += -x13*M[5] + x20*M[1] + x26*M[3] + x29*M[6] + x31*M[7] - x33*M[4] + x39*M[9] + x43*M[16] + x45*M[11] + x46*M[12] + x47*M[13] + x52*M[15] + x55*M[10] + x57*M[8] + x59*M[14] + x61*M[18] + x64*M[17] + x68*y + x9*M[2];
#pragma omp atomic
L[3] += -x13*M[4] + x30*M[3] + x31*M[6] + x40*M[9] + x43*M[15] + x45*M[10] + x47*M[12] + x53*M[2] + x56*M[5] + x57*M[7] + x58*M[8] + x59*M[13] + x60*M[14] + x61*M[17] + x64*M[16] + x65*M[11] + x67*M[18] + x68*z + x9*M[1];
#pragma omp atomic
L[4] += x24*M[0] + x26*M[1] + x30*M[2] + x39*M[4] + x40*M[5] + x45*M[7] + x50*M[3] + x55*M[6] + x65*M[8];
#pragma omp atomic
L[5] += -x13*M[2] + x26*M[0] - x33*M[1] + x39*M[3] + x45*M[5] + x46*M[6] + x47*M[7] + x55*M[4] + x59*M[8];
#pragma omp atomic
L[6] += -x13*M[1] + x30*M[0] + x40*M[3] + x45*M[4] + x47*M[6] + x56*M[2] + x59*M[7] + x60*M[8] + x65*M[5];
#pragma omp atomic
L[7] += x29*M[1] + x31*M[2] - x33*M[0] + x43*M[7] + x46*M[4] + x47*M[5] + x52*M[6] + x55*M[3] + x64*M[8];
#pragma omp atomic
L[8] += -x12*x2 + x31*M[1] + x43*M[6] + x45*M[3] + x47*M[4] + x57*M[2] + x59*M[5] + x61*M[8] + x64*M[7];
#pragma omp atomic
L[9] += x56*M[0] + x57*M[1] + x58*M[2] + x59*M[4] + x60*M[5] + x61*M[7] + x64*M[6] + x65*M[3] + x67*M[8];
#pragma omp atomic
L[10] += x39*M[1] + x40*M[2] + x50*M[0];
#pragma omp atomic
L[11] += x39*M[0] + x45*M[2] + x55*M[1];
#pragma omp atomic
L[12] += x40*M[0] + x45*M[1] + x65*M[2];
#pragma omp atomic
L[13] += x46*M[1] + x47*M[2] + x55*M[0];
#pragma omp atomic
L[14] += x45*M[0] + x47*M[1] + x59*M[2];
#pragma omp atomic
L[15] += x59*M[1] + x60*M[2] + x65*M[0];
#pragma omp atomic
L[16] += x43*M[2] + x46*M[0] + x52*M[1];
#pragma omp atomic
L[17] += x43*M[1] + x47*M[0] + x64*M[2];
#pragma omp atomic
L[18] += x59*M[0] + x61*M[2] + x64*M[1];
#pragma omp atomic
L[19] += x60*M[0] + x61*M[1] + x67*M[2];

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

void L2P_4(double x, double y, double z, double * L, double * F) {
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
x18 = x16 + x6;
x19 = x18*y;
x20 = x1*x10;
x21 = -x20;
x22 = y*(x14 + x21);
x23 = x18*z;
x24 = z*(x21 + x6);
x25 = 1.0*x;
x26 = x25*(x20 - x6);
x27 = -45.0*x9;
x28 = pow(x2, -4.5);
x29 = 105.0*x28;
x30 = x0*x29;
x31 = x*(x27 + x30);
x32 = x31*y;
x33 = x31*z;
x34 = -x10;
x35 = x8*(x30 + x34);
x36 = x1*x29;
x37 = x27 + x36;
x38 = x37*x8;
x39 = x25*x37*y;
x40 = x25*z*(x34 + x36);
x41 = 90.0*x9;
x42 = (x*x*x*x)*x29;
x43 = x29*(y*y*y*y);
x44 = x1*x30;
x45 = 105.0*x9;
x46 = -x44 - 12.0*x5;
x47 = 120.0*x9;
#pragma omp atomic
F[0] += -x*x10*x8*M[13] - x*x4*M[0] + x17*M[9] + x19*M[10] + x22*M[15] + x23*M[11] + x24*M[16] - x26*M[12] + x32*M[20] + x33*M[21] + x35*M[23] + x38*M[30] + x39*M[25] - x4*y*M[1] - x4*z*M[2] + x40*M[26] + x6*x8*M[7] + x7*y*M[4] + x7*z*M[5] + (x11 + x12)*M[3] + (x11 + x13)*M[6] + (-x17 + x26)*M[14] + (-x19 - x22)*M[17] + (-x23 - x24)*M[18] + (-x32 - x39)*M[27] + (-x33 - x40)*M[28] + (-x35 - x38)*M[32] + (-x12 - x13 + 2.0*x3)*M[8] + (x18 + x21 + x44)*M[22] + (-x0*x41 + x14 + x42)*M[19] + (-x1*x41 + x14 + x43)*M[29] + (x0*x45 + x20 - x42 + x46)*M[24] + (x1*x45 + x15 - x43 + x46)*M[31] + (210.0*x0*x1*x28 - x0*x47 - x1*x47 + x42 + x43 + 24.0*x5)*M[33];

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
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = x1*z;
x9 = (z*z);
x10 = (x*x*x);
x11 = (1.0/6.0)*q;
x12 = (1.0/2.0)*x3;
x13 = (1.0/2.0)*x0;
x14 = (y*y*y);
x15 = (1.0/2.0)*x7;
x16 = (1.0/2.0)*x9;
x17 = (z*z*z);
x18 = (x*x*x*x);
x19 = (1.0/24.0)*q;
x20 = (1.0/6.0)*x10;
x21 = q*x7;
x22 = (1.0/4.0)*x3;
x23 = q*x9;
x24 = (1.0/6.0)*x0;
x25 = (y*y*y*y);
x26 = (1.0/6.0)*x14;
x27 = (1.0/4.0)*x9;
x28 = (1.0/6.0)*x17;
x29 = (z*z*z*z);
x30 = (1.0/120.0)*q;
x31 = (1.0/24.0)*x18;
x32 = (1.0/12.0)*x10;
x33 = (1.0/12.0)*x14;
x34 = q*x3;
x35 = (1.0/12.0)*x17;
x36 = (1.0/24.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x6;
M[6] += x4*x7;
M[7] += x8;
M[8] += x4*x9;
M[9] += -x10*x11;
M[10] += -x1*x12;
M[11] += -x12*x2;
M[12] += -x13*x7;
M[13] += -x5*z;
M[14] += -x13*x9;
M[15] += -x11*x14;
M[16] += -x15*x2;
M[17] += -x1*x16;
M[18] += -x11*x17;
M[19] += x18*x19;
M[20] += x1*x20;
M[21] += x2*x20;
M[22] += x21*x22;
M[23] += x12*x8;
M[24] += x22*x23;
M[25] += x14*x24;
M[26] += x15*x6;
M[27] += x16*x5;
M[28] += x17*x24;
M[29] += x19*x25;
M[30] += x2*x26;
M[31] += x21*x27;
M[32] += x1*x28;
M[33] += x19*x29;
M[34] += -pow(x, 5)*x30;
M[35] += -x1*x31;
M[36] += -x2*x31;
M[37] += -x21*x32;
M[38] += -x20*x8;
M[39] += -x23*x32;
M[40] += -x33*x34;
M[41] += -x2*x22*x7;
M[42] += -x1*x22*x9;
M[43] += -x34*x35;
M[44] += -x25*x36;
M[45] += -x26*x6;
M[46] += -x0*x27*x7;
M[47] += -x28*x5;
M[48] += -x29*x36;
M[49] += -x30*pow(y, 5);
M[50] += -1.0/24.0*x2*x25;
M[51] += -x23*x33;
M[52] += -x21*x35;
M[53] += -1.0/24.0*x1*x29;
M[54] += -x30*pow(z, 5);
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
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = z*M[0];
x5 = y*M[1];
x6 = y*M[2];
x7 = z*M[1];
x8 = z*M[2];
x9 = x*M[3];
x10 = (x*x);
x11 = (1.0/2.0)*x10;
x12 = x*M[4];
x13 = y*M[3];
x14 = x0*y;
x15 = x*M[5];
x16 = z*M[3];
x17 = x0*z;
x18 = x*M[6];
x19 = y*M[4];
x20 = x1*y;
x21 = (y*y);
x22 = (1.0/2.0)*M[0];
x23 = x*M[7];
x24 = y*M[5];
x25 = z*M[4];
x26 = x3*y;
x27 = x1*z;
x28 = x2*z;
x29 = x*M[8];
x30 = z*M[5];
x31 = x3*z;
x32 = (z*z);
x33 = y*M[6];
x34 = (1.0/2.0)*x21;
x35 = y*M[7];
x36 = z*M[6];
x37 = x5*z;
x38 = y*M[8];
x39 = z*M[7];
x40 = x6*z;
x41 = (1.0/2.0)*x32;
x42 = z*M[8];
x43 = x*M[9];
x44 = (1.0/6.0)*(x*x*x);
x45 = x*M[10];
x46 = y*M[9];
x47 = x9*y;
x48 = x*M[11];
x49 = x*M[12];
x50 = y*M[10];
x51 = x12*y;
x52 = x*M[13];
x53 = y*M[11];
x54 = x15*y;
x55 = x*M[14];
x56 = x*M[15];
x57 = y*M[12];
x58 = x18*y;
x59 = (y*y*y);
x60 = (1.0/6.0)*M[0];
x61 = x*M[16];
x62 = y*M[13];
x63 = x23*y;
x64 = x*M[17];
x65 = y*M[14];
x66 = x29*y;
x67 = x*M[18];
x68 = (z*z*z);
x69 = y*M[15];
x70 = (1.0/6.0)*x59;
x71 = y*M[16];
x72 = y*M[17];
x73 = y*M[18];
x74 = (1.0/6.0)*x68;
x75 = (1.0/24.0)*(x*x*x*x);
x76 = (1.0/4.0)*x10;
x77 = x76*M[0];
x78 = x21*x76;
x79 = x32*x76;
x80 = (y*y*y*y);
x81 = (1.0/24.0)*M[0];
x82 = (1.0/4.0)*x21*x32;
x83 = (z*z*z*z);
x84 = (1.0/24.0)*x80;
x85 = (1.0/24.0)*x83;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + x4 + M[5];
#pragma omp atomic
Ms[6] += x5 + M[6];
#pragma omp atomic
Ms[7] += x6 + x7 + M[7];
#pragma omp atomic
Ms[8] += x8 + M[8];
#pragma omp atomic
Ms[9] += x11*M[0] + x9 + M[9];
#pragma omp atomic
Ms[10] += x11*M[1] + x12 + x13 + x14 + M[10];
#pragma omp atomic
Ms[11] += x11*M[2] + x15 + x16 + x17 + M[11];
#pragma omp atomic
Ms[12] += x18 + x19 + x20 + x21*x22 + M[12];
#pragma omp atomic
Ms[13] += x23 + x24 + x25 + x26 + x27 + x28 + M[13];
#pragma omp atomic
Ms[14] += x22*x32 + x29 + x30 + x31 + M[14];
#pragma omp atomic
Ms[15] += x33 + x34*M[1] + M[15];
#pragma omp atomic
Ms[16] += x34*M[2] + x35 + x36 + x37 + M[16];
#pragma omp atomic
Ms[17] += x38 + x39 + x40 + x41*M[1] + M[17];
#pragma omp atomic
Ms[18] += x41*M[2] + x42 + M[18];
#pragma omp atomic
Ms[19] += x11*M[3] + x43 + x44*M[0] + M[19];
#pragma omp atomic
Ms[20] += x11*x2 + x11*M[4] + x44*M[1] + x45 + x46 + x47 + M[20];
#pragma omp atomic
Ms[21] += x11*x4 + x11*M[5] + x44*M[2] + x48 + x9*z + z*M[9] + M[21];
#pragma omp atomic
Ms[22] += x0*x34 + x11*x5 + x11*M[6] + x34*M[3] + x49 + x50 + x51 + M[22];
#pragma omp atomic
Ms[23] += x11*x6 + x11*x7 + x11*M[7] + x12*z + x13*z + x14*z + x52 + x53 + x54 + z*M[10] + M[23];
#pragma omp atomic
Ms[24] += x0*x41 + x11*x8 + x11*M[8] + x15*z + x41*M[3] + x55 + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += x1*x34 + x34*M[4] + x56 + x57 + x58 + x59*x60 + M[25];
#pragma omp atomic
Ms[26] += x18*z + x19*z + x20*z + x3*x34 + x34*x4 + x34*M[5] + x61 + x62 + x63 + z*M[12] + M[26];
#pragma omp atomic
Ms[27] += x1*x41 + x2*x41 + x23*z + x24*z + x26*z + x41*M[4] + x64 + x65 + x66 + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += x29*z + x3*x41 + x41*M[5] + x60*x68 + x67 + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += x34*M[6] + x69 + x70*M[1] + M[29];
#pragma omp atomic
Ms[30] += x33*z + x34*x7 + x34*M[7] + x70*M[2] + x71 + z*M[15] + M[30];
#pragma omp atomic
Ms[31] += x34*x8 + x34*M[8] + x35*z + x41*x5 + x41*M[6] + x72 + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += x38*z + x41*x6 + x41*M[7] + x73 + x74*M[1] + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += x41*M[8] + x74*M[2] + z*M[18] + M[33];
#pragma omp atomic
Ms[34] += x*M[19] + x11*M[9] + x44*M[3] + x75*M[0] + M[34];
#pragma omp atomic
Ms[35] += x*M[20] + x11*x13 + x11*M[10] + x2*x44 + x43*y + x44*M[4] + x75*M[1] + y*M[19] + M[35];
#pragma omp atomic
Ms[36] += x*M[21] + x11*x16 + x11*M[11] + x4*x44 + x43*z + x44*M[5] + x75*M[2] + z*M[19] + M[36];
#pragma omp atomic
Ms[37] += x*M[22] + x11*x19 + x11*M[12] + x21*x77 + x34*x9 + x34*M[9] + x44*x5 + x44*M[6] + x45*y + y*M[20] + M[37];
#pragma omp atomic
Ms[38] += x*M[23] + x11*x24 + x11*x25 + x11*x28 + x11*M[13] + x44*x6 + x44*x7 + x44*M[7] + x45*z + x46*z + x47*z + x48*y + y*M[21] + z*M[20] + M[38];
#pragma omp atomic
Ms[39] += x*M[24] + x11*x30 + x11*M[14] + x32*x77 + x41*x9 + x41*M[9] + x44*x8 + x44*M[8] + x48*z + z*M[21] + M[39];
#pragma omp atomic
Ms[40] += x*M[25] + x0*x70 + x11*x33 + x11*M[15] + x12*x34 + x34*M[10] + x49*y + x70*M[3] + x78*M[1] + y*M[22] + M[40];
#pragma omp atomic
Ms[41] += x*M[26] + x11*x35 + x11*x36 + x11*x37 + x11*M[16] + x15*x34 + x16*x34 + x17*x34 + x34*M[11] + x49*z + x50*z + x51*z + x52*y + x78*M[2] + y*M[23] + z*M[22] + M[41];
#pragma omp atomic
Ms[42] += x*M[27] + x11*x38 + x11*x39 + x11*x40 + x11*M[17] + x12*x41 + x13*x41 + x14*x41 + x41*M[10] + x52*z + x53*z + x54*z + x55*y + x79*M[1] + y*M[24] + z*M[23] + M[42];
#pragma omp atomic
Ms[43] += x*M[28] + x0*x74 + x11*x42 + x11*M[18] + x15*x41 + x41*M[11] + x55*z + x74*M[3] + x79*M[2] + z*M[24] + M[43];
#pragma omp atomic
Ms[44] += x*M[29] + x1*x70 + x18*x34 + x34*M[12] + x56*y + x70*M[4] + x80*x81 + y*M[25] + M[44];
#pragma omp atomic
Ms[45] += x*M[30] + x23*x34 + x25*x34 + x27*x34 + x3*x70 + x34*M[13] + x4*x70 + x56*z + x57*z + x58*z + x61*y + x70*M[5] + y*M[26] + z*M[25] + M[45];
#pragma omp atomic
Ms[46] += x*M[31] + x18*x41 + x19*x41 + x20*x41 + x29*x34 + x30*x34 + x31*x34 + x34*M[14] + x41*M[12] + x61*z + x62*z + x63*z + x64*y + x82*M[0] + y*M[27] + z*M[26] + M[46];
#pragma omp atomic
Ms[47] += x*M[32] + x1*x74 + x2*x74 + x23*x41 + x24*x41 + x26*x41 + x41*M[13] + x64*z + x65*z + x66*z + x67*y + x74*M[4] + y*M[28] + z*M[27] + M[47];
#pragma omp atomic
Ms[48] += x*M[33] + x29*x41 + x3*x74 + x41*M[14] + x67*z + x74*M[5] + x81*x83 + z*M[28] + M[48];
#pragma omp atomic
Ms[49] += x34*M[15] + x70*M[6] + x84*M[1] + y*M[29] + M[49];
#pragma omp atomic
Ms[50] += x34*x36 + x34*M[16] + x69*z + x7*x70 + x70*M[7] + x84*M[2] + y*M[30] + z*M[29] + M[50];
#pragma omp atomic
Ms[51] += x33*x41 + x34*x39 + x34*M[17] + x41*M[15] + x70*x8 + x70*M[8] + x71*z + x82*M[1] + y*M[31] + z*M[30] + M[51];
#pragma omp atomic
Ms[52] += x34*x42 + x34*M[18] + x35*x41 + x41*M[16] + x5*x74 + x72*z + x74*M[6] + x82*M[2] + y*M[32] + z*M[31] + M[52];
#pragma omp atomic
Ms[53] += x38*x41 + x41*M[17] + x6*x74 + x73*z + x74*M[7] + x85*M[1] + y*M[33] + z*M[32] + M[53];
#pragma omp atomic
Ms[54] += x41*M[18] + x74*M[8] + x85*M[2] + z*M[33] + M[54];

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
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = x*M[0];
x3 = y*M[1];
x4 = z*M[2];
x5 = pow(R, -5);
x6 = 3.0*x5;
x7 = x*x6;
x8 = y*z;
x9 = x6*x8;
x10 = pow(R, -7);
x11 = 15.0*x10;
x12 = x*x8;
x13 = x11*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x6;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x6;
x20 = x14 + x19;
x21 = 9.0*x5;
x22 = x11*x15;
x23 = -x22;
x24 = x*(x21 + x23);
x25 = x23 + x6;
x26 = x25*y;
x27 = x11*x18;
x28 = -x27;
x29 = y*(x21 + x28);
x30 = x25*z;
x31 = z*(x28 + x6);
x32 = 1.0*x;
x33 = x32*(x27 - x6);
x34 = 45.0*x10;
x35 = -x34;
x36 = pow(R, -9);
x37 = x15*x36;
x38 = 105.0*x37;
x39 = x35 + x38;
x40 = x*x39;
x41 = x40*y;
x42 = x40*z;
x43 = x18*x36;
x44 = 105.0*x43;
x45 = x35 + x44;
x46 = x45*x8;
x47 = -x11;
x48 = x8*(x38 + x47);
x49 = x32*x45*y;
x50 = x44 + x47;
x51 = x32*x50*z;
x52 = 315.0*x36;
x53 = 945.0/pow(R, 11);
x54 = x15*x53;
x55 = x12*(x52 - x54);
x56 = x32*x8*(x18*x53 - x52);
x57 = (x*x*x*x);
x58 = 105.0*x36;
x59 = x57*x58;
x60 = 90.0*x10;
x61 = -x15*x60 + x21 + x59;
x62 = (y*y*y*y);
x63 = x58*x62;
x64 = -x18*x60 + x21 + x63;
x65 = 2.0*x0 - x16 - x19;
x66 = -225.0*x10;
x67 = -x53*x57;
x68 = x*(1050.0*x37 + x66 + x67);
x69 = x53*x62;
x70 = -x69;
x71 = y*(1050.0*x43 + x66 + x70);
x72 = x35 + 630.0*x37 + x67;
x73 = x72*y;
x74 = x72*z;
x75 = 630.0*x43;
x76 = z*(x35 + x70 + x75);
x77 = x32*(x34 + x69 - x75);
x78 = x18*x38;
x79 = x25 + x28 + x78;
x80 = -x18*x54;
x81 = x*(x18*x52 + x39 + x80);
x82 = y*(x15*x52 + x45 + x80);
x83 = z*(x38 + x50 + x80);
x84 = -x24 + x33;
x85 = -x26 - x29;
x86 = -x30 - x31;
x87 = -x41 - x49;
x88 = -x42 - x51;
x89 = -x46 - x48;
x90 = -x55 + x56;
x91 = 105.0*x10;
x92 = -12.0*x5 - x78;
x93 = x18*x91 + x22 - x63 + x92;
x94 = x15*x91 + x27 - x59 + x92;
x95 = 120.0*x10;
x96 = -x15*x95 + 210.0*x18*x37 - x18*x95 + 24.0*x5 + x59 + x63;
x97 = -x81;
x98 = -x68 + x97;
x99 = x77 + x97;
x100 = -x82;
x101 = x100 - x71;
x102 = x100 - x73;
x103 = -x83;
x104 = x103 - x74;
x105 = x103 - x76;
x106 = x71 + x73 + 2*x82;
x107 = x74 + x76 + 2*x83;
x108 = x68 - x77 + 2*x81;
x109 = x2*x6;
#pragma omp atomic
L[0] += -x1*x2 - x1*x3 - x1*x4 + x101*M[51] + x102*M[42] + x104*M[43] + x105*M[52] + x106*M[53] + x107*M[54] + x108*M[48] - x13*M[13] + x17*M[3] + x20*M[6] + x24*M[9] + x26*M[10] + x29*M[15] + x30*M[11] + x31*M[16] - x33*M[12] + x41*M[20] + x42*M[21] + x46*M[30] + x48*M[23] + x49*M[25] + x51*M[26] + x55*M[38] - x56*M[45] + x61*M[19] + x64*M[29] + x65*M[8] + x68*M[34] + x7*y*M[4] + x7*z*M[5] + x71*M[49] + x73*M[35] + x74*M[36] + x76*M[50] - x77*M[44] + x79*M[22] + x81*M[37] + x82*M[40] + x83*M[41] + x84*M[14] + x85*M[17] + x86*M[18] + x87*M[27] + x88*M[28] + x89*M[32] + x9*M[7] + x90*M[47] + x93*M[31] + x94*M[24] + x96*M[33] + x98*M[39] + x99*M[46];
#pragma omp atomic
L[1] += x102*M[27] + x104*M[28] + x108*M[33] - x13*M[7] + x17*M[0] + x24*M[3] + x26*M[4] + x3*x7 + x30*M[5] - x33*M[6] + x4*x7 + x41*M[10] + x42*M[11] + x48*M[13] + x49*M[15] + x51*M[16] + x55*M[23] - x56*M[30] + x61*M[9] + x68*M[19] + x73*M[20] + x74*M[21] - x77*M[29] + x79*M[12] + x81*M[22] + x82*M[25] + x83*M[26] + x84*M[8] + x87*M[17] + x88*M[18] + x90*M[32] + x94*M[14] + x98*M[24] + x99*M[31];
#pragma omp atomic
L[2] += x101*M[31] + x102*M[24] + x105*M[32] + x106*M[33] + x109*y - x13*M[5] + x20*M[1] + x26*M[3] + x29*M[6] + x31*M[7] - x33*M[4] + x41*M[9] + x46*M[16] + x48*M[11] + x49*M[12] + x51*M[13] + x55*M[21] - x56*M[26] + x64*M[15] + x71*M[29] + x73*M[19] + x76*M[30] - x77*M[25] + x79*M[10] + x81*M[20] + x82*M[22] + x83*M[23] + x85*M[8] + x87*M[14] + x89*M[18] + x9*M[2] + x90*M[28] + x93*M[17] + x99*M[27];
#pragma omp atomic
L[3] += x101*M[30] + x102*M[23] + x104*M[24] + x105*M[31] + x106*M[32] + x107*M[33] + x108*M[28] + x109*z - x13*M[4] + x30*M[3] + x31*M[6] + x42*M[9] + x46*M[15] + x48*M[10] + x51*M[12] + x55*M[20] - x56*M[25] + x65*M[2] + x74*M[19] + x76*M[29] + x83*M[22] + x84*M[5] + x85*M[7] + x86*M[8] + x87*M[13] + x88*M[14] + x89*M[17] + x9*M[1] + x90*M[27] + x93*M[16] + x94*M[11] + x96*M[18] + x98*M[21] + x99*M[26];
#pragma omp atomic
L[4] += x102*M[17] + x104*M[18] + x24*M[0] + x26*M[1] + x30*M[2] + x41*M[4] + x42*M[5] + x48*M[7] + x55*M[13] + x61*M[3] + x68*M[9] + x73*M[10] + x74*M[11] + x79*M[6] + x81*M[12] + x82*M[15] + x83*M[16] + x94*M[8] + x98*M[14];
#pragma omp atomic
L[5] += x102*M[14] - x13*M[2] + x26*M[0] - x33*M[1] + x41*M[3] + x48*M[5] + x49*M[6] + x51*M[7] + x55*M[11] - x56*M[16] + x73*M[9] - x77*M[15] + x79*M[4] + x81*M[10] + x82*M[12] + x83*M[13] + x87*M[8] + x90*M[18] + x99*M[17];
#pragma omp atomic
L[6] += x102*M[13] + x104*M[14] + x108*M[18] - x13*M[1] + x30*M[0] + x42*M[3] + x48*M[4] + x51*M[6] + x55*M[10] - x56*M[15] + x74*M[9] + x83*M[12] + x84*M[2] + x87*M[7] + x88*M[8] + x90*M[17] + x94*M[5] + x98*M[11] + x99*M[16];
#pragma omp atomic
L[7] += x101*M[17] + x105*M[18] + x29*M[1] + x31*M[2] - x33*M[0] + x46*M[7] + x49*M[4] + x51*M[5] - x56*M[13] + x64*M[6] + x71*M[15] + x76*M[16] - x77*M[12] + x79*M[3] + x81*M[9] + x82*M[10] + x83*M[11] + x93*M[8] + x99*M[14];
#pragma omp atomic
L[8] += x101*M[16] + x102*M[11] + x105*M[17] + x106*M[18] - x13*M[0] + x31*M[1] + x46*M[6] + x48*M[3] + x51*M[4] + x55*M[9] - x56*M[12] + x76*M[15] + x83*M[10] + x85*M[2] + x87*M[5] + x89*M[8] + x90*M[14] + x93*M[7] + x99*M[13];
#pragma omp atomic
L[9] += x101*M[15] + x102*M[10] + x104*M[11] + x105*M[16] + x106*M[17] + x107*M[18] + x108*M[14] + x84*M[0] + x85*M[1] + x86*M[2] + x87*M[4] + x88*M[5] + x89*M[7] + x90*M[13] + x93*M[6] + x94*M[3] + x96*M[8] + x98*M[9] + x99*M[12];
#pragma omp atomic
L[10] += x41*M[1] + x42*M[2] + x55*M[7] + x61*M[0] + x68*M[3] + x73*M[4] + x74*M[5] + x81*M[6] + x98*M[8];
#pragma omp atomic
L[11] += x102*M[8] + x41*M[0] + x48*M[2] + x55*M[5] + x73*M[3] + x79*M[1] + x81*M[4] + x82*M[6] + x83*M[7];
#pragma omp atomic
L[12] += x102*M[7] + x104*M[8] + x42*M[0] + x48*M[1] + x55*M[4] + x74*M[3] + x83*M[6] + x94*M[2] + x98*M[5];
#pragma omp atomic
L[13] += x49*M[1] + x51*M[2] - x56*M[7] - x77*M[6] + x79*M[0] + x81*M[3] + x82*M[4] + x83*M[5] + x99*M[8];
#pragma omp atomic
L[14] += x102*M[5] + x48*M[0] + x51*M[1] + x55*M[3] - x56*M[6] + x83*M[4] + x87*M[2] + x90*M[8] + x99*M[7];
#pragma omp atomic
L[15] += x102*M[4] + x104*M[5] + x108*M[8] + x87*M[1] + x88*M[2] + x90*M[7] + x94*M[0] + x98*M[3] + x99*M[6];
#pragma omp atomic
L[16] += x101*M[8] + x46*M[2] + x49*M[0] - x56*M[5] + x64*M[1] + x71*M[6] + x76*M[7] - x77*M[4] + x82*M[3];
#pragma omp atomic
L[17] += x101*M[7] + x105*M[8] + x46*M[1] + x51*M[0] - x56*M[4] + x76*M[6] + x83*M[3] + x93*M[2] + x99*M[5];
#pragma omp atomic
L[18] += x101*M[6] + x102*M[3] + x105*M[7] + x106*M[8] + x87*M[0] + x89*M[2] + x90*M[5] + x93*M[1] + x99*M[4];
#pragma omp atomic
L[19] += x104*M[3] + x105*M[6] + x106*M[7] + x107*M[8] + x108*M[5] + x88*M[0] + x89*M[1] + x90*M[4] + x96*M[2];
#pragma omp atomic
L[20] += x68*M[0] + x73*M[1] + x74*M[2];
#pragma omp atomic
L[21] += x55*M[2] + x73*M[0] + x81*M[1];
#pragma omp atomic
L[22] += x55*M[1] + x74*M[0] + x98*M[2];
#pragma omp atomic
L[23] += x81*M[0] + x82*M[1] + x83*M[2];
#pragma omp atomic
L[24] += x102*M[2] + x55*M[0] + x83*M[1];
#pragma omp atomic
L[25] += x102*M[1] + x104*M[2] + x98*M[0];
#pragma omp atomic
L[26] += -x56*M[2] - x77*M[1] + x82*M[0];
#pragma omp atomic
L[27] += -x56*M[1] + x83*M[0] + x99*M[2];
#pragma omp atomic
L[28] += x102*M[0] + x90*M[2] + x99*M[1];
#pragma omp atomic
L[29] += x104*M[0] + x108*M[2] + x90*M[1];
#pragma omp atomic
L[30] += x71*M[1] + x76*M[2] - x77*M[0];
#pragma omp atomic
L[31] += x101*M[2] - x56*M[0] + x76*M[1];
#pragma omp atomic
L[32] += x101*M[1] + x105*M[2] + x99*M[0];
#pragma omp atomic
L[33] += x105*M[1] + x106*M[2] + x90*M[0];
#pragma omp atomic
L[34] += x106*M[1] + x107*M[2] + x108*M[0];

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
x19 = x17 + x6;
x20 = x19*y;
x21 = x1*x10;
x22 = -x21;
x23 = y*(x15 + x22);
x24 = x19*z;
x25 = z*(x22 + x6);
x26 = 1.0*x;
x27 = x26*(x21 - x6);
x28 = 45.0*x9;
x29 = -x28;
x30 = pow(x2, -4.5);
x31 = x0*x30;
x32 = 105.0*x31;
x33 = x29 + x32;
x34 = x*x33;
x35 = x34*y;
x36 = x34*z;
x37 = -x10;
x38 = x8*(x32 + x37);
x39 = x1*x30;
x40 = 105.0*x39;
x41 = x29 + x40;
x42 = x41*x8;
x43 = x26*x41*y;
x44 = x37 + x40;
x45 = x26*x44*z;
x46 = 315.0*x30;
x47 = 945.0*pow(x2, -5.5);
x48 = x0*x47;
x49 = x11*(x46 - x48);
x50 = x26*x8*(x1*x47 - x46);
x51 = 90.0*x9;
x52 = (x*x*x*x);
x53 = 105.0*x30;
x54 = x52*x53;
x55 = (y*y*y*y);
x56 = x53*x55;
x57 = -225.0*x9;
x58 = -x47*x52;
x59 = x*(1050.0*x31 + x57 + x58);
x60 = x29 + 630.0*x31 + x58;
x61 = x60*y;
x62 = x47*x55;
x63 = -x62;
x64 = y*(1050.0*x39 + x57 + x63);
x65 = x60*z;
x66 = 630.0*x39;
x67 = z*(x29 + x63 + x66);
x68 = x26*(x28 + x62 - x66);
x69 = x1*x32;
x70 = -x1*x48;
x71 = x*(x1*x46 + x33 + x70);
x72 = y*(x0*x46 + x41 + x70);
x73 = z*(x32 + x44 + x70);
x74 = 105.0*x9;
x75 = -12.0*x5 - x69;
x76 = 120.0*x9;
x77 = -x71;
x78 = -x72;
x79 = -x73;
#pragma omp atomic
F[0] += -x*x4*M[0] - x10*x11*M[13] + x18*M[9] + x20*M[10] + x23*M[15] + x24*M[11] + x25*M[16] - x27*M[12] + x35*M[20] + x36*M[21] + x38*M[23] - x4*y*M[1] - x4*z*M[2] + x42*M[30] + x43*M[25] + x45*M[26] + x49*M[38] - x50*M[45] + x59*M[34] + x6*x8*M[7] + x61*M[35] + x64*M[49] + x65*M[36] + x67*M[50] - x68*M[44] + x7*y*M[4] + x7*z*M[5] + x71*M[37] + x72*M[40] + x73*M[41] + (x12 + x13)*M[3] + (x12 + x14)*M[6] + (-x18 + x27)*M[14] + (-x20 - x23)*M[17] + (-x24 - x25)*M[18] + (-x35 - x43)*M[27] + (-x36 - x45)*M[28] + (-x38 - x42)*M[32] + (-x49 + x50)*M[47] + (-x59 + x77)*M[39] + (-x61 + x78)*M[42] + (-x64 + x78)*M[51] + (-x65 + x79)*M[43] + (-x67 + x79)*M[52] + (x68 + x77)*M[46] + (-x13 - x14 + 2.0*x3)*M[8] + (x19 + x22 + x69)*M[22] + (x59 - x68 + 2*x71)*M[48] + (x61 + x64 + 2*x72)*M[53] + (x65 + x67 + 2*x73)*M[54] + (-x0*x51 + x15 + x54)*M[19] + (-x1*x51 + x15 + x56)*M[29] + (x0*x74 + x21 - x54 + x75)*M[24] + (x1*x74 + x16 - x56 + x75)*M[31] + (-x0*x76 + 210.0*x1*x31 - x1*x76 + 24.0*x5 + x54 + x56)*M[33];

}

void P2M_6(double x, double y, double z, double q, double * M) {
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
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = x1*z;
x9 = (z*z);
x10 = (x*x*x);
x11 = (1.0/6.0)*q;
x12 = (1.0/2.0)*x3;
x13 = (1.0/2.0)*x0;
x14 = (y*y*y);
x15 = (1.0/2.0)*x7;
x16 = (1.0/2.0)*x9;
x17 = (z*z*z);
x18 = (x*x*x*x);
x19 = (1.0/24.0)*q;
x20 = (1.0/6.0)*x10;
x21 = q*x7;
x22 = (1.0/4.0)*x3;
x23 = q*x9;
x24 = (1.0/6.0)*x0;
x25 = (y*y*y*y);
x26 = (1.0/6.0)*x14;
x27 = (1.0/4.0)*x9;
x28 = (1.0/6.0)*x17;
x29 = (z*z*z*z);
x30 = pow(x, 5);
x31 = (1.0/120.0)*q;
x32 = (1.0/24.0)*x18;
x33 = (1.0/12.0)*x10;
x34 = (1.0/12.0)*x14;
x35 = q*x3;
x36 = x2*x7;
x37 = x1*x9;
x38 = (1.0/12.0)*x17;
x39 = (1.0/24.0)*x0;
x40 = x0*x7;
x41 = pow(y, 5);
x42 = (1.0/24.0)*x25;
x43 = (1.0/24.0)*x29;
x44 = pow(z, 5);
x45 = (1.0/720.0)*q;
x46 = (1.0/120.0)*x30;
x47 = (1.0/48.0)*x18;
x48 = (1.0/36.0)*q*x10;
x49 = (1.0/48.0)*x35;
x50 = (1.0/120.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x6;
M[6] += x4*x7;
M[7] += x8;
M[8] += x4*x9;
M[9] += -x10*x11;
M[10] += -x1*x12;
M[11] += -x12*x2;
M[12] += -x13*x7;
M[13] += -x5*z;
M[14] += -x13*x9;
M[15] += -x11*x14;
M[16] += -x15*x2;
M[17] += -x1*x16;
M[18] += -x11*x17;
M[19] += x18*x19;
M[20] += x1*x20;
M[21] += x2*x20;
M[22] += x21*x22;
M[23] += x12*x8;
M[24] += x22*x23;
M[25] += x14*x24;
M[26] += x15*x6;
M[27] += x16*x5;
M[28] += x17*x24;
M[29] += x19*x25;
M[30] += x2*x26;
M[31] += x21*x27;
M[32] += x1*x28;
M[33] += x19*x29;
M[34] += -x30*x31;
M[35] += -x1*x32;
M[36] += -x2*x32;
M[37] += -x21*x33;
M[38] += -x20*x8;
M[39] += -x23*x33;
M[40] += -x34*x35;
M[41] += -x22*x36;
M[42] += -x22*x37;
M[43] += -x35*x38;
M[44] += -x25*x39;
M[45] += -x26*x6;
M[46] += -x27*x40;
M[47] += -x28*x5;
M[48] += -x29*x39;
M[49] += -x31*x41;
M[50] += -x2*x42;
M[51] += -x23*x34;
M[52] += -x21*x38;
M[53] += -x1*x43;
M[54] += -x31*x44;
M[55] += pow(x, 6)*x45;
M[56] += x1*x46;
M[57] += x2*x46;
M[58] += x21*x47;
M[59] += x32*x8;
M[60] += x23*x47;
M[61] += x14*x48;
M[62] += x33*x36;
M[63] += x33*x37;
M[64] += x17*x48;
M[65] += x25*x49;
M[66] += x2*x3*x34;
M[67] += (1.0/8.0)*x21*x3*x9;
M[68] += x1*x3*x38;
M[69] += x29*x49;
M[70] += x41*x50;
M[71] += x42*x6;
M[72] += x0*x34*x9;
M[73] += x38*x40;
M[74] += x43*x5;
M[75] += x44*x50;
M[76] += x45*pow(y, 6);
M[77] += (1.0/120.0)*x2*x41;
M[78] += (1.0/48.0)*x23*x25;
M[79] += (1.0/36.0)*q*x14*x17;
M[80] += (1.0/48.0)*x21*x29;
M[81] += (1.0/120.0)*x1*x44;
M[82] += x45*pow(z, 6);
}
void M2M_6(double x, double y, double z, double * M, double * Ms) {
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
double x164;
double x165;
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = z*M[0];
x5 = y*M[1];
x6 = y*M[2];
x7 = z*M[1];
x8 = z*M[2];
x9 = x*M[3];
x10 = (x*x);
x11 = (1.0/2.0)*x10;
x12 = x*M[4];
x13 = y*M[3];
x14 = x0*y;
x15 = x*M[5];
x16 = z*M[3];
x17 = x0*z;
x18 = x*M[6];
x19 = y*M[4];
x20 = x1*y;
x21 = (y*y);
x22 = (1.0/2.0)*M[0];
x23 = x*M[7];
x24 = y*M[5];
x25 = z*M[4];
x26 = x3*y;
x27 = x1*z;
x28 = x2*z;
x29 = x*M[8];
x30 = z*M[5];
x31 = x3*z;
x32 = (z*z);
x33 = y*M[6];
x34 = (1.0/2.0)*x21;
x35 = y*M[7];
x36 = z*M[6];
x37 = x5*z;
x38 = y*M[8];
x39 = z*M[7];
x40 = x6*z;
x41 = (1.0/2.0)*x32;
x42 = z*M[8];
x43 = x*M[9];
x44 = (x*x*x);
x45 = (1.0/6.0)*x44;
x46 = x*M[10];
x47 = y*M[9];
x48 = x9*y;
x49 = x*M[11];
x50 = z*M[9];
x51 = x9*z;
x52 = x*M[12];
x53 = y*M[10];
x54 = x12*y;
x55 = x*M[13];
x56 = y*M[11];
x57 = z*M[10];
x58 = x15*y;
x59 = x12*z;
x60 = x13*z;
x61 = x*M[14];
x62 = z*M[11];
x63 = x15*z;
x64 = x*M[15];
x65 = y*M[12];
x66 = x18*y;
x67 = (y*y*y);
x68 = (1.0/6.0)*M[0];
x69 = x*M[16];
x70 = y*M[13];
x71 = z*M[12];
x72 = x23*y;
x73 = x18*z;
x74 = x19*z;
x75 = x*M[17];
x76 = y*M[14];
x77 = z*M[13];
x78 = x29*y;
x79 = x23*z;
x80 = x24*z;
x81 = x*M[18];
x82 = z*M[14];
x83 = x29*z;
x84 = (z*z*z);
x85 = y*M[15];
x86 = (1.0/6.0)*x67;
x87 = y*M[16];
x88 = z*M[15];
x89 = x33*z;
x90 = y*M[17];
x91 = z*M[16];
x92 = x35*z;
x93 = y*M[18];
x94 = z*M[17];
x95 = x38*z;
x96 = (1.0/6.0)*x84;
x97 = z*M[18];
x98 = x*M[19];
x99 = (1.0/24.0)*(x*x*x*x);
x100 = x*M[20];
x101 = y*M[19];
x102 = x43*y;
x103 = x*M[21];
x104 = x*M[22];
x105 = y*M[20];
x106 = x46*y;
x107 = (1.0/4.0)*x10;
x108 = x21*M[0];
x109 = x*M[23];
x110 = y*M[21];
x111 = x49*y;
x112 = x*M[24];
x113 = x107*x32;
x114 = x*M[25];
x115 = y*M[22];
x116 = x52*y;
x117 = x107*x21;
x118 = x*M[26];
x119 = y*M[23];
x120 = x55*y;
x121 = x*M[27];
x122 = y*M[24];
x123 = x61*y;
x124 = x*M[28];
x125 = x*M[29];
x126 = y*M[25];
x127 = x64*y;
x128 = (y*y*y*y);
x129 = (1.0/24.0)*M[0];
x130 = x*M[30];
x131 = y*M[26];
x132 = x69*y;
x133 = x*M[31];
x134 = y*M[27];
x135 = x75*y;
x136 = (1.0/4.0)*x32;
x137 = x*M[32];
x138 = y*M[28];
x139 = x81*y;
x140 = x*M[33];
x141 = (z*z*z*z);
x142 = y*M[29];
x143 = (1.0/24.0)*x128;
x144 = y*M[30];
x145 = y*M[31];
x146 = x136*x21;
x147 = y*M[32];
x148 = y*M[33];
x149 = (1.0/24.0)*x141;
x150 = (1.0/120.0)*pow(x, 5);
x151 = (1.0/12.0)*x44;
x152 = x151*x32;
x153 = (1.0/12.0)*x10;
x154 = x153*M[0];
x155 = x151*x21;
x156 = x153*x67;
x157 = x153*x84;
x158 = pow(y, 5);
x159 = (1.0/120.0)*M[0];
x160 = (1.0/12.0)*x32*x67;
x161 = (1.0/12.0)*x84;
x162 = pow(z, 5);
x163 = (1.0/120.0)*x158;
x164 = x161*x21;
x165 = (1.0/120.0)*x162;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + x4 + M[5];
#pragma omp atomic
Ms[6] += x5 + M[6];
#pragma omp atomic
Ms[7] += x6 + x7 + M[7];
#pragma omp atomic
Ms[8] += x8 + M[8];
#pragma omp atomic
Ms[9] += x11*M[0] + x9 + M[9];
#pragma omp atomic
Ms[10] += x11*M[1] + x12 + x13 + x14 + M[10];
#pragma omp atomic
Ms[11] += x11*M[2] + x15 + x16 + x17 + M[11];
#pragma omp atomic
Ms[12] += x18 + x19 + x20 + x21*x22 + M[12];
#pragma omp atomic
Ms[13] += x23 + x24 + x25 + x26 + x27 + x28 + M[13];
#pragma omp atomic
Ms[14] += x22*x32 + x29 + x30 + x31 + M[14];
#pragma omp atomic
Ms[15] += x33 + x34*M[1] + M[15];
#pragma omp atomic
Ms[16] += x34*M[2] + x35 + x36 + x37 + M[16];
#pragma omp atomic
Ms[17] += x38 + x39 + x40 + x41*M[1] + M[17];
#pragma omp atomic
Ms[18] += x41*M[2] + x42 + M[18];
#pragma omp atomic
Ms[19] += x11*M[3] + x43 + x45*M[0] + M[19];
#pragma omp atomic
Ms[20] += x11*x2 + x11*M[4] + x45*M[1] + x46 + x47 + x48 + M[20];
#pragma omp atomic
Ms[21] += x11*x4 + x11*M[5] + x45*M[2] + x49 + x50 + x51 + M[21];
#pragma omp atomic
Ms[22] += x0*x34 + x11*x5 + x11*M[6] + x34*M[3] + x52 + x53 + x54 + M[22];
#pragma omp atomic
Ms[23] += x11*x6 + x11*x7 + x11*M[7] + x14*z + x55 + x56 + x57 + x58 + x59 + x60 + M[23];
#pragma omp atomic
Ms[24] += x0*x41 + x11*x8 + x11*M[8] + x41*M[3] + x61 + x62 + x63 + M[24];
#pragma omp atomic
Ms[25] += x1*x34 + x34*M[4] + x64 + x65 + x66 + x67*x68 + M[25];
#pragma omp atomic
Ms[26] += x20*z + x3*x34 + x34*x4 + x34*M[5] + x69 + x70 + x71 + x72 + x73 + x74 + M[26];
#pragma omp atomic
Ms[27] += x1*x41 + x2*x41 + x26*z + x41*M[4] + x75 + x76 + x77 + x78 + x79 + x80 + M[27];
#pragma omp atomic
Ms[28] += x3*x41 + x41*M[5] + x68*x84 + x81 + x82 + x83 + M[28];
#pragma omp atomic
Ms[29] += x34*M[6] + x85 + x86*M[1] + M[29];
#pragma omp atomic
Ms[30] += x34*x7 + x34*M[7] + x86*M[2] + x87 + x88 + x89 + M[30];
#pragma omp atomic
Ms[31] += x34*x8 + x34*M[8] + x41*x5 + x41*M[6] + x90 + x91 + x92 + M[31];
#pragma omp atomic
Ms[32] += x41*x6 + x41*M[7] + x93 + x94 + x95 + x96*M[1] + M[32];
#pragma omp atomic
Ms[33] += x41*M[8] + x96*M[2] + x97 + M[33];
#pragma omp atomic
Ms[34] += x11*M[9] + x45*M[3] + x98 + x99*M[0] + M[34];
#pragma omp atomic
Ms[35] += x100 + x101 + x102 + x11*x13 + x11*M[10] + x2*x45 + x45*M[4] + x99*M[1] + M[35];
#pragma omp atomic
Ms[36] += x103 + x11*x16 + x11*M[11] + x4*x45 + x43*z + x45*M[5] + x99*M[2] + z*M[19] + M[36];
#pragma omp atomic
Ms[37] += x104 + x105 + x106 + x107*x108 + x11*x19 + x11*M[12] + x34*x9 + x34*M[9] + x45*x5 + x45*M[6] + M[37];
#pragma omp atomic
Ms[38] += x109 + x11*x24 + x11*x25 + x11*x28 + x11*M[13] + x110 + x111 + x45*x6 + x45*x7 + x45*M[7] + x46*z + x47*z + x48*z + z*M[20] + M[38];
#pragma omp atomic
Ms[39] += x11*x30 + x11*M[14] + x112 + x113*M[0] + x41*x9 + x41*M[9] + x45*x8 + x45*M[8] + x49*z + z*M[21] + M[39];
#pragma omp atomic
Ms[40] += x0*x86 + x11*x33 + x11*M[15] + x114 + x115 + x116 + x117*M[1] + x12*x34 + x34*M[10] + x86*M[3] + M[40];
#pragma omp atomic
Ms[41] += x11*x35 + x11*x36 + x11*x37 + x11*M[16] + x117*M[2] + x118 + x119 + x120 + x15*x34 + x16*x34 + x17*x34 + x34*M[11] + x52*z + x53*z + x54*z + z*M[22] + M[41];
#pragma omp atomic
Ms[42] += x11*x38 + x11*x39 + x11*x40 + x11*M[17] + x113*M[1] + x12*x41 + x121 + x122 + x123 + x13*x41 + x14*x41 + x41*M[10] + x55*z + x56*z + x58*z + z*M[23] + M[42];
#pragma omp atomic
Ms[43] += x0*x96 + x11*x42 + x11*M[18] + x113*M[2] + x124 + x15*x41 + x41*M[11] + x61*z + x96*M[3] + z*M[24] + M[43];
#pragma omp atomic
Ms[44] += x1*x86 + x125 + x126 + x127 + x128*x129 + x18*x34 + x34*M[12] + x86*M[4] + M[44];
#pragma omp atomic
Ms[45] += x130 + x131 + x132 + x23*x34 + x25*x34 + x27*x34 + x3*x86 + x34*M[13] + x4*x86 + x64*z + x65*z + x66*z + x86*M[5] + z*M[25] + M[45];
#pragma omp atomic
Ms[46] += x108*x136 + x133 + x134 + x135 + x18*x41 + x19*x41 + x20*x41 + x29*x34 + x30*x34 + x31*x34 + x34*M[14] + x41*M[12] + x69*z + x70*z + x72*z + z*M[26] + M[46];
#pragma omp atomic
Ms[47] += x1*x96 + x137 + x138 + x139 + x2*x96 + x23*x41 + x24*x41 + x26*x41 + x41*M[13] + x75*z + x76*z + x78*z + x96*M[4] + z*M[27] + M[47];
#pragma omp atomic
Ms[48] += x129*x141 + x140 + x29*x41 + x3*x96 + x41*M[14] + x81*z + x96*M[5] + z*M[28] + M[48];
#pragma omp atomic
Ms[49] += x142 + x143*M[1] + x34*M[15] + x86*M[6] + M[49];
#pragma omp atomic
Ms[50] += x143*M[2] + x144 + x34*x36 + x34*M[16] + x7*x86 + x85*z + x86*M[7] + z*M[29] + M[50];
#pragma omp atomic
Ms[51] += x145 + x146*M[1] + x33*x41 + x34*x39 + x34*M[17] + x41*M[15] + x8*x86 + x86*M[8] + x87*z + z*M[30] + M[51];
#pragma omp atomic
Ms[52] += x146*M[2] + x147 + x34*x42 + x34*M[18] + x35*x41 + x41*M[16] + x5*x96 + x90*z + x96*M[6] + z*M[31] + M[52];
#pragma omp atomic
Ms[53] += x148 + x149*M[1] + x38*x41 + x41*M[17] + x6*x96 + x93*z + x96*M[7] + z*M[32] + M[53];
#pragma omp atomic
Ms[54] += x149*M[2] + x41*M[18] + x96*M[8] + z*M[33] + M[54];
#pragma omp atomic
Ms[55] += x*M[34] + x11*M[19] + x150*M[0] + x45*M[9] + x99*M[3] + M[55];
#pragma omp atomic
Ms[56] += x*M[35] + x11*x47 + x11*M[20] + x13*x45 + x150*M[1] + x2*x99 + x45*M[10] + x98*y + x99*M[4] + y*M[34] + M[56];
#pragma omp atomic
Ms[57] += x*M[36] + x11*x50 + x11*M[21] + x150*M[2] + x16*x45 + x4*x99 + x45*M[11] + x98*z + x99*M[5] + z*M[34] + M[57];
#pragma omp atomic
Ms[58] += x*M[37] + x100*y + x108*x151 + x11*x53 + x11*M[22] + x117*M[3] + x19*x45 + x34*x43 + x34*M[19] + x45*M[12] + x5*x99 + x99*M[6] + y*M[35] + M[58];
#pragma omp atomic
Ms[59] += x*M[38] + x100*z + x101*z + x102*z + x103*y + x11*x56 + x11*x57 + x11*x60 + x11*M[23] + x24*x45 + x25*x45 + x28*x45 + x45*M[13] + x6*x99 + x7*x99 + x99*M[7] + y*M[36] + z*M[35] + M[59];
#pragma omp atomic
Ms[60] += x*M[39] + x103*z + x11*x62 + x11*M[24] + x113*M[3] + x152*M[0] + x30*x45 + x41*x43 + x41*M[19] + x45*M[14] + x8*x99 + x99*M[8] + z*M[36] + M[60];
#pragma omp atomic
Ms[61] += x*M[40] + x104*y + x11*x65 + x11*M[25] + x117*M[4] + x154*x67 + x155*M[1] + x33*x45 + x34*x46 + x34*M[20] + x45*M[15] + x86*x9 + x86*M[9] + y*M[37] + M[61];
#pragma omp atomic
Ms[62] += x*M[41] + x104*z + x105*z + x106*z + x109*y + x11*x70 + x11*x71 + x11*x74 + x11*M[26] + x117*x4 + x117*M[5] + x155*M[2] + x34*x49 + x34*x50 + x34*x51 + x34*M[21] + x35*x45 + x36*x45 + x37*x45 + x45*M[16] + y*M[38] + z*M[37] + M[62];
#pragma omp atomic
Ms[63] += x*M[42] + x109*z + x11*x76 + x11*x77 + x11*x80 + x11*M[27] + x110*z + x111*z + x112*y + x113*x2 + x113*M[4] + x152*M[1] + x38*x45 + x39*x45 + x40*x45 + x41*x46 + x41*x47 + x41*x48 + x41*M[20] + x45*M[17] + y*M[39] + z*M[38] + M[63];
#pragma omp atomic
Ms[64] += x*M[43] + x11*x82 + x11*M[28] + x112*z + x113*M[5] + x152*M[2] + x154*x84 + x41*x49 + x41*M[21] + x42*x45 + x45*M[18] + x9*x96 + x96*M[9] + z*M[39] + M[64];
#pragma omp atomic
Ms[65] += x*M[44] + x0*x143 + x11*x85 + x11*M[29] + x114*y + x117*M[6] + x12*x86 + x143*M[3] + x156*M[1] + x34*x52 + x34*M[22] + x86*M[10] + y*M[40] + M[65];
#pragma omp atomic
Ms[66] += x*M[45] + x11*x87 + x11*x88 + x11*x89 + x11*M[30] + x114*z + x115*z + x116*z + x117*x7 + x117*M[7] + x118*y + x15*x86 + x156*M[2] + x16*x86 + x17*x86 + x34*x55 + x34*x57 + x34*x59 + x34*M[23] + x86*M[11] + y*M[41] + z*M[40] + M[66];
#pragma omp atomic
Ms[67] += x*M[46] + x0*x146 + x11*x90 + x11*x91 + x11*x92 + x11*M[31] + x113*x5 + x113*M[6] + x117*x8 + x117*M[8] + x118*z + x119*z + x120*z + x121*y + x146*M[3] + x34*x61 + x34*x62 + x34*x63 + x34*M[24] + x41*x52 + x41*x53 + x41*x54 + x41*M[22] + y*M[42] + z*M[41] + M[67];
#pragma omp atomic
Ms[68] += x*M[47] + x11*x93 + x11*x94 + x11*x95 + x11*M[32] + x113*x6 + x113*M[7] + x12*x96 + x121*z + x122*z + x123*z + x124*y + x13*x96 + x14*x96 + x157*M[1] + x41*x55 + x41*x56 + x41*x58 + x41*M[23] + x96*M[10] + y*M[43] + z*M[42] + M[68];
#pragma omp atomic
Ms[69] += x*M[48] + x0*x149 + x11*x97 + x11*M[33] + x113*M[8] + x124*z + x149*M[3] + x15*x96 + x157*M[2] + x41*x61 + x41*M[24] + x96*M[11] + z*M[43] + M[69];
#pragma omp atomic
Ms[70] += x*M[49] + x1*x143 + x125*y + x143*M[4] + x158*x159 + x18*x86 + x34*x64 + x34*M[25] + x86*M[12] + y*M[44] + M[70];
#pragma omp atomic
Ms[71] += x*M[50] + x125*z + x126*z + x127*z + x130*y + x143*x3 + x143*x4 + x143*M[5] + x23*x86 + x25*x86 + x27*x86 + x34*x69 + x34*x71 + x34*x73 + x34*M[26] + x86*M[13] + y*M[45] + z*M[44] + M[71];
#pragma omp atomic
Ms[72] += x*M[51] + x1*x146 + x130*z + x131*z + x132*z + x133*y + x146*M[4] + x160*M[0] + x29*x86 + x30*x86 + x31*x86 + x34*x75 + x34*x77 + x34*x79 + x34*M[27] + x41*x64 + x41*x65 + x41*x66 + x41*M[25] + x86*M[14] + y*M[46] + z*M[45] + M[72];
#pragma omp atomic
Ms[73] += x*M[52] + x108*x161 + x133*z + x134*z + x135*z + x137*y + x146*x3 + x146*M[5] + x18*x96 + x19*x96 + x20*x96 + x34*x81 + x34*x82 + x34*x83 + x34*M[28] + x41*x69 + x41*x70 + x41*x72 + x41*M[26] + x96*M[12] + y*M[47] + z*M[46] + M[73];
#pragma omp atomic
Ms[74] += x*M[53] + x1*x149 + x137*z + x138*z + x139*z + x140*y + x149*x2 + x149*M[4] + x23*x96 + x24*x96 + x26*x96 + x41*x75 + x41*x76 + x41*x78 + x41*M[27] + x96*M[13] + y*M[48] + z*M[47] + M[74];
#pragma omp atomic
Ms[75] += x*M[54] + x140*z + x149*x3 + x149*M[5] + x159*x162 + x29*x96 + x41*x81 + x41*M[28] + x96*M[14] + z*M[48] + M[75];
#pragma omp atomic
Ms[76] += x143*M[6] + x163*M[1] + x34*M[29] + x86*M[15] + y*M[49] + M[76];
#pragma omp atomic
Ms[77] += x142*z + x143*x7 + x143*M[7] + x163*M[2] + x34*x88 + x34*M[30] + x36*x86 + x86*M[16] + y*M[50] + z*M[49] + M[77];
#pragma omp atomic
Ms[78] += x143*x8 + x143*M[8] + x144*z + x146*M[6] + x160*M[1] + x34*x91 + x34*M[31] + x39*x86 + x41*x85 + x41*M[29] + x86*M[17] + y*M[51] + z*M[50] + M[78];
#pragma omp atomic
Ms[79] += x145*z + x146*M[7] + x160*M[2] + x164*M[1] + x33*x96 + x34*x94 + x34*M[32] + x41*x87 + x41*M[30] + x42*x86 + x86*M[18] + x96*M[15] + y*M[52] + z*M[51] + M[79];
#pragma omp atomic
Ms[80] += x146*M[8] + x147*z + x149*x5 + x149*M[6] + x164*M[2] + x34*x97 + x34*M[33] + x35*x96 + x41*x90 + x41*M[31] + x96*M[16] + y*M[53] + z*M[52] + M[80];
#pragma omp atomic
Ms[81] += x148*z + x149*x6 + x149*M[7] + x165*M[1] + x38*x96 + x41*x93 + x41*M[32] + x96*M[17] + y*M[54] + z*M[53] + M[81];
#pragma omp atomic
Ms[82] += x149*M[8] + x165*M[2] + x41*M[33] + x96*M[18] + z*M[54] + M[82];

}

void M2L_6(double x, double y, double z, double * M, double * L) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = pow(R, -5);
x3 = 3.0*x2;
x4 = x*y;
x5 = x3*x4;
x6 = x*z;
x7 = x3*x6;
x8 = y*z;
x9 = x3*x8;
x10 = pow(R, -7);
x11 = 15.0*x10;
x12 = x4*z;
x13 = x11*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x3;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x3;
x20 = x14 + x19;
x21 = 9.0*x2;
x22 = x11*x15;
x23 = -x22;
x24 = x*(x21 + x23);
x25 = x23 + x3;
x26 = x25*y;
x27 = x11*x18;
x28 = -x27;
x29 = y*(x21 + x28);
x30 = x25*z;
x31 = z*(x28 + x3);
x32 = 1.0*x;
x33 = x32*(x27 - x3);
x34 = 45.0*x10;
x35 = -x34;
x36 = pow(R, -9);
x37 = x15*x36;
x38 = 105.0*x37;
x39 = x35 + x38;
x40 = x39*x4;
x41 = x39*x6;
x42 = x18*x36;
x43 = 105.0*x42;
x44 = x35 + x43;
x45 = x44*x8;
x46 = -x11;
x47 = x8*(x38 + x46);
x48 = x32*y;
x49 = x44*x48;
x50 = x43 + x46;
x51 = x32*z;
x52 = x50*x51;
x53 = 315.0*x36;
x54 = pow(R, -11);
x55 = 945.0*x54;
x56 = x15*x55;
x57 = x53 - x56;
x58 = x12*x57;
x59 = x18*x55;
x60 = x32*x8*(-x53 + x59);
x61 = (x*x*x*x);
x62 = 105.0*x36;
x63 = x61*x62;
x64 = 90.0*x10;
x65 = -x15*x64 + x21 + x63;
x66 = (y*y*y*y);
x67 = x62*x66;
x68 = -x18*x64 + x21 + x67;
x69 = 2.0*x0 - x16 - x19;
x70 = -225.0*x10;
x71 = x55*x61;
x72 = -x71;
x73 = x*(1050.0*x37 + x70 + x72);
x74 = x55*x66;
x75 = -x74;
x76 = y*(1050.0*x42 + x70 + x75);
x77 = x35 + 630.0*x37 + x72;
x78 = x77*y;
x79 = x77*z;
x80 = 630.0*x42;
x81 = x35 + x75 + x80;
x82 = x81*z;
x83 = x32*(x34 + x74 - x80);
x84 = 1575.0*x36;
x85 = pow(R, -13);
x86 = 10395.0*x85;
x87 = x61*x86;
x88 = x15*x54;
x89 = x84 + x87 - 9450.0*x88;
x90 = x4*x89;
x91 = x6*x89;
x92 = 5670.0*x88;
x93 = x8*(x53 + x87 - x92);
x94 = x66*x86;
x95 = x18*x54;
x96 = x84 + x94 - 9450.0*x95;
x97 = x8*x96;
x98 = x48*x96;
x99 = x51*(x53 + x94 - 5670.0*x95);
x100 = pow(x, 6)*x86;
x101 = 14175.0*x54;
x102 = x100 - x101*x61 + 4725.0*x37 + x70;
x103 = x86*pow(y, 6);
x104 = -x101*x66 + x103 + 4725.0*x42 + x70;
x105 = x18*x38;
x106 = x105 + x25 + x28;
x107 = x18*x53;
x108 = -x18*x56;
x109 = x*(x107 + x108 + x39);
x110 = x15*x53;
x111 = y*(x108 + x110 + x44);
x112 = z*(x108 + x38 + x50);
x113 = -x24 + x33;
x114 = -x26 - x29;
x115 = -x30 - x31;
x116 = 945.0*x36;
x117 = -2835.0*x88;
x118 = x15*x18*x86;
x119 = x118 - 2835.0*x95;
x120 = x4*(x116 + x117 + x119);
x121 = x6*(x119 + x57);
x122 = x8*(x117 + x118 + x53 - x59);
x123 = -x40 - x49;
x124 = -x41 - x52;
x125 = -x45 - x47;
x126 = -x58 + x60;
x127 = 105.0*x10;
x128 = -x105 - 12.0*x2;
x129 = x127*x18 + x128 + x22 - x67;
x130 = x127*x15 + x128 + x27 - x63;
x131 = 120.0*x10;
x132 = -x131*x15 - x131*x18 + 210.0*x18*x37 + 24.0*x2 + x63 + x67;
x133 = x15*x94;
x134 = x18*x92;
x135 = -x134;
x136 = x110 + x133 + x135 + x81;
x137 = x18*x87;
x138 = x107 + x135 + x137 + x77;
x139 = 15120.0*x54;
x140 = -x100;
x141 = -x137;
x142 = 270.0*x10 + x134;
x143 = -x107 + x139*x61 + x140 + x141 + x142 - 5355.0*x37;
x144 = -x103;
x145 = -x133;
x146 = -x110 + x139*x66 + x142 + x144 + x145 - 5355.0*x42;
x147 = -x109;
x148 = x147 - x73;
x149 = x147 + x83;
x150 = -x111;
x151 = x150 - x76;
x152 = x150 - x78;
x153 = -x112;
x154 = x153 - x79;
x155 = x153 - x82;
x156 = -x120;
x157 = x156 - x90;
x158 = x156 - x98;
x159 = -x121;
x160 = x159 - x99;
x161 = x159 - x91;
x162 = -x122;
x163 = x162 - x93;
x164 = x162 - x97;
x165 = x18*x88;
x166 = -x116*x15 - x116*x18 + x141 + x145 + 11340.0*x165 + x64 + x71 + x74;
x167 = 16065.0*x54;
x168 = 20790.0*x85;
x169 = x15*x66;
x170 = -360.0*x10 - 17010.0*x165;
x171 = x103 + x137 - x167*x66 + x168*x169 + x170 + 1260.0*x37 + 6300.0*x42 + x72;
x172 = x18*x61;
x173 = x100 + x133 - x167*x61 + x168*x172 + x170 + 6300.0*x37 + 1260.0*x42 + x75;
x174 = 2*x111 + x76 + x78;
x175 = 2*x112 + x79 + x82;
x176 = 2*x109 + x73 - x83;
x177 = 17010.0*x54;
x178 = 31185.0*x85;
x179 = 720.0*x10 + x140 + x144 + 34020.0*x165 - x169*x178 - x172*x178 + x177*x61 + x177*x66 - 7560.0*x37 - 7560.0*x42;
x180 = 2*x122 + x93 + x97;
x181 = 2*x120 + x90 + x98;
x182 = 2*x121 + x91 + x99;
#pragma omp atomic
L[0] += -x*x1*M[0] - x1*y*M[1] - x1*z*M[2] + x102*M[55] + x104*M[76] + x106*M[22] + x109*M[37] + x111*M[40] + x112*M[41] + x113*M[14] + x114*M[17] + x115*M[18] + x120*M[61] + x121*M[62] + x122*M[66] + x123*M[27] + x124*M[28] + x125*M[32] + x126*M[47] + x129*M[31] - x13*M[13] + x130*M[24] + x132*M[33] + x136*M[65] + x138*M[58] + x143*M[60] + x146*M[78] + x148*M[39] + x149*M[46] + x151*M[51] + x152*M[42] + x154*M[43] + x155*M[52] + x157*M[63] + x158*M[72] + x160*M[73] + x161*M[64] + x163*M[68] + x164*M[79] + x166*M[67] + x17*M[3] + x171*M[80] + x173*M[69] + x174*M[53] + x175*M[54] + x176*M[48] + x179*M[82] + x180*M[81] + x181*M[74] + x182*M[75] + x20*M[6] + x24*M[9] + x26*M[10] + x29*M[15] + x30*M[11] + x31*M[16] - x33*M[12] + x40*M[20] + x41*M[21] + x45*M[30] + x47*M[23] + x49*M[25] + x5*M[4] + x52*M[26] + x58*M[38] - x60*M[45] + x65*M[19] + x68*M[29] + x69*M[8] + x7*M[5] + x73*M[34] + x76*M[49] + x78*M[35] + x79*M[36] + x82*M[50] - x83*M[44] + x9*M[7] + x90*M[56] + x91*M[57] + x93*M[59] + x97*M[77] + x98*M[70] + x99*M[71];
#pragma omp atomic
L[1] += x102*M[34] + x106*M[12] + x109*M[22] + x111*M[25] + x112*M[26] + x113*M[8] + x120*M[40] + x121*M[41] + x122*M[45] + x123*M[17] + x124*M[18] + x126*M[32] - x13*M[7] + x130*M[14] + x136*M[44] + x138*M[37] + x143*M[39] + x148*M[24] + x149*M[31] + x152*M[27] + x154*M[28] + x157*M[42] + x158*M[51] + x160*M[52] + x161*M[43] + x163*M[47] + x166*M[46] + x17*M[0] + x173*M[48] + x176*M[33] + x181*M[53] + x182*M[54] + x24*M[3] + x26*M[4] + x30*M[5] - x33*M[6] + x40*M[10] + x41*M[11] + x47*M[13] + x49*M[15] + x5*M[1] + x52*M[16] + x58*M[23] - x60*M[30] + x65*M[9] + x7*M[2] + x73*M[19] + x78*M[20] + x79*M[21] - x83*M[29] + x90*M[35] + x91*M[36] + x93*M[38] + x98*M[49] + x99*M[50];
#pragma omp atomic
L[2] += x104*M[49] + x106*M[10] + x109*M[20] + x111*M[22] + x112*M[23] + x114*M[8] + x120*M[37] + x121*M[38] + x122*M[41] + x123*M[14] + x125*M[18] + x126*M[28] + x129*M[17] - x13*M[5] + x136*M[40] + x138*M[35] + x146*M[51] + x149*M[27] + x151*M[31] + x152*M[24] + x155*M[32] + x157*M[39] + x158*M[46] + x160*M[47] + x163*M[43] + x164*M[52] + x166*M[42] + x171*M[53] + x174*M[33] + x180*M[54] + x181*M[48] + x20*M[1] + x26*M[3] + x29*M[6] + x31*M[7] - x33*M[4] + x40*M[9] + x45*M[16] + x47*M[11] + x49*M[12] + x5*M[0] + x52*M[13] + x58*M[21] - x60*M[26] + x68*M[15] + x76*M[29] + x78*M[19] + x82*M[30] - x83*M[25] + x9*M[2] + x90*M[34] + x93*M[36] + x97*M[50] + x98*M[44] + x99*M[45];
#pragma omp atomic
L[3] += x112*M[22] + x113*M[5] + x114*M[7] + x115*M[8] + x121*M[37] + x122*M[40] + x123*M[13] + x124*M[14] + x125*M[17] + x126*M[27] + x129*M[16] - x13*M[4] + x130*M[11] + x132*M[18] + x143*M[36] + x146*M[50] + x148*M[21] + x149*M[26] + x151*M[30] + x152*M[23] + x154*M[24] + x155*M[31] + x157*M[38] + x158*M[45] + x160*M[46] + x161*M[39] + x163*M[42] + x164*M[51] + x166*M[41] + x171*M[52] + x173*M[43] + x174*M[32] + x175*M[33] + x176*M[28] + x179*M[54] + x180*M[53] + x181*M[47] + x182*M[48] + x30*M[3] + x31*M[6] + x41*M[9] + x45*M[15] + x47*M[10] + x52*M[12] + x58*M[20] - x60*M[25] + x69*M[2] + x7*M[0] + x79*M[19] + x82*M[29] + x9*M[1] + x91*M[34] + x93*M[35] + x97*M[49] + x99*M[44];
#pragma omp atomic
L[4] += x102*M[19] + x106*M[6] + x109*M[12] + x111*M[15] + x112*M[16] + x120*M[25] + x121*M[26] + x122*M[30] + x130*M[8] + x136*M[29] + x138*M[22] + x143*M[24] + x148*M[14] + x152*M[17] + x154*M[18] + x157*M[27] + x161*M[28] + x163*M[32] + x166*M[31] + x173*M[33] + x24*M[0] + x26*M[1] + x30*M[2] + x40*M[4] + x41*M[5] + x47*M[7] + x58*M[13] + x65*M[3] + x73*M[9] + x78*M[10] + x79*M[11] + x90*M[20] + x91*M[21] + x93*M[23];
#pragma omp atomic
L[5] += x106*M[4] + x109*M[10] + x111*M[12] + x112*M[13] + x120*M[22] + x121*M[23] + x122*M[26] + x123*M[8] + x126*M[18] - x13*M[2] + x136*M[25] + x138*M[20] + x149*M[17] + x152*M[14] + x157*M[24] + x158*M[31] + x160*M[32] + x163*M[28] + x166*M[27] + x181*M[33] + x26*M[0] - x33*M[1] + x40*M[3] + x47*M[5] + x49*M[6] + x52*M[7] + x58*M[11] - x60*M[16] + x78*M[9] - x83*M[15] + x90*M[19] + x93*M[21] + x98*M[29] + x99*M[30];
#pragma omp atomic
L[6] += x112*M[12] + x113*M[2] + x121*M[22] + x122*M[25] + x123*M[7] + x124*M[8] + x126*M[17] - x13*M[1] + x130*M[5] + x143*M[21] + x148*M[11] + x149*M[16] + x152*M[13] + x154*M[14] + x157*M[23] + x158*M[30] + x160*M[31] + x161*M[24] + x163*M[27] + x166*M[26] + x173*M[28] + x176*M[18] + x181*M[32] + x182*M[33] + x30*M[0] + x41*M[3] + x47*M[4] + x52*M[6] + x58*M[10] - x60*M[15] + x79*M[9] + x91*M[19] + x93*M[20] + x99*M[29];
#pragma omp atomic
L[7] += x104*M[29] + x106*M[3] + x109*M[9] + x111*M[10] + x112*M[11] + x120*M[20] + x121*M[21] + x122*M[23] + x129*M[8] + x136*M[22] + x138*M[19] + x146*M[31] + x149*M[14] + x151*M[17] + x155*M[18] + x158*M[27] + x160*M[28] + x164*M[32] + x166*M[24] + x171*M[33] + x29*M[1] + x31*M[2] - x33*M[0] + x45*M[7] + x49*M[4] + x52*M[5] - x60*M[13] + x68*M[6] + x76*M[15] + x82*M[16] - x83*M[12] + x97*M[30] + x98*M[25] + x99*M[26];
#pragma omp atomic
L[8] += x112*M[10] + x114*M[2] + x121*M[20] + x122*M[22] + x123*M[5] + x125*M[8] + x126*M[14] + x129*M[7] - x13*M[0] + x146*M[30] + x149*M[13] + x151*M[16] + x152*M[11] + x155*M[17] + x157*M[21] + x158*M[26] + x160*M[27] + x163*M[24] + x164*M[31] + x166*M[23] + x171*M[32] + x174*M[18] + x180*M[33] + x181*M[28] + x31*M[1] + x45*M[6] + x47*M[3] + x52*M[4] + x58*M[9] - x60*M[12] + x82*M[15] + x93*M[19] + x97*M[29] + x99*M[25];
#pragma omp atomic
L[9] += x113*M[0] + x114*M[1] + x115*M[2] + x123*M[4] + x124*M[5] + x125*M[7] + x126*M[13] + x129*M[6] + x130*M[3] + x132*M[8] + x143*M[19] + x146*M[29] + x148*M[9] + x149*M[12] + x151*M[15] + x152*M[10] + x154*M[11] + x155*M[16] + x157*M[20] + x158*M[25] + x160*M[26] + x161*M[21] + x163*M[23] + x164*M[30] + x166*M[22] + x171*M[31] + x173*M[24] + x174*M[17] + x175*M[18] + x176*M[14] + x179*M[33] + x180*M[32] + x181*M[27] + x182*M[28];
#pragma omp atomic
L[10] += x102*M[9] + x109*M[6] + x120*M[15] + x121*M[16] + x138*M[12] + x143*M[14] + x148*M[8] + x157*M[17] + x161*M[18] + x40*M[1] + x41*M[2] + x58*M[7] + x65*M[0] + x73*M[3] + x78*M[4] + x79*M[5] + x90*M[10] + x91*M[11] + x93*M[13];
#pragma omp atomic
L[11] += x106*M[1] + x109*M[4] + x111*M[6] + x112*M[7] + x120*M[12] + x121*M[13] + x122*M[16] + x136*M[15] + x138*M[10] + x152*M[8] + x157*M[14] + x163*M[18] + x166*M[17] + x40*M[0] + x47*M[2] + x58*M[5] + x78*M[3] + x90*M[9] + x93*M[11];
#pragma omp atomic
L[12] += x112*M[6] + x121*M[12] + x122*M[15] + x130*M[2] + x143*M[11] + x148*M[5] + x152*M[7] + x154*M[8] + x157*M[13] + x161*M[14] + x163*M[17] + x166*M[16] + x173*M[18] + x41*M[0] + x47*M[1] + x58*M[4] + x79*M[3] + x91*M[9] + x93*M[10];
#pragma omp atomic
L[13] += x106*M[0] + x109*M[3] + x111*M[4] + x112*M[5] + x120*M[10] + x121*M[11] + x122*M[13] + x136*M[12] + x138*M[9] + x149*M[8] + x158*M[17] + x160*M[18] + x166*M[14] + x49*M[1] + x52*M[2] - x60*M[7] - x83*M[6] + x98*M[15] + x99*M[16];
#pragma omp atomic
L[14] += x112*M[4] + x121*M[10] + x122*M[12] + x123*M[2] + x126*M[8] + x149*M[7] + x152*M[5] + x157*M[11] + x158*M[16] + x160*M[17] + x163*M[14] + x166*M[13] + x181*M[18] + x47*M[0] + x52*M[1] + x58*M[3] - x60*M[6] + x93*M[9] + x99*M[15];
#pragma omp atomic
L[15] += x123*M[1] + x124*M[2] + x126*M[7] + x130*M[0] + x143*M[9] + x148*M[3] + x149*M[6] + x152*M[4] + x154*M[5] + x157*M[10] + x158*M[15] + x160*M[16] + x161*M[11] + x163*M[13] + x166*M[12] + x173*M[14] + x176*M[8] + x181*M[17] + x182*M[18];
#pragma omp atomic
L[16] += x104*M[15] + x111*M[3] + x120*M[9] + x122*M[11] + x136*M[10] + x146*M[17] + x151*M[8] + x158*M[14] + x164*M[18] + x45*M[2] + x49*M[0] - x60*M[5] + x68*M[1] + x76*M[6] + x82*M[7] - x83*M[4] + x97*M[16] + x98*M[12] + x99*M[13];
#pragma omp atomic
L[17] += x112*M[3] + x121*M[9] + x122*M[10] + x129*M[2] + x146*M[16] + x149*M[5] + x151*M[7] + x155*M[8] + x158*M[13] + x160*M[14] + x164*M[17] + x166*M[11] + x171*M[18] + x45*M[1] + x52*M[0] - x60*M[4] + x82*M[6] + x97*M[15] + x99*M[12];
#pragma omp atomic
L[18] += x123*M[0] + x125*M[2] + x126*M[5] + x129*M[1] + x146*M[15] + x149*M[4] + x151*M[6] + x152*M[3] + x155*M[7] + x157*M[9] + x158*M[12] + x160*M[13] + x163*M[11] + x164*M[16] + x166*M[10] + x171*M[17] + x174*M[8] + x180*M[18] + x181*M[14];
#pragma omp atomic
L[19] += x124*M[0] + x125*M[1] + x126*M[4] + x132*M[2] + x154*M[3] + x155*M[6] + x160*M[12] + x161*M[9] + x163*M[10] + x164*M[15] + x171*M[16] + x173*M[11] + x174*M[7] + x175*M[8] + x176*M[5] + x179*M[18] + x180*M[17] + x181*M[13] + x182*M[14];
#pragma omp atomic
L[20] += x102*M[3] + x138*M[6] + x143*M[8] + x73*M[0] + x78*M[1] + x79*M[2] + x90*M[4] + x91*M[5] + x93*M[7];
#pragma omp atomic
L[21] += x109*M[1] + x120*M[6] + x121*M[7] + x138*M[4] + x157*M[8] + x58*M[2] + x78*M[0] + x90*M[3] + x93*M[5];
#pragma omp atomic
L[22] += x121*M[6] + x143*M[5] + x148*M[2] + x157*M[7] + x161*M[8] + x58*M[1] + x79*M[0] + x91*M[3] + x93*M[4];
#pragma omp atomic
L[23] += x109*M[0] + x111*M[1] + x112*M[2] + x120*M[4] + x121*M[5] + x122*M[7] + x136*M[6] + x138*M[3] + x166*M[8];
#pragma omp atomic
L[24] += x112*M[1] + x121*M[4] + x122*M[6] + x152*M[2] + x157*M[5] + x163*M[8] + x166*M[7] + x58*M[0] + x93*M[3];
#pragma omp atomic
L[25] += x143*M[3] + x148*M[0] + x152*M[1] + x154*M[2] + x157*M[4] + x161*M[5] + x163*M[7] + x166*M[6] + x173*M[8];
#pragma omp atomic
L[26] += x111*M[0] + x120*M[3] + x122*M[5] + x136*M[4] + x158*M[8] - x60*M[2] - x83*M[1] + x98*M[6] + x99*M[7];
#pragma omp atomic
L[27] += x112*M[0] + x121*M[3] + x122*M[4] + x149*M[2] + x158*M[7] + x160*M[8] + x166*M[5] - x60*M[1] + x99*M[6];
#pragma omp atomic
L[28] += x126*M[2] + x149*M[1] + x152*M[0] + x157*M[3] + x158*M[6] + x160*M[7] + x163*M[5] + x166*M[4] + x181*M[8];
#pragma omp atomic
L[29] += x126*M[1] + x154*M[0] + x160*M[6] + x161*M[3] + x163*M[4] + x173*M[5] + x176*M[2] + x181*M[7] + x182*M[8];
#pragma omp atomic
L[30] += x104*M[6] + x136*M[3] + x146*M[8] + x76*M[1] + x82*M[2] - x83*M[0] + x97*M[7] + x98*M[4] + x99*M[5];
#pragma omp atomic
L[31] += x122*M[3] + x146*M[7] + x151*M[2] + x158*M[5] + x164*M[8] - x60*M[0] + x82*M[1] + x97*M[6] + x99*M[4];
#pragma omp atomic
L[32] += x146*M[6] + x149*M[0] + x151*M[1] + x155*M[2] + x158*M[4] + x160*M[5] + x164*M[7] + x166*M[3] + x171*M[8];
#pragma omp atomic
L[33] += x126*M[0] + x155*M[1] + x160*M[4] + x163*M[3] + x164*M[6] + x171*M[7] + x174*M[2] + x180*M[8] + x181*M[5];
#pragma omp atomic
L[34] += x171*M[6] + x173*M[3] + x174*M[1] + x175*M[2] + x176*M[0] + x179*M[8] + x180*M[7] + x181*M[4] + x182*M[5];
#pragma omp atomic
L[35] += x102*M[0] + x90*M[1] + x91*M[2];
#pragma omp atomic
L[36] += x138*M[1] + x90*M[0] + x93*M[2];
#pragma omp atomic
L[37] += x143*M[2] + x91*M[0] + x93*M[1];
#pragma omp atomic
L[38] += x120*M[1] + x121*M[2] + x138*M[0];
#pragma omp atomic
L[39] += x121*M[1] + x157*M[2] + x93*M[0];
#pragma omp atomic
L[40] += x143*M[0] + x157*M[1] + x161*M[2];
#pragma omp atomic
L[41] += x120*M[0] + x122*M[2] + x136*M[1];
#pragma omp atomic
L[42] += x121*M[0] + x122*M[1] + x166*M[2];
#pragma omp atomic
L[43] += x157*M[0] + x163*M[2] + x166*M[1];
#pragma omp atomic
L[44] += x161*M[0] + x163*M[1] + x173*M[2];
#pragma omp atomic
L[45] += x136*M[0] + x98*M[1] + x99*M[2];
#pragma omp atomic
L[46] += x122*M[0] + x158*M[2] + x99*M[1];
#pragma omp atomic
L[47] += x158*M[1] + x160*M[2] + x166*M[0];
#pragma omp atomic
L[48] += x160*M[1] + x163*M[0] + x181*M[2];
#pragma omp atomic
L[49] += x173*M[0] + x181*M[1] + x182*M[2];
#pragma omp atomic
L[50] += x104*M[1] + x97*M[2] + x98*M[0];
#pragma omp atomic
L[51] += x146*M[2] + x97*M[1] + x99*M[0];
#pragma omp atomic
L[52] += x146*M[1] + x158*M[0] + x164*M[2];
#pragma omp atomic
L[53] += x160*M[0] + x164*M[1] + x171*M[2];
#pragma omp atomic
L[54] += x171*M[1] + x180*M[2] + x181*M[0];
#pragma omp atomic
L[55] += x179*M[2] + x180*M[1] + x182*M[0];

}

void L2L_6(double x, double y, double z, double * L, double * Ls) {
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
Ls[0] += (1.0/120.0)*pow(x, 5)*L[35] + x*x0 + x*x1 + x*x4 + x*L[1] + (1.0/12.0)*x10*x17*L[53] + x10*x54*L[38] + x11*x20 + x11*x35 + x11*x43 + x11*L[7] + (1.0/12.0)*x12*x15*L[52] + x12*x52*L[41] + x13*x21 + x13*x36 + x13*x45 + x13*L[16] + x14*x22 + x14*x37 + x14*L[30] + x15*x54*L[40] + x16*x23 + x16*x32 + x16*x39 + x16*L[9] + x17*x52*L[44] + x18*x24 + x18*x33 + x18*x41 + x18*L[19] + x19*x25 + x19*x34 + x19*L[34] + x2*y + x26*x6 + x27*x6 + x28*x8 + x29*x8 + x30*x9 + x31*x9 + x47*x6 + x49*x8 + x51*x58 + x51*L[23] + x53*x57 + x53*L[25] + x55*x56 + x55*L[32] + x6*L[4] + x8*L[10] + x9*L[20] + (1.0/120.0)*pow(y, 5)*L[50] + y*L[2] + (1.0/120.0)*pow(z, 5)*L[55] + z*L[3] + L[0];
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

void L2P_6(double x, double y, double z, double * L, double * F) {
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
F[0] += (1.0/120.0)*pow(x, 5)*L[35] + x*x11*L[26] + x*x12*L[45] + x*x14*L[15] + x*x16*L[29] + x*x17*L[49] + x*x23*L[47] + x*x9*L[13] + x*L[1] + x0*x14*L[28] + x0*x16*L[48] + x0*z*L[14] + x0*L[5] + x1*x11*L[46] + x1*x9*L[27] + x1*L[6] + (1.0/12.0)*x10*x13*L[52] + x10*x20*L[41] + x11*z*L[31] + x11*L[16] + x12*z*L[51] + x12*L[30] + x13*x22*L[40] + x14*y*L[18] + x14*L[9] + x15*x20*L[44] + (1.0/12.0)*x15*x8*L[53] + x16*y*L[33] + x16*L[19] + x17*y*L[54] + x17*L[34] + x19*z*L[42] + x19*L[23] + x2*x4*L[24] + x2*x6*L[39] + x2*L[8] + x21*y*L[43] + x21*L[25] + x22*x8*L[38] + x23*L[32] + x4*y*L[11] + x4*z*L[12] + x4*L[4] + x6*y*L[21] + x6*z*L[22] + x6*L[10] + x7*y*L[36] + x7*z*L[37] + x7*L[20] + x9*z*L[17] + x9*L[7] + (1.0/120.0)*pow(y, 5)*L[50] + y*L[2] + (1.0/120.0)*pow(z, 5)*L[55] + z*L[3] + L[0];

}

void M2P_6(double x, double y, double z, double * M, double * F) {
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
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*y;
x8 = x*z;
x9 = y*z;
x10 = pow(x2, -3.5);
x11 = 15.0*x10;
x12 = x7*z;
x13 = -x4;
x14 = x0*x6;
x15 = x1*x6;
x16 = 9.0*x5;
x17 = x0*x11;
x18 = -x17;
x19 = x*(x16 + x18);
x20 = x18 + x6;
x21 = x20*y;
x22 = x1*x11;
x23 = -x22;
x24 = y*(x16 + x23);
x25 = x20*z;
x26 = z*(x23 + x6);
x27 = 1.0*x;
x28 = x27*(x22 - x6);
x29 = 45.0*x10;
x30 = -x29;
x31 = pow(x2, -4.5);
x32 = x0*x31;
x33 = 105.0*x32;
x34 = x30 + x33;
x35 = x34*x7;
x36 = x34*x8;
x37 = -x11;
x38 = x9*(x33 + x37);
x39 = x1*x31;
x40 = 105.0*x39;
x41 = x30 + x40;
x42 = x41*x9;
x43 = x27*y;
x44 = x41*x43;
x45 = x37 + x40;
x46 = x27*z;
x47 = x45*x46;
x48 = 315.0*x31;
x49 = pow(x2, -5.5);
x50 = 945.0*x49;
x51 = x0*x50;
x52 = x48 - x51;
x53 = x12*x52;
x54 = x1*x50;
x55 = x27*x9*(-x48 + x54);
x56 = 90.0*x10;
x57 = (x*x*x*x);
x58 = 105.0*x31;
x59 = x57*x58;
x60 = (y*y*y*y);
x61 = x58*x60;
x62 = -225.0*x10;
x63 = x50*x57;
x64 = -x63;
x65 = x*(1050.0*x32 + x62 + x64);
x66 = x30 + 630.0*x32 + x64;
x67 = x66*y;
x68 = x50*x60;
x69 = -x68;
x70 = y*(1050.0*x39 + x62 + x69);
x71 = x66*z;
x72 = 630.0*x39;
x73 = x30 + x69 + x72;
x74 = x73*z;
x75 = x27*(x29 + x68 - x72);
x76 = 1575.0*x31;
x77 = x0*x49;
x78 = pow(x2, -6.5);
x79 = 10395.0*x78;
x80 = x57*x79;
x81 = x76 - 9450.0*x77 + x80;
x82 = x7*x81;
x83 = x8*x81;
x84 = 5670.0*x77;
x85 = x9*(x48 + x80 - x84);
x86 = x1*x49;
x87 = x60*x79;
x88 = x76 - 9450.0*x86 + x87;
x89 = x88*x9;
x90 = x43*x88;
x91 = x46*(x48 - 5670.0*x86 + x87);
x92 = 14175.0*x49;
x93 = pow(x, 6)*x79;
x94 = x79*pow(y, 6);
x95 = x1*x33;
x96 = x1*x48;
x97 = -x1*x51;
x98 = x*(x34 + x96 + x97);
x99 = x0*x48;
x100 = y*(x41 + x97 + x99);
x101 = z*(x33 + x45 + x97);
x102 = 945.0*x31;
x103 = -2835.0*x77;
x104 = x0*x1*x79;
x105 = x104 - 2835.0*x86;
x106 = x7*(x102 + x103 + x105);
x107 = x8*(x105 + x52);
x108 = x9*(x103 + x104 + x48 - x54);
x109 = 105.0*x10;
x110 = -12.0*x5 - x95;
x111 = 120.0*x10;
x112 = x0*x87;
x113 = x1*x84;
x114 = -x113;
x115 = x1*x80;
x116 = 15120.0*x49;
x117 = -x93;
x118 = -x115;
x119 = 270.0*x10 + x113;
x120 = -x94;
x121 = -x112;
x122 = -x98;
x123 = -x100;
x124 = -x101;
x125 = -x106;
x126 = -x107;
x127 = -x108;
x128 = x1*x77;
x129 = 16065.0*x49;
x130 = 20790.0*x78;
x131 = x1*x57;
x132 = -360.0*x10 - 17010.0*x128;
x133 = x0*x60;
x134 = 17010.0*x49;
x135 = 31185.0*x78;
#pragma omp atomic
F[0] += -x*x4*M[0] + x100*M[40] + x101*M[41] + x106*M[61] + x107*M[62] + x108*M[66] - x11*x12*M[13] + x19*M[9] + x21*M[10] + x24*M[15] + x25*M[11] + x26*M[16] - x28*M[12] + x35*M[20] + x36*M[21] + x38*M[23] - x4*y*M[1] - x4*z*M[2] + x42*M[30] + x44*M[25] + x47*M[26] + x53*M[38] - x55*M[45] + x6*x7*M[4] + x6*x8*M[5] + x6*x9*M[7] + x65*M[34] + x67*M[35] + x70*M[49] + x71*M[36] + x74*M[50] - x75*M[44] + x82*M[56] + x83*M[57] + x85*M[59] + x89*M[77] + x90*M[70] + x91*M[71] + x98*M[37] + (x122 - x65)*M[39] + (x122 + x75)*M[46] + (x123 - x67)*M[42] + (x123 - x70)*M[51] + (x124 - x71)*M[43] + (x124 - x74)*M[52] + (x125 - x82)*M[63] + (x125 - x90)*M[72] + (x126 - x83)*M[64] + (x126 - x91)*M[73] + (x127 - x85)*M[68] + (x127 - x89)*M[79] + (x13 + x14)*M[3] + (x13 + x15)*M[6] + (-x19 + x28)*M[14] + (-x21 - x24)*M[17] + (-x25 - x26)*M[18] + (-x35 - x44)*M[27] + (-x36 - x47)*M[28] + (-x38 - x42)*M[32] + (-x53 + x55)*M[47] + (2*x100 + x67 + x70)*M[53] + (2*x101 + x71 + x74)*M[54] + (2*x106 + x82 + x90)*M[74] + (2*x107 + x83 + x91)*M[75] + (2*x108 + x85 + x89)*M[81] + (-x14 - x15 + 2.0*x3)*M[8] + (x20 + x23 + x95)*M[22] + (x65 - x75 + 2*x98)*M[48] + (-x0*x56 + x16 + x59)*M[19] + (-x1*x56 + x16 + x61)*M[29] + (x112 + x114 + x73 + x99)*M[65] + (x114 + x115 + x66 + x96)*M[58] + (4725.0*x32 - x57*x92 + x62 + x93)*M[55] + (4725.0*x39 - x60*x92 + x62 + x94)*M[76] + (x0*x109 + x110 + x22 - x59)*M[24] + (x1*x109 + x110 + x17 - x61)*M[31] + (-x0*x111 - x1*x111 + 210.0*x1*x32 + 24.0*x5 + x59 + x61)*M[33] + (x116*x57 + x117 + x118 + x119 - 5355.0*x32 - x96)*M[60] + (x116*x60 + x119 + x120 + x121 - 5355.0*x39 - x99)*M[78] + (x112 - x129*x57 + x130*x131 + x132 + 6300.0*x32 + 1260.0*x39 + x69 + x93)*M[69] + (x115 - x129*x60 + x130*x133 + x132 + 1260.0*x32 + 6300.0*x39 + x64 + x94)*M[80] + (-x0*x102 - x1*x102 + x118 + x121 + 11340.0*x128 + x56 + x63 + x68)*M[67] + (720.0*x10 + x117 + x120 + 34020.0*x128 - x131*x135 - x133*x135 + x134*x57 + x134*x60 - 7560.0*x32 - 7560.0*x39)*M[82];

}

void P2M_7(double x, double y, double z, double q, double * M) {
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
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = x1*z;
x9 = (z*z);
x10 = (x*x*x);
x11 = (1.0/6.0)*q;
x12 = (1.0/2.0)*x3;
x13 = (1.0/2.0)*x0;
x14 = (y*y*y);
x15 = (1.0/2.0)*x7;
x16 = (1.0/2.0)*x9;
x17 = (z*z*z);
x18 = (x*x*x*x);
x19 = (1.0/24.0)*q;
x20 = (1.0/6.0)*x10;
x21 = q*x7;
x22 = (1.0/4.0)*x3;
x23 = q*x9;
x24 = (1.0/6.0)*x0;
x25 = (y*y*y*y);
x26 = (1.0/6.0)*x14;
x27 = (1.0/4.0)*x9;
x28 = (1.0/6.0)*x17;
x29 = (z*z*z*z);
x30 = pow(x, 5);
x31 = (1.0/120.0)*q;
x32 = (1.0/24.0)*x18;
x33 = (1.0/12.0)*x10;
x34 = (1.0/12.0)*x14;
x35 = q*x3;
x36 = x2*x7;
x37 = x1*x9;
x38 = (1.0/12.0)*x17;
x39 = (1.0/24.0)*x0;
x40 = x0*x7;
x41 = pow(y, 5);
x42 = (1.0/24.0)*x25;
x43 = (1.0/24.0)*x29;
x44 = pow(z, 5);
x45 = pow(x, 6);
x46 = (1.0/720.0)*q;
x47 = (1.0/120.0)*x30;
x48 = (1.0/48.0)*x18;
x49 = q*x14;
x50 = (1.0/36.0)*x10;
x51 = q*x17;
x52 = (1.0/48.0)*x35;
x53 = x2*x3;
x54 = x3*x9;
x55 = x1*x3;
x56 = (1.0/120.0)*x0;
x57 = x0*x9;
x58 = pow(y, 6);
x59 = (1.0/120.0)*x41;
x60 = (1.0/48.0)*x25;
x61 = (1.0/36.0)*x17;
x62 = (1.0/48.0)*x29;
x63 = (1.0/120.0)*x44;
x64 = pow(z, 6);
x65 = (1.0/5040.0)*q;
x66 = (1.0/720.0)*x45;
x67 = (1.0/240.0)*x30;
x68 = (1.0/144.0)*x18;
x69 = (1.0/144.0)*x25;
x70 = q*x10;
x71 = x19*x7;
x72 = (1.0/144.0)*x29;
x73 = (1.0/240.0)*x35;
x74 = (1.0/720.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x6;
M[6] += x4*x7;
M[7] += x8;
M[8] += x4*x9;
M[9] += -x10*x11;
M[10] += -x1*x12;
M[11] += -x12*x2;
M[12] += -x13*x7;
M[13] += -x5*z;
M[14] += -x13*x9;
M[15] += -x11*x14;
M[16] += -x15*x2;
M[17] += -x1*x16;
M[18] += -x11*x17;
M[19] += x18*x19;
M[20] += x1*x20;
M[21] += x2*x20;
M[22] += x21*x22;
M[23] += x12*x8;
M[24] += x22*x23;
M[25] += x14*x24;
M[26] += x15*x6;
M[27] += x16*x5;
M[28] += x17*x24;
M[29] += x19*x25;
M[30] += x2*x26;
M[31] += x21*x27;
M[32] += x1*x28;
M[33] += x19*x29;
M[34] += -x30*x31;
M[35] += -x1*x32;
M[36] += -x2*x32;
M[37] += -x21*x33;
M[38] += -x20*x8;
M[39] += -x23*x33;
M[40] += -x34*x35;
M[41] += -x22*x36;
M[42] += -x22*x37;
M[43] += -x35*x38;
M[44] += -x25*x39;
M[45] += -x26*x6;
M[46] += -x27*x40;
M[47] += -x28*x5;
M[48] += -x29*x39;
M[49] += -x31*x41;
M[50] += -x2*x42;
M[51] += -x23*x34;
M[52] += -x21*x38;
M[53] += -x1*x43;
M[54] += -x31*x44;
M[55] += x45*x46;
M[56] += x1*x47;
M[57] += x2*x47;
M[58] += x21*x48;
M[59] += x32*x8;
M[60] += x23*x48;
M[61] += x49*x50;
M[62] += x33*x36;
M[63] += x33*x37;
M[64] += x50*x51;
M[65] += x25*x52;
M[66] += x34*x53;
M[67] += (1.0/8.0)*x21*x54;
M[68] += x38*x55;
M[69] += x29*x52;
M[70] += x41*x56;
M[71] += x42*x6;
M[72] += x34*x57;
M[73] += x38*x40;
M[74] += x43*x5;
M[75] += x44*x56;
M[76] += x46*x58;
M[77] += x2*x59;
M[78] += x23*x60;
M[79] += x49*x61;
M[80] += x21*x62;
M[81] += x1*x63;
M[82] += x46*x64;
M[83] += -pow(x, 7)*x65;
M[84] += -x1*x66;
M[85] += -x2*x66;
M[86] += -x21*x67;
M[87] += -x47*x8;
M[88] += -x23*x67;
M[89] += -x49*x68;
M[90] += -x36*x48;
M[91] += -x37*x48;
M[92] += -x51*x68;
M[93] += -x69*x70;
M[94] += -x14*x2*x50;
M[95] += -x10*x71*x9;
M[96] += -x1*x17*x50;
M[97] += -x70*x72;
M[98] += -x41*x73;
M[99] += -x53*x60;
M[100] += -x14*x19*x54;
M[101] += -x17*x3*x71;
M[102] += -x55*x62;
M[103] += -x44*x73;
M[104] += -x58*x74;
M[105] += -x59*x6;
M[106] += -x57*x60;
M[107] += -x0*x14*x61;
M[108] += -x40*x62;
M[109] += -x5*x63;
M[110] += -x64*x74;
M[111] += -x65*pow(y, 7);
M[112] += -1.0/720.0*x2*x58;
M[113] += -1.0/240.0*x23*x41;
M[114] += -x51*x69;
M[115] += -x49*x72;
M[116] += -1.0/240.0*x21*x44;
M[117] += -1.0/720.0*x1*x64;
M[118] += -x65*pow(z, 7);
}
void M2M_7(double x, double y, double z, double * M, double * Ms) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = z*M[0];
x5 = y*M[1];
x6 = y*M[2];
x7 = z*M[1];
x8 = z*M[2];
x9 = x*M[3];
x10 = (x*x);
x11 = (1.0/2.0)*x10;
x12 = x*M[4];
x13 = y*M[3];
x14 = x0*y;
x15 = x*M[5];
x16 = z*M[3];
x17 = x0*z;
x18 = x*M[6];
x19 = y*M[4];
x20 = x1*y;
x21 = (y*y);
x22 = (1.0/2.0)*M[0];
x23 = x*M[7];
x24 = y*M[5];
x25 = z*M[4];
x26 = x3*y;
x27 = x1*z;
x28 = x2*z;
x29 = x*M[8];
x30 = z*M[5];
x31 = x3*z;
x32 = (z*z);
x33 = y*M[6];
x34 = (1.0/2.0)*x21;
x35 = y*M[7];
x36 = z*M[6];
x37 = x5*z;
x38 = y*M[8];
x39 = z*M[7];
x40 = x6*z;
x41 = (1.0/2.0)*x32;
x42 = z*M[8];
x43 = x*M[9];
x44 = (x*x*x);
x45 = (1.0/6.0)*x44;
x46 = x*M[10];
x47 = y*M[9];
x48 = x9*y;
x49 = x*M[11];
x50 = z*M[9];
x51 = x9*z;
x52 = x*M[12];
x53 = y*M[10];
x54 = x12*y;
x55 = x*M[13];
x56 = y*M[11];
x57 = z*M[10];
x58 = x15*y;
x59 = x12*z;
x60 = x13*z;
x61 = x*M[14];
x62 = z*M[11];
x63 = x15*z;
x64 = x*M[15];
x65 = y*M[12];
x66 = x18*y;
x67 = (y*y*y);
x68 = (1.0/6.0)*M[0];
x69 = x*M[16];
x70 = y*M[13];
x71 = z*M[12];
x72 = x23*y;
x73 = x18*z;
x74 = x19*z;
x75 = x*M[17];
x76 = y*M[14];
x77 = z*M[13];
x78 = x29*y;
x79 = x23*z;
x80 = x24*z;
x81 = x*M[18];
x82 = z*M[14];
x83 = x29*z;
x84 = (z*z*z);
x85 = y*M[15];
x86 = (1.0/6.0)*x67;
x87 = y*M[16];
x88 = z*M[15];
x89 = x33*z;
x90 = y*M[17];
x91 = z*M[16];
x92 = x35*z;
x93 = y*M[18];
x94 = z*M[17];
x95 = x38*z;
x96 = (1.0/6.0)*x84;
x97 = z*M[18];
x98 = x*M[19];
x99 = (x*x*x*x);
x100 = (1.0/24.0)*x99;
x101 = x*M[20];
x102 = y*M[19];
x103 = x43*y;
x104 = x*M[21];
x105 = z*M[19];
x106 = x43*z;
x107 = x*M[22];
x108 = y*M[20];
x109 = x46*y;
x110 = (1.0/4.0)*x10;
x111 = x21*M[0];
x112 = x*M[23];
x113 = y*M[21];
x114 = z*M[20];
x115 = x49*y;
x116 = x46*z;
x117 = x47*z;
x118 = x*M[24];
x119 = z*M[21];
x120 = x49*z;
x121 = x110*x32;
x122 = x*M[25];
x123 = y*M[22];
x124 = x52*y;
x125 = x110*x21;
x126 = x*M[26];
x127 = y*M[23];
x128 = z*M[22];
x129 = x55*y;
x130 = x52*z;
x131 = x53*z;
x132 = x*M[27];
x133 = y*M[24];
x134 = z*M[23];
x135 = x61*y;
x136 = x55*z;
x137 = x56*z;
x138 = x*M[28];
x139 = z*M[24];
x140 = x61*z;
x141 = x*M[29];
x142 = y*M[25];
x143 = x64*y;
x144 = (y*y*y*y);
x145 = (1.0/24.0)*M[0];
x146 = x*M[30];
x147 = y*M[26];
x148 = z*M[25];
x149 = x69*y;
x150 = x64*z;
x151 = x65*z;
x152 = x*M[31];
x153 = y*M[27];
x154 = z*M[26];
x155 = x75*y;
x156 = x69*z;
x157 = x70*z;
x158 = (1.0/4.0)*x32;
x159 = x*M[32];
x160 = y*M[28];
x161 = z*M[27];
x162 = x81*y;
x163 = x75*z;
x164 = x76*z;
x165 = x*M[33];
x166 = z*M[28];
x167 = x81*z;
x168 = (z*z*z*z);
x169 = y*M[29];
x170 = (1.0/24.0)*x144;
x171 = y*M[30];
x172 = z*M[29];
x173 = x85*z;
x174 = y*M[31];
x175 = z*M[30];
x176 = x87*z;
x177 = x158*x21;
x178 = y*M[32];
x179 = z*M[31];
x180 = x90*z;
x181 = y*M[33];
x182 = z*M[32];
x183 = x93*z;
x184 = (1.0/24.0)*x168;
x185 = z*M[33];
x186 = x*M[34];
x187 = (1.0/120.0)*pow(x, 5);
x188 = x*M[35];
x189 = y*M[34];
x190 = x98*y;
x191 = x*M[36];
x192 = x*M[37];
x193 = y*M[35];
x194 = x101*y;
x195 = (1.0/12.0)*x44;
x196 = x*M[38];
x197 = y*M[36];
x198 = x104*y;
x199 = x*M[39];
x200 = x195*x32;
x201 = x*M[40];
x202 = y*M[37];
x203 = x107*y;
x204 = (1.0/12.0)*x10;
x205 = x67*M[0];
x206 = x195*x21;
x207 = x*M[41];
x208 = y*M[38];
x209 = x112*y;
x210 = x*M[42];
x211 = y*M[39];
x212 = x118*y;
x213 = x*M[43];
x214 = x204*x84;
x215 = x*M[44];
x216 = y*M[40];
x217 = x122*y;
x218 = x204*x67;
x219 = x*M[45];
x220 = y*M[41];
x221 = x126*y;
x222 = x*M[46];
x223 = y*M[42];
x224 = x132*y;
x225 = x*M[47];
x226 = y*M[43];
x227 = x138*y;
x228 = x*M[48];
x229 = x*M[49];
x230 = y*M[44];
x231 = x141*y;
x232 = pow(y, 5);
x233 = (1.0/120.0)*M[0];
x234 = x*M[50];
x235 = y*M[45];
x236 = x146*y;
x237 = x*M[51];
x238 = y*M[46];
x239 = x152*y;
x240 = (1.0/12.0)*x32;
x241 = x*M[52];
x242 = y*M[47];
x243 = x159*y;
x244 = (1.0/12.0)*x84;
x245 = x*M[53];
x246 = y*M[48];
x247 = x165*y;
x248 = x*M[54];
x249 = pow(z, 5);
x250 = y*M[49];
x251 = (1.0/120.0)*x232;
x252 = y*M[50];
x253 = y*M[51];
x254 = x240*x67;
x255 = y*M[52];
x256 = x21*x244;
x257 = y*M[53];
x258 = y*M[54];
x259 = (1.0/120.0)*x249;
x260 = (1.0/720.0)*pow(x, 6);
x261 = (1.0/48.0)*x99;
x262 = x261*x32;
x263 = (1.0/36.0)*x44;
x264 = x21*x261;
x265 = x263*x84;
x266 = (1.0/48.0)*x10;
x267 = x266*M[0];
x268 = x263*x67;
x269 = (1.0/8.0)*x10*x32;
x270 = x144*x266;
x271 = x21*x269;
x272 = x168*x266;
x273 = pow(y, 6);
x274 = (1.0/720.0)*M[0];
x275 = (1.0/48.0)*x144*x32;
x276 = (1.0/36.0)*x84;
x277 = (1.0/48.0)*x168;
x278 = pow(z, 6);
x279 = (1.0/720.0)*x273;
x280 = x276*x67;
x281 = x21*x277;
x282 = (1.0/720.0)*x278;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + x4 + M[5];
#pragma omp atomic
Ms[6] += x5 + M[6];
#pragma omp atomic
Ms[7] += x6 + x7 + M[7];
#pragma omp atomic
Ms[8] += x8 + M[8];
#pragma omp atomic
Ms[9] += x11*M[0] + x9 + M[9];
#pragma omp atomic
Ms[10] += x11*M[1] + x12 + x13 + x14 + M[10];
#pragma omp atomic
Ms[11] += x11*M[2] + x15 + x16 + x17 + M[11];
#pragma omp atomic
Ms[12] += x18 + x19 + x20 + x21*x22 + M[12];
#pragma omp atomic
Ms[13] += x23 + x24 + x25 + x26 + x27 + x28 + M[13];
#pragma omp atomic
Ms[14] += x22*x32 + x29 + x30 + x31 + M[14];
#pragma omp atomic
Ms[15] += x33 + x34*M[1] + M[15];
#pragma omp atomic
Ms[16] += x34*M[2] + x35 + x36 + x37 + M[16];
#pragma omp atomic
Ms[17] += x38 + x39 + x40 + x41*M[1] + M[17];
#pragma omp atomic
Ms[18] += x41*M[2] + x42 + M[18];
#pragma omp atomic
Ms[19] += x11*M[3] + x43 + x45*M[0] + M[19];
#pragma omp atomic
Ms[20] += x11*x2 + x11*M[4] + x45*M[1] + x46 + x47 + x48 + M[20];
#pragma omp atomic
Ms[21] += x11*x4 + x11*M[5] + x45*M[2] + x49 + x50 + x51 + M[21];
#pragma omp atomic
Ms[22] += x0*x34 + x11*x5 + x11*M[6] + x34*M[3] + x52 + x53 + x54 + M[22];
#pragma omp atomic
Ms[23] += x11*x6 + x11*x7 + x11*M[7] + x14*z + x55 + x56 + x57 + x58 + x59 + x60 + M[23];
#pragma omp atomic
Ms[24] += x0*x41 + x11*x8 + x11*M[8] + x41*M[3] + x61 + x62 + x63 + M[24];
#pragma omp atomic
Ms[25] += x1*x34 + x34*M[4] + x64 + x65 + x66 + x67*x68 + M[25];
#pragma omp atomic
Ms[26] += x20*z + x3*x34 + x34*x4 + x34*M[5] + x69 + x70 + x71 + x72 + x73 + x74 + M[26];
#pragma omp atomic
Ms[27] += x1*x41 + x2*x41 + x26*z + x41*M[4] + x75 + x76 + x77 + x78 + x79 + x80 + M[27];
#pragma omp atomic
Ms[28] += x3*x41 + x41*M[5] + x68*x84 + x81 + x82 + x83 + M[28];
#pragma omp atomic
Ms[29] += x34*M[6] + x85 + x86*M[1] + M[29];
#pragma omp atomic
Ms[30] += x34*x7 + x34*M[7] + x86*M[2] + x87 + x88 + x89 + M[30];
#pragma omp atomic
Ms[31] += x34*x8 + x34*M[8] + x41*x5 + x41*M[6] + x90 + x91 + x92 + M[31];
#pragma omp atomic
Ms[32] += x41*x6 + x41*M[7] + x93 + x94 + x95 + x96*M[1] + M[32];
#pragma omp atomic
Ms[33] += x41*M[8] + x96*M[2] + x97 + M[33];
#pragma omp atomic
Ms[34] += x100*M[0] + x11*M[9] + x45*M[3] + x98 + M[34];
#pragma omp atomic
Ms[35] += x100*M[1] + x101 + x102 + x103 + x11*x13 + x11*M[10] + x2*x45 + x45*M[4] + M[35];
#pragma omp atomic
Ms[36] += x100*M[2] + x104 + x105 + x106 + x11*x16 + x11*M[11] + x4*x45 + x45*M[5] + M[36];
#pragma omp atomic
Ms[37] += x107 + x108 + x109 + x11*x19 + x11*M[12] + x110*x111 + x34*x9 + x34*M[9] + x45*x5 + x45*M[6] + M[37];
#pragma omp atomic
Ms[38] += x11*x24 + x11*x25 + x11*x28 + x11*M[13] + x112 + x113 + x114 + x115 + x116 + x117 + x45*x6 + x45*x7 + x45*M[7] + x48*z + M[38];
#pragma omp atomic
Ms[39] += x11*x30 + x11*M[14] + x118 + x119 + x120 + x121*M[0] + x41*x9 + x41*M[9] + x45*x8 + x45*M[8] + M[39];
#pragma omp atomic
Ms[40] += x0*x86 + x11*x33 + x11*M[15] + x12*x34 + x122 + x123 + x124 + x125*M[1] + x34*M[10] + x86*M[3] + M[40];
#pragma omp atomic
Ms[41] += x11*x35 + x11*x36 + x11*x37 + x11*M[16] + x125*M[2] + x126 + x127 + x128 + x129 + x130 + x131 + x15*x34 + x16*x34 + x17*x34 + x34*M[11] + x54*z + M[41];
#pragma omp atomic
Ms[42] += x11*x38 + x11*x39 + x11*x40 + x11*M[17] + x12*x41 + x121*M[1] + x13*x41 + x132 + x133 + x134 + x135 + x136 + x137 + x14*x41 + x41*M[10] + x58*z + M[42];
#pragma omp atomic
Ms[43] += x0*x96 + x11*x42 + x11*M[18] + x121*M[2] + x138 + x139 + x140 + x15*x41 + x41*M[11] + x96*M[3] + M[43];
#pragma omp atomic
Ms[44] += x1*x86 + x141 + x142 + x143 + x144*x145 + x18*x34 + x34*M[12] + x86*M[4] + M[44];
#pragma omp atomic
Ms[45] += x146 + x147 + x148 + x149 + x150 + x151 + x23*x34 + x25*x34 + x27*x34 + x3*x86 + x34*M[13] + x4*x86 + x66*z + x86*M[5] + M[45];
#pragma omp atomic
Ms[46] += x111*x158 + x152 + x153 + x154 + x155 + x156 + x157 + x18*x41 + x19*x41 + x20*x41 + x29*x34 + x30*x34 + x31*x34 + x34*M[14] + x41*M[12] + x72*z + M[46];
#pragma omp atomic
Ms[47] += x1*x96 + x159 + x160 + x161 + x162 + x163 + x164 + x2*x96 + x23*x41 + x24*x41 + x26*x41 + x41*M[13] + x78*z + x96*M[4] + M[47];
#pragma omp atomic
Ms[48] += x145*x168 + x165 + x166 + x167 + x29*x41 + x3*x96 + x41*M[14] + x96*M[5] + M[48];
#pragma omp atomic
Ms[49] += x169 + x170*M[1] + x34*M[15] + x86*M[6] + M[49];
#pragma omp atomic
Ms[50] += x170*M[2] + x171 + x172 + x173 + x34*x36 + x34*M[16] + x7*x86 + x86*M[7] + M[50];
#pragma omp atomic
Ms[51] += x174 + x175 + x176 + x177*M[1] + x33*x41 + x34*x39 + x34*M[17] + x41*M[15] + x8*x86 + x86*M[8] + M[51];
#pragma omp atomic
Ms[52] += x177*M[2] + x178 + x179 + x180 + x34*x42 + x34*M[18] + x35*x41 + x41*M[16] + x5*x96 + x96*M[6] + M[52];
#pragma omp atomic
Ms[53] += x181 + x182 + x183 + x184*M[1] + x38*x41 + x41*M[17] + x6*x96 + x96*M[7] + M[53];
#pragma omp atomic
Ms[54] += x184*M[2] + x185 + x41*M[18] + x96*M[8] + M[54];
#pragma omp atomic
Ms[55] += x100*M[3] + x11*M[19] + x186 + x187*M[0] + x45*M[9] + M[55];
#pragma omp atomic
Ms[56] += x100*x2 + x100*M[4] + x11*x47 + x11*M[20] + x13*x45 + x187*M[1] + x188 + x189 + x190 + x45*M[10] + M[56];
#pragma omp atomic
Ms[57] += x100*x4 + x100*M[5] + x11*x50 + x11*M[21] + x16*x45 + x187*M[2] + x191 + x45*M[11] + x98*z + z*M[34] + M[57];
#pragma omp atomic
Ms[58] += x100*x5 + x100*M[6] + x11*x53 + x11*M[22] + x111*x195 + x125*M[3] + x19*x45 + x192 + x193 + x194 + x34*x43 + x34*M[19] + x45*M[12] + M[58];
#pragma omp atomic
Ms[59] += x100*x6 + x100*x7 + x100*M[7] + x101*z + x102*z + x103*z + x11*x56 + x11*x57 + x11*x60 + x11*M[23] + x196 + x197 + x198 + x24*x45 + x25*x45 + x28*x45 + x45*M[13] + z*M[35] + M[59];
#pragma omp atomic
Ms[60] += x100*x8 + x100*M[8] + x104*z + x11*x62 + x11*M[24] + x121*M[3] + x199 + x200*M[0] + x30*x45 + x41*x43 + x41*M[19] + x45*M[14] + z*M[36] + M[60];
#pragma omp atomic
Ms[61] += x11*x65 + x11*M[25] + x125*M[4] + x201 + x202 + x203 + x204*x205 + x206*M[1] + x33*x45 + x34*x46 + x34*M[20] + x45*M[15] + x86*x9 + x86*M[9] + M[61];
#pragma omp atomic
Ms[62] += x107*z + x108*z + x109*z + x11*x70 + x11*x71 + x11*x74 + x11*M[26] + x125*x4 + x125*M[5] + x206*M[2] + x207 + x208 + x209 + x34*x49 + x34*x50 + x34*x51 + x34*M[21] + x35*x45 + x36*x45 + x37*x45 + x45*M[16] + z*M[37] + M[62];
#pragma omp atomic
Ms[63] += x11*x76 + x11*x77 + x11*x80 + x11*M[27] + x112*z + x113*z + x115*z + x121*x2 + x121*M[4] + x200*M[1] + x210 + x211 + x212 + x38*x45 + x39*x45 + x40*x45 + x41*x46 + x41*x47 + x41*x48 + x41*M[20] + x45*M[17] + z*M[38] + M[63];
#pragma omp atomic
Ms[64] += x11*x82 + x11*M[28] + x118*z + x121*M[5] + x200*M[2] + x213 + x214*M[0] + x41*x49 + x41*M[21] + x42*x45 + x45*M[18] + x9*x96 + x96*M[9] + z*M[39] + M[64];
#pragma omp atomic
Ms[65] += x0*x170 + x11*x85 + x11*M[29] + x12*x86 + x125*M[6] + x170*M[3] + x215 + x216 + x217 + x218*M[1] + x34*x52 + x34*M[22] + x86*M[10] + M[65];
#pragma omp atomic
Ms[66] += x11*x87 + x11*x88 + x11*x89 + x11*M[30] + x122*z + x123*z + x124*z + x125*x7 + x125*M[7] + x15*x86 + x16*x86 + x17*x86 + x218*M[2] + x219 + x220 + x221 + x34*x55 + x34*x57 + x34*x59 + x34*M[23] + x86*M[11] + z*M[40] + M[66];
#pragma omp atomic
Ms[67] += x0*x177 + x11*x90 + x11*x91 + x11*x92 + x11*M[31] + x121*x5 + x121*M[6] + x125*x8 + x125*M[8] + x126*z + x127*z + x129*z + x177*M[3] + x222 + x223 + x224 + x34*x61 + x34*x62 + x34*x63 + x34*M[24] + x41*x52 + x41*x53 + x41*x54 + x41*M[22] + z*M[41] + M[67];
#pragma omp atomic
Ms[68] += x11*x93 + x11*x94 + x11*x95 + x11*M[32] + x12*x96 + x121*x6 + x121*M[7] + x13*x96 + x132*z + x133*z + x135*z + x14*x96 + x214*M[1] + x225 + x226 + x227 + x41*x55 + x41*x56 + x41*x58 + x41*M[23] + x96*M[10] + z*M[42] + M[68];
#pragma omp atomic
Ms[69] += x0*x184 + x11*x97 + x11*M[33] + x121*M[8] + x138*z + x15*x96 + x184*M[3] + x214*M[2] + x228 + x41*x61 + x41*M[24] + x96*M[11] + z*M[43] + M[69];
#pragma omp atomic
Ms[70] += x1*x170 + x170*M[4] + x18*x86 + x229 + x230 + x231 + x232*x233 + x34*x64 + x34*M[25] + x86*M[12] + M[70];
#pragma omp atomic
Ms[71] += x141*z + x142*z + x143*z + x170*x3 + x170*x4 + x170*M[5] + x23*x86 + x234 + x235 + x236 + x25*x86 + x27*x86 + x34*x69 + x34*x71 + x34*x73 + x34*M[26] + x86*M[13] + z*M[44] + M[71];
#pragma omp atomic
Ms[72] += x1*x177 + x146*z + x147*z + x149*z + x177*M[4] + x205*x240 + x237 + x238 + x239 + x29*x86 + x30*x86 + x31*x86 + x34*x75 + x34*x77 + x34*x79 + x34*M[27] + x41*x64 + x41*x65 + x41*x66 + x41*M[25] + x86*M[14] + z*M[45] + M[72];
#pragma omp atomic
Ms[73] += x111*x244 + x152*z + x153*z + x155*z + x177*x3 + x177*M[5] + x18*x96 + x19*x96 + x20*x96 + x241 + x242 + x243 + x34*x81 + x34*x82 + x34*x83 + x34*M[28] + x41*x69 + x41*x70 + x41*x72 + x41*M[26] + x96*M[12] + z*M[46] + M[73];
#pragma omp atomic
Ms[74] += x1*x184 + x159*z + x160*z + x162*z + x184*x2 + x184*M[4] + x23*x96 + x24*x96 + x245 + x246 + x247 + x26*x96 + x41*x75 + x41*x76 + x41*x78 + x41*M[27] + x96*M[13] + z*M[47] + M[74];
#pragma omp atomic
Ms[75] += x165*z + x184*x3 + x184*M[5] + x233*x249 + x248 + x29*x96 + x41*x81 + x41*M[28] + x96*M[14] + z*M[48] + M[75];
#pragma omp atomic
Ms[76] += x170*M[6] + x250 + x251*M[1] + x34*M[29] + x86*M[15] + M[76];
#pragma omp atomic
Ms[77] += x169*z + x170*x7 + x170*M[7] + x251*M[2] + x252 + x34*x88 + x34*M[30] + x36*x86 + x86*M[16] + z*M[49] + M[77];
#pragma omp atomic
Ms[78] += x170*x8 + x170*M[8] + x171*z + x177*M[6] + x253 + x254*M[1] + x34*x91 + x34*M[31] + x39*x86 + x41*x85 + x41*M[29] + x86*M[17] + z*M[50] + M[78];
#pragma omp atomic
Ms[79] += x174*z + x177*M[7] + x254*M[2] + x255 + x256*M[1] + x33*x96 + x34*x94 + x34*M[32] + x41*x87 + x41*M[30] + x42*x86 + x86*M[18] + x96*M[15] + z*M[51] + M[79];
#pragma omp atomic
Ms[80] += x177*M[8] + x178*z + x184*x5 + x184*M[6] + x256*M[2] + x257 + x34*x97 + x34*M[33] + x35*x96 + x41*x90 + x41*M[31] + x96*M[16] + z*M[52] + M[80];
#pragma omp atomic
Ms[81] += x181*z + x184*x6 + x184*M[7] + x258 + x259*M[1] + x38*x96 + x41*x93 + x41*M[32] + x96*M[17] + z*M[53] + M[81];
#pragma omp atomic
Ms[82] += x184*M[8] + x259*M[2] + x41*M[33] + x96*M[18] + z*M[54] + M[82];
#pragma omp atomic
Ms[83] += x*M[55] + x100*M[9] + x11*M[34] + x187*M[3] + x260*M[0] + x45*M[19] + M[83];
#pragma omp atomic
Ms[84] += x*M[56] + x100*x13 + x100*M[10] + x102*x11 + x11*M[35] + x186*y + x187*x2 + x187*M[4] + x260*M[1] + x45*x47 + x45*M[20] + y*M[55] + M[84];
#pragma omp atomic
Ms[85] += x*M[57] + x100*x16 + x100*M[11] + x105*x11 + x11*M[36] + x186*z + x187*x4 + x187*M[5] + x260*M[2] + x45*x50 + x45*M[21] + z*M[55] + M[85];
#pragma omp atomic
Ms[86] += x*M[58] + x100*x19 + x100*M[12] + x108*x11 + x11*M[37] + x111*x261 + x125*M[9] + x187*x5 + x187*M[6] + x188*y + x206*M[3] + x34*x98 + x34*M[34] + x45*x53 + x45*M[22] + y*M[56] + M[86];
#pragma omp atomic
Ms[87] += x*M[59] + x100*x24 + x100*x25 + x100*x28 + x100*M[13] + x11*x113 + x11*x114 + x11*x117 + x11*M[38] + x187*x6 + x187*x7 + x187*M[7] + x188*z + x189*z + x190*z + x191*y + x45*x56 + x45*x57 + x45*x60 + x45*M[23] + y*M[57] + z*M[56] + M[87];
#pragma omp atomic
Ms[88] += x*M[60] + x100*x30 + x100*M[14] + x11*x119 + x11*M[39] + x121*M[9] + x187*x8 + x187*M[8] + x191*z + x200*M[3] + x262*M[0] + x41*x98 + x41*M[34] + x45*x62 + x45*M[24] + z*M[57] + M[88];
#pragma omp atomic
Ms[89] += x*M[61] + x100*x33 + x100*M[15] + x101*x34 + x11*x123 + x11*M[40] + x125*M[10] + x192*y + x205*x263 + x206*M[4] + x218*M[3] + x264*M[1] + x34*M[35] + x43*x86 + x45*x65 + x45*M[25] + x86*M[19] + y*M[58] + M[89];
#pragma omp atomic
Ms[90] += x*M[62] + x100*x35 + x100*x36 + x100*x37 + x100*M[16] + x104*x34 + x105*x34 + x106*x34 + x11*x127 + x11*x128 + x11*x131 + x11*M[41] + x125*x16 + x125*M[11] + x192*z + x193*z + x194*z + x196*y + x206*x4 + x206*M[5] + x264*M[2] + x34*M[36] + x45*x70 + x45*x71 + x45*x74 + x45*M[26] + y*M[59] + z*M[58] + M[90];
#pragma omp atomic
Ms[91] += x*M[63] + x100*x38 + x100*x39 + x100*x40 + x100*M[17] + x101*x41 + x102*x41 + x103*x41 + x11*x133 + x11*x134 + x11*x137 + x11*M[42] + x121*x13 + x121*M[10] + x196*z + x197*z + x198*z + x199*y + x2*x200 + x200*M[4] + x262*M[1] + x41*M[35] + x45*x76 + x45*x77 + x45*x80 + x45*M[27] + y*M[60] + z*M[59] + M[91];
#pragma omp atomic
Ms[92] += x*M[64] + x100*x42 + x100*M[18] + x104*x41 + x11*x139 + x11*M[43] + x121*M[11] + x199*z + x200*M[5] + x214*M[3] + x262*M[2] + x265*M[0] + x41*M[36] + x43*x96 + x45*x82 + x45*M[28] + x96*M[19] + z*M[60] + M[92];
#pragma omp atomic
Ms[93] += x*M[65] + x107*x34 + x11*x142 + x11*M[44] + x125*M[12] + x144*x267 + x170*x9 + x170*M[9] + x201*y + x206*M[6] + x218*M[4] + x268*M[1] + x34*M[37] + x45*x85 + x45*M[29] + x46*x86 + x86*M[20] + y*M[61] + M[93];
#pragma omp atomic
Ms[94] += x*M[66] + x11*x147 + x11*x148 + x11*x151 + x11*M[45] + x112*x34 + x114*x34 + x116*x34 + x125*x25 + x125*M[13] + x201*z + x202*z + x203*z + x206*x7 + x206*M[7] + x207*y + x218*x4 + x218*M[5] + x268*M[2] + x34*M[38] + x45*x87 + x45*x88 + x45*x89 + x45*M[30] + x49*x86 + x50*x86 + x51*x86 + x86*M[21] + y*M[62] + z*M[61] + M[94];
#pragma omp atomic
Ms[95] += x*M[67] + x107*x41 + x108*x41 + x109*x41 + x11*x153 + x11*x154 + x11*x157 + x11*M[46] + x111*x269 + x118*x34 + x119*x34 + x120*x34 + x121*x19 + x121*M[12] + x125*x30 + x125*M[14] + x177*x9 + x177*M[9] + x200*x5 + x200*M[6] + x206*x8 + x206*M[8] + x207*z + x208*z + x209*z + x210*y + x34*M[39] + x41*M[37] + x45*x90 + x45*x91 + x45*x92 + x45*M[31] + y*M[63] + z*M[62] + M[95];
#pragma omp atomic
Ms[96] += x*M[68] + x11*x160 + x11*x161 + x11*x164 + x11*M[47] + x112*x41 + x113*x41 + x115*x41 + x121*x24 + x121*M[13] + x2*x214 + x200*x6 + x200*M[7] + x210*z + x211*z + x212*z + x213*y + x214*M[4] + x265*M[1] + x41*M[38] + x45*x93 + x45*x94 + x45*x95 + x45*M[32] + x46*x96 + x47*x96 + x48*x96 + x96*M[20] + y*M[64] + z*M[63] + M[96];
#pragma omp atomic
Ms[97] += x*M[69] + x11*x166 + x11*M[48] + x118*x41 + x121*M[14] + x168*x267 + x184*x9 + x184*M[9] + x200*M[8] + x213*z + x214*M[5] + x265*M[2] + x41*M[39] + x45*x97 + x45*M[33] + x49*x96 + x96*M[21] + z*M[64] + M[97];
#pragma omp atomic
Ms[98] += x*M[70] + x0*x251 + x11*x169 + x11*M[49] + x12*x170 + x122*x34 + x125*M[15] + x170*M[10] + x215*y + x218*M[6] + x251*M[3] + x270*M[1] + x34*M[40] + x52*x86 + x86*M[22] + y*M[65] + M[98];
#pragma omp atomic
Ms[99] += x*M[71] + x11*x171 + x11*x172 + x11*x173 + x11*M[50] + x125*x36 + x125*M[16] + x126*x34 + x128*x34 + x130*x34 + x15*x170 + x16*x170 + x17*x170 + x170*M[11] + x215*z + x216*z + x217*z + x218*x7 + x218*M[7] + x219*y + x270*M[2] + x34*M[41] + x55*x86 + x57*x86 + x59*x86 + x86*M[23] + y*M[66] + z*M[65] + M[99];
#pragma omp atomic
Ms[100] += x*M[72] + x0*x254 + x11*x174 + x11*x175 + x11*x176 + x11*M[51] + x12*x177 + x121*x33 + x121*M[15] + x122*x41 + x123*x41 + x124*x41 + x125*x39 + x125*M[17] + x132*x34 + x134*x34 + x136*x34 + x177*M[10] + x218*x8 + x218*M[8] + x219*z + x220*z + x221*z + x222*y + x254*M[3] + x271*M[1] + x34*M[42] + x41*M[40] + x61*x86 + x62*x86 + x63*x86 + x86*M[24] + y*M[67] + z*M[66] + M[100];
#pragma omp atomic
Ms[101] += x*M[73] + x0*x256 + x11*x178 + x11*x179 + x11*x180 + x11*M[52] + x121*x35 + x121*M[16] + x125*x42 + x125*M[18] + x126*x41 + x127*x41 + x129*x41 + x138*x34 + x139*x34 + x140*x34 + x15*x177 + x177*M[11] + x214*x5 + x214*M[6] + x222*z + x223*z + x224*z + x225*y + x256*M[3] + x271*M[2] + x34*M[43] + x41*M[41] + x52*x96 + x53*x96 + x54*x96 + x96*M[22] + y*M[68] + z*M[67] + M[101];
#pragma omp atomic
Ms[102] += x*M[74] + x11*x181 + x11*x182 + x11*x183 + x11*M[53] + x12*x184 + x121*x38 + x121*M[17] + x13*x184 + x132*x41 + x133*x41 + x135*x41 + x14*x184 + x184*M[10] + x214*x6 + x214*M[7] + x225*z + x226*z + x227*z + x228*y + x272*M[1] + x41*M[42] + x55*x96 + x56*x96 + x58*x96 + x96*M[23] + y*M[69] + z*M[68] + M[102];
#pragma omp atomic
Ms[103] += x*M[75] + x0*x259 + x11*x185 + x11*M[54] + x121*M[18] + x138*x41 + x15*x184 + x184*M[11] + x214*M[8] + x228*z + x259*M[3] + x272*M[2] + x41*M[43] + x61*x96 + x96*M[24] + z*M[69] + M[103];
#pragma omp atomic
Ms[104] += x*M[76] + x1*x251 + x141*x34 + x170*x18 + x170*M[12] + x229*y + x251*M[4] + x273*x274 + x34*M[44] + x64*x86 + x86*M[25] + y*M[70] + M[104];
#pragma omp atomic
Ms[105] += x*M[77] + x146*x34 + x148*x34 + x150*x34 + x170*x23 + x170*x25 + x170*x27 + x170*M[13] + x229*z + x230*z + x231*z + x234*y + x251*x3 + x251*x4 + x251*M[5] + x34*M[45] + x69*x86 + x71*x86 + x73*x86 + x86*M[26] + y*M[71] + z*M[70] + M[105];
#pragma omp atomic
Ms[106] += x*M[78] + x1*x254 + x141*x41 + x142*x41 + x143*x41 + x152*x34 + x154*x34 + x156*x34 + x170*x29 + x170*x30 + x170*x31 + x170*M[14] + x177*x18 + x177*M[12] + x234*z + x235*z + x236*z + x237*y + x254*M[4] + x275*M[0] + x34*M[46] + x41*M[44] + x75*x86 + x77*x86 + x79*x86 + x86*M[27] + y*M[72] + z*M[71] + M[106];
#pragma omp atomic
Ms[107] += x*M[79] + x1*x256 + x146*x41 + x147*x41 + x149*x41 + x159*x34 + x161*x34 + x163*x34 + x177*x23 + x177*M[13] + x205*x276 + x237*z + x238*z + x239*z + x241*y + x254*x3 + x254*M[5] + x256*M[4] + x34*M[47] + x41*M[45] + x64*x96 + x65*x96 + x66*x96 + x81*x86 + x82*x86 + x83*x86 + x86*M[28] + x96*M[25] + y*M[73] + z*M[72] + M[107];
#pragma omp atomic
Ms[108] += x*M[80] + x111*x277 + x152*x41 + x153*x41 + x155*x41 + x165*x34 + x166*x34 + x167*x34 + x177*x29 + x177*M[14] + x18*x184 + x184*x19 + x184*x20 + x184*M[12] + x241*z + x242*z + x243*z + x245*y + x256*x3 + x256*M[5] + x34*M[48] + x41*M[46] + x69*x96 + x70*x96 + x72*x96 + x96*M[26] + y*M[74] + z*M[73] + M[108];
#pragma omp atomic
Ms[109] += x*M[81] + x1*x259 + x159*x41 + x160*x41 + x162*x41 + x184*x23 + x184*x24 + x184*x26 + x184*M[13] + x2*x259 + x245*z + x246*z + x247*z + x248*y + x259*M[4] + x41*M[47] + x75*x96 + x76*x96 + x78*x96 + x96*M[27] + y*M[75] + z*M[74] + M[109];
#pragma omp atomic
Ms[110] += x*M[82] + x165*x41 + x184*x29 + x184*M[14] + x248*z + x259*x3 + x259*M[5] + x274*x278 + x41*M[48] + x81*x96 + x96*M[28] + z*M[75] + M[110];
#pragma omp atomic
Ms[111] += x170*M[15] + x251*M[6] + x279*M[1] + x34*M[49] + x86*M[29] + y*M[76] + M[111];
#pragma omp atomic
Ms[112] += x170*x36 + x170*M[16] + x172*x34 + x250*z + x251*x7 + x251*M[7] + x279*M[2] + x34*M[50] + x86*x88 + x86*M[30] + y*M[77] + z*M[76] + M[112];
#pragma omp atomic
Ms[113] += x169*x41 + x170*x39 + x170*M[17] + x175*x34 + x177*M[15] + x251*x8 + x251*M[8] + x252*z + x254*M[6] + x275*M[1] + x34*M[51] + x41*M[49] + x86*x91 + x86*M[31] + y*M[78] + z*M[77] + M[113];
#pragma omp atomic
Ms[114] += x170*x42 + x170*M[18] + x171*x41 + x177*M[16] + x179*x34 + x253*z + x254*M[7] + x256*M[6] + x275*M[2] + x280*M[1] + x34*M[52] + x41*M[50] + x85*x96 + x86*x94 + x86*M[32] + x96*M[29] + y*M[79] + z*M[78] + M[114];
#pragma omp atomic
Ms[115] += x174*x41 + x177*M[17] + x182*x34 + x184*x33 + x184*M[15] + x254*M[8] + x255*z + x256*M[7] + x280*M[2] + x281*M[1] + x34*M[53] + x41*M[51] + x86*x97 + x86*M[33] + x87*x96 + x96*M[30] + y*M[80] + z*M[79] + M[115];
#pragma omp atomic
Ms[116] += x177*M[18] + x178*x41 + x184*x35 + x184*M[16] + x185*x34 + x256*M[8] + x257*z + x259*x5 + x259*M[6] + x281*M[2] + x34*M[54] + x41*M[52] + x90*x96 + x96*M[31] + y*M[81] + z*M[80] + M[116];
#pragma omp atomic
Ms[117] += x181*x41 + x184*x38 + x184*M[17] + x258*z + x259*x6 + x259*M[7] + x282*M[1] + x41*M[53] + x93*x96 + x96*M[32] + y*M[82] + z*M[81] + M[117];
#pragma omp atomic
Ms[118] += x184*M[18] + x259*M[8] + x282*M[2] + x41*M[54] + x96*M[33] + z*M[82] + M[118];

}

void M2L_7(double x, double y, double z, double * M, double * L) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = pow(R, -5);
x3 = 3.0*x2;
x4 = x*y;
x5 = x3*x4;
x6 = x*z;
x7 = x3*x6;
x8 = y*z;
x9 = x3*x8;
x10 = pow(R, -7);
x11 = 15.0*x10;
x12 = x4*z;
x13 = x11*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x3;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x3;
x20 = x14 + x19;
x21 = 9.0*x2;
x22 = x11*x15;
x23 = -x22;
x24 = x*(x21 + x23);
x25 = x23 + x3;
x26 = x25*y;
x27 = x11*x18;
x28 = -x27;
x29 = y*(x21 + x28);
x30 = x25*z;
x31 = z*(x28 + x3);
x32 = 1.0*x;
x33 = x32*(x27 - x3);
x34 = 45.0*x10;
x35 = -x34;
x36 = pow(R, -9);
x37 = x15*x36;
x38 = 105.0*x37;
x39 = x35 + x38;
x40 = x39*x4;
x41 = x39*x6;
x42 = x18*x36;
x43 = 105.0*x42;
x44 = x35 + x43;
x45 = x44*x8;
x46 = -x11;
x47 = x8*(x38 + x46);
x48 = x32*y;
x49 = x44*x48;
x50 = x43 + x46;
x51 = x32*z;
x52 = x50*x51;
x53 = 315.0*x36;
x54 = pow(R, -11);
x55 = 945.0*x54;
x56 = x15*x55;
x57 = x53 - x56;
x58 = x12*x57;
x59 = x18*x55;
x60 = x32*x8;
x61 = x60*(-x53 + x59);
x62 = (x*x*x*x);
x63 = 105.0*x36;
x64 = x62*x63;
x65 = 90.0*x10;
x66 = -x15*x65 + x21 + x64;
x67 = (y*y*y*y);
x68 = x63*x67;
x69 = -x18*x65 + x21 + x68;
x70 = 2.0*x0 - x16 - x19;
x71 = -225.0*x10;
x72 = x55*x62;
x73 = -x72;
x74 = x*(1050.0*x37 + x71 + x73);
x75 = x55*x67;
x76 = -x75;
x77 = y*(1050.0*x42 + x71 + x76);
x78 = x35 + 630.0*x37 + x73;
x79 = x78*y;
x80 = x78*z;
x81 = 630.0*x42;
x82 = x35 + x76 + x81;
x83 = x82*z;
x84 = x32*(x34 + x75 - x81);
x85 = 1575.0*x36;
x86 = pow(R, -13);
x87 = x62*x86;
x88 = 10395.0*x87;
x89 = x15*x54;
x90 = x85 + x88 - 9450.0*x89;
x91 = x4*x90;
x92 = x6*x90;
x93 = 5670.0*x89;
x94 = x53 + x88 - x93;
x95 = x8*x94;
x96 = x67*x86;
x97 = 10395.0*x96;
x98 = x18*x54;
x99 = x85 + x97 - 9450.0*x98;
x100 = x8*x99;
x101 = x48*x99;
x102 = x53 + x97 - 5670.0*x98;
x103 = x102*x51;
x104 = 14175.0*x54;
x105 = 135135.0/pow(R, 15);
x106 = x105*x62;
x107 = x15*x86;
x108 = 103950.0*x107;
x109 = x12*(-x104 - x106 + x108);
x110 = x105*x67;
x111 = x18*x86;
x112 = x60*(x104 + x110 - 103950.0*x111);
x113 = pow(x, 6);
x114 = 10395.0*x86;
x115 = x113*x114;
x116 = -x104*x62 + x115 + 4725.0*x37 + x71;
x117 = pow(y, 6);
x118 = x114*x117;
x119 = -x104*x67 + x118 + 4725.0*x42 + x71;
x120 = 11025.0*x36;
x121 = -x105*x113;
x122 = x*(x120 + x121 + 218295.0*x87 - 99225.0*x89);
x123 = x121 + x85 + 155925.0*x87 - 42525.0*x89;
x124 = x123*y;
x125 = x105*x117;
x126 = -x125;
x127 = y*(x120 + x126 + 218295.0*x96 - 99225.0*x98);
x128 = x123*z;
x129 = 155925.0*x96;
x130 = 42525.0*x98;
x131 = z*(x126 + x129 - x130 + x85);
x132 = x32*(x125 - x129 + x130 - x85);
x133 = x18*x38;
x134 = x133 + x25 + x28;
x135 = x18*x53;
x136 = -x18*x56;
x137 = x*(x135 + x136 + x39);
x138 = x15*x53;
x139 = y*(x136 + x138 + x44);
x140 = z*(x136 + x38 + x50);
x141 = -x24 + x33;
x142 = -x26 - x29;
x143 = -x30 - x31;
x144 = -2835.0*x98;
x145 = x107*x18;
x146 = 10395.0*x145;
x147 = x144 + x146;
x148 = 945.0*x36;
x149 = -2835.0*x89;
x150 = x148 + x149;
x151 = x4*(x147 + x150);
x152 = x6*(x147 + x57);
x153 = x8*(x146 + x149 + x53 - x59);
x154 = x12*(-x105*x15*x18 + 31185.0*x107 + 31185.0*x111 - 8505.0*x54);
x155 = -x40 - x49;
x156 = -x41 - x52;
x157 = -x45 - x47;
x158 = -x58 + x61;
x159 = 105.0*x10;
x160 = -x133 - 12.0*x2;
x161 = x159*x18 + x160 + x22 - x68;
x162 = x15*x159 + x160 + x27 - x64;
x163 = 120.0*x10;
x164 = -x15*x163 - x163*x18 + 210.0*x18*x37 + 24.0*x2 + x64 + x68;
x165 = x15*x97;
x166 = x18*x93;
x167 = -x166;
x168 = x138 + x165 + x167 + x82;
x169 = x18*x88;
x170 = x135 + x167 + x169 + x78;
x171 = 31185.0*x96;
x172 = 62370.0*x145;
x173 = -x110*x15;
x174 = x172 + x173;
x175 = x*(x150 + x171 + x174 - 17010.0*x98);
x176 = x108*x18;
x177 = -x106*x18;
x178 = x*(-x104*x18 + x176 + x177 + x90);
x179 = 31185.0*x87;
x180 = 17010.0*x89;
x181 = x144 + x172 + x177;
x182 = y*(x148 + x179 - x180 + x181);
x183 = y*(-x104*x15 + x173 + x176 + x99);
x184 = z*(x181 + x94);
x185 = z*(x102 + x149 + x174);
x186 = 15120.0*x54;
x187 = -x115;
x188 = -x169;
x189 = 270.0*x10 + x166;
x190 = -x135 + x186*x62 + x187 + x188 + x189 - 5355.0*x37;
x191 = -x118;
x192 = -x165;
x193 = -x138 + x186*x67 + x189 + x191 + x192 - 5355.0*x42;
x194 = -x137;
x195 = x194 - x74;
x196 = x194 + x84;
x197 = -x139;
x198 = x197 - x77;
x199 = x197 - x79;
x200 = -x140;
x201 = x200 - x80;
x202 = x200 - x83;
x203 = -x151;
x204 = x203 - x91;
x205 = -x101 + x203;
x206 = -x152;
x207 = -x103 + x206;
x208 = x206 - x92;
x209 = -x153;
x210 = x209 - x95;
x211 = -x100 + x209;
x212 = -x154;
x213 = -x109 + x212;
x214 = x112 + x212;
x215 = x18*x89;
x216 = -x148*x15 - x148*x18 + x188 + x192 + 11340.0*x215 + x65 + x72 + x75;
x217 = 16065.0*x54;
x218 = -360.0*x10 - x18*x180;
x219 = x118 + 20790.0*x15*x96 + x169 - x217*x67 + x218 + 1260.0*x37 + 6300.0*x42 + x73;
x220 = x115 + x165 + 20790.0*x18*x87 - x217*x62 + x218 + 6300.0*x37 + 1260.0*x42 + x76;
x221 = 2*x139 + x77 + x79;
x222 = 2*x140 + x80 + x83;
x223 = 2*x137 + x74 - x84;
x224 = 17010.0*x54;
x225 = 720.0*x10 - x15*x171 - x179*x18 + x187 + x191 + 34020.0*x215 + x224*x62 + x224*x67 - 7560.0*x37 - 7560.0*x42;
x226 = x100 + 2*x153 + x95;
x227 = -x175;
x228 = x132 + x227;
x229 = -x122;
x230 = -x178;
x231 = x229 + x230;
x232 = -x124;
x233 = -x182;
x234 = x232 + x233;
x235 = -x127;
x236 = -x183;
x237 = x235 + x236;
x238 = -x128;
x239 = -x184;
x240 = x238 + x239;
x241 = -x131;
x242 = -x185;
x243 = x241 + x242;
x244 = x101 + 2*x151 + x91;
x245 = x103 + 2*x152 + x92;
x246 = x109 - x112 + 2*x154;
x247 = x227 + x230;
x248 = x233 + x236;
x249 = x239 + x242;
x250 = x122 + x175 + 2*x178;
x251 = x124 + 2*x182 + x183;
x252 = x127 + x182 + 2*x183;
x253 = x128 + 2*x184 + x185;
x254 = x131 + x184 + 2*x185;
x255 = -x132 + 2*x175 + x178;
x256 = x132 - 3*x175 - 3*x178 + x229;
x257 = -3*x182 - 3*x183 + x232 + x235;
x258 = -3*x184 - 3*x185 + x238 + x241;
#pragma omp atomic
L[0] += -x*x1*M[0] - x1*y*M[1] - x1*z*M[2] + x100*M[77] + x101*M[70] + x103*M[71] + x109*M[87] - x112*M[105] + x116*M[55] + x119*M[76] + x122*M[83] + x124*M[84] + x127*M[111] + x128*M[85] - x13*M[13] + x131*M[112] - x132*M[104] + x134*M[22] + x137*M[37] + x139*M[40] + x140*M[41] + x141*M[14] + x142*M[17] + x143*M[18] + x151*M[61] + x152*M[62] + x153*M[66] + x154*M[94] + x155*M[27] + x156*M[28] + x157*M[32] + x158*M[47] + x161*M[31] + x162*M[24] + x164*M[33] + x168*M[65] + x17*M[3] + x170*M[58] + x175*M[93] + x178*M[86] + x182*M[89] + x183*M[98] + x184*M[90] + x185*M[99] + x190*M[60] + x193*M[78] + x195*M[39] + x196*M[46] + x198*M[51] + x199*M[42] + x20*M[6] + x201*M[43] + x202*M[52] + x204*M[63] + x205*M[72] + x207*M[73] + x208*M[64] + x210*M[68] + x211*M[79] + x213*M[96] + x214*M[107] + x216*M[67] + x219*M[80] + x220*M[69] + x221*M[53] + x222*M[54] + x223*M[48] + x225*M[82] + x226*M[81] + x228*M[106] + x231*M[88] + x234*M[91] + x237*M[113] + x24*M[9] + x240*M[92] + x243*M[114] + x244*M[74] + x245*M[75] + x246*M[109] + x247*M[95] + x248*M[100] + x249*M[101] + x250*M[97] + x251*M[102] + x252*M[115] + x253*M[103] + x254*M[116] + x255*M[108] + x256*M[110] + x257*M[117] + x258*M[118] + x26*M[10] + x29*M[15] + x30*M[11] + x31*M[16] - x33*M[12] + x40*M[20] + x41*M[21] + x45*M[30] + x47*M[23] + x49*M[25] + x5*M[4] + x52*M[26] + x58*M[38] - x61*M[45] + x66*M[19] + x69*M[29] + x7*M[5] + x70*M[8] + x74*M[34] + x77*M[49] + x79*M[35] + x80*M[36] + x83*M[50] - x84*M[44] + x9*M[7] + x91*M[56] + x92*M[57] + x95*M[59];
#pragma omp atomic
L[1] += x101*M[49] + x103*M[50] + x109*M[59] - x112*M[77] + x116*M[34] + x122*M[55] + x124*M[56] + x128*M[57] - x13*M[7] - x132*M[76] + x134*M[12] + x137*M[22] + x139*M[25] + x140*M[26] + x141*M[8] + x151*M[40] + x152*M[41] + x153*M[45] + x154*M[66] + x155*M[17] + x156*M[18] + x158*M[32] + x162*M[14] + x168*M[44] + x17*M[0] + x170*M[37] + x175*M[65] + x178*M[58] + x182*M[61] + x183*M[70] + x184*M[62] + x185*M[71] + x190*M[39] + x195*M[24] + x196*M[31] + x199*M[27] + x201*M[28] + x204*M[42] + x205*M[51] + x207*M[52] + x208*M[43] + x210*M[47] + x213*M[68] + x214*M[79] + x216*M[46] + x220*M[48] + x223*M[33] + x228*M[78] + x231*M[60] + x234*M[63] + x24*M[3] + x240*M[64] + x244*M[53] + x245*M[54] + x246*M[81] + x247*M[67] + x248*M[72] + x249*M[73] + x250*M[69] + x251*M[74] + x253*M[75] + x255*M[80] + x256*M[82] + x26*M[4] + x30*M[5] - x33*M[6] + x40*M[10] + x41*M[11] + x47*M[13] + x49*M[15] + x5*M[1] + x52*M[16] + x58*M[23] - x61*M[30] + x66*M[9] + x7*M[2] + x74*M[19] + x79*M[20] + x80*M[21] - x84*M[29] + x91*M[35] + x92*M[36] + x95*M[38];
#pragma omp atomic
L[2] += x100*M[50] + x101*M[44] + x103*M[45] + x109*M[57] - x112*M[71] + x119*M[49] + x124*M[55] + x127*M[76] - x13*M[5] + x131*M[77] - x132*M[70] + x134*M[10] + x137*M[20] + x139*M[22] + x140*M[23] + x142*M[8] + x151*M[37] + x152*M[38] + x153*M[41] + x154*M[62] + x155*M[14] + x157*M[18] + x158*M[28] + x161*M[17] + x168*M[40] + x170*M[35] + x175*M[61] + x178*M[56] + x182*M[58] + x183*M[65] + x184*M[59] + x185*M[66] + x193*M[51] + x196*M[27] + x198*M[31] + x199*M[24] + x20*M[1] + x202*M[32] + x204*M[39] + x205*M[46] + x207*M[47] + x210*M[43] + x211*M[52] + x213*M[64] + x214*M[73] + x216*M[42] + x219*M[53] + x221*M[33] + x226*M[54] + x228*M[72] + x234*M[60] + x237*M[78] + x243*M[79] + x244*M[48] + x246*M[75] + x247*M[63] + x248*M[67] + x249*M[68] + x251*M[69] + x252*M[80] + x254*M[81] + x255*M[74] + x257*M[82] + x26*M[3] + x29*M[6] + x31*M[7] - x33*M[4] + x40*M[9] + x45*M[16] + x47*M[11] + x49*M[12] + x5*M[0] + x52*M[13] + x58*M[21] - x61*M[26] + x69*M[15] + x77*M[29] + x79*M[19] + x83*M[30] - x84*M[25] + x9*M[2] + x91*M[34] + x95*M[36];
#pragma omp atomic
L[3] += x100*M[49] + x103*M[44] + x109*M[56] - x112*M[70] + x128*M[55] - x13*M[4] + x131*M[76] + x140*M[22] + x141*M[5] + x142*M[7] + x143*M[8] + x152*M[37] + x153*M[40] + x154*M[61] + x155*M[13] + x156*M[14] + x157*M[17] + x158*M[27] + x161*M[16] + x162*M[11] + x164*M[18] + x184*M[58] + x185*M[65] + x190*M[36] + x193*M[50] + x195*M[21] + x196*M[26] + x198*M[30] + x199*M[23] + x201*M[24] + x202*M[31] + x204*M[38] + x205*M[45] + x207*M[46] + x208*M[39] + x210*M[42] + x211*M[51] + x213*M[63] + x214*M[72] + x216*M[41] + x219*M[52] + x220*M[43] + x221*M[32] + x222*M[33] + x223*M[28] + x225*M[54] + x226*M[53] + x228*M[71] + x231*M[57] + x234*M[59] + x237*M[77] + x240*M[60] + x243*M[78] + x244*M[47] + x245*M[48] + x246*M[74] + x247*M[62] + x248*M[66] + x249*M[67] + x250*M[64] + x251*M[68] + x252*M[79] + x253*M[69] + x254*M[80] + x255*M[73] + x256*M[75] + x257*M[81] + x258*M[82] + x30*M[3] + x31*M[6] + x41*M[9] + x45*M[15] + x47*M[10] + x52*M[12] + x58*M[20] - x61*M[25] + x7*M[0] + x70*M[2] + x80*M[19] + x83*M[29] + x9*M[1] + x92*M[34] + x95*M[35];
#pragma omp atomic
L[4] += x109*M[38] + x116*M[19] + x122*M[34] + x124*M[35] + x128*M[36] + x134*M[6] + x137*M[12] + x139*M[15] + x140*M[16] + x151*M[25] + x152*M[26] + x153*M[30] + x154*M[45] + x162*M[8] + x168*M[29] + x170*M[22] + x175*M[44] + x178*M[37] + x182*M[40] + x183*M[49] + x184*M[41] + x185*M[50] + x190*M[24] + x195*M[14] + x199*M[17] + x201*M[18] + x204*M[27] + x208*M[28] + x210*M[32] + x213*M[47] + x216*M[31] + x220*M[33] + x231*M[39] + x234*M[42] + x24*M[0] + x240*M[43] + x247*M[46] + x248*M[51] + x249*M[52] + x250*M[48] + x251*M[53] + x253*M[54] + x26*M[1] + x30*M[2] + x40*M[4] + x41*M[5] + x47*M[7] + x58*M[13] + x66*M[3] + x74*M[9] + x79*M[10] + x80*M[11] + x91*M[20] + x92*M[21] + x95*M[23];
#pragma omp atomic
L[5] += x101*M[29] + x103*M[30] + x109*M[36] - x112*M[50] + x124*M[34] - x13*M[2] - x132*M[49] + x134*M[4] + x137*M[10] + x139*M[12] + x140*M[13] + x151*M[22] + x152*M[23] + x153*M[26] + x154*M[41] + x155*M[8] + x158*M[18] + x168*M[25] + x170*M[20] + x175*M[40] + x178*M[35] + x182*M[37] + x183*M[44] + x184*M[38] + x185*M[45] + x196*M[17] + x199*M[14] + x204*M[24] + x205*M[31] + x207*M[32] + x210*M[28] + x213*M[43] + x214*M[52] + x216*M[27] + x228*M[51] + x234*M[39] + x244*M[33] + x246*M[54] + x247*M[42] + x248*M[46] + x249*M[47] + x251*M[48] + x255*M[53] + x26*M[0] - x33*M[1] + x40*M[3] + x47*M[5] + x49*M[6] + x52*M[7] + x58*M[11] - x61*M[16] + x79*M[9] - x84*M[15] + x91*M[19] + x95*M[21];
#pragma omp atomic
L[6] += x103*M[29] + x109*M[35] - x112*M[49] + x128*M[34] - x13*M[1] + x140*M[12] + x141*M[2] + x152*M[22] + x153*M[25] + x154*M[40] + x155*M[7] + x156*M[8] + x158*M[17] + x162*M[5] + x184*M[37] + x185*M[44] + x190*M[21] + x195*M[11] + x196*M[16] + x199*M[13] + x201*M[14] + x204*M[23] + x205*M[30] + x207*M[31] + x208*M[24] + x210*M[27] + x213*M[42] + x214*M[51] + x216*M[26] + x220*M[28] + x223*M[18] + x228*M[50] + x231*M[36] + x234*M[38] + x240*M[39] + x244*M[32] + x245*M[33] + x246*M[53] + x247*M[41] + x248*M[45] + x249*M[46] + x250*M[43] + x251*M[47] + x253*M[48] + x255*M[52] + x256*M[54] + x30*M[0] + x41*M[3] + x47*M[4] + x52*M[6] + x58*M[10] - x61*M[15] + x80*M[9] + x92*M[19] + x95*M[20];
#pragma omp atomic
L[7] += x100*M[30] + x101*M[25] + x103*M[26] - x112*M[45] + x119*M[29] + x127*M[49] + x131*M[50] - x132*M[44] + x134*M[3] + x137*M[9] + x139*M[10] + x140*M[11] + x151*M[20] + x152*M[21] + x153*M[23] + x154*M[38] + x161*M[8] + x168*M[22] + x170*M[19] + x175*M[37] + x178*M[34] + x182*M[35] + x183*M[40] + x184*M[36] + x185*M[41] + x193*M[31] + x196*M[14] + x198*M[17] + x202*M[18] + x205*M[27] + x207*M[28] + x211*M[32] + x214*M[47] + x216*M[24] + x219*M[33] + x228*M[46] + x237*M[51] + x243*M[52] + x247*M[39] + x248*M[42] + x249*M[43] + x252*M[53] + x254*M[54] + x255*M[48] + x29*M[1] + x31*M[2] - x33*M[0] + x45*M[7] + x49*M[4] + x52*M[5] - x61*M[13] + x69*M[6] + x77*M[15] + x83*M[16] - x84*M[12];
#pragma omp atomic
L[8] += x100*M[29] + x103*M[25] + x109*M[34] - x112*M[44] - x13*M[0] + x131*M[49] + x140*M[10] + x142*M[2] + x152*M[20] + x153*M[22] + x154*M[37] + x155*M[5] + x157*M[8] + x158*M[14] + x161*M[7] + x184*M[35] + x185*M[40] + x193*M[30] + x196*M[13] + x198*M[16] + x199*M[11] + x202*M[17] + x204*M[21] + x205*M[26] + x207*M[27] + x210*M[24] + x211*M[31] + x213*M[39] + x214*M[46] + x216*M[23] + x219*M[32] + x221*M[18] + x226*M[33] + x228*M[45] + x234*M[36] + x237*M[50] + x243*M[51] + x244*M[28] + x246*M[48] + x247*M[38] + x248*M[41] + x249*M[42] + x251*M[43] + x252*M[52] + x254*M[53] + x255*M[47] + x257*M[54] + x31*M[1] + x45*M[6] + x47*M[3] + x52*M[4] + x58*M[9] - x61*M[12] + x83*M[15] + x95*M[19];
#pragma omp atomic
L[9] += x141*M[0] + x142*M[1] + x143*M[2] + x155*M[4] + x156*M[5] + x157*M[7] + x158*M[13] + x161*M[6] + x162*M[3] + x164*M[8] + x190*M[19] + x193*M[29] + x195*M[9] + x196*M[12] + x198*M[15] + x199*M[10] + x201*M[11] + x202*M[16] + x204*M[20] + x205*M[25] + x207*M[26] + x208*M[21] + x210*M[23] + x211*M[30] + x213*M[38] + x214*M[45] + x216*M[22] + x219*M[31] + x220*M[24] + x221*M[17] + x222*M[18] + x223*M[14] + x225*M[33] + x226*M[32] + x228*M[44] + x231*M[34] + x234*M[35] + x237*M[49] + x240*M[36] + x243*M[50] + x244*M[27] + x245*M[28] + x246*M[47] + x247*M[37] + x248*M[40] + x249*M[41] + x250*M[39] + x251*M[42] + x252*M[51] + x253*M[43] + x254*M[52] + x255*M[46] + x256*M[48] + x257*M[53] + x258*M[54];
#pragma omp atomic
L[10] += x109*M[23] + x116*M[9] + x122*M[19] + x124*M[20] + x128*M[21] + x137*M[6] + x151*M[15] + x152*M[16] + x154*M[30] + x170*M[12] + x175*M[29] + x178*M[22] + x182*M[25] + x184*M[26] + x190*M[14] + x195*M[8] + x204*M[17] + x208*M[18] + x213*M[32] + x231*M[24] + x234*M[27] + x240*M[28] + x247*M[31] + x250*M[33] + x40*M[1] + x41*M[2] + x58*M[7] + x66*M[0] + x74*M[3] + x79*M[4] + x80*M[5] + x91*M[10] + x92*M[11] + x95*M[13];
#pragma omp atomic
L[11] += x109*M[21] + x124*M[19] + x134*M[1] + x137*M[4] + x139*M[6] + x140*M[7] + x151*M[12] + x152*M[13] + x153*M[16] + x154*M[26] + x168*M[15] + x170*M[10] + x175*M[25] + x178*M[20] + x182*M[22] + x183*M[29] + x184*M[23] + x185*M[30] + x199*M[8] + x204*M[14] + x210*M[18] + x213*M[28] + x216*M[17] + x234*M[24] + x247*M[27] + x248*M[31] + x249*M[32] + x251*M[33] + x40*M[0] + x47*M[2] + x58*M[5] + x79*M[3] + x91*M[9] + x95*M[11];
#pragma omp atomic
L[12] += x109*M[20] + x128*M[19] + x140*M[6] + x152*M[12] + x153*M[15] + x154*M[25] + x162*M[2] + x184*M[22] + x185*M[29] + x190*M[11] + x195*M[5] + x199*M[7] + x201*M[8] + x204*M[13] + x208*M[14] + x210*M[17] + x213*M[27] + x216*M[16] + x220*M[18] + x231*M[21] + x234*M[23] + x240*M[24] + x247*M[26] + x248*M[30] + x249*M[31] + x250*M[28] + x251*M[32] + x253*M[33] + x41*M[0] + x47*M[1] + x58*M[4] + x80*M[3] + x92*M[9] + x95*M[10];
#pragma omp atomic
L[13] += x101*M[15] + x103*M[16] - x112*M[30] - x132*M[29] + x134*M[0] + x137*M[3] + x139*M[4] + x140*M[5] + x151*M[10] + x152*M[11] + x153*M[13] + x154*M[23] + x168*M[12] + x170*M[9] + x175*M[22] + x178*M[19] + x182*M[20] + x183*M[25] + x184*M[21] + x185*M[26] + x196*M[8] + x205*M[17] + x207*M[18] + x214*M[32] + x216*M[14] + x228*M[31] + x247*M[24] + x248*M[27] + x249*M[28] + x255*M[33] + x49*M[1] + x52*M[2] - x61*M[7] - x84*M[6];
#pragma omp atomic
L[14] += x103*M[15] + x109*M[19] - x112*M[29] + x140*M[4] + x152*M[10] + x153*M[12] + x154*M[22] + x155*M[2] + x158*M[8] + x184*M[20] + x185*M[25] + x196*M[7] + x199*M[5] + x204*M[11] + x205*M[16] + x207*M[17] + x210*M[14] + x213*M[24] + x214*M[31] + x216*M[13] + x228*M[30] + x234*M[21] + x244*M[18] + x246*M[33] + x247*M[23] + x248*M[26] + x249*M[27] + x251*M[28] + x255*M[32] + x47*M[0] + x52*M[1] + x58*M[3] - x61*M[6] + x95*M[9];
#pragma omp atomic
L[15] += x155*M[1] + x156*M[2] + x158*M[7] + x162*M[0] + x190*M[9] + x195*M[3] + x196*M[6] + x199*M[4] + x201*M[5] + x204*M[10] + x205*M[15] + x207*M[16] + x208*M[11] + x210*M[13] + x213*M[23] + x214*M[30] + x216*M[12] + x220*M[14] + x223*M[8] + x228*M[29] + x231*M[19] + x234*M[20] + x240*M[21] + x244*M[17] + x245*M[18] + x246*M[32] + x247*M[22] + x248*M[25] + x249*M[26] + x250*M[24] + x251*M[27] + x253*M[28] + x255*M[31] + x256*M[33];
#pragma omp atomic
L[16] += x100*M[16] + x101*M[12] + x103*M[13] - x112*M[26] + x119*M[15] + x127*M[29] + x131*M[30] - x132*M[25] + x139*M[3] + x151*M[9] + x153*M[11] + x154*M[21] + x168*M[10] + x175*M[20] + x182*M[19] + x183*M[22] + x185*M[23] + x193*M[17] + x198*M[8] + x205*M[14] + x211*M[18] + x214*M[28] + x228*M[27] + x237*M[31] + x243*M[32] + x248*M[24] + x252*M[33] + x45*M[2] + x49*M[0] - x61*M[5] + x69*M[1] + x77*M[6] + x83*M[7] - x84*M[4];
#pragma omp atomic
L[17] += x100*M[15] + x103*M[12] - x112*M[25] + x131*M[29] + x140*M[3] + x152*M[9] + x153*M[10] + x154*M[20] + x161*M[2] + x184*M[19] + x185*M[22] + x193*M[16] + x196*M[5] + x198*M[7] + x202*M[8] + x205*M[13] + x207*M[14] + x211*M[17] + x214*M[27] + x216*M[11] + x219*M[18] + x228*M[26] + x237*M[30] + x243*M[31] + x247*M[21] + x248*M[23] + x249*M[24] + x252*M[32] + x254*M[33] + x255*M[28] + x45*M[1] + x52*M[0] - x61*M[4] + x83*M[6];
#pragma omp atomic
L[18] += x155*M[0] + x157*M[2] + x158*M[5] + x161*M[1] + x193*M[15] + x196*M[4] + x198*M[6] + x199*M[3] + x202*M[7] + x204*M[9] + x205*M[12] + x207*M[13] + x210*M[11] + x211*M[16] + x213*M[21] + x214*M[26] + x216*M[10] + x219*M[17] + x221*M[8] + x226*M[18] + x228*M[25] + x234*M[19] + x237*M[29] + x243*M[30] + x244*M[14] + x246*M[28] + x247*M[20] + x248*M[22] + x249*M[23] + x251*M[24] + x252*M[31] + x254*M[32] + x255*M[27] + x257*M[33];
#pragma omp atomic
L[19] += x156*M[0] + x157*M[1] + x158*M[4] + x164*M[2] + x201*M[3] + x202*M[6] + x207*M[12] + x208*M[9] + x210*M[10] + x211*M[15] + x213*M[20] + x214*M[25] + x219*M[16] + x220*M[11] + x221*M[7] + x222*M[8] + x223*M[5] + x225*M[18] + x226*M[17] + x240*M[19] + x243*M[29] + x244*M[13] + x245*M[14] + x246*M[27] + x249*M[22] + x250*M[21] + x251*M[23] + x252*M[30] + x253*M[24] + x254*M[31] + x255*M[26] + x256*M[28] + x257*M[32] + x258*M[33];
#pragma omp atomic
L[20] += x109*M[13] + x116*M[3] + x122*M[9] + x124*M[10] + x128*M[11] + x170*M[6] + x178*M[12] + x182*M[15] + x184*M[16] + x190*M[8] + x231*M[14] + x234*M[17] + x240*M[18] + x74*M[0] + x79*M[1] + x80*M[2] + x91*M[4] + x92*M[5] + x95*M[7];
#pragma omp atomic
L[21] += x109*M[11] + x124*M[9] + x137*M[1] + x151*M[6] + x152*M[7] + x154*M[16] + x170*M[4] + x175*M[15] + x178*M[10] + x182*M[12] + x184*M[13] + x204*M[8] + x213*M[18] + x234*M[14] + x247*M[17] + x58*M[2] + x79*M[0] + x91*M[3] + x95*M[5];
#pragma omp atomic
L[22] += x109*M[10] + x128*M[9] + x152*M[6] + x154*M[15] + x184*M[12] + x190*M[5] + x195*M[2] + x204*M[7] + x208*M[8] + x213*M[17] + x231*M[11] + x234*M[13] + x240*M[14] + x247*M[16] + x250*M[18] + x58*M[1] + x80*M[0] + x92*M[3] + x95*M[4];
#pragma omp atomic
L[23] += x137*M[0] + x139*M[1] + x140*M[2] + x151*M[4] + x152*M[5] + x153*M[7] + x154*M[13] + x168*M[6] + x170*M[3] + x175*M[12] + x178*M[9] + x182*M[10] + x183*M[15] + x184*M[11] + x185*M[16] + x216*M[8] + x247*M[14] + x248*M[17] + x249*M[18];
#pragma omp atomic
L[24] += x109*M[9] + x140*M[1] + x152*M[4] + x153*M[6] + x154*M[12] + x184*M[10] + x185*M[15] + x199*M[2] + x204*M[5] + x210*M[8] + x213*M[14] + x216*M[7] + x234*M[11] + x247*M[13] + x248*M[16] + x249*M[17] + x251*M[18] + x58*M[0] + x95*M[3];
#pragma omp atomic
L[25] += x190*M[3] + x195*M[0] + x199*M[1] + x201*M[2] + x204*M[4] + x208*M[5] + x210*M[7] + x213*M[13] + x216*M[6] + x220*M[8] + x231*M[9] + x234*M[10] + x240*M[11] + x247*M[12] + x248*M[15] + x249*M[16] + x250*M[14] + x251*M[17] + x253*M[18];
#pragma omp atomic
L[26] += x101*M[6] + x103*M[7] - x112*M[16] - x132*M[15] + x139*M[0] + x151*M[3] + x153*M[5] + x154*M[11] + x168*M[4] + x175*M[10] + x182*M[9] + x183*M[12] + x185*M[13] + x205*M[8] + x214*M[18] + x228*M[17] + x248*M[14] - x61*M[2] - x84*M[1];
#pragma omp atomic
L[27] += x103*M[6] - x112*M[15] + x140*M[0] + x152*M[3] + x153*M[4] + x154*M[10] + x184*M[9] + x185*M[12] + x196*M[2] + x205*M[7] + x207*M[8] + x214*M[17] + x216*M[5] + x228*M[16] + x247*M[11] + x248*M[13] + x249*M[14] + x255*M[18] - x61*M[1];
#pragma omp atomic
L[28] += x158*M[2] + x196*M[1] + x199*M[0] + x204*M[3] + x205*M[6] + x207*M[7] + x210*M[5] + x213*M[11] + x214*M[16] + x216*M[4] + x228*M[15] + x234*M[9] + x244*M[8] + x246*M[18] + x247*M[10] + x248*M[12] + x249*M[13] + x251*M[14] + x255*M[17];
#pragma omp atomic
L[29] += x158*M[1] + x201*M[0] + x207*M[6] + x208*M[3] + x210*M[4] + x213*M[10] + x214*M[15] + x220*M[5] + x223*M[2] + x240*M[9] + x244*M[7] + x245*M[8] + x246*M[17] + x249*M[12] + x250*M[11] + x251*M[13] + x253*M[14] + x255*M[16] + x256*M[18];
#pragma omp atomic
L[30] += x100*M[7] + x101*M[4] + x103*M[5] - x112*M[13] + x119*M[6] + x127*M[15] + x131*M[16] - x132*M[12] + x168*M[3] + x175*M[9] + x183*M[10] + x185*M[11] + x193*M[8] + x228*M[14] + x237*M[17] + x243*M[18] + x77*M[1] + x83*M[2] - x84*M[0];
#pragma omp atomic
L[31] += x100*M[6] + x103*M[4] - x112*M[12] + x131*M[15] + x153*M[3] + x154*M[9] + x185*M[10] + x193*M[7] + x198*M[2] + x205*M[5] + x211*M[8] + x214*M[14] + x228*M[13] + x237*M[16] + x243*M[17] + x248*M[11] + x252*M[18] - x61*M[0] + x83*M[1];
#pragma omp atomic
L[32] += x193*M[6] + x196*M[0] + x198*M[1] + x202*M[2] + x205*M[4] + x207*M[5] + x211*M[7] + x214*M[13] + x216*M[3] + x219*M[8] + x228*M[12] + x237*M[15] + x243*M[16] + x247*M[9] + x248*M[10] + x249*M[11] + x252*M[17] + x254*M[18] + x255*M[14];
#pragma omp atomic
L[33] += x158*M[0] + x202*M[1] + x207*M[4] + x210*M[3] + x211*M[6] + x213*M[9] + x214*M[12] + x219*M[7] + x221*M[2] + x226*M[8] + x243*M[15] + x244*M[5] + x246*M[14] + x249*M[10] + x251*M[11] + x252*M[16] + x254*M[17] + x255*M[13] + x257*M[18];
#pragma omp atomic
L[34] += x219*M[6] + x220*M[3] + x221*M[1] + x222*M[2] + x223*M[0] + x225*M[8] + x226*M[7] + x244*M[4] + x245*M[5] + x246*M[13] + x250*M[9] + x251*M[10] + x252*M[15] + x253*M[11] + x254*M[16] + x255*M[12] + x256*M[14] + x257*M[17] + x258*M[18];
#pragma omp atomic
L[35] += x109*M[7] + x116*M[0] + x122*M[3] + x124*M[4] + x128*M[5] + x178*M[6] + x231*M[8] + x91*M[1] + x92*M[2];
#pragma omp atomic
L[36] += x109*M[5] + x124*M[3] + x170*M[1] + x178*M[4] + x182*M[6] + x184*M[7] + x234*M[8] + x91*M[0] + x95*M[2];
#pragma omp atomic
L[37] += x109*M[4] + x128*M[3] + x184*M[6] + x190*M[2] + x231*M[5] + x234*M[7] + x240*M[8] + x92*M[0] + x95*M[1];
#pragma omp atomic
L[38] += x151*M[1] + x152*M[2] + x154*M[7] + x170*M[0] + x175*M[6] + x178*M[3] + x182*M[4] + x184*M[5] + x247*M[8];
#pragma omp atomic
L[39] += x109*M[3] + x152*M[1] + x154*M[6] + x184*M[4] + x204*M[2] + x213*M[8] + x234*M[5] + x247*M[7] + x95*M[0];
#pragma omp atomic
L[40] += x190*M[0] + x204*M[1] + x208*M[2] + x213*M[7] + x231*M[3] + x234*M[4] + x240*M[5] + x247*M[6] + x250*M[8];
#pragma omp atomic
L[41] += x151*M[0] + x153*M[2] + x154*M[5] + x168*M[1] + x175*M[4] + x182*M[3] + x183*M[6] + x185*M[7] + x248*M[8];
#pragma omp atomic
L[42] += x152*M[0] + x153*M[1] + x154*M[4] + x184*M[3] + x185*M[6] + x216*M[2] + x247*M[5] + x248*M[7] + x249*M[8];
#pragma omp atomic
L[43] += x204*M[0] + x210*M[2] + x213*M[5] + x216*M[1] + x234*M[3] + x247*M[4] + x248*M[6] + x249*M[7] + x251*M[8];
#pragma omp atomic
L[44] += x208*M[0] + x210*M[1] + x213*M[4] + x220*M[2] + x240*M[3] + x249*M[6] + x250*M[5] + x251*M[7] + x253*M[8];
#pragma omp atomic
L[45] += x101*M[1] + x103*M[2] - x112*M[7] - x132*M[6] + x168*M[0] + x175*M[3] + x183*M[4] + x185*M[5] + x228*M[8];
#pragma omp atomic
L[46] += x103*M[1] - x112*M[6] + x153*M[0] + x154*M[3] + x185*M[4] + x205*M[2] + x214*M[8] + x228*M[7] + x248*M[5];
#pragma omp atomic
L[47] += x205*M[1] + x207*M[2] + x214*M[7] + x216*M[0] + x228*M[6] + x247*M[3] + x248*M[4] + x249*M[5] + x255*M[8];
#pragma omp atomic
L[48] += x207*M[1] + x210*M[0] + x213*M[3] + x214*M[6] + x244*M[2] + x246*M[8] + x249*M[4] + x251*M[5] + x255*M[7];
#pragma omp atomic
L[49] += x220*M[0] + x244*M[1] + x245*M[2] + x246*M[7] + x250*M[3] + x251*M[4] + x253*M[5] + x255*M[6] + x256*M[8];
#pragma omp atomic
L[50] += x100*M[2] + x101*M[0] - x112*M[5] + x119*M[1] + x127*M[6] + x131*M[7] - x132*M[4] + x183*M[3] + x237*M[8];
#pragma omp atomic
L[51] += x100*M[1] + x103*M[0] - x112*M[4] + x131*M[6] + x185*M[3] + x193*M[2] + x228*M[5] + x237*M[7] + x243*M[8];
#pragma omp atomic
L[52] += x193*M[1] + x205*M[0] + x211*M[2] + x214*M[5] + x228*M[4] + x237*M[6] + x243*M[7] + x248*M[3] + x252*M[8];
#pragma omp atomic
L[53] += x207*M[0] + x211*M[1] + x214*M[4] + x219*M[2] + x243*M[6] + x249*M[3] + x252*M[7] + x254*M[8] + x255*M[5];
#pragma omp atomic
L[54] += x219*M[1] + x226*M[2] + x244*M[0] + x246*M[5] + x251*M[3] + x252*M[6] + x254*M[7] + x255*M[4] + x257*M[8];
#pragma omp atomic
L[55] += x225*M[2] + x226*M[1] + x245*M[0] + x246*M[4] + x253*M[3] + x254*M[6] + x256*M[5] + x257*M[7] + x258*M[8];
#pragma omp atomic
L[56] += x122*M[0] + x124*M[1] + x128*M[2];
#pragma omp atomic
L[57] += x109*M[2] + x124*M[0] + x178*M[1];
#pragma omp atomic
L[58] += x109*M[1] + x128*M[0] + x231*M[2];
#pragma omp atomic
L[59] += x178*M[0] + x182*M[1] + x184*M[2];
#pragma omp atomic
L[60] += x109*M[0] + x184*M[1] + x234*M[2];
#pragma omp atomic
L[61] += x231*M[0] + x234*M[1] + x240*M[2];
#pragma omp atomic
L[62] += x154*M[2] + x175*M[1] + x182*M[0];
#pragma omp atomic
L[63] += x154*M[1] + x184*M[0] + x247*M[2];
#pragma omp atomic
L[64] += x213*M[2] + x234*M[0] + x247*M[1];
#pragma omp atomic
L[65] += x213*M[1] + x240*M[0] + x250*M[2];
#pragma omp atomic
L[66] += x175*M[0] + x183*M[1] + x185*M[2];
#pragma omp atomic
L[67] += x154*M[0] + x185*M[1] + x248*M[2];
#pragma omp atomic
L[68] += x247*M[0] + x248*M[1] + x249*M[2];
#pragma omp atomic
L[69] += x213*M[0] + x249*M[1] + x251*M[2];
#pragma omp atomic
L[70] += x250*M[0] + x251*M[1] + x253*M[2];
#pragma omp atomic
L[71] += -x112*M[2] - x132*M[1] + x183*M[0];
#pragma omp atomic
L[72] += -x112*M[1] + x185*M[0] + x228*M[2];
#pragma omp atomic
L[73] += x214*M[2] + x228*M[1] + x248*M[0];
#pragma omp atomic
L[74] += x214*M[1] + x249*M[0] + x255*M[2];
#pragma omp atomic
L[75] += x246*M[2] + x251*M[0] + x255*M[1];
#pragma omp atomic
L[76] += x246*M[1] + x253*M[0] + x256*M[2];
#pragma omp atomic
L[77] += x127*M[1] + x131*M[2] - x132*M[0];
#pragma omp atomic
L[78] += -x112*M[0] + x131*M[1] + x237*M[2];
#pragma omp atomic
L[79] += x228*M[0] + x237*M[1] + x243*M[2];
#pragma omp atomic
L[80] += x214*M[0] + x243*M[1] + x252*M[2];
#pragma omp atomic
L[81] += x252*M[1] + x254*M[2] + x255*M[0];
#pragma omp atomic
L[82] += x246*M[0] + x254*M[1] + x257*M[2];
#pragma omp atomic
L[83] += x256*M[0] + x257*M[1] + x258*M[2];

}

void L2L_7(double x, double y, double z, double * L, double * Ls) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
x3 = z*L[14];
x4 = x3*y;
x5 = (x*x);
x6 = (1.0/2.0)*x5;
x7 = (x*x*x);
x8 = (1.0/6.0)*x7;
x9 = (x*x*x*x);
x10 = (1.0/24.0)*x9;
x11 = (1.0/120.0)*pow(x, 5);
x12 = (y*y);
x13 = (1.0/2.0)*x12;
x14 = (y*y*y);
x15 = (1.0/6.0)*x14;
x16 = (y*y*y*y);
x17 = (1.0/24.0)*x16;
x18 = (1.0/120.0)*pow(y, 5);
x19 = (z*z);
x20 = (1.0/2.0)*x19;
x21 = (z*z*z);
x22 = (1.0/6.0)*x21;
x23 = (z*z*z*z);
x24 = (1.0/24.0)*x23;
x25 = (1.0/120.0)*pow(z, 5);
x26 = x*L[13];
x27 = x*L[26];
x28 = x*L[45];
x29 = x*L[71];
x30 = x*L[15];
x31 = x*L[29];
x32 = x*L[49];
x33 = x*L[76];
x34 = y*L[11];
x35 = z*L[12];
x36 = y*L[21];
x37 = z*L[22];
x38 = y*L[36];
x39 = z*L[37];
x40 = y*L[57];
x41 = z*L[58];
x42 = y*L[18];
x43 = y*L[33];
x44 = y*L[54];
x45 = y*L[82];
x46 = z*L[17];
x47 = z*L[31];
x48 = z*L[51];
x49 = z*L[78];
x50 = y*L[28];
x51 = x*x50;
x52 = y*L[48];
x53 = x*x52;
x54 = y*L[75];
x55 = x*x54;
x56 = z*L[27];
x57 = x*x56;
x58 = z*L[46];
x59 = x*x58;
x60 = z*L[72];
x61 = x*x60;
x62 = z*L[24];
x63 = x62*y;
x64 = z*L[39];
x65 = x64*y;
x66 = z*L[60];
x67 = x66*y;
x68 = (1.0/4.0)*x5;
x69 = x12*x68;
x70 = (1.0/12.0)*x5;
x71 = x14*x70;
x72 = (1.0/48.0)*x5;
x73 = x19*x68;
x74 = x21*x70;
x75 = (1.0/12.0)*x7;
x76 = x12*x75;
x77 = (1.0/36.0)*x7;
x78 = x19*x75;
x79 = (1.0/48.0)*x9;
x80 = x12*x19;
x81 = (1.0/4.0)*x80;
x82 = (1.0/12.0)*x12*x21;
x83 = (1.0/12.0)*x14*x19;
x84 = x*L[47];
x85 = x*L[74];
x86 = x*L[73];
x87 = y*L[43];
x88 = y*L[69];
x89 = z*L[42];
x90 = z*L[67];
x91 = y*L[64];
x92 = z*L[63];
x93 = x*L[23];
x94 = x*L[41];
x95 = x*L[66];
x96 = x*L[25];
x97 = x*L[44];
x98 = x*L[70];
x99 = x*x87;
x100 = x*x88;
x101 = x*x89;
x102 = x*x90;
x103 = x*L[68];
x104 = y*L[13];
x105 = x56*y;
x106 = x*L[28];
x107 = x*L[48];
x108 = x*L[75];
x109 = y*L[23];
x110 = y*L[38];
x111 = y*L[59];
x112 = y*L[32];
x113 = y*L[53];
x114 = y*L[81];
x115 = y*L[47];
x116 = x*x115;
x117 = y*L[74];
x118 = x*x117;
x119 = x89*y;
x120 = x92*y;
x121 = y*L[68];
x122 = y*L[14];
x123 = z*L[15];
x124 = z*L[18];
x125 = z*L[28];
x126 = x125*y;
x127 = x*L[27];
x128 = x*L[46];
x129 = x*L[72];
x130 = y*L[24];
x131 = z*L[25];
x132 = y*L[39];
x133 = z*L[40];
x134 = y*L[60];
x135 = z*L[61];
x136 = z*L[32];
x137 = z*L[52];
x138 = z*L[79];
x139 = z*L[47];
x140 = x*x139;
x141 = z*L[73];
x142 = x*x141;
x143 = z*L[43];
x144 = x143*y;
x145 = z*L[64];
x146 = x145*y;
x147 = z*L[68];
x148 = x*L[38];
x149 = x*L[62];
x150 = x*L[40];
x151 = x*L[65];
x152 = x*x91;
x153 = x*x92;
x154 = x*L[43];
x155 = x*L[69];
x156 = x*x121;
x157 = x*L[42];
x158 = x*L[67];
x159 = x*x147;
x160 = y*L[26];
x161 = x58*y;
x162 = y*L[41];
x163 = y*L[62];
x164 = y*L[52];
x165 = y*L[80];
x166 = y*L[73];
x167 = x*x166;
x168 = x90*y;
x169 = y*L[27];
x170 = x139*y;
x171 = y*L[42];
x172 = y*L[63];
x173 = x147*y;
x174 = z*L[29];
x175 = z*L[33];
x176 = z*L[48];
x177 = x176*y;
x178 = z*L[44];
x179 = z*L[65];
x180 = z*L[53];
x181 = z*L[80];
x182 = z*L[74];
x183 = x*x182;
x184 = z*L[69];
x185 = x184*y;
x186 = x*L[59];
x187 = x*L[61];
x188 = x*L[64];
x189 = x*L[63];
x190 = y*L[45];
x191 = x60*y;
x192 = y*L[66];
x193 = y*L[79];
x194 = y*L[46];
x195 = x141*y;
x196 = y*L[67];
x197 = x182*y;
x198 = z*L[49];
x199 = z*L[54];
x200 = z*L[75];
x201 = x200*y;
x202 = z*L[70];
x203 = z*L[81];
x204 = y*L[71];
x205 = y*L[72];
x206 = z*L[76];
x207 = z*L[82];
#pragma omp atomic
Ls[0] += (1.0/720.0)*pow(x, 6)*L[56] + x*x0 + x*x1 + x*x4 + x*L[1] + x10*x38 + x10*x39 + x10*x67 + x10*L[20] + x11*x40 + x11*x41 + x11*L[35] + (1.0/48.0)*x12*x23*L[81] + x12*x79*L[59] + x13*x26 + x13*x46 + x13*x57 + x13*L[7] + (1.0/36.0)*x14*x21*L[80] + x14*x77*L[62] + x15*x27 + x15*x47 + x15*x59 + x15*L[16] + (1.0/48.0)*x16*x19*L[79] + x16*x72*L[66] + x17*x28 + x17*x48 + x17*x61 + x17*L[30] + x18*x29 + x18*x49 + x18*L[50] + x19*x79*L[61] + x2*y + x20*x30 + x20*x42 + x20*x51 + x20*L[9] + x21*x77*L[65] + x22*x31 + x22*x43 + x22*x53 + x22*L[19] + x23*x72*L[70] + x24*x32 + x24*x44 + x24*x55 + x24*L[34] + x25*x33 + x25*x45 + x25*L[55] + x34*x6 + x35*x6 + x36*x8 + x37*x8 + (1.0/8.0)*x5*x80*L[68] + x6*x63 + x6*L[4] + x65*x8 + x69*x89 + x69*L[23] + x71*x90 + x71*L[41] + x73*x87 + x73*L[25] + x74*x88 + x74*L[44] + x76*x92 + x76*L[38] + x78*x91 + x78*L[40] + x8*L[10] + x81*x84 + x81*L[32] + x82*x85 + x82*L[53] + x83*x86 + x83*L[52] + (1.0/720.0)*pow(y, 6)*L[77] + y*L[2] + (1.0/720.0)*pow(z, 6)*L[83] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*x34 + x*x35 + x*x63 + x*L[4] + x0 + x1 + x10*x40 + x10*x41 + x10*L[35] + x100*x22 + x101*x13 + x102*x15 + x103*x81 + x11*L[56] + x13*x56 + x13*x93 + x13*L[13] + x15*x58 + x15*x94 + x15*L[26] + x17*x60 + x17*x95 + x17*L[45] + x18*L[71] + x20*x50 + x20*x96 + x20*x99 + x20*L[15] + x22*x52 + x22*x97 + x22*L[29] + x24*x54 + x24*x98 + x24*L[49] + x25*L[76] + x36*x6 + x37*x6 + x38*x8 + x39*x8 + x4 + x6*x65 + x6*L[10] + x67*x8 + x69*x92 + x69*L[38] + x71*L[62] + x73*x91 + x73*L[40] + x74*L[65] + x76*L[59] + x78*L[61] + x8*L[20] + x81*L[47] + x82*L[74] + x83*L[73] + L[1];
#pragma omp atomic
Ls[2] += x*x104 + x*x105 + x*x3 + x*L[5] + x10*x111 + x10*x66 + x10*L[36] + x106*x20 + x107*x22 + x108*x24 + x109*x6 + x11*L[57] + x110*x8 + x112*x20 + x113*x22 + x114*x24 + x116*x20 + x118*x22 + x119*x6 + x120*x8 + x121*x73 + x13*x27 + x13*x47 + x13*x59 + x13*L[16] + x15*x28 + x15*x48 + x15*x61 + x15*L[30] + x17*x29 + x17*x49 + x17*L[50] + x18*L[77] + x2 + x20*L[18] + x22*L[33] + x24*L[54] + x25*L[82] + x46*y + x6*x62 + x6*L[11] + x64*x8 + x69*x90 + x69*L[41] + x71*L[66] + x73*L[43] + x74*L[69] + x76*L[62] + x78*L[64] + x8*L[21] + x81*x86 + x81*L[52] + x82*L[80] + x83*L[79] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*x122 + x*x123 + x*x126 + x*L[6] + x10*x134 + x10*x135 + x10*L[37] + x11*L[58] + x124*y + x127*x13 + x128*x15 + x129*x17 + x13*x136 + x13*x140 + x13*L[17] + x130*x6 + x131*x6 + x132*x8 + x133*x8 + x137*x15 + x138*x17 + x142*x15 + x144*x6 + x146*x8 + x147*x69 + x15*L[31] + x17*L[51] + x18*L[78] + x20*x31 + x20*x43 + x20*x53 + x20*L[19] + x22*x32 + x22*x44 + x22*x55 + x22*L[34] + x24*x33 + x24*x45 + x24*L[55] + x25*L[83] + x6*L[12] + x69*L[42] + x71*L[67] + x73*x88 + x73*L[44] + x74*L[70] + x76*L[63] + x78*L[65] + x8*L[22] + x81*x85 + x81*L[53] + x82*L[81] + x83*L[80] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += x*x36 + x*x37 + x*x65 + x*L[10] + x10*L[56] + x13*x148 + x13*x153 + x13*x89 + x13*L[23] + x149*x15 + x15*x90 + x15*L[41] + x150*x20 + x151*x22 + x152*x20 + x17*L[66] + x20*x87 + x20*L[25] + x22*x88 + x22*L[44] + x24*L[70] + x34 + x35 + x38*x6 + x39*x6 + x40*x8 + x41*x8 + x6*x67 + x6*L[20] + x63 + x69*L[59] + x73*L[61] + x8*L[35] + x81*L[68] + L[4];
#pragma omp atomic
Ls[5] += x*x109 + x*x119 + x*x62 + x*L[11] + x10*L[57] + x102*x13 + x104 + x105 + x110*x6 + x111*x8 + x115*x20 + x117*x22 + x120*x6 + x13*x58 + x13*x94 + x13*L[26] + x15*x60 + x15*x95 + x15*L[45] + x154*x20 + x155*x22 + x156*x20 + x17*L[71] + x20*L[28] + x22*L[48] + x24*L[75] + x3 + x6*x64 + x6*L[21] + x66*x8 + x69*L[62] + x73*L[64] + x8*L[36] + x81*L[73] + L[5];
#pragma omp atomic
Ls[6] += x*x130 + x*x131 + x*x144 + x*L[12] + x10*L[58] + x100*x20 + x122 + x123 + x126 + x13*x139 + x13*x157 + x13*x159 + x13*L[27] + x132*x6 + x133*x6 + x134*x8 + x135*x8 + x141*x15 + x146*x6 + x15*x158 + x15*L[46] + x17*L[72] + x20*x52 + x20*x97 + x20*L[29] + x22*x54 + x22*x98 + x22*L[49] + x24*L[76] + x6*L[22] + x69*L[63] + x73*L[65] + x8*L[37] + x81*L[74] + L[6];
#pragma omp atomic
Ls[7] += x*x160 + x*x161 + x10*L[59] + x13*x28 + x13*x48 + x13*x61 + x13*L[30] + x15*x29 + x15*x49 + x15*L[50] + x162*x6 + x163*x8 + x164*x20 + x165*x22 + x167*x20 + x168*x6 + x17*L[77] + x20*x84 + x20*L[32] + x22*x85 + x22*L[53] + x24*L[81] + x26 + x46 + x47*y + x57 + x6*x89 + x6*L[23] + x69*L[66] + x73*L[68] + x8*x92 + x8*L[38] + x81*L[79] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += x*x125 + x*x169 + x*x170 + x*L[14] + x10*L[60] + x107*x20 + x108*x22 + x113*x20 + x114*x22 + x118*x20 + x124 + x128*x13 + x129*x15 + x13*x137 + x13*x142 + x13*L[31] + x136*y + x138*x15 + x143*x6 + x145*x8 + x15*L[51] + x17*L[78] + x171*x6 + x172*x8 + x173*x6 + x20*L[33] + x22*L[54] + x24*L[82] + x6*L[24] + x69*L[67] + x73*L[69] + x8*L[39] + x81*L[80] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += x*x174 + x*x177 + x10*L[61] + x13*x180 + x13*x183 + x13*x84 + x13*L[32] + x15*x181 + x15*x86 + x15*L[52] + x17*L[79] + x175*y + x178*x6 + x179*x8 + x185*x6 + x20*x32 + x20*x44 + x20*x55 + x20*L[34] + x22*x33 + x22*x45 + x22*L[55] + x24*L[83] + x30 + x42 + x51 + x6*x87 + x6*L[25] + x69*L[68] + x73*L[70] + x8*x91 + x8*L[40] + x81*L[81] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += x*x38 + x*x39 + x*x67 + x*L[20] + x13*x186 + x13*x92 + x13*L[38] + x15*L[62] + x187*x20 + x20*x91 + x20*L[40] + x22*L[65] + x36 + x37 + x40*x6 + x41*x6 + x6*L[35] + x65 + x8*L[56] + L[10];
#pragma omp atomic
Ls[11] += x*x110 + x*x120 + x*x64 + x*L[21] + x109 + x111*x6 + x119 + x121*x20 + x13*x149 + x13*x90 + x13*L[41] + x15*L[66] + x188*x20 + x20*L[43] + x22*L[69] + x6*x66 + x6*L[36] + x62 + x8*L[57] + L[11];
#pragma omp atomic
Ls[12] += x*x132 + x*x133 + x*x146 + x*L[22] + x13*x147 + x13*x189 + x13*L[42] + x130 + x131 + x134*x6 + x135*x6 + x144 + x15*L[67] + x151*x20 + x20*x88 + x20*L[44] + x22*L[70] + x6*L[37] + x8*L[58] + L[12];
#pragma omp atomic
Ls[13] += x*x162 + x*x168 + x101 + x103*x20 + x13*x60 + x13*x95 + x13*L[45] + x15*L[71] + x160 + x161 + x163*x6 + x166*x20 + x20*L[47] + x22*L[74] + x56 + x6*x92 + x6*L[38] + x8*L[59] + x93 + L[13];
#pragma omp atomic
Ls[14] += x*x143 + x*x171 + x*x173 + x*L[24] + x117*x20 + x125 + x13*x141 + x13*x158 + x13*L[46] + x145*x6 + x15*L[72] + x155*x20 + x169 + x170 + x172*x6 + x20*L[48] + x22*L[75] + x6*L[39] + x8*L[60] + L[14];
#pragma omp atomic
Ls[15] += x*x178 + x*x185 + x103*x13 + x13*x182 + x13*L[47] + x15*L[73] + x174 + x177 + x179*x6 + x20*x54 + x20*x98 + x20*L[49] + x22*L[76] + x50 + x6*x91 + x6*L[40] + x8*L[61] + x96 + x99 + L[15];
#pragma omp atomic
Ls[16] += x*x190 + x*x191 + x13*x29 + x13*x49 + x13*L[50] + x15*L[77] + x192*x6 + x193*x20 + x20*x86 + x20*L[52] + x22*L[80] + x27 + x47 + x48*y + x59 + x6*x90 + x6*L[41] + x8*L[62] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += x*x194 + x*x195 + x127 + x129*x13 + x13*x138 + x13*L[51] + x136 + x137*y + x140 + x147*x6 + x15*L[78] + x165*x20 + x196*x6 + x20*x85 + x20*L[53] + x22*L[81] + x6*L[42] + x8*L[63] + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += x*x176 + x*x197 + x106 + x108*x20 + x112 + x114*x20 + x116 + x121*x6 + x13*x181 + x13*x86 + x13*L[52] + x15*L[79] + x175 + x180*y + x184*x6 + x20*L[54] + x22*L[82] + x6*L[43] + x8*L[64] + L[18];
#pragma omp atomic
Ls[19] += x*x198 + x*x201 + x13*x203 + x13*x85 + x13*L[53] + x15*L[80] + x199*y + x20*x33 + x20*x45 + x20*L[55] + x202*x6 + x22*L[83] + x31 + x43 + x53 + x6*x88 + x6*L[44] + x8*L[65] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += x*x40 + x*x41 + x*L[35] + x13*L[59] + x20*L[61] + x38 + x39 + x6*L[56] + x67 + L[20];
#pragma omp atomic
Ls[21] += x*x111 + x*x66 + x*L[36] + x110 + x120 + x13*L[62] + x20*L[64] + x6*L[57] + x64 + L[21];
#pragma omp atomic
Ls[22] += x*x134 + x*x135 + x*L[37] + x13*L[63] + x132 + x133 + x146 + x20*L[65] + x6*L[58] + L[22];
#pragma omp atomic
Ls[23] += x*x163 + x13*L[66] + x148 + x153 + x162 + x168 + x20*L[68] + x6*L[59] + x89 + L[23];
#pragma omp atomic
Ls[24] += x*x145 + x*x172 + x*L[39] + x13*L[67] + x143 + x171 + x173 + x20*L[69] + x6*L[60] + L[24];
#pragma omp atomic
Ls[25] += x*x179 + x13*L[68] + x150 + x152 + x178 + x185 + x20*L[70] + x6*L[61] + x87 + L[25];
#pragma omp atomic
Ls[26] += x*x192 + x102 + x13*L[71] + x190 + x191 + x20*L[73] + x58 + x6*L[62] + x94 + L[26];
#pragma omp atomic
Ls[27] += x*x196 + x13*L[72] + x139 + x157 + x159 + x194 + x195 + x20*L[74] + x6*L[63] + L[27];
#pragma omp atomic
Ls[28] += x*x184 + x115 + x13*L[73] + x154 + x156 + x176 + x197 + x20*L[75] + x6*L[64] + L[28];
#pragma omp atomic
Ls[29] += x*x202 + x100 + x13*L[74] + x198 + x20*L[76] + x201 + x52 + x6*L[65] + x97 + L[29];
#pragma omp atomic
Ls[30] += x*x204 + x13*L[77] + x20*L[79] + x28 + x48 + x49*y + x6*L[66] + x61 + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += x*x205 + x128 + x13*L[78] + x137 + x138*y + x142 + x20*L[80] + x6*L[67] + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += x13*L[79] + x164 + x167 + x180 + x181*y + x183 + x20*L[81] + x6*L[68] + x84 + L[32];
#pragma omp atomic
Ls[33] += x*x200 + x107 + x113 + x118 + x13*L[80] + x199 + x20*L[82] + x203*y + x6*L[69] + L[33];
#pragma omp atomic
Ls[34] += x*x206 + x13*L[81] + x20*L[83] + x207*y + x32 + x44 + x55 + x6*L[70] + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += x*L[56] + x40 + x41 + L[35];
#pragma omp atomic
Ls[36] += x*L[57] + x111 + x66 + L[36];
#pragma omp atomic
Ls[37] += x*L[58] + x134 + x135 + L[37];
#pragma omp atomic
Ls[38] += x163 + x186 + x92 + L[38];
#pragma omp atomic
Ls[39] += x*L[60] + x145 + x172 + L[39];
#pragma omp atomic
Ls[40] += x179 + x187 + x91 + L[40];
#pragma omp atomic
Ls[41] += x149 + x192 + x90 + L[41];
#pragma omp atomic
Ls[42] += x147 + x189 + x196 + L[42];
#pragma omp atomic
Ls[43] += x121 + x184 + x188 + L[43];
#pragma omp atomic
Ls[44] += x151 + x202 + x88 + L[44];
#pragma omp atomic
Ls[45] += x204 + x60 + x95 + L[45];
#pragma omp atomic
Ls[46] += x141 + x158 + x205 + L[46];
#pragma omp atomic
Ls[47] += x103 + x166 + x182 + L[47];
#pragma omp atomic
Ls[48] += x117 + x155 + x200 + L[48];
#pragma omp atomic
Ls[49] += x206 + x54 + x98 + L[49];
#pragma omp atomic
Ls[50] += x29 + x49 + y*L[77] + L[50];
#pragma omp atomic
Ls[51] += x129 + x138 + y*L[78] + L[51];
#pragma omp atomic
Ls[52] += x181 + x193 + x86 + L[52];
#pragma omp atomic
Ls[53] += x165 + x203 + x85 + L[53];
#pragma omp atomic
Ls[54] += x108 + x114 + x207 + L[54];
#pragma omp atomic
Ls[55] += x33 + x45 + z*L[83] + L[55];
#pragma omp atomic
Ls[56] += L[56];
#pragma omp atomic
Ls[57] += L[57];
#pragma omp atomic
Ls[58] += L[58];
#pragma omp atomic
Ls[59] += L[59];
#pragma omp atomic
Ls[60] += L[60];
#pragma omp atomic
Ls[61] += L[61];
#pragma omp atomic
Ls[62] += L[62];
#pragma omp atomic
Ls[63] += L[63];
#pragma omp atomic
Ls[64] += L[64];
#pragma omp atomic
Ls[65] += L[65];
#pragma omp atomic
Ls[66] += L[66];
#pragma omp atomic
Ls[67] += L[67];
#pragma omp atomic
Ls[68] += L[68];
#pragma omp atomic
Ls[69] += L[69];
#pragma omp atomic
Ls[70] += L[70];
#pragma omp atomic
Ls[71] += L[71];
#pragma omp atomic
Ls[72] += L[72];
#pragma omp atomic
Ls[73] += L[73];
#pragma omp atomic
Ls[74] += L[74];
#pragma omp atomic
Ls[75] += L[75];
#pragma omp atomic
Ls[76] += L[76];
#pragma omp atomic
Ls[77] += L[77];
#pragma omp atomic
Ls[78] += L[78];
#pragma omp atomic
Ls[79] += L[79];
#pragma omp atomic
Ls[80] += L[80];
#pragma omp atomic
Ls[81] += L[81];
#pragma omp atomic
Ls[82] += L[82];
#pragma omp atomic
Ls[83] += L[83];

}

void L2P_7(double x, double y, double z, double * L, double * F) {
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
x0 = x*y;
x1 = x*z;
x2 = y*z;
x3 = (x*x);
x4 = (1.0/2.0)*x3;
x5 = (x*x*x);
x6 = (1.0/6.0)*x5;
x7 = (x*x*x*x);
x8 = (1.0/24.0)*x7;
x9 = (1.0/120.0)*pow(x, 5);
x10 = (y*y);
x11 = (1.0/2.0)*x10;
x12 = (y*y*y);
x13 = (1.0/6.0)*x12;
x14 = (y*y*y*y);
x15 = (1.0/24.0)*x14;
x16 = (1.0/120.0)*pow(y, 5);
x17 = (z*z);
x18 = (1.0/2.0)*x17;
x19 = (z*z*z);
x20 = (1.0/6.0)*x19;
x21 = (z*z*z*z);
x22 = (1.0/24.0)*x21;
x23 = (1.0/120.0)*pow(z, 5);
x24 = (1.0/4.0)*x3;
x25 = x10*x24;
x26 = (1.0/12.0)*x3;
x27 = x12*x26;
x28 = (1.0/48.0)*x3;
x29 = x17*x24;
x30 = x19*x26;
x31 = (1.0/12.0)*x5;
x32 = x10*x31;
x33 = (1.0/36.0)*x5;
x34 = x17*x31;
x35 = (1.0/48.0)*x7;
x36 = x10*x17;
x37 = (1.0/4.0)*x36;
x38 = (1.0/12.0)*x10*x19;
x39 = (1.0/12.0)*x12*x17;
#pragma omp atomic
F[0] += (1.0/720.0)*pow(x, 6)*L[56] + x*x11*L[13] + x*x13*L[26] + x*x15*L[45] + x*x16*L[71] + x*x18*L[15] + x*x20*L[29] + x*x22*L[49] + x*x23*L[76] + x*x37*L[47] + x*x38*L[74] + x*x39*L[73] + x*L[1] + x0*x18*L[28] + x0*x20*L[48] + x0*x22*L[75] + x0*z*L[14] + x0*L[5] + x1*x11*L[27] + x1*x13*L[46] + x1*x15*L[72] + x1*L[6] + (1.0/48.0)*x10*x21*L[81] + x10*x35*L[59] + x11*z*L[17] + x11*L[7] + (1.0/36.0)*x12*x19*L[80] + x12*x33*L[62] + x13*z*L[31] + x13*L[16] + (1.0/48.0)*x14*x17*L[79] + x14*x28*L[66] + x15*z*L[51] + x15*L[30] + x16*z*L[78] + x16*L[50] + x17*x35*L[61] + x18*y*L[18] + x18*L[9] + x19*x33*L[65] + x2*x4*L[24] + x2*x6*L[39] + x2*x8*L[60] + x2*L[8] + x20*y*L[33] + x20*L[19] + x21*x28*L[70] + x22*y*L[54] + x22*L[34] + x23*y*L[82] + x23*L[55] + x25*z*L[42] + x25*L[23] + x27*z*L[67] + x27*L[41] + x29*y*L[43] + x29*L[25] + (1.0/8.0)*x3*x36*L[68] + x30*y*L[69] + x30*L[44] + x32*z*L[63] + x32*L[38] + x34*y*L[64] + x34*L[40] + x37*L[32] + x38*L[53] + x39*L[52] + x4*y*L[11] + x4*z*L[12] + x4*L[4] + x6*y*L[21] + x6*z*L[22] + x6*L[10] + x8*y*L[36] + x8*z*L[37] + x8*L[20] + x9*y*L[57] + x9*z*L[58] + x9*L[35] + (1.0/720.0)*pow(y, 6)*L[77] + y*L[2] + (1.0/720.0)*pow(z, 6)*L[83] + z*L[3] + L[0];

}

void M2P_7(double x, double y, double z, double * M, double * F) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*y;
x8 = x*z;
x9 = y*z;
x10 = pow(x2, -3.5);
x11 = 15.0*x10;
x12 = x7*z;
x13 = -x4;
x14 = x0*x6;
x15 = x1*x6;
x16 = 9.0*x5;
x17 = x0*x11;
x18 = -x17;
x19 = x*(x16 + x18);
x20 = x18 + x6;
x21 = x20*y;
x22 = x1*x11;
x23 = -x22;
x24 = y*(x16 + x23);
x25 = x20*z;
x26 = z*(x23 + x6);
x27 = 1.0*x;
x28 = x27*(x22 - x6);
x29 = 45.0*x10;
x30 = -x29;
x31 = pow(x2, -4.5);
x32 = x0*x31;
x33 = 105.0*x32;
x34 = x30 + x33;
x35 = x34*x7;
x36 = x34*x8;
x37 = -x11;
x38 = x9*(x33 + x37);
x39 = x1*x31;
x40 = 105.0*x39;
x41 = x30 + x40;
x42 = x41*x9;
x43 = x27*y;
x44 = x41*x43;
x45 = x37 + x40;
x46 = x27*z;
x47 = x45*x46;
x48 = 315.0*x31;
x49 = pow(x2, -5.5);
x50 = 945.0*x49;
x51 = x0*x50;
x52 = x48 - x51;
x53 = x12*x52;
x54 = x1*x50;
x55 = x27*x9;
x56 = x55*(-x48 + x54);
x57 = 90.0*x10;
x58 = (x*x*x*x);
x59 = 105.0*x31;
x60 = x58*x59;
x61 = (y*y*y*y);
x62 = x59*x61;
x63 = -225.0*x10;
x64 = x50*x58;
x65 = -x64;
x66 = x*(1050.0*x32 + x63 + x65);
x67 = x30 + 630.0*x32 + x65;
x68 = x67*y;
x69 = x50*x61;
x70 = -x69;
x71 = y*(1050.0*x39 + x63 + x70);
x72 = x67*z;
x73 = 630.0*x39;
x74 = x30 + x70 + x73;
x75 = x74*z;
x76 = x27*(x29 + x69 - x73);
x77 = 1575.0*x31;
x78 = x0*x49;
x79 = pow(x2, -6.5);
x80 = x58*x79;
x81 = 10395.0*x80;
x82 = x77 - 9450.0*x78 + x81;
x83 = x7*x82;
x84 = x8*x82;
x85 = 5670.0*x78;
x86 = x48 + x81 - x85;
x87 = x86*x9;
x88 = x1*x49;
x89 = x61*x79;
x90 = 10395.0*x89;
x91 = x77 - 9450.0*x88 + x90;
x92 = x9*x91;
x93 = x43*x91;
x94 = x48 - 5670.0*x88 + x90;
x95 = x46*x94;
x96 = 14175.0*x49;
x97 = x0*x79;
x98 = 103950.0*x97;
x99 = 135135.0*pow(x2, -7.5);
x100 = x58*x99;
x101 = x12*(-x100 - x96 + x98);
x102 = x1*x79;
x103 = x61*x99;
x104 = x55*(-103950.0*x102 + x103 + x96);
x105 = pow(x, 6);
x106 = 10395.0*x79;
x107 = x105*x106;
x108 = pow(y, 6);
x109 = x106*x108;
x110 = 11025.0*x31;
x111 = -x105*x99;
x112 = x*(x110 + x111 - 99225.0*x78 + 218295.0*x80);
x113 = x111 + x77 - 42525.0*x78 + 155925.0*x80;
x114 = x113*y;
x115 = x108*x99;
x116 = -x115;
x117 = y*(x110 + x116 - 99225.0*x88 + 218295.0*x89);
x118 = x113*z;
x119 = 42525.0*x88;
x120 = 155925.0*x89;
x121 = z*(x116 - x119 + x120 + x77);
x122 = x27*(x115 + x119 - x120 - x77);
x123 = x1*x33;
x124 = x1*x48;
x125 = -x1*x51;
x126 = x*(x124 + x125 + x34);
x127 = x0*x48;
x128 = y*(x125 + x127 + x41);
x129 = z*(x125 + x33 + x45);
x130 = -2835.0*x88;
x131 = x1*x97;
x132 = 10395.0*x131;
x133 = x130 + x132;
x134 = 945.0*x31;
x135 = -2835.0*x78;
x136 = x134 + x135;
x137 = x7*(x133 + x136);
x138 = x8*(x133 + x52);
x139 = x9*(x132 + x135 + x48 - x54);
x140 = x12*(-x0*x1*x99 + 31185.0*x102 - 8505.0*x49 + 31185.0*x97);
x141 = 105.0*x10;
x142 = -x123 - 12.0*x5;
x143 = 120.0*x10;
x144 = x0*x90;
x145 = x1*x85;
x146 = -x145;
x147 = x1*x81;
x148 = 31185.0*x89;
x149 = 62370.0*x131;
x150 = -x0*x103;
x151 = x149 + x150;
x152 = x*(x136 + x148 + x151 - 17010.0*x88);
x153 = x1*x98;
x154 = -x1*x100;
x155 = x*(-x1*x96 + x153 + x154 + x82);
x156 = y*(-x0*x96 + x150 + x153 + x91);
x157 = 31185.0*x80;
x158 = 17010.0*x78;
x159 = x130 + x149 + x154;
x160 = y*(x134 + x157 - x158 + x159);
x161 = z*(x135 + x151 + x94);
x162 = z*(x159 + x86);
x163 = 15120.0*x49;
x164 = -x107;
x165 = -x147;
x166 = 270.0*x10 + x145;
x167 = -x109;
x168 = -x144;
x169 = -x126;
x170 = -x128;
x171 = -x129;
x172 = -x137;
x173 = -x138;
x174 = -x139;
x175 = -x140;
x176 = x1*x78;
x177 = 16065.0*x49;
x178 = -x1*x158 - 360.0*x10;
x179 = 17010.0*x49;
x180 = -x112;
x181 = -x155;
x182 = -x152;
x183 = -x114;
x184 = -x160;
x185 = -x117;
x186 = -x156;
x187 = -x118;
x188 = -x162;
x189 = -x121;
x190 = -x161;
#pragma omp atomic
F[0] += -x*x4*M[0] + x101*M[87] - x104*M[105] - x11*x12*M[13] + x112*M[83] + x114*M[84] + x117*M[111] + x118*M[85] + x121*M[112] - x122*M[104] + x126*M[37] + x128*M[40] + x129*M[41] + x137*M[61] + x138*M[62] + x139*M[66] + x140*M[94] + x152*M[93] + x155*M[86] + x156*M[98] + x160*M[89] + x161*M[99] + x162*M[90] + x19*M[9] + x21*M[10] + x24*M[15] + x25*M[11] + x26*M[16] - x28*M[12] + x35*M[20] + x36*M[21] + x38*M[23] - x4*y*M[1] - x4*z*M[2] + x42*M[30] + x44*M[25] + x47*M[26] + x53*M[38] - x56*M[45] + x6*x7*M[4] + x6*x8*M[5] + x6*x9*M[7] + x66*M[34] + x68*M[35] + x71*M[49] + x72*M[36] + x75*M[50] - x76*M[44] + x83*M[56] + x84*M[57] + x87*M[59] + x92*M[77] + x93*M[70] + x95*M[71] + (-x101 + x175)*M[96] + (x104 + x175)*M[107] + (x122 + x182)*M[106] + (x13 + x14)*M[3] + (x13 + x15)*M[6] + (x169 - x66)*M[39] + (x169 + x76)*M[46] + (x170 - x68)*M[42] + (x170 - x71)*M[51] + (x171 - x72)*M[43] + (x171 - x75)*M[52] + (x172 - x83)*M[63] + (x172 - x93)*M[72] + (x173 - x84)*M[64] + (x173 - x95)*M[73] + (x174 - x87)*M[68] + (x174 - x92)*M[79] + (x180 + x181)*M[88] + (x181 + x182)*M[95] + (x183 + x184)*M[91] + (x184 + x186)*M[100] + (x185 + x186)*M[113] + (x187 + x188)*M[92] + (x188 + x190)*M[101] + (x189 + x190)*M[114] + (-x19 + x28)*M[14] + (-x21 - x24)*M[17] + (-x25 - x26)*M[18] + (-x35 - x44)*M[27] + (-x36 - x47)*M[28] + (-x38 - x42)*M[32] + (-x53 + x56)*M[47] + (x101 - x104 + 2*x140)*M[109] + (x112 + x152 + 2*x155)*M[97] + (x114 + x156 + 2*x160)*M[102] + (x117 + 2*x156 + x160)*M[115] + (x118 + x161 + 2*x162)*M[103] + (x121 + 2*x161 + x162)*M[116] + (-x122 + 2*x152 + x155)*M[108] + (x123 + x20 + x23)*M[22] + (2*x126 + x66 - x76)*M[48] + (2*x128 + x68 + x71)*M[53] + (2*x129 + x72 + x75)*M[54] + (2*x137 + x83 + x93)*M[74] + (2*x138 + x84 + x95)*M[75] + (2*x139 + x87 + x92)*M[81] + (-x14 - x15 + 2.0*x3)*M[8] + (-x0*x57 + x16 + x60)*M[19] + (-x1*x57 + x16 + x62)*M[29] + (x107 + 4725.0*x32 - x58*x96 + x63)*M[55] + (x109 + 4725.0*x39 - x61*x96 + x63)*M[76] + (x122 - 3*x152 - 3*x155 + x180)*M[110] + (x124 + x146 + x147 + x67)*M[58] + (x127 + x144 + x146 + x74)*M[65] + (-3*x156 - 3*x160 + x183 + x185)*M[117] + (-3*x161 - 3*x162 + x187 + x189)*M[118] + (x0*x141 + x142 + x22 - x60)*M[24] + (x1*x141 + x142 + x17 - x62)*M[31] + (-x124 + x163*x58 + x164 + x165 + x166 - 5355.0*x32)*M[60] + (-x127 + x163*x61 + x166 + x167 + x168 - 5355.0*x39)*M[78] + (-x0*x143 - x1*x143 + 210.0*x1*x32 + 24.0*x5 + x60 + x62)*M[33] + (-x0*x134 - x1*x134 + x165 + x168 + 11340.0*x176 + x57 + x64 + x69)*M[67] + (20790.0*x0*x89 + x109 + x147 - x177*x61 + x178 + 1260.0*x32 + 6300.0*x39 + x65)*M[80] + (20790.0*x1*x80 + x107 + x144 - x177*x58 + x178 + 6300.0*x32 + 1260.0*x39 + x70)*M[69] + (-x0*x148 - x1*x157 + 720.0*x10 + x164 + x167 + 34020.0*x176 + x179*x58 + x179*x61 - 7560.0*x32 - 7560.0*x39)*M[82];

}

void P2M_8(double x, double y, double z, double q, double * M) {
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
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = x1*z;
x9 = (z*z);
x10 = (x*x*x);
x11 = (1.0/6.0)*q;
x12 = (1.0/2.0)*x3;
x13 = (1.0/2.0)*x0;
x14 = (y*y*y);
x15 = (1.0/2.0)*x7;
x16 = (1.0/2.0)*x9;
x17 = (z*z*z);
x18 = (x*x*x*x);
x19 = (1.0/24.0)*q;
x20 = (1.0/6.0)*x10;
x21 = q*x7;
x22 = (1.0/4.0)*x3;
x23 = q*x9;
x24 = (1.0/6.0)*x0;
x25 = (y*y*y*y);
x26 = (1.0/6.0)*x14;
x27 = (1.0/4.0)*x9;
x28 = (1.0/6.0)*x17;
x29 = (z*z*z*z);
x30 = pow(x, 5);
x31 = (1.0/120.0)*q;
x32 = (1.0/24.0)*x18;
x33 = (1.0/12.0)*x10;
x34 = (1.0/12.0)*x14;
x35 = q*x3;
x36 = x2*x7;
x37 = x1*x9;
x38 = (1.0/12.0)*x17;
x39 = (1.0/24.0)*x0;
x40 = x0*x7;
x41 = pow(y, 5);
x42 = (1.0/24.0)*x25;
x43 = (1.0/24.0)*x29;
x44 = pow(z, 5);
x45 = pow(x, 6);
x46 = (1.0/720.0)*q;
x47 = (1.0/120.0)*x30;
x48 = (1.0/48.0)*x18;
x49 = q*x14;
x50 = (1.0/36.0)*x10;
x51 = q*x17;
x52 = (1.0/48.0)*x35;
x53 = x2*x3;
x54 = x3*x9;
x55 = x1*x3;
x56 = (1.0/120.0)*x0;
x57 = x0*x9;
x58 = pow(y, 6);
x59 = (1.0/120.0)*x41;
x60 = (1.0/48.0)*x25;
x61 = (1.0/36.0)*x17;
x62 = (1.0/48.0)*x29;
x63 = (1.0/120.0)*x44;
x64 = pow(z, 6);
x65 = pow(x, 7);
x66 = (1.0/5040.0)*q;
x67 = (1.0/720.0)*x45;
x68 = (1.0/240.0)*x30;
x69 = (1.0/144.0)*x18;
x70 = (1.0/144.0)*x25;
x71 = q*x10;
x72 = x14*x2;
x73 = x19*x7;
x74 = x1*x17;
x75 = (1.0/144.0)*x29;
x76 = (1.0/240.0)*x35;
x77 = (1.0/720.0)*x0;
x78 = x0*x14;
x79 = pow(y, 7);
x80 = (1.0/720.0)*x58;
x81 = (1.0/240.0)*x41;
x82 = (1.0/240.0)*x44;
x83 = (1.0/720.0)*x64;
x84 = pow(z, 7);
x85 = (1.0/40320.0)*q;
x86 = (1.0/5040.0)*x65;
x87 = (1.0/1440.0)*x45;
x88 = x30*x46;
x89 = (1.0/576.0)*q*x18;
x90 = (1.0/96.0)*x21;
x91 = x10*x46;
x92 = (1.0/72.0)*x10;
x93 = (1.0/1440.0)*x35;
x94 = (1.0/5040.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x6;
M[6] += x4*x7;
M[7] += x8;
M[8] += x4*x9;
M[9] += -x10*x11;
M[10] += -x1*x12;
M[11] += -x12*x2;
M[12] += -x13*x7;
M[13] += -x5*z;
M[14] += -x13*x9;
M[15] += -x11*x14;
M[16] += -x15*x2;
M[17] += -x1*x16;
M[18] += -x11*x17;
M[19] += x18*x19;
M[20] += x1*x20;
M[21] += x2*x20;
M[22] += x21*x22;
M[23] += x12*x8;
M[24] += x22*x23;
M[25] += x14*x24;
M[26] += x15*x6;
M[27] += x16*x5;
M[28] += x17*x24;
M[29] += x19*x25;
M[30] += x2*x26;
M[31] += x21*x27;
M[32] += x1*x28;
M[33] += x19*x29;
M[34] += -x30*x31;
M[35] += -x1*x32;
M[36] += -x2*x32;
M[37] += -x21*x33;
M[38] += -x20*x8;
M[39] += -x23*x33;
M[40] += -x34*x35;
M[41] += -x22*x36;
M[42] += -x22*x37;
M[43] += -x35*x38;
M[44] += -x25*x39;
M[45] += -x26*x6;
M[46] += -x27*x40;
M[47] += -x28*x5;
M[48] += -x29*x39;
M[49] += -x31*x41;
M[50] += -x2*x42;
M[51] += -x23*x34;
M[52] += -x21*x38;
M[53] += -x1*x43;
M[54] += -x31*x44;
M[55] += x45*x46;
M[56] += x1*x47;
M[57] += x2*x47;
M[58] += x21*x48;
M[59] += x32*x8;
M[60] += x23*x48;
M[61] += x49*x50;
M[62] += x33*x36;
M[63] += x33*x37;
M[64] += x50*x51;
M[65] += x25*x52;
M[66] += x34*x53;
M[67] += (1.0/8.0)*x21*x54;
M[68] += x38*x55;
M[69] += x29*x52;
M[70] += x41*x56;
M[71] += x42*x6;
M[72] += x34*x57;
M[73] += x38*x40;
M[74] += x43*x5;
M[75] += x44*x56;
M[76] += x46*x58;
M[77] += x2*x59;
M[78] += x23*x60;
M[79] += x49*x61;
M[80] += x21*x62;
M[81] += x1*x63;
M[82] += x46*x64;
M[83] += -x65*x66;
M[84] += -x1*x67;
M[85] += -x2*x67;
M[86] += -x21*x68;
M[87] += -x47*x8;
M[88] += -x23*x68;
M[89] += -x49*x69;
M[90] += -x36*x48;
M[91] += -x37*x48;
M[92] += -x51*x69;
M[93] += -x70*x71;
M[94] += -x50*x72;
M[95] += -x10*x73*x9;
M[96] += -x50*x74;
M[97] += -x71*x75;
M[98] += -x41*x76;
M[99] += -x53*x60;
M[100] += -x14*x19*x54;
M[101] += -x17*x3*x73;
M[102] += -x55*x62;
M[103] += -x44*x76;
M[104] += -x58*x77;
M[105] += -x59*x6;
M[106] += -x57*x60;
M[107] += -x61*x78;
M[108] += -x40*x62;
M[109] += -x5*x63;
M[110] += -x64*x77;
M[111] += -x66*x79;
M[112] += -x2*x80;
M[113] += -x23*x81;
M[114] += -x51*x70;
M[115] += -x49*x75;
M[116] += -x21*x82;
M[117] += -x1*x83;
M[118] += -x66*x84;
M[119] += pow(x, 8)*x85;
M[120] += x1*x86;
M[121] += x2*x86;
M[122] += x21*x87;
M[123] += x67*x8;
M[124] += x23*x87;
M[125] += x14*x88;
M[126] += x36*x68;
M[127] += x37*x68;
M[128] += x17*x88;
M[129] += x25*x89;
M[130] += x69*x72;
M[131] += x18*x9*x90;
M[132] += x69*x74;
M[133] += x29*x89;
M[134] += x41*x91;
M[135] += x10*x2*x70;
M[136] += x14*x23*x92;
M[137] += x17*x21*x92;
M[138] += x1*x10*x75;
M[139] += x44*x91;
M[140] += x58*x93;
M[141] += x53*x81;
M[142] += (1.0/96.0)*x23*x25*x3;
M[143] += (1.0/72.0)*x14*x17*x35;
M[144] += x29*x3*x90;
M[145] += x55*x82;
M[146] += x64*x93;
M[147] += x79*x94;
M[148] += x6*x80;
M[149] += x57*x81;
M[150] += x0*x17*x70;
M[151] += x75*x78;
M[152] += x40*x82;
M[153] += x5*x83;
M[154] += x84*x94;
M[155] += x85*pow(y, 8);
M[156] += (1.0/5040.0)*x2*x79;
M[157] += (1.0/1440.0)*x23*x58;
M[158] += x17*x41*x46;
M[159] += (1.0/576.0)*q*x25*x29;
M[160] += x14*x44*x46;
M[161] += (1.0/1440.0)*x21*x64;
M[162] += (1.0/5040.0)*x1*x84;
M[163] += x85*pow(z, 8);
}
void M2M_8(double x, double y, double z, double * M, double * Ms) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = z*M[0];
x5 = y*M[1];
x6 = y*M[2];
x7 = z*M[1];
x8 = z*M[2];
x9 = x*M[3];
x10 = (x*x);
x11 = (1.0/2.0)*x10;
x12 = x*M[4];
x13 = y*M[3];
x14 = x0*y;
x15 = x*M[5];
x16 = z*M[3];
x17 = x0*z;
x18 = x*M[6];
x19 = y*M[4];
x20 = x1*y;
x21 = (y*y);
x22 = (1.0/2.0)*M[0];
x23 = x*M[7];
x24 = y*M[5];
x25 = z*M[4];
x26 = x3*y;
x27 = x1*z;
x28 = x2*z;
x29 = x*M[8];
x30 = z*M[5];
x31 = x3*z;
x32 = (z*z);
x33 = y*M[6];
x34 = (1.0/2.0)*x21;
x35 = y*M[7];
x36 = z*M[6];
x37 = x5*z;
x38 = y*M[8];
x39 = z*M[7];
x40 = x6*z;
x41 = (1.0/2.0)*x32;
x42 = z*M[8];
x43 = x*M[9];
x44 = (x*x*x);
x45 = (1.0/6.0)*x44;
x46 = x*M[10];
x47 = y*M[9];
x48 = x9*y;
x49 = x*M[11];
x50 = z*M[9];
x51 = x9*z;
x52 = x*M[12];
x53 = y*M[10];
x54 = x12*y;
x55 = x*M[13];
x56 = y*M[11];
x57 = z*M[10];
x58 = x15*y;
x59 = x12*z;
x60 = x13*z;
x61 = x*M[14];
x62 = z*M[11];
x63 = x15*z;
x64 = x*M[15];
x65 = y*M[12];
x66 = x18*y;
x67 = (y*y*y);
x68 = (1.0/6.0)*M[0];
x69 = x*M[16];
x70 = y*M[13];
x71 = z*M[12];
x72 = x23*y;
x73 = x18*z;
x74 = x19*z;
x75 = x*M[17];
x76 = y*M[14];
x77 = z*M[13];
x78 = x29*y;
x79 = x23*z;
x80 = x24*z;
x81 = x*M[18];
x82 = z*M[14];
x83 = x29*z;
x84 = (z*z*z);
x85 = y*M[15];
x86 = (1.0/6.0)*x67;
x87 = y*M[16];
x88 = z*M[15];
x89 = x33*z;
x90 = y*M[17];
x91 = z*M[16];
x92 = x35*z;
x93 = y*M[18];
x94 = z*M[17];
x95 = x38*z;
x96 = (1.0/6.0)*x84;
x97 = z*M[18];
x98 = x*M[19];
x99 = (x*x*x*x);
x100 = (1.0/24.0)*x99;
x101 = x*M[20];
x102 = y*M[19];
x103 = x43*y;
x104 = x*M[21];
x105 = z*M[19];
x106 = x43*z;
x107 = x*M[22];
x108 = y*M[20];
x109 = x46*y;
x110 = (1.0/4.0)*x10;
x111 = x21*M[0];
x112 = x*M[23];
x113 = y*M[21];
x114 = z*M[20];
x115 = x49*y;
x116 = x46*z;
x117 = x47*z;
x118 = x*M[24];
x119 = z*M[21];
x120 = x49*z;
x121 = x110*x32;
x122 = x*M[25];
x123 = y*M[22];
x124 = x52*y;
x125 = x110*x21;
x126 = x*M[26];
x127 = y*M[23];
x128 = z*M[22];
x129 = x55*y;
x130 = x52*z;
x131 = x53*z;
x132 = x*M[27];
x133 = y*M[24];
x134 = z*M[23];
x135 = x61*y;
x136 = x55*z;
x137 = x56*z;
x138 = x*M[28];
x139 = z*M[24];
x140 = x61*z;
x141 = x*M[29];
x142 = y*M[25];
x143 = x64*y;
x144 = (y*y*y*y);
x145 = (1.0/24.0)*M[0];
x146 = x*M[30];
x147 = y*M[26];
x148 = z*M[25];
x149 = x69*y;
x150 = x64*z;
x151 = x65*z;
x152 = x*M[31];
x153 = y*M[27];
x154 = z*M[26];
x155 = x75*y;
x156 = x69*z;
x157 = x70*z;
x158 = (1.0/4.0)*x32;
x159 = x*M[32];
x160 = y*M[28];
x161 = z*M[27];
x162 = x81*y;
x163 = x75*z;
x164 = x76*z;
x165 = x*M[33];
x166 = z*M[28];
x167 = x81*z;
x168 = (z*z*z*z);
x169 = y*M[29];
x170 = (1.0/24.0)*x144;
x171 = y*M[30];
x172 = z*M[29];
x173 = x85*z;
x174 = y*M[31];
x175 = z*M[30];
x176 = x87*z;
x177 = x158*x21;
x178 = y*M[32];
x179 = z*M[31];
x180 = x90*z;
x181 = y*M[33];
x182 = z*M[32];
x183 = x93*z;
x184 = (1.0/24.0)*x168;
x185 = z*M[33];
x186 = x*M[34];
x187 = pow(x, 5);
x188 = (1.0/120.0)*x187;
x189 = x*M[35];
x190 = y*M[34];
x191 = x98*y;
x192 = x*M[36];
x193 = z*M[34];
x194 = x98*z;
x195 = x*M[37];
x196 = y*M[35];
x197 = x101*y;
x198 = (1.0/12.0)*x44;
x199 = x*M[38];
x200 = y*M[36];
x201 = z*M[35];
x202 = x104*y;
x203 = x101*z;
x204 = x102*z;
x205 = x*M[39];
x206 = z*M[36];
x207 = x104*z;
x208 = x198*x32;
x209 = x*M[40];
x210 = y*M[37];
x211 = x107*y;
x212 = (1.0/12.0)*x10;
x213 = x67*M[0];
x214 = x198*x21;
x215 = x*M[41];
x216 = y*M[38];
x217 = z*M[37];
x218 = x112*y;
x219 = x107*z;
x220 = x108*z;
x221 = x*M[42];
x222 = y*M[39];
x223 = z*M[38];
x224 = x118*y;
x225 = x112*z;
x226 = x113*z;
x227 = x*M[43];
x228 = z*M[39];
x229 = x118*z;
x230 = x212*x84;
x231 = x*M[44];
x232 = y*M[40];
x233 = x122*y;
x234 = x212*x67;
x235 = x*M[45];
x236 = y*M[41];
x237 = z*M[40];
x238 = x126*y;
x239 = x122*z;
x240 = x123*z;
x241 = x*M[46];
x242 = y*M[42];
x243 = z*M[41];
x244 = x132*y;
x245 = x126*z;
x246 = x127*z;
x247 = x*M[47];
x248 = y*M[43];
x249 = z*M[42];
x250 = x138*y;
x251 = x132*z;
x252 = x133*z;
x253 = x*M[48];
x254 = z*M[43];
x255 = x138*z;
x256 = x*M[49];
x257 = y*M[44];
x258 = x141*y;
x259 = pow(y, 5);
x260 = (1.0/120.0)*M[0];
x261 = x*M[50];
x262 = y*M[45];
x263 = z*M[44];
x264 = x146*y;
x265 = x141*z;
x266 = x142*z;
x267 = x*M[51];
x268 = y*M[46];
x269 = z*M[45];
x270 = x152*y;
x271 = x146*z;
x272 = x147*z;
x273 = (1.0/12.0)*x32;
x274 = x*M[52];
x275 = y*M[47];
x276 = z*M[46];
x277 = x159*y;
x278 = x152*z;
x279 = x153*z;
x280 = (1.0/12.0)*x84;
x281 = x*M[53];
x282 = y*M[48];
x283 = z*M[47];
x284 = x165*y;
x285 = x159*z;
x286 = x160*z;
x287 = x*M[54];
x288 = z*M[48];
x289 = x165*z;
x290 = pow(z, 5);
x291 = y*M[49];
x292 = (1.0/120.0)*x259;
x293 = y*M[50];
x294 = z*M[49];
x295 = x169*z;
x296 = y*M[51];
x297 = z*M[50];
x298 = x171*z;
x299 = x273*x67;
x300 = y*M[52];
x301 = z*M[51];
x302 = x174*z;
x303 = x21*x280;
x304 = y*M[53];
x305 = z*M[52];
x306 = x178*z;
x307 = y*M[54];
x308 = z*M[53];
x309 = x181*z;
x310 = (1.0/120.0)*x290;
x311 = z*M[54];
x312 = x*M[55];
x313 = (1.0/720.0)*pow(x, 6);
x314 = x*M[56];
x315 = y*M[55];
x316 = x186*y;
x317 = x*M[57];
x318 = x*M[58];
x319 = y*M[56];
x320 = x189*y;
x321 = (1.0/48.0)*x99;
x322 = x*M[59];
x323 = y*M[57];
x324 = x192*y;
x325 = x*M[60];
x326 = x32*x321;
x327 = x*M[61];
x328 = y*M[58];
x329 = x195*y;
x330 = (1.0/36.0)*x44;
x331 = x21*x321;
x332 = x*M[62];
x333 = y*M[59];
x334 = x199*y;
x335 = x*M[63];
x336 = y*M[60];
x337 = x205*y;
x338 = x*M[64];
x339 = x330*x84;
x340 = x*M[65];
x341 = y*M[61];
x342 = x209*y;
x343 = (1.0/48.0)*x10;
x344 = x144*M[0];
x345 = x330*x67;
x346 = x*M[66];
x347 = y*M[62];
x348 = x215*y;
x349 = x*M[67];
x350 = y*M[63];
x351 = x221*y;
x352 = x10*x32;
x353 = (1.0/8.0)*x352;
x354 = x*M[68];
x355 = y*M[64];
x356 = x227*y;
x357 = x*M[69];
x358 = x168*x343;
x359 = x*M[70];
x360 = y*M[65];
x361 = x231*y;
x362 = x144*x343;
x363 = x*M[71];
x364 = y*M[66];
x365 = x235*y;
x366 = x*M[72];
x367 = y*M[67];
x368 = x241*y;
x369 = x21*x353;
x370 = x*M[73];
x371 = y*M[68];
x372 = x247*y;
x373 = x*M[74];
x374 = y*M[69];
x375 = x253*y;
x376 = x*M[75];
x377 = x*M[76];
x378 = y*M[70];
x379 = x256*y;
x380 = pow(y, 6);
x381 = (1.0/720.0)*M[0];
x382 = x*M[77];
x383 = y*M[71];
x384 = x261*y;
x385 = x*M[78];
x386 = y*M[72];
x387 = x267*y;
x388 = (1.0/48.0)*x32;
x389 = x*M[79];
x390 = y*M[73];
x391 = x274*y;
x392 = (1.0/36.0)*x84;
x393 = x*M[80];
x394 = y*M[74];
x395 = x281*y;
x396 = (1.0/48.0)*x168;
x397 = x*M[81];
x398 = y*M[75];
x399 = x287*y;
x400 = x*M[82];
x401 = pow(z, 6);
x402 = y*M[76];
x403 = (1.0/720.0)*x380;
x404 = y*M[77];
x405 = y*M[78];
x406 = x144*x388;
x407 = y*M[79];
x408 = x392*x67;
x409 = y*M[80];
x410 = x21*x396;
x411 = y*M[81];
x412 = y*M[82];
x413 = (1.0/720.0)*x401;
x414 = (1.0/5040.0)*pow(x, 7);
x415 = (1.0/240.0)*x187;
x416 = x32*x415;
x417 = (1.0/144.0)*x99;
x418 = x21*x415;
x419 = x417*x84;
x420 = (1.0/144.0)*x44;
x421 = x417*x67;
x422 = x145*x21;
x423 = x32*x44;
x424 = x168*x420;
x425 = (1.0/240.0)*x10;
x426 = x425*M[0];
x427 = x144*x420;
x428 = x352*x67;
x429 = (1.0/24.0)*M[1];
x430 = x21*x423;
x431 = x10*x84;
x432 = (1.0/24.0)*M[2];
x433 = x259*x425;
x434 = x21*x431;
x435 = x290*x425;
x436 = pow(y, 7);
x437 = (1.0/5040.0)*M[0];
x438 = (1.0/240.0)*x259*x32;
x439 = (1.0/144.0)*x84;
x440 = (1.0/144.0)*x168;
x441 = (1.0/240.0)*x290;
x442 = pow(z, 7);
x443 = (1.0/5040.0)*x436;
x444 = x144*x439;
x445 = x440*x67;
x446 = x21*x441;
x447 = (1.0/5040.0)*x442;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + x4 + M[5];
#pragma omp atomic
Ms[6] += x5 + M[6];
#pragma omp atomic
Ms[7] += x6 + x7 + M[7];
#pragma omp atomic
Ms[8] += x8 + M[8];
#pragma omp atomic
Ms[9] += x11*M[0] + x9 + M[9];
#pragma omp atomic
Ms[10] += x11*M[1] + x12 + x13 + x14 + M[10];
#pragma omp atomic
Ms[11] += x11*M[2] + x15 + x16 + x17 + M[11];
#pragma omp atomic
Ms[12] += x18 + x19 + x20 + x21*x22 + M[12];
#pragma omp atomic
Ms[13] += x23 + x24 + x25 + x26 + x27 + x28 + M[13];
#pragma omp atomic
Ms[14] += x22*x32 + x29 + x30 + x31 + M[14];
#pragma omp atomic
Ms[15] += x33 + x34*M[1] + M[15];
#pragma omp atomic
Ms[16] += x34*M[2] + x35 + x36 + x37 + M[16];
#pragma omp atomic
Ms[17] += x38 + x39 + x40 + x41*M[1] + M[17];
#pragma omp atomic
Ms[18] += x41*M[2] + x42 + M[18];
#pragma omp atomic
Ms[19] += x11*M[3] + x43 + x45*M[0] + M[19];
#pragma omp atomic
Ms[20] += x11*x2 + x11*M[4] + x45*M[1] + x46 + x47 + x48 + M[20];
#pragma omp atomic
Ms[21] += x11*x4 + x11*M[5] + x45*M[2] + x49 + x50 + x51 + M[21];
#pragma omp atomic
Ms[22] += x0*x34 + x11*x5 + x11*M[6] + x34*M[3] + x52 + x53 + x54 + M[22];
#pragma omp atomic
Ms[23] += x11*x6 + x11*x7 + x11*M[7] + x14*z + x55 + x56 + x57 + x58 + x59 + x60 + M[23];
#pragma omp atomic
Ms[24] += x0*x41 + x11*x8 + x11*M[8] + x41*M[3] + x61 + x62 + x63 + M[24];
#pragma omp atomic
Ms[25] += x1*x34 + x34*M[4] + x64 + x65 + x66 + x67*x68 + M[25];
#pragma omp atomic
Ms[26] += x20*z + x3*x34 + x34*x4 + x34*M[5] + x69 + x70 + x71 + x72 + x73 + x74 + M[26];
#pragma omp atomic
Ms[27] += x1*x41 + x2*x41 + x26*z + x41*M[4] + x75 + x76 + x77 + x78 + x79 + x80 + M[27];
#pragma omp atomic
Ms[28] += x3*x41 + x41*M[5] + x68*x84 + x81 + x82 + x83 + M[28];
#pragma omp atomic
Ms[29] += x34*M[6] + x85 + x86*M[1] + M[29];
#pragma omp atomic
Ms[30] += x34*x7 + x34*M[7] + x86*M[2] + x87 + x88 + x89 + M[30];
#pragma omp atomic
Ms[31] += x34*x8 + x34*M[8] + x41*x5 + x41*M[6] + x90 + x91 + x92 + M[31];
#pragma omp atomic
Ms[32] += x41*x6 + x41*M[7] + x93 + x94 + x95 + x96*M[1] + M[32];
#pragma omp atomic
Ms[33] += x41*M[8] + x96*M[2] + x97 + M[33];
#pragma omp atomic
Ms[34] += x100*M[0] + x11*M[9] + x45*M[3] + x98 + M[34];
#pragma omp atomic
Ms[35] += x100*M[1] + x101 + x102 + x103 + x11*x13 + x11*M[10] + x2*x45 + x45*M[4] + M[35];
#pragma omp atomic
Ms[36] += x100*M[2] + x104 + x105 + x106 + x11*x16 + x11*M[11] + x4*x45 + x45*M[5] + M[36];
#pragma omp atomic
Ms[37] += x107 + x108 + x109 + x11*x19 + x11*M[12] + x110*x111 + x34*x9 + x34*M[9] + x45*x5 + x45*M[6] + M[37];
#pragma omp atomic
Ms[38] += x11*x24 + x11*x25 + x11*x28 + x11*M[13] + x112 + x113 + x114 + x115 + x116 + x117 + x45*x6 + x45*x7 + x45*M[7] + x48*z + M[38];
#pragma omp atomic
Ms[39] += x11*x30 + x11*M[14] + x118 + x119 + x120 + x121*M[0] + x41*x9 + x41*M[9] + x45*x8 + x45*M[8] + M[39];
#pragma omp atomic
Ms[40] += x0*x86 + x11*x33 + x11*M[15] + x12*x34 + x122 + x123 + x124 + x125*M[1] + x34*M[10] + x86*M[3] + M[40];
#pragma omp atomic
Ms[41] += x11*x35 + x11*x36 + x11*x37 + x11*M[16] + x125*M[2] + x126 + x127 + x128 + x129 + x130 + x131 + x15*x34 + x16*x34 + x17*x34 + x34*M[11] + x54*z + M[41];
#pragma omp atomic
Ms[42] += x11*x38 + x11*x39 + x11*x40 + x11*M[17] + x12*x41 + x121*M[1] + x13*x41 + x132 + x133 + x134 + x135 + x136 + x137 + x14*x41 + x41*M[10] + x58*z + M[42];
#pragma omp atomic
Ms[43] += x0*x96 + x11*x42 + x11*M[18] + x121*M[2] + x138 + x139 + x140 + x15*x41 + x41*M[11] + x96*M[3] + M[43];
#pragma omp atomic
Ms[44] += x1*x86 + x141 + x142 + x143 + x144*x145 + x18*x34 + x34*M[12] + x86*M[4] + M[44];
#pragma omp atomic
Ms[45] += x146 + x147 + x148 + x149 + x150 + x151 + x23*x34 + x25*x34 + x27*x34 + x3*x86 + x34*M[13] + x4*x86 + x66*z + x86*M[5] + M[45];
#pragma omp atomic
Ms[46] += x111*x158 + x152 + x153 + x154 + x155 + x156 + x157 + x18*x41 + x19*x41 + x20*x41 + x29*x34 + x30*x34 + x31*x34 + x34*M[14] + x41*M[12] + x72*z + M[46];
#pragma omp atomic
Ms[47] += x1*x96 + x159 + x160 + x161 + x162 + x163 + x164 + x2*x96 + x23*x41 + x24*x41 + x26*x41 + x41*M[13] + x78*z + x96*M[4] + M[47];
#pragma omp atomic
Ms[48] += x145*x168 + x165 + x166 + x167 + x29*x41 + x3*x96 + x41*M[14] + x96*M[5] + M[48];
#pragma omp atomic
Ms[49] += x169 + x170*M[1] + x34*M[15] + x86*M[6] + M[49];
#pragma omp atomic
Ms[50] += x170*M[2] + x171 + x172 + x173 + x34*x36 + x34*M[16] + x7*x86 + x86*M[7] + M[50];
#pragma omp atomic
Ms[51] += x174 + x175 + x176 + x177*M[1] + x33*x41 + x34*x39 + x34*M[17] + x41*M[15] + x8*x86 + x86*M[8] + M[51];
#pragma omp atomic
Ms[52] += x177*M[2] + x178 + x179 + x180 + x34*x42 + x34*M[18] + x35*x41 + x41*M[16] + x5*x96 + x96*M[6] + M[52];
#pragma omp atomic
Ms[53] += x181 + x182 + x183 + x184*M[1] + x38*x41 + x41*M[17] + x6*x96 + x96*M[7] + M[53];
#pragma omp atomic
Ms[54] += x184*M[2] + x185 + x41*M[18] + x96*M[8] + M[54];
#pragma omp atomic
Ms[55] += x100*M[3] + x11*M[19] + x186 + x188*M[0] + x45*M[9] + M[55];
#pragma omp atomic
Ms[56] += x100*x2 + x100*M[4] + x11*x47 + x11*M[20] + x13*x45 + x188*M[1] + x189 + x190 + x191 + x45*M[10] + M[56];
#pragma omp atomic
Ms[57] += x100*x4 + x100*M[5] + x11*x50 + x11*M[21] + x16*x45 + x188*M[2] + x192 + x193 + x194 + x45*M[11] + M[57];
#pragma omp atomic
Ms[58] += x100*x5 + x100*M[6] + x11*x53 + x11*M[22] + x111*x198 + x125*M[3] + x19*x45 + x195 + x196 + x197 + x34*x43 + x34*M[19] + x45*M[12] + M[58];
#pragma omp atomic
Ms[59] += x100*x6 + x100*x7 + x100*M[7] + x103*z + x11*x56 + x11*x57 + x11*x60 + x11*M[23] + x199 + x200 + x201 + x202 + x203 + x204 + x24*x45 + x25*x45 + x28*x45 + x45*M[13] + M[59];
#pragma omp atomic
Ms[60] += x100*x8 + x100*M[8] + x11*x62 + x11*M[24] + x121*M[3] + x205 + x206 + x207 + x208*M[0] + x30*x45 + x41*x43 + x41*M[19] + x45*M[14] + M[60];
#pragma omp atomic
Ms[61] += x11*x65 + x11*M[25] + x125*M[4] + x209 + x210 + x211 + x212*x213 + x214*M[1] + x33*x45 + x34*x46 + x34*M[20] + x45*M[15] + x86*x9 + x86*M[9] + M[61];
#pragma omp atomic
Ms[62] += x109*z + x11*x70 + x11*x71 + x11*x74 + x11*M[26] + x125*x4 + x125*M[5] + x214*M[2] + x215 + x216 + x217 + x218 + x219 + x220 + x34*x49 + x34*x50 + x34*x51 + x34*M[21] + x35*x45 + x36*x45 + x37*x45 + x45*M[16] + M[62];
#pragma omp atomic
Ms[63] += x11*x76 + x11*x77 + x11*x80 + x11*M[27] + x115*z + x121*x2 + x121*M[4] + x208*M[1] + x221 + x222 + x223 + x224 + x225 + x226 + x38*x45 + x39*x45 + x40*x45 + x41*x46 + x41*x47 + x41*x48 + x41*M[20] + x45*M[17] + M[63];
#pragma omp atomic
Ms[64] += x11*x82 + x11*M[28] + x121*M[5] + x208*M[2] + x227 + x228 + x229 + x230*M[0] + x41*x49 + x41*M[21] + x42*x45 + x45*M[18] + x9*x96 + x96*M[9] + M[64];
#pragma omp atomic
Ms[65] += x0*x170 + x11*x85 + x11*M[29] + x12*x86 + x125*M[6] + x170*M[3] + x231 + x232 + x233 + x234*M[1] + x34*x52 + x34*M[22] + x86*M[10] + M[65];
#pragma omp atomic
Ms[66] += x11*x87 + x11*x88 + x11*x89 + x11*M[30] + x124*z + x125*x7 + x125*M[7] + x15*x86 + x16*x86 + x17*x86 + x234*M[2] + x235 + x236 + x237 + x238 + x239 + x240 + x34*x55 + x34*x57 + x34*x59 + x34*M[23] + x86*M[11] + M[66];
#pragma omp atomic
Ms[67] += x0*x177 + x11*x90 + x11*x91 + x11*x92 + x11*M[31] + x121*x5 + x121*M[6] + x125*x8 + x125*M[8] + x129*z + x177*M[3] + x241 + x242 + x243 + x244 + x245 + x246 + x34*x61 + x34*x62 + x34*x63 + x34*M[24] + x41*x52 + x41*x53 + x41*x54 + x41*M[22] + M[67];
#pragma omp atomic
Ms[68] += x11*x93 + x11*x94 + x11*x95 + x11*M[32] + x12*x96 + x121*x6 + x121*M[7] + x13*x96 + x135*z + x14*x96 + x230*M[1] + x247 + x248 + x249 + x250 + x251 + x252 + x41*x55 + x41*x56 + x41*x58 + x41*M[23] + x96*M[10] + M[68];
#pragma omp atomic
Ms[69] += x0*x184 + x11*x97 + x11*M[33] + x121*M[8] + x15*x96 + x184*M[3] + x230*M[2] + x253 + x254 + x255 + x41*x61 + x41*M[24] + x96*M[11] + M[69];
#pragma omp atomic
Ms[70] += x1*x170 + x170*M[4] + x18*x86 + x256 + x257 + x258 + x259*x260 + x34*x64 + x34*M[25] + x86*M[12] + M[70];
#pragma omp atomic
Ms[71] += x143*z + x170*x3 + x170*x4 + x170*M[5] + x23*x86 + x25*x86 + x261 + x262 + x263 + x264 + x265 + x266 + x27*x86 + x34*x69 + x34*x71 + x34*x73 + x34*M[26] + x86*M[13] + M[71];
#pragma omp atomic
Ms[72] += x1*x177 + x149*z + x177*M[4] + x213*x273 + x267 + x268 + x269 + x270 + x271 + x272 + x29*x86 + x30*x86 + x31*x86 + x34*x75 + x34*x77 + x34*x79 + x34*M[27] + x41*x64 + x41*x65 + x41*x66 + x41*M[25] + x86*M[14] + M[72];
#pragma omp atomic
Ms[73] += x111*x280 + x155*z + x177*x3 + x177*M[5] + x18*x96 + x19*x96 + x20*x96 + x274 + x275 + x276 + x277 + x278 + x279 + x34*x81 + x34*x82 + x34*x83 + x34*M[28] + x41*x69 + x41*x70 + x41*x72 + x41*M[26] + x96*M[12] + M[73];
#pragma omp atomic
Ms[74] += x1*x184 + x162*z + x184*x2 + x184*M[4] + x23*x96 + x24*x96 + x26*x96 + x281 + x282 + x283 + x284 + x285 + x286 + x41*x75 + x41*x76 + x41*x78 + x41*M[27] + x96*M[13] + M[74];
#pragma omp atomic
Ms[75] += x184*x3 + x184*M[5] + x260*x290 + x287 + x288 + x289 + x29*x96 + x41*x81 + x41*M[28] + x96*M[14] + M[75];
#pragma omp atomic
Ms[76] += x170*M[6] + x291 + x292*M[1] + x34*M[29] + x86*M[15] + M[76];
#pragma omp atomic
Ms[77] += x170*x7 + x170*M[7] + x292*M[2] + x293 + x294 + x295 + x34*x88 + x34*M[30] + x36*x86 + x86*M[16] + M[77];
#pragma omp atomic
Ms[78] += x170*x8 + x170*M[8] + x177*M[6] + x296 + x297 + x298 + x299*M[1] + x34*x91 + x34*M[31] + x39*x86 + x41*x85 + x41*M[29] + x86*M[17] + M[78];
#pragma omp atomic
Ms[79] += x177*M[7] + x299*M[2] + x300 + x301 + x302 + x303*M[1] + x33*x96 + x34*x94 + x34*M[32] + x41*x87 + x41*M[30] + x42*x86 + x86*M[18] + x96*M[15] + M[79];
#pragma omp atomic
Ms[80] += x177*M[8] + x184*x5 + x184*M[6] + x303*M[2] + x304 + x305 + x306 + x34*x97 + x34*M[33] + x35*x96 + x41*x90 + x41*M[31] + x96*M[16] + M[80];
#pragma omp atomic
Ms[81] += x184*x6 + x184*M[7] + x307 + x308 + x309 + x310*M[1] + x38*x96 + x41*x93 + x41*M[32] + x96*M[17] + M[81];
#pragma omp atomic
Ms[82] += x184*M[8] + x310*M[2] + x311 + x41*M[33] + x96*M[18] + M[82];
#pragma omp atomic
Ms[83] += x100*M[9] + x11*M[34] + x188*M[3] + x312 + x313*M[0] + x45*M[19] + M[83];
#pragma omp atomic
Ms[84] += x100*x13 + x100*M[10] + x102*x11 + x11*M[35] + x188*x2 + x188*M[4] + x313*M[1] + x314 + x315 + x316 + x45*x47 + x45*M[20] + M[84];
#pragma omp atomic
Ms[85] += x100*x16 + x100*M[11] + x105*x11 + x11*M[36] + x186*z + x188*x4 + x188*M[5] + x313*M[2] + x317 + x45*x50 + x45*M[21] + z*M[55] + M[85];
#pragma omp atomic
Ms[86] += x100*x19 + x100*M[12] + x108*x11 + x11*M[37] + x111*x321 + x125*M[9] + x188*x5 + x188*M[6] + x214*M[3] + x318 + x319 + x320 + x34*x98 + x34*M[34] + x45*x53 + x45*M[22] + M[86];
#pragma omp atomic
Ms[87] += x100*x24 + x100*x25 + x100*x28 + x100*M[13] + x11*x113 + x11*x114 + x11*x117 + x11*M[38] + x188*x6 + x188*x7 + x188*M[7] + x189*z + x190*z + x191*z + x322 + x323 + x324 + x45*x56 + x45*x57 + x45*x60 + x45*M[23] + z*M[56] + M[87];
#pragma omp atomic
Ms[88] += x100*x30 + x100*M[14] + x11*x119 + x11*M[39] + x121*M[9] + x188*x8 + x188*M[8] + x192*z + x208*M[3] + x325 + x326*M[0] + x41*x98 + x41*M[34] + x45*x62 + x45*M[24] + z*M[57] + M[88];
#pragma omp atomic
Ms[89] += x100*x33 + x100*M[15] + x101*x34 + x11*x123 + x11*M[40] + x125*M[10] + x213*x330 + x214*M[4] + x234*M[3] + x327 + x328 + x329 + x331*M[1] + x34*M[35] + x43*x86 + x45*x65 + x45*M[25] + x86*M[19] + M[89];
#pragma omp atomic
Ms[90] += x100*x35 + x100*x36 + x100*x37 + x100*M[16] + x104*x34 + x105*x34 + x106*x34 + x11*x127 + x11*x128 + x11*x131 + x11*M[41] + x125*x16 + x125*M[11] + x195*z + x196*z + x197*z + x214*x4 + x214*M[5] + x331*M[2] + x332 + x333 + x334 + x34*M[36] + x45*x70 + x45*x71 + x45*x74 + x45*M[26] + z*M[58] + M[90];
#pragma omp atomic
Ms[91] += x100*x38 + x100*x39 + x100*x40 + x100*M[17] + x101*x41 + x102*x41 + x103*x41 + x11*x133 + x11*x134 + x11*x137 + x11*M[42] + x121*x13 + x121*M[10] + x199*z + x2*x208 + x200*z + x202*z + x208*M[4] + x326*M[1] + x335 + x336 + x337 + x41*M[35] + x45*x76 + x45*x77 + x45*x80 + x45*M[27] + z*M[59] + M[91];
#pragma omp atomic
Ms[92] += x100*x42 + x100*M[18] + x104*x41 + x11*x139 + x11*M[43] + x121*M[11] + x205*z + x208*M[5] + x230*M[3] + x326*M[2] + x338 + x339*M[0] + x41*M[36] + x43*x96 + x45*x82 + x45*M[28] + x96*M[19] + z*M[60] + M[92];
#pragma omp atomic
Ms[93] += x107*x34 + x11*x142 + x11*M[44] + x125*M[12] + x170*x9 + x170*M[9] + x214*M[6] + x234*M[4] + x34*M[37] + x340 + x341 + x342 + x343*x344 + x345*M[1] + x45*x85 + x45*M[29] + x46*x86 + x86*M[20] + M[93];
#pragma omp atomic
Ms[94] += x11*x147 + x11*x148 + x11*x151 + x11*M[45] + x112*x34 + x114*x34 + x116*x34 + x125*x25 + x125*M[13] + x209*z + x210*z + x211*z + x214*x7 + x214*M[7] + x234*x4 + x234*M[5] + x34*M[38] + x345*M[2] + x346 + x347 + x348 + x45*x87 + x45*x88 + x45*x89 + x45*M[30] + x49*x86 + x50*x86 + x51*x86 + x86*M[21] + z*M[61] + M[94];
#pragma omp atomic
Ms[95] += x107*x41 + x108*x41 + x109*x41 + x11*x153 + x11*x154 + x11*x157 + x11*M[46] + x111*x353 + x118*x34 + x119*x34 + x120*x34 + x121*x19 + x121*M[12] + x125*x30 + x125*M[14] + x177*x9 + x177*M[9] + x208*x5 + x208*M[6] + x214*x8 + x214*M[8] + x215*z + x216*z + x218*z + x34*M[39] + x349 + x350 + x351 + x41*M[37] + x45*x90 + x45*x91 + x45*x92 + x45*M[31] + z*M[62] + M[95];
#pragma omp atomic
Ms[96] += x11*x160 + x11*x161 + x11*x164 + x11*M[47] + x112*x41 + x113*x41 + x115*x41 + x121*x24 + x121*M[13] + x2*x230 + x208*x6 + x208*M[7] + x221*z + x222*z + x224*z + x230*M[4] + x339*M[1] + x354 + x355 + x356 + x41*M[38] + x45*x93 + x45*x94 + x45*x95 + x45*M[32] + x46*x96 + x47*x96 + x48*x96 + x96*M[20] + z*M[63] + M[96];
#pragma omp atomic
Ms[97] += x11*x166 + x11*M[48] + x118*x41 + x121*M[14] + x184*x9 + x184*M[9] + x208*M[8] + x227*z + x230*M[5] + x339*M[2] + x357 + x358*M[0] + x41*M[39] + x45*x97 + x45*M[33] + x49*x96 + x96*M[21] + z*M[64] + M[97];
#pragma omp atomic
Ms[98] += x0*x292 + x11*x169 + x11*M[49] + x12*x170 + x122*x34 + x125*M[15] + x170*M[10] + x234*M[6] + x292*M[3] + x34*M[40] + x359 + x360 + x361 + x362*M[1] + x52*x86 + x86*M[22] + M[98];
#pragma omp atomic
Ms[99] += x11*x171 + x11*x172 + x11*x173 + x11*M[50] + x125*x36 + x125*M[16] + x126*x34 + x128*x34 + x130*x34 + x15*x170 + x16*x170 + x17*x170 + x170*M[11] + x231*z + x232*z + x233*z + x234*x7 + x234*M[7] + x34*M[41] + x362*M[2] + x363 + x364 + x365 + x55*x86 + x57*x86 + x59*x86 + x86*M[23] + z*M[65] + M[99];
#pragma omp atomic
Ms[100] += x0*x299 + x11*x174 + x11*x175 + x11*x176 + x11*M[51] + x12*x177 + x121*x33 + x121*M[15] + x122*x41 + x123*x41 + x124*x41 + x125*x39 + x125*M[17] + x132*x34 + x134*x34 + x136*x34 + x177*M[10] + x234*x8 + x234*M[8] + x235*z + x236*z + x238*z + x299*M[3] + x34*M[42] + x366 + x367 + x368 + x369*M[1] + x41*M[40] + x61*x86 + x62*x86 + x63*x86 + x86*M[24] + z*M[66] + M[100];
#pragma omp atomic
Ms[101] += x0*x303 + x11*x178 + x11*x179 + x11*x180 + x11*M[52] + x121*x35 + x121*M[16] + x125*x42 + x125*M[18] + x126*x41 + x127*x41 + x129*x41 + x138*x34 + x139*x34 + x140*x34 + x15*x177 + x177*M[11] + x230*x5 + x230*M[6] + x241*z + x242*z + x244*z + x303*M[3] + x34*M[43] + x369*M[2] + x370 + x371 + x372 + x41*M[41] + x52*x96 + x53*x96 + x54*x96 + x96*M[22] + z*M[67] + M[101];
#pragma omp atomic
Ms[102] += x11*x181 + x11*x182 + x11*x183 + x11*M[53] + x12*x184 + x121*x38 + x121*M[17] + x13*x184 + x132*x41 + x133*x41 + x135*x41 + x14*x184 + x184*M[10] + x230*x6 + x230*M[7] + x247*z + x248*z + x250*z + x358*M[1] + x373 + x374 + x375 + x41*M[42] + x55*x96 + x56*x96 + x58*x96 + x96*M[23] + z*M[68] + M[102];
#pragma omp atomic
Ms[103] += x0*x310 + x11*x185 + x11*M[54] + x121*M[18] + x138*x41 + x15*x184 + x184*M[11] + x230*M[8] + x253*z + x310*M[3] + x358*M[2] + x376 + x41*M[43] + x61*x96 + x96*M[24] + z*M[69] + M[103];
#pragma omp atomic
Ms[104] += x1*x292 + x141*x34 + x170*x18 + x170*M[12] + x292*M[4] + x34*M[44] + x377 + x378 + x379 + x380*x381 + x64*x86 + x86*M[25] + M[104];
#pragma omp atomic
Ms[105] += x146*x34 + x148*x34 + x150*x34 + x170*x23 + x170*x25 + x170*x27 + x170*M[13] + x256*z + x257*z + x258*z + x292*x3 + x292*x4 + x292*M[5] + x34*M[45] + x382 + x383 + x384 + x69*x86 + x71*x86 + x73*x86 + x86*M[26] + z*M[70] + M[105];
#pragma omp atomic
Ms[106] += x1*x299 + x141*x41 + x142*x41 + x143*x41 + x152*x34 + x154*x34 + x156*x34 + x170*x29 + x170*x30 + x170*x31 + x170*M[14] + x177*x18 + x177*M[12] + x261*z + x262*z + x264*z + x299*M[4] + x34*M[46] + x344*x388 + x385 + x386 + x387 + x41*M[44] + x75*x86 + x77*x86 + x79*x86 + x86*M[27] + z*M[71] + M[106];
#pragma omp atomic
Ms[107] += x1*x303 + x146*x41 + x147*x41 + x149*x41 + x159*x34 + x161*x34 + x163*x34 + x177*x23 + x177*M[13] + x213*x392 + x267*z + x268*z + x270*z + x299*x3 + x299*M[5] + x303*M[4] + x34*M[47] + x389 + x390 + x391 + x41*M[45] + x64*x96 + x65*x96 + x66*x96 + x81*x86 + x82*x86 + x83*x86 + x86*M[28] + x96*M[25] + z*M[72] + M[107];
#pragma omp atomic
Ms[108] += x111*x396 + x152*x41 + x153*x41 + x155*x41 + x165*x34 + x166*x34 + x167*x34 + x177*x29 + x177*M[14] + x18*x184 + x184*x19 + x184*x20 + x184*M[12] + x274*z + x275*z + x277*z + x3*x303 + x303*M[5] + x34*M[48] + x393 + x394 + x395 + x41*M[46] + x69*x96 + x70*x96 + x72*x96 + x96*M[26] + z*M[73] + M[108];
#pragma omp atomic
Ms[109] += x1*x310 + x159*x41 + x160*x41 + x162*x41 + x184*x23 + x184*x24 + x184*x26 + x184*M[13] + x2*x310 + x281*z + x282*z + x284*z + x310*M[4] + x397 + x398 + x399 + x41*M[47] + x75*x96 + x76*x96 + x78*x96 + x96*M[27] + z*M[74] + M[109];
#pragma omp atomic
Ms[110] += x165*x41 + x184*x29 + x184*M[14] + x287*z + x3*x310 + x310*M[5] + x381*x401 + x400 + x41*M[48] + x81*x96 + x96*M[28] + z*M[75] + M[110];
#pragma omp atomic
Ms[111] += x170*M[15] + x292*M[6] + x34*M[49] + x402 + x403*M[1] + x86*M[29] + M[111];
#pragma omp atomic
Ms[112] += x170*x36 + x170*M[16] + x172*x34 + x291*z + x292*x7 + x292*M[7] + x34*M[50] + x403*M[2] + x404 + x86*x88 + x86*M[30] + z*M[76] + M[112];
#pragma omp atomic
Ms[113] += x169*x41 + x170*x39 + x170*M[17] + x175*x34 + x177*M[15] + x292*x8 + x292*M[8] + x293*z + x299*M[6] + x34*M[51] + x405 + x406*M[1] + x41*M[49] + x86*x91 + x86*M[31] + z*M[77] + M[113];
#pragma omp atomic
Ms[114] += x170*x42 + x170*M[18] + x171*x41 + x177*M[16] + x179*x34 + x296*z + x299*M[7] + x303*M[6] + x34*M[52] + x406*M[2] + x407 + x408*M[1] + x41*M[50] + x85*x96 + x86*x94 + x86*M[32] + x96*M[29] + z*M[78] + M[114];
#pragma omp atomic
Ms[115] += x174*x41 + x177*M[17] + x182*x34 + x184*x33 + x184*M[15] + x299*M[8] + x300*z + x303*M[7] + x34*M[53] + x408*M[2] + x409 + x41*M[51] + x410*M[1] + x86*x97 + x86*M[33] + x87*x96 + x96*M[30] + z*M[79] + M[115];
#pragma omp atomic
Ms[116] += x177*M[18] + x178*x41 + x184*x35 + x184*M[16] + x185*x34 + x303*M[8] + x304*z + x310*x5 + x310*M[6] + x34*M[54] + x41*M[52] + x410*M[2] + x411 + x90*x96 + x96*M[31] + z*M[80] + M[116];
#pragma omp atomic
Ms[117] += x181*x41 + x184*x38 + x184*M[17] + x307*z + x310*x6 + x310*M[7] + x41*M[53] + x412 + x413*M[1] + x93*x96 + x96*M[32] + z*M[81] + M[117];
#pragma omp atomic
Ms[118] += x184*M[18] + x310*M[8] + x41*M[54] + x413*M[2] + x96*M[33] + z*M[82] + M[118];
#pragma omp atomic
Ms[119] += x*M[83] + x100*M[19] + x11*M[55] + x188*M[9] + x313*M[3] + x414*M[0] + x45*M[34] + M[119];
#pragma omp atomic
Ms[120] += x*M[84] + x100*x47 + x100*M[20] + x102*x45 + x11*x190 + x11*M[56] + x13*x188 + x188*M[10] + x2*x313 + x312*y + x313*M[4] + x414*M[1] + x45*M[35] + y*M[83] + M[120];
#pragma omp atomic
Ms[121] += x*M[85] + x100*x50 + x100*M[21] + x105*x45 + x11*x193 + x11*M[57] + x16*x188 + x188*M[11] + x312*z + x313*x4 + x313*M[5] + x414*M[2] + x45*M[36] + z*M[83] + M[121];
#pragma omp atomic
Ms[122] += x*M[86] + x100*x53 + x100*M[22] + x108*x45 + x11*x196 + x11*M[58] + x111*x415 + x125*M[19] + x186*x34 + x188*x19 + x188*M[12] + x214*M[9] + x313*x5 + x313*M[6] + x314*y + x331*M[3] + x34*M[55] + x45*M[37] + y*M[84] + M[122];
#pragma omp atomic
Ms[123] += x*M[87] + x100*x56 + x100*x57 + x100*x60 + x100*M[23] + x11*x200 + x11*x201 + x11*x204 + x11*M[59] + x113*x45 + x114*x45 + x117*x45 + x188*x24 + x188*x25 + x188*x28 + x188*M[13] + x313*x6 + x313*x7 + x313*M[7] + x314*z + x315*z + x316*z + x317*y + x45*M[38] + y*M[85] + z*M[84] + M[123];
#pragma omp atomic
Ms[124] += x*M[88] + x100*x62 + x100*M[24] + x11*x206 + x11*M[60] + x119*x45 + x121*M[19] + x186*x41 + x188*x30 + x188*M[14] + x208*M[9] + x313*x8 + x313*M[8] + x317*z + x326*M[3] + x41*M[55] + x416*M[0] + x45*M[39] + z*M[85] + M[124];
#pragma omp atomic
Ms[125] += x*M[89] + x100*x65 + x100*M[25] + x11*x210 + x11*M[61] + x123*x45 + x125*M[20] + x188*x33 + x188*M[15] + x189*x34 + x213*x417 + x214*M[10] + x234*M[9] + x318*y + x331*M[4] + x34*M[56] + x345*M[3] + x418*M[1] + x45*M[40] + x86*x98 + x86*M[34] + y*M[86] + M[125];
#pragma omp atomic
Ms[126] += x*M[90] + x100*x70 + x100*x71 + x100*x74 + x100*M[26] + x11*x216 + x11*x217 + x11*x220 + x11*M[62] + x125*x50 + x125*M[21] + x127*x45 + x128*x45 + x131*x45 + x16*x214 + x188*x35 + x188*x36 + x188*x37 + x188*M[16] + x192*x34 + x193*x34 + x194*x34 + x214*M[11] + x318*z + x319*z + x320*z + x322*y + x331*x4 + x331*M[5] + x34*M[57] + x418*M[2] + x45*M[41] + y*M[87] + z*M[86] + M[126];
#pragma omp atomic
Ms[127] += x*M[91] + x100*x76 + x100*x77 + x100*x80 + x100*M[27] + x11*x222 + x11*x223 + x11*x226 + x11*M[63] + x121*x47 + x121*M[20] + x13*x208 + x133*x45 + x134*x45 + x137*x45 + x188*x38 + x188*x39 + x188*x40 + x188*M[17] + x189*x41 + x190*x41 + x191*x41 + x2*x326 + x208*M[10] + x322*z + x323*z + x324*z + x325*y + x326*M[4] + x41*M[56] + x416*M[1] + x45*M[42] + y*M[88] + z*M[87] + M[127];
#pragma omp atomic
Ms[128] += x*M[92] + x100*x82 + x100*M[28] + x11*x228 + x11*M[64] + x121*M[21] + x139*x45 + x188*x42 + x188*M[18] + x192*x41 + x208*M[11] + x230*M[9] + x325*z + x326*M[5] + x339*M[3] + x41*M[57] + x416*M[2] + x419*M[0] + x45*M[43] + x96*x98 + x96*M[34] + z*M[88] + M[128];
#pragma omp atomic
Ms[129] += x*M[93] + x100*x85 + x100*M[29] + x101*x86 + x11*x232 + x11*M[65] + x125*M[22] + x142*x45 + x170*x43 + x170*M[19] + x195*x34 + x214*M[12] + x234*M[10] + x327*y + x331*M[6] + x34*M[58] + x344*x420 + x345*M[4] + x362*M[3] + x421*M[1] + x45*M[44] + x86*M[35] + y*M[89] + M[129];
#pragma omp atomic
Ms[130] += x*M[94] + x100*x87 + x100*x88 + x100*x89 + x100*M[30] + x104*x86 + x105*x86 + x106*x86 + x11*x236 + x11*x237 + x11*x240 + x11*M[66] + x125*x57 + x125*M[23] + x147*x45 + x148*x45 + x151*x45 + x16*x234 + x199*x34 + x201*x34 + x203*x34 + x214*x25 + x214*M[13] + x234*M[11] + x327*z + x328*z + x329*z + x331*x7 + x331*M[7] + x332*y + x34*M[59] + x345*x4 + x345*M[5] + x421*M[2] + x45*M[45] + x86*M[36] + y*M[90] + z*M[89] + M[130];
#pragma omp atomic
Ms[131] += x*M[95] + x100*x90 + x100*x91 + x100*x92 + x100*M[31] + x11*x242 + x11*x243 + x11*x246 + x11*M[67] + x121*x53 + x121*M[22] + x125*x62 + x125*M[24] + x153*x45 + x154*x45 + x157*x45 + x177*x43 + x177*M[19] + x19*x208 + x195*x41 + x196*x41 + x197*x41 + x205*x34 + x206*x34 + x207*x34 + x208*M[12] + x214*x30 + x214*M[14] + x326*x5 + x326*M[6] + x331*x8 + x331*M[8] + x332*z + x333*z + x334*z + x335*y + x34*M[60] + x369*M[3] + x41*M[58] + x422*x423 + x45*M[46] + y*M[91] + z*M[90] + M[131];
#pragma omp atomic
Ms[132] += x*M[96] + x100*x93 + x100*x94 + x100*x95 + x100*M[32] + x101*x96 + x102*x96 + x103*x96 + x11*x248 + x11*x249 + x11*x252 + x11*M[68] + x121*x56 + x121*M[23] + x13*x230 + x160*x45 + x161*x45 + x164*x45 + x199*x41 + x2*x339 + x200*x41 + x202*x41 + x208*x24 + x208*M[13] + x230*M[10] + x326*x6 + x326*M[7] + x335*z + x336*z + x337*z + x338*y + x339*M[4] + x41*M[59] + x419*M[1] + x45*M[47] + x96*M[35] + y*M[92] + z*M[91] + M[132];
#pragma omp atomic
Ms[133] += x*M[97] + x100*x97 + x100*M[33] + x104*x96 + x11*x254 + x11*M[69] + x121*M[24] + x166*x45 + x184*x43 + x184*M[19] + x205*x41 + x208*M[14] + x230*M[11] + x326*M[8] + x338*z + x339*M[5] + x358*M[3] + x41*M[60] + x419*M[2] + x424*M[0] + x45*M[48] + x96*M[36] + z*M[92] + M[133];
#pragma omp atomic
Ms[134] += x*M[98] + x107*x86 + x11*x257 + x11*M[70] + x125*M[25] + x169*x45 + x170*x46 + x170*M[20] + x209*x34 + x214*M[15] + x234*M[12] + x259*x426 + x292*x9 + x292*M[9] + x34*M[61] + x340*y + x345*M[6] + x362*M[4] + x427*M[1] + x45*M[49] + x86*M[37] + y*M[93] + M[134];
#pragma omp atomic
Ms[135] += x*M[99] + x11*x262 + x11*x263 + x11*x266 + x11*M[71] + x112*x86 + x114*x86 + x116*x86 + x125*x71 + x125*M[26] + x170*x49 + x170*x50 + x170*x51 + x170*M[21] + x171*x45 + x172*x45 + x173*x45 + x214*x36 + x214*M[16] + x215*x34 + x217*x34 + x219*x34 + x234*x25 + x234*M[13] + x34*M[62] + x340*z + x341*z + x342*z + x345*x7 + x345*M[7] + x346*y + x362*x4 + x362*M[5] + x427*M[2] + x45*M[50] + x86*M[38] + y*M[94] + z*M[93] + M[135];
#pragma omp atomic
Ms[136] += x*M[100] + x11*x268 + x11*x269 + x11*x272 + x11*M[72] + x118*x86 + x119*x86 + x120*x86 + x121*x65 + x121*M[25] + x125*x77 + x125*M[27] + x145*x428 + x174*x45 + x175*x45 + x176*x45 + x177*x46 + x177*M[20] + x208*x33 + x208*M[15] + x209*x41 + x210*x41 + x211*x41 + x214*x39 + x214*M[17] + x221*x34 + x223*x34 + x225*x34 + x234*x30 + x234*M[14] + x299*x9 + x299*M[9] + x34*M[63] + x345*x8 + x345*M[8] + x346*z + x347*z + x348*z + x349*y + x369*M[4] + x41*M[61] + x429*x430 + x45*M[51] + x86*M[39] + y*M[95] + z*M[94] + M[136];
#pragma omp atomic
Ms[137] += x*M[101] + x107*x96 + x108*x96 + x109*x96 + x11*x275 + x11*x276 + x11*x279 + x11*M[73] + x121*x70 + x121*M[26] + x125*x82 + x125*M[28] + x177*x49 + x177*M[21] + x178*x45 + x179*x45 + x180*x45 + x19*x230 + x208*x35 + x208*M[16] + x214*x42 + x214*M[18] + x215*x41 + x216*x41 + x218*x41 + x227*x34 + x228*x34 + x229*x34 + x230*M[12] + x303*x9 + x303*M[9] + x339*x5 + x339*M[6] + x34*M[64] + x349*z + x350*z + x351*z + x354*y + x369*M[5] + x41*M[62] + x422*x431 + x430*x432 + x45*M[52] + x96*M[37] + y*M[96] + z*M[95] + M[137];
#pragma omp atomic
Ms[138] += x*M[102] + x11*x282 + x11*x283 + x11*x286 + x11*M[74] + x112*x96 + x113*x96 + x115*x96 + x121*x76 + x121*M[27] + x181*x45 + x182*x45 + x183*x45 + x184*x46 + x184*x47 + x184*x48 + x184*M[20] + x2*x358 + x208*x38 + x208*M[17] + x221*x41 + x222*x41 + x224*x41 + x230*x24 + x230*M[13] + x339*x6 + x339*M[7] + x354*z + x355*z + x356*z + x357*y + x358*M[4] + x41*M[63] + x424*M[1] + x45*M[53] + x96*M[38] + y*M[97] + z*M[96] + M[138];
#pragma omp atomic
Ms[139] += x*M[103] + x11*x288 + x11*M[75] + x118*x96 + x121*M[28] + x184*x49 + x184*M[21] + x185*x45 + x208*M[18] + x227*x41 + x230*M[14] + x290*x426 + x310*x9 + x310*M[9] + x339*M[8] + x357*z + x358*M[5] + x41*M[64] + x424*M[2] + x45*M[54] + x96*M[39] + z*M[97] + M[139];
#pragma omp atomic
Ms[140] += x*M[104] + x0*x403 + x11*x291 + x11*M[76] + x12*x292 + x122*x86 + x125*M[29] + x170*x52 + x170*M[22] + x231*x34 + x234*M[15] + x292*M[10] + x34*M[65] + x359*y + x362*M[6] + x403*M[3] + x433*M[1] + x86*M[40] + y*M[98] + M[140];
#pragma omp atomic
Ms[141] += x*M[105] + x11*x293 + x11*x294 + x11*x295 + x11*M[77] + x125*x88 + x125*M[30] + x126*x86 + x128*x86 + x130*x86 + x15*x292 + x16*x292 + x17*x292 + x170*x55 + x170*x57 + x170*x59 + x170*M[23] + x234*x36 + x234*M[16] + x235*x34 + x237*x34 + x239*x34 + x292*M[11] + x34*M[66] + x359*z + x360*z + x361*z + x362*x7 + x362*M[7] + x363*y + x433*M[2] + x86*M[41] + y*M[99] + z*M[98] + M[141];
#pragma omp atomic
Ms[142] += x*M[106] + x0*x406 + x11*x296 + x11*x297 + x11*x298 + x11*M[78] + x12*x299 + x121*x85 + x121*M[29] + x125*x91 + x125*M[31] + x132*x86 + x134*x86 + x136*x86 + x170*x61 + x170*x62 + x170*x63 + x170*M[24] + x177*x52 + x177*M[22] + x231*x41 + x232*x41 + x233*x41 + x234*x39 + x234*M[17] + x241*x34 + x243*x34 + x245*x34 + x299*M[10] + x34*M[67] + x362*x8 + x362*M[8] + x363*z + x364*z + x365*z + x366*y + x369*M[6] + x406*M[3] + x41*M[65] + x428*x429 + x86*M[42] + y*M[100] + z*M[99] + M[142];
#pragma omp atomic
Ms[143] += x*M[107] + x0*x408 + x11*x300 + x11*x301 + x11*x302 + x11*M[79] + x12*x303 + x121*x87 + x121*M[30] + x122*x96 + x123*x96 + x124*x96 + x125*x94 + x125*M[32] + x138*x86 + x139*x86 + x140*x86 + x15*x299 + x177*x55 + x177*M[23] + x230*x33 + x230*M[15] + x234*x42 + x234*M[18] + x235*x41 + x236*x41 + x238*x41 + x247*x34 + x249*x34 + x251*x34 + x299*M[11] + x303*M[10] + x34*M[68] + x366*z + x367*z + x368*z + x369*M[7] + x370*y + x408*M[3] + x41*M[66] + x428*x432 + x429*x434 + x86*M[43] + x96*M[40] + y*M[101] + z*M[100] + M[143];
#pragma omp atomic
Ms[144] += x*M[108] + x0*x410 + x11*x304 + x11*x305 + x11*x306 + x11*M[80] + x121*x90 + x121*M[31] + x125*x97 + x125*M[33] + x126*x96 + x127*x96 + x129*x96 + x15*x303 + x177*x61 + x177*M[24] + x184*x52 + x184*x53 + x184*x54 + x184*M[22] + x230*x35 + x230*M[16] + x241*x41 + x242*x41 + x244*x41 + x253*x34 + x254*x34 + x255*x34 + x303*M[11] + x34*M[69] + x358*x5 + x358*M[6] + x369*M[8] + x370*z + x371*z + x372*z + x373*y + x41*M[67] + x410*M[3] + x432*x434 + x96*M[41] + y*M[102] + z*M[101] + M[144];
#pragma omp atomic
Ms[145] += x*M[109] + x11*x307 + x11*x308 + x11*x309 + x11*M[81] + x12*x310 + x121*x93 + x121*M[32] + x13*x310 + x132*x96 + x133*x96 + x135*x96 + x14*x310 + x184*x55 + x184*x56 + x184*x58 + x184*M[23] + x230*x38 + x230*M[17] + x247*x41 + x248*x41 + x250*x41 + x310*M[10] + x358*x6 + x358*M[7] + x373*z + x374*z + x375*z + x376*y + x41*M[68] + x435*M[1] + x96*M[42] + y*M[103] + z*M[102] + M[145];
#pragma omp atomic
Ms[146] += x*M[110] + x0*x413 + x11*x311 + x11*M[82] + x121*M[33] + x138*x96 + x15*x310 + x184*x61 + x184*M[24] + x230*M[18] + x253*x41 + x310*M[11] + x358*M[8] + x376*z + x41*M[69] + x413*M[3] + x435*M[2] + x96*M[43] + z*M[103] + M[146];
#pragma omp atomic
Ms[147] += x*M[111] + x1*x403 + x141*x86 + x170*x64 + x170*M[25] + x18*x292 + x256*x34 + x292*M[12] + x34*M[70] + x377*y + x403*M[4] + x436*x437 + x86*M[44] + y*M[104] + M[147];
#pragma omp atomic
Ms[148] += x*M[112] + x146*x86 + x148*x86 + x150*x86 + x170*x69 + x170*x71 + x170*x73 + x170*M[26] + x23*x292 + x25*x292 + x261*x34 + x263*x34 + x265*x34 + x27*x292 + x292*M[13] + x3*x403 + x34*M[71] + x377*z + x378*z + x379*z + x382*y + x4*x403 + x403*M[5] + x86*M[45] + y*M[105] + z*M[104] + M[148];
#pragma omp atomic
Ms[149] += x*M[113] + x1*x406 + x152*x86 + x154*x86 + x156*x86 + x170*x75 + x170*x77 + x170*x79 + x170*M[27] + x177*x64 + x177*M[25] + x18*x299 + x256*x41 + x257*x41 + x258*x41 + x267*x34 + x269*x34 + x271*x34 + x29*x292 + x292*x30 + x292*x31 + x292*M[14] + x299*M[12] + x34*M[72] + x382*z + x383*z + x384*z + x385*y + x406*M[4] + x41*M[70] + x438*M[0] + x86*M[46] + y*M[106] + z*M[105] + M[149];
#pragma omp atomic
Ms[150] += x*M[114] + x1*x408 + x141*x96 + x142*x96 + x143*x96 + x159*x86 + x161*x86 + x163*x86 + x170*x81 + x170*x82 + x170*x83 + x170*M[28] + x177*x69 + x177*M[26] + x18*x303 + x23*x299 + x261*x41 + x262*x41 + x264*x41 + x274*x34 + x276*x34 + x278*x34 + x299*M[13] + x3*x406 + x303*M[12] + x34*M[73] + x344*x439 + x385*z + x386*z + x387*z + x389*y + x406*M[5] + x408*M[4] + x41*M[71] + x86*M[47] + x96*M[44] + y*M[107] + z*M[106] + M[150];
#pragma omp atomic
Ms[151] += x*M[115] + x1*x410 + x146*x96 + x147*x96 + x149*x96 + x165*x86 + x166*x86 + x167*x86 + x177*x75 + x177*M[27] + x184*x64 + x184*x65 + x184*x66 + x184*M[25] + x213*x440 + x23*x303 + x267*x41 + x268*x41 + x270*x41 + x281*x34 + x283*x34 + x285*x34 + x29*x299 + x299*M[14] + x3*x408 + x303*M[13] + x34*M[74] + x389*z + x390*z + x391*z + x393*y + x408*M[5] + x41*M[72] + x410*M[4] + x86*M[48] + x96*M[45] + y*M[108] + z*M[107] + M[151];
#pragma omp atomic
Ms[152] += x*M[116] + x111*x441 + x152*x96 + x153*x96 + x155*x96 + x177*x81 + x177*M[28] + x18*x310 + x184*x69 + x184*x70 + x184*x72 + x184*M[26] + x19*x310 + x20*x310 + x274*x41 + x275*x41 + x277*x41 + x287*x34 + x288*x34 + x289*x34 + x29*x303 + x3*x410 + x303*M[14] + x310*M[12] + x34*M[75] + x393*z + x394*z + x395*z + x397*y + x41*M[73] + x410*M[5] + x96*M[46] + y*M[109] + z*M[108] + M[152];
#pragma omp atomic
Ms[153] += x*M[117] + x1*x413 + x159*x96 + x160*x96 + x162*x96 + x184*x75 + x184*x76 + x184*x78 + x184*M[27] + x2*x413 + x23*x310 + x24*x310 + x26*x310 + x281*x41 + x282*x41 + x284*x41 + x310*M[13] + x397*z + x398*z + x399*z + x400*y + x41*M[74] + x413*M[4] + x96*M[47] + y*M[110] + z*M[109] + M[153];
#pragma omp atomic
Ms[154] += x*M[118] + x165*x96 + x184*x81 + x184*M[28] + x287*x41 + x29*x310 + x3*x413 + x310*M[14] + x400*z + x41*M[75] + x413*M[5] + x437*x442 + x96*M[48] + z*M[110] + M[154];
#pragma omp atomic
Ms[155] += x170*M[29] + x292*M[15] + x34*M[76] + x403*M[6] + x443*M[1] + x86*M[49] + y*M[111] + M[155];
#pragma omp atomic
Ms[156] += x170*x88 + x170*M[30] + x172*x86 + x292*x36 + x292*M[16] + x294*x34 + x34*M[77] + x402*z + x403*x7 + x403*M[7] + x443*M[2] + x86*M[50] + y*M[112] + z*M[111] + M[156];
#pragma omp atomic
Ms[157] += x170*x91 + x170*M[31] + x175*x86 + x177*M[29] + x291*x41 + x292*x39 + x292*M[17] + x297*x34 + x299*M[15] + x34*M[78] + x403*x8 + x403*M[8] + x404*z + x406*M[6] + x41*M[76] + x438*M[1] + x86*M[51] + y*M[113] + z*M[112] + M[157];
#pragma omp atomic
Ms[158] += x169*x96 + x170*x94 + x170*M[32] + x177*M[30] + x179*x86 + x292*x42 + x292*M[18] + x293*x41 + x299*M[16] + x301*x34 + x303*M[15] + x34*M[79] + x405*z + x406*M[7] + x408*M[6] + x41*M[77] + x438*M[2] + x444*M[1] + x86*M[52] + x96*M[49] + y*M[114] + z*M[113] + M[158];
#pragma omp atomic
Ms[159] += x170*x97 + x170*M[33] + x171*x96 + x177*M[31] + x182*x86 + x184*x85 + x184*M[29] + x296*x41 + x299*M[17] + x303*M[16] + x305*x34 + x34*M[80] + x406*M[8] + x407*z + x408*M[7] + x41*M[78] + x410*M[6] + x444*M[2] + x445*M[1] + x86*M[53] + x96*M[50] + y*M[115] + z*M[114] + M[159];
#pragma omp atomic
Ms[160] += x174*x96 + x177*M[32] + x184*x87 + x184*M[30] + x185*x86 + x299*M[18] + x300*x41 + x303*M[17] + x308*x34 + x310*x33 + x310*M[15] + x34*M[81] + x408*M[8] + x409*z + x41*M[79] + x410*M[7] + x445*M[2] + x446*M[1] + x86*M[54] + x96*M[51] + y*M[116] + z*M[115] + M[160];
#pragma omp atomic
Ms[161] += x177*M[33] + x178*x96 + x184*x90 + x184*M[31] + x303*M[18] + x304*x41 + x310*x35 + x310*M[16] + x311*x34 + x34*M[82] + x41*M[80] + x410*M[8] + x411*z + x413*x5 + x413*M[6] + x446*M[2] + x96*M[52] + y*M[117] + z*M[116] + M[161];
#pragma omp atomic
Ms[162] += x181*x96 + x184*x93 + x184*M[32] + x307*x41 + x310*x38 + x310*M[17] + x41*M[81] + x412*z + x413*x6 + x413*M[7] + x447*M[1] + x96*M[53] + y*M[118] + z*M[117] + M[162];
#pragma omp atomic
Ms[163] += x184*M[33] + x310*M[18] + x41*M[82] + x413*M[8] + x447*M[2] + x96*M[54] + z*M[118] + M[163];

}

void M2L_8(double x, double y, double z, double * M, double * L) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = pow(R, -5);
x3 = 3.0*x2;
x4 = x*y;
x5 = x3*x4;
x6 = x*z;
x7 = x3*x6;
x8 = y*z;
x9 = x3*x8;
x10 = pow(R, -7);
x11 = 15.0*x10;
x12 = x4*z;
x13 = x11*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x3;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x3;
x20 = x14 + x19;
x21 = 9.0*x2;
x22 = x11*x15;
x23 = -x22;
x24 = x*(x21 + x23);
x25 = x23 + x3;
x26 = x25*y;
x27 = x11*x18;
x28 = -x27;
x29 = y*(x21 + x28);
x30 = x25*z;
x31 = z*(x28 + x3);
x32 = 1.0*x;
x33 = x32*(x27 - x3);
x34 = 45.0*x10;
x35 = -x34;
x36 = pow(R, -9);
x37 = x15*x36;
x38 = 105.0*x37;
x39 = x35 + x38;
x40 = x39*x4;
x41 = x39*x6;
x42 = x18*x36;
x43 = 105.0*x42;
x44 = x35 + x43;
x45 = x44*x8;
x46 = -x11;
x47 = x8*(x38 + x46);
x48 = x32*y;
x49 = x44*x48;
x50 = x43 + x46;
x51 = x32*z;
x52 = x50*x51;
x53 = 315.0*x36;
x54 = pow(R, -11);
x55 = 945.0*x54;
x56 = x15*x55;
x57 = x53 - x56;
x58 = x12*x57;
x59 = x18*x55;
x60 = x32*x8;
x61 = x60*(-x53 + x59);
x62 = (x*x*x*x);
x63 = 105.0*x36;
x64 = x62*x63;
x65 = 90.0*x10;
x66 = -x15*x65 + x21 + x64;
x67 = (y*y*y*y);
x68 = x63*x67;
x69 = -x18*x65 + x21 + x68;
x70 = 2.0*x0 - x16 - x19;
x71 = -225.0*x10;
x72 = x55*x62;
x73 = -x72;
x74 = x*(1050.0*x37 + x71 + x73);
x75 = x55*x67;
x76 = -x75;
x77 = y*(1050.0*x42 + x71 + x76);
x78 = x35 + 630.0*x37 + x73;
x79 = x78*y;
x80 = x78*z;
x81 = 630.0*x42;
x82 = x35 + x76 + x81;
x83 = x82*z;
x84 = x32*(x34 + x75 - x81);
x85 = 1575.0*x36;
x86 = pow(R, -13);
x87 = x62*x86;
x88 = 10395.0*x87;
x89 = x15*x54;
x90 = x85 + x88 - 9450.0*x89;
x91 = x4*x90;
x92 = x6*x90;
x93 = 5670.0*x89;
x94 = x53 + x88 - x93;
x95 = x8*x94;
x96 = x67*x86;
x97 = 10395.0*x96;
x98 = x18*x54;
x99 = x85 + x97 - 9450.0*x98;
x100 = x8*x99;
x101 = x48*x99;
x102 = x53 + x97 - 5670.0*x98;
x103 = x102*x51;
x104 = 14175.0*x54;
x105 = -x104;
x106 = pow(R, -15);
x107 = 135135.0*x106;
x108 = x107*x62;
x109 = x15*x86;
x110 = 103950.0*x109;
x111 = x105 - x108 + x110;
x112 = x111*x12;
x113 = x107*x67;
x114 = x18*x86;
x115 = 103950.0*x114;
x116 = x60*(x104 + x113 - x115);
x117 = pow(x, 6);
x118 = 10395.0*x86;
x119 = x117*x118;
x120 = -x104*x62 + x119 + 4725.0*x37 + x71;
x121 = pow(y, 6);
x122 = x118*x121;
x123 = -x104*x67 + x122 + 4725.0*x42 + x71;
x124 = 11025.0*x36;
x125 = 99225.0*x54;
x126 = x107*x117;
x127 = -x126;
x128 = x127 + 218295.0*x87;
x129 = x*(x124 - x125*x15 + x128);
x130 = 42525.0*x54;
x131 = x127 - x130*x15 + x85 + 155925.0*x87;
x132 = x131*y;
x133 = x107*x121;
x134 = -x133;
x135 = x134 + 218295.0*x96;
x136 = y*(x124 - x125*x18 + x135);
x137 = x131*z;
x138 = 155925.0*x96;
x139 = x130*x18;
x140 = x134 + x138 - x139 + x85;
x141 = x140*z;
x142 = x32*(x133 - x138 + x139 - x85);
x143 = -x125;
x144 = pow(R, -17);
x145 = 2027025.0*x144;
x146 = x117*x145;
x147 = x106*x62;
x148 = 2837835.0*x147;
x149 = 1091475.0*x109 + x143 + x146 - x148;
x150 = x149*x4;
x151 = x149*x6;
x152 = x121*x145;
x153 = x106*x67;
x154 = 2837835.0*x153;
x155 = 1091475.0*x114 + x143 + x152 - x154;
x156 = x155*x8;
x157 = 2027025.0*x147;
x158 = 467775.0*x109;
x159 = x8*(x105 + x146 - x157 + x158);
x160 = x18*x38;
x161 = x160 + x25 + x28;
x162 = x155*x48;
x163 = 2027025.0*x153;
x164 = x51*(x105 + 467775.0*x114 + x152 - x163);
x165 = x18*x53;
x166 = -x18*x56;
x167 = x*(x165 + x166 + x39);
x168 = x15*x53;
x169 = y*(x166 + x168 + x44);
x170 = z*(x166 + x38 + x50);
x171 = -x24 + x33;
x172 = -x26 - x29;
x173 = -x30 - x31;
x174 = -2835.0*x98;
x175 = x109*x18;
x176 = 10395.0*x175;
x177 = x174 + x176;
x178 = 945.0*x36;
x179 = -2835.0*x89;
x180 = x178 + x179;
x181 = x4*(x177 + x180);
x182 = x6*(x177 + x57);
x183 = x8*(x176 + x179 + x53 - x59);
x184 = 31185.0*x114;
x185 = x15*x18;
x186 = -8505.0*x54;
x187 = 31185.0*x109 + x186;
x188 = x12*(-x107*x185 + x184 + x187);
x189 = -x40 - x49;
x190 = -x41 - x52;
x191 = -x45 - x47;
x192 = pow(x, 8)*x145;
x193 = 3783780.0*x106;
x194 = -x117*x193 + x124 + x192 + 2182950.0*x87 - 396900.0*x89;
x195 = x145*pow(y, 8);
x196 = -x121*x193 + x124 + x195 + 2182950.0*x96 - 396900.0*x98;
x197 = -x58 + x61;
x198 = 105.0*x10;
x199 = -x160 - 12.0*x2;
x200 = x18*x198 + x199 + x22 - x68;
x201 = x15*x198 + x199 + x27 - x64;
x202 = 120.0*x10;
x203 = -x15*x202 - x18*x202 + 210.0*x18*x37 + 24.0*x2 + x64 + x68;
x204 = x15*x97;
x205 = x18*x93;
x206 = -x205;
x207 = x168 + x204 + x206 + x82;
x208 = x18*x88;
x209 = x165 + x206 + x208 + x78;
x210 = 62370.0*x175;
x211 = -x113*x15;
x212 = x210 + x211;
x213 = 31185.0*x96;
x214 = x213 - 17010.0*x98;
x215 = x*(x180 + x212 + x214);
x216 = x104*x18;
x217 = -x216;
x218 = x110*x18;
x219 = -x108*x18;
x220 = x*(x217 + x218 + x219 + x90);
x221 = x174 + x210 + x219;
x222 = 31185.0*x87;
x223 = 17010.0*x89;
x224 = x178 + x222 - x223;
x225 = y*(x221 + x224);
x226 = x104*x15;
x227 = -x226;
x228 = y*(x211 + x218 + x227 + x99);
x229 = z*(x221 + x94);
x230 = z*(x102 + x179 + x212);
x231 = 155925.0*x109;
x232 = x106*x185;
x233 = -1351350.0*x232;
x234 = -x130 + x233;
x235 = x145*x15*x67;
x236 = -405405.0*x153 + x235;
x237 = x4*(311850.0*x114 + x231 + x234 + x236);
x238 = 155925.0*x114;
x239 = x145*x62;
x240 = x18*x239;
x241 = -405405.0*x147 + x240;
x242 = x4*(311850.0*x109 + x234 + x238 + x241);
x243 = x6*(x111 + x233 + x238 + x240);
x244 = -810810.0*x232;
x245 = x6*(187110.0*x114 + x187 + x236 + x244);
x246 = x8*(x105 - x113 + x115 + x231 + x233 + x235);
x247 = x8*(187110.0*x109 + x184 + x186 + x241 + x244);
x248 = 15120.0*x54;
x249 = -x119;
x250 = -x208;
x251 = 270.0*x10 + x205;
x252 = -x165 + x248*x62 + x249 + x250 + x251 - 5355.0*x37;
x253 = -x122;
x254 = -x204;
x255 = -x168 + x248*x67 + x251 + x253 + x254 - 5355.0*x42;
x256 = -x167;
x257 = x256 - x74;
x258 = x256 + x84;
x259 = -x169;
x260 = x259 - x77;
x261 = x259 - x79;
x262 = -x170;
x263 = x262 - x80;
x264 = x262 - x83;
x265 = -x181;
x266 = x265 - x91;
x267 = -x101 + x265;
x268 = -x182;
x269 = -x103 + x268;
x270 = x268 - x92;
x271 = -x183;
x272 = x271 - x95;
x273 = -x100 + x271;
x274 = -x188;
x275 = -x112 + x274;
x276 = x116 + x274;
x277 = x158*x18;
x278 = x146*x18;
x279 = x157*x18;
x280 = x131 + x217 + x277 + x278 - x279;
x281 = x15*x152;
x282 = x15*x163;
x283 = x140 + x227 + x277 + x281 - x282;
x284 = x18*x89;
x285 = -x15*x178 - x178*x18 + x250 + x254 + 11340.0*x284 + x65 + x72 + x75;
x286 = 3918915.0*x106;
x287 = -x195;
x288 = -x281;
x289 = -x277 - 12600.0*x36;
x290 = x121*x286 + x226 + x282 + x287 + x288 + x289 - 2338875.0*x96 + 439425.0*x98;
x291 = -x192;
x292 = -x278;
x293 = x117*x286 + x216 + x279 + x289 + x291 + x292 - 2338875.0*x87 + 439425.0*x89;
x294 = 16065.0*x54;
x295 = -360.0*x10 - x18*x223;
x296 = x122 + 20790.0*x15*x96 + x208 - x294*x67 + x295 + 1260.0*x37 + 6300.0*x42 + x73;
x297 = x119 + 20790.0*x18*x87 + x204 - x294*x62 + x295 + 6300.0*x37 + 1260.0*x42 + x76;
x298 = x15*x153;
x299 = 810810.0*x298;
x300 = -x299;
x301 = x239*x67;
x302 = x147*x18;
x303 = 810810.0*x302;
x304 = x301 - x303;
x305 = 374220.0*x175 + x214 + x224 + x300 + x304;
x306 = 2*x169 + x77 + x79;
x307 = 2*x170 + x80 + x83;
x308 = 2*x167 + x74 - x84;
x309 = 17010.0*x54;
x310 = 720.0*x10 - x15*x213 - x18*x222 + x249 + x253 + 34020.0*x284 + x309*x62 + x309*x67 - 7560.0*x37 - 7560.0*x42;
x311 = x100 + 2*x183 + x95;
x312 = -x215;
x313 = x142 + x312;
x314 = -x129;
x315 = -x220;
x316 = x314 + x315;
x317 = -x132;
x318 = -x225;
x319 = x317 + x318;
x320 = -x136;
x321 = -x228;
x322 = x320 + x321;
x323 = -x137;
x324 = -x229;
x325 = x323 + x324;
x326 = -x141;
x327 = -x230;
x328 = x326 + x327;
x329 = x101 + 2*x181 + x91;
x330 = x103 + 2*x182 + x92;
x331 = -x150;
x332 = -x242;
x333 = x331 + x332;
x334 = -x162;
x335 = -x237;
x336 = x334 + x335;
x337 = -x151;
x338 = -x243;
x339 = x337 + x338;
x340 = -x164;
x341 = -x245;
x342 = x340 + x341;
x343 = -x156;
x344 = -x246;
x345 = x343 + x344;
x346 = -x159;
x347 = -x247;
x348 = x346 + x347;
x349 = x112 - x116 + 2*x188;
x350 = x133 + x288;
x351 = -841995.0*x175 - x301 - 2520.0*x36;
x352 = x15*x154 - x222 + x303 + x350 + x351 + 31185.0*x89 - 187110.0*x96 + 59535.0*x98;
x353 = x126 + x292;
x354 = x148*x18 - x213 + x299 + x351 + x353 - 187110.0*x87 + 59535.0*x89 + 31185.0*x98;
x355 = 4054050.0*x106;
x356 = 4054050.0*x144;
x357 = x117*x18;
x358 = 1309770.0*x175 + 15120.0*x36;
x359 = -x117*x355 + x192 + x213 + x300 + x301 - 4864860.0*x302 + x356*x357 + x358 + 2525985.0*x87 - 498960.0*x89 - 45360.0*x98;
x360 = x121*x15;
x361 = -x121*x355 + x195 + x222 - 4864860.0*x298 + x304 + x356*x360 + x358 - 45360.0*x89 + 2525985.0*x96 - 498960.0*x98;
x362 = x312 + x315;
x363 = x318 + x321;
x364 = x324 + x327;
x365 = x332 + x335;
x366 = x338 + x341;
x367 = x344 + x347;
x368 = x62*x67;
x369 = x128 + x135 + 1683990.0*x175 + x278 + x281 - 3648645.0*x298 - 3648645.0*x302 + x356*x368 + 5040.0*x36 - 90720.0*x89 - 90720.0*x98;
x370 = 4189185.0*x106;
x371 = 6081075.0*x144;
x372 = -2993760.0*x175 - 20160.0*x36 - x368*x371;
x373 = x121*x370 + x287 + 8513505.0*x298 + 4459455.0*x302 + x353 - x360*x371 + x372 - 249480.0*x87 + 136080.0*x89 - 2744280.0*x96 + 589680.0*x98;
x374 = x117*x370 + x291 + 4459455.0*x298 + 8513505.0*x302 + x350 - x357*x371 + x372 - 2744280.0*x87 + 589680.0*x89 - 249480.0*x96 + 136080.0*x98;
x375 = 4324320.0*x106;
x376 = 8108100.0*x144;
x377 = -x117*x375 - x121*x375 + 12162150.0*x144*x368 + 5987520.0*x175 + x192 + x195 - 12972960.0*x298 - 12972960.0*x302 + x357*x376 + 40320.0*x36 + x360*x376 + 2993760.0*x87 - 725760.0*x89 + 2993760.0*x96 - 725760.0*x98;
x378 = x129 + x215 + 2*x220;
x379 = x132 + 2*x225 + x228;
x380 = x136 + x225 + 2*x228;
x381 = x137 + 2*x229 + x230;
x382 = x141 + x229 + 2*x230;
x383 = -x142 + 2*x215 + x220;
x384 = x150 + x237 + 2*x242;
x385 = x151 + 2*x243 + x245;
x386 = x156 + 2*x246 + x247;
x387 = x159 + x246 + 2*x247;
x388 = x162 + 2*x237 + x242;
x389 = x164 + x243 + 2*x245;
x390 = x142 - 3*x215 - 3*x220 + x314;
x391 = -3*x225 - 3*x228 + x317 + x320;
x392 = -3*x229 - 3*x230 + x323 + x326;
x393 = -3*x237 - 3*x242 + x331 + x334;
x394 = -3*x243 - 3*x245 + x337 + x340;
x395 = -3*x246 - 3*x247 + x343 + x346;
#pragma omp atomic
L[0] += -x*x1*M[0] - x1*y*M[1] - x1*z*M[2] + x100*M[77] + x101*M[70] + x103*M[71] + x112*M[87] - x116*M[105] + x120*M[55] + x123*M[76] + x129*M[83] - x13*M[13] + x132*M[84] + x136*M[111] + x137*M[85] + x141*M[112] - x142*M[104] + x150*M[120] + x151*M[121] + x156*M[156] + x159*M[123] + x161*M[22] + x162*M[147] + x164*M[148] + x167*M[37] + x169*M[40] + x17*M[3] + x170*M[41] + x171*M[14] + x172*M[17] + x173*M[18] + x181*M[61] + x182*M[62] + x183*M[66] + x188*M[94] + x189*M[27] + x190*M[28] + x191*M[32] + x194*M[119] + x196*M[155] + x197*M[47] + x20*M[6] + x200*M[31] + x201*M[24] + x203*M[33] + x207*M[65] + x209*M[58] + x215*M[93] + x220*M[86] + x225*M[89] + x228*M[98] + x229*M[90] + x230*M[99] + x237*M[134] + x24*M[9] + x242*M[125] + x243*M[126] + x245*M[135] + x246*M[141] + x247*M[130] + x252*M[60] + x255*M[78] + x257*M[39] + x258*M[46] + x26*M[10] + x260*M[51] + x261*M[42] + x263*M[43] + x264*M[52] + x266*M[63] + x267*M[72] + x269*M[73] + x270*M[64] + x272*M[68] + x273*M[79] + x275*M[96] + x276*M[107] + x280*M[122] + x283*M[140] + x285*M[67] + x29*M[15] + x290*M[157] + x293*M[124] + x296*M[80] + x297*M[69] + x30*M[11] + x305*M[129] + x306*M[53] + x307*M[54] + x308*M[48] + x31*M[16] + x310*M[82] + x311*M[81] + x313*M[106] + x316*M[88] + x319*M[91] + x322*M[113] + x325*M[92] + x328*M[114] + x329*M[74] - x33*M[12] + x330*M[75] + x333*M[127] + x336*M[149] + x339*M[128] + x342*M[150] + x345*M[158] + x348*M[132] + x349*M[109] + x352*M[142] + x354*M[131] + x359*M[133] + x361*M[159] + x362*M[95] + x363*M[100] + x364*M[101] + x365*M[136] + x366*M[137] + x367*M[143] + x369*M[144] + x373*M[161] + x374*M[146] + x377*M[163] + x378*M[97] + x379*M[102] + x380*M[115] + x381*M[103] + x382*M[116] + x383*M[108] + x384*M[138] + x385*M[139] + x386*M[160] + x387*M[145] + x388*M[151] + x389*M[152] + x390*M[110] + x391*M[117] + x392*M[118] + x393*M[153] + x394*M[154] + x395*M[162] + x40*M[20] + x41*M[21] + x45*M[30] + x47*M[23] + x49*M[25] + x5*M[4] + x52*M[26] + x58*M[38] - x61*M[45] + x66*M[19] + x69*M[29] + x7*M[5] + x70*M[8] + x74*M[34] + x77*M[49] + x79*M[35] + x80*M[36] + x83*M[50] - x84*M[44] + x9*M[7] + x91*M[56] + x92*M[57] + x95*M[59];
#pragma omp atomic
L[1] += x101*M[49] + x103*M[50] + x112*M[59] - x116*M[77] + x120*M[34] + x129*M[55] - x13*M[7] + x132*M[56] + x137*M[57] - x142*M[76] + x150*M[84] + x151*M[85] + x159*M[87] + x161*M[12] + x162*M[111] + x164*M[112] + x167*M[22] + x169*M[25] + x17*M[0] + x170*M[26] + x171*M[8] + x181*M[40] + x182*M[41] + x183*M[45] + x188*M[66] + x189*M[17] + x190*M[18] + x194*M[83] + x197*M[32] + x201*M[14] + x207*M[44] + x209*M[37] + x215*M[65] + x220*M[58] + x225*M[61] + x228*M[70] + x229*M[62] + x230*M[71] + x237*M[98] + x24*M[3] + x242*M[89] + x243*M[90] + x245*M[99] + x246*M[105] + x247*M[94] + x252*M[39] + x257*M[24] + x258*M[31] + x26*M[4] + x261*M[27] + x263*M[28] + x266*M[42] + x267*M[51] + x269*M[52] + x270*M[43] + x272*M[47] + x275*M[68] + x276*M[79] + x280*M[86] + x283*M[104] + x285*M[46] + x293*M[88] + x297*M[48] + x30*M[5] + x305*M[93] + x308*M[33] + x313*M[78] + x316*M[60] + x319*M[63] + x325*M[64] + x329*M[53] - x33*M[6] + x330*M[54] + x333*M[91] + x336*M[113] + x339*M[92] + x342*M[114] + x348*M[96] + x349*M[81] + x352*M[106] + x354*M[95] + x359*M[97] + x362*M[67] + x363*M[72] + x364*M[73] + x365*M[100] + x366*M[101] + x367*M[107] + x369*M[108] + x374*M[110] + x378*M[69] + x379*M[74] + x381*M[75] + x383*M[80] + x384*M[102] + x385*M[103] + x387*M[109] + x388*M[115] + x389*M[116] + x390*M[82] + x393*M[117] + x394*M[118] + x40*M[10] + x41*M[11] + x47*M[13] + x49*M[15] + x5*M[1] + x52*M[16] + x58*M[23] - x61*M[30] + x66*M[9] + x7*M[2] + x74*M[19] + x79*M[20] + x80*M[21] - x84*M[29] + x91*M[35] + x92*M[36] + x95*M[38];
#pragma omp atomic
L[2] += x100*M[50] + x101*M[44] + x103*M[45] + x112*M[57] - x116*M[71] + x123*M[49] - x13*M[5] + x132*M[55] + x136*M[76] + x141*M[77] - x142*M[70] + x150*M[83] + x156*M[112] + x159*M[85] + x161*M[10] + x162*M[104] + x164*M[105] + x167*M[20] + x169*M[22] + x170*M[23] + x172*M[8] + x181*M[37] + x182*M[38] + x183*M[41] + x188*M[62] + x189*M[14] + x191*M[18] + x196*M[111] + x197*M[28] + x20*M[1] + x200*M[17] + x207*M[40] + x209*M[35] + x215*M[61] + x220*M[56] + x225*M[58] + x228*M[65] + x229*M[59] + x230*M[66] + x237*M[93] + x242*M[86] + x243*M[87] + x245*M[94] + x246*M[99] + x247*M[90] + x255*M[51] + x258*M[27] + x26*M[3] + x260*M[31] + x261*M[24] + x264*M[32] + x266*M[39] + x267*M[46] + x269*M[47] + x272*M[43] + x273*M[52] + x275*M[64] + x276*M[73] + x280*M[84] + x283*M[98] + x285*M[42] + x29*M[6] + x290*M[113] + x296*M[53] + x305*M[89] + x306*M[33] + x31*M[7] + x311*M[54] + x313*M[72] + x319*M[60] + x322*M[78] + x328*M[79] + x329*M[48] - x33*M[4] + x333*M[88] + x336*M[106] + x342*M[107] + x345*M[114] + x348*M[92] + x349*M[75] + x352*M[100] + x354*M[91] + x361*M[115] + x362*M[63] + x363*M[67] + x364*M[68] + x365*M[95] + x366*M[96] + x367*M[101] + x369*M[102] + x373*M[117] + x379*M[69] + x380*M[80] + x382*M[81] + x383*M[74] + x384*M[97] + x386*M[116] + x387*M[103] + x388*M[108] + x389*M[109] + x391*M[82] + x393*M[110] + x395*M[118] + x40*M[9] + x45*M[16] + x47*M[11] + x49*M[12] + x5*M[0] + x52*M[13] + x58*M[21] - x61*M[26] + x69*M[15] + x77*M[29] + x79*M[19] + x83*M[30] - x84*M[25] + x9*M[2] + x91*M[34] + x95*M[36];
#pragma omp atomic
L[3] += x100*M[49] + x103*M[44] + x112*M[56] - x116*M[70] - x13*M[4] + x137*M[55] + x141*M[76] + x151*M[83] + x156*M[111] + x159*M[84] + x164*M[104] + x170*M[22] + x171*M[5] + x172*M[7] + x173*M[8] + x182*M[37] + x183*M[40] + x188*M[61] + x189*M[13] + x190*M[14] + x191*M[17] + x197*M[27] + x200*M[16] + x201*M[11] + x203*M[18] + x229*M[58] + x230*M[65] + x243*M[86] + x245*M[93] + x246*M[98] + x247*M[89] + x252*M[36] + x255*M[50] + x257*M[21] + x258*M[26] + x260*M[30] + x261*M[23] + x263*M[24] + x264*M[31] + x266*M[38] + x267*M[45] + x269*M[46] + x270*M[39] + x272*M[42] + x273*M[51] + x275*M[63] + x276*M[72] + x285*M[41] + x290*M[112] + x293*M[85] + x296*M[52] + x297*M[43] + x30*M[3] + x306*M[32] + x307*M[33] + x308*M[28] + x31*M[6] + x310*M[54] + x311*M[53] + x313*M[71] + x316*M[57] + x319*M[59] + x322*M[77] + x325*M[60] + x328*M[78] + x329*M[47] + x330*M[48] + x333*M[87] + x336*M[105] + x339*M[88] + x342*M[106] + x345*M[113] + x348*M[91] + x349*M[74] + x352*M[99] + x354*M[90] + x359*M[92] + x361*M[114] + x362*M[62] + x363*M[66] + x364*M[67] + x365*M[94] + x366*M[95] + x367*M[100] + x369*M[101] + x373*M[116] + x374*M[103] + x377*M[118] + x378*M[64] + x379*M[68] + x380*M[79] + x381*M[69] + x382*M[80] + x383*M[73] + x384*M[96] + x385*M[97] + x386*M[115] + x387*M[102] + x388*M[107] + x389*M[108] + x390*M[75] + x391*M[81] + x392*M[82] + x393*M[109] + x394*M[110] + x395*M[117] + x41*M[9] + x45*M[15] + x47*M[10] + x52*M[12] + x58*M[20] - x61*M[25] + x7*M[0] + x70*M[2] + x80*M[19] + x83*M[29] + x9*M[1] + x92*M[34] + x95*M[35];
#pragma omp atomic
L[4] += x112*M[38] + x120*M[19] + x129*M[34] + x132*M[35] + x137*M[36] + x150*M[56] + x151*M[57] + x159*M[59] + x161*M[6] + x167*M[12] + x169*M[15] + x170*M[16] + x181*M[25] + x182*M[26] + x183*M[30] + x188*M[45] + x194*M[55] + x201*M[8] + x207*M[29] + x209*M[22] + x215*M[44] + x220*M[37] + x225*M[40] + x228*M[49] + x229*M[41] + x230*M[50] + x237*M[70] + x24*M[0] + x242*M[61] + x243*M[62] + x245*M[71] + x246*M[77] + x247*M[66] + x252*M[24] + x257*M[14] + x26*M[1] + x261*M[17] + x263*M[18] + x266*M[27] + x270*M[28] + x272*M[32] + x275*M[47] + x280*M[58] + x283*M[76] + x285*M[31] + x293*M[60] + x297*M[33] + x30*M[2] + x305*M[65] + x316*M[39] + x319*M[42] + x325*M[43] + x333*M[63] + x339*M[64] + x348*M[68] + x352*M[78] + x354*M[67] + x359*M[69] + x362*M[46] + x363*M[51] + x364*M[52] + x365*M[72] + x366*M[73] + x367*M[79] + x369*M[80] + x374*M[82] + x378*M[48] + x379*M[53] + x381*M[54] + x384*M[74] + x385*M[75] + x387*M[81] + x40*M[4] + x41*M[5] + x47*M[7] + x58*M[13] + x66*M[3] + x74*M[9] + x79*M[10] + x80*M[11] + x91*M[20] + x92*M[21] + x95*M[23];
#pragma omp atomic
L[5] += x101*M[29] + x103*M[30] + x112*M[36] - x116*M[50] - x13*M[2] + x132*M[34] - x142*M[49] + x150*M[55] + x159*M[57] + x161*M[4] + x162*M[76] + x164*M[77] + x167*M[10] + x169*M[12] + x170*M[13] + x181*M[22] + x182*M[23] + x183*M[26] + x188*M[41] + x189*M[8] + x197*M[18] + x207*M[25] + x209*M[20] + x215*M[40] + x220*M[35] + x225*M[37] + x228*M[44] + x229*M[38] + x230*M[45] + x237*M[65] + x242*M[58] + x243*M[59] + x245*M[66] + x246*M[71] + x247*M[62] + x258*M[17] + x26*M[0] + x261*M[14] + x266*M[24] + x267*M[31] + x269*M[32] + x272*M[28] + x275*M[43] + x276*M[52] + x280*M[56] + x283*M[70] + x285*M[27] + x305*M[61] + x313*M[51] + x319*M[39] + x329*M[33] - x33*M[1] + x333*M[60] + x336*M[78] + x342*M[79] + x348*M[64] + x349*M[54] + x352*M[72] + x354*M[63] + x362*M[42] + x363*M[46] + x364*M[47] + x365*M[67] + x366*M[68] + x367*M[73] + x369*M[74] + x379*M[48] + x383*M[53] + x384*M[69] + x387*M[75] + x388*M[80] + x389*M[81] + x393*M[82] + x40*M[3] + x47*M[5] + x49*M[6] + x52*M[7] + x58*M[11] - x61*M[16] + x79*M[9] - x84*M[15] + x91*M[19] + x95*M[21];
#pragma omp atomic
L[6] += x103*M[29] + x112*M[35] - x116*M[49] - x13*M[1] + x137*M[34] + x151*M[55] + x159*M[56] + x164*M[76] + x170*M[12] + x171*M[2] + x182*M[22] + x183*M[25] + x188*M[40] + x189*M[7] + x190*M[8] + x197*M[17] + x201*M[5] + x229*M[37] + x230*M[44] + x243*M[58] + x245*M[65] + x246*M[70] + x247*M[61] + x252*M[21] + x257*M[11] + x258*M[16] + x261*M[13] + x263*M[14] + x266*M[23] + x267*M[30] + x269*M[31] + x270*M[24] + x272*M[27] + x275*M[42] + x276*M[51] + x285*M[26] + x293*M[57] + x297*M[28] + x30*M[0] + x308*M[18] + x313*M[50] + x316*M[36] + x319*M[38] + x325*M[39] + x329*M[32] + x330*M[33] + x333*M[59] + x336*M[77] + x339*M[60] + x342*M[78] + x348*M[63] + x349*M[53] + x352*M[71] + x354*M[62] + x359*M[64] + x362*M[41] + x363*M[45] + x364*M[46] + x365*M[66] + x366*M[67] + x367*M[72] + x369*M[73] + x374*M[75] + x378*M[43] + x379*M[47] + x381*M[48] + x383*M[52] + x384*M[68] + x385*M[69] + x387*M[74] + x388*M[79] + x389*M[80] + x390*M[54] + x393*M[81] + x394*M[82] + x41*M[3] + x47*M[4] + x52*M[6] + x58*M[10] - x61*M[15] + x80*M[9] + x92*M[19] + x95*M[20];
#pragma omp atomic
L[7] += x100*M[30] + x101*M[25] + x103*M[26] - x116*M[45] + x123*M[29] + x136*M[49] + x141*M[50] - x142*M[44] + x156*M[77] + x161*M[3] + x162*M[70] + x164*M[71] + x167*M[9] + x169*M[10] + x170*M[11] + x181*M[20] + x182*M[21] + x183*M[23] + x188*M[38] + x196*M[76] + x200*M[8] + x207*M[22] + x209*M[19] + x215*M[37] + x220*M[34] + x225*M[35] + x228*M[40] + x229*M[36] + x230*M[41] + x237*M[61] + x242*M[56] + x243*M[57] + x245*M[62] + x246*M[66] + x247*M[59] + x255*M[31] + x258*M[14] + x260*M[17] + x264*M[18] + x267*M[27] + x269*M[28] + x273*M[32] + x276*M[47] + x280*M[55] + x283*M[65] + x285*M[24] + x29*M[1] + x290*M[78] + x296*M[33] + x305*M[58] + x31*M[2] + x313*M[46] + x322*M[51] + x328*M[52] - x33*M[0] + x336*M[72] + x342*M[73] + x345*M[79] + x352*M[67] + x354*M[60] + x361*M[80] + x362*M[39] + x363*M[42] + x364*M[43] + x365*M[63] + x366*M[64] + x367*M[68] + x369*M[69] + x373*M[82] + x380*M[53] + x382*M[54] + x383*M[48] + x386*M[81] + x388*M[74] + x389*M[75] + x45*M[7] + x49*M[4] + x52*M[5] - x61*M[13] + x69*M[6] + x77*M[15] + x83*M[16] - x84*M[12];
#pragma omp atomic
L[8] += x100*M[29] + x103*M[25] + x112*M[34] - x116*M[44] - x13*M[0] + x141*M[49] + x156*M[76] + x159*M[55] + x164*M[70] + x170*M[10] + x172*M[2] + x182*M[20] + x183*M[22] + x188*M[37] + x189*M[5] + x191*M[8] + x197*M[14] + x200*M[7] + x229*M[35] + x230*M[40] + x243*M[56] + x245*M[61] + x246*M[65] + x247*M[58] + x255*M[30] + x258*M[13] + x260*M[16] + x261*M[11] + x264*M[17] + x266*M[21] + x267*M[26] + x269*M[27] + x272*M[24] + x273*M[31] + x275*M[39] + x276*M[46] + x285*M[23] + x290*M[77] + x296*M[32] + x306*M[18] + x31*M[1] + x311*M[33] + x313*M[45] + x319*M[36] + x322*M[50] + x328*M[51] + x329*M[28] + x333*M[57] + x336*M[71] + x342*M[72] + x345*M[78] + x348*M[60] + x349*M[48] + x352*M[66] + x354*M[59] + x361*M[79] + x362*M[38] + x363*M[41] + x364*M[42] + x365*M[62] + x366*M[63] + x367*M[67] + x369*M[68] + x373*M[81] + x379*M[43] + x380*M[52] + x382*M[53] + x383*M[47] + x384*M[64] + x386*M[80] + x387*M[69] + x388*M[73] + x389*M[74] + x391*M[54] + x393*M[75] + x395*M[82] + x45*M[6] + x47*M[3] + x52*M[4] + x58*M[9] - x61*M[12] + x83*M[15] + x95*M[19];
#pragma omp atomic
L[9] += x171*M[0] + x172*M[1] + x173*M[2] + x189*M[4] + x190*M[5] + x191*M[7] + x197*M[13] + x200*M[6] + x201*M[3] + x203*M[8] + x252*M[19] + x255*M[29] + x257*M[9] + x258*M[12] + x260*M[15] + x261*M[10] + x263*M[11] + x264*M[16] + x266*M[20] + x267*M[25] + x269*M[26] + x270*M[21] + x272*M[23] + x273*M[30] + x275*M[38] + x276*M[45] + x285*M[22] + x290*M[76] + x293*M[55] + x296*M[31] + x297*M[24] + x306*M[17] + x307*M[18] + x308*M[14] + x310*M[33] + x311*M[32] + x313*M[44] + x316*M[34] + x319*M[35] + x322*M[49] + x325*M[36] + x328*M[50] + x329*M[27] + x330*M[28] + x333*M[56] + x336*M[70] + x339*M[57] + x342*M[71] + x345*M[77] + x348*M[59] + x349*M[47] + x352*M[65] + x354*M[58] + x359*M[60] + x361*M[78] + x362*M[37] + x363*M[40] + x364*M[41] + x365*M[61] + x366*M[62] + x367*M[66] + x369*M[67] + x373*M[80] + x374*M[69] + x377*M[82] + x378*M[39] + x379*M[42] + x380*M[51] + x381*M[43] + x382*M[52] + x383*M[46] + x384*M[63] + x385*M[64] + x386*M[79] + x387*M[68] + x388*M[72] + x389*M[73] + x390*M[48] + x391*M[53] + x392*M[54] + x393*M[74] + x394*M[75] + x395*M[81];
#pragma omp atomic
L[10] += x112*M[23] + x120*M[9] + x129*M[19] + x132*M[20] + x137*M[21] + x150*M[35] + x151*M[36] + x159*M[38] + x167*M[6] + x181*M[15] + x182*M[16] + x188*M[30] + x194*M[34] + x209*M[12] + x215*M[29] + x220*M[22] + x225*M[25] + x229*M[26] + x237*M[49] + x242*M[40] + x243*M[41] + x245*M[50] + x247*M[45] + x252*M[14] + x257*M[8] + x266*M[17] + x270*M[18] + x275*M[32] + x280*M[37] + x293*M[39] + x305*M[44] + x316*M[24] + x319*M[27] + x325*M[28] + x333*M[42] + x339*M[43] + x348*M[47] + x354*M[46] + x359*M[48] + x362*M[31] + x365*M[51] + x366*M[52] + x378*M[33] + x384*M[53] + x385*M[54] + x40*M[1] + x41*M[2] + x58*M[7] + x66*M[0] + x74*M[3] + x79*M[4] + x80*M[5] + x91*M[10] + x92*M[11] + x95*M[13];
#pragma omp atomic
L[11] += x112*M[21] + x132*M[19] + x150*M[34] + x159*M[36] + x161*M[1] + x167*M[4] + x169*M[6] + x170*M[7] + x181*M[12] + x182*M[13] + x183*M[16] + x188*M[26] + x207*M[15] + x209*M[10] + x215*M[25] + x220*M[20] + x225*M[22] + x228*M[29] + x229*M[23] + x230*M[30] + x237*M[44] + x242*M[37] + x243*M[38] + x245*M[45] + x246*M[50] + x247*M[41] + x261*M[8] + x266*M[14] + x272*M[18] + x275*M[28] + x280*M[35] + x283*M[49] + x285*M[17] + x305*M[40] + x319*M[24] + x333*M[39] + x348*M[43] + x352*M[51] + x354*M[42] + x362*M[27] + x363*M[31] + x364*M[32] + x365*M[46] + x366*M[47] + x367*M[52] + x369*M[53] + x379*M[33] + x384*M[48] + x387*M[54] + x40*M[0] + x47*M[2] + x58*M[5] + x79*M[3] + x91*M[9] + x95*M[11];
#pragma omp atomic
L[12] += x112*M[20] + x137*M[19] + x151*M[34] + x159*M[35] + x170*M[6] + x182*M[12] + x183*M[15] + x188*M[25] + x201*M[2] + x229*M[22] + x230*M[29] + x243*M[37] + x245*M[44] + x246*M[49] + x247*M[40] + x252*M[11] + x257*M[5] + x261*M[7] + x263*M[8] + x266*M[13] + x270*M[14] + x272*M[17] + x275*M[27] + x285*M[16] + x293*M[36] + x297*M[18] + x316*M[21] + x319*M[23] + x325*M[24] + x333*M[38] + x339*M[39] + x348*M[42] + x352*M[50] + x354*M[41] + x359*M[43] + x362*M[26] + x363*M[30] + x364*M[31] + x365*M[45] + x366*M[46] + x367*M[51] + x369*M[52] + x374*M[54] + x378*M[28] + x379*M[32] + x381*M[33] + x384*M[47] + x385*M[48] + x387*M[53] + x41*M[0] + x47*M[1] + x58*M[4] + x80*M[3] + x92*M[9] + x95*M[10];
#pragma omp atomic
L[13] += x101*M[15] + x103*M[16] - x116*M[30] - x142*M[29] + x161*M[0] + x162*M[49] + x164*M[50] + x167*M[3] + x169*M[4] + x170*M[5] + x181*M[10] + x182*M[11] + x183*M[13] + x188*M[23] + x207*M[12] + x209*M[9] + x215*M[22] + x220*M[19] + x225*M[20] + x228*M[25] + x229*M[21] + x230*M[26] + x237*M[40] + x242*M[35] + x243*M[36] + x245*M[41] + x246*M[45] + x247*M[38] + x258*M[8] + x267*M[17] + x269*M[18] + x276*M[32] + x280*M[34] + x283*M[44] + x285*M[14] + x305*M[37] + x313*M[31] + x336*M[51] + x342*M[52] + x352*M[46] + x354*M[39] + x362*M[24] + x363*M[27] + x364*M[28] + x365*M[42] + x366*M[43] + x367*M[47] + x369*M[48] + x383*M[33] + x388*M[53] + x389*M[54] + x49*M[1] + x52*M[2] - x61*M[7] - x84*M[6];
#pragma omp atomic
L[14] += x103*M[15] + x112*M[19] - x116*M[29] + x159*M[34] + x164*M[49] + x170*M[4] + x182*M[10] + x183*M[12] + x188*M[22] + x189*M[2] + x197*M[8] + x229*M[20] + x230*M[25] + x243*M[35] + x245*M[40] + x246*M[44] + x247*M[37] + x258*M[7] + x261*M[5] + x266*M[11] + x267*M[16] + x269*M[17] + x272*M[14] + x275*M[24] + x276*M[31] + x285*M[13] + x313*M[30] + x319*M[21] + x329*M[18] + x333*M[36] + x336*M[50] + x342*M[51] + x348*M[39] + x349*M[33] + x352*M[45] + x354*M[38] + x362*M[23] + x363*M[26] + x364*M[27] + x365*M[41] + x366*M[42] + x367*M[46] + x369*M[47] + x379*M[28] + x383*M[32] + x384*M[43] + x387*M[48] + x388*M[52] + x389*M[53] + x393*M[54] + x47*M[0] + x52*M[1] + x58*M[3] - x61*M[6] + x95*M[9];
#pragma omp atomic
L[15] += x189*M[1] + x190*M[2] + x197*M[7] + x201*M[0] + x252*M[9] + x257*M[3] + x258*M[6] + x261*M[4] + x263*M[5] + x266*M[10] + x267*M[15] + x269*M[16] + x270*M[11] + x272*M[13] + x275*M[23] + x276*M[30] + x285*M[12] + x293*M[34] + x297*M[14] + x308*M[8] + x313*M[29] + x316*M[19] + x319*M[20] + x325*M[21] + x329*M[17] + x330*M[18] + x333*M[35] + x336*M[49] + x339*M[36] + x342*M[50] + x348*M[38] + x349*M[32] + x352*M[44] + x354*M[37] + x359*M[39] + x362*M[22] + x363*M[25] + x364*M[26] + x365*M[40] + x366*M[41] + x367*M[45] + x369*M[46] + x374*M[48] + x378*M[24] + x379*M[27] + x381*M[28] + x383*M[31] + x384*M[42] + x385*M[43] + x387*M[47] + x388*M[51] + x389*M[52] + x390*M[33] + x393*M[53] + x394*M[54];
#pragma omp atomic
L[16] += x100*M[16] + x101*M[12] + x103*M[13] - x116*M[26] + x123*M[15] + x136*M[29] + x141*M[30] - x142*M[25] + x156*M[50] + x162*M[44] + x164*M[45] + x169*M[3] + x181*M[9] + x183*M[11] + x188*M[21] + x196*M[49] + x207*M[10] + x215*M[20] + x225*M[19] + x228*M[22] + x230*M[23] + x237*M[37] + x242*M[34] + x245*M[38] + x246*M[41] + x247*M[36] + x255*M[17] + x260*M[8] + x267*M[14] + x273*M[18] + x276*M[28] + x283*M[40] + x290*M[51] + x305*M[35] + x313*M[27] + x322*M[31] + x328*M[32] + x336*M[46] + x342*M[47] + x345*M[52] + x352*M[42] + x361*M[53] + x363*M[24] + x365*M[39] + x367*M[43] + x380*M[33] + x386*M[54] + x388*M[48] + x45*M[2] + x49*M[0] - x61*M[5] + x69*M[1] + x77*M[6] + x83*M[7] - x84*M[4];
#pragma omp atomic
L[17] += x100*M[15] + x103*M[12] - x116*M[25] + x141*M[29] + x156*M[49] + x164*M[44] + x170*M[3] + x182*M[9] + x183*M[10] + x188*M[20] + x200*M[2] + x229*M[19] + x230*M[22] + x243*M[34] + x245*M[37] + x246*M[40] + x247*M[35] + x255*M[16] + x258*M[5] + x260*M[7] + x264*M[8] + x267*M[13] + x269*M[14] + x273*M[17] + x276*M[27] + x285*M[11] + x290*M[50] + x296*M[18] + x313*M[26] + x322*M[30] + x328*M[31] + x336*M[45] + x342*M[46] + x345*M[51] + x352*M[41] + x354*M[36] + x361*M[52] + x362*M[21] + x363*M[23] + x364*M[24] + x365*M[38] + x366*M[39] + x367*M[42] + x369*M[43] + x373*M[54] + x380*M[32] + x382*M[33] + x383*M[28] + x386*M[53] + x388*M[47] + x389*M[48] + x45*M[1] + x52*M[0] - x61*M[4] + x83*M[6];
#pragma omp atomic
L[18] += x189*M[0] + x191*M[2] + x197*M[5] + x200*M[1] + x255*M[15] + x258*M[4] + x260*M[6] + x261*M[3] + x264*M[7] + x266*M[9] + x267*M[12] + x269*M[13] + x272*M[11] + x273*M[16] + x275*M[21] + x276*M[26] + x285*M[10] + x290*M[49] + x296*M[17] + x306*M[8] + x311*M[18] + x313*M[25] + x319*M[19] + x322*M[29] + x328*M[30] + x329*M[14] + x333*M[34] + x336*M[44] + x342*M[45] + x345*M[50] + x348*M[36] + x349*M[28] + x352*M[40] + x354*M[35] + x361*M[51] + x362*M[20] + x363*M[22] + x364*M[23] + x365*M[37] + x366*M[38] + x367*M[41] + x369*M[42] + x373*M[53] + x379*M[24] + x380*M[31] + x382*M[32] + x383*M[27] + x384*M[39] + x386*M[52] + x387*M[43] + x388*M[46] + x389*M[47] + x391*M[33] + x393*M[48] + x395*M[54];
#pragma omp atomic
L[19] += x190*M[0] + x191*M[1] + x197*M[4] + x203*M[2] + x263*M[3] + x264*M[6] + x269*M[12] + x270*M[9] + x272*M[10] + x273*M[15] + x275*M[20] + x276*M[25] + x296*M[16] + x297*M[11] + x306*M[7] + x307*M[8] + x308*M[5] + x310*M[18] + x311*M[17] + x325*M[19] + x328*M[29] + x329*M[13] + x330*M[14] + x339*M[34] + x342*M[44] + x345*M[49] + x348*M[35] + x349*M[27] + x359*M[36] + x361*M[50] + x364*M[22] + x366*M[37] + x367*M[40] + x369*M[41] + x373*M[52] + x374*M[43] + x377*M[54] + x378*M[21] + x379*M[23] + x380*M[30] + x381*M[24] + x382*M[31] + x383*M[26] + x384*M[38] + x385*M[39] + x386*M[51] + x387*M[42] + x388*M[45] + x389*M[46] + x390*M[28] + x391*M[32] + x392*M[33] + x393*M[47] + x394*M[48] + x395*M[53];
#pragma omp atomic
L[20] += x112*M[13] + x120*M[3] + x129*M[9] + x132*M[10] + x137*M[11] + x150*M[20] + x151*M[21] + x159*M[23] + x194*M[19] + x209*M[6] + x220*M[12] + x225*M[15] + x229*M[16] + x242*M[25] + x243*M[26] + x247*M[30] + x252*M[8] + x280*M[22] + x293*M[24] + x305*M[29] + x316*M[14] + x319*M[17] + x325*M[18] + x333*M[27] + x339*M[28] + x348*M[32] + x354*M[31] + x359*M[33] + x74*M[0] + x79*M[1] + x80*M[2] + x91*M[4] + x92*M[5] + x95*M[7];
#pragma omp atomic
L[21] += x112*M[11] + x132*M[9] + x150*M[19] + x159*M[21] + x167*M[1] + x181*M[6] + x182*M[7] + x188*M[16] + x209*M[4] + x215*M[15] + x220*M[10] + x225*M[12] + x229*M[13] + x237*M[29] + x242*M[22] + x243*M[23] + x245*M[30] + x247*M[26] + x266*M[8] + x275*M[18] + x280*M[20] + x305*M[25] + x319*M[14] + x333*M[24] + x348*M[28] + x354*M[27] + x362*M[17] + x365*M[31] + x366*M[32] + x384*M[33] + x58*M[2] + x79*M[0] + x91*M[3] + x95*M[5];
#pragma omp atomic
L[22] += x112*M[10] + x137*M[9] + x151*M[19] + x159*M[20] + x182*M[6] + x188*M[15] + x229*M[12] + x243*M[22] + x245*M[29] + x247*M[25] + x252*M[5] + x257*M[2] + x266*M[7] + x270*M[8] + x275*M[17] + x293*M[21] + x316*M[11] + x319*M[13] + x325*M[14] + x333*M[23] + x339*M[24] + x348*M[27] + x354*M[26] + x359*M[28] + x362*M[16] + x365*M[30] + x366*M[31] + x378*M[18] + x384*M[32] + x385*M[33] + x58*M[1] + x80*M[0] + x92*M[3] + x95*M[4];
#pragma omp atomic
L[23] += x167*M[0] + x169*M[1] + x170*M[2] + x181*M[4] + x182*M[5] + x183*M[7] + x188*M[13] + x207*M[6] + x209*M[3] + x215*M[12] + x220*M[9] + x225*M[10] + x228*M[15] + x229*M[11] + x230*M[16] + x237*M[25] + x242*M[20] + x243*M[21] + x245*M[26] + x246*M[30] + x247*M[23] + x280*M[19] + x283*M[29] + x285*M[8] + x305*M[22] + x352*M[31] + x354*M[24] + x362*M[14] + x363*M[17] + x364*M[18] + x365*M[27] + x366*M[28] + x367*M[32] + x369*M[33];
#pragma omp atomic
L[24] += x112*M[9] + x159*M[19] + x170*M[1] + x182*M[4] + x183*M[6] + x188*M[12] + x229*M[10] + x230*M[15] + x243*M[20] + x245*M[25] + x246*M[29] + x247*M[22] + x261*M[2] + x266*M[5] + x272*M[8] + x275*M[14] + x285*M[7] + x319*M[11] + x333*M[21] + x348*M[24] + x352*M[30] + x354*M[23] + x362*M[13] + x363*M[16] + x364*M[17] + x365*M[26] + x366*M[27] + x367*M[31] + x369*M[32] + x379*M[18] + x384*M[28] + x387*M[33] + x58*M[0] + x95*M[3];
#pragma omp atomic
L[25] += x252*M[3] + x257*M[0] + x261*M[1] + x263*M[2] + x266*M[4] + x270*M[5] + x272*M[7] + x275*M[13] + x285*M[6] + x293*M[19] + x297*M[8] + x316*M[9] + x319*M[10] + x325*M[11] + x333*M[20] + x339*M[21] + x348*M[23] + x352*M[29] + x354*M[22] + x359*M[24] + x362*M[12] + x363*M[15] + x364*M[16] + x365*M[25] + x366*M[26] + x367*M[30] + x369*M[31] + x374*M[33] + x378*M[14] + x379*M[17] + x381*M[18] + x384*M[27] + x385*M[28] + x387*M[32];
#pragma omp atomic
L[26] += x101*M[6] + x103*M[7] - x116*M[16] - x142*M[15] + x162*M[29] + x164*M[30] + x169*M[0] + x181*M[3] + x183*M[5] + x188*M[11] + x207*M[4] + x215*M[10] + x225*M[9] + x228*M[12] + x230*M[13] + x237*M[22] + x242*M[19] + x245*M[23] + x246*M[26] + x247*M[21] + x267*M[8] + x276*M[18] + x283*M[25] + x305*M[20] + x313*M[17] + x336*M[31] + x342*M[32] + x352*M[27] + x363*M[14] + x365*M[24] + x367*M[28] + x388*M[33] - x61*M[2] - x84*M[1];
#pragma omp atomic
L[27] += x103*M[6] - x116*M[15] + x164*M[29] + x170*M[0] + x182*M[3] + x183*M[4] + x188*M[10] + x229*M[9] + x230*M[12] + x243*M[19] + x245*M[22] + x246*M[25] + x247*M[20] + x258*M[2] + x267*M[7] + x269*M[8] + x276*M[17] + x285*M[5] + x313*M[16] + x336*M[30] + x342*M[31] + x352*M[26] + x354*M[21] + x362*M[11] + x363*M[13] + x364*M[14] + x365*M[23] + x366*M[24] + x367*M[27] + x369*M[28] + x383*M[18] + x388*M[32] + x389*M[33] - x61*M[1];
#pragma omp atomic
L[28] += x197*M[2] + x258*M[1] + x261*M[0] + x266*M[3] + x267*M[6] + x269*M[7] + x272*M[5] + x275*M[11] + x276*M[16] + x285*M[4] + x313*M[15] + x319*M[9] + x329*M[8] + x333*M[19] + x336*M[29] + x342*M[30] + x348*M[21] + x349*M[18] + x352*M[25] + x354*M[20] + x362*M[10] + x363*M[12] + x364*M[13] + x365*M[22] + x366*M[23] + x367*M[26] + x369*M[27] + x379*M[14] + x383*M[17] + x384*M[24] + x387*M[28] + x388*M[31] + x389*M[32] + x393*M[33];
#pragma omp atomic
L[29] += x197*M[1] + x263*M[0] + x269*M[6] + x270*M[3] + x272*M[4] + x275*M[10] + x276*M[15] + x297*M[5] + x308*M[2] + x325*M[9] + x329*M[7] + x330*M[8] + x339*M[19] + x342*M[29] + x348*M[20] + x349*M[17] + x359*M[21] + x364*M[12] + x366*M[22] + x367*M[25] + x369*M[26] + x374*M[28] + x378*M[11] + x379*M[13] + x381*M[14] + x383*M[16] + x384*M[23] + x385*M[24] + x387*M[27] + x388*M[30] + x389*M[31] + x390*M[18] + x393*M[32] + x394*M[33];
#pragma omp atomic
L[30] += x100*M[7] + x101*M[4] + x103*M[5] - x116*M[13] + x123*M[6] + x136*M[15] + x141*M[16] - x142*M[12] + x156*M[30] + x162*M[25] + x164*M[26] + x196*M[29] + x207*M[3] + x215*M[9] + x228*M[10] + x230*M[11] + x237*M[20] + x245*M[21] + x246*M[23] + x255*M[8] + x283*M[22] + x290*M[31] + x305*M[19] + x313*M[14] + x322*M[17] + x328*M[18] + x336*M[27] + x342*M[28] + x345*M[32] + x352*M[24] + x361*M[33] + x77*M[1] + x83*M[2] - x84*M[0];
#pragma omp atomic
L[31] += x100*M[6] + x103*M[4] - x116*M[12] + x141*M[15] + x156*M[29] + x164*M[25] + x183*M[3] + x188*M[9] + x230*M[10] + x245*M[20] + x246*M[22] + x247*M[19] + x255*M[7] + x260*M[2] + x267*M[5] + x273*M[8] + x276*M[14] + x290*M[30] + x313*M[13] + x322*M[16] + x328*M[17] + x336*M[26] + x342*M[27] + x345*M[31] + x352*M[23] + x361*M[32] + x363*M[11] + x365*M[21] + x367*M[24] + x380*M[18] + x386*M[33] + x388*M[28] - x61*M[0] + x83*M[1];
#pragma omp atomic
L[32] += x255*M[6] + x258*M[0] + x260*M[1] + x264*M[2] + x267*M[4] + x269*M[5] + x273*M[7] + x276*M[13] + x285*M[3] + x290*M[29] + x296*M[8] + x313*M[12] + x322*M[15] + x328*M[16] + x336*M[25] + x342*M[26] + x345*M[30] + x352*M[22] + x354*M[19] + x361*M[31] + x362*M[9] + x363*M[10] + x364*M[11] + x365*M[20] + x366*M[21] + x367*M[23] + x369*M[24] + x373*M[33] + x380*M[17] + x382*M[18] + x383*M[14] + x386*M[32] + x388*M[27] + x389*M[28];
#pragma omp atomic
L[33] += x197*M[0] + x264*M[1] + x269*M[4] + x272*M[3] + x273*M[6] + x275*M[9] + x276*M[12] + x296*M[7] + x306*M[2] + x311*M[8] + x328*M[15] + x329*M[5] + x342*M[25] + x345*M[29] + x348*M[19] + x349*M[14] + x361*M[30] + x364*M[10] + x366*M[20] + x367*M[22] + x369*M[23] + x373*M[32] + x379*M[11] + x380*M[16] + x382*M[17] + x383*M[13] + x384*M[21] + x386*M[31] + x387*M[24] + x388*M[26] + x389*M[27] + x391*M[18] + x393*M[28] + x395*M[33];
#pragma omp atomic
L[34] += x296*M[6] + x297*M[3] + x306*M[1] + x307*M[2] + x308*M[0] + x310*M[8] + x311*M[7] + x329*M[4] + x330*M[5] + x349*M[13] + x359*M[19] + x361*M[29] + x369*M[22] + x373*M[31] + x374*M[24] + x377*M[33] + x378*M[9] + x379*M[10] + x380*M[15] + x381*M[11] + x382*M[16] + x383*M[12] + x384*M[20] + x385*M[21] + x386*M[30] + x387*M[23] + x388*M[25] + x389*M[26] + x390*M[14] + x391*M[17] + x392*M[18] + x393*M[27] + x394*M[28] + x395*M[32];
#pragma omp atomic
L[35] += x112*M[7] + x120*M[0] + x129*M[3] + x132*M[4] + x137*M[5] + x150*M[10] + x151*M[11] + x159*M[13] + x194*M[9] + x220*M[6] + x242*M[15] + x243*M[16] + x280*M[12] + x293*M[14] + x316*M[8] + x333*M[17] + x339*M[18] + x91*M[1] + x92*M[2];
#pragma omp atomic
L[36] += x112*M[5] + x132*M[3] + x150*M[9] + x159*M[11] + x209*M[1] + x220*M[4] + x225*M[6] + x229*M[7] + x242*M[12] + x243*M[13] + x247*M[16] + x280*M[10] + x305*M[15] + x319*M[8] + x333*M[14] + x348*M[18] + x354*M[17] + x91*M[0] + x95*M[2];
#pragma omp atomic
L[37] += x112*M[4] + x137*M[3] + x151*M[9] + x159*M[10] + x229*M[6] + x243*M[12] + x247*M[15] + x252*M[2] + x293*M[11] + x316*M[5] + x319*M[7] + x325*M[8] + x333*M[13] + x339*M[14] + x348*M[17] + x354*M[16] + x359*M[18] + x92*M[0] + x95*M[1];
#pragma omp atomic
L[38] += x181*M[1] + x182*M[2] + x188*M[7] + x209*M[0] + x215*M[6] + x220*M[3] + x225*M[4] + x229*M[5] + x237*M[15] + x242*M[10] + x243*M[11] + x245*M[16] + x247*M[13] + x280*M[9] + x305*M[12] + x354*M[14] + x362*M[8] + x365*M[17] + x366*M[18];
#pragma omp atomic
L[39] += x112*M[3] + x159*M[9] + x182*M[1] + x188*M[6] + x229*M[4] + x243*M[10] + x245*M[15] + x247*M[12] + x266*M[2] + x275*M[8] + x319*M[5] + x333*M[11] + x348*M[14] + x354*M[13] + x362*M[7] + x365*M[16] + x366*M[17] + x384*M[18] + x95*M[0];
#pragma omp atomic
L[40] += x252*M[0] + x266*M[1] + x270*M[2] + x275*M[7] + x293*M[9] + x316*M[3] + x319*M[4] + x325*M[5] + x333*M[10] + x339*M[11] + x348*M[13] + x354*M[12] + x359*M[14] + x362*M[6] + x365*M[15] + x366*M[16] + x378*M[8] + x384*M[17] + x385*M[18];
#pragma omp atomic
L[41] += x181*M[0] + x183*M[2] + x188*M[5] + x207*M[1] + x215*M[4] + x225*M[3] + x228*M[6] + x230*M[7] + x237*M[12] + x242*M[9] + x245*M[13] + x246*M[16] + x247*M[11] + x283*M[15] + x305*M[10] + x352*M[17] + x363*M[8] + x365*M[14] + x367*M[18];
#pragma omp atomic
L[42] += x182*M[0] + x183*M[1] + x188*M[4] + x229*M[3] + x230*M[6] + x243*M[9] + x245*M[12] + x246*M[15] + x247*M[10] + x285*M[2] + x352*M[16] + x354*M[11] + x362*M[5] + x363*M[7] + x364*M[8] + x365*M[13] + x366*M[14] + x367*M[17] + x369*M[18];
#pragma omp atomic
L[43] += x266*M[0] + x272*M[2] + x275*M[5] + x285*M[1] + x319*M[3] + x333*M[9] + x348*M[11] + x352*M[15] + x354*M[10] + x362*M[4] + x363*M[6] + x364*M[7] + x365*M[12] + x366*M[13] + x367*M[16] + x369*M[17] + x379*M[8] + x384*M[14] + x387*M[18];
#pragma omp atomic
L[44] += x270*M[0] + x272*M[1] + x275*M[4] + x297*M[2] + x325*M[3] + x339*M[9] + x348*M[10] + x359*M[11] + x364*M[6] + x366*M[12] + x367*M[15] + x369*M[16] + x374*M[18] + x378*M[5] + x379*M[7] + x381*M[8] + x384*M[13] + x385*M[14] + x387*M[17];
#pragma omp atomic
L[45] += x101*M[1] + x103*M[2] - x116*M[7] - x142*M[6] + x162*M[15] + x164*M[16] + x207*M[0] + x215*M[3] + x228*M[4] + x230*M[5] + x237*M[10] + x245*M[11] + x246*M[13] + x283*M[12] + x305*M[9] + x313*M[8] + x336*M[17] + x342*M[18] + x352*M[14];
#pragma omp atomic
L[46] += x103*M[1] - x116*M[6] + x164*M[15] + x183*M[0] + x188*M[3] + x230*M[4] + x245*M[10] + x246*M[12] + x247*M[9] + x267*M[2] + x276*M[8] + x313*M[7] + x336*M[16] + x342*M[17] + x352*M[13] + x363*M[5] + x365*M[11] + x367*M[14] + x388*M[18];
#pragma omp atomic
L[47] += x267*M[1] + x269*M[2] + x276*M[7] + x285*M[0] + x313*M[6] + x336*M[15] + x342*M[16] + x352*M[12] + x354*M[9] + x362*M[3] + x363*M[4] + x364*M[5] + x365*M[10] + x366*M[11] + x367*M[13] + x369*M[14] + x383*M[8] + x388*M[17] + x389*M[18];
#pragma omp atomic
L[48] += x269*M[1] + x272*M[0] + x275*M[3] + x276*M[6] + x329*M[2] + x342*M[15] + x348*M[9] + x349*M[8] + x364*M[4] + x366*M[10] + x367*M[12] + x369*M[13] + x379*M[5] + x383*M[7] + x384*M[11] + x387*M[14] + x388*M[16] + x389*M[17] + x393*M[18];
#pragma omp atomic
L[49] += x297*M[0] + x329*M[1] + x330*M[2] + x349*M[7] + x359*M[9] + x369*M[12] + x374*M[14] + x378*M[3] + x379*M[4] + x381*M[5] + x383*M[6] + x384*M[10] + x385*M[11] + x387*M[13] + x388*M[15] + x389*M[16] + x390*M[8] + x393*M[17] + x394*M[18];
#pragma omp atomic
L[50] += x100*M[2] + x101*M[0] - x116*M[5] + x123*M[1] + x136*M[6] + x141*M[7] - x142*M[4] + x156*M[16] + x162*M[12] + x164*M[13] + x196*M[15] + x228*M[3] + x237*M[9] + x246*M[11] + x283*M[10] + x290*M[17] + x322*M[8] + x336*M[14] + x345*M[18];
#pragma omp atomic
L[51] += x100*M[1] + x103*M[0] - x116*M[4] + x141*M[6] + x156*M[15] + x164*M[12] + x230*M[3] + x245*M[9] + x246*M[10] + x255*M[2] + x290*M[16] + x313*M[5] + x322*M[7] + x328*M[8] + x336*M[13] + x342*M[14] + x345*M[17] + x352*M[11] + x361*M[18];
#pragma omp atomic
L[52] += x255*M[1] + x267*M[0] + x273*M[2] + x276*M[5] + x290*M[15] + x313*M[4] + x322*M[6] + x328*M[7] + x336*M[12] + x342*M[13] + x345*M[16] + x352*M[10] + x361*M[17] + x363*M[3] + x365*M[9] + x367*M[11] + x380*M[8] + x386*M[18] + x388*M[14];
#pragma omp atomic
L[53] += x269*M[0] + x273*M[1] + x276*M[4] + x296*M[2] + x328*M[6] + x342*M[12] + x345*M[15] + x361*M[16] + x364*M[3] + x366*M[9] + x367*M[10] + x369*M[11] + x373*M[18] + x380*M[7] + x382*M[8] + x383*M[5] + x386*M[17] + x388*M[13] + x389*M[14];
#pragma omp atomic
L[54] += x296*M[1] + x311*M[2] + x329*M[0] + x349*M[5] + x361*M[15] + x369*M[10] + x373*M[17] + x379*M[3] + x380*M[6] + x382*M[7] + x383*M[4] + x384*M[9] + x386*M[16] + x387*M[11] + x388*M[12] + x389*M[13] + x391*M[8] + x393*M[14] + x395*M[18];
#pragma omp atomic
L[55] += x310*M[2] + x311*M[1] + x330*M[0] + x349*M[4] + x373*M[16] + x374*M[11] + x377*M[18] + x381*M[3] + x382*M[6] + x385*M[9] + x386*M[15] + x387*M[10] + x389*M[12] + x390*M[5] + x391*M[7] + x392*M[8] + x393*M[13] + x394*M[14] + x395*M[17];
#pragma omp atomic
L[56] += x129*M[0] + x132*M[1] + x137*M[2] + x150*M[4] + x151*M[5] + x159*M[7] + x194*M[3] + x280*M[6] + x293*M[8];
#pragma omp atomic
L[57] += x112*M[2] + x132*M[0] + x150*M[3] + x159*M[5] + x220*M[1] + x242*M[6] + x243*M[7] + x280*M[4] + x333*M[8];
#pragma omp atomic
L[58] += x112*M[1] + x137*M[0] + x151*M[3] + x159*M[4] + x243*M[6] + x293*M[5] + x316*M[2] + x333*M[7] + x339*M[8];
#pragma omp atomic
L[59] += x220*M[0] + x225*M[1] + x229*M[2] + x242*M[4] + x243*M[5] + x247*M[7] + x280*M[3] + x305*M[6] + x354*M[8];
#pragma omp atomic
L[60] += x112*M[0] + x159*M[3] + x229*M[1] + x243*M[4] + x247*M[6] + x319*M[2] + x333*M[5] + x348*M[8] + x354*M[7];
#pragma omp atomic
L[61] += x293*M[3] + x316*M[0] + x319*M[1] + x325*M[2] + x333*M[4] + x339*M[5] + x348*M[7] + x354*M[6] + x359*M[8];
#pragma omp atomic
L[62] += x188*M[2] + x215*M[1] + x225*M[0] + x237*M[6] + x242*M[3] + x245*M[7] + x247*M[5] + x305*M[4] + x365*M[8];
#pragma omp atomic
L[63] += x188*M[1] + x229*M[0] + x243*M[3] + x245*M[6] + x247*M[4] + x354*M[5] + x362*M[2] + x365*M[7] + x366*M[8];
#pragma omp atomic
L[64] += x275*M[2] + x319*M[0] + x333*M[3] + x348*M[5] + x354*M[4] + x362*M[1] + x365*M[6] + x366*M[7] + x384*M[8];
#pragma omp atomic
L[65] += x275*M[1] + x325*M[0] + x339*M[3] + x348*M[4] + x359*M[5] + x366*M[6] + x378*M[2] + x384*M[7] + x385*M[8];
#pragma omp atomic
L[66] += x215*M[0] + x228*M[1] + x230*M[2] + x237*M[4] + x245*M[5] + x246*M[7] + x283*M[6] + x305*M[3] + x352*M[8];
#pragma omp atomic
L[67] += x188*M[0] + x230*M[1] + x245*M[4] + x246*M[6] + x247*M[3] + x352*M[7] + x363*M[2] + x365*M[5] + x367*M[8];
#pragma omp atomic
L[68] += x352*M[6] + x354*M[3] + x362*M[0] + x363*M[1] + x364*M[2] + x365*M[4] + x366*M[5] + x367*M[7] + x369*M[8];
#pragma omp atomic
L[69] += x275*M[0] + x348*M[3] + x364*M[1] + x366*M[4] + x367*M[6] + x369*M[7] + x379*M[2] + x384*M[5] + x387*M[8];
#pragma omp atomic
L[70] += x359*M[3] + x369*M[6] + x374*M[8] + x378*M[0] + x379*M[1] + x381*M[2] + x384*M[4] + x385*M[5] + x387*M[7];
#pragma omp atomic
L[71] += -x116*M[2] - x142*M[1] + x162*M[6] + x164*M[7] + x228*M[0] + x237*M[3] + x246*M[5] + x283*M[4] + x336*M[8];
#pragma omp atomic
L[72] += -x116*M[1] + x164*M[6] + x230*M[0] + x245*M[3] + x246*M[4] + x313*M[2] + x336*M[7] + x342*M[8] + x352*M[5];
#pragma omp atomic
L[73] += x276*M[2] + x313*M[1] + x336*M[6] + x342*M[7] + x352*M[4] + x363*M[0] + x365*M[3] + x367*M[5] + x388*M[8];
#pragma omp atomic
L[74] += x276*M[1] + x342*M[6] + x364*M[0] + x366*M[3] + x367*M[4] + x369*M[5] + x383*M[2] + x388*M[7] + x389*M[8];
#pragma omp atomic
L[75] += x349*M[2] + x369*M[4] + x379*M[0] + x383*M[1] + x384*M[3] + x387*M[5] + x388*M[6] + x389*M[7] + x393*M[8];
#pragma omp atomic
L[76] += x349*M[1] + x374*M[5] + x381*M[0] + x385*M[3] + x387*M[4] + x389*M[6] + x390*M[2] + x393*M[7] + x394*M[8];
#pragma omp atomic
L[77] += x136*M[1] + x141*M[2] - x142*M[0] + x156*M[7] + x162*M[4] + x164*M[5] + x196*M[6] + x283*M[3] + x290*M[8];
#pragma omp atomic
L[78] += -x116*M[0] + x141*M[1] + x156*M[6] + x164*M[4] + x246*M[3] + x290*M[7] + x322*M[2] + x336*M[5] + x345*M[8];
#pragma omp atomic
L[79] += x290*M[6] + x313*M[0] + x322*M[1] + x328*M[2] + x336*M[4] + x342*M[5] + x345*M[7] + x352*M[3] + x361*M[8];
#pragma omp atomic
L[80] += x276*M[0] + x328*M[1] + x342*M[4] + x345*M[6] + x361*M[7] + x367*M[3] + x380*M[2] + x386*M[8] + x388*M[5];
#pragma omp atomic
L[81] += x361*M[6] + x369*M[3] + x373*M[8] + x380*M[1] + x382*M[2] + x383*M[0] + x386*M[7] + x388*M[4] + x389*M[5];
#pragma omp atomic
L[82] += x349*M[0] + x373*M[7] + x382*M[1] + x386*M[6] + x387*M[3] + x389*M[4] + x391*M[2] + x393*M[5] + x395*M[8];
#pragma omp atomic
L[83] += x373*M[6] + x374*M[3] + x377*M[8] + x390*M[0] + x391*M[1] + x392*M[2] + x393*M[4] + x394*M[5] + x395*M[7];
#pragma omp atomic
L[84] += x150*M[1] + x151*M[2] + x194*M[0];
#pragma omp atomic
L[85] += x150*M[0] + x159*M[2] + x280*M[1];
#pragma omp atomic
L[86] += x151*M[0] + x159*M[1] + x293*M[2];
#pragma omp atomic
L[87] += x242*M[1] + x243*M[2] + x280*M[0];
#pragma omp atomic
L[88] += x159*M[0] + x243*M[1] + x333*M[2];
#pragma omp atomic
L[89] += x293*M[0] + x333*M[1] + x339*M[2];
#pragma omp atomic
L[90] += x242*M[0] + x247*M[2] + x305*M[1];
#pragma omp atomic
L[91] += x243*M[0] + x247*M[1] + x354*M[2];
#pragma omp atomic
L[92] += x333*M[0] + x348*M[2] + x354*M[1];
#pragma omp atomic
L[93] += x339*M[0] + x348*M[1] + x359*M[2];
#pragma omp atomic
L[94] += x237*M[1] + x245*M[2] + x305*M[0];
#pragma omp atomic
L[95] += x245*M[1] + x247*M[0] + x365*M[2];
#pragma omp atomic
L[96] += x354*M[0] + x365*M[1] + x366*M[2];
#pragma omp atomic
L[97] += x348*M[0] + x366*M[1] + x384*M[2];
#pragma omp atomic
L[98] += x359*M[0] + x384*M[1] + x385*M[2];
#pragma omp atomic
L[99] += x237*M[0] + x246*M[2] + x283*M[1];
#pragma omp atomic
L[100] += x245*M[0] + x246*M[1] + x352*M[2];
#pragma omp atomic
L[101] += x352*M[1] + x365*M[0] + x367*M[2];
#pragma omp atomic
L[102] += x366*M[0] + x367*M[1] + x369*M[2];
#pragma omp atomic
L[103] += x369*M[1] + x384*M[0] + x387*M[2];
#pragma omp atomic
L[104] += x374*M[2] + x385*M[0] + x387*M[1];
#pragma omp atomic
L[105] += x162*M[1] + x164*M[2] + x283*M[0];
#pragma omp atomic
L[106] += x164*M[1] + x246*M[0] + x336*M[2];
#pragma omp atomic
L[107] += x336*M[1] + x342*M[2] + x352*M[0];
#pragma omp atomic
L[108] += x342*M[1] + x367*M[0] + x388*M[2];
#pragma omp atomic
L[109] += x369*M[0] + x388*M[1] + x389*M[2];
#pragma omp atomic
L[110] += x387*M[0] + x389*M[1] + x393*M[2];
#pragma omp atomic
L[111] += x374*M[0] + x393*M[1] + x394*M[2];
#pragma omp atomic
L[112] += x156*M[2] + x162*M[0] + x196*M[1];
#pragma omp atomic
L[113] += x156*M[1] + x164*M[0] + x290*M[2];
#pragma omp atomic
L[114] += x290*M[1] + x336*M[0] + x345*M[2];
#pragma omp atomic
L[115] += x342*M[0] + x345*M[1] + x361*M[2];
#pragma omp atomic
L[116] += x361*M[1] + x386*M[2] + x388*M[0];
#pragma omp atomic
L[117] += x373*M[2] + x386*M[1] + x389*M[0];
#pragma omp atomic
L[118] += x373*M[1] + x393*M[0] + x395*M[2];
#pragma omp atomic
L[119] += x377*M[2] + x394*M[0] + x395*M[1];

}

void L2L_8(double x, double y, double z, double * L, double * Ls) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
x3 = z*L[14];
x4 = x3*y;
x5 = (x*x);
x6 = (1.0/2.0)*x5;
x7 = (x*x*x);
x8 = (1.0/6.0)*x7;
x9 = (x*x*x*x);
x10 = (1.0/24.0)*x9;
x11 = pow(x, 5);
x12 = (1.0/120.0)*x11;
x13 = (1.0/720.0)*pow(x, 6);
x14 = (y*y);
x15 = (1.0/2.0)*x14;
x16 = (y*y*y);
x17 = (1.0/6.0)*x16;
x18 = (y*y*y*y);
x19 = (1.0/24.0)*x18;
x20 = pow(y, 5);
x21 = (1.0/120.0)*x20;
x22 = (1.0/720.0)*pow(y, 6);
x23 = (z*z);
x24 = (1.0/2.0)*x23;
x25 = (z*z*z);
x26 = (1.0/6.0)*x25;
x27 = (z*z*z*z);
x28 = (1.0/24.0)*x27;
x29 = pow(z, 5);
x30 = (1.0/120.0)*x29;
x31 = (1.0/720.0)*pow(z, 6);
x32 = x*L[13];
x33 = x*L[26];
x34 = x*L[45];
x35 = x*L[71];
x36 = x*L[105];
x37 = x*L[15];
x38 = x*L[29];
x39 = x*L[49];
x40 = x*L[76];
x41 = x*L[111];
x42 = y*L[11];
x43 = z*L[12];
x44 = y*L[21];
x45 = z*L[22];
x46 = y*L[36];
x47 = z*L[37];
x48 = y*L[57];
x49 = z*L[58];
x50 = y*L[85];
x51 = z*L[86];
x52 = y*L[18];
x53 = y*L[33];
x54 = y*L[54];
x55 = y*L[82];
x56 = y*L[118];
x57 = z*L[17];
x58 = z*L[31];
x59 = z*L[51];
x60 = z*L[78];
x61 = z*L[113];
x62 = y*L[28];
x63 = x*x62;
x64 = y*L[48];
x65 = x*x64;
x66 = y*L[75];
x67 = x*x66;
x68 = y*L[110];
x69 = x*x68;
x70 = z*L[27];
x71 = x*x70;
x72 = z*L[46];
x73 = x*x72;
x74 = z*L[72];
x75 = x*x74;
x76 = z*L[106];
x77 = x*x76;
x78 = z*L[24];
x79 = x78*y;
x80 = z*L[39];
x81 = x80*y;
x82 = z*L[60];
x83 = x82*y;
x84 = z*L[88];
x85 = x84*y;
x86 = (1.0/4.0)*x5;
x87 = x14*x86;
x88 = (1.0/12.0)*x5;
x89 = x16*x88;
x90 = (1.0/48.0)*x5;
x91 = x18*x90;
x92 = (1.0/240.0)*x5;
x93 = x23*x86;
x94 = x25*x88;
x95 = x27*x90;
x96 = (1.0/12.0)*x7;
x97 = x14*x96;
x98 = (1.0/36.0)*x7;
x99 = x16*x98;
x100 = (1.0/144.0)*x7;
x101 = x23*x96;
x102 = x25*x98;
x103 = (1.0/48.0)*x9;
x104 = x103*x14;
x105 = (1.0/144.0)*x9;
x106 = x103*x23;
x107 = (1.0/240.0)*x11;
x108 = x14*x23;
x109 = (1.0/4.0)*x108;
x110 = x14*x25;
x111 = (1.0/12.0)*x110;
x112 = (1.0/48.0)*x14*x27;
x113 = x16*x23;
x114 = (1.0/12.0)*x113;
x115 = (1.0/36.0)*x16*x25;
x116 = (1.0/48.0)*x18*x23;
x117 = x*L[47];
x118 = x*L[74];
x119 = x*L[109];
x120 = x*L[73];
x121 = x*L[108];
x122 = x*L[107];
x123 = y*L[43];
x124 = y*L[69];
x125 = y*L[103];
x126 = z*L[42];
x127 = z*L[67];
x128 = z*L[100];
x129 = y*L[64];
x130 = y*L[97];
x131 = z*L[63];
x132 = z*L[95];
x133 = y*L[92];
x134 = z*L[91];
x135 = (1.0/8.0)*x108*x5;
x136 = (1.0/24.0)*x5;
x137 = x*L[23];
x138 = x*L[41];
x139 = x*L[66];
x140 = x*L[99];
x141 = x*L[25];
x142 = x*L[44];
x143 = x*L[70];
x144 = x*L[104];
x145 = x*x123;
x146 = x*x124;
x147 = x*x125;
x148 = x*x126;
x149 = x*x127;
x150 = x*x128;
x151 = x*L[68];
x152 = x*L[102];
x153 = x*L[101];
x154 = y*L[13];
x155 = x70*y;
x156 = x*L[28];
x157 = x*L[48];
x158 = x*L[75];
x159 = x*L[110];
x160 = y*L[23];
x161 = y*L[38];
x162 = y*L[59];
x163 = y*L[87];
x164 = y*L[32];
x165 = y*L[53];
x166 = y*L[81];
x167 = y*L[117];
x168 = y*L[47];
x169 = x*x168;
x170 = y*L[74];
x171 = x*x170;
x172 = y*L[109];
x173 = x*x172;
x174 = x126*y;
x175 = x131*y;
x176 = x134*y;
x177 = y*L[68];
x178 = y*L[102];
x179 = y*L[96];
x180 = y*L[14];
x181 = z*L[15];
x182 = z*L[18];
x183 = z*L[28];
x184 = x183*y;
x185 = x*L[27];
x186 = x*L[46];
x187 = x*L[72];
x188 = x*L[106];
x189 = y*L[24];
x190 = z*L[25];
x191 = y*L[39];
x192 = z*L[40];
x193 = y*L[60];
x194 = z*L[61];
x195 = y*L[88];
x196 = z*L[89];
x197 = z*L[32];
x198 = z*L[52];
x199 = z*L[79];
x200 = z*L[114];
x201 = z*L[47];
x202 = x*x201;
x203 = z*L[73];
x204 = x*x203;
x205 = z*L[107];
x206 = x*x205;
x207 = z*L[43];
x208 = x207*y;
x209 = z*L[64];
x210 = x209*y;
x211 = z*L[92];
x212 = x211*y;
x213 = z*L[68];
x214 = z*L[101];
x215 = z*L[96];
x216 = x*L[38];
x217 = x*L[62];
x218 = x*L[94];
x219 = x*L[40];
x220 = x*L[65];
x221 = x*L[98];
x222 = x*x129;
x223 = x*x130;
x224 = x*x131;
x225 = x*x132;
x226 = x*L[96];
x227 = x*L[43];
x228 = x*L[69];
x229 = x*L[103];
x230 = x*x177;
x231 = x*x178;
x232 = x*L[42];
x233 = x*L[67];
x234 = x*L[100];
x235 = x*x213;
x236 = x*x214;
x237 = y*L[26];
x238 = x72*y;
x239 = y*L[41];
x240 = y*L[62];
x241 = y*L[90];
x242 = y*L[52];
x243 = y*L[80];
x244 = y*L[116];
x245 = y*L[73];
x246 = x*x245;
x247 = y*L[108];
x248 = x*x247;
x249 = x127*y;
x250 = x132*y;
x251 = y*L[101];
x252 = y*L[27];
x253 = x201*y;
x254 = y*L[42];
x255 = y*L[63];
x256 = y*L[91];
x257 = x213*y;
x258 = x215*y;
x259 = z*L[29];
x260 = z*L[33];
x261 = z*L[48];
x262 = x261*y;
x263 = z*L[44];
x264 = z*L[65];
x265 = z*L[93];
x266 = z*L[53];
x267 = z*L[80];
x268 = z*L[115];
x269 = z*L[74];
x270 = x*x269;
x271 = z*L[108];
x272 = x*x271;
x273 = z*L[69];
x274 = x273*y;
x275 = z*L[97];
x276 = x275*y;
x277 = z*L[102];
x278 = x*L[59];
x279 = x*L[90];
x280 = x*L[61];
x281 = x*L[93];
x282 = x*x133;
x283 = x*x134;
x284 = x*L[64];
x285 = x*L[97];
x286 = x*x179;
x287 = x*L[63];
x288 = x*L[95];
x289 = x*x215;
x290 = x*x251;
x291 = x*x277;
x292 = y*L[45];
x293 = x74*y;
x294 = y*L[66];
x295 = y*L[94];
x296 = y*L[79];
x297 = y*L[115];
x298 = y*L[107];
x299 = x*x298;
x300 = x128*y;
x301 = y*L[46];
x302 = x203*y;
x303 = y*L[67];
x304 = y*L[95];
x305 = x214*y;
x306 = x269*y;
x307 = x277*y;
x308 = z*L[49];
x309 = z*L[54];
x310 = z*L[75];
x311 = x310*y;
x312 = z*L[70];
x313 = z*L[98];
x314 = z*L[81];
x315 = z*L[116];
x316 = z*L[109];
x317 = x*x316;
x318 = z*L[103];
x319 = x318*y;
x320 = x*L[87];
x321 = x*L[89];
x322 = x*L[92];
x323 = x*L[91];
x324 = y*L[71];
x325 = x76*y;
x326 = y*L[99];
x327 = y*L[114];
x328 = y*L[72];
x329 = x205*y;
x330 = y*L[100];
x331 = x271*y;
x332 = x316*y;
x333 = z*L[76];
x334 = z*L[82];
x335 = z*L[110];
x336 = x335*y;
x337 = z*L[104];
x338 = z*L[117];
x339 = y*L[105];
x340 = y*L[106];
x341 = z*L[111];
x342 = z*L[118];
#pragma omp atomic
Ls[0] += (1.0/5040.0)*pow(x, 7)*L[84] + x*x0 + x*x1 + x*x4 + x*L[1] + x10*x46 + x10*x47 + x10*x83 + x10*L[20] + x100*x18*L[94] + x100*x27*L[98] + x101*x129 + x101*L[40] + x102*x130 + x102*L[65] + x104*x134 + x104*L[59] + x105*x16*L[90] + x105*x25*L[93] + x106*x133 + x106*L[61] + x107*x14*L[87] + x107*x23*L[89] + (1.0/24.0)*x108*x7*L[96] + x109*x117 + x109*L[32] + x110*x136*L[102] + x111*x118 + x111*L[53] + x112*x119 + x112*L[81] + x113*x136*L[101] + x114*x120 + x114*L[52] + x115*x121 + x115*L[80] + x116*x122 + x116*L[79] + x12*x48 + x12*x49 + x12*x85 + x12*L[35] + x123*x93 + x124*x94 + x125*x95 + x126*x87 + x127*x89 + x128*x91 + x13*x50 + x13*x51 + x13*L[56] + x131*x97 + x132*x99 + x135*L[68] + (1.0/240.0)*x14*x29*L[117] + x15*x32 + x15*x57 + x15*x71 + x15*L[7] + (1.0/144.0)*x16*x27*L[116] + x17*x33 + x17*x58 + x17*x73 + x17*L[16] + (1.0/144.0)*x18*x25*L[115] + x19*x34 + x19*x59 + x19*x75 + x19*L[30] + x2*y + (1.0/240.0)*x20*x23*L[114] + x20*x92*L[99] + x21*x35 + x21*x60 + x21*x77 + x21*L[50] + x22*x36 + x22*x61 + x22*L[77] + x24*x37 + x24*x52 + x24*x63 + x24*L[9] + x26*x38 + x26*x53 + x26*x65 + x26*L[19] + x28*x39 + x28*x54 + x28*x67 + x28*L[34] + x29*x92*L[104] + x30*x40 + x30*x55 + x30*x69 + x30*L[55] + x31*x41 + x31*x56 + x31*L[83] + x42*x6 + x43*x6 + x44*x8 + x45*x8 + x6*x79 + x6*L[4] + x8*x81 + x8*L[10] + x87*L[23] + x89*L[41] + x91*L[66] + x93*L[25] + x94*L[44] + x95*L[70] + x97*L[38] + x99*L[62] + (1.0/5040.0)*pow(y, 7)*L[112] + y*L[2] + (1.0/5040.0)*pow(z, 7)*L[119] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*x42 + x*x43 + x*x79 + x*L[4] + x0 + x1 + x10*x48 + x10*x49 + x10*x85 + x10*L[35] + x101*x133 + x101*L[61] + x102*L[93] + x104*L[87] + x106*L[89] + x109*x151 + x109*L[47] + x111*x152 + x111*L[74] + x112*L[109] + x114*x153 + x114*L[73] + x115*L[108] + x116*L[107] + x12*x50 + x12*x51 + x12*L[56] + x129*x93 + x13*L[84] + x130*x94 + x131*x87 + x132*x89 + x134*x97 + x135*L[96] + x137*x15 + x138*x17 + x139*x19 + x140*x21 + x141*x24 + x142*x26 + x143*x28 + x144*x30 + x145*x24 + x146*x26 + x147*x28 + x148*x15 + x149*x17 + x15*x70 + x15*L[13] + x150*x19 + x17*x72 + x17*L[26] + x19*x74 + x19*L[45] + x21*x76 + x21*L[71] + x22*L[105] + x24*x62 + x24*L[15] + x26*x64 + x26*L[29] + x28*x66 + x28*L[49] + x30*x68 + x30*L[76] + x31*L[111] + x4 + x44*x6 + x45*x6 + x46*x8 + x47*x8 + x6*x81 + x6*L[10] + x8*x83 + x8*L[20] + x87*L[38] + x89*L[62] + x91*L[94] + x93*L[40] + x94*L[65] + x95*L[98] + x97*L[59] + x99*L[90] + L[1];
#pragma omp atomic
Ls[2] += x*x154 + x*x155 + x*x3 + x*L[5] + x10*x162 + x10*x176 + x10*x82 + x10*L[36] + x101*x179 + x101*L[64] + x102*L[97] + x104*L[90] + x106*L[92] + x109*x120 + x109*L[52] + x111*x121 + x111*L[80] + x112*L[116] + x114*x122 + x114*L[79] + x115*L[115] + x116*L[114] + x12*x163 + x12*x84 + x12*L[57] + x127*x87 + x128*x89 + x13*L[85] + x132*x97 + x135*L[101] + x15*x33 + x15*x58 + x15*x73 + x15*L[16] + x156*x24 + x157*x26 + x158*x28 + x159*x30 + x160*x6 + x161*x8 + x164*x24 + x165*x26 + x166*x28 + x167*x30 + x169*x24 + x17*x34 + x17*x59 + x17*x75 + x17*L[30] + x171*x26 + x173*x28 + x174*x6 + x175*x8 + x177*x93 + x178*x94 + x19*x35 + x19*x60 + x19*x77 + x19*L[50] + x2 + x21*x36 + x21*x61 + x21*L[77] + x22*L[112] + x24*L[18] + x26*L[33] + x28*L[54] + x30*L[82] + x31*L[118] + x57*y + x6*x78 + x6*L[11] + x8*x80 + x8*L[21] + x87*L[41] + x89*L[66] + x91*L[99] + x93*L[43] + x94*L[69] + x95*L[103] + x97*L[62] + x99*L[94] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*x180 + x*x181 + x*x184 + x*L[6] + x10*x193 + x10*x194 + x10*x212 + x10*L[37] + x101*x130 + x101*L[65] + x102*L[98] + x104*L[91] + x106*L[93] + x109*x118 + x109*L[53] + x111*x119 + x111*L[81] + x112*L[117] + x114*x121 + x114*L[80] + x115*L[116] + x116*L[115] + x12*x195 + x12*x196 + x12*L[58] + x124*x93 + x125*x94 + x13*L[86] + x135*L[102] + x15*x185 + x15*x197 + x15*x202 + x15*L[17] + x17*x186 + x17*x198 + x17*x204 + x17*L[31] + x182*y + x187*x19 + x188*x21 + x189*x6 + x19*x199 + x19*x206 + x19*L[51] + x190*x6 + x191*x8 + x192*x8 + x200*x21 + x208*x6 + x21*L[78] + x210*x8 + x213*x87 + x214*x89 + x215*x97 + x22*L[113] + x24*x38 + x24*x53 + x24*x65 + x24*L[19] + x26*x39 + x26*x54 + x26*x67 + x26*L[34] + x28*x40 + x28*x55 + x28*x69 + x28*L[55] + x30*x41 + x30*x56 + x30*L[83] + x31*L[119] + x6*L[12] + x8*L[22] + x87*L[42] + x89*L[67] + x91*L[100] + x93*L[44] + x94*L[70] + x95*L[104] + x97*L[63] + x99*L[95] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += x*x44 + x*x45 + x*x81 + x*L[10] + x10*x50 + x10*x51 + x10*L[56] + x101*L[89] + x109*x226 + x109*L[68] + x111*L[102] + x114*L[101] + x12*L[84] + x123*x24 + x124*x26 + x125*x28 + x126*x15 + x127*x17 + x128*x19 + x133*x93 + x134*x87 + x15*x216 + x15*x224 + x15*L[23] + x17*x217 + x17*x225 + x17*L[41] + x19*x218 + x19*L[66] + x21*L[99] + x219*x24 + x220*x26 + x221*x28 + x222*x24 + x223*x26 + x24*L[25] + x26*L[44] + x28*L[70] + x30*L[104] + x42 + x43 + x46*x6 + x47*x6 + x48*x8 + x49*x8 + x6*x83 + x6*L[20] + x79 + x8*x85 + x8*L[35] + x87*L[59] + x89*L[90] + x93*L[61] + x94*L[93] + x97*L[87] + L[4];
#pragma omp atomic
Ls[5] += x*x160 + x*x174 + x*x78 + x*L[11] + x10*x163 + x10*x84 + x10*L[57] + x101*L[92] + x109*x153 + x109*L[73] + x111*L[108] + x114*L[107] + x12*L[85] + x132*x87 + x138*x15 + x139*x17 + x140*x19 + x149*x15 + x15*x72 + x15*L[26] + x150*x17 + x154 + x155 + x161*x6 + x162*x8 + x168*x24 + x17*x74 + x17*L[45] + x170*x26 + x172*x28 + x175*x6 + x176*x8 + x179*x93 + x19*x76 + x19*L[71] + x21*L[105] + x227*x24 + x228*x26 + x229*x28 + x230*x24 + x231*x26 + x24*L[28] + x26*L[48] + x28*L[75] + x3 + x30*L[110] + x6*x80 + x6*L[21] + x8*x82 + x8*L[36] + x87*L[62] + x89*L[94] + x93*L[64] + x94*L[97] + x97*L[90] + L[5];
#pragma omp atomic
Ls[6] += x*x189 + x*x190 + x*x208 + x*L[12] + x10*x195 + x10*x196 + x10*L[58] + x101*L[93] + x109*x152 + x109*L[74] + x111*L[109] + x114*L[108] + x12*L[86] + x130*x93 + x142*x24 + x143*x26 + x144*x28 + x146*x24 + x147*x26 + x15*x201 + x15*x232 + x15*x235 + x15*L[27] + x17*x203 + x17*x233 + x17*x236 + x17*L[46] + x180 + x181 + x184 + x19*x205 + x19*x234 + x19*L[72] + x191*x6 + x192*x6 + x193*x8 + x194*x8 + x21*L[106] + x210*x6 + x212*x8 + x215*x87 + x24*x64 + x24*L[29] + x26*x66 + x26*L[49] + x28*x68 + x28*L[76] + x30*L[111] + x6*L[22] + x8*L[37] + x87*L[63] + x89*L[95] + x93*L[65] + x94*L[98] + x97*L[91] + L[6];
#pragma omp atomic
Ls[7] += x*x237 + x*x238 + x10*x134 + x10*x241 + x10*L[59] + x101*L[96] + x109*x122 + x109*L[79] + x111*L[115] + x114*L[114] + x117*x24 + x118*x26 + x119*x28 + x12*L[87] + x126*x6 + x128*x87 + x131*x8 + x15*x34 + x15*x59 + x15*x75 + x15*L[30] + x17*x35 + x17*x60 + x17*x77 + x17*L[50] + x19*x36 + x19*x61 + x19*L[77] + x21*L[112] + x239*x6 + x24*x242 + x24*x246 + x24*L[32] + x240*x8 + x243*x26 + x244*x28 + x248*x26 + x249*x6 + x250*x8 + x251*x93 + x26*L[53] + x28*L[81] + x30*L[117] + x32 + x57 + x58*y + x6*L[23] + x71 + x8*L[38] + x87*L[66] + x89*L[99] + x93*L[68] + x94*L[102] + x97*L[94] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += x*x183 + x*x252 + x*x253 + x*L[14] + x10*x211 + x10*x256 + x10*L[60] + x101*L[97] + x109*x121 + x109*L[80] + x111*L[116] + x114*L[115] + x12*L[88] + x15*x186 + x15*x198 + x15*x204 + x15*L[31] + x157*x24 + x158*x26 + x159*x28 + x165*x24 + x166*x26 + x167*x28 + x17*x187 + x17*x199 + x17*x206 + x17*L[51] + x171*x24 + x173*x26 + x178*x93 + x182 + x188*x19 + x19*x200 + x19*L[78] + x197*y + x207*x6 + x209*x8 + x21*L[113] + x214*x87 + x24*L[33] + x254*x6 + x255*x8 + x257*x6 + x258*x8 + x26*L[54] + x28*L[82] + x30*L[118] + x6*L[24] + x8*L[39] + x87*L[67] + x89*L[100] + x93*L[69] + x94*L[103] + x97*L[95] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += x*x259 + x*x262 + x10*x133 + x10*x265 + x10*L[61] + x101*L[98] + x109*x119 + x109*L[81] + x111*L[117] + x114*L[116] + x117*x15 + x12*L[89] + x120*x17 + x122*x19 + x123*x6 + x125*x93 + x129*x8 + x15*x266 + x15*x270 + x15*L[32] + x17*x267 + x17*x272 + x17*L[52] + x19*x268 + x19*L[79] + x21*L[114] + x24*x39 + x24*x54 + x24*x67 + x24*L[34] + x26*x40 + x26*x55 + x26*x69 + x26*L[55] + x260*y + x263*x6 + x264*x8 + x274*x6 + x276*x8 + x277*x87 + x28*x41 + x28*x56 + x28*L[83] + x30*L[119] + x37 + x52 + x6*L[25] + x63 + x8*L[40] + x87*L[68] + x89*L[101] + x93*L[70] + x94*L[104] + x97*L[96] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += x*x46 + x*x47 + x*x83 + x*L[20] + x10*L[84] + x109*L[96] + x129*x24 + x130*x26 + x131*x15 + x132*x17 + x15*x278 + x15*x283 + x15*L[38] + x17*x279 + x17*L[62] + x19*L[94] + x24*x280 + x24*x282 + x24*L[40] + x26*x281 + x26*L[65] + x28*L[98] + x44 + x45 + x48*x6 + x49*x6 + x50*x8 + x51*x8 + x6*x85 + x6*L[35] + x8*L[56] + x81 + x87*L[87] + x93*L[89] + L[10];
#pragma omp atomic
Ls[11] += x*x161 + x*x175 + x*x80 + x*L[21] + x10*L[85] + x109*L[101] + x127*x15 + x128*x17 + x15*x217 + x15*x225 + x15*L[41] + x160 + x162*x6 + x163*x8 + x17*x218 + x17*L[66] + x174 + x176*x6 + x177*x24 + x178*x26 + x19*L[99] + x24*x284 + x24*x286 + x24*L[43] + x26*x285 + x26*L[69] + x28*L[103] + x6*x82 + x6*L[36] + x78 + x8*x84 + x8*L[57] + x87*L[90] + x93*L[92] + L[11];
#pragma omp atomic
Ls[12] += x*x191 + x*x192 + x*x210 + x*L[22] + x10*L[86] + x109*L[102] + x124*x24 + x125*x26 + x15*x213 + x15*x287 + x15*x289 + x15*L[42] + x17*x214 + x17*x288 + x17*L[67] + x189 + x19*L[100] + x190 + x193*x6 + x194*x6 + x195*x8 + x196*x8 + x208 + x212*x6 + x220*x24 + x221*x26 + x223*x24 + x24*L[44] + x26*L[70] + x28*L[104] + x6*L[37] + x8*L[58] + x87*L[91] + x93*L[93] + L[12];
#pragma omp atomic
Ls[13] += x*x239 + x*x249 + x10*L[87] + x109*L[107] + x131*x6 + x134*x8 + x137 + x139*x15 + x140*x17 + x148 + x15*x150 + x15*x74 + x15*L[45] + x151*x24 + x152*x26 + x17*x76 + x17*L[71] + x19*L[105] + x237 + x238 + x24*x245 + x24*x290 + x24*L[47] + x240*x6 + x241*x8 + x247*x26 + x250*x6 + x26*L[74] + x28*L[109] + x6*L[38] + x70 + x8*L[59] + x87*L[94] + x93*L[96] + L[13];
#pragma omp atomic
Ls[14] += x*x207 + x*x254 + x*x257 + x*L[24] + x10*L[88] + x109*L[108] + x15*x203 + x15*x233 + x15*x236 + x15*L[46] + x17*x205 + x17*x234 + x17*L[72] + x170*x24 + x172*x26 + x183 + x19*L[106] + x209*x6 + x211*x8 + x228*x24 + x229*x26 + x231*x24 + x24*L[48] + x252 + x253 + x255*x6 + x256*x8 + x258*x6 + x26*L[75] + x28*L[110] + x6*L[39] + x8*L[60] + x87*L[95] + x93*L[97] + L[14];
#pragma omp atomic
Ls[15] += x*x263 + x*x274 + x10*L[89] + x109*L[109] + x129*x6 + x133*x8 + x141 + x143*x24 + x144*x26 + x145 + x147*x24 + x15*x151 + x15*x269 + x15*x291 + x15*L[47] + x153*x17 + x17*x271 + x17*L[73] + x19*L[107] + x24*x66 + x24*L[49] + x259 + x26*x68 + x26*L[76] + x262 + x264*x6 + x265*x8 + x276*x6 + x28*L[111] + x6*L[40] + x62 + x8*L[61] + x87*L[96] + x93*L[98] + L[15];
#pragma omp atomic
Ls[16] += x*x292 + x*x293 + x10*L[90] + x109*L[114] + x120*x24 + x121*x26 + x127*x6 + x132*x8 + x15*x35 + x15*x60 + x15*x77 + x15*L[50] + x17*x36 + x17*x61 + x17*L[77] + x19*L[112] + x24*x296 + x24*x299 + x24*L[52] + x26*x297 + x26*L[80] + x28*L[116] + x294*x6 + x295*x8 + x300*x6 + x33 + x58 + x59*y + x6*L[41] + x73 + x8*L[62] + x87*L[99] + x93*L[101] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += x*x301 + x*x302 + x10*L[91] + x109*L[115] + x118*x24 + x119*x26 + x15*x187 + x15*x199 + x15*x206 + x15*L[51] + x17*x188 + x17*x200 + x17*L[78] + x185 + x19*L[113] + x197 + x198*y + x202 + x213*x6 + x215*x8 + x24*x243 + x24*x248 + x24*L[53] + x244*x26 + x26*L[81] + x28*L[117] + x303*x6 + x304*x8 + x305*x6 + x6*L[42] + x8*L[63] + x87*L[100] + x93*L[102] + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += x*x261 + x*x306 + x10*L[92] + x109*L[116] + x120*x15 + x122*x17 + x15*x267 + x15*x272 + x15*L[52] + x156 + x158*x24 + x159*x26 + x164 + x166*x24 + x167*x26 + x169 + x17*x268 + x17*L[79] + x173*x24 + x177*x6 + x179*x8 + x19*L[114] + x24*L[54] + x26*L[82] + x260 + x266*y + x273*x6 + x275*x8 + x28*L[118] + x307*x6 + x6*L[43] + x8*L[64] + x87*L[101] + x93*L[103] + L[18];
#pragma omp atomic
Ls[19] += x*x308 + x*x311 + x10*L[93] + x109*L[117] + x118*x15 + x121*x17 + x124*x6 + x130*x8 + x15*x314 + x15*x317 + x15*L[53] + x17*x315 + x17*L[80] + x19*L[115] + x24*x40 + x24*x55 + x24*x69 + x24*L[55] + x26*x41 + x26*x56 + x26*L[83] + x28*L[119] + x309*y + x312*x6 + x313*x8 + x319*x6 + x38 + x53 + x6*L[44] + x65 + x8*L[65] + x87*L[102] + x93*L[104] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += x*x48 + x*x49 + x*x85 + x*L[35] + x133*x24 + x134*x15 + x15*x320 + x15*L[59] + x17*L[90] + x24*x321 + x24*L[61] + x26*L[93] + x46 + x47 + x50*x6 + x51*x6 + x6*L[56] + x8*L[84] + x83 + L[20];
#pragma omp atomic
Ls[21] += x*x162 + x*x176 + x*x82 + x*L[36] + x132*x15 + x15*x279 + x15*L[62] + x161 + x163*x6 + x17*L[94] + x175 + x179*x24 + x24*x322 + x24*L[64] + x26*L[97] + x6*x84 + x6*L[57] + x8*L[85] + x80 + L[21];
#pragma omp atomic
Ls[22] += x*x193 + x*x194 + x*x212 + x*L[37] + x130*x24 + x15*x215 + x15*x323 + x15*L[63] + x17*L[95] + x191 + x192 + x195*x6 + x196*x6 + x210 + x24*x281 + x24*L[65] + x26*L[98] + x6*L[58] + x8*L[86] + L[22];
#pragma omp atomic
Ls[23] += x*x240 + x*x250 + x126 + x128*x15 + x134*x6 + x15*x218 + x15*L[66] + x17*L[99] + x216 + x224 + x226*x24 + x239 + x24*x251 + x24*L[68] + x241*x6 + x249 + x26*L[102] + x6*L[59] + x8*L[87] + L[23];
#pragma omp atomic
Ls[24] += x*x209 + x*x255 + x*x258 + x*L[39] + x15*x214 + x15*x288 + x15*L[67] + x17*L[100] + x178*x24 + x207 + x211*x6 + x24*x285 + x24*L[69] + x254 + x256*x6 + x257 + x26*L[103] + x6*L[60] + x8*L[88] + L[24];
#pragma omp atomic
Ls[25] += x*x264 + x*x276 + x123 + x125*x24 + x133*x6 + x15*x226 + x15*x277 + x15*L[68] + x17*L[101] + x219 + x221*x24 + x222 + x24*L[70] + x26*L[104] + x263 + x265*x6 + x274 + x6*L[61] + x8*L[89] + L[25];
#pragma omp atomic
Ls[26] += x*x294 + x*x300 + x132*x6 + x138 + x140*x15 + x149 + x15*x76 + x15*L[71] + x153*x24 + x17*L[105] + x24*x298 + x24*L[73] + x26*L[108] + x292 + x293 + x295*x6 + x6*L[62] + x72 + x8*L[90] + L[26];
#pragma omp atomic
Ls[27] += x*x303 + x*x305 + x15*x205 + x15*x234 + x15*L[72] + x152*x24 + x17*L[106] + x201 + x215*x6 + x232 + x235 + x24*x247 + x24*L[74] + x26*L[109] + x301 + x302 + x304*x6 + x6*L[63] + x8*L[91] + L[27];
#pragma omp atomic
Ls[28] += x*x273 + x*x307 + x15*x153 + x15*x271 + x15*L[73] + x168 + x17*L[107] + x172*x24 + x179*x6 + x227 + x229*x24 + x230 + x24*L[75] + x26*L[110] + x261 + x275*x6 + x306 + x6*L[64] + x8*L[92] + L[28];
#pragma omp atomic
Ls[29] += x*x312 + x*x319 + x130*x6 + x142 + x144*x24 + x146 + x15*x152 + x15*x316 + x15*L[74] + x17*L[108] + x24*x68 + x24*L[76] + x26*L[111] + x308 + x311 + x313*x6 + x6*L[65] + x64 + x8*L[93] + L[29];
#pragma omp atomic
Ls[30] += x*x324 + x*x325 + x122*x24 + x128*x6 + x15*x36 + x15*x61 + x15*L[77] + x17*L[112] + x24*x327 + x24*L[79] + x26*L[115] + x326*x6 + x34 + x59 + x6*L[66] + x60*y + x75 + x8*L[94] + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += x*x328 + x*x329 + x121*x24 + x15*x188 + x15*x200 + x15*L[78] + x17*L[113] + x186 + x198 + x199*y + x204 + x214*x6 + x24*x297 + x24*L[80] + x26*L[116] + x330*x6 + x6*L[67] + x8*L[95] + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += x*x331 + x117 + x119*x24 + x122*x15 + x15*x268 + x15*L[79] + x17*L[114] + x24*x244 + x24*L[81] + x242 + x246 + x251*x6 + x26*L[117] + x266 + x267*y + x270 + x277*x6 + x6*L[68] + x8*L[96] + L[32];
#pragma omp atomic
Ls[33] += x*x310 + x*x332 + x121*x15 + x15*x315 + x15*L[80] + x157 + x159*x24 + x165 + x167*x24 + x17*L[115] + x171 + x178*x6 + x24*L[82] + x26*L[118] + x309 + x314*y + x318*x6 + x6*L[69] + x8*L[97] + L[33];
#pragma omp atomic
Ls[34] += x*x333 + x*x336 + x119*x15 + x125*x6 + x15*x338 + x15*L[81] + x17*L[116] + x24*x41 + x24*x56 + x24*L[83] + x26*L[119] + x334*y + x337*x6 + x39 + x54 + x6*L[70] + x67 + x8*L[98] + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += x*x50 + x*x51 + x*L[56] + x15*L[87] + x24*L[89] + x48 + x49 + x6*L[84] + x85 + L[35];
#pragma omp atomic
Ls[36] += x*x163 + x*x84 + x*L[57] + x15*L[90] + x162 + x176 + x24*L[92] + x6*L[85] + x82 + L[36];
#pragma omp atomic
Ls[37] += x*x195 + x*x196 + x*L[58] + x15*L[91] + x193 + x194 + x212 + x24*L[93] + x6*L[86] + L[37];
#pragma omp atomic
Ls[38] += x*x241 + x131 + x15*L[94] + x24*L[96] + x240 + x250 + x278 + x283 + x6*L[87] + L[38];
#pragma omp atomic
Ls[39] += x*x211 + x*x256 + x*L[60] + x15*L[95] + x209 + x24*L[97] + x255 + x258 + x6*L[88] + L[39];
#pragma omp atomic
Ls[40] += x*x265 + x129 + x15*L[96] + x24*L[98] + x264 + x276 + x280 + x282 + x6*L[89] + L[40];
#pragma omp atomic
Ls[41] += x*x295 + x127 + x15*L[99] + x217 + x225 + x24*L[101] + x294 + x300 + x6*L[90] + L[41];
#pragma omp atomic
Ls[42] += x*x304 + x15*L[100] + x213 + x24*L[102] + x287 + x289 + x303 + x305 + x6*L[91] + L[42];
#pragma omp atomic
Ls[43] += x*x275 + x15*L[101] + x177 + x24*L[103] + x273 + x284 + x286 + x307 + x6*L[92] + L[43];
#pragma omp atomic
Ls[44] += x*x313 + x124 + x15*L[102] + x220 + x223 + x24*L[104] + x312 + x319 + x6*L[93] + L[44];
#pragma omp atomic
Ls[45] += x*x326 + x139 + x15*L[105] + x150 + x24*L[107] + x324 + x325 + x6*L[94] + x74 + L[45];
#pragma omp atomic
Ls[46] += x*x330 + x15*L[106] + x203 + x233 + x236 + x24*L[108] + x328 + x329 + x6*L[95] + L[46];
#pragma omp atomic
Ls[47] += x15*L[107] + x151 + x24*L[109] + x245 + x269 + x290 + x291 + x331 + x6*L[96] + L[47];
#pragma omp atomic
Ls[48] += x*x318 + x15*L[108] + x170 + x228 + x231 + x24*L[110] + x310 + x332 + x6*L[97] + L[48];
#pragma omp atomic
Ls[49] += x*x337 + x143 + x147 + x15*L[109] + x24*L[111] + x333 + x336 + x6*L[98] + x66 + L[49];
#pragma omp atomic
Ls[50] += x*x339 + x15*L[112] + x24*L[114] + x35 + x6*L[99] + x60 + x61*y + x77 + y*L[77] + L[50];
#pragma omp atomic
Ls[51] += x*x340 + x15*L[113] + x187 + x199 + x200*y + x206 + x24*L[115] + x6*L[100] + y*L[78] + L[51];
#pragma omp atomic
Ls[52] += x120 + x15*L[114] + x24*L[116] + x267 + x268*y + x272 + x296 + x299 + x6*L[101] + L[52];
#pragma omp atomic
Ls[53] += x118 + x15*L[115] + x24*L[117] + x243 + x248 + x314 + x315*y + x317 + x6*L[102] + L[53];
#pragma omp atomic
Ls[54] += x*x335 + x15*L[116] + x158 + x166 + x173 + x24*L[118] + x334 + x338*y + x6*L[103] + L[54];
#pragma omp atomic
Ls[55] += x*x341 + x15*L[117] + x24*L[119] + x342*y + x40 + x55 + x6*L[104] + x69 + z*L[83] + L[55];
#pragma omp atomic
Ls[56] += x*L[84] + x50 + x51 + L[56];
#pragma omp atomic
Ls[57] += x*L[85] + x163 + x84 + L[57];
#pragma omp atomic
Ls[58] += x*L[86] + x195 + x196 + L[58];
#pragma omp atomic
Ls[59] += x134 + x241 + x320 + L[59];
#pragma omp atomic
Ls[60] += x*L[88] + x211 + x256 + L[60];
#pragma omp atomic
Ls[61] += x133 + x265 + x321 + L[61];
#pragma omp atomic
Ls[62] += x132 + x279 + x295 + L[62];
#pragma omp atomic
Ls[63] += x215 + x304 + x323 + L[63];
#pragma omp atomic
Ls[64] += x179 + x275 + x322 + L[64];
#pragma omp atomic
Ls[65] += x130 + x281 + x313 + L[65];
#pragma omp atomic
Ls[66] += x128 + x218 + x326 + L[66];
#pragma omp atomic
Ls[67] += x214 + x288 + x330 + L[67];
#pragma omp atomic
Ls[68] += x226 + x251 + x277 + L[68];
#pragma omp atomic
Ls[69] += x178 + x285 + x318 + L[69];
#pragma omp atomic
Ls[70] += x125 + x221 + x337 + L[70];
#pragma omp atomic
Ls[71] += x140 + x339 + x76 + L[71];
#pragma omp atomic
Ls[72] += x205 + x234 + x340 + L[72];
#pragma omp atomic
Ls[73] += x153 + x271 + x298 + L[73];
#pragma omp atomic
Ls[74] += x152 + x247 + x316 + L[74];
#pragma omp atomic
Ls[75] += x172 + x229 + x335 + L[75];
#pragma omp atomic
Ls[76] += x144 + x341 + x68 + L[76];
#pragma omp atomic
Ls[77] += x36 + x61 + y*L[112] + L[77];
#pragma omp atomic
Ls[78] += x188 + x200 + y*L[113] + L[78];
#pragma omp atomic
Ls[79] += x122 + x268 + x327 + L[79];
#pragma omp atomic
Ls[80] += x121 + x297 + x315 + L[80];
#pragma omp atomic
Ls[81] += x119 + x244 + x338 + L[81];
#pragma omp atomic
Ls[82] += x159 + x167 + x342 + L[82];
#pragma omp atomic
Ls[83] += x41 + x56 + z*L[119] + L[83];
#pragma omp atomic
Ls[84] += L[84];
#pragma omp atomic
Ls[85] += L[85];
#pragma omp atomic
Ls[86] += L[86];
#pragma omp atomic
Ls[87] += L[87];
#pragma omp atomic
Ls[88] += L[88];
#pragma omp atomic
Ls[89] += L[89];
#pragma omp atomic
Ls[90] += L[90];
#pragma omp atomic
Ls[91] += L[91];
#pragma omp atomic
Ls[92] += L[92];
#pragma omp atomic
Ls[93] += L[93];
#pragma omp atomic
Ls[94] += L[94];
#pragma omp atomic
Ls[95] += L[95];
#pragma omp atomic
Ls[96] += L[96];
#pragma omp atomic
Ls[97] += L[97];
#pragma omp atomic
Ls[98] += L[98];
#pragma omp atomic
Ls[99] += L[99];
#pragma omp atomic
Ls[100] += L[100];
#pragma omp atomic
Ls[101] += L[101];
#pragma omp atomic
Ls[102] += L[102];
#pragma omp atomic
Ls[103] += L[103];
#pragma omp atomic
Ls[104] += L[104];
#pragma omp atomic
Ls[105] += L[105];
#pragma omp atomic
Ls[106] += L[106];
#pragma omp atomic
Ls[107] += L[107];
#pragma omp atomic
Ls[108] += L[108];
#pragma omp atomic
Ls[109] += L[109];
#pragma omp atomic
Ls[110] += L[110];
#pragma omp atomic
Ls[111] += L[111];
#pragma omp atomic
Ls[112] += L[112];
#pragma omp atomic
Ls[113] += L[113];
#pragma omp atomic
Ls[114] += L[114];
#pragma omp atomic
Ls[115] += L[115];
#pragma omp atomic
Ls[116] += L[116];
#pragma omp atomic
Ls[117] += L[117];
#pragma omp atomic
Ls[118] += L[118];
#pragma omp atomic
Ls[119] += L[119];

}

void L2P_8(double x, double y, double z, double * L, double * F) {
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
x0 = x*y;
x1 = x*z;
x2 = y*z;
x3 = (x*x);
x4 = (1.0/2.0)*x3;
x5 = (x*x*x);
x6 = (1.0/6.0)*x5;
x7 = (x*x*x*x);
x8 = (1.0/24.0)*x7;
x9 = pow(x, 5);
x10 = (1.0/120.0)*x9;
x11 = (1.0/720.0)*pow(x, 6);
x12 = (y*y);
x13 = (1.0/2.0)*x12;
x14 = (y*y*y);
x15 = (1.0/6.0)*x14;
x16 = (y*y*y*y);
x17 = (1.0/24.0)*x16;
x18 = pow(y, 5);
x19 = (1.0/120.0)*x18;
x20 = (1.0/720.0)*pow(y, 6);
x21 = (z*z);
x22 = (1.0/2.0)*x21;
x23 = (z*z*z);
x24 = (1.0/6.0)*x23;
x25 = (z*z*z*z);
x26 = (1.0/24.0)*x25;
x27 = pow(z, 5);
x28 = (1.0/120.0)*x27;
x29 = (1.0/720.0)*pow(z, 6);
x30 = (1.0/4.0)*x3;
x31 = x12*x30;
x32 = (1.0/12.0)*x3;
x33 = x14*x32;
x34 = (1.0/48.0)*x3;
x35 = x16*x34;
x36 = (1.0/240.0)*x3;
x37 = x21*x30;
x38 = x23*x32;
x39 = x25*x34;
x40 = (1.0/12.0)*x5;
x41 = x12*x40;
x42 = (1.0/36.0)*x5;
x43 = x14*x42;
x44 = (1.0/144.0)*x5;
x45 = x21*x40;
x46 = x23*x42;
x47 = (1.0/48.0)*x7;
x48 = x12*x47;
x49 = (1.0/144.0)*x7;
x50 = x21*x47;
x51 = (1.0/240.0)*x9;
x52 = x12*x21;
x53 = (1.0/4.0)*x52;
x54 = x12*x23;
x55 = (1.0/12.0)*x54;
x56 = (1.0/48.0)*x12*x25;
x57 = x14*x21;
x58 = (1.0/12.0)*x57;
x59 = (1.0/36.0)*x14*x23;
x60 = (1.0/48.0)*x16*x21;
x61 = (1.0/24.0)*x3;
#pragma omp atomic
F[0] += (1.0/5040.0)*pow(x, 7)*L[84] + x*x13*L[13] + x*x15*L[26] + x*x17*L[45] + x*x19*L[71] + x*x20*L[105] + x*x22*L[15] + x*x24*L[29] + x*x26*L[49] + x*x28*L[76] + x*x29*L[111] + x*x53*L[47] + x*x55*L[74] + x*x56*L[109] + x*x58*L[73] + x*x59*L[108] + x*x60*L[107] + x*L[1] + x0*x22*L[28] + x0*x24*L[48] + x0*x26*L[75] + x0*x28*L[110] + x0*z*L[14] + x0*L[5] + x1*x13*L[27] + x1*x15*L[46] + x1*x17*L[72] + x1*x19*L[106] + x1*L[6] + x10*x2*L[88] + x10*y*L[57] + x10*z*L[58] + x10*L[35] + x11*y*L[85] + x11*z*L[86] + x11*L[56] + (1.0/240.0)*x12*x27*L[117] + x12*x51*L[87] + x13*z*L[17] + x13*L[7] + (1.0/144.0)*x14*x25*L[116] + x14*x49*L[90] + x15*z*L[31] + x15*L[16] + (1.0/144.0)*x16*x23*L[115] + x16*x44*L[94] + x17*z*L[51] + x17*L[30] + (1.0/240.0)*x18*x21*L[114] + x18*x36*L[99] + x19*z*L[78] + x19*L[50] + x2*x4*L[24] + x2*x6*L[39] + x2*x8*L[60] + x2*L[8] + x20*z*L[113] + x20*L[77] + x21*x51*L[89] + x22*y*L[18] + x22*L[9] + x23*x49*L[93] + x24*y*L[33] + x24*L[19] + x25*x44*L[98] + x26*y*L[54] + x26*L[34] + x27*x36*L[104] + x28*y*L[82] + x28*L[55] + x29*y*L[118] + x29*L[83] + (1.0/8.0)*x3*x52*L[68] + x31*z*L[42] + x31*L[23] + x33*z*L[67] + x33*L[41] + x35*z*L[100] + x35*L[66] + x37*y*L[43] + x37*L[25] + x38*y*L[69] + x38*L[44] + x39*y*L[103] + x39*L[70] + x4*y*L[11] + x4*z*L[12] + x4*L[4] + x41*z*L[63] + x41*L[38] + x43*z*L[95] + x43*L[62] + x45*y*L[64] + x45*L[40] + x46*y*L[97] + x46*L[65] + x48*z*L[91] + x48*L[59] + (1.0/24.0)*x5*x52*L[96] + x50*y*L[92] + x50*L[61] + x53*L[32] + x54*x61*L[102] + x55*L[53] + x56*L[81] + x57*x61*L[101] + x58*L[52] + x59*L[80] + x6*y*L[21] + x6*z*L[22] + x6*L[10] + x60*L[79] + x8*y*L[36] + x8*z*L[37] + x8*L[20] + (1.0/5040.0)*pow(y, 7)*L[112] + y*L[2] + (1.0/5040.0)*pow(z, 7)*L[119] + z*L[3] + L[0];

}

void M2P_8(double x, double y, double z, double * M, double * F) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*y;
x8 = x*z;
x9 = y*z;
x10 = pow(x2, -3.5);
x11 = 15.0*x10;
x12 = x7*z;
x13 = -x4;
x14 = x0*x6;
x15 = x1*x6;
x16 = 9.0*x5;
x17 = x0*x11;
x18 = -x17;
x19 = x*(x16 + x18);
x20 = x18 + x6;
x21 = x20*y;
x22 = x1*x11;
x23 = -x22;
x24 = y*(x16 + x23);
x25 = x20*z;
x26 = z*(x23 + x6);
x27 = 1.0*x;
x28 = x27*(x22 - x6);
x29 = 45.0*x10;
x30 = -x29;
x31 = pow(x2, -4.5);
x32 = x0*x31;
x33 = 105.0*x32;
x34 = x30 + x33;
x35 = x34*x7;
x36 = x34*x8;
x37 = -x11;
x38 = x9*(x33 + x37);
x39 = x1*x31;
x40 = 105.0*x39;
x41 = x30 + x40;
x42 = x41*x9;
x43 = x27*y;
x44 = x41*x43;
x45 = x37 + x40;
x46 = x27*z;
x47 = x45*x46;
x48 = 315.0*x31;
x49 = pow(x2, -5.5);
x50 = 945.0*x49;
x51 = x0*x50;
x52 = x48 - x51;
x53 = x12*x52;
x54 = x1*x50;
x55 = x27*x9;
x56 = x55*(-x48 + x54);
x57 = 90.0*x10;
x58 = (x*x*x*x);
x59 = 105.0*x31;
x60 = x58*x59;
x61 = (y*y*y*y);
x62 = x59*x61;
x63 = -225.0*x10;
x64 = x50*x58;
x65 = -x64;
x66 = x*(1050.0*x32 + x63 + x65);
x67 = x30 + 630.0*x32 + x65;
x68 = x67*y;
x69 = x50*x61;
x70 = -x69;
x71 = y*(1050.0*x39 + x63 + x70);
x72 = x67*z;
x73 = 630.0*x39;
x74 = x30 + x70 + x73;
x75 = x74*z;
x76 = x27*(x29 + x69 - x73);
x77 = 1575.0*x31;
x78 = x0*x49;
x79 = pow(x2, -6.5);
x80 = x58*x79;
x81 = 10395.0*x80;
x82 = x77 - 9450.0*x78 + x81;
x83 = x7*x82;
x84 = x8*x82;
x85 = 5670.0*x78;
x86 = x48 + x81 - x85;
x87 = x86*x9;
x88 = x1*x49;
x89 = x61*x79;
x90 = 10395.0*x89;
x91 = x77 - 9450.0*x88 + x90;
x92 = x9*x91;
x93 = x43*x91;
x94 = x48 - 5670.0*x88 + x90;
x95 = x46*x94;
x96 = 14175.0*x49;
x97 = -x96;
x98 = x0*x79;
x99 = 103950.0*x98;
x100 = pow(x2, -7.5);
x101 = 135135.0*x100;
x102 = x101*x58;
x103 = -x102 + x97 + x99;
x104 = x103*x12;
x105 = x1*x79;
x106 = 103950.0*x105;
x107 = x101*x61;
x108 = x55*(-x106 + x107 + x96);
x109 = pow(x, 6);
x110 = 10395.0*x79;
x111 = x109*x110;
x112 = pow(y, 6);
x113 = x110*x112;
x114 = 11025.0*x31;
x115 = 99225.0*x49;
x116 = x101*x109;
x117 = -x116;
x118 = x117 + 218295.0*x80;
x119 = x*(-x0*x115 + x114 + x118);
x120 = 42525.0*x49;
x121 = -x0*x120 + x117 + x77 + 155925.0*x80;
x122 = x121*y;
x123 = x101*x112;
x124 = -x123;
x125 = x124 + 218295.0*x89;
x126 = y*(-x1*x115 + x114 + x125);
x127 = x121*z;
x128 = x1*x120;
x129 = 155925.0*x89;
x130 = x124 - x128 + x129 + x77;
x131 = x130*z;
x132 = x27*(x123 + x128 - x129 - x77);
x133 = -x115;
x134 = x100*x58;
x135 = 2837835.0*x134;
x136 = pow(x2, -8.5);
x137 = 2027025.0*x136;
x138 = x109*x137;
x139 = x133 - x135 + x138 + 1091475.0*x98;
x140 = x139*x7;
x141 = x139*x8;
x142 = 467775.0*x98;
x143 = 2027025.0*x134;
x144 = x9*(x138 + x142 - x143 + x97);
x145 = x100*x61;
x146 = 2837835.0*x145;
x147 = x112*x137;
x148 = 1091475.0*x105 + x133 - x146 + x147;
x149 = x148*x9;
x150 = x1*x33;
x151 = x148*x43;
x152 = 2027025.0*x145;
x153 = x46*(467775.0*x105 + x147 - x152 + x97);
x154 = x1*x48;
x155 = -x1*x51;
x156 = x*(x154 + x155 + x34);
x157 = x0*x48;
x158 = y*(x155 + x157 + x41);
x159 = z*(x155 + x33 + x45);
x160 = -2835.0*x88;
x161 = x1*x98;
x162 = 10395.0*x161;
x163 = x160 + x162;
x164 = 945.0*x31;
x165 = -2835.0*x78;
x166 = x164 + x165;
x167 = x7*(x163 + x166);
x168 = x8*(x163 + x52);
x169 = x9*(x162 + x165 + x48 - x54);
x170 = 31185.0*x105;
x171 = x0*x1;
x172 = -8505.0*x49;
x173 = x172 + 31185.0*x98;
x174 = x12*(-x101*x171 + x170 + x173);
x175 = 3783780.0*x100;
x176 = pow(x, 8)*x137;
x177 = x137*pow(y, 8);
x178 = 105.0*x10;
x179 = -x150 - 12.0*x5;
x180 = 120.0*x10;
x181 = x0*x90;
x182 = x1*x85;
x183 = -x182;
x184 = x1*x81;
x185 = 62370.0*x161;
x186 = -x0*x107;
x187 = x185 + x186;
x188 = 31185.0*x89;
x189 = x188 - 17010.0*x88;
x190 = x*(x166 + x187 + x189);
x191 = x1*x96;
x192 = -x191;
x193 = x1*x99;
x194 = -x1*x102;
x195 = x*(x192 + x193 + x194 + x82);
x196 = x0*x96;
x197 = -x196;
x198 = y*(x186 + x193 + x197 + x91);
x199 = x160 + x185 + x194;
x200 = 31185.0*x80;
x201 = 17010.0*x78;
x202 = x164 + x200 - x201;
x203 = y*(x199 + x202);
x204 = z*(x165 + x187 + x94);
x205 = z*(x199 + x86);
x206 = 155925.0*x98;
x207 = x100*x171;
x208 = -1351350.0*x207;
x209 = -x120 + x208;
x210 = x0*x137*x61;
x211 = -405405.0*x145 + x210;
x212 = x7*(311850.0*x105 + x206 + x209 + x211);
x213 = 155925.0*x105;
x214 = x137*x58;
x215 = x1*x214;
x216 = -405405.0*x134 + x215;
x217 = x7*(x209 + x213 + x216 + 311850.0*x98);
x218 = -810810.0*x207;
x219 = x8*(187110.0*x105 + x173 + x211 + x218);
x220 = x8*(x103 + x208 + x213 + x215);
x221 = x9*(x106 - x107 + x206 + x208 + x210 + x97);
x222 = x9*(x170 + x172 + x216 + x218 + 187110.0*x98);
x223 = 15120.0*x49;
x224 = -x111;
x225 = -x184;
x226 = 270.0*x10 + x182;
x227 = -x113;
x228 = -x181;
x229 = -x156;
x230 = -x158;
x231 = -x159;
x232 = -x167;
x233 = -x168;
x234 = -x169;
x235 = -x174;
x236 = x1*x142;
x237 = x0*x147;
x238 = x0*x152;
x239 = x1*x78;
x240 = x1*x138;
x241 = x1*x143;
x242 = 16065.0*x49;
x243 = -x1*x201 - 360.0*x10;
x244 = 3918915.0*x100;
x245 = -x176;
x246 = -x240;
x247 = -x236 - 12600.0*x31;
x248 = -x177;
x249 = -x237;
x250 = x0*x145;
x251 = 810810.0*x250;
x252 = -x251;
x253 = x214*x61;
x254 = x1*x134;
x255 = 810810.0*x254;
x256 = x253 - x255;
x257 = 17010.0*x49;
x258 = -x119;
x259 = -x195;
x260 = -x190;
x261 = -x122;
x262 = -x203;
x263 = -x126;
x264 = -x198;
x265 = -x127;
x266 = -x205;
x267 = -x131;
x268 = -x204;
x269 = -x140;
x270 = -x217;
x271 = -x151;
x272 = -x212;
x273 = -x141;
x274 = -x220;
x275 = -x153;
x276 = -x219;
x277 = -x144;
x278 = -x222;
x279 = -x149;
x280 = -x221;
x281 = x123 + x249;
x282 = -841995.0*x161 - x253 - 2520.0*x31;
x283 = x116 + x246;
x284 = 4054050.0*x100;
x285 = 4054050.0*x136;
x286 = x1*x109;
x287 = 1309770.0*x161 + 15120.0*x31;
x288 = x0*x112;
x289 = x58*x61;
x290 = 4189185.0*x100;
x291 = 6081075.0*x136;
x292 = -2993760.0*x161 - x289*x291 - 20160.0*x31;
x293 = 4324320.0*x100;
x294 = 8108100.0*x136;
#pragma omp atomic
F[0] += -x*x4*M[0] + x104*M[87] - x108*M[105] - x11*x12*M[13] + x119*M[83] + x122*M[84] + x126*M[111] + x127*M[85] + x131*M[112] - x132*M[104] + x140*M[120] + x141*M[121] + x144*M[123] + x149*M[156] + x151*M[147] + x153*M[148] + x156*M[37] + x158*M[40] + x159*M[41] + x167*M[61] + x168*M[62] + x169*M[66] + x174*M[94] + x19*M[9] + x190*M[93] + x195*M[86] + x198*M[98] + x203*M[89] + x204*M[99] + x205*M[90] + x21*M[10] + x212*M[134] + x217*M[125] + x219*M[135] + x220*M[126] + x221*M[141] + x222*M[130] + x24*M[15] + x25*M[11] + x26*M[16] - x28*M[12] + x35*M[20] + x36*M[21] + x38*M[23] - x4*y*M[1] - x4*z*M[2] + x42*M[30] + x44*M[25] + x47*M[26] + x53*M[38] - x56*M[45] + x6*x7*M[4] + x6*x8*M[5] + x6*x9*M[7] + x66*M[34] + x68*M[35] + x71*M[49] + x72*M[36] + x75*M[50] - x76*M[44] + x83*M[56] + x84*M[57] + x87*M[59] + x92*M[77] + x93*M[70] + x95*M[71] + (-x104 + x235)*M[96] + (x108 + x235)*M[107] + (x13 + x14)*M[3] + (x13 + x15)*M[6] + (x132 + x260)*M[106] + (-x19 + x28)*M[14] + (-x21 - x24)*M[17] + (x229 - x66)*M[39] + (x229 + x76)*M[46] + (x230 - x68)*M[42] + (x230 - x71)*M[51] + (x231 - x72)*M[43] + (x231 - x75)*M[52] + (x232 - x83)*M[63] + (x232 - x93)*M[72] + (x233 - x84)*M[64] + (x233 - x95)*M[73] + (x234 - x87)*M[68] + (x234 - x92)*M[79] + (-x25 - x26)*M[18] + (x258 + x259)*M[88] + (x259 + x260)*M[95] + (x261 + x262)*M[91] + (x262 + x264)*M[100] + (x263 + x264)*M[113] + (x265 + x266)*M[92] + (x266 + x268)*M[101] + (x267 + x268)*M[114] + (x269 + x270)*M[127] + (x270 + x272)*M[136] + (x271 + x272)*M[149] + (x273 + x274)*M[128] + (x274 + x276)*M[137] + (x275 + x276)*M[150] + (x277 + x278)*M[132] + (x278 + x280)*M[143] + (x279 + x280)*M[158] + (-x35 - x44)*M[27] + (-x36 - x47)*M[28] + (-x38 - x42)*M[32] + (-x53 + x56)*M[47] + (x104 - x108 + 2*x174)*M[109] + (x119 + x190 + 2*x195)*M[97] + (x122 + x198 + 2*x203)*M[102] + (x126 + 2*x198 + x203)*M[115] + (x127 + x204 + 2*x205)*M[103] + (x131 + 2*x204 + x205)*M[116] + (-x132 + 2*x190 + x195)*M[108] + (-x14 - x15 + 2.0*x3)*M[8] + (x140 + x212 + 2*x217)*M[138] + (x141 + x219 + 2*x220)*M[139] + (x144 + x221 + 2*x222)*M[145] + (x149 + 2*x221 + x222)*M[160] + (x150 + x20 + x23)*M[22] + (x151 + 2*x212 + x217)*M[151] + (x153 + 2*x219 + x220)*M[152] + (2*x156 + x66 - x76)*M[48] + (2*x158 + x68 + x71)*M[53] + (2*x159 + x72 + x75)*M[54] + (2*x167 + x83 + x93)*M[74] + (2*x168 + x84 + x95)*M[75] + (2*x169 + x87 + x92)*M[81] + (-x0*x57 + x16 + x60)*M[19] + (-x1*x57 + x16 + x62)*M[29] + (x111 + 4725.0*x32 - x58*x96 + x63)*M[55] + (x113 + 4725.0*x39 - x61*x96 + x63)*M[76] + (x132 - 3*x190 - 3*x195 + x258)*M[110] + (x154 + x183 + x184 + x67)*M[58] + (x157 + x181 + x183 + x74)*M[65] + (-3*x198 - 3*x203 + x261 + x263)*M[117] + (-3*x204 - 3*x205 + x265 + x267)*M[118] + (-3*x212 - 3*x217 + x269 + x271)*M[153] + (-3*x219 - 3*x220 + x273 + x275)*M[154] + (-3*x221 - 3*x222 + x277 + x279)*M[162] + (x0*x178 + x179 + x22 - x60)*M[24] + (x1*x178 + x17 + x179 - x62)*M[31] + (x121 + x192 + x236 + x240 - x241)*M[122] + (x130 + x197 + x236 + x237 - x238)*M[140] + (374220.0*x161 + x189 + x202 + x252 + x256)*M[129] + (-x109*x175 + x114 + x176 - 396900.0*x78 + 2182950.0*x80)*M[119] + (-x112*x175 + x114 + x177 - 396900.0*x88 + 2182950.0*x89)*M[155] + (-x154 + x223*x58 + x224 + x225 + x226 - 5355.0*x32)*M[60] + (-x157 + x223*x61 + x226 + x227 + x228 - 5355.0*x39)*M[78] + (-x0*x180 - x1*x180 + 210.0*x1*x32 + 24.0*x5 + x60 + x62)*M[33] + (x0*x146 - x200 + x255 + x281 + x282 + 31185.0*x78 + 59535.0*x88 - 187110.0*x89)*M[142] + (-x0*x164 - x1*x164 + x225 + x228 + 11340.0*x239 + x57 + x64 + x69)*M[67] + (20790.0*x0*x89 + x113 + x184 - x242*x61 + x243 + 1260.0*x32 + 6300.0*x39 + x65)*M[80] + (x1*x135 - x188 + x251 + x282 + x283 + 59535.0*x78 - 187110.0*x80 + 31185.0*x88)*M[131] + (20790.0*x1*x80 + x111 + x181 - x242*x58 + x243 + 6300.0*x32 + 1260.0*x39 + x70)*M[69] + (x109*x244 + x191 + x241 + x245 + x246 + x247 + 439425.0*x78 - 2338875.0*x80)*M[124] + (x112*x244 + x196 + x238 + x247 + x248 + x249 + 439425.0*x88 - 2338875.0*x89)*M[157] + (-x0*x188 - x1*x200 + 720.0*x10 + x224 + x227 + 34020.0*x239 + x257*x58 + x257*x61 - 7560.0*x32 - 7560.0*x39)*M[82] + (-x112*x284 + x177 + x200 - 4864860.0*x250 + x256 + x285*x288 + x287 - 45360.0*x78 - 498960.0*x88 + 2525985.0*x89)*M[159] + (x118 + x125 + 1683990.0*x161 + x237 + x240 - 3648645.0*x250 - 3648645.0*x254 + x285*x289 + 5040.0*x31 - 90720.0*x78 - 90720.0*x88)*M[144] + (-x109*x284 + x176 + x188 + x252 + x253 - 4864860.0*x254 + x285*x286 + x287 - 498960.0*x78 + 2525985.0*x80 - 45360.0*x88)*M[133] + (x109*x290 + x245 + 4459455.0*x250 + 8513505.0*x254 + x281 - x286*x291 + x292 + 589680.0*x78 - 2744280.0*x80 + 136080.0*x88 - 249480.0*x89)*M[146] + (x112*x290 + x248 + 8513505.0*x250 + 4459455.0*x254 + x283 - x288*x291 + x292 + 136080.0*x78 - 249480.0*x80 + 589680.0*x88 - 2744280.0*x89)*M[161] + (-x109*x293 - x112*x293 + 12162150.0*x136*x289 + 5987520.0*x161 + x176 + x177 - 12972960.0*x250 - 12972960.0*x254 + x286*x294 + x288*x294 + 40320.0*x31 - 725760.0*x78 + 2993760.0*x80 - 725760.0*x88 + 2993760.0*x89)*M[163];

}

void P2M_9(double x, double y, double z, double q, double * M) {
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
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = x1*z;
x9 = (z*z);
x10 = (x*x*x);
x11 = (1.0/6.0)*q;
x12 = (1.0/2.0)*x3;
x13 = (1.0/2.0)*x0;
x14 = (y*y*y);
x15 = (1.0/2.0)*x7;
x16 = (1.0/2.0)*x9;
x17 = (z*z*z);
x18 = (x*x*x*x);
x19 = (1.0/24.0)*q;
x20 = (1.0/6.0)*x10;
x21 = q*x7;
x22 = (1.0/4.0)*x3;
x23 = q*x9;
x24 = (1.0/6.0)*x0;
x25 = (y*y*y*y);
x26 = (1.0/6.0)*x14;
x27 = (1.0/4.0)*x9;
x28 = (1.0/6.0)*x17;
x29 = (z*z*z*z);
x30 = pow(x, 5);
x31 = (1.0/120.0)*q;
x32 = (1.0/24.0)*x18;
x33 = (1.0/12.0)*x10;
x34 = (1.0/12.0)*x14;
x35 = q*x3;
x36 = x2*x7;
x37 = x1*x9;
x38 = (1.0/12.0)*x17;
x39 = (1.0/24.0)*x0;
x40 = x0*x7;
x41 = pow(y, 5);
x42 = (1.0/24.0)*x25;
x43 = (1.0/24.0)*x29;
x44 = pow(z, 5);
x45 = pow(x, 6);
x46 = (1.0/720.0)*q;
x47 = (1.0/120.0)*x30;
x48 = (1.0/48.0)*x18;
x49 = q*x14;
x50 = (1.0/36.0)*x10;
x51 = q*x17;
x52 = (1.0/48.0)*x35;
x53 = x2*x3;
x54 = x3*x9;
x55 = x1*x3;
x56 = (1.0/120.0)*x0;
x57 = x0*x9;
x58 = pow(y, 6);
x59 = (1.0/120.0)*x41;
x60 = (1.0/48.0)*x25;
x61 = (1.0/36.0)*x17;
x62 = (1.0/48.0)*x29;
x63 = (1.0/120.0)*x44;
x64 = pow(z, 6);
x65 = pow(x, 7);
x66 = (1.0/5040.0)*q;
x67 = (1.0/720.0)*x45;
x68 = (1.0/240.0)*x30;
x69 = (1.0/144.0)*x18;
x70 = (1.0/144.0)*x25;
x71 = q*x10;
x72 = x14*x2;
x73 = x19*x7;
x74 = x1*x17;
x75 = (1.0/144.0)*x29;
x76 = (1.0/240.0)*x35;
x77 = (1.0/720.0)*x0;
x78 = x0*x14;
x79 = pow(y, 7);
x80 = (1.0/720.0)*x58;
x81 = (1.0/240.0)*x41;
x82 = (1.0/240.0)*x44;
x83 = (1.0/720.0)*x64;
x84 = pow(z, 7);
x85 = pow(x, 8);
x86 = (1.0/40320.0)*q;
x87 = (1.0/5040.0)*x65;
x88 = (1.0/1440.0)*x45;
x89 = x30*x46;
x90 = q*x25;
x91 = (1.0/576.0)*x18;
x92 = (1.0/96.0)*x21;
x93 = q*x29;
x94 = x10*x46;
x95 = x10*x2;
x96 = (1.0/72.0)*x10;
x97 = x14*x23;
x98 = x17*x21;
x99 = x1*x10;
x100 = (1.0/1440.0)*x35;
x101 = x23*x25;
x102 = x17*x35;
x103 = (1.0/5040.0)*x0;
x104 = pow(y, 8);
x105 = (1.0/5040.0)*x79;
x106 = (1.0/1440.0)*x58;
x107 = x17*x41;
x108 = (1.0/576.0)*x29;
x109 = x14*x44;
x110 = (1.0/1440.0)*x64;
x111 = (1.0/5040.0)*x84;
x112 = pow(z, 8);
x113 = (1.0/362880.0)*q;
x114 = (1.0/40320.0)*x85;
x115 = (1.0/10080.0)*x65;
x116 = (1.0/4320.0)*x45;
x117 = (1.0/2880.0)*x30;
x118 = (1.0/720.0)*x30;
x119 = (1.0/480.0)*x21;
x120 = (1.0/2880.0)*x41;
x121 = q*x18;
x122 = (1.0/288.0)*x18;
x123 = (1.0/2880.0)*x44;
x124 = (1.0/4320.0)*x71;
x125 = (1.0/288.0)*x10;
x126 = (1.0/10080.0)*x35;
x127 = (1.0/40320.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x6;
M[6] += x4*x7;
M[7] += x8;
M[8] += x4*x9;
M[9] += -x10*x11;
M[10] += -x1*x12;
M[11] += -x12*x2;
M[12] += -x13*x7;
M[13] += -x5*z;
M[14] += -x13*x9;
M[15] += -x11*x14;
M[16] += -x15*x2;
M[17] += -x1*x16;
M[18] += -x11*x17;
M[19] += x18*x19;
M[20] += x1*x20;
M[21] += x2*x20;
M[22] += x21*x22;
M[23] += x12*x8;
M[24] += x22*x23;
M[25] += x14*x24;
M[26] += x15*x6;
M[27] += x16*x5;
M[28] += x17*x24;
M[29] += x19*x25;
M[30] += x2*x26;
M[31] += x21*x27;
M[32] += x1*x28;
M[33] += x19*x29;
M[34] += -x30*x31;
M[35] += -x1*x32;
M[36] += -x2*x32;
M[37] += -x21*x33;
M[38] += -x20*x8;
M[39] += -x23*x33;
M[40] += -x34*x35;
M[41] += -x22*x36;
M[42] += -x22*x37;
M[43] += -x35*x38;
M[44] += -x25*x39;
M[45] += -x26*x6;
M[46] += -x27*x40;
M[47] += -x28*x5;
M[48] += -x29*x39;
M[49] += -x31*x41;
M[50] += -x2*x42;
M[51] += -x23*x34;
M[52] += -x21*x38;
M[53] += -x1*x43;
M[54] += -x31*x44;
M[55] += x45*x46;
M[56] += x1*x47;
M[57] += x2*x47;
M[58] += x21*x48;
M[59] += x32*x8;
M[60] += x23*x48;
M[61] += x49*x50;
M[62] += x33*x36;
M[63] += x33*x37;
M[64] += x50*x51;
M[65] += x25*x52;
M[66] += x34*x53;
M[67] += (1.0/8.0)*x21*x54;
M[68] += x38*x55;
M[69] += x29*x52;
M[70] += x41*x56;
M[71] += x42*x6;
M[72] += x34*x57;
M[73] += x38*x40;
M[74] += x43*x5;
M[75] += x44*x56;
M[76] += x46*x58;
M[77] += x2*x59;
M[78] += x23*x60;
M[79] += x49*x61;
M[80] += x21*x62;
M[81] += x1*x63;
M[82] += x46*x64;
M[83] += -x65*x66;
M[84] += -x1*x67;
M[85] += -x2*x67;
M[86] += -x21*x68;
M[87] += -x47*x8;
M[88] += -x23*x68;
M[89] += -x49*x69;
M[90] += -x36*x48;
M[91] += -x37*x48;
M[92] += -x51*x69;
M[93] += -x70*x71;
M[94] += -x50*x72;
M[95] += -x10*x73*x9;
M[96] += -x50*x74;
M[97] += -x71*x75;
M[98] += -x41*x76;
M[99] += -x53*x60;
M[100] += -x14*x19*x54;
M[101] += -x17*x3*x73;
M[102] += -x55*x62;
M[103] += -x44*x76;
M[104] += -x58*x77;
M[105] += -x59*x6;
M[106] += -x57*x60;
M[107] += -x61*x78;
M[108] += -x40*x62;
M[109] += -x5*x63;
M[110] += -x64*x77;
M[111] += -x66*x79;
M[112] += -x2*x80;
M[113] += -x23*x81;
M[114] += -x51*x70;
M[115] += -x49*x75;
M[116] += -x21*x82;
M[117] += -x1*x83;
M[118] += -x66*x84;
M[119] += x85*x86;
M[120] += x1*x87;
M[121] += x2*x87;
M[122] += x21*x88;
M[123] += x67*x8;
M[124] += x23*x88;
M[125] += x14*x89;
M[126] += x36*x68;
M[127] += x37*x68;
M[128] += x17*x89;
M[129] += x90*x91;
M[130] += x69*x72;
M[131] += x18*x9*x92;
M[132] += x69*x74;
M[133] += x91*x93;
M[134] += x41*x94;
M[135] += x70*x95;
M[136] += x96*x97;
M[137] += x96*x98;
M[138] += x75*x99;
M[139] += x44*x94;
M[140] += x100*x58;
M[141] += x53*x81;
M[142] += (1.0/96.0)*x101*x3;
M[143] += (1.0/72.0)*x102*x14;
M[144] += x29*x3*x92;
M[145] += x55*x82;
M[146] += x100*x64;
M[147] += x103*x79;
M[148] += x6*x80;
M[149] += x57*x81;
M[150] += x0*x17*x70;
M[151] += x75*x78;
M[152] += x40*x82;
M[153] += x5*x83;
M[154] += x103*x84;
M[155] += x104*x86;
M[156] += x105*x2;
M[157] += x106*x23;
M[158] += x107*x46;
M[159] += x108*x90;
M[160] += x109*x46;
M[161] += x110*x21;
M[162] += x1*x111;
M[163] += x112*x86;
M[164] += -pow(x, 9)*x113;
M[165] += -x1*x114;
M[166] += -x114*x2;
M[167] += -x115*x21;
M[168] += -x8*x87;
M[169] += -x115*x23;
M[170] += -x116*x49;
M[171] += -x36*x88;
M[172] += -x37*x88;
M[173] += -x116*x51;
M[174] += -x117*x90;
M[175] += -x118*x72;
M[176] += -x119*x30*x9;
M[177] += -x118*x74;
M[178] += -x117*x93;
M[179] += -x120*x121;
M[180] += -x2*x25*x91;
M[181] += -x122*x97;
M[182] += -x122*x98;
M[183] += -x1*x29*x91;
M[184] += -x121*x123;
M[185] += -x124*x58;
M[186] += -1.0/720.0*x41*x95;
M[187] += -x101*x125;
M[188] += -1.0/216.0*x10*x17*x49;
M[189] += -x125*x21*x29;
M[190] += -1.0/720.0*x44*x99;
M[191] += -x124*x64;
M[192] += -x126*x79;
M[193] += -x106*x53;
M[194] += -1.0/480.0*x23*x3*x41;
M[195] += -1.0/288.0*x102*x25;
M[196] += -1.0/288.0*x14*x29*x35;
M[197] += -x119*x3*x44;
M[198] += -x110*x55;
M[199] += -x126*x84;
M[200] += -x104*x127;
M[201] += -x105*x6;
M[202] += -x106*x57;
M[203] += -x107*x77;
M[204] += -x0*x108*x25;
M[205] += -x109*x77;
M[206] += -x110*x40;
M[207] += -x111*x5;
M[208] += -x112*x127;
M[209] += -x113*pow(y, 9);
M[210] += -1.0/40320.0*x104*x2;
M[211] += -1.0/10080.0*x23*x79;
M[212] += -1.0/4320.0*x51*x58;
M[213] += -x120*x93;
M[214] += -x123*x90;
M[215] += -1.0/4320.0*x49*x64;
M[216] += -1.0/10080.0*x21*x84;
M[217] += -1.0/40320.0*x1*x112;
M[218] += -x113*pow(z, 9);
}
void M2M_9(double x, double y, double z, double * M, double * Ms) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
double x448;
double x449;
double x450;
double x451;
double x452;
double x453;
double x454;
double x455;
double x456;
double x457;
double x458;
double x459;
double x460;
double x461;
double x462;
double x463;
double x464;
double x465;
double x466;
double x467;
double x468;
double x469;
double x470;
double x471;
double x472;
double x473;
double x474;
double x475;
double x476;
double x477;
double x478;
double x479;
double x480;
double x481;
double x482;
double x483;
double x484;
double x485;
double x486;
double x487;
double x488;
double x489;
double x490;
double x491;
double x492;
double x493;
double x494;
double x495;
double x496;
double x497;
double x498;
double x499;
double x500;
double x501;
double x502;
double x503;
double x504;
double x505;
double x506;
double x507;
double x508;
double x509;
double x510;
double x511;
double x512;
double x513;
double x514;
double x515;
double x516;
double x517;
double x518;
double x519;
double x520;
double x521;
double x522;
double x523;
double x524;
double x525;
double x526;
double x527;
double x528;
double x529;
double x530;
double x531;
double x532;
double x533;
double x534;
double x535;
double x536;
double x537;
double x538;
double x539;
double x540;
double x541;
double x542;
double x543;
double x544;
double x545;
double x546;
double x547;
double x548;
double x549;
double x550;
double x551;
double x552;
double x553;
double x554;
double x555;
double x556;
double x557;
double x558;
double x559;
double x560;
double x561;
double x562;
double x563;
double x564;
double x565;
double x566;
double x567;
double x568;
double x569;
double x570;
double x571;
double x572;
double x573;
double x574;
double x575;
double x576;
double x577;
double x578;
double x579;
double x580;
double x581;
double x582;
double x583;
double x584;
double x585;
double x586;
double x587;
double x588;
double x589;
double x590;
double x591;
double x592;
double x593;
double x594;
double x595;
double x596;
double x597;
double x598;
double x599;
double x600;
double x601;
double x602;
double x603;
double x604;
double x605;
double x606;
double x607;
double x608;
double x609;
double x610;
double x611;
double x612;
double x613;
double x614;
double x615;
double x616;
double x617;
double x618;
double x619;
double x620;
double x621;
double x622;
double x623;
double x624;
double x625;
double x626;
double x627;
double x628;
double x629;
double x630;
double x631;
double x632;
double x633;
double x634;
double x635;
double x636;
double x637;
double x638;
double x639;
double x640;
double x641;
double x642;
double x643;
double x644;
double x645;
double x646;
double x647;
double x648;
double x649;
double x650;
double x651;
double x652;
double x653;
double x654;
double x655;
double x656;
double x657;
double x658;
double x659;
double x660;
double x661;
double x662;
double x663;
double x664;
double x665;
double x666;
double x667;
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = z*M[0];
x5 = y*M[1];
x6 = y*M[2];
x7 = z*M[1];
x8 = z*M[2];
x9 = x*M[3];
x10 = (x*x);
x11 = (1.0/2.0)*x10;
x12 = x*M[4];
x13 = y*M[3];
x14 = x0*y;
x15 = x*M[5];
x16 = z*M[3];
x17 = x0*z;
x18 = x*M[6];
x19 = y*M[4];
x20 = x1*y;
x21 = (y*y);
x22 = (1.0/2.0)*M[0];
x23 = x*M[7];
x24 = y*M[5];
x25 = z*M[4];
x26 = x3*y;
x27 = x1*z;
x28 = x2*z;
x29 = x*M[8];
x30 = z*M[5];
x31 = x3*z;
x32 = (z*z);
x33 = y*M[6];
x34 = (1.0/2.0)*x21;
x35 = y*M[7];
x36 = z*M[6];
x37 = x5*z;
x38 = y*M[8];
x39 = z*M[7];
x40 = x6*z;
x41 = (1.0/2.0)*x32;
x42 = z*M[8];
x43 = x*M[9];
x44 = (x*x*x);
x45 = (1.0/6.0)*x44;
x46 = x*M[10];
x47 = y*M[9];
x48 = x9*y;
x49 = x*M[11];
x50 = z*M[9];
x51 = x9*z;
x52 = x*M[12];
x53 = y*M[10];
x54 = x12*y;
x55 = x*M[13];
x56 = y*M[11];
x57 = z*M[10];
x58 = x15*y;
x59 = x12*z;
x60 = x13*z;
x61 = x*M[14];
x62 = z*M[11];
x63 = x15*z;
x64 = x*M[15];
x65 = y*M[12];
x66 = x18*y;
x67 = (y*y*y);
x68 = (1.0/6.0)*M[0];
x69 = x*M[16];
x70 = y*M[13];
x71 = z*M[12];
x72 = x23*y;
x73 = x18*z;
x74 = x19*z;
x75 = x*M[17];
x76 = y*M[14];
x77 = z*M[13];
x78 = x29*y;
x79 = x23*z;
x80 = x24*z;
x81 = x*M[18];
x82 = z*M[14];
x83 = x29*z;
x84 = (z*z*z);
x85 = y*M[15];
x86 = (1.0/6.0)*x67;
x87 = y*M[16];
x88 = z*M[15];
x89 = x33*z;
x90 = y*M[17];
x91 = z*M[16];
x92 = x35*z;
x93 = y*M[18];
x94 = z*M[17];
x95 = x38*z;
x96 = (1.0/6.0)*x84;
x97 = z*M[18];
x98 = x*M[19];
x99 = (x*x*x*x);
x100 = (1.0/24.0)*x99;
x101 = x*M[20];
x102 = y*M[19];
x103 = x43*y;
x104 = x*M[21];
x105 = z*M[19];
x106 = x43*z;
x107 = x*M[22];
x108 = y*M[20];
x109 = x46*y;
x110 = (1.0/4.0)*x10;
x111 = x21*M[0];
x112 = x*M[23];
x113 = y*M[21];
x114 = z*M[20];
x115 = x49*y;
x116 = x46*z;
x117 = x47*z;
x118 = x*M[24];
x119 = z*M[21];
x120 = x49*z;
x121 = x110*x32;
x122 = x*M[25];
x123 = y*M[22];
x124 = x52*y;
x125 = x110*x21;
x126 = x*M[26];
x127 = y*M[23];
x128 = z*M[22];
x129 = x55*y;
x130 = x52*z;
x131 = x53*z;
x132 = x*M[27];
x133 = y*M[24];
x134 = z*M[23];
x135 = x61*y;
x136 = x55*z;
x137 = x56*z;
x138 = x*M[28];
x139 = z*M[24];
x140 = x61*z;
x141 = x*M[29];
x142 = y*M[25];
x143 = x64*y;
x144 = (y*y*y*y);
x145 = (1.0/24.0)*M[0];
x146 = x*M[30];
x147 = y*M[26];
x148 = z*M[25];
x149 = x69*y;
x150 = x64*z;
x151 = x65*z;
x152 = x*M[31];
x153 = y*M[27];
x154 = z*M[26];
x155 = x75*y;
x156 = x69*z;
x157 = x70*z;
x158 = (1.0/4.0)*x32;
x159 = x*M[32];
x160 = y*M[28];
x161 = z*M[27];
x162 = x81*y;
x163 = x75*z;
x164 = x76*z;
x165 = x*M[33];
x166 = z*M[28];
x167 = x81*z;
x168 = (z*z*z*z);
x169 = y*M[29];
x170 = (1.0/24.0)*x144;
x171 = y*M[30];
x172 = z*M[29];
x173 = x85*z;
x174 = y*M[31];
x175 = z*M[30];
x176 = x87*z;
x177 = x158*x21;
x178 = y*M[32];
x179 = z*M[31];
x180 = x90*z;
x181 = y*M[33];
x182 = z*M[32];
x183 = x93*z;
x184 = (1.0/24.0)*x168;
x185 = z*M[33];
x186 = x*M[34];
x187 = pow(x, 5);
x188 = (1.0/120.0)*x187;
x189 = x*M[35];
x190 = y*M[34];
x191 = x98*y;
x192 = x*M[36];
x193 = z*M[34];
x194 = x98*z;
x195 = x*M[37];
x196 = y*M[35];
x197 = x101*y;
x198 = (1.0/12.0)*x44;
x199 = x*M[38];
x200 = y*M[36];
x201 = z*M[35];
x202 = x104*y;
x203 = x101*z;
x204 = x102*z;
x205 = x*M[39];
x206 = z*M[36];
x207 = x104*z;
x208 = x198*x32;
x209 = x*M[40];
x210 = y*M[37];
x211 = x107*y;
x212 = (1.0/12.0)*x10;
x213 = x67*M[0];
x214 = x198*x21;
x215 = x*M[41];
x216 = y*M[38];
x217 = z*M[37];
x218 = x112*y;
x219 = x107*z;
x220 = x108*z;
x221 = x*M[42];
x222 = y*M[39];
x223 = z*M[38];
x224 = x118*y;
x225 = x112*z;
x226 = x113*z;
x227 = x*M[43];
x228 = z*M[39];
x229 = x118*z;
x230 = x212*x84;
x231 = x*M[44];
x232 = y*M[40];
x233 = x122*y;
x234 = x212*x67;
x235 = x*M[45];
x236 = y*M[41];
x237 = z*M[40];
x238 = x126*y;
x239 = x122*z;
x240 = x123*z;
x241 = x*M[46];
x242 = y*M[42];
x243 = z*M[41];
x244 = x132*y;
x245 = x126*z;
x246 = x127*z;
x247 = x*M[47];
x248 = y*M[43];
x249 = z*M[42];
x250 = x138*y;
x251 = x132*z;
x252 = x133*z;
x253 = x*M[48];
x254 = z*M[43];
x255 = x138*z;
x256 = x*M[49];
x257 = y*M[44];
x258 = x141*y;
x259 = pow(y, 5);
x260 = (1.0/120.0)*M[0];
x261 = x*M[50];
x262 = y*M[45];
x263 = z*M[44];
x264 = x146*y;
x265 = x141*z;
x266 = x142*z;
x267 = x*M[51];
x268 = y*M[46];
x269 = z*M[45];
x270 = x152*y;
x271 = x146*z;
x272 = x147*z;
x273 = (1.0/12.0)*x32;
x274 = x*M[52];
x275 = y*M[47];
x276 = z*M[46];
x277 = x159*y;
x278 = x152*z;
x279 = x153*z;
x280 = (1.0/12.0)*x84;
x281 = x*M[53];
x282 = y*M[48];
x283 = z*M[47];
x284 = x165*y;
x285 = x159*z;
x286 = x160*z;
x287 = x*M[54];
x288 = z*M[48];
x289 = x165*z;
x290 = pow(z, 5);
x291 = y*M[49];
x292 = (1.0/120.0)*x259;
x293 = y*M[50];
x294 = z*M[49];
x295 = x169*z;
x296 = y*M[51];
x297 = z*M[50];
x298 = x171*z;
x299 = x273*x67;
x300 = y*M[52];
x301 = z*M[51];
x302 = x174*z;
x303 = x21*x280;
x304 = y*M[53];
x305 = z*M[52];
x306 = x178*z;
x307 = y*M[54];
x308 = z*M[53];
x309 = x181*z;
x310 = (1.0/120.0)*x290;
x311 = z*M[54];
x312 = x*M[55];
x313 = pow(x, 6);
x314 = (1.0/720.0)*x313;
x315 = x*M[56];
x316 = y*M[55];
x317 = x186*y;
x318 = x*M[57];
x319 = z*M[55];
x320 = x186*z;
x321 = x*M[58];
x322 = y*M[56];
x323 = x189*y;
x324 = (1.0/48.0)*x99;
x325 = x*M[59];
x326 = y*M[57];
x327 = z*M[56];
x328 = x192*y;
x329 = x189*z;
x330 = x190*z;
x331 = x*M[60];
x332 = z*M[57];
x333 = x192*z;
x334 = x32*x324;
x335 = x*M[61];
x336 = y*M[58];
x337 = x195*y;
x338 = (1.0/36.0)*x44;
x339 = x21*x324;
x340 = x*M[62];
x341 = y*M[59];
x342 = z*M[58];
x343 = x199*y;
x344 = x195*z;
x345 = x196*z;
x346 = x*M[63];
x347 = y*M[60];
x348 = z*M[59];
x349 = x205*y;
x350 = x199*z;
x351 = x200*z;
x352 = x*M[64];
x353 = z*M[60];
x354 = x205*z;
x355 = x338*x84;
x356 = x*M[65];
x357 = y*M[61];
x358 = x209*y;
x359 = (1.0/48.0)*x10;
x360 = x144*M[0];
x361 = x338*x67;
x362 = x*M[66];
x363 = y*M[62];
x364 = z*M[61];
x365 = x215*y;
x366 = x209*z;
x367 = x210*z;
x368 = x*M[67];
x369 = y*M[63];
x370 = z*M[62];
x371 = x221*y;
x372 = x215*z;
x373 = x216*z;
x374 = x10*x32;
x375 = (1.0/8.0)*x374;
x376 = x*M[68];
x377 = y*M[64];
x378 = z*M[63];
x379 = x227*y;
x380 = x221*z;
x381 = x222*z;
x382 = x*M[69];
x383 = z*M[64];
x384 = x227*z;
x385 = x168*x359;
x386 = x*M[70];
x387 = y*M[65];
x388 = x231*y;
x389 = x144*x359;
x390 = x*M[71];
x391 = y*M[66];
x392 = z*M[65];
x393 = x235*y;
x394 = x231*z;
x395 = x232*z;
x396 = x*M[72];
x397 = y*M[67];
x398 = z*M[66];
x399 = x241*y;
x400 = x235*z;
x401 = x236*z;
x402 = x21*x375;
x403 = x*M[73];
x404 = y*M[68];
x405 = z*M[67];
x406 = x247*y;
x407 = x241*z;
x408 = x242*z;
x409 = x*M[74];
x410 = y*M[69];
x411 = z*M[68];
x412 = x253*y;
x413 = x247*z;
x414 = x248*z;
x415 = x*M[75];
x416 = z*M[69];
x417 = x253*z;
x418 = x*M[76];
x419 = y*M[70];
x420 = x256*y;
x421 = pow(y, 6);
x422 = (1.0/720.0)*M[0];
x423 = x*M[77];
x424 = y*M[71];
x425 = z*M[70];
x426 = x261*y;
x427 = x256*z;
x428 = x257*z;
x429 = x*M[78];
x430 = y*M[72];
x431 = z*M[71];
x432 = x267*y;
x433 = x261*z;
x434 = x262*z;
x435 = (1.0/48.0)*x32;
x436 = x*M[79];
x437 = y*M[73];
x438 = z*M[72];
x439 = x274*y;
x440 = x267*z;
x441 = x268*z;
x442 = (1.0/36.0)*x84;
x443 = x*M[80];
x444 = y*M[74];
x445 = z*M[73];
x446 = x281*y;
x447 = x274*z;
x448 = x275*z;
x449 = (1.0/48.0)*x168;
x450 = x*M[81];
x451 = y*M[75];
x452 = z*M[74];
x453 = x287*y;
x454 = x281*z;
x455 = x282*z;
x456 = x*M[82];
x457 = z*M[75];
x458 = x287*z;
x459 = pow(z, 6);
x460 = y*M[76];
x461 = (1.0/720.0)*x421;
x462 = y*M[77];
x463 = z*M[76];
x464 = x291*z;
x465 = y*M[78];
x466 = z*M[77];
x467 = x293*z;
x468 = x144*x435;
x469 = y*M[79];
x470 = z*M[78];
x471 = x296*z;
x472 = x442*x67;
x473 = y*M[80];
x474 = z*M[79];
x475 = x300*z;
x476 = x21*x449;
x477 = y*M[81];
x478 = z*M[80];
x479 = x304*z;
x480 = y*M[82];
x481 = z*M[81];
x482 = x307*z;
x483 = (1.0/720.0)*x459;
x484 = z*M[82];
x485 = x*M[83];
x486 = (1.0/5040.0)*pow(x, 7);
x487 = x*M[84];
x488 = y*M[83];
x489 = x312*y;
x490 = x*M[85];
x491 = x*M[86];
x492 = y*M[84];
x493 = x315*y;
x494 = (1.0/240.0)*x187;
x495 = x*M[87];
x496 = y*M[85];
x497 = x318*y;
x498 = x*M[88];
x499 = x32*x494;
x500 = x*M[89];
x501 = y*M[86];
x502 = x321*y;
x503 = (1.0/144.0)*x99;
x504 = x21*x494;
x505 = x*M[90];
x506 = y*M[87];
x507 = x325*y;
x508 = x*M[91];
x509 = y*M[88];
x510 = x331*y;
x511 = x*M[92];
x512 = x503*x84;
x513 = x*M[93];
x514 = y*M[89];
x515 = x335*y;
x516 = (1.0/144.0)*x44;
x517 = x503*x67;
x518 = x*M[94];
x519 = y*M[90];
x520 = x340*y;
x521 = x*M[95];
x522 = y*M[91];
x523 = x346*y;
x524 = x145*x21;
x525 = x32*x44;
x526 = x*M[96];
x527 = y*M[92];
x528 = x352*y;
x529 = x*M[97];
x530 = x168*x516;
x531 = x*M[98];
x532 = y*M[93];
x533 = x356*y;
x534 = (1.0/240.0)*x10;
x535 = x534*M[0];
x536 = x144*x516;
x537 = x*M[99];
x538 = y*M[94];
x539 = x362*y;
x540 = x*M[100];
x541 = y*M[95];
x542 = x368*y;
x543 = x374*x67;
x544 = (1.0/24.0)*x21;
x545 = x525*x544;
x546 = x*M[101];
x547 = y*M[96];
x548 = x376*y;
x549 = x10*x84;
x550 = x*M[102];
x551 = y*M[97];
x552 = x382*y;
x553 = x*M[103];
x554 = x*M[104];
x555 = y*M[98];
x556 = x386*y;
x557 = x259*x534;
x558 = x*M[105];
x559 = y*M[99];
x560 = x390*y;
x561 = x*M[106];
x562 = y*M[100];
x563 = x396*y;
x564 = x67*M[1];
x565 = x*M[107];
x566 = y*M[101];
x567 = x403*y;
x568 = x84*M[1];
x569 = (1.0/24.0)*x543;
x570 = x*M[108];
x571 = y*M[102];
x572 = x409*y;
x573 = x544*x549;
x574 = x*M[109];
x575 = y*M[103];
x576 = x415*y;
x577 = x290*x534;
x578 = x*M[110];
x579 = x*M[111];
x580 = y*M[104];
x581 = x418*y;
x582 = pow(y, 7);
x583 = (1.0/5040.0)*M[0];
x584 = x*M[112];
x585 = y*M[105];
x586 = x423*y;
x587 = x*M[113];
x588 = y*M[106];
x589 = x429*y;
x590 = x32*M[0];
x591 = (1.0/240.0)*x259;
x592 = x*M[114];
x593 = y*M[107];
x594 = x436*y;
x595 = (1.0/144.0)*x84;
x596 = x*M[115];
x597 = y*M[108];
x598 = x443*y;
x599 = (1.0/144.0)*x168;
x600 = x*M[116];
x601 = y*M[109];
x602 = x450*y;
x603 = (1.0/240.0)*x290;
x604 = x*M[117];
x605 = y*M[110];
x606 = x456*y;
x607 = x*M[118];
x608 = pow(z, 7);
x609 = y*M[111];
x610 = (1.0/5040.0)*x582;
x611 = y*M[112];
x612 = y*M[113];
x613 = x32*x591;
x614 = y*M[114];
x615 = x144*x595;
x616 = y*M[115];
x617 = x599*x67;
x618 = y*M[116];
x619 = x21*x603;
x620 = y*M[117];
x621 = y*M[118];
x622 = (1.0/5040.0)*x608;
x623 = (1.0/40320.0)*pow(x, 8);
x624 = (1.0/1440.0)*x313;
x625 = x187*x422;
x626 = x21*x624;
x627 = x32*x624;
x628 = (1.0/576.0)*x99;
x629 = (1.0/720.0)*x187;
x630 = x629*M[2];
x631 = (1.0/96.0)*x111;
x632 = x32*x99;
x633 = x168*x628;
x634 = x422*x44;
x635 = x144*x628;
x636 = (1.0/72.0)*x525;
x637 = (1.0/96.0)*M[1];
x638 = x21*x632;
x639 = (1.0/72.0)*x44;
x640 = x639*x84;
x641 = (1.0/96.0)*M[2];
x642 = (1.0/1440.0)*x10;
x643 = x642*M[0];
x644 = (1.0/720.0)*x44;
x645 = x259*x644;
x646 = (1.0/72.0)*x549;
x647 = x67*M[2];
x648 = x10*x168;
x649 = x290*x644;
x650 = x421*x642;
x651 = x144*x374;
x652 = x21*x648;
x653 = x459*x642;
x654 = pow(y, 8);
x655 = (1.0/40320.0)*M[0];
x656 = (1.0/1440.0)*x421;
x657 = x259*x84;
x658 = (1.0/576.0)*x168;
x659 = x290*x67;
x660 = (1.0/1440.0)*x459;
x661 = pow(z, 8);
x662 = (1.0/40320.0)*x654;
x663 = x32*x656;
x664 = x144*x658;
x665 = (1.0/720.0)*M[2];
x666 = x21*x660;
x667 = (1.0/40320.0)*x661;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + x4 + M[5];
#pragma omp atomic
Ms[6] += x5 + M[6];
#pragma omp atomic
Ms[7] += x6 + x7 + M[7];
#pragma omp atomic
Ms[8] += x8 + M[8];
#pragma omp atomic
Ms[9] += x11*M[0] + x9 + M[9];
#pragma omp atomic
Ms[10] += x11*M[1] + x12 + x13 + x14 + M[10];
#pragma omp atomic
Ms[11] += x11*M[2] + x15 + x16 + x17 + M[11];
#pragma omp atomic
Ms[12] += x18 + x19 + x20 + x21*x22 + M[12];
#pragma omp atomic
Ms[13] += x23 + x24 + x25 + x26 + x27 + x28 + M[13];
#pragma omp atomic
Ms[14] += x22*x32 + x29 + x30 + x31 + M[14];
#pragma omp atomic
Ms[15] += x33 + x34*M[1] + M[15];
#pragma omp atomic
Ms[16] += x34*M[2] + x35 + x36 + x37 + M[16];
#pragma omp atomic
Ms[17] += x38 + x39 + x40 + x41*M[1] + M[17];
#pragma omp atomic
Ms[18] += x41*M[2] + x42 + M[18];
#pragma omp atomic
Ms[19] += x11*M[3] + x43 + x45*M[0] + M[19];
#pragma omp atomic
Ms[20] += x11*x2 + x11*M[4] + x45*M[1] + x46 + x47 + x48 + M[20];
#pragma omp atomic
Ms[21] += x11*x4 + x11*M[5] + x45*M[2] + x49 + x50 + x51 + M[21];
#pragma omp atomic
Ms[22] += x0*x34 + x11*x5 + x11*M[6] + x34*M[3] + x52 + x53 + x54 + M[22];
#pragma omp atomic
Ms[23] += x11*x6 + x11*x7 + x11*M[7] + x14*z + x55 + x56 + x57 + x58 + x59 + x60 + M[23];
#pragma omp atomic
Ms[24] += x0*x41 + x11*x8 + x11*M[8] + x41*M[3] + x61 + x62 + x63 + M[24];
#pragma omp atomic
Ms[25] += x1*x34 + x34*M[4] + x64 + x65 + x66 + x67*x68 + M[25];
#pragma omp atomic
Ms[26] += x20*z + x3*x34 + x34*x4 + x34*M[5] + x69 + x70 + x71 + x72 + x73 + x74 + M[26];
#pragma omp atomic
Ms[27] += x1*x41 + x2*x41 + x26*z + x41*M[4] + x75 + x76 + x77 + x78 + x79 + x80 + M[27];
#pragma omp atomic
Ms[28] += x3*x41 + x41*M[5] + x68*x84 + x81 + x82 + x83 + M[28];
#pragma omp atomic
Ms[29] += x34*M[6] + x85 + x86*M[1] + M[29];
#pragma omp atomic
Ms[30] += x34*x7 + x34*M[7] + x86*M[2] + x87 + x88 + x89 + M[30];
#pragma omp atomic
Ms[31] += x34*x8 + x34*M[8] + x41*x5 + x41*M[6] + x90 + x91 + x92 + M[31];
#pragma omp atomic
Ms[32] += x41*x6 + x41*M[7] + x93 + x94 + x95 + x96*M[1] + M[32];
#pragma omp atomic
Ms[33] += x41*M[8] + x96*M[2] + x97 + M[33];
#pragma omp atomic
Ms[34] += x100*M[0] + x11*M[9] + x45*M[3] + x98 + M[34];
#pragma omp atomic
Ms[35] += x100*M[1] + x101 + x102 + x103 + x11*x13 + x11*M[10] + x2*x45 + x45*M[4] + M[35];
#pragma omp atomic
Ms[36] += x100*M[2] + x104 + x105 + x106 + x11*x16 + x11*M[11] + x4*x45 + x45*M[5] + M[36];
#pragma omp atomic
Ms[37] += x107 + x108 + x109 + x11*x19 + x11*M[12] + x110*x111 + x34*x9 + x34*M[9] + x45*x5 + x45*M[6] + M[37];
#pragma omp atomic
Ms[38] += x11*x24 + x11*x25 + x11*x28 + x11*M[13] + x112 + x113 + x114 + x115 + x116 + x117 + x45*x6 + x45*x7 + x45*M[7] + x48*z + M[38];
#pragma omp atomic
Ms[39] += x11*x30 + x11*M[14] + x118 + x119 + x120 + x121*M[0] + x41*x9 + x41*M[9] + x45*x8 + x45*M[8] + M[39];
#pragma omp atomic
Ms[40] += x0*x86 + x11*x33 + x11*M[15] + x12*x34 + x122 + x123 + x124 + x125*M[1] + x34*M[10] + x86*M[3] + M[40];
#pragma omp atomic
Ms[41] += x11*x35 + x11*x36 + x11*x37 + x11*M[16] + x125*M[2] + x126 + x127 + x128 + x129 + x130 + x131 + x15*x34 + x16*x34 + x17*x34 + x34*M[11] + x54*z + M[41];
#pragma omp atomic
Ms[42] += x11*x38 + x11*x39 + x11*x40 + x11*M[17] + x12*x41 + x121*M[1] + x13*x41 + x132 + x133 + x134 + x135 + x136 + x137 + x14*x41 + x41*M[10] + x58*z + M[42];
#pragma omp atomic
Ms[43] += x0*x96 + x11*x42 + x11*M[18] + x121*M[2] + x138 + x139 + x140 + x15*x41 + x41*M[11] + x96*M[3] + M[43];
#pragma omp atomic
Ms[44] += x1*x86 + x141 + x142 + x143 + x144*x145 + x18*x34 + x34*M[12] + x86*M[4] + M[44];
#pragma omp atomic
Ms[45] += x146 + x147 + x148 + x149 + x150 + x151 + x23*x34 + x25*x34 + x27*x34 + x3*x86 + x34*M[13] + x4*x86 + x66*z + x86*M[5] + M[45];
#pragma omp atomic
Ms[46] += x111*x158 + x152 + x153 + x154 + x155 + x156 + x157 + x18*x41 + x19*x41 + x20*x41 + x29*x34 + x30*x34 + x31*x34 + x34*M[14] + x41*M[12] + x72*z + M[46];
#pragma omp atomic
Ms[47] += x1*x96 + x159 + x160 + x161 + x162 + x163 + x164 + x2*x96 + x23*x41 + x24*x41 + x26*x41 + x41*M[13] + x78*z + x96*M[4] + M[47];
#pragma omp atomic
Ms[48] += x145*x168 + x165 + x166 + x167 + x29*x41 + x3*x96 + x41*M[14] + x96*M[5] + M[48];
#pragma omp atomic
Ms[49] += x169 + x170*M[1] + x34*M[15] + x86*M[6] + M[49];
#pragma omp atomic
Ms[50] += x170*M[2] + x171 + x172 + x173 + x34*x36 + x34*M[16] + x7*x86 + x86*M[7] + M[50];
#pragma omp atomic
Ms[51] += x174 + x175 + x176 + x177*M[1] + x33*x41 + x34*x39 + x34*M[17] + x41*M[15] + x8*x86 + x86*M[8] + M[51];
#pragma omp atomic
Ms[52] += x177*M[2] + x178 + x179 + x180 + x34*x42 + x34*M[18] + x35*x41 + x41*M[16] + x5*x96 + x96*M[6] + M[52];
#pragma omp atomic
Ms[53] += x181 + x182 + x183 + x184*M[1] + x38*x41 + x41*M[17] + x6*x96 + x96*M[7] + M[53];
#pragma omp atomic
Ms[54] += x184*M[2] + x185 + x41*M[18] + x96*M[8] + M[54];
#pragma omp atomic
Ms[55] += x100*M[3] + x11*M[19] + x186 + x188*M[0] + x45*M[9] + M[55];
#pragma omp atomic
Ms[56] += x100*x2 + x100*M[4] + x11*x47 + x11*M[20] + x13*x45 + x188*M[1] + x189 + x190 + x191 + x45*M[10] + M[56];
#pragma omp atomic
Ms[57] += x100*x4 + x100*M[5] + x11*x50 + x11*M[21] + x16*x45 + x188*M[2] + x192 + x193 + x194 + x45*M[11] + M[57];
#pragma omp atomic
Ms[58] += x100*x5 + x100*M[6] + x11*x53 + x11*M[22] + x111*x198 + x125*M[3] + x19*x45 + x195 + x196 + x197 + x34*x43 + x34*M[19] + x45*M[12] + M[58];
#pragma omp atomic
Ms[59] += x100*x6 + x100*x7 + x100*M[7] + x103*z + x11*x56 + x11*x57 + x11*x60 + x11*M[23] + x199 + x200 + x201 + x202 + x203 + x204 + x24*x45 + x25*x45 + x28*x45 + x45*M[13] + M[59];
#pragma omp atomic
Ms[60] += x100*x8 + x100*M[8] + x11*x62 + x11*M[24] + x121*M[3] + x205 + x206 + x207 + x208*M[0] + x30*x45 + x41*x43 + x41*M[19] + x45*M[14] + M[60];
#pragma omp atomic
Ms[61] += x11*x65 + x11*M[25] + x125*M[4] + x209 + x210 + x211 + x212*x213 + x214*M[1] + x33*x45 + x34*x46 + x34*M[20] + x45*M[15] + x86*x9 + x86*M[9] + M[61];
#pragma omp atomic
Ms[62] += x109*z + x11*x70 + x11*x71 + x11*x74 + x11*M[26] + x125*x4 + x125*M[5] + x214*M[2] + x215 + x216 + x217 + x218 + x219 + x220 + x34*x49 + x34*x50 + x34*x51 + x34*M[21] + x35*x45 + x36*x45 + x37*x45 + x45*M[16] + M[62];
#pragma omp atomic
Ms[63] += x11*x76 + x11*x77 + x11*x80 + x11*M[27] + x115*z + x121*x2 + x121*M[4] + x208*M[1] + x221 + x222 + x223 + x224 + x225 + x226 + x38*x45 + x39*x45 + x40*x45 + x41*x46 + x41*x47 + x41*x48 + x41*M[20] + x45*M[17] + M[63];
#pragma omp atomic
Ms[64] += x11*x82 + x11*M[28] + x121*M[5] + x208*M[2] + x227 + x228 + x229 + x230*M[0] + x41*x49 + x41*M[21] + x42*x45 + x45*M[18] + x9*x96 + x96*M[9] + M[64];
#pragma omp atomic
Ms[65] += x0*x170 + x11*x85 + x11*M[29] + x12*x86 + x125*M[6] + x170*M[3] + x231 + x232 + x233 + x234*M[1] + x34*x52 + x34*M[22] + x86*M[10] + M[65];
#pragma omp atomic
Ms[66] += x11*x87 + x11*x88 + x11*x89 + x11*M[30] + x124*z + x125*x7 + x125*M[7] + x15*x86 + x16*x86 + x17*x86 + x234*M[2] + x235 + x236 + x237 + x238 + x239 + x240 + x34*x55 + x34*x57 + x34*x59 + x34*M[23] + x86*M[11] + M[66];
#pragma omp atomic
Ms[67] += x0*x177 + x11*x90 + x11*x91 + x11*x92 + x11*M[31] + x121*x5 + x121*M[6] + x125*x8 + x125*M[8] + x129*z + x177*M[3] + x241 + x242 + x243 + x244 + x245 + x246 + x34*x61 + x34*x62 + x34*x63 + x34*M[24] + x41*x52 + x41*x53 + x41*x54 + x41*M[22] + M[67];
#pragma omp atomic
Ms[68] += x11*x93 + x11*x94 + x11*x95 + x11*M[32] + x12*x96 + x121*x6 + x121*M[7] + x13*x96 + x135*z + x14*x96 + x230*M[1] + x247 + x248 + x249 + x250 + x251 + x252 + x41*x55 + x41*x56 + x41*x58 + x41*M[23] + x96*M[10] + M[68];
#pragma omp atomic
Ms[69] += x0*x184 + x11*x97 + x11*M[33] + x121*M[8] + x15*x96 + x184*M[3] + x230*M[2] + x253 + x254 + x255 + x41*x61 + x41*M[24] + x96*M[11] + M[69];
#pragma omp atomic
Ms[70] += x1*x170 + x170*M[4] + x18*x86 + x256 + x257 + x258 + x259*x260 + x34*x64 + x34*M[25] + x86*M[12] + M[70];
#pragma omp atomic
Ms[71] += x143*z + x170*x3 + x170*x4 + x170*M[5] + x23*x86 + x25*x86 + x261 + x262 + x263 + x264 + x265 + x266 + x27*x86 + x34*x69 + x34*x71 + x34*x73 + x34*M[26] + x86*M[13] + M[71];
#pragma omp atomic
Ms[72] += x1*x177 + x149*z + x177*M[4] + x213*x273 + x267 + x268 + x269 + x270 + x271 + x272 + x29*x86 + x30*x86 + x31*x86 + x34*x75 + x34*x77 + x34*x79 + x34*M[27] + x41*x64 + x41*x65 + x41*x66 + x41*M[25] + x86*M[14] + M[72];
#pragma omp atomic
Ms[73] += x111*x280 + x155*z + x177*x3 + x177*M[5] + x18*x96 + x19*x96 + x20*x96 + x274 + x275 + x276 + x277 + x278 + x279 + x34*x81 + x34*x82 + x34*x83 + x34*M[28] + x41*x69 + x41*x70 + x41*x72 + x41*M[26] + x96*M[12] + M[73];
#pragma omp atomic
Ms[74] += x1*x184 + x162*z + x184*x2 + x184*M[4] + x23*x96 + x24*x96 + x26*x96 + x281 + x282 + x283 + x284 + x285 + x286 + x41*x75 + x41*x76 + x41*x78 + x41*M[27] + x96*M[13] + M[74];
#pragma omp atomic
Ms[75] += x184*x3 + x184*M[5] + x260*x290 + x287 + x288 + x289 + x29*x96 + x41*x81 + x41*M[28] + x96*M[14] + M[75];
#pragma omp atomic
Ms[76] += x170*M[6] + x291 + x292*M[1] + x34*M[29] + x86*M[15] + M[76];
#pragma omp atomic
Ms[77] += x170*x7 + x170*M[7] + x292*M[2] + x293 + x294 + x295 + x34*x88 + x34*M[30] + x36*x86 + x86*M[16] + M[77];
#pragma omp atomic
Ms[78] += x170*x8 + x170*M[8] + x177*M[6] + x296 + x297 + x298 + x299*M[1] + x34*x91 + x34*M[31] + x39*x86 + x41*x85 + x41*M[29] + x86*M[17] + M[78];
#pragma omp atomic
Ms[79] += x177*M[7] + x299*M[2] + x300 + x301 + x302 + x303*M[1] + x33*x96 + x34*x94 + x34*M[32] + x41*x87 + x41*M[30] + x42*x86 + x86*M[18] + x96*M[15] + M[79];
#pragma omp atomic
Ms[80] += x177*M[8] + x184*x5 + x184*M[6] + x303*M[2] + x304 + x305 + x306 + x34*x97 + x34*M[33] + x35*x96 + x41*x90 + x41*M[31] + x96*M[16] + M[80];
#pragma omp atomic
Ms[81] += x184*x6 + x184*M[7] + x307 + x308 + x309 + x310*M[1] + x38*x96 + x41*x93 + x41*M[32] + x96*M[17] + M[81];
#pragma omp atomic
Ms[82] += x184*M[8] + x310*M[2] + x311 + x41*M[33] + x96*M[18] + M[82];
#pragma omp atomic
Ms[83] += x100*M[9] + x11*M[34] + x188*M[3] + x312 + x314*M[0] + x45*M[19] + M[83];
#pragma omp atomic
Ms[84] += x100*x13 + x100*M[10] + x102*x11 + x11*M[35] + x188*x2 + x188*M[4] + x314*M[1] + x315 + x316 + x317 + x45*x47 + x45*M[20] + M[84];
#pragma omp atomic
Ms[85] += x100*x16 + x100*M[11] + x105*x11 + x11*M[36] + x188*x4 + x188*M[5] + x314*M[2] + x318 + x319 + x320 + x45*x50 + x45*M[21] + M[85];
#pragma omp atomic
Ms[86] += x100*x19 + x100*M[12] + x108*x11 + x11*M[37] + x111*x324 + x125*M[9] + x188*x5 + x188*M[6] + x214*M[3] + x321 + x322 + x323 + x34*x98 + x34*M[34] + x45*x53 + x45*M[22] + M[86];
#pragma omp atomic
Ms[87] += x100*x24 + x100*x25 + x100*x28 + x100*M[13] + x11*x113 + x11*x114 + x11*x117 + x11*M[38] + x188*x6 + x188*x7 + x188*M[7] + x191*z + x325 + x326 + x327 + x328 + x329 + x330 + x45*x56 + x45*x57 + x45*x60 + x45*M[23] + M[87];
#pragma omp atomic
Ms[88] += x100*x30 + x100*M[14] + x11*x119 + x11*M[39] + x121*M[9] + x188*x8 + x188*M[8] + x208*M[3] + x331 + x332 + x333 + x334*M[0] + x41*x98 + x41*M[34] + x45*x62 + x45*M[24] + M[88];
#pragma omp atomic
Ms[89] += x100*x33 + x100*M[15] + x101*x34 + x11*x123 + x11*M[40] + x125*M[10] + x213*x338 + x214*M[4] + x234*M[3] + x335 + x336 + x337 + x339*M[1] + x34*M[35] + x43*x86 + x45*x65 + x45*M[25] + x86*M[19] + M[89];
#pragma omp atomic
Ms[90] += x100*x35 + x100*x36 + x100*x37 + x100*M[16] + x104*x34 + x105*x34 + x106*x34 + x11*x127 + x11*x128 + x11*x131 + x11*M[41] + x125*x16 + x125*M[11] + x197*z + x214*x4 + x214*M[5] + x339*M[2] + x34*M[36] + x340 + x341 + x342 + x343 + x344 + x345 + x45*x70 + x45*x71 + x45*x74 + x45*M[26] + M[90];
#pragma omp atomic
Ms[91] += x100*x38 + x100*x39 + x100*x40 + x100*M[17] + x101*x41 + x102*x41 + x103*x41 + x11*x133 + x11*x134 + x11*x137 + x11*M[42] + x121*x13 + x121*M[10] + x2*x208 + x202*z + x208*M[4] + x334*M[1] + x346 + x347 + x348 + x349 + x350 + x351 + x41*M[35] + x45*x76 + x45*x77 + x45*x80 + x45*M[27] + M[91];
#pragma omp atomic
Ms[92] += x100*x42 + x100*M[18] + x104*x41 + x11*x139 + x11*M[43] + x121*M[11] + x208*M[5] + x230*M[3] + x334*M[2] + x352 + x353 + x354 + x355*M[0] + x41*M[36] + x43*x96 + x45*x82 + x45*M[28] + x96*M[19] + M[92];
#pragma omp atomic
Ms[93] += x107*x34 + x11*x142 + x11*M[44] + x125*M[12] + x170*x9 + x170*M[9] + x214*M[6] + x234*M[4] + x34*M[37] + x356 + x357 + x358 + x359*x360 + x361*M[1] + x45*x85 + x45*M[29] + x46*x86 + x86*M[20] + M[93];
#pragma omp atomic
Ms[94] += x11*x147 + x11*x148 + x11*x151 + x11*M[45] + x112*x34 + x114*x34 + x116*x34 + x125*x25 + x125*M[13] + x211*z + x214*x7 + x214*M[7] + x234*x4 + x234*M[5] + x34*M[38] + x361*M[2] + x362 + x363 + x364 + x365 + x366 + x367 + x45*x87 + x45*x88 + x45*x89 + x45*M[30] + x49*x86 + x50*x86 + x51*x86 + x86*M[21] + M[94];
#pragma omp atomic
Ms[95] += x107*x41 + x108*x41 + x109*x41 + x11*x153 + x11*x154 + x11*x157 + x11*M[46] + x111*x375 + x118*x34 + x119*x34 + x120*x34 + x121*x19 + x121*M[12] + x125*x30 + x125*M[14] + x177*x9 + x177*M[9] + x208*x5 + x208*M[6] + x214*x8 + x214*M[8] + x218*z + x34*M[39] + x368 + x369 + x370 + x371 + x372 + x373 + x41*M[37] + x45*x90 + x45*x91 + x45*x92 + x45*M[31] + M[95];
#pragma omp atomic
Ms[96] += x11*x160 + x11*x161 + x11*x164 + x11*M[47] + x112*x41 + x113*x41 + x115*x41 + x121*x24 + x121*M[13] + x2*x230 + x208*x6 + x208*M[7] + x224*z + x230*M[4] + x355*M[1] + x376 + x377 + x378 + x379 + x380 + x381 + x41*M[38] + x45*x93 + x45*x94 + x45*x95 + x45*M[32] + x46*x96 + x47*x96 + x48*x96 + x96*M[20] + M[96];
#pragma omp atomic
Ms[97] += x11*x166 + x11*M[48] + x118*x41 + x121*M[14] + x184*x9 + x184*M[9] + x208*M[8] + x230*M[5] + x355*M[2] + x382 + x383 + x384 + x385*M[0] + x41*M[39] + x45*x97 + x45*M[33] + x49*x96 + x96*M[21] + M[97];
#pragma omp atomic
Ms[98] += x0*x292 + x11*x169 + x11*M[49] + x12*x170 + x122*x34 + x125*M[15] + x170*M[10] + x234*M[6] + x292*M[3] + x34*M[40] + x386 + x387 + x388 + x389*M[1] + x52*x86 + x86*M[22] + M[98];
#pragma omp atomic
Ms[99] += x11*x171 + x11*x172 + x11*x173 + x11*M[50] + x125*x36 + x125*M[16] + x126*x34 + x128*x34 + x130*x34 + x15*x170 + x16*x170 + x17*x170 + x170*M[11] + x233*z + x234*x7 + x234*M[7] + x34*M[41] + x389*M[2] + x390 + x391 + x392 + x393 + x394 + x395 + x55*x86 + x57*x86 + x59*x86 + x86*M[23] + M[99];
#pragma omp atomic
Ms[100] += x0*x299 + x11*x174 + x11*x175 + x11*x176 + x11*M[51] + x12*x177 + x121*x33 + x121*M[15] + x122*x41 + x123*x41 + x124*x41 + x125*x39 + x125*M[17] + x132*x34 + x134*x34 + x136*x34 + x177*M[10] + x234*x8 + x234*M[8] + x238*z + x299*M[3] + x34*M[42] + x396 + x397 + x398 + x399 + x400 + x401 + x402*M[1] + x41*M[40] + x61*x86 + x62*x86 + x63*x86 + x86*M[24] + M[100];
#pragma omp atomic
Ms[101] += x0*x303 + x11*x178 + x11*x179 + x11*x180 + x11*M[52] + x121*x35 + x121*M[16] + x125*x42 + x125*M[18] + x126*x41 + x127*x41 + x129*x41 + x138*x34 + x139*x34 + x140*x34 + x15*x177 + x177*M[11] + x230*x5 + x230*M[6] + x244*z + x303*M[3] + x34*M[43] + x402*M[2] + x403 + x404 + x405 + x406 + x407 + x408 + x41*M[41] + x52*x96 + x53*x96 + x54*x96 + x96*M[22] + M[101];
#pragma omp atomic
Ms[102] += x11*x181 + x11*x182 + x11*x183 + x11*M[53] + x12*x184 + x121*x38 + x121*M[17] + x13*x184 + x132*x41 + x133*x41 + x135*x41 + x14*x184 + x184*M[10] + x230*x6 + x230*M[7] + x250*z + x385*M[1] + x409 + x41*M[42] + x410 + x411 + x412 + x413 + x414 + x55*x96 + x56*x96 + x58*x96 + x96*M[23] + M[102];
#pragma omp atomic
Ms[103] += x0*x310 + x11*x185 + x11*M[54] + x121*M[18] + x138*x41 + x15*x184 + x184*M[11] + x230*M[8] + x310*M[3] + x385*M[2] + x41*M[43] + x415 + x416 + x417 + x61*x96 + x96*M[24] + M[103];
#pragma omp atomic
Ms[104] += x1*x292 + x141*x34 + x170*x18 + x170*M[12] + x292*M[4] + x34*M[44] + x418 + x419 + x420 + x421*x422 + x64*x86 + x86*M[25] + M[104];
#pragma omp atomic
Ms[105] += x146*x34 + x148*x34 + x150*x34 + x170*x23 + x170*x25 + x170*x27 + x170*M[13] + x258*z + x292*x3 + x292*x4 + x292*M[5] + x34*M[45] + x423 + x424 + x425 + x426 + x427 + x428 + x69*x86 + x71*x86 + x73*x86 + x86*M[26] + M[105];
#pragma omp atomic
Ms[106] += x1*x299 + x141*x41 + x142*x41 + x143*x41 + x152*x34 + x154*x34 + x156*x34 + x170*x29 + x170*x30 + x170*x31 + x170*M[14] + x177*x18 + x177*M[12] + x264*z + x299*M[4] + x34*M[46] + x360*x435 + x41*M[44] + x429 + x430 + x431 + x432 + x433 + x434 + x75*x86 + x77*x86 + x79*x86 + x86*M[27] + M[106];
#pragma omp atomic
Ms[107] += x1*x303 + x146*x41 + x147*x41 + x149*x41 + x159*x34 + x161*x34 + x163*x34 + x177*x23 + x177*M[13] + x213*x442 + x270*z + x299*x3 + x299*M[5] + x303*M[4] + x34*M[47] + x41*M[45] + x436 + x437 + x438 + x439 + x440 + x441 + x64*x96 + x65*x96 + x66*x96 + x81*x86 + x82*x86 + x83*x86 + x86*M[28] + x96*M[25] + M[107];
#pragma omp atomic
Ms[108] += x111*x449 + x152*x41 + x153*x41 + x155*x41 + x165*x34 + x166*x34 + x167*x34 + x177*x29 + x177*M[14] + x18*x184 + x184*x19 + x184*x20 + x184*M[12] + x277*z + x3*x303 + x303*M[5] + x34*M[48] + x41*M[46] + x443 + x444 + x445 + x446 + x447 + x448 + x69*x96 + x70*x96 + x72*x96 + x96*M[26] + M[108];
#pragma omp atomic
Ms[109] += x1*x310 + x159*x41 + x160*x41 + x162*x41 + x184*x23 + x184*x24 + x184*x26 + x184*M[13] + x2*x310 + x284*z + x310*M[4] + x41*M[47] + x450 + x451 + x452 + x453 + x454 + x455 + x75*x96 + x76*x96 + x78*x96 + x96*M[27] + M[109];
#pragma omp atomic
Ms[110] += x165*x41 + x184*x29 + x184*M[14] + x3*x310 + x310*M[5] + x41*M[48] + x422*x459 + x456 + x457 + x458 + x81*x96 + x96*M[28] + M[110];
#pragma omp atomic
Ms[111] += x170*M[15] + x292*M[6] + x34*M[49] + x460 + x461*M[1] + x86*M[29] + M[111];
#pragma omp atomic
Ms[112] += x170*x36 + x170*M[16] + x172*x34 + x292*x7 + x292*M[7] + x34*M[50] + x461*M[2] + x462 + x463 + x464 + x86*x88 + x86*M[30] + M[112];
#pragma omp atomic
Ms[113] += x169*x41 + x170*x39 + x170*M[17] + x175*x34 + x177*M[15] + x292*x8 + x292*M[8] + x299*M[6] + x34*M[51] + x41*M[49] + x465 + x466 + x467 + x468*M[1] + x86*x91 + x86*M[31] + M[113];
#pragma omp atomic
Ms[114] += x170*x42 + x170*M[18] + x171*x41 + x177*M[16] + x179*x34 + x299*M[7] + x303*M[6] + x34*M[52] + x41*M[50] + x468*M[2] + x469 + x470 + x471 + x472*M[1] + x85*x96 + x86*x94 + x86*M[32] + x96*M[29] + M[114];
#pragma omp atomic
Ms[115] += x174*x41 + x177*M[17] + x182*x34 + x184*x33 + x184*M[15] + x299*M[8] + x303*M[7] + x34*M[53] + x41*M[51] + x472*M[2] + x473 + x474 + x475 + x476*M[1] + x86*x97 + x86*M[33] + x87*x96 + x96*M[30] + M[115];
#pragma omp atomic
Ms[116] += x177*M[18] + x178*x41 + x184*x35 + x184*M[16] + x185*x34 + x303*M[8] + x310*x5 + x310*M[6] + x34*M[54] + x41*M[52] + x476*M[2] + x477 + x478 + x479 + x90*x96 + x96*M[31] + M[116];
#pragma omp atomic
Ms[117] += x181*x41 + x184*x38 + x184*M[17] + x310*x6 + x310*M[7] + x41*M[53] + x480 + x481 + x482 + x483*M[1] + x93*x96 + x96*M[32] + M[117];
#pragma omp atomic
Ms[118] += x184*M[18] + x310*M[8] + x41*M[54] + x483*M[2] + x484 + x96*M[33] + M[118];
#pragma omp atomic
Ms[119] += x100*M[19] + x11*M[55] + x188*M[9] + x314*M[3] + x45*M[34] + x485 + x486*M[0] + M[119];
#pragma omp atomic
Ms[120] += x100*x47 + x100*M[20] + x102*x45 + x11*x190 + x11*M[56] + x13*x188 + x188*M[10] + x2*x314 + x314*M[4] + x45*M[35] + x486*M[1] + x487 + x488 + x489 + M[120];
#pragma omp atomic
Ms[121] += x100*x50 + x100*M[21] + x105*x45 + x11*x193 + x11*M[57] + x16*x188 + x188*M[11] + x312*z + x314*x4 + x314*M[5] + x45*M[36] + x486*M[2] + x490 + z*M[83] + M[121];
#pragma omp atomic
Ms[122] += x100*x53 + x100*M[22] + x108*x45 + x11*x196 + x11*M[58] + x111*x494 + x125*M[19] + x186*x34 + x188*x19 + x188*M[12] + x214*M[9] + x314*x5 + x314*M[6] + x339*M[3] + x34*M[55] + x45*M[37] + x491 + x492 + x493 + M[122];
#pragma omp atomic
Ms[123] += x100*x56 + x100*x57 + x100*x60 + x100*M[23] + x11*x200 + x11*x201 + x11*x204 + x11*M[59] + x113*x45 + x114*x45 + x117*x45 + x188*x24 + x188*x25 + x188*x28 + x188*M[13] + x314*x6 + x314*x7 + x314*M[7] + x315*z + x316*z + x317*z + x45*M[38] + x495 + x496 + x497 + z*M[84] + M[123];
#pragma omp atomic
Ms[124] += x100*x62 + x100*M[24] + x11*x206 + x11*M[60] + x119*x45 + x121*M[19] + x186*x41 + x188*x30 + x188*M[14] + x208*M[9] + x314*x8 + x314*M[8] + x318*z + x334*M[3] + x41*M[55] + x45*M[39] + x498 + x499*M[0] + z*M[85] + M[124];
#pragma omp atomic
Ms[125] += x100*x65 + x100*M[25] + x11*x210 + x11*M[61] + x123*x45 + x125*M[20] + x188*x33 + x188*M[15] + x189*x34 + x213*x503 + x214*M[10] + x234*M[9] + x339*M[4] + x34*M[56] + x361*M[3] + x45*M[40] + x500 + x501 + x502 + x504*M[1] + x86*x98 + x86*M[34] + M[125];
#pragma omp atomic
Ms[126] += x100*x70 + x100*x71 + x100*x74 + x100*M[26] + x11*x216 + x11*x217 + x11*x220 + x11*M[62] + x125*x50 + x125*M[21] + x127*x45 + x128*x45 + x131*x45 + x16*x214 + x188*x35 + x188*x36 + x188*x37 + x188*M[16] + x192*x34 + x193*x34 + x194*x34 + x214*M[11] + x321*z + x322*z + x323*z + x339*x4 + x339*M[5] + x34*M[57] + x45*M[41] + x504*M[2] + x505 + x506 + x507 + z*M[86] + M[126];
#pragma omp atomic
Ms[127] += x100*x76 + x100*x77 + x100*x80 + x100*M[27] + x11*x222 + x11*x223 + x11*x226 + x11*M[63] + x121*x47 + x121*M[20] + x13*x208 + x133*x45 + x134*x45 + x137*x45 + x188*x38 + x188*x39 + x188*x40 + x188*M[17] + x189*x41 + x190*x41 + x191*x41 + x2*x334 + x208*M[10] + x325*z + x326*z + x328*z + x334*M[4] + x41*M[56] + x45*M[42] + x499*M[1] + x508 + x509 + x510 + z*M[87] + M[127];
#pragma omp atomic
Ms[128] += x100*x82 + x100*M[28] + x11*x228 + x11*M[64] + x121*M[21] + x139*x45 + x188*x42 + x188*M[18] + x192*x41 + x208*M[11] + x230*M[9] + x331*z + x334*M[5] + x355*M[3] + x41*M[57] + x45*M[43] + x499*M[2] + x511 + x512*M[0] + x96*x98 + x96*M[34] + z*M[88] + M[128];
#pragma omp atomic
Ms[129] += x100*x85 + x100*M[29] + x101*x86 + x11*x232 + x11*M[65] + x125*M[22] + x142*x45 + x170*x43 + x170*M[19] + x195*x34 + x214*M[12] + x234*M[10] + x339*M[6] + x34*M[58] + x360*x516 + x361*M[4] + x389*M[3] + x45*M[44] + x513 + x514 + x515 + x517*M[1] + x86*M[35] + M[129];
#pragma omp atomic
Ms[130] += x100*x87 + x100*x88 + x100*x89 + x100*M[30] + x104*x86 + x105*x86 + x106*x86 + x11*x236 + x11*x237 + x11*x240 + x11*M[66] + x125*x57 + x125*M[23] + x147*x45 + x148*x45 + x151*x45 + x16*x234 + x199*x34 + x201*x34 + x203*x34 + x214*x25 + x214*M[13] + x234*M[11] + x335*z + x336*z + x337*z + x339*x7 + x339*M[7] + x34*M[59] + x361*x4 + x361*M[5] + x45*M[45] + x517*M[2] + x518 + x519 + x520 + x86*M[36] + z*M[89] + M[130];
#pragma omp atomic
Ms[131] += x100*x90 + x100*x91 + x100*x92 + x100*M[31] + x11*x242 + x11*x243 + x11*x246 + x11*M[67] + x121*x53 + x121*M[22] + x125*x62 + x125*M[24] + x153*x45 + x154*x45 + x157*x45 + x177*x43 + x177*M[19] + x19*x208 + x195*x41 + x196*x41 + x197*x41 + x205*x34 + x206*x34 + x207*x34 + x208*M[12] + x214*x30 + x214*M[14] + x334*x5 + x334*M[6] + x339*x8 + x339*M[8] + x34*M[60] + x340*z + x341*z + x343*z + x402*M[3] + x41*M[58] + x45*M[46] + x521 + x522 + x523 + x524*x525 + z*M[90] + M[131];
#pragma omp atomic
Ms[132] += x100*x93 + x100*x94 + x100*x95 + x100*M[32] + x101*x96 + x102*x96 + x103*x96 + x11*x248 + x11*x249 + x11*x252 + x11*M[68] + x121*x56 + x121*M[23] + x13*x230 + x160*x45 + x161*x45 + x164*x45 + x199*x41 + x2*x355 + x200*x41 + x202*x41 + x208*x24 + x208*M[13] + x230*M[10] + x334*x6 + x334*M[7] + x346*z + x347*z + x349*z + x355*M[4] + x41*M[59] + x45*M[47] + x512*M[1] + x526 + x527 + x528 + x96*M[35] + z*M[91] + M[132];
#pragma omp atomic
Ms[133] += x100*x97 + x100*M[33] + x104*x96 + x11*x254 + x11*M[69] + x121*M[24] + x166*x45 + x184*x43 + x184*M[19] + x205*x41 + x208*M[14] + x230*M[11] + x334*M[8] + x352*z + x355*M[5] + x385*M[3] + x41*M[60] + x45*M[48] + x512*M[2] + x529 + x530*M[0] + x96*M[36] + z*M[92] + M[133];
#pragma omp atomic
Ms[134] += x107*x86 + x11*x257 + x11*M[70] + x125*M[25] + x169*x45 + x170*x46 + x170*M[20] + x209*x34 + x214*M[15] + x234*M[12] + x259*x535 + x292*x9 + x292*M[9] + x34*M[61] + x361*M[6] + x389*M[4] + x45*M[49] + x531 + x532 + x533 + x536*M[1] + x86*M[37] + M[134];
#pragma omp atomic
Ms[135] += x11*x262 + x11*x263 + x11*x266 + x11*M[71] + x112*x86 + x114*x86 + x116*x86 + x125*x71 + x125*M[26] + x170*x49 + x170*x50 + x170*x51 + x170*M[21] + x171*x45 + x172*x45 + x173*x45 + x214*x36 + x214*M[16] + x215*x34 + x217*x34 + x219*x34 + x234*x25 + x234*M[13] + x34*M[62] + x356*z + x357*z + x358*z + x361*x7 + x361*M[7] + x389*x4 + x389*M[5] + x45*M[50] + x536*M[2] + x537 + x538 + x539 + x86*M[38] + z*M[93] + M[135];
#pragma omp atomic
Ms[136] += x11*x268 + x11*x269 + x11*x272 + x11*M[72] + x118*x86 + x119*x86 + x120*x86 + x121*x65 + x121*M[25] + x125*x77 + x125*M[27] + x145*x543 + x174*x45 + x175*x45 + x176*x45 + x177*x46 + x177*M[20] + x208*x33 + x208*M[15] + x209*x41 + x210*x41 + x211*x41 + x214*x39 + x214*M[17] + x221*x34 + x223*x34 + x225*x34 + x234*x30 + x234*M[14] + x299*x9 + x299*M[9] + x34*M[63] + x361*x8 + x361*M[8] + x362*z + x363*z + x365*z + x402*M[4] + x41*M[61] + x45*M[51] + x540 + x541 + x542 + x545*M[1] + x86*M[39] + z*M[94] + M[136];
#pragma omp atomic
Ms[137] += x107*x96 + x108*x96 + x109*x96 + x11*x275 + x11*x276 + x11*x279 + x11*M[73] + x121*x70 + x121*M[26] + x125*x82 + x125*M[28] + x177*x49 + x177*M[21] + x178*x45 + x179*x45 + x180*x45 + x19*x230 + x208*x35 + x208*M[16] + x214*x42 + x214*M[18] + x215*x41 + x216*x41 + x218*x41 + x227*x34 + x228*x34 + x229*x34 + x230*M[12] + x303*x9 + x303*M[9] + x34*M[64] + x355*x5 + x355*M[6] + x368*z + x369*z + x371*z + x402*M[5] + x41*M[62] + x45*M[52] + x524*x549 + x545*M[2] + x546 + x547 + x548 + x96*M[37] + z*M[95] + M[137];
#pragma omp atomic
Ms[138] += x11*x282 + x11*x283 + x11*x286 + x11*M[74] + x112*x96 + x113*x96 + x115*x96 + x121*x76 + x121*M[27] + x181*x45 + x182*x45 + x183*x45 + x184*x46 + x184*x47 + x184*x48 + x184*M[20] + x2*x385 + x208*x38 + x208*M[17] + x221*x41 + x222*x41 + x224*x41 + x230*x24 + x230*M[13] + x355*x6 + x355*M[7] + x376*z + x377*z + x379*z + x385*M[4] + x41*M[63] + x45*M[53] + x530*M[1] + x550 + x551 + x552 + x96*M[38] + z*M[96] + M[138];
#pragma omp atomic
Ms[139] += x11*x288 + x11*M[75] + x118*x96 + x121*M[28] + x184*x49 + x184*M[21] + x185*x45 + x208*M[18] + x227*x41 + x230*M[14] + x290*x535 + x310*x9 + x310*M[9] + x355*M[8] + x382*z + x385*M[5] + x41*M[64] + x45*M[54] + x530*M[2] + x553 + x96*M[39] + z*M[97] + M[139];
#pragma omp atomic
Ms[140] += x0*x461 + x11*x291 + x11*M[76] + x12*x292 + x122*x86 + x125*M[29] + x170*x52 + x170*M[22] + x231*x34 + x234*M[15] + x292*M[10] + x34*M[65] + x389*M[6] + x461*M[3] + x554 + x555 + x556 + x557*M[1] + x86*M[40] + M[140];
#pragma omp atomic
Ms[141] += x11*x293 + x11*x294 + x11*x295 + x11*M[77] + x125*x88 + x125*M[30] + x126*x86 + x128*x86 + x130*x86 + x15*x292 + x16*x292 + x17*x292 + x170*x55 + x170*x57 + x170*x59 + x170*M[23] + x234*x36 + x234*M[16] + x235*x34 + x237*x34 + x239*x34 + x292*M[11] + x34*M[66] + x386*z + x387*z + x388*z + x389*x7 + x389*M[7] + x557*M[2] + x558 + x559 + x560 + x86*M[41] + z*M[98] + M[141];
#pragma omp atomic
Ms[142] += x0*x468 + x11*x296 + x11*x297 + x11*x298 + x11*M[78] + x12*x299 + x121*x85 + x121*M[29] + x125*x91 + x125*M[31] + x132*x86 + x134*x86 + x136*x86 + x170*x61 + x170*x62 + x170*x63 + x170*M[24] + x177*x52 + x177*M[22] + x231*x41 + x232*x41 + x233*x41 + x234*x39 + x234*M[17] + x241*x34 + x243*x34 + x245*x34 + x299*M[10] + x34*M[67] + (1.0/24.0)*x374*x564 + x389*x8 + x389*M[8] + x390*z + x391*z + x393*z + x402*M[6] + x41*M[65] + x468*M[3] + x561 + x562 + x563 + x86*M[42] + z*M[99] + M[142];
#pragma omp atomic
Ms[143] += x0*x472 + x10*x544*x568 + x11*x300 + x11*x301 + x11*x302 + x11*M[79] + x12*x303 + x121*x87 + x121*M[30] + x122*x96 + x123*x96 + x124*x96 + x125*x94 + x125*M[32] + x138*x86 + x139*x86 + x140*x86 + x15*x299 + x177*x55 + x177*M[23] + x230*x33 + x230*M[15] + x234*x42 + x234*M[18] + x235*x41 + x236*x41 + x238*x41 + x247*x34 + x249*x34 + x251*x34 + x299*M[11] + x303*M[10] + x34*M[68] + x396*z + x397*z + x399*z + x402*M[7] + x41*M[66] + x472*M[3] + x565 + x566 + x567 + x569*M[2] + x86*M[43] + x96*M[40] + z*M[100] + M[143];
#pragma omp atomic
Ms[144] += x0*x476 + x11*x304 + x11*x305 + x11*x306 + x11*M[80] + x121*x90 + x121*M[31] + x125*x97 + x125*M[33] + x126*x96 + x127*x96 + x129*x96 + x15*x303 + x177*x61 + x177*M[24] + x184*x52 + x184*x53 + x184*x54 + x184*M[22] + x230*x35 + x230*M[16] + x241*x41 + x242*x41 + x244*x41 + x253*x34 + x254*x34 + x255*x34 + x303*M[11] + x34*M[69] + x385*x5 + x385*M[6] + x402*M[8] + x403*z + x404*z + x406*z + x41*M[67] + x476*M[3] + x570 + x571 + x572 + x573*M[2] + x96*M[41] + z*M[101] + M[144];
#pragma omp atomic
Ms[145] += x11*x307 + x11*x308 + x11*x309 + x11*M[81] + x12*x310 + x121*x93 + x121*M[32] + x13*x310 + x132*x96 + x133*x96 + x135*x96 + x14*x310 + x184*x55 + x184*x56 + x184*x58 + x184*M[23] + x230*x38 + x230*M[17] + x247*x41 + x248*x41 + x250*x41 + x310*M[10] + x385*x6 + x385*M[7] + x409*z + x41*M[68] + x410*z + x412*z + x574 + x575 + x576 + x577*M[1] + x96*M[42] + z*M[102] + M[145];
#pragma omp atomic
Ms[146] += x0*x483 + x11*x311 + x11*M[82] + x121*M[33] + x138*x96 + x15*x310 + x184*x61 + x184*M[24] + x230*M[18] + x253*x41 + x310*M[11] + x385*M[8] + x41*M[69] + x415*z + x483*M[3] + x577*M[2] + x578 + x96*M[43] + z*M[103] + M[146];
#pragma omp atomic
Ms[147] += x1*x461 + x141*x86 + x170*x64 + x170*M[25] + x18*x292 + x256*x34 + x292*M[12] + x34*M[70] + x461*M[4] + x579 + x580 + x581 + x582*x583 + x86*M[44] + M[147];
#pragma omp atomic
Ms[148] += x146*x86 + x148*x86 + x150*x86 + x170*x69 + x170*x71 + x170*x73 + x170*M[26] + x23*x292 + x25*x292 + x261*x34 + x263*x34 + x265*x34 + x27*x292 + x292*M[13] + x3*x461 + x34*M[71] + x4*x461 + x418*z + x419*z + x420*z + x461*M[5] + x584 + x585 + x586 + x86*M[45] + z*M[104] + M[148];
#pragma omp atomic
Ms[149] += x1*x468 + x152*x86 + x154*x86 + x156*x86 + x170*x75 + x170*x77 + x170*x79 + x170*M[27] + x177*x64 + x177*M[25] + x18*x299 + x256*x41 + x257*x41 + x258*x41 + x267*x34 + x269*x34 + x271*x34 + x29*x292 + x292*x30 + x292*x31 + x292*M[14] + x299*M[12] + x34*M[72] + x41*M[70] + x423*z + x424*z + x426*z + x468*M[4] + x587 + x588 + x589 + x590*x591 + x86*M[46] + z*M[105] + M[149];
#pragma omp atomic
Ms[150] += x1*x472 + x141*x96 + x142*x96 + x143*x96 + x159*x86 + x161*x86 + x163*x86 + x170*x81 + x170*x82 + x170*x83 + x170*M[28] + x177*x69 + x177*M[26] + x18*x303 + x23*x299 + x261*x41 + x262*x41 + x264*x41 + x274*x34 + x276*x34 + x278*x34 + x299*M[13] + x3*x468 + x303*M[12] + x34*M[73] + x360*x595 + x41*M[71] + x429*z + x430*z + x432*z + x468*M[5] + x472*M[4] + x592 + x593 + x594 + x86*M[47] + x96*M[44] + z*M[106] + M[150];
#pragma omp atomic
Ms[151] += x1*x476 + x146*x96 + x147*x96 + x149*x96 + x165*x86 + x166*x86 + x167*x86 + x177*x75 + x177*M[27] + x184*x64 + x184*x65 + x184*x66 + x184*M[25] + x213*x599 + x23*x303 + x267*x41 + x268*x41 + x270*x41 + x281*x34 + x283*x34 + x285*x34 + x29*x299 + x299*M[14] + x3*x472 + x303*M[13] + x34*M[74] + x41*M[72] + x436*z + x437*z + x439*z + x472*M[5] + x476*M[4] + x596 + x597 + x598 + x86*M[48] + x96*M[45] + z*M[107] + M[151];
#pragma omp atomic
Ms[152] += x111*x603 + x152*x96 + x153*x96 + x155*x96 + x177*x81 + x177*M[28] + x18*x310 + x184*x69 + x184*x70 + x184*x72 + x184*M[26] + x19*x310 + x20*x310 + x274*x41 + x275*x41 + x277*x41 + x287*x34 + x288*x34 + x289*x34 + x29*x303 + x3*x476 + x303*M[14] + x310*M[12] + x34*M[75] + x41*M[73] + x443*z + x444*z + x446*z + x476*M[5] + x600 + x601 + x602 + x96*M[46] + z*M[108] + M[152];
#pragma omp atomic
Ms[153] += x1*x483 + x159*x96 + x160*x96 + x162*x96 + x184*x75 + x184*x76 + x184*x78 + x184*M[27] + x2*x483 + x23*x310 + x24*x310 + x26*x310 + x281*x41 + x282*x41 + x284*x41 + x310*M[13] + x41*M[74] + x450*z + x451*z + x453*z + x483*M[4] + x604 + x605 + x606 + x96*M[47] + z*M[109] + M[153];
#pragma omp atomic
Ms[154] += x165*x96 + x184*x81 + x184*M[28] + x287*x41 + x29*x310 + x3*x483 + x310*M[14] + x41*M[75] + x456*z + x483*M[5] + x583*x608 + x607 + x96*M[48] + z*M[110] + M[154];
#pragma omp atomic
Ms[155] += x170*M[29] + x292*M[15] + x34*M[76] + x461*M[6] + x609 + x610*M[1] + x86*M[49] + M[155];
#pragma omp atomic
Ms[156] += x170*x88 + x170*M[30] + x172*x86 + x292*x36 + x292*M[16] + x294*x34 + x34*M[77] + x460*z + x461*x7 + x461*M[7] + x610*M[2] + x611 + x86*M[50] + z*M[111] + M[156];
#pragma omp atomic
Ms[157] += x170*x91 + x170*M[31] + x175*x86 + x177*M[29] + x291*x41 + x292*x39 + x292*M[17] + x297*x34 + x299*M[15] + x34*M[78] + x41*M[76] + x461*x8 + x461*M[8] + x462*z + x468*M[6] + x612 + x613*M[1] + x86*M[51] + z*M[112] + M[157];
#pragma omp atomic
Ms[158] += x169*x96 + x170*x94 + x170*M[32] + x177*M[30] + x179*x86 + x292*x42 + x292*M[18] + x293*x41 + x299*M[16] + x301*x34 + x303*M[15] + x34*M[79] + x41*M[77] + x465*z + x468*M[7] + x472*M[6] + x613*M[2] + x614 + x615*M[1] + x86*M[52] + x96*M[49] + z*M[113] + M[158];
#pragma omp atomic
Ms[159] += x170*x97 + x170*M[33] + x171*x96 + x177*M[31] + x182*x86 + x184*x85 + x184*M[29] + x296*x41 + x299*M[17] + x303*M[16] + x305*x34 + x34*M[80] + x41*M[78] + x468*M[8] + x469*z + x472*M[7] + x476*M[6] + x615*M[2] + x616 + x617*M[1] + x86*M[53] + x96*M[50] + z*M[114] + M[159];
#pragma omp atomic
Ms[160] += x174*x96 + x177*M[32] + x184*x87 + x184*M[30] + x185*x86 + x299*M[18] + x300*x41 + x303*M[17] + x308*x34 + x310*x33 + x310*M[15] + x34*M[81] + x41*M[79] + x472*M[8] + x473*z + x476*M[7] + x617*M[2] + x618 + x619*M[1] + x86*M[54] + x96*M[51] + z*M[115] + M[160];
#pragma omp atomic
Ms[161] += x177*M[33] + x178*x96 + x184*x90 + x184*M[31] + x303*M[18] + x304*x41 + x310*x35 + x310*M[16] + x311*x34 + x34*M[82] + x41*M[80] + x476*M[8] + x477*z + x483*x5 + x483*M[6] + x619*M[2] + x620 + x96*M[52] + z*M[116] + M[161];
#pragma omp atomic
Ms[162] += x181*x96 + x184*x93 + x184*M[32] + x307*x41 + x310*x38 + x310*M[17] + x41*M[81] + x480*z + x483*x6 + x483*M[7] + x621 + x622*M[1] + x96*M[53] + z*M[117] + M[162];
#pragma omp atomic
Ms[163] += x184*M[33] + x310*M[18] + x41*M[82] + x483*M[8] + x622*M[2] + x96*M[54] + z*M[118] + M[163];
#pragma omp atomic
Ms[164] += x*M[119] + x100*M[34] + x11*M[83] + x188*M[19] + x314*M[9] + x45*M[55] + x486*M[3] + x623*M[0] + M[164];
#pragma omp atomic
Ms[165] += x*M[120] + x100*x102 + x100*M[35] + x11*x316 + x11*M[84] + x13*x314 + x188*x47 + x188*M[20] + x190*x45 + x2*x486 + x314*M[10] + x45*M[56] + x485*y + x486*M[4] + x623*M[1] + y*M[119] + M[165];
#pragma omp atomic
Ms[166] += x*M[121] + x100*x105 + x100*M[36] + x11*x319 + x11*M[85] + x16*x314 + x188*x50 + x188*M[21] + x193*x45 + x314*M[11] + x4*x486 + x45*M[57] + x485*z + x486*M[5] + x623*M[2] + z*M[119] + M[166];
#pragma omp atomic
Ms[167] += x*M[122] + x100*x108 + x100*M[37] + x11*x322 + x11*M[86] + x111*x624 + x125*M[34] + x188*x53 + x188*M[22] + x19*x314 + x196*x45 + x214*M[19] + x312*x34 + x314*M[12] + x339*M[9] + x34*M[83] + x45*M[58] + x486*x5 + x486*M[6] + x487*y + x504*M[3] + y*M[120] + M[167];
#pragma omp atomic
Ms[168] += x*M[123] + x100*x113 + x100*x114 + x100*x117 + x100*M[38] + x11*x326 + x11*x327 + x11*x330 + x11*M[87] + x188*x56 + x188*x57 + x188*x60 + x188*M[23] + x200*x45 + x201*x45 + x204*x45 + x24*x314 + x25*x314 + x28*x314 + x314*M[13] + x45*M[59] + x486*x6 + x486*x7 + x486*M[7] + x487*z + x488*z + x489*z + x490*y + y*M[121] + z*M[120] + M[168];
#pragma omp atomic
Ms[169] += x*M[124] + x100*x119 + x100*M[39] + x11*x332 + x11*M[88] + x121*M[34] + x188*x62 + x188*M[24] + x206*x45 + x208*M[19] + x30*x314 + x312*x41 + x314*M[14] + x334*M[9] + x41*M[83] + x45*M[60] + x486*x8 + x486*M[8] + x490*z + x499*M[3] + x590*x624 + z*M[121] + M[169];
#pragma omp atomic
Ms[170] += x*M[125] + x100*x123 + x100*M[40] + x11*x336 + x11*M[89] + x125*M[35] + x186*x86 + x188*x65 + x188*M[25] + x210*x45 + x214*M[20] + x234*M[19] + x314*x33 + x314*M[15] + x315*x34 + x339*M[10] + x34*M[84] + x361*M[9] + x45*M[61] + x491*y + x504*M[4] + x517*M[3] + x625*x67 + x626*M[1] + x86*M[55] + y*M[122] + M[170];
#pragma omp atomic
Ms[171] += x*M[126] + x100*x127 + x100*x128 + x100*x131 + x100*M[41] + x105*x125 + x11*x341 + x11*x342 + x11*x345 + x11*M[90] + x125*M[36] + x16*x339 + x188*x70 + x188*x71 + x188*x74 + x188*M[26] + x214*x50 + x214*M[21] + x216*x45 + x217*x45 + x220*x45 + x314*x35 + x314*x36 + x314*x37 + x314*M[16] + x318*x34 + x319*x34 + x320*x34 + x339*M[11] + x34*M[85] + x4*x504 + x45*M[62] + x491*z + x492*z + x493*z + x495*y + x504*M[5] + x626*M[2] + y*M[123] + z*M[122] + M[171];
#pragma omp atomic
Ms[172] += x*M[127] + x100*x133 + x100*x134 + x100*x137 + x100*M[42] + x102*x121 + x11*x347 + x11*x348 + x11*x351 + x11*M[91] + x121*M[35] + x13*x334 + x188*x76 + x188*x77 + x188*x80 + x188*M[27] + x2*x499 + x208*x47 + x208*M[20] + x222*x45 + x223*x45 + x226*x45 + x314*x38 + x314*x39 + x314*x40 + x314*M[17] + x315*x41 + x316*x41 + x317*x41 + x334*M[10] + x41*M[84] + x45*M[63] + x495*z + x496*z + x497*z + x498*y + x499*M[4] + x627*M[1] + y*M[124] + z*M[123] + M[172];
#pragma omp atomic
Ms[173] += x*M[128] + x100*x139 + x100*M[43] + x11*x353 + x11*M[92] + x121*M[36] + x186*x96 + x188*x82 + x188*M[28] + x208*M[21] + x228*x45 + x230*M[19] + x314*x42 + x314*M[18] + x318*x41 + x334*M[11] + x355*M[9] + x41*M[85] + x45*M[64] + x498*z + x499*M[5] + x512*M[3] + x625*x84 + x627*M[2] + x96*M[55] + z*M[124] + M[173];
#pragma omp atomic
Ms[174] += x*M[129] + x100*x142 + x100*M[44] + x11*x357 + x11*M[93] + x125*M[37] + x170*x98 + x170*M[34] + x188*x85 + x188*M[29] + x189*x86 + x214*M[22] + x232*x45 + x234*M[20] + x321*x34 + x339*M[12] + x34*M[86] + x360*x628 + x361*M[10] + x389*M[9] + x45*M[65] + x500*y + x504*M[6] + x517*M[4] + x536*M[3] + x564*x629 + x86*M[56] + y*M[125] + M[174];
#pragma omp atomic
Ms[175] += x*M[130] + x100*x147 + x100*x148 + x100*x151 + x100*M[45] + x11*x363 + x11*x364 + x11*x367 + x11*M[94] + x114*x125 + x125*M[38] + x16*x361 + x188*x87 + x188*x88 + x188*x89 + x188*M[30] + x192*x86 + x193*x86 + x194*x86 + x214*x57 + x214*M[23] + x234*x50 + x234*M[21] + x236*x45 + x237*x45 + x240*x45 + x25*x339 + x325*x34 + x327*x34 + x329*x34 + x339*M[13] + x34*M[87] + x361*M[11] + x4*x517 + x45*M[66] + x500*z + x501*z + x502*z + x504*x7 + x504*M[7] + x505*y + x517*M[5] + x630*x67 + x86*M[57] + y*M[126] + z*M[125] + M[175];
#pragma omp atomic
Ms[176] += x*M[131] + x100*x153 + x100*x154 + x100*x157 + x100*M[46] + x108*x121 + x11*x369 + x11*x370 + x11*x373 + x11*M[95] + x119*x125 + x121*M[37] + x125*M[39] + x177*x98 + x177*M[34] + x188*x90 + x188*x91 + x188*x92 + x188*M[31] + x19*x334 + x208*x53 + x208*M[22] + x214*x62 + x214*M[24] + x242*x45 + x243*x45 + x246*x45 + x30*x339 + x321*x41 + x322*x41 + x323*x41 + x331*x34 + x332*x34 + x333*x34 + x334*M[12] + x339*M[14] + x34*M[88] + x402*M[9] + x41*M[86] + x45*M[67] + x499*x5 + x499*M[6] + x504*x8 + x504*M[8] + x505*z + x506*z + x507*z + x508*y + x545*M[3] + x631*x632 + y*M[127] + z*M[126] + M[176];
#pragma omp atomic
Ms[177] += x*M[132] + x100*x160 + x100*x161 + x100*x164 + x100*M[47] + x11*x377 + x11*x378 + x11*x381 + x11*M[96] + x113*x121 + x121*M[38] + x13*x355 + x188*x93 + x188*x94 + x188*x95 + x188*M[32] + x189*x96 + x190*x96 + x191*x96 + x2*x512 + x208*x56 + x208*M[23] + x230*x47 + x230*M[20] + x24*x334 + x248*x45 + x249*x45 + x252*x45 + x325*x41 + x326*x41 + x328*x41 + x334*M[13] + x355*M[10] + x41*M[87] + x45*M[68] + x499*x6 + x499*M[7] + x508*z + x509*z + x510*z + x511*y + x512*M[4] + x568*x629 + x96*M[56] + y*M[128] + z*M[127] + M[177];
#pragma omp atomic
Ms[178] += x*M[133] + x100*x166 + x100*M[48] + x11*x383 + x11*M[97] + x121*M[39] + x184*x98 + x184*M[34] + x188*x97 + x188*M[33] + x192*x96 + x208*M[24] + x230*M[21] + x254*x45 + x331*x41 + x334*M[14] + x355*M[11] + x385*M[9] + x41*M[88] + x45*M[69] + x499*M[8] + x511*z + x512*M[5] + x530*M[3] + x630*x84 + x633*M[0] + x96*M[57] + z*M[128] + M[178];
#pragma omp atomic
Ms[179] += x*M[134] + x100*x169 + x100*M[49] + x101*x170 + x11*x387 + x11*M[98] + x125*M[40] + x170*M[35] + x195*x86 + x214*M[25] + x234*M[22] + x257*x45 + x259*x634 + x292*x43 + x292*M[19] + x335*x34 + x339*M[15] + x34*M[89] + x361*M[12] + x389*M[10] + x45*M[70] + x513*y + x517*M[6] + x536*M[4] + x557*M[3] + x635*M[1] + x86*M[58] + y*M[129] + M[179];
#pragma omp atomic
Ms[180] += x*M[135] + x100*x171 + x100*x172 + x100*x173 + x100*M[50] + x104*x170 + x105*x170 + x106*x170 + x11*x391 + x11*x392 + x11*x395 + x11*M[99] + x125*x128 + x125*M[41] + x16*x389 + x170*M[36] + x199*x86 + x201*x86 + x203*x86 + x214*x71 + x214*M[26] + x234*x57 + x234*M[23] + x25*x361 + x262*x45 + x263*x45 + x266*x45 + x339*x36 + x339*M[16] + x34*x340 + x34*x342 + x34*x344 + x34*M[90] + x361*M[13] + x389*M[11] + x4*x536 + x45*M[71] + x513*z + x514*z + x515*z + x517*x7 + x517*M[7] + x518*y + x536*M[5] + x635*M[2] + x86*M[59] + y*M[130] + z*M[129] + M[180];
#pragma omp atomic
Ms[181] += x*M[136] + x100*x174 + x100*x175 + x100*x176 + x100*M[51] + x101*x177 + x11*x397 + x11*x398 + x11*x401 + x11*M[100] + x121*x123 + x121*M[40] + x125*x134 + x125*M[42] + x177*M[35] + x205*x86 + x206*x86 + x207*x86 + x208*x65 + x208*M[25] + x213*x636 + x214*x77 + x214*M[27] + x234*x62 + x234*M[24] + x268*x45 + x269*x45 + x272*x45 + x299*x43 + x299*M[19] + x30*x361 + x33*x334 + x334*M[15] + x335*x41 + x336*x41 + x337*x41 + x339*x39 + x339*M[17] + x34*x346 + x34*x348 + x34*x350 + x34*M[91] + x361*M[14] + x402*M[10] + x41*M[89] + x45*M[72] + x517*x8 + x517*M[8] + x518*z + x519*z + x520*z + x521*y + x545*M[4] + x569*M[3] + x637*x638 + x86*M[60] + y*M[131] + z*M[130] + M[181];
#pragma omp atomic
Ms[182] += x*M[137] + x100*x178 + x100*x179 + x100*x180 + x100*M[52] + x104*x177 + x11*x404 + x11*x405 + x11*x408 + x11*M[101] + x111*x640 + x121*x127 + x121*M[41] + x125*x139 + x125*M[43] + x177*M[36] + x19*x355 + x195*x96 + x196*x96 + x197*x96 + x208*x70 + x208*M[26] + x214*x82 + x214*M[28] + x230*x53 + x230*M[22] + x275*x45 + x276*x45 + x279*x45 + x303*x43 + x303*M[19] + x334*x35 + x334*M[16] + x339*x42 + x339*M[18] + x34*x352 + x34*x353 + x34*x354 + x34*M[92] + x340*x41 + x341*x41 + x343*x41 + x355*M[12] + x402*M[11] + x41*M[90] + x45*M[73] + x5*x512 + x512*M[6] + x521*z + x522*z + x523*z + x526*y + x545*M[5] + x573*M[3] + x638*x641 + x96*M[58] + y*M[132] + z*M[131] + M[182];
#pragma omp atomic
Ms[183] += x*M[138] + x100*x181 + x100*x182 + x100*x183 + x100*M[53] + x101*x184 + x102*x184 + x103*x184 + x11*x410 + x11*x411 + x11*x414 + x11*M[102] + x121*x133 + x121*M[42] + x13*x385 + x184*M[35] + x199*x96 + x2*x530 + x200*x96 + x202*x96 + x208*x76 + x208*M[27] + x230*x56 + x230*M[23] + x24*x355 + x282*x45 + x283*x45 + x286*x45 + x334*x38 + x334*M[17] + x346*x41 + x347*x41 + x349*x41 + x355*M[13] + x385*M[10] + x41*M[91] + x45*M[74] + x512*x6 + x512*M[7] + x526*z + x527*z + x528*z + x529*y + x530*M[4] + x633*M[1] + x96*M[59] + y*M[133] + z*M[132] + M[183];
#pragma omp atomic
Ms[184] += x*M[139] + x100*x185 + x100*M[54] + x104*x184 + x11*x416 + x11*M[103] + x121*M[43] + x184*M[36] + x205*x96 + x208*M[28] + x230*M[24] + x288*x45 + x290*x634 + x310*x43 + x310*M[19] + x334*M[18] + x352*x41 + x355*M[14] + x385*M[11] + x41*M[92] + x45*M[75] + x512*M[8] + x529*z + x530*M[5] + x577*M[3] + x633*M[2] + x96*M[60] + z*M[133] + M[184];
#pragma omp atomic
Ms[185] += x*M[140] + x107*x170 + x11*x419 + x11*M[104] + x125*M[44] + x170*M[37] + x209*x86 + x214*M[29] + x234*M[25] + x291*x45 + x292*x46 + x292*M[20] + x34*x356 + x34*M[93] + x361*M[15] + x389*M[12] + x421*x643 + x45*M[76] + x461*x9 + x461*M[9] + x531*y + x536*M[6] + x557*M[4] + x645*M[1] + x86*M[61] + y*M[134] + M[185];
#pragma omp atomic
Ms[186] += x*M[141] + x11*x424 + x11*x425 + x11*x428 + x11*M[105] + x112*x170 + x114*x170 + x116*x170 + x125*x148 + x125*M[45] + x170*M[38] + x214*x88 + x214*M[30] + x215*x86 + x217*x86 + x219*x86 + x234*x71 + x234*M[26] + x25*x389 + x292*x49 + x292*x50 + x292*x51 + x292*M[21] + x293*x45 + x294*x45 + x295*x45 + x34*x362 + x34*x364 + x34*x366 + x34*M[94] + x36*x361 + x361*M[16] + x389*M[13] + x4*x557 + x45*M[77] + x531*z + x532*z + x533*z + x536*x7 + x536*M[7] + x537*y + x557*M[5] + x645*M[2] + x86*M[62] + y*M[135] + z*M[134] + M[186];
#pragma omp atomic
Ms[187] += x*M[142] + x107*x177 + x11*x430 + x11*x431 + x11*x434 + x11*M[106] + x118*x170 + x119*x170 + x120*x170 + x121*x142 + x121*M[44] + x125*x154 + x125*M[46] + x170*M[39] + x177*M[37] + x208*x85 + x208*M[29] + x214*x91 + x214*M[31] + x221*x86 + x223*x86 + x225*x86 + x234*x77 + x234*M[27] + x296*x45 + x297*x45 + x298*x45 + x299*x46 + x299*M[20] + x30*x389 + x34*x368 + x34*x370 + x34*x372 + x34*M[95] + x356*x41 + x357*x41 + x358*x41 + (1.0/96.0)*x360*x374 + x361*x39 + x361*M[17] + x389*M[14] + x402*M[12] + x41*M[93] + x45*M[78] + x468*x9 + x468*M[9] + x536*x8 + x536*M[8] + x537*z + x538*z + x539*z + x540*y + x545*M[6] + x564*x636 + x569*M[4] + x86*M[63] + y*M[136] + z*M[135] + M[187];
#pragma omp atomic
Ms[188] += x*M[143] + x11*x437 + x11*x438 + x11*x441 + x11*M[107] + x112*x177 + x121*x147 + x121*M[45] + x125*x161 + x125*M[47] + x177*M[38] + x208*x87 + x208*M[30] + x209*x96 + x21*x568*x639 + x210*x96 + x211*x96 + x213*x646 + x214*x94 + x214*M[32] + x227*x86 + x228*x86 + x229*x86 + x230*x65 + x230*M[25] + x234*x82 + x234*M[28] + x299*x49 + x299*M[21] + x300*x45 + x301*x45 + x302*x45 + x303*x46 + x303*M[20] + x33*x355 + x34*x376 + x34*x378 + x34*x380 + x34*M[96] + x355*M[15] + x361*x42 + x361*M[18] + x362*x41 + x363*x41 + x365*x41 + x402*M[13] + x41*M[94] + x45*M[79] + x472*x9 + x472*M[9] + x540*z + x541*z + x542*z + x545*M[7] + x546*y + x569*M[5] + x573*M[4] + x636*x647 + x86*M[64] + x96*M[61] + y*M[137] + z*M[136] + M[188];
#pragma omp atomic
Ms[189] += x*M[144] + x107*x184 + x108*x184 + x109*x184 + x11*x444 + x11*x445 + x11*x448 + x11*M[108] + x118*x177 + x121*x153 + x121*M[46] + x125*x166 + x125*M[48] + x177*M[39] + x184*M[37] + x19*x385 + x208*x90 + x208*M[31] + x21*x640*M[2] + x214*x97 + x214*M[33] + x215*x96 + x216*x96 + x218*x96 + x230*x70 + x230*M[26] + x303*x49 + x303*M[21] + x304*x45 + x305*x45 + x306*x45 + x34*x382 + x34*x383 + x34*x384 + x34*M[97] + x35*x355 + x355*M[16] + x368*x41 + x369*x41 + x371*x41 + x385*M[12] + x402*M[14] + x41*M[95] + x45*M[80] + x476*x9 + x476*M[9] + x5*x530 + x530*M[6] + x545*M[8] + x546*z + x547*z + x548*z + x550*y + x573*M[5] + x631*x648 + x96*M[62] + y*M[138] + z*M[137] + M[189];
#pragma omp atomic
Ms[190] += x*M[145] + x11*x451 + x11*x452 + x11*x455 + x11*M[109] + x112*x184 + x113*x184 + x115*x184 + x121*x160 + x121*M[47] + x184*M[38] + x2*x577 + x208*x93 + x208*M[32] + x221*x96 + x222*x96 + x224*x96 + x230*x76 + x230*M[27] + x24*x385 + x307*x45 + x308*x45 + x309*x45 + x310*x46 + x310*x47 + x310*x48 + x310*M[20] + x355*x38 + x355*M[17] + x376*x41 + x377*x41 + x379*x41 + x385*M[13] + x41*M[96] + x45*M[81] + x530*x6 + x530*M[7] + x550*z + x551*z + x552*z + x553*y + x577*M[4] + x649*M[1] + x96*M[63] + y*M[139] + z*M[138] + M[190];
#pragma omp atomic
Ms[191] += x*M[146] + x11*x457 + x11*M[110] + x118*x184 + x121*M[48] + x184*M[39] + x208*M[33] + x227*x96 + x230*M[28] + x310*x49 + x310*M[21] + x311*x45 + x355*M[18] + x382*x41 + x385*M[14] + x41*M[97] + x45*M[82] + x459*x643 + x483*x9 + x483*M[9] + x530*M[8] + x553*z + x577*M[5] + x649*M[2] + x96*M[64] + z*M[139] + M[191];
#pragma omp atomic
Ms[192] += x*M[147] + x0*x610 + x11*x460 + x11*M[111] + x12*x461 + x122*x170 + x125*M[49] + x170*M[40] + x231*x86 + x234*M[29] + x292*x52 + x292*M[22] + x34*x386 + x34*M[98] + x389*M[15] + x461*M[10] + x554*y + x557*M[6] + x610*M[3] + x650*M[1] + x86*M[65] + y*M[140] + M[192];
#pragma omp atomic
Ms[193] += x*M[148] + x11*x462 + x11*x463 + x11*x464 + x11*M[112] + x125*x172 + x125*M[50] + x126*x170 + x128*x170 + x130*x170 + x15*x461 + x16*x461 + x17*x461 + x170*M[41] + x234*x88 + x234*M[30] + x235*x86 + x237*x86 + x239*x86 + x292*x55 + x292*x57 + x292*x59 + x292*M[23] + x34*x390 + x34*x392 + x34*x394 + x34*M[99] + x36*x389 + x389*M[16] + x461*M[11] + x554*z + x555*z + x556*z + x557*x7 + x557*M[7] + x558*y + x650*M[2] + x86*M[66] + y*M[141] + z*M[140] + M[193];
#pragma omp atomic
Ms[194] += x*M[149] + x0*x613 + x11*x465 + x11*x466 + x11*x467 + x11*M[113] + x12*x468 + x121*x169 + x121*M[49] + x122*x177 + x125*x175 + x125*M[51] + x132*x170 + x134*x170 + x136*x170 + x170*M[42] + x177*M[40] + x234*x91 + x234*M[31] + x241*x86 + x243*x86 + x245*x86 + x292*x61 + x292*x62 + x292*x63 + x292*M[24] + x299*x52 + x299*M[22] + x34*x396 + x34*x398 + x34*x400 + x34*M[100] + x386*x41 + x387*x41 + x388*x41 + x389*x39 + x389*M[17] + x402*M[15] + x41*M[98] + x468*M[10] + x557*x8 + x557*M[8] + x558*z + x559*z + x560*z + x561*y + x569*M[6] + x613*M[3] + x637*x651 + x86*M[67] + y*M[142] + z*M[141] + M[194];
#pragma omp atomic
Ms[195] += x*M[150] + x0*x615 + x11*x469 + x11*x470 + x11*x471 + x11*M[114] + x12*x472 + x121*x171 + x121*M[50] + x125*x179 + x125*M[52] + x126*x177 + x138*x170 + x139*x170 + x140*x170 + x15*x468 + x170*M[43] + x177*M[41] + x230*x85 + x230*M[29] + x231*x96 + x232*x96 + x233*x96 + x234*x94 + x234*M[32] + x247*x86 + x249*x86 + x251*x86 + x299*x55 + x299*M[23] + x303*x52 + x303*M[22] + x34*x403 + x34*x405 + x34*x407 + x34*M[101] + x389*x42 + x389*M[18] + x390*x41 + x391*x41 + x393*x41 + x402*M[16] + x41*M[99] + x468*M[11] + x472*M[10] + x561*z + x562*z + x563*z + x564*x646 + x565*y + x569*M[7] + x573*M[6] + x615*M[3] + x641*x651 + x86*M[68] + x96*M[65] + y*M[143] + z*M[142] + M[195];
#pragma omp atomic
Ms[196] += x*M[151] + x0*x617 + x11*x473 + x11*x474 + x11*x475 + x11*M[115] + x12*x476 + x121*x174 + x121*M[51] + x122*x184 + x123*x184 + x124*x184 + x125*x182 + x125*M[53] + x132*x177 + x15*x472 + x177*M[42] + x184*M[40] + x230*x87 + x230*M[30] + x234*x97 + x234*M[33] + x235*x96 + x236*x96 + x238*x96 + x253*x86 + x254*x86 + x255*x86 + x299*x61 + x299*M[24] + x303*x55 + x303*M[23] + x33*x385 + x34*x409 + x34*x411 + x34*x413 + x34*M[102] + x385*M[15] + x396*x41 + x397*x41 + x399*x41 + x402*M[17] + x41*M[100] + x472*M[11] + x476*M[10] + x565*z + x566*z + x567*z + x569*M[8] + x570*y + x573*M[7] + x617*M[3] + x637*x652 + x646*x647 + x86*M[69] + x96*M[66] + y*M[144] + z*M[143] + M[196];
#pragma omp atomic
Ms[197] += x*M[152] + x0*x619 + x11*x477 + x11*x478 + x11*x479 + x11*M[116] + x121*x178 + x121*M[52] + x125*x185 + x125*M[54] + x126*x184 + x127*x184 + x129*x184 + x138*x177 + x15*x476 + x177*M[43] + x184*M[41] + x230*x90 + x230*M[31] + x241*x96 + x242*x96 + x244*x96 + x303*x61 + x303*M[24] + x310*x52 + x310*x53 + x310*x54 + x310*M[22] + x34*x415 + x34*x416 + x34*x417 + x34*M[103] + x35*x385 + x385*M[16] + x402*M[18] + x403*x41 + x404*x41 + x406*x41 + x41*M[101] + x476*M[11] + x5*x577 + x570*z + x571*z + x572*z + x573*M[8] + x574*y + x577*M[6] + x619*M[3] + x641*x652 + x96*M[67] + y*M[145] + z*M[144] + M[197];
#pragma omp atomic
Ms[198] += x*M[153] + x11*x480 + x11*x481 + x11*x482 + x11*M[117] + x12*x483 + x121*x181 + x121*M[53] + x13*x483 + x132*x184 + x133*x184 + x135*x184 + x14*x483 + x184*M[42] + x230*x93 + x230*M[32] + x247*x96 + x248*x96 + x250*x96 + x310*x55 + x310*x56 + x310*x58 + x310*M[23] + x38*x385 + x385*M[17] + x409*x41 + x41*x410 + x41*x412 + x41*M[102] + x483*M[10] + x574*z + x575*z + x576*z + x577*x6 + x577*M[7] + x578*y + x653*M[1] + x96*M[68] + y*M[146] + z*M[145] + M[198];
#pragma omp atomic
Ms[199] += x*M[154] + x0*x622 + x11*x484 + x11*M[118] + x121*M[54] + x138*x184 + x15*x483 + x184*M[43] + x230*M[33] + x253*x96 + x310*x61 + x310*M[24] + x385*M[18] + x41*x415 + x41*M[103] + x483*M[11] + x577*M[8] + x578*z + x622*M[3] + x653*M[2] + x96*M[69] + z*M[146] + M[199];
#pragma omp atomic
Ms[200] += x*M[155] + x1*x610 + x141*x170 + x170*M[44] + x18*x461 + x256*x86 + x292*x64 + x292*M[25] + x34*x418 + x34*M[104] + x461*M[12] + x579*y + x610*M[4] + x654*x655 + x86*M[70] + y*M[147] + M[200];
#pragma omp atomic
Ms[201] += x*M[156] + x146*x170 + x148*x170 + x150*x170 + x170*M[45] + x23*x461 + x25*x461 + x261*x86 + x263*x86 + x265*x86 + x27*x461 + x292*x69 + x292*x71 + x292*x73 + x292*M[26] + x3*x610 + x34*x423 + x34*x425 + x34*x427 + x34*M[105] + x4*x610 + x461*M[13] + x579*z + x580*z + x581*z + x584*y + x610*M[5] + x86*M[71] + y*M[148] + z*M[147] + M[201];
#pragma omp atomic
Ms[202] += x*M[157] + x1*x613 + x141*x177 + x152*x170 + x154*x170 + x156*x170 + x170*M[46] + x177*M[44] + x18*x468 + x267*x86 + x269*x86 + x271*x86 + x29*x461 + x292*x75 + x292*x77 + x292*x79 + x292*M[27] + x299*x64 + x299*M[25] + x30*x461 + x31*x461 + x34*x429 + x34*x431 + x34*x433 + x34*M[106] + x41*x418 + x41*x419 + x41*x420 + x41*M[104] + x461*M[14] + x468*M[12] + x584*z + x585*z + x586*z + x587*y + x590*x656 + x613*M[4] + x86*M[72] + y*M[149] + z*M[148] + M[202];
#pragma omp atomic
Ms[203] += x*M[158] + x1*x615 + x146*x177 + x159*x170 + x161*x170 + x163*x170 + x170*M[47] + x177*M[45] + x18*x472 + x23*x468 + x256*x96 + x257*x96 + x258*x96 + x274*x86 + x276*x86 + x278*x86 + x292*x81 + x292*x82 + x292*x83 + x292*M[28] + x299*x69 + x299*M[26] + x3*x613 + x303*x64 + x303*M[25] + x34*x436 + x34*x438 + x34*x440 + x34*M[107] + x41*x423 + x41*x424 + x41*x426 + x41*M[105] + x422*x657 + x468*M[13] + x472*M[12] + x587*z + x588*z + x589*z + x592*y + x613*M[5] + x615*M[4] + x86*M[73] + x96*M[70] + y*M[150] + z*M[149] + M[203];
#pragma omp atomic
Ms[204] += x*M[159] + x1*x617 + x141*x184 + x142*x184 + x143*x184 + x152*x177 + x165*x170 + x166*x170 + x167*x170 + x170*M[48] + x177*M[46] + x18*x476 + x184*M[44] + x23*x472 + x261*x96 + x262*x96 + x264*x96 + x281*x86 + x283*x86 + x285*x86 + x29*x468 + x299*x75 + x299*M[27] + x3*x615 + x303*x69 + x303*M[26] + x34*x443 + x34*x445 + x34*x447 + x34*M[108] + x360*x658 + x41*x429 + x41*x430 + x41*x432 + x41*M[106] + x468*M[14] + x472*M[13] + x476*M[12] + x592*z + x593*z + x594*z + x596*y + x615*M[5] + x617*M[4] + x86*M[74] + x96*M[71] + y*M[151] + z*M[150] + M[204];
#pragma omp atomic
Ms[205] += x*M[160] + x1*x619 + x146*x184 + x147*x184 + x149*x184 + x159*x177 + x177*M[47] + x184*M[45] + x23*x476 + x267*x96 + x268*x96 + x270*x96 + x287*x86 + x288*x86 + x289*x86 + x29*x472 + x299*x81 + x299*M[28] + x3*x617 + x303*x75 + x303*M[27] + x310*x64 + x310*x65 + x310*x66 + x310*M[25] + x34*x450 + x34*x452 + x34*x454 + x34*M[109] + x41*x436 + x41*x437 + x41*x439 + x41*M[107] + x422*x659 + x472*M[14] + x476*M[13] + x596*z + x597*z + x598*z + x600*y + x617*M[5] + x619*M[4] + x86*M[75] + x96*M[72] + y*M[152] + z*M[151] + M[205];
#pragma omp atomic
Ms[206] += x*M[161] + x111*x660 + x152*x184 + x153*x184 + x155*x184 + x165*x177 + x177*M[48] + x18*x483 + x184*M[46] + x19*x483 + x20*x483 + x274*x96 + x275*x96 + x277*x96 + x29*x476 + x3*x619 + x303*x81 + x303*M[28] + x310*x69 + x310*x70 + x310*x72 + x310*M[26] + x34*x456 + x34*x457 + x34*x458 + x34*M[110] + x41*x443 + x41*x444 + x41*x446 + x41*M[108] + x476*M[14] + x483*M[12] + x600*z + x601*z + x602*z + x604*y + x619*M[5] + x96*M[73] + y*M[153] + z*M[152] + M[206];
#pragma omp atomic
Ms[207] += x*M[162] + x1*x622 + x159*x184 + x160*x184 + x162*x184 + x184*M[47] + x2*x622 + x23*x483 + x24*x483 + x26*x483 + x281*x96 + x282*x96 + x284*x96 + x310*x75 + x310*x76 + x310*x78 + x310*M[27] + x41*x450 + x41*x451 + x41*x453 + x41*M[109] + x483*M[13] + x604*z + x605*z + x606*z + x607*y + x622*M[4] + x96*M[74] + y*M[154] + z*M[153] + M[207];
#pragma omp atomic
Ms[208] += x*M[163] + x165*x184 + x184*M[48] + x287*x96 + x29*x483 + x3*x622 + x310*x81 + x310*M[28] + x41*x456 + x41*M[110] + x483*M[14] + x607*z + x622*M[5] + x655*x661 + x96*M[75] + z*M[154] + M[208];
#pragma omp atomic
Ms[209] += x170*M[49] + x292*M[29] + x34*M[111] + x461*M[15] + x610*M[6] + x662*M[1] + x86*M[76] + y*M[155] + M[209];
#pragma omp atomic
Ms[210] += x170*x172 + x170*M[50] + x292*x88 + x292*M[30] + x294*x86 + x34*x463 + x34*M[112] + x36*x461 + x461*M[16] + x609*z + x610*x7 + x610*M[7] + x662*M[2] + x86*M[77] + y*M[156] + z*M[155] + M[210];
#pragma omp atomic
Ms[211] += x170*x175 + x170*M[51] + x177*M[49] + x292*x91 + x292*M[31] + x297*x86 + x299*M[29] + x34*x466 + x34*M[113] + x39*x461 + x41*x460 + x41*M[111] + x461*M[17] + x468*M[15] + x610*x8 + x610*M[8] + x611*z + x613*M[6] + x663*M[1] + x86*M[78] + y*M[157] + z*M[156] + M[211];
#pragma omp atomic
Ms[212] += x170*x179 + x170*M[52] + x177*M[50] + (1.0/720.0)*x259*x568 + x291*x96 + x292*x94 + x292*M[32] + x299*M[30] + x301*x86 + x303*M[29] + x34*x470 + x34*M[114] + x41*x462 + x41*M[112] + x42*x461 + x461*M[18] + x468*M[16] + x472*M[15] + x612*z + x613*M[7] + x615*M[6] + x663*M[2] + x86*M[79] + x96*M[76] + y*M[158] + z*M[157] + M[212];
#pragma omp atomic
Ms[213] += x169*x184 + x170*x182 + x170*M[53] + x177*M[51] + x184*M[49] + x292*x97 + x292*M[33] + x293*x96 + x299*M[31] + x303*M[30] + x305*x86 + x34*x474 + x34*M[115] + x41*x465 + x41*M[113] + x468*M[17] + x472*M[16] + x476*M[15] + x613*M[8] + x614*z + x615*M[7] + x617*M[6] + x657*x665 + x664*M[1] + x86*M[80] + x96*M[77] + y*M[159] + z*M[158] + M[213];
#pragma omp atomic
Ms[214] += x170*x185 + x170*M[54] + x171*x184 + x177*M[52] + x184*M[50] + (1.0/720.0)*x290*x564 + x296*x96 + x299*M[32] + x303*M[31] + x308*x86 + x310*x85 + x310*M[29] + x34*x478 + x34*M[116] + x41*x469 + x41*M[114] + x468*M[18] + x472*M[17] + x476*M[16] + x615*M[8] + x616*z + x617*M[7] + x619*M[6] + x664*M[2] + x86*M[81] + x96*M[78] + y*M[160] + z*M[159] + M[214];
#pragma omp atomic
Ms[215] += x174*x184 + x177*M[53] + x184*M[51] + x299*M[33] + x300*x96 + x303*M[32] + x310*x87 + x310*M[30] + x311*x86 + x33*x483 + x34*x481 + x34*M[117] + x41*x473 + x41*M[115] + x472*M[18] + x476*M[17] + x483*M[15] + x617*M[8] + x618*z + x619*M[7] + x659*x665 + x666*M[1] + x86*M[82] + x96*M[79] + y*M[161] + z*M[160] + M[215];
#pragma omp atomic
Ms[216] += x177*M[54] + x178*x184 + x184*M[52] + x303*M[33] + x304*x96 + x310*x90 + x310*M[31] + x34*x484 + x34*M[118] + x35*x483 + x41*x477 + x41*M[116] + x476*M[18] + x483*M[16] + x5*x622 + x619*M[8] + x620*z + x622*M[6] + x666*M[2] + x96*M[80] + y*M[162] + z*M[161] + M[216];
#pragma omp atomic
Ms[217] += x181*x184 + x184*M[53] + x307*x96 + x310*x93 + x310*M[32] + x38*x483 + x41*x480 + x41*M[117] + x483*M[17] + x6*x622 + x621*z + x622*M[7] + x667*M[1] + x96*M[81] + y*M[163] + z*M[162] + M[217];
#pragma omp atomic
Ms[218] += x184*M[54] + x310*M[33] + x41*M[118] + x483*M[18] + x622*M[8] + x667*M[2] + x96*M[82] + z*M[163] + M[218];

}

void M2L_9(double x, double y, double z, double * M, double * L) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
double x448;
double x449;
double x450;
double x451;
double x452;
double x453;
double x454;
double x455;
double x456;
double x457;
double x458;
double x459;
double x460;
double x461;
double x462;
double x463;
double x464;
double x465;
double x466;
double x467;
double x468;
double x469;
double x470;
double x471;
double x472;
double x473;
double x474;
double x475;
double x476;
double x477;
double x478;
double x479;
double x480;
double x481;
double x482;
double x483;
double x484;
double x485;
double x486;
double x487;
double x488;
double x489;
double x490;
double x491;
double x492;
double x493;
double x494;
double x495;
double x496;
double x497;
double x498;
double x499;
double x500;
double x501;
double x502;
double x503;
double x504;
double x505;
double x506;
double x507;
double x508;
double x509;
double x510;
double x511;
double x512;
double x513;
double x514;
double x515;
double x516;
double x517;
double x518;
double x519;
double x520;
double x521;
double x522;
double x523;
double x524;
double x525;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = pow(R, -5);
x3 = 3.0*x2;
x4 = x*y;
x5 = x3*x4;
x6 = x*z;
x7 = x3*x6;
x8 = y*z;
x9 = x3*x8;
x10 = pow(R, -7);
x11 = 15.0*x10;
x12 = x4*z;
x13 = x11*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x3;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x3;
x20 = x14 + x19;
x21 = 9.0*x2;
x22 = x11*x15;
x23 = -x22;
x24 = x*(x21 + x23);
x25 = x23 + x3;
x26 = x25*y;
x27 = x11*x18;
x28 = -x27;
x29 = y*(x21 + x28);
x30 = x25*z;
x31 = z*(x28 + x3);
x32 = 1.0*x;
x33 = x32*(x27 - x3);
x34 = 45.0*x10;
x35 = -x34;
x36 = pow(R, -9);
x37 = x15*x36;
x38 = 105.0*x37;
x39 = x35 + x38;
x40 = x39*x4;
x41 = x39*x6;
x42 = x18*x36;
x43 = 105.0*x42;
x44 = x35 + x43;
x45 = x44*x8;
x46 = -x11;
x47 = x8*(x38 + x46);
x48 = x32*y;
x49 = x44*x48;
x50 = x43 + x46;
x51 = x32*z;
x52 = x50*x51;
x53 = 315.0*x36;
x54 = pow(R, -11);
x55 = 945.0*x54;
x56 = x15*x55;
x57 = x53 - x56;
x58 = x12*x57;
x59 = x18*x55;
x60 = x32*x8;
x61 = x60*(-x53 + x59);
x62 = (x*x*x*x);
x63 = 105.0*x36;
x64 = x62*x63;
x65 = 90.0*x10;
x66 = -x15*x65 + x21 + x64;
x67 = (y*y*y*y);
x68 = x63*x67;
x69 = -x18*x65 + x21 + x68;
x70 = 2.0*x0 - x16 - x19;
x71 = -225.0*x10;
x72 = x55*x62;
x73 = -x72;
x74 = x*(1050.0*x37 + x71 + x73);
x75 = x55*x67;
x76 = -x75;
x77 = y*(1050.0*x42 + x71 + x76);
x78 = x35 + 630.0*x37 + x73;
x79 = x78*y;
x80 = x78*z;
x81 = 630.0*x42;
x82 = x35 + x76 + x81;
x83 = x82*z;
x84 = x32*(x34 + x75 - x81);
x85 = 1575.0*x36;
x86 = pow(R, -13);
x87 = x62*x86;
x88 = 10395.0*x87;
x89 = x15*x54;
x90 = x85 + x88 - 9450.0*x89;
x91 = x4*x90;
x92 = x6*x90;
x93 = 5670.0*x89;
x94 = x53 + x88 - x93;
x95 = x8*x94;
x96 = x67*x86;
x97 = 10395.0*x96;
x98 = x18*x54;
x99 = x85 + x97 - 9450.0*x98;
x100 = x8*x99;
x101 = x48*x99;
x102 = x53 + x97 - 5670.0*x98;
x103 = x102*x51;
x104 = 14175.0*x54;
x105 = -x104;
x106 = pow(R, -15);
x107 = 135135.0*x106;
x108 = x107*x62;
x109 = x15*x86;
x110 = 103950.0*x109;
x111 = x105 - x108 + x110;
x112 = x111*x12;
x113 = x107*x67;
x114 = x18*x86;
x115 = 103950.0*x114;
x116 = x60*(x104 + x113 - x115);
x117 = pow(x, 6);
x118 = 10395.0*x86;
x119 = x117*x118;
x120 = -x104*x62 + x119 + 4725.0*x37 + x71;
x121 = pow(y, 6);
x122 = x118*x121;
x123 = -x104*x67 + x122 + 4725.0*x42 + x71;
x124 = 11025.0*x36;
x125 = 99225.0*x54;
x126 = x107*x117;
x127 = -x126;
x128 = x127 + 218295.0*x87;
x129 = x*(x124 - x125*x15 + x128);
x130 = 42525.0*x54;
x131 = x127 - x130*x15 + x85 + 155925.0*x87;
x132 = x131*y;
x133 = x107*x121;
x134 = -x133;
x135 = x134 + 218295.0*x96;
x136 = y*(x124 - x125*x18 + x135);
x137 = x131*z;
x138 = 155925.0*x96;
x139 = x130*x18;
x140 = x134 + x138 - x139 + x85;
x141 = x140*z;
x142 = x32*(x133 - x138 + x139 - x85);
x143 = pow(R, -17);
x144 = x117*x143;
x145 = 2027025.0*x144;
x146 = x106*x62;
x147 = 2837835.0*x146;
x148 = -x125;
x149 = 1091475.0*x86;
x150 = x148 + x149*x15;
x151 = x145 - x147 + x150;
x152 = x151*x4;
x153 = x151*x6;
x154 = x149*x18;
x155 = x121*x143;
x156 = 2027025.0*x155;
x157 = x106*x67;
x158 = 2837835.0*x157;
x159 = x154 + x156 - x158;
x160 = x148 + x159;
x161 = x160*x8;
x162 = 2027025.0*x146;
x163 = -x162;
x164 = 467775.0*x86;
x165 = x15*x164;
x166 = x105 + x145 + x163 + x165;
x167 = x166*x8;
x168 = x18*x38;
x169 = x168 + x25 + x28;
x170 = x160*x48;
x171 = 2027025.0*x157;
x172 = -x171;
x173 = x105 + x156 + x164*x18 + x172;
x174 = x173*x51;
x175 = 34459425.0/pow(R, 19);
x176 = x117*x175;
x177 = x143*x62;
x178 = 42567525.0*x177;
x179 = x106*x15;
x180 = 14189175.0*x179;
x181 = x12*(x149 - x176 + x178 - x180);
x182 = x18*x53;
x183 = -x18*x56;
x184 = x*(x182 + x183 + x39);
x185 = x15*x53;
x186 = y*(x183 + x185 + x44);
x187 = z*(x183 + x38 + x50);
x188 = x121*x175;
x189 = x143*x67;
x190 = 42567525.0*x189;
x191 = x106*x18;
x192 = x60*(-x149 + x188 - x190 + 14189175.0*x191);
x193 = -x24 + x33;
x194 = -x26 - x29;
x195 = -x30 - x31;
x196 = -2835.0*x98;
x197 = x109*x18;
x198 = 10395.0*x197;
x199 = x196 + x198;
x200 = 945.0*x36;
x201 = -2835.0*x89;
x202 = x200 + x201;
x203 = x4*(x199 + x202);
x204 = x6*(x199 + x57);
x205 = x8*(x198 + x201 + x53 - x59);
x206 = 31185.0*x114;
x207 = x15*x18;
x208 = -8505.0*x54;
x209 = 31185.0*x109;
x210 = x208 + x209;
x211 = x12*(-x107*x207 + x206 + x210);
x212 = -x40 - x49;
x213 = -x41 - x52;
x214 = -x45 - x47;
x215 = pow(x, 8);
x216 = 2027025.0*x143;
x217 = x215*x216;
x218 = 3783780.0*x106;
x219 = -x117*x218 + x124 + x217 + 2182950.0*x87 - 396900.0*x89;
x220 = pow(y, 8);
x221 = x216*x220;
x222 = -x121*x218 + x124 + x221 + 2182950.0*x96 - 396900.0*x98;
x223 = -x58 + x61;
x224 = -893025.0*x54;
x225 = -x175*x215;
x226 = x*(13097700.0*x109 + 72972900.0*x144 - 51081030.0*x146 + x224 + x225);
x227 = x175*x220;
x228 = -x227;
x229 = y*(13097700.0*x114 + 72972900.0*x155 - 51081030.0*x157 + x224 + x228);
x230 = 4365900.0*x109 + 56756700.0*x144 - 28378350.0*x146 + x148 + x225;
x231 = x230*y;
x232 = x230*z;
x233 = 56756700.0*x155;
x234 = 28378350.0*x157;
x235 = 4365900.0*x114;
x236 = z*(x148 + x228 + x233 - x234 + x235);
x237 = x32*(x125 + x227 - x233 + x234 - x235);
x238 = 105.0*x10;
x239 = -x168 - 12.0*x2;
x240 = x18*x238 + x22 + x239 - x68;
x241 = x15*x238 + x239 + x27 - x64;
x242 = 120.0*x10;
x243 = -x15*x242 - x18*x242 + 210.0*x18*x37 + 24.0*x2 + x64 + x68;
x244 = x109*x67;
x245 = 10395.0*x244;
x246 = x18*x93;
x247 = -x246;
x248 = x185 + x245 + x247 + x82;
x249 = x18*x88;
x250 = x182 + x247 + x249 + x78;
x251 = 62370.0*x197;
x252 = -x113*x15;
x253 = x251 + x252;
x254 = 31185.0*x96;
x255 = x254 - 17010.0*x98;
x256 = x*(x202 + x253 + x255);
x257 = x104*x18;
x258 = -x257;
x259 = x110*x18;
x260 = -x108*x18;
x261 = x*(x258 + x259 + x260 + x90);
x262 = x196 + x251 + x260;
x263 = 31185.0*x87;
x264 = 17010.0*x89;
x265 = x200 + x263 - x264;
x266 = y*(x262 + x265);
x267 = x104*x15;
x268 = -x267;
x269 = y*(x252 + x259 + x268 + x99);
x270 = z*(x262 + x94);
x271 = z*(x102 + x201 + x253);
x272 = 155925.0*x109;
x273 = 311850.0*x114;
x274 = -x130;
x275 = x179*x18;
x276 = -1351350.0*x275;
x277 = x274 + x276;
x278 = -405405.0*x157;
x279 = x15*x189;
x280 = 2027025.0*x279;
x281 = x278 + x280;
x282 = x4*(x272 + x273 + x277 + x281);
x283 = 155925.0*x114;
x284 = 311850.0*x109;
x285 = -405405.0*x146;
x286 = x177*x18;
x287 = 2027025.0*x286;
x288 = x285 + x287;
x289 = x4*(x277 + x283 + x284 + x288);
x290 = x6*(x111 + x276 + x283 + x287);
x291 = 187110.0*x114;
x292 = -810810.0*x275;
x293 = x6*(x210 + x281 + x291 + x292);
x294 = x8*(x105 - x113 + x115 + x272 + x276 + x280);
x295 = 187110.0*x109 + x208;
x296 = x8*(x206 + x288 + x292 + x295);
x297 = 6081075.0*x177;
x298 = x175*x62;
x299 = 20270250.0*x143*x207 + x164;
x300 = x12*(-4054050.0*x179 - x18*x298 - 2027025.0*x191 + x297 + x299);
x301 = x12*(-x15*x175*x67 - 2027025.0*x179 + 6081075.0*x189 - 4054050.0*x191 + x299);
x302 = 15120.0*x54;
x303 = -x119;
x304 = -x249;
x305 = 270.0*x10 + x246;
x306 = -x182 + x302*x62 + x303 + x304 + x305 - 5355.0*x37;
x307 = -x122;
x308 = -x245;
x309 = -x185 + x302*x67 + x305 + x307 + x308 - 5355.0*x42;
x310 = -x184;
x311 = x310 - x74;
x312 = x310 + x84;
x313 = -x186;
x314 = x313 - x77;
x315 = x313 - x79;
x316 = -x187;
x317 = x316 - x80;
x318 = x316 - x83;
x319 = -x203;
x320 = x319 - x91;
x321 = -x101 + x319;
x322 = -x204;
x323 = -x103 + x322;
x324 = x322 - x92;
x325 = -x205;
x326 = x325 - x95;
x327 = -x100 + x325;
x328 = -x211;
x329 = -x112 + x328;
x330 = x116 + x328;
x331 = x165*x18;
x332 = x145*x18;
x333 = x162*x18;
x334 = x131 + x258 + x331 + x332 - x333;
x335 = x15*x156;
x336 = x15*x171;
x337 = x140 + x268 + x331 + x335 - x336;
x338 = x18*x89;
x339 = -x15*x200 - x18*x200 + x304 + x308 + 11340.0*x338 + x65 + x72 + x75;
x340 = -x18*x180;
x341 = -x176*x18;
x342 = x*(x151 + x154 + x178*x18 + x340 + x341);
x343 = 6081075.0*x155;
x344 = -6081075.0*x275;
x345 = x274 + x344;
x346 = -x15*x188;
x347 = x272 + 30405375.0*x279 + x346;
x348 = x*(1403325.0*x114 - 6081075.0*x157 + x343 + x345 + x347);
x349 = y*(x15*x190 + x150 + x159 + x340 + x346);
x350 = 6081075.0*x144;
x351 = x283 + 30405375.0*x286 + x341;
x352 = y*(1403325.0*x109 - 6081075.0*x146 + x345 + x350 + x351);
x353 = z*(x173 + x344 + x347);
x354 = z*(x166 + x344 + x351);
x355 = 3918915.0*x106;
x356 = -x221;
x357 = -x335;
x358 = -x331 - 12600.0*x36;
x359 = x121*x355 + x267 + x336 + x356 + x357 + x358 - 2338875.0*x96 + 439425.0*x98;
x360 = -x217;
x361 = -x332;
x362 = x117*x355 + x257 + x333 + x358 + x360 + x361 - 2338875.0*x87 + 439425.0*x89;
x363 = 16065.0*x54;
x364 = -360.0*x10 - x18*x264;
x365 = x122 + 20790.0*x244 + x249 - x363*x67 + x364 + 1260.0*x37 + 6300.0*x42 + x73;
x366 = x119 + 20790.0*x18*x87 + x245 - x363*x62 + x364 + 6300.0*x37 + 1260.0*x42 + x76;
x367 = x15*x157;
x368 = 810810.0*x367;
x369 = -x368;
x370 = x177*x67;
x371 = 2027025.0*x370;
x372 = x146*x18;
x373 = 810810.0*x372;
x374 = x371 - x373;
x375 = 374220.0*x197 + x255 + x265 + x369 + x374;
x376 = -x298*x67;
x377 = x274 - 8108100.0*x275 + x376;
x378 = x285 + 12162150.0*x286;
x379 = x*(935550.0*x114 + x172 + 20270250.0*x279 + x284 + x377 + x378);
x380 = x278 + 12162150.0*x279;
x381 = y*(935550.0*x109 + x163 + x273 + 20270250.0*x286 + x377 + x380);
x382 = z*(-4864860.0*x275 + x291 + x295 + x376 + x378 + x380);
x383 = 2*x186 + x77 + x79;
x384 = 2*x187 + x80 + x83;
x385 = 2*x184 + x74 - x84;
x386 = 17010.0*x54;
x387 = 720.0*x10 - x18*x263 - x209*x67 + x303 + x307 + 34020.0*x338 - 7560.0*x37 + x386*x62 + x386*x67 - 7560.0*x42;
x388 = x100 + 2*x205 + x95;
x389 = -x256;
x390 = x142 + x389;
x391 = -x129;
x392 = -x261;
x393 = x391 + x392;
x394 = -x132;
x395 = -x266;
x396 = x394 + x395;
x397 = -x136;
x398 = -x269;
x399 = x397 + x398;
x400 = -x137;
x401 = -x270;
x402 = x400 + x401;
x403 = -x141;
x404 = -x271;
x405 = x403 + x404;
x406 = x101 + 2*x203 + x91;
x407 = x103 + 2*x204 + x92;
x408 = -x152;
x409 = -x289;
x410 = x408 + x409;
x411 = -x170;
x412 = -x282;
x413 = x411 + x412;
x414 = -x153;
x415 = -x290;
x416 = x414 + x415;
x417 = -x174;
x418 = -x293;
x419 = x417 + x418;
x420 = -x161;
x421 = -x294;
x422 = x420 + x421;
x423 = -x167;
x424 = -x296;
x425 = x423 + x424;
x426 = x112 - x116 + 2*x211;
x427 = -x301;
x428 = x192 + x427;
x429 = -x181;
x430 = -x300;
x431 = x429 + x430;
x432 = x133 + x357;
x433 = -841995.0*x197 - 2520.0*x36 - x371;
x434 = x15*x158 - x263 + x373 + x432 + x433 + 31185.0*x89 - 187110.0*x96 + 59535.0*x98;
x435 = x126 + x361;
x436 = x147*x18 - x254 + x368 + x433 + x435 - 187110.0*x87 + 59535.0*x89 + 31185.0*x98;
x437 = 4054050.0*x106;
x438 = x144*x18;
x439 = 1309770.0*x197 + 15120.0*x36;
x440 = -x117*x437 + x217 + x254 + x369 + x371 - 4864860.0*x372 + 4054050.0*x438 + x439 + 2525985.0*x87 - 498960.0*x89 - 45360.0*x98;
x441 = x15*x155;
x442 = -x121*x437 + x221 + x263 - 4864860.0*x367 + x374 + x439 + 4054050.0*x441 - 45360.0*x89 + 2525985.0*x96 - 498960.0*x98;
x443 = x389 + x392;
x444 = x395 + x398;
x445 = x401 + x404;
x446 = x409 + x412;
x447 = x415 + x418;
x448 = x421 + x424;
x449 = x427 + x430;
x450 = -x226;
x451 = -x342;
x452 = x450 + x451;
x453 = -x348;
x454 = x237 + x453;
x455 = -x229;
x456 = -x349;
x457 = x455 + x456;
x458 = -x231;
x459 = -x352;
x460 = x458 + x459;
x461 = -x232;
x462 = -x354;
x463 = x461 + x462;
x464 = -x236;
x465 = -x353;
x466 = x464 + x465;
x467 = x128 + x135 + 1683990.0*x197 + x332 + x335 + 5040.0*x36 - 3648645.0*x367 + 4054050.0*x370 - 3648645.0*x372 - 90720.0*x89 - 90720.0*x98;
x468 = 4189185.0*x106;
x469 = -2993760.0*x197 - x297*x67 - 20160.0*x36;
x470 = x121*x468 - x15*x343 + x356 + 8513505.0*x367 + 4459455.0*x372 + x435 + x469 - 249480.0*x87 + 136080.0*x89 - 2744280.0*x96 + 589680.0*x98;
x471 = x117*x468 - x18*x350 + x360 + 4459455.0*x367 + 8513505.0*x372 + x432 + x469 - 2744280.0*x87 + 589680.0*x89 - 249480.0*x96 + 136080.0*x98;
x472 = 4324320.0*x106;
x473 = -x117*x472 - x121*x472 + 5987520.0*x197 + x217 + x221 + 40320.0*x36 - 12972960.0*x367 + 12162150.0*x370 - 12972960.0*x372 + 8108100.0*x438 + 8108100.0*x441 + 2993760.0*x87 - 725760.0*x89 + 2993760.0*x96 - 725760.0*x98;
x474 = x129 + x256 + 2*x261;
x475 = x132 + 2*x266 + x269;
x476 = x136 + x266 + 2*x269;
x477 = x137 + 2*x270 + x271;
x478 = x141 + x270 + 2*x271;
x479 = -x142 + 2*x256 + x261;
x480 = x152 + x282 + 2*x289;
x481 = x153 + 2*x290 + x293;
x482 = x161 + 2*x294 + x296;
x483 = x167 + x294 + 2*x296;
x484 = x170 + 2*x282 + x289;
x485 = x174 + x290 + 2*x293;
x486 = x181 + 2*x300 + x301;
x487 = -x192 + x300 + 2*x301;
x488 = -x379;
x489 = x451 + x488;
x490 = x453 + x488;
x491 = -x381;
x492 = x456 + x491;
x493 = x459 + x491;
x494 = -x382;
x495 = x465 + x494;
x496 = x462 + x494;
x497 = x142 - 3*x256 - 3*x261 + x391;
x498 = -3*x266 - 3*x269 + x394 + x397;
x499 = -3*x270 - 3*x271 + x400 + x403;
x500 = -3*x282 - 3*x289 + x408 + x411;
x501 = -3*x290 - 3*x293 + x414 + x417;
x502 = -3*x294 - 3*x296 + x420 + x423;
x503 = x192 - 3*x300 - 3*x301 + x429;
x504 = x226 + 2*x342 + x379;
x505 = x229 + 2*x349 + x381;
x506 = x231 + 2*x352 + x381;
x507 = x232 + 2*x354 + x382;
x508 = x236 + 2*x353 + x382;
x509 = -x237;
x510 = 2*x348 + x379 + x509;
x511 = x342 + x348 + 2*x379;
x512 = x349 + x352 + 2*x381;
x513 = x353 + x354 + 2*x382;
x514 = -3*x379;
x515 = -3*x342 + x450 + x453 + x514;
x516 = x237 - 3*x348 + x451 + x514;
x517 = -3*x381;
x518 = -3*x349 + x455 + x459 + x517;
x519 = -3*x352 + x456 + x458 + x517;
x520 = -3*x382;
x521 = -3*x354 + x461 + x465 + x520;
x522 = -3*x353 + x462 + x464 + x520;
x523 = x229 + x231 + 4*x349 + 4*x352 + 6*x381;
x524 = x232 + x236 + 4*x353 + 4*x354 + 6*x382;
x525 = x226 + 4*x342 + 4*x348 + 6*x379 + x509;
#pragma omp atomic
L[0] += -x*x1*M[0] - x1*y*M[1] - x1*z*M[2] + x100*M[77] + x101*M[70] + x103*M[71] + x112*M[87] - x116*M[105] + x120*M[55] + x123*M[76] + x129*M[83] - x13*M[13] + x132*M[84] + x136*M[111] + x137*M[85] + x141*M[112] - x142*M[104] + x152*M[120] + x153*M[121] + x161*M[156] + x167*M[123] + x169*M[22] + x17*M[3] + x170*M[147] + x174*M[148] + x181*M[168] + x184*M[37] + x186*M[40] + x187*M[41] - x192*M[201] + x193*M[14] + x194*M[17] + x195*M[18] + x20*M[6] + x203*M[61] + x204*M[62] + x205*M[66] + x211*M[94] + x212*M[27] + x213*M[28] + x214*M[32] + x219*M[119] + x222*M[155] + x223*M[47] + x226*M[164] + x229*M[209] + x231*M[165] + x232*M[166] + x236*M[210] - x237*M[200] + x24*M[9] + x240*M[31] + x241*M[24] + x243*M[33] + x248*M[65] + x250*M[58] + x256*M[93] + x26*M[10] + x261*M[86] + x266*M[89] + x269*M[98] + x270*M[90] + x271*M[99] + x282*M[134] + x289*M[125] + x29*M[15] + x290*M[126] + x293*M[135] + x294*M[141] + x296*M[130] + x30*M[11] + x300*M[175] + x301*M[186] + x306*M[60] + x309*M[78] + x31*M[16] + x311*M[39] + x312*M[46] + x314*M[51] + x315*M[42] + x317*M[43] + x318*M[52] + x320*M[63] + x321*M[72] + x323*M[73] + x324*M[64] + x326*M[68] + x327*M[79] + x329*M[96] - x33*M[12] + x330*M[107] + x334*M[122] + x337*M[140] + x339*M[67] + x342*M[167] + x348*M[185] + x349*M[192] + x352*M[170] + x353*M[193] + x354*M[171] + x359*M[157] + x362*M[124] + x365*M[80] + x366*M[69] + x375*M[129] + x379*M[174] + x381*M[179] + x382*M[180] + x383*M[53] + x384*M[54] + x385*M[48] + x387*M[82] + x388*M[81] + x390*M[106] + x393*M[88] + x396*M[91] + x399*M[113] + x40*M[20] + x402*M[92] + x405*M[114] + x406*M[74] + x407*M[75] + x41*M[21] + x410*M[127] + x413*M[149] + x416*M[128] + x419*M[150] + x422*M[158] + x425*M[132] + x426*M[109] + x428*M[203] + x431*M[177] + x434*M[142] + x436*M[131] + x440*M[133] + x442*M[159] + x443*M[95] + x444*M[100] + x445*M[101] + x446*M[136] + x447*M[137] + x448*M[143] + x449*M[188] + x45*M[30] + x452*M[169] + x454*M[202] + x457*M[211] + x460*M[172] + x463*M[173] + x466*M[212] + x467*M[144] + x47*M[23] + x470*M[161] + x471*M[146] + x473*M[163] + x474*M[97] + x475*M[102] + x476*M[115] + x477*M[103] + x478*M[116] + x479*M[108] + x480*M[138] + x481*M[139] + x482*M[160] + x483*M[145] + x484*M[151] + x485*M[152] + x486*M[190] + x487*M[205] + x489*M[176] + x49*M[25] + x490*M[187] + x492*M[194] + x493*M[181] + x495*M[195] + x496*M[182] + x497*M[110] + x498*M[117] + x499*M[118] + x5*M[4] + x500*M[153] + x501*M[154] + x502*M[162] + x503*M[207] + x504*M[178] + x505*M[213] + x506*M[183] + x507*M[184] + x508*M[214] + x510*M[204] + x511*M[189] + x512*M[196] + x513*M[197] + x515*M[191] + x516*M[206] + x518*M[215] + x519*M[198] + x52*M[26] + x521*M[199] + x522*M[216] + x523*M[217] + x524*M[218] + x525*M[208] + x58*M[38] - x61*M[45] + x66*M[19] + x69*M[29] + x7*M[5] + x70*M[8] + x74*M[34] + x77*M[49] + x79*M[35] + x80*M[36] + x83*M[50] - x84*M[44] + x9*M[7] + x91*M[56] + x92*M[57] + x95*M[59];
#pragma omp atomic
L[1] += x101*M[49] + x103*M[50] + x112*M[59] - x116*M[77] + x120*M[34] + x129*M[55] - x13*M[7] + x132*M[56] + x137*M[57] - x142*M[76] + x152*M[84] + x153*M[85] + x167*M[87] + x169*M[12] + x17*M[0] + x170*M[111] + x174*M[112] + x181*M[123] + x184*M[22] + x186*M[25] + x187*M[26] - x192*M[156] + x193*M[8] + x203*M[40] + x204*M[41] + x205*M[45] + x211*M[66] + x212*M[17] + x213*M[18] + x219*M[83] + x223*M[32] + x226*M[119] + x231*M[120] + x232*M[121] - x237*M[155] + x24*M[3] + x241*M[14] + x248*M[44] + x250*M[37] + x256*M[65] + x26*M[4] + x261*M[58] + x266*M[61] + x269*M[70] + x270*M[62] + x271*M[71] + x282*M[98] + x289*M[89] + x290*M[90] + x293*M[99] + x294*M[105] + x296*M[94] + x30*M[5] + x300*M[130] + x301*M[141] + x306*M[39] + x311*M[24] + x312*M[31] + x315*M[27] + x317*M[28] + x320*M[42] + x321*M[51] + x323*M[52] + x324*M[43] + x326*M[47] + x329*M[68] - x33*M[6] + x330*M[79] + x334*M[86] + x337*M[104] + x339*M[46] + x342*M[122] + x348*M[140] + x349*M[147] + x352*M[125] + x353*M[148] + x354*M[126] + x362*M[88] + x366*M[48] + x375*M[93] + x379*M[129] + x381*M[134] + x382*M[135] + x385*M[33] + x390*M[78] + x393*M[60] + x396*M[63] + x40*M[10] + x402*M[64] + x406*M[53] + x407*M[54] + x41*M[11] + x410*M[91] + x413*M[113] + x416*M[92] + x419*M[114] + x425*M[96] + x426*M[81] + x428*M[158] + x431*M[132] + x434*M[106] + x436*M[95] + x440*M[97] + x443*M[67] + x444*M[72] + x445*M[73] + x446*M[100] + x447*M[101] + x448*M[107] + x449*M[143] + x452*M[124] + x454*M[157] + x460*M[127] + x463*M[128] + x467*M[108] + x47*M[13] + x471*M[110] + x474*M[69] + x475*M[74] + x477*M[75] + x479*M[80] + x480*M[102] + x481*M[103] + x483*M[109] + x484*M[115] + x485*M[116] + x486*M[145] + x487*M[160] + x489*M[131] + x49*M[15] + x490*M[142] + x492*M[149] + x493*M[136] + x495*M[150] + x496*M[137] + x497*M[82] + x5*M[1] + x500*M[117] + x501*M[118] + x503*M[162] + x504*M[133] + x506*M[138] + x507*M[139] + x510*M[159] + x511*M[144] + x512*M[151] + x513*M[152] + x515*M[146] + x516*M[161] + x519*M[153] + x52*M[16] + x521*M[154] + x525*M[163] + x58*M[23] - x61*M[30] + x66*M[9] + x7*M[2] + x74*M[19] + x79*M[20] + x80*M[21] - x84*M[29] + x91*M[35] + x92*M[36] + x95*M[38];
#pragma omp atomic
L[2] += x100*M[50] + x101*M[44] + x103*M[45] + x112*M[57] - x116*M[71] + x123*M[49] - x13*M[5] + x132*M[55] + x136*M[76] + x141*M[77] - x142*M[70] + x152*M[83] + x161*M[112] + x167*M[85] + x169*M[10] + x170*M[104] + x174*M[105] + x181*M[121] + x184*M[20] + x186*M[22] + x187*M[23] - x192*M[148] + x194*M[8] + x20*M[1] + x203*M[37] + x204*M[38] + x205*M[41] + x211*M[62] + x212*M[14] + x214*M[18] + x222*M[111] + x223*M[28] + x229*M[155] + x231*M[119] + x236*M[156] - x237*M[147] + x240*M[17] + x248*M[40] + x250*M[35] + x256*M[61] + x26*M[3] + x261*M[56] + x266*M[58] + x269*M[65] + x270*M[59] + x271*M[66] + x282*M[93] + x289*M[86] + x29*M[6] + x290*M[87] + x293*M[94] + x294*M[99] + x296*M[90] + x300*M[126] + x301*M[135] + x309*M[51] + x31*M[7] + x312*M[27] + x314*M[31] + x315*M[24] + x318*M[32] + x320*M[39] + x321*M[46] + x323*M[47] + x326*M[43] + x327*M[52] + x329*M[64] - x33*M[4] + x330*M[73] + x334*M[84] + x337*M[98] + x339*M[42] + x342*M[120] + x348*M[134] + x349*M[140] + x352*M[122] + x353*M[141] + x354*M[123] + x359*M[113] + x365*M[53] + x375*M[89] + x379*M[125] + x381*M[129] + x382*M[130] + x383*M[33] + x388*M[54] + x390*M[72] + x396*M[60] + x399*M[78] + x40*M[9] + x405*M[79] + x406*M[48] + x410*M[88] + x413*M[106] + x419*M[107] + x422*M[114] + x425*M[92] + x426*M[75] + x428*M[150] + x431*M[128] + x434*M[100] + x436*M[91] + x442*M[115] + x443*M[63] + x444*M[67] + x445*M[68] + x446*M[95] + x447*M[96] + x448*M[101] + x449*M[137] + x45*M[16] + x454*M[149] + x457*M[157] + x460*M[124] + x466*M[158] + x467*M[102] + x47*M[11] + x470*M[117] + x475*M[69] + x476*M[80] + x478*M[81] + x479*M[74] + x480*M[97] + x482*M[116] + x483*M[103] + x484*M[108] + x485*M[109] + x486*M[139] + x487*M[152] + x489*M[127] + x49*M[12] + x490*M[136] + x492*M[142] + x493*M[131] + x495*M[143] + x496*M[132] + x498*M[82] + x5*M[0] + x500*M[110] + x502*M[118] + x503*M[154] + x505*M[159] + x506*M[133] + x508*M[160] + x510*M[151] + x511*M[138] + x512*M[144] + x513*M[145] + x516*M[153] + x518*M[161] + x519*M[146] + x52*M[13] + x522*M[162] + x523*M[163] + x58*M[21] - x61*M[26] + x69*M[15] + x77*M[29] + x79*M[19] + x83*M[30] - x84*M[25] + x9*M[2] + x91*M[34] + x95*M[36];
#pragma omp atomic
L[3] += x100*M[49] + x103*M[44] + x112*M[56] - x116*M[70] - x13*M[4] + x137*M[55] + x141*M[76] + x153*M[83] + x161*M[111] + x167*M[84] + x174*M[104] + x181*M[120] + x187*M[22] - x192*M[147] + x193*M[5] + x194*M[7] + x195*M[8] + x204*M[37] + x205*M[40] + x211*M[61] + x212*M[13] + x213*M[14] + x214*M[17] + x223*M[27] + x232*M[119] + x236*M[155] + x240*M[16] + x241*M[11] + x243*M[18] + x270*M[58] + x271*M[65] + x290*M[86] + x293*M[93] + x294*M[98] + x296*M[89] + x30*M[3] + x300*M[125] + x301*M[134] + x306*M[36] + x309*M[50] + x31*M[6] + x311*M[21] + x312*M[26] + x314*M[30] + x315*M[23] + x317*M[24] + x318*M[31] + x320*M[38] + x321*M[45] + x323*M[46] + x324*M[39] + x326*M[42] + x327*M[51] + x329*M[63] + x330*M[72] + x339*M[41] + x353*M[140] + x354*M[122] + x359*M[112] + x362*M[85] + x365*M[52] + x366*M[43] + x382*M[129] + x383*M[32] + x384*M[33] + x385*M[28] + x387*M[54] + x388*M[53] + x390*M[71] + x393*M[57] + x396*M[59] + x399*M[77] + x402*M[60] + x405*M[78] + x406*M[47] + x407*M[48] + x41*M[9] + x410*M[87] + x413*M[105] + x416*M[88] + x419*M[106] + x422*M[113] + x425*M[91] + x426*M[74] + x428*M[149] + x431*M[127] + x434*M[99] + x436*M[90] + x440*M[92] + x442*M[114] + x443*M[62] + x444*M[66] + x445*M[67] + x446*M[94] + x447*M[95] + x448*M[100] + x449*M[136] + x45*M[15] + x452*M[121] + x454*M[148] + x457*M[156] + x460*M[123] + x463*M[124] + x466*M[157] + x467*M[101] + x47*M[10] + x470*M[116] + x471*M[103] + x473*M[118] + x474*M[64] + x475*M[68] + x476*M[79] + x477*M[69] + x478*M[80] + x479*M[73] + x480*M[96] + x481*M[97] + x482*M[115] + x483*M[102] + x484*M[107] + x485*M[108] + x486*M[138] + x487*M[151] + x489*M[126] + x490*M[135] + x492*M[141] + x493*M[130] + x495*M[142] + x496*M[131] + x497*M[75] + x498*M[81] + x499*M[82] + x500*M[109] + x501*M[110] + x502*M[117] + x503*M[153] + x504*M[128] + x505*M[158] + x506*M[132] + x507*M[133] + x508*M[159] + x510*M[150] + x511*M[137] + x512*M[143] + x513*M[144] + x515*M[139] + x516*M[152] + x518*M[160] + x519*M[145] + x52*M[12] + x521*M[146] + x522*M[161] + x523*M[162] + x524*M[163] + x525*M[154] + x58*M[20] - x61*M[25] + x7*M[0] + x70*M[2] + x80*M[19] + x83*M[29] + x9*M[1] + x92*M[34] + x95*M[35];
#pragma omp atomic
L[4] += x112*M[38] + x120*M[19] + x129*M[34] + x132*M[35] + x137*M[36] + x152*M[56] + x153*M[57] + x167*M[59] + x169*M[6] + x181*M[87] + x184*M[12] + x186*M[15] + x187*M[16] + x203*M[25] + x204*M[26] + x205*M[30] + x211*M[45] + x219*M[55] + x226*M[83] + x231*M[84] + x232*M[85] + x24*M[0] + x241*M[8] + x248*M[29] + x250*M[22] + x256*M[44] + x26*M[1] + x261*M[37] + x266*M[40] + x269*M[49] + x270*M[41] + x271*M[50] + x282*M[70] + x289*M[61] + x290*M[62] + x293*M[71] + x294*M[77] + x296*M[66] + x30*M[2] + x300*M[94] + x301*M[105] + x306*M[24] + x311*M[14] + x315*M[17] + x317*M[18] + x320*M[27] + x324*M[28] + x326*M[32] + x329*M[47] + x334*M[58] + x337*M[76] + x339*M[31] + x342*M[86] + x348*M[104] + x349*M[111] + x352*M[89] + x353*M[112] + x354*M[90] + x362*M[60] + x366*M[33] + x375*M[65] + x379*M[93] + x381*M[98] + x382*M[99] + x393*M[39] + x396*M[42] + x40*M[4] + x402*M[43] + x41*M[5] + x410*M[63] + x416*M[64] + x425*M[68] + x431*M[96] + x434*M[78] + x436*M[67] + x440*M[69] + x443*M[46] + x444*M[51] + x445*M[52] + x446*M[72] + x447*M[73] + x448*M[79] + x449*M[107] + x452*M[88] + x460*M[91] + x463*M[92] + x467*M[80] + x47*M[7] + x471*M[82] + x474*M[48] + x475*M[53] + x477*M[54] + x480*M[74] + x481*M[75] + x483*M[81] + x486*M[109] + x489*M[95] + x490*M[106] + x492*M[113] + x493*M[100] + x495*M[114] + x496*M[101] + x504*M[97] + x506*M[102] + x507*M[103] + x511*M[108] + x512*M[115] + x513*M[116] + x515*M[110] + x519*M[117] + x521*M[118] + x58*M[13] + x66*M[3] + x74*M[9] + x79*M[10] + x80*M[11] + x91*M[20] + x92*M[21] + x95*M[23];
#pragma omp atomic
L[5] += x101*M[29] + x103*M[30] + x112*M[36] - x116*M[50] - x13*M[2] + x132*M[34] - x142*M[49] + x152*M[55] + x167*M[57] + x169*M[4] + x170*M[76] + x174*M[77] + x181*M[85] + x184*M[10] + x186*M[12] + x187*M[13] - x192*M[112] + x203*M[22] + x204*M[23] + x205*M[26] + x211*M[41] + x212*M[8] + x223*M[18] + x231*M[83] - x237*M[111] + x248*M[25] + x250*M[20] + x256*M[40] + x26*M[0] + x261*M[35] + x266*M[37] + x269*M[44] + x270*M[38] + x271*M[45] + x282*M[65] + x289*M[58] + x290*M[59] + x293*M[66] + x294*M[71] + x296*M[62] + x300*M[90] + x301*M[99] + x312*M[17] + x315*M[14] + x320*M[24] + x321*M[31] + x323*M[32] + x326*M[28] + x329*M[43] - x33*M[1] + x330*M[52] + x334*M[56] + x337*M[70] + x339*M[27] + x342*M[84] + x348*M[98] + x349*M[104] + x352*M[86] + x353*M[105] + x354*M[87] + x375*M[61] + x379*M[89] + x381*M[93] + x382*M[94] + x390*M[51] + x396*M[39] + x40*M[3] + x406*M[33] + x410*M[60] + x413*M[78] + x419*M[79] + x425*M[64] + x426*M[54] + x428*M[114] + x431*M[92] + x434*M[72] + x436*M[63] + x443*M[42] + x444*M[46] + x445*M[47] + x446*M[67] + x447*M[68] + x448*M[73] + x449*M[101] + x454*M[113] + x460*M[88] + x467*M[74] + x47*M[5] + x475*M[48] + x479*M[53] + x480*M[69] + x483*M[75] + x484*M[80] + x485*M[81] + x486*M[103] + x487*M[116] + x489*M[91] + x49*M[6] + x490*M[100] + x492*M[106] + x493*M[95] + x495*M[107] + x496*M[96] + x500*M[82] + x503*M[118] + x506*M[97] + x510*M[115] + x511*M[102] + x512*M[108] + x513*M[109] + x516*M[117] + x519*M[110] + x52*M[7] + x58*M[11] - x61*M[16] + x79*M[9] - x84*M[15] + x91*M[19] + x95*M[21];
#pragma omp atomic
L[6] += x103*M[29] + x112*M[35] - x116*M[49] - x13*M[1] + x137*M[34] + x153*M[55] + x167*M[56] + x174*M[76] + x181*M[84] + x187*M[12] - x192*M[111] + x193*M[2] + x204*M[22] + x205*M[25] + x211*M[40] + x212*M[7] + x213*M[8] + x223*M[17] + x232*M[83] + x241*M[5] + x270*M[37] + x271*M[44] + x290*M[58] + x293*M[65] + x294*M[70] + x296*M[61] + x30*M[0] + x300*M[89] + x301*M[98] + x306*M[21] + x311*M[11] + x312*M[16] + x315*M[13] + x317*M[14] + x320*M[23] + x321*M[30] + x323*M[31] + x324*M[24] + x326*M[27] + x329*M[42] + x330*M[51] + x339*M[26] + x353*M[104] + x354*M[86] + x362*M[57] + x366*M[28] + x382*M[93] + x385*M[18] + x390*M[50] + x393*M[36] + x396*M[38] + x402*M[39] + x406*M[32] + x407*M[33] + x41*M[3] + x410*M[59] + x413*M[77] + x416*M[60] + x419*M[78] + x425*M[63] + x426*M[53] + x428*M[113] + x431*M[91] + x434*M[71] + x436*M[62] + x440*M[64] + x443*M[41] + x444*M[45] + x445*M[46] + x446*M[66] + x447*M[67] + x448*M[72] + x449*M[100] + x452*M[85] + x454*M[112] + x460*M[87] + x463*M[88] + x467*M[73] + x47*M[4] + x471*M[75] + x474*M[43] + x475*M[47] + x477*M[48] + x479*M[52] + x480*M[68] + x481*M[69] + x483*M[74] + x484*M[79] + x485*M[80] + x486*M[102] + x487*M[115] + x489*M[90] + x490*M[99] + x492*M[105] + x493*M[94] + x495*M[106] + x496*M[95] + x497*M[54] + x500*M[81] + x501*M[82] + x503*M[117] + x504*M[92] + x506*M[96] + x507*M[97] + x510*M[114] + x511*M[101] + x512*M[107] + x513*M[108] + x515*M[103] + x516*M[116] + x519*M[109] + x52*M[6] + x521*M[110] + x525*M[118] + x58*M[10] - x61*M[15] + x80*M[9] + x92*M[19] + x95*M[20];
#pragma omp atomic
L[7] += x100*M[30] + x101*M[25] + x103*M[26] - x116*M[45] + x123*M[29] + x136*M[49] + x141*M[50] - x142*M[44] + x161*M[77] + x169*M[3] + x170*M[70] + x174*M[71] + x184*M[9] + x186*M[10] + x187*M[11] - x192*M[105] + x203*M[20] + x204*M[21] + x205*M[23] + x211*M[38] + x222*M[76] + x229*M[111] + x236*M[112] - x237*M[104] + x240*M[8] + x248*M[22] + x250*M[19] + x256*M[37] + x261*M[34] + x266*M[35] + x269*M[40] + x270*M[36] + x271*M[41] + x282*M[61] + x289*M[56] + x29*M[1] + x290*M[57] + x293*M[62] + x294*M[66] + x296*M[59] + x300*M[87] + x301*M[94] + x309*M[31] + x31*M[2] + x312*M[14] + x314*M[17] + x318*M[18] + x321*M[27] + x323*M[28] + x327*M[32] - x33*M[0] + x330*M[47] + x334*M[55] + x337*M[65] + x339*M[24] + x342*M[83] + x348*M[93] + x349*M[98] + x352*M[84] + x353*M[99] + x354*M[85] + x359*M[78] + x365*M[33] + x375*M[58] + x379*M[86] + x381*M[89] + x382*M[90] + x390*M[46] + x399*M[51] + x405*M[52] + x413*M[72] + x419*M[73] + x422*M[79] + x428*M[107] + x434*M[67] + x436*M[60] + x442*M[80] + x443*M[39] + x444*M[42] + x445*M[43] + x446*M[63] + x447*M[64] + x448*M[68] + x449*M[96] + x45*M[7] + x454*M[106] + x457*M[113] + x466*M[114] + x467*M[69] + x470*M[82] + x476*M[53] + x478*M[54] + x479*M[48] + x482*M[81] + x484*M[74] + x485*M[75] + x487*M[109] + x489*M[88] + x49*M[4] + x490*M[95] + x492*M[100] + x493*M[91] + x495*M[101] + x496*M[92] + x505*M[115] + x508*M[116] + x510*M[108] + x511*M[97] + x512*M[102] + x513*M[103] + x516*M[110] + x518*M[117] + x52*M[5] + x522*M[118] - x61*M[13] + x69*M[6] + x77*M[15] + x83*M[16] - x84*M[12];
#pragma omp atomic
L[8] += x100*M[29] + x103*M[25] + x112*M[34] - x116*M[44] - x13*M[0] + x141*M[49] + x161*M[76] + x167*M[55] + x174*M[70] + x181*M[83] + x187*M[10] - x192*M[104] + x194*M[2] + x204*M[20] + x205*M[22] + x211*M[37] + x212*M[5] + x214*M[8] + x223*M[14] + x236*M[111] + x240*M[7] + x270*M[35] + x271*M[40] + x290*M[56] + x293*M[61] + x294*M[65] + x296*M[58] + x300*M[86] + x301*M[93] + x309*M[30] + x31*M[1] + x312*M[13] + x314*M[16] + x315*M[11] + x318*M[17] + x320*M[21] + x321*M[26] + x323*M[27] + x326*M[24] + x327*M[31] + x329*M[39] + x330*M[46] + x339*M[23] + x353*M[98] + x354*M[84] + x359*M[77] + x365*M[32] + x382*M[89] + x383*M[18] + x388*M[33] + x390*M[45] + x396*M[36] + x399*M[50] + x405*M[51] + x406*M[28] + x410*M[57] + x413*M[71] + x419*M[72] + x422*M[78] + x425*M[60] + x426*M[48] + x428*M[106] + x431*M[88] + x434*M[66] + x436*M[59] + x442*M[79] + x443*M[38] + x444*M[41] + x445*M[42] + x446*M[62] + x447*M[63] + x448*M[67] + x449*M[95] + x45*M[6] + x454*M[105] + x457*M[112] + x460*M[85] + x466*M[113] + x467*M[68] + x47*M[3] + x470*M[81] + x475*M[43] + x476*M[52] + x478*M[53] + x479*M[47] + x480*M[64] + x482*M[80] + x483*M[69] + x484*M[73] + x485*M[74] + x486*M[97] + x487*M[108] + x489*M[87] + x490*M[94] + x492*M[99] + x493*M[90] + x495*M[100] + x496*M[91] + x498*M[54] + x500*M[75] + x502*M[82] + x503*M[110] + x505*M[114] + x506*M[92] + x508*M[115] + x510*M[107] + x511*M[96] + x512*M[101] + x513*M[102] + x516*M[109] + x518*M[116] + x519*M[103] + x52*M[4] + x522*M[117] + x523*M[118] + x58*M[9] - x61*M[12] + x83*M[15] + x95*M[19];
#pragma omp atomic
L[9] += x193*M[0] + x194*M[1] + x195*M[2] + x212*M[4] + x213*M[5] + x214*M[7] + x223*M[13] + x240*M[6] + x241*M[3] + x243*M[8] + x306*M[19] + x309*M[29] + x311*M[9] + x312*M[12] + x314*M[15] + x315*M[10] + x317*M[11] + x318*M[16] + x320*M[20] + x321*M[25] + x323*M[26] + x324*M[21] + x326*M[23] + x327*M[30] + x329*M[38] + x330*M[45] + x339*M[22] + x359*M[76] + x362*M[55] + x365*M[31] + x366*M[24] + x383*M[17] + x384*M[18] + x385*M[14] + x387*M[33] + x388*M[32] + x390*M[44] + x393*M[34] + x396*M[35] + x399*M[49] + x402*M[36] + x405*M[50] + x406*M[27] + x407*M[28] + x410*M[56] + x413*M[70] + x416*M[57] + x419*M[71] + x422*M[77] + x425*M[59] + x426*M[47] + x428*M[105] + x431*M[87] + x434*M[65] + x436*M[58] + x440*M[60] + x442*M[78] + x443*M[37] + x444*M[40] + x445*M[41] + x446*M[61] + x447*M[62] + x448*M[66] + x449*M[94] + x452*M[83] + x454*M[104] + x457*M[111] + x460*M[84] + x463*M[85] + x466*M[112] + x467*M[67] + x470*M[80] + x471*M[69] + x473*M[82] + x474*M[39] + x475*M[42] + x476*M[51] + x477*M[43] + x478*M[52] + x479*M[46] + x480*M[63] + x481*M[64] + x482*M[79] + x483*M[68] + x484*M[72] + x485*M[73] + x486*M[96] + x487*M[107] + x489*M[86] + x490*M[93] + x492*M[98] + x493*M[89] + x495*M[99] + x496*M[90] + x497*M[48] + x498*M[53] + x499*M[54] + x500*M[74] + x501*M[75] + x502*M[81] + x503*M[109] + x504*M[88] + x505*M[113] + x506*M[91] + x507*M[92] + x508*M[114] + x510*M[106] + x511*M[95] + x512*M[100] + x513*M[101] + x515*M[97] + x516*M[108] + x518*M[115] + x519*M[102] + x521*M[103] + x522*M[116] + x523*M[117] + x524*M[118] + x525*M[110];
#pragma omp atomic
L[10] += x112*M[23] + x120*M[9] + x129*M[19] + x132*M[20] + x137*M[21] + x152*M[35] + x153*M[36] + x167*M[38] + x181*M[59] + x184*M[6] + x203*M[15] + x204*M[16] + x211*M[30] + x219*M[34] + x226*M[55] + x231*M[56] + x232*M[57] + x250*M[12] + x256*M[29] + x261*M[22] + x266*M[25] + x270*M[26] + x282*M[49] + x289*M[40] + x290*M[41] + x293*M[50] + x296*M[45] + x300*M[66] + x301*M[77] + x306*M[14] + x311*M[8] + x320*M[17] + x324*M[18] + x329*M[32] + x334*M[37] + x342*M[58] + x348*M[76] + x352*M[61] + x354*M[62] + x362*M[39] + x375*M[44] + x379*M[65] + x381*M[70] + x382*M[71] + x393*M[24] + x396*M[27] + x40*M[1] + x402*M[28] + x41*M[2] + x410*M[42] + x416*M[43] + x425*M[47] + x431*M[68] + x436*M[46] + x440*M[48] + x443*M[31] + x446*M[51] + x447*M[52] + x449*M[79] + x452*M[60] + x460*M[63] + x463*M[64] + x474*M[33] + x480*M[53] + x481*M[54] + x486*M[81] + x489*M[67] + x490*M[78] + x493*M[72] + x496*M[73] + x504*M[69] + x506*M[74] + x507*M[75] + x511*M[80] + x515*M[82] + x58*M[7] + x66*M[0] + x74*M[3] + x79*M[4] + x80*M[5] + x91*M[10] + x92*M[11] + x95*M[13];
#pragma omp atomic
L[11] += x112*M[21] + x132*M[19] + x152*M[34] + x167*M[36] + x169*M[1] + x181*M[57] + x184*M[4] + x186*M[6] + x187*M[7] + x203*M[12] + x204*M[13] + x205*M[16] + x211*M[26] + x231*M[55] + x248*M[15] + x250*M[10] + x256*M[25] + x261*M[20] + x266*M[22] + x269*M[29] + x270*M[23] + x271*M[30] + x282*M[44] + x289*M[37] + x290*M[38] + x293*M[45] + x294*M[50] + x296*M[41] + x300*M[62] + x301*M[71] + x315*M[8] + x320*M[14] + x326*M[18] + x329*M[28] + x334*M[35] + x337*M[49] + x339*M[17] + x342*M[56] + x348*M[70] + x349*M[76] + x352*M[58] + x353*M[77] + x354*M[59] + x375*M[40] + x379*M[61] + x381*M[65] + x382*M[66] + x396*M[24] + x40*M[0] + x410*M[39] + x425*M[43] + x431*M[64] + x434*M[51] + x436*M[42] + x443*M[27] + x444*M[31] + x445*M[32] + x446*M[46] + x447*M[47] + x448*M[52] + x449*M[73] + x460*M[60] + x467*M[53] + x47*M[2] + x475*M[33] + x480*M[48] + x483*M[54] + x486*M[75] + x489*M[63] + x490*M[72] + x492*M[78] + x493*M[67] + x495*M[79] + x496*M[68] + x506*M[69] + x511*M[74] + x512*M[80] + x513*M[81] + x519*M[82] + x58*M[5] + x79*M[3] + x91*M[9] + x95*M[11];
#pragma omp atomic
L[12] += x112*M[20] + x137*M[19] + x153*M[34] + x167*M[35] + x181*M[56] + x187*M[6] + x204*M[12] + x205*M[15] + x211*M[25] + x232*M[55] + x241*M[2] + x270*M[22] + x271*M[29] + x290*M[37] + x293*M[44] + x294*M[49] + x296*M[40] + x300*M[61] + x301*M[70] + x306*M[11] + x311*M[5] + x315*M[7] + x317*M[8] + x320*M[13] + x324*M[14] + x326*M[17] + x329*M[27] + x339*M[16] + x353*M[76] + x354*M[58] + x362*M[36] + x366*M[18] + x382*M[65] + x393*M[21] + x396*M[23] + x402*M[24] + x41*M[0] + x410*M[38] + x416*M[39] + x425*M[42] + x431*M[63] + x434*M[50] + x436*M[41] + x440*M[43] + x443*M[26] + x444*M[30] + x445*M[31] + x446*M[45] + x447*M[46] + x448*M[51] + x449*M[72] + x452*M[57] + x460*M[59] + x463*M[60] + x467*M[52] + x47*M[1] + x471*M[54] + x474*M[28] + x475*M[32] + x477*M[33] + x480*M[47] + x481*M[48] + x483*M[53] + x486*M[74] + x489*M[62] + x490*M[71] + x492*M[77] + x493*M[66] + x495*M[78] + x496*M[67] + x504*M[64] + x506*M[68] + x507*M[69] + x511*M[73] + x512*M[79] + x513*M[80] + x515*M[75] + x519*M[81] + x521*M[82] + x58*M[4] + x80*M[3] + x92*M[9] + x95*M[10];
#pragma omp atomic
L[13] += x101*M[15] + x103*M[16] - x116*M[30] - x142*M[29] + x169*M[0] + x170*M[49] + x174*M[50] + x184*M[3] + x186*M[4] + x187*M[5] - x192*M[77] + x203*M[10] + x204*M[11] + x205*M[13] + x211*M[23] - x237*M[76] + x248*M[12] + x250*M[9] + x256*M[22] + x261*M[19] + x266*M[20] + x269*M[25] + x270*M[21] + x271*M[26] + x282*M[40] + x289*M[35] + x290*M[36] + x293*M[41] + x294*M[45] + x296*M[38] + x300*M[59] + x301*M[66] + x312*M[8] + x321*M[17] + x323*M[18] + x330*M[32] + x334*M[34] + x337*M[44] + x339*M[14] + x342*M[55] + x348*M[65] + x349*M[70] + x352*M[56] + x353*M[71] + x354*M[57] + x375*M[37] + x379*M[58] + x381*M[61] + x382*M[62] + x390*M[31] + x413*M[51] + x419*M[52] + x428*M[79] + x434*M[46] + x436*M[39] + x443*M[24] + x444*M[27] + x445*M[28] + x446*M[42] + x447*M[43] + x448*M[47] + x449*M[68] + x454*M[78] + x467*M[48] + x479*M[33] + x484*M[53] + x485*M[54] + x487*M[81] + x489*M[60] + x49*M[1] + x490*M[67] + x492*M[72] + x493*M[63] + x495*M[73] + x496*M[64] + x510*M[80] + x511*M[69] + x512*M[74] + x513*M[75] + x516*M[82] + x52*M[2] - x61*M[7] - x84*M[6];
#pragma omp atomic
L[14] += x103*M[15] + x112*M[19] - x116*M[29] + x167*M[34] + x174*M[49] + x181*M[55] + x187*M[4] - x192*M[76] + x204*M[10] + x205*M[12] + x211*M[22] + x212*M[2] + x223*M[8] + x270*M[20] + x271*M[25] + x290*M[35] + x293*M[40] + x294*M[44] + x296*M[37] + x300*M[58] + x301*M[65] + x312*M[7] + x315*M[5] + x320*M[11] + x321*M[16] + x323*M[17] + x326*M[14] + x329*M[24] + x330*M[31] + x339*M[13] + x353*M[70] + x354*M[56] + x382*M[61] + x390*M[30] + x396*M[21] + x406*M[18] + x410*M[36] + x413*M[50] + x419*M[51] + x425*M[39] + x426*M[33] + x428*M[78] + x431*M[60] + x434*M[45] + x436*M[38] + x443*M[23] + x444*M[26] + x445*M[27] + x446*M[41] + x447*M[42] + x448*M[46] + x449*M[67] + x454*M[77] + x460*M[57] + x467*M[47] + x47*M[0] + x475*M[28] + x479*M[32] + x480*M[43] + x483*M[48] + x484*M[52] + x485*M[53] + x486*M[69] + x487*M[80] + x489*M[59] + x490*M[66] + x492*M[71] + x493*M[62] + x495*M[72] + x496*M[63] + x500*M[54] + x503*M[82] + x506*M[64] + x510*M[79] + x511*M[68] + x512*M[73] + x513*M[74] + x516*M[81] + x519*M[75] + x52*M[1] + x58*M[3] - x61*M[6] + x95*M[9];
#pragma omp atomic
L[15] += x212*M[1] + x213*M[2] + x223*M[7] + x241*M[0] + x306*M[9] + x311*M[3] + x312*M[6] + x315*M[4] + x317*M[5] + x320*M[10] + x321*M[15] + x323*M[16] + x324*M[11] + x326*M[13] + x329*M[23] + x330*M[30] + x339*M[12] + x362*M[34] + x366*M[14] + x385*M[8] + x390*M[29] + x393*M[19] + x396*M[20] + x402*M[21] + x406*M[17] + x407*M[18] + x410*M[35] + x413*M[49] + x416*M[36] + x419*M[50] + x425*M[38] + x426*M[32] + x428*M[77] + x431*M[59] + x434*M[44] + x436*M[37] + x440*M[39] + x443*M[22] + x444*M[25] + x445*M[26] + x446*M[40] + x447*M[41] + x448*M[45] + x449*M[66] + x452*M[55] + x454*M[76] + x460*M[56] + x463*M[57] + x467*M[46] + x471*M[48] + x474*M[24] + x475*M[27] + x477*M[28] + x479*M[31] + x480*M[42] + x481*M[43] + x483*M[47] + x484*M[51] + x485*M[52] + x486*M[68] + x487*M[79] + x489*M[58] + x490*M[65] + x492*M[70] + x493*M[61] + x495*M[71] + x496*M[62] + x497*M[33] + x500*M[53] + x501*M[54] + x503*M[81] + x504*M[60] + x506*M[63] + x507*M[64] + x510*M[78] + x511*M[67] + x512*M[72] + x513*M[73] + x515*M[69] + x516*M[80] + x519*M[74] + x521*M[75] + x525*M[82];
#pragma omp atomic
L[16] += x100*M[16] + x101*M[12] + x103*M[13] - x116*M[26] + x123*M[15] + x136*M[29] + x141*M[30] - x142*M[25] + x161*M[50] + x170*M[44] + x174*M[45] + x186*M[3] - x192*M[71] + x203*M[9] + x205*M[11] + x211*M[21] + x222*M[49] + x229*M[76] + x236*M[77] - x237*M[70] + x248*M[10] + x256*M[20] + x266*M[19] + x269*M[22] + x271*M[23] + x282*M[37] + x289*M[34] + x293*M[38] + x294*M[41] + x296*M[36] + x300*M[57] + x301*M[62] + x309*M[17] + x314*M[8] + x321*M[14] + x327*M[18] + x330*M[28] + x337*M[40] + x348*M[61] + x349*M[65] + x352*M[55] + x353*M[66] + x359*M[51] + x375*M[35] + x379*M[56] + x381*M[58] + x382*M[59] + x390*M[27] + x399*M[31] + x405*M[32] + x413*M[46] + x419*M[47] + x422*M[52] + x428*M[73] + x434*M[42] + x442*M[53] + x444*M[24] + x446*M[39] + x448*M[43] + x449*M[64] + x45*M[2] + x454*M[72] + x457*M[78] + x466*M[79] + x476*M[33] + x482*M[54] + x484*M[48] + x487*M[75] + x49*M[0] + x490*M[63] + x492*M[67] + x493*M[60] + x495*M[68] + x505*M[80] + x508*M[81] + x510*M[74] + x512*M[69] + x518*M[82] - x61*M[5] + x69*M[1] + x77*M[6] + x83*M[7] - x84*M[4];
#pragma omp atomic
L[17] += x100*M[15] + x103*M[12] - x116*M[25] + x141*M[29] + x161*M[49] + x174*M[44] + x187*M[3] - x192*M[70] + x204*M[9] + x205*M[10] + x211*M[20] + x236*M[76] + x240*M[2] + x270*M[19] + x271*M[22] + x290*M[34] + x293*M[37] + x294*M[40] + x296*M[35] + x300*M[56] + x301*M[61] + x309*M[16] + x312*M[5] + x314*M[7] + x318*M[8] + x321*M[13] + x323*M[14] + x327*M[17] + x330*M[27] + x339*M[11] + x353*M[65] + x354*M[55] + x359*M[50] + x365*M[18] + x382*M[58] + x390*M[26] + x399*M[30] + x405*M[31] + x413*M[45] + x419*M[46] + x422*M[51] + x428*M[72] + x434*M[41] + x436*M[36] + x442*M[52] + x443*M[21] + x444*M[23] + x445*M[24] + x446*M[38] + x447*M[39] + x448*M[42] + x449*M[63] + x45*M[1] + x454*M[71] + x457*M[77] + x466*M[78] + x467*M[43] + x470*M[54] + x476*M[32] + x478*M[33] + x479*M[28] + x482*M[53] + x484*M[47] + x485*M[48] + x487*M[74] + x489*M[57] + x490*M[62] + x492*M[66] + x493*M[59] + x495*M[67] + x496*M[60] + x505*M[79] + x508*M[80] + x510*M[73] + x511*M[64] + x512*M[68] + x513*M[69] + x516*M[75] + x518*M[81] + x52*M[0] + x522*M[82] - x61*M[4] + x83*M[6];
#pragma omp atomic
L[18] += x212*M[0] + x214*M[2] + x223*M[5] + x240*M[1] + x309*M[15] + x312*M[4] + x314*M[6] + x315*M[3] + x318*M[7] + x320*M[9] + x321*M[12] + x323*M[13] + x326*M[11] + x327*M[16] + x329*M[21] + x330*M[26] + x339*M[10] + x359*M[49] + x365*M[17] + x383*M[8] + x388*M[18] + x390*M[25] + x396*M[19] + x399*M[29] + x405*M[30] + x406*M[14] + x410*M[34] + x413*M[44] + x419*M[45] + x422*M[50] + x425*M[36] + x426*M[28] + x428*M[71] + x431*M[57] + x434*M[40] + x436*M[35] + x442*M[51] + x443*M[20] + x444*M[22] + x445*M[23] + x446*M[37] + x447*M[38] + x448*M[41] + x449*M[62] + x454*M[70] + x457*M[76] + x460*M[55] + x466*M[77] + x467*M[42] + x470*M[53] + x475*M[24] + x476*M[31] + x478*M[32] + x479*M[27] + x480*M[39] + x482*M[52] + x483*M[43] + x484*M[46] + x485*M[47] + x486*M[64] + x487*M[73] + x489*M[56] + x490*M[61] + x492*M[65] + x493*M[58] + x495*M[66] + x496*M[59] + x498*M[33] + x500*M[48] + x502*M[54] + x503*M[75] + x505*M[78] + x506*M[60] + x508*M[79] + x510*M[72] + x511*M[63] + x512*M[67] + x513*M[68] + x516*M[74] + x518*M[80] + x519*M[69] + x522*M[81] + x523*M[82];
#pragma omp atomic
L[19] += x213*M[0] + x214*M[1] + x223*M[4] + x243*M[2] + x317*M[3] + x318*M[6] + x323*M[12] + x324*M[9] + x326*M[10] + x327*M[15] + x329*M[20] + x330*M[25] + x365*M[16] + x366*M[11] + x383*M[7] + x384*M[8] + x385*M[5] + x387*M[18] + x388*M[17] + x402*M[19] + x405*M[29] + x406*M[13] + x407*M[14] + x416*M[34] + x419*M[44] + x422*M[49] + x425*M[35] + x426*M[27] + x428*M[70] + x431*M[56] + x440*M[36] + x442*M[50] + x445*M[22] + x447*M[37] + x448*M[40] + x449*M[61] + x463*M[55] + x466*M[76] + x467*M[41] + x470*M[52] + x471*M[43] + x473*M[54] + x474*M[21] + x475*M[23] + x476*M[30] + x477*M[24] + x478*M[31] + x479*M[26] + x480*M[38] + x481*M[39] + x482*M[51] + x483*M[42] + x484*M[45] + x485*M[46] + x486*M[63] + x487*M[72] + x495*M[65] + x496*M[58] + x497*M[28] + x498*M[32] + x499*M[33] + x500*M[47] + x501*M[48] + x502*M[53] + x503*M[74] + x504*M[57] + x505*M[77] + x506*M[59] + x507*M[60] + x508*M[78] + x510*M[71] + x511*M[62] + x512*M[66] + x513*M[67] + x515*M[64] + x516*M[73] + x518*M[79] + x519*M[68] + x521*M[69] + x522*M[80] + x523*M[81] + x524*M[82] + x525*M[75];
#pragma omp atomic
L[20] += x112*M[13] + x120*M[3] + x129*M[9] + x132*M[10] + x137*M[11] + x152*M[20] + x153*M[21] + x167*M[23] + x181*M[38] + x219*M[19] + x226*M[34] + x231*M[35] + x232*M[36] + x250*M[6] + x261*M[12] + x266*M[15] + x270*M[16] + x289*M[25] + x290*M[26] + x296*M[30] + x300*M[45] + x306*M[8] + x334*M[22] + x342*M[37] + x352*M[40] + x354*M[41] + x362*M[24] + x375*M[29] + x379*M[44] + x381*M[49] + x382*M[50] + x393*M[14] + x396*M[17] + x402*M[18] + x410*M[27] + x416*M[28] + x425*M[32] + x431*M[47] + x436*M[31] + x440*M[33] + x452*M[39] + x460*M[42] + x463*M[43] + x489*M[46] + x493*M[51] + x496*M[52] + x504*M[48] + x506*M[53] + x507*M[54] + x74*M[0] + x79*M[1] + x80*M[2] + x91*M[4] + x92*M[5] + x95*M[7];
#pragma omp atomic
L[21] += x112*M[11] + x132*M[9] + x152*M[19] + x167*M[21] + x181*M[36] + x184*M[1] + x203*M[6] + x204*M[7] + x211*M[16] + x231*M[34] + x250*M[4] + x256*M[15] + x261*M[10] + x266*M[12] + x270*M[13] + x282*M[29] + x289*M[22] + x290*M[23] + x293*M[30] + x296*M[26] + x300*M[41] + x301*M[50] + x320*M[8] + x329*M[18] + x334*M[20] + x342*M[35] + x348*M[49] + x352*M[37] + x354*M[38] + x375*M[25] + x379*M[40] + x381*M[44] + x382*M[45] + x396*M[14] + x410*M[24] + x425*M[28] + x431*M[43] + x436*M[27] + x443*M[17] + x446*M[31] + x447*M[32] + x449*M[52] + x460*M[39] + x480*M[33] + x486*M[54] + x489*M[42] + x490*M[51] + x493*M[46] + x496*M[47] + x506*M[48] + x511*M[53] + x58*M[2] + x79*M[0] + x91*M[3] + x95*M[5];
#pragma omp atomic
L[22] += x112*M[10] + x137*M[9] + x153*M[19] + x167*M[20] + x181*M[35] + x204*M[6] + x211*M[15] + x232*M[34] + x270*M[12] + x290*M[22] + x293*M[29] + x296*M[25] + x300*M[40] + x301*M[49] + x306*M[5] + x311*M[2] + x320*M[7] + x324*M[8] + x329*M[17] + x354*M[37] + x362*M[21] + x382*M[44] + x393*M[11] + x396*M[13] + x402*M[14] + x410*M[23] + x416*M[24] + x425*M[27] + x431*M[42] + x436*M[26] + x440*M[28] + x443*M[16] + x446*M[30] + x447*M[31] + x449*M[51] + x452*M[36] + x460*M[38] + x463*M[39] + x474*M[18] + x480*M[32] + x481*M[33] + x486*M[53] + x489*M[41] + x490*M[50] + x493*M[45] + x496*M[46] + x504*M[43] + x506*M[47] + x507*M[48] + x511*M[52] + x515*M[54] + x58*M[1] + x80*M[0] + x92*M[3] + x95*M[4];
#pragma omp atomic
L[23] += x184*M[0] + x186*M[1] + x187*M[2] + x203*M[4] + x204*M[5] + x205*M[7] + x211*M[13] + x248*M[6] + x250*M[3] + x256*M[12] + x261*M[9] + x266*M[10] + x269*M[15] + x270*M[11] + x271*M[16] + x282*M[25] + x289*M[20] + x290*M[21] + x293*M[26] + x294*M[30] + x296*M[23] + x300*M[38] + x301*M[45] + x334*M[19] + x337*M[29] + x339*M[8] + x342*M[34] + x348*M[44] + x349*M[49] + x352*M[35] + x353*M[50] + x354*M[36] + x375*M[22] + x379*M[37] + x381*M[40] + x382*M[41] + x434*M[31] + x436*M[24] + x443*M[14] + x444*M[17] + x445*M[18] + x446*M[27] + x447*M[28] + x448*M[32] + x449*M[47] + x467*M[33] + x489*M[39] + x490*M[46] + x492*M[51] + x493*M[42] + x495*M[52] + x496*M[43] + x511*M[48] + x512*M[53] + x513*M[54];
#pragma omp atomic
L[24] += x112*M[9] + x167*M[19] + x181*M[34] + x187*M[1] + x204*M[4] + x205*M[6] + x211*M[12] + x270*M[10] + x271*M[15] + x290*M[20] + x293*M[25] + x294*M[29] + x296*M[22] + x300*M[37] + x301*M[44] + x315*M[2] + x320*M[5] + x326*M[8] + x329*M[14] + x339*M[7] + x353*M[49] + x354*M[35] + x382*M[40] + x396*M[11] + x410*M[21] + x425*M[24] + x431*M[39] + x434*M[30] + x436*M[23] + x443*M[13] + x444*M[16] + x445*M[17] + x446*M[26] + x447*M[27] + x448*M[31] + x449*M[46] + x460*M[36] + x467*M[32] + x475*M[18] + x480*M[28] + x483*M[33] + x486*M[48] + x489*M[38] + x490*M[45] + x492*M[50] + x493*M[41] + x495*M[51] + x496*M[42] + x506*M[43] + x511*M[47] + x512*M[52] + x513*M[53] + x519*M[54] + x58*M[0] + x95*M[3];
#pragma omp atomic
L[25] += x306*M[3] + x311*M[0] + x315*M[1] + x317*M[2] + x320*M[4] + x324*M[5] + x326*M[7] + x329*M[13] + x339*M[6] + x362*M[19] + x366*M[8] + x393*M[9] + x396*M[10] + x402*M[11] + x410*M[20] + x416*M[21] + x425*M[23] + x431*M[38] + x434*M[29] + x436*M[22] + x440*M[24] + x443*M[12] + x444*M[15] + x445*M[16] + x446*M[25] + x447*M[26] + x448*M[30] + x449*M[45] + x452*M[34] + x460*M[35] + x463*M[36] + x467*M[31] + x471*M[33] + x474*M[14] + x475*M[17] + x477*M[18] + x480*M[27] + x481*M[28] + x483*M[32] + x486*M[47] + x489*M[37] + x490*M[44] + x492*M[49] + x493*M[40] + x495*M[50] + x496*M[41] + x504*M[39] + x506*M[42] + x507*M[43] + x511*M[46] + x512*M[51] + x513*M[52] + x515*M[48] + x519*M[53] + x521*M[54];
#pragma omp atomic
L[26] += x101*M[6] + x103*M[7] - x116*M[16] - x142*M[15] + x170*M[29] + x174*M[30] + x186*M[0] - x192*M[50] + x203*M[3] + x205*M[5] + x211*M[11] - x237*M[49] + x248*M[4] + x256*M[10] + x266*M[9] + x269*M[12] + x271*M[13] + x282*M[22] + x289*M[19] + x293*M[23] + x294*M[26] + x296*M[21] + x300*M[36] + x301*M[41] + x321*M[8] + x330*M[18] + x337*M[25] + x348*M[40] + x349*M[44] + x352*M[34] + x353*M[45] + x375*M[20] + x379*M[35] + x381*M[37] + x382*M[38] + x390*M[17] + x413*M[31] + x419*M[32] + x428*M[52] + x434*M[27] + x444*M[14] + x446*M[24] + x448*M[28] + x449*M[43] + x454*M[51] + x484*M[33] + x487*M[54] + x490*M[42] + x492*M[46] + x493*M[39] + x495*M[47] + x510*M[53] + x512*M[48] - x61*M[2] - x84*M[1];
#pragma omp atomic
L[27] += x103*M[6] - x116*M[15] + x174*M[29] + x187*M[0] - x192*M[49] + x204*M[3] + x205*M[4] + x211*M[10] + x270*M[9] + x271*M[12] + x290*M[19] + x293*M[22] + x294*M[25] + x296*M[20] + x300*M[35] + x301*M[40] + x312*M[2] + x321*M[7] + x323*M[8] + x330*M[17] + x339*M[5] + x353*M[44] + x354*M[34] + x382*M[37] + x390*M[16] + x413*M[30] + x419*M[31] + x428*M[51] + x434*M[26] + x436*M[21] + x443*M[11] + x444*M[13] + x445*M[14] + x446*M[23] + x447*M[24] + x448*M[27] + x449*M[42] + x454*M[50] + x467*M[28] + x479*M[18] + x484*M[32] + x485*M[33] + x487*M[53] + x489*M[36] + x490*M[41] + x492*M[45] + x493*M[38] + x495*M[46] + x496*M[39] + x510*M[52] + x511*M[43] + x512*M[47] + x513*M[48] + x516*M[54] - x61*M[1];
#pragma omp atomic
L[28] += x223*M[2] + x312*M[1] + x315*M[0] + x320*M[3] + x321*M[6] + x323*M[7] + x326*M[5] + x329*M[11] + x330*M[16] + x339*M[4] + x390*M[15] + x396*M[9] + x406*M[8] + x410*M[19] + x413*M[29] + x419*M[30] + x425*M[21] + x426*M[18] + x428*M[50] + x431*M[36] + x434*M[25] + x436*M[20] + x443*M[10] + x444*M[12] + x445*M[13] + x446*M[22] + x447*M[23] + x448*M[26] + x449*M[41] + x454*M[49] + x460*M[34] + x467*M[27] + x475*M[14] + x479*M[17] + x480*M[24] + x483*M[28] + x484*M[31] + x485*M[32] + x486*M[43] + x487*M[52] + x489*M[35] + x490*M[40] + x492*M[44] + x493*M[37] + x495*M[45] + x496*M[38] + x500*M[33] + x503*M[54] + x506*M[39] + x510*M[51] + x511*M[42] + x512*M[46] + x513*M[47] + x516*M[53] + x519*M[48];
#pragma omp atomic
L[29] += x223*M[1] + x317*M[0] + x323*M[6] + x324*M[3] + x326*M[4] + x329*M[10] + x330*M[15] + x366*M[5] + x385*M[2] + x402*M[9] + x406*M[7] + x407*M[8] + x416*M[19] + x419*M[29] + x425*M[20] + x426*M[17] + x428*M[49] + x431*M[35] + x440*M[21] + x445*M[12] + x447*M[22] + x448*M[25] + x449*M[40] + x463*M[34] + x467*M[26] + x471*M[28] + x474*M[11] + x475*M[13] + x477*M[14] + x479*M[16] + x480*M[23] + x481*M[24] + x483*M[27] + x484*M[30] + x485*M[31] + x486*M[42] + x487*M[51] + x495*M[44] + x496*M[37] + x497*M[18] + x500*M[32] + x501*M[33] + x503*M[53] + x504*M[36] + x506*M[38] + x507*M[39] + x510*M[50] + x511*M[41] + x512*M[45] + x513*M[46] + x515*M[43] + x516*M[52] + x519*M[47] + x521*M[48] + x525*M[54];
#pragma omp atomic
L[30] += x100*M[7] + x101*M[4] + x103*M[5] - x116*M[13] + x123*M[6] + x136*M[15] + x141*M[16] - x142*M[12] + x161*M[30] + x170*M[25] + x174*M[26] - x192*M[45] + x222*M[29] + x229*M[49] + x236*M[50] - x237*M[44] + x248*M[3] + x256*M[9] + x269*M[10] + x271*M[11] + x282*M[20] + x293*M[21] + x294*M[23] + x301*M[38] + x309*M[8] + x337*M[22] + x348*M[37] + x349*M[40] + x353*M[41] + x359*M[31] + x375*M[19] + x379*M[34] + x381*M[35] + x382*M[36] + x390*M[14] + x399*M[17] + x405*M[18] + x413*M[27] + x419*M[28] + x422*M[32] + x428*M[47] + x434*M[24] + x442*M[33] + x454*M[46] + x457*M[51] + x466*M[52] + x490*M[39] + x492*M[42] + x495*M[43] + x505*M[53] + x508*M[54] + x510*M[48] + x77*M[1] + x83*M[2] - x84*M[0];
#pragma omp atomic
L[31] += x100*M[6] + x103*M[4] - x116*M[12] + x141*M[15] + x161*M[29] + x174*M[25] - x192*M[44] + x205*M[3] + x211*M[9] + x236*M[49] + x271*M[10] + x293*M[20] + x294*M[22] + x296*M[19] + x300*M[34] + x301*M[37] + x309*M[7] + x314*M[2] + x321*M[5] + x327*M[8] + x330*M[14] + x353*M[40] + x359*M[30] + x382*M[35] + x390*M[13] + x399*M[16] + x405*M[17] + x413*M[26] + x419*M[27] + x422*M[31] + x428*M[46] + x434*M[23] + x442*M[32] + x444*M[11] + x446*M[21] + x448*M[24] + x449*M[39] + x454*M[45] + x457*M[50] + x466*M[51] + x476*M[18] + x482*M[33] + x484*M[28] + x487*M[48] + x490*M[38] + x492*M[41] + x493*M[36] + x495*M[42] + x505*M[52] + x508*M[53] + x510*M[47] + x512*M[43] + x518*M[54] - x61*M[0] + x83*M[1];
#pragma omp atomic
L[32] += x309*M[6] + x312*M[0] + x314*M[1] + x318*M[2] + x321*M[4] + x323*M[5] + x327*M[7] + x330*M[13] + x339*M[3] + x359*M[29] + x365*M[8] + x390*M[12] + x399*M[15] + x405*M[16] + x413*M[25] + x419*M[26] + x422*M[30] + x428*M[45] + x434*M[22] + x436*M[19] + x442*M[31] + x443*M[9] + x444*M[10] + x445*M[11] + x446*M[20] + x447*M[21] + x448*M[23] + x449*M[38] + x454*M[44] + x457*M[49] + x466*M[50] + x467*M[24] + x470*M[33] + x476*M[17] + x478*M[18] + x479*M[14] + x482*M[32] + x484*M[27] + x485*M[28] + x487*M[47] + x489*M[34] + x490*M[37] + x492*M[40] + x493*M[35] + x495*M[41] + x496*M[36] + x505*M[51] + x508*M[52] + x510*M[46] + x511*M[39] + x512*M[42] + x513*M[43] + x516*M[48] + x518*M[53] + x522*M[54];
#pragma omp atomic
L[33] += x223*M[0] + x318*M[1] + x323*M[4] + x326*M[3] + x327*M[6] + x329*M[9] + x330*M[12] + x365*M[7] + x383*M[2] + x388*M[8] + x405*M[15] + x406*M[5] + x419*M[25] + x422*M[29] + x425*M[19] + x426*M[14] + x428*M[44] + x431*M[34] + x442*M[30] + x445*M[10] + x447*M[20] + x448*M[22] + x449*M[37] + x466*M[49] + x467*M[23] + x470*M[32] + x475*M[11] + x476*M[16] + x478*M[17] + x479*M[13] + x480*M[21] + x482*M[31] + x483*M[24] + x484*M[26] + x485*M[27] + x486*M[39] + x487*M[46] + x495*M[40] + x496*M[35] + x498*M[18] + x500*M[28] + x502*M[33] + x503*M[48] + x505*M[50] + x506*M[36] + x508*M[51] + x510*M[45] + x511*M[38] + x512*M[41] + x513*M[42] + x516*M[47] + x518*M[52] + x519*M[43] + x522*M[53] + x523*M[54];
#pragma omp atomic
L[34] += x365*M[6] + x366*M[3] + x383*M[1] + x384*M[2] + x385*M[0] + x387*M[8] + x388*M[7] + x406*M[4] + x407*M[5] + x426*M[13] + x440*M[19] + x442*M[29] + x467*M[22] + x470*M[31] + x471*M[24] + x473*M[33] + x474*M[9] + x475*M[10] + x476*M[15] + x477*M[11] + x478*M[16] + x479*M[12] + x480*M[20] + x481*M[21] + x482*M[30] + x483*M[23] + x484*M[25] + x485*M[26] + x486*M[38] + x487*M[45] + x497*M[14] + x498*M[17] + x499*M[18] + x500*M[27] + x501*M[28] + x502*M[32] + x503*M[47] + x504*M[34] + x505*M[49] + x506*M[35] + x507*M[36] + x508*M[50] + x510*M[44] + x511*M[37] + x512*M[40] + x513*M[41] + x515*M[39] + x516*M[46] + x518*M[51] + x519*M[42] + x521*M[43] + x522*M[52] + x523*M[53] + x524*M[54] + x525*M[48];
#pragma omp atomic
L[35] += x112*M[7] + x120*M[0] + x129*M[3] + x132*M[4] + x137*M[5] + x152*M[10] + x153*M[11] + x167*M[13] + x181*M[23] + x219*M[9] + x226*M[19] + x231*M[20] + x232*M[21] + x261*M[6] + x289*M[15] + x290*M[16] + x300*M[30] + x334*M[12] + x342*M[22] + x352*M[25] + x354*M[26] + x362*M[14] + x379*M[29] + x393*M[8] + x410*M[17] + x416*M[18] + x431*M[32] + x452*M[24] + x460*M[27] + x463*M[28] + x489*M[31] + x504*M[33] + x91*M[1] + x92*M[2];
#pragma omp atomic
L[36] += x112*M[5] + x132*M[3] + x152*M[9] + x167*M[11] + x181*M[21] + x231*M[19] + x250*M[1] + x261*M[4] + x266*M[6] + x270*M[7] + x289*M[12] + x290*M[13] + x296*M[16] + x300*M[26] + x334*M[10] + x342*M[20] + x352*M[22] + x354*M[23] + x375*M[15] + x379*M[25] + x381*M[29] + x382*M[30] + x396*M[8] + x410*M[14] + x425*M[18] + x431*M[28] + x436*M[17] + x460*M[24] + x489*M[27] + x493*M[31] + x496*M[32] + x506*M[33] + x91*M[0] + x95*M[2];
#pragma omp atomic
L[37] += x112*M[4] + x137*M[3] + x153*M[9] + x167*M[10] + x181*M[20] + x232*M[19] + x270*M[6] + x290*M[12] + x296*M[15] + x300*M[25] + x306*M[2] + x354*M[22] + x362*M[11] + x382*M[29] + x393*M[5] + x396*M[7] + x402*M[8] + x410*M[13] + x416*M[14] + x425*M[17] + x431*M[27] + x436*M[16] + x440*M[18] + x452*M[21] + x460*M[23] + x463*M[24] + x489*M[26] + x493*M[30] + x496*M[31] + x504*M[28] + x506*M[32] + x507*M[33] + x92*M[0] + x95*M[1];
#pragma omp atomic
L[38] += x203*M[1] + x204*M[2] + x211*M[7] + x250*M[0] + x256*M[6] + x261*M[3] + x266*M[4] + x270*M[5] + x282*M[15] + x289*M[10] + x290*M[11] + x293*M[16] + x296*M[13] + x300*M[23] + x301*M[30] + x334*M[9] + x342*M[19] + x348*M[29] + x352*M[20] + x354*M[21] + x375*M[12] + x379*M[22] + x381*M[25] + x382*M[26] + x436*M[14] + x443*M[8] + x446*M[17] + x447*M[18] + x449*M[32] + x489*M[24] + x490*M[31] + x493*M[27] + x496*M[28] + x511*M[33];
#pragma omp atomic
L[39] += x112*M[3] + x167*M[9] + x181*M[19] + x204*M[1] + x211*M[6] + x270*M[4] + x290*M[10] + x293*M[15] + x296*M[12] + x300*M[22] + x301*M[29] + x320*M[2] + x329*M[8] + x354*M[20] + x382*M[25] + x396*M[5] + x410*M[11] + x425*M[14] + x431*M[24] + x436*M[13] + x443*M[7] + x446*M[16] + x447*M[17] + x449*M[31] + x460*M[21] + x480*M[18] + x486*M[33] + x489*M[23] + x490*M[30] + x493*M[26] + x496*M[27] + x506*M[28] + x511*M[32] + x95*M[0];
#pragma omp atomic
L[40] += x306*M[0] + x320*M[1] + x324*M[2] + x329*M[7] + x362*M[9] + x393*M[3] + x396*M[4] + x402*M[5] + x410*M[10] + x416*M[11] + x425*M[13] + x431*M[23] + x436*M[12] + x440*M[14] + x443*M[6] + x446*M[15] + x447*M[16] + x449*M[30] + x452*M[19] + x460*M[20] + x463*M[21] + x474*M[8] + x480*M[17] + x481*M[18] + x486*M[32] + x489*M[22] + x490*M[29] + x493*M[25] + x496*M[26] + x504*M[24] + x506*M[27] + x507*M[28] + x511*M[31] + x515*M[33];
#pragma omp atomic
L[41] += x203*M[0] + x205*M[2] + x211*M[5] + x248*M[1] + x256*M[4] + x266*M[3] + x269*M[6] + x271*M[7] + x282*M[12] + x289*M[9] + x293*M[13] + x294*M[16] + x296*M[11] + x300*M[21] + x301*M[26] + x337*M[15] + x348*M[25] + x349*M[29] + x352*M[19] + x353*M[30] + x375*M[10] + x379*M[20] + x381*M[22] + x382*M[23] + x434*M[17] + x444*M[8] + x446*M[14] + x448*M[18] + x449*M[28] + x490*M[27] + x492*M[31] + x493*M[24] + x495*M[32] + x512*M[33];
#pragma omp atomic
L[42] += x204*M[0] + x205*M[1] + x211*M[4] + x270*M[3] + x271*M[6] + x290*M[9] + x293*M[12] + x294*M[15] + x296*M[10] + x300*M[20] + x301*M[25] + x339*M[2] + x353*M[29] + x354*M[19] + x382*M[22] + x434*M[16] + x436*M[11] + x443*M[5] + x444*M[7] + x445*M[8] + x446*M[13] + x447*M[14] + x448*M[17] + x449*M[27] + x467*M[18] + x489*M[21] + x490*M[26] + x492*M[30] + x493*M[23] + x495*M[31] + x496*M[24] + x511*M[28] + x512*M[32] + x513*M[33];
#pragma omp atomic
L[43] += x320*M[0] + x326*M[2] + x329*M[5] + x339*M[1] + x396*M[3] + x410*M[9] + x425*M[11] + x431*M[21] + x434*M[15] + x436*M[10] + x443*M[4] + x444*M[6] + x445*M[7] + x446*M[12] + x447*M[13] + x448*M[16] + x449*M[26] + x460*M[19] + x467*M[17] + x475*M[8] + x480*M[14] + x483*M[18] + x486*M[28] + x489*M[20] + x490*M[25] + x492*M[29] + x493*M[22] + x495*M[30] + x496*M[23] + x506*M[24] + x511*M[27] + x512*M[31] + x513*M[32] + x519*M[33];
#pragma omp atomic
L[44] += x324*M[0] + x326*M[1] + x329*M[4] + x366*M[2] + x402*M[3] + x416*M[9] + x425*M[10] + x431*M[20] + x440*M[11] + x445*M[6] + x447*M[12] + x448*M[15] + x449*M[25] + x463*M[19] + x467*M[16] + x471*M[18] + x474*M[5] + x475*M[7] + x477*M[8] + x480*M[13] + x481*M[14] + x483*M[17] + x486*M[27] + x495*M[29] + x496*M[22] + x504*M[21] + x506*M[23] + x507*M[24] + x511*M[26] + x512*M[30] + x513*M[31] + x515*M[28] + x519*M[32] + x521*M[33];
#pragma omp atomic
L[45] += x101*M[1] + x103*M[2] - x116*M[7] - x142*M[6] + x170*M[15] + x174*M[16] - x192*M[30] - x237*M[29] + x248*M[0] + x256*M[3] + x269*M[4] + x271*M[5] + x282*M[10] + x293*M[11] + x294*M[13] + x301*M[23] + x337*M[12] + x348*M[22] + x349*M[25] + x353*M[26] + x375*M[9] + x379*M[19] + x381*M[20] + x382*M[21] + x390*M[8] + x413*M[17] + x419*M[18] + x428*M[32] + x434*M[14] + x454*M[31] + x490*M[24] + x492*M[27] + x495*M[28] + x510*M[33];
#pragma omp atomic
L[46] += x103*M[1] - x116*M[6] + x174*M[15] - x192*M[29] + x205*M[0] + x211*M[3] + x271*M[4] + x293*M[10] + x294*M[12] + x296*M[9] + x300*M[19] + x301*M[22] + x321*M[2] + x330*M[8] + x353*M[25] + x382*M[20] + x390*M[7] + x413*M[16] + x419*M[17] + x428*M[31] + x434*M[13] + x444*M[5] + x446*M[11] + x448*M[14] + x449*M[24] + x454*M[30] + x484*M[18] + x487*M[33] + x490*M[23] + x492*M[26] + x493*M[21] + x495*M[27] + x510*M[32] + x512*M[28];
#pragma omp atomic
L[47] += x321*M[1] + x323*M[2] + x330*M[7] + x339*M[0] + x390*M[6] + x413*M[15] + x419*M[16] + x428*M[30] + x434*M[12] + x436*M[9] + x443*M[3] + x444*M[4] + x445*M[5] + x446*M[10] + x447*M[11] + x448*M[13] + x449*M[23] + x454*M[29] + x467*M[14] + x479*M[8] + x484*M[17] + x485*M[18] + x487*M[32] + x489*M[19] + x490*M[22] + x492*M[25] + x493*M[20] + x495*M[26] + x496*M[21] + x510*M[31] + x511*M[24] + x512*M[27] + x513*M[28] + x516*M[33];
#pragma omp atomic
L[48] += x323*M[1] + x326*M[0] + x329*M[3] + x330*M[6] + x406*M[2] + x419*M[15] + x425*M[9] + x426*M[8] + x428*M[29] + x431*M[19] + x445*M[4] + x447*M[10] + x448*M[12] + x449*M[22] + x467*M[13] + x475*M[5] + x479*M[7] + x480*M[11] + x483*M[14] + x484*M[16] + x485*M[17] + x486*M[24] + x487*M[31] + x495*M[25] + x496*M[20] + x500*M[18] + x503*M[33] + x506*M[21] + x510*M[30] + x511*M[23] + x512*M[26] + x513*M[27] + x516*M[32] + x519*M[28];
#pragma omp atomic
L[49] += x366*M[0] + x406*M[1] + x407*M[2] + x426*M[7] + x440*M[9] + x467*M[12] + x471*M[14] + x474*M[3] + x475*M[4] + x477*M[5] + x479*M[6] + x480*M[10] + x481*M[11] + x483*M[13] + x484*M[15] + x485*M[16] + x486*M[23] + x487*M[30] + x497*M[8] + x500*M[17] + x501*M[18] + x503*M[32] + x504*M[19] + x506*M[20] + x507*M[21] + x510*M[29] + x511*M[22] + x512*M[25] + x513*M[26] + x515*M[24] + x516*M[31] + x519*M[27] + x521*M[28] + x525*M[33];
#pragma omp atomic
L[50] += x100*M[2] + x101*M[0] - x116*M[5] + x123*M[1] + x136*M[6] + x141*M[7] - x142*M[4] + x161*M[16] + x170*M[12] + x174*M[13] - x192*M[26] + x222*M[15] + x229*M[29] + x236*M[30] - x237*M[25] + x269*M[3] + x282*M[9] + x294*M[11] + x301*M[21] + x337*M[10] + x348*M[20] + x349*M[22] + x353*M[23] + x359*M[17] + x381*M[19] + x399*M[8] + x413*M[14] + x422*M[18] + x428*M[28] + x454*M[27] + x457*M[31] + x466*M[32] + x492*M[24] + x505*M[33];
#pragma omp atomic
L[51] += x100*M[1] + x103*M[0] - x116*M[4] + x141*M[6] + x161*M[15] + x174*M[12] - x192*M[25] + x236*M[29] + x271*M[3] + x293*M[9] + x294*M[10] + x301*M[20] + x309*M[2] + x353*M[22] + x359*M[16] + x382*M[19] + x390*M[5] + x399*M[7] + x405*M[8] + x413*M[13] + x419*M[14] + x422*M[17] + x428*M[27] + x434*M[11] + x442*M[18] + x454*M[26] + x457*M[30] + x466*M[31] + x490*M[21] + x492*M[23] + x495*M[24] + x505*M[32] + x508*M[33] + x510*M[28];
#pragma omp atomic
L[52] += x309*M[1] + x321*M[0] + x327*M[2] + x330*M[5] + x359*M[15] + x390*M[4] + x399*M[6] + x405*M[7] + x413*M[12] + x419*M[13] + x422*M[16] + x428*M[26] + x434*M[10] + x442*M[17] + x444*M[3] + x446*M[9] + x448*M[11] + x449*M[21] + x454*M[25] + x457*M[29] + x466*M[30] + x476*M[8] + x482*M[18] + x484*M[14] + x487*M[28] + x490*M[20] + x492*M[22] + x493*M[19] + x495*M[23] + x505*M[31] + x508*M[32] + x510*M[27] + x512*M[24] + x518*M[33];
#pragma omp atomic
L[53] += x323*M[0] + x327*M[1] + x330*M[4] + x365*M[2] + x405*M[6] + x419*M[12] + x422*M[15] + x428*M[25] + x442*M[16] + x445*M[3] + x447*M[9] + x448*M[10] + x449*M[20] + x466*M[29] + x467*M[11] + x470*M[18] + x476*M[7] + x478*M[8] + x479*M[5] + x482*M[17] + x484*M[13] + x485*M[14] + x487*M[27] + x495*M[22] + x496*M[19] + x505*M[30] + x508*M[31] + x510*M[26] + x511*M[21] + x512*M[23] + x513*M[24] + x516*M[28] + x518*M[32] + x522*M[33];
#pragma omp atomic
L[54] += x365*M[1] + x388*M[2] + x406*M[0] + x426*M[5] + x442*M[15] + x467*M[10] + x470*M[17] + x475*M[3] + x476*M[6] + x478*M[7] + x479*M[4] + x480*M[9] + x482*M[16] + x483*M[11] + x484*M[12] + x485*M[13] + x486*M[21] + x487*M[26] + x498*M[8] + x500*M[14] + x502*M[18] + x503*M[28] + x505*M[29] + x506*M[19] + x508*M[30] + x510*M[25] + x511*M[20] + x512*M[22] + x513*M[23] + x516*M[27] + x518*M[31] + x519*M[24] + x522*M[32] + x523*M[33];
#pragma omp atomic
L[55] += x387*M[2] + x388*M[1] + x407*M[0] + x426*M[4] + x470*M[16] + x471*M[11] + x473*M[18] + x477*M[3] + x478*M[6] + x481*M[9] + x482*M[15] + x483*M[10] + x485*M[12] + x486*M[20] + x487*M[25] + x497*M[5] + x498*M[7] + x499*M[8] + x500*M[13] + x501*M[14] + x502*M[17] + x503*M[27] + x507*M[19] + x508*M[29] + x513*M[22] + x515*M[21] + x516*M[26] + x518*M[30] + x519*M[23] + x521*M[24] + x522*M[31] + x523*M[32] + x524*M[33] + x525*M[28];
#pragma omp atomic
L[56] += x129*M[0] + x132*M[1] + x137*M[2] + x152*M[4] + x153*M[5] + x167*M[7] + x181*M[13] + x219*M[3] + x226*M[9] + x231*M[10] + x232*M[11] + x334*M[6] + x342*M[12] + x352*M[15] + x354*M[16] + x362*M[8] + x452*M[14] + x460*M[17] + x463*M[18];
#pragma omp atomic
L[57] += x112*M[2] + x132*M[0] + x152*M[3] + x167*M[5] + x181*M[11] + x231*M[9] + x261*M[1] + x289*M[6] + x290*M[7] + x300*M[16] + x334*M[4] + x342*M[10] + x352*M[12] + x354*M[13] + x379*M[15] + x410*M[8] + x431*M[18] + x460*M[14] + x489*M[17];
#pragma omp atomic
L[58] += x112*M[1] + x137*M[0] + x153*M[3] + x167*M[4] + x181*M[10] + x232*M[9] + x290*M[6] + x300*M[15] + x354*M[12] + x362*M[5] + x393*M[2] + x410*M[7] + x416*M[8] + x431*M[17] + x452*M[11] + x460*M[13] + x463*M[14] + x489*M[16] + x504*M[18];
#pragma omp atomic
L[59] += x261*M[0] + x266*M[1] + x270*M[2] + x289*M[4] + x290*M[5] + x296*M[7] + x300*M[13] + x334*M[3] + x342*M[9] + x352*M[10] + x354*M[11] + x375*M[6] + x379*M[12] + x381*M[15] + x382*M[16] + x436*M[8] + x489*M[14] + x493*M[17] + x496*M[18];
#pragma omp atomic
L[60] += x112*M[0] + x167*M[3] + x181*M[9] + x270*M[1] + x290*M[4] + x296*M[6] + x300*M[12] + x354*M[10] + x382*M[15] + x396*M[2] + x410*M[5] + x425*M[8] + x431*M[14] + x436*M[7] + x460*M[11] + x489*M[13] + x493*M[16] + x496*M[17] + x506*M[18];
#pragma omp atomic
L[61] += x362*M[3] + x393*M[0] + x396*M[1] + x402*M[2] + x410*M[4] + x416*M[5] + x425*M[7] + x431*M[13] + x436*M[6] + x440*M[8] + x452*M[9] + x460*M[10] + x463*M[11] + x489*M[12] + x493*M[15] + x496*M[16] + x504*M[14] + x506*M[17] + x507*M[18];
#pragma omp atomic
L[62] += x211*M[2] + x256*M[1] + x266*M[0] + x282*M[6] + x289*M[3] + x293*M[7] + x296*M[5] + x300*M[11] + x301*M[16] + x348*M[15] + x352*M[9] + x375*M[4] + x379*M[10] + x381*M[12] + x382*M[13] + x446*M[8] + x449*M[18] + x490*M[17] + x493*M[14];
#pragma omp atomic
L[63] += x211*M[1] + x270*M[0] + x290*M[3] + x293*M[6] + x296*M[4] + x300*M[10] + x301*M[15] + x354*M[9] + x382*M[12] + x436*M[5] + x443*M[2] + x446*M[7] + x447*M[8] + x449*M[17] + x489*M[11] + x490*M[16] + x493*M[13] + x496*M[14] + x511*M[18];
#pragma omp atomic
L[64] += x329*M[2] + x396*M[0] + x410*M[3] + x425*M[5] + x431*M[11] + x436*M[4] + x443*M[1] + x446*M[6] + x447*M[7] + x449*M[16] + x460*M[9] + x480*M[8] + x486*M[18] + x489*M[10] + x490*M[15] + x493*M[12] + x496*M[13] + x506*M[14] + x511*M[17];
#pragma omp atomic
L[65] += x329*M[1] + x402*M[0] + x416*M[3] + x425*M[4] + x431*M[10] + x440*M[5] + x447*M[6] + x449*M[15] + x463*M[9] + x474*M[2] + x480*M[7] + x481*M[8] + x486*M[17] + x496*M[12] + x504*M[11] + x506*M[13] + x507*M[14] + x511*M[16] + x515*M[18];
#pragma omp atomic
L[66] += x256*M[0] + x269*M[1] + x271*M[2] + x282*M[4] + x293*M[5] + x294*M[7] + x301*M[13] + x337*M[6] + x348*M[12] + x349*M[15] + x353*M[16] + x375*M[3] + x379*M[9] + x381*M[10] + x382*M[11] + x434*M[8] + x490*M[14] + x492*M[17] + x495*M[18];
#pragma omp atomic
L[67] += x211*M[0] + x271*M[1] + x293*M[4] + x294*M[6] + x296*M[3] + x300*M[9] + x301*M[12] + x353*M[15] + x382*M[10] + x434*M[7] + x444*M[2] + x446*M[5] + x448*M[8] + x449*M[14] + x490*M[13] + x492*M[16] + x493*M[11] + x495*M[17] + x512*M[18];
#pragma omp atomic
L[68] += x434*M[6] + x436*M[3] + x443*M[0] + x444*M[1] + x445*M[2] + x446*M[4] + x447*M[5] + x448*M[7] + x449*M[13] + x467*M[8] + x489*M[9] + x490*M[12] + x492*M[15] + x493*M[10] + x495*M[16] + x496*M[11] + x511*M[14] + x512*M[17] + x513*M[18];
#pragma omp atomic
L[69] += x329*M[0] + x425*M[3] + x431*M[9] + x445*M[1] + x447*M[4] + x448*M[6] + x449*M[12] + x467*M[7] + x475*M[2] + x480*M[5] + x483*M[8] + x486*M[14] + x495*M[15] + x496*M[10] + x506*M[11] + x511*M[13] + x512*M[16] + x513*M[17] + x519*M[18];
#pragma omp atomic
L[70] += x440*M[3] + x467*M[6] + x471*M[8] + x474*M[0] + x475*M[1] + x477*M[2] + x480*M[4] + x481*M[5] + x483*M[7] + x486*M[13] + x504*M[9] + x506*M[10] + x507*M[11] + x511*M[12] + x512*M[15] + x513*M[16] + x515*M[14] + x519*M[17] + x521*M[18];
#pragma omp atomic
L[71] += -x116*M[2] - x142*M[1] + x170*M[6] + x174*M[7] - x192*M[16] - x237*M[15] + x269*M[0] + x282*M[3] + x294*M[5] + x301*M[11] + x337*M[4] + x348*M[10] + x349*M[12] + x353*M[13] + x381*M[9] + x413*M[8] + x428*M[18] + x454*M[17] + x492*M[14];
#pragma omp atomic
L[72] += -x116*M[1] + x174*M[6] - x192*M[15] + x271*M[0] + x293*M[3] + x294*M[4] + x301*M[10] + x353*M[12] + x382*M[9] + x390*M[2] + x413*M[7] + x419*M[8] + x428*M[17] + x434*M[5] + x454*M[16] + x490*M[11] + x492*M[13] + x495*M[14] + x510*M[18];
#pragma omp atomic
L[73] += x330*M[2] + x390*M[1] + x413*M[6] + x419*M[7] + x428*M[16] + x434*M[4] + x444*M[0] + x446*M[3] + x448*M[5] + x449*M[11] + x454*M[15] + x484*M[8] + x487*M[18] + x490*M[10] + x492*M[12] + x493*M[9] + x495*M[13] + x510*M[17] + x512*M[14];
#pragma omp atomic
L[74] += x330*M[1] + x419*M[6] + x428*M[15] + x445*M[0] + x447*M[3] + x448*M[4] + x449*M[10] + x467*M[5] + x479*M[2] + x484*M[7] + x485*M[8] + x487*M[17] + x495*M[12] + x496*M[9] + x510*M[16] + x511*M[11] + x512*M[13] + x513*M[14] + x516*M[18];
#pragma omp atomic
L[75] += x426*M[2] + x467*M[4] + x475*M[0] + x479*M[1] + x480*M[3] + x483*M[5] + x484*M[6] + x485*M[7] + x486*M[11] + x487*M[16] + x500*M[8] + x503*M[18] + x506*M[9] + x510*M[15] + x511*M[10] + x512*M[12] + x513*M[13] + x516*M[17] + x519*M[14];
#pragma omp atomic
L[76] += x426*M[1] + x471*M[5] + x477*M[0] + x481*M[3] + x483*M[4] + x485*M[6] + x486*M[10] + x487*M[15] + x497*M[2] + x500*M[7] + x501*M[8] + x503*M[17] + x507*M[9] + x513*M[12] + x515*M[11] + x516*M[16] + x519*M[13] + x521*M[14] + x525*M[18];
#pragma omp atomic
L[77] += x136*M[1] + x141*M[2] - x142*M[0] + x161*M[7] + x170*M[4] + x174*M[5] - x192*M[13] + x222*M[6] + x229*M[15] + x236*M[16] - x237*M[12] + x337*M[3] + x348*M[9] + x349*M[10] + x353*M[11] + x359*M[8] + x454*M[14] + x457*M[17] + x466*M[18];
#pragma omp atomic
L[78] += -x116*M[0] + x141*M[1] + x161*M[6] + x174*M[4] - x192*M[12] + x236*M[15] + x294*M[3] + x301*M[9] + x353*M[10] + x359*M[7] + x399*M[2] + x413*M[5] + x422*M[8] + x428*M[14] + x454*M[13] + x457*M[16] + x466*M[17] + x492*M[11] + x505*M[18];
#pragma omp atomic
L[79] += x359*M[6] + x390*M[0] + x399*M[1] + x405*M[2] + x413*M[4] + x419*M[5] + x422*M[7] + x428*M[13] + x434*M[3] + x442*M[8] + x454*M[12] + x457*M[15] + x466*M[16] + x490*M[9] + x492*M[10] + x495*M[11] + x505*M[17] + x508*M[18] + x510*M[14];
#pragma omp atomic
L[80] += x330*M[0] + x405*M[1] + x419*M[4] + x422*M[6] + x428*M[12] + x442*M[7] + x448*M[3] + x449*M[9] + x466*M[15] + x476*M[2] + x482*M[8] + x484*M[5] + x487*M[14] + x495*M[10] + x505*M[16] + x508*M[17] + x510*M[13] + x512*M[11] + x518*M[18];
#pragma omp atomic
L[81] += x442*M[6] + x467*M[3] + x470*M[8] + x476*M[1] + x478*M[2] + x479*M[0] + x482*M[7] + x484*M[4] + x485*M[5] + x487*M[13] + x505*M[15] + x508*M[16] + x510*M[12] + x511*M[9] + x512*M[10] + x513*M[11] + x516*M[14] + x518*M[17] + x522*M[18];
#pragma omp atomic
L[82] += x426*M[0] + x470*M[7] + x478*M[1] + x482*M[6] + x483*M[3] + x485*M[4] + x486*M[9] + x487*M[12] + x498*M[2] + x500*M[5] + x502*M[8] + x503*M[14] + x508*M[15] + x513*M[10] + x516*M[13] + x518*M[16] + x519*M[11] + x522*M[17] + x523*M[18];
#pragma omp atomic
L[83] += x470*M[6] + x471*M[3] + x473*M[8] + x497*M[0] + x498*M[1] + x499*M[2] + x500*M[4] + x501*M[5] + x502*M[7] + x503*M[13] + x515*M[9] + x516*M[12] + x518*M[15] + x519*M[10] + x521*M[11] + x522*M[16] + x523*M[17] + x524*M[18] + x525*M[14];
#pragma omp atomic
L[84] += x152*M[1] + x153*M[2] + x181*M[7] + x219*M[0] + x226*M[3] + x231*M[4] + x232*M[5] + x342*M[6] + x452*M[8];
#pragma omp atomic
L[85] += x152*M[0] + x167*M[2] + x181*M[5] + x231*M[3] + x334*M[1] + x342*M[4] + x352*M[6] + x354*M[7] + x460*M[8];
#pragma omp atomic
L[86] += x153*M[0] + x167*M[1] + x181*M[4] + x232*M[3] + x354*M[6] + x362*M[2] + x452*M[5] + x460*M[7] + x463*M[8];
#pragma omp atomic
L[87] += x289*M[1] + x290*M[2] + x300*M[7] + x334*M[0] + x342*M[3] + x352*M[4] + x354*M[5] + x379*M[6] + x489*M[8];
#pragma omp atomic
L[88] += x167*M[0] + x181*M[3] + x290*M[1] + x300*M[6] + x354*M[4] + x410*M[2] + x431*M[8] + x460*M[5] + x489*M[7];
#pragma omp atomic
L[89] += x362*M[0] + x410*M[1] + x416*M[2] + x431*M[7] + x452*M[3] + x460*M[4] + x463*M[5] + x489*M[6] + x504*M[8];
#pragma omp atomic
L[90] += x289*M[0] + x296*M[2] + x300*M[5] + x352*M[3] + x375*M[1] + x379*M[4] + x381*M[6] + x382*M[7] + x493*M[8];
#pragma omp atomic
L[91] += x290*M[0] + x296*M[1] + x300*M[4] + x354*M[3] + x382*M[6] + x436*M[2] + x489*M[5] + x493*M[7] + x496*M[8];
#pragma omp atomic
L[92] += x410*M[0] + x425*M[2] + x431*M[5] + x436*M[1] + x460*M[3] + x489*M[4] + x493*M[6] + x496*M[7] + x506*M[8];
#pragma omp atomic
L[93] += x416*M[0] + x425*M[1] + x431*M[4] + x440*M[2] + x463*M[3] + x496*M[6] + x504*M[5] + x506*M[7] + x507*M[8];
#pragma omp atomic
L[94] += x282*M[1] + x293*M[2] + x301*M[7] + x348*M[6] + x375*M[0] + x379*M[3] + x381*M[4] + x382*M[5] + x490*M[8];
#pragma omp atomic
L[95] += x293*M[1] + x296*M[0] + x300*M[3] + x301*M[6] + x382*M[4] + x446*M[2] + x449*M[8] + x490*M[7] + x493*M[5];
#pragma omp atomic
L[96] += x436*M[0] + x446*M[1] + x447*M[2] + x449*M[7] + x489*M[3] + x490*M[6] + x493*M[4] + x496*M[5] + x511*M[8];
#pragma omp atomic
L[97] += x425*M[0] + x431*M[3] + x447*M[1] + x449*M[6] + x480*M[2] + x486*M[8] + x496*M[4] + x506*M[5] + x511*M[7];
#pragma omp atomic
L[98] += x440*M[0] + x480*M[1] + x481*M[2] + x486*M[7] + x504*M[3] + x506*M[4] + x507*M[5] + x511*M[6] + x515*M[8];
#pragma omp atomic
L[99] += x282*M[0] + x294*M[2] + x301*M[5] + x337*M[1] + x348*M[4] + x349*M[6] + x353*M[7] + x381*M[3] + x492*M[8];
#pragma omp atomic
L[100] += x293*M[0] + x294*M[1] + x301*M[4] + x353*M[6] + x382*M[3] + x434*M[2] + x490*M[5] + x492*M[7] + x495*M[8];
#pragma omp atomic
L[101] += x434*M[1] + x446*M[0] + x448*M[2] + x449*M[5] + x490*M[4] + x492*M[6] + x493*M[3] + x495*M[7] + x512*M[8];
#pragma omp atomic
L[102] += x447*M[0] + x448*M[1] + x449*M[4] + x467*M[2] + x495*M[6] + x496*M[3] + x511*M[5] + x512*M[7] + x513*M[8];
#pragma omp atomic
L[103] += x467*M[1] + x480*M[0] + x483*M[2] + x486*M[5] + x506*M[3] + x511*M[4] + x512*M[6] + x513*M[7] + x519*M[8];
#pragma omp atomic
L[104] += x471*M[2] + x481*M[0] + x483*M[1] + x486*M[4] + x507*M[3] + x513*M[6] + x515*M[5] + x519*M[7] + x521*M[8];
#pragma omp atomic
L[105] += x170*M[1] + x174*M[2] - x192*M[7] - x237*M[6] + x337*M[0] + x348*M[3] + x349*M[4] + x353*M[5] + x454*M[8];
#pragma omp atomic
L[106] += x174*M[1] - x192*M[6] + x294*M[0] + x301*M[3] + x353*M[4] + x413*M[2] + x428*M[8] + x454*M[7] + x492*M[5];
#pragma omp atomic
L[107] += x413*M[1] + x419*M[2] + x428*M[7] + x434*M[0] + x454*M[6] + x490*M[3] + x492*M[4] + x495*M[5] + x510*M[8];
#pragma omp atomic
L[108] += x419*M[1] + x428*M[6] + x448*M[0] + x449*M[3] + x484*M[2] + x487*M[8] + x495*M[4] + x510*M[7] + x512*M[5];
#pragma omp atomic
L[109] += x467*M[0] + x484*M[1] + x485*M[2] + x487*M[7] + x510*M[6] + x511*M[3] + x512*M[4] + x513*M[5] + x516*M[8];
#pragma omp atomic
L[110] += x483*M[0] + x485*M[1] + x486*M[3] + x487*M[6] + x500*M[2] + x503*M[8] + x513*M[4] + x516*M[7] + x519*M[5];
#pragma omp atomic
L[111] += x471*M[0] + x500*M[1] + x501*M[2] + x503*M[7] + x515*M[3] + x516*M[6] + x519*M[4] + x521*M[5] + x525*M[8];
#pragma omp atomic
L[112] += x161*M[2] + x170*M[0] - x192*M[5] + x222*M[1] + x229*M[6] + x236*M[7] - x237*M[4] + x349*M[3] + x457*M[8];
#pragma omp atomic
L[113] += x161*M[1] + x174*M[0] - x192*M[4] + x236*M[6] + x353*M[3] + x359*M[2] + x454*M[5] + x457*M[7] + x466*M[8];
#pragma omp atomic
L[114] += x359*M[1] + x413*M[0] + x422*M[2] + x428*M[5] + x454*M[4] + x457*M[6] + x466*M[7] + x492*M[3] + x505*M[8];
#pragma omp atomic
L[115] += x419*M[0] + x422*M[1] + x428*M[4] + x442*M[2] + x466*M[6] + x495*M[3] + x505*M[7] + x508*M[8] + x510*M[5];
#pragma omp atomic
L[116] += x442*M[1] + x482*M[2] + x484*M[0] + x487*M[5] + x505*M[6] + x508*M[7] + x510*M[4] + x512*M[3] + x518*M[8];
#pragma omp atomic
L[117] += x470*M[2] + x482*M[1] + x485*M[0] + x487*M[4] + x508*M[6] + x513*M[3] + x516*M[5] + x518*M[7] + x522*M[8];
#pragma omp atomic
L[118] += x470*M[1] + x500*M[0] + x502*M[2] + x503*M[5] + x516*M[4] + x518*M[6] + x519*M[3] + x522*M[7] + x523*M[8];
#pragma omp atomic
L[119] += x473*M[2] + x501*M[0] + x502*M[1] + x503*M[4] + x521*M[3] + x522*M[6] + x523*M[7] + x524*M[8] + x525*M[5];
#pragma omp atomic
L[120] += x226*M[0] + x231*M[1] + x232*M[2];
#pragma omp atomic
L[121] += x181*M[2] + x231*M[0] + x342*M[1];
#pragma omp atomic
L[122] += x181*M[1] + x232*M[0] + x452*M[2];
#pragma omp atomic
L[123] += x342*M[0] + x352*M[1] + x354*M[2];
#pragma omp atomic
L[124] += x181*M[0] + x354*M[1] + x460*M[2];
#pragma omp atomic
L[125] += x452*M[0] + x460*M[1] + x463*M[2];
#pragma omp atomic
L[126] += x300*M[2] + x352*M[0] + x379*M[1];
#pragma omp atomic
L[127] += x300*M[1] + x354*M[0] + x489*M[2];
#pragma omp atomic
L[128] += x431*M[2] + x460*M[0] + x489*M[1];
#pragma omp atomic
L[129] += x431*M[1] + x463*M[0] + x504*M[2];
#pragma omp atomic
L[130] += x379*M[0] + x381*M[1] + x382*M[2];
#pragma omp atomic
L[131] += x300*M[0] + x382*M[1] + x493*M[2];
#pragma omp atomic
L[132] += x489*M[0] + x493*M[1] + x496*M[2];
#pragma omp atomic
L[133] += x431*M[0] + x496*M[1] + x506*M[2];
#pragma omp atomic
L[134] += x504*M[0] + x506*M[1] + x507*M[2];
#pragma omp atomic
L[135] += x301*M[2] + x348*M[1] + x381*M[0];
#pragma omp atomic
L[136] += x301*M[1] + x382*M[0] + x490*M[2];
#pragma omp atomic
L[137] += x449*M[2] + x490*M[1] + x493*M[0];
#pragma omp atomic
L[138] += x449*M[1] + x496*M[0] + x511*M[2];
#pragma omp atomic
L[139] += x486*M[2] + x506*M[0] + x511*M[1];
#pragma omp atomic
L[140] += x486*M[1] + x507*M[0] + x515*M[2];
#pragma omp atomic
L[141] += x348*M[0] + x349*M[1] + x353*M[2];
#pragma omp atomic
L[142] += x301*M[0] + x353*M[1] + x492*M[2];
#pragma omp atomic
L[143] += x490*M[0] + x492*M[1] + x495*M[2];
#pragma omp atomic
L[144] += x449*M[0] + x495*M[1] + x512*M[2];
#pragma omp atomic
L[145] += x511*M[0] + x512*M[1] + x513*M[2];
#pragma omp atomic
L[146] += x486*M[0] + x513*M[1] + x519*M[2];
#pragma omp atomic
L[147] += x515*M[0] + x519*M[1] + x521*M[2];
#pragma omp atomic
L[148] += -x192*M[2] - x237*M[1] + x349*M[0];
#pragma omp atomic
L[149] += -x192*M[1] + x353*M[0] + x454*M[2];
#pragma omp atomic
L[150] += x428*M[2] + x454*M[1] + x492*M[0];
#pragma omp atomic
L[151] += x428*M[1] + x495*M[0] + x510*M[2];
#pragma omp atomic
L[152] += x487*M[2] + x510*M[1] + x512*M[0];
#pragma omp atomic
L[153] += x487*M[1] + x513*M[0] + x516*M[2];
#pragma omp atomic
L[154] += x503*M[2] + x516*M[1] + x519*M[0];
#pragma omp atomic
L[155] += x503*M[1] + x521*M[0] + x525*M[2];
#pragma omp atomic
L[156] += x229*M[1] + x236*M[2] - x237*M[0];
#pragma omp atomic
L[157] += -x192*M[0] + x236*M[1] + x457*M[2];
#pragma omp atomic
L[158] += x454*M[0] + x457*M[1] + x466*M[2];
#pragma omp atomic
L[159] += x428*M[0] + x466*M[1] + x505*M[2];
#pragma omp atomic
L[160] += x505*M[1] + x508*M[2] + x510*M[0];
#pragma omp atomic
L[161] += x487*M[0] + x508*M[1] + x518*M[2];
#pragma omp atomic
L[162] += x516*M[0] + x518*M[1] + x522*M[2];
#pragma omp atomic
L[163] += x503*M[0] + x522*M[1] + x523*M[2];
#pragma omp atomic
L[164] += x523*M[1] + x524*M[2] + x525*M[0];

}

void L2L_9(double x, double y, double z, double * L, double * Ls) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
double x448;
double x449;
double x450;
double x451;
double x452;
double x453;
double x454;
double x455;
double x456;
double x457;
double x458;
double x459;
double x460;
double x461;
double x462;
double x463;
double x464;
double x465;
double x466;
double x467;
double x468;
double x469;
double x470;
double x471;
double x472;
double x473;
double x474;
double x475;
double x476;
double x477;
double x478;
double x479;
double x480;
double x481;
double x482;
double x483;
double x484;
double x485;
double x486;
double x487;
double x488;
double x489;
double x490;
double x491;
double x492;
double x493;
double x494;
double x495;
double x496;
double x497;
double x498;
double x499;
double x500;
double x501;
double x502;
double x503;
double x504;
double x505;
double x506;
double x507;
double x508;
double x509;
double x510;
double x511;
double x512;
double x513;
double x514;
double x515;
double x516;
double x517;
double x518;
double x519;
double x520;
double x521;
double x522;
double x523;
double x524;
double x525;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
x3 = z*L[14];
x4 = x3*y;
x5 = (x*x);
x6 = (1.0/2.0)*x5;
x7 = (x*x*x);
x8 = (1.0/6.0)*x7;
x9 = (x*x*x*x);
x10 = (1.0/24.0)*x9;
x11 = pow(x, 5);
x12 = (1.0/120.0)*x11;
x13 = pow(x, 6);
x14 = (1.0/720.0)*x13;
x15 = (1.0/5040.0)*pow(x, 7);
x16 = (y*y);
x17 = (1.0/2.0)*x16;
x18 = (y*y*y);
x19 = (1.0/6.0)*x18;
x20 = (y*y*y*y);
x21 = (1.0/24.0)*x20;
x22 = pow(y, 5);
x23 = (1.0/120.0)*x22;
x24 = pow(y, 6);
x25 = (1.0/720.0)*x24;
x26 = (1.0/5040.0)*pow(y, 7);
x27 = (z*z);
x28 = (1.0/2.0)*x27;
x29 = (z*z*z);
x30 = (1.0/6.0)*x29;
x31 = (z*z*z*z);
x32 = (1.0/24.0)*x31;
x33 = pow(z, 5);
x34 = (1.0/120.0)*x33;
x35 = pow(z, 6);
x36 = (1.0/720.0)*x35;
x37 = (1.0/5040.0)*pow(z, 7);
x38 = x*L[13];
x39 = x*L[26];
x40 = x*L[45];
x41 = x*L[71];
x42 = x*L[105];
x43 = x*L[148];
x44 = x*L[15];
x45 = x*L[29];
x46 = x*L[49];
x47 = x*L[76];
x48 = x*L[111];
x49 = x*L[155];
x50 = y*L[11];
x51 = z*L[12];
x52 = y*L[21];
x53 = z*L[22];
x54 = y*L[36];
x55 = z*L[37];
x56 = y*L[57];
x57 = z*L[58];
x58 = y*L[85];
x59 = z*L[86];
x60 = y*L[121];
x61 = z*L[122];
x62 = y*L[18];
x63 = y*L[33];
x64 = y*L[54];
x65 = y*L[82];
x66 = y*L[118];
x67 = y*L[163];
x68 = z*L[17];
x69 = z*L[31];
x70 = z*L[51];
x71 = z*L[78];
x72 = z*L[113];
x73 = z*L[157];
x74 = y*L[28];
x75 = x*x74;
x76 = y*L[48];
x77 = x*x76;
x78 = y*L[75];
x79 = x*x78;
x80 = y*L[110];
x81 = x*x80;
x82 = y*L[154];
x83 = x*x82;
x84 = z*L[27];
x85 = x*x84;
x86 = z*L[46];
x87 = x*x86;
x88 = z*L[72];
x89 = x*x88;
x90 = z*L[106];
x91 = x*x90;
x92 = z*L[149];
x93 = x*x92;
x94 = z*L[24];
x95 = x94*y;
x96 = z*L[39];
x97 = x96*y;
x98 = z*L[60];
x99 = x98*y;
x100 = z*L[88];
x101 = x100*y;
x102 = z*L[124];
x103 = x102*y;
x104 = (1.0/4.0)*x5;
x105 = x104*x16;
x106 = (1.0/12.0)*x5;
x107 = x106*x18;
x108 = (1.0/48.0)*x5;
x109 = x108*x20;
x110 = (1.0/240.0)*x5;
x111 = x110*x22;
x112 = (1.0/1440.0)*x5;
x113 = x104*x27;
x114 = x106*x29;
x115 = x108*x31;
x116 = x110*x33;
x117 = (1.0/12.0)*x7;
x118 = x117*x16;
x119 = (1.0/36.0)*x7;
x120 = x119*x18;
x121 = (1.0/144.0)*x7;
x122 = x121*x20;
x123 = (1.0/720.0)*x7;
x124 = x117*x27;
x125 = x119*x29;
x126 = x121*x31;
x127 = (1.0/48.0)*x9;
x128 = x127*x16;
x129 = (1.0/144.0)*x9;
x130 = x129*x18;
x131 = (1.0/576.0)*x9;
x132 = x127*x27;
x133 = x129*x29;
x134 = (1.0/240.0)*x11;
x135 = x134*x16;
x136 = (1.0/720.0)*x11;
x137 = x134*x27;
x138 = (1.0/1440.0)*x13;
x139 = x16*x27;
x140 = (1.0/4.0)*x139;
x141 = x16*x29;
x142 = (1.0/12.0)*x141;
x143 = x16*x31;
x144 = (1.0/48.0)*x143;
x145 = (1.0/240.0)*x16*x33;
x146 = x18*x27;
x147 = (1.0/12.0)*x146;
x148 = x18*x29;
x149 = (1.0/36.0)*x148;
x150 = (1.0/144.0)*x18*x31;
x151 = x20*x27;
x152 = (1.0/48.0)*x151;
x153 = (1.0/144.0)*x20*x29;
x154 = (1.0/240.0)*x22*x27;
x155 = x*L[47];
x156 = x*L[74];
x157 = x*L[109];
x158 = x*L[153];
x159 = x*L[73];
x160 = x*L[108];
x161 = x*L[152];
x162 = x*L[107];
x163 = x*L[151];
x164 = x*L[150];
x165 = y*L[43];
x166 = y*L[69];
x167 = y*L[103];
x168 = y*L[146];
x169 = z*L[42];
x170 = z*L[67];
x171 = z*L[100];
x172 = z*L[142];
x173 = y*L[64];
x174 = y*L[97];
x175 = y*L[139];
x176 = z*L[63];
x177 = z*L[95];
x178 = z*L[136];
x179 = y*L[92];
x180 = y*L[133];
x181 = z*L[91];
x182 = z*L[131];
x183 = y*L[128];
x184 = z*L[127];
x185 = (1.0/8.0)*x139*x5;
x186 = (1.0/24.0)*x5;
x187 = x141*x186;
x188 = (1.0/96.0)*x5;
x189 = x146*x186;
x190 = (1.0/24.0)*x139*x7;
x191 = (1.0/72.0)*x7;
x192 = x*L[23];
x193 = x*L[41];
x194 = x*L[66];
x195 = x*L[99];
x196 = x*L[141];
x197 = x*L[25];
x198 = x*L[44];
x199 = x*L[70];
x200 = x*L[104];
x201 = x*L[147];
x202 = x*x165;
x203 = x*x166;
x204 = x*x167;
x205 = x*x168;
x206 = x*x169;
x207 = x*x170;
x208 = x*x171;
x209 = x*x172;
x210 = x*L[68];
x211 = x*L[102];
x212 = x*L[145];
x213 = x*L[101];
x214 = x*L[144];
x215 = x*L[143];
x216 = y*L[13];
x217 = x84*y;
x218 = x*L[28];
x219 = x*L[48];
x220 = x*L[75];
x221 = x*L[110];
x222 = x*L[154];
x223 = y*L[23];
x224 = y*L[38];
x225 = y*L[59];
x226 = y*L[87];
x227 = y*L[123];
x228 = y*L[32];
x229 = y*L[53];
x230 = y*L[81];
x231 = y*L[117];
x232 = y*L[162];
x233 = y*L[47];
x234 = x*x233;
x235 = y*L[74];
x236 = x*x235;
x237 = y*L[109];
x238 = x*x237;
x239 = y*L[153];
x240 = x*x239;
x241 = x169*y;
x242 = x176*y;
x243 = x181*y;
x244 = x184*y;
x245 = y*L[68];
x246 = y*L[102];
x247 = y*L[145];
x248 = y*L[96];
x249 = y*L[138];
x250 = y*L[132];
x251 = y*L[14];
x252 = z*L[15];
x253 = z*L[18];
x254 = z*L[28];
x255 = x254*y;
x256 = x*L[27];
x257 = x*L[46];
x258 = x*L[72];
x259 = x*L[106];
x260 = x*L[149];
x261 = y*L[24];
x262 = z*L[25];
x263 = y*L[39];
x264 = z*L[40];
x265 = y*L[60];
x266 = z*L[61];
x267 = y*L[88];
x268 = z*L[89];
x269 = y*L[124];
x270 = z*L[125];
x271 = z*L[32];
x272 = z*L[52];
x273 = z*L[79];
x274 = z*L[114];
x275 = z*L[158];
x276 = z*L[47];
x277 = x*x276;
x278 = z*L[73];
x279 = x*x278;
x280 = z*L[107];
x281 = x*x280;
x282 = z*L[150];
x283 = x*x282;
x284 = z*L[43];
x285 = x284*y;
x286 = z*L[64];
x287 = x286*y;
x288 = z*L[92];
x289 = x288*y;
x290 = z*L[128];
x291 = x290*y;
x292 = z*L[68];
x293 = z*L[101];
x294 = z*L[143];
x295 = z*L[96];
x296 = z*L[137];
x297 = z*L[132];
x298 = x*L[38];
x299 = x*L[62];
x300 = x*L[94];
x301 = x*L[135];
x302 = x*L[40];
x303 = x*L[65];
x304 = x*L[98];
x305 = x*L[140];
x306 = x*x173;
x307 = x*x174;
x308 = x*x175;
x309 = x*x176;
x310 = x*x177;
x311 = x*x178;
x312 = x*L[96];
x313 = x*L[138];
x314 = x*L[137];
x315 = x*L[43];
x316 = x*L[69];
x317 = x*L[103];
x318 = x*L[146];
x319 = x*x245;
x320 = x*x246;
x321 = x*x247;
x322 = x*L[42];
x323 = x*L[67];
x324 = x*L[100];
x325 = x*L[142];
x326 = x*x292;
x327 = x*x293;
x328 = x*x294;
x329 = y*L[26];
x330 = x86*y;
x331 = y*L[41];
x332 = y*L[62];
x333 = y*L[90];
x334 = y*L[126];
x335 = y*L[52];
x336 = y*L[80];
x337 = y*L[116];
x338 = y*L[161];
x339 = y*L[73];
x340 = x*x339;
x341 = y*L[108];
x342 = x*x341;
x343 = y*L[152];
x344 = x*x343;
x345 = x170*y;
x346 = x177*y;
x347 = x182*y;
x348 = y*L[101];
x349 = y*L[144];
x350 = y*L[137];
x351 = y*L[27];
x352 = x276*y;
x353 = y*L[42];
x354 = y*L[63];
x355 = y*L[91];
x356 = y*L[127];
x357 = x292*y;
x358 = x295*y;
x359 = x297*y;
x360 = z*L[29];
x361 = z*L[33];
x362 = z*L[48];
x363 = x362*y;
x364 = z*L[44];
x365 = z*L[65];
x366 = z*L[93];
x367 = z*L[129];
x368 = z*L[53];
x369 = z*L[80];
x370 = z*L[115];
x371 = z*L[159];
x372 = z*L[74];
x373 = x*x372;
x374 = z*L[108];
x375 = x*x374;
x376 = z*L[151];
x377 = x*x376;
x378 = z*L[69];
x379 = x378*y;
x380 = z*L[97];
x381 = x380*y;
x382 = z*L[133];
x383 = x382*y;
x384 = z*L[102];
x385 = z*L[144];
x386 = z*L[138];
x387 = x*L[59];
x388 = x*L[90];
x389 = x*L[130];
x390 = x*L[61];
x391 = x*L[93];
x392 = x*L[134];
x393 = x*x179;
x394 = x*x180;
x395 = x*x181;
x396 = x*x182;
x397 = x*L[132];
x398 = x*L[64];
x399 = x*L[97];
x400 = x*L[139];
x401 = x*x248;
x402 = x*x249;
x403 = x*L[63];
x404 = x*L[95];
x405 = x*L[136];
x406 = x*x295;
x407 = x*x296;
x408 = x*x348;
x409 = x*x349;
x410 = x*x384;
x411 = x*x385;
x412 = y*L[45];
x413 = x88*y;
x414 = y*L[66];
x415 = y*L[94];
x416 = y*L[130];
x417 = y*L[79];
x418 = y*L[115];
x419 = y*L[160];
x420 = y*L[107];
x421 = x*x420;
x422 = y*L[151];
x423 = x*x422;
x424 = x171*y;
x425 = x178*y;
x426 = y*L[143];
x427 = y*L[46];
x428 = x278*y;
x429 = y*L[67];
x430 = y*L[95];
x431 = y*L[131];
x432 = x293*y;
x433 = x296*y;
x434 = x372*y;
x435 = x384*y;
x436 = x386*y;
x437 = z*L[49];
x438 = z*L[54];
x439 = z*L[75];
x440 = x439*y;
x441 = z*L[70];
x442 = z*L[98];
x443 = z*L[134];
x444 = z*L[81];
x445 = z*L[116];
x446 = z*L[160];
x447 = z*L[109];
x448 = x*x447;
x449 = z*L[152];
x450 = x*x449;
x451 = z*L[103];
x452 = x451*y;
x453 = z*L[139];
x454 = x453*y;
x455 = z*L[145];
x456 = x*L[87];
x457 = x*L[126];
x458 = x*L[89];
x459 = x*L[129];
x460 = x*x183;
x461 = x*x184;
x462 = x*L[92];
x463 = x*L[133];
x464 = x*x250;
x465 = x*L[91];
x466 = x*L[131];
x467 = x*x297;
x468 = x*x350;
x469 = x*x386;
x470 = x*x426;
x471 = x*x455;
x472 = y*L[71];
x473 = x90*y;
x474 = y*L[99];
x475 = y*L[135];
x476 = y*L[114];
x477 = y*L[159];
x478 = y*L[150];
x479 = x*x478;
x480 = x172*y;
x481 = y*L[72];
x482 = x280*y;
x483 = y*L[100];
x484 = y*L[136];
x485 = x294*y;
x486 = x374*y;
x487 = x385*y;
x488 = x447*y;
x489 = x455*y;
x490 = z*L[76];
x491 = z*L[82];
x492 = z*L[110];
x493 = x492*y;
x494 = z*L[104];
x495 = z*L[140];
x496 = z*L[117];
x497 = z*L[161];
x498 = z*L[153];
x499 = x*x498;
x500 = z*L[146];
x501 = x500*y;
x502 = x*L[123];
x503 = x*L[125];
x504 = x*L[128];
x505 = x*L[127];
x506 = y*L[105];
x507 = x92*y;
x508 = y*L[141];
x509 = y*L[158];
x510 = y*L[106];
x511 = x282*y;
x512 = y*L[142];
x513 = x376*y;
x514 = x449*y;
x515 = x498*y;
x516 = z*L[111];
x517 = z*L[118];
x518 = z*L[154];
x519 = x518*y;
x520 = z*L[147];
x521 = z*L[162];
x522 = y*L[148];
x523 = y*L[149];
x524 = z*L[155];
x525 = z*L[163];
#pragma omp atomic
Ls[0] += (1.0/40320.0)*pow(x, 8)*L[120] + x*x0 + x*x1 + x*x4 + x*L[1] + x10*x54 + x10*x55 + x10*x99 + x10*L[20] + x101*x12 + x103*x14 + x105*x169 + x105*L[23] + x107*x170 + x107*L[41] + x109*x171 + x109*L[66] + x111*x172 + x111*L[99] + x112*x24*L[141] + x112*x35*L[147] + x113*x165 + x113*L[25] + x114*x166 + x114*L[44] + x115*x167 + x115*L[70] + x116*x168 + x116*L[104] + x118*x176 + x118*L[38] + x12*x56 + x12*x57 + x12*L[35] + x120*x177 + x120*L[62] + x122*x178 + x122*L[94] + x123*x22*L[135] + x123*x33*L[140] + x124*x173 + x124*L[40] + x125*x174 + x125*L[65] + x126*x175 + x126*L[98] + x128*x181 + x128*L[59] + x130*x182 + x130*L[90] + x131*x20*L[130] + x131*x31*L[134] + x132*x179 + x132*L[61] + x133*x180 + x133*L[93] + x135*x184 + x135*L[87] + x136*x18*L[126] + x136*x29*L[129] + x137*x183 + x137*L[89] + x138*x16*L[123] + x138*x27*L[125] + (1.0/96.0)*x139*x9*L[132] + x14*x58 + x14*x59 + x14*L[56] + x140*x155 + x140*L[32] + x141*x191*L[138] + x142*x156 + x142*L[53] + x143*x188*L[145] + x144*x157 + x144*L[81] + x145*x158 + x145*L[117] + x146*x191*L[137] + x147*x159 + x147*L[52] + (1.0/72.0)*x148*x5*L[144] + x149*x160 + x149*L[80] + x15*x60 + x15*x61 + x15*L[84] + x150*x161 + x150*L[116] + x151*x188*L[143] + x152*x162 + x152*L[79] + x153*x163 + x153*L[115] + x154*x164 + x154*L[114] + (1.0/1440.0)*x16*x35*L[162] + x17*x38 + x17*x68 + x17*x85 + x17*L[7] + (1.0/720.0)*x18*x33*L[161] + x185*L[68] + x187*L[102] + x189*L[101] + x19*x39 + x19*x69 + x19*x87 + x19*L[16] + x190*L[96] + x2*y + (1.0/576.0)*x20*x31*L[160] + x21*x40 + x21*x70 + x21*x89 + x21*L[30] + (1.0/720.0)*x22*x29*L[159] + x23*x41 + x23*x71 + x23*x91 + x23*L[50] + (1.0/1440.0)*x24*x27*L[158] + x25*x42 + x25*x72 + x25*x93 + x25*L[77] + x26*x43 + x26*x73 + x26*L[112] + x28*x44 + x28*x62 + x28*x75 + x28*L[9] + x30*x45 + x30*x63 + x30*x77 + x30*L[19] + x32*x46 + x32*x64 + x32*x79 + x32*L[34] + x34*x47 + x34*x65 + x34*x81 + x34*L[55] + x36*x48 + x36*x66 + x36*x83 + x36*L[83] + x37*x49 + x37*x67 + x37*L[119] + x50*x6 + x51*x6 + x52*x8 + x53*x8 + x6*x95 + x6*L[4] + x8*x97 + x8*L[10] + (1.0/40320.0)*pow(y, 8)*L[156] + y*L[2] + (1.0/40320.0)*pow(z, 8)*L[164] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*x50 + x*x51 + x*x95 + x*L[4] + x0 + x1 + x10*x101 + x10*x56 + x10*x57 + x10*L[35] + x103*x12 + x105*x176 + x105*L[38] + x107*x177 + x107*L[62] + x109*x178 + x109*L[94] + x111*L[135] + x113*x173 + x113*L[40] + x114*x174 + x114*L[65] + x115*x175 + x115*L[98] + x116*L[140] + x118*x181 + x118*L[59] + x12*x58 + x12*x59 + x12*L[56] + x120*x182 + x120*L[90] + x122*L[130] + x124*x179 + x124*L[61] + x125*x180 + x125*L[93] + x126*L[134] + x128*x184 + x128*L[87] + x130*L[126] + x132*x183 + x132*L[89] + x133*L[129] + x135*L[123] + x137*L[125] + x14*x60 + x14*x61 + x14*L[84] + x140*x210 + x140*L[47] + x142*x211 + x142*L[74] + x144*x212 + x144*L[109] + x145*L[153] + x147*x213 + x147*L[73] + x149*x214 + x149*L[108] + x15*L[120] + x150*L[152] + x152*x215 + x152*L[107] + x153*L[151] + x154*L[150] + x17*x192 + x17*x206 + x17*x84 + x17*L[13] + x185*L[96] + x187*L[138] + x189*L[137] + x19*x193 + x19*x207 + x19*x86 + x19*L[26] + x190*L[132] + x194*x21 + x195*x23 + x196*x25 + x197*x28 + x198*x30 + x199*x32 + x200*x34 + x201*x36 + x202*x28 + x203*x30 + x204*x32 + x205*x34 + x208*x21 + x209*x23 + x21*x88 + x21*L[45] + x23*x90 + x23*L[71] + x25*x92 + x25*L[105] + x26*L[148] + x28*x74 + x28*L[15] + x30*x76 + x30*L[29] + x32*x78 + x32*L[49] + x34*x80 + x34*L[76] + x36*x82 + x36*L[111] + x37*L[155] + x4 + x52*x6 + x53*x6 + x54*x8 + x55*x8 + x6*x97 + x6*L[10] + x8*x99 + x8*L[20] + L[1];
#pragma omp atomic
Ls[2] += x*x216 + x*x217 + x*x3 + x*L[5] + x10*x225 + x10*x243 + x10*x98 + x10*L[36] + x100*x12 + x102*x14 + x105*x170 + x105*L[41] + x107*x171 + x107*L[66] + x109*x172 + x109*L[99] + x111*L[141] + x113*x245 + x113*L[43] + x114*x246 + x114*L[69] + x115*x247 + x115*L[103] + x116*L[146] + x118*x177 + x118*L[62] + x12*x226 + x12*x244 + x12*L[57] + x120*x178 + x120*L[94] + x122*L[135] + x124*x248 + x124*L[64] + x125*x249 + x125*L[97] + x126*L[139] + x128*x182 + x128*L[90] + x130*L[130] + x132*x250 + x132*L[92] + x133*L[133] + x135*L[126] + x137*L[128] + x14*x227 + x14*L[85] + x140*x159 + x140*L[52] + x142*x160 + x142*L[80] + x144*x161 + x144*L[116] + x145*L[161] + x147*x162 + x147*L[79] + x149*x163 + x149*L[115] + x15*L[121] + x150*L[160] + x152*x164 + x152*L[114] + x153*L[159] + x154*L[158] + x17*x39 + x17*x69 + x17*x87 + x17*L[16] + x185*L[101] + x187*L[144] + x189*L[143] + x19*x40 + x19*x70 + x19*x89 + x19*L[30] + x190*L[137] + x2 + x21*x41 + x21*x71 + x21*x91 + x21*L[50] + x218*x28 + x219*x30 + x220*x32 + x221*x34 + x222*x36 + x223*x6 + x224*x8 + x228*x28 + x229*x30 + x23*x42 + x23*x72 + x23*x93 + x23*L[77] + x230*x32 + x231*x34 + x232*x36 + x234*x28 + x236*x30 + x238*x32 + x240*x34 + x241*x6 + x242*x8 + x25*x43 + x25*x73 + x25*L[112] + x26*L[156] + x28*L[18] + x30*L[33] + x32*L[54] + x34*L[82] + x36*L[118] + x37*L[163] + x6*x94 + x6*L[11] + x68*y + x8*x96 + x8*L[21] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*x251 + x*x252 + x*x255 + x*L[6] + x10*x265 + x10*x266 + x10*x289 + x10*L[37] + x105*x292 + x105*L[42] + x107*x293 + x107*L[67] + x109*x294 + x109*L[100] + x111*L[142] + x113*x166 + x113*L[44] + x114*x167 + x114*L[70] + x115*x168 + x115*L[104] + x116*L[147] + x118*x295 + x118*L[63] + x12*x267 + x12*x268 + x12*x291 + x12*L[58] + x120*x296 + x120*L[95] + x122*L[136] + x124*x174 + x124*L[65] + x125*x175 + x125*L[98] + x126*L[140] + x128*x297 + x128*L[91] + x130*L[131] + x132*x180 + x132*L[93] + x133*L[134] + x135*L[127] + x137*L[129] + x14*x269 + x14*x270 + x14*L[86] + x140*x156 + x140*L[53] + x142*x157 + x142*L[81] + x144*x158 + x144*L[117] + x145*L[162] + x147*x160 + x147*L[80] + x149*x161 + x149*L[116] + x15*L[122] + x150*L[161] + x152*x163 + x152*L[115] + x153*L[160] + x154*L[159] + x17*x256 + x17*x271 + x17*x277 + x17*L[17] + x185*L[102] + x187*L[145] + x189*L[144] + x19*x257 + x19*x272 + x19*x279 + x19*L[31] + x190*L[138] + x21*x258 + x21*x273 + x21*x281 + x21*L[51] + x23*x259 + x23*x274 + x23*x283 + x23*L[78] + x25*x260 + x25*x275 + x25*L[113] + x253*y + x26*L[157] + x261*x6 + x262*x6 + x263*x8 + x264*x8 + x28*x45 + x28*x63 + x28*x77 + x28*L[19] + x285*x6 + x287*x8 + x30*x46 + x30*x64 + x30*x79 + x30*L[34] + x32*x47 + x32*x65 + x32*x81 + x32*L[55] + x34*x48 + x34*x66 + x34*x83 + x34*L[83] + x36*x49 + x36*x67 + x36*L[119] + x37*L[164] + x6*L[12] + x8*L[22] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += x*x52 + x*x53 + x*x97 + x*L[10] + x10*x103 + x10*x58 + x10*x59 + x10*L[56] + x101*x8 + x105*x181 + x105*L[59] + x107*x182 + x107*L[90] + x109*L[130] + x113*x179 + x113*L[61] + x114*x180 + x114*L[93] + x115*L[134] + x118*x184 + x118*L[87] + x12*x60 + x12*x61 + x12*L[84] + x120*L[126] + x124*x183 + x124*L[89] + x125*L[129] + x128*L[123] + x132*L[125] + x14*L[120] + x140*x312 + x140*L[68] + x142*x313 + x142*L[102] + x144*L[145] + x147*x314 + x147*L[101] + x149*L[144] + x152*L[143] + x165*x28 + x166*x30 + x167*x32 + x168*x34 + x169*x17 + x17*x298 + x17*x309 + x17*L[23] + x170*x19 + x171*x21 + x172*x23 + x185*L[132] + x19*x299 + x19*x310 + x19*L[41] + x21*x300 + x21*x311 + x21*L[66] + x23*x301 + x23*L[99] + x25*L[141] + x28*x302 + x28*x306 + x28*L[25] + x30*x303 + x30*x307 + x30*L[44] + x304*x32 + x305*x34 + x308*x32 + x32*L[70] + x34*L[104] + x36*L[147] + x50 + x51 + x54*x6 + x55*x6 + x56*x8 + x57*x8 + x6*x99 + x6*L[20] + x8*L[35] + x95 + L[4];
#pragma omp atomic
Ls[5] += x*x223 + x*x241 + x*x94 + x*L[11] + x10*x100 + x10*x226 + x10*x244 + x10*L[57] + x102*x12 + x105*x177 + x105*L[62] + x107*x178 + x107*L[94] + x109*L[135] + x113*x248 + x113*L[64] + x114*x249 + x114*L[97] + x115*L[139] + x118*x182 + x118*L[90] + x12*x227 + x12*L[85] + x120*L[130] + x124*x250 + x124*L[92] + x125*L[133] + x128*L[126] + x132*L[128] + x14*L[121] + x140*x213 + x140*L[73] + x142*x214 + x142*L[108] + x144*L[152] + x147*x215 + x147*L[107] + x149*L[151] + x152*L[150] + x17*x193 + x17*x207 + x17*x86 + x17*L[26] + x185*L[137] + x19*x194 + x19*x208 + x19*x88 + x19*L[45] + x195*x21 + x196*x23 + x209*x21 + x21*x90 + x21*L[71] + x216 + x217 + x224*x6 + x225*x8 + x23*x92 + x23*L[105] + x233*x28 + x235*x30 + x237*x32 + x239*x34 + x242*x6 + x243*x8 + x25*L[148] + x28*x315 + x28*x319 + x28*L[28] + x3 + x30*x316 + x30*x320 + x30*L[48] + x317*x32 + x318*x34 + x32*x321 + x32*L[75] + x34*L[110] + x36*L[154] + x6*x96 + x6*L[21] + x8*x98 + x8*L[36] + L[5];
#pragma omp atomic
Ls[6] += x*x261 + x*x262 + x*x285 + x*L[12] + x10*x267 + x10*x268 + x10*x291 + x10*L[58] + x105*x295 + x105*L[63] + x107*x296 + x107*L[95] + x109*L[136] + x113*x174 + x113*L[65] + x114*x175 + x114*L[98] + x115*L[140] + x118*x297 + x118*L[91] + x12*x269 + x12*x270 + x12*L[86] + x120*L[131] + x124*x180 + x124*L[93] + x125*L[134] + x128*L[127] + x132*L[129] + x14*L[122] + x140*x211 + x140*L[74] + x142*x212 + x142*L[109] + x144*L[153] + x147*x214 + x147*L[108] + x149*L[152] + x152*L[151] + x17*x276 + x17*x322 + x17*x326 + x17*L[27] + x185*L[138] + x19*x278 + x19*x323 + x19*x327 + x19*L[46] + x198*x28 + x199*x30 + x200*x32 + x201*x34 + x203*x28 + x204*x30 + x205*x32 + x21*x280 + x21*x324 + x21*x328 + x21*L[72] + x23*x282 + x23*x325 + x23*L[106] + x25*L[149] + x251 + x252 + x255 + x263*x6 + x264*x6 + x265*x8 + x266*x8 + x28*x76 + x28*L[29] + x287*x6 + x289*x8 + x30*x78 + x30*L[49] + x32*x80 + x32*L[76] + x34*x82 + x34*L[111] + x36*L[155] + x6*L[22] + x8*L[37] + L[6];
#pragma omp atomic
Ls[7] += x*x329 + x*x330 + x10*x181 + x10*x333 + x10*x347 + x10*L[59] + x105*x171 + x105*L[66] + x107*x172 + x107*L[99] + x109*L[141] + x113*x348 + x113*L[68] + x114*x349 + x114*L[102] + x115*L[145] + x118*x178 + x118*L[94] + x12*x184 + x12*x334 + x12*L[87] + x120*L[135] + x124*x350 + x124*L[96] + x125*L[138] + x128*L[130] + x132*L[132] + x14*L[123] + x140*x162 + x140*L[79] + x142*x163 + x142*L[115] + x144*L[160] + x147*x164 + x147*L[114] + x149*L[159] + x152*L[158] + x155*x28 + x156*x30 + x157*x32 + x158*x34 + x169*x6 + x17*x40 + x17*x70 + x17*x89 + x17*L[30] + x176*x8 + x185*L[143] + x19*x41 + x19*x71 + x19*x91 + x19*L[50] + x21*x42 + x21*x72 + x21*x93 + x21*L[77] + x23*x43 + x23*x73 + x23*L[112] + x25*L[156] + x28*x335 + x28*x340 + x28*L[32] + x30*x336 + x30*x342 + x30*L[53] + x32*x337 + x32*x344 + x32*L[81] + x331*x6 + x332*x8 + x338*x34 + x34*L[117] + x345*x6 + x346*x8 + x36*L[162] + x38 + x6*L[23] + x68 + x69*y + x8*L[38] + x85 + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += x*x254 + x*x351 + x*x352 + x*L[14] + x10*x288 + x10*x355 + x10*x359 + x10*L[60] + x105*x293 + x105*L[67] + x107*x294 + x107*L[100] + x109*L[142] + x113*x246 + x113*L[69] + x114*x247 + x114*L[103] + x115*L[146] + x118*x296 + x118*L[95] + x12*x290 + x12*x356 + x12*L[88] + x120*L[136] + x124*x249 + x124*L[97] + x125*L[139] + x128*L[131] + x132*L[133] + x14*L[124] + x140*x160 + x140*L[80] + x142*x161 + x142*L[116] + x144*L[161] + x147*x163 + x147*L[115] + x149*L[160] + x152*L[159] + x17*x257 + x17*x272 + x17*x279 + x17*L[31] + x185*L[144] + x19*x258 + x19*x273 + x19*x281 + x19*L[51] + x21*x259 + x21*x274 + x21*x283 + x21*L[78] + x219*x28 + x220*x30 + x221*x32 + x222*x34 + x229*x28 + x23*x260 + x23*x275 + x23*L[113] + x230*x30 + x231*x32 + x232*x34 + x236*x28 + x238*x30 + x240*x32 + x25*L[157] + x253 + x271*y + x28*L[33] + x284*x6 + x286*x8 + x30*L[54] + x32*L[82] + x34*L[118] + x353*x6 + x354*x8 + x357*x6 + x358*x8 + x36*L[163] + x6*L[24] + x8*L[39] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += x*x360 + x*x363 + x10*x179 + x10*x366 + x10*x383 + x10*L[61] + x105*x384 + x105*L[68] + x107*x385 + x107*L[101] + x109*L[143] + x113*x167 + x113*L[70] + x114*x168 + x114*L[104] + x115*L[147] + x118*x386 + x118*L[96] + x12*x183 + x12*x367 + x12*L[89] + x120*L[137] + x124*x175 + x124*L[98] + x125*L[140] + x128*L[132] + x132*L[134] + x14*L[125] + x140*x157 + x140*L[81] + x142*x158 + x142*L[117] + x144*L[162] + x147*x161 + x147*L[116] + x149*L[161] + x152*L[160] + x155*x17 + x159*x19 + x162*x21 + x164*x23 + x165*x6 + x17*x368 + x17*x373 + x17*L[32] + x173*x8 + x185*L[145] + x19*x369 + x19*x375 + x19*L[52] + x21*x370 + x21*x377 + x21*L[79] + x23*x371 + x23*L[114] + x25*L[158] + x28*x46 + x28*x64 + x28*x79 + x28*L[34] + x30*x47 + x30*x65 + x30*x81 + x30*L[55] + x32*x48 + x32*x66 + x32*x83 + x32*L[83] + x34*x49 + x34*x67 + x34*L[119] + x36*L[164] + x361*y + x364*x6 + x365*x8 + x379*x6 + x381*x8 + x44 + x6*L[25] + x62 + x75 + x8*L[40] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += x*x54 + x*x55 + x*x99 + x*L[20] + x10*x60 + x10*x61 + x10*L[84] + x101*x6 + x103*x8 + x105*x184 + x105*L[87] + x107*L[126] + x113*x183 + x113*L[89] + x114*L[129] + x118*L[123] + x12*L[120] + x124*L[125] + x140*x397 + x140*L[96] + x142*L[138] + x147*L[137] + x17*x176 + x17*x387 + x17*x395 + x17*L[38] + x173*x28 + x174*x30 + x175*x32 + x177*x19 + x178*x21 + x19*x388 + x19*x396 + x19*L[62] + x21*x389 + x21*L[94] + x23*L[135] + x28*x390 + x28*x393 + x28*L[40] + x30*x391 + x30*x394 + x30*L[65] + x32*x392 + x32*L[98] + x34*L[140] + x52 + x53 + x56*x6 + x57*x6 + x58*x8 + x59*x8 + x6*L[35] + x8*L[56] + x97 + L[10];
#pragma omp atomic
Ls[11] += x*x224 + x*x242 + x*x96 + x*L[21] + x10*x102 + x10*x227 + x10*L[85] + x100*x8 + x105*x182 + x105*L[90] + x107*L[130] + x113*x250 + x113*L[92] + x114*L[133] + x118*L[126] + x12*L[121] + x124*L[128] + x140*x314 + x140*L[101] + x142*L[144] + x147*L[143] + x17*x170 + x17*x299 + x17*x310 + x17*L[41] + x171*x19 + x172*x21 + x19*x300 + x19*x311 + x19*L[66] + x21*x301 + x21*L[99] + x223 + x225*x6 + x226*x8 + x23*L[141] + x241 + x243*x6 + x244*x8 + x245*x28 + x246*x30 + x247*x32 + x28*x398 + x28*x401 + x28*L[43] + x30*x399 + x30*x402 + x30*L[69] + x32*x400 + x32*L[103] + x34*L[146] + x6*x98 + x6*L[36] + x8*L[57] + x94 + L[11];
#pragma omp atomic
Ls[12] += x*x263 + x*x264 + x*x287 + x*L[22] + x10*x269 + x10*x270 + x10*L[86] + x105*x297 + x105*L[91] + x107*L[131] + x113*x180 + x113*L[93] + x114*L[134] + x118*L[127] + x12*L[122] + x124*L[129] + x140*x313 + x140*L[102] + x142*L[145] + x147*L[144] + x166*x28 + x167*x30 + x168*x32 + x17*x292 + x17*x403 + x17*x406 + x17*L[42] + x19*x293 + x19*x404 + x19*x407 + x19*L[67] + x21*x294 + x21*x405 + x21*L[100] + x23*L[142] + x261 + x262 + x265*x6 + x266*x6 + x267*x8 + x268*x8 + x28*x303 + x28*x307 + x28*L[44] + x285 + x289*x6 + x291*x8 + x30*x304 + x30*x308 + x30*L[70] + x305*x32 + x32*L[104] + x34*L[147] + x6*L[37] + x8*L[58] + L[12];
#pragma omp atomic
Ls[13] += x*x331 + x*x345 + x10*x184 + x10*x334 + x10*L[87] + x105*x178 + x105*L[94] + x107*L[135] + x113*x350 + x113*L[96] + x114*L[138] + x118*L[130] + x12*L[123] + x124*L[132] + x140*x215 + x140*L[107] + x142*L[151] + x147*L[150] + x17*x194 + x17*x208 + x17*x88 + x17*L[45] + x176*x6 + x181*x8 + x19*x195 + x19*x209 + x19*x90 + x19*L[71] + x192 + x196*x21 + x206 + x21*x92 + x21*L[105] + x210*x28 + x211*x30 + x212*x32 + x23*L[148] + x28*x339 + x28*x408 + x28*L[47] + x30*x341 + x30*x409 + x30*L[74] + x32*x343 + x32*L[109] + x329 + x330 + x332*x6 + x333*x8 + x34*L[153] + x346*x6 + x347*x8 + x6*L[38] + x8*L[59] + x84 + L[13];
#pragma omp atomic
Ls[14] += x*x284 + x*x353 + x*x357 + x*L[24] + x10*x290 + x10*x356 + x10*L[88] + x105*x296 + x105*L[95] + x107*L[136] + x113*x249 + x113*L[97] + x114*L[139] + x118*L[131] + x12*L[124] + x124*L[133] + x140*x214 + x140*L[108] + x142*L[152] + x147*L[151] + x17*x278 + x17*x323 + x17*x327 + x17*L[46] + x19*x280 + x19*x324 + x19*x328 + x19*L[72] + x21*x282 + x21*x325 + x21*L[106] + x23*L[149] + x235*x28 + x237*x30 + x239*x32 + x254 + x28*x316 + x28*x320 + x28*L[48] + x286*x6 + x288*x8 + x30*x317 + x30*x321 + x30*L[75] + x318*x32 + x32*L[110] + x34*L[154] + x351 + x352 + x354*x6 + x355*x8 + x358*x6 + x359*x8 + x6*L[39] + x8*L[60] + L[14];
#pragma omp atomic
Ls[15] += x*x364 + x*x379 + x10*x183 + x10*x367 + x10*L[89] + x105*x386 + x105*L[96] + x107*L[137] + x113*x175 + x113*L[98] + x114*L[140] + x118*L[132] + x12*L[125] + x124*L[134] + x140*x212 + x140*L[109] + x142*L[153] + x147*L[152] + x17*x210 + x17*x372 + x17*x410 + x17*L[47] + x173*x6 + x179*x8 + x19*x213 + x19*x374 + x19*x411 + x19*L[73] + x197 + x199*x28 + x200*x30 + x201*x32 + x202 + x204*x28 + x205*x30 + x21*x215 + x21*x376 + x21*L[107] + x23*L[150] + x28*x78 + x28*L[49] + x30*x80 + x30*L[76] + x32*x82 + x32*L[111] + x34*L[155] + x360 + x363 + x365*x6 + x366*x8 + x381*x6 + x383*x8 + x6*L[40] + x74 + x8*L[61] + L[15];
#pragma omp atomic
Ls[16] += x*x412 + x*x413 + x10*x182 + x10*x416 + x10*L[90] + x105*x172 + x105*L[99] + x107*L[141] + x113*x426 + x113*L[101] + x114*L[144] + x118*L[135] + x12*L[126] + x124*L[137] + x140*x164 + x140*L[114] + x142*L[159] + x147*L[158] + x159*x28 + x160*x30 + x161*x32 + x17*x41 + x17*x71 + x17*x91 + x17*L[50] + x170*x6 + x177*x8 + x19*x42 + x19*x72 + x19*x93 + x19*L[77] + x21*x43 + x21*x73 + x21*L[112] + x23*L[156] + x28*x417 + x28*x421 + x28*L[52] + x30*x418 + x30*x423 + x30*L[80] + x32*x419 + x32*L[116] + x34*L[161] + x39 + x414*x6 + x415*x8 + x424*x6 + x425*x8 + x6*L[41] + x69 + x70*y + x8*L[62] + x87 + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += x*x427 + x*x428 + x10*x297 + x10*x431 + x10*L[91] + x105*x294 + x105*L[100] + x107*L[142] + x113*x349 + x113*L[102] + x114*L[145] + x118*L[136] + x12*L[127] + x124*L[138] + x140*x163 + x140*L[115] + x142*L[160] + x147*L[159] + x156*x28 + x157*x30 + x158*x32 + x17*x258 + x17*x273 + x17*x281 + x17*L[51] + x19*x259 + x19*x274 + x19*x283 + x19*L[78] + x21*x260 + x21*x275 + x21*L[113] + x23*L[157] + x256 + x271 + x272*y + x277 + x28*x336 + x28*x342 + x28*L[53] + x292*x6 + x295*x8 + x30*x337 + x30*x344 + x30*L[81] + x32*x338 + x32*L[117] + x34*L[162] + x429*x6 + x430*x8 + x432*x6 + x433*x8 + x6*L[42] + x8*L[63] + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += x*x362 + x*x434 + x10*x250 + x10*x382 + x10*L[92] + x105*x385 + x105*L[101] + x107*L[143] + x113*x247 + x113*L[103] + x114*L[146] + x118*L[137] + x12*L[128] + x124*L[139] + x140*x161 + x140*L[116] + x142*L[161] + x147*L[160] + x159*x17 + x162*x19 + x164*x21 + x17*x369 + x17*x375 + x17*L[52] + x19*x370 + x19*x377 + x19*L[79] + x21*x371 + x21*L[114] + x218 + x220*x28 + x221*x30 + x222*x32 + x228 + x23*L[158] + x230*x28 + x231*x30 + x232*x32 + x234 + x238*x28 + x240*x30 + x245*x6 + x248*x8 + x28*L[54] + x30*L[82] + x32*L[118] + x34*L[163] + x361 + x368*y + x378*x6 + x380*x8 + x435*x6 + x436*x8 + x6*L[43] + x8*L[64] + L[18];
#pragma omp atomic
Ls[19] += x*x437 + x*x440 + x10*x180 + x10*x443 + x10*L[93] + x105*x455 + x105*L[102] + x107*L[144] + x113*x168 + x113*L[104] + x114*L[147] + x118*L[138] + x12*L[129] + x124*L[140] + x140*x158 + x140*L[117] + x142*L[162] + x147*L[161] + x156*x17 + x160*x19 + x163*x21 + x166*x6 + x17*x444 + x17*x448 + x17*L[53] + x174*x8 + x19*x445 + x19*x450 + x19*L[80] + x21*x446 + x21*L[115] + x23*L[159] + x28*x47 + x28*x65 + x28*x81 + x28*L[55] + x30*x48 + x30*x66 + x30*x83 + x30*L[83] + x32*x49 + x32*x67 + x32*L[119] + x34*L[164] + x438*y + x441*x6 + x442*x8 + x45 + x452*x6 + x454*x8 + x6*L[44] + x63 + x77 + x8*L[65] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += x*x101 + x*x56 + x*x57 + x*L[35] + x10*L[120] + x103*x6 + x105*L[123] + x113*L[125] + x140*L[132] + x17*x181 + x17*x456 + x17*x461 + x17*L[59] + x179*x28 + x180*x30 + x182*x19 + x19*x457 + x19*L[90] + x21*L[130] + x28*x458 + x28*x460 + x28*L[61] + x30*x459 + x30*L[93] + x32*L[134] + x54 + x55 + x58*x6 + x59*x6 + x6*L[56] + x60*x8 + x61*x8 + x8*L[84] + x99 + L[20];
#pragma omp atomic
Ls[21] += x*x225 + x*x243 + x*x98 + x*L[36] + x10*L[121] + x100*x6 + x102*x8 + x105*L[126] + x113*L[128] + x140*L[137] + x17*x177 + x17*x388 + x17*x396 + x17*L[62] + x178*x19 + x19*x389 + x19*L[94] + x21*L[135] + x224 + x226*x6 + x227*x8 + x242 + x244*x6 + x248*x28 + x249*x30 + x28*x462 + x28*x464 + x28*L[64] + x30*x463 + x30*L[97] + x32*L[139] + x6*L[57] + x8*L[85] + x96 + L[21];
#pragma omp atomic
Ls[22] += x*x265 + x*x266 + x*x289 + x*L[37] + x10*L[122] + x105*L[127] + x113*L[129] + x140*L[138] + x17*x295 + x17*x465 + x17*x467 + x17*L[63] + x174*x28 + x175*x30 + x19*x296 + x19*x466 + x19*L[95] + x21*L[136] + x263 + x264 + x267*x6 + x268*x6 + x269*x8 + x270*x8 + x28*x391 + x28*x394 + x28*L[65] + x287 + x291*x6 + x30*x392 + x30*L[98] + x32*L[140] + x6*L[58] + x8*L[86] + L[22];
#pragma omp atomic
Ls[23] += x*x332 + x*x346 + x10*L[123] + x105*L[130] + x113*L[132] + x140*L[143] + x169 + x17*x171 + x17*x300 + x17*x311 + x17*L[66] + x172*x19 + x181*x6 + x184*x8 + x19*x301 + x19*L[99] + x21*L[141] + x28*x312 + x28*x348 + x28*x468 + x28*L[68] + x298 + x30*x313 + x30*x349 + x30*L[102] + x309 + x32*L[145] + x331 + x333*x6 + x334*x8 + x345 + x347*x6 + x6*L[59] + x8*L[87] + L[23];
#pragma omp atomic
Ls[24] += x*x286 + x*x354 + x*x358 + x*L[39] + x10*L[124] + x105*L[131] + x113*L[133] + x140*L[144] + x17*x293 + x17*x404 + x17*x407 + x17*L[67] + x19*x294 + x19*x405 + x19*L[100] + x21*L[142] + x246*x28 + x247*x30 + x28*x399 + x28*x402 + x28*L[69] + x284 + x288*x6 + x290*x8 + x30*x400 + x30*L[103] + x32*L[146] + x353 + x355*x6 + x356*x8 + x357 + x359*x6 + x6*L[60] + x8*L[88] + L[24];
#pragma omp atomic
Ls[25] += x*x365 + x*x381 + x10*L[125] + x105*L[132] + x113*L[134] + x140*L[145] + x165 + x167*x28 + x168*x30 + x17*x312 + x17*x384 + x17*x469 + x17*L[68] + x179*x6 + x183*x8 + x19*x314 + x19*x385 + x19*L[101] + x21*L[143] + x28*x304 + x28*x308 + x28*L[70] + x30*x305 + x30*L[104] + x302 + x306 + x32*L[147] + x364 + x366*x6 + x367*x8 + x379 + x383*x6 + x6*L[61] + x8*L[89] + L[25];
#pragma omp atomic
Ls[26] += x*x414 + x*x424 + x10*L[126] + x105*L[135] + x113*L[137] + x140*L[150] + x17*x195 + x17*x209 + x17*x90 + x17*L[71] + x177*x6 + x182*x8 + x19*x196 + x19*x92 + x19*L[105] + x193 + x207 + x21*L[148] + x213*x28 + x214*x30 + x28*x420 + x28*x470 + x28*L[73] + x30*x422 + x30*L[108] + x32*L[152] + x412 + x413 + x415*x6 + x416*x8 + x425*x6 + x6*L[62] + x8*L[90] + x86 + L[26];
#pragma omp atomic
Ls[27] += x*x429 + x*x432 + x10*L[127] + x105*L[136] + x113*L[138] + x140*L[151] + x17*x280 + x17*x324 + x17*x328 + x17*L[72] + x19*x282 + x19*x325 + x19*L[106] + x21*L[149] + x211*x28 + x212*x30 + x276 + x28*x341 + x28*x409 + x28*L[74] + x295*x6 + x297*x8 + x30*x343 + x30*L[109] + x32*L[153] + x322 + x326 + x427 + x428 + x430*x6 + x431*x8 + x433*x6 + x6*L[63] + x8*L[91] + L[27];
#pragma omp atomic
Ls[28] += x*x378 + x*x435 + x10*L[128] + x105*L[137] + x113*L[139] + x140*L[152] + x17*x213 + x17*x374 + x17*x411 + x17*L[73] + x19*x215 + x19*x376 + x19*L[107] + x21*L[150] + x233 + x237*x28 + x239*x30 + x248*x6 + x250*x8 + x28*x317 + x28*x321 + x28*L[75] + x30*x318 + x30*L[110] + x315 + x319 + x32*L[154] + x362 + x380*x6 + x382*x8 + x434 + x436*x6 + x6*L[64] + x8*L[92] + L[28];
#pragma omp atomic
Ls[29] += x*x441 + x*x452 + x10*L[129] + x105*L[138] + x113*L[140] + x140*L[153] + x17*x211 + x17*x447 + x17*x471 + x17*L[74] + x174*x6 + x180*x8 + x19*x214 + x19*x449 + x19*L[108] + x198 + x200*x28 + x201*x30 + x203 + x205*x28 + x21*L[151] + x28*x80 + x28*L[76] + x30*x82 + x30*L[111] + x32*L[155] + x437 + x440 + x442*x6 + x443*x8 + x454*x6 + x6*L[65] + x76 + x8*L[93] + L[29];
#pragma omp atomic
Ls[30] += x*x472 + x*x473 + x10*L[130] + x105*L[141] + x113*L[143] + x140*L[158] + x162*x28 + x163*x30 + x17*x42 + x17*x72 + x17*x93 + x17*L[77] + x171*x6 + x178*x8 + x19*x43 + x19*x73 + x19*L[112] + x21*L[156] + x28*x476 + x28*x479 + x28*L[79] + x30*x477 + x30*L[115] + x32*L[160] + x40 + x474*x6 + x475*x8 + x480*x6 + x6*L[66] + x70 + x71*y + x8*L[94] + x89 + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += x*x481 + x*x482 + x10*L[131] + x105*L[142] + x113*L[144] + x140*L[159] + x160*x28 + x161*x30 + x17*x259 + x17*x274 + x17*x283 + x17*L[78] + x19*x260 + x19*x275 + x19*L[113] + x21*L[157] + x257 + x272 + x273*y + x279 + x28*x418 + x28*x423 + x28*L[80] + x293*x6 + x296*x8 + x30*x419 + x30*L[116] + x32*L[161] + x483*x6 + x484*x8 + x485*x6 + x6*L[67] + x8*L[95] + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += x*x486 + x10*L[132] + x105*L[143] + x113*L[145] + x140*L[160] + x155 + x157*x28 + x158*x30 + x162*x17 + x164*x19 + x17*x370 + x17*x377 + x17*L[79] + x19*x371 + x19*L[114] + x21*L[158] + x28*x337 + x28*x344 + x28*L[81] + x30*x338 + x30*L[117] + x32*L[162] + x335 + x340 + x348*x6 + x350*x8 + x368 + x369*y + x373 + x384*x6 + x386*x8 + x487*x6 + x6*L[68] + x8*L[96] + L[32];
#pragma omp atomic
Ls[33] += x*x439 + x*x488 + x10*L[133] + x105*L[144] + x113*L[146] + x140*L[161] + x160*x17 + x163*x19 + x17*x445 + x17*x450 + x17*L[80] + x19*x446 + x19*L[115] + x21*L[159] + x219 + x221*x28 + x222*x30 + x229 + x231*x28 + x232*x30 + x236 + x240*x28 + x246*x6 + x249*x8 + x28*L[82] + x30*L[118] + x32*L[163] + x438 + x444*y + x451*x6 + x453*x8 + x489*x6 + x6*L[69] + x8*L[97] + L[33];
#pragma omp atomic
Ls[34] += x*x490 + x*x493 + x10*L[134] + x105*L[145] + x113*L[147] + x140*L[162] + x157*x17 + x161*x19 + x167*x6 + x17*x496 + x17*x499 + x17*L[81] + x175*x8 + x19*x497 + x19*L[116] + x21*L[160] + x28*x48 + x28*x66 + x28*x83 + x28*L[83] + x30*x49 + x30*x67 + x30*L[119] + x32*L[164] + x46 + x491*y + x494*x6 + x495*x8 + x501*x6 + x6*L[70] + x64 + x79 + x8*L[98] + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += x*x103 + x*x58 + x*x59 + x*L[56] + x101 + x17*x184 + x17*x502 + x17*L[87] + x183*x28 + x19*L[126] + x28*x503 + x28*L[89] + x30*L[129] + x56 + x57 + x6*x60 + x6*x61 + x6*L[84] + x8*L[120] + L[35];
#pragma omp atomic
Ls[36] += x*x100 + x*x226 + x*x244 + x*L[57] + x102*x6 + x17*x182 + x17*x457 + x17*L[90] + x19*L[130] + x225 + x227*x6 + x243 + x250*x28 + x28*x504 + x28*L[92] + x30*L[133] + x6*L[85] + x8*L[121] + x98 + L[36];
#pragma omp atomic
Ls[37] += x*x267 + x*x268 + x*x291 + x*L[58] + x17*x297 + x17*x505 + x17*L[91] + x180*x28 + x19*L[131] + x265 + x266 + x269*x6 + x270*x6 + x28*x459 + x28*L[93] + x289 + x30*L[134] + x6*L[86] + x8*L[122] + L[37];
#pragma omp atomic
Ls[38] += x*x333 + x*x347 + x17*x178 + x17*x389 + x17*L[94] + x176 + x184*x6 + x19*L[135] + x28*x350 + x28*x397 + x28*L[96] + x30*L[138] + x332 + x334*x6 + x346 + x387 + x395 + x6*L[87] + x8*L[123] + L[38];
#pragma omp atomic
Ls[39] += x*x288 + x*x355 + x*x359 + x*L[60] + x17*x296 + x17*x466 + x17*L[95] + x19*L[136] + x249*x28 + x28*x463 + x28*L[97] + x286 + x290*x6 + x30*L[139] + x354 + x356*x6 + x358 + x6*L[88] + x8*L[124] + L[39];
#pragma omp atomic
Ls[40] += x*x366 + x*x383 + x17*x386 + x17*x397 + x17*L[96] + x173 + x175*x28 + x183*x6 + x19*L[137] + x28*x392 + x28*L[98] + x30*L[140] + x365 + x367*x6 + x381 + x390 + x393 + x6*L[89] + x8*L[125] + L[40];
#pragma omp atomic
Ls[41] += x*x415 + x*x425 + x17*x172 + x17*x301 + x17*L[99] + x170 + x182*x6 + x19*L[141] + x28*x314 + x28*x426 + x28*L[101] + x299 + x30*L[144] + x310 + x414 + x416*x6 + x424 + x6*L[90] + x8*L[126] + L[41];
#pragma omp atomic
Ls[42] += x*x430 + x*x433 + x17*x294 + x17*x405 + x17*L[100] + x19*L[142] + x28*x313 + x28*x349 + x28*L[102] + x292 + x297*x6 + x30*L[145] + x403 + x406 + x429 + x431*x6 + x432 + x6*L[91] + x8*L[127] + L[42];
#pragma omp atomic
Ls[43] += x*x380 + x*x436 + x17*x314 + x17*x385 + x17*L[101] + x19*L[143] + x245 + x247*x28 + x250*x6 + x28*x400 + x28*L[103] + x30*L[146] + x378 + x382*x6 + x398 + x401 + x435 + x6*L[92] + x8*L[128] + L[43];
#pragma omp atomic
Ls[44] += x*x442 + x*x454 + x166 + x168*x28 + x17*x313 + x17*x455 + x17*L[102] + x180*x6 + x19*L[144] + x28*x305 + x28*L[104] + x30*L[147] + x303 + x307 + x441 + x443*x6 + x452 + x6*L[93] + x8*L[129] + L[44];
#pragma omp atomic
Ls[45] += x*x474 + x*x480 + x17*x196 + x17*x92 + x17*L[105] + x178*x6 + x19*L[148] + x194 + x208 + x215*x28 + x28*x478 + x28*L[107] + x30*L[151] + x472 + x473 + x475*x6 + x6*L[94] + x8*L[130] + x88 + L[45];
#pragma omp atomic
Ls[46] += x*x483 + x*x485 + x17*x282 + x17*x325 + x17*L[106] + x19*L[149] + x214*x28 + x278 + x28*x422 + x28*L[108] + x296*x6 + x30*L[152] + x323 + x327 + x481 + x482 + x484*x6 + x6*L[95] + x8*L[131] + L[46];
#pragma omp atomic
Ls[47] += x*x487 + x17*x215 + x17*x376 + x17*L[107] + x19*L[150] + x210 + x212*x28 + x28*x343 + x28*L[109] + x30*L[153] + x339 + x350*x6 + x372 + x386*x6 + x408 + x410 + x486 + x6*L[96] + x8*L[132] + L[47];
#pragma omp atomic
Ls[48] += x*x451 + x*x489 + x17*x214 + x17*x449 + x17*L[108] + x19*L[151] + x235 + x239*x28 + x249*x6 + x28*x318 + x28*L[110] + x30*L[154] + x316 + x320 + x439 + x453*x6 + x488 + x6*L[97] + x8*L[133] + L[48];
#pragma omp atomic
Ls[49] += x*x494 + x*x501 + x17*x212 + x17*x498 + x17*L[109] + x175*x6 + x19*L[152] + x199 + x201*x28 + x204 + x28*x82 + x28*L[111] + x30*L[155] + x490 + x493 + x495*x6 + x6*L[98] + x78 + x8*L[134] + L[49];
#pragma omp atomic
Ls[50] += x*x506 + x*x507 + x164*x28 + x17*x43 + x17*x73 + x17*L[112] + x172*x6 + x19*L[156] + x28*x509 + x28*L[114] + x30*L[159] + x41 + x508*x6 + x6*L[99] + x71 + x72*y + x8*L[135] + x91 + y*L[77] + L[50];
#pragma omp atomic
Ls[51] += x*x510 + x*x511 + x163*x28 + x17*x260 + x17*x275 + x17*L[113] + x19*L[157] + x258 + x273 + x274*y + x28*x477 + x28*L[115] + x281 + x294*x6 + x30*L[160] + x512*x6 + x6*L[100] + x8*L[136] + y*L[78] + L[51];
#pragma omp atomic
Ls[52] += x*x513 + x159 + x161*x28 + x164*x17 + x17*x371 + x17*L[114] + x19*L[158] + x28*x419 + x28*L[116] + x30*L[161] + x369 + x370*y + x375 + x385*x6 + x417 + x421 + x426*x6 + x6*L[101] + x8*L[137] + L[52];
#pragma omp atomic
Ls[53] += x*x514 + x156 + x158*x28 + x163*x17 + x17*x446 + x17*L[115] + x19*L[159] + x28*x338 + x28*L[117] + x30*L[162] + x336 + x342 + x349*x6 + x444 + x445*y + x448 + x455*x6 + x6*L[102] + x8*L[138] + L[53];
#pragma omp atomic
Ls[54] += x*x492 + x*x515 + x161*x17 + x17*x497 + x17*L[116] + x19*L[160] + x220 + x222*x28 + x230 + x232*x28 + x238 + x247*x6 + x28*L[118] + x30*L[163] + x491 + x496*y + x500*x6 + x6*L[103] + x8*L[139] + L[54];
#pragma omp atomic
Ls[55] += x*x516 + x*x519 + x158*x17 + x168*x6 + x17*x521 + x17*L[117] + x19*L[161] + x28*x49 + x28*x67 + x28*L[119] + x30*L[164] + x47 + x517*y + x520*x6 + x6*L[104] + x65 + x8*L[140] + x81 + z*L[83] + L[55];
#pragma omp atomic
Ls[56] += x*x60 + x*x61 + x*L[84] + x103 + x17*L[123] + x28*L[125] + x58 + x59 + x6*L[120] + L[56];
#pragma omp atomic
Ls[57] += x*x102 + x*x227 + x*L[85] + x100 + x17*L[126] + x226 + x244 + x28*L[128] + x6*L[121] + L[57];
#pragma omp atomic
Ls[58] += x*x269 + x*x270 + x*L[86] + x17*L[127] + x267 + x268 + x28*L[129] + x291 + x6*L[122] + L[58];
#pragma omp atomic
Ls[59] += x*x334 + x17*L[130] + x181 + x28*L[132] + x333 + x347 + x456 + x461 + x6*L[123] + L[59];
#pragma omp atomic
Ls[60] += x*x290 + x*x356 + x*L[88] + x17*L[131] + x28*L[133] + x288 + x355 + x359 + x6*L[124] + L[60];
#pragma omp atomic
Ls[61] += x*x367 + x17*L[132] + x179 + x28*L[134] + x366 + x383 + x458 + x460 + x6*L[125] + L[61];
#pragma omp atomic
Ls[62] += x*x416 + x17*L[135] + x177 + x28*L[137] + x388 + x396 + x415 + x425 + x6*L[126] + L[62];
#pragma omp atomic
Ls[63] += x*x431 + x17*L[136] + x28*L[138] + x295 + x430 + x433 + x465 + x467 + x6*L[127] + L[63];
#pragma omp atomic
Ls[64] += x*x382 + x17*L[137] + x248 + x28*L[139] + x380 + x436 + x462 + x464 + x6*L[128] + L[64];
#pragma omp atomic
Ls[65] += x*x443 + x17*L[138] + x174 + x28*L[140] + x391 + x394 + x442 + x454 + x6*L[129] + L[65];
#pragma omp atomic
Ls[66] += x*x475 + x17*L[141] + x171 + x28*L[143] + x300 + x311 + x474 + x480 + x6*L[130] + L[66];
#pragma omp atomic
Ls[67] += x*x484 + x17*L[142] + x28*L[144] + x293 + x404 + x407 + x483 + x485 + x6*L[131] + L[67];
#pragma omp atomic
Ls[68] += x17*L[143] + x28*L[145] + x312 + x348 + x384 + x468 + x469 + x487 + x6*L[132] + L[68];
#pragma omp atomic
Ls[69] += x*x453 + x17*L[144] + x246 + x28*L[146] + x399 + x402 + x451 + x489 + x6*L[133] + L[69];
#pragma omp atomic
Ls[70] += x*x495 + x167 + x17*L[145] + x28*L[147] + x304 + x308 + x494 + x501 + x6*L[134] + L[70];
#pragma omp atomic
Ls[71] += x*x508 + x17*L[148] + x195 + x209 + x28*L[150] + x506 + x507 + x6*L[135] + x90 + L[71];
#pragma omp atomic
Ls[72] += x*x512 + x17*L[149] + x28*L[151] + x280 + x324 + x328 + x510 + x511 + x6*L[136] + L[72];
#pragma omp atomic
Ls[73] += x17*L[150] + x213 + x28*L[152] + x374 + x411 + x420 + x470 + x513 + x6*L[137] + L[73];
#pragma omp atomic
Ls[74] += x17*L[151] + x211 + x28*L[153] + x341 + x409 + x447 + x471 + x514 + x6*L[138] + L[74];
#pragma omp atomic
Ls[75] += x*x500 + x17*L[152] + x237 + x28*L[154] + x317 + x321 + x492 + x515 + x6*L[139] + L[75];
#pragma omp atomic
Ls[76] += x*x520 + x17*L[153] + x200 + x205 + x28*L[155] + x516 + x519 + x6*L[140] + x80 + L[76];
#pragma omp atomic
Ls[77] += x*x522 + x17*L[156] + x28*L[158] + x42 + x6*L[141] + x72 + x73*y + x93 + y*L[112] + L[77];
#pragma omp atomic
Ls[78] += x*x523 + x17*L[157] + x259 + x274 + x275*y + x28*L[159] + x283 + x6*L[142] + y*L[113] + L[78];
#pragma omp atomic
Ls[79] += x162 + x17*L[158] + x28*L[160] + x370 + x371*y + x377 + x476 + x479 + x6*L[143] + L[79];
#pragma omp atomic
Ls[80] += x160 + x17*L[159] + x28*L[161] + x418 + x423 + x445 + x446*y + x450 + x6*L[144] + L[80];
#pragma omp atomic
Ls[81] += x157 + x17*L[160] + x28*L[162] + x337 + x344 + x496 + x497*y + x499 + x6*L[145] + L[81];
#pragma omp atomic
Ls[82] += x*x518 + x17*L[161] + x221 + x231 + x240 + x28*L[163] + x517 + x521*y + x6*L[146] + L[82];
#pragma omp atomic
Ls[83] += x*x524 + x17*L[162] + x28*L[164] + x48 + x525*y + x6*L[147] + x66 + x83 + z*L[119] + L[83];
#pragma omp atomic
Ls[84] += x*L[120] + x60 + x61 + L[84];
#pragma omp atomic
Ls[85] += x*L[121] + x102 + x227 + L[85];
#pragma omp atomic
Ls[86] += x*L[122] + x269 + x270 + L[86];
#pragma omp atomic
Ls[87] += x184 + x334 + x502 + L[87];
#pragma omp atomic
Ls[88] += x*L[124] + x290 + x356 + L[88];
#pragma omp atomic
Ls[89] += x183 + x367 + x503 + L[89];
#pragma omp atomic
Ls[90] += x182 + x416 + x457 + L[90];
#pragma omp atomic
Ls[91] += x297 + x431 + x505 + L[91];
#pragma omp atomic
Ls[92] += x250 + x382 + x504 + L[92];
#pragma omp atomic
Ls[93] += x180 + x443 + x459 + L[93];
#pragma omp atomic
Ls[94] += x178 + x389 + x475 + L[94];
#pragma omp atomic
Ls[95] += x296 + x466 + x484 + L[95];
#pragma omp atomic
Ls[96] += x350 + x386 + x397 + L[96];
#pragma omp atomic
Ls[97] += x249 + x453 + x463 + L[97];
#pragma omp atomic
Ls[98] += x175 + x392 + x495 + L[98];
#pragma omp atomic
Ls[99] += x172 + x301 + x508 + L[99];
#pragma omp atomic
Ls[100] += x294 + x405 + x512 + L[100];
#pragma omp atomic
Ls[101] += x314 + x385 + x426 + L[101];
#pragma omp atomic
Ls[102] += x313 + x349 + x455 + L[102];
#pragma omp atomic
Ls[103] += x247 + x400 + x500 + L[103];
#pragma omp atomic
Ls[104] += x168 + x305 + x520 + L[104];
#pragma omp atomic
Ls[105] += x196 + x522 + x92 + L[105];
#pragma omp atomic
Ls[106] += x282 + x325 + x523 + L[106];
#pragma omp atomic
Ls[107] += x215 + x376 + x478 + L[107];
#pragma omp atomic
Ls[108] += x214 + x422 + x449 + L[108];
#pragma omp atomic
Ls[109] += x212 + x343 + x498 + L[109];
#pragma omp atomic
Ls[110] += x239 + x318 + x518 + L[110];
#pragma omp atomic
Ls[111] += x201 + x524 + x82 + L[111];
#pragma omp atomic
Ls[112] += x43 + x73 + y*L[156] + L[112];
#pragma omp atomic
Ls[113] += x260 + x275 + y*L[157] + L[113];
#pragma omp atomic
Ls[114] += x164 + x371 + x509 + L[114];
#pragma omp atomic
Ls[115] += x163 + x446 + x477 + L[115];
#pragma omp atomic
Ls[116] += x161 + x419 + x497 + L[116];
#pragma omp atomic
Ls[117] += x158 + x338 + x521 + L[117];
#pragma omp atomic
Ls[118] += x222 + x232 + x525 + L[118];
#pragma omp atomic
Ls[119] += x49 + x67 + z*L[164] + L[119];
#pragma omp atomic
Ls[120] += L[120];
#pragma omp atomic
Ls[121] += L[121];
#pragma omp atomic
Ls[122] += L[122];
#pragma omp atomic
Ls[123] += L[123];
#pragma omp atomic
Ls[124] += L[124];
#pragma omp atomic
Ls[125] += L[125];
#pragma omp atomic
Ls[126] += L[126];
#pragma omp atomic
Ls[127] += L[127];
#pragma omp atomic
Ls[128] += L[128];
#pragma omp atomic
Ls[129] += L[129];
#pragma omp atomic
Ls[130] += L[130];
#pragma omp atomic
Ls[131] += L[131];
#pragma omp atomic
Ls[132] += L[132];
#pragma omp atomic
Ls[133] += L[133];
#pragma omp atomic
Ls[134] += L[134];
#pragma omp atomic
Ls[135] += L[135];
#pragma omp atomic
Ls[136] += L[136];
#pragma omp atomic
Ls[137] += L[137];
#pragma omp atomic
Ls[138] += L[138];
#pragma omp atomic
Ls[139] += L[139];
#pragma omp atomic
Ls[140] += L[140];
#pragma omp atomic
Ls[141] += L[141];
#pragma omp atomic
Ls[142] += L[142];
#pragma omp atomic
Ls[143] += L[143];
#pragma omp atomic
Ls[144] += L[144];
#pragma omp atomic
Ls[145] += L[145];
#pragma omp atomic
Ls[146] += L[146];
#pragma omp atomic
Ls[147] += L[147];
#pragma omp atomic
Ls[148] += L[148];
#pragma omp atomic
Ls[149] += L[149];
#pragma omp atomic
Ls[150] += L[150];
#pragma omp atomic
Ls[151] += L[151];
#pragma omp atomic
Ls[152] += L[152];
#pragma omp atomic
Ls[153] += L[153];
#pragma omp atomic
Ls[154] += L[154];
#pragma omp atomic
Ls[155] += L[155];
#pragma omp atomic
Ls[156] += L[156];
#pragma omp atomic
Ls[157] += L[157];
#pragma omp atomic
Ls[158] += L[158];
#pragma omp atomic
Ls[159] += L[159];
#pragma omp atomic
Ls[160] += L[160];
#pragma omp atomic
Ls[161] += L[161];
#pragma omp atomic
Ls[162] += L[162];
#pragma omp atomic
Ls[163] += L[163];
#pragma omp atomic
Ls[164] += L[164];

}

void L2P_9(double x, double y, double z, double * L, double * F) {
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
x0 = x*y;
x1 = x*z;
x2 = y*z;
x3 = (x*x);
x4 = (1.0/2.0)*x3;
x5 = (x*x*x);
x6 = (1.0/6.0)*x5;
x7 = (x*x*x*x);
x8 = (1.0/24.0)*x7;
x9 = pow(x, 5);
x10 = (1.0/120.0)*x9;
x11 = pow(x, 6);
x12 = (1.0/720.0)*x11;
x13 = (1.0/5040.0)*pow(x, 7);
x14 = (y*y);
x15 = (1.0/2.0)*x14;
x16 = (y*y*y);
x17 = (1.0/6.0)*x16;
x18 = (y*y*y*y);
x19 = (1.0/24.0)*x18;
x20 = pow(y, 5);
x21 = (1.0/120.0)*x20;
x22 = pow(y, 6);
x23 = (1.0/720.0)*x22;
x24 = (1.0/5040.0)*pow(y, 7);
x25 = (z*z);
x26 = (1.0/2.0)*x25;
x27 = (z*z*z);
x28 = (1.0/6.0)*x27;
x29 = (z*z*z*z);
x30 = (1.0/24.0)*x29;
x31 = pow(z, 5);
x32 = (1.0/120.0)*x31;
x33 = pow(z, 6);
x34 = (1.0/720.0)*x33;
x35 = (1.0/5040.0)*pow(z, 7);
x36 = (1.0/4.0)*x3;
x37 = x14*x36;
x38 = (1.0/12.0)*x3;
x39 = x16*x38;
x40 = (1.0/48.0)*x3;
x41 = x18*x40;
x42 = (1.0/240.0)*x3;
x43 = x20*x42;
x44 = (1.0/1440.0)*x3;
x45 = x25*x36;
x46 = x27*x38;
x47 = x29*x40;
x48 = x31*x42;
x49 = (1.0/12.0)*x5;
x50 = x14*x49;
x51 = (1.0/36.0)*x5;
x52 = x16*x51;
x53 = (1.0/144.0)*x5;
x54 = x18*x53;
x55 = (1.0/720.0)*x5;
x56 = x25*x49;
x57 = x27*x51;
x58 = x29*x53;
x59 = (1.0/48.0)*x7;
x60 = x14*x59;
x61 = (1.0/144.0)*x7;
x62 = x16*x61;
x63 = (1.0/576.0)*x7;
x64 = x25*x59;
x65 = x27*x61;
x66 = (1.0/240.0)*x9;
x67 = x14*x66;
x68 = (1.0/720.0)*x9;
x69 = x25*x66;
x70 = (1.0/1440.0)*x11;
x71 = x14*x25;
x72 = (1.0/4.0)*x71;
x73 = x14*x27;
x74 = (1.0/12.0)*x73;
x75 = x14*x29;
x76 = (1.0/48.0)*x75;
x77 = (1.0/240.0)*x14*x31;
x78 = x16*x25;
x79 = (1.0/12.0)*x78;
x80 = x16*x27;
x81 = (1.0/36.0)*x80;
x82 = (1.0/144.0)*x16*x29;
x83 = x18*x25;
x84 = (1.0/48.0)*x83;
x85 = (1.0/144.0)*x18*x27;
x86 = (1.0/240.0)*x20*x25;
x87 = (1.0/24.0)*x3;
x88 = (1.0/96.0)*x3;
x89 = (1.0/72.0)*x5;
#pragma omp atomic
F[0] += (1.0/40320.0)*pow(x, 8)*L[120] + x*x15*L[13] + x*x17*L[26] + x*x19*L[45] + x*x21*L[71] + x*x23*L[105] + x*x24*L[148] + x*x26*L[15] + x*x28*L[29] + x*x30*L[49] + x*x32*L[76] + x*x34*L[111] + x*x35*L[155] + x*x72*L[47] + x*x74*L[74] + x*x76*L[109] + x*x77*L[153] + x*x79*L[73] + x*x81*L[108] + x*x82*L[152] + x*x84*L[107] + x*x85*L[151] + x*x86*L[150] + x*L[1] + x0*x26*L[28] + x0*x28*L[48] + x0*x30*L[75] + x0*x32*L[110] + x0*x34*L[154] + x0*z*L[14] + x0*L[5] + x1*x15*L[27] + x1*x17*L[46] + x1*x19*L[72] + x1*x21*L[106] + x1*x23*L[149] + x1*L[6] + x10*x2*L[88] + x10*y*L[57] + x10*z*L[58] + x10*L[35] + x12*x2*L[124] + x12*y*L[85] + x12*z*L[86] + x12*L[56] + x13*y*L[121] + x13*z*L[122] + x13*L[84] + (1.0/1440.0)*x14*x33*L[162] + x14*x70*L[123] + x15*z*L[17] + x15*L[7] + (1.0/720.0)*x16*x31*L[161] + x16*x68*L[126] + x17*z*L[31] + x17*L[16] + (1.0/576.0)*x18*x29*L[160] + x18*x63*L[130] + x19*z*L[51] + x19*L[30] + x2*x4*L[24] + x2*x6*L[39] + x2*x8*L[60] + x2*L[8] + (1.0/720.0)*x20*x27*L[159] + x20*x55*L[135] + x21*z*L[78] + x21*L[50] + (1.0/1440.0)*x22*x25*L[158] + x22*x44*L[141] + x23*z*L[113] + x23*L[77] + x24*z*L[157] + x24*L[112] + x25*x70*L[125] + x26*y*L[18] + x26*L[9] + x27*x68*L[129] + x28*y*L[33] + x28*L[19] + x29*x63*L[134] + (1.0/8.0)*x3*x71*L[68] + (1.0/72.0)*x3*x80*L[144] + x30*y*L[54] + x30*L[34] + x31*x55*L[140] + x32*y*L[82] + x32*L[55] + x33*x44*L[147] + x34*y*L[118] + x34*L[83] + x35*y*L[163] + x35*L[119] + x37*z*L[42] + x37*L[23] + x39*z*L[67] + x39*L[41] + x4*y*L[11] + x4*z*L[12] + x4*L[4] + x41*z*L[100] + x41*L[66] + x43*z*L[142] + x43*L[99] + x45*y*L[43] + x45*L[25] + x46*y*L[69] + x46*L[44] + x47*y*L[103] + x47*L[70] + x48*y*L[146] + x48*L[104] + (1.0/24.0)*x5*x71*L[96] + x50*z*L[63] + x50*L[38] + x52*z*L[95] + x52*L[62] + x54*z*L[136] + x54*L[94] + x56*y*L[64] + x56*L[40] + x57*y*L[97] + x57*L[65] + x58*y*L[139] + x58*L[98] + x6*y*L[21] + x6*z*L[22] + x6*L[10] + x60*z*L[91] + x60*L[59] + x62*z*L[131] + x62*L[90] + x64*y*L[92] + x64*L[61] + x65*y*L[133] + x65*L[93] + x67*z*L[127] + x67*L[87] + x69*y*L[128] + x69*L[89] + (1.0/96.0)*x7*x71*L[132] + x72*L[32] + x73*x87*L[102] + x73*x89*L[138] + x74*L[53] + x75*x88*L[145] + x76*L[81] + x77*L[117] + x78*x87*L[101] + x78*x89*L[137] + x79*L[52] + x8*y*L[36] + x8*z*L[37] + x8*L[20] + x81*L[80] + x82*L[116] + x83*x88*L[143] + x84*L[79] + x85*L[115] + x86*L[114] + (1.0/40320.0)*pow(y, 8)*L[156] + y*L[2] + (1.0/40320.0)*pow(z, 8)*L[164] + z*L[3] + L[0];

}

void M2P_9(double x, double y, double z, double * M, double * F) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*y;
x8 = x*z;
x9 = y*z;
x10 = pow(x2, -3.5);
x11 = 15.0*x10;
x12 = x7*z;
x13 = -x4;
x14 = x0*x6;
x15 = x1*x6;
x16 = 9.0*x5;
x17 = x0*x11;
x18 = -x17;
x19 = x*(x16 + x18);
x20 = x18 + x6;
x21 = x20*y;
x22 = x1*x11;
x23 = -x22;
x24 = y*(x16 + x23);
x25 = x20*z;
x26 = z*(x23 + x6);
x27 = 1.0*x;
x28 = x27*(x22 - x6);
x29 = 45.0*x10;
x30 = -x29;
x31 = pow(x2, -4.5);
x32 = x0*x31;
x33 = 105.0*x32;
x34 = x30 + x33;
x35 = x34*x7;
x36 = x34*x8;
x37 = -x11;
x38 = x9*(x33 + x37);
x39 = x1*x31;
x40 = 105.0*x39;
x41 = x30 + x40;
x42 = x41*x9;
x43 = x27*y;
x44 = x41*x43;
x45 = x37 + x40;
x46 = x27*z;
x47 = x45*x46;
x48 = 315.0*x31;
x49 = pow(x2, -5.5);
x50 = 945.0*x49;
x51 = x0*x50;
x52 = x48 - x51;
x53 = x12*x52;
x54 = x1*x50;
x55 = x27*x9;
x56 = x55*(-x48 + x54);
x57 = 90.0*x10;
x58 = (x*x*x*x);
x59 = 105.0*x31;
x60 = x58*x59;
x61 = (y*y*y*y);
x62 = x59*x61;
x63 = -225.0*x10;
x64 = x50*x58;
x65 = -x64;
x66 = x*(1050.0*x32 + x63 + x65);
x67 = x30 + 630.0*x32 + x65;
x68 = x67*y;
x69 = x50*x61;
x70 = -x69;
x71 = y*(1050.0*x39 + x63 + x70);
x72 = x67*z;
x73 = 630.0*x39;
x74 = x30 + x70 + x73;
x75 = x74*z;
x76 = x27*(x29 + x69 - x73);
x77 = 1575.0*x31;
x78 = x0*x49;
x79 = pow(x2, -6.5);
x80 = x58*x79;
x81 = 10395.0*x80;
x82 = x77 - 9450.0*x78 + x81;
x83 = x7*x82;
x84 = x8*x82;
x85 = 5670.0*x78;
x86 = x48 + x81 - x85;
x87 = x86*x9;
x88 = x1*x49;
x89 = x61*x79;
x90 = 10395.0*x89;
x91 = x77 - 9450.0*x88 + x90;
x92 = x9*x91;
x93 = x43*x91;
x94 = x48 - 5670.0*x88 + x90;
x95 = x46*x94;
x96 = 14175.0*x49;
x97 = -x96;
x98 = x0*x79;
x99 = 103950.0*x98;
x100 = pow(x2, -7.5);
x101 = 135135.0*x100;
x102 = x101*x58;
x103 = -x102 + x97 + x99;
x104 = x103*x12;
x105 = x1*x79;
x106 = 103950.0*x105;
x107 = x101*x61;
x108 = x55*(-x106 + x107 + x96);
x109 = pow(x, 6);
x110 = 10395.0*x79;
x111 = x109*x110;
x112 = pow(y, 6);
x113 = x110*x112;
x114 = 11025.0*x31;
x115 = 99225.0*x49;
x116 = x101*x109;
x117 = -x116;
x118 = x117 + 218295.0*x80;
x119 = x*(-x0*x115 + x114 + x118);
x120 = 42525.0*x49;
x121 = -x0*x120 + x117 + x77 + 155925.0*x80;
x122 = x121*y;
x123 = x101*x112;
x124 = -x123;
x125 = x124 + 218295.0*x89;
x126 = y*(-x1*x115 + x114 + x125);
x127 = x121*z;
x128 = x1*x120;
x129 = 155925.0*x89;
x130 = x124 - x128 + x129 + x77;
x131 = x130*z;
x132 = x27*(x123 + x128 - x129 - x77);
x133 = pow(x2, -8.5);
x134 = x109*x133;
x135 = 2027025.0*x134;
x136 = x100*x58;
x137 = 2837835.0*x136;
x138 = -x115;
x139 = 1091475.0*x79;
x140 = x0*x139 + x138;
x141 = x135 - x137 + x140;
x142 = x141*x7;
x143 = x141*x8;
x144 = 467775.0*x79;
x145 = x0*x144;
x146 = 2027025.0*x136;
x147 = -x146;
x148 = x135 + x145 + x147 + x97;
x149 = x148*x9;
x150 = x1*x139;
x151 = x112*x133;
x152 = 2027025.0*x151;
x153 = x100*x61;
x154 = 2837835.0*x153;
x155 = x150 + x152 - x154;
x156 = x138 + x155;
x157 = x156*x9;
x158 = x1*x33;
x159 = x156*x43;
x160 = 2027025.0*x153;
x161 = -x160;
x162 = x1*x144 + x152 + x161 + x97;
x163 = x162*x46;
x164 = x0*x100;
x165 = 14189175.0*x164;
x166 = x133*x58;
x167 = 42567525.0*x166;
x168 = 34459425.0*pow(x2, -9.5);
x169 = x109*x168;
x170 = x12*(x139 - x165 + x167 - x169);
x171 = x1*x48;
x172 = -x1*x51;
x173 = x*(x171 + x172 + x34);
x174 = x0*x48;
x175 = y*(x172 + x174 + x41);
x176 = z*(x172 + x33 + x45);
x177 = x1*x100;
x178 = x133*x61;
x179 = 42567525.0*x178;
x180 = x112*x168;
x181 = x55*(-x139 + 14189175.0*x177 - x179 + x180);
x182 = -2835.0*x88;
x183 = x1*x98;
x184 = 10395.0*x183;
x185 = x182 + x184;
x186 = 945.0*x31;
x187 = -2835.0*x78;
x188 = x186 + x187;
x189 = x7*(x185 + x188);
x190 = x8*(x185 + x52);
x191 = x9*(x184 + x187 + x48 - x54);
x192 = 31185.0*x105;
x193 = x0*x1;
x194 = -8505.0*x49;
x195 = 31185.0*x98;
x196 = x194 + x195;
x197 = x12*(-x101*x193 + x192 + x196);
x198 = 3783780.0*x100;
x199 = pow(x, 8);
x200 = 2027025.0*x133;
x201 = x199*x200;
x202 = pow(y, 8);
x203 = x200*x202;
x204 = -893025.0*x49;
x205 = -x168*x199;
x206 = x*(72972900.0*x134 - 51081030.0*x136 + x204 + x205 + 13097700.0*x98);
x207 = 56756700.0*x134 - 28378350.0*x136 + x138 + x205 + 4365900.0*x98;
x208 = x207*y;
x209 = x168*x202;
x210 = -x209;
x211 = y*(13097700.0*x105 + 72972900.0*x151 - 51081030.0*x153 + x204 + x210);
x212 = x207*z;
x213 = 4365900.0*x105;
x214 = 28378350.0*x153;
x215 = 56756700.0*x151;
x216 = z*(x138 + x210 + x213 - x214 + x215);
x217 = x27*(x115 + x209 - x213 + x214 - x215);
x218 = 105.0*x10;
x219 = -x158 - 12.0*x5;
x220 = 120.0*x10;
x221 = x61*x98;
x222 = 10395.0*x221;
x223 = x1*x85;
x224 = -x223;
x225 = x1*x81;
x226 = 62370.0*x183;
x227 = -x0*x107;
x228 = x226 + x227;
x229 = 31185.0*x89;
x230 = x229 - 17010.0*x88;
x231 = x*(x188 + x228 + x230);
x232 = x1*x96;
x233 = -x232;
x234 = x1*x99;
x235 = -x1*x102;
x236 = x*(x233 + x234 + x235 + x82);
x237 = x0*x96;
x238 = -x237;
x239 = y*(x227 + x234 + x238 + x91);
x240 = x182 + x226 + x235;
x241 = 31185.0*x80;
x242 = 17010.0*x78;
x243 = x186 + x241 - x242;
x244 = y*(x240 + x243);
x245 = z*(x187 + x228 + x94);
x246 = z*(x240 + x86);
x247 = 155925.0*x98;
x248 = 311850.0*x105;
x249 = -x120;
x250 = x1*x164;
x251 = -1351350.0*x250;
x252 = x249 + x251;
x253 = -405405.0*x153;
x254 = x0*x178;
x255 = 2027025.0*x254;
x256 = x253 + x255;
x257 = x7*(x247 + x248 + x252 + x256);
x258 = 155925.0*x105;
x259 = 311850.0*x98;
x260 = -405405.0*x136;
x261 = x1*x166;
x262 = 2027025.0*x261;
x263 = x260 + x262;
x264 = x7*(x252 + x258 + x259 + x263);
x265 = 187110.0*x105;
x266 = -810810.0*x250;
x267 = x8*(x196 + x256 + x265 + x266);
x268 = x8*(x103 + x251 + x258 + x262);
x269 = x9*(x106 - x107 + x247 + x251 + x255 + x97);
x270 = x194 + 187110.0*x98;
x271 = x9*(x192 + x263 + x266 + x270);
x272 = 20270250.0*x133*x193 + x144;
x273 = x12*(-x0*x168*x61 - 2027025.0*x164 - 4054050.0*x177 + 6081075.0*x178 + x272);
x274 = 6081075.0*x166;
x275 = x168*x58;
x276 = x12*(-x1*x275 - 4054050.0*x164 - 2027025.0*x177 + x272 + x274);
x277 = 15120.0*x49;
x278 = -x111;
x279 = -x225;
x280 = 270.0*x10 + x223;
x281 = -x113;
x282 = -x222;
x283 = -x173;
x284 = -x175;
x285 = -x176;
x286 = -x189;
x287 = -x190;
x288 = -x191;
x289 = -x197;
x290 = x1*x145;
x291 = x0*x152;
x292 = x0*x160;
x293 = x1*x78;
x294 = x1*x135;
x295 = x1*x146;
x296 = 6081075.0*x151;
x297 = -6081075.0*x250;
x298 = x249 + x297;
x299 = -x0*x180;
x300 = x247 + 30405375.0*x254 + x299;
x301 = x*(1403325.0*x105 - 6081075.0*x153 + x296 + x298 + x300);
x302 = -x1*x165;
x303 = -x1*x169;
x304 = x*(x1*x167 + x141 + x150 + x302 + x303);
x305 = y*(x0*x179 + x140 + x155 + x299 + x302);
x306 = 6081075.0*x134;
x307 = x258 + 30405375.0*x261 + x303;
x308 = y*(-6081075.0*x136 + x298 + x306 + x307 + 1403325.0*x98);
x309 = z*(x162 + x297 + x300);
x310 = z*(x148 + x297 + x307);
x311 = 16065.0*x49;
x312 = -x1*x242 - 360.0*x10;
x313 = 3918915.0*x100;
x314 = -x201;
x315 = -x294;
x316 = -x290 - 12600.0*x31;
x317 = -x203;
x318 = -x291;
x319 = x0*x153;
x320 = 810810.0*x319;
x321 = -x320;
x322 = x166*x61;
x323 = 2027025.0*x322;
x324 = x1*x136;
x325 = 810810.0*x324;
x326 = x323 - x325;
x327 = -x275*x61;
x328 = x249 - 8108100.0*x250 + x327;
x329 = x260 + 12162150.0*x261;
x330 = x*(935550.0*x105 + x161 + 20270250.0*x254 + x259 + x328 + x329);
x331 = x253 + 12162150.0*x254;
x332 = y*(x147 + x248 + 20270250.0*x261 + x328 + x331 + 935550.0*x98);
x333 = z*(-4864860.0*x250 + x265 + x270 + x327 + x329 + x331);
x334 = 17010.0*x49;
x335 = -x119;
x336 = -x236;
x337 = -x231;
x338 = -x122;
x339 = -x244;
x340 = -x126;
x341 = -x239;
x342 = -x127;
x343 = -x246;
x344 = -x131;
x345 = -x245;
x346 = -x142;
x347 = -x264;
x348 = -x159;
x349 = -x257;
x350 = -x143;
x351 = -x268;
x352 = -x163;
x353 = -x267;
x354 = -x149;
x355 = -x271;
x356 = -x157;
x357 = -x269;
x358 = -x170;
x359 = -x276;
x360 = -x273;
x361 = x123 + x318;
x362 = -841995.0*x183 - 2520.0*x31 - x323;
x363 = x116 + x315;
x364 = 4054050.0*x100;
x365 = x1*x134;
x366 = 1309770.0*x183 + 15120.0*x31;
x367 = x0*x151;
x368 = -x206;
x369 = -x304;
x370 = -x301;
x371 = -x208;
x372 = -x308;
x373 = -x211;
x374 = -x305;
x375 = -x212;
x376 = -x310;
x377 = -x216;
x378 = -x309;
x379 = 4189185.0*x100;
x380 = -2993760.0*x183 - x274*x61 - 20160.0*x31;
x381 = 4324320.0*x100;
x382 = -x330;
x383 = -x332;
x384 = -x333;
x385 = -x217;
x386 = -3*x330;
x387 = -3*x332;
x388 = -3*x333;
#pragma omp atomic
F[0] += -x*x4*M[0] + x104*M[87] - x108*M[105] - x11*x12*M[13] + x119*M[83] + x122*M[84] + x126*M[111] + x127*M[85] + x131*M[112] - x132*M[104] + x142*M[120] + x143*M[121] + x149*M[123] + x157*M[156] + x159*M[147] + x163*M[148] + x170*M[168] + x173*M[37] + x175*M[40] + x176*M[41] - x181*M[201] + x189*M[61] + x19*M[9] + x190*M[62] + x191*M[66] + x197*M[94] + x206*M[164] + x208*M[165] + x21*M[10] + x211*M[209] + x212*M[166] + x216*M[210] - x217*M[200] + x231*M[93] + x236*M[86] + x239*M[98] + x24*M[15] + x244*M[89] + x245*M[99] + x246*M[90] + x25*M[11] + x257*M[134] + x26*M[16] + x264*M[125] + x267*M[135] + x268*M[126] + x269*M[141] + x271*M[130] + x273*M[186] + x276*M[175] - x28*M[12] + x301*M[185] + x304*M[167] + x305*M[192] + x308*M[170] + x309*M[193] + x310*M[171] + x330*M[174] + x332*M[179] + x333*M[180] + x35*M[20] + x36*M[21] + x38*M[23] - x4*y*M[1] - x4*z*M[2] + x42*M[30] + x44*M[25] + x47*M[26] + x53*M[38] - x56*M[45] + x6*x7*M[4] + x6*x8*M[5] + x6*x9*M[7] + x66*M[34] + x68*M[35] + x71*M[49] + x72*M[36] + x75*M[50] - x76*M[44] + x83*M[56] + x84*M[57] + x87*M[59] + x92*M[77] + x93*M[70] + x95*M[71] + (-x104 + x289)*M[96] + (x108 + x289)*M[107] + (x13 + x14)*M[3] + (x13 + x15)*M[6] + (x132 + x337)*M[106] + (x181 + x360)*M[203] + (-x19 + x28)*M[14] + (-x21 - x24)*M[17] + (x217 + x370)*M[202] + (-x25 - x26)*M[18] + (x283 - x66)*M[39] + (x283 + x76)*M[46] + (x284 - x68)*M[42] + (x284 - x71)*M[51] + (x285 - x72)*M[43] + (x285 - x75)*M[52] + (x286 - x83)*M[63] + (x286 - x93)*M[72] + (x287 - x84)*M[64] + (x287 - x95)*M[73] + (x288 - x87)*M[68] + (x288 - x92)*M[79] + (x335 + x336)*M[88] + (x336 + x337)*M[95] + (x338 + x339)*M[91] + (x339 + x341)*M[100] + (x340 + x341)*M[113] + (x342 + x343)*M[92] + (x343 + x345)*M[101] + (x344 + x345)*M[114] + (x346 + x347)*M[127] + (x347 + x349)*M[136] + (x348 + x349)*M[149] + (-x35 - x44)*M[27] + (x350 + x351)*M[128] + (x351 + x353)*M[137] + (x352 + x353)*M[150] + (x354 + x355)*M[132] + (x355 + x357)*M[143] + (x356 + x357)*M[158] + (x358 + x359)*M[177] + (x359 + x360)*M[188] + (-x36 - x47)*M[28] + (x368 + x369)*M[169] + (x369 + x382)*M[176] + (x370 + x382)*M[187] + (x371 + x372)*M[172] + (x372 + x383)*M[181] + (x373 + x374)*M[211] + (x374 + x383)*M[194] + (x375 + x376)*M[173] + (x376 + x384)*M[182] + (x377 + x378)*M[212] + (x378 + x384)*M[195] + (-x38 - x42)*M[32] + (-x53 + x56)*M[47] + (x104 - x108 + 2*x197)*M[109] + (x119 + x231 + 2*x236)*M[97] + (x122 + x239 + 2*x244)*M[102] + (x126 + 2*x239 + x244)*M[115] + (x127 + x245 + 2*x246)*M[103] + (x131 + 2*x245 + x246)*M[116] + (-x132 + 2*x231 + x236)*M[108] + (-x14 - x15 + 2.0*x3)*M[8] + (x142 + x257 + 2*x264)*M[138] + (x143 + x267 + 2*x268)*M[139] + (x149 + x269 + 2*x271)*M[145] + (x157 + 2*x269 + x271)*M[160] + (x158 + x20 + x23)*M[22] + (x159 + 2*x257 + x264)*M[151] + (x163 + 2*x267 + x268)*M[152] + (x170 + x273 + 2*x276)*M[190] + (2*x173 + x66 - x76)*M[48] + (2*x175 + x68 + x71)*M[53] + (2*x176 + x72 + x75)*M[54] + (-x181 + 2*x273 + x276)*M[205] + (2*x189 + x83 + x93)*M[74] + (2*x190 + x84 + x95)*M[75] + (2*x191 + x87 + x92)*M[81] + (x206 + 2*x304 + x330)*M[178] + (x208 + 2*x308 + x332)*M[183] + (x211 + 2*x305 + x332)*M[213] + (x212 + 2*x310 + x333)*M[184] + (x216 + 2*x309 + x333)*M[214] + (x301 + x304 + 2*x330)*M[189] + (2*x301 + x330 + x385)*M[204] + (x305 + x308 + 2*x332)*M[196] + (x309 + x310 + 2*x333)*M[197] + (-x0*x57 + x16 + x60)*M[19] + (-x1*x57 + x16 + x62)*M[29] + (x111 + 4725.0*x32 - x58*x96 + x63)*M[55] + (x113 + 4725.0*x39 - x61*x96 + x63)*M[76] + (x132 - 3*x231 - 3*x236 + x335)*M[110] + (x171 + x224 + x225 + x67)*M[58] + (x174 + x222 + x224 + x74)*M[65] + (x181 - 3*x273 - 3*x276 + x358)*M[207] + (x217 - 3*x301 + x369 + x386)*M[206] + (-3*x239 - 3*x244 + x338 + x340)*M[117] + (-3*x245 - 3*x246 + x342 + x344)*M[118] + (-3*x257 - 3*x264 + x346 + x348)*M[153] + (-3*x267 - 3*x268 + x350 + x352)*M[154] + (-3*x269 - 3*x271 + x354 + x356)*M[162] + (-3*x304 + x368 + x370 + x386)*M[191] + (-3*x305 + x372 + x373 + x387)*M[215] + (-3*x308 + x371 + x374 + x387)*M[198] + (-3*x309 + x376 + x377 + x388)*M[216] + (-3*x310 + x375 + x378 + x388)*M[199] + (x0*x218 + x219 + x22 - x60)*M[24] + (x1*x218 + x17 + x219 - x62)*M[31] + (x121 + x233 + x290 + x294 - x295)*M[122] + (x130 + x238 + x290 + x291 - x292)*M[140] + (374220.0*x183 + x230 + x243 + x321 + x326)*M[129] + (x206 + 4*x301 + 4*x304 + 6*x330 + x385)*M[208] + (x208 + x211 + 4*x305 + 4*x308 + 6*x332)*M[217] + (x212 + x216 + 4*x309 + 4*x310 + 6*x333)*M[218] + (-x109*x198 + x114 + x201 - 396900.0*x78 + 2182950.0*x80)*M[119] + (-x112*x198 + x114 + x203 - 396900.0*x88 + 2182950.0*x89)*M[155] + (-x171 + x277*x58 + x278 + x279 + x280 - 5355.0*x32)*M[60] + (-x174 + x277*x61 + x280 + x281 + x282 - 5355.0*x39)*M[78] + (-x0*x220 - x1*x220 + 210.0*x1*x32 + 24.0*x5 + x60 + x62)*M[33] + (x113 + 20790.0*x221 + x225 - x311*x61 + x312 + 1260.0*x32 + 6300.0*x39 + x65)*M[80] + (x0*x154 - x241 + x325 + x361 + x362 + 31185.0*x78 + 59535.0*x88 - 187110.0*x89)*M[142] + (-x0*x186 - x1*x186 + x279 + x282 + 11340.0*x293 + x57 + x64 + x69)*M[67] + (x1*x137 - x229 + x320 + x362 + x363 + 59535.0*x78 - 187110.0*x80 + 31185.0*x88)*M[131] + (20790.0*x1*x80 + x111 + x222 - x311*x58 + x312 + 6300.0*x32 + 1260.0*x39 + x70)*M[69] + (x109*x313 + x232 + x295 + x314 + x315 + x316 + 439425.0*x78 - 2338875.0*x80)*M[124] + (x112*x313 + x237 + x292 + x316 + x317 + x318 + 439425.0*x88 - 2338875.0*x89)*M[157] + (-x1*x241 + 720.0*x10 - x195*x61 + x278 + x281 + 34020.0*x293 - 7560.0*x32 + x334*x58 + x334*x61 - 7560.0*x39)*M[82] + (-x112*x364 + x203 + x241 - 4864860.0*x319 + x326 + x366 + 4054050.0*x367 - 45360.0*x78 - 498960.0*x88 + 2525985.0*x89)*M[159] + (x118 + x125 + 1683990.0*x183 + x291 + x294 + 5040.0*x31 - 3648645.0*x319 + 4054050.0*x322 - 3648645.0*x324 - 90720.0*x78 - 90720.0*x88)*M[144] + (-x0*x296 + x112*x379 + x317 + 8513505.0*x319 + 4459455.0*x324 + x363 + x380 + 136080.0*x78 - 249480.0*x80 + 589680.0*x88 - 2744280.0*x89)*M[161] + (-x1*x306 + x109*x379 + x314 + 4459455.0*x319 + 8513505.0*x324 + x361 + x380 + 589680.0*x78 - 2744280.0*x80 + 136080.0*x88 - 249480.0*x89)*M[146] + (-x109*x364 + x201 + x229 + x321 + x323 - 4864860.0*x324 + 4054050.0*x365 + x366 - 498960.0*x78 + 2525985.0*x80 - 45360.0*x88)*M[133] + (-x109*x381 - x112*x381 + 5987520.0*x183 + x201 + x203 + 40320.0*x31 - 12972960.0*x319 + 12162150.0*x322 - 12972960.0*x324 + 8108100.0*x365 + 8108100.0*x367 - 725760.0*x78 + 2993760.0*x80 - 725760.0*x88 + 2993760.0*x89)*M[163];

}

void P2M_10(double x, double y, double z, double q, double * M) {
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
x0 = q*x;
x1 = q*y;
x2 = q*z;
x3 = (x*x);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = x0*z;
x7 = (y*y);
x8 = x1*z;
x9 = (z*z);
x10 = (x*x*x);
x11 = (1.0/6.0)*q;
x12 = (1.0/2.0)*x3;
x13 = (1.0/2.0)*x0;
x14 = (y*y*y);
x15 = (1.0/2.0)*x7;
x16 = (1.0/2.0)*x9;
x17 = (z*z*z);
x18 = (x*x*x*x);
x19 = (1.0/24.0)*q;
x20 = (1.0/6.0)*x10;
x21 = q*x7;
x22 = (1.0/4.0)*x3;
x23 = q*x9;
x24 = (1.0/6.0)*x0;
x25 = (y*y*y*y);
x26 = (1.0/6.0)*x14;
x27 = (1.0/4.0)*x9;
x28 = (1.0/6.0)*x17;
x29 = (z*z*z*z);
x30 = pow(x, 5);
x31 = (1.0/120.0)*q;
x32 = (1.0/24.0)*x18;
x33 = (1.0/12.0)*x10;
x34 = (1.0/12.0)*x14;
x35 = q*x3;
x36 = x2*x7;
x37 = x1*x9;
x38 = (1.0/12.0)*x17;
x39 = (1.0/24.0)*x0;
x40 = x0*x7;
x41 = pow(y, 5);
x42 = (1.0/24.0)*x25;
x43 = (1.0/24.0)*x29;
x44 = pow(z, 5);
x45 = pow(x, 6);
x46 = (1.0/720.0)*q;
x47 = (1.0/120.0)*x30;
x48 = (1.0/48.0)*x18;
x49 = q*x14;
x50 = (1.0/36.0)*x10;
x51 = q*x17;
x52 = (1.0/48.0)*x35;
x53 = x2*x3;
x54 = x3*x9;
x55 = x1*x3;
x56 = (1.0/120.0)*x0;
x57 = x0*x9;
x58 = pow(y, 6);
x59 = (1.0/120.0)*x41;
x60 = (1.0/48.0)*x25;
x61 = (1.0/36.0)*x17;
x62 = (1.0/48.0)*x29;
x63 = (1.0/120.0)*x44;
x64 = pow(z, 6);
x65 = pow(x, 7);
x66 = (1.0/5040.0)*q;
x67 = (1.0/720.0)*x45;
x68 = (1.0/240.0)*x30;
x69 = (1.0/144.0)*x18;
x70 = (1.0/144.0)*x25;
x71 = q*x10;
x72 = x14*x2;
x73 = x19*x7;
x74 = x1*x17;
x75 = (1.0/144.0)*x29;
x76 = (1.0/240.0)*x35;
x77 = (1.0/720.0)*x0;
x78 = x0*x14;
x79 = pow(y, 7);
x80 = (1.0/720.0)*x58;
x81 = (1.0/240.0)*x41;
x82 = (1.0/240.0)*x44;
x83 = (1.0/720.0)*x64;
x84 = pow(z, 7);
x85 = pow(x, 8);
x86 = (1.0/40320.0)*q;
x87 = (1.0/5040.0)*x65;
x88 = (1.0/1440.0)*x45;
x89 = x30*x46;
x90 = q*x25;
x91 = (1.0/576.0)*x18;
x92 = (1.0/96.0)*x21;
x93 = q*x29;
x94 = x10*x46;
x95 = x10*x2;
x96 = (1.0/72.0)*x10;
x97 = x14*x23;
x98 = x17*x21;
x99 = x1*x10;
x100 = (1.0/1440.0)*x35;
x101 = x23*x25;
x102 = x17*x35;
x103 = (1.0/5040.0)*x0;
x104 = x0*x17;
x105 = pow(y, 8);
x106 = (1.0/5040.0)*x79;
x107 = (1.0/1440.0)*x58;
x108 = x17*x41;
x109 = (1.0/576.0)*x29;
x110 = x14*x44;
x111 = (1.0/1440.0)*x64;
x112 = (1.0/5040.0)*x84;
x113 = pow(z, 8);
x114 = pow(x, 9);
x115 = (1.0/362880.0)*q;
x116 = (1.0/40320.0)*x85;
x117 = (1.0/10080.0)*x65;
x118 = (1.0/4320.0)*x45;
x119 = (1.0/2880.0)*x30;
x120 = (1.0/720.0)*x30;
x121 = (1.0/480.0)*x21;
x122 = (1.0/2880.0)*x41;
x123 = q*x18;
x124 = x2*x25;
x125 = (1.0/288.0)*x18;
x126 = x1*x29;
x127 = (1.0/2880.0)*x44;
x128 = (1.0/4320.0)*x71;
x129 = (1.0/288.0)*x10;
x130 = x17*x49;
x131 = x21*x29;
x132 = (1.0/10080.0)*x35;
x133 = x23*x41;
x134 = x29*x35;
x135 = (1.0/40320.0)*x0;
x136 = x0*x25;
x137 = pow(y, 9);
x138 = (1.0/40320.0)*x105;
x139 = (1.0/10080.0)*x79;
x140 = (1.0/4320.0)*x58;
x141 = (1.0/4320.0)*x64;
x142 = (1.0/10080.0)*x84;
x143 = (1.0/40320.0)*x113;
x144 = pow(z, 9);
x145 = (1.0/3628800.0)*q;
x146 = (1.0/362880.0)*x114;
x147 = (1.0/80640.0)*x85;
x148 = (1.0/30240.0)*x65;
x149 = (1.0/17280.0)*x45;
x150 = (1.0/2880.0)*x21;
x151 = (1.0/14400.0)*q*x30;
x152 = (1.0/1440.0)*x30;
x153 = (1.0/17280.0)*x123;
x154 = (1.0/1152.0)*x18;
x155 = (1.0/30240.0)*x71;
x156 = (1.0/1440.0)*x10;
x157 = (1.0/864.0)*x10;
x158 = (1.0/80640.0)*x35;
x159 = (1.0/362880.0)*x0;
M[0] += -x0;
M[1] += -x1;
M[2] += -x2;
M[3] += x3*x4;
M[4] += x5;
M[5] += x6;
M[6] += x4*x7;
M[7] += x8;
M[8] += x4*x9;
M[9] += -x10*x11;
M[10] += -x1*x12;
M[11] += -x12*x2;
M[12] += -x13*x7;
M[13] += -x5*z;
M[14] += -x13*x9;
M[15] += -x11*x14;
M[16] += -x15*x2;
M[17] += -x1*x16;
M[18] += -x11*x17;
M[19] += x18*x19;
M[20] += x1*x20;
M[21] += x2*x20;
M[22] += x21*x22;
M[23] += x12*x8;
M[24] += x22*x23;
M[25] += x14*x24;
M[26] += x15*x6;
M[27] += x16*x5;
M[28] += x17*x24;
M[29] += x19*x25;
M[30] += x2*x26;
M[31] += x21*x27;
M[32] += x1*x28;
M[33] += x19*x29;
M[34] += -x30*x31;
M[35] += -x1*x32;
M[36] += -x2*x32;
M[37] += -x21*x33;
M[38] += -x20*x8;
M[39] += -x23*x33;
M[40] += -x34*x35;
M[41] += -x22*x36;
M[42] += -x22*x37;
M[43] += -x35*x38;
M[44] += -x25*x39;
M[45] += -x26*x6;
M[46] += -x27*x40;
M[47] += -x28*x5;
M[48] += -x29*x39;
M[49] += -x31*x41;
M[50] += -x2*x42;
M[51] += -x23*x34;
M[52] += -x21*x38;
M[53] += -x1*x43;
M[54] += -x31*x44;
M[55] += x45*x46;
M[56] += x1*x47;
M[57] += x2*x47;
M[58] += x21*x48;
M[59] += x32*x8;
M[60] += x23*x48;
M[61] += x49*x50;
M[62] += x33*x36;
M[63] += x33*x37;
M[64] += x50*x51;
M[65] += x25*x52;
M[66] += x34*x53;
M[67] += (1.0/8.0)*x21*x54;
M[68] += x38*x55;
M[69] += x29*x52;
M[70] += x41*x56;
M[71] += x42*x6;
M[72] += x34*x57;
M[73] += x38*x40;
M[74] += x43*x5;
M[75] += x44*x56;
M[76] += x46*x58;
M[77] += x2*x59;
M[78] += x23*x60;
M[79] += x49*x61;
M[80] += x21*x62;
M[81] += x1*x63;
M[82] += x46*x64;
M[83] += -x65*x66;
M[84] += -x1*x67;
M[85] += -x2*x67;
M[86] += -x21*x68;
M[87] += -x47*x8;
M[88] += -x23*x68;
M[89] += -x49*x69;
M[90] += -x36*x48;
M[91] += -x37*x48;
M[92] += -x51*x69;
M[93] += -x70*x71;
M[94] += -x50*x72;
M[95] += -x10*x73*x9;
M[96] += -x50*x74;
M[97] += -x71*x75;
M[98] += -x41*x76;
M[99] += -x53*x60;
M[100] += -x14*x19*x54;
M[101] += -x17*x3*x73;
M[102] += -x55*x62;
M[103] += -x44*x76;
M[104] += -x58*x77;
M[105] += -x59*x6;
M[106] += -x57*x60;
M[107] += -x61*x78;
M[108] += -x40*x62;
M[109] += -x5*x63;
M[110] += -x64*x77;
M[111] += -x66*x79;
M[112] += -x2*x80;
M[113] += -x23*x81;
M[114] += -x51*x70;
M[115] += -x49*x75;
M[116] += -x21*x82;
M[117] += -x1*x83;
M[118] += -x66*x84;
M[119] += x85*x86;
M[120] += x1*x87;
M[121] += x2*x87;
M[122] += x21*x88;
M[123] += x67*x8;
M[124] += x23*x88;
M[125] += x14*x89;
M[126] += x36*x68;
M[127] += x37*x68;
M[128] += x17*x89;
M[129] += x90*x91;
M[130] += x69*x72;
M[131] += x18*x9*x92;
M[132] += x69*x74;
M[133] += x91*x93;
M[134] += x41*x94;
M[135] += x70*x95;
M[136] += x96*x97;
M[137] += x96*x98;
M[138] += x75*x99;
M[139] += x44*x94;
M[140] += x100*x58;
M[141] += x53*x81;
M[142] += (1.0/96.0)*x101*x3;
M[143] += (1.0/72.0)*x102*x14;
M[144] += x29*x3*x92;
M[145] += x55*x82;
M[146] += x100*x64;
M[147] += x103*x79;
M[148] += x6*x80;
M[149] += x57*x81;
M[150] += x104*x70;
M[151] += x75*x78;
M[152] += x40*x82;
M[153] += x5*x83;
M[154] += x103*x84;
M[155] += x105*x86;
M[156] += x106*x2;
M[157] += x107*x23;
M[158] += x108*x46;
M[159] += x109*x90;
M[160] += x110*x46;
M[161] += x111*x21;
M[162] += x1*x112;
M[163] += x113*x86;
M[164] += -x114*x115;
M[165] += -x1*x116;
M[166] += -x116*x2;
M[167] += -x117*x21;
M[168] += -x8*x87;
M[169] += -x117*x23;
M[170] += -x118*x49;
M[171] += -x36*x88;
M[172] += -x37*x88;
M[173] += -x118*x51;
M[174] += -x119*x90;
M[175] += -x120*x72;
M[176] += -x121*x30*x9;
M[177] += -x120*x74;
M[178] += -x119*x93;
M[179] += -x122*x123;
M[180] += -x124*x91;
M[181] += -x125*x97;
M[182] += -x125*x98;
M[183] += -x126*x91;
M[184] += -x123*x127;
M[185] += -x128*x58;
M[186] += -1.0/720.0*x41*x95;
M[187] += -x101*x129;
M[188] += -1.0/216.0*x10*x130;
M[189] += -x129*x131;
M[190] += -1.0/720.0*x44*x99;
M[191] += -x128*x64;
M[192] += -x132*x79;
M[193] += -x107*x53;
M[194] += -1.0/480.0*x133*x3;
M[195] += -1.0/288.0*x102*x25;
M[196] += -1.0/288.0*x134*x14;
M[197] += -x121*x3*x44;
M[198] += -x111*x55;
M[199] += -x132*x84;
M[200] += -x105*x135;
M[201] += -x106*x6;
M[202] += -x107*x57;
M[203] += -x108*x77;
M[204] += -x109*x136;
M[205] += -x110*x77;
M[206] += -x111*x40;
M[207] += -x112*x5;
M[208] += -x113*x135;
M[209] += -x115*x137;
M[210] += -x138*x2;
M[211] += -x139*x23;
M[212] += -x140*x51;
M[213] += -x122*x93;
M[214] += -x127*x90;
M[215] += -x141*x49;
M[216] += -x142*x21;
M[217] += -x1*x143;
M[218] += -x115*x144;
M[219] += pow(x, 10)*x145;
M[220] += x1*x146;
M[221] += x146*x2;
M[222] += x147*x21;
M[223] += x116*x8;
M[224] += x147*x23;
M[225] += x148*x49;
M[226] += x117*x36;
M[227] += x117*x37;
M[228] += x148*x51;
M[229] += x149*x90;
M[230] += x118*x72;
M[231] += x150*x45*x9;
M[232] += x118*x74;
M[233] += x149*x93;
M[234] += x151*x41;
M[235] += x119*x124;
M[236] += x152*x97;
M[237] += x152*x98;
M[238] += x119*x126;
M[239] += x151*x44;
M[240] += x153*x58;
M[241] += x122*x18*x2;
M[242] += x101*x154;
M[243] += (1.0/864.0)*x130*x18;
M[244] += x131*x154;
M[245] += x1*x127*x18;
M[246] += x153*x64;
M[247] += x155*x79;
M[248] += x140*x95;
M[249] += x133*x156;
M[250] += x157*x25*x51;
M[251] += x157*x29*x49;
M[252] += x156*x21*x44;
M[253] += x141*x99;
M[254] += x155*x84;
M[255] += x105*x158;
M[256] += x139*x53;
M[257] += (1.0/2880.0)*x23*x3*x58;
M[258] += x100*x108;
M[259] += (1.0/1152.0)*x134*x25;
M[260] += x100*x110;
M[261] += x150*x3*x64;
M[262] += x142*x55;
M[263] += x113*x158;
M[264] += x137*x159;
M[265] += x138*x6;
M[266] += x139*x57;
M[267] += x104*x140;
M[268] += x0*x122*x29;
M[269] += x127*x136;
M[270] += x141*x78;
M[271] += x142*x40;
M[272] += x143*x5;
M[273] += x144*x159;
M[274] += x145*pow(y, 10);
M[275] += (1.0/362880.0)*x137*x2;
M[276] += (1.0/80640.0)*x105*x23;
M[277] += (1.0/30240.0)*x51*x79;
M[278] += (1.0/17280.0)*x58*x93;
M[279] += (1.0/14400.0)*q*x41*x44;
M[280] += (1.0/17280.0)*x64*x90;
M[281] += (1.0/30240.0)*x49*x84;
M[282] += (1.0/80640.0)*x113*x21;
M[283] += (1.0/362880.0)*x1*x144;
M[284] += x145*pow(z, 10);
}
void M2M_10(double x, double y, double z, double * M, double * Ms) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
double x448;
double x449;
double x450;
double x451;
double x452;
double x453;
double x454;
double x455;
double x456;
double x457;
double x458;
double x459;
double x460;
double x461;
double x462;
double x463;
double x464;
double x465;
double x466;
double x467;
double x468;
double x469;
double x470;
double x471;
double x472;
double x473;
double x474;
double x475;
double x476;
double x477;
double x478;
double x479;
double x480;
double x481;
double x482;
double x483;
double x484;
double x485;
double x486;
double x487;
double x488;
double x489;
double x490;
double x491;
double x492;
double x493;
double x494;
double x495;
double x496;
double x497;
double x498;
double x499;
double x500;
double x501;
double x502;
double x503;
double x504;
double x505;
double x506;
double x507;
double x508;
double x509;
double x510;
double x511;
double x512;
double x513;
double x514;
double x515;
double x516;
double x517;
double x518;
double x519;
double x520;
double x521;
double x522;
double x523;
double x524;
double x525;
double x526;
double x527;
double x528;
double x529;
double x530;
double x531;
double x532;
double x533;
double x534;
double x535;
double x536;
double x537;
double x538;
double x539;
double x540;
double x541;
double x542;
double x543;
double x544;
double x545;
double x546;
double x547;
double x548;
double x549;
double x550;
double x551;
double x552;
double x553;
double x554;
double x555;
double x556;
double x557;
double x558;
double x559;
double x560;
double x561;
double x562;
double x563;
double x564;
double x565;
double x566;
double x567;
double x568;
double x569;
double x570;
double x571;
double x572;
double x573;
double x574;
double x575;
double x576;
double x577;
double x578;
double x579;
double x580;
double x581;
double x582;
double x583;
double x584;
double x585;
double x586;
double x587;
double x588;
double x589;
double x590;
double x591;
double x592;
double x593;
double x594;
double x595;
double x596;
double x597;
double x598;
double x599;
double x600;
double x601;
double x602;
double x603;
double x604;
double x605;
double x606;
double x607;
double x608;
double x609;
double x610;
double x611;
double x612;
double x613;
double x614;
double x615;
double x616;
double x617;
double x618;
double x619;
double x620;
double x621;
double x622;
double x623;
double x624;
double x625;
double x626;
double x627;
double x628;
double x629;
double x630;
double x631;
double x632;
double x633;
double x634;
double x635;
double x636;
double x637;
double x638;
double x639;
double x640;
double x641;
double x642;
double x643;
double x644;
double x645;
double x646;
double x647;
double x648;
double x649;
double x650;
double x651;
double x652;
double x653;
double x654;
double x655;
double x656;
double x657;
double x658;
double x659;
double x660;
double x661;
double x662;
double x663;
double x664;
double x665;
double x666;
double x667;
double x668;
double x669;
double x670;
double x671;
double x672;
double x673;
double x674;
double x675;
double x676;
double x677;
double x678;
double x679;
double x680;
double x681;
double x682;
double x683;
double x684;
double x685;
double x686;
double x687;
double x688;
double x689;
double x690;
double x691;
double x692;
double x693;
double x694;
double x695;
double x696;
double x697;
double x698;
double x699;
double x700;
double x701;
double x702;
double x703;
double x704;
double x705;
double x706;
double x707;
double x708;
double x709;
double x710;
double x711;
double x712;
double x713;
double x714;
double x715;
double x716;
double x717;
double x718;
double x719;
double x720;
double x721;
double x722;
double x723;
double x724;
double x725;
double x726;
double x727;
double x728;
double x729;
double x730;
double x731;
double x732;
double x733;
double x734;
double x735;
double x736;
double x737;
double x738;
double x739;
double x740;
double x741;
double x742;
double x743;
double x744;
double x745;
double x746;
double x747;
double x748;
double x749;
double x750;
double x751;
double x752;
double x753;
double x754;
double x755;
double x756;
double x757;
double x758;
double x759;
double x760;
double x761;
double x762;
double x763;
double x764;
double x765;
double x766;
double x767;
double x768;
double x769;
double x770;
double x771;
double x772;
double x773;
double x774;
double x775;
double x776;
double x777;
double x778;
double x779;
double x780;
double x781;
double x782;
double x783;
double x784;
double x785;
double x786;
double x787;
double x788;
double x789;
double x790;
double x791;
double x792;
double x793;
double x794;
double x795;
double x796;
double x797;
double x798;
double x799;
double x800;
double x801;
double x802;
double x803;
double x804;
double x805;
double x806;
double x807;
double x808;
double x809;
double x810;
double x811;
double x812;
double x813;
double x814;
double x815;
double x816;
double x817;
double x818;
double x819;
double x820;
double x821;
double x822;
double x823;
double x824;
double x825;
double x826;
double x827;
double x828;
double x829;
double x830;
double x831;
double x832;
double x833;
double x834;
double x835;
double x836;
double x837;
double x838;
double x839;
double x840;
double x841;
double x842;
double x843;
double x844;
double x845;
double x846;
double x847;
double x848;
double x849;
double x850;
double x851;
double x852;
double x853;
double x854;
double x855;
double x856;
double x857;
double x858;
double x859;
double x860;
double x861;
double x862;
double x863;
double x864;
double x865;
double x866;
double x867;
double x868;
double x869;
double x870;
double x871;
double x872;
double x873;
double x874;
double x875;
double x876;
double x877;
double x878;
double x879;
double x880;
double x881;
double x882;
double x883;
double x884;
double x885;
double x886;
double x887;
double x888;
double x889;
double x890;
double x891;
double x892;
double x893;
double x894;
double x895;
double x896;
double x897;
double x898;
double x899;
double x900;
double x901;
double x902;
double x903;
double x904;
double x905;
double x906;
double x907;
double x908;
double x909;
double x910;
double x911;
double x912;
double x913;
double x914;
double x915;
double x916;
double x917;
double x918;
double x919;
double x920;
double x921;
double x922;
double x923;
double x924;
double x925;
double x926;
double x927;
double x928;
double x929;
double x930;
double x931;
double x932;
double x933;
double x934;
double x935;
double x936;
double x937;
double x938;
double x939;
double x940;
double x941;
double x942;
double x943;
double x944;
double x945;
double x946;
double x947;
double x948;
double x949;
double x950;
double x951;
double x952;
double x953;
x0 = x*M[0];
x1 = x*M[1];
x2 = y*M[0];
x3 = x*M[2];
x4 = z*M[0];
x5 = y*M[1];
x6 = y*M[2];
x7 = z*M[1];
x8 = z*M[2];
x9 = x*M[3];
x10 = (x*x);
x11 = (1.0/2.0)*x10;
x12 = x*M[4];
x13 = y*M[3];
x14 = x0*y;
x15 = x*M[5];
x16 = z*M[3];
x17 = x0*z;
x18 = x*M[6];
x19 = y*M[4];
x20 = x1*y;
x21 = (y*y);
x22 = (1.0/2.0)*M[0];
x23 = x*M[7];
x24 = y*M[5];
x25 = z*M[4];
x26 = x3*y;
x27 = x1*z;
x28 = x2*z;
x29 = x*M[8];
x30 = z*M[5];
x31 = x3*z;
x32 = (z*z);
x33 = y*M[6];
x34 = (1.0/2.0)*x21;
x35 = y*M[7];
x36 = z*M[6];
x37 = x5*z;
x38 = y*M[8];
x39 = z*M[7];
x40 = x6*z;
x41 = (1.0/2.0)*x32;
x42 = z*M[8];
x43 = x*M[9];
x44 = (x*x*x);
x45 = (1.0/6.0)*x44;
x46 = x*M[10];
x47 = y*M[9];
x48 = x9*y;
x49 = x*M[11];
x50 = z*M[9];
x51 = x9*z;
x52 = x*M[12];
x53 = y*M[10];
x54 = x12*y;
x55 = x*M[13];
x56 = y*M[11];
x57 = z*M[10];
x58 = x15*y;
x59 = x12*z;
x60 = x13*z;
x61 = x*M[14];
x62 = z*M[11];
x63 = x15*z;
x64 = x*M[15];
x65 = y*M[12];
x66 = x18*y;
x67 = (y*y*y);
x68 = (1.0/6.0)*M[0];
x69 = x*M[16];
x70 = y*M[13];
x71 = z*M[12];
x72 = x23*y;
x73 = x18*z;
x74 = x19*z;
x75 = x*M[17];
x76 = y*M[14];
x77 = z*M[13];
x78 = x29*y;
x79 = x23*z;
x80 = x24*z;
x81 = x*M[18];
x82 = z*M[14];
x83 = x29*z;
x84 = (z*z*z);
x85 = y*M[15];
x86 = (1.0/6.0)*x67;
x87 = y*M[16];
x88 = z*M[15];
x89 = x33*z;
x90 = y*M[17];
x91 = z*M[16];
x92 = x35*z;
x93 = y*M[18];
x94 = z*M[17];
x95 = x38*z;
x96 = (1.0/6.0)*x84;
x97 = z*M[18];
x98 = x*M[19];
x99 = (x*x*x*x);
x100 = (1.0/24.0)*x99;
x101 = x*M[20];
x102 = y*M[19];
x103 = x43*y;
x104 = x*M[21];
x105 = z*M[19];
x106 = x43*z;
x107 = x*M[22];
x108 = y*M[20];
x109 = x46*y;
x110 = (1.0/4.0)*x10;
x111 = x21*M[0];
x112 = x*M[23];
x113 = y*M[21];
x114 = z*M[20];
x115 = x49*y;
x116 = x46*z;
x117 = x47*z;
x118 = x*M[24];
x119 = z*M[21];
x120 = x49*z;
x121 = x110*x32;
x122 = x*M[25];
x123 = y*M[22];
x124 = x52*y;
x125 = x110*x21;
x126 = x*M[26];
x127 = y*M[23];
x128 = z*M[22];
x129 = x55*y;
x130 = x52*z;
x131 = x53*z;
x132 = x*M[27];
x133 = y*M[24];
x134 = z*M[23];
x135 = x61*y;
x136 = x55*z;
x137 = x56*z;
x138 = x*M[28];
x139 = z*M[24];
x140 = x61*z;
x141 = x*M[29];
x142 = y*M[25];
x143 = x64*y;
x144 = (y*y*y*y);
x145 = (1.0/24.0)*M[0];
x146 = x*M[30];
x147 = y*M[26];
x148 = z*M[25];
x149 = x69*y;
x150 = x64*z;
x151 = x65*z;
x152 = x*M[31];
x153 = y*M[27];
x154 = z*M[26];
x155 = x75*y;
x156 = x69*z;
x157 = x70*z;
x158 = (1.0/4.0)*x32;
x159 = x*M[32];
x160 = y*M[28];
x161 = z*M[27];
x162 = x81*y;
x163 = x75*z;
x164 = x76*z;
x165 = x*M[33];
x166 = z*M[28];
x167 = x81*z;
x168 = (z*z*z*z);
x169 = y*M[29];
x170 = (1.0/24.0)*x144;
x171 = y*M[30];
x172 = z*M[29];
x173 = x85*z;
x174 = y*M[31];
x175 = z*M[30];
x176 = x87*z;
x177 = x158*x21;
x178 = y*M[32];
x179 = z*M[31];
x180 = x90*z;
x181 = y*M[33];
x182 = z*M[32];
x183 = x93*z;
x184 = (1.0/24.0)*x168;
x185 = z*M[33];
x186 = x*M[34];
x187 = pow(x, 5);
x188 = (1.0/120.0)*x187;
x189 = x*M[35];
x190 = y*M[34];
x191 = x98*y;
x192 = x*M[36];
x193 = z*M[34];
x194 = x98*z;
x195 = x*M[37];
x196 = y*M[35];
x197 = x101*y;
x198 = (1.0/12.0)*x44;
x199 = x*M[38];
x200 = y*M[36];
x201 = z*M[35];
x202 = x104*y;
x203 = x101*z;
x204 = x102*z;
x205 = x*M[39];
x206 = z*M[36];
x207 = x104*z;
x208 = x198*x32;
x209 = x*M[40];
x210 = y*M[37];
x211 = x107*y;
x212 = (1.0/12.0)*x10;
x213 = x67*M[0];
x214 = x198*x21;
x215 = x*M[41];
x216 = y*M[38];
x217 = z*M[37];
x218 = x112*y;
x219 = x107*z;
x220 = x108*z;
x221 = x*M[42];
x222 = y*M[39];
x223 = z*M[38];
x224 = x118*y;
x225 = x112*z;
x226 = x113*z;
x227 = x*M[43];
x228 = z*M[39];
x229 = x118*z;
x230 = x212*x84;
x231 = x*M[44];
x232 = y*M[40];
x233 = x122*y;
x234 = x212*x67;
x235 = x*M[45];
x236 = y*M[41];
x237 = z*M[40];
x238 = x126*y;
x239 = x122*z;
x240 = x123*z;
x241 = x*M[46];
x242 = y*M[42];
x243 = z*M[41];
x244 = x132*y;
x245 = x126*z;
x246 = x127*z;
x247 = x*M[47];
x248 = y*M[43];
x249 = z*M[42];
x250 = x138*y;
x251 = x132*z;
x252 = x133*z;
x253 = x*M[48];
x254 = z*M[43];
x255 = x138*z;
x256 = x*M[49];
x257 = y*M[44];
x258 = x141*y;
x259 = pow(y, 5);
x260 = (1.0/120.0)*M[0];
x261 = x*M[50];
x262 = y*M[45];
x263 = z*M[44];
x264 = x146*y;
x265 = x141*z;
x266 = x142*z;
x267 = x*M[51];
x268 = y*M[46];
x269 = z*M[45];
x270 = x152*y;
x271 = x146*z;
x272 = x147*z;
x273 = (1.0/12.0)*x32;
x274 = x*M[52];
x275 = y*M[47];
x276 = z*M[46];
x277 = x159*y;
x278 = x152*z;
x279 = x153*z;
x280 = (1.0/12.0)*x84;
x281 = x*M[53];
x282 = y*M[48];
x283 = z*M[47];
x284 = x165*y;
x285 = x159*z;
x286 = x160*z;
x287 = x*M[54];
x288 = z*M[48];
x289 = x165*z;
x290 = pow(z, 5);
x291 = y*M[49];
x292 = (1.0/120.0)*x259;
x293 = y*M[50];
x294 = z*M[49];
x295 = x169*z;
x296 = y*M[51];
x297 = z*M[50];
x298 = x171*z;
x299 = x273*x67;
x300 = y*M[52];
x301 = z*M[51];
x302 = x174*z;
x303 = x21*x280;
x304 = y*M[53];
x305 = z*M[52];
x306 = x178*z;
x307 = y*M[54];
x308 = z*M[53];
x309 = x181*z;
x310 = (1.0/120.0)*x290;
x311 = z*M[54];
x312 = x*M[55];
x313 = pow(x, 6);
x314 = (1.0/720.0)*x313;
x315 = x*M[56];
x316 = y*M[55];
x317 = x186*y;
x318 = x*M[57];
x319 = z*M[55];
x320 = x186*z;
x321 = x*M[58];
x322 = y*M[56];
x323 = x189*y;
x324 = (1.0/48.0)*x99;
x325 = x*M[59];
x326 = y*M[57];
x327 = z*M[56];
x328 = x192*y;
x329 = x189*z;
x330 = x190*z;
x331 = x*M[60];
x332 = z*M[57];
x333 = x192*z;
x334 = x32*x324;
x335 = x*M[61];
x336 = y*M[58];
x337 = x195*y;
x338 = (1.0/36.0)*x44;
x339 = x21*x324;
x340 = x*M[62];
x341 = y*M[59];
x342 = z*M[58];
x343 = x199*y;
x344 = x195*z;
x345 = x196*z;
x346 = x*M[63];
x347 = y*M[60];
x348 = z*M[59];
x349 = x205*y;
x350 = x199*z;
x351 = x200*z;
x352 = x*M[64];
x353 = z*M[60];
x354 = x205*z;
x355 = x338*x84;
x356 = x*M[65];
x357 = y*M[61];
x358 = x209*y;
x359 = (1.0/48.0)*x10;
x360 = x144*M[0];
x361 = x338*x67;
x362 = x*M[66];
x363 = y*M[62];
x364 = z*M[61];
x365 = x215*y;
x366 = x209*z;
x367 = x210*z;
x368 = x*M[67];
x369 = y*M[63];
x370 = z*M[62];
x371 = x221*y;
x372 = x215*z;
x373 = x216*z;
x374 = x10*x32;
x375 = (1.0/8.0)*x374;
x376 = x*M[68];
x377 = y*M[64];
x378 = z*M[63];
x379 = x227*y;
x380 = x221*z;
x381 = x222*z;
x382 = x*M[69];
x383 = z*M[64];
x384 = x227*z;
x385 = x168*x359;
x386 = x*M[70];
x387 = y*M[65];
x388 = x231*y;
x389 = x144*x359;
x390 = x*M[71];
x391 = y*M[66];
x392 = z*M[65];
x393 = x235*y;
x394 = x231*z;
x395 = x232*z;
x396 = x*M[72];
x397 = y*M[67];
x398 = z*M[66];
x399 = x241*y;
x400 = x235*z;
x401 = x236*z;
x402 = x21*x375;
x403 = x*M[73];
x404 = y*M[68];
x405 = z*M[67];
x406 = x247*y;
x407 = x241*z;
x408 = x242*z;
x409 = x*M[74];
x410 = y*M[69];
x411 = z*M[68];
x412 = x253*y;
x413 = x247*z;
x414 = x248*z;
x415 = x*M[75];
x416 = z*M[69];
x417 = x253*z;
x418 = x*M[76];
x419 = y*M[70];
x420 = x256*y;
x421 = pow(y, 6);
x422 = (1.0/720.0)*M[0];
x423 = x*M[77];
x424 = y*M[71];
x425 = z*M[70];
x426 = x261*y;
x427 = x256*z;
x428 = x257*z;
x429 = x*M[78];
x430 = y*M[72];
x431 = z*M[71];
x432 = x267*y;
x433 = x261*z;
x434 = x262*z;
x435 = (1.0/48.0)*x32;
x436 = x*M[79];
x437 = y*M[73];
x438 = z*M[72];
x439 = x274*y;
x440 = x267*z;
x441 = x268*z;
x442 = (1.0/36.0)*x84;
x443 = x*M[80];
x444 = y*M[74];
x445 = z*M[73];
x446 = x281*y;
x447 = x274*z;
x448 = x275*z;
x449 = (1.0/48.0)*x168;
x450 = x*M[81];
x451 = y*M[75];
x452 = z*M[74];
x453 = x287*y;
x454 = x281*z;
x455 = x282*z;
x456 = x*M[82];
x457 = z*M[75];
x458 = x287*z;
x459 = pow(z, 6);
x460 = y*M[76];
x461 = (1.0/720.0)*x421;
x462 = y*M[77];
x463 = z*M[76];
x464 = x291*z;
x465 = y*M[78];
x466 = z*M[77];
x467 = x293*z;
x468 = x144*x435;
x469 = y*M[79];
x470 = z*M[78];
x471 = x296*z;
x472 = x442*x67;
x473 = y*M[80];
x474 = z*M[79];
x475 = x300*z;
x476 = x21*x449;
x477 = y*M[81];
x478 = z*M[80];
x479 = x304*z;
x480 = y*M[82];
x481 = z*M[81];
x482 = x307*z;
x483 = (1.0/720.0)*x459;
x484 = z*M[82];
x485 = x*M[83];
x486 = pow(x, 7);
x487 = (1.0/5040.0)*x486;
x488 = x*M[84];
x489 = y*M[83];
x490 = x312*y;
x491 = x*M[85];
x492 = z*M[83];
x493 = x312*z;
x494 = x*M[86];
x495 = y*M[84];
x496 = x315*y;
x497 = (1.0/240.0)*x187;
x498 = x*M[87];
x499 = y*M[85];
x500 = z*M[84];
x501 = x318*y;
x502 = x315*z;
x503 = x316*z;
x504 = x*M[88];
x505 = z*M[85];
x506 = x318*z;
x507 = x32*x497;
x508 = x*M[89];
x509 = y*M[86];
x510 = x321*y;
x511 = (1.0/144.0)*x99;
x512 = x21*x497;
x513 = x*M[90];
x514 = y*M[87];
x515 = z*M[86];
x516 = x325*y;
x517 = x321*z;
x518 = x322*z;
x519 = x*M[91];
x520 = y*M[88];
x521 = z*M[87];
x522 = x331*y;
x523 = x325*z;
x524 = x326*z;
x525 = x*M[92];
x526 = z*M[88];
x527 = x331*z;
x528 = x511*x84;
x529 = x*M[93];
x530 = y*M[89];
x531 = x335*y;
x532 = (1.0/144.0)*x44;
x533 = x511*x67;
x534 = x*M[94];
x535 = y*M[90];
x536 = z*M[89];
x537 = x340*y;
x538 = x335*z;
x539 = x336*z;
x540 = x*M[95];
x541 = y*M[91];
x542 = z*M[90];
x543 = x346*y;
x544 = x340*z;
x545 = x341*z;
x546 = x145*x21;
x547 = x32*x44;
x548 = x*M[96];
x549 = y*M[92];
x550 = z*M[91];
x551 = x352*y;
x552 = x346*z;
x553 = x347*z;
x554 = x*M[97];
x555 = z*M[92];
x556 = x352*z;
x557 = x168*x532;
x558 = x*M[98];
x559 = y*M[93];
x560 = x356*y;
x561 = (1.0/240.0)*x10;
x562 = x259*M[0];
x563 = x144*x532;
x564 = x*M[99];
x565 = y*M[94];
x566 = z*M[93];
x567 = x362*y;
x568 = x356*z;
x569 = x357*z;
x570 = x*M[100];
x571 = y*M[95];
x572 = z*M[94];
x573 = x368*y;
x574 = x362*z;
x575 = x363*z;
x576 = x374*x67;
x577 = (1.0/24.0)*x21;
x578 = x547*M[1];
x579 = x*M[101];
x580 = y*M[96];
x581 = z*M[95];
x582 = x376*y;
x583 = x368*z;
x584 = x369*z;
x585 = x10*x84;
x586 = x577*M[2];
x587 = x*M[102];
x588 = y*M[97];
x589 = z*M[96];
x590 = x382*y;
x591 = x376*z;
x592 = x377*z;
x593 = x*M[103];
x594 = z*M[97];
x595 = x382*z;
x596 = x290*x561;
x597 = x*M[104];
x598 = y*M[98];
x599 = x386*y;
x600 = x259*x561;
x601 = x*M[105];
x602 = y*M[99];
x603 = z*M[98];
x604 = x390*y;
x605 = x386*z;
x606 = x387*z;
x607 = x*M[106];
x608 = y*M[100];
x609 = z*M[99];
x610 = x396*y;
x611 = x390*z;
x612 = x391*z;
x613 = x67*M[1];
x614 = (1.0/24.0)*x374;
x615 = x*M[107];
x616 = y*M[101];
x617 = z*M[100];
x618 = x403*y;
x619 = x396*z;
x620 = x397*z;
x621 = x84*M[1];
x622 = x10*x621;
x623 = x67*M[2];
x624 = x*M[108];
x625 = y*M[102];
x626 = z*M[101];
x627 = x409*y;
x628 = x403*z;
x629 = x404*z;
x630 = x*M[109];
x631 = y*M[103];
x632 = z*M[102];
x633 = x415*y;
x634 = x409*z;
x635 = x410*z;
x636 = x*M[110];
x637 = z*M[103];
x638 = x415*z;
x639 = x*M[111];
x640 = y*M[104];
x641 = x418*y;
x642 = pow(y, 7);
x643 = (1.0/5040.0)*M[0];
x644 = x*M[112];
x645 = y*M[105];
x646 = z*M[104];
x647 = x423*y;
x648 = x418*z;
x649 = x419*z;
x650 = x*M[113];
x651 = y*M[106];
x652 = z*M[105];
x653 = x429*y;
x654 = x423*z;
x655 = x424*z;
x656 = (1.0/240.0)*x32;
x657 = x*M[114];
x658 = y*M[107];
x659 = z*M[106];
x660 = x436*y;
x661 = x429*z;
x662 = x430*z;
x663 = (1.0/144.0)*x84;
x664 = x*M[115];
x665 = y*M[108];
x666 = z*M[107];
x667 = x443*y;
x668 = x436*z;
x669 = x437*z;
x670 = (1.0/144.0)*x168;
x671 = x*M[116];
x672 = y*M[109];
x673 = z*M[108];
x674 = x450*y;
x675 = x443*z;
x676 = x444*z;
x677 = (1.0/240.0)*x290;
x678 = x*M[117];
x679 = y*M[110];
x680 = z*M[109];
x681 = x456*y;
x682 = x450*z;
x683 = x451*z;
x684 = x*M[118];
x685 = z*M[110];
x686 = x456*z;
x687 = pow(z, 7);
x688 = y*M[111];
x689 = (1.0/5040.0)*x642;
x690 = y*M[112];
x691 = z*M[111];
x692 = x460*z;
x693 = y*M[113];
x694 = z*M[112];
x695 = x462*z;
x696 = x259*x656;
x697 = y*M[114];
x698 = z*M[113];
x699 = x465*z;
x700 = x144*x663;
x701 = y*M[115];
x702 = z*M[114];
x703 = x469*z;
x704 = x67*x670;
x705 = y*M[116];
x706 = z*M[115];
x707 = x473*z;
x708 = x21*x677;
x709 = y*M[117];
x710 = z*M[116];
x711 = x477*z;
x712 = y*M[118];
x713 = z*M[117];
x714 = x480*z;
x715 = (1.0/5040.0)*x687;
x716 = z*M[118];
x717 = x*M[119];
x718 = (1.0/40320.0)*pow(x, 8);
x719 = x*M[120];
x720 = y*M[119];
x721 = x485*y;
x722 = x*M[121];
x723 = x*M[122];
x724 = y*M[120];
x725 = x488*y;
x726 = (1.0/1440.0)*x313;
x727 = x*M[123];
x728 = y*M[121];
x729 = x491*y;
x730 = x*M[124];
x731 = x32*x726;
x732 = x*M[125];
x733 = y*M[122];
x734 = x494*y;
x735 = x187*x422;
x736 = x21*x726;
x737 = x*M[126];
x738 = y*M[123];
x739 = x498*y;
x740 = x*M[127];
x741 = y*M[124];
x742 = x504*y;
x743 = x*M[128];
x744 = x*M[129];
x745 = y*M[125];
x746 = x508*y;
x747 = (1.0/576.0)*x99;
x748 = (1.0/720.0)*x187;
x749 = x*M[130];
x750 = y*M[126];
x751 = x513*y;
x752 = x748*M[2];
x753 = x*M[131];
x754 = y*M[127];
x755 = x519*y;
x756 = x547*x577;
x757 = (1.0/96.0)*x111;
x758 = x32*x99;
x759 = x*M[132];
x760 = y*M[128];
x761 = x525*y;
x762 = x*M[133];
x763 = x168*x747;
x764 = x*M[134];
x765 = y*M[129];
x766 = x529*y;
x767 = x422*x44;
x768 = x144*x747;
x769 = x*M[135];
x770 = y*M[130];
x771 = x534*y;
x772 = x*M[136];
x773 = y*M[131];
x774 = x540*y;
x775 = (1.0/24.0)*x576;
x776 = (1.0/72.0)*x547;
x777 = (1.0/96.0)*M[1];
x778 = x21*x758;
x779 = x*M[137];
x780 = y*M[132];
x781 = x548*y;
x782 = x577*x585;
x783 = (1.0/72.0)*x44;
x784 = x111*x84;
x785 = (1.0/96.0)*M[2];
x786 = x*M[138];
x787 = y*M[133];
x788 = x554*y;
x789 = x*M[139];
x790 = x*M[140];
x791 = y*M[134];
x792 = x558*y;
x793 = (1.0/1440.0)*x10;
x794 = x421*M[0];
x795 = (1.0/720.0)*x44;
x796 = x259*M[1];
x797 = x*M[141];
x798 = y*M[135];
x799 = x564*y;
x800 = x795*M[2];
x801 = x*M[142];
x802 = y*M[136];
x803 = x570*y;
x804 = (1.0/96.0)*x374;
x805 = x*M[143];
x806 = y*M[137];
x807 = x579*y;
x808 = (1.0/72.0)*x585;
x809 = x21*x783;
x810 = x*M[144];
x811 = y*M[138];
x812 = x587*y;
x813 = x10*x168;
x814 = x809*x84;
x815 = x*M[145];
x816 = y*M[139];
x817 = x593*y;
x818 = x290*x795;
x819 = x*M[146];
x820 = x459*x793;
x821 = x*M[147];
x822 = y*M[140];
x823 = x597*y;
x824 = x421*x793;
x825 = x*M[148];
x826 = y*M[141];
x827 = x601*y;
x828 = x*M[149];
x829 = y*M[142];
x830 = x607*y;
x831 = x144*x374;
x832 = x*M[150];
x833 = y*M[143];
x834 = x615*y;
x835 = x*M[151];
x836 = y*M[144];
x837 = x624*y;
x838 = x21*x813;
x839 = x*M[152];
x840 = y*M[145];
x841 = x630*y;
x842 = x*M[153];
x843 = y*M[146];
x844 = x636*y;
x845 = x*M[154];
x846 = x*M[155];
x847 = y*M[147];
x848 = x639*y;
x849 = pow(y, 8);
x850 = (1.0/40320.0)*M[0];
x851 = x*M[156];
x852 = y*M[148];
x853 = x644*y;
x854 = x*M[157];
x855 = y*M[149];
x856 = x650*y;
x857 = (1.0/1440.0)*x32;
x858 = x*M[158];
x859 = y*M[150];
x860 = x657*y;
x861 = x259*x84;
x862 = x*M[159];
x863 = y*M[151];
x864 = x664*y;
x865 = (1.0/576.0)*x168;
x866 = x*M[160];
x867 = y*M[152];
x868 = x671*y;
x869 = x290*x67;
x870 = x*M[161];
x871 = y*M[153];
x872 = x678*y;
x873 = (1.0/1440.0)*x459;
x874 = x*M[162];
x875 = y*M[154];
x876 = x684*y;
x877 = x*M[163];
x878 = pow(z, 8);
x879 = y*M[155];
x880 = (1.0/40320.0)*x849;
x881 = y*M[156];
x882 = y*M[157];
x883 = x421*x857;
x884 = y*M[158];
x885 = y*M[159];
x886 = x144*x865;
x887 = (1.0/720.0)*M[2];
x888 = y*M[160];
x889 = y*M[161];
x890 = x21*x873;
x891 = y*M[162];
x892 = y*M[163];
x893 = (1.0/40320.0)*x878;
x894 = (1.0/362880.0)*pow(x, 9);
x895 = (1.0/10080.0)*x486;
x896 = x32*x895;
x897 = x748*M[3];
x898 = (1.0/4320.0)*x313;
x899 = x21*x895;
x900 = x84*x898;
x901 = x67*x748;
x902 = (1.0/2880.0)*x187;
x903 = (1.0/96.0)*x778;
x904 = (1.0/480.0)*x111;
x905 = x187*x32;
x906 = x748*x84;
x907 = x168*x902;
x908 = x259*x795;
x909 = (1.0/2880.0)*x99;
x910 = x144*x902;
x911 = x67*x776;
x912 = (1.0/288.0)*x758;
x913 = x21*M[1];
x914 = (1.0/480.0)*x905;
x915 = (1.0/288.0)*x99;
x916 = x21*M[2];
x917 = x290*x909;
x918 = (1.0/4320.0)*x44;
x919 = x259*M[2];
x920 = x144*x804;
x921 = (1.0/288.0)*x360;
x922 = x67*x808;
x923 = (1.0/216.0)*x44*x84;
x924 = (1.0/96.0)*x838;
x925 = (1.0/288.0)*x168*x44;
x926 = x459*x918;
x927 = (1.0/10080.0)*x10;
x928 = x927*M[0];
x929 = x421*x918;
x930 = (1.0/480.0)*x374;
x931 = (1.0/288.0)*x144;
x932 = x931*M[2];
x933 = (1.0/288.0)*x813;
x934 = x10*x290;
x935 = x642*x927;
x936 = (1.0/720.0)*M[3];
x937 = (1.0/720.0)*x861;
x938 = (1.0/720.0)*x869;
x939 = (1.0/480.0)*x934;
x940 = x687*x927;
x941 = pow(y, 9);
x942 = (1.0/362880.0)*M[0];
x943 = (1.0/10080.0)*x32*x642;
x944 = (1.0/4320.0)*x84;
x945 = (1.0/2880.0)*x168;
x946 = (1.0/2880.0)*x290;
x947 = (1.0/4320.0)*x459;
x948 = (1.0/10080.0)*x687;
x949 = pow(z, 9);
x950 = (1.0/362880.0)*x941;
x951 = x144*x946;
x952 = x21*x948;
x953 = (1.0/362880.0)*x949;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += M[1];
#pragma omp atomic
Ms[2] += M[2];
#pragma omp atomic
Ms[3] += x0 + M[3];
#pragma omp atomic
Ms[4] += x1 + x2 + M[4];
#pragma omp atomic
Ms[5] += x3 + x4 + M[5];
#pragma omp atomic
Ms[6] += x5 + M[6];
#pragma omp atomic
Ms[7] += x6 + x7 + M[7];
#pragma omp atomic
Ms[8] += x8 + M[8];
#pragma omp atomic
Ms[9] += x11*M[0] + x9 + M[9];
#pragma omp atomic
Ms[10] += x11*M[1] + x12 + x13 + x14 + M[10];
#pragma omp atomic
Ms[11] += x11*M[2] + x15 + x16 + x17 + M[11];
#pragma omp atomic
Ms[12] += x18 + x19 + x20 + x21*x22 + M[12];
#pragma omp atomic
Ms[13] += x23 + x24 + x25 + x26 + x27 + x28 + M[13];
#pragma omp atomic
Ms[14] += x22*x32 + x29 + x30 + x31 + M[14];
#pragma omp atomic
Ms[15] += x33 + x34*M[1] + M[15];
#pragma omp atomic
Ms[16] += x34*M[2] + x35 + x36 + x37 + M[16];
#pragma omp atomic
Ms[17] += x38 + x39 + x40 + x41*M[1] + M[17];
#pragma omp atomic
Ms[18] += x41*M[2] + x42 + M[18];
#pragma omp atomic
Ms[19] += x11*M[3] + x43 + x45*M[0] + M[19];
#pragma omp atomic
Ms[20] += x11*x2 + x11*M[4] + x45*M[1] + x46 + x47 + x48 + M[20];
#pragma omp atomic
Ms[21] += x11*x4 + x11*M[5] + x45*M[2] + x49 + x50 + x51 + M[21];
#pragma omp atomic
Ms[22] += x0*x34 + x11*x5 + x11*M[6] + x34*M[3] + x52 + x53 + x54 + M[22];
#pragma omp atomic
Ms[23] += x11*x6 + x11*x7 + x11*M[7] + x14*z + x55 + x56 + x57 + x58 + x59 + x60 + M[23];
#pragma omp atomic
Ms[24] += x0*x41 + x11*x8 + x11*M[8] + x41*M[3] + x61 + x62 + x63 + M[24];
#pragma omp atomic
Ms[25] += x1*x34 + x34*M[4] + x64 + x65 + x66 + x67*x68 + M[25];
#pragma omp atomic
Ms[26] += x20*z + x3*x34 + x34*x4 + x34*M[5] + x69 + x70 + x71 + x72 + x73 + x74 + M[26];
#pragma omp atomic
Ms[27] += x1*x41 + x2*x41 + x26*z + x41*M[4] + x75 + x76 + x77 + x78 + x79 + x80 + M[27];
#pragma omp atomic
Ms[28] += x3*x41 + x41*M[5] + x68*x84 + x81 + x82 + x83 + M[28];
#pragma omp atomic
Ms[29] += x34*M[6] + x85 + x86*M[1] + M[29];
#pragma omp atomic
Ms[30] += x34*x7 + x34*M[7] + x86*M[2] + x87 + x88 + x89 + M[30];
#pragma omp atomic
Ms[31] += x34*x8 + x34*M[8] + x41*x5 + x41*M[6] + x90 + x91 + x92 + M[31];
#pragma omp atomic
Ms[32] += x41*x6 + x41*M[7] + x93 + x94 + x95 + x96*M[1] + M[32];
#pragma omp atomic
Ms[33] += x41*M[8] + x96*M[2] + x97 + M[33];
#pragma omp atomic
Ms[34] += x100*M[0] + x11*M[9] + x45*M[3] + x98 + M[34];
#pragma omp atomic
Ms[35] += x100*M[1] + x101 + x102 + x103 + x11*x13 + x11*M[10] + x2*x45 + x45*M[4] + M[35];
#pragma omp atomic
Ms[36] += x100*M[2] + x104 + x105 + x106 + x11*x16 + x11*M[11] + x4*x45 + x45*M[5] + M[36];
#pragma omp atomic
Ms[37] += x107 + x108 + x109 + x11*x19 + x11*M[12] + x110*x111 + x34*x9 + x34*M[9] + x45*x5 + x45*M[6] + M[37];
#pragma omp atomic
Ms[38] += x11*x24 + x11*x25 + x11*x28 + x11*M[13] + x112 + x113 + x114 + x115 + x116 + x117 + x45*x6 + x45*x7 + x45*M[7] + x48*z + M[38];
#pragma omp atomic
Ms[39] += x11*x30 + x11*M[14] + x118 + x119 + x120 + x121*M[0] + x41*x9 + x41*M[9] + x45*x8 + x45*M[8] + M[39];
#pragma omp atomic
Ms[40] += x0*x86 + x11*x33 + x11*M[15] + x12*x34 + x122 + x123 + x124 + x125*M[1] + x34*M[10] + x86*M[3] + M[40];
#pragma omp atomic
Ms[41] += x11*x35 + x11*x36 + x11*x37 + x11*M[16] + x125*M[2] + x126 + x127 + x128 + x129 + x130 + x131 + x15*x34 + x16*x34 + x17*x34 + x34*M[11] + x54*z + M[41];
#pragma omp atomic
Ms[42] += x11*x38 + x11*x39 + x11*x40 + x11*M[17] + x12*x41 + x121*M[1] + x13*x41 + x132 + x133 + x134 + x135 + x136 + x137 + x14*x41 + x41*M[10] + x58*z + M[42];
#pragma omp atomic
Ms[43] += x0*x96 + x11*x42 + x11*M[18] + x121*M[2] + x138 + x139 + x140 + x15*x41 + x41*M[11] + x96*M[3] + M[43];
#pragma omp atomic
Ms[44] += x1*x86 + x141 + x142 + x143 + x144*x145 + x18*x34 + x34*M[12] + x86*M[4] + M[44];
#pragma omp atomic
Ms[45] += x146 + x147 + x148 + x149 + x150 + x151 + x23*x34 + x25*x34 + x27*x34 + x3*x86 + x34*M[13] + x4*x86 + x66*z + x86*M[5] + M[45];
#pragma omp atomic
Ms[46] += x111*x158 + x152 + x153 + x154 + x155 + x156 + x157 + x18*x41 + x19*x41 + x20*x41 + x29*x34 + x30*x34 + x31*x34 + x34*M[14] + x41*M[12] + x72*z + M[46];
#pragma omp atomic
Ms[47] += x1*x96 + x159 + x160 + x161 + x162 + x163 + x164 + x2*x96 + x23*x41 + x24*x41 + x26*x41 + x41*M[13] + x78*z + x96*M[4] + M[47];
#pragma omp atomic
Ms[48] += x145*x168 + x165 + x166 + x167 + x29*x41 + x3*x96 + x41*M[14] + x96*M[5] + M[48];
#pragma omp atomic
Ms[49] += x169 + x170*M[1] + x34*M[15] + x86*M[6] + M[49];
#pragma omp atomic
Ms[50] += x170*M[2] + x171 + x172 + x173 + x34*x36 + x34*M[16] + x7*x86 + x86*M[7] + M[50];
#pragma omp atomic
Ms[51] += x174 + x175 + x176 + x177*M[1] + x33*x41 + x34*x39 + x34*M[17] + x41*M[15] + x8*x86 + x86*M[8] + M[51];
#pragma omp atomic
Ms[52] += x177*M[2] + x178 + x179 + x180 + x34*x42 + x34*M[18] + x35*x41 + x41*M[16] + x5*x96 + x96*M[6] + M[52];
#pragma omp atomic
Ms[53] += x181 + x182 + x183 + x184*M[1] + x38*x41 + x41*M[17] + x6*x96 + x96*M[7] + M[53];
#pragma omp atomic
Ms[54] += x184*M[2] + x185 + x41*M[18] + x96*M[8] + M[54];
#pragma omp atomic
Ms[55] += x100*M[3] + x11*M[19] + x186 + x188*M[0] + x45*M[9] + M[55];
#pragma omp atomic
Ms[56] += x100*x2 + x100*M[4] + x11*x47 + x11*M[20] + x13*x45 + x188*M[1] + x189 + x190 + x191 + x45*M[10] + M[56];
#pragma omp atomic
Ms[57] += x100*x4 + x100*M[5] + x11*x50 + x11*M[21] + x16*x45 + x188*M[2] + x192 + x193 + x194 + x45*M[11] + M[57];
#pragma omp atomic
Ms[58] += x100*x5 + x100*M[6] + x11*x53 + x11*M[22] + x111*x198 + x125*M[3] + x19*x45 + x195 + x196 + x197 + x34*x43 + x34*M[19] + x45*M[12] + M[58];
#pragma omp atomic
Ms[59] += x100*x6 + x100*x7 + x100*M[7] + x103*z + x11*x56 + x11*x57 + x11*x60 + x11*M[23] + x199 + x200 + x201 + x202 + x203 + x204 + x24*x45 + x25*x45 + x28*x45 + x45*M[13] + M[59];
#pragma omp atomic
Ms[60] += x100*x8 + x100*M[8] + x11*x62 + x11*M[24] + x121*M[3] + x205 + x206 + x207 + x208*M[0] + x30*x45 + x41*x43 + x41*M[19] + x45*M[14] + M[60];
#pragma omp atomic
Ms[61] += x11*x65 + x11*M[25] + x125*M[4] + x209 + x210 + x211 + x212*x213 + x214*M[1] + x33*x45 + x34*x46 + x34*M[20] + x45*M[15] + x86*x9 + x86*M[9] + M[61];
#pragma omp atomic
Ms[62] += x109*z + x11*x70 + x11*x71 + x11*x74 + x11*M[26] + x125*x4 + x125*M[5] + x214*M[2] + x215 + x216 + x217 + x218 + x219 + x220 + x34*x49 + x34*x50 + x34*x51 + x34*M[21] + x35*x45 + x36*x45 + x37*x45 + x45*M[16] + M[62];
#pragma omp atomic
Ms[63] += x11*x76 + x11*x77 + x11*x80 + x11*M[27] + x115*z + x121*x2 + x121*M[4] + x208*M[1] + x221 + x222 + x223 + x224 + x225 + x226 + x38*x45 + x39*x45 + x40*x45 + x41*x46 + x41*x47 + x41*x48 + x41*M[20] + x45*M[17] + M[63];
#pragma omp atomic
Ms[64] += x11*x82 + x11*M[28] + x121*M[5] + x208*M[2] + x227 + x228 + x229 + x230*M[0] + x41*x49 + x41*M[21] + x42*x45 + x45*M[18] + x9*x96 + x96*M[9] + M[64];
#pragma omp atomic
Ms[65] += x0*x170 + x11*x85 + x11*M[29] + x12*x86 + x125*M[6] + x170*M[3] + x231 + x232 + x233 + x234*M[1] + x34*x52 + x34*M[22] + x86*M[10] + M[65];
#pragma omp atomic
Ms[66] += x11*x87 + x11*x88 + x11*x89 + x11*M[30] + x124*z + x125*x7 + x125*M[7] + x15*x86 + x16*x86 + x17*x86 + x234*M[2] + x235 + x236 + x237 + x238 + x239 + x240 + x34*x55 + x34*x57 + x34*x59 + x34*M[23] + x86*M[11] + M[66];
#pragma omp atomic
Ms[67] += x0*x177 + x11*x90 + x11*x91 + x11*x92 + x11*M[31] + x121*x5 + x121*M[6] + x125*x8 + x125*M[8] + x129*z + x177*M[3] + x241 + x242 + x243 + x244 + x245 + x246 + x34*x61 + x34*x62 + x34*x63 + x34*M[24] + x41*x52 + x41*x53 + x41*x54 + x41*M[22] + M[67];
#pragma omp atomic
Ms[68] += x11*x93 + x11*x94 + x11*x95 + x11*M[32] + x12*x96 + x121*x6 + x121*M[7] + x13*x96 + x135*z + x14*x96 + x230*M[1] + x247 + x248 + x249 + x250 + x251 + x252 + x41*x55 + x41*x56 + x41*x58 + x41*M[23] + x96*M[10] + M[68];
#pragma omp atomic
Ms[69] += x0*x184 + x11*x97 + x11*M[33] + x121*M[8] + x15*x96 + x184*M[3] + x230*M[2] + x253 + x254 + x255 + x41*x61 + x41*M[24] + x96*M[11] + M[69];
#pragma omp atomic
Ms[70] += x1*x170 + x170*M[4] + x18*x86 + x256 + x257 + x258 + x259*x260 + x34*x64 + x34*M[25] + x86*M[12] + M[70];
#pragma omp atomic
Ms[71] += x143*z + x170*x3 + x170*x4 + x170*M[5] + x23*x86 + x25*x86 + x261 + x262 + x263 + x264 + x265 + x266 + x27*x86 + x34*x69 + x34*x71 + x34*x73 + x34*M[26] + x86*M[13] + M[71];
#pragma omp atomic
Ms[72] += x1*x177 + x149*z + x177*M[4] + x213*x273 + x267 + x268 + x269 + x270 + x271 + x272 + x29*x86 + x30*x86 + x31*x86 + x34*x75 + x34*x77 + x34*x79 + x34*M[27] + x41*x64 + x41*x65 + x41*x66 + x41*M[25] + x86*M[14] + M[72];
#pragma omp atomic
Ms[73] += x111*x280 + x155*z + x177*x3 + x177*M[5] + x18*x96 + x19*x96 + x20*x96 + x274 + x275 + x276 + x277 + x278 + x279 + x34*x81 + x34*x82 + x34*x83 + x34*M[28] + x41*x69 + x41*x70 + x41*x72 + x41*M[26] + x96*M[12] + M[73];
#pragma omp atomic
Ms[74] += x1*x184 + x162*z + x184*x2 + x184*M[4] + x23*x96 + x24*x96 + x26*x96 + x281 + x282 + x283 + x284 + x285 + x286 + x41*x75 + x41*x76 + x41*x78 + x41*M[27] + x96*M[13] + M[74];
#pragma omp atomic
Ms[75] += x184*x3 + x184*M[5] + x260*x290 + x287 + x288 + x289 + x29*x96 + x41*x81 + x41*M[28] + x96*M[14] + M[75];
#pragma omp atomic
Ms[76] += x170*M[6] + x291 + x292*M[1] + x34*M[29] + x86*M[15] + M[76];
#pragma omp atomic
Ms[77] += x170*x7 + x170*M[7] + x292*M[2] + x293 + x294 + x295 + x34*x88 + x34*M[30] + x36*x86 + x86*M[16] + M[77];
#pragma omp atomic
Ms[78] += x170*x8 + x170*M[8] + x177*M[6] + x296 + x297 + x298 + x299*M[1] + x34*x91 + x34*M[31] + x39*x86 + x41*x85 + x41*M[29] + x86*M[17] + M[78];
#pragma omp atomic
Ms[79] += x177*M[7] + x299*M[2] + x300 + x301 + x302 + x303*M[1] + x33*x96 + x34*x94 + x34*M[32] + x41*x87 + x41*M[30] + x42*x86 + x86*M[18] + x96*M[15] + M[79];
#pragma omp atomic
Ms[80] += x177*M[8] + x184*x5 + x184*M[6] + x303*M[2] + x304 + x305 + x306 + x34*x97 + x34*M[33] + x35*x96 + x41*x90 + x41*M[31] + x96*M[16] + M[80];
#pragma omp atomic
Ms[81] += x184*x6 + x184*M[7] + x307 + x308 + x309 + x310*M[1] + x38*x96 + x41*x93 + x41*M[32] + x96*M[17] + M[81];
#pragma omp atomic
Ms[82] += x184*M[8] + x310*M[2] + x311 + x41*M[33] + x96*M[18] + M[82];
#pragma omp atomic
Ms[83] += x100*M[9] + x11*M[34] + x188*M[3] + x312 + x314*M[0] + x45*M[19] + M[83];
#pragma omp atomic
Ms[84] += x100*x13 + x100*M[10] + x102*x11 + x11*M[35] + x188*x2 + x188*M[4] + x314*M[1] + x315 + x316 + x317 + x45*x47 + x45*M[20] + M[84];
#pragma omp atomic
Ms[85] += x100*x16 + x100*M[11] + x105*x11 + x11*M[36] + x188*x4 + x188*M[5] + x314*M[2] + x318 + x319 + x320 + x45*x50 + x45*M[21] + M[85];
#pragma omp atomic
Ms[86] += x100*x19 + x100*M[12] + x108*x11 + x11*M[37] + x111*x324 + x125*M[9] + x188*x5 + x188*M[6] + x214*M[3] + x321 + x322 + x323 + x34*x98 + x34*M[34] + x45*x53 + x45*M[22] + M[86];
#pragma omp atomic
Ms[87] += x100*x24 + x100*x25 + x100*x28 + x100*M[13] + x11*x113 + x11*x114 + x11*x117 + x11*M[38] + x188*x6 + x188*x7 + x188*M[7] + x191*z + x325 + x326 + x327 + x328 + x329 + x330 + x45*x56 + x45*x57 + x45*x60 + x45*M[23] + M[87];
#pragma omp atomic
Ms[88] += x100*x30 + x100*M[14] + x11*x119 + x11*M[39] + x121*M[9] + x188*x8 + x188*M[8] + x208*M[3] + x331 + x332 + x333 + x334*M[0] + x41*x98 + x41*M[34] + x45*x62 + x45*M[24] + M[88];
#pragma omp atomic
Ms[89] += x100*x33 + x100*M[15] + x101*x34 + x11*x123 + x11*M[40] + x125*M[10] + x213*x338 + x214*M[4] + x234*M[3] + x335 + x336 + x337 + x339*M[1] + x34*M[35] + x43*x86 + x45*x65 + x45*M[25] + x86*M[19] + M[89];
#pragma omp atomic
Ms[90] += x100*x35 + x100*x36 + x100*x37 + x100*M[16] + x104*x34 + x105*x34 + x106*x34 + x11*x127 + x11*x128 + x11*x131 + x11*M[41] + x125*x16 + x125*M[11] + x197*z + x214*x4 + x214*M[5] + x339*M[2] + x34*M[36] + x340 + x341 + x342 + x343 + x344 + x345 + x45*x70 + x45*x71 + x45*x74 + x45*M[26] + M[90];
#pragma omp atomic
Ms[91] += x100*x38 + x100*x39 + x100*x40 + x100*M[17] + x101*x41 + x102*x41 + x103*x41 + x11*x133 + x11*x134 + x11*x137 + x11*M[42] + x121*x13 + x121*M[10] + x2*x208 + x202*z + x208*M[4] + x334*M[1] + x346 + x347 + x348 + x349 + x350 + x351 + x41*M[35] + x45*x76 + x45*x77 + x45*x80 + x45*M[27] + M[91];
#pragma omp atomic
Ms[92] += x100*x42 + x100*M[18] + x104*x41 + x11*x139 + x11*M[43] + x121*M[11] + x208*M[5] + x230*M[3] + x334*M[2] + x352 + x353 + x354 + x355*M[0] + x41*M[36] + x43*x96 + x45*x82 + x45*M[28] + x96*M[19] + M[92];
#pragma omp atomic
Ms[93] += x107*x34 + x11*x142 + x11*M[44] + x125*M[12] + x170*x9 + x170*M[9] + x214*M[6] + x234*M[4] + x34*M[37] + x356 + x357 + x358 + x359*x360 + x361*M[1] + x45*x85 + x45*M[29] + x46*x86 + x86*M[20] + M[93];
#pragma omp atomic
Ms[94] += x11*x147 + x11*x148 + x11*x151 + x11*M[45] + x112*x34 + x114*x34 + x116*x34 + x125*x25 + x125*M[13] + x211*z + x214*x7 + x214*M[7] + x234*x4 + x234*M[5] + x34*M[38] + x361*M[2] + x362 + x363 + x364 + x365 + x366 + x367 + x45*x87 + x45*x88 + x45*x89 + x45*M[30] + x49*x86 + x50*x86 + x51*x86 + x86*M[21] + M[94];
#pragma omp atomic
Ms[95] += x107*x41 + x108*x41 + x109*x41 + x11*x153 + x11*x154 + x11*x157 + x11*M[46] + x111*x375 + x118*x34 + x119*x34 + x120*x34 + x121*x19 + x121*M[12] + x125*x30 + x125*M[14] + x177*x9 + x177*M[9] + x208*x5 + x208*M[6] + x214*x8 + x214*M[8] + x218*z + x34*M[39] + x368 + x369 + x370 + x371 + x372 + x373 + x41*M[37] + x45*x90 + x45*x91 + x45*x92 + x45*M[31] + M[95];
#pragma omp atomic
Ms[96] += x11*x160 + x11*x161 + x11*x164 + x11*M[47] + x112*x41 + x113*x41 + x115*x41 + x121*x24 + x121*M[13] + x2*x230 + x208*x6 + x208*M[7] + x224*z + x230*M[4] + x355*M[1] + x376 + x377 + x378 + x379 + x380 + x381 + x41*M[38] + x45*x93 + x45*x94 + x45*x95 + x45*M[32] + x46*x96 + x47*x96 + x48*x96 + x96*M[20] + M[96];
#pragma omp atomic
Ms[97] += x11*x166 + x11*M[48] + x118*x41 + x121*M[14] + x184*x9 + x184*M[9] + x208*M[8] + x230*M[5] + x355*M[2] + x382 + x383 + x384 + x385*M[0] + x41*M[39] + x45*x97 + x45*M[33] + x49*x96 + x96*M[21] + M[97];
#pragma omp atomic
Ms[98] += x0*x292 + x11*x169 + x11*M[49] + x12*x170 + x122*x34 + x125*M[15] + x170*M[10] + x234*M[6] + x292*M[3] + x34*M[40] + x386 + x387 + x388 + x389*M[1] + x52*x86 + x86*M[22] + M[98];
#pragma omp atomic
Ms[99] += x11*x171 + x11*x172 + x11*x173 + x11*M[50] + x125*x36 + x125*M[16] + x126*x34 + x128*x34 + x130*x34 + x15*x170 + x16*x170 + x17*x170 + x170*M[11] + x233*z + x234*x7 + x234*M[7] + x34*M[41] + x389*M[2] + x390 + x391 + x392 + x393 + x394 + x395 + x55*x86 + x57*x86 + x59*x86 + x86*M[23] + M[99];
#pragma omp atomic
Ms[100] += x0*x299 + x11*x174 + x11*x175 + x11*x176 + x11*M[51] + x12*x177 + x121*x33 + x121*M[15] + x122*x41 + x123*x41 + x124*x41 + x125*x39 + x125*M[17] + x132*x34 + x134*x34 + x136*x34 + x177*M[10] + x234*x8 + x234*M[8] + x238*z + x299*M[3] + x34*M[42] + x396 + x397 + x398 + x399 + x400 + x401 + x402*M[1] + x41*M[40] + x61*x86 + x62*x86 + x63*x86 + x86*M[24] + M[100];
#pragma omp atomic
Ms[101] += x0*x303 + x11*x178 + x11*x179 + x11*x180 + x11*M[52] + x121*x35 + x121*M[16] + x125*x42 + x125*M[18] + x126*x41 + x127*x41 + x129*x41 + x138*x34 + x139*x34 + x140*x34 + x15*x177 + x177*M[11] + x230*x5 + x230*M[6] + x244*z + x303*M[3] + x34*M[43] + x402*M[2] + x403 + x404 + x405 + x406 + x407 + x408 + x41*M[41] + x52*x96 + x53*x96 + x54*x96 + x96*M[22] + M[101];
#pragma omp atomic
Ms[102] += x11*x181 + x11*x182 + x11*x183 + x11*M[53] + x12*x184 + x121*x38 + x121*M[17] + x13*x184 + x132*x41 + x133*x41 + x135*x41 + x14*x184 + x184*M[10] + x230*x6 + x230*M[7] + x250*z + x385*M[1] + x409 + x41*M[42] + x410 + x411 + x412 + x413 + x414 + x55*x96 + x56*x96 + x58*x96 + x96*M[23] + M[102];
#pragma omp atomic
Ms[103] += x0*x310 + x11*x185 + x11*M[54] + x121*M[18] + x138*x41 + x15*x184 + x184*M[11] + x230*M[8] + x310*M[3] + x385*M[2] + x41*M[43] + x415 + x416 + x417 + x61*x96 + x96*M[24] + M[103];
#pragma omp atomic
Ms[104] += x1*x292 + x141*x34 + x170*x18 + x170*M[12] + x292*M[4] + x34*M[44] + x418 + x419 + x420 + x421*x422 + x64*x86 + x86*M[25] + M[104];
#pragma omp atomic
Ms[105] += x146*x34 + x148*x34 + x150*x34 + x170*x23 + x170*x25 + x170*x27 + x170*M[13] + x258*z + x292*x3 + x292*x4 + x292*M[5] + x34*M[45] + x423 + x424 + x425 + x426 + x427 + x428 + x69*x86 + x71*x86 + x73*x86 + x86*M[26] + M[105];
#pragma omp atomic
Ms[106] += x1*x299 + x141*x41 + x142*x41 + x143*x41 + x152*x34 + x154*x34 + x156*x34 + x170*x29 + x170*x30 + x170*x31 + x170*M[14] + x177*x18 + x177*M[12] + x264*z + x299*M[4] + x34*M[46] + x360*x435 + x41*M[44] + x429 + x430 + x431 + x432 + x433 + x434 + x75*x86 + x77*x86 + x79*x86 + x86*M[27] + M[106];
#pragma omp atomic
Ms[107] += x1*x303 + x146*x41 + x147*x41 + x149*x41 + x159*x34 + x161*x34 + x163*x34 + x177*x23 + x177*M[13] + x213*x442 + x270*z + x299*x3 + x299*M[5] + x303*M[4] + x34*M[47] + x41*M[45] + x436 + x437 + x438 + x439 + x440 + x441 + x64*x96 + x65*x96 + x66*x96 + x81*x86 + x82*x86 + x83*x86 + x86*M[28] + x96*M[25] + M[107];
#pragma omp atomic
Ms[108] += x111*x449 + x152*x41 + x153*x41 + x155*x41 + x165*x34 + x166*x34 + x167*x34 + x177*x29 + x177*M[14] + x18*x184 + x184*x19 + x184*x20 + x184*M[12] + x277*z + x3*x303 + x303*M[5] + x34*M[48] + x41*M[46] + x443 + x444 + x445 + x446 + x447 + x448 + x69*x96 + x70*x96 + x72*x96 + x96*M[26] + M[108];
#pragma omp atomic
Ms[109] += x1*x310 + x159*x41 + x160*x41 + x162*x41 + x184*x23 + x184*x24 + x184*x26 + x184*M[13] + x2*x310 + x284*z + x310*M[4] + x41*M[47] + x450 + x451 + x452 + x453 + x454 + x455 + x75*x96 + x76*x96 + x78*x96 + x96*M[27] + M[109];
#pragma omp atomic
Ms[110] += x165*x41 + x184*x29 + x184*M[14] + x3*x310 + x310*M[5] + x41*M[48] + x422*x459 + x456 + x457 + x458 + x81*x96 + x96*M[28] + M[110];
#pragma omp atomic
Ms[111] += x170*M[15] + x292*M[6] + x34*M[49] + x460 + x461*M[1] + x86*M[29] + M[111];
#pragma omp atomic
Ms[112] += x170*x36 + x170*M[16] + x172*x34 + x292*x7 + x292*M[7] + x34*M[50] + x461*M[2] + x462 + x463 + x464 + x86*x88 + x86*M[30] + M[112];
#pragma omp atomic
Ms[113] += x169*x41 + x170*x39 + x170*M[17] + x175*x34 + x177*M[15] + x292*x8 + x292*M[8] + x299*M[6] + x34*M[51] + x41*M[49] + x465 + x466 + x467 + x468*M[1] + x86*x91 + x86*M[31] + M[113];
#pragma omp atomic
Ms[114] += x170*x42 + x170*M[18] + x171*x41 + x177*M[16] + x179*x34 + x299*M[7] + x303*M[6] + x34*M[52] + x41*M[50] + x468*M[2] + x469 + x470 + x471 + x472*M[1] + x85*x96 + x86*x94 + x86*M[32] + x96*M[29] + M[114];
#pragma omp atomic
Ms[115] += x174*x41 + x177*M[17] + x182*x34 + x184*x33 + x184*M[15] + x299*M[8] + x303*M[7] + x34*M[53] + x41*M[51] + x472*M[2] + x473 + x474 + x475 + x476*M[1] + x86*x97 + x86*M[33] + x87*x96 + x96*M[30] + M[115];
#pragma omp atomic
Ms[116] += x177*M[18] + x178*x41 + x184*x35 + x184*M[16] + x185*x34 + x303*M[8] + x310*x5 + x310*M[6] + x34*M[54] + x41*M[52] + x476*M[2] + x477 + x478 + x479 + x90*x96 + x96*M[31] + M[116];
#pragma omp atomic
Ms[117] += x181*x41 + x184*x38 + x184*M[17] + x310*x6 + x310*M[7] + x41*M[53] + x480 + x481 + x482 + x483*M[1] + x93*x96 + x96*M[32] + M[117];
#pragma omp atomic
Ms[118] += x184*M[18] + x310*M[8] + x41*M[54] + x483*M[2] + x484 + x96*M[33] + M[118];
#pragma omp atomic
Ms[119] += x100*M[19] + x11*M[55] + x188*M[9] + x314*M[3] + x45*M[34] + x485 + x487*M[0] + M[119];
#pragma omp atomic
Ms[120] += x100*x47 + x100*M[20] + x102*x45 + x11*x190 + x11*M[56] + x13*x188 + x188*M[10] + x2*x314 + x314*M[4] + x45*M[35] + x487*M[1] + x488 + x489 + x490 + M[120];
#pragma omp atomic
Ms[121] += x100*x50 + x100*M[21] + x105*x45 + x11*x193 + x11*M[57] + x16*x188 + x188*M[11] + x314*x4 + x314*M[5] + x45*M[36] + x487*M[2] + x491 + x492 + x493 + M[121];
#pragma omp atomic
Ms[122] += x100*x53 + x100*M[22] + x108*x45 + x11*x196 + x11*M[58] + x111*x497 + x125*M[19] + x186*x34 + x188*x19 + x188*M[12] + x214*M[9] + x314*x5 + x314*M[6] + x339*M[3] + x34*M[55] + x45*M[37] + x494 + x495 + x496 + M[122];
#pragma omp atomic
Ms[123] += x100*x56 + x100*x57 + x100*x60 + x100*M[23] + x11*x200 + x11*x201 + x11*x204 + x11*M[59] + x113*x45 + x114*x45 + x117*x45 + x188*x24 + x188*x25 + x188*x28 + x188*M[13] + x314*x6 + x314*x7 + x314*M[7] + x317*z + x45*M[38] + x498 + x499 + x500 + x501 + x502 + x503 + M[123];
#pragma omp atomic
Ms[124] += x100*x62 + x100*M[24] + x11*x206 + x11*M[60] + x119*x45 + x121*M[19] + x186*x41 + x188*x30 + x188*M[14] + x208*M[9] + x314*x8 + x314*M[8] + x334*M[3] + x41*M[55] + x45*M[39] + x504 + x505 + x506 + x507*M[0] + M[124];
#pragma omp atomic
Ms[125] += x100*x65 + x100*M[25] + x11*x210 + x11*M[61] + x123*x45 + x125*M[20] + x188*x33 + x188*M[15] + x189*x34 + x213*x511 + x214*M[10] + x234*M[9] + x339*M[4] + x34*M[56] + x361*M[3] + x45*M[40] + x508 + x509 + x510 + x512*M[1] + x86*x98 + x86*M[34] + M[125];
#pragma omp atomic
Ms[126] += x100*x70 + x100*x71 + x100*x74 + x100*M[26] + x11*x216 + x11*x217 + x11*x220 + x11*M[62] + x125*x50 + x125*M[21] + x127*x45 + x128*x45 + x131*x45 + x16*x214 + x188*x35 + x188*x36 + x188*x37 + x188*M[16] + x192*x34 + x193*x34 + x194*x34 + x214*M[11] + x323*z + x339*x4 + x339*M[5] + x34*M[57] + x45*M[41] + x512*M[2] + x513 + x514 + x515 + x516 + x517 + x518 + M[126];
#pragma omp atomic
Ms[127] += x100*x76 + x100*x77 + x100*x80 + x100*M[27] + x11*x222 + x11*x223 + x11*x226 + x11*M[63] + x121*x47 + x121*M[20] + x13*x208 + x133*x45 + x134*x45 + x137*x45 + x188*x38 + x188*x39 + x188*x40 + x188*M[17] + x189*x41 + x190*x41 + x191*x41 + x2*x334 + x208*M[10] + x328*z + x334*M[4] + x41*M[56] + x45*M[42] + x507*M[1] + x519 + x520 + x521 + x522 + x523 + x524 + M[127];
#pragma omp atomic
Ms[128] += x100*x82 + x100*M[28] + x11*x228 + x11*M[64] + x121*M[21] + x139*x45 + x188*x42 + x188*M[18] + x192*x41 + x208*M[11] + x230*M[9] + x334*M[5] + x355*M[3] + x41*M[57] + x45*M[43] + x507*M[2] + x525 + x526 + x527 + x528*M[0] + x96*x98 + x96*M[34] + M[128];
#pragma omp atomic
Ms[129] += x100*x85 + x100*M[29] + x101*x86 + x11*x232 + x11*M[65] + x125*M[22] + x142*x45 + x170*x43 + x170*M[19] + x195*x34 + x214*M[12] + x234*M[10] + x339*M[6] + x34*M[58] + x360*x532 + x361*M[4] + x389*M[3] + x45*M[44] + x529 + x530 + x531 + x533*M[1] + x86*M[35] + M[129];
#pragma omp atomic
Ms[130] += x100*x87 + x100*x88 + x100*x89 + x100*M[30] + x104*x86 + x105*x86 + x106*x86 + x11*x236 + x11*x237 + x11*x240 + x11*M[66] + x125*x57 + x125*M[23] + x147*x45 + x148*x45 + x151*x45 + x16*x234 + x199*x34 + x201*x34 + x203*x34 + x214*x25 + x214*M[13] + x234*M[11] + x337*z + x339*x7 + x339*M[7] + x34*M[59] + x361*x4 + x361*M[5] + x45*M[45] + x533*M[2] + x534 + x535 + x536 + x537 + x538 + x539 + x86*M[36] + M[130];
#pragma omp atomic
Ms[131] += x100*x90 + x100*x91 + x100*x92 + x100*M[31] + x11*x242 + x11*x243 + x11*x246 + x11*M[67] + x121*x53 + x121*M[22] + x125*x62 + x125*M[24] + x153*x45 + x154*x45 + x157*x45 + x177*x43 + x177*M[19] + x19*x208 + x195*x41 + x196*x41 + x197*x41 + x205*x34 + x206*x34 + x207*x34 + x208*M[12] + x214*x30 + x214*M[14] + x334*x5 + x334*M[6] + x339*x8 + x339*M[8] + x34*M[60] + x343*z + x402*M[3] + x41*M[58] + x45*M[46] + x540 + x541 + x542 + x543 + x544 + x545 + x546*x547 + M[131];
#pragma omp atomic
Ms[132] += x100*x93 + x100*x94 + x100*x95 + x100*M[32] + x101*x96 + x102*x96 + x103*x96 + x11*x248 + x11*x249 + x11*x252 + x11*M[68] + x121*x56 + x121*M[23] + x13*x230 + x160*x45 + x161*x45 + x164*x45 + x199*x41 + x2*x355 + x200*x41 + x202*x41 + x208*x24 + x208*M[13] + x230*M[10] + x334*x6 + x334*M[7] + x349*z + x355*M[4] + x41*M[59] + x45*M[47] + x528*M[1] + x548 + x549 + x550 + x551 + x552 + x553 + x96*M[35] + M[132];
#pragma omp atomic
Ms[133] += x100*x97 + x100*M[33] + x104*x96 + x11*x254 + x11*M[69] + x121*M[24] + x166*x45 + x184*x43 + x184*M[19] + x205*x41 + x208*M[14] + x230*M[11] + x334*M[8] + x355*M[5] + x385*M[3] + x41*M[60] + x45*M[48] + x528*M[2] + x554 + x555 + x556 + x557*M[0] + x96*M[36] + M[133];
#pragma omp atomic
Ms[134] += x107*x86 + x11*x257 + x11*M[70] + x125*M[25] + x169*x45 + x170*x46 + x170*M[20] + x209*x34 + x214*M[15] + x234*M[12] + x292*x9 + x292*M[9] + x34*M[61] + x361*M[6] + x389*M[4] + x45*M[49] + x558 + x559 + x560 + x561*x562 + x563*M[1] + x86*M[37] + M[134];
#pragma omp atomic
Ms[135] += x11*x262 + x11*x263 + x11*x266 + x11*M[71] + x112*x86 + x114*x86 + x116*x86 + x125*x71 + x125*M[26] + x170*x49 + x170*x50 + x170*x51 + x170*M[21] + x171*x45 + x172*x45 + x173*x45 + x214*x36 + x214*M[16] + x215*x34 + x217*x34 + x219*x34 + x234*x25 + x234*M[13] + x34*M[62] + x358*z + x361*x7 + x361*M[7] + x389*x4 + x389*M[5] + x45*M[50] + x563*M[2] + x564 + x565 + x566 + x567 + x568 + x569 + x86*M[38] + M[135];
#pragma omp atomic
Ms[136] += x11*x268 + x11*x269 + x11*x272 + x11*M[72] + x118*x86 + x119*x86 + x120*x86 + x121*x65 + x121*M[25] + x125*x77 + x125*M[27] + x145*x576 + x174*x45 + x175*x45 + x176*x45 + x177*x46 + x177*M[20] + x208*x33 + x208*M[15] + x209*x41 + x210*x41 + x211*x41 + x214*x39 + x214*M[17] + x221*x34 + x223*x34 + x225*x34 + x234*x30 + x234*M[14] + x299*x9 + x299*M[9] + x34*M[63] + x361*x8 + x361*M[8] + x365*z + x402*M[4] + x41*M[61] + x45*M[51] + x570 + x571 + x572 + x573 + x574 + x575 + x577*x578 + x86*M[39] + M[136];
#pragma omp atomic
Ms[137] += x107*x96 + x108*x96 + x109*x96 + x11*x275 + x11*x276 + x11*x279 + x11*M[73] + x121*x70 + x121*M[26] + x125*x82 + x125*M[28] + x177*x49 + x177*M[21] + x178*x45 + x179*x45 + x180*x45 + x19*x230 + x208*x35 + x208*M[16] + x214*x42 + x214*M[18] + x215*x41 + x216*x41 + x218*x41 + x227*x34 + x228*x34 + x229*x34 + x230*M[12] + x303*x9 + x303*M[9] + x34*M[64] + x355*x5 + x355*M[6] + x371*z + x402*M[5] + x41*M[62] + x45*M[52] + x546*x585 + x547*x586 + x579 + x580 + x581 + x582 + x583 + x584 + x96*M[37] + M[137];
#pragma omp atomic
Ms[138] += x11*x282 + x11*x283 + x11*x286 + x11*M[74] + x112*x96 + x113*x96 + x115*x96 + x121*x76 + x121*M[27] + x181*x45 + x182*x45 + x183*x45 + x184*x46 + x184*x47 + x184*x48 + x184*M[20] + x2*x385 + x208*x38 + x208*M[17] + x221*x41 + x222*x41 + x224*x41 + x230*x24 + x230*M[13] + x355*x6 + x355*M[7] + x379*z + x385*M[4] + x41*M[63] + x45*M[53] + x557*M[1] + x587 + x588 + x589 + x590 + x591 + x592 + x96*M[38] + M[138];
#pragma omp atomic
Ms[139] += x11*x288 + x11*M[75] + x118*x96 + x121*M[28] + x184*x49 + x184*M[21] + x185*x45 + x208*M[18] + x227*x41 + x230*M[14] + x310*x9 + x310*M[9] + x355*M[8] + x385*M[5] + x41*M[64] + x45*M[54] + x557*M[2] + x593 + x594 + x595 + x596*M[0] + x96*M[39] + M[139];
#pragma omp atomic
Ms[140] += x0*x461 + x11*x291 + x11*M[76] + x12*x292 + x122*x86 + x125*M[29] + x170*x52 + x170*M[22] + x231*x34 + x234*M[15] + x292*M[10] + x34*M[65] + x389*M[6] + x461*M[3] + x597 + x598 + x599 + x600*M[1] + x86*M[40] + M[140];
#pragma omp atomic
Ms[141] += x11*x293 + x11*x294 + x11*x295 + x11*M[77] + x125*x88 + x125*M[30] + x126*x86 + x128*x86 + x130*x86 + x15*x292 + x16*x292 + x17*x292 + x170*x55 + x170*x57 + x170*x59 + x170*M[23] + x234*x36 + x234*M[16] + x235*x34 + x237*x34 + x239*x34 + x292*M[11] + x34*M[66] + x388*z + x389*x7 + x389*M[7] + x600*M[2] + x601 + x602 + x603 + x604 + x605 + x606 + x86*M[41] + M[141];
#pragma omp atomic
Ms[142] += x0*x468 + x11*x296 + x11*x297 + x11*x298 + x11*M[78] + x12*x299 + x121*x85 + x121*M[29] + x125*x91 + x125*M[31] + x132*x86 + x134*x86 + x136*x86 + x170*x61 + x170*x62 + x170*x63 + x170*M[24] + x177*x52 + x177*M[22] + x231*x41 + x232*x41 + x233*x41 + x234*x39 + x234*M[17] + x241*x34 + x243*x34 + x245*x34 + x299*M[10] + x34*M[67] + x389*x8 + x389*M[8] + x393*z + x402*M[6] + x41*M[65] + x468*M[3] + x607 + x608 + x609 + x610 + x611 + x612 + x613*x614 + x86*M[42] + M[142];
#pragma omp atomic
Ms[143] += x0*x472 + x11*x300 + x11*x301 + x11*x302 + x11*M[79] + x12*x303 + x121*x87 + x121*M[30] + x122*x96 + x123*x96 + x124*x96 + x125*x94 + x125*M[32] + x138*x86 + x139*x86 + x140*x86 + x15*x299 + x177*x55 + x177*M[23] + x230*x33 + x230*M[15] + x234*x42 + x234*M[18] + x235*x41 + x236*x41 + x238*x41 + x247*x34 + x249*x34 + x251*x34 + x299*M[11] + x303*M[10] + x34*M[68] + x399*z + x402*M[7] + x41*M[66] + x472*M[3] + x577*x622 + x614*x623 + x615 + x616 + x617 + x618 + x619 + x620 + x86*M[43] + x96*M[40] + M[143];
#pragma omp atomic
Ms[144] += x0*x476 + x11*x304 + x11*x305 + x11*x306 + x11*M[80] + x121*x90 + x121*M[31] + x125*x97 + x125*M[33] + x126*x96 + x127*x96 + x129*x96 + x15*x303 + x177*x61 + x177*M[24] + x184*x52 + x184*x53 + x184*x54 + x184*M[22] + x230*x35 + x230*M[16] + x241*x41 + x242*x41 + x244*x41 + x253*x34 + x254*x34 + x255*x34 + x303*M[11] + x34*M[69] + x385*x5 + x385*M[6] + x402*M[8] + x406*z + x41*M[67] + x476*M[3] + x585*x586 + x624 + x625 + x626 + x627 + x628 + x629 + x96*M[41] + M[144];
#pragma omp atomic
Ms[145] += x11*x307 + x11*x308 + x11*x309 + x11*M[81] + x12*x310 + x121*x93 + x121*M[32] + x13*x310 + x132*x96 + x133*x96 + x135*x96 + x14*x310 + x184*x55 + x184*x56 + x184*x58 + x184*M[23] + x230*x38 + x230*M[17] + x247*x41 + x248*x41 + x250*x41 + x310*M[10] + x385*x6 + x385*M[7] + x41*M[68] + x412*z + x596*M[1] + x630 + x631 + x632 + x633 + x634 + x635 + x96*M[42] + M[145];
#pragma omp atomic
Ms[146] += x0*x483 + x11*x311 + x11*M[82] + x121*M[33] + x138*x96 + x15*x310 + x184*x61 + x184*M[24] + x230*M[18] + x253*x41 + x310*M[11] + x385*M[8] + x41*M[69] + x483*M[3] + x596*M[2] + x636 + x637 + x638 + x96*M[43] + M[146];
#pragma omp atomic
Ms[147] += x1*x461 + x141*x86 + x170*x64 + x170*M[25] + x18*x292 + x256*x34 + x292*M[12] + x34*M[70] + x461*M[4] + x639 + x640 + x641 + x642*x643 + x86*M[44] + M[147];
#pragma omp atomic
Ms[148] += x146*x86 + x148*x86 + x150*x86 + x170*x69 + x170*x71 + x170*x73 + x170*M[26] + x23*x292 + x25*x292 + x261*x34 + x263*x34 + x265*x34 + x27*x292 + x292*M[13] + x3*x461 + x34*M[71] + x4*x461 + x420*z + x461*M[5] + x644 + x645 + x646 + x647 + x648 + x649 + x86*M[45] + M[148];
#pragma omp atomic
Ms[149] += x1*x468 + x152*x86 + x154*x86 + x156*x86 + x170*x75 + x170*x77 + x170*x79 + x170*M[27] + x177*x64 + x177*M[25] + x18*x299 + x256*x41 + x257*x41 + x258*x41 + x267*x34 + x269*x34 + x271*x34 + x29*x292 + x292*x30 + x292*x31 + x292*M[14] + x299*M[12] + x34*M[72] + x41*M[70] + x426*z + x468*M[4] + x562*x656 + x650 + x651 + x652 + x653 + x654 + x655 + x86*M[46] + M[149];
#pragma omp atomic
Ms[150] += x1*x472 + x141*x96 + x142*x96 + x143*x96 + x159*x86 + x161*x86 + x163*x86 + x170*x81 + x170*x82 + x170*x83 + x170*M[28] + x177*x69 + x177*M[26] + x18*x303 + x23*x299 + x261*x41 + x262*x41 + x264*x41 + x274*x34 + x276*x34 + x278*x34 + x299*M[13] + x3*x468 + x303*M[12] + x34*M[73] + x360*x663 + x41*M[71] + x432*z + x468*M[5] + x472*M[4] + x657 + x658 + x659 + x660 + x661 + x662 + x86*M[47] + x96*M[44] + M[150];
#pragma omp atomic
Ms[151] += x1*x476 + x146*x96 + x147*x96 + x149*x96 + x165*x86 + x166*x86 + x167*x86 + x177*x75 + x177*M[27] + x184*x64 + x184*x65 + x184*x66 + x184*M[25] + x213*x670 + x23*x303 + x267*x41 + x268*x41 + x270*x41 + x281*x34 + x283*x34 + x285*x34 + x29*x299 + x299*M[14] + x3*x472 + x303*M[13] + x34*M[74] + x41*M[72] + x439*z + x472*M[5] + x476*M[4] + x664 + x665 + x666 + x667 + x668 + x669 + x86*M[48] + x96*M[45] + M[151];
#pragma omp atomic
Ms[152] += x111*x677 + x152*x96 + x153*x96 + x155*x96 + x177*x81 + x177*M[28] + x18*x310 + x184*x69 + x184*x70 + x184*x72 + x184*M[26] + x19*x310 + x20*x310 + x274*x41 + x275*x41 + x277*x41 + x287*x34 + x288*x34 + x289*x34 + x29*x303 + x3*x476 + x303*M[14] + x310*M[12] + x34*M[75] + x41*M[73] + x446*z + x476*M[5] + x671 + x672 + x673 + x674 + x675 + x676 + x96*M[46] + M[152];
#pragma omp atomic
Ms[153] += x1*x483 + x159*x96 + x160*x96 + x162*x96 + x184*x75 + x184*x76 + x184*x78 + x184*M[27] + x2*x483 + x23*x310 + x24*x310 + x26*x310 + x281*x41 + x282*x41 + x284*x41 + x310*M[13] + x41*M[74] + x453*z + x483*M[4] + x678 + x679 + x680 + x681 + x682 + x683 + x96*M[47] + M[153];
#pragma omp atomic
Ms[154] += x165*x96 + x184*x81 + x184*M[28] + x287*x41 + x29*x310 + x3*x483 + x310*M[14] + x41*M[75] + x483*M[5] + x643*x687 + x684 + x685 + x686 + x96*M[48] + M[154];
#pragma omp atomic
Ms[155] += x170*M[29] + x292*M[15] + x34*M[76] + x461*M[6] + x688 + x689*M[1] + x86*M[49] + M[155];
#pragma omp atomic
Ms[156] += x170*x88 + x170*M[30] + x172*x86 + x292*x36 + x292*M[16] + x294*x34 + x34*M[77] + x461*x7 + x461*M[7] + x689*M[2] + x690 + x691 + x692 + x86*M[50] + M[156];
#pragma omp atomic
Ms[157] += x170*x91 + x170*M[31] + x175*x86 + x177*M[29] + x291*x41 + x292*x39 + x292*M[17] + x297*x34 + x299*M[15] + x34*M[78] + x41*M[76] + x461*x8 + x461*M[8] + x468*M[6] + x693 + x694 + x695 + x696*M[1] + x86*M[51] + M[157];
#pragma omp atomic
Ms[158] += x169*x96 + x170*x94 + x170*M[32] + x177*M[30] + x179*x86 + x292*x42 + x292*M[18] + x293*x41 + x299*M[16] + x301*x34 + x303*M[15] + x34*M[79] + x41*M[77] + x468*M[7] + x472*M[6] + x696*M[2] + x697 + x698 + x699 + x700*M[1] + x86*M[52] + x96*M[49] + M[158];
#pragma omp atomic
Ms[159] += x170*x97 + x170*M[33] + x171*x96 + x177*M[31] + x182*x86 + x184*x85 + x184*M[29] + x296*x41 + x299*M[17] + x303*M[16] + x305*x34 + x34*M[80] + x41*M[78] + x468*M[8] + x472*M[7] + x476*M[6] + x700*M[2] + x701 + x702 + x703 + x704*M[1] + x86*M[53] + x96*M[50] + M[159];
#pragma omp atomic
Ms[160] += x174*x96 + x177*M[32] + x184*x87 + x184*M[30] + x185*x86 + x299*M[18] + x300*x41 + x303*M[17] + x308*x34 + x310*x33 + x310*M[15] + x34*M[81] + x41*M[79] + x472*M[8] + x476*M[7] + x704*M[2] + x705 + x706 + x707 + x708*M[1] + x86*M[54] + x96*M[51] + M[160];
#pragma omp atomic
Ms[161] += x177*M[33] + x178*x96 + x184*x90 + x184*M[31] + x303*M[18] + x304*x41 + x310*x35 + x310*M[16] + x311*x34 + x34*M[82] + x41*M[80] + x476*M[8] + x483*x5 + x483*M[6] + x708*M[2] + x709 + x710 + x711 + x96*M[52] + M[161];
#pragma omp atomic
Ms[162] += x181*x96 + x184*x93 + x184*M[32] + x307*x41 + x310*x38 + x310*M[17] + x41*M[81] + x483*x6 + x483*M[7] + x712 + x713 + x714 + x715*M[1] + x96*M[53] + M[162];
#pragma omp atomic
Ms[163] += x184*M[33] + x310*M[18] + x41*M[82] + x483*M[8] + x715*M[2] + x716 + x96*M[54] + M[163];
#pragma omp atomic
Ms[164] += x100*M[34] + x11*M[83] + x188*M[19] + x314*M[9] + x45*M[55] + x487*M[3] + x717 + x718*M[0] + M[164];
#pragma omp atomic
Ms[165] += x100*x102 + x100*M[35] + x11*x316 + x11*M[84] + x13*x314 + x188*x47 + x188*M[20] + x190*x45 + x2*x487 + x314*M[10] + x45*M[56] + x487*M[4] + x718*M[1] + x719 + x720 + x721 + M[165];
#pragma omp atomic
Ms[166] += x100*x105 + x100*M[36] + x11*x319 + x11*M[85] + x16*x314 + x188*x50 + x188*M[21] + x193*x45 + x314*M[11] + x4*x487 + x45*M[57] + x485*z + x487*M[5] + x718*M[2] + x722 + z*M[119] + M[166];
#pragma omp atomic
Ms[167] += x100*x108 + x100*M[37] + x11*x322 + x11*M[86] + x111*x726 + x125*M[34] + x188*x53 + x188*M[22] + x19*x314 + x196*x45 + x214*M[19] + x312*x34 + x314*M[12] + x339*M[9] + x34*M[83] + x45*M[58] + x487*x5 + x487*M[6] + x512*M[3] + x723 + x724 + x725 + M[167];
#pragma omp atomic
Ms[168] += x100*x113 + x100*x114 + x100*x117 + x100*M[38] + x11*x326 + x11*x327 + x11*x330 + x11*M[87] + x188*x56 + x188*x57 + x188*x60 + x188*M[23] + x200*x45 + x201*x45 + x204*x45 + x24*x314 + x25*x314 + x28*x314 + x314*M[13] + x45*M[59] + x487*x6 + x487*x7 + x487*M[7] + x488*z + x489*z + x490*z + x727 + x728 + x729 + z*M[120] + M[168];
#pragma omp atomic
Ms[169] += x100*x119 + x100*M[39] + x11*x332 + x11*M[88] + x121*M[34] + x188*x62 + x188*M[24] + x206*x45 + x208*M[19] + x30*x314 + x312*x41 + x314*M[14] + x334*M[9] + x41*M[83] + x45*M[60] + x487*x8 + x487*M[8] + x491*z + x507*M[3] + x730 + x731*M[0] + z*M[121] + M[169];
#pragma omp atomic
Ms[170] += x100*x123 + x100*M[40] + x11*x336 + x11*M[89] + x125*M[35] + x186*x86 + x188*x65 + x188*M[25] + x210*x45 + x214*M[20] + x234*M[19] + x314*x33 + x314*M[15] + x315*x34 + x339*M[10] + x34*M[84] + x361*M[9] + x45*M[61] + x512*M[4] + x533*M[3] + x67*x735 + x732 + x733 + x734 + x736*M[1] + x86*M[55] + M[170];
#pragma omp atomic
Ms[171] += x100*x127 + x100*x128 + x100*x131 + x100*M[41] + x105*x125 + x11*x341 + x11*x342 + x11*x345 + x11*M[90] + x125*M[36] + x16*x339 + x188*x70 + x188*x71 + x188*x74 + x188*M[26] + x214*x50 + x214*M[21] + x216*x45 + x217*x45 + x220*x45 + x314*x35 + x314*x36 + x314*x37 + x314*M[16] + x318*x34 + x319*x34 + x320*x34 + x339*M[11] + x34*M[85] + x4*x512 + x45*M[62] + x494*z + x495*z + x496*z + x512*M[5] + x736*M[2] + x737 + x738 + x739 + z*M[122] + M[171];
#pragma omp atomic
Ms[172] += x100*x133 + x100*x134 + x100*x137 + x100*M[42] + x102*x121 + x11*x347 + x11*x348 + x11*x351 + x11*M[91] + x121*M[35] + x13*x334 + x188*x76 + x188*x77 + x188*x80 + x188*M[27] + x2*x507 + x208*x47 + x208*M[20] + x222*x45 + x223*x45 + x226*x45 + x314*x38 + x314*x39 + x314*x40 + x314*M[17] + x315*x41 + x316*x41 + x317*x41 + x334*M[10] + x41*M[84] + x45*M[63] + x498*z + x499*z + x501*z + x507*M[4] + x731*M[1] + x740 + x741 + x742 + z*M[123] + M[172];
#pragma omp atomic
Ms[173] += x100*x139 + x100*M[43] + x11*x353 + x11*M[92] + x121*M[36] + x186*x96 + x188*x82 + x188*M[28] + x208*M[21] + x228*x45 + x230*M[19] + x314*x42 + x314*M[18] + x318*x41 + x334*M[11] + x355*M[9] + x41*M[85] + x45*M[64] + x504*z + x507*M[5] + x528*M[3] + x731*M[2] + x735*x84 + x743 + x96*M[55] + z*M[124] + M[173];
#pragma omp atomic
Ms[174] += x100*x142 + x100*M[44] + x11*x357 + x11*M[93] + x125*M[37] + x170*x98 + x170*M[34] + x188*x85 + x188*M[29] + x189*x86 + x214*M[22] + x232*x45 + x234*M[20] + x321*x34 + x339*M[12] + x34*M[86] + x360*x747 + x361*M[10] + x389*M[9] + x45*M[65] + x512*M[6] + x533*M[4] + x563*M[3] + x613*x748 + x744 + x745 + x746 + x86*M[56] + M[174];
#pragma omp atomic
Ms[175] += x100*x147 + x100*x148 + x100*x151 + x100*M[45] + x11*x363 + x11*x364 + x11*x367 + x11*M[94] + x114*x125 + x125*M[38] + x16*x361 + x188*x87 + x188*x88 + x188*x89 + x188*M[30] + x192*x86 + x193*x86 + x194*x86 + x214*x57 + x214*M[23] + x234*x50 + x234*M[21] + x236*x45 + x237*x45 + x240*x45 + x25*x339 + x325*x34 + x327*x34 + x329*x34 + x339*M[13] + x34*M[87] + x361*M[11] + x4*x533 + x45*M[66] + x508*z + x509*z + x510*z + x512*x7 + x512*M[7] + x533*M[5] + x67*x752 + x749 + x750 + x751 + x86*M[57] + z*M[125] + M[175];
#pragma omp atomic
Ms[176] += x100*x153 + x100*x154 + x100*x157 + x100*M[46] + x108*x121 + x11*x369 + x11*x370 + x11*x373 + x11*M[95] + x119*x125 + x121*M[37] + x125*M[39] + x177*x98 + x177*M[34] + x188*x90 + x188*x91 + x188*x92 + x188*M[31] + x19*x334 + x208*x53 + x208*M[22] + x214*x62 + x214*M[24] + x242*x45 + x243*x45 + x246*x45 + x30*x339 + x321*x41 + x322*x41 + x323*x41 + x331*x34 + x332*x34 + x333*x34 + x334*M[12] + x339*M[14] + x34*M[88] + x402*M[9] + x41*M[86] + x45*M[67] + x5*x507 + x507*M[6] + x512*x8 + x512*M[8] + x513*z + x514*z + x516*z + x753 + x754 + x755 + x756*M[3] + x757*x758 + z*M[126] + M[176];
#pragma omp atomic
Ms[177] += x100*x160 + x100*x161 + x100*x164 + x100*M[47] + x11*x377 + x11*x378 + x11*x381 + x11*M[96] + x113*x121 + x121*M[38] + x13*x355 + x188*x93 + x188*x94 + x188*x95 + x188*M[32] + x189*x96 + x190*x96 + x191*x96 + x2*x528 + x208*x56 + x208*M[23] + x230*x47 + x230*M[20] + x24*x334 + x248*x45 + x249*x45 + x252*x45 + x325*x41 + x326*x41 + x328*x41 + x334*M[13] + x355*M[10] + x41*M[87] + x45*M[68] + x507*x6 + x507*M[7] + x519*z + x520*z + x522*z + x528*M[4] + x621*x748 + x759 + x760 + x761 + x96*M[56] + z*M[127] + M[177];
#pragma omp atomic
Ms[178] += x100*x166 + x100*M[48] + x11*x383 + x11*M[97] + x121*M[39] + x184*x98 + x184*M[34] + x188*x97 + x188*M[33] + x192*x96 + x208*M[24] + x230*M[21] + x254*x45 + x331*x41 + x334*M[14] + x355*M[11] + x385*M[9] + x41*M[88] + x45*M[69] + x507*M[8] + x525*z + x528*M[5] + x557*M[3] + x752*x84 + x762 + x763*M[0] + x96*M[57] + z*M[128] + M[178];
#pragma omp atomic
Ms[179] += x100*x169 + x100*M[49] + x101*x170 + x11*x387 + x11*M[98] + x125*M[40] + x170*M[35] + x195*x86 + x214*M[25] + x234*M[22] + x257*x45 + x259*x767 + x292*x43 + x292*M[19] + x335*x34 + x339*M[15] + x34*M[89] + x361*M[12] + x389*M[10] + x45*M[70] + x533*M[6] + x563*M[4] + x600*M[3] + x764 + x765 + x766 + x768*M[1] + x86*M[58] + M[179];
#pragma omp atomic
Ms[180] += x100*x171 + x100*x172 + x100*x173 + x100*M[50] + x104*x170 + x105*x170 + x106*x170 + x11*x391 + x11*x392 + x11*x395 + x11*M[99] + x125*x128 + x125*M[41] + x16*x389 + x170*M[36] + x199*x86 + x201*x86 + x203*x86 + x214*x71 + x214*M[26] + x234*x57 + x234*M[23] + x25*x361 + x262*x45 + x263*x45 + x266*x45 + x339*x36 + x339*M[16] + x34*x340 + x34*x342 + x34*x344 + x34*M[90] + x361*M[13] + x389*M[11] + x4*x563 + x45*M[71] + x529*z + x530*z + x531*z + x533*x7 + x533*M[7] + x563*M[5] + x768*M[2] + x769 + x770 + x771 + x86*M[59] + z*M[129] + M[180];
#pragma omp atomic
Ms[181] += x100*x174 + x100*x175 + x100*x176 + x100*M[51] + x101*x177 + x11*x397 + x11*x398 + x11*x401 + x11*M[100] + x121*x123 + x121*M[40] + x125*x134 + x125*M[42] + x177*M[35] + x205*x86 + x206*x86 + x207*x86 + x208*x65 + x208*M[25] + x213*x776 + x214*x77 + x214*M[27] + x234*x62 + x234*M[24] + x268*x45 + x269*x45 + x272*x45 + x299*x43 + x299*M[19] + x30*x361 + x33*x334 + x334*M[15] + x335*x41 + x336*x41 + x337*x41 + x339*x39 + x339*M[17] + x34*x346 + x34*x348 + x34*x350 + x34*M[91] + x361*M[14] + x402*M[10] + x41*M[89] + x45*M[72] + x533*x8 + x533*M[8] + x534*z + x535*z + x537*z + x756*M[4] + x772 + x773 + x774 + x775*M[3] + x777*x778 + x86*M[60] + z*M[130] + M[181];
#pragma omp atomic
Ms[182] += x100*x178 + x100*x179 + x100*x180 + x100*M[52] + x104*x177 + x11*x404 + x11*x405 + x11*x408 + x11*M[101] + x121*x127 + x121*M[41] + x125*x139 + x125*M[43] + x177*M[36] + x19*x355 + x195*x96 + x196*x96 + x197*x96 + x208*x70 + x208*M[26] + x214*x82 + x214*M[28] + x230*x53 + x230*M[22] + x275*x45 + x276*x45 + x279*x45 + x303*x43 + x303*M[19] + x334*x35 + x334*M[16] + x339*x42 + x339*M[18] + x34*x352 + x34*x353 + x34*x354 + x34*M[92] + x340*x41 + x341*x41 + x343*x41 + x355*M[12] + x402*M[11] + x41*M[90] + x45*M[73] + x5*x528 + x528*M[6] + x540*z + x541*z + x543*z + x756*M[5] + x778*x785 + x779 + x780 + x781 + x782*M[3] + x783*x784 + x96*M[58] + z*M[131] + M[182];
#pragma omp atomic
Ms[183] += x100*x181 + x100*x182 + x100*x183 + x100*M[53] + x101*x184 + x102*x184 + x103*x184 + x11*x410 + x11*x411 + x11*x414 + x11*M[102] + x121*x133 + x121*M[42] + x13*x385 + x184*M[35] + x199*x96 + x2*x557 + x200*x96 + x202*x96 + x208*x76 + x208*M[27] + x230*x56 + x230*M[23] + x24*x355 + x282*x45 + x283*x45 + x286*x45 + x334*x38 + x334*M[17] + x346*x41 + x347*x41 + x349*x41 + x355*M[13] + x385*M[10] + x41*M[91] + x45*M[74] + x528*x6 + x528*M[7] + x548*z + x549*z + x551*z + x557*M[4] + x763*M[1] + x786 + x787 + x788 + x96*M[59] + z*M[132] + M[183];
#pragma omp atomic
Ms[184] += x100*x185 + x100*M[54] + x104*x184 + x11*x416 + x11*M[103] + x121*M[43] + x184*M[36] + x205*x96 + x208*M[28] + x230*M[24] + x288*x45 + x290*x767 + x310*x43 + x310*M[19] + x334*M[18] + x352*x41 + x355*M[14] + x385*M[11] + x41*M[92] + x45*M[75] + x528*M[8] + x554*z + x557*M[5] + x596*M[3] + x763*M[2] + x789 + x96*M[60] + z*M[133] + M[184];
#pragma omp atomic
Ms[185] += x107*x170 + x11*x419 + x11*M[104] + x125*M[44] + x170*M[37] + x209*x86 + x214*M[29] + x234*M[25] + x291*x45 + x292*x46 + x292*M[20] + x34*x356 + x34*M[93] + x361*M[15] + x389*M[12] + x45*M[76] + x461*x9 + x461*M[9] + x563*M[6] + x600*M[4] + x790 + x791 + x792 + x793*x794 + x795*x796 + x86*M[61] + M[185];
#pragma omp atomic
Ms[186] += x11*x424 + x11*x425 + x11*x428 + x11*M[105] + x112*x170 + x114*x170 + x116*x170 + x125*x148 + x125*M[45] + x170*M[38] + x214*x88 + x214*M[30] + x215*x86 + x217*x86 + x219*x86 + x234*x71 + x234*M[26] + x25*x389 + x259*x800 + x292*x49 + x292*x50 + x292*x51 + x292*M[21] + x293*x45 + x294*x45 + x295*x45 + x34*x362 + x34*x364 + x34*x366 + x34*M[94] + x36*x361 + x361*M[16] + x389*M[13] + x4*x600 + x45*M[77] + x558*z + x559*z + x560*z + x563*x7 + x563*M[7] + x600*M[5] + x797 + x798 + x799 + x86*M[62] + z*M[134] + M[186];
#pragma omp atomic
Ms[187] += x107*x177 + x11*x430 + x11*x431 + x11*x434 + x11*M[106] + x118*x170 + x119*x170 + x120*x170 + x121*x142 + x121*M[44] + x125*x154 + x125*M[46] + x170*M[39] + x177*M[37] + x208*x85 + x208*M[29] + x214*x91 + x214*M[31] + x221*x86 + x223*x86 + x225*x86 + x234*x77 + x234*M[27] + x296*x45 + x297*x45 + x298*x45 + x299*x46 + x299*M[20] + x30*x389 + x34*x368 + x34*x370 + x34*x372 + x34*M[95] + x356*x41 + x357*x41 + x358*x41 + x360*x804 + x361*x39 + x361*M[17] + x389*M[14] + x402*M[12] + x41*M[93] + x45*M[78] + x468*x9 + x468*M[9] + x563*x8 + x563*M[8] + x564*z + x565*z + x567*z + x613*x776 + x756*M[6] + x775*M[4] + x801 + x802 + x803 + x86*M[63] + z*M[135] + M[187];
#pragma omp atomic
Ms[188] += x11*x437 + x11*x438 + x11*x441 + x11*M[107] + x112*x177 + x121*x147 + x121*M[45] + x125*x161 + x125*M[47] + x177*M[38] + x208*x87 + x208*M[30] + x209*x96 + x210*x96 + x211*x96 + x213*x808 + x214*x94 + x214*M[32] + x227*x86 + x228*x86 + x229*x86 + x230*x65 + x230*M[25] + x234*x82 + x234*M[28] + x299*x49 + x299*M[21] + x300*x45 + x301*x45 + x302*x45 + x303*x46 + x303*M[20] + x33*x355 + x34*x376 + x34*x378 + x34*x380 + x34*M[96] + x355*M[15] + x361*x42 + x361*M[18] + x362*x41 + x363*x41 + x365*x41 + x402*M[13] + x41*M[94] + x45*M[79] + x472*x9 + x472*M[9] + x570*z + x571*z + x573*z + x621*x809 + x623*x776 + x756*M[7] + x775*M[5] + x782*M[4] + x805 + x806 + x807 + x86*M[64] + x96*M[61] + z*M[136] + M[188];
#pragma omp atomic
Ms[189] += x107*x184 + x108*x184 + x109*x184 + x11*x444 + x11*x445 + x11*x448 + x11*M[108] + x118*x177 + x121*x153 + x121*M[46] + x125*x166 + x125*M[48] + x177*M[39] + x184*M[37] + x19*x385 + x208*x90 + x208*M[31] + x214*x97 + x214*M[33] + x215*x96 + x216*x96 + x218*x96 + x230*x70 + x230*M[26] + x303*x49 + x303*M[21] + x304*x45 + x305*x45 + x306*x45 + x34*x382 + x34*x383 + x34*x384 + x34*M[97] + x35*x355 + x355*M[16] + x368*x41 + x369*x41 + x371*x41 + x385*M[12] + x402*M[14] + x41*M[95] + x45*M[80] + x476*x9 + x476*M[9] + x5*x557 + x557*M[6] + x579*z + x580*z + x582*z + x756*M[8] + x757*x813 + x782*M[5] + x810 + x811 + x812 + x814*M[2] + x96*M[62] + z*M[137] + M[189];
#pragma omp atomic
Ms[190] += x11*x451 + x11*x452 + x11*x455 + x11*M[109] + x112*x184 + x113*x184 + x115*x184 + x121*x160 + x121*M[47] + x184*M[38] + x2*x596 + x208*x93 + x208*M[32] + x221*x96 + x222*x96 + x224*x96 + x230*x76 + x230*M[27] + x24*x385 + x307*x45 + x308*x45 + x309*x45 + x310*x46 + x310*x47 + x310*x48 + x310*M[20] + x355*x38 + x355*M[17] + x376*x41 + x377*x41 + x379*x41 + x385*M[13] + x41*M[96] + x45*M[81] + x557*x6 + x557*M[7] + x587*z + x588*z + x590*z + x596*M[4] + x815 + x816 + x817 + x818*M[1] + x96*M[63] + z*M[138] + M[190];
#pragma omp atomic
Ms[191] += x11*x457 + x11*M[110] + x118*x184 + x121*M[48] + x184*M[39] + x208*M[33] + x227*x96 + x230*M[28] + x290*x800 + x310*x49 + x310*M[21] + x311*x45 + x355*M[18] + x382*x41 + x385*M[14] + x41*M[97] + x45*M[82] + x483*x9 + x483*M[9] + x557*M[8] + x593*z + x596*M[5] + x819 + x820*M[0] + x96*M[64] + z*M[139] + M[191];
#pragma omp atomic
Ms[192] += x0*x689 + x11*x460 + x11*M[111] + x12*x461 + x122*x170 + x125*M[49] + x170*M[40] + x231*x86 + x234*M[29] + x292*x52 + x292*M[22] + x34*x386 + x34*M[98] + x389*M[15] + x461*M[10] + x600*M[6] + x689*M[3] + x821 + x822 + x823 + x824*M[1] + x86*M[65] + M[192];
#pragma omp atomic
Ms[193] += x11*x462 + x11*x463 + x11*x464 + x11*M[112] + x125*x172 + x125*M[50] + x126*x170 + x128*x170 + x130*x170 + x15*x461 + x16*x461 + x17*x461 + x170*M[41] + x234*x88 + x234*M[30] + x235*x86 + x237*x86 + x239*x86 + x292*x55 + x292*x57 + x292*x59 + x292*M[23] + x34*x390 + x34*x392 + x34*x394 + x34*M[99] + x36*x389 + x389*M[16] + x461*M[11] + x597*z + x598*z + x599*z + x600*x7 + x600*M[7] + x824*M[2] + x825 + x826 + x827 + x86*M[66] + z*M[140] + M[193];
#pragma omp atomic
Ms[194] += x0*x696 + x11*x465 + x11*x466 + x11*x467 + x11*M[113] + x12*x468 + x121*x169 + x121*M[49] + x122*x177 + x125*x175 + x125*M[51] + x132*x170 + x134*x170 + x136*x170 + x170*M[42] + x177*M[40] + x234*x91 + x234*M[31] + x241*x86 + x243*x86 + x245*x86 + x292*x61 + x292*x62 + x292*x63 + x292*M[24] + x299*x52 + x299*M[22] + x34*x396 + x34*x398 + x34*x400 + x34*M[100] + x386*x41 + x387*x41 + x388*x41 + x389*x39 + x389*M[17] + x402*M[15] + x41*M[98] + x468*M[10] + x600*x8 + x600*M[8] + x601*z + x602*z + x604*z + x696*M[3] + x775*M[6] + x777*x831 + x828 + x829 + x830 + x86*M[67] + z*M[141] + M[194];
#pragma omp atomic
Ms[195] += x0*x700 + x11*x469 + x11*x470 + x11*x471 + x11*M[114] + x12*x472 + x121*x171 + x121*M[50] + x125*x179 + x125*M[52] + x126*x177 + x138*x170 + x139*x170 + x140*x170 + x15*x468 + x170*M[43] + x177*M[41] + x230*x85 + x230*M[29] + x231*x96 + x232*x96 + x233*x96 + x234*x94 + x234*M[32] + x247*x86 + x249*x86 + x251*x86 + x299*x55 + x299*M[23] + x303*x52 + x303*M[22] + x34*x403 + x34*x405 + x34*x407 + x34*M[101] + x389*x42 + x389*M[18] + x390*x41 + x391*x41 + x393*x41 + x402*M[16] + x41*M[99] + x468*M[11] + x472*M[10] + x607*z + x608*z + x610*z + x613*x808 + x700*M[3] + x775*M[7] + x782*M[6] + x785*x831 + x832 + x833 + x834 + x86*M[68] + x96*M[65] + z*M[142] + M[195];
#pragma omp atomic
Ms[196] += x0*x704 + x11*x473 + x11*x474 + x11*x475 + x11*M[115] + x12*x476 + x121*x174 + x121*M[51] + x122*x184 + x123*x184 + x124*x184 + x125*x182 + x125*M[53] + x132*x177 + x15*x472 + x177*M[42] + x184*M[40] + x230*x87 + x230*M[30] + x234*x97 + x234*M[33] + x235*x96 + x236*x96 + x238*x96 + x253*x86 + x254*x86 + x255*x86 + x299*x61 + x299*M[24] + x303*x55 + x303*M[23] + x33*x385 + x34*x409 + x34*x411 + x34*x413 + x34*M[102] + x385*M[15] + x396*x41 + x397*x41 + x399*x41 + x402*M[17] + x41*M[100] + x472*M[11] + x476*M[10] + x615*z + x616*z + x618*z + x623*x808 + x704*M[3] + x775*M[8] + x777*x838 + x782*M[7] + x835 + x836 + x837 + x86*M[69] + x96*M[66] + z*M[143] + M[196];
#pragma omp atomic
Ms[197] += x0*x708 + x11*x477 + x11*x478 + x11*x479 + x11*M[116] + x121*x178 + x121*M[52] + x125*x185 + x125*M[54] + x126*x184 + x127*x184 + x129*x184 + x138*x177 + x15*x476 + x177*M[43] + x184*M[41] + x230*x90 + x230*M[31] + x241*x96 + x242*x96 + x244*x96 + x303*x61 + x303*M[24] + x310*x52 + x310*x53 + x310*x54 + x310*M[22] + x34*x415 + x34*x416 + x34*x417 + x34*M[103] + x35*x385 + x385*M[16] + x402*M[18] + x403*x41 + x404*x41 + x406*x41 + x41*M[101] + x476*M[11] + x5*x596 + x596*M[6] + x624*z + x625*z + x627*z + x708*M[3] + x782*M[8] + x785*x838 + x839 + x840 + x841 + x96*M[67] + z*M[144] + M[197];
#pragma omp atomic
Ms[198] += x11*x480 + x11*x481 + x11*x482 + x11*M[117] + x12*x483 + x121*x181 + x121*M[53] + x13*x483 + x132*x184 + x133*x184 + x135*x184 + x14*x483 + x184*M[42] + x230*x93 + x230*M[32] + x247*x96 + x248*x96 + x250*x96 + x310*x55 + x310*x56 + x310*x58 + x310*M[23] + x38*x385 + x385*M[17] + x409*x41 + x41*x410 + x41*x412 + x41*M[102] + x483*M[10] + x596*x6 + x596*M[7] + x630*z + x631*z + x633*z + x820*M[1] + x842 + x843 + x844 + x96*M[68] + z*M[145] + M[198];
#pragma omp atomic
Ms[199] += x0*x715 + x11*x484 + x11*M[118] + x121*M[54] + x138*x184 + x15*x483 + x184*M[43] + x230*M[33] + x253*x96 + x310*x61 + x310*M[24] + x385*M[18] + x41*x415 + x41*M[103] + x483*M[11] + x596*M[8] + x636*z + x715*M[3] + x820*M[2] + x845 + x96*M[69] + z*M[146] + M[199];
#pragma omp atomic
Ms[200] += x1*x689 + x141*x170 + x170*M[44] + x18*x461 + x256*x86 + x292*x64 + x292*M[25] + x34*x418 + x34*M[104] + x461*M[12] + x689*M[4] + x846 + x847 + x848 + x849*x850 + x86*M[70] + M[200];
#pragma omp atomic
Ms[201] += x146*x170 + x148*x170 + x150*x170 + x170*M[45] + x23*x461 + x25*x461 + x261*x86 + x263*x86 + x265*x86 + x27*x461 + x292*x69 + x292*x71 + x292*x73 + x292*M[26] + x3*x689 + x34*x423 + x34*x425 + x34*x427 + x34*M[105] + x4*x689 + x461*M[13] + x639*z + x640*z + x641*z + x689*M[5] + x851 + x852 + x853 + x86*M[71] + z*M[147] + M[201];
#pragma omp atomic
Ms[202] += x1*x696 + x141*x177 + x152*x170 + x154*x170 + x156*x170 + x170*M[46] + x177*M[44] + x18*x468 + x267*x86 + x269*x86 + x271*x86 + x29*x461 + x292*x75 + x292*x77 + x292*x79 + x292*M[27] + x299*x64 + x299*M[25] + x30*x461 + x31*x461 + x34*x429 + x34*x431 + x34*x433 + x34*M[106] + x41*x418 + x41*x419 + x41*x420 + x41*M[104] + x461*M[14] + x468*M[12] + x644*z + x645*z + x647*z + x696*M[4] + x794*x857 + x854 + x855 + x856 + x86*M[72] + z*M[148] + M[202];
#pragma omp atomic
Ms[203] += x1*x700 + x146*x177 + x159*x170 + x161*x170 + x163*x170 + x170*M[47] + x177*M[45] + x18*x472 + x23*x468 + x256*x96 + x257*x96 + x258*x96 + x274*x86 + x276*x86 + x278*x86 + x292*x81 + x292*x82 + x292*x83 + x292*M[28] + x299*x69 + x299*M[26] + x3*x696 + x303*x64 + x303*M[25] + x34*x436 + x34*x438 + x34*x440 + x34*M[107] + x41*x423 + x41*x424 + x41*x426 + x41*M[105] + x422*x861 + x468*M[13] + x472*M[12] + x650*z + x651*z + x653*z + x696*M[5] + x700*M[4] + x858 + x859 + x86*M[73] + x860 + x96*M[70] + z*M[149] + M[203];
#pragma omp atomic
Ms[204] += x1*x704 + x141*x184 + x142*x184 + x143*x184 + x152*x177 + x165*x170 + x166*x170 + x167*x170 + x170*M[48] + x177*M[46] + x18*x476 + x184*M[44] + x23*x472 + x261*x96 + x262*x96 + x264*x96 + x281*x86 + x283*x86 + x285*x86 + x29*x468 + x299*x75 + x299*M[27] + x3*x700 + x303*x69 + x303*M[26] + x34*x443 + x34*x445 + x34*x447 + x34*M[108] + x360*x865 + x41*x429 + x41*x430 + x41*x432 + x41*M[106] + x468*M[14] + x472*M[13] + x476*M[12] + x657*z + x658*z + x660*z + x700*M[5] + x704*M[4] + x86*M[74] + x862 + x863 + x864 + x96*M[71] + z*M[150] + M[204];
#pragma omp atomic
Ms[205] += x1*x708 + x146*x184 + x147*x184 + x149*x184 + x159*x177 + x177*M[47] + x184*M[45] + x23*x476 + x267*x96 + x268*x96 + x270*x96 + x287*x86 + x288*x86 + x289*x86 + x29*x472 + x299*x81 + x299*M[28] + x3*x704 + x303*x75 + x303*M[27] + x310*x64 + x310*x65 + x310*x66 + x310*M[25] + x34*x450 + x34*x452 + x34*x454 + x34*M[109] + x41*x436 + x41*x437 + x41*x439 + x41*M[107] + x422*x869 + x472*M[14] + x476*M[13] + x664*z + x665*z + x667*z + x704*M[5] + x708*M[4] + x86*M[75] + x866 + x867 + x868 + x96*M[72] + z*M[151] + M[205];
#pragma omp atomic
Ms[206] += x111*x873 + x152*x184 + x153*x184 + x155*x184 + x165*x177 + x177*M[48] + x18*x483 + x184*M[46] + x19*x483 + x20*x483 + x274*x96 + x275*x96 + x277*x96 + x29*x476 + x3*x708 + x303*x81 + x303*M[28] + x310*x69 + x310*x70 + x310*x72 + x310*M[26] + x34*x456 + x34*x457 + x34*x458 + x34*M[110] + x41*x443 + x41*x444 + x41*x446 + x41*M[108] + x476*M[14] + x483*M[12] + x671*z + x672*z + x674*z + x708*M[5] + x870 + x871 + x872 + x96*M[73] + z*M[152] + M[206];
#pragma omp atomic
Ms[207] += x1*x715 + x159*x184 + x160*x184 + x162*x184 + x184*M[47] + x2*x715 + x23*x483 + x24*x483 + x26*x483 + x281*x96 + x282*x96 + x284*x96 + x310*x75 + x310*x76 + x310*x78 + x310*M[27] + x41*x450 + x41*x451 + x41*x453 + x41*M[109] + x483*M[13] + x678*z + x679*z + x681*z + x715*M[4] + x874 + x875 + x876 + x96*M[74] + z*M[153] + M[207];
#pragma omp atomic
Ms[208] += x165*x184 + x184*M[48] + x287*x96 + x29*x483 + x3*x715 + x310*x81 + x310*M[28] + x41*x456 + x41*M[110] + x483*M[14] + x684*z + x715*M[5] + x850*x878 + x877 + x96*M[75] + z*M[154] + M[208];
#pragma omp atomic
Ms[209] += x170*M[49] + x292*M[29] + x34*M[111] + x461*M[15] + x689*M[6] + x86*M[76] + x879 + x880*M[1] + M[209];
#pragma omp atomic
Ms[210] += x170*x172 + x170*M[50] + x292*x88 + x292*M[30] + x294*x86 + x34*x463 + x34*M[112] + x36*x461 + x461*M[16] + x688*z + x689*x7 + x689*M[7] + x86*M[77] + x880*M[2] + x881 + z*M[155] + M[210];
#pragma omp atomic
Ms[211] += x170*x175 + x170*M[51] + x177*M[49] + x292*x91 + x292*M[31] + x297*x86 + x299*M[29] + x34*x466 + x34*M[113] + x39*x461 + x41*x460 + x41*M[111] + x461*M[17] + x468*M[15] + x689*x8 + x689*M[8] + x690*z + x696*M[6] + x86*M[78] + x882 + x883*M[1] + z*M[156] + M[211];
#pragma omp atomic
Ms[212] += x170*x179 + x170*M[52] + x177*M[50] + (1.0/720.0)*x259*x621 + x291*x96 + x292*x94 + x292*M[32] + x299*M[30] + x301*x86 + x303*M[29] + x34*x470 + x34*M[114] + x41*x462 + x41*M[112] + x42*x461 + x461*M[18] + x468*M[16] + x472*M[15] + x693*z + x696*M[7] + x700*M[6] + x86*M[79] + x883*M[2] + x884 + x96*M[76] + z*M[157] + M[212];
#pragma omp atomic
Ms[213] += x169*x184 + x170*x182 + x170*M[53] + x177*M[51] + x184*M[49] + x292*x97 + x292*M[33] + x293*x96 + x299*M[31] + x303*M[30] + x305*x86 + x34*x474 + x34*M[115] + x41*x465 + x41*M[113] + x468*M[17] + x472*M[16] + x476*M[15] + x696*M[8] + x697*z + x700*M[7] + x704*M[6] + x86*M[80] + x861*x887 + x885 + x886*M[1] + x96*M[77] + z*M[158] + M[213];
#pragma omp atomic
Ms[214] += x170*x185 + x170*M[54] + x171*x184 + x177*M[52] + x184*M[50] + (1.0/720.0)*x290*x613 + x296*x96 + x299*M[32] + x303*M[31] + x308*x86 + x310*x85 + x310*M[29] + x34*x478 + x34*M[116] + x41*x469 + x41*M[114] + x468*M[18] + x472*M[17] + x476*M[16] + x700*M[8] + x701*z + x704*M[7] + x708*M[6] + x86*M[81] + x886*M[2] + x888 + x96*M[78] + z*M[159] + M[214];
#pragma omp atomic
Ms[215] += x174*x184 + x177*M[53] + x184*M[51] + x299*M[33] + x300*x96 + x303*M[32] + x310*x87 + x310*M[30] + x311*x86 + x33*x483 + x34*x481 + x34*M[117] + x41*x473 + x41*M[115] + x472*M[18] + x476*M[17] + x483*M[15] + x704*M[8] + x705*z + x708*M[7] + x86*M[82] + x869*x887 + x889 + x890*M[1] + x96*M[79] + z*M[160] + M[215];
#pragma omp atomic
Ms[216] += x177*M[54] + x178*x184 + x184*M[52] + x303*M[33] + x304*x96 + x310*x90 + x310*M[31] + x34*x484 + x34*M[118] + x35*x483 + x41*x477 + x41*M[116] + x476*M[18] + x483*M[16] + x5*x715 + x708*M[8] + x709*z + x715*M[6] + x890*M[2] + x891 + x96*M[80] + z*M[161] + M[216];
#pragma omp atomic
Ms[217] += x181*x184 + x184*M[53] + x307*x96 + x310*x93 + x310*M[32] + x38*x483 + x41*x480 + x41*M[117] + x483*M[17] + x6*x715 + x712*z + x715*M[7] + x892 + x893*M[1] + x96*M[81] + z*M[162] + M[217];
#pragma omp atomic
Ms[218] += x184*M[54] + x310*M[33] + x41*M[118] + x483*M[18] + x715*M[8] + x893*M[2] + x96*M[82] + z*M[163] + M[218];
#pragma omp atomic
Ms[219] += x*M[164] + x100*M[55] + x11*M[119] + x188*M[34] + x314*M[19] + x45*M[83] + x487*M[9] + x718*M[3] + x894*M[0] + M[219];
#pragma omp atomic
Ms[220] += x*M[165] + x100*x190 + x100*M[56] + x102*x188 + x11*x489 + x11*M[120] + x13*x487 + x188*M[35] + x2*x718 + x314*x47 + x314*M[20] + x316*x45 + x45*M[84] + x487*M[10] + x717*y + x718*M[4] + x894*M[1] + y*M[164] + M[220];
#pragma omp atomic
Ms[221] += x*M[166] + x100*x193 + x100*M[57] + x105*x188 + x11*x492 + x11*M[121] + x16*x487 + x188*M[36] + x314*x50 + x314*M[21] + x319*x45 + x4*x718 + x45*M[85] + x487*M[11] + x717*z + x718*M[5] + x894*M[2] + z*M[164] + M[221];
#pragma omp atomic
Ms[222] += x*M[167] + x100*x196 + x100*M[58] + x108*x188 + x11*x495 + x11*M[122] + x111*x895 + x125*M[55] + x188*M[37] + x19*x487 + x214*M[34] + x314*x53 + x314*M[22] + x322*x45 + x339*M[19] + x34*x485 + x34*M[119] + x45*M[86] + x487*M[12] + x5*x718 + x512*M[9] + x718*M[6] + x719*y + x736*M[3] + y*M[165] + M[222];
#pragma omp atomic
Ms[223] += x*M[168] + x100*x200 + x100*x201 + x100*x204 + x100*M[59] + x11*x499 + x11*x500 + x11*x503 + x11*M[123] + x113*x188 + x114*x188 + x117*x188 + x188*M[38] + x24*x487 + x25*x487 + x28*x487 + x314*x56 + x314*x57 + x314*x60 + x314*M[23] + x326*x45 + x327*x45 + x330*x45 + x45*M[87] + x487*M[13] + x6*x718 + x7*x718 + x718*M[7] + x719*z + x720*z + x721*z + x722*y + y*M[166] + z*M[165] + M[223];
#pragma omp atomic
Ms[224] += x*M[169] + x100*x206 + x100*M[60] + x11*x505 + x11*M[124] + x119*x188 + x121*M[55] + x188*M[39] + x208*M[34] + x30*x487 + x314*x62 + x314*M[24] + x332*x45 + x334*M[19] + x41*x485 + x41*M[119] + x45*M[88] + x487*M[14] + x507*M[9] + x718*x8 + x718*M[8] + x722*z + x731*M[3] + x896*M[0] + z*M[166] + M[224];
#pragma omp atomic
Ms[225] += x*M[170] + x100*x210 + x100*M[61] + x11*x509 + x11*M[125] + x123*x188 + x125*M[56] + x188*M[40] + x213*x898 + x214*M[35] + x234*M[34] + x312*x86 + x314*x65 + x314*M[25] + x33*x487 + x336*x45 + x339*M[20] + x34*x488 + x34*M[120] + x361*M[19] + x45*M[89] + x487*M[15] + x512*M[10] + x533*M[9] + x67*x897 + x723*y + x736*M[4] + x86*M[83] + x899*M[1] + y*M[167] + M[225];
#pragma omp atomic
Ms[226] += x*M[171] + x100*x216 + x100*x217 + x100*x220 + x100*M[62] + x105*x214 + x11*x514 + x11*x515 + x11*x518 + x11*M[126] + x125*x193 + x125*M[57] + x127*x188 + x128*x188 + x131*x188 + x16*x512 + x188*M[41] + x214*M[36] + x314*x70 + x314*x71 + x314*x74 + x314*M[26] + x339*x50 + x339*M[21] + x34*x491 + x34*x492 + x34*x493 + x34*M[121] + x341*x45 + x342*x45 + x345*x45 + x35*x487 + x36*x487 + x37*x487 + x4*x736 + x45*M[90] + x487*M[16] + x512*M[11] + x723*z + x724*z + x725*z + x727*y + x736*M[5] + x899*M[2] + y*M[168] + z*M[167] + M[226];
#pragma omp atomic
Ms[227] += x*M[172] + x100*x222 + x100*x223 + x100*x226 + x100*M[63] + x102*x208 + x11*x520 + x11*x521 + x11*x524 + x11*M[127] + x121*x190 + x121*M[56] + x13*x507 + x133*x188 + x134*x188 + x137*x188 + x188*M[42] + x2*x731 + x208*M[35] + x314*x76 + x314*x77 + x314*x80 + x314*M[27] + x334*x47 + x334*M[20] + x347*x45 + x348*x45 + x351*x45 + x38*x487 + x39*x487 + x40*x487 + x41*x488 + x41*x489 + x41*x490 + x41*M[120] + x45*M[91] + x487*M[17] + x507*M[10] + x727*z + x728*z + x729*z + x730*y + x731*M[4] + x896*M[1] + y*M[169] + z*M[168] + M[227];
#pragma omp atomic
Ms[228] += x*M[173] + x100*x228 + x100*M[64] + x11*x526 + x11*M[128] + x121*M[57] + x139*x188 + x188*M[43] + x208*M[36] + x230*M[34] + x312*x96 + x314*x82 + x314*M[28] + x334*M[21] + x353*x45 + x355*M[19] + x41*x491 + x41*M[121] + x42*x487 + x45*M[92] + x487*M[18] + x507*M[11] + x528*M[9] + x730*z + x731*M[5] + x84*x897 + x896*M[2] + x900*M[0] + x96*M[83] + z*M[169] + M[228];
#pragma omp atomic
Ms[229] += x*M[174] + x100*x232 + x100*M[65] + x11*x530 + x11*M[129] + x125*M[58] + x142*x188 + x170*x186 + x170*M[55] + x188*M[44] + x214*M[37] + x234*M[35] + x314*x85 + x314*M[29] + x315*x86 + x339*M[22] + x34*x494 + x34*M[122] + x357*x45 + x360*x902 + x361*M[20] + x389*M[19] + x45*M[93] + x512*M[12] + x533*M[10] + x563*M[9] + x613*x898 + x732*y + x736*M[6] + x768*M[3] + x86*M[84] + x901*M[4] + y*M[170] + M[229];
#pragma omp atomic
Ms[230] += x*M[175] + x100*x236 + x100*x237 + x100*x240 + x100*M[66] + x105*x234 + x11*x535 + x11*x536 + x11*x539 + x11*M[130] + x114*x214 + x125*x201 + x125*M[59] + x147*x188 + x148*x188 + x151*x188 + x16*x533 + x188*M[45] + x214*M[38] + x234*M[36] + x25*x512 + x314*x87 + x314*x88 + x314*x89 + x314*M[30] + x318*x86 + x319*x86 + x320*x86 + x339*x57 + x339*M[23] + x34*x498 + x34*x500 + x34*x502 + x34*M[123] + x361*x50 + x361*M[21] + x363*x45 + x364*x45 + x367*x45 + x4*x901 + x45*M[94] + x512*M[13] + x533*M[11] + x623*x898 + x7*x736 + x732*z + x733*z + x734*z + x736*M[7] + x737*y + x86*M[85] + x901*M[5] + y*M[171] + z*M[170] + M[230];
#pragma omp atomic
Ms[231] += x*M[176] + x100*x242 + x100*x243 + x100*x246 + x100*M[67] + x108*x208 + x11*x541 + x11*x542 + x11*x545 + x11*M[131] + x119*x214 + x121*x196 + x121*M[58] + x125*x206 + x125*M[60] + x153*x188 + x154*x188 + x157*x188 + x177*x186 + x177*M[55] + x188*M[46] + x19*x507 + x208*M[37] + x214*M[39] + x30*x512 + x314*x90 + x314*x91 + x314*x92 + x314*M[31] + x334*x53 + x334*M[22] + x339*x62 + x339*M[24] + x34*x504 + x34*x505 + x34*x506 + x34*M[124] + x369*x45 + x370*x45 + x373*x45 + x402*M[19] + x41*x494 + x41*x495 + x41*x496 + x41*M[122] + x45*M[95] + x5*x731 + x507*M[12] + x512*M[14] + x731*M[6] + x736*x8 + x736*M[8] + x737*z + x738*z + x739*z + x740*y + x756*M[9] + x903*M[3] + x904*x905 + y*M[172] + z*M[171] + M[231];
#pragma omp atomic
Ms[232] += x*M[177] + x100*x248 + x100*x249 + x100*x252 + x100*M[68] + x102*x230 + x11*x549 + x11*x550 + x11*x553 + x11*M[132] + x113*x208 + x121*x200 + x121*M[59] + x13*x528 + x160*x188 + x161*x188 + x164*x188 + x188*M[47] + x2*x906 + x208*M[38] + x230*M[35] + x24*x507 + x314*x93 + x314*x94 + x314*x95 + x314*M[32] + x315*x96 + x316*x96 + x317*x96 + x334*x56 + x334*M[23] + x355*x47 + x355*M[20] + x377*x45 + x378*x45 + x381*x45 + x41*x498 + x41*x499 + x41*x501 + x41*M[123] + x45*M[96] + x507*M[13] + x528*M[10] + x6*x731 + x621*x898 + x731*M[7] + x740*z + x741*z + x742*z + x743*y + x906*M[4] + x96*M[84] + y*M[173] + z*M[172] + M[232];
#pragma omp atomic
Ms[233] += x*M[178] + x100*x254 + x100*M[69] + x11*x555 + x11*M[133] + x121*M[60] + x166*x188 + x184*x186 + x184*M[55] + x188*M[48] + x208*M[39] + x230*M[36] + x314*x97 + x314*M[33] + x318*x96 + x334*M[24] + x355*M[21] + x383*x45 + x385*M[19] + x41*x504 + x41*M[124] + x45*M[97] + x507*M[14] + x528*M[11] + x557*M[9] + x731*M[8] + x743*z + x763*M[3] + x900*M[2] + x906*M[5] + x907*M[0] + x96*M[85] + z*M[173] + M[233];
#pragma omp atomic
Ms[234] += x*M[179] + x100*x257 + x100*M[70] + x11*x559 + x11*M[134] + x125*M[61] + x169*x188 + x170*x189 + x170*M[56] + x188*M[49] + x214*M[40] + x234*M[37] + x292*x98 + x292*M[34] + x321*x86 + x339*M[25] + x34*x508 + x34*M[125] + x361*M[22] + x387*x45 + x389*M[20] + x45*M[98] + x512*M[15] + x533*M[12] + x562*x909 + x563*M[10] + x600*M[9] + x744*y + x768*M[4] + x86*M[86] + x901*M[6] + x908*M[3] + x910*M[1] + y*M[174] + M[234];
#pragma omp atomic
Ms[235] += x*M[180] + x100*x262 + x100*x263 + x100*x266 + x100*M[71] + x11*x565 + x11*x566 + x11*x569 + x11*M[135] + x114*x234 + x125*x217 + x125*M[62] + x128*x214 + x16*x563 + x170*x192 + x170*x193 + x170*x194 + x170*M[57] + x171*x188 + x172*x188 + x173*x188 + x188*M[50] + x214*M[41] + x234*M[38] + x25*x533 + x325*x86 + x327*x86 + x329*x86 + x339*x71 + x339*M[26] + x34*x513 + x34*x515 + x34*x517 + x34*M[126] + x36*x512 + x361*x57 + x361*M[23] + x389*x50 + x389*M[21] + x391*x45 + x392*x45 + x395*x45 + x4*x768 + x45*M[99] + x512*M[16] + x533*M[13] + x563*M[11] + x7*x901 + x744*z + x745*z + x746*z + x749*y + x768*M[5] + x86*M[87] + x901*M[7] + x910*M[2] + y*M[175] + z*M[174] + M[235];
#pragma omp atomic
Ms[236] += x*M[181] + x100*x268 + x100*x269 + x100*x272 + x100*M[72] + x11*x571 + x11*x572 + x11*x575 + x11*M[136] + x119*x234 + x121*x210 + x121*M[61] + x123*x208 + x125*x223 + x125*M[63] + x134*x214 + x174*x188 + x175*x188 + x176*x188 + x177*x189 + x177*M[56] + x188*M[51] + x208*M[40] + x213*x912 + x214*M[42] + x234*M[39] + x299*x98 + x299*M[34] + x30*x533 + x33*x507 + x331*x86 + x332*x86 + x333*x86 + x334*x65 + x334*M[25] + x339*x77 + x339*M[27] + x34*x519 + x34*x521 + x34*x523 + x34*M[127] + x361*x62 + x361*M[24] + x39*x512 + x397*x45 + x398*x45 + x401*x45 + x402*M[20] + x41*x508 + x41*x509 + x41*x510 + x41*M[125] + x45*M[100] + x507*M[15] + x512*M[17] + x533*M[14] + x749*z + x750*z + x751*z + x753*y + x756*M[10] + x775*M[9] + x8*x901 + x86*M[88] + x901*M[8] + x903*M[4] + x911*M[3] + x913*x914 + y*M[176] + z*M[175] + M[236];
#pragma omp atomic
Ms[237] += x*M[182] + x100*x275 + x100*x276 + x100*x279 + x100*M[73] + x108*x230 + x11*x580 + x11*x581 + x11*x584 + x11*M[137] + x121*x216 + x121*M[62] + x125*x228 + x125*M[64] + x127*x208 + x139*x214 + x177*x192 + x177*M[57] + x178*x188 + x179*x188 + x180*x188 + x188*M[52] + x19*x528 + x208*M[41] + x214*M[43] + x230*M[37] + x303*x98 + x303*M[34] + x321*x96 + x322*x96 + x323*x96 + x334*x70 + x334*M[26] + x339*x82 + x339*M[28] + x34*x525 + x34*x526 + x34*x527 + x34*M[128] + x35*x507 + x355*x53 + x355*M[22] + x402*M[21] + x404*x45 + x405*x45 + x408*x45 + x41*x513 + x41*x514 + x41*x516 + x41*M[126] + x42*x512 + x45*M[101] + x5*x906 + x507*M[16] + x512*M[18] + x528*M[12] + x753*z + x754*z + x755*z + x756*M[11] + x759*y + x782*M[9] + x784*x915 + x814*M[3] + x903*M[5] + x906*M[6] + x914*x916 + x96*M[86] + y*M[177] + z*M[176] + M[237];
#pragma omp atomic
Ms[238] += x*M[183] + x100*x282 + x100*x283 + x100*x286 + x100*M[74] + x11*x588 + x11*x589 + x11*x592 + x11*M[138] + x113*x230 + x121*x222 + x121*M[63] + x13*x557 + x133*x208 + x181*x188 + x182*x188 + x183*x188 + x184*x189 + x184*x190 + x184*x191 + x184*M[56] + x188*M[53] + x2*x763 + x208*M[42] + x230*M[38] + x24*x528 + x325*x96 + x326*x96 + x328*x96 + x334*x76 + x334*M[27] + x355*x56 + x355*M[23] + x38*x507 + x385*x47 + x385*M[20] + x41*x519 + x41*x520 + x41*x522 + x41*M[127] + x410*x45 + x411*x45 + x414*x45 + x45*M[102] + x507*M[17] + x528*M[13] + x557*M[10] + x6*x906 + x759*z + x760*z + x761*z + x762*y + x763*M[4] + x906*M[7] + x907*M[1] + x96*M[87] + y*M[178] + z*M[177] + M[238];
#pragma omp atomic
Ms[239] += x*M[184] + x100*x288 + x100*M[75] + x11*x594 + x11*M[139] + x121*M[64] + x184*x192 + x184*M[57] + x185*x188 + x188*M[54] + x208*M[43] + x230*M[39] + x310*x98 + x310*M[34] + x331*x96 + x334*M[28] + x355*M[24] + x385*M[21] + x41*x525 + x41*M[128] + x416*x45 + x45*M[103] + x507*M[18] + x528*M[14] + x557*M[11] + x596*M[9] + x762*z + x763*M[5] + x818*M[3] + x906*M[8] + x907*M[2] + x917*M[0] + x96*M[88] + z*M[178] + M[239];
#pragma omp atomic
Ms[240] += x*M[185] + x100*x291 + x100*M[76] + x101*x292 + x11*x598 + x11*M[140] + x125*M[65] + x170*x195 + x170*M[58] + x214*M[44] + x234*M[40] + x292*M[35] + x335*x86 + x339*M[29] + x34*x529 + x34*M[129] + x361*M[25] + x389*M[22] + x419*x45 + x43*x461 + x45*M[104] + x461*M[19] + x533*M[15] + x563*M[12] + x600*M[10] + x764*y + x768*M[6] + x794*x918 + x796*x909 + x824*M[3] + x86*M[89] + x908*M[4] + y*M[179] + M[240];
#pragma omp atomic
Ms[241] += x*M[186] + x100*x293 + x100*x294 + x100*x295 + x100*M[77] + x104*x292 + x105*x292 + x106*x292 + x11*x602 + x11*x603 + x11*x606 + x11*M[141] + x125*x237 + x125*M[66] + x128*x234 + x148*x214 + x16*x600 + x170*x199 + x170*x201 + x170*x203 + x170*M[59] + x214*M[45] + x234*M[41] + x25*x563 + x292*M[36] + x339*x88 + x339*M[30] + x34*x534 + x34*x536 + x34*x538 + x34*M[130] + x340*x86 + x342*x86 + x344*x86 + x36*x533 + x361*x71 + x361*M[26] + x389*x57 + x389*M[23] + x4*x908 + x424*x45 + x425*x45 + x428*x45 + x45*M[105] + x533*M[16] + x563*M[13] + x600*M[11] + x7*x768 + x764*z + x765*z + x766*z + x768*M[7] + x769*y + x86*M[90] + x908*M[5] + x909*x919 + y*M[180] + z*M[179] + M[241];
#pragma omp atomic
Ms[242] += x*M[187] + x100*x296 + x100*x297 + x100*x298 + x100*M[78] + x101*x299 + x11*x608 + x11*x609 + x11*x612 + x11*M[142] + x121*x232 + x121*M[65] + x125*x243 + x125*M[67] + x134*x234 + x142*x208 + x154*x214 + x170*x205 + x170*x206 + x170*x207 + x170*M[60] + x177*x195 + x177*M[58] + x208*M[44] + x214*M[46] + x234*M[42] + x299*M[35] + x30*x563 + x334*x85 + x334*M[29] + x339*x91 + x339*M[31] + x34*x540 + x34*x542 + x34*x544 + x34*M[131] + x346*x86 + x348*x86 + x350*x86 + x361*x77 + x361*M[27] + x389*x62 + x389*M[24] + x39*x533 + x402*M[22] + x41*x529 + x41*x530 + x41*x531 + x41*M[129] + x43*x468 + x430*x45 + x431*x45 + x434*x45 + x45*M[106] + x468*M[19] + x533*M[17] + x547*x921 + x563*M[14] + x613*x912 + x756*M[12] + x768*x8 + x768*M[8] + x769*z + x770*z + x771*z + x772*y + x775*M[10] + x86*M[91] + x903*M[6] + x911*M[4] + x920*M[3] + y*M[181] + z*M[180] + M[242];
#pragma omp atomic
Ms[243] += x*M[188] + x100*x300 + x100*x301 + x100*x302 + x100*M[79] + x101*x303 + x104*x299 + x11*x616 + x11*x617 + x11*x620 + x11*M[143] + x121*x236 + x121*M[66] + x123*x230 + x125*x249 + x125*M[68] + x139*x234 + x147*x208 + x161*x214 + x177*x199 + x177*M[59] + x208*M[45] + x21*x621*x915 + x213*x923 + x214*M[47] + x230*M[40] + x234*M[43] + x299*M[36] + x303*M[35] + x33*x528 + x334*x87 + x334*M[30] + x335*x96 + x336*x96 + x337*x96 + x339*x94 + x339*M[32] + x34*x548 + x34*x550 + x34*x552 + x34*M[132] + x352*x86 + x353*x86 + x354*x86 + x355*x65 + x355*M[25] + x361*x82 + x361*M[28] + x402*M[23] + x41*x534 + x41*x535 + x41*x537 + x41*M[130] + x42*x533 + x43*x472 + x437*x45 + x438*x45 + x441*x45 + x45*M[107] + x472*M[19] + x528*M[15] + x533*M[18] + x623*x912 + x756*M[13] + x772*z + x773*z + x774*z + x775*M[11] + x779*y + x782*M[10] + x814*M[4] + x86*M[92] + x903*M[7] + x911*M[5] + x922*M[3] + x96*M[89] + y*M[182] + z*M[181] + M[243];
#pragma omp atomic
Ms[244] += x*M[189] + x100*x304 + x100*x305 + x100*x306 + x100*M[80] + x104*x303 + x11*x625 + x11*x626 + x11*x629 + x11*M[144] + x111*x925 + x121*x242 + x121*M[67] + x125*x254 + x125*M[69] + x127*x230 + x153*x208 + x166*x214 + x177*x205 + x177*M[60] + x184*x195 + x184*x196 + x184*x197 + x184*M[58] + x19*x557 + x208*M[46] + x214*M[48] + x230*M[41] + x303*M[36] + x334*x90 + x334*M[31] + x339*x97 + x339*M[33] + x34*x554 + x34*x555 + x34*x556 + x34*M[133] + x340*x96 + x341*x96 + x343*x96 + x35*x528 + x355*x70 + x355*M[26] + x385*x53 + x385*M[22] + x402*M[24] + x41*x540 + x41*x541 + x41*x543 + x41*M[131] + x43*x476 + x444*x45 + x445*x45 + x448*x45 + x45*M[108] + x476*M[19] + x5*x763 + x528*M[16] + x557*M[12] + x756*M[14] + x763*M[6] + x779*z + x780*z + x781*z + x782*M[11] + x786*y + x814*M[5] + x84*x915*x916 + x903*M[8] + x924*M[3] + x96*M[90] + y*M[183] + z*M[182] + M[244];
#pragma omp atomic
Ms[245] += x*M[190] + x100*x307 + x100*x308 + x100*x309 + x100*M[81] + x101*x310 + x102*x310 + x103*x310 + x11*x631 + x11*x632 + x11*x635 + x11*M[145] + x121*x248 + x121*M[68] + x13*x596 + x133*x230 + x160*x208 + x184*x199 + x184*x200 + x184*x202 + x184*M[59] + x2*x818 + x208*M[47] + x230*M[42] + x24*x557 + x310*M[35] + x334*x93 + x334*M[32] + x346*x96 + x347*x96 + x349*x96 + x355*x76 + x355*M[27] + x38*x528 + x385*x56 + x385*M[23] + x41*x548 + x41*x549 + x41*x551 + x41*M[132] + x45*x451 + x45*x452 + x45*x455 + x45*M[109] + x528*M[17] + x557*M[13] + x596*M[10] + x6*x763 + x763*M[7] + x786*z + x787*z + x788*z + x789*y + x818*M[4] + x917*M[1] + x96*M[91] + y*M[184] + z*M[183] + M[245];
#pragma omp atomic
Ms[246] += x*M[191] + x100*x311 + x100*M[82] + x104*x310 + x11*x637 + x11*M[146] + x121*M[69] + x184*x205 + x184*M[60] + x208*M[48] + x230*M[43] + x310*M[36] + x334*M[33] + x352*x96 + x355*M[28] + x385*M[24] + x41*x554 + x41*M[133] + x43*x483 + x45*x457 + x45*M[110] + x483*M[19] + x528*M[18] + x557*M[14] + x596*M[11] + x763*M[8] + x789*z + x818*M[5] + x820*M[3] + x917*M[2] + x926*M[0] + x96*M[92] + z*M[184] + M[246];
#pragma omp atomic
Ms[247] += x*M[192] + x107*x292 + x11*x640 + x11*M[147] + x125*M[70] + x170*x209 + x170*M[61] + x214*M[49] + x234*M[44] + x292*M[37] + x34*x558 + x34*M[134] + x356*x86 + x361*M[29] + x389*M[25] + x45*x460 + x45*M[111] + x46*x461 + x461*M[20] + x563*M[15] + x600*M[12] + x642*x928 + x689*x9 + x689*M[9] + x790*y + x824*M[4] + x86*M[93] + x908*M[6] + x929*M[1] + y*M[185] + M[247];
#pragma omp atomic
Ms[248] += x*M[193] + x11*x645 + x11*x646 + x11*x649 + x11*M[148] + x112*x292 + x114*x292 + x116*x292 + x125*x263 + x125*M[71] + x148*x234 + x170*x215 + x170*x217 + x170*x219 + x170*M[62] + x172*x214 + x214*M[50] + x234*M[45] + x25*x600 + x292*M[38] + x34*x564 + x34*x566 + x34*x568 + x34*M[135] + x36*x563 + x361*x88 + x361*M[30] + x362*x86 + x364*x86 + x366*x86 + x389*x71 + x389*M[26] + x4*x824 + x45*x462 + x45*x463 + x45*x464 + x45*M[112] + x461*x49 + x461*x50 + x461*x51 + x461*M[21] + x563*M[16] + x600*M[13] + x7*x908 + x790*z + x791*z + x792*z + x797*y + x824*M[5] + x86*M[94] + x908*M[7] + x929*M[2] + y*M[186] + z*M[185] + M[248];
#pragma omp atomic
Ms[249] += x*M[194] + x107*x299 + x11*x651 + x11*x652 + x11*x655 + x11*M[149] + x118*x292 + x119*x292 + x120*x292 + x121*x257 + x121*M[70] + x125*x269 + x125*M[72] + x154*x234 + x169*x208 + x170*x221 + x170*x223 + x170*x225 + x170*M[63] + x175*x214 + x177*x209 + x177*M[61] + x208*M[49] + x214*M[51] + x234*M[46] + x292*M[39] + x299*M[37] + x30*x600 + x34*x570 + x34*x572 + x34*x574 + x34*M[136] + x361*x91 + x361*M[31] + x368*x86 + x370*x86 + x372*x86 + x389*x77 + x389*M[27] + x39*x563 + x402*M[25] + x41*x558 + x41*x559 + x41*x560 + x41*M[134] + x45*x465 + x45*x466 + x45*x467 + x45*M[113] + x46*x468 + x468*M[20] + x562*x930 + x563*M[17] + x578*x931 + x600*M[14] + x696*x9 + x696*M[9] + x756*M[15] + x775*M[12] + x797*z + x798*z + x799*z + x8*x908 + x801*y + x86*M[95] + x908*M[8] + x911*M[6] + x920*M[4] + y*M[187] + z*M[186] + M[249];
#pragma omp atomic
Ms[250] += x*M[195] + x107*x303 + x11*x658 + x11*x659 + x11*x662 + x11*M[150] + x112*x299 + x121*x262 + x121*M[71] + x125*x276 + x125*M[73] + x142*x230 + x161*x234 + x170*x227 + x170*x228 + x170*x229 + x170*M[64] + x171*x208 + x177*x215 + x177*M[62] + x179*x214 + x208*M[50] + x214*M[52] + x230*M[44] + x234*M[47] + x299*M[38] + x303*M[37] + x34*x579 + x34*x581 + x34*x583 + x34*M[137] + x355*x85 + x355*M[29] + x356*x96 + x357*x96 + x358*x96 + x361*x94 + x361*M[32] + x376*x86 + x378*x86 + x380*x86 + x389*x82 + x389*M[28] + x402*M[26] + x41*x564 + x41*x565 + x41*x567 + x41*M[135] + x42*x563 + x45*x469 + x45*x470 + x45*x471 + x45*M[114] + x46*x472 + x468*x49 + x468*M[21] + x472*M[20] + x547*x932 + x563*M[18] + x585*x921 + x613*x923 + x700*x9 + x700*M[9] + x756*M[16] + x775*M[13] + x782*M[12] + x801*z + x802*z + x803*z + x805*y + x814*M[6] + x86*M[96] + x911*M[7] + x920*M[5] + x922*M[4] + x96*M[93] + y*M[188] + z*M[187] + M[250];
#pragma omp atomic
Ms[251] += x*M[196] + x11*x665 + x11*x666 + x11*x669 + x11*M[151] + x112*x303 + x118*x299 + x121*x268 + x121*M[72] + x125*x283 + x125*M[74] + x147*x230 + x166*x234 + x174*x208 + x177*x221 + x177*M[63] + x182*x214 + x184*x209 + x184*x210 + x184*x211 + x184*M[61] + x208*M[51] + x213*x933 + x214*M[53] + x230*M[45] + x234*M[48] + x299*M[39] + x303*M[38] + x33*x557 + x34*x587 + x34*x589 + x34*x591 + x34*M[138] + x355*x87 + x355*M[30] + x361*x97 + x361*M[33] + x362*x96 + x363*x96 + x365*x96 + x382*x86 + x383*x86 + x384*x86 + x385*x65 + x385*M[25] + x402*M[27] + x41*x570 + x41*x571 + x41*x573 + x41*M[136] + x45*x473 + x45*x474 + x45*x475 + x45*M[115] + x46*x476 + x472*x49 + x472*M[21] + x476*M[20] + x557*M[15] + x623*x923 + x704*x9 + x704*M[9] + x756*M[17] + x775*M[14] + x782*M[13] + x805*z + x806*z + x807*z + x810*y + x814*M[7] + x86*M[97] + x911*M[8] + x913*x925 + x922*M[5] + x924*M[4] + x96*M[94] + y*M[189] + z*M[188] + M[251];
#pragma omp atomic
Ms[252] += x*M[197] + x107*x310 + x108*x310 + x109*x310 + x11*x672 + x11*x673 + x11*x676 + x11*M[152] + x118*x303 + x121*x275 + x121*M[73] + x125*x288 + x125*M[75] + x153*x230 + x177*x227 + x177*M[64] + x178*x208 + x184*x215 + x184*x216 + x184*x218 + x184*M[62] + x185*x214 + x19*x596 + x208*M[52] + x214*M[54] + x230*M[46] + x303*M[39] + x310*M[37] + x34*x593 + x34*x594 + x34*x595 + x34*M[139] + x35*x557 + x355*x90 + x355*M[31] + x368*x96 + x369*x96 + x371*x96 + x385*x70 + x385*M[26] + x402*M[28] + x41*x579 + x41*x580 + x41*x582 + x41*M[137] + x45*x477 + x45*x478 + x45*x479 + x45*M[116] + x476*x49 + x476*M[21] + x5*x818 + x557*M[16] + x596*M[12] + x708*x9 + x708*M[9] + x756*M[18] + x782*M[14] + x810*z + x811*z + x812*z + x814*M[8] + x815*y + x818*M[6] + x904*x934 + x916*x925 + x924*M[5] + x96*M[95] + y*M[190] + z*M[189] + M[252];
#pragma omp atomic
Ms[253] += x*M[198] + x11*x679 + x11*x680 + x11*x683 + x11*M[153] + x112*x310 + x113*x310 + x115*x310 + x121*x282 + x121*M[74] + x160*x230 + x181*x208 + x184*x221 + x184*x222 + x184*x224 + x184*M[63] + x2*x820 + x208*M[53] + x230*M[47] + x24*x596 + x310*M[38] + x355*x93 + x355*M[32] + x376*x96 + x377*x96 + x379*x96 + x38*x557 + x385*x76 + x385*M[27] + x41*x587 + x41*x588 + x41*x590 + x41*M[138] + x45*x480 + x45*x481 + x45*x482 + x45*M[117] + x46*x483 + x47*x483 + x48*x483 + x483*M[20] + x557*M[17] + x596*M[13] + x6*x818 + x815*z + x816*z + x817*z + x818*M[7] + x819*y + x820*M[4] + x926*M[1] + x96*M[96] + y*M[191] + z*M[190] + M[253];
#pragma omp atomic
Ms[254] += x*M[199] + x11*x685 + x11*M[154] + x118*x310 + x121*M[75] + x184*x227 + x184*M[64] + x208*M[54] + x230*M[48] + x310*M[39] + x355*M[33] + x382*x96 + x385*M[28] + x41*x593 + x41*M[139] + x45*x484 + x45*M[118] + x483*x49 + x483*M[21] + x557*M[18] + x596*M[14] + x687*x928 + x715*x9 + x715*M[9] + x818*M[8] + x819*z + x820*M[5] + x926*M[2] + x96*M[97] + z*M[191] + M[254];
#pragma omp atomic
Ms[255] += x*M[200] + x0*x880 + x11*x688 + x11*M[155] + x12*x689 + x122*x292 + x125*M[76] + x170*x231 + x170*M[65] + x234*M[49] + x292*M[40] + x34*x597 + x34*M[140] + x386*x86 + x389*M[29] + x461*x52 + x461*M[22] + x600*M[15] + x689*M[10] + x821*y + x824*M[6] + x86*M[98] + x880*M[3] + x935*M[1] + y*M[192] + M[255];
#pragma omp atomic
Ms[256] += x*M[201] + x11*x690 + x11*x691 + x11*x692 + x11*M[156] + x125*x294 + x125*M[77] + x126*x292 + x128*x292 + x130*x292 + x15*x689 + x16*x689 + x17*x689 + x170*x235 + x170*x237 + x170*x239 + x170*M[66] + x172*x234 + x234*M[50] + x292*M[41] + x34*x601 + x34*x603 + x34*x605 + x34*M[141] + x36*x600 + x389*x88 + x389*M[30] + x390*x86 + x392*x86 + x394*x86 + x461*x55 + x461*x57 + x461*x59 + x461*M[23] + x600*M[16] + x689*M[11] + x7*x824 + x821*z + x822*z + x823*z + x824*M[7] + x825*y + x86*M[99] + x935*M[2] + y*M[193] + z*M[192] + M[256];
#pragma omp atomic
Ms[257] += x*M[202] + x0*x883 + x11*x693 + x11*x694 + x11*x695 + x11*M[157] + x12*x696 + x121*x291 + x121*M[76] + x122*x299 + x125*x297 + x125*M[78] + x132*x292 + x134*x292 + x136*x292 + x170*x241 + x170*x243 + x170*x245 + x170*M[67] + x175*x234 + x177*x231 + x177*M[65] + x234*M[51] + x292*M[42] + x299*M[40] + x34*x607 + x34*x609 + x34*x611 + x34*M[142] + x389*x91 + x389*M[31] + x39*x600 + x396*x86 + x398*x86 + x400*x86 + x402*M[29] + x41*x597 + x41*x598 + x41*x599 + x41*M[140] + x461*x61 + x461*x62 + x461*x63 + x461*M[24] + x468*x52 + x468*M[22] + x600*M[17] + x696*M[10] + x775*M[15] + x796*x930 + x8*x824 + x824*M[8] + x825*z + x826*z + x827*z + x828*y + x86*M[100] + x883*M[3] + x920*M[6] + y*M[194] + z*M[193] + M[257];
#pragma omp atomic
Ms[258] += x*M[203] + x0*x937 + x11*x697 + x11*x698 + x11*x699 + x11*M[158] + x12*x700 + x121*x293 + x121*M[77] + x122*x303 + x125*x301 + x125*M[79] + x126*x299 + x138*x292 + x139*x292 + x140*x292 + x15*x696 + x169*x230 + x170*x247 + x170*x249 + x170*x251 + x170*M[68] + x177*x235 + x177*M[66] + x179*x234 + x230*M[49] + x234*M[52] + x292*M[43] + x299*M[41] + x303*M[40] + x34*x615 + x34*x617 + x34*x619 + x34*M[143] + x386*x96 + x387*x96 + x388*x96 + x389*x94 + x389*M[32] + x402*M[30] + x403*x86 + x405*x86 + x407*x86 + x41*x601 + x41*x602 + x41*x604 + x41*M[141] + x42*x600 + x468*x55 + x468*M[23] + x472*x52 + x472*M[22] + x600*M[18] + x622*x931 + x696*M[11] + x700*M[10] + x775*M[16] + x782*M[15] + x828*z + x829*z + x830*z + x832*y + x86*M[101] + x861*x936 + x919*x930 + x920*M[7] + x922*M[6] + x96*M[98] + y*M[195] + z*M[194] + M[258];
#pragma omp atomic
Ms[259] += x*M[204] + x0*x886 + x11*x701 + x11*x702 + x11*x703 + x11*M[159] + x12*x704 + x121*x296 + x121*M[78] + x125*x305 + x125*M[80] + x126*x303 + x132*x299 + x15*x700 + x170*x253 + x170*x254 + x170*x255 + x170*M[69] + x171*x230 + x177*x241 + x177*M[67] + x182*x234 + x184*x231 + x184*x232 + x184*x233 + x184*M[65] + x230*M[50] + x234*M[53] + x299*M[42] + x303*M[41] + x34*x624 + x34*x626 + x34*x628 + x34*M[144] + x385*x85 + x385*M[29] + x389*x97 + x389*M[33] + x390*x96 + x391*x96 + x393*x96 + x402*M[31] + x409*x86 + x41*x607 + x41*x608 + x41*x610 + x41*M[142] + x411*x86 + x413*x86 + x468*x61 + x468*M[24] + x472*x55 + x472*M[23] + x476*x52 + x476*M[22] + x585*x932 + x613*x933 + x700*M[11] + x704*M[10] + x775*M[17] + x782*M[16] + x832*z + x833*z + x834*z + x835*y + x86*M[102] + x886*M[3] + x920*M[8] + x922*M[7] + x924*M[6] + x96*M[99] + y*M[196] + z*M[195] + M[259];
#pragma omp atomic
Ms[260] += x*M[205] + x0*x938 + x11*x705 + x11*x706 + x11*x707 + x11*M[160] + x12*x708 + x121*x300 + x121*M[79] + x122*x310 + x123*x310 + x124*x310 + x125*x308 + x125*M[81] + x132*x303 + x138*x299 + x15*x704 + x174*x230 + x177*x247 + x177*M[68] + x184*x235 + x184*x236 + x184*x238 + x184*M[66] + x185*x234 + x230*M[51] + x234*M[54] + x299*M[43] + x303*M[42] + x310*M[40] + x33*x596 + x34*x630 + x34*x632 + x34*x634 + x34*M[145] + x385*x87 + x385*M[30] + x396*x96 + x397*x96 + x399*x96 + x402*M[32] + x41*x615 + x41*x616 + x41*x618 + x41*M[143] + x415*x86 + x416*x86 + x417*x86 + x472*x61 + x472*M[24] + x476*x55 + x476*M[23] + x596*M[15] + x623*x933 + x704*M[11] + x708*M[10] + x775*M[18] + x782*M[17] + x835*z + x836*z + x837*z + x839*y + x86*M[103] + x869*x936 + x913*x939 + x922*M[8] + x924*M[7] + x96*M[100] + y*M[197] + z*M[196] + M[260];
#pragma omp atomic
Ms[261] += x*M[206] + x0*x890 + x11*x709 + x11*x710 + x11*x711 + x11*M[161] + x121*x304 + x121*M[80] + x125*x311 + x125*M[82] + x126*x310 + x127*x310 + x129*x310 + x138*x303 + x15*x708 + x177*x253 + x177*M[69] + x178*x230 + x184*x241 + x184*x242 + x184*x244 + x184*M[67] + x230*M[52] + x303*M[43] + x310*M[41] + x34*x636 + x34*x637 + x34*x638 + x34*M[146] + x35*x596 + x385*x90 + x385*M[31] + x402*M[33] + x403*x96 + x404*x96 + x406*x96 + x41*x624 + x41*x625 + x41*x627 + x41*M[144] + x476*x61 + x476*M[24] + x483*x52 + x483*x53 + x483*x54 + x483*M[22] + x5*x820 + x596*M[16] + x708*M[11] + x782*M[18] + x820*M[6] + x839*z + x840*z + x841*z + x842*y + x890*M[3] + x916*x939 + x924*M[8] + x96*M[101] + y*M[198] + z*M[197] + M[261];
#pragma omp atomic
Ms[262] += x*M[207] + x11*x712 + x11*x713 + x11*x714 + x11*M[162] + x12*x715 + x121*x307 + x121*M[81] + x13*x715 + x132*x310 + x133*x310 + x135*x310 + x14*x715 + x181*x230 + x184*x247 + x184*x248 + x184*x250 + x184*M[68] + x230*M[53] + x310*M[42] + x38*x596 + x385*x93 + x385*M[32] + x409*x96 + x41*x630 + x41*x631 + x41*x633 + x41*M[145] + x410*x96 + x412*x96 + x483*x55 + x483*x56 + x483*x58 + x483*M[23] + x596*M[17] + x6*x820 + x715*M[10] + x820*M[7] + x842*z + x843*z + x844*z + x845*y + x940*M[1] + x96*M[102] + y*M[199] + z*M[198] + M[262];
#pragma omp atomic
Ms[263] += x*M[208] + x0*x893 + x11*x716 + x11*M[163] + x121*M[82] + x138*x310 + x15*x715 + x184*x253 + x184*M[69] + x230*M[54] + x310*M[43] + x385*M[33] + x41*x636 + x41*M[146] + x415*x96 + x483*x61 + x483*M[24] + x596*M[18] + x715*M[11] + x820*M[8] + x845*z + x893*M[3] + x940*M[2] + x96*M[103] + z*M[199] + M[263];
#pragma omp atomic
Ms[264] += x*M[209] + x1*x880 + x141*x292 + x170*x256 + x170*M[70] + x18*x689 + x292*M[44] + x34*x639 + x34*M[147] + x418*x86 + x461*x64 + x461*M[25] + x689*M[12] + x846*y + x86*M[104] + x880*M[4] + x941*x942 + y*M[200] + M[264];
#pragma omp atomic
Ms[265] += x*M[210] + x146*x292 + x148*x292 + x150*x292 + x170*x261 + x170*x263 + x170*x265 + x170*M[71] + x23*x689 + x25*x689 + x27*x689 + x292*M[45] + x3*x880 + x34*x644 + x34*x646 + x34*x648 + x34*M[148] + x4*x880 + x423*x86 + x425*x86 + x427*x86 + x461*x69 + x461*x71 + x461*x73 + x461*M[26] + x689*M[13] + x846*z + x847*z + x848*z + x851*y + x86*M[105] + x880*M[5] + y*M[201] + z*M[200] + M[265];
#pragma omp atomic
Ms[266] += x*M[211] + x1*x883 + x141*x299 + x152*x292 + x154*x292 + x156*x292 + x170*x267 + x170*x269 + x170*x271 + x170*M[72] + x177*x256 + x177*M[70] + x18*x696 + x29*x689 + x292*M[46] + x299*M[44] + x30*x689 + x31*x689 + x34*x650 + x34*x652 + x34*x654 + x34*M[149] + x41*x639 + x41*x640 + x41*x641 + x41*M[147] + x429*x86 + x431*x86 + x433*x86 + x461*x75 + x461*x77 + x461*x79 + x461*M[27] + x468*x64 + x468*M[25] + x689*M[14] + x696*M[12] + x851*z + x852*z + x853*z + x854*y + x86*M[106] + x883*M[4] + x943*M[0] + y*M[202] + z*M[201] + M[266];
#pragma omp atomic
Ms[267] += x*M[212] + x1*x937 + x141*x303 + x146*x299 + x159*x292 + x161*x292 + x163*x292 + x170*x274 + x170*x276 + x170*x278 + x170*M[73] + x177*x261 + x177*M[71] + x18*x700 + x23*x696 + x292*M[47] + x299*M[45] + x3*x883 + x303*M[44] + x34*x657 + x34*x659 + x34*x661 + x34*M[150] + x41*x644 + x41*x645 + x41*x647 + x41*M[148] + x418*x96 + x419*x96 + x420*x96 + x436*x86 + x438*x86 + x440*x86 + x461*x81 + x461*x82 + x461*x83 + x461*M[28] + x468*x69 + x468*M[26] + x472*x64 + x472*M[25] + x696*M[13] + x700*M[12] + x794*x944 + x854*z + x855*z + x856*z + x858*y + x86*M[107] + x883*M[5] + x937*M[4] + x96*M[104] + y*M[203] + z*M[202] + M[267];
#pragma omp atomic
Ms[268] += x*M[213] + x1*x886 + x146*x303 + x152*x299 + x165*x292 + x166*x292 + x167*x292 + x170*x281 + x170*x283 + x170*x285 + x170*M[74] + x177*x267 + x177*M[72] + x18*x704 + x184*x256 + x184*x257 + x184*x258 + x184*M[70] + x23*x700 + x29*x696 + x292*M[48] + x299*M[46] + x3*x937 + x303*M[45] + x34*x664 + x34*x666 + x34*x668 + x34*M[151] + x41*x650 + x41*x651 + x41*x653 + x41*M[149] + x423*x96 + x424*x96 + x426*x96 + x443*x86 + x445*x86 + x447*x86 + x468*x75 + x468*M[27] + x472*x69 + x472*M[26] + x476*x64 + x476*M[25] + x562*x945 + x696*M[14] + x700*M[13] + x704*M[12] + x858*z + x859*z + x86*M[108] + x860*z + x862*y + x886*M[4] + x937*M[5] + x96*M[105] + y*M[204] + z*M[203] + M[268];
#pragma omp atomic
Ms[269] += x*M[214] + x1*x938 + x141*x310 + x142*x310 + x143*x310 + x152*x303 + x159*x299 + x170*x287 + x170*x288 + x170*x289 + x170*M[75] + x177*x274 + x177*M[73] + x18*x708 + x184*x261 + x184*x262 + x184*x264 + x184*M[71] + x23*x704 + x29*x700 + x299*M[47] + x3*x886 + x303*M[46] + x310*M[44] + x34*x671 + x34*x673 + x34*x675 + x34*M[152] + x360*x946 + x41*x657 + x41*x658 + x41*x660 + x41*M[150] + x429*x96 + x430*x96 + x432*x96 + x450*x86 + x452*x86 + x454*x86 + x468*x81 + x468*M[28] + x472*x75 + x472*M[27] + x476*x69 + x476*M[26] + x700*M[14] + x704*M[13] + x708*M[12] + x86*M[109] + x862*z + x863*z + x864*z + x866*y + x886*M[5] + x938*M[4] + x96*M[106] + y*M[205] + z*M[204] + M[269];
#pragma omp atomic
Ms[270] += x*M[215] + x1*x890 + x146*x310 + x147*x310 + x149*x310 + x159*x303 + x165*x299 + x177*x281 + x177*M[74] + x184*x267 + x184*x268 + x184*x270 + x184*M[72] + x213*x947 + x23*x708 + x29*x704 + x299*M[48] + x3*x938 + x303*M[47] + x310*M[45] + x34*x678 + x34*x680 + x34*x682 + x34*M[153] + x41*x664 + x41*x665 + x41*x667 + x41*M[151] + x436*x96 + x437*x96 + x439*x96 + x456*x86 + x457*x86 + x458*x86 + x472*x81 + x472*M[28] + x476*x75 + x476*M[27] + x483*x64 + x483*x65 + x483*x66 + x483*M[25] + x704*M[14] + x708*M[13] + x86*M[110] + x866*z + x867*z + x868*z + x870*y + x890*M[4] + x938*M[5] + x96*M[107] + y*M[206] + z*M[205] + M[270];
#pragma omp atomic
Ms[271] += x*M[216] + x111*x948 + x152*x310 + x153*x310 + x155*x310 + x165*x303 + x177*x287 + x177*M[75] + x18*x715 + x184*x274 + x184*x275 + x184*x277 + x184*M[73] + x19*x715 + x20*x715 + x29*x708 + x3*x890 + x303*M[48] + x310*M[46] + x34*x684 + x34*x685 + x34*x686 + x34*M[154] + x41*x671 + x41*x672 + x41*x674 + x41*M[152] + x443*x96 + x444*x96 + x446*x96 + x476*x81 + x476*M[28] + x483*x69 + x483*x70 + x483*x72 + x483*M[26] + x708*M[14] + x715*M[12] + x870*z + x871*z + x872*z + x874*y + x890*M[5] + x96*M[108] + y*M[207] + z*M[206] + M[271];
#pragma omp atomic
Ms[272] += x*M[217] + x1*x893 + x159*x310 + x160*x310 + x162*x310 + x184*x281 + x184*x282 + x184*x284 + x184*M[74] + x2*x893 + x23*x715 + x24*x715 + x26*x715 + x310*M[47] + x41*x678 + x41*x679 + x41*x681 + x41*M[153] + x450*x96 + x451*x96 + x453*x96 + x483*x75 + x483*x76 + x483*x78 + x483*M[27] + x715*M[13] + x874*z + x875*z + x876*z + x877*y + x893*M[4] + x96*M[109] + y*M[208] + z*M[207] + M[272];
#pragma omp atomic
Ms[273] += x*M[218] + x165*x310 + x184*x287 + x184*M[75] + x29*x715 + x3*x893 + x310*M[48] + x41*x684 + x41*M[154] + x456*x96 + x483*x81 + x483*M[28] + x715*M[14] + x877*z + x893*M[5] + x942*x949 + x96*M[110] + z*M[208] + M[273];
#pragma omp atomic
Ms[274] += x170*M[76] + x292*M[49] + x34*M[155] + x461*M[29] + x689*M[15] + x86*M[111] + x880*M[6] + x950*M[1] + y*M[209] + M[274];
#pragma omp atomic
Ms[275] += x170*x294 + x170*M[77] + x172*x292 + x292*M[50] + x34*x691 + x34*M[156] + x36*x689 + x461*x88 + x461*M[30] + x463*x86 + x689*M[16] + x7*x880 + x86*M[112] + x879*z + x880*M[7] + x950*M[2] + y*M[210] + z*M[209] + M[275];
#pragma omp atomic
Ms[276] += x170*x297 + x170*M[78] + x175*x292 + x177*M[76] + x292*M[51] + x299*M[49] + x34*x694 + x34*M[157] + x39*x689 + x41*x688 + x41*M[155] + x461*x91 + x461*M[31] + x466*x86 + x468*M[29] + x689*M[17] + x696*M[15] + x8*x880 + x86*M[113] + x880*M[8] + x881*z + x883*M[6] + x943*M[1] + y*M[211] + z*M[210] + M[276];
#pragma omp atomic
Ms[277] += x170*x301 + x170*M[79] + x177*M[77] + x179*x292 + x292*M[52] + x299*M[50] + x303*M[49] + x34*x698 + x34*M[158] + x41*x690 + x41*M[156] + x42*x689 + (1.0/4320.0)*x421*x621 + x460*x96 + x461*x94 + x461*M[32] + x468*M[30] + x470*x86 + x472*M[29] + x689*M[18] + x696*M[16] + x700*M[15] + x86*M[114] + x882*z + x883*M[7] + x937*M[6] + x943*M[2] + x96*M[111] + y*M[212] + z*M[211] + M[277];
#pragma omp atomic
Ms[278] += x170*x305 + x170*M[80] + x177*M[78] + x182*x292 + x184*x291 + x184*M[76] + x292*M[53] + x299*M[51] + x303*M[50] + x34*x702 + x34*M[159] + x41*x693 + x41*M[157] + x421*x944*M[2] + x461*x97 + x461*M[33] + x462*x96 + x468*M[31] + x472*M[30] + x474*x86 + x476*M[29] + x696*M[17] + x700*M[16] + x704*M[15] + x796*x945 + x86*M[115] + x883*M[8] + x884*z + x886*M[6] + x937*M[7] + x96*M[112] + y*M[213] + z*M[212] + M[278];
#pragma omp atomic
Ms[279] += x169*x310 + x170*x308 + x170*M[81] + x177*M[79] + x184*x293 + x184*M[77] + x185*x292 + x292*M[54] + x299*M[52] + x303*M[51] + x310*M[49] + x34*x706 + x34*M[160] + x41*x697 + x41*M[158] + x465*x96 + x468*M[32] + x472*M[31] + x476*M[30] + x478*x86 + x696*M[18] + x700*M[17] + x704*M[16] + x708*M[15] + x86*M[116] + x885*z + x886*M[7] + x919*x945 + x937*M[8] + x938*M[6] + x951*M[1] + x96*M[113] + y*M[214] + z*M[213] + M[279];
#pragma omp atomic
Ms[280] += x170*x311 + x170*M[82] + x171*x310 + x177*M[80] + x184*x296 + x184*M[78] + x299*M[53] + x303*M[52] + x310*M[50] + x34*x710 + x34*M[161] + x41*x701 + x41*M[159] + x468*M[33] + x469*x96 + x472*M[32] + x476*M[31] + x481*x86 + x483*x85 + x483*M[29] + x613*x947 + x700*M[18] + x704*M[17] + x708*M[16] + x86*M[117] + x886*M[8] + x888*z + x890*M[6] + x938*M[7] + x951*M[2] + x96*M[114] + y*M[215] + z*M[214] + M[280];
#pragma omp atomic
Ms[281] += x174*x310 + x177*M[81] + x184*x300 + x184*M[79] + x299*M[54] + x303*M[53] + x310*M[51] + x33*x715 + x34*x713 + x34*M[162] + x41*x705 + x41*M[160] + x472*M[33] + x473*x96 + x476*M[32] + x483*x87 + x483*M[30] + x484*x86 + x623*x947 + x704*M[18] + x708*M[17] + x715*M[15] + x86*M[118] + x889*z + x890*M[7] + x938*M[8] + x952*M[1] + x96*M[115] + y*M[216] + z*M[215] + M[281];
#pragma omp atomic
Ms[282] += x177*M[82] + x178*x310 + x184*x304 + x184*M[80] + x303*M[54] + x310*M[52] + x34*x716 + x34*M[163] + x35*x715 + x41*x709 + x41*M[161] + x476*M[33] + x477*x96 + x483*x90 + x483*M[31] + x5*x893 + x708*M[18] + x715*M[16] + x890*M[8] + x891*z + x893*M[6] + x952*M[2] + x96*M[116] + y*M[217] + z*M[216] + M[282];
#pragma omp atomic
Ms[283] += x181*x310 + x184*x307 + x184*M[81] + x310*M[53] + x38*x715 + x41*x712 + x41*M[162] + x480*x96 + x483*x93 + x483*M[32] + x6*x893 + x715*M[17] + x892*z + x893*M[7] + x953*M[1] + x96*M[117] + y*M[218] + z*M[217] + M[283];
#pragma omp atomic
Ms[284] += x184*M[82] + x310*M[54] + x41*M[163] + x483*M[33] + x715*M[18] + x893*M[8] + x953*M[2] + x96*M[118] + z*M[218] + M[284];

}

void M2L_10(double x, double y, double z, double * M, double * L) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
double x448;
double x449;
double x450;
double x451;
double x452;
double x453;
double x454;
double x455;
double x456;
double x457;
double x458;
double x459;
double x460;
double x461;
double x462;
double x463;
double x464;
double x465;
double x466;
double x467;
double x468;
double x469;
double x470;
double x471;
double x472;
double x473;
double x474;
double x475;
double x476;
double x477;
double x478;
double x479;
double x480;
double x481;
double x482;
double x483;
double x484;
double x485;
double x486;
double x487;
double x488;
double x489;
double x490;
double x491;
double x492;
double x493;
double x494;
double x495;
double x496;
double x497;
double x498;
double x499;
double x500;
double x501;
double x502;
double x503;
double x504;
double x505;
double x506;
double x507;
double x508;
double x509;
double x510;
double x511;
double x512;
double x513;
double x514;
double x515;
double x516;
double x517;
double x518;
double x519;
double x520;
double x521;
double x522;
double x523;
double x524;
double x525;
double x526;
double x527;
double x528;
double x529;
double x530;
double x531;
double x532;
double x533;
double x534;
double x535;
double x536;
double x537;
double x538;
double x539;
double x540;
double x541;
double x542;
double x543;
double x544;
double x545;
double x546;
double x547;
double x548;
double x549;
double x550;
double x551;
double x552;
double x553;
double x554;
double x555;
double x556;
double x557;
double x558;
double x559;
double x560;
double x561;
double x562;
double x563;
double x564;
double x565;
double x566;
double x567;
double x568;
double x569;
double x570;
double x571;
double x572;
double x573;
double x574;
double x575;
double x576;
double x577;
double x578;
double x579;
double x580;
double x581;
double x582;
double x583;
double x584;
double x585;
double x586;
double x587;
double x588;
double x589;
double x590;
double x591;
double x592;
double x593;
double x594;
double x595;
double x596;
double x597;
double x598;
double x599;
double x600;
double x601;
double x602;
double x603;
double x604;
double x605;
double x606;
double x607;
double x608;
double x609;
double x610;
double x611;
double x612;
double x613;
double x614;
double x615;
double x616;
double x617;
double x618;
double x619;
double x620;
double x621;
double x622;
double x623;
double x624;
double x625;
double x626;
double x627;
double x628;
double x629;
double x630;
double x631;
double x632;
double x633;
double x634;
double x635;
double x636;
double x637;
double x638;
double x639;
double x640;
double x641;
double x642;
double x643;
double x644;
double x645;
double x646;
double x647;
double x648;
double x649;
double x650;
double x651;
double x652;
double x653;
double x654;
double x655;
double x656;
double x657;
double x658;
double x659;
double x660;
double x661;
double x662;
double x663;
double x664;
double x665;
double x666;
double x667;
double x668;
double x669;
double x670;
double x671;
double x672;
double x673;
double x674;
double x675;
double x676;
double x677;
double x678;
double x679;
double x680;
double x681;
double x682;
double x683;
double x684;
double x685;
double x686;
double x687;
double x688;
double x689;
double x690;
double x691;
double x692;
double x693;
double x694;
double x695;
double x696;
double x697;
double x698;
double x699;
double x700;
double x701;
double x702;
double x703;
double x704;
double x705;
double x706;
double x707;
double x708;
double x709;
double x710;
double x711;
double x712;
double x713;
double x714;
double x715;
double x716;
double x717;
double x718;
double x719;
double x720;
double x721;
double x722;
double x723;
double x724;
double x725;
double x726;
double x727;
double x728;
double x729;
double x730;
double x731;
x0 = (1 / (R*R*R));
x1 = 1.0*x0;
x2 = pow(R, -5);
x3 = 3.0*x2;
x4 = x*y;
x5 = x3*x4;
x6 = x*z;
x7 = x3*x6;
x8 = y*z;
x9 = x3*x8;
x10 = pow(R, -7);
x11 = 15.0*x10;
x12 = x4*z;
x13 = x11*x12;
x14 = -x1;
x15 = (x*x);
x16 = x15*x3;
x17 = x14 + x16;
x18 = (y*y);
x19 = x18*x3;
x20 = x14 + x19;
x21 = 9.0*x2;
x22 = x11*x15;
x23 = -x22;
x24 = x*(x21 + x23);
x25 = x23 + x3;
x26 = x25*y;
x27 = x11*x18;
x28 = -x27;
x29 = y*(x21 + x28);
x30 = x25*z;
x31 = z*(x28 + x3);
x32 = 1.0*x;
x33 = x32*(x27 - x3);
x34 = 45.0*x10;
x35 = -x34;
x36 = pow(R, -9);
x37 = x15*x36;
x38 = 105.0*x37;
x39 = x35 + x38;
x40 = x39*x4;
x41 = x39*x6;
x42 = x18*x36;
x43 = 105.0*x42;
x44 = x35 + x43;
x45 = x44*x8;
x46 = -x11;
x47 = x8*(x38 + x46);
x48 = x32*y;
x49 = x44*x48;
x50 = x43 + x46;
x51 = x32*z;
x52 = x50*x51;
x53 = 315.0*x36;
x54 = pow(R, -11);
x55 = 945.0*x54;
x56 = x15*x55;
x57 = x53 - x56;
x58 = x12*x57;
x59 = x18*x55;
x60 = x32*x8;
x61 = x60*(-x53 + x59);
x62 = (x*x*x*x);
x63 = 105.0*x36;
x64 = x62*x63;
x65 = 90.0*x10;
x66 = -x15*x65 + x21 + x64;
x67 = (y*y*y*y);
x68 = x63*x67;
x69 = -x18*x65 + x21 + x68;
x70 = 2.0*x0 - x16 - x19;
x71 = -225.0*x10;
x72 = x55*x62;
x73 = -x72;
x74 = x*(1050.0*x37 + x71 + x73);
x75 = x55*x67;
x76 = -x75;
x77 = y*(1050.0*x42 + x71 + x76);
x78 = x35 + 630.0*x37 + x73;
x79 = x78*y;
x80 = x78*z;
x81 = 630.0*x42;
x82 = x35 + x76 + x81;
x83 = x82*z;
x84 = x32*(x34 + x75 - x81);
x85 = 1575.0*x36;
x86 = pow(R, -13);
x87 = x62*x86;
x88 = 10395.0*x87;
x89 = x15*x54;
x90 = x85 + x88 - 9450.0*x89;
x91 = x4*x90;
x92 = x6*x90;
x93 = 5670.0*x89;
x94 = x53 + x88 - x93;
x95 = x8*x94;
x96 = x67*x86;
x97 = 10395.0*x96;
x98 = x18*x54;
x99 = x85 + x97 - 9450.0*x98;
x100 = x8*x99;
x101 = x48*x99;
x102 = x53 + x97 - 5670.0*x98;
x103 = x102*x51;
x104 = 14175.0*x54;
x105 = -x104;
x106 = pow(R, -15);
x107 = 135135.0*x106;
x108 = x107*x62;
x109 = x15*x86;
x110 = 103950.0*x109;
x111 = x105 - x108 + x110;
x112 = x111*x12;
x113 = x107*x67;
x114 = x18*x86;
x115 = 103950.0*x114;
x116 = x60*(x104 + x113 - x115);
x117 = pow(x, 6);
x118 = 10395.0*x86;
x119 = x117*x118;
x120 = -x104*x62 + x119 + 4725.0*x37 + x71;
x121 = pow(y, 6);
x122 = x118*x121;
x123 = -x104*x67 + x122 + 4725.0*x42 + x71;
x124 = 11025.0*x36;
x125 = 99225.0*x54;
x126 = x107*x117;
x127 = -x126;
x128 = x127 + 218295.0*x87;
x129 = x*(x124 - x125*x15 + x128);
x130 = 42525.0*x54;
x131 = x127 - x130*x15 + x85 + 155925.0*x87;
x132 = x131*y;
x133 = x107*x121;
x134 = -x133;
x135 = x134 + 218295.0*x96;
x136 = y*(x124 - x125*x18 + x135);
x137 = x131*z;
x138 = 155925.0*x96;
x139 = x130*x18;
x140 = x134 + x138 - x139 + x85;
x141 = x140*z;
x142 = x32*(x133 - x138 + x139 - x85);
x143 = pow(R, -17);
x144 = x117*x143;
x145 = 2027025.0*x144;
x146 = x106*x62;
x147 = 2837835.0*x146;
x148 = -x125;
x149 = 1091475.0*x86;
x150 = x149*x15;
x151 = x148 + x150;
x152 = x145 - x147 + x151;
x153 = x152*x4;
x154 = x152*x6;
x155 = x149*x18;
x156 = x148 + x155;
x157 = x121*x143;
x158 = 2027025.0*x157;
x159 = x106*x67;
x160 = 2837835.0*x159;
x161 = x158 - x160;
x162 = x156 + x161;
x163 = x162*x8;
x164 = 2027025.0*x146;
x165 = -x164;
x166 = 467775.0*x86;
x167 = x15*x166;
x168 = x105 + x145 + x165 + x167;
x169 = x168*x8;
x170 = x18*x38;
x171 = x170 + x25 + x28;
x172 = x162*x48;
x173 = 2027025.0*x159;
x174 = -x173;
x175 = x105 + x158 + x166*x18 + x174;
x176 = x175*x51;
x177 = x143*x62;
x178 = 42567525.0*x177;
x179 = pow(R, -19);
x180 = 34459425.0*x179;
x181 = x117*x180;
x182 = x106*x15;
x183 = 14189175.0*x182;
x184 = -x183;
x185 = x149 + x184;
x186 = x178 - x181 + x185;
x187 = x12*x186;
x188 = x18*x53;
x189 = -x18*x56;
x190 = x*(x188 + x189 + x39);
x191 = x15*x53;
x192 = y*(x189 + x191 + x44);
x193 = z*(x189 + x38 + x50);
x194 = x121*x180;
x195 = x143*x67;
x196 = 42567525.0*x195;
x197 = x106*x18;
x198 = 14189175.0*x197;
x199 = x60*(-x149 + x194 - x196 + x198);
x200 = -x24 + x33;
x201 = -x26 - x29;
x202 = -x30 - x31;
x203 = -2835.0*x98;
x204 = x109*x18;
x205 = 10395.0*x204;
x206 = x203 + x205;
x207 = 945.0*x36;
x208 = -2835.0*x89;
x209 = x207 + x208;
x210 = x4*(x206 + x209);
x211 = x6*(x206 + x57);
x212 = x8*(x205 + x208 + x53 - x59);
x213 = 31185.0*x114;
x214 = x15*x18;
x215 = -8505.0*x54;
x216 = 31185.0*x109;
x217 = x215 + x216;
x218 = x12*(-x107*x214 + x213 + x217);
x219 = -x40 - x49;
x220 = -x41 - x52;
x221 = -x45 - x47;
x222 = pow(x, 8);
x223 = 2027025.0*x143;
x224 = x222*x223;
x225 = 3783780.0*x106;
x226 = -x117*x225 + x124 + x224 + 2182950.0*x87 - 396900.0*x89;
x227 = pow(y, 8);
x228 = x223*x227;
x229 = -x121*x225 + x124 + x228 + 2182950.0*x96 - 396900.0*x98;
x230 = -x58 + x61;
x231 = -893025.0*x54;
x232 = x180*x222;
x233 = -x232;
x234 = x*(13097700.0*x109 + 72972900.0*x144 - 51081030.0*x146 + x231 + x233);
x235 = x180*x227;
x236 = -x235;
x237 = y*(13097700.0*x114 + 72972900.0*x157 - 51081030.0*x159 + x231 + x236);
x238 = 4365900.0*x109;
x239 = 56756700.0*x144 - 28378350.0*x146 + x233 + x238;
x240 = x148 + x239;
x241 = x240*y;
x242 = x240*z;
x243 = 4365900.0*x114;
x244 = 56756700.0*x157;
x245 = 28378350.0*x159;
x246 = x236 + x243 + x244 - x245;
x247 = z*(x148 + x246);
x248 = x32*(x125 + x235 - x243 - x244 + x245);
x249 = 9823275.0*x86;
x250 = pow(R, -21);
x251 = 654729075.0*x250;
x252 = x222*x251;
x253 = x117*x179;
x254 = 766215450.0*x177 - 170270100.0*x182 + x249 + x252 - 1240539300.0*x253;
x255 = x254*x4;
x256 = x254*x6;
x257 = 964863900.0*x253;
x258 = 425675250.0*x177;
x259 = 56756700.0*x182;
x260 = x8*(x149 + x252 - x257 + x258 - x259);
x261 = x227*x251;
x262 = x121*x179;
x263 = 766215450.0*x195 - 170270100.0*x197 + x249 + x261 - 1240539300.0*x262;
x264 = x263*x8;
x265 = 105.0*x10;
x266 = -x170 - 12.0*x2;
x267 = x18*x265 + x22 + x266 - x68;
x268 = x15*x265 + x266 + x27 - x64;
x269 = x263*x48;
x270 = 964863900.0*x262;
x271 = 425675250.0*x195;
x272 = x51*(x149 - 56756700.0*x197 + x261 - x270 + x271);
x273 = pow(x, 10)*x251;
x274 = 1550674125.0*x179;
x275 = 49116375.0*x109 + 1277025750.0*x144 - 425675250.0*x146 - x222*x274 + x231 + x273;
x276 = x251*pow(y, 10);
x277 = 49116375.0*x114 + 1277025750.0*x157 - 425675250.0*x159 - x227*x274 + x231 + x276;
x278 = 120.0*x10;
x279 = -x15*x278 - x18*x278 + 210.0*x18*x37 + 24.0*x2 + x64 + x68;
x280 = x109*x67;
x281 = 10395.0*x280;
x282 = x18*x93;
x283 = -x282;
x284 = x191 + x281 + x283 + x82;
x285 = x18*x88;
x286 = x188 + x283 + x285 + x78;
x287 = 62370.0*x204;
x288 = -x113*x15;
x289 = x287 + x288;
x290 = 31185.0*x96;
x291 = x290 - 17010.0*x98;
x292 = x*(x209 + x289 + x291);
x293 = x104*x18;
x294 = -x293;
x295 = x110*x18;
x296 = -x108*x18;
x297 = x*(x294 + x295 + x296 + x90);
x298 = x203 + x287 + x296;
x299 = 31185.0*x87;
x300 = 17010.0*x89;
x301 = x207 + x299 - x300;
x302 = y*(x298 + x301);
x303 = x104*x15;
x304 = -x303;
x305 = y*(x288 + x295 + x304 + x99);
x306 = z*(x298 + x94);
x307 = z*(x102 + x208 + x289);
x308 = 155925.0*x109;
x309 = 311850.0*x114;
x310 = -x130;
x311 = x18*x182;
x312 = -1351350.0*x311;
x313 = x310 + x312;
x314 = -405405.0*x159;
x315 = x15*x195;
x316 = 2027025.0*x315;
x317 = x314 + x316;
x318 = x4*(x308 + x309 + x313 + x317);
x319 = 155925.0*x114;
x320 = 311850.0*x109;
x321 = -405405.0*x146;
x322 = x177*x18;
x323 = 2027025.0*x322;
x324 = x321 + x323;
x325 = x4*(x313 + x319 + x320 + x324);
x326 = x6*(x111 + x312 + x319 + x323);
x327 = 187110.0*x114;
x328 = -810810.0*x311;
x329 = x6*(x217 + x317 + x327 + x328);
x330 = x8*(x105 - x113 + x115 + x308 + x312 + x316);
x331 = 187110.0*x109 + x215;
x332 = x8*(x213 + x324 + x328 + x331);
x333 = -2027025.0*x197;
x334 = x18*x62;
x335 = x143*x214;
x336 = x166 + 20270250.0*x335;
x337 = 6081075.0*x177;
x338 = -4054050.0*x182 + x337;
x339 = x12*(-x180*x334 + x333 + x336 + x338);
x340 = -2027025.0*x182;
x341 = x15*x67;
x342 = 6081075.0*x195 - 4054050.0*x197;
x343 = x12*(-x180*x341 + x336 + x340 + x342);
x344 = 15120.0*x54;
x345 = -x119;
x346 = -x285;
x347 = 270.0*x10 + x282;
x348 = -x188 + x344*x62 + x345 + x346 + x347 - 5355.0*x37;
x349 = -x122;
x350 = -x281;
x351 = -x191 + x344*x67 + x347 + x349 + x350 - 5355.0*x42;
x352 = -x190;
x353 = x352 - x74;
x354 = x352 + x84;
x355 = -x192;
x356 = x355 - x77;
x357 = x355 - x79;
x358 = -x193;
x359 = x358 - x80;
x360 = x358 - x83;
x361 = -x210;
x362 = x361 - x91;
x363 = -x101 + x361;
x364 = -x211;
x365 = -x103 + x364;
x366 = x364 - x92;
x367 = -x212;
x368 = x367 - x95;
x369 = -x100 + x367;
x370 = -x218;
x371 = -x112 + x370;
x372 = x116 + x370;
x373 = x167*x18;
x374 = x145*x18;
x375 = x164*x18;
x376 = x131 + x294 + x373 + x374 - x375;
x377 = x15*x158;
x378 = x15*x173;
x379 = x140 + x304 + x373 + x377 - x378;
x380 = x18*x89;
x381 = -x15*x207 - x18*x207 + x346 + x350 + 11340.0*x380 + x65 + x72 + x75;
x382 = -x18*x181;
x383 = x155 - x18*x183;
x384 = x*(x152 + x178*x18 + x382 + x383);
x385 = -6081075.0*x311;
x386 = x310 + x385;
x387 = 30405375.0*x195;
x388 = -x15*x194;
x389 = x15*x387 + x308 + x388;
x390 = 6081075.0*x157;
x391 = 1403325.0*x114 - 6081075.0*x159 + x390;
x392 = x*(x386 + x389 + x391);
x393 = y*(x15*x196 + x151 + x161 + x383 + x388);
x394 = 30405375.0*x177;
x395 = x18*x394 + x319 + x382;
x396 = 6081075.0*x144;
x397 = 1403325.0*x109 - 6081075.0*x146 + x396;
x398 = y*(x386 + x395 + x397);
x399 = z*(x175 + x385 + x389);
x400 = z*(x168 + x385 + x395);
x401 = 212837625.0*x335;
x402 = x401 + 3274425.0*x86;
x403 = x117*x251;
x404 = x18*x403;
x405 = -103378275.0*x253 + x404;
x406 = -x198;
x407 = 723647925.0*x179;
x408 = -x334*x407 + x406;
x409 = x4*(127702575.0*x177 - 42567525.0*x182 + x402 + x405 + x408);
x410 = -x341*x407;
x411 = x121*x251;
x412 = x15*x411;
x413 = -103378275.0*x262 + x412;
x414 = x4*(x184 + 127702575.0*x195 - 42567525.0*x197 + x402 + x410 + x413);
x415 = 91216125.0*x195;
x416 = 516891375.0*x179;
x417 = x166 + 91216125.0*x335;
x418 = x6*(-18243225.0*x197 + x340 - x341*x416 + x413 + x415 + x417);
x419 = x6*(x186 + x401 + x404 + x408);
x420 = 91216125.0*x177;
x421 = x8*(-18243225.0*x182 + x333 - x334*x416 + x405 + x417 + x420);
x422 = x8*(x185 - x194 + x196 + x401 + x406 + x410 + x412);
x423 = 3918915.0*x106;
x424 = -x228;
x425 = 2338875.0*x86;
x426 = -x377;
x427 = -12600.0*x36 - x373;
x428 = x121*x423 + x303 + x378 + x424 - x425*x67 + x426 + x427 + 439425.0*x98;
x429 = -x224;
x430 = -x374;
x431 = x117*x423 + x293 + x375 - x425*x62 + x427 + x429 + x430 + 439425.0*x89;
x432 = 16065.0*x54;
x433 = -360.0*x10 - x18*x300;
x434 = x122 + 20790.0*x280 + x285 + 1260.0*x37 + 6300.0*x42 - x432*x67 + x433 + x73;
x435 = x119 + 20790.0*x18*x87 + x281 + 6300.0*x37 + 1260.0*x42 - x432*x62 + x433 + x76;
x436 = x15*x159;
x437 = 810810.0*x436;
x438 = -x437;
x439 = x177*x67;
x440 = 2027025.0*x439;
x441 = x146*x18;
x442 = 810810.0*x441;
x443 = x440 - x442;
x444 = 374220.0*x204 + x291 + x301 + x438 + x443;
x445 = x62*x67;
x446 = -x180*x445;
x447 = x310 - 8108100.0*x311 + x446;
x448 = x321 + 12162150.0*x322;
x449 = 935550.0*x114 + x174;
x450 = x*(20270250.0*x315 + x320 + x447 + x448 + x449);
x451 = x314 + 12162150.0*x315;
x452 = 935550.0*x109 + x165;
x453 = y*(x309 + 20270250.0*x322 + x447 + x451 + x452);
x454 = z*(-4864860.0*x311 + x327 + x331 + x446 + x448 + x451);
x455 = x251*x445;
x456 = 344594250.0*x179;
x457 = -x341*x456 + x387 + x455;
x458 = -x334*x456 + x394;
x459 = x4*(-20270250.0*x182 - 20270250.0*x197 + 202702500.0*x335 + x425 + x457 + x458);
x460 = 206756550.0*x179;
x461 = x166 + 121621500.0*x335;
x462 = x6*(-12162150.0*x197 - x334*x460 + x338 + x457 + x461);
x463 = x8*(-12162150.0*x182 - x341*x460 + x342 + x455 + x458 + x461);
x464 = 2*x192 + x77 + x79;
x465 = 2*x193 + x80 + x83;
x466 = 2*x190 + x74 - x84;
x467 = 17010.0*x54;
x468 = 720.0*x10 - x18*x299 - x216*x67 + x345 + x349 - 7560.0*x37 + 34020.0*x380 - 7560.0*x42 + x467*x62 + x467*x67;
x469 = x100 + 2*x212 + x95;
x470 = -x292;
x471 = x142 + x470;
x472 = -x129;
x473 = -x297;
x474 = x472 + x473;
x475 = -x132;
x476 = -x302;
x477 = x475 + x476;
x478 = -x136;
x479 = -x305;
x480 = x478 + x479;
x481 = -x137;
x482 = -x306;
x483 = x481 + x482;
x484 = -x141;
x485 = -x307;
x486 = x484 + x485;
x487 = x101 + 2*x210 + x91;
x488 = x103 + 2*x211 + x92;
x489 = x15*x271;
x490 = x15*x261;
x491 = x18*x259;
x492 = -x491;
x493 = x15*x270;
x494 = x151 + x246 + x489 + x490 + x492 - x493;
x495 = x18*x258;
x496 = x18*x252;
x497 = x18*x257;
x498 = x156 + x239 + x492 + x495 + x496 - x497;
x499 = -x153;
x500 = -x325;
x501 = x499 + x500;
x502 = -x172;
x503 = -x318;
x504 = x502 + x503;
x505 = -x154;
x506 = -x326;
x507 = x505 + x506;
x508 = -x176;
x509 = -x329;
x510 = x508 + x509;
x511 = -x163;
x512 = -x330;
x513 = x511 + x512;
x514 = -x169;
x515 = -x332;
x516 = x514 + x515;
x517 = x112 - x116 + 2*x218;
x518 = -x343;
x519 = x199 + x518;
x520 = -x187;
x521 = -x339;
x522 = x520 + x521;
x523 = 1585133550.0*x179;
x524 = -x273;
x525 = -x496;
x526 = x491 + 992250.0*x54;
x527 = -53482275.0*x109 - 1333782450.0*x144 + 454053600.0*x146 - x155 + x222*x523 - x495 + x497 + x524 + x525 + x526;
x528 = -x276;
x529 = -x490;
x530 = -53482275.0*x114 - x150 - 1333782450.0*x157 + 454053600.0*x159 + x227*x523 - x489 + x493 + x526 + x528 + x529;
x531 = x133 + x426;
x532 = -841995.0*x204 - 2520.0*x36 - x440;
x533 = x15*x160 - x299 + x442 + x531 + x532 + 31185.0*x89 - 187110.0*x96 + 59535.0*x98;
x534 = x126 + x430;
x535 = x147*x18 - x290 + x437 + x532 + x534 - 187110.0*x87 + 59535.0*x89 + 31185.0*x98;
x536 = 4054050.0*x106;
x537 = x144*x18;
x538 = 1309770.0*x204 + 15120.0*x36;
x539 = -x117*x536 + x224 + x290 + x438 + x440 - 4864860.0*x441 + 4054050.0*x537 + x538 + 2525985.0*x87 - 498960.0*x89 - 45360.0*x98;
x540 = x15*x157;
x541 = -x121*x536 + x228 + x299 - 4864860.0*x436 + x443 + x538 + 4054050.0*x540 - 45360.0*x89 + 2525985.0*x96 - 498960.0*x98;
x542 = x470 + x473;
x543 = x476 + x479;
x544 = x482 + x485;
x545 = x18*x420;
x546 = x411*x62;
x547 = x15*x262;
x548 = 206756550.0*x547;
x549 = x546 - x548;
x550 = x416*x445;
x551 = -x550;
x552 = x310 - 36486450.0*x311 + x551;
x553 = 182432250.0*x315 + x391 + x452 + x545 + x549 + x552;
x554 = x18*x253;
x555 = 206756550.0*x554;
x556 = -x555;
x557 = x15*x415;
x558 = x403*x67;
x559 = x557 + x558;
x560 = 182432250.0*x322 + x397 + x449 + x552 + x556 + x559;
x561 = x500 + x503;
x562 = x506 + x509;
x563 = x512 + x515;
x564 = x518 + x521;
x565 = -x234;
x566 = -x384;
x567 = x565 + x566;
x568 = -x392;
x569 = x248 + x568;
x570 = -x237;
x571 = -x393;
x572 = x570 + x571;
x573 = -x241;
x574 = -x398;
x575 = x573 + x574;
x576 = -x242;
x577 = -x400;
x578 = x576 + x577;
x579 = -x247;
x580 = -x399;
x581 = x579 + x580;
x582 = -x255;
x583 = -x409;
x584 = x582 + x583;
x585 = -x269;
x586 = -x414;
x587 = x585 + x586;
x588 = -x272;
x589 = -x418;
x590 = x588 + x589;
x591 = -x256;
x592 = -x419;
x593 = x591 + x592;
x594 = -x260;
x595 = -x421;
x596 = x594 + x595;
x597 = -x264;
x598 = -x422;
x599 = x597 + x598;
x600 = x128 + x135 + 1683990.0*x204 + 5040.0*x36 + x374 + x377 - 3648645.0*x436 + 4054050.0*x439 - 3648645.0*x441 - 90720.0*x89 - 90720.0*x98;
x601 = 4189185.0*x106;
x602 = -2993760.0*x204 - x337*x67 - 20160.0*x36;
x603 = x121*x601 - x15*x390 + x424 + 8513505.0*x436 + 4459455.0*x441 + x534 + x602 - 249480.0*x87 + 136080.0*x89 - 2744280.0*x96 + 589680.0*x98;
x604 = x117*x601 - x18*x396 + x429 + 4459455.0*x436 + 8513505.0*x441 + x531 + x602 - 2744280.0*x87 + 589680.0*x89 - 249480.0*x96 + 136080.0*x98;
x605 = -x558;
x606 = 93243150.0*x311 + 141750.0*x54 + x550;
x607 = x232 + x525;
x608 = -5769225.0*x109 - 2027025.0*x114 - 62837775.0*x144 + 34459425.0*x146 + x173 - 608107500.0*x322 + 1171620450.0*x554 - x557 + x605 + x606 + x607;
x609 = -x546;
x610 = x235 + x529;
x611 = -2027025.0*x109 - 5769225.0*x114 - 62837775.0*x157 + 34459425.0*x159 + x164 - 608107500.0*x315 - x545 + 1171620450.0*x547 + x606 + x609 + x610;
x612 = 4324320.0*x106;
x613 = -x117*x612 - x121*x612 + 5987520.0*x204 + x224 + x228 + 40320.0*x36 - 12972960.0*x436 + 12162150.0*x439 - 12972960.0*x441 + 8108100.0*x537 + 8108100.0*x540 + 2993760.0*x87 - 725760.0*x89 + 2993760.0*x96 - 725760.0*x98;
x614 = 1619592975.0*x179;
x615 = 1309458150.0*x250;
x616 = x15*x227;
x617 = -149999850.0*x311 - 1134000.0*x54 + x551;
x618 = 3118500.0*x109 + 59251500.0*x114 + 1396620225.0*x157 - 488513025.0*x159 + x165 - x227*x614 + x276 + 1033782750.0*x315 + x545 + x546 - 2136484350.0*x547 + x615*x616 + x617;
x619 = x18*x222;
x620 = 59251500.0*x109 + 3118500.0*x114 + 1396620225.0*x144 - 488513025.0*x146 + x174 - x222*x614 + x273 + 1033782750.0*x322 - 2136484350.0*x554 + x559 + x615*x619 + x617;
x621 = x129 + x292 + 2*x297;
x622 = x132 + 2*x302 + x305;
x623 = x136 + x302 + 2*x305;
x624 = x137 + 2*x306 + x307;
x625 = x141 + x306 + 2*x307;
x626 = x179*x445;
x627 = -x396 + x555 + x605;
x628 = -x390 + x548 + x609;
x629 = 8108100.0*x146 - x15*x425 + 8108100.0*x159 - x18*x425 + 72972900.0*x311 - 273648375.0*x315 - 273648375.0*x322 + 85050.0*x54 + 1033782750.0*x626 + x627 + x628;
x630 = -x142 + 2*x292 + x297;
x631 = x153 + x318 + 2*x325;
x632 = x154 + 2*x326 + x329;
x633 = x163 + 2*x330 + x332;
x634 = x169 + x330 + 2*x332;
x635 = x172 + 2*x318 + x325;
x636 = x176 + x326 + 2*x329;
x637 = x187 + 2*x339 + x343;
x638 = -x199 + x339 + 2*x343;
x639 = -x450;
x640 = x566 + x639;
x641 = x568 + x639;
x642 = -x453;
x643 = x571 + x642;
x644 = x574 + x642;
x645 = -x454;
x646 = x580 + x645;
x647 = x577 + x645;
x648 = x121*x62;
x649 = x236 + x490;
x650 = -x274*x445 - 166216050.0*x311 - 226800.0*x54;
x651 = 8108100.0*x114 - 10135125.0*x146 + 68918850.0*x157 - 42567525.0*x159 + x238 + 881755875.0*x315 + 364864500.0*x322 + x396 - 1378377000.0*x547 + x556 + x558 + x615*x648 + x649 + x650;
x652 = x117*x67;
x653 = x233 + x496;
x654 = 8108100.0*x109 + 68918850.0*x144 - 42567525.0*x146 - 10135125.0*x159 + x243 + 364864500.0*x315 + 881755875.0*x322 + x390 + x549 - 1378377000.0*x554 + x615*x652 + x650 + x653;
x655 = -x459;
x656 = x583 + x655;
x657 = x586 + x655;
x658 = -x462;
x659 = x589 + x658;
x660 = x592 + x658;
x661 = -x463;
x662 = x595 + x661;
x663 = x598 + x661;
x664 = 1654052400.0*x179;
x665 = 1964187225.0*x250;
x666 = -x652*x665;
x667 = 316215900.0*x311 + 1360800.0*x54 + 2067565500.0*x626;
x668 = -67359600.0*x109 - 7484400.0*x114 - 1465539075.0*x144 + 531080550.0*x146 + 12162150.0*x159 + x222*x664 - 456080625.0*x315 - 1915538625.0*x322 + x524 + 3514861350.0*x554 - x619*x665 + x628 + x666 + x667;
x669 = -x648*x665;
x670 = -7484400.0*x109 - 67359600.0*x114 + 12162150.0*x146 - 1465539075.0*x157 + 531080550.0*x159 + x227*x664 - 1915538625.0*x315 - 456080625.0*x322 + x528 + 3514861350.0*x547 - x616*x665 + x627 + x667 + x669;
x671 = x142 - 3*x292 - 3*x297 + x472;
x672 = -3*x302 - 3*x305 + x475 + x478;
x673 = -3*x306 - 3*x307 + x481 + x484;
x674 = -12474000.0*x109 - 12474000.0*x114 - 74999925.0*x144 + 52702650.0*x146 - 74999925.0*x157 + 52702650.0*x159 + 332432100.0*x311 - 1246620375.0*x315 - 1246620375.0*x322 + 453600.0*x54 + 1585133550.0*x547 + 1585133550.0*x554 + x607 + x610 + 3101348250.0*x626 + x666 + x669;
x675 = -3*x318 - 3*x325 + x499 + x502;
x676 = -3*x326 - 3*x329 + x505 + x508;
x677 = -3*x330 - 3*x332 + x511 + x514;
x678 = x199 - 3*x339 - 3*x343 + x520;
x679 = 1688511825.0*x179;
x680 = 2618916300.0*x250;
x681 = 3928374450.0*x250;
x682 = -648648000.0*x311 - 1814400.0*x54 - 5168913750.0*x626;
x683 = 19958400.0*x109 + 79833600.0*x114 + 81081000.0*x144 - 64864800.0*x146 + 1540539000.0*x157 - 583783200.0*x159 - x227*x679 + x276 + 3162159000.0*x315 + 1702701000.0*x322 - 5099994900.0*x547 - 1791890100.0*x554 + x616*x680 + x648*x681 + x652*x680 + x653 + x682;
x684 = 79833600.0*x109 + 19958400.0*x114 + 1540539000.0*x144 - 583783200.0*x146 + 81081000.0*x157 - 64864800.0*x159 - x222*x679 + x273 + 1702701000.0*x315 + 3162159000.0*x322 - 1791890100.0*x547 - 5099994900.0*x554 + x619*x680 + x648*x680 + x649 + x652*x681 + x682;
x685 = 1722971250.0*x179;
x686 = 3273645375.0*x250;
x687 = 6547290750.0*x250;
x688 = -99792000.0*x109 - 99792000.0*x114 - 1621620000.0*x144 + 648648000.0*x146 - 1621620000.0*x157 + 648648000.0*x159 + x222*x685 + x227*x685 + 1297296000.0*x311 - 4864860000.0*x315 - 4864860000.0*x322 + x524 + x528 + 3628800.0*x54 + 6891885000.0*x547 + 6891885000.0*x554 - x616*x686 - x619*x686 + 10337827500.0*x626 - x648*x687 - x652*x687;
x689 = x234 + 2*x384 + x450;
x690 = x237 + 2*x393 + x453;
x691 = x241 + 2*x398 + x453;
x692 = x242 + 2*x400 + x454;
x693 = x247 + 2*x399 + x454;
x694 = -x248;
x695 = 2*x392 + x450 + x694;
x696 = x255 + 2*x409 + x459;
x697 = x256 + 2*x419 + x462;
x698 = x260 + 2*x421 + x463;
x699 = x264 + 2*x422 + x463;
x700 = x269 + 2*x414 + x459;
x701 = x272 + 2*x418 + x462;
x702 = x384 + x392 + 2*x450;
x703 = x393 + x398 + 2*x453;
x704 = x399 + x400 + 2*x454;
x705 = x409 + x414 + 2*x459;
x706 = x418 + x419 + 2*x462;
x707 = x421 + x422 + 2*x463;
x708 = -3*x450;
x709 = -3*x384 + x565 + x568 + x708;
x710 = x248 - 3*x392 + x566 + x708;
x711 = -3*x453;
x712 = -3*x393 + x570 + x574 + x711;
x713 = -3*x398 + x571 + x573 + x711;
x714 = -3*x454;
x715 = -3*x400 + x576 + x580 + x714;
x716 = -3*x399 + x577 + x579 + x714;
x717 = -3*x459;
x718 = -3*x409 + x582 + x586 + x717;
x719 = -3*x414 + x583 + x585 + x717;
x720 = -3*x462;
x721 = -3*x418 + x588 + x592 + x720;
x722 = -3*x419 + x589 + x591 + x720;
x723 = -3*x463;
x724 = -3*x421 + x594 + x598 + x723;
x725 = -3*x422 + x595 + x597 + x723;
x726 = x237 + x241 + 4*x393 + 4*x398 + 6*x453;
x727 = x242 + x247 + 4*x399 + 4*x400 + 6*x454;
x728 = x234 + 4*x384 + 4*x392 + 6*x450 + x694;
x729 = x260 + x264 + 4*x421 + 4*x422 + 6*x463;
x730 = x255 + x269 + 4*x409 + 4*x414 + 6*x459;
x731 = x256 + x272 + 4*x418 + 4*x419 + 6*x462;
#pragma omp atomic
L[0] += -x*x1*M[0] - x1*y*M[1] - x1*z*M[2] + x100*M[77] + x101*M[70] + x103*M[71] + x112*M[87] - x116*M[105] + x120*M[55] + x123*M[76] + x129*M[83] - x13*M[13] + x132*M[84] + x136*M[111] + x137*M[85] + x141*M[112] - x142*M[104] + x153*M[120] + x154*M[121] + x163*M[156] + x169*M[123] + x17*M[3] + x171*M[22] + x172*M[147] + x176*M[148] + x187*M[168] + x190*M[37] + x192*M[40] + x193*M[41] - x199*M[201] + x20*M[6] + x200*M[14] + x201*M[17] + x202*M[18] + x210*M[61] + x211*M[62] + x212*M[66] + x218*M[94] + x219*M[27] + x220*M[28] + x221*M[32] + x226*M[119] + x229*M[155] + x230*M[47] + x234*M[164] + x237*M[209] + x24*M[9] + x241*M[165] + x242*M[166] + x247*M[210] - x248*M[200] + x255*M[220] + x256*M[221] + x26*M[10] + x260*M[223] + x264*M[275] + x267*M[31] + x268*M[24] + x269*M[264] + x272*M[265] + x275*M[219] + x277*M[274] + x279*M[33] + x284*M[65] + x286*M[58] + x29*M[15] + x292*M[93] + x297*M[86] + x30*M[11] + x302*M[89] + x305*M[98] + x306*M[90] + x307*M[99] + x31*M[16] + x318*M[134] + x325*M[125] + x326*M[126] + x329*M[135] - x33*M[12] + x330*M[141] + x332*M[130] + x339*M[175] + x343*M[186] + x348*M[60] + x351*M[78] + x353*M[39] + x354*M[46] + x356*M[51] + x357*M[42] + x359*M[43] + x360*M[52] + x362*M[63] + x363*M[72] + x365*M[73] + x366*M[64] + x368*M[68] + x369*M[79] + x371*M[96] + x372*M[107] + x376*M[122] + x379*M[140] + x381*M[67] + x384*M[167] + x392*M[185] + x393*M[192] + x398*M[170] + x399*M[193] + x40*M[20] + x400*M[171] + x409*M[225] + x41*M[21] + x414*M[247] + x418*M[248] + x419*M[226] + x421*M[230] + x422*M[256] + x428*M[157] + x431*M[124] + x434*M[80] + x435*M[69] + x444*M[129] + x45*M[30] + x450*M[174] + x453*M[179] + x454*M[180] + x459*M[234] + x462*M[235] + x463*M[241] + x464*M[53] + x465*M[54] + x466*M[48] + x468*M[82] + x469*M[81] + x47*M[23] + x471*M[106] + x474*M[88] + x477*M[91] + x480*M[113] + x483*M[92] + x486*M[114] + x487*M[74] + x488*M[75] + x49*M[25] + x494*M[255] + x498*M[222] + x5*M[4] + x501*M[127] + x504*M[149] + x507*M[128] + x510*M[150] + x513*M[158] + x516*M[132] + x517*M[109] + x519*M[203] + x52*M[26] + x522*M[177] + x527*M[224] + x530*M[276] + x533*M[142] + x535*M[131] + x539*M[133] + x541*M[159] + x542*M[95] + x543*M[100] + x544*M[101] + x553*M[240] + x560*M[229] + x561*M[136] + x562*M[137] + x563*M[143] + x564*M[188] + x567*M[169] + x569*M[202] + x572*M[211] + x575*M[172] + x578*M[173] + x58*M[38] + x581*M[212] + x584*M[227] + x587*M[266] + x590*M[267] + x593*M[228] + x596*M[232] + x599*M[277] + x600*M[144] + x603*M[161] + x604*M[146] + x608*M[231] - x61*M[45] + x611*M[257] + x613*M[163] + x618*M[278] + x620*M[233] + x621*M[97] + x622*M[102] + x623*M[115] + x624*M[103] + x625*M[116] + x629*M[242] + x630*M[108] + x631*M[138] + x632*M[139] + x633*M[160] + x634*M[145] + x635*M[151] + x636*M[152] + x637*M[190] + x638*M[205] + x640*M[176] + x641*M[187] + x643*M[194] + x644*M[181] + x646*M[195] + x647*M[182] + x651*M[259] + x654*M[244] + x656*M[236] + x657*M[249] + x659*M[250] + x66*M[19] + x660*M[237] + x662*M[243] + x663*M[258] + x668*M[246] + x670*M[280] + x671*M[110] + x672*M[117] + x673*M[118] + x674*M[261] + x675*M[153] + x676*M[154] + x677*M[162] + x678*M[207] + x683*M[282] + x684*M[263] + x688*M[284] + x689*M[178] + x69*M[29] + x690*M[213] + x691*M[183] + x692*M[184] + x693*M[214] + x695*M[204] + x696*M[238] + x697*M[239] + x698*M[245] + x699*M[279] + x7*M[5] + x70*M[8] + x700*M[268] + x701*M[269] + x702*M[189] + x703*M[196] + x704*M[197] + x705*M[251] + x706*M[252] + x707*M[260] + x709*M[191] + x710*M[206] + x712*M[215] + x713*M[198] + x715*M[199] + x716*M[216] + x718*M[253] + x719*M[270] + x721*M[271] + x722*M[254] + x724*M[262] + x725*M[281] + x726*M[217] + x727*M[218] + x728*M[208] + x729*M[283] + x730*M[272] + x731*M[273] + x74*M[34] + x77*M[49] + x79*M[35] + x80*M[36] + x83*M[50] - x84*M[44] + x9*M[7] + x91*M[56] + x92*M[57] + x95*M[59];
#pragma omp atomic
L[1] += x101*M[49] + x103*M[50] + x112*M[59] - x116*M[77] + x120*M[34] + x129*M[55] - x13*M[7] + x132*M[56] + x137*M[57] - x142*M[76] + x153*M[84] + x154*M[85] + x169*M[87] + x17*M[0] + x171*M[12] + x172*M[111] + x176*M[112] + x187*M[123] + x190*M[22] + x192*M[25] + x193*M[26] - x199*M[156] + x200*M[8] + x210*M[40] + x211*M[41] + x212*M[45] + x218*M[66] + x219*M[17] + x220*M[18] + x226*M[83] + x230*M[32] + x234*M[119] + x24*M[3] + x241*M[120] + x242*M[121] - x248*M[155] + x255*M[165] + x256*M[166] + x26*M[4] + x260*M[168] + x268*M[14] + x269*M[209] + x272*M[210] + x275*M[164] + x284*M[44] + x286*M[37] + x292*M[65] + x297*M[58] + x30*M[5] + x302*M[61] + x305*M[70] + x306*M[62] + x307*M[71] + x318*M[98] + x325*M[89] + x326*M[90] + x329*M[99] - x33*M[6] + x330*M[105] + x332*M[94] + x339*M[130] + x343*M[141] + x348*M[39] + x353*M[24] + x354*M[31] + x357*M[27] + x359*M[28] + x362*M[42] + x363*M[51] + x365*M[52] + x366*M[43] + x368*M[47] + x371*M[68] + x372*M[79] + x376*M[86] + x379*M[104] + x381*M[46] + x384*M[122] + x392*M[140] + x393*M[147] + x398*M[125] + x399*M[148] + x40*M[10] + x400*M[126] + x409*M[170] + x41*M[11] + x414*M[192] + x418*M[193] + x419*M[171] + x421*M[175] + x422*M[201] + x431*M[88] + x435*M[48] + x444*M[93] + x450*M[129] + x453*M[134] + x454*M[135] + x459*M[179] + x462*M[180] + x463*M[186] + x466*M[33] + x47*M[13] + x471*M[78] + x474*M[60] + x477*M[63] + x483*M[64] + x487*M[53] + x488*M[54] + x49*M[15] + x494*M[200] + x498*M[167] + x5*M[1] + x501*M[91] + x504*M[113] + x507*M[92] + x510*M[114] + x516*M[96] + x517*M[81] + x519*M[158] + x52*M[16] + x522*M[132] + x527*M[169] + x533*M[106] + x535*M[95] + x539*M[97] + x542*M[67] + x543*M[72] + x544*M[73] + x553*M[185] + x560*M[174] + x561*M[100] + x562*M[101] + x563*M[107] + x564*M[143] + x567*M[124] + x569*M[157] + x575*M[127] + x578*M[128] + x58*M[23] + x584*M[172] + x587*M[211] + x590*M[212] + x593*M[173] + x596*M[177] + x600*M[108] + x604*M[110] + x608*M[176] - x61*M[30] + x611*M[202] + x620*M[178] + x621*M[69] + x622*M[74] + x624*M[75] + x629*M[187] + x630*M[80] + x631*M[102] + x632*M[103] + x634*M[109] + x635*M[115] + x636*M[116] + x637*M[145] + x638*M[160] + x640*M[131] + x641*M[142] + x643*M[149] + x644*M[136] + x646*M[150] + x647*M[137] + x651*M[204] + x654*M[189] + x656*M[181] + x657*M[194] + x659*M[195] + x66*M[9] + x660*M[182] + x662*M[188] + x663*M[203] + x668*M[191] + x671*M[82] + x674*M[206] + x675*M[117] + x676*M[118] + x678*M[162] + x684*M[208] + x689*M[133] + x691*M[138] + x692*M[139] + x695*M[159] + x696*M[183] + x697*M[184] + x698*M[190] + x7*M[2] + x700*M[213] + x701*M[214] + x702*M[144] + x703*M[151] + x704*M[152] + x705*M[196] + x706*M[197] + x707*M[205] + x709*M[146] + x710*M[161] + x713*M[153] + x715*M[154] + x718*M[198] + x719*M[215] + x721*M[216] + x722*M[199] + x724*M[207] + x728*M[163] + x730*M[217] + x731*M[218] + x74*M[19] + x79*M[20] + x80*M[21] - x84*M[29] + x91*M[35] + x92*M[36] + x95*M[38];
#pragma omp atomic
L[2] += x100*M[50] + x101*M[44] + x103*M[45] + x112*M[57] - x116*M[71] + x123*M[49] - x13*M[5] + x132*M[55] + x136*M[76] + x141*M[77] - x142*M[70] + x153*M[83] + x163*M[112] + x169*M[85] + x171*M[10] + x172*M[104] + x176*M[105] + x187*M[121] + x190*M[20] + x192*M[22] + x193*M[23] - x199*M[148] + x20*M[1] + x201*M[8] + x210*M[37] + x211*M[38] + x212*M[41] + x218*M[62] + x219*M[14] + x221*M[18] + x229*M[111] + x230*M[28] + x237*M[155] + x241*M[119] + x247*M[156] - x248*M[147] + x255*M[164] + x26*M[3] + x260*M[166] + x264*M[210] + x267*M[17] + x269*M[200] + x272*M[201] + x277*M[209] + x284*M[40] + x286*M[35] + x29*M[6] + x292*M[61] + x297*M[56] + x302*M[58] + x305*M[65] + x306*M[59] + x307*M[66] + x31*M[7] + x318*M[93] + x325*M[86] + x326*M[87] + x329*M[94] - x33*M[4] + x330*M[99] + x332*M[90] + x339*M[126] + x343*M[135] + x351*M[51] + x354*M[27] + x356*M[31] + x357*M[24] + x360*M[32] + x362*M[39] + x363*M[46] + x365*M[47] + x368*M[43] + x369*M[52] + x371*M[64] + x372*M[73] + x376*M[84] + x379*M[98] + x381*M[42] + x384*M[120] + x392*M[134] + x393*M[140] + x398*M[122] + x399*M[141] + x40*M[9] + x400*M[123] + x409*M[167] + x414*M[185] + x418*M[186] + x419*M[168] + x421*M[171] + x422*M[193] + x428*M[113] + x434*M[53] + x444*M[89] + x45*M[16] + x450*M[125] + x453*M[129] + x454*M[130] + x459*M[174] + x462*M[175] + x463*M[180] + x464*M[33] + x469*M[54] + x47*M[11] + x471*M[72] + x477*M[60] + x480*M[78] + x486*M[79] + x487*M[48] + x49*M[12] + x494*M[192] + x498*M[165] + x5*M[0] + x501*M[88] + x504*M[106] + x510*M[107] + x513*M[114] + x516*M[92] + x517*M[75] + x519*M[150] + x52*M[13] + x522*M[128] + x530*M[211] + x533*M[100] + x535*M[91] + x541*M[115] + x542*M[63] + x543*M[67] + x544*M[68] + x553*M[179] + x560*M[170] + x561*M[95] + x562*M[96] + x563*M[101] + x564*M[137] + x569*M[149] + x572*M[157] + x575*M[124] + x58*M[21] + x581*M[158] + x584*M[169] + x587*M[202] + x590*M[203] + x596*M[173] + x599*M[212] + x600*M[102] + x603*M[117] + x608*M[172] - x61*M[26] + x611*M[194] + x618*M[213] + x622*M[69] + x623*M[80] + x625*M[81] + x629*M[181] + x630*M[74] + x631*M[97] + x633*M[116] + x634*M[103] + x635*M[108] + x636*M[109] + x637*M[139] + x638*M[152] + x640*M[127] + x641*M[136] + x643*M[142] + x644*M[131] + x646*M[143] + x647*M[132] + x651*M[196] + x654*M[183] + x656*M[176] + x657*M[187] + x659*M[188] + x660*M[177] + x662*M[182] + x663*M[195] + x670*M[215] + x672*M[82] + x674*M[198] + x675*M[110] + x677*M[118] + x678*M[154] + x683*M[217] + x69*M[15] + x690*M[159] + x691*M[133] + x693*M[160] + x695*M[151] + x696*M[178] + x698*M[184] + x699*M[214] + x700*M[204] + x701*M[205] + x702*M[138] + x703*M[144] + x704*M[145] + x705*M[189] + x706*M[190] + x707*M[197] + x710*M[153] + x712*M[161] + x713*M[146] + x716*M[162] + x718*M[191] + x719*M[206] + x721*M[207] + x724*M[199] + x725*M[216] + x726*M[163] + x729*M[218] + x730*M[208] + x77*M[29] + x79*M[19] + x83*M[30] - x84*M[25] + x9*M[2] + x91*M[34] + x95*M[36];
#pragma omp atomic
L[3] += x100*M[49] + x103*M[44] + x112*M[56] - x116*M[70] - x13*M[4] + x137*M[55] + x141*M[76] + x154*M[83] + x163*M[111] + x169*M[84] + x176*M[104] + x187*M[120] + x193*M[22] - x199*M[147] + x200*M[5] + x201*M[7] + x202*M[8] + x211*M[37] + x212*M[40] + x218*M[61] + x219*M[13] + x220*M[14] + x221*M[17] + x230*M[27] + x242*M[119] + x247*M[155] + x256*M[164] + x260*M[165] + x264*M[209] + x267*M[16] + x268*M[11] + x272*M[200] + x279*M[18] + x30*M[3] + x306*M[58] + x307*M[65] + x31*M[6] + x326*M[86] + x329*M[93] + x330*M[98] + x332*M[89] + x339*M[125] + x343*M[134] + x348*M[36] + x351*M[50] + x353*M[21] + x354*M[26] + x356*M[30] + x357*M[23] + x359*M[24] + x360*M[31] + x362*M[38] + x363*M[45] + x365*M[46] + x366*M[39] + x368*M[42] + x369*M[51] + x371*M[63] + x372*M[72] + x381*M[41] + x399*M[140] + x400*M[122] + x41*M[9] + x418*M[185] + x419*M[167] + x421*M[170] + x422*M[192] + x428*M[112] + x431*M[85] + x434*M[52] + x435*M[43] + x45*M[15] + x454*M[129] + x462*M[174] + x463*M[179] + x464*M[32] + x465*M[33] + x466*M[28] + x468*M[54] + x469*M[53] + x47*M[10] + x471*M[71] + x474*M[57] + x477*M[59] + x480*M[77] + x483*M[60] + x486*M[78] + x487*M[47] + x488*M[48] + x501*M[87] + x504*M[105] + x507*M[88] + x510*M[106] + x513*M[113] + x516*M[91] + x517*M[74] + x519*M[149] + x52*M[12] + x522*M[127] + x527*M[166] + x530*M[210] + x533*M[99] + x535*M[90] + x539*M[92] + x541*M[114] + x542*M[62] + x543*M[66] + x544*M[67] + x561*M[94] + x562*M[95] + x563*M[100] + x564*M[136] + x567*M[121] + x569*M[148] + x572*M[156] + x575*M[123] + x578*M[124] + x58*M[20] + x581*M[157] + x584*M[168] + x587*M[201] + x590*M[202] + x593*M[169] + x596*M[172] + x599*M[211] + x600*M[101] + x603*M[116] + x604*M[103] + x608*M[171] - x61*M[25] + x611*M[193] + x613*M[118] + x618*M[212] + x620*M[173] + x621*M[64] + x622*M[68] + x623*M[79] + x624*M[69] + x625*M[80] + x629*M[180] + x630*M[73] + x631*M[96] + x632*M[97] + x633*M[115] + x634*M[102] + x635*M[107] + x636*M[108] + x637*M[138] + x638*M[151] + x640*M[126] + x641*M[135] + x643*M[141] + x644*M[130] + x646*M[142] + x647*M[131] + x651*M[195] + x654*M[182] + x656*M[175] + x657*M[186] + x659*M[187] + x660*M[176] + x662*M[181] + x663*M[194] + x668*M[184] + x670*M[214] + x671*M[75] + x672*M[81] + x673*M[82] + x674*M[197] + x675*M[109] + x676*M[110] + x677*M[117] + x678*M[153] + x683*M[216] + x684*M[199] + x688*M[218] + x689*M[128] + x690*M[158] + x691*M[132] + x692*M[133] + x693*M[159] + x695*M[150] + x696*M[177] + x697*M[178] + x698*M[183] + x699*M[213] + x7*M[0] + x70*M[2] + x700*M[203] + x701*M[204] + x702*M[137] + x703*M[143] + x704*M[144] + x705*M[188] + x706*M[189] + x707*M[196] + x709*M[139] + x710*M[152] + x712*M[160] + x713*M[145] + x715*M[146] + x716*M[161] + x718*M[190] + x719*M[205] + x721*M[206] + x722*M[191] + x724*M[198] + x725*M[215] + x726*M[162] + x727*M[163] + x728*M[154] + x729*M[217] + x730*M[207] + x731*M[208] + x80*M[19] + x83*M[29] + x9*M[1] + x92*M[34] + x95*M[35];
#pragma omp atomic
L[4] += x112*M[38] + x120*M[19] + x129*M[34] + x132*M[35] + x137*M[36] + x153*M[56] + x154*M[57] + x169*M[59] + x171*M[6] + x187*M[87] + x190*M[12] + x192*M[15] + x193*M[16] + x210*M[25] + x211*M[26] + x212*M[30] + x218*M[45] + x226*M[55] + x234*M[83] + x24*M[0] + x241*M[84] + x242*M[85] + x255*M[120] + x256*M[121] + x26*M[1] + x260*M[123] + x268*M[8] + x275*M[119] + x284*M[29] + x286*M[22] + x292*M[44] + x297*M[37] + x30*M[2] + x302*M[40] + x305*M[49] + x306*M[41] + x307*M[50] + x318*M[70] + x325*M[61] + x326*M[62] + x329*M[71] + x330*M[77] + x332*M[66] + x339*M[94] + x343*M[105] + x348*M[24] + x353*M[14] + x357*M[17] + x359*M[18] + x362*M[27] + x366*M[28] + x368*M[32] + x371*M[47] + x376*M[58] + x379*M[76] + x381*M[31] + x384*M[86] + x392*M[104] + x393*M[111] + x398*M[89] + x399*M[112] + x40*M[4] + x400*M[90] + x409*M[125] + x41*M[5] + x414*M[147] + x418*M[148] + x419*M[126] + x421*M[130] + x422*M[156] + x431*M[60] + x435*M[33] + x444*M[65] + x450*M[93] + x453*M[98] + x454*M[99] + x459*M[134] + x462*M[135] + x463*M[141] + x47*M[7] + x474*M[39] + x477*M[42] + x483*M[43] + x494*M[155] + x498*M[122] + x501*M[63] + x507*M[64] + x516*M[68] + x522*M[96] + x527*M[124] + x533*M[78] + x535*M[67] + x539*M[69] + x542*M[46] + x543*M[51] + x544*M[52] + x553*M[140] + x560*M[129] + x561*M[72] + x562*M[73] + x563*M[79] + x564*M[107] + x567*M[88] + x575*M[91] + x578*M[92] + x58*M[13] + x584*M[127] + x593*M[128] + x596*M[132] + x600*M[80] + x604*M[82] + x608*M[131] + x611*M[157] + x620*M[133] + x621*M[48] + x622*M[53] + x624*M[54] + x629*M[142] + x631*M[74] + x632*M[75] + x634*M[81] + x637*M[109] + x640*M[95] + x641*M[106] + x643*M[113] + x644*M[100] + x646*M[114] + x647*M[101] + x651*M[159] + x654*M[144] + x656*M[136] + x657*M[149] + x659*M[150] + x66*M[3] + x660*M[137] + x662*M[143] + x663*M[158] + x668*M[146] + x674*M[161] + x684*M[163] + x689*M[97] + x691*M[102] + x692*M[103] + x696*M[138] + x697*M[139] + x698*M[145] + x702*M[108] + x703*M[115] + x704*M[116] + x705*M[151] + x706*M[152] + x707*M[160] + x709*M[110] + x713*M[117] + x715*M[118] + x718*M[153] + x722*M[154] + x724*M[162] + x74*M[9] + x79*M[10] + x80*M[11] + x91*M[20] + x92*M[21] + x95*M[23];
#pragma omp atomic
L[5] += x101*M[29] + x103*M[30] + x112*M[36] - x116*M[50] - x13*M[2] + x132*M[34] - x142*M[49] + x153*M[55] + x169*M[57] + x171*M[4] + x172*M[76] + x176*M[77] + x187*M[85] + x190*M[10] + x192*M[12] + x193*M[13] - x199*M[112] + x210*M[22] + x211*M[23] + x212*M[26] + x218*M[41] + x219*M[8] + x230*M[18] + x241*M[83] - x248*M[111] + x255*M[119] + x26*M[0] + x260*M[121] + x269*M[155] + x272*M[156] + x284*M[25] + x286*M[20] + x292*M[40] + x297*M[35] + x302*M[37] + x305*M[44] + x306*M[38] + x307*M[45] + x318*M[65] + x325*M[58] + x326*M[59] + x329*M[66] - x33*M[1] + x330*M[71] + x332*M[62] + x339*M[90] + x343*M[99] + x354*M[17] + x357*M[14] + x362*M[24] + x363*M[31] + x365*M[32] + x368*M[28] + x371*M[43] + x372*M[52] + x376*M[56] + x379*M[70] + x381*M[27] + x384*M[84] + x392*M[98] + x393*M[104] + x398*M[86] + x399*M[105] + x40*M[3] + x400*M[87] + x409*M[122] + x414*M[140] + x418*M[141] + x419*M[123] + x421*M[126] + x422*M[148] + x444*M[61] + x450*M[89] + x453*M[93] + x454*M[94] + x459*M[129] + x462*M[130] + x463*M[135] + x47*M[5] + x471*M[51] + x477*M[39] + x487*M[33] + x49*M[6] + x494*M[147] + x498*M[120] + x501*M[60] + x504*M[78] + x510*M[79] + x516*M[64] + x517*M[54] + x519*M[114] + x52*M[7] + x522*M[92] + x533*M[72] + x535*M[63] + x542*M[42] + x543*M[46] + x544*M[47] + x553*M[134] + x560*M[125] + x561*M[67] + x562*M[68] + x563*M[73] + x564*M[101] + x569*M[113] + x575*M[88] + x58*M[11] + x584*M[124] + x587*M[157] + x590*M[158] + x596*M[128] + x600*M[74] + x608*M[127] - x61*M[16] + x611*M[149] + x622*M[48] + x629*M[136] + x630*M[53] + x631*M[69] + x634*M[75] + x635*M[80] + x636*M[81] + x637*M[103] + x638*M[116] + x640*M[91] + x641*M[100] + x643*M[106] + x644*M[95] + x646*M[107] + x647*M[96] + x651*M[151] + x654*M[138] + x656*M[131] + x657*M[142] + x659*M[143] + x660*M[132] + x662*M[137] + x663*M[150] + x674*M[153] + x675*M[82] + x678*M[118] + x691*M[97] + x695*M[115] + x696*M[133] + x698*M[139] + x700*M[159] + x701*M[160] + x702*M[102] + x703*M[108] + x704*M[109] + x705*M[144] + x706*M[145] + x707*M[152] + x710*M[117] + x713*M[110] + x718*M[146] + x719*M[161] + x721*M[162] + x724*M[154] + x730*M[163] + x79*M[9] - x84*M[15] + x91*M[19] + x95*M[21];
#pragma omp atomic
L[6] += x103*M[29] + x112*M[35] - x116*M[49] - x13*M[1] + x137*M[34] + x154*M[55] + x169*M[56] + x176*M[76] + x187*M[84] + x193*M[12] - x199*M[111] + x200*M[2] + x211*M[22] + x212*M[25] + x218*M[40] + x219*M[7] + x220*M[8] + x230*M[17] + x242*M[83] + x256*M[119] + x260*M[120] + x268*M[5] + x272*M[155] + x30*M[0] + x306*M[37] + x307*M[44] + x326*M[58] + x329*M[65] + x330*M[70] + x332*M[61] + x339*M[89] + x343*M[98] + x348*M[21] + x353*M[11] + x354*M[16] + x357*M[13] + x359*M[14] + x362*M[23] + x363*M[30] + x365*M[31] + x366*M[24] + x368*M[27] + x371*M[42] + x372*M[51] + x381*M[26] + x399*M[104] + x400*M[86] + x41*M[3] + x418*M[140] + x419*M[122] + x421*M[125] + x422*M[147] + x431*M[57] + x435*M[28] + x454*M[93] + x462*M[129] + x463*M[134] + x466*M[18] + x47*M[4] + x471*M[50] + x474*M[36] + x477*M[38] + x483*M[39] + x487*M[32] + x488*M[33] + x501*M[59] + x504*M[77] + x507*M[60] + x510*M[78] + x516*M[63] + x517*M[53] + x519*M[113] + x52*M[6] + x522*M[91] + x527*M[121] + x533*M[71] + x535*M[62] + x539*M[64] + x542*M[41] + x543*M[45] + x544*M[46] + x561*M[66] + x562*M[67] + x563*M[72] + x564*M[100] + x567*M[85] + x569*M[112] + x575*M[87] + x578*M[88] + x58*M[10] + x584*M[123] + x587*M[156] + x590*M[157] + x593*M[124] + x596*M[127] + x600*M[73] + x604*M[75] + x608*M[126] - x61*M[15] + x611*M[148] + x620*M[128] + x621*M[43] + x622*M[47] + x624*M[48] + x629*M[135] + x630*M[52] + x631*M[68] + x632*M[69] + x634*M[74] + x635*M[79] + x636*M[80] + x637*M[102] + x638*M[115] + x640*M[90] + x641*M[99] + x643*M[105] + x644*M[94] + x646*M[106] + x647*M[95] + x651*M[150] + x654*M[137] + x656*M[130] + x657*M[141] + x659*M[142] + x660*M[131] + x662*M[136] + x663*M[149] + x668*M[139] + x671*M[54] + x674*M[152] + x675*M[81] + x676*M[82] + x678*M[117] + x684*M[154] + x689*M[92] + x691*M[96] + x692*M[97] + x695*M[114] + x696*M[132] + x697*M[133] + x698*M[138] + x700*M[158] + x701*M[159] + x702*M[101] + x703*M[107] + x704*M[108] + x705*M[143] + x706*M[144] + x707*M[151] + x709*M[103] + x710*M[116] + x713*M[109] + x715*M[110] + x718*M[145] + x719*M[160] + x721*M[161] + x722*M[146] + x724*M[153] + x728*M[118] + x730*M[162] + x731*M[163] + x80*M[9] + x92*M[19] + x95*M[20];
#pragma omp atomic
L[7] += x100*M[30] + x101*M[25] + x103*M[26] - x116*M[45] + x123*M[29] + x136*M[49] + x141*M[50] - x142*M[44] + x163*M[77] + x171*M[3] + x172*M[70] + x176*M[71] + x190*M[9] + x192*M[10] + x193*M[11] - x199*M[105] + x210*M[20] + x211*M[21] + x212*M[23] + x218*M[38] + x229*M[76] + x237*M[111] + x247*M[112] - x248*M[104] + x264*M[156] + x267*M[8] + x269*M[147] + x272*M[148] + x277*M[155] + x284*M[22] + x286*M[19] + x29*M[1] + x292*M[37] + x297*M[34] + x302*M[35] + x305*M[40] + x306*M[36] + x307*M[41] + x31*M[2] + x318*M[61] + x325*M[56] + x326*M[57] + x329*M[62] - x33*M[0] + x330*M[66] + x332*M[59] + x339*M[87] + x343*M[94] + x351*M[31] + x354*M[14] + x356*M[17] + x360*M[18] + x363*M[27] + x365*M[28] + x369*M[32] + x372*M[47] + x376*M[55] + x379*M[65] + x381*M[24] + x384*M[83] + x392*M[93] + x393*M[98] + x398*M[84] + x399*M[99] + x400*M[85] + x409*M[120] + x414*M[134] + x418*M[135] + x419*M[121] + x421*M[123] + x422*M[141] + x428*M[78] + x434*M[33] + x444*M[58] + x45*M[7] + x450*M[86] + x453*M[89] + x454*M[90] + x459*M[125] + x462*M[126] + x463*M[130] + x471*M[46] + x480*M[51] + x486*M[52] + x49*M[4] + x494*M[140] + x498*M[119] + x504*M[72] + x510*M[73] + x513*M[79] + x519*M[107] + x52*M[5] + x530*M[157] + x533*M[67] + x535*M[60] + x541*M[80] + x542*M[39] + x543*M[42] + x544*M[43] + x553*M[129] + x560*M[122] + x561*M[63] + x562*M[64] + x563*M[68] + x564*M[96] + x569*M[106] + x572*M[113] + x581*M[114] + x587*M[149] + x590*M[150] + x599*M[158] + x600*M[69] + x603*M[82] + x608*M[124] - x61*M[13] + x611*M[142] + x618*M[159] + x623*M[53] + x625*M[54] + x629*M[131] + x630*M[48] + x633*M[81] + x635*M[74] + x636*M[75] + x638*M[109] + x640*M[88] + x641*M[95] + x643*M[100] + x644*M[91] + x646*M[101] + x647*M[92] + x651*M[144] + x654*M[133] + x656*M[127] + x657*M[136] + x659*M[137] + x660*M[128] + x662*M[132] + x663*M[143] + x670*M[161] + x674*M[146] + x683*M[163] + x69*M[6] + x690*M[115] + x693*M[116] + x695*M[108] + x699*M[160] + x700*M[151] + x701*M[152] + x702*M[97] + x703*M[102] + x704*M[103] + x705*M[138] + x706*M[139] + x707*M[145] + x710*M[110] + x712*M[117] + x716*M[118] + x719*M[153] + x721*M[154] + x725*M[162] + x77*M[15] + x83*M[16] - x84*M[12];
#pragma omp atomic
L[8] += x100*M[29] + x103*M[25] + x112*M[34] - x116*M[44] - x13*M[0] + x141*M[49] + x163*M[76] + x169*M[55] + x176*M[70] + x187*M[83] + x193*M[10] - x199*M[104] + x201*M[2] + x211*M[20] + x212*M[22] + x218*M[37] + x219*M[5] + x221*M[8] + x230*M[14] + x247*M[111] + x260*M[119] + x264*M[155] + x267*M[7] + x272*M[147] + x306*M[35] + x307*M[40] + x31*M[1] + x326*M[56] + x329*M[61] + x330*M[65] + x332*M[58] + x339*M[86] + x343*M[93] + x351*M[30] + x354*M[13] + x356*M[16] + x357*M[11] + x360*M[17] + x362*M[21] + x363*M[26] + x365*M[27] + x368*M[24] + x369*M[31] + x371*M[39] + x372*M[46] + x381*M[23] + x399*M[98] + x400*M[84] + x418*M[134] + x419*M[120] + x421*M[122] + x422*M[140] + x428*M[77] + x434*M[32] + x45*M[6] + x454*M[89] + x462*M[125] + x463*M[129] + x464*M[18] + x469*M[33] + x47*M[3] + x471*M[45] + x477*M[36] + x480*M[50] + x486*M[51] + x487*M[28] + x501*M[57] + x504*M[71] + x510*M[72] + x513*M[78] + x516*M[60] + x517*M[48] + x519*M[106] + x52*M[4] + x522*M[88] + x530*M[156] + x533*M[66] + x535*M[59] + x541*M[79] + x542*M[38] + x543*M[41] + x544*M[42] + x561*M[62] + x562*M[63] + x563*M[67] + x564*M[95] + x569*M[105] + x572*M[112] + x575*M[85] + x58*M[9] + x581*M[113] + x584*M[121] + x587*M[148] + x590*M[149] + x596*M[124] + x599*M[157] + x600*M[68] + x603*M[81] + x608*M[123] - x61*M[12] + x611*M[141] + x618*M[158] + x622*M[43] + x623*M[52] + x625*M[53] + x629*M[130] + x630*M[47] + x631*M[64] + x633*M[80] + x634*M[69] + x635*M[73] + x636*M[74] + x637*M[97] + x638*M[108] + x640*M[87] + x641*M[94] + x643*M[99] + x644*M[90] + x646*M[100] + x647*M[91] + x651*M[143] + x654*M[132] + x656*M[126] + x657*M[135] + x659*M[136] + x660*M[127] + x662*M[131] + x663*M[142] + x670*M[160] + x672*M[54] + x674*M[145] + x675*M[75] + x677*M[82] + x678*M[110] + x683*M[162] + x690*M[114] + x691*M[92] + x693*M[115] + x695*M[107] + x696*M[128] + x698*M[133] + x699*M[159] + x700*M[150] + x701*M[151] + x702*M[96] + x703*M[101] + x704*M[102] + x705*M[137] + x706*M[138] + x707*M[144] + x710*M[109] + x712*M[116] + x713*M[103] + x716*M[117] + x718*M[139] + x719*M[152] + x721*M[153] + x724*M[146] + x725*M[161] + x726*M[118] + x729*M[163] + x730*M[154] + x83*M[15] + x95*M[19];
#pragma omp atomic
L[9] += x200*M[0] + x201*M[1] + x202*M[2] + x219*M[4] + x220*M[5] + x221*M[7] + x230*M[13] + x267*M[6] + x268*M[3] + x279*M[8] + x348*M[19] + x351*M[29] + x353*M[9] + x354*M[12] + x356*M[15] + x357*M[10] + x359*M[11] + x360*M[16] + x362*M[20] + x363*M[25] + x365*M[26] + x366*M[21] + x368*M[23] + x369*M[30] + x371*M[38] + x372*M[45] + x381*M[22] + x428*M[76] + x431*M[55] + x434*M[31] + x435*M[24] + x464*M[17] + x465*M[18] + x466*M[14] + x468*M[33] + x469*M[32] + x471*M[44] + x474*M[34] + x477*M[35] + x480*M[49] + x483*M[36] + x486*M[50] + x487*M[27] + x488*M[28] + x501*M[56] + x504*M[70] + x507*M[57] + x510*M[71] + x513*M[77] + x516*M[59] + x517*M[47] + x519*M[105] + x522*M[87] + x527*M[119] + x530*M[155] + x533*M[65] + x535*M[58] + x539*M[60] + x541*M[78] + x542*M[37] + x543*M[40] + x544*M[41] + x561*M[61] + x562*M[62] + x563*M[66] + x564*M[94] + x567*M[83] + x569*M[104] + x572*M[111] + x575*M[84] + x578*M[85] + x581*M[112] + x584*M[120] + x587*M[147] + x590*M[148] + x593*M[121] + x596*M[123] + x599*M[156] + x600*M[67] + x603*M[80] + x604*M[69] + x608*M[122] + x611*M[140] + x613*M[82] + x618*M[157] + x620*M[124] + x621*M[39] + x622*M[42] + x623*M[51] + x624*M[43] + x625*M[52] + x629*M[129] + x630*M[46] + x631*M[63] + x632*M[64] + x633*M[79] + x634*M[68] + x635*M[72] + x636*M[73] + x637*M[96] + x638*M[107] + x640*M[86] + x641*M[93] + x643*M[98] + x644*M[89] + x646*M[99] + x647*M[90] + x651*M[142] + x654*M[131] + x656*M[125] + x657*M[134] + x659*M[135] + x660*M[126] + x662*M[130] + x663*M[141] + x668*M[133] + x670*M[159] + x671*M[48] + x672*M[53] + x673*M[54] + x674*M[144] + x675*M[74] + x676*M[75] + x677*M[81] + x678*M[109] + x683*M[161] + x684*M[146] + x688*M[163] + x689*M[88] + x690*M[113] + x691*M[91] + x692*M[92] + x693*M[114] + x695*M[106] + x696*M[127] + x697*M[128] + x698*M[132] + x699*M[158] + x700*M[149] + x701*M[150] + x702*M[95] + x703*M[100] + x704*M[101] + x705*M[136] + x706*M[137] + x707*M[143] + x709*M[97] + x710*M[108] + x712*M[115] + x713*M[102] + x715*M[103] + x716*M[116] + x718*M[138] + x719*M[151] + x721*M[152] + x722*M[139] + x724*M[145] + x725*M[160] + x726*M[117] + x727*M[118] + x728*M[110] + x729*M[162] + x730*M[153] + x731*M[154];
#pragma omp atomic
L[10] += x112*M[23] + x120*M[9] + x129*M[19] + x132*M[20] + x137*M[21] + x153*M[35] + x154*M[36] + x169*M[38] + x187*M[59] + x190*M[6] + x210*M[15] + x211*M[16] + x218*M[30] + x226*M[34] + x234*M[55] + x241*M[56] + x242*M[57] + x255*M[84] + x256*M[85] + x260*M[87] + x275*M[83] + x286*M[12] + x292*M[29] + x297*M[22] + x302*M[25] + x306*M[26] + x318*M[49] + x325*M[40] + x326*M[41] + x329*M[50] + x332*M[45] + x339*M[66] + x343*M[77] + x348*M[14] + x353*M[8] + x362*M[17] + x366*M[18] + x371*M[32] + x376*M[37] + x384*M[58] + x392*M[76] + x398*M[61] + x40*M[1] + x400*M[62] + x409*M[89] + x41*M[2] + x414*M[111] + x418*M[112] + x419*M[90] + x421*M[94] + x431*M[39] + x444*M[44] + x450*M[65] + x453*M[70] + x454*M[71] + x459*M[98] + x462*M[99] + x463*M[105] + x474*M[24] + x477*M[27] + x483*M[28] + x498*M[86] + x501*M[42] + x507*M[43] + x516*M[47] + x522*M[68] + x527*M[88] + x535*M[46] + x539*M[48] + x542*M[31] + x553*M[104] + x560*M[93] + x561*M[51] + x562*M[52] + x564*M[79] + x567*M[60] + x575*M[63] + x578*M[64] + x58*M[7] + x584*M[91] + x593*M[92] + x596*M[96] + x608*M[95] + x620*M[97] + x621*M[33] + x629*M[106] + x631*M[53] + x632*M[54] + x637*M[81] + x640*M[67] + x641*M[78] + x644*M[72] + x647*M[73] + x654*M[108] + x656*M[100] + x657*M[113] + x659*M[114] + x66*M[0] + x660*M[101] + x662*M[107] + x668*M[110] + x689*M[69] + x691*M[74] + x692*M[75] + x696*M[102] + x697*M[103] + x698*M[109] + x702*M[80] + x705*M[115] + x706*M[116] + x709*M[82] + x718*M[117] + x722*M[118] + x74*M[3] + x79*M[4] + x80*M[5] + x91*M[10] + x92*M[11] + x95*M[13];
#pragma omp atomic
L[11] += x112*M[21] + x132*M[19] + x153*M[34] + x169*M[36] + x171*M[1] + x187*M[57] + x190*M[4] + x192*M[6] + x193*M[7] + x210*M[12] + x211*M[13] + x212*M[16] + x218*M[26] + x241*M[55] + x255*M[83] + x260*M[85] + x284*M[15] + x286*M[10] + x292*M[25] + x297*M[20] + x302*M[22] + x305*M[29] + x306*M[23] + x307*M[30] + x318*M[44] + x325*M[37] + x326*M[38] + x329*M[45] + x330*M[50] + x332*M[41] + x339*M[62] + x343*M[71] + x357*M[8] + x362*M[14] + x368*M[18] + x371*M[28] + x376*M[35] + x379*M[49] + x381*M[17] + x384*M[56] + x392*M[70] + x393*M[76] + x398*M[58] + x399*M[77] + x40*M[0] + x400*M[59] + x409*M[86] + x414*M[104] + x418*M[105] + x419*M[87] + x421*M[90] + x422*M[112] + x444*M[40] + x450*M[61] + x453*M[65] + x454*M[66] + x459*M[93] + x462*M[94] + x463*M[99] + x47*M[2] + x477*M[24] + x494*M[111] + x498*M[84] + x501*M[39] + x516*M[43] + x522*M[64] + x533*M[51] + x535*M[42] + x542*M[27] + x543*M[31] + x544*M[32] + x553*M[98] + x560*M[89] + x561*M[46] + x562*M[47] + x563*M[52] + x564*M[73] + x575*M[60] + x58*M[5] + x584*M[88] + x596*M[92] + x600*M[53] + x608*M[91] + x611*M[113] + x622*M[33] + x629*M[100] + x631*M[48] + x634*M[54] + x637*M[75] + x640*M[63] + x641*M[72] + x643*M[78] + x644*M[67] + x646*M[79] + x647*M[68] + x651*M[115] + x654*M[102] + x656*M[95] + x657*M[106] + x659*M[107] + x660*M[96] + x662*M[101] + x663*M[114] + x674*M[117] + x691*M[69] + x696*M[97] + x698*M[103] + x702*M[74] + x703*M[80] + x704*M[81] + x705*M[108] + x706*M[109] + x707*M[116] + x713*M[82] + x718*M[110] + x724*M[118] + x79*M[3] + x91*M[9] + x95*M[11];
#pragma omp atomic
L[12] += x112*M[20] + x137*M[19] + x154*M[34] + x169*M[35] + x187*M[56] + x193*M[6] + x211*M[12] + x212*M[15] + x218*M[25] + x242*M[55] + x256*M[83] + x260*M[84] + x268*M[2] + x306*M[22] + x307*M[29] + x326*M[37] + x329*M[44] + x330*M[49] + x332*M[40] + x339*M[61] + x343*M[70] + x348*M[11] + x353*M[5] + x357*M[7] + x359*M[8] + x362*M[13] + x366*M[14] + x368*M[17] + x371*M[27] + x381*M[16] + x399*M[76] + x400*M[58] + x41*M[0] + x418*M[104] + x419*M[86] + x421*M[89] + x422*M[111] + x431*M[36] + x435*M[18] + x454*M[65] + x462*M[93] + x463*M[98] + x47*M[1] + x474*M[21] + x477*M[23] + x483*M[24] + x501*M[38] + x507*M[39] + x516*M[42] + x522*M[63] + x527*M[85] + x533*M[50] + x535*M[41] + x539*M[43] + x542*M[26] + x543*M[30] + x544*M[31] + x561*M[45] + x562*M[46] + x563*M[51] + x564*M[72] + x567*M[57] + x575*M[59] + x578*M[60] + x58*M[4] + x584*M[87] + x593*M[88] + x596*M[91] + x600*M[52] + x604*M[54] + x608*M[90] + x611*M[112] + x620*M[92] + x621*M[28] + x622*M[32] + x624*M[33] + x629*M[99] + x631*M[47] + x632*M[48] + x634*M[53] + x637*M[74] + x640*M[62] + x641*M[71] + x643*M[77] + x644*M[66] + x646*M[78] + x647*M[67] + x651*M[114] + x654*M[101] + x656*M[94] + x657*M[105] + x659*M[106] + x660*M[95] + x662*M[100] + x663*M[113] + x668*M[103] + x674*M[116] + x684*M[118] + x689*M[64] + x691*M[68] + x692*M[69] + x696*M[96] + x697*M[97] + x698*M[102] + x702*M[73] + x703*M[79] + x704*M[80] + x705*M[107] + x706*M[108] + x707*M[115] + x709*M[75] + x713*M[81] + x715*M[82] + x718*M[109] + x722*M[110] + x724*M[117] + x80*M[3] + x92*M[9] + x95*M[10];
#pragma omp atomic
L[13] += x101*M[15] + x103*M[16] - x116*M[30] - x142*M[29] + x171*M[0] + x172*M[49] + x176*M[50] + x190*M[3] + x192*M[4] + x193*M[5] - x199*M[77] + x210*M[10] + x211*M[11] + x212*M[13] + x218*M[23] - x248*M[76] + x269*M[111] + x272*M[112] + x284*M[12] + x286*M[9] + x292*M[22] + x297*M[19] + x302*M[20] + x305*M[25] + x306*M[21] + x307*M[26] + x318*M[40] + x325*M[35] + x326*M[36] + x329*M[41] + x330*M[45] + x332*M[38] + x339*M[59] + x343*M[66] + x354*M[8] + x363*M[17] + x365*M[18] + x372*M[32] + x376*M[34] + x379*M[44] + x381*M[14] + x384*M[55] + x392*M[65] + x393*M[70] + x398*M[56] + x399*M[71] + x400*M[57] + x409*M[84] + x414*M[98] + x418*M[99] + x419*M[85] + x421*M[87] + x422*M[105] + x444*M[37] + x450*M[58] + x453*M[61] + x454*M[62] + x459*M[89] + x462*M[90] + x463*M[94] + x471*M[31] + x49*M[1] + x494*M[104] + x498*M[83] + x504*M[51] + x510*M[52] + x519*M[79] + x52*M[2] + x533*M[46] + x535*M[39] + x542*M[24] + x543*M[27] + x544*M[28] + x553*M[93] + x560*M[86] + x561*M[42] + x562*M[43] + x563*M[47] + x564*M[68] + x569*M[78] + x587*M[113] + x590*M[114] + x600*M[48] + x608*M[88] - x61*M[7] + x611*M[106] + x629*M[95] + x630*M[33] + x635*M[53] + x636*M[54] + x638*M[81] + x640*M[60] + x641*M[67] + x643*M[72] + x644*M[63] + x646*M[73] + x647*M[64] + x651*M[108] + x654*M[97] + x656*M[91] + x657*M[100] + x659*M[101] + x660*M[92] + x662*M[96] + x663*M[107] + x674*M[110] + x695*M[80] + x700*M[115] + x701*M[116] + x702*M[69] + x703*M[74] + x704*M[75] + x705*M[102] + x706*M[103] + x707*M[109] + x710*M[82] + x719*M[117] + x721*M[118] - x84*M[6];
#pragma omp atomic
L[14] += x103*M[15] + x112*M[19] - x116*M[29] + x169*M[34] + x176*M[49] + x187*M[55] + x193*M[4] - x199*M[76] + x211*M[10] + x212*M[12] + x218*M[22] + x219*M[2] + x230*M[8] + x260*M[83] + x272*M[111] + x306*M[20] + x307*M[25] + x326*M[35] + x329*M[40] + x330*M[44] + x332*M[37] + x339*M[58] + x343*M[65] + x354*M[7] + x357*M[5] + x362*M[11] + x363*M[16] + x365*M[17] + x368*M[14] + x371*M[24] + x372*M[31] + x381*M[13] + x399*M[70] + x400*M[56] + x418*M[98] + x419*M[84] + x421*M[86] + x422*M[104] + x454*M[61] + x462*M[89] + x463*M[93] + x47*M[0] + x471*M[30] + x477*M[21] + x487*M[18] + x501*M[36] + x504*M[50] + x510*M[51] + x516*M[39] + x517*M[33] + x519*M[78] + x52*M[1] + x522*M[60] + x533*M[45] + x535*M[38] + x542*M[23] + x543*M[26] + x544*M[27] + x561*M[41] + x562*M[42] + x563*M[46] + x564*M[67] + x569*M[77] + x575*M[57] + x58*M[3] + x584*M[85] + x587*M[112] + x590*M[113] + x596*M[88] + x600*M[47] + x608*M[87] - x61*M[6] + x611*M[105] + x622*M[28] + x629*M[94] + x630*M[32] + x631*M[43] + x634*M[48] + x635*M[52] + x636*M[53] + x637*M[69] + x638*M[80] + x640*M[59] + x641*M[66] + x643*M[71] + x644*M[62] + x646*M[72] + x647*M[63] + x651*M[107] + x654*M[96] + x656*M[90] + x657*M[99] + x659*M[100] + x660*M[91] + x662*M[95] + x663*M[106] + x674*M[109] + x675*M[54] + x678*M[82] + x691*M[64] + x695*M[79] + x696*M[92] + x698*M[97] + x700*M[114] + x701*M[115] + x702*M[68] + x703*M[73] + x704*M[74] + x705*M[101] + x706*M[102] + x707*M[108] + x710*M[81] + x713*M[75] + x718*M[103] + x719*M[116] + x721*M[117] + x724*M[110] + x730*M[118] + x95*M[9];
#pragma omp atomic
L[15] += x219*M[1] + x220*M[2] + x230*M[7] + x268*M[0] + x348*M[9] + x353*M[3] + x354*M[6] + x357*M[4] + x359*M[5] + x362*M[10] + x363*M[15] + x365*M[16] + x366*M[11] + x368*M[13] + x371*M[23] + x372*M[30] + x381*M[12] + x431*M[34] + x435*M[14] + x466*M[8] + x471*M[29] + x474*M[19] + x477*M[20] + x483*M[21] + x487*M[17] + x488*M[18] + x501*M[35] + x504*M[49] + x507*M[36] + x510*M[50] + x516*M[38] + x517*M[32] + x519*M[77] + x522*M[59] + x527*M[83] + x533*M[44] + x535*M[37] + x539*M[39] + x542*M[22] + x543*M[25] + x544*M[26] + x561*M[40] + x562*M[41] + x563*M[45] + x564*M[66] + x567*M[55] + x569*M[76] + x575*M[56] + x578*M[57] + x584*M[84] + x587*M[111] + x590*M[112] + x593*M[85] + x596*M[87] + x600*M[46] + x604*M[48] + x608*M[86] + x611*M[104] + x620*M[88] + x621*M[24] + x622*M[27] + x624*M[28] + x629*M[93] + x630*M[31] + x631*M[42] + x632*M[43] + x634*M[47] + x635*M[51] + x636*M[52] + x637*M[68] + x638*M[79] + x640*M[58] + x641*M[65] + x643*M[70] + x644*M[61] + x646*M[71] + x647*M[62] + x651*M[106] + x654*M[95] + x656*M[89] + x657*M[98] + x659*M[99] + x660*M[90] + x662*M[94] + x663*M[105] + x668*M[97] + x671*M[33] + x674*M[108] + x675*M[53] + x676*M[54] + x678*M[81] + x684*M[110] + x689*M[60] + x691*M[63] + x692*M[64] + x695*M[78] + x696*M[91] + x697*M[92] + x698*M[96] + x700*M[113] + x701*M[114] + x702*M[67] + x703*M[72] + x704*M[73] + x705*M[100] + x706*M[101] + x707*M[107] + x709*M[69] + x710*M[80] + x713*M[74] + x715*M[75] + x718*M[102] + x719*M[115] + x721*M[116] + x722*M[103] + x724*M[109] + x728*M[82] + x730*M[117] + x731*M[118];
#pragma omp atomic
L[16] += x100*M[16] + x101*M[12] + x103*M[13] - x116*M[26] + x123*M[15] + x136*M[29] + x141*M[30] - x142*M[25] + x163*M[50] + x172*M[44] + x176*M[45] + x192*M[3] - x199*M[71] + x210*M[9] + x212*M[11] + x218*M[21] + x229*M[49] + x237*M[76] + x247*M[77] - x248*M[70] + x264*M[112] + x269*M[104] + x272*M[105] + x277*M[111] + x284*M[10] + x292*M[20] + x302*M[19] + x305*M[22] + x307*M[23] + x318*M[37] + x325*M[34] + x329*M[38] + x330*M[41] + x332*M[36] + x339*M[57] + x343*M[62] + x351*M[17] + x356*M[8] + x363*M[14] + x369*M[18] + x372*M[28] + x379*M[40] + x392*M[61] + x393*M[65] + x398*M[55] + x399*M[66] + x409*M[83] + x414*M[93] + x418*M[94] + x421*M[85] + x422*M[99] + x428*M[51] + x444*M[35] + x45*M[2] + x450*M[56] + x453*M[58] + x454*M[59] + x459*M[86] + x462*M[87] + x463*M[90] + x471*M[27] + x480*M[31] + x486*M[32] + x49*M[0] + x494*M[98] + x504*M[46] + x510*M[47] + x513*M[52] + x519*M[73] + x530*M[113] + x533*M[42] + x541*M[53] + x543*M[24] + x553*M[89] + x560*M[84] + x561*M[39] + x563*M[43] + x564*M[64] + x569*M[72] + x572*M[78] + x581*M[79] + x587*M[106] + x590*M[107] + x599*M[114] - x61*M[5] + x611*M[100] + x618*M[115] + x623*M[33] + x629*M[91] + x633*M[54] + x635*M[48] + x638*M[75] + x641*M[63] + x643*M[67] + x644*M[60] + x646*M[68] + x651*M[102] + x656*M[88] + x657*M[95] + x659*M[96] + x662*M[92] + x663*M[101] + x670*M[117] + x69*M[1] + x690*M[80] + x693*M[81] + x695*M[74] + x699*M[116] + x700*M[108] + x701*M[109] + x703*M[69] + x705*M[97] + x707*M[103] + x712*M[82] + x719*M[110] + x725*M[118] + x77*M[6] + x83*M[7] - x84*M[4];
#pragma omp atomic
L[17] += x100*M[15] + x103*M[12] - x116*M[25] + x141*M[29] + x163*M[49] + x176*M[44] + x193*M[3] - x199*M[70] + x211*M[9] + x212*M[10] + x218*M[20] + x247*M[76] + x264*M[111] + x267*M[2] + x272*M[104] + x306*M[19] + x307*M[22] + x326*M[34] + x329*M[37] + x330*M[40] + x332*M[35] + x339*M[56] + x343*M[61] + x351*M[16] + x354*M[5] + x356*M[7] + x360*M[8] + x363*M[13] + x365*M[14] + x369*M[17] + x372*M[27] + x381*M[11] + x399*M[65] + x400*M[55] + x418*M[93] + x419*M[83] + x421*M[84] + x422*M[98] + x428*M[50] + x434*M[18] + x45*M[1] + x454*M[58] + x462*M[86] + x463*M[89] + x471*M[26] + x480*M[30] + x486*M[31] + x504*M[45] + x510*M[46] + x513*M[51] + x519*M[72] + x52*M[0] + x530*M[112] + x533*M[41] + x535*M[36] + x541*M[52] + x542*M[21] + x543*M[23] + x544*M[24] + x561*M[38] + x562*M[39] + x563*M[42] + x564*M[63] + x569*M[71] + x572*M[77] + x581*M[78] + x587*M[105] + x590*M[106] + x599*M[113] + x600*M[43] + x603*M[54] + x608*M[85] - x61*M[4] + x611*M[99] + x618*M[114] + x623*M[32] + x625*M[33] + x629*M[90] + x630*M[28] + x633*M[53] + x635*M[47] + x636*M[48] + x638*M[74] + x640*M[57] + x641*M[62] + x643*M[66] + x644*M[59] + x646*M[67] + x647*M[60] + x651*M[101] + x654*M[92] + x656*M[87] + x657*M[94] + x659*M[95] + x660*M[88] + x662*M[91] + x663*M[100] + x670*M[116] + x674*M[103] + x683*M[118] + x690*M[79] + x693*M[80] + x695*M[73] + x699*M[115] + x700*M[107] + x701*M[108] + x702*M[64] + x703*M[68] + x704*M[69] + x705*M[96] + x706*M[97] + x707*M[102] + x710*M[75] + x712*M[81] + x716*M[82] + x719*M[109] + x721*M[110] + x725*M[117] + x83*M[6];
#pragma omp atomic
L[18] += x219*M[0] + x221*M[2] + x230*M[5] + x267*M[1] + x351*M[15] + x354*M[4] + x356*M[6] + x357*M[3] + x360*M[7] + x362*M[9] + x363*M[12] + x365*M[13] + x368*M[11] + x369*M[16] + x371*M[21] + x372*M[26] + x381*M[10] + x428*M[49] + x434*M[17] + x464*M[8] + x469*M[18] + x471*M[25] + x477*M[19] + x480*M[29] + x486*M[30] + x487*M[14] + x501*M[34] + x504*M[44] + x510*M[45] + x513*M[50] + x516*M[36] + x517*M[28] + x519*M[71] + x522*M[57] + x530*M[111] + x533*M[40] + x535*M[35] + x541*M[51] + x542*M[20] + x543*M[22] + x544*M[23] + x561*M[37] + x562*M[38] + x563*M[41] + x564*M[62] + x569*M[70] + x572*M[76] + x575*M[55] + x581*M[77] + x584*M[83] + x587*M[104] + x590*M[105] + x596*M[85] + x599*M[112] + x600*M[42] + x603*M[53] + x608*M[84] + x611*M[98] + x618*M[113] + x622*M[24] + x623*M[31] + x625*M[32] + x629*M[89] + x630*M[27] + x631*M[39] + x633*M[52] + x634*M[43] + x635*M[46] + x636*M[47] + x637*M[64] + x638*M[73] + x640*M[56] + x641*M[61] + x643*M[65] + x644*M[58] + x646*M[66] + x647*M[59] + x651*M[100] + x654*M[91] + x656*M[86] + x657*M[93] + x659*M[94] + x660*M[87] + x662*M[90] + x663*M[99] + x670*M[115] + x672*M[33] + x674*M[102] + x675*M[48] + x677*M[54] + x678*M[75] + x683*M[117] + x690*M[78] + x691*M[60] + x693*M[79] + x695*M[72] + x696*M[88] + x698*M[92] + x699*M[114] + x700*M[106] + x701*M[107] + x702*M[63] + x703*M[67] + x704*M[68] + x705*M[95] + x706*M[96] + x707*M[101] + x710*M[74] + x712*M[80] + x713*M[69] + x716*M[81] + x718*M[97] + x719*M[108] + x721*M[109] + x724*M[103] + x725*M[116] + x726*M[82] + x729*M[118] + x730*M[110];
#pragma omp atomic
L[19] += x220*M[0] + x221*M[1] + x230*M[4] + x279*M[2] + x359*M[3] + x360*M[6] + x365*M[12] + x366*M[9] + x368*M[10] + x369*M[15] + x371*M[20] + x372*M[25] + x434*M[16] + x435*M[11] + x464*M[7] + x465*M[8] + x466*M[5] + x468*M[18] + x469*M[17] + x483*M[19] + x486*M[29] + x487*M[13] + x488*M[14] + x507*M[34] + x510*M[44] + x513*M[49] + x516*M[35] + x517*M[27] + x519*M[70] + x522*M[56] + x539*M[36] + x541*M[50] + x544*M[22] + x562*M[37] + x563*M[40] + x564*M[61] + x578*M[55] + x581*M[76] + x590*M[104] + x593*M[83] + x596*M[84] + x599*M[111] + x600*M[41] + x603*M[52] + x604*M[43] + x613*M[54] + x618*M[112] + x620*M[85] + x621*M[21] + x622*M[23] + x623*M[30] + x624*M[24] + x625*M[31] + x630*M[26] + x631*M[38] + x632*M[39] + x633*M[51] + x634*M[42] + x635*M[45] + x636*M[46] + x637*M[63] + x638*M[72] + x646*M[65] + x647*M[58] + x651*M[99] + x654*M[90] + x659*M[93] + x660*M[86] + x662*M[89] + x663*M[98] + x668*M[92] + x670*M[114] + x671*M[28] + x672*M[32] + x673*M[33] + x674*M[101] + x675*M[47] + x676*M[48] + x677*M[53] + x678*M[74] + x683*M[116] + x684*M[103] + x688*M[118] + x689*M[57] + x690*M[77] + x691*M[59] + x692*M[60] + x693*M[78] + x695*M[71] + x696*M[87] + x697*M[88] + x698*M[91] + x699*M[113] + x700*M[105] + x701*M[106] + x702*M[62] + x703*M[66] + x704*M[67] + x705*M[94] + x706*M[95] + x707*M[100] + x709*M[64] + x710*M[73] + x712*M[79] + x713*M[68] + x715*M[69] + x716*M[80] + x718*M[96] + x719*M[107] + x721*M[108] + x722*M[97] + x724*M[102] + x725*M[115] + x726*M[81] + x727*M[82] + x728*M[75] + x729*M[117] + x730*M[109] + x731*M[110];
#pragma omp atomic
L[20] += x112*M[13] + x120*M[3] + x129*M[9] + x132*M[10] + x137*M[11] + x153*M[20] + x154*M[21] + x169*M[23] + x187*M[38] + x226*M[19] + x234*M[34] + x241*M[35] + x242*M[36] + x255*M[56] + x256*M[57] + x260*M[59] + x275*M[55] + x286*M[6] + x297*M[12] + x302*M[15] + x306*M[16] + x325*M[25] + x326*M[26] + x332*M[30] + x339*M[45] + x348*M[8] + x376*M[22] + x384*M[37] + x398*M[40] + x400*M[41] + x409*M[61] + x419*M[62] + x421*M[66] + x431*M[24] + x444*M[29] + x450*M[44] + x453*M[49] + x454*M[50] + x459*M[70] + x462*M[71] + x463*M[77] + x474*M[14] + x477*M[17] + x483*M[18] + x498*M[58] + x501*M[27] + x507*M[28] + x516*M[32] + x522*M[47] + x527*M[60] + x535*M[31] + x539*M[33] + x553*M[76] + x560*M[65] + x567*M[39] + x575*M[42] + x578*M[43] + x584*M[63] + x593*M[64] + x596*M[68] + x608*M[67] + x620*M[69] + x629*M[78] + x640*M[46] + x644*M[51] + x647*M[52] + x654*M[80] + x656*M[72] + x660*M[73] + x662*M[79] + x668*M[82] + x689*M[48] + x691*M[53] + x692*M[54] + x696*M[74] + x697*M[75] + x698*M[81] + x74*M[0] + x79*M[1] + x80*M[2] + x91*M[4] + x92*M[5] + x95*M[7];
#pragma omp atomic
L[21] += x112*M[11] + x132*M[9] + x153*M[19] + x169*M[21] + x187*M[36] + x190*M[1] + x210*M[6] + x211*M[7] + x218*M[16] + x241*M[34] + x255*M[55] + x260*M[57] + x286*M[4] + x292*M[15] + x297*M[10] + x302*M[12] + x306*M[13] + x318*M[29] + x325*M[22] + x326*M[23] + x329*M[30] + x332*M[26] + x339*M[41] + x343*M[50] + x362*M[8] + x371*M[18] + x376*M[20] + x384*M[35] + x392*M[49] + x398*M[37] + x400*M[38] + x409*M[58] + x414*M[76] + x418*M[77] + x419*M[59] + x421*M[62] + x444*M[25] + x450*M[40] + x453*M[44] + x454*M[45] + x459*M[65] + x462*M[66] + x463*M[71] + x477*M[14] + x498*M[56] + x501*M[24] + x516*M[28] + x522*M[43] + x535*M[27] + x542*M[17] + x553*M[70] + x560*M[61] + x561*M[31] + x562*M[32] + x564*M[52] + x575*M[39] + x58*M[2] + x584*M[60] + x596*M[64] + x608*M[63] + x629*M[72] + x631*M[33] + x637*M[54] + x640*M[42] + x641*M[51] + x644*M[46] + x647*M[47] + x654*M[74] + x656*M[67] + x657*M[78] + x659*M[79] + x660*M[68] + x662*M[73] + x691*M[48] + x696*M[69] + x698*M[75] + x702*M[53] + x705*M[80] + x706*M[81] + x718*M[82] + x79*M[0] + x91*M[3] + x95*M[5];
#pragma omp atomic
L[22] += x112*M[10] + x137*M[9] + x154*M[19] + x169*M[20] + x187*M[35] + x211*M[6] + x218*M[15] + x242*M[34] + x256*M[55] + x260*M[56] + x306*M[12] + x326*M[22] + x329*M[29] + x332*M[25] + x339*M[40] + x343*M[49] + x348*M[5] + x353*M[2] + x362*M[7] + x366*M[8] + x371*M[17] + x400*M[37] + x418*M[76] + x419*M[58] + x421*M[61] + x431*M[21] + x454*M[44] + x462*M[65] + x463*M[70] + x474*M[11] + x477*M[13] + x483*M[14] + x501*M[23] + x507*M[24] + x516*M[27] + x522*M[42] + x527*M[57] + x535*M[26] + x539*M[28] + x542*M[16] + x561*M[30] + x562*M[31] + x564*M[51] + x567*M[36] + x575*M[38] + x578*M[39] + x58*M[1] + x584*M[59] + x593*M[60] + x596*M[63] + x608*M[62] + x620*M[64] + x621*M[18] + x629*M[71] + x631*M[32] + x632*M[33] + x637*M[53] + x640*M[41] + x641*M[50] + x644*M[45] + x647*M[46] + x654*M[73] + x656*M[66] + x657*M[77] + x659*M[78] + x660*M[67] + x662*M[72] + x668*M[75] + x689*M[43] + x691*M[47] + x692*M[48] + x696*M[68] + x697*M[69] + x698*M[74] + x702*M[52] + x705*M[79] + x706*M[80] + x709*M[54] + x718*M[81] + x722*M[82] + x80*M[0] + x92*M[3] + x95*M[4];
#pragma omp atomic
L[23] += x190*M[0] + x192*M[1] + x193*M[2] + x210*M[4] + x211*M[5] + x212*M[7] + x218*M[13] + x284*M[6] + x286*M[3] + x292*M[12] + x297*M[9] + x302*M[10] + x305*M[15] + x306*M[11] + x307*M[16] + x318*M[25] + x325*M[20] + x326*M[21] + x329*M[26] + x330*M[30] + x332*M[23] + x339*M[38] + x343*M[45] + x376*M[19] + x379*M[29] + x381*M[8] + x384*M[34] + x392*M[44] + x393*M[49] + x398*M[35] + x399*M[50] + x400*M[36] + x409*M[56] + x414*M[70] + x418*M[71] + x419*M[57] + x421*M[59] + x422*M[77] + x444*M[22] + x450*M[37] + x453*M[40] + x454*M[41] + x459*M[61] + x462*M[62] + x463*M[66] + x494*M[76] + x498*M[55] + x533*M[31] + x535*M[24] + x542*M[14] + x543*M[17] + x544*M[18] + x553*M[65] + x560*M[58] + x561*M[27] + x562*M[28] + x563*M[32] + x564*M[47] + x600*M[33] + x608*M[60] + x611*M[78] + x629*M[67] + x640*M[39] + x641*M[46] + x643*M[51] + x644*M[42] + x646*M[52] + x647*M[43] + x651*M[80] + x654*M[69] + x656*M[63] + x657*M[72] + x659*M[73] + x660*M[64] + x662*M[68] + x663*M[79] + x674*M[82] + x702*M[48] + x703*M[53] + x704*M[54] + x705*M[74] + x706*M[75] + x707*M[81];
#pragma omp atomic
L[24] += x112*M[9] + x169*M[19] + x187*M[34] + x193*M[1] + x211*M[4] + x212*M[6] + x218*M[12] + x260*M[55] + x306*M[10] + x307*M[15] + x326*M[20] + x329*M[25] + x330*M[29] + x332*M[22] + x339*M[37] + x343*M[44] + x357*M[2] + x362*M[5] + x368*M[8] + x371*M[14] + x381*M[7] + x399*M[49] + x400*M[35] + x418*M[70] + x419*M[56] + x421*M[58] + x422*M[76] + x454*M[40] + x462*M[61] + x463*M[65] + x477*M[11] + x501*M[21] + x516*M[24] + x522*M[39] + x533*M[30] + x535*M[23] + x542*M[13] + x543*M[16] + x544*M[17] + x561*M[26] + x562*M[27] + x563*M[31] + x564*M[46] + x575*M[36] + x58*M[0] + x584*M[57] + x596*M[60] + x600*M[32] + x608*M[59] + x611*M[77] + x622*M[18] + x629*M[66] + x631*M[28] + x634*M[33] + x637*M[48] + x640*M[38] + x641*M[45] + x643*M[50] + x644*M[41] + x646*M[51] + x647*M[42] + x651*M[79] + x654*M[68] + x656*M[62] + x657*M[71] + x659*M[72] + x660*M[63] + x662*M[67] + x663*M[78] + x674*M[81] + x691*M[43] + x696*M[64] + x698*M[69] + x702*M[47] + x703*M[52] + x704*M[53] + x705*M[73] + x706*M[74] + x707*M[80] + x713*M[54] + x718*M[75] + x724*M[82] + x95*M[3];
#pragma omp atomic
L[25] += x348*M[3] + x353*M[0] + x357*M[1] + x359*M[2] + x362*M[4] + x366*M[5] + x368*M[7] + x371*M[13] + x381*M[6] + x431*M[19] + x435*M[8] + x474*M[9] + x477*M[10] + x483*M[11] + x501*M[20] + x507*M[21] + x516*M[23] + x522*M[38] + x527*M[55] + x533*M[29] + x535*M[22] + x539*M[24] + x542*M[12] + x543*M[15] + x544*M[16] + x561*M[25] + x562*M[26] + x563*M[30] + x564*M[45] + x567*M[34] + x575*M[35] + x578*M[36] + x584*M[56] + x593*M[57] + x596*M[59] + x600*M[31] + x604*M[33] + x608*M[58] + x611*M[76] + x620*M[60] + x621*M[14] + x622*M[17] + x624*M[18] + x629*M[65] + x631*M[27] + x632*M[28] + x634*M[32] + x637*M[47] + x640*M[37] + x641*M[44] + x643*M[49] + x644*M[40] + x646*M[50] + x647*M[41] + x651*M[78] + x654*M[67] + x656*M[61] + x657*M[70] + x659*M[71] + x660*M[62] + x662*M[66] + x663*M[77] + x668*M[69] + x674*M[80] + x684*M[82] + x689*M[39] + x691*M[42] + x692*M[43] + x696*M[63] + x697*M[64] + x698*M[68] + x702*M[46] + x703*M[51] + x704*M[52] + x705*M[72] + x706*M[73] + x707*M[79] + x709*M[48] + x713*M[53] + x715*M[54] + x718*M[74] + x722*M[75] + x724*M[81];
#pragma omp atomic
L[26] += x101*M[6] + x103*M[7] - x116*M[16] - x142*M[15] + x172*M[29] + x176*M[30] + x192*M[0] - x199*M[50] + x210*M[3] + x212*M[5] + x218*M[11] - x248*M[49] + x269*M[76] + x272*M[77] + x284*M[4] + x292*M[10] + x302*M[9] + x305*M[12] + x307*M[13] + x318*M[22] + x325*M[19] + x329*M[23] + x330*M[26] + x332*M[21] + x339*M[36] + x343*M[41] + x363*M[8] + x372*M[18] + x379*M[25] + x392*M[40] + x393*M[44] + x398*M[34] + x399*M[45] + x409*M[55] + x414*M[65] + x418*M[66] + x421*M[57] + x422*M[71] + x444*M[20] + x450*M[35] + x453*M[37] + x454*M[38] + x459*M[58] + x462*M[59] + x463*M[62] + x471*M[17] + x494*M[70] + x504*M[31] + x510*M[32] + x519*M[52] + x533*M[27] + x543*M[14] + x553*M[61] + x560*M[56] + x561*M[24] + x563*M[28] + x564*M[43] + x569*M[51] + x587*M[78] + x590*M[79] - x61*M[2] + x611*M[72] + x629*M[63] + x635*M[33] + x638*M[54] + x641*M[42] + x643*M[46] + x644*M[39] + x646*M[47] + x651*M[74] + x656*M[60] + x657*M[67] + x659*M[68] + x662*M[64] + x663*M[73] + x695*M[53] + x700*M[80] + x701*M[81] + x703*M[48] + x705*M[69] + x707*M[75] + x719*M[82] - x84*M[1];
#pragma omp atomic
L[27] += x103*M[6] - x116*M[15] + x176*M[29] + x193*M[0] - x199*M[49] + x211*M[3] + x212*M[4] + x218*M[10] + x272*M[76] + x306*M[9] + x307*M[12] + x326*M[19] + x329*M[22] + x330*M[25] + x332*M[20] + x339*M[35] + x343*M[40] + x354*M[2] + x363*M[7] + x365*M[8] + x372*M[17] + x381*M[5] + x399*M[44] + x400*M[34] + x418*M[65] + x419*M[55] + x421*M[56] + x422*M[70] + x454*M[37] + x462*M[58] + x463*M[61] + x471*M[16] + x504*M[30] + x510*M[31] + x519*M[51] + x533*M[26] + x535*M[21] + x542*M[11] + x543*M[13] + x544*M[14] + x561*M[23] + x562*M[24] + x563*M[27] + x564*M[42] + x569*M[50] + x587*M[77] + x590*M[78] + x600*M[28] + x608*M[57] - x61*M[1] + x611*M[71] + x629*M[62] + x630*M[18] + x635*M[32] + x636*M[33] + x638*M[53] + x640*M[36] + x641*M[41] + x643*M[45] + x644*M[38] + x646*M[46] + x647*M[39] + x651*M[73] + x654*M[64] + x656*M[59] + x657*M[66] + x659*M[67] + x660*M[60] + x662*M[63] + x663*M[72] + x674*M[75] + x695*M[52] + x700*M[79] + x701*M[80] + x702*M[43] + x703*M[47] + x704*M[48] + x705*M[68] + x706*M[69] + x707*M[74] + x710*M[54] + x719*M[81] + x721*M[82];
#pragma omp atomic
L[28] += x230*M[2] + x354*M[1] + x357*M[0] + x362*M[3] + x363*M[6] + x365*M[7] + x368*M[5] + x371*M[11] + x372*M[16] + x381*M[4] + x471*M[15] + x477*M[9] + x487*M[8] + x501*M[19] + x504*M[29] + x510*M[30] + x516*M[21] + x517*M[18] + x519*M[50] + x522*M[36] + x533*M[25] + x535*M[20] + x542*M[10] + x543*M[12] + x544*M[13] + x561*M[22] + x562*M[23] + x563*M[26] + x564*M[41] + x569*M[49] + x575*M[34] + x584*M[55] + x587*M[76] + x590*M[77] + x596*M[57] + x600*M[27] + x608*M[56] + x611*M[70] + x622*M[14] + x629*M[61] + x630*M[17] + x631*M[24] + x634*M[28] + x635*M[31] + x636*M[32] + x637*M[43] + x638*M[52] + x640*M[35] + x641*M[40] + x643*M[44] + x644*M[37] + x646*M[45] + x647*M[38] + x651*M[72] + x654*M[63] + x656*M[58] + x657*M[65] + x659*M[66] + x660*M[59] + x662*M[62] + x663*M[71] + x674*M[74] + x675*M[33] + x678*M[54] + x691*M[39] + x695*M[51] + x696*M[60] + x698*M[64] + x700*M[78] + x701*M[79] + x702*M[42] + x703*M[46] + x704*M[47] + x705*M[67] + x706*M[68] + x707*M[73] + x710*M[53] + x713*M[48] + x718*M[69] + x719*M[80] + x721*M[81] + x724*M[75] + x730*M[82];
#pragma omp atomic
L[29] += x230*M[1] + x359*M[0] + x365*M[6] + x366*M[3] + x368*M[4] + x371*M[10] + x372*M[15] + x435*M[5] + x466*M[2] + x483*M[9] + x487*M[7] + x488*M[8] + x507*M[19] + x510*M[29] + x516*M[20] + x517*M[17] + x519*M[49] + x522*M[35] + x539*M[21] + x544*M[12] + x562*M[22] + x563*M[25] + x564*M[40] + x578*M[34] + x590*M[76] + x593*M[55] + x596*M[56] + x600*M[26] + x604*M[28] + x620*M[57] + x621*M[11] + x622*M[13] + x624*M[14] + x630*M[16] + x631*M[23] + x632*M[24] + x634*M[27] + x635*M[30] + x636*M[31] + x637*M[42] + x638*M[51] + x646*M[44] + x647*M[37] + x651*M[71] + x654*M[62] + x659*M[65] + x660*M[58] + x662*M[61] + x663*M[70] + x668*M[64] + x671*M[18] + x674*M[73] + x675*M[32] + x676*M[33] + x678*M[53] + x684*M[75] + x689*M[36] + x691*M[38] + x692*M[39] + x695*M[50] + x696*M[59] + x697*M[60] + x698*M[63] + x700*M[77] + x701*M[78] + x702*M[41] + x703*M[45] + x704*M[46] + x705*M[66] + x706*M[67] + x707*M[72] + x709*M[43] + x710*M[52] + x713*M[47] + x715*M[48] + x718*M[68] + x719*M[79] + x721*M[80] + x722*M[69] + x724*M[74] + x728*M[54] + x730*M[81] + x731*M[82];
#pragma omp atomic
L[30] += x100*M[7] + x101*M[4] + x103*M[5] - x116*M[13] + x123*M[6] + x136*M[15] + x141*M[16] - x142*M[12] + x163*M[30] + x172*M[25] + x176*M[26] - x199*M[45] + x229*M[29] + x237*M[49] + x247*M[50] - x248*M[44] + x264*M[77] + x269*M[70] + x272*M[71] + x277*M[76] + x284*M[3] + x292*M[9] + x305*M[10] + x307*M[11] + x318*M[20] + x329*M[21] + x330*M[23] + x343*M[38] + x351*M[8] + x379*M[22] + x392*M[37] + x393*M[40] + x399*M[41] + x414*M[61] + x418*M[62] + x422*M[66] + x428*M[31] + x444*M[19] + x450*M[34] + x453*M[35] + x454*M[36] + x459*M[56] + x462*M[57] + x463*M[59] + x471*M[14] + x480*M[17] + x486*M[18] + x494*M[65] + x504*M[27] + x510*M[28] + x513*M[32] + x519*M[47] + x530*M[78] + x533*M[24] + x541*M[33] + x553*M[58] + x560*M[55] + x569*M[46] + x572*M[51] + x581*M[52] + x587*M[72] + x590*M[73] + x599*M[79] + x611*M[67] + x618*M[80] + x629*M[60] + x641*M[39] + x643*M[42] + x646*M[43] + x651*M[69] + x657*M[63] + x659*M[64] + x663*M[68] + x670*M[82] + x690*M[53] + x693*M[54] + x695*M[48] + x699*M[81] + x700*M[74] + x701*M[75] + x77*M[1] + x83*M[2] - x84*M[0];
#pragma omp atomic
L[31] += x100*M[6] + x103*M[4] - x116*M[12] + x141*M[15] + x163*M[29] + x176*M[25] - x199*M[44] + x212*M[3] + x218*M[9] + x247*M[49] + x264*M[76] + x272*M[70] + x307*M[10] + x329*M[20] + x330*M[22] + x332*M[19] + x339*M[34] + x343*M[37] + x351*M[7] + x356*M[2] + x363*M[5] + x369*M[8] + x372*M[14] + x399*M[40] + x418*M[61] + x421*M[55] + x422*M[65] + x428*M[30] + x454*M[35] + x462*M[56] + x463*M[58] + x471*M[13] + x480*M[16] + x486*M[17] + x504*M[26] + x510*M[27] + x513*M[31] + x519*M[46] + x530*M[77] + x533*M[23] + x541*M[32] + x543*M[11] + x561*M[21] + x563*M[24] + x564*M[39] + x569*M[45] + x572*M[50] + x581*M[51] + x587*M[71] + x590*M[72] + x599*M[78] - x61*M[0] + x611*M[66] + x618*M[79] + x623*M[18] + x629*M[59] + x633*M[33] + x635*M[28] + x638*M[48] + x641*M[38] + x643*M[41] + x644*M[36] + x646*M[42] + x651*M[68] + x656*M[57] + x657*M[62] + x659*M[63] + x662*M[60] + x663*M[67] + x670*M[81] + x690*M[52] + x693*M[53] + x695*M[47] + x699*M[80] + x700*M[73] + x701*M[74] + x703*M[43] + x705*M[64] + x707*M[69] + x712*M[54] + x719*M[75] + x725*M[82] + x83*M[1];
#pragma omp atomic
L[32] += x351*M[6] + x354*M[0] + x356*M[1] + x360*M[2] + x363*M[4] + x365*M[5] + x369*M[7] + x372*M[13] + x381*M[3] + x428*M[29] + x434*M[8] + x471*M[12] + x480*M[15] + x486*M[16] + x504*M[25] + x510*M[26] + x513*M[30] + x519*M[45] + x530*M[76] + x533*M[22] + x535*M[19] + x541*M[31] + x542*M[9] + x543*M[10] + x544*M[11] + x561*M[20] + x562*M[21] + x563*M[23] + x564*M[38] + x569*M[44] + x572*M[49] + x581*M[50] + x587*M[70] + x590*M[71] + x599*M[77] + x600*M[24] + x603*M[33] + x608*M[55] + x611*M[65] + x618*M[78] + x623*M[17] + x625*M[18] + x629*M[58] + x630*M[14] + x633*M[32] + x635*M[27] + x636*M[28] + x638*M[47] + x640*M[34] + x641*M[37] + x643*M[40] + x644*M[35] + x646*M[41] + x647*M[36] + x651*M[67] + x654*M[60] + x656*M[56] + x657*M[61] + x659*M[62] + x660*M[57] + x662*M[59] + x663*M[66] + x670*M[80] + x674*M[69] + x683*M[82] + x690*M[51] + x693*M[52] + x695*M[46] + x699*M[79] + x700*M[72] + x701*M[73] + x702*M[39] + x703*M[42] + x704*M[43] + x705*M[63] + x706*M[64] + x707*M[68] + x710*M[48] + x712*M[53] + x716*M[54] + x719*M[74] + x721*M[75] + x725*M[81];
#pragma omp atomic
L[33] += x230*M[0] + x360*M[1] + x365*M[4] + x368*M[3] + x369*M[6] + x371*M[9] + x372*M[12] + x434*M[7] + x464*M[2] + x469*M[8] + x486*M[15] + x487*M[5] + x510*M[25] + x513*M[29] + x516*M[19] + x517*M[14] + x519*M[44] + x522*M[34] + x541*M[30] + x544*M[10] + x562*M[20] + x563*M[22] + x564*M[37] + x581*M[49] + x590*M[70] + x596*M[55] + x599*M[76] + x600*M[23] + x603*M[32] + x618*M[77] + x622*M[11] + x623*M[16] + x625*M[17] + x630*M[13] + x631*M[21] + x633*M[31] + x634*M[24] + x635*M[26] + x636*M[27] + x637*M[39] + x638*M[46] + x646*M[40] + x647*M[35] + x651*M[66] + x654*M[59] + x659*M[61] + x660*M[56] + x662*M[58] + x663*M[65] + x670*M[79] + x672*M[18] + x674*M[68] + x675*M[28] + x677*M[33] + x678*M[48] + x683*M[81] + x690*M[50] + x691*M[36] + x693*M[51] + x695*M[45] + x696*M[57] + x698*M[60] + x699*M[78] + x700*M[71] + x701*M[72] + x702*M[38] + x703*M[41] + x704*M[42] + x705*M[62] + x706*M[63] + x707*M[67] + x710*M[47] + x712*M[52] + x713*M[43] + x716*M[53] + x718*M[64] + x719*M[73] + x721*M[74] + x724*M[69] + x725*M[80] + x726*M[54] + x729*M[82] + x730*M[75];
#pragma omp atomic
L[34] += x434*M[6] + x435*M[3] + x464*M[1] + x465*M[2] + x466*M[0] + x468*M[8] + x469*M[7] + x487*M[4] + x488*M[5] + x517*M[13] + x539*M[19] + x541*M[29] + x600*M[22] + x603*M[31] + x604*M[24] + x613*M[33] + x618*M[76] + x620*M[55] + x621*M[9] + x622*M[10] + x623*M[15] + x624*M[11] + x625*M[16] + x630*M[12] + x631*M[20] + x632*M[21] + x633*M[30] + x634*M[23] + x635*M[25] + x636*M[26] + x637*M[38] + x638*M[45] + x651*M[65] + x654*M[58] + x668*M[60] + x670*M[78] + x671*M[14] + x672*M[17] + x673*M[18] + x674*M[67] + x675*M[27] + x676*M[28] + x677*M[32] + x678*M[47] + x683*M[80] + x684*M[69] + x688*M[82] + x689*M[34] + x690*M[49] + x691*M[35] + x692*M[36] + x693*M[50] + x695*M[44] + x696*M[56] + x697*M[57] + x698*M[59] + x699*M[77] + x700*M[70] + x701*M[71] + x702*M[37] + x703*M[40] + x704*M[41] + x705*M[61] + x706*M[62] + x707*M[66] + x709*M[39] + x710*M[46] + x712*M[51] + x713*M[42] + x715*M[43] + x716*M[52] + x718*M[63] + x719*M[72] + x721*M[73] + x722*M[64] + x724*M[68] + x725*M[79] + x726*M[53] + x727*M[54] + x728*M[48] + x729*M[81] + x730*M[74] + x731*M[75];
#pragma omp atomic
L[35] += x112*M[7] + x120*M[0] + x129*M[3] + x132*M[4] + x137*M[5] + x153*M[10] + x154*M[11] + x169*M[13] + x187*M[23] + x226*M[9] + x234*M[19] + x241*M[20] + x242*M[21] + x255*M[35] + x256*M[36] + x260*M[38] + x275*M[34] + x297*M[6] + x325*M[15] + x326*M[16] + x339*M[30] + x376*M[12] + x384*M[22] + x398*M[25] + x400*M[26] + x409*M[40] + x419*M[41] + x421*M[45] + x431*M[14] + x450*M[29] + x459*M[49] + x462*M[50] + x474*M[8] + x498*M[37] + x501*M[17] + x507*M[18] + x522*M[32] + x527*M[39] + x560*M[44] + x567*M[24] + x575*M[27] + x578*M[28] + x584*M[42] + x593*M[43] + x596*M[47] + x608*M[46] + x620*M[48] + x640*M[31] + x656*M[51] + x660*M[52] + x689*M[33] + x696*M[53] + x697*M[54] + x91*M[1] + x92*M[2];
#pragma omp atomic
L[36] += x112*M[5] + x132*M[3] + x153*M[9] + x169*M[11] + x187*M[21] + x241*M[19] + x255*M[34] + x260*M[36] + x286*M[1] + x297*M[4] + x302*M[6] + x306*M[7] + x325*M[12] + x326*M[13] + x332*M[16] + x339*M[26] + x376*M[10] + x384*M[20] + x398*M[22] + x400*M[23] + x409*M[37] + x419*M[38] + x421*M[41] + x444*M[15] + x450*M[25] + x453*M[29] + x454*M[30] + x459*M[44] + x462*M[45] + x463*M[50] + x477*M[8] + x498*M[35] + x501*M[14] + x516*M[18] + x522*M[28] + x535*M[17] + x553*M[49] + x560*M[40] + x575*M[24] + x584*M[39] + x596*M[43] + x608*M[42] + x629*M[51] + x640*M[27] + x644*M[31] + x647*M[32] + x654*M[53] + x656*M[46] + x660*M[47] + x662*M[52] + x691*M[33] + x696*M[48] + x698*M[54] + x91*M[0] + x95*M[2];
#pragma omp atomic
L[37] += x112*M[4] + x137*M[3] + x154*M[9] + x169*M[10] + x187*M[20] + x242*M[19] + x256*M[34] + x260*M[35] + x306*M[6] + x326*M[12] + x332*M[15] + x339*M[25] + x348*M[2] + x400*M[22] + x419*M[37] + x421*M[40] + x431*M[11] + x454*M[29] + x462*M[44] + x463*M[49] + x474*M[5] + x477*M[7] + x483*M[8] + x501*M[13] + x507*M[14] + x516*M[17] + x522*M[27] + x527*M[36] + x535*M[16] + x539*M[18] + x567*M[21] + x575*M[23] + x578*M[24] + x584*M[38] + x593*M[39] + x596*M[42] + x608*M[41] + x620*M[43] + x629*M[50] + x640*M[26] + x644*M[30] + x647*M[31] + x654*M[52] + x656*M[45] + x660*M[46] + x662*M[51] + x668*M[54] + x689*M[28] + x691*M[32] + x692*M[33] + x696*M[47] + x697*M[48] + x698*M[53] + x92*M[0] + x95*M[1];
#pragma omp atomic
L[38] += x210*M[1] + x211*M[2] + x218*M[7] + x286*M[0] + x292*M[6] + x297*M[3] + x302*M[4] + x306*M[5] + x318*M[15] + x325*M[10] + x326*M[11] + x329*M[16] + x332*M[13] + x339*M[23] + x343*M[30] + x376*M[9] + x384*M[19] + x392*M[29] + x398*M[20] + x400*M[21] + x409*M[35] + x414*M[49] + x418*M[50] + x419*M[36] + x421*M[38] + x444*M[12] + x450*M[22] + x453*M[25] + x454*M[26] + x459*M[40] + x462*M[41] + x463*M[45] + x498*M[34] + x535*M[14] + x542*M[8] + x553*M[44] + x560*M[37] + x561*M[17] + x562*M[18] + x564*M[32] + x608*M[39] + x629*M[46] + x640*M[24] + x641*M[31] + x644*M[27] + x647*M[28] + x654*M[48] + x656*M[42] + x657*M[51] + x659*M[52] + x660*M[43] + x662*M[47] + x702*M[33] + x705*M[53] + x706*M[54];
#pragma omp atomic
L[39] += x112*M[3] + x169*M[9] + x187*M[19] + x211*M[1] + x218*M[6] + x260*M[34] + x306*M[4] + x326*M[10] + x329*M[15] + x332*M[12] + x339*M[22] + x343*M[29] + x362*M[2] + x371*M[8] + x400*M[20] + x418*M[49] + x419*M[35] + x421*M[37] + x454*M[25] + x462*M[40] + x463*M[44] + x477*M[5] + x501*M[11] + x516*M[14] + x522*M[24] + x535*M[13] + x542*M[7] + x561*M[16] + x562*M[17] + x564*M[31] + x575*M[21] + x584*M[36] + x596*M[39] + x608*M[38] + x629*M[45] + x631*M[18] + x637*M[33] + x640*M[23] + x641*M[30] + x644*M[26] + x647*M[27] + x654*M[47] + x656*M[41] + x657*M[50] + x659*M[51] + x660*M[42] + x662*M[46] + x691*M[28] + x696*M[43] + x698*M[48] + x702*M[32] + x705*M[52] + x706*M[53] + x718*M[54] + x95*M[0];
#pragma omp atomic
L[40] += x348*M[0] + x362*M[1] + x366*M[2] + x371*M[7] + x431*M[9] + x474*M[3] + x477*M[4] + x483*M[5] + x501*M[10] + x507*M[11] + x516*M[13] + x522*M[23] + x527*M[34] + x535*M[12] + x539*M[14] + x542*M[6] + x561*M[15] + x562*M[16] + x564*M[30] + x567*M[19] + x575*M[20] + x578*M[21] + x584*M[35] + x593*M[36] + x596*M[38] + x608*M[37] + x620*M[39] + x621*M[8] + x629*M[44] + x631*M[17] + x632*M[18] + x637*M[32] + x640*M[22] + x641*M[29] + x644*M[25] + x647*M[26] + x654*M[46] + x656*M[40] + x657*M[49] + x659*M[50] + x660*M[41] + x662*M[45] + x668*M[48] + x689*M[24] + x691*M[27] + x692*M[28] + x696*M[42] + x697*M[43] + x698*M[47] + x702*M[31] + x705*M[51] + x706*M[52] + x709*M[33] + x718*M[53] + x722*M[54];
#pragma omp atomic
L[41] += x210*M[0] + x212*M[2] + x218*M[5] + x284*M[1] + x292*M[4] + x302*M[3] + x305*M[6] + x307*M[7] + x318*M[12] + x325*M[9] + x329*M[13] + x330*M[16] + x332*M[11] + x339*M[21] + x343*M[26] + x379*M[15] + x392*M[25] + x393*M[29] + x398*M[19] + x399*M[30] + x409*M[34] + x414*M[44] + x418*M[45] + x421*M[36] + x422*M[50] + x444*M[10] + x450*M[20] + x453*M[22] + x454*M[23] + x459*M[37] + x462*M[38] + x463*M[41] + x494*M[49] + x533*M[17] + x543*M[8] + x553*M[40] + x560*M[35] + x561*M[14] + x563*M[18] + x564*M[28] + x611*M[51] + x629*M[42] + x641*M[27] + x643*M[31] + x644*M[24] + x646*M[32] + x651*M[53] + x656*M[39] + x657*M[46] + x659*M[47] + x662*M[43] + x663*M[52] + x703*M[33] + x705*M[48] + x707*M[54];
#pragma omp atomic
L[42] += x211*M[0] + x212*M[1] + x218*M[4] + x306*M[3] + x307*M[6] + x326*M[9] + x329*M[12] + x330*M[15] + x332*M[10] + x339*M[20] + x343*M[25] + x381*M[2] + x399*M[29] + x400*M[19] + x418*M[44] + x419*M[34] + x421*M[35] + x422*M[49] + x454*M[22] + x462*M[37] + x463*M[40] + x533*M[16] + x535*M[11] + x542*M[5] + x543*M[7] + x544*M[8] + x561*M[13] + x562*M[14] + x563*M[17] + x564*M[27] + x600*M[18] + x608*M[36] + x611*M[50] + x629*M[41] + x640*M[21] + x641*M[26] + x643*M[30] + x644*M[23] + x646*M[31] + x647*M[24] + x651*M[52] + x654*M[43] + x656*M[38] + x657*M[45] + x659*M[46] + x660*M[39] + x662*M[42] + x663*M[51] + x674*M[54] + x702*M[28] + x703*M[32] + x704*M[33] + x705*M[47] + x706*M[48] + x707*M[53];
#pragma omp atomic
L[43] += x362*M[0] + x368*M[2] + x371*M[5] + x381*M[1] + x477*M[3] + x501*M[9] + x516*M[11] + x522*M[21] + x533*M[15] + x535*M[10] + x542*M[4] + x543*M[6] + x544*M[7] + x561*M[12] + x562*M[13] + x563*M[16] + x564*M[26] + x575*M[19] + x584*M[34] + x596*M[36] + x600*M[17] + x608*M[35] + x611*M[49] + x622*M[8] + x629*M[40] + x631*M[14] + x634*M[18] + x637*M[28] + x640*M[20] + x641*M[25] + x643*M[29] + x644*M[22] + x646*M[30] + x647*M[23] + x651*M[51] + x654*M[42] + x656*M[37] + x657*M[44] + x659*M[45] + x660*M[38] + x662*M[41] + x663*M[50] + x674*M[53] + x691*M[24] + x696*M[39] + x698*M[43] + x702*M[27] + x703*M[31] + x704*M[32] + x705*M[46] + x706*M[47] + x707*M[52] + x713*M[33] + x718*M[48] + x724*M[54];
#pragma omp atomic
L[44] += x366*M[0] + x368*M[1] + x371*M[4] + x435*M[2] + x483*M[3] + x507*M[9] + x516*M[10] + x522*M[20] + x539*M[11] + x544*M[6] + x562*M[12] + x563*M[15] + x564*M[25] + x578*M[19] + x593*M[34] + x596*M[35] + x600*M[16] + x604*M[18] + x620*M[36] + x621*M[5] + x622*M[7] + x624*M[8] + x631*M[13] + x632*M[14] + x634*M[17] + x637*M[27] + x646*M[29] + x647*M[22] + x651*M[50] + x654*M[41] + x659*M[44] + x660*M[37] + x662*M[40] + x663*M[49] + x668*M[43] + x674*M[52] + x684*M[54] + x689*M[21] + x691*M[23] + x692*M[24] + x696*M[38] + x697*M[39] + x698*M[42] + x702*M[26] + x703*M[30] + x704*M[31] + x705*M[45] + x706*M[46] + x707*M[51] + x709*M[28] + x713*M[32] + x715*M[33] + x718*M[47] + x722*M[48] + x724*M[53];
#pragma omp atomic
L[45] += x101*M[1] + x103*M[2] - x116*M[7] - x142*M[6] + x172*M[15] + x176*M[16] - x199*M[30] - x248*M[29] + x269*M[49] + x272*M[50] + x284*M[0] + x292*M[3] + x305*M[4] + x307*M[5] + x318*M[10] + x329*M[11] + x330*M[13] + x343*M[23] + x379*M[12] + x392*M[22] + x393*M[25] + x399*M[26] + x414*M[40] + x418*M[41] + x422*M[45] + x444*M[9] + x450*M[19] + x453*M[20] + x454*M[21] + x459*M[35] + x462*M[36] + x463*M[38] + x471*M[8] + x494*M[44] + x504*M[17] + x510*M[18] + x519*M[32] + x533*M[14] + x553*M[37] + x560*M[34] + x569*M[31] + x587*M[51] + x590*M[52] + x611*M[46] + x629*M[39] + x641*M[24] + x643*M[27] + x646*M[28] + x651*M[48] + x657*M[42] + x659*M[43] + x663*M[47] + x695*M[33] + x700*M[53] + x701*M[54];
#pragma omp atomic
L[46] += x103*M[1] - x116*M[6] + x176*M[15] - x199*M[29] + x212*M[0] + x218*M[3] + x272*M[49] + x307*M[4] + x329*M[10] + x330*M[12] + x332*M[9] + x339*M[19] + x343*M[22] + x363*M[2] + x372*M[8] + x399*M[25] + x418*M[40] + x421*M[34] + x422*M[44] + x454*M[20] + x462*M[35] + x463*M[37] + x471*M[7] + x504*M[16] + x510*M[17] + x519*M[31] + x533*M[13] + x543*M[5] + x561*M[11] + x563*M[14] + x564*M[24] + x569*M[30] + x587*M[50] + x590*M[51] + x611*M[45] + x629*M[38] + x635*M[18] + x638*M[33] + x641*M[23] + x643*M[26] + x644*M[21] + x646*M[27] + x651*M[47] + x656*M[36] + x657*M[41] + x659*M[42] + x662*M[39] + x663*M[46] + x695*M[32] + x700*M[52] + x701*M[53] + x703*M[28] + x705*M[43] + x707*M[48] + x719*M[54];
#pragma omp atomic
L[47] += x363*M[1] + x365*M[2] + x372*M[7] + x381*M[0] + x471*M[6] + x504*M[15] + x510*M[16] + x519*M[30] + x533*M[12] + x535*M[9] + x542*M[3] + x543*M[4] + x544*M[5] + x561*M[10] + x562*M[11] + x563*M[13] + x564*M[23] + x569*M[29] + x587*M[49] + x590*M[50] + x600*M[14] + x608*M[34] + x611*M[44] + x629*M[37] + x630*M[8] + x635*M[17] + x636*M[18] + x638*M[32] + x640*M[19] + x641*M[22] + x643*M[25] + x644*M[20] + x646*M[26] + x647*M[21] + x651*M[46] + x654*M[39] + x656*M[35] + x657*M[40] + x659*M[41] + x660*M[36] + x662*M[38] + x663*M[45] + x674*M[48] + x695*M[31] + x700*M[51] + x701*M[52] + x702*M[24] + x703*M[27] + x704*M[28] + x705*M[42] + x706*M[43] + x707*M[47] + x710*M[33] + x719*M[53] + x721*M[54];
#pragma omp atomic
L[48] += x365*M[1] + x368*M[0] + x371*M[3] + x372*M[6] + x487*M[2] + x510*M[15] + x516*M[9] + x517*M[8] + x519*M[29] + x522*M[19] + x544*M[4] + x562*M[10] + x563*M[12] + x564*M[22] + x590*M[49] + x596*M[34] + x600*M[13] + x622*M[5] + x630*M[7] + x631*M[11] + x634*M[14] + x635*M[16] + x636*M[17] + x637*M[24] + x638*M[31] + x646*M[25] + x647*M[20] + x651*M[45] + x654*M[38] + x659*M[40] + x660*M[35] + x662*M[37] + x663*M[44] + x674*M[47] + x675*M[18] + x678*M[33] + x691*M[21] + x695*M[30] + x696*M[36] + x698*M[39] + x700*M[50] + x701*M[51] + x702*M[23] + x703*M[26] + x704*M[27] + x705*M[41] + x706*M[42] + x707*M[46] + x710*M[32] + x713*M[28] + x718*M[43] + x719*M[52] + x721*M[53] + x724*M[48] + x730*M[54];
#pragma omp atomic
L[49] += x435*M[0] + x487*M[1] + x488*M[2] + x517*M[7] + x539*M[9] + x600*M[12] + x604*M[14] + x620*M[34] + x621*M[3] + x622*M[4] + x624*M[5] + x630*M[6] + x631*M[10] + x632*M[11] + x634*M[13] + x635*M[15] + x636*M[16] + x637*M[23] + x638*M[30] + x651*M[44] + x654*M[37] + x668*M[39] + x671*M[8] + x674*M[46] + x675*M[17] + x676*M[18] + x678*M[32] + x684*M[48] + x689*M[19] + x691*M[20] + x692*M[21] + x695*M[29] + x696*M[35] + x697*M[36] + x698*M[38] + x700*M[49] + x701*M[50] + x702*M[22] + x703*M[25] + x704*M[26] + x705*M[40] + x706*M[41] + x707*M[45] + x709*M[24] + x710*M[31] + x713*M[27] + x715*M[28] + x718*M[42] + x719*M[51] + x721*M[52] + x722*M[43] + x724*M[47] + x728*M[33] + x730*M[53] + x731*M[54];
#pragma omp atomic
L[50] += x100*M[2] + x101*M[0] - x116*M[5] + x123*M[1] + x136*M[6] + x141*M[7] - x142*M[4] + x163*M[16] + x172*M[12] + x176*M[13] - x199*M[26] + x229*M[15] + x237*M[29] + x247*M[30] - x248*M[25] + x264*M[50] + x269*M[44] + x272*M[45] + x277*M[49] + x305*M[3] + x318*M[9] + x330*M[11] + x343*M[21] + x379*M[10] + x392*M[20] + x393*M[22] + x399*M[23] + x414*M[37] + x418*M[38] + x422*M[41] + x428*M[17] + x453*M[19] + x459*M[34] + x463*M[36] + x480*M[8] + x494*M[40] + x504*M[14] + x513*M[18] + x519*M[28] + x530*M[51] + x553*M[35] + x569*M[27] + x572*M[31] + x581*M[32] + x587*M[46] + x590*M[47] + x599*M[52] + x611*M[42] + x618*M[53] + x643*M[24] + x657*M[39] + x663*M[43] + x690*M[33] + x699*M[54] + x700*M[48];
#pragma omp atomic
L[51] += x100*M[1] + x103*M[0] - x116*M[4] + x141*M[6] + x163*M[15] + x176*M[12] - x199*M[25] + x247*M[29] + x264*M[49] + x272*M[44] + x307*M[3] + x329*M[9] + x330*M[10] + x343*M[20] + x351*M[2] + x399*M[22] + x418*M[37] + x422*M[40] + x428*M[16] + x454*M[19] + x462*M[34] + x463*M[35] + x471*M[5] + x480*M[7] + x486*M[8] + x504*M[13] + x510*M[14] + x513*M[17] + x519*M[27] + x530*M[50] + x533*M[11] + x541*M[18] + x569*M[26] + x572*M[30] + x581*M[31] + x587*M[45] + x590*M[46] + x599*M[51] + x611*M[41] + x618*M[52] + x629*M[36] + x641*M[21] + x643*M[23] + x646*M[24] + x651*M[43] + x657*M[38] + x659*M[39] + x663*M[42] + x670*M[54] + x690*M[32] + x693*M[33] + x695*M[28] + x699*M[53] + x700*M[47] + x701*M[48];
#pragma omp atomic
L[52] += x351*M[1] + x363*M[0] + x369*M[2] + x372*M[5] + x428*M[15] + x471*M[4] + x480*M[6] + x486*M[7] + x504*M[12] + x510*M[13] + x513*M[16] + x519*M[26] + x530*M[49] + x533*M[10] + x541*M[17] + x543*M[3] + x561*M[9] + x563*M[11] + x564*M[21] + x569*M[25] + x572*M[29] + x581*M[30] + x587*M[44] + x590*M[45] + x599*M[50] + x611*M[40] + x618*M[51] + x623*M[8] + x629*M[35] + x633*M[18] + x635*M[14] + x638*M[28] + x641*M[20] + x643*M[22] + x644*M[19] + x646*M[23] + x651*M[42] + x656*M[34] + x657*M[37] + x659*M[38] + x662*M[36] + x663*M[41] + x670*M[53] + x690*M[31] + x693*M[32] + x695*M[27] + x699*M[52] + x700*M[46] + x701*M[47] + x703*M[24] + x705*M[39] + x707*M[43] + x712*M[33] + x719*M[48] + x725*M[54];
#pragma omp atomic
L[53] += x365*M[0] + x369*M[1] + x372*M[4] + x434*M[2] + x486*M[6] + x510*M[12] + x513*M[15] + x519*M[25] + x541*M[16] + x544*M[3] + x562*M[9] + x563*M[10] + x564*M[20] + x581*M[29] + x590*M[44] + x599*M[49] + x600*M[11] + x603*M[18] + x618*M[50] + x623*M[7] + x625*M[8] + x630*M[5] + x633*M[17] + x635*M[13] + x636*M[14] + x638*M[27] + x646*M[22] + x647*M[19] + x651*M[41] + x654*M[36] + x659*M[37] + x660*M[34] + x662*M[35] + x663*M[40] + x670*M[52] + x674*M[43] + x683*M[54] + x690*M[30] + x693*M[31] + x695*M[26] + x699*M[51] + x700*M[45] + x701*M[46] + x702*M[21] + x703*M[23] + x704*M[24] + x705*M[38] + x706*M[39] + x707*M[42] + x710*M[28] + x712*M[32] + x716*M[33] + x719*M[47] + x721*M[48] + x725*M[53];
#pragma omp atomic
L[54] += x434*M[1] + x469*M[2] + x487*M[0] + x517*M[5] + x541*M[15] + x600*M[10] + x603*M[17] + x618*M[49] + x622*M[3] + x623*M[6] + x625*M[7] + x630*M[4] + x631*M[9] + x633*M[16] + x634*M[11] + x635*M[12] + x636*M[13] + x637*M[21] + x638*M[26] + x651*M[40] + x654*M[35] + x670*M[51] + x672*M[8] + x674*M[42] + x675*M[14] + x677*M[18] + x678*M[28] + x683*M[53] + x690*M[29] + x691*M[19] + x693*M[30] + x695*M[25] + x696*M[34] + x698*M[36] + x699*M[50] + x700*M[44] + x701*M[45] + x702*M[20] + x703*M[22] + x704*M[23] + x705*M[37] + x706*M[38] + x707*M[41] + x710*M[27] + x712*M[31] + x713*M[24] + x716*M[32] + x718*M[39] + x719*M[46] + x721*M[47] + x724*M[43] + x725*M[52] + x726*M[33] + x729*M[54] + x730*M[48];
#pragma omp atomic
L[55] += x468*M[2] + x469*M[1] + x488*M[0] + x517*M[4] + x603*M[16] + x604*M[11] + x613*M[18] + x624*M[3] + x625*M[6] + x632*M[9] + x633*M[15] + x634*M[10] + x636*M[12] + x637*M[20] + x638*M[25] + x668*M[36] + x670*M[50] + x671*M[5] + x672*M[7] + x673*M[8] + x674*M[41] + x675*M[13] + x676*M[14] + x677*M[17] + x678*M[27] + x683*M[52] + x684*M[43] + x688*M[54] + x692*M[19] + x693*M[29] + x697*M[34] + x698*M[35] + x699*M[49] + x701*M[44] + x704*M[22] + x706*M[37] + x707*M[40] + x709*M[21] + x710*M[26] + x712*M[30] + x713*M[23] + x715*M[24] + x716*M[31] + x718*M[38] + x719*M[45] + x721*M[46] + x722*M[39] + x724*M[42] + x725*M[51] + x726*M[32] + x727*M[33] + x728*M[28] + x729*M[53] + x730*M[47] + x731*M[48];
#pragma omp atomic
L[56] += x129*M[0] + x132*M[1] + x137*M[2] + x153*M[4] + x154*M[5] + x169*M[7] + x187*M[13] + x226*M[3] + x234*M[9] + x241*M[10] + x242*M[11] + x255*M[20] + x256*M[21] + x260*M[23] + x275*M[19] + x376*M[6] + x384*M[12] + x398*M[15] + x400*M[16] + x409*M[25] + x419*M[26] + x421*M[30] + x431*M[8] + x498*M[22] + x527*M[24] + x560*M[29] + x567*M[14] + x575*M[17] + x578*M[18] + x584*M[27] + x593*M[28] + x596*M[32] + x608*M[31] + x620*M[33];
#pragma omp atomic
L[57] += x112*M[2] + x132*M[0] + x153*M[3] + x169*M[5] + x187*M[11] + x241*M[9] + x255*M[19] + x260*M[21] + x297*M[1] + x325*M[6] + x326*M[7] + x339*M[16] + x376*M[4] + x384*M[10] + x398*M[12] + x400*M[13] + x409*M[22] + x419*M[23] + x421*M[26] + x450*M[15] + x459*M[29] + x462*M[30] + x498*M[20] + x501*M[8] + x522*M[18] + x560*M[25] + x575*M[14] + x584*M[24] + x596*M[28] + x608*M[27] + x640*M[17] + x656*M[31] + x660*M[32] + x696*M[33];
#pragma omp atomic
L[58] += x112*M[1] + x137*M[0] + x154*M[3] + x169*M[4] + x187*M[10] + x242*M[9] + x256*M[19] + x260*M[20] + x326*M[6] + x339*M[15] + x400*M[12] + x419*M[22] + x421*M[25] + x431*M[5] + x462*M[29] + x474*M[2] + x501*M[7] + x507*M[8] + x522*M[17] + x527*M[21] + x567*M[11] + x575*M[13] + x578*M[14] + x584*M[23] + x593*M[24] + x596*M[27] + x608*M[26] + x620*M[28] + x640*M[16] + x656*M[30] + x660*M[31] + x689*M[18] + x696*M[32] + x697*M[33];
#pragma omp atomic
L[59] += x297*M[0] + x302*M[1] + x306*M[2] + x325*M[4] + x326*M[5] + x332*M[7] + x339*M[13] + x376*M[3] + x384*M[9] + x398*M[10] + x400*M[11] + x409*M[20] + x419*M[21] + x421*M[23] + x444*M[6] + x450*M[12] + x453*M[15] + x454*M[16] + x459*M[25] + x462*M[26] + x463*M[30] + x498*M[19] + x535*M[8] + x553*M[29] + x560*M[22] + x608*M[24] + x629*M[31] + x640*M[14] + x644*M[17] + x647*M[18] + x654*M[33] + x656*M[27] + x660*M[28] + x662*M[32];
#pragma omp atomic
L[60] += x112*M[0] + x169*M[3] + x187*M[9] + x260*M[19] + x306*M[1] + x326*M[4] + x332*M[6] + x339*M[12] + x400*M[10] + x419*M[20] + x421*M[22] + x454*M[15] + x462*M[25] + x463*M[29] + x477*M[2] + x501*M[5] + x516*M[8] + x522*M[14] + x535*M[7] + x575*M[11] + x584*M[21] + x596*M[24] + x608*M[23] + x629*M[30] + x640*M[13] + x644*M[16] + x647*M[17] + x654*M[32] + x656*M[26] + x660*M[27] + x662*M[31] + x691*M[18] + x696*M[28] + x698*M[33];
#pragma omp atomic
L[61] += x431*M[3] + x474*M[0] + x477*M[1] + x483*M[2] + x501*M[4] + x507*M[5] + x516*M[7] + x522*M[13] + x527*M[19] + x535*M[6] + x539*M[8] + x567*M[9] + x575*M[10] + x578*M[11] + x584*M[20] + x593*M[21] + x596*M[23] + x608*M[22] + x620*M[24] + x629*M[29] + x640*M[12] + x644*M[15] + x647*M[16] + x654*M[31] + x656*M[25] + x660*M[26] + x662*M[30] + x668*M[33] + x689*M[14] + x691*M[17] + x692*M[18] + x696*M[27] + x697*M[28] + x698*M[32];
#pragma omp atomic
L[62] += x218*M[2] + x292*M[1] + x302*M[0] + x318*M[6] + x325*M[3] + x329*M[7] + x332*M[5] + x339*M[11] + x343*M[16] + x392*M[15] + x398*M[9] + x409*M[19] + x414*M[29] + x418*M[30] + x421*M[21] + x444*M[4] + x450*M[10] + x453*M[12] + x454*M[13] + x459*M[22] + x462*M[23] + x463*M[26] + x553*M[25] + x560*M[20] + x561*M[8] + x564*M[18] + x629*M[27] + x641*M[17] + x644*M[14] + x656*M[24] + x657*M[31] + x659*M[32] + x662*M[28] + x705*M[33];
#pragma omp atomic
L[63] += x218*M[1] + x306*M[0] + x326*M[3] + x329*M[6] + x332*M[4] + x339*M[10] + x343*M[15] + x400*M[9] + x418*M[29] + x419*M[19] + x421*M[20] + x454*M[12] + x462*M[22] + x463*M[25] + x535*M[5] + x542*M[2] + x561*M[7] + x562*M[8] + x564*M[17] + x608*M[21] + x629*M[26] + x640*M[11] + x641*M[16] + x644*M[13] + x647*M[14] + x654*M[28] + x656*M[23] + x657*M[30] + x659*M[31] + x660*M[24] + x662*M[27] + x702*M[18] + x705*M[32] + x706*M[33];
#pragma omp atomic
L[64] += x371*M[2] + x477*M[0] + x501*M[3] + x516*M[5] + x522*M[11] + x535*M[4] + x542*M[1] + x561*M[6] + x562*M[7] + x564*M[16] + x575*M[9] + x584*M[19] + x596*M[21] + x608*M[20] + x629*M[25] + x631*M[8] + x637*M[18] + x640*M[10] + x641*M[15] + x644*M[12] + x647*M[13] + x654*M[27] + x656*M[22] + x657*M[29] + x659*M[30] + x660*M[23] + x662*M[26] + x691*M[14] + x696*M[24] + x698*M[28] + x702*M[17] + x705*M[31] + x706*M[32] + x718*M[33];
#pragma omp atomic
L[65] += x371*M[1] + x483*M[0] + x507*M[3] + x516*M[4] + x522*M[10] + x539*M[5] + x562*M[6] + x564*M[15] + x578*M[9] + x593*M[19] + x596*M[20] + x620*M[21] + x621*M[2] + x631*M[7] + x632*M[8] + x637*M[17] + x647*M[12] + x654*M[26] + x659*M[29] + x660*M[22] + x662*M[25] + x668*M[28] + x689*M[11] + x691*M[13] + x692*M[14] + x696*M[23] + x697*M[24] + x698*M[27] + x702*M[16] + x705*M[30] + x706*M[31] + x709*M[18] + x718*M[32] + x722*M[33];
#pragma omp atomic
L[66] += x292*M[0] + x305*M[1] + x307*M[2] + x318*M[4] + x329*M[5] + x330*M[7] + x343*M[13] + x379*M[6] + x392*M[12] + x393*M[15] + x399*M[16] + x414*M[25] + x418*M[26] + x422*M[30] + x444*M[3] + x450*M[9] + x453*M[10] + x454*M[11] + x459*M[20] + x462*M[21] + x463*M[23] + x494*M[29] + x533*M[8] + x553*M[22] + x560*M[19] + x611*M[31] + x629*M[24] + x641*M[14] + x643*M[17] + x646*M[18] + x651*M[33] + x657*M[27] + x659*M[28] + x663*M[32];
#pragma omp atomic
L[67] += x218*M[0] + x307*M[1] + x329*M[4] + x330*M[6] + x332*M[3] + x339*M[9] + x343*M[12] + x399*M[15] + x418*M[25] + x421*M[19] + x422*M[29] + x454*M[10] + x462*M[20] + x463*M[22] + x533*M[7] + x543*M[2] + x561*M[5] + x563*M[8] + x564*M[14] + x611*M[30] + x629*M[23] + x641*M[13] + x643*M[16] + x644*M[11] + x646*M[17] + x651*M[32] + x656*M[21] + x657*M[26] + x659*M[27] + x662*M[24] + x663*M[31] + x703*M[18] + x705*M[28] + x707*M[33];
#pragma omp atomic
L[68] += x533*M[6] + x535*M[3] + x542*M[0] + x543*M[1] + x544*M[2] + x561*M[4] + x562*M[5] + x563*M[7] + x564*M[13] + x600*M[8] + x608*M[19] + x611*M[29] + x629*M[22] + x640*M[9] + x641*M[12] + x643*M[15] + x644*M[10] + x646*M[16] + x647*M[11] + x651*M[31] + x654*M[24] + x656*M[20] + x657*M[25] + x659*M[26] + x660*M[21] + x662*M[23] + x663*M[30] + x674*M[33] + x702*M[14] + x703*M[17] + x704*M[18] + x705*M[27] + x706*M[28] + x707*M[32];
#pragma omp atomic
L[69] += x371*M[0] + x516*M[3] + x522*M[9] + x544*M[1] + x562*M[4] + x563*M[6] + x564*M[12] + x596*M[19] + x600*M[7] + x622*M[2] + x631*M[5] + x634*M[8] + x637*M[14] + x646*M[15] + x647*M[10] + x651*M[30] + x654*M[23] + x659*M[25] + x660*M[20] + x662*M[22] + x663*M[29] + x674*M[32] + x691*M[11] + x696*M[21] + x698*M[24] + x702*M[13] + x703*M[16] + x704*M[17] + x705*M[26] + x706*M[27] + x707*M[31] + x713*M[18] + x718*M[28] + x724*M[33];
#pragma omp atomic
L[70] += x539*M[3] + x600*M[6] + x604*M[8] + x620*M[19] + x621*M[0] + x622*M[1] + x624*M[2] + x631*M[4] + x632*M[5] + x634*M[7] + x637*M[13] + x651*M[29] + x654*M[22] + x668*M[24] + x674*M[31] + x684*M[33] + x689*M[9] + x691*M[10] + x692*M[11] + x696*M[20] + x697*M[21] + x698*M[23] + x702*M[12] + x703*M[15] + x704*M[16] + x705*M[25] + x706*M[26] + x707*M[30] + x709*M[14] + x713*M[17] + x715*M[18] + x718*M[27] + x722*M[28] + x724*M[32];
#pragma omp atomic
L[71] += -x116*M[2] - x142*M[1] + x172*M[6] + x176*M[7] - x199*M[16] - x248*M[15] + x269*M[29] + x272*M[30] + x305*M[0] + x318*M[3] + x330*M[5] + x343*M[11] + x379*M[4] + x392*M[10] + x393*M[12] + x399*M[13] + x414*M[22] + x418*M[23] + x422*M[26] + x453*M[9] + x459*M[19] + x463*M[21] + x494*M[25] + x504*M[8] + x519*M[18] + x553*M[20] + x569*M[17] + x587*M[31] + x590*M[32] + x611*M[27] + x643*M[14] + x657*M[24] + x663*M[28] + x700*M[33];
#pragma omp atomic
L[72] += -x116*M[1] + x176*M[6] - x199*M[15] + x272*M[29] + x307*M[0] + x329*M[3] + x330*M[4] + x343*M[10] + x399*M[12] + x418*M[22] + x422*M[25] + x454*M[9] + x462*M[19] + x463*M[20] + x471*M[2] + x504*M[7] + x510*M[8] + x519*M[17] + x533*M[5] + x569*M[16] + x587*M[30] + x590*M[31] + x611*M[26] + x629*M[21] + x641*M[11] + x643*M[13] + x646*M[14] + x651*M[28] + x657*M[23] + x659*M[24] + x663*M[27] + x695*M[18] + x700*M[32] + x701*M[33];
#pragma omp atomic
L[73] += x372*M[2] + x471*M[1] + x504*M[6] + x510*M[7] + x519*M[16] + x533*M[4] + x543*M[0] + x561*M[3] + x563*M[5] + x564*M[11] + x569*M[15] + x587*M[29] + x590*M[30] + x611*M[25] + x629*M[20] + x635*M[8] + x638*M[18] + x641*M[10] + x643*M[12] + x644*M[9] + x646*M[13] + x651*M[27] + x656*M[19] + x657*M[22] + x659*M[23] + x662*M[21] + x663*M[26] + x695*M[17] + x700*M[31] + x701*M[32] + x703*M[14] + x705*M[24] + x707*M[28] + x719*M[33];
#pragma omp atomic
L[74] += x372*M[1] + x510*M[6] + x519*M[15] + x544*M[0] + x562*M[3] + x563*M[4] + x564*M[10] + x590*M[29] + x600*M[5] + x630*M[2] + x635*M[7] + x636*M[8] + x638*M[17] + x646*M[12] + x647*M[9] + x651*M[26] + x654*M[21] + x659*M[22] + x660*M[19] + x662*M[20] + x663*M[25] + x674*M[28] + x695*M[16] + x700*M[30] + x701*M[31] + x702*M[11] + x703*M[13] + x704*M[14] + x705*M[23] + x706*M[24] + x707*M[27] + x710*M[18] + x719*M[32] + x721*M[33];
#pragma omp atomic
L[75] += x517*M[2] + x600*M[4] + x622*M[0] + x630*M[1] + x631*M[3] + x634*M[5] + x635*M[6] + x636*M[7] + x637*M[11] + x638*M[16] + x651*M[25] + x654*M[20] + x674*M[27] + x675*M[8] + x678*M[18] + x691*M[9] + x695*M[15] + x696*M[19] + x698*M[21] + x700*M[29] + x701*M[30] + x702*M[10] + x703*M[12] + x704*M[13] + x705*M[22] + x706*M[23] + x707*M[26] + x710*M[17] + x713*M[14] + x718*M[24] + x719*M[31] + x721*M[32] + x724*M[28] + x730*M[33];
#pragma omp atomic
L[76] += x517*M[1] + x604*M[5] + x624*M[0] + x632*M[3] + x634*M[4] + x636*M[6] + x637*M[10] + x638*M[15] + x668*M[21] + x671*M[2] + x674*M[26] + x675*M[7] + x676*M[8] + x678*M[17] + x684*M[28] + x692*M[9] + x697*M[19] + x698*M[20] + x701*M[29] + x704*M[12] + x706*M[22] + x707*M[25] + x709*M[11] + x710*M[16] + x713*M[13] + x715*M[14] + x718*M[23] + x719*M[30] + x721*M[31] + x722*M[24] + x724*M[27] + x728*M[18] + x730*M[32] + x731*M[33];
#pragma omp atomic
L[77] += x136*M[1] + x141*M[2] - x142*M[0] + x163*M[7] + x172*M[4] + x176*M[5] - x199*M[13] + x229*M[6] + x237*M[15] + x247*M[16] - x248*M[12] + x264*M[30] + x269*M[25] + x272*M[26] + x277*M[29] + x379*M[3] + x392*M[9] + x393*M[10] + x399*M[11] + x414*M[20] + x418*M[21] + x422*M[23] + x428*M[8] + x494*M[22] + x530*M[31] + x553*M[19] + x569*M[14] + x572*M[17] + x581*M[18] + x587*M[27] + x590*M[28] + x599*M[32] + x611*M[24] + x618*M[33];
#pragma omp atomic
L[78] += -x116*M[0] + x141*M[1] + x163*M[6] + x176*M[4] - x199*M[12] + x247*M[15] + x264*M[29] + x272*M[25] + x330*M[3] + x343*M[9] + x399*M[10] + x418*M[20] + x422*M[22] + x428*M[7] + x463*M[19] + x480*M[2] + x504*M[5] + x513*M[8] + x519*M[14] + x530*M[30] + x569*M[13] + x572*M[16] + x581*M[17] + x587*M[26] + x590*M[27] + x599*M[31] + x611*M[23] + x618*M[32] + x643*M[11] + x657*M[21] + x663*M[24] + x690*M[18] + x699*M[33] + x700*M[28];
#pragma omp atomic
L[79] += x428*M[6] + x471*M[0] + x480*M[1] + x486*M[2] + x504*M[4] + x510*M[5] + x513*M[7] + x519*M[13] + x530*M[29] + x533*M[3] + x541*M[8] + x569*M[12] + x572*M[15] + x581*M[16] + x587*M[25] + x590*M[26] + x599*M[30] + x611*M[22] + x618*M[31] + x629*M[19] + x641*M[9] + x643*M[10] + x646*M[11] + x651*M[24] + x657*M[20] + x659*M[21] + x663*M[23] + x670*M[33] + x690*M[17] + x693*M[18] + x695*M[14] + x699*M[32] + x700*M[27] + x701*M[28];
#pragma omp atomic
L[80] += x372*M[0] + x486*M[1] + x510*M[4] + x513*M[6] + x519*M[12] + x541*M[7] + x563*M[3] + x564*M[9] + x581*M[15] + x590*M[25] + x599*M[29] + x618*M[30] + x623*M[2] + x633*M[8] + x635*M[5] + x638*M[14] + x646*M[10] + x651*M[23] + x659*M[20] + x662*M[19] + x663*M[22] + x670*M[32] + x690*M[16] + x693*M[17] + x695*M[13] + x699*M[31] + x700*M[26] + x701*M[27] + x703*M[11] + x705*M[21] + x707*M[24] + x712*M[18] + x719*M[28] + x725*M[33];
#pragma omp atomic
L[81] += x541*M[6] + x600*M[3] + x603*M[8] + x618*M[29] + x623*M[1] + x625*M[2] + x630*M[0] + x633*M[7] + x635*M[4] + x636*M[5] + x638*M[13] + x651*M[22] + x654*M[19] + x670*M[31] + x674*M[24] + x683*M[33] + x690*M[15] + x693*M[16] + x695*M[12] + x699*M[30] + x700*M[25] + x701*M[26] + x702*M[9] + x703*M[10] + x704*M[11] + x705*M[20] + x706*M[21] + x707*M[23] + x710*M[14] + x712*M[17] + x716*M[18] + x719*M[27] + x721*M[28] + x725*M[32];
#pragma omp atomic
L[82] += x517*M[0] + x603*M[7] + x625*M[1] + x633*M[6] + x634*M[3] + x636*M[4] + x637*M[9] + x638*M[12] + x670*M[30] + x672*M[2] + x674*M[23] + x675*M[5] + x677*M[8] + x678*M[14] + x683*M[32] + x693*M[15] + x698*M[19] + x699*M[29] + x701*M[25] + x704*M[10] + x706*M[20] + x707*M[22] + x710*M[13] + x712*M[16] + x713*M[11] + x716*M[17] + x718*M[21] + x719*M[26] + x721*M[27] + x724*M[24] + x725*M[31] + x726*M[18] + x729*M[33] + x730*M[28];
#pragma omp atomic
L[83] += x603*M[6] + x604*M[3] + x613*M[8] + x668*M[19] + x670*M[29] + x671*M[0] + x672*M[1] + x673*M[2] + x674*M[22] + x675*M[4] + x676*M[5] + x677*M[7] + x678*M[13] + x683*M[31] + x684*M[24] + x688*M[33] + x709*M[9] + x710*M[12] + x712*M[15] + x713*M[10] + x715*M[11] + x716*M[16] + x718*M[20] + x719*M[25] + x721*M[26] + x722*M[21] + x724*M[23] + x725*M[30] + x726*M[17] + x727*M[18] + x728*M[14] + x729*M[32] + x730*M[27] + x731*M[28];
#pragma omp atomic
L[84] += x153*M[1] + x154*M[2] + x187*M[7] + x226*M[0] + x234*M[3] + x241*M[4] + x242*M[5] + x255*M[10] + x256*M[11] + x260*M[13] + x275*M[9] + x384*M[6] + x409*M[15] + x419*M[16] + x498*M[12] + x527*M[14] + x567*M[8] + x584*M[17] + x593*M[18];
#pragma omp atomic
L[85] += x153*M[0] + x169*M[2] + x187*M[5] + x241*M[3] + x255*M[9] + x260*M[11] + x376*M[1] + x384*M[4] + x398*M[6] + x400*M[7] + x409*M[12] + x419*M[13] + x421*M[16] + x498*M[10] + x560*M[15] + x575*M[8] + x584*M[14] + x596*M[18] + x608*M[17];
#pragma omp atomic
L[86] += x154*M[0] + x169*M[1] + x187*M[4] + x242*M[3] + x256*M[9] + x260*M[10] + x400*M[6] + x419*M[12] + x421*M[15] + x431*M[2] + x527*M[11] + x567*M[5] + x575*M[7] + x578*M[8] + x584*M[13] + x593*M[14] + x596*M[17] + x608*M[16] + x620*M[18];
#pragma omp atomic
L[87] += x325*M[1] + x326*M[2] + x339*M[7] + x376*M[0] + x384*M[3] + x398*M[4] + x400*M[5] + x409*M[10] + x419*M[11] + x421*M[13] + x450*M[6] + x459*M[15] + x462*M[16] + x498*M[9] + x560*M[12] + x608*M[14] + x640*M[8] + x656*M[17] + x660*M[18];
#pragma omp atomic
L[88] += x169*M[0] + x187*M[3] + x260*M[9] + x326*M[1] + x339*M[6] + x400*M[4] + x419*M[10] + x421*M[12] + x462*M[15] + x501*M[2] + x522*M[8] + x575*M[5] + x584*M[11] + x596*M[14] + x608*M[13] + x640*M[7] + x656*M[16] + x660*M[17] + x696*M[18];
#pragma omp atomic
L[89] += x431*M[0] + x501*M[1] + x507*M[2] + x522*M[7] + x527*M[9] + x567*M[3] + x575*M[4] + x578*M[5] + x584*M[10] + x593*M[11] + x596*M[13] + x608*M[12] + x620*M[14] + x640*M[6] + x656*M[15] + x660*M[16] + x689*M[8] + x696*M[17] + x697*M[18];
#pragma omp atomic
L[90] += x325*M[0] + x332*M[2] + x339*M[5] + x398*M[3] + x409*M[9] + x421*M[11] + x444*M[1] + x450*M[4] + x453*M[6] + x454*M[7] + x459*M[12] + x462*M[13] + x463*M[16] + x553*M[15] + x560*M[10] + x629*M[17] + x644*M[8] + x656*M[14] + x662*M[18];
#pragma omp atomic
L[91] += x326*M[0] + x332*M[1] + x339*M[4] + x400*M[3] + x419*M[9] + x421*M[10] + x454*M[6] + x462*M[12] + x463*M[15] + x535*M[2] + x608*M[11] + x629*M[16] + x640*M[5] + x644*M[7] + x647*M[8] + x654*M[18] + x656*M[13] + x660*M[14] + x662*M[17];
#pragma omp atomic
L[92] += x501*M[0] + x516*M[2] + x522*M[5] + x535*M[1] + x575*M[3] + x584*M[9] + x596*M[11] + x608*M[10] + x629*M[15] + x640*M[4] + x644*M[6] + x647*M[7] + x654*M[17] + x656*M[12] + x660*M[13] + x662*M[16] + x691*M[8] + x696*M[14] + x698*M[18];
#pragma omp atomic
L[93] += x507*M[0] + x516*M[1] + x522*M[4] + x539*M[2] + x578*M[3] + x593*M[9] + x596*M[10] + x620*M[11] + x647*M[6] + x654*M[16] + x660*M[12] + x662*M[15] + x668*M[18] + x689*M[5] + x691*M[7] + x692*M[8] + x696*M[13] + x697*M[14] + x698*M[17];
#pragma omp atomic
L[94] += x318*M[1] + x329*M[2] + x343*M[7] + x392*M[6] + x414*M[15] + x418*M[16] + x444*M[0] + x450*M[3] + x453*M[4] + x454*M[5] + x459*M[10] + x462*M[11] + x463*M[13] + x553*M[12] + x560*M[9] + x629*M[14] + x641*M[8] + x657*M[17] + x659*M[18];
#pragma omp atomic
L[95] += x329*M[1] + x332*M[0] + x339*M[3] + x343*M[6] + x418*M[15] + x421*M[9] + x454*M[4] + x462*M[10] + x463*M[12] + x561*M[2] + x564*M[8] + x629*M[13] + x641*M[7] + x644*M[5] + x656*M[11] + x657*M[16] + x659*M[17] + x662*M[14] + x705*M[18];
#pragma omp atomic
L[96] += x535*M[0] + x561*M[1] + x562*M[2] + x564*M[7] + x608*M[9] + x629*M[12] + x640*M[3] + x641*M[6] + x644*M[4] + x647*M[5] + x654*M[14] + x656*M[10] + x657*M[15] + x659*M[16] + x660*M[11] + x662*M[13] + x702*M[8] + x705*M[17] + x706*M[18];
#pragma omp atomic
L[97] += x516*M[0] + x522*M[3] + x562*M[1] + x564*M[6] + x596*M[9] + x631*M[2] + x637*M[8] + x647*M[4] + x654*M[13] + x659*M[15] + x660*M[10] + x662*M[12] + x691*M[5] + x696*M[11] + x698*M[14] + x702*M[7] + x705*M[16] + x706*M[17] + x718*M[18];
#pragma omp atomic
L[98] += x539*M[0] + x620*M[9] + x631*M[1] + x632*M[2] + x637*M[7] + x654*M[12] + x668*M[14] + x689*M[3] + x691*M[4] + x692*M[5] + x696*M[10] + x697*M[11] + x698*M[13] + x702*M[6] + x705*M[15] + x706*M[16] + x709*M[8] + x718*M[17] + x722*M[18];
#pragma omp atomic
L[99] += x318*M[0] + x330*M[2] + x343*M[5] + x379*M[1] + x392*M[4] + x393*M[6] + x399*M[7] + x414*M[12] + x418*M[13] + x422*M[16] + x453*M[3] + x459*M[9] + x463*M[11] + x494*M[15] + x553*M[10] + x611*M[17] + x643*M[8] + x657*M[14] + x663*M[18];
#pragma omp atomic
L[100] += x329*M[0] + x330*M[1] + x343*M[4] + x399*M[6] + x418*M[12] + x422*M[15] + x454*M[3] + x462*M[9] + x463*M[10] + x533*M[2] + x611*M[16] + x629*M[11] + x641*M[5] + x643*M[7] + x646*M[8] + x651*M[18] + x657*M[13] + x659*M[14] + x663*M[17];
#pragma omp atomic
L[101] += x533*M[1] + x561*M[0] + x563*M[2] + x564*M[5] + x611*M[15] + x629*M[10] + x641*M[4] + x643*M[6] + x644*M[3] + x646*M[7] + x651*M[17] + x656*M[9] + x657*M[12] + x659*M[13] + x662*M[11] + x663*M[16] + x703*M[8] + x705*M[14] + x707*M[18];
#pragma omp atomic
L[102] += x562*M[0] + x563*M[1] + x564*M[4] + x600*M[2] + x646*M[6] + x647*M[3] + x651*M[16] + x654*M[11] + x659*M[12] + x660*M[9] + x662*M[10] + x663*M[15] + x674*M[18] + x702*M[5] + x703*M[7] + x704*M[8] + x705*M[13] + x706*M[14] + x707*M[17];
#pragma omp atomic
L[103] += x600*M[1] + x631*M[0] + x634*M[2] + x637*M[5] + x651*M[15] + x654*M[10] + x674*M[17] + x691*M[3] + x696*M[9] + x698*M[11] + x702*M[4] + x703*M[6] + x704*M[7] + x705*M[12] + x706*M[13] + x707*M[16] + x713*M[8] + x718*M[14] + x724*M[18];
#pragma omp atomic
L[104] += x604*M[2] + x632*M[0] + x634*M[1] + x637*M[4] + x668*M[11] + x674*M[16] + x684*M[18] + x692*M[3] + x697*M[9] + x698*M[10] + x704*M[6] + x706*M[12] + x707*M[15] + x709*M[5] + x713*M[7] + x715*M[8] + x718*M[13] + x722*M[14] + x724*M[17];
#pragma omp atomic
L[105] += x172*M[1] + x176*M[2] - x199*M[7] - x248*M[6] + x269*M[15] + x272*M[16] + x379*M[0] + x392*M[3] + x393*M[4] + x399*M[5] + x414*M[10] + x418*M[11] + x422*M[13] + x494*M[12] + x553*M[9] + x569*M[8] + x587*M[17] + x590*M[18] + x611*M[14];
#pragma omp atomic
L[106] += x176*M[1] - x199*M[6] + x272*M[15] + x330*M[0] + x343*M[3] + x399*M[4] + x418*M[10] + x422*M[12] + x463*M[9] + x504*M[2] + x519*M[8] + x569*M[7] + x587*M[16] + x590*M[17] + x611*M[13] + x643*M[5] + x657*M[11] + x663*M[14] + x700*M[18];
#pragma omp atomic
L[107] += x504*M[1] + x510*M[2] + x519*M[7] + x533*M[0] + x569*M[6] + x587*M[15] + x590*M[16] + x611*M[12] + x629*M[9] + x641*M[3] + x643*M[4] + x646*M[5] + x651*M[14] + x657*M[10] + x659*M[11] + x663*M[13] + x695*M[8] + x700*M[17] + x701*M[18];
#pragma omp atomic
L[108] += x510*M[1] + x519*M[6] + x563*M[0] + x564*M[3] + x590*M[15] + x635*M[2] + x638*M[8] + x646*M[4] + x651*M[13] + x659*M[10] + x662*M[9] + x663*M[12] + x695*M[7] + x700*M[16] + x701*M[17] + x703*M[5] + x705*M[11] + x707*M[14] + x719*M[18];
#pragma omp atomic
L[109] += x600*M[0] + x635*M[1] + x636*M[2] + x638*M[7] + x651*M[12] + x654*M[9] + x674*M[14] + x695*M[6] + x700*M[15] + x701*M[16] + x702*M[3] + x703*M[4] + x704*M[5] + x705*M[10] + x706*M[11] + x707*M[13] + x710*M[8] + x719*M[17] + x721*M[18];
#pragma omp atomic
L[110] += x634*M[0] + x636*M[1] + x637*M[3] + x638*M[6] + x674*M[13] + x675*M[2] + x678*M[8] + x698*M[9] + x701*M[15] + x704*M[4] + x706*M[10] + x707*M[12] + x710*M[7] + x713*M[5] + x718*M[11] + x719*M[16] + x721*M[17] + x724*M[14] + x730*M[18];
#pragma omp atomic
L[111] += x604*M[0] + x668*M[9] + x674*M[12] + x675*M[1] + x676*M[2] + x678*M[7] + x684*M[14] + x709*M[3] + x710*M[6] + x713*M[4] + x715*M[5] + x718*M[10] + x719*M[15] + x721*M[16] + x722*M[11] + x724*M[13] + x728*M[8] + x730*M[17] + x731*M[18];
#pragma omp atomic
L[112] += x163*M[2] + x172*M[0] - x199*M[5] + x229*M[1] + x237*M[6] + x247*M[7] - x248*M[4] + x264*M[16] + x269*M[12] + x272*M[13] + x277*M[15] + x393*M[3] + x414*M[9] + x422*M[11] + x494*M[10] + x530*M[17] + x572*M[8] + x587*M[14] + x599*M[18];
#pragma omp atomic
L[113] += x163*M[1] + x176*M[0] - x199*M[4] + x247*M[6] + x264*M[15] + x272*M[12] + x399*M[3] + x418*M[9] + x422*M[10] + x428*M[2] + x530*M[16] + x569*M[5] + x572*M[7] + x581*M[8] + x587*M[13] + x590*M[14] + x599*M[17] + x611*M[11] + x618*M[18];
#pragma omp atomic
L[114] += x428*M[1] + x504*M[0] + x513*M[2] + x519*M[5] + x530*M[15] + x569*M[4] + x572*M[6] + x581*M[7] + x587*M[12] + x590*M[13] + x599*M[16] + x611*M[10] + x618*M[17] + x643*M[3] + x657*M[9] + x663*M[11] + x690*M[8] + x699*M[18] + x700*M[14];
#pragma omp atomic
L[115] += x510*M[0] + x513*M[1] + x519*M[4] + x541*M[2] + x581*M[6] + x590*M[12] + x599*M[15] + x618*M[16] + x646*M[3] + x651*M[11] + x659*M[9] + x663*M[10] + x670*M[18] + x690*M[7] + x693*M[8] + x695*M[5] + x699*M[17] + x700*M[13] + x701*M[14];
#pragma omp atomic
L[116] += x541*M[1] + x618*M[15] + x633*M[2] + x635*M[0] + x638*M[5] + x651*M[10] + x670*M[17] + x690*M[6] + x693*M[7] + x695*M[4] + x699*M[16] + x700*M[12] + x701*M[13] + x703*M[3] + x705*M[9] + x707*M[11] + x712*M[8] + x719*M[14] + x725*M[18];
#pragma omp atomic
L[117] += x603*M[2] + x633*M[1] + x636*M[0] + x638*M[4] + x670*M[16] + x674*M[11] + x683*M[18] + x693*M[6] + x699*M[15] + x701*M[12] + x704*M[3] + x706*M[9] + x707*M[10] + x710*M[5] + x712*M[7] + x716*M[8] + x719*M[13] + x721*M[14] + x725*M[17];
#pragma omp atomic
L[118] += x603*M[1] + x670*M[15] + x674*M[10] + x675*M[0] + x677*M[2] + x678*M[5] + x683*M[17] + x710*M[4] + x712*M[6] + x713*M[3] + x716*M[7] + x718*M[9] + x719*M[12] + x721*M[13] + x724*M[11] + x725*M[16] + x726*M[8] + x729*M[18] + x730*M[14];
#pragma omp atomic
L[119] += x613*M[2] + x676*M[0] + x677*M[1] + x678*M[4] + x683*M[16] + x684*M[11] + x688*M[18] + x715*M[3] + x716*M[6] + x721*M[12] + x722*M[9] + x724*M[10] + x725*M[15] + x726*M[7] + x727*M[8] + x728*M[5] + x729*M[17] + x730*M[13] + x731*M[14];
#pragma omp atomic
L[120] += x234*M[0] + x241*M[1] + x242*M[2] + x255*M[4] + x256*M[5] + x260*M[7] + x275*M[3] + x498*M[6] + x527*M[8];
#pragma omp atomic
L[121] += x187*M[2] + x241*M[0] + x255*M[3] + x260*M[5] + x384*M[1] + x409*M[6] + x419*M[7] + x498*M[4] + x584*M[8];
#pragma omp atomic
L[122] += x187*M[1] + x242*M[0] + x256*M[3] + x260*M[4] + x419*M[6] + x527*M[5] + x567*M[2] + x584*M[7] + x593*M[8];
#pragma omp atomic
L[123] += x384*M[0] + x398*M[1] + x400*M[2] + x409*M[4] + x419*M[5] + x421*M[7] + x498*M[3] + x560*M[6] + x608*M[8];
#pragma omp atomic
L[124] += x187*M[0] + x260*M[3] + x400*M[1] + x419*M[4] + x421*M[6] + x575*M[2] + x584*M[5] + x596*M[8] + x608*M[7];
#pragma omp atomic
L[125] += x527*M[3] + x567*M[0] + x575*M[1] + x578*M[2] + x584*M[4] + x593*M[5] + x596*M[7] + x608*M[6] + x620*M[8];
#pragma omp atomic
L[126] += x339*M[2] + x398*M[0] + x409*M[3] + x421*M[5] + x450*M[1] + x459*M[6] + x462*M[7] + x560*M[4] + x656*M[8];
#pragma omp atomic
L[127] += x339*M[1] + x400*M[0] + x419*M[3] + x421*M[4] + x462*M[6] + x608*M[5] + x640*M[2] + x656*M[7] + x660*M[8];
#pragma omp atomic
L[128] += x522*M[2] + x575*M[0] + x584*M[3] + x596*M[5] + x608*M[4] + x640*M[1] + x656*M[6] + x660*M[7] + x696*M[8];
#pragma omp atomic
L[129] += x522*M[1] + x578*M[0] + x593*M[3] + x596*M[4] + x620*M[5] + x660*M[6] + x689*M[2] + x696*M[7] + x697*M[8];
#pragma omp atomic
L[130] += x450*M[0] + x453*M[1] + x454*M[2] + x459*M[4] + x462*M[5] + x463*M[7] + x553*M[6] + x560*M[3] + x629*M[8];
#pragma omp atomic
L[131] += x339*M[0] + x421*M[3] + x454*M[1] + x462*M[4] + x463*M[6] + x629*M[7] + x644*M[2] + x656*M[5] + x662*M[8];
#pragma omp atomic
L[132] += x608*M[3] + x629*M[6] + x640*M[0] + x644*M[1] + x647*M[2] + x654*M[8] + x656*M[4] + x660*M[5] + x662*M[7];
#pragma omp atomic
L[133] += x522*M[0] + x596*M[3] + x647*M[1] + x654*M[7] + x660*M[4] + x662*M[6] + x691*M[2] + x696*M[5] + x698*M[8];
#pragma omp atomic
L[134] += x620*M[3] + x654*M[6] + x668*M[8] + x689*M[0] + x691*M[1] + x692*M[2] + x696*M[4] + x697*M[5] + x698*M[7];
#pragma omp atomic
L[135] += x343*M[2] + x392*M[1] + x414*M[6] + x418*M[7] + x453*M[0] + x459*M[3] + x463*M[5] + x553*M[4] + x657*M[8];
#pragma omp atomic
L[136] += x343*M[1] + x418*M[6] + x454*M[0] + x462*M[3] + x463*M[4] + x629*M[5] + x641*M[2] + x657*M[7] + x659*M[8];
#pragma omp atomic
L[137] += x564*M[2] + x629*M[4] + x641*M[1] + x644*M[0] + x656*M[3] + x657*M[6] + x659*M[7] + x662*M[5] + x705*M[8];
#pragma omp atomic
L[138] += x564*M[1] + x647*M[0] + x654*M[5] + x659*M[6] + x660*M[3] + x662*M[4] + x702*M[2] + x705*M[7] + x706*M[8];
#pragma omp atomic
L[139] += x637*M[2] + x654*M[4] + x691*M[0] + x696*M[3] + x698*M[5] + x702*M[1] + x705*M[6] + x706*M[7] + x718*M[8];
#pragma omp atomic
L[140] += x637*M[1] + x668*M[5] + x692*M[0] + x697*M[3] + x698*M[4] + x706*M[6] + x709*M[2] + x718*M[7] + x722*M[8];
#pragma omp atomic
L[141] += x392*M[0] + x393*M[1] + x399*M[2] + x414*M[4] + x418*M[5] + x422*M[7] + x494*M[6] + x553*M[3] + x611*M[8];
#pragma omp atomic
L[142] += x343*M[0] + x399*M[1] + x418*M[4] + x422*M[6] + x463*M[3] + x611*M[7] + x643*M[2] + x657*M[5] + x663*M[8];
#pragma omp atomic
L[143] += x611*M[6] + x629*M[3] + x641*M[0] + x643*M[1] + x646*M[2] + x651*M[8] + x657*M[4] + x659*M[5] + x663*M[7];
#pragma omp atomic
L[144] += x564*M[0] + x646*M[1] + x651*M[7] + x659*M[4] + x662*M[3] + x663*M[6] + x703*M[2] + x705*M[5] + x707*M[8];
#pragma omp atomic
L[145] += x651*M[6] + x654*M[3] + x674*M[8] + x702*M[0] + x703*M[1] + x704*M[2] + x705*M[4] + x706*M[5] + x707*M[7];
#pragma omp atomic
L[146] += x637*M[0] + x674*M[7] + x698*M[3] + x704*M[1] + x706*M[4] + x707*M[6] + x713*M[2] + x718*M[5] + x724*M[8];
#pragma omp atomic
L[147] += x668*M[3] + x674*M[6] + x684*M[8] + x709*M[0] + x713*M[1] + x715*M[2] + x718*M[4] + x722*M[5] + x724*M[7];
#pragma omp atomic
L[148] += -x199*M[2] - x248*M[1] + x269*M[6] + x272*M[7] + x393*M[0] + x414*M[3] + x422*M[5] + x494*M[4] + x587*M[8];
#pragma omp atomic
L[149] += -x199*M[1] + x272*M[6] + x399*M[0] + x418*M[3] + x422*M[4] + x569*M[2] + x587*M[7] + x590*M[8] + x611*M[5];
#pragma omp atomic
L[150] += x519*M[2] + x569*M[1] + x587*M[6] + x590*M[7] + x611*M[4] + x643*M[0] + x657*M[3] + x663*M[5] + x700*M[8];
#pragma omp atomic
L[151] += x519*M[1] + x590*M[6] + x646*M[0] + x651*M[5] + x659*M[3] + x663*M[4] + x695*M[2] + x700*M[7] + x701*M[8];
#pragma omp atomic
L[152] += x638*M[2] + x651*M[4] + x695*M[1] + x700*M[6] + x701*M[7] + x703*M[0] + x705*M[3] + x707*M[5] + x719*M[8];
#pragma omp atomic
L[153] += x638*M[1] + x674*M[5] + x701*M[6] + x704*M[0] + x706*M[3] + x707*M[4] + x710*M[2] + x719*M[7] + x721*M[8];
#pragma omp atomic
L[154] += x674*M[4] + x678*M[2] + x710*M[1] + x713*M[0] + x718*M[3] + x719*M[6] + x721*M[7] + x724*M[5] + x730*M[8];
#pragma omp atomic
L[155] += x678*M[1] + x684*M[5] + x715*M[0] + x721*M[6] + x722*M[3] + x724*M[4] + x728*M[2] + x730*M[7] + x731*M[8];
#pragma omp atomic
L[156] += x237*M[1] + x247*M[2] - x248*M[0] + x264*M[7] + x269*M[4] + x272*M[5] + x277*M[6] + x494*M[3] + x530*M[8];
#pragma omp atomic
L[157] += -x199*M[0] + x247*M[1] + x264*M[6] + x272*M[4] + x422*M[3] + x530*M[7] + x572*M[2] + x587*M[5] + x599*M[8];
#pragma omp atomic
L[158] += x530*M[6] + x569*M[0] + x572*M[1] + x581*M[2] + x587*M[4] + x590*M[5] + x599*M[7] + x611*M[3] + x618*M[8];
#pragma omp atomic
L[159] += x519*M[0] + x581*M[1] + x590*M[4] + x599*M[6] + x618*M[7] + x663*M[3] + x690*M[2] + x699*M[8] + x700*M[5];
#pragma omp atomic
L[160] += x618*M[6] + x651*M[3] + x670*M[8] + x690*M[1] + x693*M[2] + x695*M[0] + x699*M[7] + x700*M[4] + x701*M[5];
#pragma omp atomic
L[161] += x638*M[0] + x670*M[7] + x693*M[1] + x699*M[6] + x701*M[4] + x707*M[3] + x712*M[2] + x719*M[5] + x725*M[8];
#pragma omp atomic
L[162] += x670*M[6] + x674*M[3] + x683*M[8] + x710*M[0] + x712*M[1] + x716*M[2] + x719*M[4] + x721*M[5] + x725*M[7];
#pragma omp atomic
L[163] += x678*M[0] + x683*M[7] + x716*M[1] + x721*M[4] + x724*M[3] + x725*M[6] + x726*M[2] + x729*M[8] + x730*M[5];
#pragma omp atomic
L[164] += x683*M[6] + x684*M[3] + x688*M[8] + x726*M[1] + x727*M[2] + x728*M[0] + x729*M[7] + x730*M[4] + x731*M[5];
#pragma omp atomic
L[165] += x255*M[1] + x256*M[2] + x275*M[0];
#pragma omp atomic
L[166] += x255*M[0] + x260*M[2] + x498*M[1];
#pragma omp atomic
L[167] += x256*M[0] + x260*M[1] + x527*M[2];
#pragma omp atomic
L[168] += x409*M[1] + x419*M[2] + x498*M[0];
#pragma omp atomic
L[169] += x260*M[0] + x419*M[1] + x584*M[2];
#pragma omp atomic
L[170] += x527*M[0] + x584*M[1] + x593*M[2];
#pragma omp atomic
L[171] += x409*M[0] + x421*M[2] + x560*M[1];
#pragma omp atomic
L[172] += x419*M[0] + x421*M[1] + x608*M[2];
#pragma omp atomic
L[173] += x584*M[0] + x596*M[2] + x608*M[1];
#pragma omp atomic
L[174] += x593*M[0] + x596*M[1] + x620*M[2];
#pragma omp atomic
L[175] += x459*M[1] + x462*M[2] + x560*M[0];
#pragma omp atomic
L[176] += x421*M[0] + x462*M[1] + x656*M[2];
#pragma omp atomic
L[177] += x608*M[0] + x656*M[1] + x660*M[2];
#pragma omp atomic
L[178] += x596*M[0] + x660*M[1] + x696*M[2];
#pragma omp atomic
L[179] += x620*M[0] + x696*M[1] + x697*M[2];
#pragma omp atomic
L[180] += x459*M[0] + x463*M[2] + x553*M[1];
#pragma omp atomic
L[181] += x462*M[0] + x463*M[1] + x629*M[2];
#pragma omp atomic
L[182] += x629*M[1] + x656*M[0] + x662*M[2];
#pragma omp atomic
L[183] += x654*M[2] + x660*M[0] + x662*M[1];
#pragma omp atomic
L[184] += x654*M[1] + x696*M[0] + x698*M[2];
#pragma omp atomic
L[185] += x668*M[2] + x697*M[0] + x698*M[1];
#pragma omp atomic
L[186] += x414*M[1] + x418*M[2] + x553*M[0];
#pragma omp atomic
L[187] += x418*M[1] + x463*M[0] + x657*M[2];
#pragma omp atomic
L[188] += x629*M[0] + x657*M[1] + x659*M[2];
#pragma omp atomic
L[189] += x659*M[1] + x662*M[0] + x705*M[2];
#pragma omp atomic
L[190] += x654*M[0] + x705*M[1] + x706*M[2];
#pragma omp atomic
L[191] += x698*M[0] + x706*M[1] + x718*M[2];
#pragma omp atomic
L[192] += x668*M[0] + x718*M[1] + x722*M[2];
#pragma omp atomic
L[193] += x414*M[0] + x422*M[2] + x494*M[1];
#pragma omp atomic
L[194] += x418*M[0] + x422*M[1] + x611*M[2];
#pragma omp atomic
L[195] += x611*M[1] + x657*M[0] + x663*M[2];
#pragma omp atomic
L[196] += x651*M[2] + x659*M[0] + x663*M[1];
#pragma omp atomic
L[197] += x651*M[1] + x705*M[0] + x707*M[2];
#pragma omp atomic
L[198] += x674*M[2] + x706*M[0] + x707*M[1];
#pragma omp atomic
L[199] += x674*M[1] + x718*M[0] + x724*M[2];
#pragma omp atomic
L[200] += x684*M[2] + x722*M[0] + x724*M[1];
#pragma omp atomic
L[201] += x269*M[1] + x272*M[2] + x494*M[0];
#pragma omp atomic
L[202] += x272*M[1] + x422*M[0] + x587*M[2];
#pragma omp atomic
L[203] += x587*M[1] + x590*M[2] + x611*M[0];
#pragma omp atomic
L[204] += x590*M[1] + x663*M[0] + x700*M[2];
#pragma omp atomic
L[205] += x651*M[0] + x700*M[1] + x701*M[2];
#pragma omp atomic
L[206] += x701*M[1] + x707*M[0] + x719*M[2];
#pragma omp atomic
L[207] += x674*M[0] + x719*M[1] + x721*M[2];
#pragma omp atomic
L[208] += x721*M[1] + x724*M[0] + x730*M[2];
#pragma omp atomic
L[209] += x684*M[0] + x730*M[1] + x731*M[2];
#pragma omp atomic
L[210] += x264*M[2] + x269*M[0] + x277*M[1];
#pragma omp atomic
L[211] += x264*M[1] + x272*M[0] + x530*M[2];
#pragma omp atomic
L[212] += x530*M[1] + x587*M[0] + x599*M[2];
#pragma omp atomic
L[213] += x590*M[0] + x599*M[1] + x618*M[2];
#pragma omp atomic
L[214] += x618*M[1] + x699*M[2] + x700*M[0];
#pragma omp atomic
L[215] += x670*M[2] + x699*M[1] + x701*M[0];
#pragma omp atomic
L[216] += x670*M[1] + x719*M[0] + x725*M[2];
#pragma omp atomic
L[217] += x683*M[2] + x721*M[0] + x725*M[1];
#pragma omp atomic
L[218] += x683*M[1] + x729*M[2] + x730*M[0];
#pragma omp atomic
L[219] += x688*M[2] + x729*M[1] + x731*M[0];

}

void L2L_10(double x, double y, double z, double * L, double * Ls) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
double x448;
double x449;
double x450;
double x451;
double x452;
double x453;
double x454;
double x455;
double x456;
double x457;
double x458;
double x459;
double x460;
double x461;
double x462;
double x463;
double x464;
double x465;
double x466;
double x467;
double x468;
double x469;
double x470;
double x471;
double x472;
double x473;
double x474;
double x475;
double x476;
double x477;
double x478;
double x479;
double x480;
double x481;
double x482;
double x483;
double x484;
double x485;
double x486;
double x487;
double x488;
double x489;
double x490;
double x491;
double x492;
double x493;
double x494;
double x495;
double x496;
double x497;
double x498;
double x499;
double x500;
double x501;
double x502;
double x503;
double x504;
double x505;
double x506;
double x507;
double x508;
double x509;
double x510;
double x511;
double x512;
double x513;
double x514;
double x515;
double x516;
double x517;
double x518;
double x519;
double x520;
double x521;
double x522;
double x523;
double x524;
double x525;
double x526;
double x527;
double x528;
double x529;
double x530;
double x531;
double x532;
double x533;
double x534;
double x535;
double x536;
double x537;
double x538;
double x539;
double x540;
double x541;
double x542;
double x543;
double x544;
double x545;
double x546;
double x547;
double x548;
double x549;
double x550;
double x551;
double x552;
double x553;
double x554;
double x555;
double x556;
double x557;
double x558;
double x559;
double x560;
double x561;
double x562;
double x563;
double x564;
double x565;
double x566;
double x567;
double x568;
double x569;
double x570;
double x571;
double x572;
double x573;
double x574;
double x575;
double x576;
double x577;
double x578;
double x579;
double x580;
double x581;
double x582;
double x583;
double x584;
double x585;
double x586;
double x587;
double x588;
double x589;
double x590;
double x591;
double x592;
double x593;
double x594;
double x595;
double x596;
double x597;
double x598;
double x599;
double x600;
double x601;
double x602;
double x603;
double x604;
double x605;
double x606;
double x607;
double x608;
double x609;
double x610;
double x611;
double x612;
double x613;
double x614;
double x615;
double x616;
double x617;
double x618;
double x619;
double x620;
double x621;
double x622;
double x623;
double x624;
double x625;
double x626;
double x627;
double x628;
double x629;
double x630;
double x631;
double x632;
double x633;
double x634;
double x635;
double x636;
double x637;
double x638;
double x639;
double x640;
double x641;
double x642;
double x643;
double x644;
double x645;
double x646;
double x647;
double x648;
double x649;
double x650;
double x651;
double x652;
double x653;
double x654;
double x655;
double x656;
double x657;
double x658;
double x659;
double x660;
double x661;
double x662;
double x663;
double x664;
double x665;
double x666;
double x667;
double x668;
double x669;
double x670;
double x671;
double x672;
double x673;
double x674;
double x675;
double x676;
double x677;
double x678;
double x679;
double x680;
double x681;
double x682;
double x683;
double x684;
double x685;
double x686;
double x687;
double x688;
double x689;
double x690;
double x691;
double x692;
double x693;
double x694;
double x695;
double x696;
double x697;
double x698;
double x699;
double x700;
double x701;
double x702;
double x703;
double x704;
double x705;
double x706;
double x707;
double x708;
double x709;
double x710;
double x711;
double x712;
double x713;
double x714;
double x715;
double x716;
double x717;
double x718;
double x719;
double x720;
double x721;
double x722;
double x723;
double x724;
double x725;
double x726;
double x727;
double x728;
double x729;
double x730;
double x731;
double x732;
double x733;
double x734;
double x735;
double x736;
double x737;
double x738;
double x739;
double x740;
double x741;
double x742;
double x743;
double x744;
double x745;
double x746;
double x747;
double x748;
double x749;
double x750;
double x751;
double x752;
double x753;
double x754;
double x755;
double x756;
double x757;
double x758;
double x759;
double x760;
double x761;
double x762;
double x763;
double x764;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
x3 = z*L[14];
x4 = x3*y;
x5 = (x*x);
x6 = (1.0/2.0)*x5;
x7 = (x*x*x);
x8 = (1.0/6.0)*x7;
x9 = (x*x*x*x);
x10 = (1.0/24.0)*x9;
x11 = pow(x, 5);
x12 = (1.0/120.0)*x11;
x13 = pow(x, 6);
x14 = (1.0/720.0)*x13;
x15 = pow(x, 7);
x16 = (1.0/5040.0)*x15;
x17 = (1.0/40320.0)*pow(x, 8);
x18 = (y*y);
x19 = (1.0/2.0)*x18;
x20 = (y*y*y);
x21 = (1.0/6.0)*x20;
x22 = (y*y*y*y);
x23 = (1.0/24.0)*x22;
x24 = pow(y, 5);
x25 = (1.0/120.0)*x24;
x26 = pow(y, 6);
x27 = (1.0/720.0)*x26;
x28 = pow(y, 7);
x29 = (1.0/5040.0)*x28;
x30 = (1.0/40320.0)*pow(y, 8);
x31 = (z*z);
x32 = (1.0/2.0)*x31;
x33 = (z*z*z);
x34 = (1.0/6.0)*x33;
x35 = (z*z*z*z);
x36 = (1.0/24.0)*x35;
x37 = pow(z, 5);
x38 = (1.0/120.0)*x37;
x39 = pow(z, 6);
x40 = (1.0/720.0)*x39;
x41 = pow(z, 7);
x42 = (1.0/5040.0)*x41;
x43 = (1.0/40320.0)*pow(z, 8);
x44 = x*L[13];
x45 = x*L[26];
x46 = x*L[45];
x47 = x*L[71];
x48 = x*L[105];
x49 = x*L[148];
x50 = x*L[201];
x51 = x*L[15];
x52 = x*L[29];
x53 = x*L[49];
x54 = x*L[76];
x55 = x*L[111];
x56 = x*L[155];
x57 = x*L[209];
x58 = y*L[11];
x59 = z*L[12];
x60 = y*L[21];
x61 = z*L[22];
x62 = y*L[36];
x63 = z*L[37];
x64 = y*L[57];
x65 = z*L[58];
x66 = y*L[85];
x67 = z*L[86];
x68 = y*L[121];
x69 = z*L[122];
x70 = y*L[166];
x71 = z*L[167];
x72 = y*L[18];
x73 = y*L[33];
x74 = y*L[54];
x75 = y*L[82];
x76 = y*L[118];
x77 = y*L[163];
x78 = y*L[218];
x79 = z*L[17];
x80 = z*L[31];
x81 = z*L[51];
x82 = z*L[78];
x83 = z*L[113];
x84 = z*L[157];
x85 = z*L[211];
x86 = y*L[28];
x87 = x*x86;
x88 = y*L[48];
x89 = x*x88;
x90 = y*L[75];
x91 = x*x90;
x92 = y*L[110];
x93 = x*x92;
x94 = y*L[154];
x95 = x*x94;
x96 = y*L[208];
x97 = x*x96;
x98 = z*L[27];
x99 = x*x98;
x100 = z*L[46];
x101 = x*x100;
x102 = z*L[72];
x103 = x*x102;
x104 = z*L[106];
x105 = x*x104;
x106 = z*L[149];
x107 = x*x106;
x108 = z*L[202];
x109 = x*x108;
x110 = z*L[24];
x111 = x110*y;
x112 = z*L[39];
x113 = x112*y;
x114 = z*L[60];
x115 = x114*y;
x116 = z*L[88];
x117 = x116*y;
x118 = z*L[124];
x119 = x118*y;
x120 = z*L[169];
x121 = x120*y;
x122 = (1.0/4.0)*x5;
x123 = x122*x18;
x124 = (1.0/12.0)*x5;
x125 = x124*x20;
x126 = (1.0/48.0)*x5;
x127 = x126*x22;
x128 = (1.0/240.0)*x5;
x129 = x128*x24;
x130 = (1.0/1440.0)*x5;
x131 = x130*x26;
x132 = (1.0/10080.0)*x5;
x133 = x122*x31;
x134 = x124*x33;
x135 = x126*x35;
x136 = x128*x37;
x137 = x130*x39;
x138 = (1.0/12.0)*x7;
x139 = x138*x18;
x140 = (1.0/36.0)*x7;
x141 = x140*x20;
x142 = (1.0/144.0)*x7;
x143 = x142*x22;
x144 = (1.0/720.0)*x7;
x145 = x144*x24;
x146 = (1.0/4320.0)*x7;
x147 = x138*x31;
x148 = x140*x33;
x149 = x142*x35;
x150 = x144*x37;
x151 = (1.0/48.0)*x9;
x152 = x151*x18;
x153 = (1.0/144.0)*x9;
x154 = x153*x20;
x155 = (1.0/576.0)*x9;
x156 = x155*x22;
x157 = (1.0/2880.0)*x9;
x158 = x151*x31;
x159 = x153*x33;
x160 = x155*x35;
x161 = (1.0/240.0)*x11;
x162 = x161*x18;
x163 = (1.0/720.0)*x11;
x164 = x163*x20;
x165 = (1.0/2880.0)*x11;
x166 = x161*x31;
x167 = x163*x33;
x168 = (1.0/1440.0)*x13;
x169 = x168*x18;
x170 = (1.0/4320.0)*x13;
x171 = x168*x31;
x172 = (1.0/10080.0)*x15;
x173 = x18*x31;
x174 = (1.0/4.0)*x173;
x175 = x18*x33;
x176 = (1.0/12.0)*x175;
x177 = x18*x35;
x178 = (1.0/48.0)*x177;
x179 = x18*x37;
x180 = (1.0/240.0)*x179;
x181 = (1.0/1440.0)*x18*x39;
x182 = x20*x31;
x183 = (1.0/12.0)*x182;
x184 = x20*x33;
x185 = (1.0/36.0)*x184;
x186 = x20*x35;
x187 = (1.0/144.0)*x186;
x188 = (1.0/720.0)*x20*x37;
x189 = x22*x31;
x190 = (1.0/48.0)*x189;
x191 = x22*x33;
x192 = (1.0/144.0)*x191;
x193 = (1.0/576.0)*x22*x35;
x194 = x24*x31;
x195 = (1.0/240.0)*x194;
x196 = (1.0/720.0)*x24*x33;
x197 = (1.0/1440.0)*x26*x31;
x198 = x*L[47];
x199 = x*L[74];
x200 = x*L[109];
x201 = x*L[153];
x202 = x*L[207];
x203 = x*L[73];
x204 = x*L[108];
x205 = x*L[152];
x206 = x*L[206];
x207 = x*L[107];
x208 = x*L[151];
x209 = x*L[205];
x210 = x*L[150];
x211 = x*L[204];
x212 = x*L[203];
x213 = y*L[43];
x214 = y*L[69];
x215 = y*L[103];
x216 = y*L[146];
x217 = y*L[199];
x218 = z*L[42];
x219 = z*L[67];
x220 = z*L[100];
x221 = z*L[142];
x222 = z*L[194];
x223 = y*L[64];
x224 = y*L[97];
x225 = y*L[139];
x226 = y*L[191];
x227 = z*L[63];
x228 = z*L[95];
x229 = z*L[136];
x230 = z*L[187];
x231 = y*L[92];
x232 = y*L[133];
x233 = y*L[184];
x234 = z*L[91];
x235 = z*L[131];
x236 = z*L[181];
x237 = y*L[128];
x238 = y*L[178];
x239 = z*L[127];
x240 = z*L[176];
x241 = y*L[173];
x242 = z*L[172];
x243 = (1.0/8.0)*x173*x5;
x244 = (1.0/24.0)*x5;
x245 = x175*x244;
x246 = (1.0/96.0)*x5;
x247 = x177*x246;
x248 = (1.0/480.0)*x5;
x249 = x182*x244;
x250 = (1.0/72.0)*x184*x5;
x251 = (1.0/288.0)*x5;
x252 = x189*x246;
x253 = (1.0/24.0)*x173*x7;
x254 = (1.0/72.0)*x7;
x255 = x175*x254;
x256 = (1.0/288.0)*x7;
x257 = x182*x254;
x258 = (1.0/96.0)*x173*x9;
x259 = (1.0/288.0)*x9;
x260 = x*L[23];
x261 = x*L[41];
x262 = x*L[66];
x263 = x*L[99];
x264 = x*L[141];
x265 = x*L[193];
x266 = x*L[25];
x267 = x*L[44];
x268 = x*L[70];
x269 = x*L[104];
x270 = x*L[147];
x271 = x*L[200];
x272 = x*x213;
x273 = x*x214;
x274 = x*x215;
x275 = x*x216;
x276 = x*x217;
x277 = x*x218;
x278 = x*x219;
x279 = x*x220;
x280 = x*x221;
x281 = x*x222;
x282 = x*L[68];
x283 = x*L[102];
x284 = x*L[145];
x285 = x*L[198];
x286 = x*L[101];
x287 = x*L[144];
x288 = x*L[197];
x289 = x*L[143];
x290 = x*L[196];
x291 = x*L[195];
x292 = y*L[13];
x293 = x98*y;
x294 = x*L[28];
x295 = x*L[48];
x296 = x*L[75];
x297 = x*L[110];
x298 = x*L[154];
x299 = x*L[208];
x300 = y*L[23];
x301 = y*L[38];
x302 = y*L[59];
x303 = y*L[87];
x304 = y*L[123];
x305 = y*L[168];
x306 = y*L[32];
x307 = y*L[53];
x308 = y*L[81];
x309 = y*L[117];
x310 = y*L[162];
x311 = y*L[217];
x312 = y*L[47];
x313 = x*x312;
x314 = y*L[74];
x315 = x*x314;
x316 = y*L[109];
x317 = x*x316;
x318 = y*L[153];
x319 = x*x318;
x320 = y*L[207];
x321 = x*x320;
x322 = x218*y;
x323 = x227*y;
x324 = x234*y;
x325 = x239*y;
x326 = x242*y;
x327 = y*L[68];
x328 = y*L[102];
x329 = y*L[145];
x330 = y*L[198];
x331 = y*L[96];
x332 = y*L[138];
x333 = y*L[190];
x334 = y*L[132];
x335 = y*L[183];
x336 = y*L[177];
x337 = y*L[14];
x338 = z*L[15];
x339 = z*L[18];
x340 = z*L[28];
x341 = x340*y;
x342 = x*L[27];
x343 = x*L[46];
x344 = x*L[72];
x345 = x*L[106];
x346 = x*L[149];
x347 = x*L[202];
x348 = y*L[24];
x349 = z*L[25];
x350 = y*L[39];
x351 = z*L[40];
x352 = y*L[60];
x353 = z*L[61];
x354 = y*L[88];
x355 = z*L[89];
x356 = y*L[124];
x357 = z*L[125];
x358 = y*L[169];
x359 = z*L[170];
x360 = z*L[32];
x361 = z*L[52];
x362 = z*L[79];
x363 = z*L[114];
x364 = z*L[158];
x365 = z*L[212];
x366 = z*L[47];
x367 = x*x366;
x368 = z*L[73];
x369 = x*x368;
x370 = z*L[107];
x371 = x*x370;
x372 = z*L[150];
x373 = x*x372;
x374 = z*L[203];
x375 = x*x374;
x376 = z*L[43];
x377 = x376*y;
x378 = z*L[64];
x379 = x378*y;
x380 = z*L[92];
x381 = x380*y;
x382 = z*L[128];
x383 = x382*y;
x384 = z*L[173];
x385 = x384*y;
x386 = z*L[68];
x387 = z*L[101];
x388 = z*L[143];
x389 = z*L[195];
x390 = z*L[96];
x391 = z*L[137];
x392 = z*L[188];
x393 = z*L[132];
x394 = z*L[182];
x395 = z*L[177];
x396 = x*L[38];
x397 = x*L[62];
x398 = x*L[94];
x399 = x*L[135];
x400 = x*L[186];
x401 = x*L[40];
x402 = x*L[65];
x403 = x*L[98];
x404 = x*L[140];
x405 = x*L[192];
x406 = x*x223;
x407 = x*x224;
x408 = x*x225;
x409 = x*x226;
x410 = x*x227;
x411 = x*x228;
x412 = x*x229;
x413 = x*x230;
x414 = x*L[96];
x415 = x*L[138];
x416 = x*L[190];
x417 = x*L[137];
x418 = x*L[189];
x419 = x*L[188];
x420 = x*L[43];
x421 = x*L[69];
x422 = x*L[103];
x423 = x*L[146];
x424 = x*L[199];
x425 = x*x327;
x426 = x*x328;
x427 = x*x329;
x428 = x*x330;
x429 = x*L[42];
x430 = x*L[67];
x431 = x*L[100];
x432 = x*L[142];
x433 = x*L[194];
x434 = x*x386;
x435 = x*x387;
x436 = x*x388;
x437 = x*x389;
x438 = y*L[26];
x439 = x100*y;
x440 = y*L[41];
x441 = y*L[62];
x442 = y*L[90];
x443 = y*L[126];
x444 = y*L[171];
x445 = y*L[52];
x446 = y*L[80];
x447 = y*L[116];
x448 = y*L[161];
x449 = y*L[216];
x450 = y*L[73];
x451 = x*x450;
x452 = y*L[108];
x453 = x*x452;
x454 = y*L[152];
x455 = x*x454;
x456 = y*L[206];
x457 = x*x456;
x458 = x219*y;
x459 = x228*y;
x460 = x235*y;
x461 = x240*y;
x462 = y*L[101];
x463 = y*L[144];
x464 = y*L[197];
x465 = y*L[137];
x466 = y*L[189];
x467 = y*L[182];
x468 = y*L[27];
x469 = x366*y;
x470 = y*L[42];
x471 = y*L[63];
x472 = y*L[91];
x473 = y*L[127];
x474 = y*L[172];
x475 = x386*y;
x476 = x390*y;
x477 = x393*y;
x478 = x395*y;
x479 = z*L[29];
x480 = z*L[33];
x481 = z*L[48];
x482 = x481*y;
x483 = z*L[44];
x484 = z*L[65];
x485 = z*L[93];
x486 = z*L[129];
x487 = z*L[174];
x488 = z*L[53];
x489 = z*L[80];
x490 = z*L[115];
x491 = z*L[159];
x492 = z*L[213];
x493 = z*L[74];
x494 = x*x493;
x495 = z*L[108];
x496 = x*x495;
x497 = z*L[151];
x498 = x*x497;
x499 = z*L[204];
x500 = x*x499;
x501 = z*L[69];
x502 = x501*y;
x503 = z*L[97];
x504 = x503*y;
x505 = z*L[133];
x506 = x505*y;
x507 = z*L[178];
x508 = x507*y;
x509 = z*L[102];
x510 = z*L[144];
x511 = z*L[196];
x512 = z*L[138];
x513 = z*L[189];
x514 = z*L[183];
x515 = x*L[59];
x516 = x*L[90];
x517 = x*L[130];
x518 = x*L[180];
x519 = x*L[61];
x520 = x*L[93];
x521 = x*L[134];
x522 = x*L[185];
x523 = x*x231;
x524 = x*x232;
x525 = x*x233;
x526 = x*x234;
x527 = x*x235;
x528 = x*x236;
x529 = x*L[132];
x530 = x*L[183];
x531 = x*L[182];
x532 = x*L[64];
x533 = x*L[97];
x534 = x*L[139];
x535 = x*L[191];
x536 = x*x331;
x537 = x*x332;
x538 = x*x333;
x539 = x*L[63];
x540 = x*L[95];
x541 = x*L[136];
x542 = x*L[187];
x543 = x*x390;
x544 = x*x391;
x545 = x*x392;
x546 = x*x462;
x547 = x*x463;
x548 = x*x464;
x549 = x*x509;
x550 = x*x510;
x551 = x*x511;
x552 = y*L[45];
x553 = x102*y;
x554 = y*L[66];
x555 = y*L[94];
x556 = y*L[130];
x557 = y*L[175];
x558 = y*L[79];
x559 = y*L[115];
x560 = y*L[160];
x561 = y*L[215];
x562 = y*L[107];
x563 = x*x562;
x564 = y*L[151];
x565 = x*x564;
x566 = y*L[205];
x567 = x*x566;
x568 = x220*y;
x569 = x229*y;
x570 = x236*y;
x571 = y*L[143];
x572 = y*L[196];
x573 = y*L[188];
x574 = y*L[46];
x575 = x368*y;
x576 = y*L[67];
x577 = y*L[95];
x578 = y*L[131];
x579 = y*L[176];
x580 = x387*y;
x581 = x391*y;
x582 = x394*y;
x583 = x493*y;
x584 = x509*y;
x585 = x512*y;
x586 = x514*y;
x587 = z*L[49];
x588 = z*L[54];
x589 = z*L[75];
x590 = x589*y;
x591 = z*L[70];
x592 = z*L[98];
x593 = z*L[134];
x594 = z*L[179];
x595 = z*L[81];
x596 = z*L[116];
x597 = z*L[160];
x598 = z*L[214];
x599 = z*L[109];
x600 = x*x599;
x601 = z*L[152];
x602 = x*x601;
x603 = z*L[205];
x604 = x*x603;
x605 = z*L[103];
x606 = x605*y;
x607 = z*L[139];
x608 = x607*y;
x609 = z*L[184];
x610 = x609*y;
x611 = z*L[145];
x612 = z*L[197];
x613 = z*L[190];
x614 = x*L[87];
x615 = x*L[126];
x616 = x*L[175];
x617 = x*L[89];
x618 = x*L[129];
x619 = x*L[179];
x620 = x*x237;
x621 = x*x238;
x622 = x*x239;
x623 = x*x240;
x624 = x*L[177];
x625 = x*L[92];
x626 = x*L[133];
x627 = x*L[184];
x628 = x*x334;
x629 = x*x335;
x630 = x*L[91];
x631 = x*L[131];
x632 = x*L[181];
x633 = x*x393;
x634 = x*x394;
x635 = x*x465;
x636 = x*x466;
x637 = x*x512;
x638 = x*x513;
x639 = x*x571;
x640 = x*x572;
x641 = x*x611;
x642 = x*x612;
x643 = y*L[71];
x644 = x104*y;
x645 = y*L[99];
x646 = y*L[135];
x647 = y*L[180];
x648 = y*L[114];
x649 = y*L[159];
x650 = y*L[214];
x651 = y*L[150];
x652 = x*x651;
x653 = y*L[204];
x654 = x*x653;
x655 = x221*y;
x656 = x230*y;
x657 = y*L[195];
x658 = y*L[72];
x659 = x370*y;
x660 = y*L[100];
x661 = y*L[136];
x662 = y*L[181];
x663 = x388*y;
x664 = x392*y;
x665 = x495*y;
x666 = x510*y;
x667 = x513*y;
x668 = x599*y;
x669 = x611*y;
x670 = x613*y;
x671 = z*L[76];
x672 = z*L[82];
x673 = z*L[110];
x674 = x673*y;
x675 = z*L[104];
x676 = z*L[140];
x677 = z*L[185];
x678 = z*L[117];
x679 = z*L[161];
x680 = z*L[215];
x681 = z*L[153];
x682 = x*x681;
x683 = z*L[206];
x684 = x*x683;
x685 = z*L[146];
x686 = x685*y;
x687 = z*L[191];
x688 = x687*y;
x689 = z*L[198];
x690 = x*L[123];
x691 = x*L[171];
x692 = x*L[125];
x693 = x*L[174];
x694 = x*x241;
x695 = x*x242;
x696 = x*L[128];
x697 = x*L[178];
x698 = x*x336;
x699 = x*L[127];
x700 = x*L[176];
x701 = x*x395;
x702 = x*x467;
x703 = x*x514;
x704 = x*x573;
x705 = x*x613;
x706 = x*x657;
x707 = x*x689;
x708 = y*L[105];
x709 = x106*y;
x710 = y*L[141];
x711 = y*L[186];
x712 = y*L[158];
x713 = y*L[213];
x714 = y*L[203];
x715 = x*x714;
x716 = x222*y;
x717 = y*L[106];
x718 = x372*y;
x719 = y*L[142];
x720 = y*L[187];
x721 = x389*y;
x722 = x497*y;
x723 = x511*y;
x724 = x601*y;
x725 = x612*y;
x726 = x681*y;
x727 = x689*y;
x728 = z*L[111];
x729 = z*L[118];
x730 = z*L[154];
x731 = x730*y;
x732 = z*L[147];
x733 = z*L[192];
x734 = z*L[162];
x735 = z*L[216];
x736 = z*L[207];
x737 = x*x736;
x738 = z*L[199];
x739 = x738*y;
x740 = x*L[168];
x741 = x*L[170];
x742 = x*L[173];
x743 = x*L[172];
x744 = y*L[148];
x745 = x108*y;
x746 = y*L[193];
x747 = y*L[212];
x748 = y*L[149];
x749 = x374*y;
x750 = y*L[194];
x751 = x499*y;
x752 = x603*y;
x753 = x683*y;
x754 = x736*y;
x755 = z*L[155];
x756 = z*L[163];
x757 = z*L[208];
x758 = x757*y;
x759 = z*L[200];
x760 = z*L[217];
x761 = y*L[201];
x762 = y*L[202];
x763 = z*L[209];
x764 = z*L[218];
#pragma omp atomic
Ls[0] += (1.0/362880.0)*pow(x, 9)*L[165] + x*x0 + x*x1 + x*x4 + x*L[1] + x10*x115 + x10*x62 + x10*x63 + x10*L[20] + x101*x21 + x103*x23 + x105*x25 + x107*x27 + x109*x29 + (1.0/480.0)*x11*x173*L[177] + x111*x6 + x113*x8 + x117*x12 + x119*x14 + x12*x64 + x12*x65 + x12*L[35] + x121*x16 + x123*x218 + x123*L[23] + x125*x219 + x125*L[41] + x127*x220 + x127*L[66] + x129*x221 + x129*L[99] + x131*x222 + x131*L[141] + x132*x28*L[193] + x132*x41*L[200] + x133*x213 + x133*L[25] + x134*x214 + x134*L[44] + x135*x215 + x135*L[70] + x136*x216 + x136*L[104] + x137*x217 + x137*L[147] + x139*x227 + x139*L[38] + x14*x66 + x14*x67 + x14*L[56] + x141*x228 + x141*L[62] + x143*x229 + x143*L[94] + x145*x230 + x145*L[135] + x146*x26*L[186] + x146*x39*L[192] + x147*x223 + x147*L[40] + x148*x224 + x148*L[65] + x149*x225 + x149*L[98] + x150*x226 + x150*L[140] + x152*x234 + x152*L[59] + x154*x235 + x154*L[90] + x156*x236 + x156*L[130] + x157*x24*L[180] + x157*x37*L[185] + x158*x231 + x158*L[61] + x159*x232 + x159*L[93] + x16*x68 + x16*x69 + x16*L[84] + x160*x233 + x160*L[134] + x162*x239 + x162*L[87] + x164*x240 + x164*L[126] + x165*x22*L[175] + x165*x35*L[179] + x166*x237 + x166*L[89] + x167*x238 + x167*L[129] + x169*x242 + x169*L[123] + x17*x70 + x17*x71 + x17*L[120] + x170*x20*L[171] + x170*x33*L[174] + x171*x241 + x171*L[125] + x172*x18*L[168] + x172*x31*L[170] + x174*x198 + x174*L[32] + x175*x259*L[183] + x176*x199 + x176*L[53] + x177*x256*L[190] + x178*x200 + x178*L[81] + x179*x248*L[198] + (1.0/10080.0)*x18*x41*L[217] + x180*x201 + x180*L[117] + x181*x202 + x181*L[162] + x182*x259*L[182] + x183*x203 + x183*L[52] + (1.0/216.0)*x184*x7*L[189] + x185*x204 + x185*L[80] + x186*x251*L[197] + x187*x205 + x187*L[116] + x188*x206 + x188*L[161] + x189*x256*L[188] + x19*x44 + x19*x79 + x19*x99 + x19*L[7] + x190*x207 + x190*L[79] + x191*x251*L[196] + x192*x208 + x192*L[115] + x193*x209 + x193*L[160] + x194*x248*L[195] + x195*x210 + x195*L[114] + x196*x211 + x196*L[159] + x197*x212 + x197*L[158] + x2*y + (1.0/4320.0)*x20*x39*L[216] + x21*x45 + x21*x80 + x21*L[16] + (1.0/2880.0)*x22*x37*L[215] + x23*x46 + x23*x81 + x23*L[30] + (1.0/2880.0)*x24*x35*L[214] + x243*L[68] + x245*L[102] + x247*L[145] + x249*L[101] + x25*x47 + x25*x82 + x25*L[50] + x250*L[144] + x252*L[143] + x253*L[96] + x255*L[138] + x257*L[137] + x258*L[132] + (1.0/4320.0)*x26*x33*L[213] + x27*x48 + x27*x83 + x27*L[77] + (1.0/10080.0)*x28*x31*L[212] + x29*x49 + x29*x84 + x29*L[112] + x30*x50 + x30*x85 + x30*L[156] + x32*x51 + x32*x72 + x32*x87 + x32*L[9] + x34*x52 + x34*x73 + x34*x89 + x34*L[19] + x36*x53 + x36*x74 + x36*x91 + x36*L[34] + x38*x54 + x38*x75 + x38*x93 + x38*L[55] + x40*x55 + x40*x76 + x40*x95 + x40*L[83] + x42*x56 + x42*x77 + x42*x97 + x42*L[119] + x43*x57 + x43*x78 + x43*L[164] + x58*x6 + x59*x6 + x6*L[4] + x60*x8 + x61*x8 + x8*L[10] + (1.0/362880.0)*pow(y, 9)*L[210] + y*L[2] + (1.0/362880.0)*pow(z, 9)*L[219] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += x*x111 + x*x58 + x*x59 + x*L[4] + x0 + x1 + x10*x117 + x10*x64 + x10*x65 + x10*L[35] + x100*x21 + x102*x23 + x104*x25 + x106*x27 + x108*x29 + x113*x6 + x115*x8 + x119*x12 + x12*x66 + x12*x67 + x12*L[56] + x121*x14 + x123*x227 + x123*L[38] + x125*x228 + x125*L[62] + x127*x229 + x127*L[94] + x129*x230 + x129*L[135] + x131*L[186] + x133*x223 + x133*L[40] + x134*x224 + x134*L[65] + x135*x225 + x135*L[98] + x136*x226 + x136*L[140] + x137*L[192] + x139*x234 + x139*L[59] + x14*x68 + x14*x69 + x14*L[84] + x141*x235 + x141*L[90] + x143*x236 + x143*L[130] + x145*L[180] + x147*x231 + x147*L[61] + x148*x232 + x148*L[93] + x149*x233 + x149*L[134] + x150*L[185] + x152*x239 + x152*L[87] + x154*x240 + x154*L[126] + x156*L[175] + x158*x237 + x158*L[89] + x159*x238 + x159*L[129] + x16*x70 + x16*x71 + x16*L[120] + x160*L[179] + x162*x242 + x162*L[123] + x164*L[171] + x166*x241 + x166*L[125] + x167*L[174] + x169*L[168] + x17*L[165] + x171*L[170] + x174*x282 + x174*L[47] + x176*x283 + x176*L[74] + x178*x284 + x178*L[109] + x180*x285 + x180*L[153] + x181*L[207] + x183*x286 + x183*L[73] + x185*x287 + x185*L[108] + x187*x288 + x187*L[152] + x188*L[206] + x19*x260 + x19*x277 + x19*x98 + x19*L[13] + x190*x289 + x190*L[107] + x192*x290 + x192*L[151] + x193*L[205] + x195*x291 + x195*L[150] + x196*L[204] + x197*L[203] + x21*x261 + x21*x278 + x21*L[26] + x23*x262 + x23*x279 + x23*L[45] + x243*L[96] + x245*L[138] + x247*L[190] + x249*L[137] + x25*x263 + x25*x280 + x25*L[71] + x250*L[189] + x252*L[188] + x253*L[132] + x255*L[183] + x257*L[182] + x258*L[177] + x264*x27 + x265*x29 + x266*x32 + x267*x34 + x268*x36 + x269*x38 + x27*x281 + x27*L[105] + x270*x40 + x271*x42 + x272*x32 + x273*x34 + x274*x36 + x275*x38 + x276*x40 + x29*L[148] + x30*L[201] + x32*x86 + x32*L[15] + x34*x88 + x34*L[29] + x36*x90 + x36*L[49] + x38*x92 + x38*L[76] + x4 + x40*x94 + x40*L[111] + x42*x96 + x42*L[155] + x43*L[209] + x6*x60 + x6*x61 + x6*L[10] + x62*x8 + x63*x8 + x8*L[20] + L[1];
#pragma omp atomic
Ls[2] += x*x292 + x*x293 + x*x3 + x*L[5] + x10*x114 + x10*x302 + x10*x324 + x10*L[36] + x101*x19 + x103*x21 + x105*x23 + x107*x25 + x109*x27 + x110*x6 + x112*x8 + x116*x12 + x118*x14 + x12*x303 + x12*x325 + x12*L[57] + x120*x16 + x123*x219 + x123*L[41] + x125*x220 + x125*L[66] + x127*x221 + x127*L[99] + x129*x222 + x129*L[141] + x131*L[193] + x133*x327 + x133*L[43] + x134*x328 + x134*L[69] + x135*x329 + x135*L[103] + x136*x330 + x136*L[146] + x137*L[199] + x139*x228 + x139*L[62] + x14*x304 + x14*x326 + x14*L[85] + x141*x229 + x141*L[94] + x143*x230 + x143*L[135] + x145*L[186] + x147*x331 + x147*L[64] + x148*x332 + x148*L[97] + x149*x333 + x149*L[139] + x150*L[191] + x152*x235 + x152*L[90] + x154*x236 + x154*L[130] + x156*L[180] + x158*x334 + x158*L[92] + x159*x335 + x159*L[133] + x16*x305 + x16*L[121] + x160*L[184] + x162*x240 + x162*L[126] + x164*L[175] + x166*x336 + x166*L[128] + x167*L[178] + x169*L[171] + x17*L[166] + x171*L[173] + x174*x203 + x174*L[52] + x176*x204 + x176*L[80] + x178*x205 + x178*L[116] + x180*x206 + x180*L[161] + x181*L[216] + x183*x207 + x183*L[79] + x185*x208 + x185*L[115] + x187*x209 + x187*L[160] + x188*L[215] + x19*x45 + x19*x80 + x19*L[16] + x190*x210 + x190*L[114] + x192*x211 + x192*L[159] + x193*L[214] + x195*x212 + x195*L[158] + x196*L[213] + x197*L[212] + x2 + x21*x46 + x21*x81 + x21*L[30] + x23*x47 + x23*x82 + x23*L[50] + x243*L[101] + x245*L[144] + x247*L[197] + x249*L[143] + x25*x48 + x25*x83 + x25*L[77] + x250*L[196] + x252*L[195] + x253*L[137] + x255*L[189] + x257*L[188] + x258*L[182] + x27*x49 + x27*x84 + x27*L[112] + x29*x50 + x29*x85 + x29*L[156] + x294*x32 + x295*x34 + x296*x36 + x297*x38 + x298*x40 + x299*x42 + x30*L[210] + x300*x6 + x301*x8 + x306*x32 + x307*x34 + x308*x36 + x309*x38 + x310*x40 + x311*x42 + x313*x32 + x315*x34 + x317*x36 + x319*x38 + x32*L[18] + x321*x40 + x322*x6 + x323*x8 + x34*L[33] + x36*L[54] + x38*L[82] + x40*L[118] + x42*L[163] + x43*L[218] + x6*L[11] + x79*y + x8*L[21] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*x337 + x*x338 + x*x341 + x*L[6] + x10*x352 + x10*x353 + x10*x381 + x10*L[37] + x12*x354 + x12*x355 + x12*x383 + x12*L[58] + x123*x386 + x123*L[42] + x125*x387 + x125*L[67] + x127*x388 + x127*L[100] + x129*x389 + x129*L[142] + x131*L[194] + x133*x214 + x133*L[44] + x134*x215 + x134*L[70] + x135*x216 + x135*L[104] + x136*x217 + x136*L[147] + x137*L[200] + x139*x390 + x139*L[63] + x14*x356 + x14*x357 + x14*x385 + x14*L[86] + x141*x391 + x141*L[95] + x143*x392 + x143*L[136] + x145*L[187] + x147*x224 + x147*L[65] + x148*x225 + x148*L[98] + x149*x226 + x149*L[140] + x150*L[192] + x152*x393 + x152*L[91] + x154*x394 + x154*L[131] + x156*L[181] + x158*x232 + x158*L[93] + x159*x233 + x159*L[134] + x16*x358 + x16*x359 + x16*L[122] + x160*L[185] + x162*x395 + x162*L[127] + x164*L[176] + x166*x238 + x166*L[129] + x167*L[179] + x169*L[172] + x17*L[167] + x171*L[174] + x174*x199 + x174*L[53] + x176*x200 + x176*L[81] + x178*x201 + x178*L[117] + x180*x202 + x180*L[162] + x181*L[217] + x183*x204 + x183*L[80] + x185*x205 + x185*L[116] + x187*x206 + x187*L[161] + x188*L[216] + x19*x342 + x19*x360 + x19*x367 + x19*L[17] + x190*x208 + x190*L[115] + x192*x209 + x192*L[160] + x193*L[215] + x195*x211 + x195*L[159] + x196*L[214] + x197*L[213] + x21*x343 + x21*x361 + x21*x369 + x21*L[31] + x23*x344 + x23*x362 + x23*x371 + x23*L[51] + x243*L[102] + x245*L[145] + x247*L[198] + x249*L[144] + x25*x345 + x25*x363 + x25*x373 + x25*L[78] + x250*L[197] + x252*L[196] + x253*L[138] + x255*L[190] + x257*L[189] + x258*L[183] + x27*x346 + x27*x364 + x27*x375 + x27*L[113] + x29*x347 + x29*x365 + x29*L[157] + x30*L[211] + x32*x52 + x32*x73 + x32*x89 + x32*L[19] + x339*y + x34*x53 + x34*x74 + x34*x91 + x34*L[34] + x348*x6 + x349*x6 + x350*x8 + x351*x8 + x36*x54 + x36*x75 + x36*x93 + x36*L[55] + x377*x6 + x379*x8 + x38*x55 + x38*x76 + x38*x95 + x38*L[83] + x40*x56 + x40*x77 + x40*x97 + x40*L[119] + x42*x57 + x42*x78 + x42*L[164] + x43*L[219] + x6*L[12] + x8*L[22] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += x*x113 + x*x60 + x*x61 + x*L[10] + x10*x119 + x10*x66 + x10*x67 + x10*L[56] + x111 + x115*x6 + x117*x8 + x12*x121 + x12*x68 + x12*x69 + x12*L[84] + x123*x234 + x123*L[59] + x125*x235 + x125*L[90] + x127*x236 + x127*L[130] + x129*L[180] + x133*x231 + x133*L[61] + x134*x232 + x134*L[93] + x135*x233 + x135*L[134] + x136*L[185] + x139*x239 + x139*L[87] + x14*x70 + x14*x71 + x14*L[120] + x141*x240 + x141*L[126] + x143*L[175] + x147*x237 + x147*L[89] + x148*x238 + x148*L[129] + x149*L[179] + x152*x242 + x152*L[123] + x154*L[171] + x158*x241 + x158*L[125] + x159*L[174] + x16*L[165] + x162*L[168] + x166*L[170] + x174*x414 + x174*L[68] + x176*x415 + x176*L[102] + x178*x416 + x178*L[145] + x180*L[198] + x183*x417 + x183*L[101] + x185*x418 + x185*L[144] + x187*L[197] + x19*x218 + x19*x396 + x19*x410 + x19*L[23] + x190*x419 + x190*L[143] + x192*L[196] + x195*L[195] + x21*x219 + x21*x397 + x21*x411 + x21*L[41] + x213*x32 + x214*x34 + x215*x36 + x216*x38 + x217*x40 + x220*x23 + x221*x25 + x222*x27 + x23*x398 + x23*x412 + x23*L[66] + x243*L[132] + x245*L[183] + x249*L[182] + x25*x399 + x25*x413 + x25*L[99] + x253*L[177] + x27*x400 + x27*L[141] + x29*L[193] + x32*x401 + x32*x406 + x32*L[25] + x34*x402 + x34*x407 + x34*L[44] + x36*x403 + x36*x408 + x36*L[70] + x38*x404 + x38*x409 + x38*L[104] + x40*x405 + x40*L[147] + x42*L[200] + x58 + x59 + x6*x62 + x6*x63 + x6*L[20] + x64*x8 + x65*x8 + x8*L[35] + L[4];
#pragma omp atomic
Ls[5] += x*x110 + x*x300 + x*x322 + x*L[11] + x10*x116 + x10*x303 + x10*x325 + x10*L[57] + x100*x19 + x102*x21 + x104*x23 + x106*x25 + x108*x27 + x112*x6 + x114*x8 + x118*x12 + x12*x304 + x12*x326 + x12*L[85] + x120*x14 + x123*x228 + x123*L[62] + x125*x229 + x125*L[94] + x127*x230 + x127*L[135] + x129*L[186] + x133*x331 + x133*L[64] + x134*x332 + x134*L[97] + x135*x333 + x135*L[139] + x136*L[191] + x139*x235 + x139*L[90] + x14*x305 + x14*L[121] + x141*x236 + x141*L[130] + x143*L[180] + x147*x334 + x147*L[92] + x148*x335 + x148*L[133] + x149*L[184] + x152*x240 + x152*L[126] + x154*L[175] + x158*x336 + x158*L[128] + x159*L[178] + x16*L[166] + x162*L[171] + x166*L[173] + x174*x286 + x174*L[73] + x176*x287 + x176*L[108] + x178*x288 + x178*L[152] + x180*L[206] + x183*x289 + x183*L[107] + x185*x290 + x185*L[151] + x187*L[205] + x19*x261 + x19*x278 + x19*L[26] + x190*x291 + x190*L[150] + x192*L[204] + x195*L[203] + x21*x262 + x21*x279 + x21*L[45] + x23*x263 + x23*x280 + x23*L[71] + x243*L[137] + x245*L[189] + x249*L[188] + x25*x264 + x25*x281 + x25*L[105] + x253*L[182] + x265*x27 + x27*L[148] + x29*L[201] + x292 + x293 + x3 + x301*x6 + x302*x8 + x312*x32 + x314*x34 + x316*x36 + x318*x38 + x32*x420 + x32*x425 + x32*L[28] + x320*x40 + x323*x6 + x324*x8 + x34*x421 + x34*x426 + x34*L[48] + x36*x422 + x36*x427 + x36*L[75] + x38*x423 + x38*x428 + x38*L[110] + x40*x424 + x40*L[154] + x42*L[208] + x6*L[21] + x8*L[36] + L[5];
#pragma omp atomic
Ls[6] += x*x348 + x*x349 + x*x377 + x*L[12] + x10*x354 + x10*x355 + x10*x383 + x10*L[58] + x12*x356 + x12*x357 + x12*x385 + x12*L[86] + x123*x390 + x123*L[63] + x125*x391 + x125*L[95] + x127*x392 + x127*L[136] + x129*L[187] + x133*x224 + x133*L[65] + x134*x225 + x134*L[98] + x135*x226 + x135*L[140] + x136*L[192] + x139*x393 + x139*L[91] + x14*x358 + x14*x359 + x14*L[122] + x141*x394 + x141*L[131] + x143*L[181] + x147*x232 + x147*L[93] + x148*x233 + x148*L[134] + x149*L[185] + x152*x395 + x152*L[127] + x154*L[176] + x158*x238 + x158*L[129] + x159*L[179] + x16*L[167] + x162*L[172] + x166*L[174] + x174*x283 + x174*L[74] + x176*x284 + x176*L[109] + x178*x285 + x178*L[153] + x180*L[207] + x183*x287 + x183*L[108] + x185*x288 + x185*L[152] + x187*L[206] + x19*x366 + x19*x429 + x19*x434 + x19*L[27] + x190*x290 + x190*L[151] + x192*L[205] + x195*L[204] + x21*x368 + x21*x430 + x21*x435 + x21*L[46] + x23*x370 + x23*x431 + x23*x436 + x23*L[72] + x243*L[138] + x245*L[190] + x249*L[189] + x25*x372 + x25*x432 + x25*x437 + x25*L[106] + x253*L[183] + x267*x32 + x268*x34 + x269*x36 + x27*x374 + x27*x433 + x27*L[149] + x270*x38 + x271*x40 + x273*x32 + x274*x34 + x275*x36 + x276*x38 + x29*L[202] + x32*x88 + x32*L[29] + x337 + x338 + x34*x90 + x34*L[49] + x341 + x350*x6 + x351*x6 + x352*x8 + x353*x8 + x36*x92 + x36*L[76] + x379*x6 + x38*x94 + x38*L[111] + x381*x8 + x40*x96 + x40*L[155] + x42*L[209] + x6*L[22] + x8*L[37] + L[6];
#pragma omp atomic
Ls[7] += x*x438 + x*x439 + x10*x234 + x10*x442 + x10*x460 + x10*L[59] + x103*x19 + x105*x21 + x107*x23 + x109*x25 + x12*x239 + x12*x443 + x12*x461 + x12*L[87] + x123*x220 + x123*L[66] + x125*x221 + x125*L[99] + x127*x222 + x127*L[141] + x129*L[193] + x133*x462 + x133*L[68] + x134*x463 + x134*L[102] + x135*x464 + x135*L[145] + x136*L[198] + x139*x229 + x139*L[94] + x14*x242 + x14*x444 + x14*L[123] + x141*x230 + x141*L[135] + x143*L[186] + x147*x465 + x147*L[96] + x148*x466 + x148*L[138] + x149*L[190] + x152*x236 + x152*L[130] + x154*L[180] + x158*x467 + x158*L[132] + x159*L[183] + x16*L[168] + x162*L[175] + x166*L[177] + x174*x207 + x174*L[79] + x176*x208 + x176*L[115] + x178*x209 + x178*L[160] + x180*L[215] + x183*x210 + x183*L[114] + x185*x211 + x185*L[159] + x187*L[214] + x19*x46 + x19*x81 + x19*L[30] + x190*x212 + x190*L[158] + x192*L[213] + x195*L[212] + x198*x32 + x199*x34 + x200*x36 + x201*x38 + x202*x40 + x21*x47 + x21*x82 + x21*L[50] + x218*x6 + x227*x8 + x23*x48 + x23*x83 + x23*L[77] + x243*L[143] + x245*L[196] + x249*L[195] + x25*x49 + x25*x84 + x25*L[112] + x253*L[188] + x27*x50 + x27*x85 + x27*L[156] + x29*L[210] + x32*x445 + x32*x451 + x32*L[32] + x34*x446 + x34*x453 + x34*L[53] + x36*x447 + x36*x455 + x36*L[81] + x38*x448 + x38*x457 + x38*L[117] + x40*x449 + x40*L[162] + x42*L[217] + x44 + x440*x6 + x441*x8 + x458*x6 + x459*x8 + x6*L[23] + x79 + x8*L[38] + x80*y + x99 + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += x*x340 + x*x468 + x*x469 + x*L[14] + x10*x380 + x10*x472 + x10*x477 + x10*L[60] + x12*x382 + x12*x473 + x12*x478 + x12*L[88] + x123*x387 + x123*L[67] + x125*x388 + x125*L[100] + x127*x389 + x127*L[142] + x129*L[194] + x133*x328 + x133*L[69] + x134*x329 + x134*L[103] + x135*x330 + x135*L[146] + x136*L[199] + x139*x391 + x139*L[95] + x14*x384 + x14*x474 + x14*L[124] + x141*x392 + x141*L[136] + x143*L[187] + x147*x332 + x147*L[97] + x148*x333 + x148*L[139] + x149*L[191] + x152*x394 + x152*L[131] + x154*L[181] + x158*x335 + x158*L[133] + x159*L[184] + x16*L[169] + x162*L[176] + x166*L[178] + x174*x204 + x174*L[80] + x176*x205 + x176*L[116] + x178*x206 + x178*L[161] + x180*L[216] + x183*x208 + x183*L[115] + x185*x209 + x185*L[160] + x187*L[215] + x19*x343 + x19*x361 + x19*x369 + x19*L[31] + x190*x211 + x190*L[159] + x192*L[214] + x195*L[213] + x21*x344 + x21*x362 + x21*x371 + x21*L[51] + x23*x345 + x23*x363 + x23*x373 + x23*L[78] + x243*L[144] + x245*L[197] + x249*L[196] + x25*x346 + x25*x364 + x25*x375 + x25*L[113] + x253*L[189] + x27*x347 + x27*x365 + x27*L[157] + x29*L[211] + x295*x32 + x296*x34 + x297*x36 + x298*x38 + x299*x40 + x307*x32 + x308*x34 + x309*x36 + x310*x38 + x311*x40 + x315*x32 + x317*x34 + x319*x36 + x32*L[33] + x321*x38 + x339 + x34*L[54] + x36*L[82] + x360*y + x376*x6 + x378*x8 + x38*L[118] + x40*L[163] + x42*L[218] + x470*x6 + x471*x8 + x475*x6 + x476*x8 + x6*L[24] + x8*L[39] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += x*x479 + x*x482 + x10*x231 + x10*x485 + x10*x506 + x10*L[61] + x12*x237 + x12*x486 + x12*x508 + x12*L[89] + x123*x509 + x123*L[68] + x125*x510 + x125*L[101] + x127*x511 + x127*L[143] + x129*L[195] + x133*x215 + x133*L[70] + x134*x216 + x134*L[104] + x135*x217 + x135*L[147] + x136*L[200] + x139*x512 + x139*L[96] + x14*x241 + x14*x487 + x14*L[125] + x141*x513 + x141*L[137] + x143*L[188] + x147*x225 + x147*L[98] + x148*x226 + x148*L[140] + x149*L[192] + x152*x514 + x152*L[132] + x154*L[182] + x158*x233 + x158*L[134] + x159*L[185] + x16*L[170] + x162*L[177] + x166*L[179] + x174*x200 + x174*L[81] + x176*x201 + x176*L[117] + x178*x202 + x178*L[162] + x180*L[217] + x183*x205 + x183*L[116] + x185*x206 + x185*L[161] + x187*L[216] + x19*x198 + x19*x488 + x19*x494 + x19*L[32] + x190*x209 + x190*L[160] + x192*L[215] + x195*L[214] + x203*x21 + x207*x23 + x21*x489 + x21*x496 + x21*L[52] + x210*x25 + x212*x27 + x213*x6 + x223*x8 + x23*x490 + x23*x498 + x23*L[79] + x243*L[145] + x245*L[198] + x249*L[197] + x25*x491 + x25*x500 + x25*L[114] + x253*L[190] + x27*x492 + x27*L[158] + x29*L[212] + x32*x53 + x32*x74 + x32*x91 + x32*L[34] + x34*x54 + x34*x75 + x34*x93 + x34*L[55] + x36*x55 + x36*x76 + x36*x95 + x36*L[83] + x38*x56 + x38*x77 + x38*x97 + x38*L[119] + x40*x57 + x40*x78 + x40*L[164] + x42*L[219] + x480*y + x483*x6 + x484*x8 + x502*x6 + x504*x8 + x51 + x6*L[25] + x72 + x8*L[40] + x87 + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += x*x115 + x*x62 + x*x63 + x*L[20] + x10*x121 + x10*x68 + x10*x69 + x10*L[84] + x113 + x117*x6 + x119*x8 + x12*x70 + x12*x71 + x12*L[120] + x123*x239 + x123*L[87] + x125*x240 + x125*L[126] + x127*L[175] + x133*x237 + x133*L[89] + x134*x238 + x134*L[129] + x135*L[179] + x139*x242 + x139*L[123] + x14*L[165] + x141*L[171] + x147*x241 + x147*L[125] + x148*L[174] + x152*L[168] + x158*L[170] + x174*x529 + x174*L[96] + x176*x530 + x176*L[138] + x178*L[190] + x183*x531 + x183*L[137] + x185*L[189] + x19*x227 + x19*x515 + x19*x526 + x19*L[38] + x190*L[188] + x21*x228 + x21*x516 + x21*x527 + x21*L[62] + x223*x32 + x224*x34 + x225*x36 + x226*x38 + x229*x23 + x23*x517 + x23*x528 + x23*L[94] + x230*x25 + x243*L[177] + x25*x518 + x25*L[135] + x27*L[186] + x32*x519 + x32*x523 + x32*L[40] + x34*x520 + x34*x524 + x34*L[65] + x36*x521 + x36*x525 + x36*L[98] + x38*x522 + x38*L[140] + x40*L[192] + x6*x64 + x6*x65 + x6*L[35] + x60 + x61 + x66*x8 + x67*x8 + x8*L[56] + L[10];
#pragma omp atomic
Ls[11] += x*x112 + x*x301 + x*x323 + x*L[21] + x10*x118 + x10*x304 + x10*x326 + x10*L[85] + x110 + x114*x6 + x116*x8 + x12*x120 + x12*x305 + x12*L[121] + x123*x235 + x123*L[90] + x125*x236 + x125*L[130] + x127*L[180] + x133*x334 + x133*L[92] + x134*x335 + x134*L[133] + x135*L[184] + x139*x240 + x139*L[126] + x14*L[166] + x141*L[175] + x147*x336 + x147*L[128] + x148*L[178] + x152*L[171] + x158*L[173] + x174*x417 + x174*L[101] + x176*x418 + x176*L[144] + x178*L[197] + x183*x419 + x183*L[143] + x185*L[196] + x19*x219 + x19*x397 + x19*x411 + x19*L[41] + x190*L[195] + x21*x220 + x21*x398 + x21*x412 + x21*L[66] + x221*x23 + x222*x25 + x23*x399 + x23*x413 + x23*L[99] + x243*L[182] + x25*x400 + x25*L[141] + x27*L[193] + x300 + x302*x6 + x303*x8 + x32*x327 + x32*x532 + x32*x536 + x32*L[43] + x322 + x324*x6 + x325*x8 + x328*x34 + x329*x36 + x330*x38 + x34*x533 + x34*x537 + x34*L[69] + x36*x534 + x36*x538 + x36*L[103] + x38*x535 + x38*L[146] + x40*L[199] + x6*L[36] + x8*L[57] + L[11];
#pragma omp atomic
Ls[12] += x*x350 + x*x351 + x*x379 + x*L[22] + x10*x356 + x10*x357 + x10*x385 + x10*L[86] + x12*x358 + x12*x359 + x12*L[122] + x123*x393 + x123*L[91] + x125*x394 + x125*L[131] + x127*L[181] + x133*x232 + x133*L[93] + x134*x233 + x134*L[134] + x135*L[185] + x139*x395 + x139*L[127] + x14*L[167] + x141*L[176] + x147*x238 + x147*L[129] + x148*L[179] + x152*L[172] + x158*L[174] + x174*x415 + x174*L[102] + x176*x416 + x176*L[145] + x178*L[198] + x183*x418 + x183*L[144] + x185*L[197] + x19*x386 + x19*x539 + x19*x543 + x19*L[42] + x190*L[196] + x21*x387 + x21*x540 + x21*x544 + x21*L[67] + x214*x32 + x215*x34 + x216*x36 + x217*x38 + x23*x388 + x23*x541 + x23*x545 + x23*L[100] + x243*L[183] + x25*x389 + x25*x542 + x25*L[142] + x27*L[194] + x32*x402 + x32*x407 + x32*L[44] + x34*x403 + x34*x408 + x34*L[70] + x348 + x349 + x352*x6 + x353*x6 + x354*x8 + x355*x8 + x36*x404 + x36*x409 + x36*L[104] + x377 + x38*x405 + x38*L[147] + x381*x6 + x383*x8 + x40*L[200] + x6*L[37] + x8*L[58] + L[12];
#pragma omp atomic
Ls[13] += x*x440 + x*x458 + x10*x239 + x10*x443 + x10*x461 + x10*L[87] + x102*x19 + x104*x21 + x106*x23 + x108*x25 + x12*x242 + x12*x444 + x12*L[123] + x123*x229 + x123*L[94] + x125*x230 + x125*L[135] + x127*L[186] + x133*x465 + x133*L[96] + x134*x466 + x134*L[138] + x135*L[190] + x139*x236 + x139*L[130] + x14*L[168] + x141*L[180] + x147*x467 + x147*L[132] + x148*L[183] + x152*L[175] + x158*L[177] + x174*x289 + x174*L[107] + x176*x290 + x176*L[151] + x178*L[205] + x183*x291 + x183*L[150] + x185*L[204] + x19*x262 + x19*x279 + x19*L[45] + x190*L[203] + x21*x263 + x21*x280 + x21*L[71] + x227*x6 + x23*x264 + x23*x281 + x23*L[105] + x234*x8 + x243*L[188] + x25*x265 + x25*L[148] + x260 + x27*L[201] + x277 + x282*x32 + x283*x34 + x284*x36 + x285*x38 + x32*x450 + x32*x546 + x32*L[47] + x34*x452 + x34*x547 + x34*L[74] + x36*x454 + x36*x548 + x36*L[109] + x38*x456 + x38*L[153] + x40*L[207] + x438 + x439 + x441*x6 + x442*x8 + x459*x6 + x460*x8 + x6*L[38] + x8*L[59] + x98 + L[13];
#pragma omp atomic
Ls[14] += x*x376 + x*x470 + x*x475 + x*L[24] + x10*x382 + x10*x473 + x10*x478 + x10*L[88] + x12*x384 + x12*x474 + x12*L[124] + x123*x391 + x123*L[95] + x125*x392 + x125*L[136] + x127*L[187] + x133*x332 + x133*L[97] + x134*x333 + x134*L[139] + x135*L[191] + x139*x394 + x139*L[131] + x14*L[169] + x141*L[181] + x147*x335 + x147*L[133] + x148*L[184] + x152*L[176] + x158*L[178] + x174*x287 + x174*L[108] + x176*x288 + x176*L[152] + x178*L[206] + x183*x290 + x183*L[151] + x185*L[205] + x19*x368 + x19*x430 + x19*x435 + x19*L[46] + x190*L[204] + x21*x370 + x21*x431 + x21*x436 + x21*L[72] + x23*x372 + x23*x432 + x23*x437 + x23*L[106] + x243*L[189] + x25*x374 + x25*x433 + x25*L[149] + x27*L[202] + x314*x32 + x316*x34 + x318*x36 + x32*x421 + x32*x426 + x32*L[48] + x320*x38 + x34*x422 + x34*x427 + x34*L[75] + x340 + x36*x423 + x36*x428 + x36*L[110] + x378*x6 + x38*x424 + x38*L[154] + x380*x8 + x40*L[208] + x468 + x469 + x471*x6 + x472*x8 + x476*x6 + x477*x8 + x6*L[39] + x8*L[60] + L[14];
#pragma omp atomic
Ls[15] += x*x483 + x*x502 + x10*x237 + x10*x486 + x10*x508 + x10*L[89] + x12*x241 + x12*x487 + x12*L[125] + x123*x512 + x123*L[96] + x125*x513 + x125*L[137] + x127*L[188] + x133*x225 + x133*L[98] + x134*x226 + x134*L[140] + x135*L[192] + x139*x514 + x139*L[132] + x14*L[170] + x141*L[182] + x147*x233 + x147*L[134] + x148*L[185] + x152*L[177] + x158*L[179] + x174*x284 + x174*L[109] + x176*x285 + x176*L[153] + x178*L[207] + x183*x288 + x183*L[152] + x185*L[206] + x19*x282 + x19*x493 + x19*x549 + x19*L[47] + x190*L[205] + x21*x286 + x21*x495 + x21*x550 + x21*L[73] + x223*x6 + x23*x289 + x23*x497 + x23*x551 + x23*L[107] + x231*x8 + x243*L[190] + x25*x291 + x25*x499 + x25*L[150] + x266 + x268*x32 + x269*x34 + x27*L[203] + x270*x36 + x271*x38 + x272 + x274*x32 + x275*x34 + x276*x36 + x32*x90 + x32*L[49] + x34*x92 + x34*L[76] + x36*x94 + x36*L[111] + x38*x96 + x38*L[155] + x40*L[209] + x479 + x482 + x484*x6 + x485*x8 + x504*x6 + x506*x8 + x6*L[40] + x8*L[61] + x86 + L[15];
#pragma omp atomic
Ls[16] += x*x552 + x*x553 + x10*x235 + x10*x556 + x10*x570 + x10*L[90] + x101 + x105*x19 + x107*x21 + x109*x23 + x12*x240 + x12*x557 + x12*L[126] + x123*x221 + x123*L[99] + x125*x222 + x125*L[141] + x127*L[193] + x133*x571 + x133*L[101] + x134*x572 + x134*L[144] + x135*L[197] + x139*x230 + x139*L[135] + x14*L[171] + x141*L[186] + x147*x573 + x147*L[137] + x148*L[189] + x152*L[180] + x158*L[182] + x174*x210 + x174*L[114] + x176*x211 + x176*L[159] + x178*L[214] + x183*x212 + x183*L[158] + x185*L[213] + x19*x47 + x19*x82 + x19*L[50] + x190*L[212] + x203*x32 + x204*x34 + x205*x36 + x206*x38 + x21*x48 + x21*x83 + x21*L[77] + x219*x6 + x228*x8 + x23*x49 + x23*x84 + x23*L[112] + x243*L[195] + x25*x50 + x25*x85 + x25*L[156] + x27*L[210] + x32*x558 + x32*x563 + x32*L[52] + x34*x559 + x34*x565 + x34*L[80] + x36*x560 + x36*x567 + x36*L[116] + x38*x561 + x38*L[161] + x40*L[216] + x45 + x554*x6 + x555*x8 + x568*x6 + x569*x8 + x6*L[41] + x8*L[62] + x80 + x81*y + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += x*x574 + x*x575 + x10*x393 + x10*x578 + x10*x582 + x10*L[91] + x12*x395 + x12*x579 + x12*L[127] + x123*x388 + x123*L[100] + x125*x389 + x125*L[142] + x127*L[194] + x133*x463 + x133*L[102] + x134*x464 + x134*L[145] + x135*L[198] + x139*x392 + x139*L[136] + x14*L[172] + x141*L[187] + x147*x466 + x147*L[138] + x148*L[190] + x152*L[181] + x158*L[183] + x174*x208 + x174*L[115] + x176*x209 + x176*L[160] + x178*L[215] + x183*x211 + x183*L[159] + x185*L[214] + x19*x344 + x19*x362 + x19*x371 + x19*L[51] + x190*L[213] + x199*x32 + x200*x34 + x201*x36 + x202*x38 + x21*x345 + x21*x363 + x21*x373 + x21*L[78] + x23*x346 + x23*x364 + x23*x375 + x23*L[113] + x243*L[196] + x25*x347 + x25*x365 + x25*L[157] + x27*L[211] + x32*x446 + x32*x453 + x32*L[53] + x34*x447 + x34*x455 + x34*L[81] + x342 + x36*x448 + x36*x457 + x36*L[117] + x360 + x361*y + x367 + x38*x449 + x38*L[162] + x386*x6 + x390*x8 + x40*L[217] + x576*x6 + x577*x8 + x580*x6 + x581*x8 + x6*L[42] + x8*L[63] + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += x*x481 + x*x583 + x10*x334 + x10*x505 + x10*x586 + x10*L[92] + x12*x336 + x12*x507 + x12*L[128] + x123*x510 + x123*L[101] + x125*x511 + x125*L[143] + x127*L[195] + x133*x329 + x133*L[103] + x134*x330 + x134*L[146] + x135*L[199] + x139*x513 + x139*L[137] + x14*L[173] + x141*L[188] + x147*x333 + x147*L[139] + x148*L[191] + x152*L[182] + x158*L[184] + x174*x205 + x174*L[116] + x176*x206 + x176*L[161] + x178*L[216] + x183*x209 + x183*L[160] + x185*L[215] + x19*x203 + x19*x489 + x19*x496 + x19*L[52] + x190*L[214] + x207*x21 + x21*x490 + x21*x498 + x21*L[79] + x210*x23 + x212*x25 + x23*x491 + x23*x500 + x23*L[114] + x243*L[197] + x25*x492 + x25*L[158] + x27*L[212] + x294 + x296*x32 + x297*x34 + x298*x36 + x299*x38 + x306 + x308*x32 + x309*x34 + x310*x36 + x311*x38 + x313 + x317*x32 + x319*x34 + x32*L[54] + x321*x36 + x327*x6 + x331*x8 + x34*L[82] + x36*L[118] + x38*L[163] + x40*L[218] + x480 + x488*y + x501*x6 + x503*x8 + x584*x6 + x585*x8 + x6*L[43] + x8*L[64] + L[18];
#pragma omp atomic
Ls[19] += x*x587 + x*x590 + x10*x232 + x10*x593 + x10*x610 + x10*L[93] + x12*x238 + x12*x594 + x12*L[129] + x123*x611 + x123*L[102] + x125*x612 + x125*L[144] + x127*L[196] + x133*x216 + x133*L[104] + x134*x217 + x134*L[147] + x135*L[200] + x139*x613 + x139*L[138] + x14*L[174] + x141*L[189] + x147*x226 + x147*L[140] + x148*L[192] + x152*L[183] + x158*L[185] + x174*x201 + x174*L[117] + x176*x202 + x176*L[162] + x178*L[217] + x183*x206 + x183*L[161] + x185*L[216] + x19*x199 + x19*x595 + x19*x600 + x19*L[53] + x190*L[215] + x204*x21 + x208*x23 + x21*x596 + x21*x602 + x21*L[80] + x211*x25 + x214*x6 + x224*x8 + x23*x597 + x23*x604 + x23*L[115] + x243*L[198] + x25*x598 + x25*L[159] + x27*L[213] + x32*x54 + x32*x75 + x32*x93 + x32*L[55] + x34*x55 + x34*x76 + x34*x95 + x34*L[83] + x36*x56 + x36*x77 + x36*x97 + x36*L[119] + x38*x57 + x38*x78 + x38*L[164] + x40*L[219] + x52 + x588*y + x591*x6 + x592*x8 + x6*x606 + x6*L[44] + x608*x8 + x73 + x8*L[65] + x89 + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += x*x117 + x*x64 + x*x65 + x*L[35] + x10*x70 + x10*x71 + x10*L[120] + x115 + x119*x6 + x12*L[165] + x121*x8 + x123*x242 + x123*L[123] + x125*L[171] + x133*x241 + x133*L[125] + x134*L[174] + x139*L[168] + x147*L[170] + x174*x624 + x174*L[132] + x176*L[183] + x183*L[182] + x19*x234 + x19*x614 + x19*x622 + x19*L[59] + x21*x235 + x21*x615 + x21*x623 + x21*L[90] + x23*x236 + x23*x616 + x23*L[130] + x231*x32 + x232*x34 + x233*x36 + x25*L[180] + x32*x617 + x32*x620 + x32*L[61] + x34*x618 + x34*x621 + x34*L[93] + x36*x619 + x36*L[134] + x38*L[185] + x6*x66 + x6*x67 + x6*L[56] + x62 + x63 + x68*x8 + x69*x8 + x8*L[84] + L[20];
#pragma omp atomic
Ls[21] += x*x114 + x*x302 + x*x324 + x*L[36] + x10*x120 + x10*x305 + x10*L[121] + x112 + x116*x6 + x118*x8 + x12*L[166] + x123*x240 + x123*L[126] + x125*L[175] + x133*x336 + x133*L[128] + x134*L[178] + x139*L[171] + x147*L[173] + x174*x531 + x174*L[137] + x176*L[189] + x183*L[188] + x19*x228 + x19*x516 + x19*x527 + x19*L[62] + x21*x229 + x21*x517 + x21*x528 + x21*L[94] + x23*x230 + x23*x518 + x23*L[135] + x25*L[186] + x301 + x303*x6 + x304*x8 + x32*x331 + x32*x625 + x32*x628 + x32*L[64] + x323 + x325*x6 + x326*x8 + x332*x34 + x333*x36 + x34*x626 + x34*x629 + x34*L[97] + x36*x627 + x36*L[139] + x38*L[191] + x6*L[57] + x8*L[85] + L[21];
#pragma omp atomic
Ls[22] += x*x352 + x*x353 + x*x381 + x*L[37] + x10*x358 + x10*x359 + x10*L[122] + x12*L[167] + x123*x395 + x123*L[127] + x125*L[176] + x133*x238 + x133*L[129] + x134*L[179] + x139*L[172] + x147*L[174] + x174*x530 + x174*L[138] + x176*L[190] + x183*L[189] + x19*x390 + x19*x630 + x19*x633 + x19*L[63] + x21*x391 + x21*x631 + x21*x634 + x21*L[95] + x224*x32 + x225*x34 + x226*x36 + x23*x392 + x23*x632 + x23*L[136] + x25*L[187] + x32*x520 + x32*x524 + x32*L[65] + x34*x521 + x34*x525 + x34*L[98] + x350 + x351 + x354*x6 + x355*x6 + x356*x8 + x357*x8 + x36*x522 + x36*L[140] + x379 + x38*L[192] + x383*x6 + x385*x8 + x6*L[58] + x8*L[86] + L[22];
#pragma omp atomic
Ls[23] += x*x441 + x*x459 + x10*x242 + x10*x444 + x10*L[123] + x12*L[168] + x123*x236 + x123*L[130] + x125*L[180] + x133*x467 + x133*L[132] + x134*L[183] + x139*L[175] + x147*L[177] + x174*x419 + x174*L[143] + x176*L[196] + x183*L[195] + x19*x220 + x19*x398 + x19*x412 + x19*L[66] + x21*x221 + x21*x399 + x21*x413 + x21*L[99] + x218 + x222*x23 + x23*x400 + x23*L[141] + x234*x6 + x239*x8 + x25*L[193] + x32*x414 + x32*x462 + x32*x635 + x32*L[68] + x34*x415 + x34*x463 + x34*x636 + x34*L[102] + x36*x416 + x36*x464 + x36*L[145] + x38*L[198] + x396 + x410 + x440 + x442*x6 + x443*x8 + x458 + x460*x6 + x461*x8 + x6*L[59] + x8*L[87] + L[23];
#pragma omp atomic
Ls[24] += x*x378 + x*x471 + x*x476 + x*L[39] + x10*x384 + x10*x474 + x10*L[124] + x12*L[169] + x123*x394 + x123*L[131] + x125*L[181] + x133*x335 + x133*L[133] + x134*L[184] + x139*L[176] + x147*L[178] + x174*x418 + x174*L[144] + x176*L[197] + x183*L[196] + x19*x387 + x19*x540 + x19*x544 + x19*L[67] + x21*x388 + x21*x541 + x21*x545 + x21*L[100] + x23*x389 + x23*x542 + x23*L[142] + x25*L[194] + x32*x328 + x32*x533 + x32*x537 + x32*L[69] + x329*x34 + x330*x36 + x34*x534 + x34*x538 + x34*L[103] + x36*x535 + x36*L[146] + x376 + x38*L[199] + x380*x6 + x382*x8 + x470 + x472*x6 + x473*x8 + x475 + x477*x6 + x478*x8 + x6*L[60] + x8*L[88] + L[24];
#pragma omp atomic
Ls[25] += x*x484 + x*x504 + x10*x241 + x10*x487 + x10*L[125] + x12*L[170] + x123*x514 + x123*L[132] + x125*L[182] + x133*x233 + x133*L[134] + x134*L[185] + x139*L[177] + x147*L[179] + x174*x416 + x174*L[145] + x176*L[198] + x183*L[197] + x19*x414 + x19*x509 + x19*x637 + x19*L[68] + x21*x417 + x21*x510 + x21*x638 + x21*L[101] + x213 + x215*x32 + x216*x34 + x217*x36 + x23*x419 + x23*x511 + x23*L[143] + x231*x6 + x237*x8 + x25*L[195] + x32*x403 + x32*x408 + x32*L[70] + x34*x404 + x34*x409 + x34*L[104] + x36*x405 + x36*L[147] + x38*L[200] + x401 + x406 + x483 + x485*x6 + x486*x8 + x502 + x506*x6 + x508*x8 + x6*L[61] + x8*L[89] + L[25];
#pragma omp atomic
Ls[26] += x*x554 + x*x568 + x10*x240 + x10*x557 + x10*L[126] + x100 + x104*x19 + x106*x21 + x108*x23 + x12*L[171] + x123*x230 + x123*L[135] + x125*L[186] + x133*x573 + x133*L[137] + x134*L[189] + x139*L[180] + x147*L[182] + x174*x291 + x174*L[150] + x176*L[204] + x183*L[203] + x19*x263 + x19*x280 + x19*L[71] + x21*x264 + x21*x281 + x21*L[105] + x228*x6 + x23*x265 + x23*L[148] + x235*x8 + x25*L[201] + x261 + x278 + x286*x32 + x287*x34 + x288*x36 + x32*x562 + x32*x639 + x32*L[73] + x34*x564 + x34*x640 + x34*L[108] + x36*x566 + x36*L[152] + x38*L[206] + x552 + x553 + x555*x6 + x556*x8 + x569*x6 + x570*x8 + x6*L[62] + x8*L[90] + L[26];
#pragma omp atomic
Ls[27] += x*x576 + x*x580 + x10*x395 + x10*x579 + x10*L[127] + x12*L[172] + x123*x392 + x123*L[136] + x125*L[187] + x133*x466 + x133*L[138] + x134*L[190] + x139*L[181] + x147*L[183] + x174*x290 + x174*L[151] + x176*L[205] + x183*L[204] + x19*x370 + x19*x431 + x19*x436 + x19*L[72] + x21*x372 + x21*x432 + x21*x437 + x21*L[106] + x23*x374 + x23*x433 + x23*L[149] + x25*L[202] + x283*x32 + x284*x34 + x285*x36 + x32*x452 + x32*x547 + x32*L[74] + x34*x454 + x34*x548 + x34*L[109] + x36*x456 + x36*L[153] + x366 + x38*L[207] + x390*x6 + x393*x8 + x429 + x434 + x574 + x575 + x577*x6 + x578*x8 + x581*x6 + x582*x8 + x6*L[63] + x8*L[91] + L[27];
#pragma omp atomic
Ls[28] += x*x501 + x*x584 + x10*x336 + x10*x507 + x10*L[128] + x12*L[173] + x123*x513 + x123*L[137] + x125*L[188] + x133*x333 + x133*L[139] + x134*L[191] + x139*L[182] + x147*L[184] + x174*x288 + x174*L[152] + x176*L[206] + x183*L[205] + x19*x286 + x19*x495 + x19*x550 + x19*L[73] + x21*x289 + x21*x497 + x21*x551 + x21*L[107] + x23*x291 + x23*x499 + x23*L[150] + x25*L[203] + x312 + x316*x32 + x318*x34 + x32*x422 + x32*x427 + x32*L[75] + x320*x36 + x331*x6 + x334*x8 + x34*x423 + x34*x428 + x34*L[110] + x36*x424 + x36*L[154] + x38*L[208] + x420 + x425 + x481 + x503*x6 + x505*x8 + x583 + x585*x6 + x586*x8 + x6*L[64] + x8*L[92] + L[28];
#pragma omp atomic
Ls[29] += x*x591 + x*x606 + x10*x238 + x10*x594 + x10*L[129] + x12*L[174] + x123*x613 + x123*L[138] + x125*L[189] + x133*x226 + x133*L[140] + x134*L[192] + x139*L[183] + x147*L[185] + x174*x285 + x174*L[153] + x176*L[207] + x183*L[206] + x19*x283 + x19*x599 + x19*x641 + x19*L[74] + x21*x287 + x21*x601 + x21*x642 + x21*L[108] + x224*x6 + x23*x290 + x23*x603 + x23*L[151] + x232*x8 + x25*L[204] + x267 + x269*x32 + x270*x34 + x271*x36 + x273 + x275*x32 + x276*x34 + x32*x92 + x32*L[76] + x34*x94 + x34*L[111] + x36*x96 + x36*L[155] + x38*L[209] + x587 + x590 + x592*x6 + x593*x8 + x6*x608 + x6*L[65] + x610*x8 + x8*L[93] + x88 + L[29];
#pragma omp atomic
Ls[30] += x*x643 + x*x644 + x10*x236 + x10*x647 + x10*L[130] + x103 + x107*x19 + x109*x21 + x12*L[175] + x123*x222 + x123*L[141] + x125*L[193] + x133*x657 + x133*L[143] + x134*L[196] + x139*L[186] + x147*L[188] + x174*x212 + x174*L[158] + x176*L[213] + x183*L[212] + x19*x48 + x19*x83 + x19*L[77] + x207*x32 + x208*x34 + x209*x36 + x21*x49 + x21*x84 + x21*L[112] + x220*x6 + x229*x8 + x23*x50 + x23*x85 + x23*L[156] + x25*L[210] + x32*x648 + x32*x652 + x32*L[79] + x34*x649 + x34*x654 + x34*L[115] + x36*x650 + x36*L[160] + x38*L[215] + x46 + x6*x645 + x6*x655 + x6*L[66] + x646*x8 + x656*x8 + x8*L[94] + x81 + x82*y + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += x*x658 + x*x659 + x10*x394 + x10*x662 + x10*L[131] + x12*L[176] + x123*x389 + x123*L[142] + x125*L[194] + x133*x572 + x133*L[144] + x134*L[197] + x139*L[187] + x147*L[189] + x174*x211 + x174*L[159] + x176*L[214] + x183*L[213] + x19*x345 + x19*x363 + x19*x373 + x19*L[78] + x204*x32 + x205*x34 + x206*x36 + x21*x346 + x21*x364 + x21*x375 + x21*L[113] + x23*x347 + x23*x365 + x23*L[157] + x25*L[211] + x32*x559 + x32*x565 + x32*L[80] + x34*x560 + x34*x567 + x34*L[116] + x343 + x36*x561 + x36*L[161] + x361 + x362*y + x369 + x38*L[216] + x387*x6 + x391*x8 + x6*x660 + x6*x663 + x6*L[67] + x661*x8 + x664*x8 + x8*L[95] + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += x*x665 + x10*x467 + x10*x514 + x10*L[132] + x12*L[177] + x123*x511 + x123*L[143] + x125*L[195] + x133*x464 + x133*L[145] + x134*L[198] + x139*L[188] + x147*L[190] + x174*x209 + x174*L[160] + x176*L[215] + x183*L[214] + x19*x207 + x19*x490 + x19*x498 + x19*L[79] + x198 + x200*x32 + x201*x34 + x202*x36 + x21*x210 + x21*x491 + x21*x500 + x21*L[114] + x212*x23 + x23*x492 + x23*L[158] + x25*L[212] + x32*x447 + x32*x455 + x32*L[81] + x34*x448 + x34*x457 + x34*L[117] + x36*x449 + x36*L[162] + x38*L[217] + x445 + x451 + x462*x6 + x465*x8 + x488 + x489*y + x494 + x509*x6 + x512*x8 + x6*x666 + x6*L[68] + x667*x8 + x8*L[96] + L[32];
#pragma omp atomic
Ls[33] += x*x589 + x*x668 + x10*x335 + x10*x609 + x10*L[133] + x12*L[178] + x123*x612 + x123*L[144] + x125*L[196] + x133*x330 + x133*L[146] + x134*L[199] + x139*L[189] + x147*L[191] + x174*x206 + x174*L[161] + x176*L[216] + x183*L[215] + x19*x204 + x19*x596 + x19*x602 + x19*L[80] + x208*x21 + x21*x597 + x21*x604 + x21*L[115] + x211*x23 + x23*x598 + x23*L[159] + x25*L[213] + x295 + x297*x32 + x298*x34 + x299*x36 + x307 + x309*x32 + x310*x34 + x311*x36 + x315 + x319*x32 + x32*L[82] + x321*x34 + x328*x6 + x332*x8 + x34*L[118] + x36*L[163] + x38*L[218] + x588 + x595*y + x6*x605 + x6*x669 + x6*L[69] + x607*x8 + x670*x8 + x8*L[97] + L[33];
#pragma omp atomic
Ls[34] += x*x671 + x*x674 + x10*x233 + x10*x677 + x10*L[134] + x12*L[179] + x123*x689 + x123*L[145] + x125*L[197] + x133*x217 + x133*L[147] + x134*L[200] + x139*L[190] + x147*L[192] + x174*x202 + x174*L[162] + x176*L[217] + x183*L[216] + x19*x200 + x19*x678 + x19*x682 + x19*L[81] + x205*x21 + x209*x23 + x21*x679 + x21*x684 + x21*L[116] + x215*x6 + x225*x8 + x23*x680 + x23*L[160] + x25*L[214] + x32*x55 + x32*x76 + x32*x95 + x32*L[83] + x34*x56 + x34*x77 + x34*x97 + x34*L[119] + x36*x57 + x36*x78 + x36*L[164] + x38*L[219] + x53 + x6*x675 + x6*x686 + x6*L[70] + x672*y + x676*x8 + x688*x8 + x74 + x8*L[98] + x91 + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += x*x119 + x*x66 + x*x67 + x*L[56] + x10*L[165] + x117 + x121*x6 + x123*L[168] + x133*L[170] + x174*L[177] + x19*x239 + x19*x690 + x19*x695 + x19*L[87] + x21*x240 + x21*x691 + x21*L[126] + x23*L[175] + x237*x32 + x238*x34 + x32*x692 + x32*x694 + x32*L[89] + x34*x693 + x34*L[129] + x36*L[179] + x6*x68 + x6*x69 + x6*L[84] + x64 + x65 + x70*x8 + x71*x8 + x8*L[120] + L[35];
#pragma omp atomic
Ls[36] += x*x116 + x*x303 + x*x325 + x*L[57] + x10*L[166] + x114 + x118*x6 + x120*x8 + x123*L[171] + x133*L[173] + x174*L[182] + x19*x235 + x19*x615 + x19*x623 + x19*L[90] + x21*x236 + x21*x616 + x21*L[130] + x23*L[180] + x302 + x304*x6 + x305*x8 + x32*x334 + x32*x696 + x32*x698 + x32*L[92] + x324 + x326*x6 + x335*x34 + x34*x697 + x34*L[133] + x36*L[184] + x6*L[85] + x8*L[121] + L[36];
#pragma omp atomic
Ls[37] += x*x354 + x*x355 + x*x383 + x*L[58] + x10*L[167] + x123*L[172] + x133*L[174] + x174*L[183] + x19*x393 + x19*x699 + x19*x701 + x19*L[91] + x21*x394 + x21*x700 + x21*L[131] + x23*L[181] + x232*x32 + x233*x34 + x32*x618 + x32*x621 + x32*L[93] + x34*x619 + x34*L[134] + x352 + x353 + x356*x6 + x357*x6 + x358*x8 + x359*x8 + x36*L[185] + x381 + x385*x6 + x6*L[86] + x8*L[122] + L[37];
#pragma omp atomic
Ls[38] += x*x442 + x*x460 + x10*L[168] + x123*L[175] + x133*L[177] + x174*L[188] + x19*x229 + x19*x517 + x19*x528 + x19*L[94] + x21*x230 + x21*x518 + x21*L[135] + x227 + x23*L[186] + x239*x6 + x242*x8 + x32*x465 + x32*x529 + x32*x702 + x32*L[96] + x34*x466 + x34*x530 + x34*L[138] + x36*L[190] + x441 + x443*x6 + x444*x8 + x459 + x461*x6 + x515 + x526 + x6*L[87] + x8*L[123] + L[38];
#pragma omp atomic
Ls[39] += x*x380 + x*x472 + x*x477 + x*L[60] + x10*L[169] + x123*L[176] + x133*L[178] + x174*L[189] + x19*x391 + x19*x631 + x19*x634 + x19*L[95] + x21*x392 + x21*x632 + x21*L[136] + x23*L[187] + x32*x332 + x32*x626 + x32*x629 + x32*L[97] + x333*x34 + x34*x627 + x34*L[139] + x36*L[191] + x378 + x382*x6 + x384*x8 + x471 + x473*x6 + x474*x8 + x476 + x478*x6 + x6*L[88] + x8*L[124] + L[39];
#pragma omp atomic
Ls[40] += x*x485 + x*x506 + x10*L[170] + x123*L[177] + x133*L[179] + x174*L[190] + x19*x512 + x19*x529 + x19*x703 + x19*L[96] + x21*x513 + x21*x531 + x21*L[137] + x223 + x225*x32 + x226*x34 + x23*L[188] + x237*x6 + x241*x8 + x32*x521 + x32*x525 + x32*L[98] + x34*x522 + x34*L[140] + x36*L[192] + x484 + x486*x6 + x487*x8 + x504 + x508*x6 + x519 + x523 + x6*L[89] + x8*L[125] + L[40];
#pragma omp atomic
Ls[41] += x*x555 + x*x569 + x10*L[171] + x123*L[180] + x133*L[182] + x174*L[195] + x19*x221 + x19*x399 + x19*x413 + x19*L[99] + x21*x222 + x21*x400 + x21*L[141] + x219 + x23*L[193] + x235*x6 + x240*x8 + x32*x417 + x32*x571 + x32*x704 + x32*L[101] + x34*x418 + x34*x572 + x34*L[144] + x36*L[197] + x397 + x411 + x554 + x556*x6 + x557*x8 + x568 + x570*x6 + x6*L[90] + x8*L[126] + L[41];
#pragma omp atomic
Ls[42] += x*x577 + x*x581 + x10*L[172] + x123*L[181] + x133*L[183] + x174*L[196] + x19*x388 + x19*x541 + x19*x545 + x19*L[100] + x21*x389 + x21*x542 + x21*L[142] + x23*L[194] + x32*x415 + x32*x463 + x32*x636 + x32*L[102] + x34*x416 + x34*x464 + x34*L[145] + x36*L[198] + x386 + x393*x6 + x395*x8 + x539 + x543 + x576 + x578*x6 + x579*x8 + x580 + x582*x6 + x6*L[91] + x8*L[127] + L[42];
#pragma omp atomic
Ls[43] += x*x503 + x*x585 + x10*L[173] + x123*L[182] + x133*L[184] + x174*L[197] + x19*x417 + x19*x510 + x19*x638 + x19*L[101] + x21*x419 + x21*x511 + x21*L[143] + x23*L[195] + x32*x329 + x32*x534 + x32*x538 + x32*L[103] + x327 + x330*x34 + x334*x6 + x336*x8 + x34*x535 + x34*L[146] + x36*L[199] + x501 + x505*x6 + x507*x8 + x532 + x536 + x584 + x586*x6 + x6*L[92] + x8*L[128] + L[43];
#pragma omp atomic
Ls[44] += x*x592 + x*x608 + x10*L[174] + x123*L[183] + x133*L[185] + x174*L[198] + x19*x415 + x19*x611 + x19*x705 + x19*L[102] + x21*x418 + x21*x612 + x21*L[144] + x214 + x216*x32 + x217*x34 + x23*L[196] + x232*x6 + x238*x8 + x32*x404 + x32*x409 + x32*L[104] + x34*x405 + x34*L[147] + x36*L[200] + x402 + x407 + x591 + x593*x6 + x594*x8 + x6*x610 + x6*L[93] + x606 + x8*L[129] + L[44];
#pragma omp atomic
Ls[45] += x*x645 + x*x655 + x10*L[175] + x102 + x106*x19 + x108*x21 + x123*L[186] + x133*L[188] + x174*L[203] + x19*x264 + x19*x281 + x19*L[105] + x21*x265 + x21*L[148] + x229*x6 + x23*L[201] + x236*x8 + x262 + x279 + x289*x32 + x290*x34 + x32*x651 + x32*x706 + x32*L[107] + x34*x653 + x34*L[151] + x36*L[205] + x6*x646 + x6*x656 + x6*L[94] + x643 + x644 + x647*x8 + x8*L[130] + L[45];
#pragma omp atomic
Ls[46] += x*x660 + x*x663 + x10*L[176] + x123*L[187] + x133*L[189] + x174*L[204] + x19*x372 + x19*x432 + x19*x437 + x19*L[106] + x21*x374 + x21*x433 + x21*L[149] + x23*L[202] + x287*x32 + x288*x34 + x32*x564 + x32*x640 + x32*L[108] + x34*x566 + x34*L[152] + x36*L[206] + x368 + x391*x6 + x394*x8 + x430 + x435 + x6*x661 + x6*x664 + x6*L[95] + x658 + x659 + x662*x8 + x8*L[131] + L[46];
#pragma omp atomic
Ls[47] += x*x666 + x10*L[177] + x123*L[188] + x133*L[190] + x174*L[205] + x19*x289 + x19*x497 + x19*x551 + x19*L[107] + x21*x291 + x21*x499 + x21*L[150] + x23*L[203] + x282 + x284*x32 + x285*x34 + x32*x454 + x32*x548 + x32*L[109] + x34*x456 + x34*L[153] + x36*L[207] + x450 + x465*x6 + x467*x8 + x493 + x512*x6 + x514*x8 + x546 + x549 + x6*x667 + x6*L[96] + x665 + x8*L[132] + L[47];
#pragma omp atomic
Ls[48] += x*x605 + x*x669 + x10*L[178] + x123*L[189] + x133*L[191] + x174*L[206] + x19*x287 + x19*x601 + x19*x642 + x19*L[108] + x21*x290 + x21*x603 + x21*L[151] + x23*L[204] + x314 + x318*x32 + x32*x423 + x32*x428 + x32*L[110] + x320*x34 + x332*x6 + x335*x8 + x34*x424 + x34*L[154] + x36*L[208] + x421 + x426 + x589 + x6*x607 + x6*x670 + x6*L[97] + x609*x8 + x668 + x8*L[133] + L[48];
#pragma omp atomic
Ls[49] += x*x675 + x*x686 + x10*L[179] + x123*L[190] + x133*L[192] + x174*L[207] + x19*x284 + x19*x681 + x19*x707 + x19*L[109] + x21*x288 + x21*x683 + x21*L[152] + x225*x6 + x23*L[205] + x233*x8 + x268 + x270*x32 + x271*x34 + x274 + x276*x32 + x32*x94 + x32*L[111] + x34*x96 + x34*L[155] + x36*L[209] + x6*x676 + x6*x688 + x6*L[98] + x671 + x674 + x677*x8 + x8*L[134] + x90 + L[49];
#pragma omp atomic
Ls[50] += x*x708 + x*x709 + x10*L[180] + x105 + x109*x19 + x123*L[193] + x133*L[195] + x174*L[212] + x19*x49 + x19*x84 + x19*L[112] + x21*x50 + x21*x85 + x21*L[156] + x210*x32 + x211*x34 + x221*x6 + x23*L[210] + x230*x8 + x32*x712 + x32*x715 + x32*L[114] + x34*x713 + x34*L[159] + x36*L[214] + x47 + x6*x710 + x6*x716 + x6*L[99] + x711*x8 + x8*L[135] + x82 + x83*y + y*L[77] + L[50];
#pragma omp atomic
Ls[51] += x*x717 + x*x718 + x10*L[181] + x123*L[194] + x133*L[196] + x174*L[213] + x19*x346 + x19*x364 + x19*x375 + x19*L[113] + x208*x32 + x209*x34 + x21*x347 + x21*x365 + x21*L[157] + x23*L[211] + x32*x649 + x32*x654 + x32*L[115] + x34*x650 + x34*L[160] + x344 + x36*L[215] + x362 + x363*y + x371 + x388*x6 + x392*x8 + x6*x719 + x6*x721 + x6*L[100] + x720*x8 + x8*L[136] + y*L[78] + L[51];
#pragma omp atomic
Ls[52] += x*x722 + x10*L[182] + x123*L[195] + x133*L[197] + x174*L[214] + x19*x210 + x19*x491 + x19*x500 + x19*L[114] + x203 + x205*x32 + x206*x34 + x21*x212 + x21*x492 + x21*L[158] + x23*L[212] + x32*x560 + x32*x567 + x32*L[116] + x34*x561 + x34*L[161] + x36*L[216] + x489 + x490*y + x496 + x510*x6 + x513*x8 + x558 + x563 + x571*x6 + x573*x8 + x6*x723 + x6*L[101] + x8*L[137] + L[52];
#pragma omp atomic
Ls[53] += x*x724 + x10*L[183] + x123*L[196] + x133*L[198] + x174*L[215] + x19*x208 + x19*x597 + x19*x604 + x19*L[115] + x199 + x201*x32 + x202*x34 + x21*x211 + x21*x598 + x21*L[159] + x23*L[213] + x32*x448 + x32*x457 + x32*L[117] + x34*x449 + x34*L[162] + x36*L[217] + x446 + x453 + x463*x6 + x466*x8 + x595 + x596*y + x6*x611 + x6*x725 + x6*L[102] + x600 + x613*x8 + x8*L[138] + L[53];
#pragma omp atomic
Ls[54] += x*x673 + x*x726 + x10*L[184] + x123*L[197] + x133*L[199] + x174*L[216] + x19*x205 + x19*x679 + x19*x684 + x19*L[116] + x209*x21 + x21*x680 + x21*L[160] + x23*L[214] + x296 + x298*x32 + x299*x34 + x308 + x310*x32 + x311*x34 + x317 + x32*x321 + x32*L[118] + x329*x6 + x333*x8 + x34*L[163] + x36*L[218] + x6*x685 + x6*x727 + x6*L[103] + x672 + x678*y + x687*x8 + x8*L[139] + L[54];
#pragma omp atomic
Ls[55] += x*x728 + x*x731 + x10*L[185] + x123*L[198] + x133*L[200] + x174*L[217] + x19*x201 + x19*x734 + x19*x737 + x19*L[117] + x206*x21 + x21*x735 + x21*L[161] + x216*x6 + x226*x8 + x23*L[215] + x32*x56 + x32*x77 + x32*x97 + x32*L[119] + x34*x57 + x34*x78 + x34*L[164] + x36*L[219] + x54 + x6*x732 + x6*x739 + x6*L[104] + x729*y + x733*x8 + x75 + x8*L[140] + x93 + z*L[83] + L[55];
#pragma omp atomic
Ls[56] += x*x121 + x*x68 + x*x69 + x*L[84] + x119 + x19*x242 + x19*x740 + x19*L[123] + x21*L[171] + x241*x32 + x32*x741 + x32*L[125] + x34*L[174] + x6*x70 + x6*x71 + x6*L[120] + x66 + x67 + x8*L[165] + L[56];
#pragma omp atomic
Ls[57] += x*x118 + x*x304 + x*x326 + x*L[85] + x116 + x120*x6 + x19*x240 + x19*x691 + x19*L[126] + x21*L[175] + x303 + x305*x6 + x32*x336 + x32*x742 + x32*L[128] + x325 + x34*L[178] + x6*L[121] + x8*L[166] + L[57];
#pragma omp atomic
Ls[58] += x*x356 + x*x357 + x*x385 + x*L[86] + x19*x395 + x19*x743 + x19*L[127] + x21*L[176] + x238*x32 + x32*x693 + x32*L[129] + x34*L[179] + x354 + x355 + x358*x6 + x359*x6 + x383 + x6*L[122] + x8*L[167] + L[58];
#pragma omp atomic
Ls[59] += x*x443 + x*x461 + x19*x236 + x19*x616 + x19*L[130] + x21*L[180] + x234 + x242*x6 + x32*x467 + x32*x624 + x32*L[132] + x34*L[183] + x442 + x444*x6 + x460 + x6*L[123] + x614 + x622 + x8*L[168] + L[59];
#pragma omp atomic
Ls[60] += x*x382 + x*x473 + x*x478 + x*L[88] + x19*x394 + x19*x700 + x19*L[131] + x21*L[181] + x32*x335 + x32*x697 + x32*L[133] + x34*L[184] + x380 + x384*x6 + x472 + x474*x6 + x477 + x6*L[124] + x8*L[169] + L[60];
#pragma omp atomic
Ls[61] += x*x486 + x*x508 + x19*x514 + x19*x624 + x19*L[132] + x21*L[182] + x231 + x233*x32 + x241*x6 + x32*x619 + x32*L[134] + x34*L[185] + x485 + x487*x6 + x506 + x6*L[125] + x617 + x620 + x8*L[170] + L[61];
#pragma omp atomic
Ls[62] += x*x556 + x*x570 + x19*x230 + x19*x518 + x19*L[135] + x21*L[186] + x228 + x240*x6 + x32*x531 + x32*x573 + x32*L[137] + x34*L[189] + x516 + x527 + x555 + x557*x6 + x569 + x6*L[126] + x8*L[171] + L[62];
#pragma omp atomic
Ls[63] += x*x578 + x*x582 + x19*x392 + x19*x632 + x19*L[136] + x21*L[187] + x32*x466 + x32*x530 + x32*L[138] + x34*L[190] + x390 + x395*x6 + x577 + x579*x6 + x581 + x6*L[127] + x630 + x633 + x8*L[172] + L[63];
#pragma omp atomic
Ls[64] += x*x505 + x*x586 + x19*x513 + x19*x531 + x19*L[137] + x21*L[188] + x32*x333 + x32*x627 + x32*L[139] + x331 + x336*x6 + x34*L[191] + x503 + x507*x6 + x585 + x6*L[128] + x625 + x628 + x8*L[173] + L[64];
#pragma omp atomic
Ls[65] += x*x593 + x*x610 + x19*x530 + x19*x613 + x19*L[138] + x21*L[189] + x224 + x226*x32 + x238*x6 + x32*x522 + x32*L[140] + x34*L[192] + x520 + x524 + x592 + x594*x6 + x6*L[129] + x608 + x8*L[174] + L[65];
#pragma omp atomic
Ls[66] += x*x646 + x*x656 + x19*x222 + x19*x400 + x19*L[141] + x21*L[193] + x220 + x236*x6 + x32*x419 + x32*x657 + x32*L[143] + x34*L[196] + x398 + x412 + x6*x647 + x6*L[130] + x645 + x655 + x8*L[175] + L[66];
#pragma omp atomic
Ls[67] += x*x661 + x*x664 + x19*x389 + x19*x542 + x19*L[142] + x21*L[194] + x32*x418 + x32*x572 + x32*L[144] + x34*L[197] + x387 + x394*x6 + x540 + x544 + x6*x662 + x6*L[131] + x660 + x663 + x8*L[176] + L[67];
#pragma omp atomic
Ls[68] += x*x667 + x19*x419 + x19*x511 + x19*L[143] + x21*L[195] + x32*x416 + x32*x464 + x32*L[145] + x34*L[198] + x414 + x462 + x467*x6 + x509 + x514*x6 + x6*L[132] + x635 + x637 + x666 + x8*L[177] + L[68];
#pragma omp atomic
Ls[69] += x*x607 + x*x670 + x19*x418 + x19*x612 + x19*L[144] + x21*L[196] + x32*x330 + x32*x535 + x32*L[146] + x328 + x335*x6 + x34*L[199] + x533 + x537 + x6*x609 + x6*L[133] + x605 + x669 + x8*L[178] + L[69];
#pragma omp atomic
Ls[70] += x*x676 + x*x688 + x19*x416 + x19*x689 + x19*L[145] + x21*L[197] + x215 + x217*x32 + x233*x6 + x32*x405 + x32*L[147] + x34*L[200] + x403 + x408 + x6*x677 + x6*L[134] + x675 + x686 + x8*L[179] + L[70];
#pragma omp atomic
Ls[71] += x*x710 + x*x716 + x104 + x108*x19 + x19*x265 + x19*L[148] + x21*L[201] + x230*x6 + x263 + x280 + x291*x32 + x32*x714 + x32*L[150] + x34*L[204] + x6*x711 + x6*L[135] + x708 + x709 + x8*L[180] + L[71];
#pragma omp atomic
Ls[72] += x*x719 + x*x721 + x19*x374 + x19*x433 + x19*L[149] + x21*L[202] + x290*x32 + x32*x653 + x32*L[151] + x34*L[205] + x370 + x392*x6 + x431 + x436 + x6*x720 + x6*L[136] + x717 + x718 + x8*L[181] + L[72];
#pragma omp atomic
Ls[73] += x*x723 + x19*x291 + x19*x499 + x19*L[150] + x21*L[203] + x286 + x288*x32 + x32*x566 + x32*L[152] + x34*L[206] + x495 + x513*x6 + x550 + x562 + x573*x6 + x6*L[137] + x639 + x722 + x8*L[182] + L[73];
#pragma omp atomic
Ls[74] += x*x725 + x19*x290 + x19*x603 + x19*L[151] + x21*L[204] + x283 + x285*x32 + x32*x456 + x32*L[153] + x34*L[207] + x452 + x466*x6 + x547 + x599 + x6*x613 + x6*L[138] + x641 + x724 + x8*L[183] + L[74];
#pragma omp atomic
Ls[75] += x*x685 + x*x727 + x19*x288 + x19*x683 + x19*L[152] + x21*L[205] + x316 + x32*x320 + x32*x424 + x32*L[154] + x333*x6 + x34*L[208] + x422 + x427 + x6*x687 + x6*L[139] + x673 + x726 + x8*L[184] + L[75];
#pragma omp atomic
Ls[76] += x*x732 + x*x739 + x19*x285 + x19*x736 + x19*L[153] + x21*L[206] + x226*x6 + x269 + x271*x32 + x275 + x32*x96 + x32*L[155] + x34*L[209] + x6*x733 + x6*L[140] + x728 + x731 + x8*L[185] + x92 + L[76];
#pragma omp atomic
Ls[77] += x*x744 + x*x745 + x107 + x19*x50 + x19*x85 + x19*L[156] + x21*L[210] + x212*x32 + x222*x6 + x32*x747 + x32*L[158] + x34*L[213] + x48 + x6*x746 + x6*L[141] + x8*L[186] + x83 + x84*y + y*L[112] + L[77];
#pragma omp atomic
Ls[78] += x*x748 + x*x749 + x19*x347 + x19*x365 + x19*L[157] + x21*L[211] + x211*x32 + x32*x713 + x32*L[159] + x34*L[214] + x345 + x363 + x364*y + x373 + x389*x6 + x6*x750 + x6*L[142] + x8*L[187] + y*L[113] + L[78];
#pragma omp atomic
Ls[79] += x*x751 + x19*x212 + x19*x492 + x19*L[158] + x207 + x209*x32 + x21*L[212] + x32*x650 + x32*L[160] + x34*L[215] + x490 + x491*y + x498 + x511*x6 + x6*x657 + x6*L[143] + x648 + x652 + x8*L[188] + L[79];
#pragma omp atomic
Ls[80] += x*x752 + x19*x211 + x19*x598 + x19*L[159] + x204 + x206*x32 + x21*L[213] + x32*x561 + x32*L[161] + x34*L[216] + x559 + x565 + x572*x6 + x596 + x597*y + x6*x612 + x6*L[144] + x602 + x8*L[189] + L[80];
#pragma omp atomic
Ls[81] += x*x753 + x19*x209 + x19*x680 + x19*L[160] + x200 + x202*x32 + x21*L[214] + x32*x449 + x32*L[162] + x34*L[217] + x447 + x455 + x464*x6 + x6*x689 + x6*L[145] + x678 + x679*y + x682 + x8*L[190] + L[81];
#pragma omp atomic
Ls[82] += x*x730 + x*x754 + x19*x206 + x19*x735 + x19*L[161] + x21*L[215] + x297 + x299*x32 + x309 + x311*x32 + x319 + x32*L[163] + x330*x6 + x34*L[218] + x6*x738 + x6*L[146] + x729 + x734*y + x8*L[191] + L[82];
#pragma omp atomic
Ls[83] += x*x755 + x*x758 + x19*x202 + x19*x760 + x19*L[162] + x21*L[216] + x217*x6 + x32*x57 + x32*x78 + x32*L[164] + x34*L[219] + x55 + x6*x759 + x6*L[147] + x756*y + x76 + x8*L[192] + x95 + z*L[119] + L[83];
#pragma omp atomic
Ls[84] += x*x70 + x*x71 + x*L[120] + x121 + x19*L[168] + x32*L[170] + x6*L[165] + x68 + x69 + L[84];
#pragma omp atomic
Ls[85] += x*x120 + x*x305 + x*L[121] + x118 + x19*L[171] + x304 + x32*L[173] + x326 + x6*L[166] + L[85];
#pragma omp atomic
Ls[86] += x*x358 + x*x359 + x*L[122] + x19*L[172] + x32*L[174] + x356 + x357 + x385 + x6*L[167] + L[86];
#pragma omp atomic
Ls[87] += x*x444 + x19*L[175] + x239 + x32*L[177] + x443 + x461 + x6*L[168] + x690 + x695 + L[87];
#pragma omp atomic
Ls[88] += x*x384 + x*x474 + x*L[124] + x19*L[176] + x32*L[178] + x382 + x473 + x478 + x6*L[169] + L[88];
#pragma omp atomic
Ls[89] += x*x487 + x19*L[177] + x237 + x32*L[179] + x486 + x508 + x6*L[170] + x692 + x694 + L[89];
#pragma omp atomic
Ls[90] += x*x557 + x19*L[180] + x235 + x32*L[182] + x556 + x570 + x6*L[171] + x615 + x623 + L[90];
#pragma omp atomic
Ls[91] += x*x579 + x19*L[181] + x32*L[183] + x393 + x578 + x582 + x6*L[172] + x699 + x701 + L[91];
#pragma omp atomic
Ls[92] += x*x507 + x19*L[182] + x32*L[184] + x334 + x505 + x586 + x6*L[173] + x696 + x698 + L[92];
#pragma omp atomic
Ls[93] += x*x594 + x19*L[183] + x232 + x32*L[185] + x593 + x6*L[174] + x610 + x618 + x621 + L[93];
#pragma omp atomic
Ls[94] += x*x647 + x19*L[186] + x229 + x32*L[188] + x517 + x528 + x6*L[175] + x646 + x656 + L[94];
#pragma omp atomic
Ls[95] += x*x662 + x19*L[187] + x32*L[189] + x391 + x6*L[176] + x631 + x634 + x661 + x664 + L[95];
#pragma omp atomic
Ls[96] += x19*L[188] + x32*L[190] + x465 + x512 + x529 + x6*L[177] + x667 + x702 + x703 + L[96];
#pragma omp atomic
Ls[97] += x*x609 + x19*L[189] + x32*L[191] + x332 + x6*L[178] + x607 + x626 + x629 + x670 + L[97];
#pragma omp atomic
Ls[98] += x*x677 + x19*L[190] + x225 + x32*L[192] + x521 + x525 + x6*L[179] + x676 + x688 + L[98];
#pragma omp atomic
Ls[99] += x*x711 + x19*L[193] + x221 + x32*L[195] + x399 + x413 + x6*L[180] + x710 + x716 + L[99];
#pragma omp atomic
Ls[100] += x*x720 + x19*L[194] + x32*L[196] + x388 + x541 + x545 + x6*L[181] + x719 + x721 + L[100];
#pragma omp atomic
Ls[101] += x19*L[195] + x32*L[197] + x417 + x510 + x571 + x6*L[182] + x638 + x704 + x723 + L[101];
#pragma omp atomic
Ls[102] += x19*L[196] + x32*L[198] + x415 + x463 + x6*L[183] + x611 + x636 + x705 + x725 + L[102];
#pragma omp atomic
Ls[103] += x*x687 + x19*L[197] + x32*L[199] + x329 + x534 + x538 + x6*L[184] + x685 + x727 + L[103];
#pragma omp atomic
Ls[104] += x*x733 + x19*L[198] + x216 + x32*L[200] + x404 + x409 + x6*L[185] + x732 + x739 + L[104];
#pragma omp atomic
Ls[105] += x*x746 + x106 + x19*L[201] + x264 + x281 + x32*L[203] + x6*L[186] + x744 + x745 + L[105];
#pragma omp atomic
Ls[106] += x*x750 + x19*L[202] + x32*L[204] + x372 + x432 + x437 + x6*L[187] + x748 + x749 + L[106];
#pragma omp atomic
Ls[107] += x19*L[203] + x289 + x32*L[205] + x497 + x551 + x6*L[188] + x651 + x706 + x751 + L[107];
#pragma omp atomic
Ls[108] += x19*L[204] + x287 + x32*L[206] + x564 + x6*L[189] + x601 + x640 + x642 + x752 + L[108];
#pragma omp atomic
Ls[109] += x19*L[205] + x284 + x32*L[207] + x454 + x548 + x6*L[190] + x681 + x707 + x753 + L[109];
#pragma omp atomic
Ls[110] += x*x738 + x19*L[206] + x318 + x32*L[208] + x423 + x428 + x6*L[191] + x730 + x754 + L[110];
#pragma omp atomic
Ls[111] += x*x759 + x19*L[207] + x270 + x276 + x32*L[209] + x6*L[192] + x755 + x758 + x94 + L[111];
#pragma omp atomic
Ls[112] += x*x761 + x109 + x19*L[210] + x32*L[212] + x49 + x6*L[193] + x84 + x85*y + y*L[156] + L[112];
#pragma omp atomic
Ls[113] += x*x762 + x19*L[211] + x32*L[213] + x346 + x364 + x365*y + x375 + x6*L[194] + y*L[157] + L[113];
#pragma omp atomic
Ls[114] += x19*L[212] + x210 + x32*L[214] + x491 + x492*y + x500 + x6*L[195] + x712 + x715 + L[114];
#pragma omp atomic
Ls[115] += x19*L[213] + x208 + x32*L[215] + x597 + x598*y + x6*L[196] + x604 + x649 + x654 + L[115];
#pragma omp atomic
Ls[116] += x19*L[214] + x205 + x32*L[216] + x560 + x567 + x6*L[197] + x679 + x680*y + x684 + L[116];
#pragma omp atomic
Ls[117] += x19*L[215] + x201 + x32*L[217] + x448 + x457 + x6*L[198] + x734 + x735*y + x737 + L[117];
#pragma omp atomic
Ls[118] += x*x757 + x19*L[216] + x298 + x310 + x32*L[218] + x321 + x6*L[199] + x756 + x760*y + L[118];
#pragma omp atomic
Ls[119] += x*x763 + x19*L[217] + x32*L[219] + x56 + x6*L[200] + x764*y + x77 + x97 + z*L[164] + L[119];
#pragma omp atomic
Ls[120] += x*L[165] + x70 + x71 + L[120];
#pragma omp atomic
Ls[121] += x*L[166] + x120 + x305 + L[121];
#pragma omp atomic
Ls[122] += x*L[167] + x358 + x359 + L[122];
#pragma omp atomic
Ls[123] += x242 + x444 + x740 + L[123];
#pragma omp atomic
Ls[124] += x*L[169] + x384 + x474 + L[124];
#pragma omp atomic
Ls[125] += x241 + x487 + x741 + L[125];
#pragma omp atomic
Ls[126] += x240 + x557 + x691 + L[126];
#pragma omp atomic
Ls[127] += x395 + x579 + x743 + L[127];
#pragma omp atomic
Ls[128] += x336 + x507 + x742 + L[128];
#pragma omp atomic
Ls[129] += x238 + x594 + x693 + L[129];
#pragma omp atomic
Ls[130] += x236 + x616 + x647 + L[130];
#pragma omp atomic
Ls[131] += x394 + x662 + x700 + L[131];
#pragma omp atomic
Ls[132] += x467 + x514 + x624 + L[132];
#pragma omp atomic
Ls[133] += x335 + x609 + x697 + L[133];
#pragma omp atomic
Ls[134] += x233 + x619 + x677 + L[134];
#pragma omp atomic
Ls[135] += x230 + x518 + x711 + L[135];
#pragma omp atomic
Ls[136] += x392 + x632 + x720 + L[136];
#pragma omp atomic
Ls[137] += x513 + x531 + x573 + L[137];
#pragma omp atomic
Ls[138] += x466 + x530 + x613 + L[138];
#pragma omp atomic
Ls[139] += x333 + x627 + x687 + L[139];
#pragma omp atomic
Ls[140] += x226 + x522 + x733 + L[140];
#pragma omp atomic
Ls[141] += x222 + x400 + x746 + L[141];
#pragma omp atomic
Ls[142] += x389 + x542 + x750 + L[142];
#pragma omp atomic
Ls[143] += x419 + x511 + x657 + L[143];
#pragma omp atomic
Ls[144] += x418 + x572 + x612 + L[144];
#pragma omp atomic
Ls[145] += x416 + x464 + x689 + L[145];
#pragma omp atomic
Ls[146] += x330 + x535 + x738 + L[146];
#pragma omp atomic
Ls[147] += x217 + x405 + x759 + L[147];
#pragma omp atomic
Ls[148] += x108 + x265 + x761 + L[148];
#pragma omp atomic
Ls[149] += x374 + x433 + x762 + L[149];
#pragma omp atomic
Ls[150] += x291 + x499 + x714 + L[150];
#pragma omp atomic
Ls[151] += x290 + x603 + x653 + L[151];
#pragma omp atomic
Ls[152] += x288 + x566 + x683 + L[152];
#pragma omp atomic
Ls[153] += x285 + x456 + x736 + L[153];
#pragma omp atomic
Ls[154] += x320 + x424 + x757 + L[154];
#pragma omp atomic
Ls[155] += x271 + x763 + x96 + L[155];
#pragma omp atomic
Ls[156] += x50 + x85 + y*L[210] + L[156];
#pragma omp atomic
Ls[157] += x347 + x365 + y*L[211] + L[157];
#pragma omp atomic
Ls[158] += x212 + x492 + x747 + L[158];
#pragma omp atomic
Ls[159] += x211 + x598 + x713 + L[159];
#pragma omp atomic
Ls[160] += x209 + x650 + x680 + L[160];
#pragma omp atomic
Ls[161] += x206 + x561 + x735 + L[161];
#pragma omp atomic
Ls[162] += x202 + x449 + x760 + L[162];
#pragma omp atomic
Ls[163] += x299 + x311 + x764 + L[163];
#pragma omp atomic
Ls[164] += x57 + x78 + z*L[219] + L[164];
#pragma omp atomic
Ls[165] += L[165];
#pragma omp atomic
Ls[166] += L[166];
#pragma omp atomic
Ls[167] += L[167];
#pragma omp atomic
Ls[168] += L[168];
#pragma omp atomic
Ls[169] += L[169];
#pragma omp atomic
Ls[170] += L[170];
#pragma omp atomic
Ls[171] += L[171];
#pragma omp atomic
Ls[172] += L[172];
#pragma omp atomic
Ls[173] += L[173];
#pragma omp atomic
Ls[174] += L[174];
#pragma omp atomic
Ls[175] += L[175];
#pragma omp atomic
Ls[176] += L[176];
#pragma omp atomic
Ls[177] += L[177];
#pragma omp atomic
Ls[178] += L[178];
#pragma omp atomic
Ls[179] += L[179];
#pragma omp atomic
Ls[180] += L[180];
#pragma omp atomic
Ls[181] += L[181];
#pragma omp atomic
Ls[182] += L[182];
#pragma omp atomic
Ls[183] += L[183];
#pragma omp atomic
Ls[184] += L[184];
#pragma omp atomic
Ls[185] += L[185];
#pragma omp atomic
Ls[186] += L[186];
#pragma omp atomic
Ls[187] += L[187];
#pragma omp atomic
Ls[188] += L[188];
#pragma omp atomic
Ls[189] += L[189];
#pragma omp atomic
Ls[190] += L[190];
#pragma omp atomic
Ls[191] += L[191];
#pragma omp atomic
Ls[192] += L[192];
#pragma omp atomic
Ls[193] += L[193];
#pragma omp atomic
Ls[194] += L[194];
#pragma omp atomic
Ls[195] += L[195];
#pragma omp atomic
Ls[196] += L[196];
#pragma omp atomic
Ls[197] += L[197];
#pragma omp atomic
Ls[198] += L[198];
#pragma omp atomic
Ls[199] += L[199];
#pragma omp atomic
Ls[200] += L[200];
#pragma omp atomic
Ls[201] += L[201];
#pragma omp atomic
Ls[202] += L[202];
#pragma omp atomic
Ls[203] += L[203];
#pragma omp atomic
Ls[204] += L[204];
#pragma omp atomic
Ls[205] += L[205];
#pragma omp atomic
Ls[206] += L[206];
#pragma omp atomic
Ls[207] += L[207];
#pragma omp atomic
Ls[208] += L[208];
#pragma omp atomic
Ls[209] += L[209];
#pragma omp atomic
Ls[210] += L[210];
#pragma omp atomic
Ls[211] += L[211];
#pragma omp atomic
Ls[212] += L[212];
#pragma omp atomic
Ls[213] += L[213];
#pragma omp atomic
Ls[214] += L[214];
#pragma omp atomic
Ls[215] += L[215];
#pragma omp atomic
Ls[216] += L[216];
#pragma omp atomic
Ls[217] += L[217];
#pragma omp atomic
Ls[218] += L[218];
#pragma omp atomic
Ls[219] += L[219];

}

void L2P_10(double x, double y, double z, double * L, double * F) {
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
x0 = x*y;
x1 = x*z;
x2 = y*z;
x3 = (x*x);
x4 = (1.0/2.0)*x3;
x5 = (x*x*x);
x6 = (1.0/6.0)*x5;
x7 = (x*x*x*x);
x8 = (1.0/24.0)*x7;
x9 = pow(x, 5);
x10 = (1.0/120.0)*x9;
x11 = pow(x, 6);
x12 = (1.0/720.0)*x11;
x13 = pow(x, 7);
x14 = (1.0/5040.0)*x13;
x15 = (1.0/40320.0)*pow(x, 8);
x16 = (y*y);
x17 = (1.0/2.0)*x16;
x18 = (y*y*y);
x19 = (1.0/6.0)*x18;
x20 = (y*y*y*y);
x21 = (1.0/24.0)*x20;
x22 = pow(y, 5);
x23 = (1.0/120.0)*x22;
x24 = pow(y, 6);
x25 = (1.0/720.0)*x24;
x26 = pow(y, 7);
x27 = (1.0/5040.0)*x26;
x28 = (1.0/40320.0)*pow(y, 8);
x29 = (z*z);
x30 = (1.0/2.0)*x29;
x31 = (z*z*z);
x32 = (1.0/6.0)*x31;
x33 = (z*z*z*z);
x34 = (1.0/24.0)*x33;
x35 = pow(z, 5);
x36 = (1.0/120.0)*x35;
x37 = pow(z, 6);
x38 = (1.0/720.0)*x37;
x39 = pow(z, 7);
x40 = (1.0/5040.0)*x39;
x41 = (1.0/40320.0)*pow(z, 8);
x42 = (1.0/4.0)*x3;
x43 = x16*x42;
x44 = (1.0/12.0)*x3;
x45 = x18*x44;
x46 = (1.0/48.0)*x3;
x47 = x20*x46;
x48 = (1.0/240.0)*x3;
x49 = x22*x48;
x50 = (1.0/1440.0)*x3;
x51 = x24*x50;
x52 = (1.0/10080.0)*x3;
x53 = x29*x42;
x54 = x31*x44;
x55 = x33*x46;
x56 = x35*x48;
x57 = x37*x50;
x58 = (1.0/12.0)*x5;
x59 = x16*x58;
x60 = (1.0/36.0)*x5;
x61 = x18*x60;
x62 = (1.0/144.0)*x5;
x63 = x20*x62;
x64 = (1.0/720.0)*x5;
x65 = x22*x64;
x66 = (1.0/4320.0)*x5;
x67 = x29*x58;
x68 = x31*x60;
x69 = x33*x62;
x70 = x35*x64;
x71 = (1.0/48.0)*x7;
x72 = x16*x71;
x73 = (1.0/144.0)*x7;
x74 = x18*x73;
x75 = (1.0/576.0)*x7;
x76 = x20*x75;
x77 = (1.0/2880.0)*x7;
x78 = x29*x71;
x79 = x31*x73;
x80 = x33*x75;
x81 = (1.0/240.0)*x9;
x82 = x16*x81;
x83 = (1.0/720.0)*x9;
x84 = x18*x83;
x85 = (1.0/2880.0)*x9;
x86 = x29*x81;
x87 = x31*x83;
x88 = (1.0/1440.0)*x11;
x89 = x16*x88;
x90 = (1.0/4320.0)*x11;
x91 = x29*x88;
x92 = (1.0/10080.0)*x13;
x93 = x16*x29;
x94 = (1.0/4.0)*x93;
x95 = x16*x31;
x96 = (1.0/12.0)*x95;
x97 = x16*x33;
x98 = (1.0/48.0)*x97;
x99 = x16*x35;
x100 = (1.0/240.0)*x99;
x101 = (1.0/1440.0)*x16*x37;
x102 = x18*x29;
x103 = (1.0/12.0)*x102;
x104 = x18*x31;
x105 = (1.0/36.0)*x104;
x106 = x18*x33;
x107 = (1.0/144.0)*x106;
x108 = (1.0/720.0)*x18*x35;
x109 = x20*x29;
x110 = (1.0/48.0)*x109;
x111 = x20*x31;
x112 = (1.0/144.0)*x111;
x113 = (1.0/576.0)*x20*x33;
x114 = x22*x29;
x115 = (1.0/240.0)*x114;
x116 = (1.0/720.0)*x22*x31;
x117 = (1.0/1440.0)*x24*x29;
x118 = (1.0/24.0)*x3;
x119 = (1.0/96.0)*x3;
x120 = (1.0/480.0)*x3;
x121 = (1.0/288.0)*x3;
x122 = (1.0/72.0)*x5;
x123 = (1.0/288.0)*x5;
x124 = (1.0/288.0)*x7;
#pragma omp atomic
F[0] += (1.0/362880.0)*pow(x, 9)*L[165] + x*x100*L[153] + x*x101*L[207] + x*x103*L[73] + x*x105*L[108] + x*x107*L[152] + x*x108*L[206] + x*x110*L[107] + x*x112*L[151] + x*x113*L[205] + x*x115*L[150] + x*x116*L[204] + x*x117*L[203] + x*x17*L[13] + x*x19*L[26] + x*x21*L[45] + x*x23*L[71] + x*x25*L[105] + x*x27*L[148] + x*x28*L[201] + x*x30*L[15] + x*x32*L[29] + x*x34*L[49] + x*x36*L[76] + x*x38*L[111] + x*x40*L[155] + x*x41*L[209] + x*x94*L[47] + x*x96*L[74] + x*x98*L[109] + x*L[1] + x0*x30*L[28] + x0*x32*L[48] + x0*x34*L[75] + x0*x36*L[110] + x0*x38*L[154] + x0*x40*L[208] + x0*z*L[14] + x0*L[5] + x1*x17*L[27] + x1*x19*L[46] + x1*x21*L[72] + x1*x23*L[106] + x1*x25*L[149] + x1*x27*L[202] + x1*L[6] + x10*x2*L[88] + x10*y*L[57] + x10*z*L[58] + x10*L[35] + x100*L[117] + x101*L[162] + x102*x118*L[101] + x102*x122*L[137] + x102*x124*L[182] + x103*L[52] + (1.0/72.0)*x104*x3*L[144] + (1.0/216.0)*x104*x5*L[189] + x105*L[80] + x106*x121*L[197] + x107*L[116] + x108*L[161] + x109*x119*L[143] + x109*x123*L[188] + x110*L[79] + x111*x121*L[196] + x112*L[115] + x113*L[160] + x114*x120*L[195] + x115*L[114] + x116*L[159] + x117*L[158] + x118*x95*L[102] + x119*x97*L[145] + x12*x2*L[124] + x12*y*L[85] + x12*z*L[86] + x12*L[56] + x120*x99*L[198] + x122*x95*L[138] + x123*x97*L[190] + x124*x95*L[183] + x14*x2*L[169] + x14*y*L[121] + x14*z*L[122] + x14*L[84] + x15*y*L[166] + x15*z*L[167] + x15*L[120] + (1.0/10080.0)*x16*x39*L[217] + x16*x92*L[168] + x17*z*L[17] + x17*L[7] + (1.0/4320.0)*x18*x37*L[216] + x18*x90*L[171] + x19*z*L[31] + x19*L[16] + x2*x4*L[24] + x2*x6*L[39] + x2*x8*L[60] + x2*L[8] + (1.0/2880.0)*x20*x35*L[215] + x20*x85*L[175] + x21*z*L[51] + x21*L[30] + (1.0/2880.0)*x22*x33*L[214] + x22*x77*L[180] + x23*z*L[78] + x23*L[50] + (1.0/4320.0)*x24*x31*L[213] + x24*x66*L[186] + x25*z*L[113] + x25*L[77] + (1.0/10080.0)*x26*x29*L[212] + x26*x52*L[193] + x27*z*L[157] + x27*L[112] + x28*z*L[211] + x28*L[156] + x29*x92*L[170] + (1.0/8.0)*x3*x93*L[68] + x30*y*L[18] + x30*L[9] + x31*x90*L[174] + x32*y*L[33] + x32*L[19] + x33*x85*L[179] + x34*y*L[54] + x34*L[34] + x35*x77*L[185] + x36*y*L[82] + x36*L[55] + x37*x66*L[192] + x38*y*L[118] + x38*L[83] + x39*x52*L[200] + x4*y*L[11] + x4*z*L[12] + x4*L[4] + x40*y*L[163] + x40*L[119] + x41*y*L[218] + x41*L[164] + x43*z*L[42] + x43*L[23] + x45*z*L[67] + x45*L[41] + x47*z*L[100] + x47*L[66] + x49*z*L[142] + x49*L[99] + (1.0/24.0)*x5*x93*L[96] + x51*z*L[194] + x51*L[141] + x53*y*L[43] + x53*L[25] + x54*y*L[69] + x54*L[44] + x55*y*L[103] + x55*L[70] + x56*y*L[146] + x56*L[104] + x57*y*L[199] + x57*L[147] + x59*z*L[63] + x59*L[38] + x6*y*L[21] + x6*z*L[22] + x6*L[10] + x61*z*L[95] + x61*L[62] + x63*z*L[136] + x63*L[94] + x65*z*L[187] + x65*L[135] + x67*y*L[64] + x67*L[40] + x68*y*L[97] + x68*L[65] + x69*y*L[139] + x69*L[98] + (1.0/96.0)*x7*x93*L[132] + x70*y*L[191] + x70*L[140] + x72*z*L[91] + x72*L[59] + x74*z*L[131] + x74*L[90] + x76*z*L[181] + x76*L[130] + x78*y*L[92] + x78*L[61] + x79*y*L[133] + x79*L[93] + x8*y*L[36] + x8*z*L[37] + x8*L[20] + x80*y*L[184] + x80*L[134] + x82*z*L[127] + x82*L[87] + x84*z*L[176] + x84*L[126] + x86*y*L[128] + x86*L[89] + x87*y*L[178] + x87*L[129] + x89*z*L[172] + x89*L[123] + (1.0/480.0)*x9*x93*L[177] + x91*y*L[173] + x91*L[125] + x94*L[32] + x96*L[53] + x98*L[81] + (1.0/362880.0)*pow(y, 9)*L[210] + y*L[2] + (1.0/362880.0)*pow(z, 9)*L[219] + z*L[3] + L[0];

}

void M2P_10(double x, double y, double z, double * M, double * F) {
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
double x164;
double x165;
double x166;
double x167;
double x168;
double x169;
double x170;
double x171;
double x172;
double x173;
double x174;
double x175;
double x176;
double x177;
double x178;
double x179;
double x180;
double x181;
double x182;
double x183;
double x184;
double x185;
double x186;
double x187;
double x188;
double x189;
double x190;
double x191;
double x192;
double x193;
double x194;
double x195;
double x196;
double x197;
double x198;
double x199;
double x200;
double x201;
double x202;
double x203;
double x204;
double x205;
double x206;
double x207;
double x208;
double x209;
double x210;
double x211;
double x212;
double x213;
double x214;
double x215;
double x216;
double x217;
double x218;
double x219;
double x220;
double x221;
double x222;
double x223;
double x224;
double x225;
double x226;
double x227;
double x228;
double x229;
double x230;
double x231;
double x232;
double x233;
double x234;
double x235;
double x236;
double x237;
double x238;
double x239;
double x240;
double x241;
double x242;
double x243;
double x244;
double x245;
double x246;
double x247;
double x248;
double x249;
double x250;
double x251;
double x252;
double x253;
double x254;
double x255;
double x256;
double x257;
double x258;
double x259;
double x260;
double x261;
double x262;
double x263;
double x264;
double x265;
double x266;
double x267;
double x268;
double x269;
double x270;
double x271;
double x272;
double x273;
double x274;
double x275;
double x276;
double x277;
double x278;
double x279;
double x280;
double x281;
double x282;
double x283;
double x284;
double x285;
double x286;
double x287;
double x288;
double x289;
double x290;
double x291;
double x292;
double x293;
double x294;
double x295;
double x296;
double x297;
double x298;
double x299;
double x300;
double x301;
double x302;
double x303;
double x304;
double x305;
double x306;
double x307;
double x308;
double x309;
double x310;
double x311;
double x312;
double x313;
double x314;
double x315;
double x316;
double x317;
double x318;
double x319;
double x320;
double x321;
double x322;
double x323;
double x324;
double x325;
double x326;
double x327;
double x328;
double x329;
double x330;
double x331;
double x332;
double x333;
double x334;
double x335;
double x336;
double x337;
double x338;
double x339;
double x340;
double x341;
double x342;
double x343;
double x344;
double x345;
double x346;
double x347;
double x348;
double x349;
double x350;
double x351;
double x352;
double x353;
double x354;
double x355;
double x356;
double x357;
double x358;
double x359;
double x360;
double x361;
double x362;
double x363;
double x364;
double x365;
double x366;
double x367;
double x368;
double x369;
double x370;
double x371;
double x372;
double x373;
double x374;
double x375;
double x376;
double x377;
double x378;
double x379;
double x380;
double x381;
double x382;
double x383;
double x384;
double x385;
double x386;
double x387;
double x388;
double x389;
double x390;
double x391;
double x392;
double x393;
double x394;
double x395;
double x396;
double x397;
double x398;
double x399;
double x400;
double x401;
double x402;
double x403;
double x404;
double x405;
double x406;
double x407;
double x408;
double x409;
double x410;
double x411;
double x412;
double x413;
double x414;
double x415;
double x416;
double x417;
double x418;
double x419;
double x420;
double x421;
double x422;
double x423;
double x424;
double x425;
double x426;
double x427;
double x428;
double x429;
double x430;
double x431;
double x432;
double x433;
double x434;
double x435;
double x436;
double x437;
double x438;
double x439;
double x440;
double x441;
double x442;
double x443;
double x444;
double x445;
double x446;
double x447;
double x448;
double x449;
double x450;
double x451;
double x452;
double x453;
double x454;
double x455;
double x456;
double x457;
double x458;
double x459;
double x460;
double x461;
double x462;
double x463;
double x464;
double x465;
double x466;
double x467;
double x468;
double x469;
double x470;
double x471;
double x472;
double x473;
double x474;
double x475;
double x476;
double x477;
double x478;
double x479;
double x480;
double x481;
double x482;
double x483;
double x484;
double x485;
double x486;
double x487;
double x488;
double x489;
double x490;
double x491;
double x492;
double x493;
double x494;
double x495;
double x496;
double x497;
double x498;
double x499;
double x500;
double x501;
double x502;
double x503;
double x504;
double x505;
double x506;
double x507;
double x508;
double x509;
double x510;
double x511;
double x512;
double x513;
double x514;
double x515;
double x516;
double x517;
double x518;
double x519;
double x520;
double x521;
double x522;
double x523;
double x524;
double x525;
double x526;
double x527;
double x528;
double x529;
double x530;
double x531;
double x532;
double x533;
double x534;
double x535;
double x536;
double x537;
double x538;
double x539;
double x540;
double x541;
double x542;
double x543;
x0 = (x*x);
x1 = (y*y);
x2 = x0 + x1 + (z*z);
x3 = pow(x2, -1.5);
x4 = 1.0*x3;
x5 = pow(x2, -2.5);
x6 = 3.0*x5;
x7 = x*y;
x8 = x*z;
x9 = y*z;
x10 = pow(x2, -3.5);
x11 = 15.0*x10;
x12 = x7*z;
x13 = -x4;
x14 = x0*x6;
x15 = x1*x6;
x16 = 9.0*x5;
x17 = x0*x11;
x18 = -x17;
x19 = x*(x16 + x18);
x20 = x18 + x6;
x21 = x20*y;
x22 = x1*x11;
x23 = -x22;
x24 = y*(x16 + x23);
x25 = x20*z;
x26 = z*(x23 + x6);
x27 = 1.0*x;
x28 = x27*(x22 - x6);
x29 = 45.0*x10;
x30 = -x29;
x31 = pow(x2, -4.5);
x32 = x0*x31;
x33 = 105.0*x32;
x34 = x30 + x33;
x35 = x34*x7;
x36 = x34*x8;
x37 = -x11;
x38 = x9*(x33 + x37);
x39 = x1*x31;
x40 = 105.0*x39;
x41 = x30 + x40;
x42 = x41*x9;
x43 = x27*y;
x44 = x41*x43;
x45 = x37 + x40;
x46 = x27*z;
x47 = x45*x46;
x48 = 315.0*x31;
x49 = pow(x2, -5.5);
x50 = 945.0*x49;
x51 = x0*x50;
x52 = x48 - x51;
x53 = x12*x52;
x54 = x1*x50;
x55 = x27*x9;
x56 = x55*(-x48 + x54);
x57 = 90.0*x10;
x58 = (x*x*x*x);
x59 = 105.0*x31;
x60 = x58*x59;
x61 = (y*y*y*y);
x62 = x59*x61;
x63 = -225.0*x10;
x64 = x50*x58;
x65 = -x64;
x66 = x*(1050.0*x32 + x63 + x65);
x67 = x30 + 630.0*x32 + x65;
x68 = x67*y;
x69 = x50*x61;
x70 = -x69;
x71 = y*(1050.0*x39 + x63 + x70);
x72 = x67*z;
x73 = 630.0*x39;
x74 = x30 + x70 + x73;
x75 = x74*z;
x76 = x27*(x29 + x69 - x73);
x77 = 1575.0*x31;
x78 = x0*x49;
x79 = pow(x2, -6.5);
x80 = x58*x79;
x81 = 10395.0*x80;
x82 = x77 - 9450.0*x78 + x81;
x83 = x7*x82;
x84 = x8*x82;
x85 = 5670.0*x78;
x86 = x48 + x81 - x85;
x87 = x86*x9;
x88 = x1*x49;
x89 = x61*x79;
x90 = 10395.0*x89;
x91 = x77 - 9450.0*x88 + x90;
x92 = x9*x91;
x93 = x43*x91;
x94 = x48 - 5670.0*x88 + x90;
x95 = x46*x94;
x96 = 14175.0*x49;
x97 = -x96;
x98 = x0*x79;
x99 = 103950.0*x98;
x100 = pow(x2, -7.5);
x101 = 135135.0*x100;
x102 = x101*x58;
x103 = -x102 + x97 + x99;
x104 = x103*x12;
x105 = x1*x79;
x106 = 103950.0*x105;
x107 = x101*x61;
x108 = x55*(-x106 + x107 + x96);
x109 = pow(x, 6);
x110 = 10395.0*x79;
x111 = x109*x110;
x112 = pow(y, 6);
x113 = x110*x112;
x114 = 11025.0*x31;
x115 = 99225.0*x49;
x116 = x101*x109;
x117 = -x116;
x118 = x117 + 218295.0*x80;
x119 = x*(-x0*x115 + x114 + x118);
x120 = 42525.0*x49;
x121 = -x0*x120 + x117 + x77 + 155925.0*x80;
x122 = x121*y;
x123 = x101*x112;
x124 = -x123;
x125 = x124 + 218295.0*x89;
x126 = y*(-x1*x115 + x114 + x125);
x127 = x121*z;
x128 = x1*x120;
x129 = 155925.0*x89;
x130 = x124 - x128 + x129 + x77;
x131 = x130*z;
x132 = x27*(x123 + x128 - x129 - x77);
x133 = pow(x2, -8.5);
x134 = x109*x133;
x135 = 2027025.0*x134;
x136 = x100*x58;
x137 = 2837835.0*x136;
x138 = -x115;
x139 = 1091475.0*x79;
x140 = x0*x139;
x141 = x138 + x140;
x142 = x135 - x137 + x141;
x143 = x142*x7;
x144 = x142*x8;
x145 = 467775.0*x79;
x146 = x0*x145;
x147 = 2027025.0*x136;
x148 = -x147;
x149 = x135 + x146 + x148 + x97;
x150 = x149*x9;
x151 = x1*x139;
x152 = x138 + x151;
x153 = x112*x133;
x154 = 2027025.0*x153;
x155 = x100*x61;
x156 = 2837835.0*x155;
x157 = x154 - x156;
x158 = x152 + x157;
x159 = x158*x9;
x160 = x1*x33;
x161 = x158*x43;
x162 = 2027025.0*x155;
x163 = -x162;
x164 = x1*x145 + x154 + x163 + x97;
x165 = x164*x46;
x166 = x133*x58;
x167 = 42567525.0*x166;
x168 = pow(x2, -9.5);
x169 = 34459425.0*x168;
x170 = x109*x169;
x171 = x0*x100;
x172 = 14189175.0*x171;
x173 = -x172;
x174 = x139 + x173;
x175 = x167 - x170 + x174;
x176 = x12*x175;
x177 = x1*x48;
x178 = -x1*x51;
x179 = x*(x177 + x178 + x34);
x180 = x0*x48;
x181 = y*(x178 + x180 + x41);
x182 = z*(x178 + x33 + x45);
x183 = x1*x100;
x184 = 14189175.0*x183;
x185 = x133*x61;
x186 = 42567525.0*x185;
x187 = x112*x169;
x188 = x55*(-x139 + x184 - x186 + x187);
x189 = -2835.0*x88;
x190 = x1*x98;
x191 = 10395.0*x190;
x192 = x189 + x191;
x193 = 945.0*x31;
x194 = -2835.0*x78;
x195 = x193 + x194;
x196 = x7*(x192 + x195);
x197 = x8*(x192 + x52);
x198 = x9*(x191 + x194 + x48 - x54);
x199 = 31185.0*x105;
x200 = x0*x1;
x201 = -8505.0*x49;
x202 = 31185.0*x98;
x203 = x201 + x202;
x204 = x12*(-x101*x200 + x199 + x203);
x205 = 3783780.0*x100;
x206 = pow(x, 8);
x207 = 2027025.0*x133;
x208 = x206*x207;
x209 = pow(y, 8);
x210 = x207*x209;
x211 = -893025.0*x49;
x212 = x169*x206;
x213 = -x212;
x214 = x*(72972900.0*x134 - 51081030.0*x136 + x211 + x213 + 13097700.0*x98);
x215 = 4365900.0*x98;
x216 = 56756700.0*x134 - 28378350.0*x136 + x213 + x215;
x217 = x138 + x216;
x218 = x217*y;
x219 = x169*x209;
x220 = -x219;
x221 = y*(13097700.0*x105 + 72972900.0*x153 - 51081030.0*x155 + x211 + x220);
x222 = x217*z;
x223 = 4365900.0*x105;
x224 = 56756700.0*x153;
x225 = 28378350.0*x155;
x226 = x220 + x223 + x224 - x225;
x227 = z*(x138 + x226);
x228 = x27*(x115 + x219 - x223 - x224 + x225);
x229 = 9823275.0*x79;
x230 = x109*x168;
x231 = pow(x2, -10.5);
x232 = 654729075.0*x231;
x233 = x206*x232;
x234 = 766215450.0*x166 - 170270100.0*x171 + x229 - 1240539300.0*x230 + x233;
x235 = x234*x7;
x236 = x234*x8;
x237 = 56756700.0*x171;
x238 = 425675250.0*x166;
x239 = 964863900.0*x230;
x240 = x9*(x139 + x233 - x237 + x238 - x239);
x241 = x112*x168;
x242 = x209*x232;
x243 = -170270100.0*x183 + 766215450.0*x185 + x229 - 1240539300.0*x241 + x242;
x244 = x243*x9;
x245 = 105.0*x10;
x246 = -x160 - 12.0*x5;
x247 = x243*x43;
x248 = 425675250.0*x185;
x249 = 964863900.0*x241;
x250 = x46*(x139 - 56756700.0*x183 + x242 + x248 - x249);
x251 = 1550674125.0*x168;
x252 = pow(x, 10)*x232;
x253 = x232*pow(y, 10);
x254 = 120.0*x10;
x255 = x61*x98;
x256 = 10395.0*x255;
x257 = x1*x85;
x258 = -x257;
x259 = x1*x81;
x260 = 62370.0*x190;
x261 = -x0*x107;
x262 = x260 + x261;
x263 = 31185.0*x89;
x264 = x263 - 17010.0*x88;
x265 = x*(x195 + x262 + x264);
x266 = x1*x96;
x267 = -x266;
x268 = x1*x99;
x269 = -x1*x102;
x270 = x*(x267 + x268 + x269 + x82);
x271 = x0*x96;
x272 = -x271;
x273 = y*(x261 + x268 + x272 + x91);
x274 = x189 + x260 + x269;
x275 = 31185.0*x80;
x276 = 17010.0*x78;
x277 = x193 + x275 - x276;
x278 = y*(x274 + x277);
x279 = z*(x194 + x262 + x94);
x280 = z*(x274 + x86);
x281 = 155925.0*x98;
x282 = 311850.0*x105;
x283 = -x120;
x284 = x1*x171;
x285 = -1351350.0*x284;
x286 = x283 + x285;
x287 = -405405.0*x155;
x288 = x0*x185;
x289 = 2027025.0*x288;
x290 = x287 + x289;
x291 = x7*(x281 + x282 + x286 + x290);
x292 = 155925.0*x105;
x293 = 311850.0*x98;
x294 = -405405.0*x136;
x295 = x1*x166;
x296 = 2027025.0*x295;
x297 = x294 + x296;
x298 = x7*(x286 + x292 + x293 + x297);
x299 = 187110.0*x105;
x300 = -810810.0*x284;
x301 = x8*(x203 + x290 + x299 + x300);
x302 = x8*(x103 + x285 + x292 + x296);
x303 = x9*(x106 - x107 + x281 + x285 + x289 + x97);
x304 = x201 + 187110.0*x98;
x305 = x9*(x199 + x297 + x300 + x304);
x306 = -2027025.0*x171;
x307 = x0*x61;
x308 = x133*x200;
x309 = x145 + 20270250.0*x308;
x310 = -4054050.0*x183 + 6081075.0*x185;
x311 = x12*(-x169*x307 + x306 + x309 + x310);
x312 = -2027025.0*x183;
x313 = x1*x58;
x314 = 6081075.0*x166;
x315 = -4054050.0*x171 + x314;
x316 = x12*(-x169*x313 + x309 + x312 + x315);
x317 = 15120.0*x49;
x318 = -x111;
x319 = -x259;
x320 = 270.0*x10 + x257;
x321 = -x113;
x322 = -x256;
x323 = -x179;
x324 = -x181;
x325 = -x182;
x326 = -x196;
x327 = -x197;
x328 = -x198;
x329 = -x204;
x330 = x1*x146;
x331 = x0*x154;
x332 = x0*x162;
x333 = x1*x78;
x334 = x1*x135;
x335 = x1*x147;
x336 = -6081075.0*x284;
x337 = x283 + x336;
x338 = 30405375.0*x185;
x339 = -x0*x187;
x340 = x0*x338 + x281 + x339;
x341 = 6081075.0*x153;
x342 = 1403325.0*x105 - 6081075.0*x155 + x341;
x343 = x*(x337 + x340 + x342);
x344 = -x1*x170;
x345 = -x1*x172 + x151;
x346 = x*(x1*x167 + x142 + x344 + x345);
x347 = y*(x0*x186 + x141 + x157 + x339 + x345);
x348 = 30405375.0*x166;
x349 = x1*x348 + x292 + x344;
x350 = 6081075.0*x134;
x351 = -6081075.0*x136 + x350 + 1403325.0*x98;
x352 = y*(x337 + x349 + x351);
x353 = z*(x164 + x336 + x340);
x354 = z*(x149 + x336 + x349);
x355 = 723647925.0*x168;
x356 = -x307*x355;
x357 = 212837625.0*x308;
x358 = x357 + 3274425.0*x79;
x359 = x112*x232;
x360 = x0*x359;
x361 = -103378275.0*x241 + x360;
x362 = x7*(x173 - 42567525.0*x183 + 127702575.0*x185 + x356 + x358 + x361);
x363 = x109*x232;
x364 = x1*x363;
x365 = -103378275.0*x230 + x364;
x366 = -x184;
x367 = -x313*x355 + x366;
x368 = x7*(127702575.0*x166 - 42567525.0*x171 + x358 + x365 + x367);
x369 = 91216125.0*x185;
x370 = 516891375.0*x168;
x371 = x145 + 91216125.0*x308;
x372 = x8*(-18243225.0*x183 + x306 - x307*x370 + x361 + x369 + x371);
x373 = x8*(x175 + x357 + x364 + x367);
x374 = x9*(x174 + x186 - x187 + x356 + x357 + x360 + x366);
x375 = 91216125.0*x166;
x376 = x9*(-18243225.0*x171 + x312 - x313*x370 + x365 + x371 + x375);
x377 = 16065.0*x49;
x378 = -x1*x276 - 360.0*x10;
x379 = 3918915.0*x100;
x380 = -x208;
x381 = 2338875.0*x79;
x382 = -x334;
x383 = -12600.0*x31 - x330;
x384 = -x210;
x385 = -x331;
x386 = x0*x155;
x387 = 810810.0*x386;
x388 = -x387;
x389 = x166*x61;
x390 = 2027025.0*x389;
x391 = x1*x136;
x392 = 810810.0*x391;
x393 = x390 - x392;
x394 = x58*x61;
x395 = -x169*x394;
x396 = x283 - 8108100.0*x284 + x395;
x397 = x294 + 12162150.0*x295;
x398 = 935550.0*x105 + x163;
x399 = x*(20270250.0*x288 + x293 + x396 + x397 + x398);
x400 = x287 + 12162150.0*x288;
x401 = x148 + 935550.0*x98;
x402 = y*(x282 + 20270250.0*x295 + x396 + x400 + x401);
x403 = z*(-4864860.0*x284 + x299 + x304 + x395 + x397 + x400);
x404 = x232*x394;
x405 = 344594250.0*x168;
x406 = -x307*x405 + x338 + x404;
x407 = -x313*x405 + x348;
x408 = x7*(-20270250.0*x171 - 20270250.0*x183 + 202702500.0*x308 + x381 + x406 + x407);
x409 = 206756550.0*x168;
x410 = x145 + 121621500.0*x308;
x411 = x8*(-12162150.0*x183 - x313*x409 + x315 + x406 + x410);
x412 = x9*(-12162150.0*x171 - x307*x409 + x310 + x404 + x407 + x410);
x413 = 17010.0*x49;
x414 = -x119;
x415 = -x270;
x416 = -x265;
x417 = -x122;
x418 = -x278;
x419 = -x126;
x420 = -x273;
x421 = -x127;
x422 = -x280;
x423 = -x131;
x424 = -x279;
x425 = x0*x248;
x426 = x0*x242;
x427 = x1*x237;
x428 = -x427;
x429 = x0*x249;
x430 = x1*x238;
x431 = x1*x233;
x432 = x1*x239;
x433 = -x143;
x434 = -x298;
x435 = -x161;
x436 = -x291;
x437 = -x144;
x438 = -x302;
x439 = -x165;
x440 = -x301;
x441 = -x150;
x442 = -x305;
x443 = -x159;
x444 = -x303;
x445 = -x176;
x446 = -x316;
x447 = -x311;
x448 = 1585133550.0*x168;
x449 = -x252;
x450 = -x431;
x451 = x427 + 992250.0*x49;
x452 = -x253;
x453 = -x426;
x454 = x123 + x385;
x455 = -841995.0*x190 - 2520.0*x31 - x390;
x456 = x116 + x382;
x457 = 4054050.0*x100;
x458 = x1*x134;
x459 = 1309770.0*x190 + 15120.0*x31;
x460 = x0*x153;
x461 = x1*x375;
x462 = x359*x58;
x463 = x0*x241;
x464 = 206756550.0*x463;
x465 = x462 - x464;
x466 = x370*x394;
x467 = -x466;
x468 = x283 - 36486450.0*x284 + x467;
x469 = x1*x230;
x470 = 206756550.0*x469;
x471 = -x470;
x472 = x0*x369;
x473 = x363*x61;
x474 = x472 + x473;
x475 = -x214;
x476 = -x346;
x477 = -x343;
x478 = -x218;
x479 = -x352;
x480 = -x221;
x481 = -x347;
x482 = -x222;
x483 = -x354;
x484 = -x227;
x485 = -x353;
x486 = -x235;
x487 = -x368;
x488 = -x247;
x489 = -x362;
x490 = -x236;
x491 = -x373;
x492 = -x250;
x493 = -x372;
x494 = -x240;
x495 = -x376;
x496 = -x244;
x497 = -x374;
x498 = 4189185.0*x100;
x499 = -2993760.0*x190 - 20160.0*x31 - x314*x61;
x500 = -x462;
x501 = x219 + x453;
x502 = 93243150.0*x284 + x466 + 141750.0*x49;
x503 = -x473;
x504 = x212 + x450;
x505 = 4324320.0*x100;
x506 = 1619592975.0*x168;
x507 = 1309458150.0*x231;
x508 = x1*x206;
x509 = -149999850.0*x284 + x467 - 1134000.0*x49;
x510 = x0*x209;
x511 = x168*x394;
x512 = -x350 + x470 + x503;
x513 = -x341 + x464 + x500;
x514 = -x399;
x515 = -x402;
x516 = -x403;
x517 = x112*x58;
x518 = x220 + x426;
x519 = -x251*x394 - 166216050.0*x284 - 226800.0*x49;
x520 = x109*x61;
x521 = x213 + x431;
x522 = -x408;
x523 = -x411;
x524 = -x412;
x525 = 1654052400.0*x168;
x526 = 1964187225.0*x231;
x527 = -x520*x526;
x528 = 316215900.0*x284 + 1360800.0*x49 + 2067565500.0*x511;
x529 = -x517*x526;
x530 = 1688511825.0*x168;
x531 = 2618916300.0*x231;
x532 = 3928374450.0*x231;
x533 = -648648000.0*x284 - 1814400.0*x49 - 5168913750.0*x511;
x534 = 1722971250.0*x168;
x535 = 3273645375.0*x231;
x536 = 6547290750.0*x231;
x537 = -x228;
x538 = -3*x399;
x539 = -3*x402;
x540 = -3*x403;
x541 = -3*x408;
x542 = -3*x411;
x543 = -3*x412;
#pragma omp atomic
F[0] += -x*x4*M[0] + x104*M[87] - x108*M[105] - x11*x12*M[13] + x119*M[83] + x122*M[84] + x126*M[111] + x127*M[85] + x131*M[112] - x132*M[104] + x143*M[120] + x144*M[121] + x150*M[123] + x159*M[156] + x161*M[147] + x165*M[148] + x176*M[168] + x179*M[37] + x181*M[40] + x182*M[41] - x188*M[201] + x19*M[9] + x196*M[61] + x197*M[62] + x198*M[66] + x204*M[94] + x21*M[10] + x214*M[164] + x218*M[165] + x221*M[209] + x222*M[166] + x227*M[210] - x228*M[200] + x235*M[220] + x236*M[221] + x24*M[15] + x240*M[223] + x244*M[275] + x247*M[264] + x25*M[11] + x250*M[265] + x26*M[16] + x265*M[93] + x270*M[86] + x273*M[98] + x278*M[89] + x279*M[99] - x28*M[12] + x280*M[90] + x291*M[134] + x298*M[125] + x301*M[135] + x302*M[126] + x303*M[141] + x305*M[130] + x311*M[186] + x316*M[175] + x343*M[185] + x346*M[167] + x347*M[192] + x35*M[20] + x352*M[170] + x353*M[193] + x354*M[171] + x36*M[21] + x362*M[247] + x368*M[225] + x372*M[248] + x373*M[226] + x374*M[256] + x376*M[230] + x38*M[23] + x399*M[174] - x4*y*M[1] - x4*z*M[2] + x402*M[179] + x403*M[180] + x408*M[234] + x411*M[235] + x412*M[241] + x42*M[30] + x44*M[25] + x47*M[26] + x53*M[38] - x56*M[45] + x6*x7*M[4] + x6*x8*M[5] + x6*x9*M[7] + x66*M[34] + x68*M[35] + x71*M[49] + x72*M[36] + x75*M[50] - x76*M[44] + x83*M[56] + x84*M[57] + x87*M[59] + x92*M[77] + x93*M[70] + x95*M[71] + (-x104 + x329)*M[96] + (x108 + x329)*M[107] + (x13 + x14)*M[3] + (x13 + x15)*M[6] + (x132 + x416)*M[106] + (x188 + x447)*M[203] + (-x19 + x28)*M[14] + (-x21 - x24)*M[17] + (x228 + x477)*M[202] + (-x25 - x26)*M[18] + (x323 - x66)*M[39] + (x323 + x76)*M[46] + (x324 - x68)*M[42] + (x324 - x71)*M[51] + (x325 - x72)*M[43] + (x325 - x75)*M[52] + (x326 - x83)*M[63] + (x326 - x93)*M[72] + (x327 - x84)*M[64] + (x327 - x95)*M[73] + (x328 - x87)*M[68] + (x328 - x92)*M[79] + (-x35 - x44)*M[27] + (-x36 - x47)*M[28] + (-x38 - x42)*M[32] + (x414 + x415)*M[88] + (x415 + x416)*M[95] + (x417 + x418)*M[91] + (x418 + x420)*M[100] + (x419 + x420)*M[113] + (x421 + x422)*M[92] + (x422 + x424)*M[101] + (x423 + x424)*M[114] + (x433 + x434)*M[127] + (x434 + x436)*M[136] + (x435 + x436)*M[149] + (x437 + x438)*M[128] + (x438 + x440)*M[137] + (x439 + x440)*M[150] + (x441 + x442)*M[132] + (x442 + x444)*M[143] + (x443 + x444)*M[158] + (x445 + x446)*M[177] + (x446 + x447)*M[188] + (x475 + x476)*M[169] + (x476 + x514)*M[176] + (x477 + x514)*M[187] + (x478 + x479)*M[172] + (x479 + x515)*M[181] + (x480 + x481)*M[211] + (x481 + x515)*M[194] + (x482 + x483)*M[173] + (x483 + x516)*M[182] + (x484 + x485)*M[212] + (x485 + x516)*M[195] + (x486 + x487)*M[227] + (x487 + x522)*M[236] + (x488 + x489)*M[266] + (x489 + x522)*M[249] + (x490 + x491)*M[228] + (x491 + x523)*M[237] + (x492 + x493)*M[267] + (x493 + x523)*M[250] + (x494 + x495)*M[232] + (x495 + x524)*M[243] + (x496 + x497)*M[277] + (x497 + x524)*M[258] + (-x53 + x56)*M[47] + (x104 - x108 + 2*x204)*M[109] + (x119 + x265 + 2*x270)*M[97] + (x122 + x273 + 2*x278)*M[102] + (x126 + 2*x273 + x278)*M[115] + (x127 + x279 + 2*x280)*M[103] + (x131 + 2*x279 + x280)*M[116] + (-x132 + 2*x265 + x270)*M[108] + (-x14 - x15 + 2.0*x3)*M[8] + (x143 + x291 + 2*x298)*M[138] + (x144 + x301 + 2*x302)*M[139] + (x150 + x303 + 2*x305)*M[145] + (x159 + 2*x303 + x305)*M[160] + (x160 + x20 + x23)*M[22] + (x161 + 2*x291 + x298)*M[151] + (x165 + 2*x301 + x302)*M[152] + (x176 + x311 + 2*x316)*M[190] + (2*x179 + x66 - x76)*M[48] + (2*x181 + x68 + x71)*M[53] + (2*x182 + x72 + x75)*M[54] + (-x188 + 2*x311 + x316)*M[205] + (2*x196 + x83 + x93)*M[74] + (2*x197 + x84 + x95)*M[75] + (2*x198 + x87 + x92)*M[81] + (x214 + 2*x346 + x399)*M[178] + (x218 + 2*x352 + x402)*M[183] + (x221 + 2*x347 + x402)*M[213] + (x222 + 2*x354 + x403)*M[184] + (x227 + 2*x353 + x403)*M[214] + (x235 + 2*x368 + x408)*M[238] + (x236 + 2*x373 + x411)*M[239] + (x240 + 2*x376 + x412)*M[245] + (x244 + 2*x374 + x412)*M[279] + (x247 + 2*x362 + x408)*M[268] + (x250 + 2*x372 + x411)*M[269] + (x343 + x346 + 2*x399)*M[189] + (2*x343 + x399 + x537)*M[204] + (x347 + x352 + 2*x402)*M[196] + (x353 + x354 + 2*x403)*M[197] + (x362 + x368 + 2*x408)*M[251] + (x372 + x373 + 2*x411)*M[252] + (x374 + x376 + 2*x412)*M[260] + (-x0*x57 + x16 + x60)*M[19] + (-x1*x57 + x16 + x62)*M[29] + (x111 + 4725.0*x32 - x58*x96 + x63)*M[55] + (x113 + 4725.0*x39 - x61*x96 + x63)*M[76] + (x132 - 3*x265 - 3*x270 + x414)*M[110] + (x177 + x258 + x259 + x67)*M[58] + (x180 + x256 + x258 + x74)*M[65] + (x188 - 3*x311 - 3*x316 + x445)*M[207] + (x228 - 3*x343 + x476 + x538)*M[206] + (-3*x273 - 3*x278 + x417 + x419)*M[117] + (-3*x279 - 3*x280 + x421 + x423)*M[118] + (-3*x291 - 3*x298 + x433 + x435)*M[153] + (-3*x301 - 3*x302 + x437 + x439)*M[154] + (-3*x303 - 3*x305 + x441 + x443)*M[162] + (-3*x346 + x475 + x477 + x538)*M[191] + (-3*x347 + x479 + x480 + x539)*M[215] + (-3*x352 + x478 + x481 + x539)*M[198] + (-3*x353 + x483 + x484 + x540)*M[216] + (-3*x354 + x482 + x485 + x540)*M[199] + (-3*x362 + x487 + x488 + x541)*M[270] + (-3*x368 + x486 + x489 + x541)*M[253] + (-3*x372 + x491 + x492 + x542)*M[271] + (-3*x373 + x490 + x493 + x542)*M[254] + (-3*x374 + x495 + x496 + x543)*M[281] + (-3*x376 + x494 + x497 + x543)*M[262] + (x0*x245 + x22 + x246 - x60)*M[24] + (x1*x245 + x17 + x246 - x62)*M[31] + (x121 + x267 + x330 + x334 - x335)*M[122] + (x130 + x272 + x330 + x331 - x332)*M[140] + (374220.0*x190 + x264 + x277 + x388 + x393)*M[129] + (x214 + 4*x343 + 4*x346 + 6*x399 + x537)*M[208] + (x218 + x221 + 4*x347 + 4*x352 + 6*x402)*M[217] + (x222 + x227 + 4*x353 + 4*x354 + 6*x403)*M[218] + (x235 + x247 + 4*x362 + 4*x368 + 6*x408)*M[272] + (x236 + x250 + 4*x372 + 4*x373 + 6*x411)*M[273] + (x240 + x244 + 4*x374 + 4*x376 + 6*x412)*M[283] + (-x109*x205 + x114 + x208 - 396900.0*x78 + 2182950.0*x80)*M[119] + (-x112*x205 + x114 + x210 - 396900.0*x88 + 2182950.0*x89)*M[155] + (49116375.0*x105 + 1277025750.0*x153 - 425675250.0*x155 - x209*x251 + x211 + x253)*M[274] + (1277025750.0*x134 - 425675250.0*x136 - x206*x251 + x211 + x252 + 49116375.0*x98)*M[219] + (x141 + x226 + x425 + x426 + x428 - x429)*M[255] + (x152 + x216 + x428 + x430 + x431 - x432)*M[222] + (-x177 + x317*x58 + x318 + x319 - 5355.0*x32 + x320)*M[60] + (-x180 + x317*x61 + x320 + x321 + x322 - 5355.0*x39)*M[78] + (182432250.0*x288 + x342 + x401 + x461 + x465 + x468)*M[240] + (182432250.0*x295 + x351 + x398 + x468 + x471 + x474)*M[229] + (-x0*x254 - x1*x254 + 210.0*x1*x32 + 24.0*x5 + x60 + x62)*M[33] + (x113 + 20790.0*x255 + x259 + 1260.0*x32 - x377*x61 + x378 + 6300.0*x39 + x65)*M[80] + (x0*x156 - x275 + x392 + x454 + x455 + 31185.0*x78 + 59535.0*x88 - 187110.0*x89)*M[142] + (-x0*x193 - x1*x193 + x319 + x322 + 11340.0*x333 + x57 + x64 + x69)*M[67] + (x1*x137 - x263 + x387 + x455 + x456 + 59535.0*x78 - 187110.0*x80 + 31185.0*x88)*M[131] + (20790.0*x1*x80 + x111 + x256 + 6300.0*x32 - x377*x58 + x378 + 1260.0*x39 + x70)*M[69] + (x109*x379 + x266 + x335 + x380 - x381*x58 + x382 + x383 + 439425.0*x78)*M[124] + (x112*x379 + x271 + x332 - x381*x61 + x383 + x384 + x385 + 439425.0*x88)*M[157] + (-53482275.0*x105 - x140 - 1333782450.0*x153 + 454053600.0*x155 + x209*x448 - x425 + x429 + x451 + x452 + x453)*M[276] + (-1333782450.0*x134 + 454053600.0*x136 - x151 + x206*x448 - x430 + x432 + x449 + x450 + x451 - 53482275.0*x98)*M[224] + (-x1*x275 + 720.0*x10 - x202*x61 + x318 - 7560.0*x32 + x321 + 34020.0*x333 - 7560.0*x39 + x413*x58 + x413*x61)*M[82] + (-x112*x457 + x210 + x275 - 4864860.0*x386 + x393 + x459 + 4054050.0*x460 - 45360.0*x78 - 498960.0*x88 + 2525985.0*x89)*M[159] + (-5769225.0*x105 + x147 - 62837775.0*x153 + 34459425.0*x155 - 608107500.0*x288 - x461 + 1171620450.0*x463 + x500 + x501 + x502 - 2027025.0*x98)*M[257] + (-2027025.0*x105 - 62837775.0*x134 + 34459425.0*x136 + x162 - 608107500.0*x295 + 1171620450.0*x469 - x472 + x502 + x503 + x504 - 5769225.0*x98)*M[231] + (x118 + x125 + 1683990.0*x190 + 5040.0*x31 + x331 + x334 - 3648645.0*x386 + 4054050.0*x389 - 3648645.0*x391 - 90720.0*x78 - 90720.0*x88)*M[144] + (-x0*x341 + x112*x498 + x384 + 8513505.0*x386 + 4459455.0*x391 + x456 + x499 + 136080.0*x78 - 249480.0*x80 + 589680.0*x88 - 2744280.0*x89)*M[161] + (-x0*x381 - x1*x381 + 8108100.0*x136 + 8108100.0*x155 + 72972900.0*x284 - 273648375.0*x288 - 273648375.0*x295 + 85050.0*x49 + 1033782750.0*x511 + x512 + x513)*M[242] + (-x1*x350 + x109*x498 + x380 + 4459455.0*x386 + 8513505.0*x391 + x454 + x499 + 589680.0*x78 - 2744280.0*x80 + 136080.0*x88 - 249480.0*x89)*M[146] + (-x109*x457 + x208 + x263 + x388 + x390 - 4864860.0*x391 + 4054050.0*x458 + x459 - 498960.0*x78 + 2525985.0*x80 - 45360.0*x88)*M[133] + (3118500.0*x105 + 1396620225.0*x134 - 488513025.0*x136 + x163 - x206*x506 + x252 + 1033782750.0*x295 - 2136484350.0*x469 + x474 + x507*x508 + x509 + 59251500.0*x98)*M[233] + (59251500.0*x105 + x148 + 1396620225.0*x153 - 488513025.0*x155 - x209*x506 + x253 + 1033782750.0*x288 + x461 + x462 - 2136484350.0*x463 + x507*x510 + x509 + 3118500.0*x98)*M[278] + (68918850.0*x134 - 42567525.0*x136 - 10135125.0*x155 + x223 + 364864500.0*x288 + 881755875.0*x295 + x341 + x465 - 1378377000.0*x469 + x507*x520 + x519 + x521 + 8108100.0*x98)*M[244] + (-67359600.0*x105 + 12162150.0*x136 - 1465539075.0*x153 + 531080550.0*x155 + x209*x525 - 1915538625.0*x288 - 456080625.0*x295 + x452 + 3514861350.0*x463 - x510*x526 + x512 + x528 + x529 - 7484400.0*x98)*M[280] + (-7484400.0*x105 - 1465539075.0*x134 + 531080550.0*x136 + 12162150.0*x155 + x206*x525 - 456080625.0*x288 - 1915538625.0*x295 + x449 + 3514861350.0*x469 - x508*x526 + x513 + x527 + x528 - 67359600.0*x98)*M[246] + (8108100.0*x105 - 10135125.0*x136 + 68918850.0*x153 - 42567525.0*x155 + x215 + 881755875.0*x288 + 364864500.0*x295 + x350 - 1378377000.0*x463 + x471 + x473 + x507*x517 + x518 + x519)*M[259] + (-x109*x505 - x112*x505 + 5987520.0*x190 + x208 + x210 + 40320.0*x31 - 12972960.0*x386 + 12162150.0*x389 - 12972960.0*x391 + 8108100.0*x458 + 8108100.0*x460 - 725760.0*x78 + 2993760.0*x80 - 725760.0*x88 + 2993760.0*x89)*M[163] + (-12474000.0*x105 - 74999925.0*x134 + 52702650.0*x136 - 74999925.0*x153 + 52702650.0*x155 + 332432100.0*x284 - 1246620375.0*x288 - 1246620375.0*x295 + 1585133550.0*x463 + 1585133550.0*x469 + 453600.0*x49 + x501 + x504 + 3101348250.0*x511 + x527 + x529 - 12474000.0*x98)*M[261] + (19958400.0*x105 + 1540539000.0*x134 - 583783200.0*x136 + 81081000.0*x153 - 64864800.0*x155 - x206*x530 + x252 + 1702701000.0*x288 + 3162159000.0*x295 - 1791890100.0*x463 - 5099994900.0*x469 + x508*x531 + x517*x531 + x518 + x520*x532 + x533 + 79833600.0*x98)*M[263] + (79833600.0*x105 + 81081000.0*x134 - 64864800.0*x136 + 1540539000.0*x153 - 583783200.0*x155 - x209*x530 + x253 + 3162159000.0*x288 + 1702701000.0*x295 - 5099994900.0*x463 - 1791890100.0*x469 + x510*x531 + x517*x532 + x520*x531 + x521 + x533 + 19958400.0*x98)*M[282] + (-99792000.0*x105 - 1621620000.0*x134 + 648648000.0*x136 - 1621620000.0*x153 + 648648000.0*x155 + x206*x534 + x209*x534 + 1297296000.0*x284 - 4864860000.0*x288 - 4864860000.0*x295 + x449 + x452 + 6891885000.0*x463 + 6891885000.0*x469 + 3628800.0*x49 - x508*x535 - x510*x535 + 10337827500.0*x511 - x517*x536 - x520*x536 - 99792000.0*x98)*M[284];

}

void P2M(double x, double y, double z, double q, double * M, int order) {
switch (order) {
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
  case 6:
    P2M_6(x, y, z, q, M);
    break;
  case 7:
    P2M_7(x, y, z, q, M);
    break;
  case 8:
    P2M_8(x, y, z, q, M);
    break;
  case 9:
    P2M_9(x, y, z, q, M);
    break;
  case 10:
    P2M_10(x, y, z, q, M);
    break;
  }
}
void M2M(double x, double y, double z, double * M, double * Ms, int order) {
switch (order) {
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
  case 6:
    M2M_6(x, y, z, M, Ms);
    break;
  case 7:
    M2M_7(x, y, z, M, Ms);
    break;
  case 8:
    M2M_8(x, y, z, M, Ms);
    break;
  case 9:
    M2M_9(x, y, z, M, Ms);
    break;
  case 10:
    M2M_10(x, y, z, M, Ms);
    break;
  }
}
void M2L(double x, double y, double z, double * M, double * L, int order) {
switch (order) {
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
  case 6:
    M2L_6(x, y, z, M, L);
    break;
  case 7:
    M2L_7(x, y, z, M, L);
    break;
  case 8:
    M2L_8(x, y, z, M, L);
    break;
  case 9:
    M2L_9(x, y, z, M, L);
    break;
  case 10:
    M2L_10(x, y, z, M, L);
    break;
  }
}
void L2L(double x, double y, double z, double * L, double * Ls, int order) {
switch (order) {
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
  case 6:
    L2L_6(x, y, z, L, Ls);
    break;
  case 7:
    L2L_7(x, y, z, L, Ls);
    break;
  case 8:
    L2L_8(x, y, z, L, Ls);
    break;
  case 9:
    L2L_9(x, y, z, L, Ls);
    break;
  case 10:
    L2L_10(x, y, z, L, Ls);
    break;
  }
}
void L2P(double x, double y, double z, double * L, double * F, int order) {
switch (order) {
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
  case 6:
    L2P_6(x, y, z, L, F);
    break;
  case 7:
    L2P_7(x, y, z, L, F);
    break;
  case 8:
    L2P_8(x, y, z, L, F);
    break;
  case 9:
    L2P_9(x, y, z, L, F);
    break;
  case 10:
    L2P_10(x, y, z, L, F);
    break;
  }
}
void M2P(double x, double y, double z, double * M, double * F, int order) {
switch (order) {
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
  case 6:
    M2P_6(x, y, z, M, F);
    break;
  case 7:
    M2P_7(x, y, z, M, F);
    break;
  case 8:
    M2P_8(x, y, z, M, F);
    break;
  case 9:
    M2P_9(x, y, z, M, F);
    break;
  case 10:
    M2P_10(x, y, z, M, F);
    break;
  }
}
