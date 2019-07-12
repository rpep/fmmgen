#include "operators.h"
#include "math.h"
void P2M_1(double x, double y, double z, double q, double * M) {
M[0] += -q*x;
M[1] += -q*y;
M[2] += -q*z;
}
void M2M_1(double x, double y, double z, double * M, double * Ms) {
Ms[0] += M[0];
Ms[1] += M[1];
Ms[2] += M[2];
}

void M2L_1(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double x0;
x0 = 1.0/pow(R, 3);
L[0] += -x*x0*M[0] - x0*y*M[1] - x0*z*M[2];
}

void L2L_1(double x, double y, double z, double * L, double * Ls) {
Ls[0] += L[0];
}

void L2P_1(double x, double y, double z, double * L, double * F) {
F[0] += 0;
F[1] += 0;
F[2] += 0;
}

void M2P_1(double x, double y, double z, double * M, double * F) {
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
x0 = pow(x, 2);
x1 = pow(y, 2);
x2 = pow(z, 2);
x3 = x0 + x1 + x2;
x4 = 1.0*pow(x3, -1.5);
x5 = 3.0*pow(x3, -2.5);
x6 = x5*M[1];
x7 = x*y;
x8 = x5*M[2];
x9 = x8*z;
x10 = x5*M[0];
F[0] += -x*x9 - x0*x10 + x4*M[0] - x6*x7;
F[1] += -x1*x6 - x10*x7 + x4*M[1] - x9*y;
F[2] += -x*x10*z - x2*x8 + x4*M[2] - x6*y*z;
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
M[3] += pow(x, 2)*x2;
M[4] += x0*y;
M[5] += x0*z;
M[6] += x2*pow(y, 2);
M[7] += x1*z;
M[8] += x2*pow(z, 2);
}
void M2M_2(double x, double y, double z, double * M, double * Ms) {
Ms[0] += M[0];
Ms[1] += M[1];
Ms[2] += M[2];
Ms[3] += x*M[0] + M[3];
Ms[4] += x*M[1] + y*M[0] + M[4];
Ms[5] += x*M[2] + z*M[0] + M[5];
Ms[6] += y*M[1] + M[6];
Ms[7] += y*M[2] + z*M[1] + M[7];
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
x0 = 1.0/pow(R, 3);
x1 = x*M[0];
x2 = y*M[1];
x3 = z*M[2];
x4 = 3.0/pow(R, 5);
x5 = x*x4;
x6 = x4*y;
x7 = -x0;
x8 = pow(x, 2)*x4 + x7;
x9 = x4*pow(y, 2) + x7;
x10 = x4*pow(z, 2) + x7;
x11 = x4*z;
L[0] += -x0*x1 - x0*x2 - x0*x3 + x10*M[8] + x5*y*M[4] + x5*z*M[5] + x6*z*M[7] + x8*M[3] + x9*M[6];
L[1] += x2*x5 + x3*x5 + x8*M[0];
L[2] += x1*x6 + x3*x6 + x9*M[1];
L[3] += x1*x11 + x10*M[2] + x11*x2;
}

void L2L_2(double x, double y, double z, double * L, double * Ls) {
Ls[0] += x*L[1] + y*L[2] + z*L[3] + L[0];
Ls[1] += L[1];
Ls[2] += L[2];
Ls[3] += L[3];
}

void L2P_2(double x, double y, double z, double * L, double * F) {
F[0] += -L[1];
F[1] += -L[2];
F[2] += -L[3];
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
x0 = pow(x, 2);
x1 = pow(y, 2);
x2 = pow(z, 2);
x3 = x0 + x1 + x2;
x4 = 1.0*pow(x3, -1.5);
x5 = pow(x3, -2.5);
x6 = x5*y;
x7 = 3.0*x6;
x8 = x5*z;
x9 = 3.0*x8;
x10 = x*x5;
x11 = 3.0*x10;
x12 = x11*y;
x13 = z*M[2];
x14 = 3.0*x5;
x15 = 15.0*pow(x3, -3.5);
x16 = x*x15;
x17 = x16*y*z;
x18 = x15*y;
x19 = x0*x18;
x20 = x15*z;
x21 = x0*x20;
x22 = x1*x16;
x23 = x16*x2;
x24 = x1*x20;
x25 = x18*x2;
F[0] += -x0*x14*M[0] - x11*x13 - x12*M[1] + x17*M[7] + x19*M[4] + x21*M[5] + x4*M[0] - x7*M[4] - x9*M[5] - (x11 - x22)*M[6] - (x11 - x23)*M[8] - (-pow(x, 3)*x15 + 9.0*x10)*M[3];
F[1] += -x1*x14*M[1] - x11*M[4] - x12*M[0] - x13*x7 + x17*M[5] + x22*M[4] + x24*M[7] + x4*M[1] - x9*M[7] - (-x19 + x7)*M[3] - (-x25 + x7)*M[8] - (-x15*pow(y, 3) + 9.0*x6)*M[6];
F[2] += -x11*z*M[0] - x11*M[5] - x14*x2*M[2] + x17*M[4] + x23*M[5] + x25*M[7] + x4*M[2] - x7*z*M[1] - x7*M[7] - (-x21 + x9)*M[3] - (-x24 + x9)*M[6] - (-x15*pow(z, 3) + 9.0*x8)*M[8];
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
x3 = pow(x, 2);
x4 = (1.0/2.0)*q;
x5 = x0*y;
x6 = pow(y, 2);
x7 = pow(z, 2);
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
M[9] += -pow(x, 3)*x8;
M[10] += -x1*x9;
M[11] += -x2*x9;
M[12] += -x10*x6;
M[13] += -x5*z;
M[14] += -x10*x7;
M[15] += -x8*pow(y, 3);
M[16] += -1.0/2.0*x2*x6;
M[17] += -1.0/2.0*x1*x7;
M[18] += -x8*pow(z, 3);
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
x6 = (1.0/2.0)*pow(x, 2);
x7 = pow(y, 2);
x8 = (1.0/2.0)*M[0];
x9 = pow(z, 2);
x10 = (1.0/2.0)*x7;
x11 = (1.0/2.0)*x9;
Ms[0] += M[0];
Ms[1] += M[1];
Ms[2] += M[2];
Ms[3] += x0 + M[3];
Ms[4] += x1 + x2 + M[4];
Ms[5] += x3 + z*M[0] + M[5];
Ms[6] += x4 + M[6];
Ms[7] += x5 + z*M[1] + M[7];
Ms[8] += z*M[2] + M[8];
Ms[9] += x*M[3] + x6*M[0] + M[9];
Ms[10] += x*M[4] + x0*y + x6*M[1] + y*M[3] + M[10];
Ms[11] += x*M[5] + x0*z + x6*M[2] + z*M[3] + M[11];
Ms[12] += x*M[6] + x1*y + x7*x8 + y*M[4] + M[12];
Ms[13] += x*M[7] + x1*z + x2*z + x3*y + y*M[5] + z*M[4] + M[13];
Ms[14] += x*M[8] + x3*z + x8*x9 + z*M[5] + M[14];
Ms[15] += x10*M[1] + y*M[6] + M[15];
Ms[16] += x10*M[2] + x4*z + y*M[7] + z*M[6] + M[16];
Ms[17] += x11*M[1] + x5*z + y*M[8] + z*M[7] + M[17];
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
double x39;
double x40;
double x41;
double x42;
x0 = 1.0/pow(R, 3);
x1 = x*M[0];
x2 = y*M[1];
x3 = z*M[2];
x4 = pow(R, -5);
x5 = 3.0*x4;
x6 = x*x5;
x7 = y*M[4];
x8 = z*M[5];
x9 = x5*y;
x10 = z*M[7];
x11 = 15.0/pow(R, 7);
x12 = x*x11*y;
x13 = -x0;
x14 = pow(x, 2);
x15 = x13 + x14*x5;
x16 = pow(y, 2);
x17 = x13 + x16*x5;
x18 = pow(z, 2);
x19 = x13 + x18*x5;
x20 = 9.0*x4;
x21 = -x11*x14;
x22 = x*(x20 + x21);
x23 = x21 + x5;
x24 = x23*y;
x25 = x11*x16;
x26 = -x25;
x27 = y*(x20 + x26);
x28 = x23*z;
x29 = z*(x26 + x5);
x30 = x11*x18;
x31 = z*(x20 - x30);
x32 = -x5;
x33 = x25 + x32;
x34 = 1.0*x;
x35 = x33*x34;
x36 = x30 + x32;
x37 = x34*x36;
x38 = 1.0*x36;
x39 = x38*y;
x40 = x5*z;
x41 = x11*z;
x42 = x*x41;
L[0] += -x0*x1 - x0*x2 - x0*x3 + x10*x9 - x12*z*M[13] + x15*M[3] + x17*M[6] + x19*M[8] + x22*M[9] + x24*M[10] + x27*M[15] + x28*M[11] + x29*M[16] + x31*M[18] - x35*M[12] - x37*M[14] - x39*M[17] + x6*x7 + x6*x8;
L[1] += -x10*x12 + x15*M[0] + x2*x6 + x22*M[3] + x24*M[4] + x28*M[5] + x3*x6 - x35*M[6] - x37*M[8];
L[2] += x1*x9 - x12*x8 + x17*M[1] + x24*M[3] + x27*M[6] + x29*M[7] + x3*x9 - x35*M[4] - x39*M[8];
L[3] += x1*x40 + x19*M[2] + x2*x40 + x28*M[3] + x29*M[6] + x31*M[8] - x37*M[5] - x39*M[7] - x42*x7;
L[4] += x22*M[0] + x24*M[1] + x28*M[2];
L[5] += -x12*x3 + x24*M[0] - x35*M[1];
L[6] += -x2*x42 + x28*M[0] - x37*M[2];
L[7] += -1.0*x1*x33 + x27*M[1] + x29*M[2];
L[8] += -x1*x41*y + x29*M[1] - x39*M[2];
L[9] += -x1*x38 - x2*x38 + x31*M[2];
}

void L2L_3(double x, double y, double z, double * L, double * Ls) {
double x0;
double x1;
double x2;
x0 = y*L[5];
x1 = z*L[6];
x2 = z*L[8];
Ls[0] += (1.0/2.0)*pow(x, 2)*L[4] + x*x0 + x*x1 + x*L[1] + x2*y + (1.0/2.0)*pow(y, 2)*L[7] + y*L[2] + (1.0/2.0)*pow(z, 2)*L[9] + z*L[3] + L[0];
Ls[1] += x*L[4] + x0 + x1 + L[1];
Ls[2] += x*L[5] + x2 + y*L[7] + L[2];
Ls[3] += x*L[6] + y*L[8] + z*L[9] + L[3];
Ls[4] += L[4];
Ls[5] += L[5];
Ls[6] += L[6];
Ls[7] += L[7];
Ls[8] += L[8];
Ls[9] += L[9];
}

void L2P_3(double x, double y, double z, double * L, double * F) {
F[0] += -x*L[4] - y*L[5] - z*L[6] - L[1];
F[1] += -x*L[5] - y*L[7] - z*L[8] - L[2];
F[2] += -x*L[6] - y*L[8] - z*L[9] - L[3];
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
x0 = pow(x, 2);
x1 = pow(y, 2);
x2 = pow(z, 2);
x3 = x0 + x1 + x2;
x4 = 1.0*pow(x3, -1.5);
x5 = pow(x3, -2.5);
x6 = 3.0*x5;
x7 = x6*y;
x8 = x6*z;
x9 = x*x6;
x10 = x9*y;
x11 = z*M[2];
x12 = pow(x3, -3.5);
x13 = 15.0*x12;
x14 = x13*y;
x15 = z*M[13];
x16 = x*x13;
x17 = x16*y;
x18 = x17*z;
x19 = x0*x13;
x20 = x19*y;
x21 = x19*z;
x22 = 105.0*pow(x3, -4.5);
x23 = x22*y;
x24 = x0*x23;
x25 = 9.0*x5;
x26 = -x19;
x27 = x1*x13;
x28 = -x27;
x29 = x28 + x6;
x30 = x13*x2;
x31 = -x30;
x32 = x31 + x6;
x33 = pow(x, 3);
x34 = x*x27;
x35 = x*x30;
x36 = x*x12;
x37 = x22*x33;
x38 = x*M[9];
x39 = -45.0*x36;
x40 = x37 + x39;
x41 = y*M[10];
x42 = z*M[11];
x43 = x*x22;
x44 = x1*x43;
x45 = 1.0*x;
x46 = x45*M[12];
x47 = x2*x43;
x48 = x16 - x47;
x49 = x45*M[14];
x50 = y*M[15];
x51 = 1.0*y*M[17];
x52 = z*M[16];
x53 = z*M[18];
x54 = x27*z;
x55 = x26 + x6;
x56 = pow(y, 3);
x57 = x30*y;
x58 = x12*y;
x59 = 45.0*x58;
x60 = x22*x56;
x61 = -x59;
x62 = x2*x23;
x63 = x14 - x62;
x64 = -x14 + x24;
x65 = pow(z, 3);
x66 = x12*z;
x67 = 45.0*x66;
x68 = x22*x65;
x69 = x67 - x68;
x70 = -x67;
x71 = x22*z;
x72 = x0*x71;
x73 = x13*z;
x74 = x1*x71;
x75 = -x73;
x76 = x72 + x75;
F[0] += -x0*x6*M[0] - x10*M[1] - x11*x9 + x14*x15 - x15*x24 + x18*M[7] + x20*M[4] + x21*M[5] - x29*M[12] - x32*M[14] - x38*(-75.0*x36 + x37) + x4*M[0] - x40*x41 - x40*x42 + x46*(x16 - x44) + x48*x49 + x48*x51 - x50*(x39 + x44) - x52*(-x16 + x44) - x53*(x39 + x47) - x7*M[4] - x8*M[5] - (x25 + x26)*M[9] - (-x34 + x9)*M[6] - (-x35 + x9)*M[8] - (x*x25 - x13*x33)*M[3];
F[1] += -x1*x6*M[1] - x10*M[0] - x11*x7 + x15*x16 - x15*x44 + x18*M[5] - x32*M[17] + x34*M[4] - x38*(x24 + x61) + x4*M[1] - x41*x64 - x42*x64 + x46*(x59 - x60) + x49*x63 - x50*(-75.0*x58 + x60) + x51*x63 - x52*(x60 + x61) - x53*(x61 + x62) + x54*M[7] - x55*M[10] - x8*M[7] - x9*M[4] - (-x20 + x7)*M[3] - (x25 + x28)*M[15] - (-x57 + x7)*M[8] - (-x13*x56 + x25*y)*M[6];
F[2] += x17*M[13] + x18*M[4] - x2*x6*M[2] - x29*M[16] + x35*M[5] - x38*(x70 + x72) + x4*M[2] - x41*x76 - x42*x76 + x46*(x73 - x74) - x47*y*M[13] + x49*x69 - x50*(x70 + x74) + x51*x69 - x52*(x74 + x75) - x53*(-75.0*x66 + x68) - x55*M[11] + x57*M[7] - x7*z*M[1] - x7*M[7] - x9*z*M[0] - x9*M[5] - (-x21 + x8)*M[3] - (x25 + x31)*M[18] - (-x54 + x8)*M[6] - (-x13*x65 + x25*z)*M[8];
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
  }
}
