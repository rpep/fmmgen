// main.c
#include "stdio.h"
#include "stdlib.h"

#include "operators.h"

unsigned int Nterms(const unsigned int order) {
  int nterms = 0;
  for(unsigned int n = 0; n <= order; n++) {
    nterms += (n*(n + 1)) / 2;
  }
  return nterms;
}

int main(void) {
  // Expansion Order:
  const unsigned int order = 2;
  // Number of terms in expansion at order:
  unsigned int N = Nterms(order);
  // Single multipole array set to zero:

  double *M = (double *) calloc(N, sizeof(double));

  double x, y, z, q;
  x = 1.1;
  y = 2.4;
  z = -1.2;
  q = 1e-5;
   
  // Calculate P2M expansion up to order
  P2M(x, y, z, q, M, order);

  for(unsigned int i = 0; i < N; i++) {
    printf("M[%d] = %g\n", i, M[i]);
  }


  free(M);
  return 0;
}
