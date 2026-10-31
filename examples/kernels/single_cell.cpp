// Far-field check for generated operators: S2M -> M2M -> M2L -> L2L -> L2P
// against direct summation, for two source cells and two target cells.
//
// Conventions of the generated operators (all accumulate with +=):
//   S2M(centre - source)  M2M(parent - child)  M2L(target - source)
//   L2L(child - parent)   L2P(body - centre)   P2P(target - source)
// Output F = {potential, -grad potential}.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "yukawa_cpp.h"

static int64_t state = 12345;
static double uniform() {  // in [-1, 1); the same generator is used by single_cell.f90
  state = (state * 1103515245 + 12345) % 2147483648;
  return 2.0 * double(state) / 2147483648.0 - 1.0;
}

static size_t nterms(int p) { return size_t(p + 1) * (p + 2) * (p + 3) / 6; }

struct Pt { double x, y, z; };

static Pt around(Pt c, double r) { return {c.x + r * uniform(), c.y + r * uniform(), c.z + r * uniform()}; }
static Pt sub(Pt a, Pt b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static Pt add(Pt a, Pt b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }

int main() {
  const Pt cs = {0, 0, 0}, ct = {2.5, 0.5, -0.5};  // parent source and target cell centres
  const Pt off = {0.1, 0.1, 0.1};
  const Pt skid[2] = {add(cs, off), sub(cs, off)}, tkid[2] = {add(ct, off), sub(ct, off)};
  const int ns = 20, nt = 8;

  std::vector<Pt> src(ns), tgt(nt);
  std::vector<double> q(ns);
  std::vector<int> sc(ns), tc(nt);
  for (int i = 0; i < ns; i++) { sc[i] = i % 2; src[i] = around(skid[sc[i]], 0.1); q[i] = uniform(); }
  for (int i = 0; i < nt; i++) { tc[i] = i % 2; tgt[i] = around(tkid[tc[i]], 0.1); }

  // Direct reference through the generated P2P.
  std::vector<double> ref(4 * nt, 0.0);
  for (int i = 0; i < nt; i++)
    for (int j = 0; j < ns; j++) {
      Pt d = sub(tgt[i], src[j]);
      P2P(d.x, d.y, d.z, &q[j], &ref[4 * i]);
    }

  std::printf("kernel exp(-r)/r, %d sources, %d targets\n order  max relative error in potential\n", ns, nt);
  for (int p = FMMGEN_MINORDER; p <= FMMGEN_MAXORDER; p += 2) {
    const size_t n = nterms(p);
    std::vector<double> M[2], Mp(n, 0.0), L(n, 0.0), Lk[2];
    for (int k = 0; k < 2; k++) { M[k].assign(n, 0.0); Lk[k].assign(n, 0.0); }
    for (int j = 0; j < ns; j++) {
      Pt d = sub(skid[sc[j]], src[j]);
      S2M(d.x, d.y, d.z, &q[j], M[sc[j]].data(), p);
    }
    for (int k = 0; k < 2; k++) { Pt d = sub(cs, skid[k]); M2M(d.x, d.y, d.z, M[k].data(), Mp.data(), p); }
    { Pt d = sub(ct, cs); M2L(d.x, d.y, d.z, Mp.data(), L.data(), p); }
    for (int k = 0; k < 2; k++) { Pt d = sub(tkid[k], ct); L2L(d.x, d.y, d.z, L.data(), Lk[k].data(), p); }
    double err = 0.0, scale = 0.0;
    for (int i = 0; i < nt; i++) {
      double F[4] = {0, 0, 0, 0};
      Pt d = sub(tgt[i], tkid[tc[i]]);
      L2P(d.x, d.y, d.z, Lk[tc[i]].data(), F, p);
      err = std::fmax(err, std::fabs(F[0] - ref[4 * i]));
      scale = std::fmax(scale, std::fabs(ref[4 * i]));
    }
    std::printf(" %5d  %.3e\n", p, err / scale);
  }
  return 0;
}
