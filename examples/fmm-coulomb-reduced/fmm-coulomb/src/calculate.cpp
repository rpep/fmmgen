#include "calculate.hpp"
#include "operators.h"
#include "tree.hpp"
#include "utils.hpp"
#include <iostream>
#include <stack>

inline void P2P(double x, double y, double z, double q, double *F) {
  double R = sqrt(x * x + y * y + z * z);
  double R3 = R * R * R;
  F[0] += q / R;
  F[1] += q * x / R3;
  F[2] += q * y / R3;
  F[3] += q * z / R3;
}

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  size_t cell, size_t ncrit,
                  size_t exporder) {
#pragma omp parallel for schedule(dynamic, 32)
  for(size_t c = 0; c < cells.size(); c++) {
    if (!(cells[c].nleaf >= ncrit)) {
      /*    for (size_t octant = 0; octant < 8; octant++) {
      if (cells[cell].nchild & (1 << octant)) {
        evaluate_P2M(particles, cells, cells[cell].child[octant], ncrit, exporder);
      }
    }
  }
  else { */
      for (size_t i = 0; i < (cells[c].nleaf); i++) {
	size_t l = cells[c].leaf[i];
	double dx = particles[l].x - cells[c].x;
	double dy = particles[l].y - cells[c].y;
	double dz = particles[l].z - cells[c].z;
	P2M(dx, dy, dz, particles[l].q, cells[c].M.data(), exporder);
      }
    }
  }
}

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  size_t exporder) {
  /*
  evaluate_M2M(particles, cells)

  This function evaluates the multipole to
  multipole kernel. It does this by working up the
  tree from the leaf nodes, which is possible
  by iterating backwards through the nodes because
  of the way the tree is constructed.
  */

  for (size_t i = cells.size() - 1; i > 0; i--) {
    size_t p = cells[i].parent;
    double dx = cells[p].x - cells[i].x;
    double dy = cells[p].y - cells[i].y;
    double dz = cells[p].z - cells[i].z;
    M2M(dx, dy, dz, cells[i].M.data(), cells[p].M.data(), exporder);
  }
}


void P2P_Cells(size_t A, size_t B, std::vector<Cell> &cells,
  std::vector<Particle> &particles, double *F) {
    // A - target
    // B - source

    // P2P for the pair of cells
  //  std::cout << "    P2P_Cells("<<A<<","<<B<<")"<<std::endl;
  for (size_t p1 = 0; p1 < cells[A].nleaf; p1++) {
    size_t l1 = cells[A].leaf[p1];
    for (size_t p2 = 0; p2 < cells[B].nleaf; p2++) {
      size_t l2 = cells[B].leaf[p2];
      if (l2 != l1) {
	double dx = particles[l1].x - particles[l2].x;
	double dy = particles[l1].y - particles[l2].y;
	double dz = particles[l1].z - particles[l2].z;
	//std::cout << "      P2P("<<l1<<","<< l2 << ")" << std::endl;
	P2P(dx, dy, dz, particles[p2].q, &F[4 * l1]);
      }
    }
  }
}


void interact(std::vector<Cell> &cells, std::vector<Particle> &particles, double theta, size_t order, size_t ncrit, double *F) {
    interact_dehnen(0, 0, cells, particles, theta, order, ncrit, F);
}


void interact_dehnen(size_t A, size_t B, std::vector<Cell> &cells, std::vector<Particle> &particles, double theta, size_t order, size_t ncrit, double *F) {
  //  std::cout << "interact_dehnen("<<A<<","<<B<<")"<<std::endl;
  double dx = cells[A].x - cells[B].x;
  double dy = cells[A].y - cells[B].y;
  double dz = cells[A].z - cells[B].z;
  double R = sqrt(dx*dx + dy*dy + dz*dz);
  // If multipole acceptance criteria for FMM is met:
  //  std::cout << "  R*theta = " << R << std::endl;
  //std::cout << "  cells["<<A<<"].nchild = " << cells[A].nchild << std::endl;
  //std::cout << "  cells["<<B<<"].nchild = " << cells[B].nchild << std::endl;

  //  std::cout << "  cells["<<A<<"].rmax = " << cells[A].rmax << std::endl;
  // std::cout << "  cells["<<B<<"].rmax = " << cells[B].rmax << std::endl;
  // std::cout << "  (cells[" << A << "].rmax >= cells[" << B << "].rmax) = " << (cells[A].rmax >= cells[B].rmax) << std::endl;

  if (R*theta > (cells[A].rmax + cells[B].rmax)) {
    // then evaluate multipole to local from B's multipole to A's local
    // std::cout << "M2L(" << A << "," << B << ")" << std::endl;
    M2L(dx, dy, dz, cells[B].M.data(), cells[A].L.data(), order);
  }

  else if (cells[A].nchild == 0 && cells[B].nchild == 0) {
    if (cells[B].nleaf >= ncrit) {
      M2L(dx, dy, dz, cells[B].M.data(), cells[A].L.data(), order);
    }
    else {
      P2P_Cells(A, B, cells,particles, F);
    }
  }

  else if (cells[B].nchild == 0 || (cells[A].rmax >= cells[B].rmax && cells[A].nchild != 0)) {
    for(int oa = 0; oa < 8; oa++) {
      // For all 8 children of A, if child exists
      {
	if (cells[A].nchild & (1 << oa)) {
	  size_t a = cells[A].child[oa];
	  interact_dehnen(a, B, cells, particles, theta, order, ncrit, F);
	}
      }
    }
  }

  else {
    for(int ob = 0; ob < 8; ob++) {
      // for all 8 children of B, if child exists:
      if (cells[B].nchild & (1 << ob)) {
        int b = cells[B].child[ob];
        interact_dehnen(A, b, cells, particles, theta, order, ncrit, F);
      }
    }
  }
  //  std::cout << "Done interact_dehnen(" << A << "," << B << ")" << std::endl;
}



void evaluate_L2L(std::vector<Cell> &cells, size_t exporder) {
  // Note, cells can't have more than one parent so this is
  // inherently thread safe.
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < cells.size(); i++) {
    for (int octant = 0; octant < 8; octant++) {
      if (cells[i].nchild & (1 << octant)) {
        // for child in cell i
        size_t c = cells[i].child[octant];
        double dx = cells[c].x - cells[i].x;
        double dy = cells[c].y - cells[i].y;
        double dz = cells[c].z - cells[i].z;
        L2L(dx, dy, dz, cells[i].L.data(), cells[c].L.data(), exporder);
      }
    }
  }
}

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, size_t ncrit, size_t exporder) {
  // Each particle resides only in a single cell, so this is inherently
  // parallel.
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < cells.size(); i++) {
    if (cells[i].nleaf < ncrit) {
      // std::cout << "cell " << i << " is a leaf " << std::endl;
      for (size_t p = 0; p < cells[i].nleaf; p++) {
	size_t k = cells[i].leaf[p];
	// std::cout << "L2P from " << i << " to particle " << k << std::endl;
        double dx = particles[k].x - cells[i].x;
        double dy = particles[k].y - cells[i].y;
        double dz = particles[k].z - cells[i].z;
        L2P(dx, dy, dz, cells[i].L.data(), &F[4 * k], exporder);
      }
    }
  }
}

// void evaluate_approx(std::vector<Particle> &particles, std::vector<Cell>
// &cells,
//                      std::vector<double> &F, size_t n_crit, double
//                      theta, size_t exp_order) {
//   for (size_t i = 0; i < particles.size(); i++) {
//     evaluate_M2P_and_P2P(particles, 0, i, cells, F, n_crit, theta,
//     exp_order);
//   }
// }

void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &F) {
  #pragma omp parallel for
  for (size_t i = 0; i < particles.size(); i++) {
    for (size_t j = 0; j < particles.size(); j++) {
      if (i != j) {
        double dx = particles[i].x - particles[j].x;
        double dy = particles[i].y - particles[j].y;
        double dz = particles[i].z - particles[j].z;
        // calculation of R and R3 will be inlined by compiler
        // so no need to worry about that.
        P2P(dx, dy, dz, particles[j].q, &F[4 * i]);
      }
    }
  }
}
