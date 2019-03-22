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
                  unsigned int cell, unsigned int ncrit,
                  unsigned int exporder) {
  if (cells[cell].nleaf > ncrit) {
    for (unsigned int octant = 0; octant < 8; octant++) {
      if (cells[cell].nchild & (1 << octant)) {
        evaluate_P2M(particles, cells, cells[cell].child[octant], ncrit,
                     exporder);
      }
    }
  } else {
    for (unsigned int i = 0; i < (cells[cell].nleaf); i++) {
      unsigned int l = cells[cell].leaf[i];
      double dx = particles[l].x - cells[cell].x;
      double dy = particles[l].y - cells[cell].y;
      double dz = particles[l].z - cells[cell].z;
      P2M(dx, dy, dz, particles[l].q, cells[cell].M.data(), exporder);
    }
  }
}

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  unsigned int exporder) {
  /*
  evaluate_M2M(particles, cells)

  This function evaluates the multipole to
  multipole kernel. It does this by working up the
  tree from the leaf nodes, which is possible
  by iterating backwards through the nodes because
  of the way the tree is constructed.
  */

  for (unsigned int i = cells.size() - 1; i > 0; i--) {
    unsigned int p = cells[i].parent;
    double dx = cells[p].x - cells[i].x;
    double dy = cells[p].y - cells[i].y;
    double dz = cells[p].z - cells[i].z;
    M2M(dx, dy, dz, cells[i].M.data(), cells[p].M.data(), exporder);
  }
}

void interact(unsigned int a, unsigned int b, double *F,
              std::vector<Cell> &cells, std::vector<Particle> &particles,
              std::stack<std::pair<unsigned int, unsigned int>> &stack,
              unsigned int ncrit, double theta, unsigned int order) {

  // interact(A, B)
  //   if A and B are both leafs then
  //     call P2P kernel
  //   else
  //     if A and B satisfy MAC then
  //       call M2L kernel
  //     else
  //       push pair (A, B) to stack
  //     end if
  //   end if

  // If both cells are leaf cells
  // std::cout << "Interact(" << a << ", " << b << ")" << std::endl;
  if ((cells[a].nleaf < ncrit) && (cells[b].nleaf < ncrit)) {
    // std::cout << "Both leaves" << std::endl;
    // P2P for the pair of cells
    for (unsigned int p1 = 0; p1 < cells[a].nleaf; p1++) {
      unsigned int l1 = cells[a].leaf[p1];
      for (unsigned int p2 = 0; p2 < cells[b].nleaf; p2++) {
        unsigned int l2 = cells[b].leaf[p2];
        if (l2 != l1) {
          double dx = particles[p1].x - particles[p2].x;
          double dy = particles[p1].y - particles[p2].y;
          double dz = particles[p2].z - particles[p2].z;
          // std::cout << "P2P disabled" << std::endl;
          P2P(dx, dy, dz, particles[p2].q, &F[4 * p1]);
        }
      }
    }

  } else {
    /// std::cout << "Not leaves!" << std::endl;
    double dx = cells[a].x - cells[b].x;
    double dy = cells[a].y - cells[b].y;
    double dz = cells[a].z - cells[b].z;
    double r = sqrt(dx * dx + dy * dy + dz * dz);
    if ((cells[a].r + cells[b].r) > theta * r) {
      // std::cout << "M2L" << std::endl;
      // std::cout << "dx = " << dx;
      // std::cout << "dy = " << dy;
      // std::cout << "dz = " << dz;
      // std::cout << "dr = " << r << std::endl;

      M2L(dx, dy, dz, cells[b].M.data(), cells[a].L.data(), order);
    } else {
      // std::cout << "Adding to stack" << std::endl;
      std::pair<unsigned int, unsigned int> pair = std::make_pair(a, b);
      stack.push(pair);
    }
  }
}

void FMMDualTreeTraversal(std::vector<Particle> &particles,
                          std::vector<Cell> &cells, double *F,
                          unsigned int ncrit, double theta,
                          unsigned int exporder) {
  std::stack<std::pair<unsigned int, unsigned int>> stack;
  // push pair of root cells (A, B) to stack
  stack.push(std::make_pair<unsigned int, unsigned int>(0, 0));
  // while stack not empty
  while (stack.size() > 0) {
    std::pair<unsigned int, unsigned int> pair = stack.top();
    stack.pop();
    unsigned int A = pair.first;
    unsigned int B = pair.second;
    // std::cout << "FMM - A = " << A << " B = " << B << std::endl;
    if (cells[A].r > cells[B].r) {
      //     if target cell is larger then source cell
      //     for all children a of target cell A
      for (unsigned int octant = 0; octant < 8; octant++) {
        if (cells[A].nchild & (1 << octant)) {
          unsigned int a = cells[A].child[octant];
          interact(a, B, F, cells, particles, stack, ncrit, theta, exporder);
        }
      }
    } else {
      for (unsigned int i = 0; i < cells[B].nleaf; i++) {
        // for all children b of source cell B
        for (unsigned int octant = 0; octant < 8; octant++) {
          // If a child exists in a given octant:
          if (cells[B].nchild & (1 << octant)) {
            unsigned int b = cells[B].child[octant];
            interact(A, b, F, cells, particles, stack, ncrit, theta, exporder);
          }
        }
      }
    }
  }
}

void evaluate_L2L(std::vector<Cell> &cells, unsigned int exporder) {
  for (unsigned int i = 0; i < cells.size(); i++) {
    for (unsigned int octant = 0; octant < 8; octant++) {
      if (cells[i].nchild & (1 << octant)) {
        // for child in cell i
        unsigned int c = cells[i].child[octant];
        double dx = cells[i].x - cells[c].x;
        double dy = cells[i].y - cells[c].y;
        double dz = cells[i].z - cells[c].z;
        L2L(dx, dy, dz, cells[i].L.data(), cells[c].L.data(), exporder);
      }
    }
  }
}

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, unsigned int ncrit, unsigned int exporder) {
  for (unsigned int i = 0; i < cells.size(); i++) {
    if (cells[i].nleaf < ncrit) {
      for (unsigned int p = 0; p < cells[i].nleaf; p++) {
        double dx = cells[i].x - particles[p].x;
        double dy = cells[i].y - particles[p].y;
        double dz = cells[i].z - particles[p].z;
        L2P(dx, dy, dz, cells[i].L.data(), &F[4 * p], exporder);
      }
    }
  }
}

void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p,
                          unsigned int i, std::vector<Cell> &cells,
                          std::vector<double> &F, unsigned int n_crit,
                          double theta, unsigned int exporder) {
  // For cell i, start at cell p (initially 0).
  double dx, dy, dz, r;
  unsigned int c;
  unsigned int j;
  // If cell p is not a leaf cell:
  if (cells[p].nleaf >= n_crit) {
    // Iterate through it's children
    for (unsigned int octant = 0; octant < 8; octant++) {
      // If a child exists in a given octant:
      if (cells[p].nchild & (1 << octant)) {
        // Get the child's index
        c = cells[p].child[octant];
        // Calculate the distance from the particle to the child cell
        dx = particles[i].x - cells[c].x;
        dy = particles[i].y - cells[c].y;
        dz = particles[i].z - cells[c].z;
        r = sqrt(dx * dx + dy * dy + dz * dz);
        // Apply the Barnes-Hut criterion:
        if (cells[c].r > theta * r) {
          // If the cell is 'near':
          evaluate_M2P_and_P2P(particles, c, i, cells, F, n_crit, theta,
                               exporder);
        } else {
          // If the cell is 'far', calculate the potential and field
          // on the particle from the multipole expansion:
          M2P(dx, dy, dz, cells[c].M.data(), &F[4 * i], exporder);
        }
      }
    }
  } else {
    // loop in leaf cell's particles
    for (unsigned int l = 0; l < (cells[p].nleaf); l++) {
      // Get the particle index:
      j = cells[p].leaf[l];
      if (i != j) {
        // Calculate the interparticle distance:
        dx = particles[i].x - particles[j].x;
        dy = particles[i].y - particles[j].y;
        dz = particles[i].z - particles[j].z;
        // Compute the field:
        P2P(dx, dy, dz, particles[j].q, &F[4 * i]);
      }
    }
  }
}

void evaluate_approx(std::vector<Particle> &particles, std::vector<Cell> &cells,
                     std::vector<double> &F, unsigned int n_crit, double theta,
                     unsigned int exp_order) {
  for (unsigned int i = 0; i < particles.size(); i++) {
    evaluate_M2P_and_P2P(particles, 0, i, cells, F, n_crit, theta, exp_order);
  }
}

void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &F) {
#pragma omp parallel for
  for (unsigned int i = 0; i < particles.size(); i++) {
    for (unsigned int j = 0; j < particles.size(); j++) {
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
