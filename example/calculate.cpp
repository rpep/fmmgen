#include "calculate.hpp"
#include "operators.h"
#include "tree.hpp"
#include "utils.hpp"
#include <iostream>
#include <stack>
#include <cmath>

void M_sanity_check(const std::vector<Cell> &cells) {
	double M0 = 0;
	for(size_t c = 1; c < cells.size(); c++) {
      if (cells[c].nchild == 0) {
		    M0 += cells[c].M[0];
	    }
  }
	std::cout << "Cell 0 has M0 = " << cells[0].M[0] << std::endl;
	std::cout << "Other cells   = " << M0 << std::endl;
  if (std::abs((cells[0].M[0] - M0)/M0) > 10e-10) {
    throw std::runtime_error("M0 sanity check failed");
  }
}


// interact_dehnen (the non-lazy, immediate-evaluation traversal) and P2P_Cells
// were removed here. interact_dehnen was only ever called from itself, so it
// was unreachable; P2P_Cells lost its last caller when evaluate_P2P_lazy was
// restructured to group by target below.

void interact_dehnen_lazy(const size_t A, const size_t B,
                          const std::vector<Cell> &cells,
                          const std::vector<Particle> &particles,
			                    const double theta, const size_t order,
                          const size_t ncrit,
                          std::vector<Interaction> &M2L_list,
                          std::vector<Interaction> &P2P_list) {
  const double dx = cells[A].x - cells[B].x;
  const double dy = cells[A].y - cells[B].y;
  const double dz = cells[A].z - cells[B].z;
  const double R = sqrt(dx*dx + dy*dy + dz*dz);

  if (R*theta > (cells[A].rmax + cells[B].rmax)) {
    //if (cells[A].nleaf < ncrit && cells[B].nleaf < ncrit) {
      M2L_list.push_back(Interaction{(uint32_t)A, (uint32_t)B});
    //}
  }

  else if (cells[A].nchild == 0 && cells[B].nchild == 0) {
    if (cells[B].nleaf >= ncrit) {
      M2L_list.push_back(Interaction{(uint32_t)A, (uint32_t)B});
      M2L(dx, dy, dz, cells[B].M, cells[A].L, order);
    }
    else {
      //if (cells[A].nleaf < ncrit and cells[B].nleaf < ncrit) {
    	P2P_list.push_back(Interaction{(uint32_t)A, (uint32_t)B});
      //}
    }
  }

  else if (cells[B].nchild == 0 || (cells[A].rmax >= cells[B].rmax && cells[A].nchild != 0)) {
    for(int oa = 0; oa < 8; oa++) {
      // For all 8 children of A, if child exists
      if (cells[A].nchild & (1 << oa)) {
    	int a = cells[A].child[oa];
    	interact_dehnen_lazy(a, B, cells, particles, theta, order, ncrit, M2L_list, P2P_list);
      }
    }
  }

  else {
    for(int ob = 0; ob < 8; ob++) {
      // for all 8 children of B, if child exists:
      if (cells[B].nchild & (1 << ob)) {
        int b = cells[B].child[ob];
        interact_dehnen_lazy(A, b, cells, particles, theta, order, ncrit, M2L_list, P2P_list);
      }
    }
  }
}

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
              size_t cell, size_t ncrit, size_t exporder) {
  // Allocate thread-local temporary array once per thread
  size_t msize = Msize(exporder, FMMGEN_SOURCEORDER);
  thread_local std::vector<double> M_thread;

  // Ensure thread-local buffer is sized correctly
  if (M_thread.size() != msize) {
    M_thread.resize(msize);
  }

  double *M = M_thread.data();

  #pragma omp for
  for(size_t c = 0; c < cells.size(); c++) {
    if (cells[c].nleaf < ncrit) {
      for(size_t i = 0; i < cells[c].nleaf; i++) {
        size_t l = cells[c].leaf[i];
        // Walter dehnen's definition:
        // (-1)^m / m! (x_a - z_a)^m
        double dx = (cells[c].x - particles[l].r[0]);
        double dy = (cells[c].y - particles[l].r[1]);
        double dz = (cells[c].z - particles[l].r[2]);
        for(int k = 0; k < FMMGEN_SOURCESIZE; k++) {
          M[k] = particles[l].S[k];
        }
        M2M(dx, dy, dz, M, cells[c].M, exporder);
      }
   }
  }
}

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  const std::vector<std::vector<size_t>> &levels, size_t exporder) {
  /*
  evaluate_M2M with level-synchronous parallel traversal.

  Process tree bottom-up, parallelizing within each level.
  Cells at the same level can be processed in parallel since
  they write to different parent cells.
  */
  // Bottom-up: start from deepest level and work up to level 1 (skip root at level 0).
  //
  // Parallelise over the PARENT cells at level l-1, not over the children at
  // level l. Siblings share a parent, so distributing children across threads
  // means up to 8 threads accumulate into the same cells[p].M concurrently.
  // Iterating parents gives each thread exclusive ownership of the cell it
  // writes, which is the same pattern evaluate_L2L already uses.
  for (int l = levels.size() - 1; l > 0; l--) {
    #pragma omp for schedule(static)
    for (size_t i = 0; i < levels[l-1].size(); i++) {
      size_t p = levels[l-1][i];
      for (int octant = 0; octant < 8; octant++) {
        if (cells[p].nchild & (1 << octant)) {
          size_t c = cells[p].child[octant];
          double dx = (cells[p].x - cells[c].x);
          double dy = (cells[p].y - cells[c].y);
          double dz = (cells[p].z - cells[c].z);
          M2M(dx, dy, dz, cells[c].M, cells[p].M, exporder);
        }
      }
    }
  }
}


void evaluate_M2L_lazy(std::vector<Cell> &cells,
                       std::vector<Interaction> &M2L_list,
                       std::vector<size_t> &M2L_group, size_t order) {
    // M2L_list is grouped by target, so one thread can own an entire target.
    // That thread is the only writer of cells[A].L, which lets M2L accumulate
    // straight into it. Compared with the previous version this removes, per
    // interaction: a heap allocation, the zeroing of the temporary, and an
    // lsize-long loop of atomic adds (364 atomics per interaction at p=11,
    // roughly 660M over a typical run).
    //
    // Entries for target A live in [M2L_group[A], M2L_group[A+1]); empty
    // targets are skipped. Groups vary widely in size, hence dynamic
    // scheduling.
    const size_t ngroups = M2L_group.empty() ? 0 : M2L_group.size() - 1;
    #pragma omp for schedule(dynamic, 16)
    for (size_t A = 0; A < ngroups; A++) {
        const size_t begin = M2L_group[A], end = M2L_group[A+1];
        if (begin == end) continue;
        double *const L = cells[A].L;
        const double ax = cells[A].x, ay = cells[A].y, az = cells[A].z;
        for (size_t i = begin; i < end; i++) {
            const size_t B = M2L_list[i].source;
            M2L(ax - cells[B].x, ay - cells[B].y, az - cells[B].z,
                cells[B].M, L, order);
        }
    }

}

void evaluate_P2P_lazy(std::vector<Cell> &cells, std::vector<Particle> &particles,
                       std::vector<Interaction> &P2P_list,
                       std::vector<size_t> &P2P_group, double *F) {
   // Grouped by target, as evaluate_M2L_lazy is. Every particle belongs to
   // exactly one leaf, so the thread owning target cell A is the sole writer of
   // F for A's particles and no atomics are needed. Removing them also makes
   // the whole solver bit-reproducible: accumulation order is now fixed by the
   // (stable) counting sort rather than by which thread got there first.
   //
   // Loop nesting is deliberately source-cell outer, target-particle inner --
   // the same order as the original code. Hoisting the target particle
   // outermost so its accumulator can live in registers looks attractive, and
   // is ~1.15x faster at N=20k, but it makes every target particle re-stream
   // all of its source cells. Once that working set stops fitting in cache the
   // trade inverts: at N=100k it measured 0.79x, i.e. materially slower.
   const size_t ngroups = P2P_group.empty() ? 0 : P2P_group.size() - 1;
   #pragma omp for schedule(dynamic, 16)
   for (size_t A = 0; A < ngroups; A++) {
       const size_t begin = P2P_group[A], end = P2P_group[A+1];
       if (begin == end) continue;
       for (size_t i = begin; i < end; i++) {
           const size_t B = P2P_list[i].source;
           for (size_t p1 = 0; p1 < cells[A].nleaf; p1++) {
               const size_t l1 = cells[A].leaf[p1];
               double *const Fl = &F[FMMGEN_OUTPUTSIZE*l1];
               for (size_t p2 = 0; p2 < cells[B].nleaf; p2++) {
                   const size_t l2 = cells[B].leaf[p2];
                   if (l2 != l1) {
                       P2P(particles[l1].r[0] - particles[l2].r[0],
                           particles[l1].r[1] - particles[l2].r[1],
                           particles[l1].r[2] - particles[l2].r[2],
                           particles[l2].S, Fl);
                   }
               }
           }
       }
   }
}



void evaluate_L2L(std::vector<Cell> &cells, const std::vector<std::vector<size_t>> &levels,
                  size_t exporder) {
  /*
  evaluate_L2L with level-synchronous parallel traversal.

  Process tree top-down, parallelizing within each level.
  Cells at the same level can be processed in parallel.
  */
  // Top-down: start from root level and work down
  for (size_t l = 0; l < levels.size() - 1; l++) {
    #pragma omp for schedule(static)
    for (size_t i = 0; i < levels[l].size(); i++) {
      size_t p = levels[l][i];
      for (int octant = 0; octant < 8; octant++) {
        if (cells[p].nchild & (1 << octant)) {
          size_t c = cells[p].child[octant];
          double dx = cells[c].x - cells[p].x;
          double dy = cells[c].y - cells[p].y;
          double dz = cells[c].z - cells[p].z;
          // Each thread processes different parent cells, so no race condition
          L2L(dx, dy, dz, cells[p].L, cells[c].L, exporder);
        }
      }
    }
  }
}

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, size_t ncrit, size_t exporder) {
  #pragma omp for schedule(runtime)
  for (size_t i = 0; i < cells.size(); i++) {
    if (cells[i].nleaf < ncrit) {
      for (size_t p = 0; p < cells[i].nleaf; p++) {
	    size_t k = cells[i].leaf[p];
        double dx = particles[k].r[0] - cells[i].x;
        double dy = particles[k].r[1] - cells[i].y;
        double dz = particles[k].r[2] - cells[i].z;
        L2P(dx, dy, dz, cells[i].L, &F[FMMGEN_OUTPUTSIZE*k], exporder);
      }
    }
  }
}

void evaluate_direct(std::vector<Particle> &particles, double *F, size_t n) {
  #pragma omp parallel for schedule(runtime)
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (i != j) {
      	double dx = particles[i].r[0] - particles[j].r[0];
      	double dy = particles[i].r[1] - particles[j].r[1];
      	double dz = particles[i].r[2] - particles[j].r[2];
      	P2P(dx, dy, dz, particles[j].S, &F[FMMGEN_OUTPUTSIZE*i]);
      }
    }
  }
}

void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p, unsigned int i,
  std::vector<Cell> &cells, double *F, unsigned int n_crit, double theta,
  unsigned int exporder) {
  // For particle i, start at cell p
  double dx, dy, dz, r;
  int c;
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
        dx = particles[i].r[0] - cells[c].x;
        dy = particles[i].r[1] - cells[c].y;
        dz = particles[i].r[2] - cells[c].z;
        r = sqrt(dx*dx + dy*dy + dz*dz);
        // Apply the Barnes-Hut criterion:
        if (cells[c].r > theta * r) {
            // If the cell is 'near':
            evaluate_M2P_and_P2P(particles, c, i, cells, F, n_crit, theta, exporder);
        }
        else {
            // If the cell is 'far', calculate the potential and field
            // on the particle from the multipole expansion:
            M2P(dx, dy, dz, cells[c].M, &F[FMMGEN_OUTPUTSIZE*i], exporder);
        }
      }
    }
  }
  else {
    // loop in leaf cell's particles
    for(unsigned int l = 0; l < (cells[p].nleaf); l++) {
      // Get the particle index:
      j = cells[p].leaf[l];
      if (i != j) {
        // Calculate the interparticle distance:
        dx = particles[i].r[0] - particles[j].r[0];
        dy = particles[i].r[1] - particles[j].r[1];
        dz = particles[i].r[2] - particles[j].r[2];
        // Compute the field:
        P2P(dx, dy, dz, particles[j].S, &F[FMMGEN_OUTPUTSIZE*i]);
      }
    }
  }
}
