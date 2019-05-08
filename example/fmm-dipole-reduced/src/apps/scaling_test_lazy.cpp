#include "calculate.hpp"
#include "tree.hpp"
#include "utils.hpp"
#include "operators.h"
#include "omp.h"
#include <chrono>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

int main(int argc, char **argv) {
  // Set initial parameters by user input from
  // the command line:
  size_t Nparticles = std::stoul(argv[1]);
  size_t ncrit = std::stoul(argv[2]);
  double theta = std::stod(argv[3]);
  size_t maxorder = std::stoul(argv[4]);

  std::cout << "Scaling Test Parameters" << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << "Nparticles = " << Nparticles << std::endl;
  std::cout << "ncrit      = " << ncrit << std::endl;
  std::cout << "theta      = " << theta << std::endl;
  std::cout << "maxorder   = " << maxorder << std::endl;

  std::vector<double> F_exact(3 * Nparticles, 0.0);
  std::vector<Particle> particles;
  std::default_random_engine generator(0.0);
  std::uniform_real_distribution<double> distribution(-1, 1);

  double mux_total = 0.0;
  double muy_total = 0.0;
  double muz_total = 0.0;

  double *mu = new double[3*Nparticles];

  for (size_t i = 0; i < Nparticles; i++) {
    double mux = distribution(generator);
    double muy = distribution(generator);
    double muz = distribution(generator);

    double mod = std::sqrt(mux*mux + muy*muy + muz*muz);
    mux /= mod;
    muy /= mod;
    muz /= mod;

    mux_total += mux;
    muy_total += muy;
    muz_total += muz;
    mu[3*i+0] = mux;
    mu[3*i+1] = muy;
    mu[3*i+2] = muz;

    Particle tmp(distribution(generator), distribution(generator),
           distribution(generator), &mu[3*i]);

    particles.push_back(tmp);
  }



  //std::cout << "\n\n\n" << std::endl;


  std::cout << "Direct\n------" << std::endl;
  Timer timer1;
  evaluate_direct(particles, F_exact);
  double t1 = timer1.elapsed();
  std::cout << "Time = " << t1 << std::endl;

  for (size_t order = 2; order < maxorder; order++) {
    std::cout << "Order " << order << "\n-------" << std::endl;
    std::vector<double> F_approx(3 * Nparticles, 0.0);

    auto root = Cell(0.0, 0.0, 0.0, 1.0, 0, order, 0, ncrit);
    auto cells = build_tree(particles, root, ncrit, order);

    std::cout << "Tree built with " << cells.size() << " cells.\n\n\n" << std::endl;
    std::vector<std::pair<size_t, size_t>> M2L_Cell_list;
    std::vector<std::pair<size_t, size_t>> P2P_Cell_list;

    double *M = new double[cells.size() * (Nterms(order) - Nterms(0))];
    double *L = new double[cells.size() * Nterms(order - 1)];

    for(int i = 0; i < cells.size(); i++) {
      cells[i].M = &M[i*(Nterms(order) - Nterms(0))];
      cells[i].L = &L[i*(Nterms(order - 1))];
    }


    interact_dehnen_lazy(0, 0, cells, particles, theta, order, ncrit, M2L_Cell_list, P2P_Cell_list);

    std::sort(M2L_Cell_list.begin(), M2L_Cell_list.end(),
             [](std::pair<size_t, size_t> &left, std::pair<size_t, size_t> &right) {
  	            return left.first < right.first;
  	           }
             );

    std::vector<omp_lock_t> M2L_locks(M2L_Cell_list.size());
    for(size_t i = 0; i < M2L_locks.size(); i++) {
      omp_init_lock(&M2L_locks[i]);
    }

    std::vector<omp_lock_t> P2P_locks(P2P_Cell_list.size());
    for(size_t i = 0; i < P2P_locks.size(); i++) {
      omp_init_lock(&P2P_locks[i]);
    }

    Timer timer2;
    evaluate_P2M(particles, cells, 0, ncrit, order);
    evaluate_M2M(particles, cells, order);

    #pragma omp parallel
    {
      #pragma omp single nowait
      {
        for(size_t i = 0; i < M2L_Cell_list.size(); i++) {
          #pragma omp task firstprivate(i)
          {
            size_t B = M2L_Cell_list[i].first;
            size_t A = M2L_Cell_list[i].second;
            omp_set_lock(&M2L_locks[A]);
            double dx =	cells[A].x - cells[B].x;
            double dy =	cells[A].y - cells[B].y;
            double dz =	cells[A].z - cells[B].z;
            M2L(dx, dy, dz, cells[B].M, cells[A].L, order);
            omp_unset_lock(&M2L_locks[A]);
          }
        }

        for(size_t i = 0; i < P2P_Cell_list.size(); i++) {
          #pragma omp task firstprivate(i)
          {
            size_t A = P2P_Cell_list[i].first;
            size_t B = P2P_Cell_list[i].second;
            omp_set_lock(&P2P_locks[A]);
            P2P_Cells(A, B, cells, particles, F_approx.data());
            omp_unset_lock(&P2P_locks[A]);
          }
        }
      }
    }

    evaluate_L2L(cells, order);
    evaluate_L2P(particles, cells, F_approx.data(), ncrit, order);

    double t2 = timer2.elapsed();

    double vrel_err = 0;
    double Exrel_err = 0;
    double Eyrel_err = 0;
    double Ezrel_err = 0;

    std::string filename = "error_order_" + std::to_string(order) + ".txt";
    std::ofstream fout(filename);
    for (size_t i = 0; i < particles.size(); i++) {
      double exerr =
        (F_exact[3 * i + 0] - F_approx[3 * i + 0]) / F_exact[3 * i + 0];
      double eyerr =
        (F_exact[3 * i + 1] - F_approx[3 * i + 1]) / F_exact[3 * i + 1];
      double ezerr =
        (F_exact[3 * i + 2] - F_approx[3 * i + 2]) / F_exact[3 * i + 2];
      fout << exerr << "," << eyerr << "," << ezerr << std::endl;
      Exrel_err += sqrt(exerr * exerr);
      Eyrel_err += sqrt(eyerr * eyerr);
      Ezrel_err += sqrt(ezerr * ezerr);
    }

    Exrel_err /= particles.size();
    Eyrel_err /= particles.size();
    Ezrel_err /= particles.size();

    std::cerr << "Rel errs = " << std::setw(10) << Exrel_err;
    std::cerr << ", " << std::setw(10) << Eyrel_err;
    std::cerr << ", " << std::setw(10) << Ezrel_err << std::endl;

    std::cout << "Approx. calculation  = " << t2 << " seconds. "
  	    << std::setw(10) << t2 / t1 * 100 << "% of direct time."
  	    << std::endl;


    for(int i = 0; i < M2L_locks.size(); i++) {
      omp_destroy_lock(&M2L_locks[i]);
    }

    for(size_t i = 0; i < P2P_locks.size(); i++) {
      omp_destroy_lock(&P2P_locks[i]);
    }

    delete[] M;
    delete[] L;

  }





  return 0;
}
