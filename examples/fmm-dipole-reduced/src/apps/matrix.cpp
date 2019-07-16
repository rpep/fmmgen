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

  // Array containing r and mu for faster memory access
  double *r = new double[3*Nparticles];
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


    r[3*i+0] = distribution(generator);
    r[3*i+1] = distribution(generator);
    r[3*i+2] = distribution(generator);
    mu[3*i+0] = mux;
    mu[3*i+1] = muy;
    mu[3*i+2] = muz;
    Particle tmp(&r[3*i],
		 &mu[3*i]);

    particles.push_back(tmp);
  }



  //std::cout << "\n\n\n" << std::endl;


  std::cout << "Direct\n------" << std::endl;
  Timer timer1;
  evaluate_direct(particles, F_exact, Nparticles);
  double t1 = timer1.elapsed();
  std::cout << "Time = " << t1 << std::endl;

  size_t order = 2;
  std::cout << "Order " << order << "\n-------" << std::endl;
  std::vector<double> F_approx(3 * Nparticles, 0.0);
  
  auto root = Cell(0.0, 0.0, 0.0, 1.0, 0, order, 0, ncrit);
  auto cells = build_tree(particles, root, ncrit, order);
  
  std::cout << "Tree built with " << cells.size() << " cells.\n\n\n" << std::endl;
  std::vector<std::pair<size_t, size_t>> M2L_Cell_list;
  std::vector<std::pair<size_t, size_t>> P2P_Cell_list;
  std::vector<double> M(cells.size() * (Nterms(order) - 1), 0.0);
  std::vector<double> L(cells.size() * Nterms(order - 1), 0.0);
  
  for(size_t i = 0; i < cells.size(); i++) {
    cells[i].M = &M[i*(Nterms(order) - Nterms(0))];
    cells[i].L = &L[i*(Nterms(order - 1))];
  }

  interact_dehnen_lazy(0, 0, cells, particles, theta, order, ncrit, M2L_Cell_list, P2P_Cell_list);
  
  std::sort(M2L_Cell_list.begin(), M2L_Cell_list.end(),
	    [](std::pair<size_t, size_t> &left, std::pair<size_t, size_t> &right) {
	      return left.first < right.first;
  	           }
	    );

  std::sort(P2P_Cell_list.begin(), P2P_Cell_list.end());
  
  std::string m2lfile = "m2lfile_" + std::to_string(Nparticles) + ".txt";
  std::ofstream m2lout(m2lfile);

  std::string p2pfile = "p2pfile_" + std::to_string(Nparticles) + ".txt";
  std::ofstream p2pout(p2pfile);
  
  for (size_t i = 0; i < M2L_Cell_list.size(); i++) {
    m2lout << M2L_Cell_list[i].first << "," << M2L_Cell_list[i].second << std::endl;
  }

  for (size_t i = 0; i < P2P_Cell_list.size(); i++) {
    p2pout << P2P_Cell_list[i].first << "," << P2P_Cell_list[i].second << std::endl;
  }
  
  
  delete[] mu;
  return 0;
}
