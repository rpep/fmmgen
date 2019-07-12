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

  std::cout << "Scaling Test Parameters" << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << "Nparticles = " << Nparticles << std::endl;

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



  std::cout << "Direct\n------" << std::endl;
  auto t1 = omp_get_wtime();
  evaluate_direct(particles, F_exact, Nparticles);
  auto t2 = omp_get_wtime();
  std::cout << "Time = " << t2 - t1 << std::endl;

  delete[] mu;
  delete[] r;
  return 0;
}
