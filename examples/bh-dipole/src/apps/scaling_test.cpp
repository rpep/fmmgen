#include "calculate.hpp"
#include "tree.hpp"
#include "utils.hpp"
#include <string>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>


int main(int argc, char **argv) {
  // Set initial parameters by user input from
  // the command line:
  unsigned int Nparticles = std::stoul(argv[1]);
  unsigned int ncrit = std::stoul(argv[2]);
  double theta = std::stod(argv[3]);
  unsigned int maxorder = std::stoul(argv[4]);

  std::cout << "Scaling Test Parameters" << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << "Nparticles = " << Nparticles << std::endl;
  std::cout << "ncrit      = " << ncrit << std::endl;
  std::cout << "theta      = " << theta << std::endl;
  std::cout << "maxorder   = " << maxorder << std::endl;

  std::vector<double> F_exact(4 * Nparticles, 0.0);
  std::vector<Particle> particles;
  std::default_random_engine generator(0.0);
  std::uniform_real_distribution<double> distribution(-1, 1);
  
  double mux_total = 0.0;
  double muy_total = 0.0;
  double muz_total = 0.0;

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

    Particle tmp(distribution(generator), distribution(generator),
           distribution(generator), mux, muy, muz);

    particles.push_back(tmp);
  }
  
  std::cout << "Direct\n------" << std::endl;
  Timer timer1;
  evaluate_direct(particles, F_exact);
  double t1 = timer1.elapsed();
  std::cout << "Time = " << t1 << std::endl;

  for(unsigned int order = 0; order < maxorder; order++) {
    std::cout << "Order " << order << "\n-------" << std::endl;
    std::vector<double> F_approx(4 * Nparticles, 0.0); 
    Timer timer2;
    auto root = Cell(0.0, 0.0, 0.0, 1.0, 0, order, 0, ncrit);
    auto cells = build_tree(particles, root, ncrit, order);
    std::cout << "Tree built with " << cells.size() << " cells." << std::endl;
    evaluate_P2M(particles, cells, 0, ncrit, order);
    evaluate_M2M(particles, cells, order);
    
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < particles.size(); i++) {
      evaluate_M2P_and_P2P(particles, 0, i, cells, F_approx, ncrit, theta, order);
    }

    double t2 = timer2.elapsed();
    double vrel_err = 0;
    double Exrel_err = 0;
    double Eyrel_err = 0;
    double Ezrel_err = 0;


    std::string filename = "error_order_" + std::to_string(order) + ".txt";
    std::ofstream fout(filename);
    for (unsigned int i = 0; i < particles.size(); i++) {
      double verr  = (F_exact[4*i+0] - F_approx[4*i+0]) / F_exact[4*i+0];
      double exerr = (F_exact[4*i+1] - F_approx[4*i+1]) / F_exact[4*i+1];
      double eyerr = (F_exact[4*i+2] - F_approx[4*i+2]) / F_exact[4*i+2];
      double ezerr = (F_exact[4*i+3] - F_approx[4*i+3]) / F_exact[4*i+3];
      fout << verr << "," << exerr << "," << eyerr << "," << ezerr << std::endl;
      vrel_err += sqrt(verr * verr);
      Exrel_err += sqrt(exerr*exerr);
      Eyrel_err += sqrt(eyerr*eyerr);
      Ezrel_err += sqrt(ezerr*ezerr);
    }
    
    vrel_err /= particles.size();
    Exrel_err /= particles.size();
    Eyrel_err /= particles.size();
    Ezrel_err /= particles.size();
    
    std::cout << "Rel errs = " << std::setw(10) << vrel_err;
    std::cout << ", " << std::setw(10) << Exrel_err;
    std::cout << ", " << std::setw(10) << Eyrel_err;
    std::cout << ", " << std::setw(10) << Ezrel_err << std::endl;

    std::cout << "Approx. calculation  = " << t2 << " seconds. " << std::setw(10) << t2/t1 * 100 << "% of direct time." << std::endl;

  }

  return 0;
}
