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
  size_t type = std::stoul(argv[4]); // type == 0, FMM, type == 1, BH

  const size_t calc_direct = 1;	  
  std::cout << "Scaling Test Parameters" << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << "Nparticles = " << Nparticles << std::endl;
  std::cout << "ncrit      = " << ncrit << std::endl;
  std::cout << "theta      = " << theta << std::endl;
  std::cout << "FMMGEN_MINORDER = " << FMMGEN_MINORDER << std::endl;
  std::cout << "FMMGEN_MAXORDER = " << FMMGEN_MAXORDER << std::endl;
  std::cout << "FMMGEN_SOURCEORDER = " << FMMGEN_SOURCEORDER << std::endl;  
  std::cout << "FMMGEN_OUTPUTSIZE = " << FMMGEN_OUTPUTSIZE << std::endl;
  std::cout << "FMMGEN_SOURCESIZE = " << FMMGEN_SOURCESIZE << std::endl;
  std::cout << "FMMGEN TYPE = ";
  if (type == 0) {
    std::cout << "Fast Multipole Method (Lazy Evaluation)" << std::endl;
  }
  if (type == 1) {
    std::cout << "Barnes-Hut Method" << std::endl;
  }

  std::vector<double> F_exact(FMMGEN_OUTPUTSIZE * Nparticles, 0.0);
  std::vector<double> F_approx(FMMGEN_OUTPUTSIZE * Nparticles, 0.0);
  std::default_random_engine generator(0.0);
  std::uniform_real_distribution<double> distribution(-1, 1);

  double mux_total = 0.0;
  double muy_total = 0.0;
  double muz_total = 0.0;

  // Array containing r and source strengths
  double *r = new double[3*Nparticles];
  double *S = new double[FMMGEN_SOURCESIZE*Nparticles];
  auto filename = "particles_n_" + std::to_string(Nparticles) + ".txt";
  std::ofstream fout;
  fout.open(filename);
  for (size_t i = 0; i < Nparticles; i++) {
    for(int j = 0; j < 3; j++) {
      r[3*i+j] = distribution(generator) * 1e-9;
      fout << r[3*i+j] << ",";
    }
    for(int j = 0; j < FMMGEN_SOURCESIZE; j++) {
      S[FMMGEN_SOURCESIZE*i + j] = distribution(generator);
      fout << S[FMMGEN_SOURCESIZE*i + j] << ",";
    }
    fout << std::endl;
  }
  fout.close();


  double t_direct;
  double t_approx;

  for (size_t order = FMMGEN_MINORDER; order < FMMGEN_MAXORDER; order++) {
    Tree tree = build_tree(r, S, Nparticles, ncrit, order, theta);
    std::cout << "Tree built with " << tree.cells.size() << " cells.\n\n\n" << std::endl;
    std::cout << "Order " << order << "\n-------" << std::endl;
    std::fill(F_approx.begin(), F_approx.end(), 0);
    if (order == FMMGEN_MINORDER && calc_direct) {
      Timer timer;
      tree.compute_field_exact(F_exact.data());
      t_direct = timer.elapsed();
      std::cout << "t_direct = " << t_direct << std::endl;
    }

    Timer timer;
    if (type == 0) {
	   tree.compute_field_fmm(F_approx.data());
    }
    else if (type == 1) {
	   tree.compute_field_bh(F_approx.data());
    }
    t_approx = timer.elapsed();

    if (calc_direct) {
        double Exrel_err = 0;
        double Eyrel_err = 0;
        double Ezrel_err = 0;

        auto filename = "errors_lazy_p_" + std::to_string(order) +
                               "_n_" + std::to_string(Nparticles) +
                               "_ncrit_" + std::to_string(ncrit) +
                               "_theta_" + std::to_string(theta) + 
    			   "_type_" + std::to_string(type) + ".txt";
        std::ofstream fout(filename);

        double errs[FMMGEN_OUTPUTSIZE] = {0.0};
        for (size_t i = 0; i < Nparticles; i++) {
            for(int k = 0; k < FMMGEN_OUTPUTSIZE; k++) {
              double err = (F_exact[FMMGEN_OUTPUTSIZE * i + k] - F_approx[FMMGEN_OUTPUTSIZE * i + k]) / F_exact[FMMGEN_OUTPUTSIZE * i + k];
              fout << err << ",";
              errs[k] += sqrt(err * err);
          }
          fout << std::endl;
        }

        std::cerr << "Rel errs = " << std::scientific;
        for(int k = 0; k < FMMGEN_OUTPUTSIZE; k++) {
            std::cerr << std::setw(10) << errs[k] / Nparticles << ",";
        }
        std::cout << std::endl;
    }

    std::cout << "Approx. calculation  = " << t_approx << " seconds. " << std::endl;
    if (calc_direct) {
	    std:: cout << std::setw(10) << t_approx / t_direct * 100 << "% of direct time."
  	    << std::endl;
    }

    if (calc_direct) {
        auto filename = "field_p_" + std::to_string(order) +
                             "_n_" + std::to_string(Nparticles) +
                             "_ncrit_" + std::to_string(ncrit) +
                             "_theta_" + std::to_string(theta) + 
             "_type_" + std::to_string(type) + ".txt";
        std::ofstream fout(filename);

        for (size_t i = 0; i < Nparticles; i++) {
            for(int k = 0; k < FMMGEN_OUTPUTSIZE; k++) {
              fout << F_exact[FMMGEN_OUTPUTSIZE * i + k] << "," << F_approx[FMMGEN_OUTPUTSIZE * i + k] << ",";
            }
          fout << std::endl;
        }
    }



    






  	
  }









  delete[] r;
  delete[] S;
  return 0;
}
