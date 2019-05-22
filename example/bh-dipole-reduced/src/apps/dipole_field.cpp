#include<iostream>
#include<iomanip>
#include<fstream>
#include<chrono>
#include<random>
#include "calculate.hpp"
#include "utils.hpp"
#include "tree.hpp"
#include "operators.h"

int main(int argc, char** argv) {
  unsigned int exporder = 5;

  for(int order = 0; order < exporder; order++) {
    unsigned int ncrit = 100;
    unsigned int Nparticles = 1;
    std::vector<double> F_exact(4*Nparticles, 0.0);
    std::vector<double> F_approx(4*Nparticles, 0.0);
    
    std::vector<Particle> particles;
    std::default_random_engine generator(0.0);
    std::uniform_real_distribution<double> distribution(-1, 1);
    
    double scale = 1e-8;
    double mux = 1.0;
    double muy = 0.0;
    double muz = 0.0;
    
    double *M = new double[Nterms(order)];
    
    auto root = Cell(0.0, 0.0,  0.0, 1.0, 0, order, 0, ncrit);
    
    double xa = -1.0*scale;
    double ya = 0.1*scale;
    double za = -0.3*scale;
    
    Particle tmp1(xa, ya, za, mux, muy, muz);
    
    particles.push_back(tmp1);
    auto cells = build_tree(particles, root, ncrit, order);
    std::cout << "Tree built with " << cells.size() << " cells." << std::endl;
    evaluate_P2M(particles, cells, 0, ncrit, order);

  
    double xb = 20.1*scale;
    double yb = -0.3*scale;
    double zb = 1.0*scale;
    
    P2P(xb-xa, yb-ya, zb-za, mux, muy, muz, F_exact.data());

    double Fvals[4] = {0.0, 0.0, 0.0, 0.0};

    M2P(xb, yb, zb, cells[0].M.data(), Fvals, order);

    F_approx[0] -= Fvals[0];
    F_approx[1] -= Fvals[1];
    F_approx[2] -= Fvals[2];
    F_approx[3] -= Fvals[3];

    std::cout << "Order = " << order << std::endl;
    std::cout << "\tExact = " << F_exact[0] << ", " << F_exact[1] << ", " << F_exact[2] << ", " << F_exact[3] << std::endl;    
    std::cout << "\tApprox = " << F_approx[0] << ", " << F_approx[1] << ", " << F_approx[2] << ", " << F_approx[3] << std::endl;
   
   

    delete[] M;
  }
  return 0;
}
