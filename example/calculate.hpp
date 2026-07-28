//############################################
//#
//# Functions for running the Barnes-Hut
//# method for gravitational source particles.
//#
//# (C) Ryan Pepper, 2018
//# University of Southampton, UK
//#
//#
//###########################################
#pragma once
#include "tree.hpp"
#include "utils.hpp"
#include <iostream>
#include<omp.h>

void M_sanity_check(const std::vector<Cell> &cells);

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
		  size_t cell, size_t ncrit, size_t exporder);

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  const std::vector<std::vector<size_t>> &levels, size_t exporder);


void evaluate_L2L(std::vector<Cell> &cells, const std::vector<std::vector<size_t>> &levels,
                  size_t exporder);

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, size_t ncrit, size_t exporder);

void evaluate_direct(std::vector<Particle> &particles, double *F, size_t Nparticles);


void interact_dehnen_lazy(const size_t A, const size_t B, const std::vector<Cell> &cells, const std::vector<Particle> &particles,
			  const double theta, const size_t order, const size_t ncrit,
			  std::vector<Interaction> &M2L_list,
			  std::vector<Interaction> &P2P_list);




void evaluate_M2L_lazy(std::vector<Cell> &cells,
                    std::vector<Interaction> &M2L_list,
                    std::vector<size_t> &M2L_group, size_t order);

void evaluate_P2P_lazy(std::vector<Cell> &cells, std::vector<Particle> &particles,
                    std::vector<Interaction> &P2P_list,
                    std::vector<size_t> &P2P_group, double *F);

void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p, unsigned int i,
   std::vector<Cell> &cells, double *F, unsigned int n_crit, double theta,
   unsigned int exporder);
