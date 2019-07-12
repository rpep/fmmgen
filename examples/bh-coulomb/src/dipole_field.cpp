#include<iostream>
#include<iomanip>
#include<fstream>
#include<chrono>
#include<random>
#include "calculate.hpp"
#include "utils.hpp"
#include "tree.hpp"

class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_t> m_beg;

public:
	Timer() : m_beg(clock_t::now())
	{
	}

	void reset()
	{
		m_beg = clock_t::now();
	}

	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};



int main(int argc, char** argv) {
  unsigned int order = 5;
  unsigned int ncrit = 100;
  unsigned int Nparticles = 2;
  std::vector<double> F_exact(4*Nparticles, 0.0);
  std::vector<double> F_approx(4*Nparticles, 0.0);

  std::vector<Particle> particles;
  std::default_random_engine generator(0.0);
  std::uniform_real_distribution<double> distribution(-1, 1);

  auto root = Cell(0.0, 0.0, 0.0, 10.0, 0, order, 0, ncrit);
  std::cout << "Multipole arrays have " << root.M.size() << " elements." << std::endl;
  Particle tmp1(0, 0, 1, -1.0);
  Particle tmp2(0, 0, -1, 1.0);
  particles.push_back(tmp1);
  particles.push_back(tmp2);
  auto cells = build_tree(particles, root, ncrit, order);
  std::cout << "Tree built with " << cells.size() << " cells." << std::endl;
  evaluate_P2M(particles, cells, 0, ncrit, order);

  for(int i = 0; i < cells[0].M.size(); i++) {
    std::cout << "M[" << i << "] = " << cells[0].M[i] << std::endl;
  }

  return 0;
}
