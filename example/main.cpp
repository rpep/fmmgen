#include "calculate.hpp"
#include "tree.hpp"
#include "utils.hpp"
#include "variant.hpp"
#include "operators.h"
#include "omp.h"
#include <chrono>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "args.hxx"

#ifdef __linux__
#include <sched.h>
#include <sstream>
#endif

/*! \brief Pin one thread per physical core.

    Worth up to 2.3x. Left to the operating system, threads are placed two to a
    core on the SMT siblings while physical cores sit idle, and the consequence
    is that more threads can be SLOWER: measured on a 2x48-core Sapphire Rapids
    node, 192 unbound threads took 36.3 ms against 26.1 ms for 48, and the first
    parallel region of the solve inflated from 0.18 ms to 7.8 ms.

    This cannot be done through omp.h. OMP_PLACES and OMP_PROC_BIND are
    read-only ICVs fixed before the program starts -- libgomp parses them in a
    constructor that runs before main, so calling setenv here is too late
    (verified: omp_get_proc_bind still reports false afterwards) -- and there is
    deliberately no omp_set_proc_bind in the API. The proc_bind(spread) clause
    exists but is ignored when binding is disabled, which is the default. So the
    placement is done here explicitly.

    Setting OMP_PROC_BIND or OMP_PLACES in the environment disables this and
    hands the decision back to the OpenMP runtime -- including OMP_PROC_BIND=false,
    which is an explicit request for no binding at all.
*/
static void bind_threads_to_cores() {
#ifdef __linux__
  // Any explicit setting wins, whatever it says. Testing omp_get_proc_bind()
  // alone would not honour OMP_PROC_BIND=false, since that reads back the same
  // as "unset".
  if (std::getenv("OMP_PROC_BIND") || std::getenv("OMP_PLACES")) {
    return;
  }
  if (omp_get_proc_bind() != omp_proc_bind_false) {
    return;   // the runtime is already binding; do not fight it
  }

  cpu_set_t allowed;
  CPU_ZERO(&allowed);
  if (sched_getaffinity(0, sizeof allowed, &allowed) != 0) {
    return;
  }

  // One representative CPU per physical core, then the remaining SMT siblings.
  // The representative is the lowest-numbered sibling, read from sysfs rather
  // than assumed: this node enumerates all 96 physical cores before any sibling,
  // but round-robin enumeration is common elsewhere and assuming either one
  // silently recreates the problem this function exists to avoid.
  std::vector<int> cores, siblings;
  for (int cpu = 0; cpu < CPU_SETSIZE; cpu++) {
    if (!CPU_ISSET(cpu, &allowed)) continue;

    std::ifstream f("/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/thread_siblings_list");
    int first = cpu;
    if (f) {
      std::string list;
      std::getline(f, list);
      std::istringstream ss(list);
      int v;
      if (ss >> v) first = v;   // lists are ascending, so the first entry is the core
    }
    (first == cpu ? cores : siblings).push_back(cpu);
  }

  if (cores.empty()) return;

  // Unless the user asked for a specific count, use one thread per physical
  // core: SMT measured worthless here (P2P does not improve at all past the
  // physical core count) and it costs a factor of two when it displaces a core.
  if (!std::getenv("OMP_NUM_THREADS")) {
    omp_set_num_threads((int)cores.size());
  }

  std::vector<int> order = cores;
  order.insert(order.end(), siblings.begin(), siblings.end());

  #pragma omp parallel
  {
    const int t = omp_get_thread_num();
    cpu_set_t one;
    CPU_ZERO(&one);
    CPU_SET(order[t % order.size()], &one);
    // Threads are drawn from libgomp's persistent pool, so the affinity set
    // here survives into every later parallel region at the same team size.
    sched_setaffinity(0, sizeof one, &one);
  }

  std::cout << "Bound " << omp_get_max_threads() << " threads across "
            << cores.size() << " physical cores." << std::endl;
#endif
}

struct Args {
  size_t Nparticles, ncrit, type;
  double theta;
  bool nodirect;
  bool compressed_set;
  int compressed;
  std::string filelabel;
  bool have_filelabel;
};

// D=2 generates particles confined to the z=0 plane (a "planar" tree, see
// tree.hpp) with positions stored 2-wide; D=3 is the ordinary octree case
// with 3-wide positions. Source strengths (FMMGEN_SOURCESIZE-wide) are always
// full, regardless of D -- a planar dipole still has 3 moment components.
template <int D>
int run(const Args &args) {
  const size_t Nparticles = args.Nparticles;
  const size_t ncrit = args.ncrit;
  const double theta = args.theta;
  const size_t type = args.type;

  std::cout << "Scaling Test Parameters" << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << "Dimension  = " << D << (D == 2 ? " (planar)" : "") << std::endl;
  std::cout << "Nparticles = " << Nparticles << std::endl;
  fmm_select(args.compressed_set ? args.compressed != 0 : false);
  std::cout << "operators  = " << fmm->name
            << (fmm_have_compressed() ? "" : " (compressed set not generated)") << std::endl;
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

  // Array containing D-wide positions and FMMGEN_SOURCESIZE-wide source
  // strengths. For D=2 this never allocates or fills a z-component -- a
  // planar run has nothing to zero out, unlike forcing r[3*i+2]=0.0 in a
  // fixed-3-wide array would.
  double *r = new double[D*Nparticles];
  double *S = new double[FMMGEN_SOURCESIZE*Nparticles];
  auto filename = "particles_n_" + std::to_string(Nparticles) + ".txt";
  std::ofstream fout;
  fout.open(filename);
  for (size_t i = 0; i < Nparticles; i++) {
    for(int j = 0; j < D; j++) {
      r[D*i+j] = distribution(generator) * 1e-9;
      fout << r[D*i+j] << ",";
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

  auto base_filename = "_n_" + std::to_string(Nparticles) +
    "_ncrit_" + std::to_string(ncrit) +
    "_theta_" + std::to_string(theta) +
    "_type_" + std::to_string(type) +
    "_d_" + std::to_string(D);
  if (args.have_filelabel) {
    base_filename += "_label_" + args.filelabel;
  }
  base_filename += ".txt";


  auto time_filename = "times" + base_filename;

  std::ofstream timeout(time_filename);

  for (size_t order = FMMGEN_MINORDER; order < FMMGEN_MAXORDER; order++) {
    auto errs_filename = "errors_p_" + std::to_string(order) + "_" + base_filename;
    auto field_filename = "field_p_" + std::to_string(order) + "_" + base_filename;

    // If you're wanting to use the library, this is the part you need to look at!
    // Warning: Check what the standard of your field of study is. Various
    // conventions apply for the sign with regards to dipole and quadrupole
    // orientation!
    Tree<D> tree = build_tree<D>(r, S, Nparticles, ncrit, order, theta);
    std::fill(F_approx.begin(), F_approx.end(), 0);
    if (order == FMMGEN_MINORDER && !args.nodirect) {
      std::cout << "Direct\n-------" << std::endl;
      Timer timer;
      tree.compute_field_exact(F_exact.data());
      t_direct = timer.elapsed();
      std::cout << "t_direct = " << t_direct << std::endl;
      timeout << "direct,"<<t_direct<<std::endl;
    }
    std::cout << "Order " << order << "\n-------" << std::endl;

    #ifdef FMMLIBDEBUG
    std::cout << "Tree built with " << tree.cells.size() << " cells.\n\n\n" << std::endl;
    #endif

    Timer timer;
    // Check the type of simulation, run the appropriate calculation:
    if (type == 0) {
	   tree.compute_field_fmm(F_approx.data());
    }
    else if (type == 1) {
	   tree.compute_field_bh(F_approx.data());
    }
    t_approx = timer.elapsed();
    timeout << order << "," << t_approx << std::endl;

    // If direct calculation is enabled, check the error:
    if (!args.nodirect) {
        std::ofstream errout(errs_filename);

        double errs[FMMGEN_OUTPUTSIZE] = {0.0};
        for (size_t i = 0; i < Nparticles; i++) {
            for(int k = 0; k < FMMGEN_OUTPUTSIZE; k++) {
              double err = (F_exact[FMMGEN_OUTPUTSIZE * i + k] - F_approx[FMMGEN_OUTPUTSIZE * i + k]) / F_exact[FMMGEN_OUTPUTSIZE * i + k];
              fout << err << ",";
              errs[k] += std::abs(err);
          }
          errout << std::endl;
        }

        std::cerr << "Rel errs = " << std::scientific;
        for(int k = 0; k < FMMGEN_OUTPUTSIZE; k++) {
            std::cerr << std::setw(16) << errs[k] / Nparticles << ",";
        }
        std::cout << std::endl;
    }

    std::cout << "Approx. calculation  = " << t_approx << " seconds. " << std::endl;

    // If direct calculation enabled, print the field to a file for checking
    if (!args.nodirect) {
      std::ofstream fieldout(field_filename);
      // Default ostream precision is 6 significant figures, which floors any
      // relative error computed from this file at ~1e-7. At high expansion
      // order the solver is well below that, so the default silently reports
      // the file's rounding rather than the method's accuracy.
      fieldout << std::setprecision(17);
      for (size_t i = 0; i < Nparticles; i++) {
	for(int k = 0; k < FMMGEN_OUTPUTSIZE; k++) {
	  fieldout << F_exact[FMMGEN_OUTPUTSIZE * i + k] << "," << F_approx[FMMGEN_OUTPUTSIZE * i + k] << ",";
	}
	fieldout << std::endl;
      }
    }
  }
  delete[] r;
  delete[] S;
  return 0;
}

int main(int argc, const char **argv) {
  bind_threads_to_cores();

  // Set initial parameters by user input from
  // the command line:
  args::ArgumentParser parser("Fmmgen Example Code.", "This shows how the program scales and prints field,\n"
                                                      "error, and particle position and source strengths\n"
                                                      "to a file");

  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Flag nodirect(parser, "Disable calculate direct", "Disables direct field calculation", {'d', "nodirect"});
  args::ValueFlag<size_t> nparticles(parser, "nparticles", "The total number of particles", {'n', "nparticles"});
  args::ValueFlag<float> thet(parser, "theta", "The opening angle parameter which controls error", {'t', "theta"});
  args::ValueFlag<size_t> nc(parser, "ncrit", "The maximum number of particles in a cell", {"n", "ncrit"});
  args::ValueFlag<size_t> typ(parser, "type", "Type of field evaluation - 0 for FMM and 1 for Barnes-Hut", {"T", "type"});
  args::ValueFlag<std::string> filelabel(parser, "label", "Label for the output files", {"l", "label"});
  args::ValueFlag<int> compressed(parser, "compress", "Use harmonic-compressed operators (0/1)", {"c", "compress"});
  args::ValueFlag<size_t> dim(parser, "dim", "Spatial dimension of particle positions: 2 (planar, z=0) or 3 (default)", {"D", "dim"});

  try
  {
      parser.ParseCLI(argc, argv);
  }
  catch (const args::Completion& e)
  {
      std::cout << e.what();
      return 0;
  }
  catch (const args::Help&)
  {
      std::cout << parser;
      return 0;
  }
  catch (const args::ParseError& e)
  {
      std::cerr << e.what() << std::endl;
      std::cerr << parser;
      return 1;
  }

  Args args_;
  args_.Nparticles = nparticles ? args::get(nparticles) : 10000;
  args_.ncrit = nc ? args::get(nc) : 64;
  args_.theta = thet ? args::get(thet) : 0.4;
  args_.type = typ ? args::get(typ) : 0;
  args_.nodirect = (bool)nodirect;
  args_.compressed_set = (bool)compressed;
  args_.compressed = compressed ? args::get(compressed) : 0;
  args_.have_filelabel = (bool)filelabel;
  args_.filelabel = filelabel ? args::get(filelabel) : "";

  if (args_.type > 1) {
    throw std::runtime_error("Type must be either 0 (Fast Multipole) or 1 (Barnes-Hut)");
  }

  const size_t D = dim ? args::get(dim) : 3;
  if (D == 2) return run<2>(args_);
  if (D == 3) return run<3>(args_);
  throw std::runtime_error("--dim must be 2 or 3");
}
