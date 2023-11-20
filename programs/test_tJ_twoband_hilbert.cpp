// -----------------------------------------------------------------------

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

// -----------------------------------------------------------------------

#ifndef CNTR_USE_OMP
#define CNTR_USE_OMP
#endif

#ifndef CNTR_USE_MPI
#define CNTR_USE_MPI
#endif

#define NCA_SOLVER_ASSERT_0 0
#define NCA_SOLVER_ASSERT_1 0

#include <cntr/cntr.hpp>
#include <cntr/utils/read_inputfile.hpp>

// -----------------------------------------------------------------------

#include "./ppsc/ppsc.hpp"
#include "./ppsc/solver.hpp"

#include "./ppsc/hilbert_spaces/two_band_fermi_tJ.hpp"

// -----------------------------------------------------------------------

using namespace std;
typedef ppsc::operator_type operator_type;

// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  ppsc::hilbert_spaces::two_band_fermi_tJ h;
  h.init();

  int norb=2;
  int spin_degeneracy=2;
  int nf = 2* 2;

  std::vector<operator_type> n_op;
  n_op.resize(nf);
  for (int orb = 0; orb < norb; orb++) {
    for (int spin = 0; spin < spin_degeneracy; spin++) {
      // std::cout << "set c[orb="<<orb<<",spin="<<spin<<"]" <<std::endl;
      h.set_operator_matrix_n(orb, spin, n_op[h.flavor(orb, spin, 0)]);
      std::cout << "density operator for: orb, spin "
	  << orb << ", " << spin << std::endl;
      std::cout << n_op[h.flavor(orb, spin, 0)].to_dense() << std::endl;
    }
  }

  

  std::cout << h;
}
