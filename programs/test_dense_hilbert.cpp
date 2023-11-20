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

#include "./ppsc/hilbert_spaces/single_band_bose.hpp"
#include "./ppsc/hilbert_spaces/single_band_bose_diag.hpp"
#include "./ppsc/hilbert_spaces/single_band_fermi_sc.hpp"
#include "./ppsc/hilbert_spaces/single_band_fermi_diag.hpp"
#include "./ppsc/hilbert_spaces/n_band_fermi_dense.hpp"

// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------
void test_single_band_boson() {
  ppsc::hilbert_spaces::single_band_bose h;
  int nmax = 4;
  h.init(nmax);
  std::cout << h;
}

// -----------------------------------------------------------------------
void test_single_band_boson_diag() {
  ppsc::hilbert_spaces::single_band_bose_diag h;
  int nmax = 4;
  h.init(nmax);
  std::cout << h;
}

// -----------------------------------------------------------------------
void test_single_band_fermi_sc() {
  ppsc::hilbert_spaces::single_band_fermi_sc h;
  h.init();
  std::cout << h;
}

// -----------------------------------------------------------------------
void test_single_band_fermi_diag() {
  ppsc::hilbert_spaces::single_band_fermi_diag h;
  h.init();
  std::cout << h;
}

// -----------------------------------------------------------------------
void test_n_band_fermi_dense(int norb) {
  ppsc::hilbert_spaces::n_band_fermi_dense h;
  h.init(norb);
  std::cout << h;
}

// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  //test_single_band_boson();
  //test_single_band_boson_diag();
  //test_single_band_fermi_sc();
  //test_single_band_fermi_diag();
  //test_n_band_fermi_dense(1);
  test_n_band_fermi_dense(2);
}
