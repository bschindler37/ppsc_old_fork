
#ifndef _PPSC_OCA_SIGMA_TSTP_NOP_HPP
#define _PPSC_OCA_SIGMA_TSTP_NOP_HPP

// -----------------------------------------------------------------------
//
// OCA real-time pseudo-particle self-energy diagram integrator
// (version without operators for scalar pseudo-particle Green's funct.)
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#include <algorithm> // std::max, std::min
// -----------------------------------------------------------------------

#include "ppsc/data_types.hpp"
#include "integration.hpp"
#include "sigma_mat_nop.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
void oca_sdiagram_tstep_refactor(
  int n, // current time step
  int sig, // sign of the total diagram, this hardcode and special for one band Hubbard w. diag ppGf
  gf_tstp_type &s_in, // output argument self-energy diagram time-step
  ppgf_type &ga_in, ppgf_type &gb_in, ppgf_type &gc_in, // the three branches of pseud-particle Green's functions
  gf_type &lam1_in, gf_type &lam2_in, // two arcs of hybridization functions
  double beta, double h,
  int kt, // integration order
  int nomp) {

  // -- Helper functions
  using mam::get_column_vector;
  using std::max;
  using std::min;

  // -- Integration routines
  using oca::integrals::convolution2b;
  using oca::integrals::matsubara_convolution_lowside;
  using oca::integrals::gregory_sum;

  // -- Dimensions
  int mmax = s_in.ntau();
  double dtau = beta / mmax;

  // number of time steps to compute (if n < kt then do all steps up to kt)
  int n1 = max(n, kt);
  int n11 = n1 + 1; // number of preceding time steps ?

  // -------------------------------------------------------------------
  // -- Matrix types, generalize according to template parameters!

  using size_type = oca::size_type;
  using value_type = oca::value_type;

  size_type nflavour = 1; // FIXME FROM INPUT!
  using matrix_type = mam::static_matrix_type<1, 1>;
  //using matrix_type = mam::dynamic_matrix_type;

  // -------------------------------------------------------------------
  // Forget whatever was in the pseudo particle self-energy timestep: s

  s_in.clear();
  oca::timeslice_array_matrix_ref<matrix_type> s(s_in);

  // -------------------------------------------------------------------

  oca::herm_matrix_matrix_ref<matrix_type> ga(ga_in);
  oca::herm_matrix_matrix_ref<matrix_type> gb(gb_in);
  oca::herm_matrix_matrix_ref<matrix_type> gc(gc_in);

  oca::herm_matrix_matrix_ref<matrix_type> lam1(lam1_in);
  oca::herm_matrix_matrix_ref<matrix_type> lam2(lam2_in);

  // -------------------------------------------------------------------
  // -- Temporary storage

  mam::matrix_matrix<matrix_type> data(n11, n11, nflavour, nflavour);
  mam::matrix_matrix<matrix_type> data2_tv(n1+1, mmax+1, nflavour, nflavour);
  mam::array_matrix<matrix_type> lam2_full_vt(mmax+1, nflavour, nflavour);
  mam::array_matrix<matrix_type> gc_full_vt(mmax+1, nflavour, nflavour);

  // =====================================================================
  // Component by component integrals
  // =====================================================================

#pragma omp parallel num_threads(nomp)
  { // omp parallel

  // -------------------------------------------------------------------
  // -- ret component (1)

  { OCA_TIMER("MATRIX: ret (1)") // ret (1)

#pragma omp for schedule(dynamic)
    for(int i1 = 0; i1 < n1+1; i1++) {
      int i2max = min(i1 + kt, n1);
      for( int i2 : range(0, i2max+1) ) {
	data(i1, i2).noalias() = ga.ret(n, i1) * gb.ret(i1, i2) * lam1.gtr(n, i2);
      }
    }

#pragma omp for schedule(dynamic)
    for(int j = n; j > -1; j--) {
      int i1min = min(j, n1 - kt);
      int dim = n1 - i1min;
      mam::matrix_matrix<matrix_type> data1(dim+1, dim+1, nflavour, nflavour);
      for( int i1 : range(i1min, n1+1) ) {
	int i2max = min(i1 + kt, n1);
	matrix_type lam2_gtr(lam2.ret(i1, j) + lam2.les(i1, j));
	for( int i2 : range(i1min, i2max+1) ) {
	  data1(i1 - i1min, i2 - i1min).noalias() = data(i1, i2) * gc.ret(i2, j) * lam2_gtr;
	}
      }
      s.ret(j).noalias() = -h*h * convolution2b(data1, j - i1min, n - i1min, kt);
    }

  } // ret (1)

  // -------------------------------------------------------------------
  // -- vt component (1)

  { OCA_TIMER("MATRIX: vt (1)") // vt (1)

#pragma omp for schedule(dynamic)
    for(int t1 = 0; t1 < n1+1; t1++) {
      int i2max = min(t1 + kt, n1);
      for( int t2 : range(0, i2max+1) ) {
	data(t1, t2).noalias() = gb.ret(t2, t1) * lam2.gtr(t2, n) * gc.ret(t1, n);
      }
    }

#pragma omp for schedule(dynamic)
    for(int m = 0; m < mmax+1; m++) {
      mam::matrix_matrix<matrix_type> data1(n1+1, n1+1, nflavour, nflavour);
      for( int i2 : range(0, n1+1) ) {
	int i1max = min(i2 + kt, n1);
	for( int i1 : range(0, i1max+1) ) {
	  data1(i2, i1).noalias() = data(i2, i1) * ga.vt(m, i1) * lam1.vt(m, i2);
	}
      }
      s.tv(mmax-m).noalias() += sig*h*h * convolution2b(data1, 0, n, kt).adjoint();
    }

  } // vt (1)

  // -------------------------------------------------------------------
  // -- les component (1)

  { OCA_TIMER("MATRIX: les (1)") // les (1)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      mam::matrix_matrix<matrix_type> data1(n1+1, n1+1, nflavour, nflavour);
      for( int i2 : range(0, n1+1) ) {
	int i1max = min(i2 + kt, n1);
	for( int i1 : range(0, i1max+1) ) {
	  data1(i2, i1).noalias() = data(i2, i1) * ga.les(j, i1) * lam1.les(j, i2);
	}
      }
      s.les(j).noalias() += -h*h * convolution2b(data1, 0, n, kt);
    }

  } // les (1)

  // -------------------------------------------------------------------
  // -- vt component (2)

  { OCA_TIMER("MATRIX: vt (2)") // vt (2)

#pragma omp for schedule(dynamic)
    for(int t = 0; t < n1+1; t++) {
      for( int tau : range(0, mmax+1) ) {
	data2_tv(t, tau).noalias() = gb.vt(tau, t) * lam2.vt(tau, n);
      }
    }

#pragma omp for schedule(dynamic)
    for(int m = 0; m < mmax+1; m++) {
      mam::array_matrix<matrix_type> f1(n1+1, nflavour, nflavour);
      for( int i2 : range(0, n1+1) ) {
	mam::array_matrix_ref<matrix_type> view = get_column_vector(i2, data2_tv);
	matrix_type z1 = matsubara_convolution_lowside(m, ga.mat, view, mmax, kt);
	f1(i2).noalias() = z1 * lam1.vt(m, i2) * gc.ret(i2, n);
      }
      s.tv(mmax-m).noalias() += -sig * dtau*h * gregory_sum(kt, n, f1).adjoint();
    }

  } // vt (2)

  // -------------------------------------------------------------------
  // -- les component (2)

  { OCA_TIMER("MATRIX: les (2)") // les (2)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      mam::array_matrix<matrix_type> f1(n1+1, nflavour, nflavour);
      mam::array_matrix<matrix_type> f2(mmax+1, nflavour, nflavour);
      for( int i2 : range(0, n1+1) ) {
	for( int m : range(0, mmax+1) ) {
	  f2(m).noalias() = ga.tv(j, m) * data2_tv(i2, m);
	}
	f1(i2).noalias() =  lam1.les(j, i2) * gc.ret(i2, n) * gregory_sum(kt, mmax, f2);
      }
      s.les(j).noalias() += value_type(0.0, -dtau*h) * gregory_sum(kt, n, f1);
    }

  } // les (2)

  // -------------------------------------------------------------------
  // -- vt component (3)

  { OCA_TIMER("MATRIX: vt (3)") // vt (3)

#pragma omp for schedule(dynamic)
    for(int tau = 0; tau < mmax+1; tau++) {
      lam2_full_vt(tau).noalias() = lam2.vt(tau, n);
      gc_full_vt(tau).noalias() = gc.vt(tau, n);
    }

#pragma omp for schedule(dynamic)
    for(int m = 0; m < mmax+1; m++) {
      matrix_type z1 = oca_sdiagram_mat_integral(
        m, ga.mat, gb.mat, gc_full_vt, lam1.mat, lam2_full_vt, mmax, kt);
      s.tv(mmax - m).noalias() += value_type(0.0, -sig*dtau*dtau) * z1.adjoint();
    }

  } // vt (3)

  // -------------------------------------------------------------------
  // -- les component (3)

  { OCA_TIMER("MATRIX: les (3)") // les (3)

#pragma omp for schedule(dynamic)
    for(int j = n; j > -1; j--) {
      int i1max = max(j, kt);
      mam::array_matrix<matrix_type> f1(i1max+1, nflavour, nflavour);
      mam::array_matrix<matrix_type> f2(n1+1, nflavour, nflavour);
      for( int i1 : range(0, i1max+1) ) {
	for( int i2 : range(0, n1+1) ) {
	  f2(i2).noalias() = gb.les(i1, i2) * gc.ret(i2, n) * lam1.les(j, i2);
	}
	f1(i1).noalias() = ga.ret(j, i1) * lam2.les(i1, n) * gregory_sum(kt, n, f2);
      }
      s.les(j).noalias() += h*h * gregory_sum(kt, j, f1);
    }

  } // les (3)

  // -------------------------------------------------------------------
  // -- les component (4)

  { OCA_TIMER("MATRIX: les (4)") // les (4)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      mam::array_matrix<matrix_type> f1(mmax+1, nflavour, nflavour);
      mam::array_matrix<matrix_type> f2(mmax+1, nflavour, nflavour);

      for( int m : range(0, mmax+1) )
	f2(m).noalias() = lam1.tv(j, m) * gc.vt(m, n);

      for( int m : range(0, mmax+1) )
	f1(m).noalias() = ga.tv(j, m) * lam2.vt(m, n)
	  * matsubara_convolution_lowside(m, gb.mat, f2, mmax, kt);

      s.les(j).noalias() += value_type(0.0, dtau*dtau)* gregory_sum(kt, mmax, f1);
    }

  } // les (4)

  // -------------------------------------------------------------------
  // -- les component (5)

  { OCA_TIMER("MATRIX: les (5)") // les (5)
#pragma omp for schedule(dynamic)
    for(int i1 = 0; i1 < n1+1; i1++) {
      int i2max = min(i1 + kt, n1);
      for( auto i2 : range(0, i2max+1) ) {
	data(i1, i2).noalias() = gb.ret(i1, i2) * gc.les(i2, n) * lam2.les(i1, n);
      }
    }

#pragma omp for schedule(dynamic)
    for(int j = n; j > -1; j--) {
      int i1max = max(j, kt);
      mam::matrix_matrix<matrix_type> data1(i1max+1, i1max+1, nflavour, nflavour);
      for( int i1 : range(0, i1max+1) ) {
	int i2max = min(i1 + kt, i1max);
	for( int i2 : range(0, i2max+1) ) {
	  data1(i1, i2).noalias() = data(i1, i2) * ga.ret(j, i1) * lam1.gtr(j, i2);
	}
      }
      s.les(j).noalias() += -h*h * convolution2b(data1, 0, j, kt);
    }

  } // les (5)

  // -------------------------------------------------------------------
  // -- les component (6)

  { OCA_TIMER("MATRIX: les (6)") // les (6)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      int i1max = max(j, kt);
      mam::array_matrix<matrix_type> f1(i1max+1, nflavour, nflavour);
      mam::array_matrix<matrix_type> f2(mmax+1, nflavour, nflavour);
      for( int i1 : range(0, i1max+1) ) {
	for( int m : range(0, mmax+1) ) {
	  f2(m).noalias() = gc.vt(m, n) * lam1.tv(j, m) * gb.tv(i1, m);
	}
	f1(i1).noalias() = ga.ret(j, i1) * lam2.les(i1, n) * gregory_sum(kt, mmax, f2);
      }
      s.les(j).noalias() += oca::value_type(0.0, dtau*h) * gregory_sum(kt, j, f1);
    }

  } // les (6)

  // -------------------------------------------------------------------
  } // omp parallel
  // -------------------------------------------------------------------


} // oca_sdiagram_tstep_refactor

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_TSTP_NOP_HPP
