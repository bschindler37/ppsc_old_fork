
#ifndef _PPSC_OCA_SIGMA_TSTP_HPP
#define _PPSC_OCA_SIGMA_TSTP_HPP

// -----------------------------------------------------------------------
//
// OCA real-time pseudo-particle self-energy diagram integrator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#include <algorithm> // std::max, std::min
// -----------------------------------------------------------------------

#include "ppsc/data_types.hpp"
#include "integration.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//
// One-crossing approximation self-energy diagram calculator
//
// Diagram structure:
//
// \Sigma = [ O1 GA O2 GB O3 GC O4 ] x L1 L2
// \Sigma(t, t') = [ O_1 G_A(t, t_1) O_2 G_B(t_1, t_2) O_3 G_C(t_2, t') O_4 ] x Lambda_1(t, t_2) Lambda_2(t_1, t')
// with t' fixed at the current time-step (n) on the contour
//
// -----------------------------------------------------------------------
//template<int SDIM, int GADIM, int GBDIM, int GCDIM>
template<class DIAGT>
void sdiagram_tstp(
  int n, // current time step
  DIAGT & d, // oca-diagram instance
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
  int mmax = d.s.ntau() - 1;
  double dtau = beta / mmax;

  // number of time steps to compute (if n < kt then do all steps up to kt)
  int n1 = max(n, kt);
  int n11 = n1 + 1; // number of preceding time steps ?

  // -------------------------------------------------------------------
  // -- Matrix types

  using size_type = oca::size_type;
  using value_type = oca::value_type;

  using sigma_matrix_type = typename DIAGT::sigma_matrix_type;
  using ga_matrix_type = typename DIAGT::ga_matrix_type;
  using gb_matrix_type = typename DIAGT::gb_matrix_type;
  using gc_matrix_type = typename DIAGT::gc_matrix_type;
  using lambda_matrix_type = typename DIAGT::lambda_matrix_type;
  using ab_matrix_type = typename DIAGT::ab_matrix_type;
  using bc_matrix_type = typename DIAGT::bc_matrix_type;
  using ac_matrix_type = typename DIAGT::ac_matrix_type;

  // -------------------------------------------------------------------
  // -- Temporary storage

  mam::matrix_matrix<ab_matrix_type> data_ab(n11, n11, d.adim, d.bdim);
  mam::matrix_matrix<bc_matrix_type> data_bc(n11, n11, d.bdim, d.cdim);

  mam::matrix_matrix<ab_matrix_type> data2_tv_ab(n1+1, mmax+1, d.adim, d.bdim);

  mam::array_matrix<lambda_matrix_type> lam2_full_vt(mmax+1, d.ldim, d.ldim);
  mam::array_matrix_ref<lambda_matrix_type> lam2_full_vt_ref(lam2_full_vt.data(), mmax+1, d.ldim, d.ldim);

  mam::array_matrix<gc_matrix_type> gc_full_vt(mmax+1, d.cdim, d.cdim);
  mam::array_matrix_ref<gc_matrix_type> gc_full_vt_ref(gc_full_vt.data(), mmax+1, d.cdim, d.cdim);

  // -------------------------------------------------------------------

  // =====================================================================
  // Component by component integrals
  // =====================================================================

#pragma omp parallel num_threads(nomp)
  { // omp parallel

  // -- Thread local temporary variables
  mam::matrix_matrix<ac_matrix_type> data_ac(n11, n11, d.adim, d.cdim); // bad scaling nomp * nt^2

  int nmax = max(mmax+1, n1+1);
  mam::array_matrix<ab_matrix_type> f_ab(nmax, d.adim, d.bdim);
  mam::array_matrix<bc_matrix_type> f_bc(nmax, d.bdim, d.cdim);
  mam::array_matrix<ac_matrix_type> f_ac(nmax, d.adim, d.cdim);

  // -------------------------------------------------------------------
  // -- ret component (1)

  { OCA_TIMER("MATRIX: ret (1)") // ret (1)

#pragma omp for schedule(dynamic)
    for(int i1 = 0; i1 < n1+1; i1++) {
      int i2max = min(i1 + kt, n1);
      ga_matrix_type ga_ret = d.ga.ret(n, i1);
      for( int i2 : range(0, i2max+1) ) {
	data_ab(i1, i2).noalias() =  d.lam1.gtr(n, i2)(0,0)
	  * (ga_ret * d.op2 * d.gb.ret(i1, i2));
      }
    }

#pragma omp for schedule(dynamic)
    for(int j = n; j > -1; j--) {
      int i1min = min(j, n1 - kt);
      int dim = n1 - i1min;
      for( int i1 : range(i1min, n1+1) ) {
	int i2max = min(i1 + kt, n1);
	lambda_matrix_type lam2_gtr = d.lam2.gtr(i1, j);
	for( int i2 : range(i1min, i2max+1) ) {
	  data_ac(i1 - i1min, i2 - i1min).noalias() = lam2_gtr(0,0)
	    * (data_ab(i1, i2) * d.op3 * d.gc.ret(i2, j));
	}
      }
      ac_matrix_type conv_ac = convolution2b(data_ac, j - i1min, n - i1min, kt);
      d.s.ret(j).noalias() = -h*h * (d.op1 * conv_ac * d.op4);
    }

  } // ret (1)

  // -------------------------------------------------------------------
  // -- vt component (1)

  { OCA_TIMER("MATRIX: vt (1)") // vt (1)

#pragma omp for schedule(dynamic)
    for(int t1 = 0; t1 < n1+1; t1++) {
      int i2max = min(t1 + kt, n1);
      gc_matrix_type gc_ret = d.gc.ret(t1, n);
      for( int t2 : range(0, i2max+1) ) {
	data_bc(t1, t2).noalias() = d.lam2.gtr(t2, n)(0,0)
	  * (d.gb.ret(t2, t1) * d.op3  * gc_ret);
      }
    }

#pragma omp for schedule(dynamic)
    for(int m = 0; m < mmax+1; m++) {
      for( int i2 : range(0, n1+1) ) {
	int i1max = min(i2 + kt, n1);
	lambda_matrix_type lam1_vt = d.lam1.vt(m, i2);
	for( int i1 : range(0, i1max+1) ) {
	  data_ac(i2, i1).noalias() = lam1_vt(0,0)
	    * (d.ga.vt(m, i1) * d.op2 * data_bc(i2, i1));
	}
      }
      ac_matrix_type conv_ac = convolution2b(data_ac, 0, n, kt);
      d.s.tv(mmax-m).noalias() += d.sig*h*h * (d.op1 * conv_ac * d.op4).adjoint();
    }

  } // vt (1)

  // -------------------------------------------------------------------
  // -- les component (1)

  { OCA_TIMER("MATRIX: les (1)") // les (1)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      for( int i2 : range(0, n1+1) ) {
	int i1max = min(i2 + kt, n1);
	lambda_matrix_type lam1_les = d.lam1.les(j, i2);
	for( int i1 : range(0, i1max+1) ) {
	  data_ac(i2, i1).noalias() = lam1_les(0,0)
	    * (d.ga.les(j, i1) * d.op2 * data_bc(i2, i1));
	}
      }
      ac_matrix_type conv_ac = convolution2b(data_ac, 0, n, kt);
      d.s.les(j).noalias() += -h*h * (d.op1 * conv_ac * d.op4);
    }

  } // les (1)

  // -------------------------------------------------------------------
  // -- vt component (2)

  { OCA_TIMER("MATRIX: vt (2)") // vt (2)

#pragma omp for schedule(dynamic)
    for( int tau = 0; tau < mmax+1; tau++) {
      lambda_matrix_type lam2_vt = d.lam2.vt(tau, n);
      for( int t : range(0, n1+1) ) {
	data2_tv_ab(t, tau).noalias() = lam2_vt(0,0) * (d.op2 * d.gb.vt(tau, t));
      }
    }

#pragma omp for schedule(dynamic)
    for(int m = 0; m < mmax+1; m++) {
      for( int i2 : range(0, n1+1) ) {
	mam::array_matrix_ref<ab_matrix_type> view_ab = get_column_vector(i2, data2_tv_ab);
	ab_matrix_type mat_conv_ab = matsubara_convolution_lowside(m, d.ga.mat, view_ab, mmax, kt);
	f_ac(i2).noalias() = d.lam1.vt(m, i2)(0,0) * (mat_conv_ab * d.op3 * d.gc.ret(i2, n));
      }
      ac_matrix_type integr_ac = gregory_sum(kt, n, f_ac);
      d.s.tv(mmax-m).noalias() += -d.sig * dtau*h * (d.op1 * integr_ac * d.op4).adjoint();
    }

  } // vt (2)

  // -------------------------------------------------------------------
  // -- les component (2)

  { OCA_TIMER("MATRIX: les (2)") // les (2)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      for( int i2 : range(0, n1+1) ) {
	for( int m : range(0, mmax+1) )
	  f_ab(m).noalias() = d.ga.tv(j, m) * data2_tv_ab(i2, m);
	ab_matrix_type integr_ab = gregory_sum(kt, mmax, f_ab);
	f_ac(i2).noalias() =  d.lam1.les(j, i2)(0,0) * (integr_ab * d.op3 * d.gc.ret(i2, n));
      }
      ac_matrix_type integr_ac = gregory_sum(kt, n, f_ac);
      d.s.les(j).noalias() += value_type(0.0, -dtau*h) * (d.op1 * integr_ac * d.op4);
    }

  } // les (2)

  // -------------------------------------------------------------------
  // -- vt component (3)

  { OCA_TIMER("MATRIX: vt (3)") // vt (3)

#pragma omp for schedule(dynamic)
    for(int tau = 0; tau < mmax+1; tau++) {
      lam2_full_vt_ref(tau).noalias() = d.lam2.vt(tau, n);
      gc_full_vt_ref(tau).noalias() = d.gc.vt(tau, n);
    }

#pragma omp for schedule(dynamic)
    for(int m = 0; m < mmax+1; m++) {
      sigma_matrix_type sigma_vt = sdiagram_mat_integral(
        m, d, d.ga.mat, d.gb.mat, gc_full_vt_ref, d.lam1.mat, lam2_full_vt_ref, mmax, kt);
      d.s.tv(mmax - m).noalias() += value_type(0.0, -d.sig*dtau*dtau) * sigma_vt.adjoint();
    }

  } // vt (3)

  // -------------------------------------------------------------------
  // -- les component (3)

  { OCA_TIMER("MATRIX: les (3)") // les (3)

#pragma omp for schedule(dynamic)
    for(int j = n; j > -1; j--) {
      int i1max = max(j, kt);
      for( int i1 : range(0, i1max+1) ) {
	for( int i2 : range(0, n1+1) )
	  f_bc(i2).noalias() =  d.lam1.les(j, i2)(0,0) * (d.gb.les(i1, i2) * d.op3 * d.gc.ret(i2, n));
	bc_matrix_type integr_bc = gregory_sum(kt, n, f_bc);
	f_ac(i1).noalias() = d.lam2.les(i1, n)(0,0) * (d.ga.ret(j, i1) * d.op2 * integr_bc);
      }
      ac_matrix_type integr_ac = gregory_sum(kt, j, f_ac);
      d.s.les(j).noalias() += h*h * (d.op1 * integr_ac * d.op4);
    }

  } // les (3)

  // -------------------------------------------------------------------
  // -- les component (4)

  { OCA_TIMER("MATRIX: les (4)") // les (4)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      for( int m : range(0, mmax+1) )
	f_bc(m).noalias() = d.lam1.tv(j, m)(0,0) * (d.op3 * d.gc.vt(m, n));
      for( int m : range(0, mmax+1) ) {
	bc_matrix_type mat_conv_bc = matsubara_convolution_lowside(m, d.gb.mat, f_bc, mmax, kt);
	f_ac(m).noalias() = d.lam2.vt(m, n)(0,0) * (d.ga.tv(j, m) * d.op2 * mat_conv_bc);
      }
      ac_matrix_type integr_ac = gregory_sum(kt, mmax, f_ac);
      d.s.les(j).noalias() += value_type(0.0, dtau*dtau) * (d.op1 * integr_ac * d.op4);
    }

  } // les (4)

  // -------------------------------------------------------------------
  // -- les component (5)

  { OCA_TIMER("MATRIX: les (5)") // les (5)
#pragma omp for schedule(dynamic)
    for(int i1 = 0; i1 < n1+1; i1++) {
      int i2max = min(i1 + kt, n1);
      lambda_matrix_type lam2_les = d.lam2.les(i1, n);
      for( auto i2 : range(0, i2max+1) ) {
	data_bc(i1, i2).noalias() = lam2_les(0,0) * (d.gb.ret(i1, i2) * d.op3 * d.gc.les(i2, n));
      }
    }

#pragma omp for schedule(dynamic)
    for(int j = n; j > -1; j--) {
      int i1max = max(j, kt);
      for( int i1 : range(0, i1max+1) ) {
	int i2max = min(i1 + kt, i1max);
	ga_matrix_type ga_ret = d.ga.ret(j, i1);
	for( int i2 : range(0, i2max+1) ) {
	  data_ac(i1, i2).noalias() = d.lam1.gtr(j, i2)(0,0) * (ga_ret * d.op2 * data_bc(i1, i2));
	}
      }
      ac_matrix_type conv_ac = convolution2b(data_ac, 0, j, kt);
      d.s.les(j).noalias() += -h*h * (d.op1 * conv_ac * d.op4);
    }

  } // les (5)

  // -------------------------------------------------------------------
  // -- les component (6)

  { OCA_TIMER("MATRIX: les (6)") // les (6)

#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      int i1max = max(j, kt);
      for( int i1 : range(0, i1max+1) ) {
	for( int m : range(0, mmax+1) ) {
	  f_bc(m).noalias() = d.lam1.tv(j, m)(0,0) * (d.gb.tv(i1, m) * d.op3 * d.gc.vt(m, n));
	}
	bc_matrix_type integr_bc = gregory_sum(kt, mmax, f_bc);
	f_ac(i1).noalias() = d.lam2.les(i1, n)(0,0) * (d.ga.ret(j, i1) * d.op2 * integr_bc);
      }
      ac_matrix_type integr_ac = gregory_sum(kt, j, f_ac);
      d.s.les(j).noalias() += oca::value_type(0.0, dtau*h) *  (d.op1 * integr_ac * d.op4);
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

#endif // _PPSC_OCA_SIGMA_TSTP_HPP
