
#ifndef _PPSC_OCA_SIGMA_MAT_NOP_HPP
#define _PPSC_OCA_SIGMA_MAT_NOP_HPP

// -----------------------------------------------------------------------
//
// Equilibrium OCA pseudo particle self energy integrator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/data_types.hpp"
#include "integration.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class ARRGA, class ARRGB, class ARRGC, class ARRL1, class ARRL2>
typename ARRGA::matrix_type oca_sdiagram_mat_integral(int m,
					   ARRGA & ga, ARRGB & gb, ARRGC & gc,
					   ARRL1 & lam1, ARRL2 & lam2,
					   int mmax, int kt) {

  using size_type = oca::size_type;
  typedef typename ARRGA::matrix_type matrix_type;

  // -- Integration routines
  using oca::integrals::gregory_sum;
  using oca::integrals::poly_interpol;

  double dtau = 1.0; // -- makes no sense to divide with one everywhere below?
  double tau = m * dtau;

  int m1 = (m >= kt ? m : kt);
  size_type nflavour = ga.get_flavour_shape()[0];

  mam::array_matrix<matrix_type> f1_mat(m1+1, nflavour, nflavour);
  mam::array_matrix<matrix_type> f2_mat(m1+1, nflavour, nflavour);

  matrix_type res_mat(matrix_type::Zero(nflavour, nflavour));

  if(m < kt) { // First kt points, do poly interp

    double dtau1 = tau / kt;

    // -- int dtau2 ...
    for( auto j : range(1, kt+1) ) { // -- Three leg vertex calc G_b L1 G_c

      double tau1 = dtau1 * j;
      double dtau2 = tau1 / kt;

      for( auto l : range(0, kt+1) ) {

	double tau2 = l * dtau2;

	matrix_type gc_interp = poly_interpol(gc, tau2/dtau, kt);
        matrix_type lam1_interp = poly_interpol(lam1, (tau-tau2)/dtau, kt);
        matrix_type gb_interp = poly_interpol(gb, (tau1-tau2)/dtau, kt);

	f2_mat(l).noalias() = gc_interp * lam1_interp * gb_interp;
      }

      matrix_type f2_int = dtau2 * gregory_sum(kt, kt, f2_mat); // auto broken here.. ?

      matrix_type lam2_interp = poly_interpol(lam2, tau1/dtau, kt);
      matrix_type ga_interp = poly_interpol(ga, (tau-tau1)/dtau, kt);

      f1_mat(j).noalias() = ga_interp * lam2_interp * f2_int;
    }

    res_mat.noalias() = dtau1 * gregory_sum(kt, kt, f1_mat);

  } else { // -- Use gregory integration for the rest

    for( auto j : range(1, m+1) ) {

      double tau1 = dtau * j;
      matrix_type z1_mat(matrix_type::Zero(nflavour, nflavour));

      // -- int dtau2 ...
      if(j < kt) {

	double dtau2=tau1/kt;

	for( auto l : range(0, kt+1) ) {

	  double tau2 = l * dtau2;
	  int i1 = floor((tau-tau2)/dtau)-kt/2;
	  if(i1 > mmax - kt) i1 = mmax - kt;
	  if(i1 < 0) i1 = 0;

	  mam::array_matrix_ref<matrix_type> lam1_i1_ref(
	    &lam1[i1][0][0], m1+1, nflavour, nflavour);

	  matrix_type lam1_interp = poly_interpol(lam1_i1_ref, (tau - tau2)/dtau - i1, kt);
	  matrix_type gb_interp = poly_interpol(gb, (tau1-tau2)/dtau, kt);
	  matrix_type gc_interp = poly_interpol(gc, tau2/dtau, kt);

	  f2_mat(l).noalias() = gb_interp * lam1_interp * gc_interp;
	}

	z1_mat.noalias() = dtau2 * gregory_sum(kt, kt, f2_mat);

      } else {

	for(int l=0;l<=j;l++) f2_mat(l).noalias() = gc(l) * lam1(m-l) * gb(j-l);

	z1_mat.noalias() = dtau * gregory_sum(kt, j, f2_mat);

      } // if/else j<kt

      f1_mat(j).noalias() = z1_mat * lam2(j) * ga(m-j);

    } // for j

    res_mat.noalias() = dtau * gregory_sum(kt, m, f1_mat);

  } // if/else m<kt

  return res_mat;
}

// -----------------------------------------------------------------------
void oca_sdiagram_mat_refactor(gf_tstp_type &s, ppgf_type &ga, ppgf_type &gb, ppgf_type &gc,
		      gf_type &lam1, gf_type &lam2, double beta, int kt, int nomp) {
  int mmax=s.ntau();
  int nprocess;
  double dtau;

  dtau = beta / mmax;
  nprocess=nomp;

  auto nflavour = ga.size1();

  assert(nflavour == 1);

  using matrix_type = mam::static_matrix_type<1, 1>;

  mam::array_matrix_ref<matrix_type> ga_ref(ga.matptr(0), ga.ntau()+1, nflavour, nflavour);
  mam::array_matrix_ref<matrix_type> gb_ref(gb.matptr(0), gb.ntau()+1, nflavour, nflavour);
  mam::array_matrix_ref<matrix_type> gc_ref(gc.matptr(0), gb.ntau()+1, nflavour, nflavour);
  mam::array_matrix_ref<matrix_type> lam1_ref(lam1.matptr(0), lam1.ntau()+1, nflavour, nflavour);
  mam::array_matrix_ref<matrix_type> lam2_ref(lam2.matptr(0), lam2.ntau()+1, nflavour, nflavour);

  mam::array_matrix_ref<matrix_type> s_ref(s.matptr(0), s.ntau()+1, nflavour, nflavour);

  s_ref(0) = matrix_type::Zero(nflavour, nflavour);

#pragma omp parallel num_threads(nprocess)
  {
#pragma omp for schedule(dynamic)
    for(int m = mmax; m > 1; m--) {
      s_ref(m) = dtau * dtau
	* oca_sdiagram_mat_integral(m, ga_ref, gb_ref, gc_ref, lam1_ref, lam2_ref, mmax, kt);
    }
  }
  return;
}

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_MAT_NOP_HPP
