
#ifndef _PPSC_OCA_SIGMA_MAT_HPP
#define _PPSC_OCA_SIGMA_MAT_HPP

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
#include "sigma_diag.hpp"
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class DIAGT>
typename DIAGT::sigma_matrix_type sdiagram_mat_integral(
  int m, DIAGT & d,
  typename DIAGT::ga_array_type & ga,
  typename DIAGT::gb_array_type & gb,
  typename DIAGT::gc_array_type & gc,
  typename DIAGT::lambda_array_type & lam1,
  typename DIAGT::lambda_array_type & lam2,
  int mmax, int kt) {

  using sigma_matrix_type = typename DIAGT::sigma_matrix_type;
  using ga_matrix_type = typename DIAGT::ga_matrix_type;
  using gb_matrix_type = typename DIAGT::gb_matrix_type;
  using gc_matrix_type = typename DIAGT::gc_matrix_type;
  using lambda_matrix_type = typename DIAGT::lambda_matrix_type;
  using bc_matrix_type = typename DIAGT::bc_matrix_type;
  using ac_matrix_type = typename DIAGT::ac_matrix_type;

  // -- Integration routines
  using mam::get_view_from;
  using oca::integrals::gregory_sum;
  using oca::integrals::poly_interpol;

  double dtau = 1.0; // -- makes no sense to divide with one everywhere below?
  double tau = m * dtau;

  int m1 = (m >= kt ? m : kt);

  mam::array_matrix<ac_matrix_type> f_ac(m1+1, d.adim, d.cdim);
  mam::array_matrix<bc_matrix_type> f_bc(m1+1, d.bdim, d.cdim);

  sigma_matrix_type res(sigma_matrix_type::Zero(d.sdim, d.sdim));

  if(m == 0) {

    return res;

  } else if(m < kt) { // First kt points, do poly interp

    double dtau1 = tau / kt;

    // -- int dtau2 ...
    for( auto j : range(1, kt+1) ) { // -- Three leg vertex calc G_b L1 G_c

      double tau1 = dtau1 * j;
      double dtau2 = tau1 / kt;

      for( auto l : range(0, kt+1) ) {

	double tau2 = l * dtau2;

	lambda_matrix_type lam1_interp = poly_interpol(lam1, (tau-tau2)/dtau, kt);
        gb_matrix_type gb_interp = poly_interpol(gb, (tau1-tau2)/dtau, kt);
	gc_matrix_type gc_interp = poly_interpol(gc, tau2/dtau, kt);

	f_bc(l).noalias() = lam1_interp(0,0) * (gb_interp * d.op3 * gc_interp);
      }

      bc_matrix_type integr_bc = dtau2 * gregory_sum(kt, kt, f_bc);
      lambda_matrix_type lam2_interp = poly_interpol(lam2, tau1/dtau, kt);
      ga_matrix_type ga_interp = poly_interpol(ga, (tau-tau1)/dtau, kt);
      f_ac(j).noalias() = lam2_interp(0,0) * (ga_interp * d.op2 * integr_bc);

    }

    ac_matrix_type integr_ac = gregory_sum(kt, kt, f_ac);
    res.noalias() = dtau1 * (d.op1 * integr_ac * d.op4);

  } else { // -- Use gregory integration for the rest

    for( auto j : range(1, m+1) ) {

      double tau1 = dtau * j;
      bc_matrix_type integr_bc(bc_matrix_type::Zero(d.bdim, d.cdim));

      // -- int dtau2 ...
      if(j < kt) {

	double dtau2=tau1/kt;

	for( auto l : range(0, kt+1) ) {

	  double tau2 = l * dtau2;
	  int i1 = floor((tau-tau2)/dtau)-kt/2;
	  if(i1 > mmax - kt) i1 = mmax - kt;
	  if(i1 < 0) i1 = 0;

	  mam::array_matrix_ref<lambda_matrix_type> lam1_i1_ref = get_view_from(i1, lam1);

	  lambda_matrix_type lam1_interp = poly_interpol(lam1_i1_ref, (tau - tau2)/dtau - i1, kt);
	  gb_matrix_type gb_interp = poly_interpol(gb, (tau1-tau2)/dtau, kt);
	  gc_matrix_type gc_interp = poly_interpol(gc, tau2/dtau, kt);

	  f_bc(l).noalias() = lam1_interp(0,0) * (gb_interp * d.op3 * gc_interp);
	}

	integr_bc.noalias() = dtau2 * gregory_sum(kt, kt, f_bc);

      } else {

	for(int l = 0; l <= j; l++) f_bc(l).noalias() = lam1(m-l)(0,0) * (gb(j-l) * d.op3 * gc(l));
	integr_bc.noalias() = dtau * gregory_sum(kt, j, f_bc);

      } // if/else j<kt

      f_ac(j).noalias() = lam2(j)(0,0) * (ga(m-j) * d.op2 * integr_bc);

    } // for j

    ac_matrix_type integr_ac = gregory_sum(kt, m, f_ac);
    res.noalias() = dtau * (d.op1 * integr_ac * d.op4);

  } // if/else m<kt

  return res;
}

// -----------------------------------------------------------------------
template<class DIAGT>
void sdiagram_mat(DIAGT & d, double beta, int kt, int nomp) {

  int mmax = d.s.ntau() - 1;
  double dtau = beta / mmax;
  int nprocess = nomp;

  {
 #pragma omp parallel for num_threads(nomp) schedule(dynamic)
    for(int m = 0; m <= mmax; m++) { // include m = 0

      d.s.mat(m).noalias() = dtau * dtau
	* sdiagram_mat_integral(m, d,
				d.ga.mat, d.gb.mat, d.gc.mat,
				d.lam1.mat, d.lam2.mat, mmax, kt);
    }
  }
}

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_MAT_HPP
