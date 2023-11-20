#ifndef _PPSC_OCA_GF_MAT_HPP
#define _PPSC_OCA_GF_MAT_HPP

// -----------------------------------------------------------------------
//
// OCA equilibrium single-particle Green's function diagram
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//#include <algorithm> // std::min, std::max
// -----------------------------------------------------------------------

#include "data_types.hpp"
#include "integration.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class DIAGT>
typename DIAGT::g_matrix_type gdiagram_mat_integral(
  int m,
  DIAGT & d,
  typename DIAGT::ga_array_type & ga,
  typename DIAGT::gb_array_type & gb,
  typename DIAGT::gc_array_type & gc,
  typename DIAGT::gd_array_type & gd,
  typename DIAGT::lambda_array_type & lam,
  int mmax, int kt, bool equilibrium=false) {

  // -- matrix types
  using g_matrix_type = typename DIAGT::g_matrix_type;

  using ga_matrix_type = typename DIAGT::ga_matrix_type;
  using gb_matrix_type = typename DIAGT::gb_matrix_type;
  using gc_matrix_type = typename DIAGT::gc_matrix_type;
  using gd_matrix_type = typename DIAGT::gd_matrix_type;

  using ba_matrix_type = typename DIAGT::ba_matrix_type;
  using ca_matrix_type = typename DIAGT::ca_matrix_type;
  using da_matrix_type = typename DIAGT::da_matrix_type;

  using dc_matrix_type = typename DIAGT::dc_matrix_type;
  
  using lambda_matrix_type = typename DIAGT::lambda_matrix_type;

  // -- Helper routines
  using mam::get_view_from;

  // -- Integration routines
  using oca::integrals::gregory_sum;
  using oca::integrals::poly_interpol;

  // -- Temporary storage
  mam::array_matrix<ba_matrix_type> f_ba(mmax+1, d.bdim, d.adim);
  mam::array_matrix<da_matrix_type> f_da(mmax+1, d.ddim, d.adim);
  mam::array_matrix<dc_matrix_type> f_dc(mmax+1, d.ddim, d.cdim);

  assert(kt > 0);

  g_matrix_type g(g_matrix_type::Zero(d.gdim, d.gdim));
  
  // ---------------------------------------------------------------------  
  // -- Complication that makes the routine usable both in equil and real time.
  mam::array_matrix<gd_matrix_type> gd_flip(0, d.ddim, d.ddim); // avoid alloc in real-time
  mam::array_matrix_ref<gd_matrix_type> gd_ref_flip(gd.data(), mmax+1, d.ddim, d.ddim);

  if(equilibrium) {
    
    // -- Re-route gd_ref_flip to point to gd_flip in equilibrium
    new (&gd_flip) mam::array_matrix<gd_matrix_type>(mmax+1, d.ddim, d.ddim);
    new (&gd_ref_flip) mam::array_matrix_ref<gd_matrix_type>(gd_flip.data(), mmax+1, d.ddim, d.ddim);
    
    for(int tau = 0; tau < mmax+1; tau++) {
      gd_ref_flip(tau) = -d.gd.sig() * gd(mmax - tau);
    }
  }
  // ---------------------------------------------------------------------  

  double dtau = 1.0; // beta = 1 ?
  double tau = m * dtau;
  
  if(m < kt) {
    double dtau1 = tau / kt;
    for(int i1 = 0; i1 <= kt; i1++) {
      double tau1 = i1 * dtau1;
      for(int i2 = m; i2 <= mmax; i2++) {
	double tau2 = i2 * dtau;
	double y = (mmax + tau1 - tau2)/dtau; 
	int l = floor(y) - kt/2; 
	if(l > mmax - kt) l = mmax - kt; // l = min(l, mmax-kt)
	if(l < 0) l = 0; // l = max(l, 0)
	
	mam::array_matrix_ref<lambda_matrix_type> lam_view = get_view_from(l, lam);
	lambda_matrix_type lam_interp = poly_interpol(lam_view, y - l, kt);
	f_dc(i2 - m).noalias() =  lam_interp(0,0) * (gd_ref_flip(i2) * d.op1 * gc(i2 - m));
      }

      dc_matrix_type integr_dc = gregory_sum(kt, mmax - m, f_dc);
      ga_matrix_type ga_interp = poly_interpol(ga, tau1/dtau, kt);
      gb_matrix_type gb_interp = poly_interpol(gb, (tau-tau1)/dtau, kt);

      f_da(i1).noalias() = integr_dc * d.op2 * gb_interp * d.op3 * ga_interp;
    }

    da_matrix_type integr_da = gregory_sum(kt, kt, f_da);
    g(0,0) = - dtau1 * (integr_da * d.op4).trace(); // statistical signs??

  }else if(m>mmax-kt){
    double dtau2=(mmax-tau)/kt;
    for(int i2=0;i2<=kt;i2++){
      double tau2=tau+i2*dtau2;
      for(int i1=0;i1<=m;i1++){
	double tau1=i1*dtau;
	double y=(mmax+tau1-tau2)/dtau;
	int l=floor(y)-kt/2;
	if(l>mmax-kt) l=mmax-kt;
	if(l<0) l=0;

	mam::array_matrix_ref<lambda_matrix_type> lam_view = get_view_from(l, lam);
	lambda_matrix_type lam_interp = poly_interpol(lam_view, y-l, kt);
	f_ba(i1).noalias() = lam_interp(0,0) * (gb(m - i1) * d.op3 * ga(i1));
      }
      
      double y=tau2/dtau;
      int l=floor(y)-kt/2;
      if(l>mmax-kt) l=mmax-kt;
      if(l<0) l=0;

      ba_matrix_type integr_ba = gregory_sum(kt, m, f_ba);
      gc_matrix_type gc_interp = poly_interpol(gc, (tau2 - tau)/dtau, kt);
      mam::array_matrix_ref<gd_matrix_type> gd_view = get_view_from(l, gd_ref_flip);
      gd_matrix_type gd_interp = poly_interpol(gd_view, y-l, kt);

      f_da(i2).noalias() = gd_interp * d.op1 * gc_interp * d.op2 * integr_ba;

    }
    da_matrix_type integr_da = gregory_sum(kt, kt, f_da);
    g(0,0) = -dtau2 * (integr_da * d.op4).trace(); // statistical signs??

  }else{
    for(int i1=0;i1<=m;i1++){
      for(int i2=m;i2<=mmax;i2++){
	f_dc(i2 - m).noalias() = lam(mmax - i2 + i1)(0,0) * (gd_ref_flip(i2) * d.op1 * gc(i2 - m));
      }
      dc_matrix_type integr_dc = gregory_sum(kt, mmax - m, f_dc);
      f_da(i1) = integr_dc * d.op2 * gb(m - i1) * d.op3 * ga(i1);
    }
    da_matrix_type integr_da = gregory_sum(kt, m, f_da);
    g(0,0) = -(integr_da * d.op4).trace(); // statistical signs??

  }

  return g;
}

// -----------------------------------------------------------------------
template<class DIAGT>
void gdiagram_mat(DIAGT & d, double beta, int kt, int nomp) {

  int mmax = d.g.ntau() - 1;
  double prefactor = (beta * beta) / (mmax * mmax);
  
  { // scope
#pragma omp parallel for num_threads(nomp)
    for(int m = mmax; m >= 0; m--) {
      d.g.mat(m).noalias() = prefactor
        * gdiagram_mat_integral(
            m, d, d.ga.mat, d.gb.mat, d.gc.mat, d.gd.mat, d.lam.mat, mmax, kt, true);
    } // end for m
  } // end scope
}
  
// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_MAT_HPP
