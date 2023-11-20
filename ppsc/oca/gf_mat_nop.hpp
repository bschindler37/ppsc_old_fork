
#ifndef _PPSC_OCA_GF_MAT_NOP_HPP
#define _PPSC_OCA_GF_MAT_NOP_HPP

// -----------------------------------------------------------------------
//
// OCA equilibrium single-particle Green's function diagram
// (version without operators for scalar pseudo-particle Green's funct.)
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
template< class GMAT, class ARRA, class ARRB, class ARRC, class ARRD, class ARRL >
GMAT oca_gdiagram_mat_integral_refactor_matrix(
  int m, int sig, ARRA & ga, ARRB & gb, ARRC & gc, ARRD & gd,
  ARRL & lam, int mmax, int kt, bool equilibrium=false) {

  assert(kt > 0);
  
  // -- Helper routines
  using mam::get_view_from;
  //using std::min;
  //using std::max;

  // -- Integration routines
  using oca::integrals::gregory_sum;
  using oca::integrals::poly_interpol;

  // -- value types
  using matrix_type = GMAT;
  oca::size_type nflavour = 1; // FIXME FOR MATRIX VERSION
  
  // -- return value
  matrix_type g(nflavour, nflavour);

  // -- Temporary storage
  mam::array_matrix<matrix_type> f1(mmax+1, nflavour, nflavour);
  mam::array_matrix<matrix_type> f2(mmax+1, nflavour, nflavour);

  // ---------------------------------------------------------------------  
  // -- Complication that makes the routine usable both in equil and real time.
  mam::array_matrix<matrix_type> gd_flip(0, nflavour, nflavour); // avoid alloc in real-time
  mam::array_matrix_ref<matrix_type> gd_ref_flip(gd.data(), mmax+1, nflavour, nflavour);

  if(equilibrium) {
    // -- Re-route gd_ref_flip to point to gd_flip in equilibrium
    new (&gd_flip) mam::array_matrix<matrix_type>(mmax+1, nflavour, nflavour);
    new (&gd_ref_flip) mam::array_matrix_ref<matrix_type>(gd_flip.data(), mmax+1, nflavour, nflavour);
    
    for(int tau = 0; tau < mmax+1; tau++) {
      gd_ref_flip(tau) = -sig * gd(mmax - tau);
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
	
	mam::array_matrix_ref<matrix_type> view = get_view_from(l, lam);
	f2(i2 - m).noalias() = poly_interpol(view, y - l, kt) * gc(i2 - m) * gd_ref_flip(i2);	
      }
      f1(i1).noalias() = gregory_sum(kt, mmax - m, f2)
	* poly_interpol(ga, tau1/dtau, kt)
	* poly_interpol(gb, (tau-tau1)/dtau, kt);
    }
    
    g = -sig*sig*dtau1 * gregory_sum(kt, kt, f1);
    
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

	mam::array_matrix_ref<matrix_type> view = get_view_from(l, lam);
	f1(i1).noalias() = poly_interpol(view, y-l, kt) * ga(i1) * gb(m - i1);
      }
      
      double y=tau2/dtau;
      int l=floor(y)-kt/2;
      if(l>mmax-kt) l=mmax-kt;
      if(l<0) l=0;

      mam::array_matrix_ref<matrix_type> view = get_view_from(l, gd_ref_flip);
      f2(i2).noalias() = gregory_sum(kt, m, f1)
	* poly_interpol(gc, (tau2 - tau)/dtau, kt)
	* poly_interpol(view, y-l, kt);
    }
    g = -sig*sig*dtau2 * gregory_sum(kt, kt, f2);

  }else{
    for(int i1=0;i1<=m;i1++){
      for(int i2=m;i2<=mmax;i2++){
	f2(i2 - m).noalias() = lam(mmax - i2 + i1) * gc(i2 - m) * gd_ref_flip(i2);
      }
      f1(i1) = gregory_sum(kt, mmax - m, f2) * ga(i1) * gb(m - i1);
    }
    g = -sig*sig * gregory_sum(kt, m, f1);

  }
  return g;
}

// -----------------------------------------------------------------------
void oca_gdiagram_mat_refactor(
  gf_tstp_type &g_in, int sig,
  ppgf_type &ga_in, ppgf_type &gb_in, ppgf_type &gc_in, ppgf_type &gd_in,
  gf_type &lam_in,
  double beta, int kt, int nomp) {

  int nprocess = nomp;
  int mmax = g_in.ntau();
  
  oca::size_type nflavour = 1;
  using matrix_type = mam::static_matrix_type<1, 1>;

  oca::herm_matrix_matrix_ref<matrix_type> ga(ga_in);
  oca::herm_matrix_matrix_ref<matrix_type> gb(gb_in);
  oca::herm_matrix_matrix_ref<matrix_type> gc(gc_in);
  oca::herm_matrix_matrix_ref<matrix_type> gd(gd_in);
  oca::herm_matrix_matrix_ref<matrix_type> lam(lam_in);

  oca::timeslice_array_matrix_ref<matrix_type> g(g_in);

  { // scope
#pragma omp parallel for num_threads(nprocess)    
    for(int m = mmax; m >= 0; m--) {
      g.mat(m).noalias() =
        oca_gdiagram_mat_integral_refactor_matrix<matrix_type>(
          m, sig, ga.mat, gb.mat, gc.mat, gd.mat, lam.mat, mmax, kt, true);
      g.mat(m) *= (beta * beta) / (mmax * mmax);
    } // end for m
  } // end scope

}
  
// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_MAT_NOP_HPP
