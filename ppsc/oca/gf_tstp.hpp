
#ifndef _PPSC_OCA_GF_TSTP_HPP
#define _PPSC_OCA_GF_TSTP_HPP

// -----------------------------------------------------------------------
//
// OCA real-time single-particle Green's function diagram
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "data_types.hpp"
#include "integration.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------
  
// -----------------------------------------------------------------------
template<class DIAGT>
void gdiagram_tstp(int n, DIAGT & d, double beta, double h, int kt, int nomp) {
    
  // -- helper routines
  using mam::get_column_vector;
  using std::max;
  using std::min;

  // -- Integration routines
  using oca::integrals::gregory_sum;
  using oca::integrals::poly_integral;
  using oca::integrals::matsubara_convolution_lowside_operator;
  using oca::integrals::matsubara_convolution_highside_operator;

  // -- Time boundaries and discretizations
  int mmax = d.g.ntau() - 1;
  int mmax1 = mmax + 1;
  int n1 = max(n, kt);
  int n11 = n1 + 1;
  double dtau = beta / mmax;

  // -- Value types
  using size_type = oca::size_type;
  using value_type = oca::value_type;

  // -- Matrix types
  using g_matrix_type = typename DIAGT::g_matrix_type;
  using lambda_matrix_type = typename DIAGT::lambda_matrix_type;

  using ga_matrix_type = typename DIAGT::ga_matrix_type;
  using gb_matrix_type = typename DIAGT::gb_matrix_type;
  using gc_matrix_type = typename DIAGT::gc_matrix_type;
  using gd_matrix_type = typename DIAGT::gd_matrix_type;

  using ba_matrix_type = typename DIAGT::ba_matrix_type;
  using ca_matrix_type = typename DIAGT::ca_matrix_type;
  using da_matrix_type = typename DIAGT::da_matrix_type;
  using dc_matrix_type = typename DIAGT::dc_matrix_type;
    
  // -- Temporary storage
  mam::matrix_matrix<ga_matrix_type> data_a(n11, mmax1, d.adim, d.adim);
  mam::matrix_matrix<gd_matrix_type> data_d(n11, mmax1, d.ddim, d.ddim);
  
  mam::array_matrix<ga_matrix_type> gavt(mmax1, d.adim, d.adim);
  mam::array_matrix<gd_matrix_type> gdtv(mmax1, d.ddim, d.ddim);

  mam::array_matrix_ref<ga_matrix_type> gavt_ref(gavt.data(), mmax1, d.adim, d.adim);
  mam::array_matrix_ref<gd_matrix_type> gdtv_ref(gdtv.data(), mmax1, d.ddim, d.ddim);
  
  // ---------------------------------------------------------------------
  // -- OpenMP spawn threads
  
#pragma omp parallel num_threads(nomp)
  { // parallell calc

    // -- Thread local temporary storage
    int nmax = max(n11, mmax1);
    mam::array_matrix<ba_matrix_type> f_ba(nmax, d.bdim, d.adim);
    mam::array_matrix<da_matrix_type> f_da(nmax, d.ddim, d.adim);
    mam::array_matrix<dc_matrix_type> f_dc(nmax, d.ddim, d.cdim);
    
    // -------------------------------------------------------------------
    // -- precomp data and data1
    
#pragma omp for schedule(dynamic)
    for (int i2 = 0; i2 < n11; i2++) { 
      for (int m : range(0, mmax1)) {
	data_a(i2, m).noalias() = d.lam.vt(m, i2)(0,0) * d.ga.vt(m, n);
      }
    }
    
#pragma omp for schedule(dynamic)
    for (int i1 = 0; i1 < n11; i1++) { 
      for (int m : range(0, mmax1)) {
	data_d(i1, m).noalias() = d.lam.tv(i1, m)(0,0) * d.gd.tv(n, m);
      }
    }

    // -------------------------------------------------------------------
    // -- pre-fetch for vt component (2)
    
#pragma omp for schedule(dynamic)
    for (int tau = 0; tau < mmax1; tau++) {
      gavt(tau) = d.ga.vt(tau, n);
      gdtv(tau) = d.gd.tv(n, tau);
    }

    // -------------------------------------------------------------------
    // -- vt component g^vt(m,n)
    // -------------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int m = 0; m < mmax1; m++) {

      // -- vt component (1)

      for (int i1 : range(0, n11)) {
        for (int i2 : range(0, n11)) {
	  f_dc(i2) = d.lam.gtr(i1, i2)(0,0) * (d.gd.ret(n, i2) * d.op1 * d.gc.tv(i2, m));
        }
	dc_matrix_type integr_dc = gregory_sum(kt, n, f_dc);
	f_da(i1) = integr_dc * d.op2 * d.gb.vt(m, i1) * d.op3 * d.ga.ret(i1, n);
      }

      da_matrix_type integr_vt1_da = -h*h * gregory_sum(kt, n, f_da);

      // -- vt component (2)
            
      g_matrix_type vt2 = value_type(0., dtau*dtau)
	* gdiagram_mat_integral(
	  m, d, gavt_ref, d.gb.mat, d.gc.mat, gdtv_ref, d.lam.mat, mmax, kt);

      // -- vt component (3)

      for (int i2 : range(0, n11)) {
	mam::array_matrix_ref<ga_matrix_type> view_a = get_column_vector(i2, data_a);
	ba_matrix_type conv_ba =
	  matsubara_convolution_lowside_operator(m, d.gb.mat, d.op3, view_a, mmax, kt);
	f_da(i2).noalias() =  d.gd.ret(n, i2) * d.op1 * d.gc.tv(i2, m) * d.op2 * conv_ba;
      }

      da_matrix_type integr_vt3_da = h*dtau * gregory_sum(kt, n, f_da);

      // -- vt component (4)
      
      for (int i1 : range(0, n11)) {
	mam::array_matrix_ref<gd_matrix_type> view_d = mam::get_column_vector(i1, data_d);
	dc_matrix_type conv_dc =
	  matsubara_convolution_highside_operator(m, view_d, d.op1, d.gc.mat, mmax, kt);
	
        f_da(i1).noalias() = conv_dc * d.op2 * d.gb.vt(m, i1) * d.op3 * d.ga.ret(i1, n);
      }

      da_matrix_type integr_vt4_da = -h*dtau * gregory_sum(kt, n, f_da);
      
      // -- vt result

      d.g.tv(mmax - m) = vt2.adjoint()
	+ ( ((integr_vt1_da + integr_vt3_da + integr_vt4_da) * d.op4).trace()
	* g_matrix_type::Identity() ).adjoint();
      
    } // for m
      
    // -------------------------------------------------------------------
    // -- les component g^les(j,n)
    // -------------------------------------------------------------------
    
#pragma omp for schedule(dynamic)
    for (int j = 0; j < n+1; j++) { // NB! n != n1

      int i2min = min(j, n1 - kt);

      for (int i2 : range(i2min, n11)) {

        // -- les component (1)
	
        for (int i1 : range(0, n11))
	  f_ba(i1).noalias() = d.lam.gtr(i1, i2)(0,0)
	    * (d.gb.les(j, i1) * d.op3 * d.ga.ret(i1, n));

	ba_matrix_type les1_ba = -h * gregory_sum(kt, n, f_ba);

	// -- les component (2)

	for (int i1 : range(0, mmax1))
	  f_ba(i1).noalias() = d.lam.vt(i1, i2)(0,0) *
	    (d.gb.tv(j, i1) * d.op3 * d.ga.vt(i1, n));

	ba_matrix_type les2_ba = value_type(0., -dtau)
	  * gregory_sum(kt, mmax, f_ba);
	
	// -- les component (3)

	int i1max = max(j, kt);
        for (int i1 : range(0, i1max+1))
	  f_ba(i1).noalias() = d.lam.les(i1, i2)(0,0)
	    * (d.gb.ret(j, i1) * d.op3 * d.ga.les(i1, n));
	
	ba_matrix_type les3_ba = h * gregory_sum(kt, j, f_ba);
	
	// -- partial result
	
	f_da(i2 - i2min).noalias() = d.gd.ret(n, i2) * d.op1 * d.gc.ret(i2, j)
	  * d.op2 * (les1_ba + les2_ba + les3_ba);
	
      } // for i2

      // -- les result

      da_matrix_type integr_da = poly_integral(f_da, j - i2min, n - i2min, kt);
      d.g.les(j).noalias() = h * (integr_da * d.op4).trace()
	* g_matrix_type::Identity();
	
    } // for j

    // -------------------------------------------------------------------
    // -- gtr component g^gtr(n,j)
    // -------------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int j = 0; j < n+1; j++) { // NB! n != n1
      
      int i1min = min(j, n1 - kt);

      for (int i1 = i1min; i1 <= n1; i1++) {

	// -- gtr component (1)

        for (int i2 : range(0, n11))
	  f_dc(i2).noalias() =
	    d.lam.les(i1, i2)(0,0) * (d.gd.les(j, i2) * d.op1 * d.gc.ret(i2, n));

	dc_matrix_type gtr1_dc = -h * gregory_sum(kt, n, f_dc);
	
        // -- gtr component (2)

        for (int i2 : range(0, mmax1))
	  f_dc(i2).noalias() =
	    d.lam.tv(i1, i2)(0,0) * (d.gd.tv(j, i2) * d.op1 * d.gc.vt(i2, n));

	dc_matrix_type gtr2_dc = value_type(0., -dtau) * gregory_sum(kt, mmax, f_dc);

        // -- gtr component (3)

	int i2max = max(j, kt);
        for (int i2 : range(0, i2max+1)) 
	  f_dc(i2).noalias() =
	    d.lam.gtr(i1, i2)(0,0) * (d.gd.ret(j, i2) * d.op1 * d.gc.les(i2, n));

	dc_matrix_type gtr3_dc = h * gregory_sum(kt, j, f_dc);

	// -- partial result

	f_da(i1 - i1min) = (gtr1_dc + gtr2_dc + gtr3_dc) * d.op2 * d.gb.ret(n, i1) * d.op3 * d.ga.ret(i1, j);
	
      } // for i1

      // -- ret result

      da_matrix_type integr_da = poly_integral(f_da, j - i1min, n - i1min, kt);
      d.g.ret(j) = h * (integr_da * d.op4).trace() * g_matrix_type::Identity() + d.g.les(j).adjoint(); 

    } // for j

  // ---------------------------------------------------------------------
  } // end omp parallell section
  // ---------------------------------------------------------------------
  
} // oca_gdiagram_tstep_refactor

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_TSTP_HPP
