
#ifndef _PPSC_OCA_GF_TSTP_NOP_HPP
#define _PPSC_OCA_GF_TSTP_NOP_HPP

// -----------------------------------------------------------------------
//
// OCA real-time single-particle Green's function diagram
// (version without operators for scalar pseudo-particle Green's funct.)
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
void oca_gdiagram_tstep_refactor(
  int n,
  gf_tstp_type &g_in, int sigd,
  ppgf_type &ga_in, ppgf_type &gb_in, ppgf_type &gc_in, ppgf_type &gd_in,
  gf_type &lam_in,
  double beta,
  double h,
  int kt,
  int nomp
				 ) {

  // -- helper routines
  using mam::get_column_vector;

  // -- Integration routines
  using oca::integrals::gregory_sum;
  using oca::integrals::poly_integral;
  using oca::integrals::matsubara_convolution_lowside;
  using oca::integrals::matsubara_convolution_highside;

  // -- Value types
  using size_type = oca::size_type;
  using value_type = oca::value_type;
  using matrix_type = mam::static_matrix_type<1, 1>;
  //using matrix_product_type = oca::matrix_product_type;

  size_type nflavour = 1;

  int nprocess = 1, pid;
  nprocess = nomp;

  int mmax = g_in.ntau();
  int mmax1 = mmax + 1;
  int n1 = (n <= kt ? kt : n);
  int n11 = n1 + 1;

  double dtau = beta / mmax;

  // -- Argument wrappers
  g_in.clear(); // clear g_in timestep
  oca::timeslice_array_matrix_ref<matrix_type> g(g_in);

  oca::herm_matrix_matrix_ref<matrix_type> ga(ga_in);
  oca::herm_matrix_matrix_ref<matrix_type> gb(gb_in);
  oca::herm_matrix_matrix_ref<matrix_type> gc(gc_in);
  oca::herm_matrix_matrix_ref<matrix_type> gd(gd_in);
  oca::herm_matrix_matrix_ref<matrix_type> lam(lam_in);

  // -- Temporary storage
  mam::matrix_matrix<matrix_type> data(n11, mmax1, nflavour, nflavour);
  mam::matrix_matrix<matrix_type> data1(n11, mmax1, nflavour, nflavour);
  
  mam::array_matrix<matrix_type> gavt(mmax1, nflavour, nflavour);
  mam::array_matrix<matrix_type> gdtv(mmax1, nflavour, nflavour);
  
  // ---------------------------------------------------------------------
  // -- OpenMP spawn threads
  
#pragma omp parallel num_threads(nprocess)
  { // parallell calc

    // -- Thread local temporary storage
    mam::array_matrix<matrix_type> f1(n11 + mmax1, nflavour, nflavour);
    mam::array_matrix<matrix_type> f2(n11 + mmax1, nflavour, nflavour);
    
    // -------------------------------------------------------------------
    // -- precomp data and data1
    
    // data[i2*mmax1+m] = ga^vt(m,n)*lam^vt(m,i2)
  
#pragma omp for schedule(dynamic)
    for (int i2 = 0; i2 <= n1; i2++) { 
      for (int m = 0; m <= mmax; m++) {
	data(i2, m).noalias() = lam.vt(m, i2) * ga.vt(m, n);
      }
    }
    
    // data1[i1*mmax1+m] = gd^tv(n,m)*lam^tv(i1,m)
    
#pragma omp for schedule(dynamic)
    for (int i1 = 0; i1 <= n1; i1++) { 
      for (int m = 0; m <= mmax; m++) {
	data1(i1, m).noalias() = gd.tv(n, m) * lam.tv(i1, m);
      }
    }
    
    // -------------------------------------------------------------------
    // -- pre-fetch for vt component (2)
    
#pragma omp for schedule(dynamic)
    for (int tau = 0; tau <= mmax; tau++) {
      gavt(tau) = ga.vt(tau, n);
      gdtv(tau) = gd.tv(n, tau);
    }
    
    // -------------------------------------------------------------------
    // -- vt component g^vt(m,n)
    // -------------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int m = 0; m <= mmax; m++) {
      
      // -- vt component (1)
      for (int i1 = 0; i1 <= n1; i1++) {
        for (int i2 = 0; i2 <= n1; i2++) {
	  f2(i2) = lam.gtr(i1, i2) * gd.ret(n, i2) * gc.tv(i2, m);
        }
	f1(i1) = gregory_sum(kt, n, f2)
	  * gb.vt(m, i1) * ga.ret(i1, n);
      }
      
      matrix_type vt1 = -h*h * gregory_sum(kt, n, f1);

      // -- vt component (2)
            
      matrix_type vt2 = value_type(0., 1.) * dtau * dtau
	* oca_gdiagram_mat_integral_refactor_matrix<matrix_type>(
          m, sigd, gavt, gb.mat, gc.mat, gdtv, lam.mat, mmax, kt);

      // -- vt component (3)

      for (int i2 = 0; i2 <= n1; i2++) {
	mam::array_matrix_ref<matrix_type> data_view = get_column_vector(i2, data);
	f2(i2).noalias() = matsubara_convolution_lowside(
          m, gb.mat, data_view, mmax, kt) * gd.ret(n, i2) * gc.tv(i2, m);
      }

      matrix_type vt3 = h * dtau * gregory_sum(kt, n, f2);

      // -- vt component (4)
      
      for (int i1 = 0; i1 <= n1; i1++) {
	mam::array_matrix_ref<matrix_type> data_view = mam::get_column_vector(i1, data1);
        f1(i1).noalias() =
	  matsubara_convolution_highside(m, gc.mat, data_view, mmax, kt)
	  * ga.ret(i1, n) * gb.vt(m, i1);
      }

      matrix_type vt4 = -h * dtau * gregory_sum(kt, n, f1);

      // -- vt result

      g.tv(mmax - m) = (vt1 + vt2 + vt3 + vt4).adjoint();
      
    } // for m

    // -------------------------------------------------------------------
    // -- les component g^les(j,n)
    // -------------------------------------------------------------------
    
#pragma omp for schedule(dynamic)
    for (int j = 0; j <= n; j++) {

      int i2min = n1 - kt;
      if (i2min > j)
        i2min = j;

      for (int i2 = i2min; i2 <= n1; i2++) {

        // -- les component (1)
	
        for (int i1 = 0; i1 <= n1; i1++)
	  f1(i1).noalias() = lam.gtr(i1, i2) * ga.ret(i1, n) * gb.les(j, i1);

	matrix_type les1 = -h * gregory_sum(kt, n, f1);

	// -- les component (2)

        for (int i1 = 0; i1 <= mmax; i1++)
	  f1(i1).noalias() = lam.vt(i1, i2) * ga.vt(i1, n) * gb.tv(j, i1);

	matrix_type les2 = -dtau * value_type(0., 1.) * gregory_sum(kt, mmax, f1);
	
	// -- les component (3)

        int i1max = (j >= kt ? j : kt); // max(j, kt)
        for (int i1 = 0; i1 <= i1max; i1++)
	  f1(i1).noalias() = lam.les(i1, i2) * ga.les(i1, n) * gb.ret(j, i1);
	
	matrix_type les3 = h * gregory_sum(kt, j, f1);
	
	// -- partial result
	
	f2(i2 - i2min).noalias() = (les1 + les2 + les3)
	  * gc.ret(i2, j) * gd.ret(n, i2);
	
      } // for i2

      // -- les result
      
      g.les(j).noalias() = h * poly_integral(f2, j - i2min, n - i2min, kt);
	
    } // for j

    // -------------------------------------------------------------------
    // -- gtr component g^gtr(n,j)
    // -------------------------------------------------------------------

#pragma omp for schedule(dynamic)
    for (int j = 0; j <= n; j++) {
      
      int i1min = n1 - kt;
      if (i1min > j)
        i1min = j;

      for (int i1 = i1min; i1 <= n1; i1++) {

	// -- gtr component (1)

        for (int i2 = 0; i2 <= n1; i2++)
	  f2(i2).noalias() =
	    lam.les(i1, i2) * gc.ret(i2, n) * gd.les(j, i2);

	matrix_type gtr1 = -h * gregory_sum(kt, n, f2);
	
        // -- gtr component (2)

        for (int i2 = 0; i2 <= mmax; i2++)
	  f2(i2).noalias() =
	    lam.tv(i1, i2) * gc.vt(i2, n) * gd.tv(j, i2);

	matrix_type gtr2 = value_type(0., -dtau) * gregory_sum(kt, mmax, f2);

        // -- gtr component (3)

        int i2max = (j >= kt ? j : kt); // max(j, kt)
        for (int i2 = 0; i2 <= i2max; i2++) 
	  f2(i2).noalias() =
	    lam.gtr(i1, i2) * gc.les(i2, n) * gd.ret(j, i2);

	matrix_type gtr3 = h * gregory_sum(kt, j, f2);

	// -- partial result

	f1(i1 - i1min) = (gtr1 + gtr2 + gtr3) * ga.ret(i1, j) * gb.ret(n, i1);
	
      } // for i1

      // -- ret result
      
      g.ret(j) = h * poly_integral(f1, j - i1min, n - i1min, kt)
	+ g.les(j).adjoint(); 

    } // for j

  // ---------------------------------------------------------------------
  } // end omp parallell section
  // ---------------------------------------------------------------------

} // oca_gdiagram_tstep_refactor

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_TSTP_NOP_HPP
 
