
#ifndef _TEST_VIEWS_HPP
#define _TEST_VIEWS_HPP

// -----------------------------------------------------------------------
//
// Views and datatypes for multidimensional arrays of matrices
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "multi_array_matrix.hpp"
#include "contour_functions.hpp"

// -----------------------------------------------------------------------  
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// -- Test by comparing views and getters for cntr::herm_matrix
// -----------------------------------------------------------------------

template<class GF>
void compare_views_on_herm_matrix(GF & G) {
//void compare_views_on_herm_matrix(cntr::herm_matrix<double> & G) {

  auto nt = G.nt()+1;
  auto ntau = G.ntau()+1;
  auto nflavour = G.size1();
  
  auto gmat = array_matrix_ref(G.matptr(0), ntau, nflavour);
  auto gles = upper_triangular_matrix_matrix_ref(G.lesptr(0, 0), nt, nflavour);
  auto gles_herm = hermitian_upper_triangular_matrix_matrix_ref(G.lesptr(0, 0), nt, nflavour);
  auto gret = lower_triangular_matrix_matrix_ref(G.retptr(0, 0), nt, nflavour);
  auto gret_herm = hermitian_lower_triangular_matrix_matrix_ref(G.retptr(0, 0), nt, nflavour);
  auto gtv = matrix_matrix_ref(G.tvptr(0, 0), nt, ntau, nflavour);

  // -- Compare tv component
  for( auto tau : irange(0, ntau) ) {
    matrix_type mat; G.get_mat(tau, mat);
    auto diff = std::abs((mat - gmat(tau)).sum()); assert(diff < 1e-9);
  }
  std::cout << "Gmat ok!" << std::endl;

  // -- Compare tv component
  for( auto t : irange(0, nt) ) {
    for( auto tau : irange(0, ntau) ) {
      matrix_type mat; G.get_tv(t, tau, mat);
      auto diff = std::abs((mat - gtv(t, tau)).sum()); assert(diff < 1e-9);      
    }
  }
  std::cout << "Gtv ok!" << std::endl;

  // -- Compare retarded component (lower triangular)
  for(size_type t1 = 0; t1 < nt; t1++) {
    for(size_type t2 = 0; t2 <= t1; t2++) {
      matrix_type mat;
      G.get_ret(t1, t2, mat);
      auto mat_view(gret(t1, t2));
      auto diff = std::abs((mat - mat_view).sum()); assert(diff < 1e-9);
    }
  }

  // -- Compare retarded component (hermitian lower triangular)
  for(size_type t1 = 0; t1 < nt; t1++) {
    for(size_type t2 = 0; t2 < nt; t2++) {
      matrix_type mat;
      G.get_ret(t1, t2, mat);
      auto mat_view(gret_herm(t1, t2));
      auto diff = std::abs((mat - mat_view).sum()); assert(diff < 1e-9);
    }
  }
  std::cout << "Gret_herm ok!" << std::endl;

    // -- Compare lesser component (upper triangular)
  for(size_type t1 = 0; t1 < nt; t1++) {
    for(size_type t2 = t1; t2 < nt; t2++) {
      matrix_type mat;
      G.get_les(t1, t2, mat);
      auto mat_view(gles(t1, t2));
      auto diff = std::abs((mat - mat_view).sum()); assert(diff < 1e-9);
    }
  }
  std::cout << "Gles ok!" << std::endl;

  // -- Compare lesser component (hermitian upper triangular)
  for(size_type t1 = 0; t1 < nt; t1++) {
    for(size_type t2 = 0; t2 < nt; t2++) {
      matrix_type mat;
      G.get_les(t1, t2, mat);
      auto mat_view(gles_herm(t1, t2));
      auto diff = std::abs((mat - mat_view).sum()); assert(diff < 1e-9);
    }
  }
  std::cout << "Gles_herm ok!" << std::endl;
  
}


template<class GF>
void compare_herm_matrix_and_view(GF & G) {

  herm_matrix_matrix_ref Gref(G);
  
  // -- Compare tv component
  for( auto tau : range(0, Gref.ntau()) ) {
    matrix_type mat; G.get_mat(tau, mat);
    auto diff = std::abs((mat - Gref.mat(tau)).sum()); assert(diff < 1e-9);
  }
  std::cout << "Gmat ok!" << std::endl;

  // -- Compare tv component
  for( auto t : range(0, Gref.nt()) ) {
    for( auto tau : range(0, Gref.ntau()) ) {
      matrix_type mat; G.get_tv(t, tau, mat);
      auto diff = std::abs((mat - Gref.tv(t, tau)).sum()); assert(diff < 1e-9);      
    }
  }
  std::cout << "Gtv ok!" << std::endl;

  // -- Compare retarded component (hermitian lower triangular)
  for( auto t1 : range(0, Gref.nt()) ) {
    for( auto t2 : range(0, Gref.nt()) ) {
      matrix_type mat; G.get_ret(t1, t2, mat);
      auto diff = std::abs((mat - Gref.ret(t1, t2)).sum()); assert(diff < 1e-9);
    }
  }
  std::cout << "Gret_herm ok!" << std::endl;

  // -- Compare lesser component (hermitian upper triangular)
  for( auto t1 : range(0, Gref.nt()) ) {
    for( auto t2 : range(0, Gref.nt()) ) {
      matrix_type mat; G.get_les(t1, t2, mat);
      auto diff = std::abs((mat - Gref.les(t1, t2)).sum()); assert(diff < 1e-9);
    }
  }
  std::cout << "Gles_herm ok!" << std::endl;
  
}
  
// -----------------------------------------------------------------------
} // end namespace oca
// -----------------------------------------------------------------------

#endif // _TEST_VIEWS_HPP
