
#ifndef _MULTI_ARRAY_MATRIX_COMMON_HPP
#define _MULTI_ARRAY_MATRIX_COMMON_HPP

// -----------------------------------------------------------------------
//
// Views and datatypes for multidimensional arrays of matrices
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------  
#include <cstring> // memcpy
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
// -----------------------------------------------------------------------  

#include "ppsc/cntr_types.hpp"
#include "ppsc/range.hpp"

// -----------------------------------------------------------------------  
namespace ppsc {
namespace mam {
// -----------------------------------------------------------------------

using namespace ppsc;  
  
template<size_type DIMS>
using array_view_type = boost::multi_array_ref<value_type, DIMS>;

template<size_type DIMS>
using sub_array_view_type = boost::detail::multi_array::sub_array<value_type, DIMS>;
  
// -----------------------------------------------------------------------

typedef size_t size_type;

// Underlying data in libcntr is in RowMajor order
// thus for the views to work out correctly we need
// EXPLICITLY to use Eigen::RowMajor

typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dynamic_matrix_type;
typedef Eigen::Map<dynamic_matrix_type> dynamic_matrix_view_type;

// NB! Static matrices with COLS == 1 are required to be ColMajor
// see, https://forum.kde.org/viewtopic.php?f=74&t=98536
template<int ROWS, int COLS>
using static_matrix_type = Eigen::Matrix<value_type, ROWS, COLS,
 (COLS == 1 && ROWS > 1) ? Eigen::ColMajor : Eigen::RowMajor>; // ColMaj if COLS==1

template<int ROWS, int COLS>
using static_matrix_view_type = Eigen::Map<static_matrix_type<ROWS, COLS> >;

typedef dynamic_matrix_type matrix_type;
typedef dynamic_matrix_view_type matrix_view_type;

// -----------------------------------------------------------------------
// -- Helper functions, array concatenation and product reduction

template<typename T, size_type LEN1, size_type LEN2>
boost::array<T, LEN1+LEN2> concatenate_arrays(boost::array<T, LEN1> arr1,
					      boost::array<T, LEN2> arr2) {
  boost::array<T, LEN1+LEN2> out_arr;
  for(size_type i = 0; i < arr1.size(); i++) out_arr[i] = arr1[i];
  for(size_type i = 0; i < arr2.size(); i++) out_arr[i + arr1.size()] = arr2[i];
  return out_arr;
}

template<typename T, size_type LEN>
T product_reduction(const boost::array<T, LEN> arr) {
  T product = 1;
  for( auto & el : arr ) { product *= el; }
  return product;
}

template<typename T, size_type LEN>
void print_array(const boost::array<T, LEN> arr) {
  for( auto & el : arr ) std::cout << el << ", ";
  std::cout << std::endl;
}

// -----------------------------------------------------------------------
} // end namespace mam
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _MULTI_ARRAY_MATRIX_COMMON_HPP
