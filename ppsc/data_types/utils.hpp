#ifndef _PPSC_MAM_UTILS_HPP
#define _PPSC_MAM_UTILS_HPP

// -----------------------------------------------------------------------
namespace ppsc {
namespace mam {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------

template<class A, class B> bool compare_matrices(A a, B b, double tol = 1e-9) {

  if(a.get_array_shape().size() != b.get_array_shape().size()) return false;
  if(a.get_array_shape() != b.get_array_shape()) return false;
  if(a.get_flavour_shape() != b.get_flavour_shape()) return false;

  return compare_matrices(a, b, a.get_array_shape()[0], a.get_array_shape()[1], tol);
}

// -----------------------------------------------------------------------
  
template<class A, class B>
bool compare_matrices(A & a, B & b, size_type size1, size_type size2, double tol = 1e-9) {
  
  for( auto t1 : range(0, size1) ) {
    for( auto t2 : range(0, size2) ) {
      auto diff = std::abs((a(t1, t2) - b(t1, t2)).sum());
      if(diff > tol) {
	std::cout << "NOT EQUAL: t1, t2 = " << t1 << ", " << t2 << std::endl;
	std::cout << a(t1, t2) << ", " << b(t1, t2) << std::endl;
	return false;
      }
    }
  }
  return true;
}

// -----------------------------------------------------------------------
  
template<class A, class B> bool compare_arrays(A & a, B & b, double tol = 1e-9) {

  if(a.get_array_shape().size() != b.get_array_shape().size()) {
    std::cout << "--> compare_arrays: array_shape.size() differs!" << std::endl;
    return false;
  }
  if(a.get_array_shape() != b.get_array_shape()) {
    std::cout << "--> compare_arrays: array_shape differs!" << std::endl;    
    return false;
  }
  if(a.get_flavour_shape() != b.get_flavour_shape()) {
    std::cout << "--> compare_arrays: flavour_shape differs!" << std::endl;    
    return false;
  }
  
  for( auto t1 : range(0, a.get_array_shape()[0]) ) {
    auto diff = (a(t1) - b(t1)).cwiseAbs().sum();
    if(diff > tol) {
      std::cout << "NOT EQUAL: t1, a, b, diff = " << t1 << std::endl;
      std::cout << a(t1) << ", " << b(t1) << ", " << diff << std::endl;
      return false;
    }
  }
  return true;
}  
  
// -----------------------------------------------------------------------
} // namespace mam
} // namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_MAM_UTILS_HPP
