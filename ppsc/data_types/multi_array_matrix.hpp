
#ifndef _MULTI_ARRAY_MATRIX_HPP
#define _MULTI_ARRAY_MATRIX_HPP

// -----------------------------------------------------------------------
//
// Views and datatypes for multidimensional arrays of matrices
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------

#include "common.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace mam {

// -----------------------------------------------------------------------
// -- multi_array_matrix_base implementation

template <size_type array_dimensions, class MAT, class Derived>
class multi_array_matrix_base {

public:
  constexpr static bool has_data = false;
  constexpr static size_type flavour_dimensions = 2;
  constexpr static size_type dimensions = array_dimensions + flavour_dimensions;

  typedef multi_array_matrix_base<array_dimensions, MAT, Derived> this_type;

  typedef std::complex<double> value_type;

  typedef MAT matrix_type;
  typedef Eigen::Map<matrix_type> matrix_view_type;

  typedef array_view_type<dimensions> view_type;
  typedef sub_array_view_type<dimensions - 1> sub_view_type;

  typedef boost::array<size_type, dimensions> shape_type;
  typedef boost::array<size_type, array_dimensions> array_shape_type;
  typedef boost::array<size_type, flavour_dimensions> flavour_shape_type;

  template <size_type D>
  using array_shape_type_template = boost::array<size_type, D>;

  multi_array_matrix_base(array_shape_type array_shape,
                          flavour_shape_type flavour_shape)
      : array_shape(array_shape), flavour_shape(flavour_shape),
        shape(concatenate_arrays(array_shape, flavour_shape)),
        total_size(product_reduction(shape)),
        total_flavour_size(product_reduction(flavour_shape)),
        view_(NULL, shape) {

    assert(MAT::RowsAtCompileTime == -1 ||
           MAT::RowsAtCompileTime == flavour_shape[0]);
    assert(MAT::ColsAtCompileTime == -1 ||
           MAT::ColsAtCompileTime == flavour_shape[1]);
  }

  view_type view() { return view_; }
  sub_view_type operator[](size_type index) { return view_[index]; }

  // -----------------------------------------------------------------------
  // -- Original get_matrix_view

  /*
  matrix_view_type get_matrix_view(array_shape_type array_index) {
    shape_type index(concatenate_arrays(array_index, flavour_shape_type({0,
  0}))); return matrix_view_type(& this->view()(index), flavour_shape[0],
  flavour_shape[1]);
  }
  */

  // -- SFINAE version(s) with 1D and 2D specializations...

  // -- 1D specialization (optimization)
  template <size_type T = array_dimensions>
  matrix_view_type
  get_matrix_view(array_shape_type_template<T> array_index,
                  typename std::enable_if<(T == 1)>::type * = 0) {
    return matrix_view_type(data() + array_index[0] * total_flavour_size,
                            flavour_shape[0], flavour_shape[1]);
  }

  // -- 2D specialization (optimization)
  template <size_type T = array_dimensions>
  matrix_view_type
  get_matrix_view(array_shape_type_template<T> array_index,
                  typename std::enable_if<(T == 2)>::type * = 0) {
    return matrix_view_type(
        data() + (array_shape[1] * array_index[0] + array_index[1]) *
                     total_flavour_size,
        flavour_shape[0], flavour_shape[1]);
  }

  // -- ND general case
  template <size_type T = array_dimensions>
  matrix_view_type
  get_matrix_view(array_shape_type array_index,
                  typename std::enable_if<(T > 2)>::type * = 0) {
    shape_type index(
        concatenate_arrays(array_index, flavour_shape_type({0, 0})));
    return matrix_view_type(&this->view()(index), flavour_shape[0],
                            flavour_shape[1]);
  }

  // -----------------------------------------------------------------------

  const array_shape_type get_array_shape() const { return array_shape; }
  const flavour_shape_type get_flavour_shape() const { return flavour_shape; }

  size_type get_total_size() { return total_size; }

  value_type *data() { return static_cast<Derived *>(this)->data(); }

  void _init_view() { new (&view_) view_type(data(), shape); }

  // -----------------------------------------------------------------------
  // -- inplace scalar multiplication
  this_type &operator*=(value_type scalar) {
    for (int i = 0; i < view_.num_elements(); i++)
      view_.data()[i] *= scalar;
    return *this;
  }

  // -----------------------------------------------------------------------
  // -- check for zero
  bool is_zero(double tol = 1e-12) {
    for (int i = 0; i < view_.num_elements(); i++)
      if (std::abs(view_.data()[i]) > tol)
        return false;
    return true;
  }

  // -----------------------------------------------------------------------
  // -- templated compare
  template <class T> bool is_equal_to(T &other_in, double tol = 1e-12) {
    this_type other(static_cast<this_type>(other_in));
    if (view_.num_elements() != other.view().num_elements())
      return false;
    for (int i = 0; i < view_.num_elements(); i++)
      if (std::abs(view_.data()[i] - other.view().data()[i]) > tol)
        return false;
    return true;
  }

  template <class T> double abs_diff(T &other_in) {
    double diff = 0.0;
    this_type other(static_cast<this_type>(other_in));

    if (view_.num_elements() != other.view().num_elements())
      return -1.0;
    if (view_.num_elements() == 0)
      return 0.0;

    for (int i = 0; i < view_.num_elements(); i++)
      diff += std::abs(view_.data()[i] - other.view().data()[i]);
    diff /= view_.num_elements();

    return diff;
  }

  // -----------------------------------------------------------------------

private:
  array_shape_type array_shape;
  flavour_shape_type flavour_shape;
  shape_type shape;
  size_type total_size;
  size_type total_flavour_size;
  view_type view_;
};

// -----------------------------------------------------------------------
// -- Local buffer

template <int array_dimensions, class MAT>
class multi_array_matrix
    : public multi_array_matrix_base<
          array_dimensions, MAT, multi_array_matrix<array_dimensions, MAT>> {
public:
  constexpr static bool has_data = true;

  typedef multi_array_matrix<array_dimensions, MAT> this_type;
  typedef multi_array_matrix_base<array_dimensions, MAT, this_type> base_type;

  typedef typename base_type::array_shape_type array_shape_type;
  typedef typename base_type::flavour_shape_type flavour_shape_type;

  typedef std::vector<value_type> buffer_type;

  multi_array_matrix(array_shape_type array_shape,
                     flavour_shape_type flavour_shape)
      : base_type(array_shape, flavour_shape),
        buffer(base_type::get_total_size()) {
    base_type::_init_view();
  }

  // -- Move ctor from other matrix type (dynamic <-> static)
  /*
  template<class OTHER_MAT>
  multi_array_matrix(multi_array_matrix<array_dimensions, OTHER_MAT> && other) :
    base_type(other.get_array_shape(), other.get_flavour_shape()),
  buffer(std::move(other.buffer))
  {

    static_assert( !((MAT::RowsAtCompileTime == OTHER_MAT::RowsAtCompileTime) &&
                     (MAT::ColsAtCompileTime == OTHER_MAT::ColsAtCompileTime)),
                   "multi_array_matrix Move ctor from templated matrix used for
  the same type. Should Not Happen." );

    constexpr int dynamic = -1;

    if( (OTHER_MAT::RowsAtCompileTime == dynamic) && (MAT::RowsAtCompileTime !=
  dynamic) ) assert( MAT::RowsAtCompileTime == other.get_flavour_shape()[0] );

    if( (OTHER_MAT::ColsAtCompileTime == dynamic) && (MAT::ColsAtCompileTime !=
  dynamic) ) assert( MAT::ColsAtCompileTime == other.get_flavour_shape()[1] );
  }
    */

  value_type *data() { return buffer.data(); }

private:
  buffer_type buffer;
};

// -----------------------------------------------------------------------
// -- Reference pointer

template <int array_dimensions, class MAT>
class multi_array_matrix_ref
    : public multi_array_matrix_base<
          array_dimensions, MAT,
          multi_array_matrix_ref<array_dimensions, MAT>> {
public:
  typedef multi_array_matrix_ref<array_dimensions, MAT> this_type;
  typedef multi_array_matrix_base<array_dimensions, MAT, this_type> base_type;

  typedef typename base_type::array_shape_type array_shape_type;
  typedef typename base_type::flavour_shape_type flavour_shape_type;

  multi_array_matrix_ref(value_type *ptr, array_shape_type array_shape,
                         flavour_shape_type flavour_shape)
      : base_type(array_shape, flavour_shape), ptr(ptr) {
    base_type::_init_view();
  }

  // -- Copy ctor from other matrix type (dynamic <-> static)
  template <class OTHER_MAT>
  multi_array_matrix_ref(
      const multi_array_matrix_ref<array_dimensions, OTHER_MAT> &other)
      : multi_array_matrix_ref(other.data(), other.get_array_shape(),
                               other.get_flavour_shape()) {

    static_assert(!((static_cast<int>(MAT::RowsAtCompileTime) ==
                     static_cast<int>(OTHER_MAT::RowsAtCompileTime)) &&
                    (static_cast<int>(MAT::ColsAtCompileTime) ==
                     static_cast<int>(OTHER_MAT::ColsAtCompileTime))),
                  "multi_array_matrix_ref template copy ctor used for the same "
                  "type. Should Not Happen.");

    constexpr int dynamic = -1;

    assert(((OTHER_MAT::RowsAtCompileTime == dynamic) &&
            (MAT::RowsAtCompileTime != dynamic))
               ? // if
               (MAT::RowsAtCompileTime ==
                other.get_flavour_shape()[0]) // then check this
               : true && "multi_array_matrix_ref: dynamic to static rows cast "
                         "failed.");

    assert(((OTHER_MAT::ColsAtCompileTime == dynamic) &&
            (MAT::ColsAtCompileTime != dynamic))
               ? // if
               (MAT::ColsAtCompileTime ==
                other.get_flavour_shape()[1]) // then check this
               : true && "multi_array_matrix_ref: dynamic to static cols cast "
                         "failed.");
  }

  // value_type * data() { return ptr; } // -- pointer accessor
  value_type *data() const { return ptr; } // -- pointer accessor

private:
  value_type *ptr;
};

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// -- array type

template <class T> class array_matrix_adaptor : public T {

public:
  typedef T base_type;
  using base_type::base_type;
  typedef typename base_type::matrix_view_type matrix_view_type;

  template <class A = T>
  array_matrix_adaptor(size_type size, size_type rows, size_type cols,
                       typename std::enable_if<A::has_data>::type * = 0)
      : base_type({{size}}, {{rows, cols}}) {}

  template <class A = T>
  array_matrix_adaptor(value_type *ptr, size_type size, size_type rows,
                       size_type cols,
                       typename std::enable_if<!A::has_data>::type * = 0)
      : base_type(ptr, {{size}}, {{rows, cols}}) {}

  /*
  template<class A = T, typename = std::enable_if<A::has_data> >
  array_matrix_adaptor(size_type size, size_type rows, size_type cols) :
    base_type({{size}}, {{rows, cols}}) {}

  template<class A = T, typename = std::enable_if<!A::has_data> >
  array_matrix_adaptor(value_type * ptr, size_type size, size_type rows,
  size_type cols) : base_type(ptr, {{size}}, {{rows, cols}}) {}
  */

  matrix_view_type operator()(size_type index) {
    return base_type::get_matrix_view({{index}});
  }
};

// -----------------------------------------------------------------------

// -- array_matrix
template <class M>
using array_matrix = array_matrix_adaptor<multi_array_matrix<1, M>>;

// -- array_matrix_ref
template <class M>
using array_matrix_ref = array_matrix_adaptor<multi_array_matrix_ref<1, M>>;

// -----------------------------------------------------------------------

template <class ARR>
array_matrix_ref<typename ARR::matrix_type> get_view_from(size_type idx,
                                                          ARR &array) {

  assert(idx >= 0);
  assert(idx < array.get_array_shape()[0]);

  static_assert(ARR::dimensions == 3, "get_view_from called with non-array");

  return array_matrix_ref<typename ARR::matrix_type>(
      &array[idx][0][0], array.get_array_shape()[0] - idx,
      array.get_flavour_shape()[0], array.get_flavour_shape()[1]);
}

// -----------------------------------------------------------------------

template <class MAT>
array_matrix_ref<typename MAT::matrix_type> get_column_vector(size_type idx,
                                                              MAT &matrix) {

  assert(idx >= 0);
  assert(idx < matrix.get_array_shape()[0]);

  static_assert(MAT::dimensions == 4,
                "get_column_vector called with non-matrix");

  return array_matrix_ref<typename MAT::matrix_type>(
      &matrix[idx][0][0][0], matrix.get_array_shape()[1],
      matrix.get_flavour_shape()[0], matrix.get_flavour_shape()[1]);
}

// -----------------------------------------------------------------------
// -- matrix type

template <class T> class matrix_matrix_adaptor : public T {

public:
  typedef T base_type;
  using base_type::base_type;
  typedef typename base_type::matrix_view_type matrix_view_type;

  template <class A = T>
  matrix_matrix_adaptor(size_type size1, size_type size2, size_type rows,
                        size_type cols,
                        typename std::enable_if<A::has_data>::type * = 0)
      : base_type({{size1, size2}}, {{rows, cols}}) {}

  template <class A = T>
  matrix_matrix_adaptor(value_type *ptr, size_type size1, size_type size2,
                        size_type rows, size_type cols,
                        typename std::enable_if<!A::has_data>::type * = 0)
      : base_type(ptr, {{size1, size2}}, {{rows, cols}}) {}

  matrix_view_type operator()(size_type i1, size_type i2) {
    return base_type::get_matrix_view({{i1, i2}});
  }
};

// -----------------------------------------------------------------------

// -- matrix_matrix
template <class M>
using matrix_matrix = matrix_matrix_adaptor<multi_array_matrix<2, M>>;

// -- matrix_matrix_ref
template <class M>
using matrix_matrix_ref = matrix_matrix_adaptor<multi_array_matrix_ref<2, M>>;

// -----------------------------------------------------------------------
// -- triangular types

template <class GETTER> class triangular_matrix_matrix_base : public GETTER {

public:
  typedef GETTER base_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::matrix_view_type matrix_view_type;

  // -- dynamic <-> static conversion
  template <class T>
  triangular_matrix_matrix_base(const T &other) : base_type(other) {}

  template <class A = base_type>
  triangular_matrix_matrix_base(
      size_type size, size_type rows, size_type cols,
      typename std::enable_if<A::has_data>::type * = 0)
      : base_type({{size * (size + 1) / 2}}, {{rows, cols}}) {}

  template <class A = base_type>
  triangular_matrix_matrix_base(
      value_type *ptr, size_type size, size_type rows, size_type cols,
      typename std::enable_if<!A::has_data>::type * = 0)
      : base_type(ptr, {{size * (size + 1) / 2}}, {{rows, cols}}) {}

  // -- Plain triangular matrices return views (no copys) & are assignable
  template <class A = GETTER>
  matrix_view_type operator()(size_type index1, size_type index2,
                              typename std::enable_if<A::is_view>::type * = 0) {
    return base_type::operator()(index1, index2);
  }

  // -- Herm types can not return view, and return a matrix copy instead
  template <class A = GETTER>
  matrix_type
  operator()(size_type index1, size_type index2,
             typename std::enable_if<!A::is_view>::type * = 0) const {
    return base_type::operator()(index1, index2);
  }
};

// -----------------------------------------------------------------------
// -- getter adaptors

template <class M, class IDX>
class plain_getter : public array_matrix_ref<M>, IDX {
public:
  typedef IDX idx_type;
  typedef array_matrix_ref<M> base_type;
  typedef typename base_type::matrix_view_type matrix_view_type;

  using base_type::base_type;

  constexpr static bool is_view = true;

  matrix_view_type operator()(size_type index1, size_type index2) {
    size_type index = idx_type::get_index(index1, index2);
    return base_type::get_matrix_view({{index}});
  }
};

template <class M, class IDX>
class hermitian_getter : public array_matrix_ref<M>, IDX {
public:
  typedef IDX idx_type;
  typedef array_matrix_ref<M> base_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::matrix_view_type matrix_view_type;

  using base_type::base_type;

  constexpr static bool is_view = false;

  matrix_type operator()(size_type index1, size_type index2) const {
    if (idx_type::allowed_range(index1, index2)) {
      size_type index = idx_type::get_index(index1, index2);

      // we do not modify anything so we are safely const...
      // return matrix_type(T::get_matrix_view({{index}}));
      return matrix_type(const_cast<hermitian_getter<M, IDX> *>(this)
                             ->base_type::get_matrix_view({{index}}));

    } else {
      size_type index = idx_type::get_index(index2, index1);

      return -matrix_type(const_cast<hermitian_getter<M, IDX> *>(this)
                              ->base_type::get_matrix_view({{index}}))
                  .adjoint();
    }
  }
};

// -----------------------------------------------------------------------

class upper_triangular_indexer {
public:
  bool allowed_range(size_type index1, size_type index2) const {
    return index1 <= index2;
  }
  size_type get_index(size_type index1, size_type index2) const {
    assert(index1 <= index2);
    return index1 + index2 * (index2 + 1) / 2;
  }
};

class lower_triangular_indexer {
public:
  bool allowed_range(size_type index1, size_type index2) const {
    return index1 >= index2;
  }
  size_type get_index(size_type index1, size_type index2) const {
    assert(index1 >= index2);
    return index1 * (index1 + 1) / 2 + index2;
  }
};

// -----------------------------------------------------------------------

// -- Upper triangular

template <class M>
using upper_triangular_matrix_matrix_ref =
    triangular_matrix_matrix_base<plain_getter<M, upper_triangular_indexer>>;

template <class M>
using hermitian_upper_triangular_matrix_matrix_ref =
    triangular_matrix_matrix_base<
        hermitian_getter<M, upper_triangular_indexer>>;

// -- Lower triangular

template <class M>
using lower_triangular_matrix_matrix_ref =
    triangular_matrix_matrix_base<plain_getter<M, lower_triangular_indexer>>;

template <class M>
using hermitian_lower_triangular_matrix_matrix_ref =
    triangular_matrix_matrix_base<
        hermitian_getter<M, lower_triangular_indexer>>;

// -----------------------------------------------------------------------
// -- vt from tv matrix view
/*
template< class T > class vt_from_tv_matrix_matrix_adaptor : public T {

public:

  typedef T base_type;
  using base_type::base_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::matrix_view_type matrix_view_type;

  // -- Only view constructor
  template<class A = T>
  vt_from_tv_matrix_matrix_adaptor(
    value_type * ptr, size_type nt, size_type ntau, size_type rows, size_type
cols, int sig, typename std::enable_if<!A::has_data>::type * = 0) :
    base_type(ptr, {{nt, ntau}}, {{rows, cols}}), sig(sig) {}

  matrix_type operator() (size_type tau, size_type t) {

    size_type ntau = base_type::get_array_shape()[1];
    assert(tau < ntau);

    size_type rev_tau = ntau - 1 - tau;
    return -sig * matrix_type(base_type::get_matrix_view({{t,
rev_tau}})).adjoint();
  }

private:
  int sig;

};
*/

template <class M> class vt_from_tv_matrix_matrix_ref {

  typedef M matrix_type;
  typedef matrix_matrix_ref<M> tv_type;

public:
  // -- Only view constructor
  vt_from_tv_matrix_matrix_ref(value_type *ptr, size_type nt, size_type ntau,
                               size_type rows, size_type cols, int sig)
      : tv(ptr, nt, ntau, rows, cols), sig(sig) {}

  template <class MOTHER>
  vt_from_tv_matrix_matrix_ref(
      const vt_from_tv_matrix_matrix_ref<MOTHER> &other)
      : tv(other.tv), sig(other.sig) {}

  matrix_type operator()(size_type tau, size_type t) {

    size_type ntau = tv.get_array_shape()[1];
    assert(tau < ntau);

    size_type rev_tau = ntau - 1 - tau;
    return -sig * tv(t, rev_tau).adjoint();
  }

  int sig;
  tv_type tv;
};

// -----------------------------------------------------------------------

// -- vt_from_tv_matrix_matrix_ref
// template<class M> using vt_from_tv_matrix_matrix_ref =
//  vt_from_tv_matrix_matrix_adaptor< multi_array_matrix_ref<2, M> >;

// -----------------------------------------------------------------------
// -- gtr from ret and les matrix view

template <class M> class gtr_from_ret_and_less_matrix_matrix_ref {

  typedef M matrix_type;

  typedef hermitian_lower_triangular_matrix_matrix_ref<M> ret_ref_type;
  typedef hermitian_upper_triangular_matrix_matrix_ref<M> les_ref_type;

public:
  gtr_from_ret_and_less_matrix_matrix_ref(value_type *ret_ptr,
                                          value_type *les_ptr, size_type nt,
                                          size_type rows, size_type cols)
      : ret(ret_ptr, nt, rows, cols), les(les_ptr, nt, rows, cols) {}

  template <class MOTHER>
  gtr_from_ret_and_less_matrix_matrix_ref(
      const gtr_from_ret_and_less_matrix_matrix_ref<MOTHER> &other)
      : ret(other.ret), les(other.les) {}

  const matrix_type operator()(size_type index1, size_type index2) {
    return ret(index1, index2) + les(index1, index2);
  }

  ret_ref_type ret;
  les_ref_type les;
};

// -----------------------------------------------------------------------
} // namespace mam
} // namespace ppsc
// -----------------------------------------------------------------------

#endif // _MULTI_ARRAY_MATRIX_HPP
