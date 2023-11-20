
#ifndef _CONTOUR_FUNCTIONS_HPP
#define _CONTOUR_FUNCTIONS_HPP

// -----------------------------------------------------------------------
//
// Views and datatypes for multidimensional arrays of matrices
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "multi_array_matrix.hpp"
#include "utils.hpp"

// -----------------------------------------------------------------------  
namespace ppsc {
namespace cntr {
// -----------------------------------------------------------------------

class contour_base {

public:

  contour_base(size_type nt, size_type ntau, size_type nflavour, int sig) :
    nt_(nt), ntau_(ntau), nflavour_(nflavour), sig_(sig) {}

  const size_type nt() const { return nt_; }
  const size_type ntau() const { return ntau_; }
  const size_type nflavour() const { return nflavour_; }
  const int sig() const { return sig_; }

protected:
  
  const size_type nt_, ntau_, nflavour_;
  const int sig_;

};

// -----------------------------------------------------------------------

template< class M >
class herm_matrix_matrix_ref : public contour_base {

  typedef ppsc::mam::array_matrix_ref<M> mat_ref_type;
  typedef ppsc::mam::hermitian_lower_triangular_matrix_matrix_ref<M> ret_ref_type;
  typedef ppsc::mam::hermitian_upper_triangular_matrix_matrix_ref<M> les_ref_type;

  typedef ppsc::mam::lower_triangular_matrix_matrix_ref<M> lt_ret_ref_type;
  typedef ppsc::mam::upper_triangular_matrix_matrix_ref<M> ut_les_ref_type;
  
  typedef ppsc::mam::gtr_from_ret_and_less_matrix_matrix_ref<M> gtr_ref_type;
  typedef ppsc::mam::matrix_matrix_ref<M> tv_ref_type;
  typedef ppsc::mam::vt_from_tv_matrix_matrix_ref<M> vt_ref_type;
  
public:

  /*
  // -- This hides default constructor...
  template<class HERM_MATRIX> herm_matrix_matrix_ref(HERM_MATRIX & G) :
    contour_base(G.nt()+1, G.ntau()+1, G.size1(), G.sig()),
    mat(G.matptr(0), ntau_, nflavour_, nflavour_),
    ret(G.retptr(0, 0), nt_, nflavour_, nflavour_),
    ret_lt(G.retptr(0, 0), nt_, nflavour_, nflavour_),
    les(G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    les_ut(G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    gtr(G.retptr(0, 0), G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    tv(G.tvptr(0, 0), nt_, ntau_, nflavour_, nflavour_),
    vt(G.tvptr(0, 0), nt_, ntau_, nflavour_, nflavour_, G.sig())
  {}
  */

  herm_matrix_matrix_ref(ppgf_type & G) :
    contour_base(G.nt()+1, G.ntau()+1, G.size1(), G.sig()),
    mat(G.matptr(0), ntau_, nflavour_, nflavour_),
    ret(G.retptr(0, 0), nt_, nflavour_, nflavour_),
    ret_lt(G.retptr(0, 0), nt_, nflavour_, nflavour_),
    les(G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    les_ut(G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    gtr(G.retptr(0, 0), G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    tv(G.tvptr(0, 0), nt_, ntau_, nflavour_, nflavour_),
    vt(G.tvptr(0, 0), nt_, ntau_, nflavour_, nflavour_, G.sig()) {

    assert( (M::RowsAtCompileTime == -1) || (nflavour_ == M::RowsAtCompileTime) );
    assert( (M::ColsAtCompileTime == -1) || (nflavour_ == M::ColsAtCompileTime) );
  }

  herm_matrix_matrix_ref(gf_type & G) :
    contour_base(G.nt()+1, G.ntau()+1, G.size1(), G.sig()),
    mat(G.matptr(0), ntau_, nflavour_, nflavour_),
    ret(G.retptr(0, 0), nt_, nflavour_, nflavour_),
    ret_lt(G.retptr(0, 0), nt_, nflavour_, nflavour_),
    les(G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    les_ut(G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    gtr(G.retptr(0, 0), G.lesptr(0, 0), nt_, nflavour_, nflavour_),
    tv(G.tvptr(0, 0), nt_, ntau_, nflavour_, nflavour_),
    vt(G.tvptr(0, 0), nt_, ntau_, nflavour_, nflavour_, G.sig()) {

    assert( (M::RowsAtCompileTime == -1) || (nflavour_ == M::RowsAtCompileTime) );
    assert( (M::ColsAtCompileTime == -1) || (nflavour_ == M::ColsAtCompileTime) );
  }

  template<class MOTHER>
  herm_matrix_matrix_ref(const herm_matrix_matrix_ref<MOTHER> & other) :
    contour_base(other.nt(), other.ntau(), other.nflavour(), other.sig()),
    mat(other.mat), ret(other.ret), ret_lt(other.ret_lt),
    les(other.les), les_ut(other.les_ut), gtr(other.gtr), tv(other.tv), vt(other.vt) {}
  
public:

  mat_ref_type mat;
  ret_ref_type ret;
  lt_ret_ref_type ret_lt;
  les_ref_type les;
  ut_les_ref_type les_ut;
  gtr_ref_type gtr;
  tv_ref_type tv;
  vt_ref_type vt;

};

// -----------------------------------------------------------------------
// -- own data and assignable... ?
  
template< class M >
class herm_matrix_matrix : public contour_base {

  typedef ppsc::mam::array_matrix<M> mat_type;
  typedef ppsc::mam::matrix_matrix<M> ret_type;
  typedef ppsc::mam::matrix_matrix<M> les_type;
  typedef ppsc::mam::matrix_matrix<M> tv_type;

  typedef ppsc::mam::hermitian_lower_triangular_matrix_matrix_ref<M> ret_ref_type;
  typedef ppsc::mam::hermitian_upper_triangular_matrix_matrix_ref<M> les_ref_type;
  typedef ppsc::mam::gtr_from_ret_and_less_matrix_matrix_ref<M> gtr_ref_type;
  typedef ppsc::mam::vt_from_tv_matrix_matrix_ref<M> vt_ref_type;
  
public:

  herm_matrix_matrix(
    size_type nt, size_type ntau, size_type nflavour, int sig) :
    contour_base(nt, ntau, nflavour, sig),
    mat(ntau_, nflavour_, nflavour_),
    ass_ret(nt_, nflavour_, nflavour_),
    ret(ass_ret.data(), nt_, nflavour_, nflavour_),
    ass_les(nt_, nflavour_, nflavour_),
    les(ass_les.data(), nt_, nflavour_, nflavour_),
    gtr(ass_ret.data(), ass_les.data(), nt_, nflavour_, nflavour_),
    tv(nt_, ntau_, nflavour_, nflavour_),
    vt(tv.data(), nt_, ntau_, nflavour_, nflavour_, sig)
  {}

public:

  mat_type mat;
  ret_type ass_ret;
  ret_ref_type ret;
  les_type ass_les;
  les_ref_type les;
  gtr_ref_type gtr;
  tv_type tv;
  vt_ref_type vt;

};
  
// -----------------------------------------------------------------------

template< class M >
class full_matrix_matrix : public contour_base {

  typedef ppsc::mam::array_matrix<M> mat_type;
  typedef ppsc::mam::matrix_matrix<M> ret_type;
  typedef ppsc::mam::matrix_matrix<M> les_type;
  typedef ppsc::mam::matrix_matrix<M> gtr_type;
  typedef ppsc::mam::matrix_matrix<M> tv_type;
  typedef ppsc::mam::matrix_matrix<M> vt_type;
  
public:

  full_matrix_matrix(size_type nt, size_type ntau, size_type nflavour, int sig) :
    contour_base(nt, ntau, nflavour, sig),
    mat(ntau_, nflavour_, nflavour_),
    ret(nt_, nt_, nflavour_, nflavour_),
    les(nt_, nt_, nflavour_, nflavour_),
    gtr(nt_, nt_, nflavour_, nflavour_),
    tv(nt_, ntau_, nflavour_, nflavour_),
    vt(ntau_, nt_, nflavour_, nflavour_)
  {}

  full_matrix_matrix(herm_matrix_matrix_ref<M> & G) :
    contour_base(G.nt(), G.ntau(), G.nflavour(), G.sig()),
    mat(ntau_, nflavour_, nflavour_),
    ret(nt_, nt_, nflavour_, nflavour_),
    les(nt_, nt_, nflavour_, nflavour_),
    gtr(nt_, nt_, nflavour_, nflavour_),
    tv(nt_, ntau_, nflavour_, nflavour_),
    vt(ntau_, nt_, nflavour_, nflavour_)
  {

    for( auto tau : range(0, ntau()) ) {
      mat(tau) = G.mat(tau);
    }

    for( auto t1 : range(0, nt()) ) {
      for( auto t2 : range(0, nt()) ) {
	ret(t1, t2) = G.ret(t1, t2);
	les(t1, t2) = G.les(t1, t2);
	gtr(t1, t2) = G.gtr(t1, t2);
      }
    }

    for( auto t : range(0, nt()) ) {
      for( auto tau : range(0, ntau()) ) {
	tv(t, tau) = G.tv(t, tau);
	vt(tau, t) = G.vt(tau, t);	  
      }
    }
   
  }

public:

  mat_type mat;
  ret_type ret;
  les_type les;
  gtr_type gtr;
  tv_type tv;
  vt_type vt;

};
  
// -----------------------------------------------------------------------
template< class GA, class GB>
bool compare_greens_functions(GA & ga, GB & gb, double tol = 1e-9) {
  assert( ga.nt() == gb.nt() );
  assert( ga.ntau() == gb.ntau() );
  assert( ga.nflavour() == gb.nflavour() );

  assert( compare_arrays(ga.mat, gb.mat, tol=tol) );
  assert( compare_matrices(ga.ret, gb.ret, ga.nt(), ga.nt(), tol=tol) );
  assert( compare_matrices(ga.les, gb.les, ga.nt(), ga.nt(), tol=tol) );
  assert( compare_matrices(ga.gtr, gb.gtr, ga.nt(), ga.nt(), tol=tol) );
  assert( compare_matrices(ga.tv, gb.tv, ga.nt(), ga.ntau(), tol=tol) );
  assert( compare_matrices(ga.vt, gb.vt, ga.ntau(), ga.nt(), tol=tol) );
}

// -----------------------------------------------------------------------

template<class M>
class full_array_matrix : public contour_base {

  typedef ppsc::mam::array_matrix<M> gtr_type;
  typedef ppsc::mam::array_matrix<M> les_type;
  typedef ppsc::mam::array_matrix<M> vt_type;
  
public:

  full_array_matrix(size_type nt, size_type ntau, size_type nflavour, int sig) :
    contour_base(nt, ntau, nflavour, sig),
    gtr(nt_, nflavour_, nflavour_),
    les(nt_, nflavour_, nflavour_),
    vt(ntau_, nflavour_, nflavour_)
  {}

public:

  gtr_type gtr;
  les_type les;
  vt_type vt;

};
  
// -----------------------------------------------------------------------

template<class M>
class timeslice_array_matrix_ref : public contour_base {

  typedef ppsc::mam::array_matrix_ref<M> ret_type;
  typedef ppsc::mam::array_matrix_ref<M> les_type;
  typedef ppsc::mam::array_matrix_ref<M> tv_type;
  typedef ppsc::mam::array_matrix_ref<M> mat_type;
  
public:

  /*
  template<class SLICE>
  timeslice_array_matrix_ref(SLICE & slice) :
    contour_base(slice.tstp_+1, slice.ntau()+1, slice.size1(), slice.sig_),
    ret(slice.retptr(0), nt_, nflavour_, nflavour_),
    les(slice.lesptr(0), nt_, nflavour_, nflavour_),
    tv(slice.tvptr(0), ntau_, nflavour_, nflavour_),
    mat(slice.matptr(0), ntau_, nflavour_, nflavour_)
  {}
  */

  timeslice_array_matrix_ref(gf_tstp_type & slice) :
    contour_base(slice.tstp_+1, slice.ntau()+1, slice.size1(), slice.sig_),
    ret(slice.retptr(0), nt_, nflavour_, nflavour_),
    les(slice.lesptr(0), nt_, nflavour_, nflavour_),
    tv(slice.tvptr(0), ntau_, nflavour_, nflavour_),
    mat(slice.matptr(0), ntau_, nflavour_, nflavour_)
  {}

  template<class MOTHER>
  timeslice_array_matrix_ref(const timeslice_array_matrix_ref<MOTHER> & other) :
    contour_base(other.nt(), other.ntau(), other.nflavour(), other.sig()),
    ret(other.ret), les(other.les), tv(other.tv), mat(other.mat) {}

  // -- element-wise inplace scalar multiplication
  timeslice_array_matrix_ref<M> & operator*=(ppsc::mam::value_type scalar) {
    tv *= scalar; ret *= scalar; les *= scalar; return *this; }

  bool is_zero(double tol=1e-12) {
    return ret.is_zero(tol) && les.is_zero(tol) && tv.is_zero(tol) && mat.is_zero(tol); }

  template<class MOTHER>
  bool is_equal_to(const timeslice_array_matrix_ref<MOTHER> & other, double tol=1e-12) {
    return ret.is_equal_to(other.ret, tol) && les.is_equal_to(other.les, tol)
      && tv.is_equal_to(other.tv, tol) && mat.is_equal_to(other.mat, tol);
  }

  template<class MOTHER>
  double abs_diff(const timeslice_array_matrix_ref<MOTHER> & other) {
    double diff = ret.abs_diff(other.ret) + les.abs_diff(other.les) +
      + tv.abs_diff(other.tv) + mat.abs_diff(other.mat);
    return diff;
  }
  
public:

  ret_type ret;
  les_type les;
  tv_type tv;
  mat_type mat;

};

// -----------------------------------------------------------------------

template<class M>
class timeslice_array_matrix : public contour_base {

  typedef ppsc::mam::array_matrix<M> ret_type;
  typedef ppsc::mam::array_matrix<M> les_type;
  typedef ppsc::mam::array_matrix<M> tv_type;
  typedef ppsc::mam::array_matrix<M> mat_type;
  
public:

  timeslice_array_matrix(size_type nt, size_type ntau, size_type nflavour, int sig) :
    contour_base(nt, ntau, nflavour, sig),
    ret(nt_, nflavour_, nflavour_),
    les(nt_, nflavour_, nflavour_),
    tv(ntau_, nflavour_, nflavour_),
    mat(ntau_, nflavour_, nflavour_)
  {}

public:

  ret_type ret;
  les_type les;
  tv_type tv;
  mat_type mat;

};
  
// -----------------------------------------------------------------------
template<class SLICET>
timeslice_array_matrix<ppsc::mam::static_matrix_type<1, 1> >
get_timeslice_component(size_type i1, size_type i2, SLICET & ts) {

  timeslice_array_matrix<ppsc::mam::static_matrix_type<1, 1> >
    component(ts.nt(), ts.ntau(), 1, ts.sig());

  for( auto t : range(0, ts.nt()) ) {
    component.ret(t)(0,0) = ts.ret(t)(i1, i2);
    component.les(t)(0,0) = ts.les(t)(i1, i2);
  }

  for( auto tau : range(0, ts.ntau()) ) {
    component.mat(tau)(0,0) = ts.mat(tau)(i1, i2);
    component.tv(tau)(0,0) = ts.tv(tau)(i1, i2);
  }

  return component;
}

// -----------------------------------------------------------------------
template<class TS1, class TS2>
bool compare_timeslices(TS1 & a, TS2 & b, double tol = 1e-9) {

  bool ret = compare_arrays(a.ret, b.ret, tol=tol);
  if(!ret) std::cerr << "ERROR: compare_timeslices ret" << std::endl;

  bool les = compare_arrays(a.les, b.les, tol=tol);
  if(!les) std::cerr << "ERROR: compare_timeslices les" << std::endl;

  bool tv  = compare_arrays(a.tv,  b.tv,  tol=tol);
  if(!tv)  std::cerr << "ERROR: compare_timeslices tv"  << std::endl;

  bool mat = compare_arrays(a.mat, b.mat, tol=tol);
  if(!mat) std::cerr << "ERROR: compare_timeslices mat" << std::endl;
  
  return ret && les && tv && mat;
}
  
// -----------------------------------------------------------------------
} // end namespace cntr
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _CONTOUR_FUNCTIONS_HPP
