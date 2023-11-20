
#ifndef _PPSC_OCA_INTEGRATION_HPP
#define _PPSC_OCA_INTEGRATION_HPP

// -----------------------------------------------------------------------
//
// Matrix valued contour integrators (used for OCA diagrams)
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/data_types.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
namespace integrals {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// -- auxiliary convolution routines:
// -- rewritten from old code with new integration object
// -- sum_{l=0}^{max(k,n)} v1[l]*w(k,n,l)

template<typename V>
typename V::matrix_type gregory_sum(
  integration::Integrator<double> & integrator, int n, V & vector) {
  //  const oca::integration::Integrator<double> & integrator, int n, V & vector) {

  using matrix_type = typename V::matrix_type;
  auto shape = vector.get_flavour_shape();
  matrix_type out(matrix_type::Zero(shape[0], shape[1]));

  int k = integrator.get_k();
  if(n < k + 1) { // -- Warmup integration
    for(int l=0; l <= k; l++)
      out += integrator.gregory_weights(n, l) * vector(l);
  } else { // -- Full range integration
    int greg_j = (n<=2*k+1 ? n-k-1 : k);
    for(int l=0; l <= k; l++)
      out += integrator.gregory_weights(n, l) * vector(l);
    for(int l=k+1; l < n-greg_j; l++)
      out += vector(l);
    for(int l1=0; l1 <= greg_j; l1++)
      out += integrator.gregory_omega(l1) * vector(n-l1);
  }
  return out;
}

// -----------------------------------------------------------------------
template<typename V>
typename V::matrix_type gregory_sum(int kt, int n, V & vector) {
  return gregory_sum(integration::I<double>(kt), n, vector);
}

// -----------------------------------------------------------------------
template<typename V>
typename V::matrix_type poly_integral(
  V & vector, int i0, int i1,
  integration::Integrator<double> & integrator) {
  //const oca::integration::Integrator<double> & integrator) {

  using mam::get_view_from;
  typedef typename V::matrix_type matrix_type;

  int kt = integrator.get_k();

  if(i1 - i0 > kt) { // -- Straigh sub inteval integration
    mam::array_matrix_ref<matrix_type> view = get_view_from(i0, vector);
    return gregory_sum(integrator, i1 - i0, view);

  } else { // -- Intepolate outside range

    auto shape = vector.get_flavour_shape();
    matrix_type out(matrix_type::Zero(shape[0], shape[1]));

    for(int j = 0; j <= kt; j++)
      out += integrator.poly_integration(i0, i1, j) * vector(j);

    return out;
  }
}

// -----------------------------------------------------------------------
template<typename V>
typename V::matrix_type poly_integral(
  V & vector,int i0, int i1, int kt) {
  return poly_integral(vector, i0, i1, integration::I<double>(kt));
  //  return poly_integral(vector, i0, i1, oca::integration::I<double>(kt));
}

// -----------------------------------------------------------------------
//  convolution_2b
// -----------------------------------------------------------------------
//
// c = int_{jh}^{nh} dt1 dt2 f(t1,t2)
// f_ij=f[i*(dim+1)+j]=f(ih,jh),
// i=0...imax, j=0...min(i+kt,imax), imax=max(kt,n)
//
// imax = max(kt, n)
// jmax = min(i + kt, imax)
// c = \sum_{i=0}^imax \sum_{j=0}^jmax(i) f_ij
//
// Do not understand the inputs
//
// dim : outer dim of f?
// n : current time step
// j : what is this ???
//
// from the integral above it looks as if we do
//
// \int_j^n di \int_?^? dj f(i, j) ?? repeated labels...?
//
// -- Integrate over [j, n] in both indices?
// -- Integrate over t2 first
// -- store in temporary variable
// -- and integrate temporary over t1

template<typename M>
typename M::matrix_type convolution2b(
  M & mat, int j, int n,
  integration::Integrator<double> & integrator)
//  const oca::integration::Integrator<double> & integrator)
{

  using mam::get_column_vector;
  typedef typename M::matrix_type matrix_type;

  int kt = integrator.get_k();

  int i1max = (n>=kt ? n : kt);
  int i1min = i1max-kt;
  if(i1min > j) i1min = j;

  auto shape = mat.get_flavour_shape();
  mam::array_matrix<matrix_type> f1(i1max+1, shape[0], shape[1]);

  for(int i = i1min; i < j; i++) {
    mam::array_matrix_ref<matrix_type> vec_view = get_column_vector(i, mat);
    f1(i) = -1.0 * poly_integral(vec_view, i-i1min, j-i1min, integrator);
  }

  f1(j) = matrix_type::Zero(shape[0], shape[1]);

  for(int i = j + 1; i <= i1max; i++) {
    mam::array_matrix_ref<matrix_type> array_view = get_column_vector(i, mat);
    f1(i) = poly_integral(array_view, j-i1min, i-i1min, integrator);
  }

  return poly_integral(f1, j - i1min, n - i1min, integrator);
}

// -----------------------------------------------------------------------
template<typename M>
typename M::matrix_type convolution2b(M & mat, int j, int n, int kt) {
  return convolution2b(mat, j, n, integration::I<double>(kt));
  //return convolution2b(mat, j, n, oca::integration::I<double>(kt));
}

// -----------------------------------------------------------------------
//
//  matsubara convolution
//
//  c(x0) = int_0^X dx a(j-x)b(x)  for x0=j, X=mmax
//  a is antiperiodic, a(x+X)=-a(x)
//  a[0..mmax]: a[0]=a(0^+), a[mmax]=a(X^-)
//
// -----------------------------------------------------------------------

// c(x0) = int_0^j dx a(j-x)b(x)  for x0=j

// -----------------------------------------------------------------------
template<class A1, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_lowside(
  int j, A1 & aa, A2 & bb, int mmax, integration::Integrator<double> &I) {

  using out_matrix_type = mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime,
						   A2::matrix_type::ColsAtCompileTime>;

  assert( j <= mmax );
  assert( aa.get_array_shape()[0] >= mmax );
  assert( bb.get_array_shape()[0] >= mmax );
  assert( aa.get_flavour_shape()[1] == bb.get_flavour_shape()[0] ) ; // Req for matrix multiplications

  int k = I.get_k();

  auto rows = aa.get_flavour_shape()[0];
  auto cols = bb.get_flavour_shape()[1];
  out_matrix_type out(out_matrix_type::Zero(rows, cols));
  out_matrix_type tmp;

  if(j > 0 && j < k) {
    for( auto l : range(0, k+1) ) {
      for( auto m : range(0, k+1) )
	{
	  tmp = aa(l) * bb(m);
	  out += I.rcorr(j,l,m) * tmp ;
	}
    }
  } else if (j > 0) {
    mam::array_matrix<out_matrix_type> f1(j+1, rows, cols);
    for( auto l : range(0, j+1) ) f1(l) = aa(j-l) * bb(l);
    out = gregory_sum(I, j, f1);
  } else {
    // this is actually happening..
  }
  return out;
}
// -----------------------------------------------------------------------
template<class A1, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_lowside(int j, A1 & aa, A2 & bb, int mmax, int kt) {
  return matsubara_convolution_lowside(j, aa, bb, mmax, integration::I<double>(kt));
}

// -----------------------------------------------------------------------
template<class A1, class M, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_lowside_operator(
  int j, A1 & aa, M & op, A2 & bb, int mmax, integration::Integrator<double> &I) {

  using out_matrix_type = mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime,
						   A2::matrix_type::ColsAtCompileTime>;

  assert( aa.get_array_shape()[0] == bb.get_array_shape()[0] );
  assert( aa.get_flavour_shape()[1] == op.rows() );
  assert( bb.get_flavour_shape()[0] == op.cols() );

  int k = I.get_k();

  auto rows = aa.get_flavour_shape()[0];
  auto cols = bb.get_flavour_shape()[1];

  out_matrix_type out(out_matrix_type::Zero(rows, cols));
  out_matrix_type tmp;
  
  if(j > 0 && j < k) {
    for( auto l : range(0, k+1) ) {
      for( auto m : range(0, k+1) ) {
	tmp = aa(l) * op * bb(m);
	out += I.rcorr(j,l,m) * tmp;
      }
    }
  } else if (j > 0) {
    mam::array_matrix<out_matrix_type> f1(j+1, rows, cols);
    for( auto l : range(0, j+1) ) f1(l) = aa(j-l) * op * bb(l);
    out = gregory_sum(I, j, f1);
  } else {
    // this is actually happening..
  }
  return out;
}
// -----------------------------------------------------------------------
template<class A1, class M, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_lowside_operator(int j, A1 & aa, M & op, A2 & bb, int mmax, int kt) {
  return matsubara_convolution_lowside_operator(j, aa, op, bb, mmax, integration::I<double>(kt));
}

// -----------------------------------------------------------------------
// c = int_j^mmax dx a(x-j)b(x)
template<class A1, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_highside(
  int j, A1 & aa, A2 & bb,int mmax, integration::Integrator<double> & I)
{
  using out_matrix_type = mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime,
						   A2::matrix_type::ColsAtCompileTime>;

  int k = I.get_k();
  int j1 = mmax - j;

  auto rows = aa.get_flavour_shape()[0];
  auto cols = bb.get_flavour_shape()[1];
  out_matrix_type out(out_matrix_type::Zero(rows, cols));
  out_matrix_type tmp;

    if(j1 > 0 && j1 < k) {
      for(int l=0;l<=k;l++){
	for(int m=0;m<=k;m++){
	  tmp = aa(m) * bb(mmax - l);
	  out +=  I.rcorr(j1, l, m) * tmp;
	}
      }
    } else if (j1 > 0) {
      mam::array_matrix<out_matrix_type> f1(j1+1, rows, cols);
      for(int l=0;l<=j1;l++)
	f1(l) = aa(l) * bb(mmax - (j1 - l));
      out = gregory_sum(I, j1, f1);
    }
  return out;
}

// -----------------------------------------------------------------------
template<class A1, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_highside(int j, A1 & aa, A2 & bb,int mmax, int kt) {
    return matsubara_convolution_highside(j, aa, bb, mmax, integration::I<double>(kt));
}

// -----------------------------------------------------------------------
// c = int_j^mmax dx a(x-j)b(x)
template<class A1, class M, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_highside_operator(
  int j, A1 & aa, M & op, A2 & bb,int mmax, integration::Integrator<double> & I)
{
  using out_matrix_type = mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime,
						   A2::matrix_type::ColsAtCompileTime>;

  assert( aa.get_array_shape()[0] == bb.get_array_shape()[0] );
  assert( aa.get_flavour_shape()[1] == op.rows() );
  assert( bb.get_flavour_shape()[0] == op.cols() );

  int k = I.get_k();
  int j1 = mmax - j;

  auto rows = aa.get_flavour_shape()[0];
  auto cols = bb.get_flavour_shape()[1];
  out_matrix_type out(out_matrix_type::Zero(rows, cols));
  out_matrix_type tmp;

    if(j1 > 0 && j1 < k) {
      for(int l=0;l<=k;l++){
	for(int m=0;m<=k;m++){
	  //out +=  I.rcorr(j1, l, m) * aa(m) * op * bb(mmax - l);
	  tmp = aa(mmax - l) * op * bb(m);
	  out += I.rcorr(j1, l, m) * tmp;
	}
      }
    } else if (j1 > 0) {
      mam::array_matrix<out_matrix_type> f1(j1+1, rows, cols);
      for(int l=0;l<=j1;l++)
	//f1(l) = aa(l) * op * bb(mmax - (j1 - l));
	f1(l) = aa(mmax - (j1 - l)) * op * bb(l);
      out = gregory_sum(I, j1, f1);
    }
  return out;
}

// -----------------------------------------------------------------------
template<class A1, class M, class A2>
mam::static_matrix_type<A1::matrix_type::RowsAtCompileTime, A2::matrix_type::ColsAtCompileTime>
  matsubara_convolution_highside_operator(int j, A1 & aa, M & op, A2 & bb,int mmax, int kt) {
  return matsubara_convolution_highside_operator(j, aa, op, bb, mmax, integration::I<double>(kt));
}

// -----------------------------------------------------------------------
// -- interpolate function f at x

template<class ARR>
typename ARR::matrix_type poly_interpol(ARR & f, double x,
					integration::Integrator<double> &I) {

  typedef typename ARR::matrix_type matrix_type;

  int kt = I.get_k();
  double x1 = 1.0;

  auto shape = f.get_flavour_shape();
  matrix_type sum(matrix_type::Zero(shape[0], shape[1]));

  for(int l=0;l<=kt;l++) {
    matrix_type fl(matrix_type::Zero(shape[0], shape[1]));
    for(int i=0;i<=kt;i++) {
      fl += I.poly_interpolation(l, i) * f(i);
    }
    sum += fl * x1;
    x1 *= x;
  }
  return sum;
}

// -----------------------------------------------------------------------
template<class ARR>
typename ARR::matrix_type poly_interpol(ARR & f, double x, int kt) {
  return poly_interpol(f, x, integration::I<double>(kt));
}

// -----------------------------------------------------------------------
} // end namespace integration
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_INTEGRATION_HPP
