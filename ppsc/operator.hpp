#ifndef _PPSC_OPERATOR_H
#define _PPSC_OPERATOR_H

#include "ppsc/ppsc.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace operators {
// -----------------------------------------------------------------------

  using ppsc::mam::dynamic_matrix_type;
  
// -----------------------------------------------------------------------
// Blockmatrix (MATRIX works with Eigen library)
// ** maps a sector s to precisely one sector to_sector_[s]
// -----------------------------------------------------------------------
template <class MATRIX> class Operator {
public:
  int ns_;
  std::vector<int> ssdim_;
  std::vector<int> to_sector_;
  std::vector<MATRIX> M_;

  // ---------------------------------------------------------------------
  void InitDiagZero(std::vector<int> &ssdim) {
    int s, dim;
    ns_ = ssdim.size();
    assert(ns_ >= 0);
    ssdim_ = ssdim;
    to_sector_.resize(ns_);
    M_.resize(ns_);
    for (s = 0; s < ns_; s++) {
      dim = ssdim_[s];
      assert(dim >= 0);
      to_sector_[s] = s;
      M_[s].resize(dim, dim);
      M_[s].setZero(); // works with Eigen
    }
  }

  // ---------------------------------------------------------------------
  void InitZero(std::vector<int> &ssdim, std::vector<int> &to_sector) {
    int s, dim1, dim2, s1;
    ns_ = ssdim.size();
    assert(ns_ >= 0);
    assert(ns_ == to_sector.size());
    ssdim_ = ssdim;
    to_sector_ = to_sector;
    M_.resize(ns_);
    for (s = 0; s < ns_; s++) {
      s1 = to_sector_[s];
      assert(-1 <= s1 && s1 <= ns_ - 1);
      if (s1 == -1) {
        M_[s].resize(0, 0);
      } else {
        dim1 = ssdim_[s];
        dim2 = ssdim_[s1];
        assert(dim1 >= 0 && dim2 >= 0);
        M_[s].resize(dim2, dim1);
        M_[s].setZero(); // works with Eigen
      }
    }
  }

  // ---------------------------------------------------------------------
  void SetZero(void) {
    for (int s = 0; s < ns_; s++)
      if (to_sector_[s] != -1)
        M_[s].setZero();
  }

  // ---------------------------------------------------------------------
  // trace A*B
  template <class SCALAR> SCALAR trace(void) {
    int s, s1;
    SCALAR res = 0.0;
    // with resize of C
    for (s = 0; s < ns_; s++) {
      s1 = to_sector_[s];
      if (s1 == s)
        res += M_[s].trace();
    }
    return res;
  }

  // ---------------------------------------------------------------------
  Operator<MATRIX> & operator*=(const double scale) {
    for(int s=0; s < ns_; s++) M_[s] *= scale;
    return *this;
  }

  // ---------------------------------------------------------------------
  Operator<MATRIX> & operator*=(const std::complex<double> scale) {
    for(int s=0; s < ns_; s++) M_[s] *= scale;
    return *this;
  }
  
  // ---------------------------------------------------------------------
  Operator<MATRIX> operator*(const double scale) {
    Operator<MATRIX> prod = *this;
    prod *= scale;
    return prod;
  }

  // ---------------------------------------------------------------------
  Operator<MATRIX> operator*(const std::complex<double> scale) {
    Operator<MATRIX> prod = *this;
    prod *= scale;
    return prod;
  }
  
  // ---------------------------------------------------------------------
  template<class HILB>
  static Operator<MATRIX> Identity(HILB & hil_) {
    Operator<MATRIX> I;
    I.InitDiagZero(hil_.ssdim_);
    for(int idx = 0; idx < I.M_.size(); idx++) {
      MATRIX m = I.M_[idx];
      I.M_[idx] = MATRIX::Identity(m.rows(), m.cols());
    }
    return I;
  }

  // ---------------------------------------------------------------------
  int total_size() {
    int dim = 0;
    for( auto sector_dim : ssdim_ ) dim += sector_dim;
    return dim;
  }

  // ---------------------------------------------------------------------
  dynamic_matrix_type to_dense() {

    // -- Allocate big dense matrix
    int tsize = this->total_size();
    dynamic_matrix_type mat(dynamic_matrix_type::Zero(tsize, tsize));

    std::vector<int> sector_start(ns_);

    // -- Calculate indices for the start of each sector
    int dim = 0;
    for( auto sidx : range(0, ns_) ) {
      sector_start[sidx] = dim;
      int sector_dim = ssdim_[sidx];
      dim += sector_dim;      
    }

    // -- Fill in the blocks of ther operator in the full dense matrix
    for( auto s1 : range(0, ns_) ) {
      int s2 = to_sector_[s1];
      
      if(s2 == -1) continue;

      // -- Block start indices (row, col)
      int s1_start = sector_start[s1];
      int s2_start = sector_start[s2];

      // -- Block sizes (row, col)
      int n1 = ssdim_[s1];
      int n2 = ssdim_[s2];

      mat.block(s2_start, s1_start, n2, n1) = M_[s1];
    }
    
    return mat;
  }
};

// -----------------------------------------------------------------------
// "Print function" for Operator<MATRIX> type
template <class MATRIX>
std::ostream &operator<<(std::ostream &os, const Operator<MATRIX> &obj) {

  os << "sectors = " << obj.M_.size() << std::endl;

  for (size_t s = 0; s < obj.M_.size(); s++) {
    os << "s1, s2 = " << s << ", " << obj.to_sector_[s] << std::endl;
    os << obj.M_[s] << std::endl;
  }
  return os;
}
// -----------------------------------------------------------------------

typedef Operator<mam::dynamic_matrix_type> cOp;

// -----------------------------------------------------------------------
// some inefficient algebra with the operators, used to construct Hamiltonian
// etc. ...
// C ->  A*C
template <class MATRIX>
void left_mult(Operator<MATRIX> &C, Operator<MATRIX> &A) {
  int ns = A.ns_, s, s1, s2;
  // with resize of C
  assert(C.ns_ == ns);
  for (s = 0; s < C.ns_; s++) {
    assert(C.ssdim_[s] == A.ssdim_[s]);
    // std::cout << "sector: " << s << std::endl;
    s1 = C.to_sector_[s];
    if (s1 == -1) {
      s2 = -1;
    } else {
      s2 = A.to_sector_[s1];
    }
    if (s2 == -1) {
      // war wohl nix!
      C.to_sector_[s] = -1;
      C.M_[s].resize(0, 0);
    } else {
      MATRIX temp;
      C.to_sector_[s] = s2;
      temp = A.M_[s1] * C.M_[s];
      C.M_[s] = temp;
    }
  }
}

// -----------------------------------------------------------------------
template <class MATRIX>
Operator<MATRIX> operator*(const double scale, const Operator<MATRIX> & A) {
  Operator<MATRIX> prod = A;
  prod *= scale;
  return prod;
}

template <class MATRIX>
Operator<MATRIX> operator*(const std::complex<double> scale,
			   const Operator<MATRIX> & A) {
  Operator<MATRIX> prod = A;
  prod *= scale;
  return prod;
}

// -----------------------------------------------------------------------
// Inplace right side multiplication
// C = A * B
  
template <class MATRIX>
Operator<MATRIX> operator*(const Operator<MATRIX> & A, const Operator<MATRIX> & B) {

  assert(A.ns_ == B.ns_);

  Operator<MATRIX> C;

  C.ns_ = A.ns_;
  C.ssdim_ = A.ssdim_;
  
  for(int sb_col = 0; sb_col < B.ns_; sb_col++) {

    int sb_row = B.to_sector_[sb_col];
    int sa_col, sa_row;
    
    if(sb_row != -1) {
      sa_col = sb_row;
      sa_row = A.to_sector_[sa_col];
    } else {
      C.to_sector_.push_back(-1);
      C.M_.push_back(MATRIX::Zero(0,0));
      continue;
    }
    
    C.to_sector_.push_back(sa_row);
    
    if(sa_row != -1) {
      C.M_.push_back( A.M_[sa_col] * B.M_[sb_col] );
    } else {
      C.M_.push_back(MATRIX::Zero(0,0));
    }
  }

  return C;
}

// -----------------------------------------------------------------------
template <class MATRIX>
Operator<MATRIX> operator+(const Operator<MATRIX> & A, const Operator<MATRIX> & B) {

  assert(A.ns_ == B.ns_);

  Operator<MATRIX> C;

  C.ns_ = A.ns_;
  C.ssdim_ = A.ssdim_;

  for(int s = 0; s < B.ns_; s++) {
    assert( A.to_sector_[s] == B.to_sector_[s] || A.to_sector_[s] == -1 || B.to_sector_[s] == -1);

    // NB! We need to treat the cases where one of the matrices is zero in the given sector.
    
    if(A.to_sector_[s] == -1 && B.to_sector_[s] == -1) {
      // A is zero, B is zero
      C.to_sector_.push_back(-1);
      C.M_.push_back( MATRIX::Zero(0,0) );

    } else if (A.to_sector_[s] != -1 && B.to_sector_[s] == -1) {
      // A has contrib, B is zero
      C.to_sector_.push_back( A.to_sector_[s] );
      C.M_.push_back( A.M_[s] );
      
    } else if (A.to_sector_[s] == -1 && B.to_sector_[s] != -1) {      
      // A is zero, B has contrib
      C.to_sector_.push_back( B.to_sector_[s] );
      C.M_.push_back( B.M_[s] );

    } else if (A.to_sector_[s] == B.to_sector_[s]) {      
      // A has contrib, B has contrib, add both together!
      C.to_sector_.push_back( B.to_sector_[s] );
      C.M_.push_back( A.M_[s] + B.M_[s] );

    } else {
      std::cerr << "--> ppsc::operators::operator+: ERROR in operator addition." << std::endl;
      exit(0);
    }
  }
  return C;
}

// -----------------------------------------------------------------------
template <class MATRIX>
Operator<MATRIX> operator-(const Operator<MATRIX> & A, const Operator<MATRIX> & B) {
  Operator<MATRIX> B_minus = B;
  B_minus *= -1;
  return A + B_minus;
}

// -----------------------------------------------------------------------
// C ->  C + alpha*A ...
// fails if C is not diagonal,
// if sector of A is pointing to -1, it is assumed to be 0
template <class MATRIX, class SCALAR>
void incr_force_diagonal(Operator<MATRIX> &C, Operator<MATRIX> &A,
                         SCALAR alpha) {
  //                         SCALAR &alpha) {
  int ns = A.ns_, s, s1;
  // with resize of C
  assert(C.ns_ == ns);
  for (s = 0; s < ns; s++) {
    assert(C.ssdim_[s] == A.ssdim_[s]);
    if (C.to_sector_[s] == -1) {
      // make this diagonal
      C.M_[s].resize(C.ssdim_[s], C.ssdim_[s]);
      C.M_[s].setZero();
      C.to_sector_[s] = s;
    }
    if (C.to_sector_[s] != s) {
      std::cerr << "--> incr_force_diagonal: "
                << "\n C not diagonal" << std::endl;
      exit(0);
    }
    s1 = A.to_sector_[s];
    if (s1 == -1) {
      // Assume sector is just 0, and ignore
    } else {
      if (s1 != s) {
        std::cerr << "--> incr_force_diagonal: "
                  << "\n A not diagonal" << std::endl;
        exit(0);
      }
      C.M_[s] += alpha * A.M_[s];
    }
  }
}

// -----------------------------------------------------------------------
// -- trace A*B
template <class MATRIX, class SCALAR>
SCALAR trace_AB(Operator<MATRIX> &A, Operator<MATRIX> &B) {
  int ns = A.ns_, s, s1, s2;
  SCALAR res = 0.0;
  // with resize of C
  assert(B.ns_ == ns);
  for (s = 0; s < ns; s++) {
    s1 = B.to_sector_[s];
    if (s1 != -1) {
      assert(A.ssdim_[s1] == B.ssdim_[s1]);
      s2 = A.to_sector_[s1];
      if (s2 == s) {
        // do not need to multiply matrices for this, but anyway
        res += (A.M_[s1] * B.M_[s]).trace();
      }
    }
  }
  return res;
}

// -----------------------------------------------------------------------
// Block-vector (VECTOR works with Eigen library)
template <class VECTOR> class Vector {
public:
  int ns_;
  std::vector<int> ssdim_;
  std::vector<VECTOR> M_;
  void InitZero(std::vector<int> &ssdim) {
    int s, dim;
    ns_ = ssdim.size();
    assert(ns_ >= 0);
    ssdim_ = ssdim;
    M_.resize(ns_);
    for (s = 0; s < ns_; s++) {
      dim = ssdim_[s];
      assert(dim >= 0);
      M_[s].resize(dim);
      M_[s].setZero(); // works with Eigen
    }
  }
  void SetZero(void) {
    for (int s = 0; s < ns_; s++)
      M_[s].setZero();
  }
  void printf(void) {
    for (int s = 0; s < ns_; s++) {
      std::printf("sector %d (dim=%d)\n", s, ssdim_[s]);
      for (int i = 0; i < ssdim_[s]; i++) {
        std::printf("%.10g %.10g , ", M_[s](i).real(), M_[s](i).imag());
      }
      std::printf("\n");
    }
    std::printf("...........\n");
  }
};

// -----------------------------------------------------------------------

typedef Vector<cdvector> cVec;
typedef Vector<dvector> dVec;

// -----------------------------------------------------------------------

// V = alpha*V + beta*Operator*V1 ;  V=V1 not allowed
// do not check dimensions
template <class SCALAR, class VECTOR, class MATRIX>
void OperatorMatrixProduct(SCALAR &alpha, Vector<VECTOR> &V, SCALAR &beta,
                           Operator<MATRIX> &Op, Vector<VECTOR> &V1) {
  int s, s1, ns = V.ns_;
  assert(V1.ns_ <= ns); // ??
  assert(Op.ns_ <= ns); // ??

  for (s = 0; s < ns; s++)
    V.M_[s] *= alpha;
  for (s = 0; s < ns; s++) {
    s1 = Op.to_sector_[s];
    // std::cout << "Op : " << s << " --> " << s1 << std::endl;
    if (s1 != -1) {
      V.M_[s1] += (beta * Op.M_[s] * V1.M_[s]);
      // if(s==0) std::cout << s << " ..... : \n" << Op.M_[s]*V1.M_[s] << " --->
      // s1= " << s1 << ": \n" << V.M_[s1] << std::endl;
    }
  }
}

// -----------------------------------------------------------------------
// alpha = Vdagger * M * V

template <class SCALAR, class VECTOR, class MATRIX>
SCALAR OperatorExpectationValue(Vector<VECTOR> &V, Operator<MATRIX> &Op) {
  int s, s1, ns = V.ns_;
  SCALAR alpha, a1;
  assert(Op.ns_ == ns);
  alpha = 0.0;
  for (s = 0; s < ns; s++) {
    s1 = Op.to_sector_[s];
    if (s1 != -1) {
      cdvector tmp = Op.M_[s] * V.M_[s];
      a1 = V.M_[s1].dot(tmp);
      alpha += a1;
    }
  }
  return alpha;
}

// -----------------------------------------------------------------------
template <class SCALAR, class VECTOR> SCALAR Norm2(Vector<VECTOR> &V) {
  int s, ns = V.ns_;
  SCALAR alpha = 0.0, a1;
  for (s = 0; s < ns; s++) {
    a1 = V.M_[s].dot(V.M_[s]);
    alpha += a1;
  }
  return alpha;
}

// -----------------------------------------------------------------------
} // namespace operators
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OPERATOR_H
