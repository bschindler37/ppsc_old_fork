
#ifndef _PPSC_OCA_SIGMA_DIAG_HPP
#define _PPSC_OCA_SIGMA_DIAG_HPP

// -----------------------------------------------------------------------
//
// OCA pseudo particle self energy diagram descriptor
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "sigma.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

using namespace ppsc;
  
// -----------------------------------------------------------------------
// -- OCA pseudo-particle self-energy diagram (with matrix operators)
//
// Calculates the double integral
//
// \Sigma_{OCA}(t, t') = \int_{t'}^t dt_1 \int_{t'}^{t_1} dt_2
//   [ G_a(t, t_1) G_b(t_1, t_2) G_c(t_2, t') ]
//     x Lambda_1(t, t_2) Lambda_2(t_1, t')
//
// where G_a, G_b, G_c are _scalar_ pseudo-particle Green's functions
// and Lambda_1 and Lambda_2 are (hardcoded? fermionic) _scalar_
// hybridization functions
//
// The integral is performed for a given self-energy time-step
// \Sigma(t,t'), where t' is fixed at the current time 'n'.
// (i.e. tmax on the contour), while t takes all possible values.
//
//   0 ---------->--------------╮
//   0 ╭------t----<------------╯t' = tmax (n)
//     |      
//     v
//     |
//  -i\beta
//
// -- Diagram structure
//
//                    Ga                                L1
// G_a(t, t') : t ====<==== t', Lambda_1(t, t') : t ----<---- t'
//
// \Sigma_{OCA}(t, t') =
//
//                        L1(t1, t')
//        L2(t, t2)  ----------<----------
//     ---------<---/---------            \
//    /            /          \            \  
//   |            |            |            | 
//  [o1]====<====[o2]====<====[o3]====<====[o4]
//   t     Ga     t1    Gb     t2    Gc     t'
//
// -----------------------------------------------------------------------
  
// -----------------------------------------------------------------------
template<int SDIM=-1, int ADIM=-1, int BDIM=-1, int CDIM=-1>
class sdiagram_configuration {

public:

  // -- Matrix types
  
  //using lambda_matrix_type = mam::static_matrix_type<1, 1>;
  using lambda_matrix_type = mam::dynamic_matrix_type;
  using dynamic_matrix_type = mam::dynamic_matrix_type;

  using sigma_matrix_type = mam::static_matrix_type<SDIM, SDIM>;

  using ga_matrix_type = mam::static_matrix_type<ADIM, ADIM>;
  using gb_matrix_type = mam::static_matrix_type<BDIM, BDIM>;
  using gc_matrix_type = mam::static_matrix_type<CDIM, CDIM>;

  using ab_matrix_type = mam::static_matrix_type<ADIM, BDIM>;
  using ac_matrix_type = mam::static_matrix_type<ADIM, CDIM>;
  using bc_matrix_type = mam::static_matrix_type<BDIM, CDIM>;
  
  using op1_matrix_type = mam::static_matrix_type<SDIM, ADIM>;
  using op2_matrix_type = mam::static_matrix_type<ADIM, BDIM>;
  using op3_matrix_type = mam::static_matrix_type<BDIM, CDIM>;
  using op4_matrix_type = mam::static_matrix_type<CDIM, SDIM>;

  // -- Contour compound types
  
  using sigma_type = ppsc::cntr::timeslice_array_matrix_ref<sigma_matrix_type>;
  
  using ga_type = ppsc::cntr::herm_matrix_matrix_ref<ga_matrix_type>;
  using gb_type = ppsc::cntr::herm_matrix_matrix_ref<gb_matrix_type>;
  using gc_type = ppsc::cntr::herm_matrix_matrix_ref<gc_matrix_type>;

  using lambda_type = ppsc::cntr::herm_matrix_matrix_ref<lambda_matrix_type>;

  // -- Array matrix types

  using ga_array_type = mam::array_matrix_ref<ga_matrix_type>;
  using gb_array_type = mam::array_matrix_ref<gb_matrix_type>;
  using gc_array_type = mam::array_matrix_ref<gc_matrix_type>;

  using lambda_array_type = mam::array_matrix_ref<lambda_matrix_type>;
  
  sdiagram_configuration(gf_tstp_type & s_in, int sig,
			 dynamic_matrix_type op1_in, dynamic_matrix_type op2_in,
			 dynamic_matrix_type op3_in, dynamic_matrix_type op4_in,
			 ppgf_type &ga_in, ppgf_type &gb_in, ppgf_type &gc_in,
                         gf_type &lam1_in, gf_type &lam2_in,
			 int sigma_sector, value_type prefactor=1.0) :
    s(s_in), sig(sig), op1(op1_in), op2(op2_in), op3(op3_in), op4(op4_in),
    ga(ga_in), gb(gb_in), gc(gc_in), lam1(lam1_in), lam2(lam2_in),
    sigma_sector(sigma_sector), prefactor(prefactor),
    sdim(s.nflavour()), adim(ga.nflavour()), bdim(gb.nflavour()), cdim(gc.nflavour())
  {static_matrix_asserts();}

  sdiagram_configuration(gf_tstp_type & s_in, int sig,
			 dynamic_matrix_type op1_in, dynamic_matrix_type op2_in,
			 dynamic_matrix_type op3_in, dynamic_matrix_type op4_in,
			 ppgf_type &ga_in, ppgf_type &gb_in, ppgf_type &gc_in,
			 lambda_type &lam1_in, lambda_type &lam2_in, // crtp data types for hybridiz
			 int sigma_sector, value_type prefactor=1.0) : 
    s(s_in), sig(sig), op1(op1_in), op2(op2_in), op3(op3_in), op4(op4_in),
    ga(ga_in), gb(gb_in), gc(gc_in), lam1(lam1_in), lam2(lam2_in),
    sigma_sector(sigma_sector), prefactor(prefactor),
    sdim(s.nflavour()), adim(ga.nflavour()), bdim(gb.nflavour()), cdim(gc.nflavour())
  {static_matrix_asserts();}

  sdiagram_configuration(const sdiagram_configuration<-1, -1, -1, -1> & d) :
    s(d.s), sig(d.sig), op1(d.op1), op2(d.op2), op3(d.op3), op4(d.op4),
    ga(d.ga), gb(d.gb), gc(d.gc), lam1(d.lam1), lam2(d.lam2),
    sigma_sector(d.sigma_sector), prefactor(d.prefactor),
    sdim(s.nflavour()), adim(ga.nflavour()), bdim(gb.nflavour()), cdim(gc.nflavour())
  {static_matrix_asserts();}

  int out_idx() { return sigma_sector; }
  int out_size() { return s.nflavour(); }

  void static_matrix_asserts();
  
  sigma_type s;

  int sig;
  
  ga_type ga;
  gb_type gb;
  gc_type gc;

  lambda_type lam1, lam2;

  op1_matrix_type op1;
  op2_matrix_type op2;
  op3_matrix_type op3;
  op4_matrix_type op4;

  int sigma_sector;
  value_type prefactor;

  const int sdim;
  const int adim;
  const int bdim;
  const int cdim;
  const int ldim = 1;
  
};

typedef std::vector<sdiagram_configuration<> > sdiagrams_type;

// -----------------------------------------------------------------------
template<int SDIM, int ADIM, int BDIM, int CDIM>
std::ostream &operator<<(std::ostream & os,
			 const sdiagram_configuration<SDIM, ADIM, BDIM, CDIM> & d) {

  typedef sdiagram_configuration<SDIM, ADIM, BDIM, CDIM> diagram_type;
  
  os << "--> oca::sdiagram_configuration: Information" << std::endl;
  os << "sdim, adim, bdim, cdim = "
     << d.sdim << ", " << d.adim << ", " << d.bdim << ", " << d.cdim << std::endl;

  os << "--> Pseudo particle matrix types" << std::endl;

  os << std::endl;  

  os << "d.s: nt, ntau, nflavour = " << d.s.nt() << ", " << d.s.ntau() << ", " << d.s.nflavour() << std::endl;
  os << "sigma_matrix_type::RowsAtCompileTime = " << diagram_type::sigma_matrix_type::RowsAtCompileTime << std::endl;
  os << "sigma_matrix_type::ColsAtCompileTime = " << diagram_type::sigma_matrix_type::ColsAtCompileTime << std::endl;
  
  os << std::endl;

  os << "d.ga: nt, ntau, nflavour = " << d.ga.nt() << ", " << d.ga.ntau() << ", " << d.ga.nflavour() << std::endl;
  os << "ga_matrix_type::RowsAtCompileTime = " << diagram_type::ga_matrix_type::RowsAtCompileTime << std::endl;
  os << "ga_matrix_type::ColsAtCompileTime = " << diagram_type::ga_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;
  
  os << "d.gb: nt, ntau, nflavour = " << d.gb.nt() << ", " << d.gb.ntau() << ", " << d.gb.nflavour() << std::endl;
  os << "gb_matrix_type::RowsAtCompileTime = " << diagram_type::gb_matrix_type::RowsAtCompileTime << std::endl;
  os << "gb_matrix_type::ColsAtCompileTime = " << diagram_type::gb_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;

  os << "d.gc: nt, ntau, nflavour = " << d.gc.nt() << ", " << d.gc.ntau() << ", " << d.gc.nflavour() << std::endl;
  os << "gc_matrix_type::RowsAtCompileTime = " << diagram_type::gc_matrix_type::RowsAtCompileTime << std::endl;
  os << "gc_matrix_type::ColsAtCompileTime = " << diagram_type::gc_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;

  os << "--> Product matrix types" << std::endl;  

  os << std::endl;  

  os << "ab_matrix_type::RowsAtCompileTime = " << diagram_type::ab_matrix_type::RowsAtCompileTime << std::endl;
  os << "ab_matrix_type::ColsAtCompileTime = " << diagram_type::ab_matrix_type::ColsAtCompileTime << std::endl;
  os << std::endl;
  
  os << "ac_matrix_type::ColsAtCompileTime = " << diagram_type::ac_matrix_type::ColsAtCompileTime << std::endl;
  os << "ac_matrix_type::ColsAtCompileTime = " << diagram_type::ac_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;

  os << "bc_matrix_type::ColsAtCompileTime = " << diagram_type::bc_matrix_type::ColsAtCompileTime << std::endl;
  os << "bc_matrix_type::ColsAtCompileTime = " << diagram_type::bc_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;  

  os << "--> Operator matrix types" << std::endl;  

  os << std::endl;

  os << "op1: rows, cols = " << d.op1.rows() << ", " << d.op1.cols() << std::endl;
  os << "op1_matrix_type::RowsAtCompileTime = " << diagram_type::op1_matrix_type::RowsAtCompileTime << std::endl;
  os << "op1_matrix_type::ColsAtCompileTime = " << diagram_type::op1_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;

  os << "op2: rows, cols = " << d.op2.rows() << ", " << d.op2.cols() << std::endl;
  os << "op2_matrix_type::RowsAtCompileTime = " << diagram_type::op2_matrix_type::RowsAtCompileTime << std::endl;
  os << "op2_matrix_type::ColsAtCompileTime = " << diagram_type::op2_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;

  os << "op3: rows, cols = " << d.op3.rows() << ", " << d.op3.cols() << std::endl;
  os << "op3_matrix_type::RowsAtCompileTime = " << diagram_type::op3_matrix_type::RowsAtCompileTime << std::endl;
  os << "op3_matrix_type::ColsAtCompileTime = " << diagram_type::op3_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;

  os << "op4: rows, cols = " << d.op4.rows() << ", " << d.op4.cols() << std::endl;
  os << "op4_matrix_type::RowsAtCompileTime = " << diagram_type::op4_matrix_type::RowsAtCompileTime << std::endl;
  os << "op4_matrix_type::ColsAtCompileTime = " << diagram_type::op4_matrix_type::ColsAtCompileTime << std::endl;

  os << std::endl;
  
  return os;
}

// -----------------------------------------------------------------------
template<int SDIM, int ADIM, int BDIM, int CDIM>
void sdiagram_configuration<SDIM, ADIM, BDIM, CDIM>::static_matrix_asserts() {

  static_assert( sigma_matrix_type::RowsAtCompileTime == sigma_matrix_type::ColsAtCompileTime, "" );
  static_assert( ga_matrix_type::RowsAtCompileTime == ga_matrix_type::ColsAtCompileTime, "" );
  static_assert( gb_matrix_type::RowsAtCompileTime == gb_matrix_type::ColsAtCompileTime, "" );
  static_assert( gc_matrix_type::RowsAtCompileTime == gc_matrix_type::ColsAtCompileTime, "" );
  
  static_assert( ab_matrix_type::RowsAtCompileTime == ga_matrix_type::ColsAtCompileTime, "" );
  static_assert( ab_matrix_type::ColsAtCompileTime == gb_matrix_type::RowsAtCompileTime, "" );

  static_assert( ac_matrix_type::RowsAtCompileTime == ga_matrix_type::ColsAtCompileTime, "" );
  static_assert( ac_matrix_type::ColsAtCompileTime == gc_matrix_type::RowsAtCompileTime, "" );

  static_assert( bc_matrix_type::RowsAtCompileTime == gb_matrix_type::ColsAtCompileTime, "" );
  static_assert( bc_matrix_type::ColsAtCompileTime == gc_matrix_type::RowsAtCompileTime, "" );
  
  static_assert( op1_matrix_type::RowsAtCompileTime == sigma_matrix_type::RowsAtCompileTime, "" );
  static_assert( op1_matrix_type::ColsAtCompileTime == ga_matrix_type::RowsAtCompileTime, "" );

  static_assert( op2_matrix_type::RowsAtCompileTime == ga_matrix_type::RowsAtCompileTime, "" );
  static_assert( op2_matrix_type::ColsAtCompileTime == gb_matrix_type::RowsAtCompileTime, "" );

  static_assert( op3_matrix_type::RowsAtCompileTime == gb_matrix_type::RowsAtCompileTime, "" );
  static_assert( op3_matrix_type::ColsAtCompileTime == gc_matrix_type::RowsAtCompileTime, "" );

  static_assert( op4_matrix_type::RowsAtCompileTime == gc_matrix_type::RowsAtCompileTime, "" );
  static_assert( op4_matrix_type::ColsAtCompileTime == sigma_matrix_type::RowsAtCompileTime, "" );
}

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_DIAG_HPP
