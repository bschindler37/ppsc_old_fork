
#ifndef _PPSC_OCA_GF_DIAG_HPP
#define _PPSC_OCA_GF_DIAG_HPP

// -----------------------------------------------------------------------
//
// OCA single-particle Green's function diagram descriptor
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "gf.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

using namespace ppsc;
  
// -----------------------------------------------------------------------
// -- OCA single-particle Green's function (with matrix operators)
//
// Calculates the integral
//
// G_{OCA}(t, t') = \int_{t'}^t dt_1 \int_{t'}^{t_1} dt_2
//   [ G_a(t_1, t) G_b(t', t_1) G_c(t_2, t') G_d(t, t_2) ] x Lambda(t_1, t_2) 
//
// where G_a, G_b, G_c, G_d are _scalar_ psedo-particle Green's functions
// and Lambda is a (hardcoded? fermionic) _scalar hybridization function
//
// The integral is performed for a given Green's function time-step
// G_{OCA}(t,t'), where t' is fixed at the current time 'n'.
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
//  G_{sp}(t, t') =
//                         Lambda(t1, t2)
//                    ---------->------------
//                   /                       \
//            Gd    |      Gc           Gb    |      Ga          
//      []====<====[o1]====<====[o2]====<====[o3]====<====[o4]
//      t'          t2           t            t1           t'
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<int ADIM=-1, int BDIM=-1, int CDIM=-1, int DDIM=-1>
class gdiagram_configuration {

public:

  // -- Matrix types
  
  //using lambda_matrix_type = mam::static_matrix_type<1, 1>;
  using lambda_matrix_type = mam::dynamic_matrix_type;
  using dynamic_matrix_type = mam::dynamic_matrix_type;

  using g_matrix_type = mam::static_matrix_type<1, 1>;

  using ga_matrix_type = mam::static_matrix_type<ADIM, ADIM>;
  using gb_matrix_type = mam::static_matrix_type<BDIM, BDIM>;
  using gc_matrix_type = mam::static_matrix_type<CDIM, CDIM>;
  using gd_matrix_type = mam::static_matrix_type<DDIM, DDIM>;

  using ba_matrix_type = mam::static_matrix_type<BDIM, ADIM>;
  using ca_matrix_type = mam::static_matrix_type<CDIM, ADIM>;
  using da_matrix_type = mam::static_matrix_type<DDIM, ADIM>;

  using cb_matrix_type = mam::static_matrix_type<CDIM, BDIM>;
  using db_matrix_type = mam::static_matrix_type<DDIM, BDIM>;
  
  using dc_matrix_type = mam::static_matrix_type<DDIM, CDIM>;
  
  using op1_matrix_type = mam::static_matrix_type<DDIM, CDIM>;
  using op2_matrix_type = mam::static_matrix_type<CDIM, BDIM>;
  using op3_matrix_type = mam::static_matrix_type<BDIM, ADIM>;
  using op4_matrix_type = mam::static_matrix_type<ADIM, DDIM>;

  // -- Contour compound types
  
  using g_type = ppsc::cntr::timeslice_array_matrix_ref<g_matrix_type>;
  
  using ga_type = ppsc::cntr::herm_matrix_matrix_ref<ga_matrix_type>;
  using gb_type = ppsc::cntr::herm_matrix_matrix_ref<gb_matrix_type>;
  using gc_type = ppsc::cntr::herm_matrix_matrix_ref<gc_matrix_type>;
  using gd_type = ppsc::cntr::herm_matrix_matrix_ref<gd_matrix_type>;

  using lambda_type = ppsc::cntr::herm_matrix_matrix_ref<lambda_matrix_type>;

  // -- Array matrix types

  using ga_array_type = mam::array_matrix_ref<ga_matrix_type>;
  using gb_array_type = mam::array_matrix_ref<gb_matrix_type>;
  using gc_array_type = mam::array_matrix_ref<gc_matrix_type>;
  using gd_array_type = mam::array_matrix_ref<gd_matrix_type>;

  using lambda_array_type = mam::array_matrix_ref<lambda_matrix_type>;
  
  gdiagram_configuration(gf_tstp_type & g_in, int sig,
			 dynamic_matrix_type op1_in, dynamic_matrix_type op2_in,
			 dynamic_matrix_type op3_in, dynamic_matrix_type op4_in,
			 ppgf_type & ga_in, ppgf_type & gb_in, ppgf_type & gc_in, ppgf_type & gd_in,
			 gf_type &lam_in,
			 int vertex_idx, value_type prefactor=1.0) :
    g(g_in), sig(sig), op1(op1_in), op2(op2_in), op3(op3_in), op4(op4_in),
    ga(ga_in), gb(gb_in), gc(gc_in), gd(gd_in), lam(lam_in),
    vertex_idx(vertex_idx), prefactor(prefactor),
    adim(ga.nflavour()), bdim(gb.nflavour()), cdim(gc.nflavour()), ddim(gd.nflavour())
  {}

  gdiagram_configuration(gf_tstp_type & g_in, int sig,
			 dynamic_matrix_type op1_in, dynamic_matrix_type op2_in,
			 dynamic_matrix_type op3_in, dynamic_matrix_type op4_in,
			 ppgf_type & ga_in, ppgf_type & gb_in, ppgf_type & gc_in, ppgf_type & gd_in,
     			 lambda_type &lam_in, // crtp data type for hybridiz
			 int vertex_idx, value_type prefactor=1.0) : 
    g(g_in), sig(sig), op1(op1_in), op2(op2_in), op3(op3_in), op4(op4_in),
    ga(ga_in), gb(gb_in), gc(gc_in), gd(gd_in), lam(lam_in),
    vertex_idx(vertex_idx), prefactor(prefactor),
    adim(ga.nflavour()), bdim(gb.nflavour()), cdim(gc.nflavour()), ddim(gd.nflavour())
  {}

  gdiagram_configuration(const gdiagram_configuration<-1, -1, -1, -1> & d) :
    g(d.g), sig(d.sig), op1(d.op1), op2(d.op2), op3(d.op3), op4(d.op4),
    ga(d.ga), gb(d.gb), gc(d.gc), gd(d.gd), lam(d.lam),
    vertex_idx(d.vertex_idx), prefactor(d.prefactor),
    adim(ga.nflavour()), bdim(gb.nflavour()), cdim(gc.nflavour()), ddim(gd.nflavour())
  {}
  
  int out_idx() { return vertex_idx; }
  int out_size() { return g.nflavour(); }

  g_type g;
  int sig;
  
  op1_matrix_type op1;
  op2_matrix_type op2;
  op3_matrix_type op3;
  op4_matrix_type op4;

  ga_type ga;
  gb_type gb;
  gc_type gc;
  gd_type gd;

  lambda_type lam;

  int vertex_idx;
  value_type prefactor;

  const int gdim = 1;
  const int adim;
  const int bdim;
  const int cdim;
  const int ddim;
  const int ldim = 1;
  
};

typedef std::vector<gdiagram_configuration<> > gdiagrams_type;

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_DIAG_HPP
