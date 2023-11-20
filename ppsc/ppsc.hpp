
#ifndef _PPSC_HPP
#define _PPSC_HPP

// -----------------------------------------------------------------------
//
// Pseudo particle strong coupling expansion utilities
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/cntr_types.hpp"
#include "ppsc/data_types.hpp"
#include "ppsc/operator.hpp"
#include "ppsc/serialization.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

using mam::range;
using mam::value_type;

typedef ppsc::operators::Operator<mam::dynamic_matrix_type> operator_type;
typedef std::vector<operator_type> operators_type;

// -----------------------------------------------------------------------
gf_tstps_type get_list_of_timesteps(int nt, int ntau, std::vector<int> size_vec);
ppgfs_type get_list_of_ppgfs(int nt, int ntau, std::vector<int> size_vec, std::vector<int> sig_vec);
functions_type get_list_of_functions(int nt, std::vector<int> size_vec);
operators_type get_list_of_operators(int nt, std::vector<int> size_vec);
  
void set_bwd_from_fwd(int t, gf_type & gf_bwd_out, gf_type & gf_fwd_in);
void set_bwd_from_fwd_single_component(int t, gf_type & gf_bwd_out, gf_type & gf_fwd_in, gf_type & gf_fwd_in_orbital_transpose);

void solve_pp_dyson(int tstp, ppgfs_type & ppGfs, ppgfs_type & ppSigmas,
		    functions_type H, double pp_mu, double beta, double h,
		    int kt, int nomp);

void set_t0_from_mat(ppgfs_type & ppGfs);

mam::dynamic_matrix_type density_matrix(int tstp, ppgf_type & ppGf_in);
void update_density_matrix(int tstp, ppgfs_type & ppGfs, operators_type & rho);

void normalize_ppgfs_mat(ppgfs_type & ppGfs, double & pp_mu, double beta);
  
double ppgfs_density(int tstp, ppgfs_type & ppGfs);

value_type expectation_value(int tstp, operators_type & rho, operator_type & op);
void extrapolate_ppgf_timestep(int tstp, ppgfs_type & ppGfs, int kt);
void density_matrix_eigenvalues(int tstp, operators_type & rho);

void zeroth_order(ppgfs_type & ppGfs, functions_type & H, double & pp_mu, double beta);
std::vector<double> get_kinetic_energy(gf_type & Gloc, gf_type & Delta, double beta, double h, int kt);
void force_matsubara_hermitian(ppgf_type & G_in);

// -----------------------------------------------------------------------
class pseudo_particle_interaction {

public:

  //typedef mam::static_matrix_type<1, 1> lam_matrix_type;
  typedef mam::dynamic_matrix_type lam_matrix_type;
  typedef cntr::herm_matrix_matrix_ref<lam_matrix_type> lam_ref_type;

  // scalar 
  pseudo_particle_interaction(gf_type & lam,
			      operator_type op1,
			      operator_type op2,
			      int sig=+1, // statistics +1 Boson, -1 Fermion
			      int dir=+1  // hybridization function direction fwd +1, bwd -1
			      ) :
    lam(lam), op1(op1), op2(op2), sig(sig), dir(dir), idx1(0), idx2(0) {
    assert( lam.size1() == 1 && lam.size2() == 1 );
  }

  pseudo_particle_interaction(gf_type & lam,
			      int idx1, int idx2,
			      operator_type op1,
			      operator_type op2,
			      int sig=+1, // statistics +1 Boson, -1 Fermion
			      int dir=+1  // hybridization function direction fwd +1, bwd -1
			      ) :
    lam(lam), op1(op1), op2(op2), sig(sig), dir(dir), idx1(idx1), idx2(idx2) {}
  
  lam_ref_type lam;
  operator_type op1, op2;
  int sig, dir;
  int idx1, idx2; // component indices lam(t, t')(idx1, idx2)
  
};

typedef pseudo_particle_interaction pp_int_type;
typedef std::vector<pp_int_type> pp_ints_type;
  
// -----------------------------------------------------------------------
class greens_function_vertex {

public:

  greens_function_vertex(int idx1, int idx2,
			 operator_type op1,
			 operator_type op2) :
    idx1(idx1), idx2(idx2), op1(op1), op2(op2) {}
  
  int idx1, idx2;
  operator_type op1, op2;
  
};

typedef greens_function_vertex gf_vert_type;
typedef std::vector<gf_vert_type> gf_verts_type;
  
// -----------------------------------------------------------------------
class hilbert_space_type {
public:

  template<class HST>
  hilbert_space_type(HST & hs) :
    ssdim_(hs.ssdim_), sig_(hs.sig_), ns_(hs.ns_), nh_(hs.nh_) {}
  //sector_sizes(hs.ssdim_), sector_sig(hs.sig_) {}

  //std::vector<int> sector_sizes;
  //std::vector<int> sector_sig;
  std::vector<int> ssdim_;
  std::vector<int> sig_;
  int ns_;
  unsigned long nh_;
  
};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HPP
