// -----------------------------------------------------------------------
//
// Pseudo particle strong coupling expansion solver
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2017)
//
// -----------------------------------------------------------------------
#pragma once

#include "ppsc/solver.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
  namespace lang_firsov {
// -----------------------------------------------------------------------
    
// -----------------------------------------------------------------------
// -- Boson exponent amplitude

// -- Eq. (28) in PRB, 88, 165108 (2013)
// -- Why are the s-factors missing?

std::complex<double> phi(std::complex<double> t, std::complex<double> tprime,
			 double beta, double w0, double g) {

  std::complex<double> I(0., 1.);
  return std::pow(g/w0, 2) / std::sinh(0.5*beta*w0) * (
      std::cosh( w0*(0.5*beta - I*(t - tprime))) - std::cosh(0.5*beta*w0)
      );
}

// -----------------------------------------------------------------------
// -- Construct the bosonic Lang-Firsov kernel function K

void phonon_kernel(ppsc::gf_type &K_in,
		   double w0, double g, double beta, double h) {

  using ppsc::range;
  std::complex<double> I(0., 1.);
  ppsc::cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> K(K_in);

  for( auto m : range(0, K.ntau()) ) {
    double tau = beta * m / (K.ntau() - 1);
    K.mat(m)(0,0) = -std::exp(phi(-I*tau, 0., beta, w0, g));
  }

  for( auto t1 : range(0, K.nt()) ) {
    for( auto t2 : range(0, t1+1) ) {
      std::complex<double> exp_phi_gtr = std::exp(phi(h*t1, h*t2, beta, w0, g));
      std::complex<double> exp_phi_les = std::exp(phi(h*t2, h*t1, beta, w0, g));
      K.ret_lt(t1, t2)(0,0) = -I * (exp_phi_gtr - exp_phi_les);
      K.les_ut(t2, t1)(0,0) = -I * std::exp(phi(h*t1, h*t2, beta, w0, g));
    }
    for( auto m : range(0, K.ntau()) ) {
      double tau = beta * m / (K.ntau() - 1);
      K.tv(t1, m)(0,0) = -I*std::exp(phi(-I*tau, h*t1, beta, w0, g));
    }
  }      
}

// -----------------------------------------------------------------------
void multiply_equal(int tstp, ppsc::gf_type &Delta_in, ppsc::gf_type &K_in) {

  using ppsc::range;
  std::complex<double> I(0., 1.);
  typedef mam::dynamic_matrix_type matrix_type;
  ppsc::cntr::herm_matrix_matrix_ref<matrix_type> K(K_in);
  ppsc::cntr::herm_matrix_matrix_ref<matrix_type> Delta(Delta_in);

  if( tstp == -1 ) {
    for( auto tau : range(0, K.ntau()) ) Delta.mat(tau) *= -K.mat(tau);
  } else { // tstp >= 0
    for( auto t2 : range(0, tstp+1) ) {

      matrix_type DeltaK_les = I * Delta.les(t2, tstp) * K.les(t2, tstp);
      matrix_type DeltaK_ret = I * Delta.gtr(tstp, t2) * K.gtr(tstp, t2)
                             - I * Delta.les(tstp, t2) * K.les(tstp, t2);
      
      Delta.ret_lt(tstp, t2) = DeltaK_ret;
      Delta.les_ut(t2, tstp) = DeltaK_les;
    }
    for( auto tau : range(0, K.ntau()) ) Delta.tv(tstp, tau) *= I * K.tv(tstp, tau);
  }   
}    

// -----------------------------------------------------------------------
void multiply_equal(int tstp, ppgf_type &G_in, ppsc::gf_type &K_in) {

  // Pseudo-particle Green's function bubble with scalar Gf

  // Differs with the normal Green's function bubble in the
  // treatment of the retarded component, that for pseudo particles
  // is given by the greater component only.
  
  using ppsc::range;
  std::complex<double> I(0., 1.);
  typedef mam::dynamic_matrix_type matrix_type;
  ppsc::cntr::herm_matrix_matrix_ref<matrix_type> K(K_in);
  ppsc::cntr::herm_matrix_matrix_ref<matrix_type> G(G_in);

  if( tstp == -1 ) {
    for( auto tau : range(0, K.ntau()) ) G.mat(tau) *= -K.mat(tau)(0,0);
  } else { // tstp >= 0
    for( auto t2 : range(0, tstp+1) ) {
      G.ret_lt(tstp, t2) = I * G.ret(tstp, t2) * K.gtr(tstp, t2)(0,0);
      G.les_ut(t2, tstp) = I * G.les(t2, tstp) * K.les(t2, tstp)(0,0);
    }
    for( auto tau : range(0, K.ntau()) ) G.tv(tstp, tau) *= I * K.tv(tstp, tau)(0,0);
  }   
}

// -----------------------------------------------------------------------
void set_timestep(int tstp, ppgf_type &A_in, ppgf_type &B_in) {

  typedef mam::dynamic_matrix_type matrix_type;
  cntr::herm_matrix_matrix_ref<matrix_type> A(A_in);
  cntr::herm_matrix_matrix_ref<matrix_type> B(B_in);

  if( tstp == -1 ) {
    for( auto tau : range(0, A.ntau()) ) A.mat(tau) = B.mat(tau);
  } else { // tstp >= 0
    for( auto t2 : range(0, tstp+1) ) {
      A.ret_lt(tstp, t2) = B.ret_lt(tstp, t2);
      A.les_ut(t2, tstp) = B.les_ut(t2, tstp);
    }
    for( auto tau : range(0, A.ntau()) ) A.tv(tstp, tau) = B.tv(tstp, tau);
  }   
}    

void set_timestep(int tstp, ppgfs_type &A, ppgfs_type &B) {
  for( auto idx : range(0, A.size()) ) set_timestep(tstp, A[idx], B[idx]);
}
    
// -----------------------------------------------------------------------
void add_timestep(int tstp, ppgf_type &A_in, ppgf_type &B_in) {

  typedef mam::dynamic_matrix_type matrix_type;
  cntr::herm_matrix_matrix_ref<matrix_type> A(A_in);
  cntr::herm_matrix_matrix_ref<matrix_type> B(B_in);

  if( tstp == -1 ) {
    for( auto tau : range(0, A.ntau()) ) A.mat(tau) += B.mat(tau);
  } else { // tstp >= 0
    for( auto t2 : range(0, tstp+1) ) {
      A.ret_lt(tstp, t2) += B.ret(tstp, t2);
      A.les_ut(t2, tstp) += B.les(t2, tstp);
    }
    for( auto tau : range(0, A.ntau()) ) A.tv(tstp, tau) += B.tv(tstp, tau);
  }   
}    
    
void add_timestep(int tstp, ppgfs_type &A, ppgfs_type &B) {
  for( auto idx : range(0, A.size()) ) add_timestep(tstp, A[idx], B[idx]);
}
    
// -----------------------------------------------------------------------
  } // end namespace lang_firsov
// -----------------------------------------------------------------------

template<class HAM>
class solver_lang_firsov : public solver<HAM> {

public:

  typedef solver<HAM> B;
  
  // ---------------------------------------------------------------------

  template<class HILB>
  solver_lang_firsov(int nt, int ntau, double beta, double h, int kt, int nomp,
		     HILB & hilbert_space, int expansion_order=1) :
    solver<HAM>::solver(nt, ntau, beta, h, kt, nomp,
			hilbert_space, expansion_order),
    ppGfs_bare(get_list_of_ppgfs(nt, ntau, hilbert_space.ssdim_, hilbert_space.sig_)),
    ppSigmas_oca(get_list_of_ppgfs(nt, ntau, hilbert_space.ssdim_, hilbert_space.sig_)),
    K(nt, ntau, 1, +1) {}

  // ---------------------------------------------------------------------

  void update_hamiltonian() {
    solver<HAM>::update_hamiltonian();

    lang_firsov::phonon_kernel(
      K, B::hamiltonian.omega, B::hamiltonian.g, B::beta, B::h);
  }

  // ---------------------------------------------------------------------

  void update_sigma(int tstp) {    

    B::nca_sdh.set_pp_self_energy(
      B::ppSigmas, tstp, B::ntau, B::beta, B::h, B::kt, B::nomp);

    if(has_oca()) {

      lang_firsov::set_timestep(tstp, ppGfs_bare, B::ppGfs); // Store the bare ppGfs
      add_phonon_kernel(tstp, B::ppGfs); // Renormalize ppGfs

      B::oca_sdh.set_pp_self_energy(
      	ppSigmas_oca, tstp, B::ntau, B::beta, B::h, B::kt, B::nomp);

      lang_firsov::set_timestep(tstp, B::ppGfs, ppGfs_bare); // Revert ppGfs to bare
      
      add_phonon_kernel(tstp, ppSigmas_oca);
      lang_firsov::add_timestep(tstp, B::ppSigmas, ppSigmas_oca);
    }
    
  }

  // ---------------------------------------------------------------------

  ppsc::gf_tstps_type get_spgf(int tstp) {

    ppsc::gf_tstps_type gf_tstp_vertex_list =
      B::nca_gdh.evaluate_to_out_idx(tstp, B::ntau, B::beta, B::h, B::kt, B::nomp);
    
    if (has_oca()) {

      // Save ppGfs to bare
      // Renormalize ppGfs
      lang_firsov::set_timestep(tstp, ppGfs_bare, B::ppGfs); // Store the bare ppGfs
      add_phonon_kernel(tstp, B::ppGfs); // Renormalize ppGfs

      // Eval OCA spgf
      ppsc::gf_tstps_type oca_gf_tstp_vertex_list =
	B::oca_gdh.evaluate_to_out_idx(tstp, B::ntau, B::beta, B::h, B::kt, B::nomp);

      // Restore ppGfs
      lang_firsov::set_timestep(tstp, B::ppGfs, ppGfs_bare); // Revert ppGfs to bare

      for( auto idx : range(0, gf_tstp_vertex_list.size()) )
	gf_tstp_vertex_list[idx].incr(oca_gf_tstp_vertex_list[idx], 1.0);
    }

    return gf_tstp_vertex_list;
  }

  // ---------------------------------------------------------------------

  void pp_step(int tstp) {
    int n1 = (tstp <= B::kt && tstp >= 0 ? 0 : tstp);
    int n2 = (tstp <= B::kt && tstp >= 0 ? B::kt : tstp);

    for(int n = n1; n <= n2; n++) update_sigma(n);
    B::solve_dyson(tstp);
    if(tstp == -1) B::normalize_ppgf();
    for(int n = n1; n <= n2; n++) B::update_density_matrix(n);
    for(int n = n1; n <= n2; n++) B::hamiltonian.update_exp_vals(n, B::rho);
  }
  
  // ---------------------------------------------------------------------

  void add_phonon_kernel(int tstp, gf_type & Delta) {
    lang_firsov::multiply_equal(tstp, Delta, K);
  }

  void add_phonon_kernel(int tstp, ppgf_type & G) {
    lang_firsov::multiply_equal(tstp, G, K);
  }

  void add_phonon_kernel(int tstp, ppgfs_type & ppGfs) {
    for( auto & ppGf : ppGfs ) add_phonon_kernel(tstp, ppGf);
  }
  
  // ---------------------------------------------------------------------

  void store(hid_t file_id, bool store_pp=false, bool store_K=false) {

    if(store_K){
      hid_t group_id = create_group(file_id, "K");
      store_herm_greens_function(group_id, K);
      close_group(group_id);
    }
    
    B::store(file_id, store_pp);
  }
  
  // ---------------------------------------------------------------------

  ppgfs_type ppGfs_bare;
  ppgfs_type ppSigmas_oca;
  
  gf_type K; // Phonon kernel function
  gf_type K_inverse; // Inverse phonon kernel function
};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------
