
#ifndef _PPSC_CPP
#define _PPSC_CPP

// -----------------------------------------------------------------------
//
// Pseudo particle strong coupling expansion utilities
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

//#include "./dyson/pp_dyson_mat.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
gf_tstps_type get_list_of_timesteps(int nt, int ntau, std::vector<int> size_vec) {
  gf_tstps_type tstps;
  for( auto idx : range(0, size_vec.size()) )
    tstps.push_back(gf_tstp_type(nt, ntau, size_vec[idx]));
  return tstps;
}

// -----------------------------------------------------------------------
ppgfs_type get_list_of_ppgfs(int nt, int ntau,
  std::vector<int> size_vec, std::vector<int> sig_vec) {
  ppgfs_type ppGfs;
  for( auto idx : range(0, size_vec.size()) )
    ppGfs.push_back(ppgf_type(nt, ntau, size_vec[idx], sig_vec[idx]));

  return ppGfs;
}

// -----------------------------------------------------------------------
functions_type get_list_of_functions(int nt, std::vector<int> size_vec) {
  functions_type functions;
  for( auto idx : range(0, size_vec.size()) )
    functions.push_back(function_type(nt, size_vec[idx]));
  return functions;
}

// -----------------------------------------------------------------------
operators_type get_list_of_operators(int nt, std::vector<int> size_vec) {
  operators_type ops(nt+2);
  for( auto idx : range(0, ops.size()) ) ops[idx].InitDiagZero(size_vec);
  return ops;
}
  
// -----------------------------------------------------------------------
void set_bwd_from_fwd_depr(int t, gf_type & gf_bwd_out, gf_type & gf_fwd_in) {

  // Construct response function for backwards propagating case
  //
  // G^M_{bwd}(\tau - 0) = G^M_{fwd}(0 - \tau) = \xi G^M_{fwd}(\beta - \tau)
  // G^R_{bwd}(t, t')    = G^A_{fwd}(t', t)    = - [G^R_{fwd}(t, t')]^+
  // G^<_{bwd}(t, t')    = G^>_{fwd}(t', t)    = G^R_{fwd}(t', t) - [G^<_{fwd}(t, t')]^+
  // G^{tv}_{bwd}(t, tau') = G^{vt}_{fwd}(tau', t) = - [G^{vt}_{fwd}(t, tau')]^+
  //
  // Currently all these relations are multiplied with an extra fermi minus sign
  // This does not work out (so well) for bosonic propagators...
  
  typedef mam::dynamic_matrix_type matrix_type;

  cntr::herm_matrix_matrix_ref<matrix_type> gf_fwd(gf_fwd_in);
  cntr::herm_matrix_matrix_ref<matrix_type> gf_bwd(gf_bwd_out);
   
  int ntau = gf_fwd.ntau();

  assert( gf_fwd.sig() == gf_bwd.sig() );
  int xi = gf_fwd.sig();

  if(t == -1) {
    for( auto tau : range(0, ntau) )
      // additional minus sign assuming fermionic
      //gf_bwd.mat(tau).noalias() = gf_fwd.mat(ntau - 1 - tau);

      // explicit sign
      gf_bwd.mat(tau).noalias() = xi * xi * gf_fwd.mat(ntau - 1 - tau);

  } else {
    
    for( auto j : range(0, t) ) {
      // additional minus sign assuming fermionic
      //gf_bwd.ret_lt(t, j).noalias() = -gf_fwd.ret(t, j).adjoint();
      //gf_bwd.les_ut(j, t).noalias() = -gf_fwd.ret(t, j) + gf_fwd.les(j, t).adjoint();

      // explicit sign
      gf_bwd.ret_lt(t, j).noalias() = -xi * (-gf_fwd.ret(t, j).adjoint());
      gf_bwd.les_ut(j, t).noalias() = -xi * (-gf_fwd.ret(t, j) + gf_fwd.les(j, t).adjoint());

    }

    // additional minus sign assuming fermionic
    //gf_bwd.ret_lt(t, t).noalias() = gf_fwd.ret(t, t);
    //gf_bwd.les_ut(t, t).noalias() = gf_fwd.ret(t, t).adjoint() - gf_fwd.les(t, t);
    
    // explicit sign
    gf_bwd.ret_lt(t, t).noalias() = -xi * gf_fwd.ret(t, t);
    gf_bwd.les_ut(t, t).noalias() = -xi * (gf_fwd.ret(t, t).adjoint() - gf_fwd.les(t, t));

    for( auto tau : range(0, ntau) ) {
      // additional minus sign assuming fermionic
      //gf_bwd.tv(t, ntau - 1 - tau).noalias() = -gf_fwd.tv(t, tau).adjoint();

      // explicit sign
      gf_bwd.tv(t, ntau - 1 - tau).noalias() = (xi * xi) * -gf_fwd.tv(t, tau).adjoint();

    }
  } // t != -1
}

// -----------------------------------------------------------------------
void set_bwd_from_fwd(int t, gf_type & gf_bwd_out, gf_type & gf_fwd_in) {

  // Construct response function for backwards propagating case
  //
  // G^M_{bwd}(\tau - 0) = G^M_{fwd}(0 - \tau) = \xi G^M_{fwd}(\beta - \tau)
  // G^R_{bwd}(t, t')    = G^A_{fwd}(t', t)    = [G^R_{fwd}(t, t')]^+
  // G^<_{bwd}(t, t')    = G^>_{fwd}(t', t)    = G^R_{fwd}(t', t) - [G^<_{fwd}(t, t')]^+
  // G^{tv}_{bwd}(t, \beta - \tau') = G^{vt}_{fwd}(tau', t) = - \xi [G^{tv}_{fwd}(t, tau')]^+
  //
  // NB! Only valid for matrix valued Green's functions
  // (and single components that are their own anti-hermitian conjugate on
  // time transpose.)

  typedef mam::dynamic_matrix_type matrix_type;

  cntr::herm_matrix_matrix_ref<matrix_type> gf_fwd(gf_fwd_in);
  cntr::herm_matrix_matrix_ref<matrix_type> gf_bwd(gf_bwd_out);

  int ntau = gf_fwd.ntau();

  assert( gf_fwd.sig() == gf_bwd.sig() );
  int xi = gf_fwd.sig();

  if(t == -1) {

    for( auto tau : range(0, ntau) )
      gf_bwd.mat(tau).noalias() = xi * gf_fwd.mat(ntau - 1 - tau);

  } else {
    
    for( auto j : range(0, t+1) ) {
      gf_bwd.ret_lt(t, j).noalias() = gf_fwd.ret(t, j).adjoint();
      gf_bwd.les_ut(j, t).noalias() = gf_fwd.gtr(t, j);
    }

    for( auto tau : range(0, ntau) )
      gf_bwd.tv(t, tau).noalias() = gf_fwd.vt(tau, t);
    
  } // t != -1
}  

// -----------------------------------------------------------------------
  void set_bwd_from_fwd_single_component(int t, gf_type & gf_bwd_out,
    gf_type & gf_fwd_in, gf_type & gf_fwd_in_orbital_transpose) {

  // Construct response function for backwards propagating case
  //
  // G^M_{bwd}(\tau - 0) = G^M_{fwd}(0 - \tau) = \xi G^M_{fwd}(\beta - \tau)
  // G^R_{bwd}(t, t')    = G^A_{fwd}(t', t)    = [G^R_{fwd}(t, t')]^+
  // G^<_{bwd}(t, t')    = G^>_{fwd}(t', t)    = G^R_{fwd}(t', t) - [G^<_{fwd}(t, t')]^+
  // G^{tv}_{bwd}(t, \beta - \tau') = G^{vt}_{fwd}(tau', t) = - \xi [G^{tv}_{fwd}(t, tau')]^+
  //
  // Special case for one-component Green's functions with provided orbital transpose.

    
  typedef mam::dynamic_matrix_type matrix_type;

  cntr::herm_matrix_matrix_ref<matrix_type> gf_fwd(gf_fwd_in);
  cntr::herm_matrix_matrix_ref<matrix_type> gf_fwd_orb_trans(gf_fwd_in_orbital_transpose);
  cntr::herm_matrix_matrix_ref<matrix_type> gf_bwd(gf_bwd_out);

  int ntau = gf_fwd.ntau();

  assert( gf_fwd.sig() == gf_bwd.sig() );
  int xi = gf_fwd.sig();

  if(t == -1) {

    for( auto tau : range(0, ntau) )
      gf_bwd.mat(tau).noalias() = xi * gf_fwd.mat(ntau - 1 - tau);

  } else {
    
    for( auto j : range(0, t+1) ) {
      gf_bwd.ret_lt(t, j).noalias() = gf_fwd_orb_trans.ret(t, j).adjoint();
      gf_bwd.les_ut(j, t).noalias() = gf_fwd.ret(t, j) - gf_fwd_orb_trans.les(j, t).adjoint();
    }

    for( auto tau : range(0, ntau) )
      gf_bwd.tv(t, tau).noalias() = gf_fwd_orb_trans.vt(tau, t);
    
  } // t != -1
}  

  
// -----------------------------------------------------------------------
void solve_pp_dyson(int tstp, ppgfs_type & ppGfs, ppgfs_type & ppSigmas,
		    functions_type H, double pp_mu, double beta, double h,
		    int kt, int nomp) {

  assert( ppGfs.size() == ppSigmas.size() );
  assert( ppGfs.size() == H.size() );

#pragma omp parallel for num_threads(nomp)
  for(int idx = 0; idx < ppGfs.size(); idx++) {
    if(tstp == -1) {

      //force_matsubara_hermitian(ppSigmas[idx]);

      ::cntr::pseudodyson_mat(
        ppGfs[idx], pp_mu, H[idx], ppSigmas[idx],
	integration::I<double>(kt), beta);

      /*
      dyson::ppsc_pseudodyson_mat(
        ppGfs[idx], pp_mu, H[idx], ppSigmas[idx],
	integration::I<double>(kt), beta);
      */
      
      //force_matsubara_hermitian(ppGfs[idx]);
      
    } else if(tstp == 0) {
      set_t0_from_mat(ppGfs);
    } else if(tstp <= kt) {
      ::cntr::pseudodyson_start(
	ppGfs[idx], pp_mu, H[idx], ppSigmas[idx],
	integration::I<double>(tstp), beta, h);
    } else {
      ::cntr::pseudodyson_timestep(
	tstp, ppGfs[idx], pp_mu, H[idx], ppSigmas[idx],
	integration::I<double>(kt), beta, h);
    }
  } // idx
  
}

// -----------------------------------------------------------------------
void set_t0_from_mat(ppgfs_type & ppGfs) {
  for( auto idx : range(0, ppGfs.size()) )
    ::cntr::set_t0_from_mat(ppGfs[idx]);
}

// -----------------------------------------------------------------------
mam::dynamic_matrix_type density_matrix(int tstp, ppgf_type & ppGf_in) {
  cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> ppGf(ppGf_in);

  // -- Seemingly the pp dyson solver gives transposed indices
  // -- here we compensate for that to follow the same convention

  /*
  if(tstp == -1) return -ppGf.mat(ppGf.ntau()-1).transpose();
  else return value_type(0., ppGf.sig()) * ppGf.les(tstp, tstp).transpose();
  */

  if(tstp == -1) return -ppGf.mat(ppGf.ntau()-1);
  else return value_type(0., ppGf.sig()) * ppGf.les(tstp, tstp);
}

// -----------------------------------------------------------------------
void update_density_matrix(int tstp, ppgfs_type & ppGfs, operators_type & rho) {

  for( auto idx : range(0, ppGfs.size()) )
    rho[tstp + 1].M_[idx] = density_matrix(tstp, ppGfs[idx]);
}

// -----------------------------------------------------------------------
double ppgfs_density(int tstp, ppgfs_type & ppGfs) {
  value_type dens = 0.0;
  for( auto idx : range(0, ppGfs.size()) ) 
    dens += (density_matrix(tstp, ppGfs[idx])).trace();
  return dens.real();
}

// -----------------------------------------------------------------------
void normalize_ppgfs_mat(ppgfs_type & ppGfs, double & pp_mu, double beta) {

  int tstp = -1;
  double dens = ppgfs_density(tstp, ppGfs);

  if(dens <= 0.0) {
    std::cerr << "--> normalize_ppgfs_mat: Can not normalize dens<=0." << std::endl;
    exit(0);
  }

  double dlam = - std::log(dens) / beta;

  for( auto idx : range(0, ppGfs.size()) ) {
    cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> ppGf(ppGfs[idx]);

    double dtau = beta / (ppGf.ntau()-1);
    
    for( auto m : range(0, ppGf.ntau()) ) {
      ppGf.mat(m) *= std::exp(dlam * (dtau * m));
    }
  }

  pp_mu += dlam;

  std::cout << "dens, pp_mu = " << dens << ", " << pp_mu << std::endl;
}

// -----------------------------------------------------------------------
value_type expectation_value(int tstp, operators_type & rho, operator_type & op){
  value_type exp_val = 0.0;
  for( auto idx : range(0, op.M_.size()) )
    if(op.to_sector_[idx] == idx)
      exp_val += (rho[tstp + 1].M_[idx] * op.M_[idx]).trace();
  return exp_val;
}

// -----------------------------------------------------------------------
void extrapolate_ppgf_timestep(int tstp, ppgfs_type & ppGfs, int kt) {
  int kt_extrap = (tstp <= kt ? tstp : kt);
  if(tstp == -1) set_t0_from_mat(ppGfs);
  for( auto idx : range(0, ppGfs.size()) )
    ::cntr::extrapolate_timestep(tstp, ppGfs[idx],
			       integration::I<double>(kt_extrap));
}

// -----------------------------------------------------------------------
void density_matrix_eigenvalues(int tstp, operators_type & rho) {

  std::cout << "--> density_matrix_eigenvalues" << std::endl;
  
  int ns = rho[tstp+1].M_.size();
  
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    std::vector<vector_type> evals(ns);
  
  for( auto idx : range(0, ns) ) {
    mam::dynamic_matrix_type rho_s = rho[tstp + 1].M_[idx];
    
    Eigen::SelfAdjointEigenSolver<mam::dynamic_matrix_type> eigh(rho_s);
    evals[idx] = eigh.eigenvalues();
    
    std::cout << "sector, evals = " << idx << " : " << evals[idx] << std::endl;
  }
  std::cout << std::endl;
}

// -----------------------------------------------------------------------
void zeroth_order(ppgfs_type & ppGfs, functions_type & H, double & pp_mu, double beta) {

  //std::cout << "--> zeroth order" << std::endl;
  
  int ns = ppGfs.size();

  double emin = 1e100;

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
  
  std::vector<vector_type> evals(ns);
  std::vector<mam::dynamic_matrix_type> evecs(ns);

  for( auto idx : range(0, ppGfs.size()) ) {
    cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> ppGf(ppGfs[idx]);
    double dtau = beta / (ppGf.ntau()-1);
    mam::dynamic_matrix_type H_sector;
    H[idx].get_value(-1, H_sector);

    Eigen::SelfAdjointEigenSolver<mam::dynamic_matrix_type> eigh(H_sector);
    evals[idx] = eigh.eigenvalues();
    evecs[idx] = eigh.eigenvectors();

    /*
    std::cout << "sector idx = " << idx << std::endl;
    std::cout << "H_sector = " << std::endl << H_sector << std::endl;
    std::cout << "evals = " << evals[idx] << std::endl;
    std::cout << "evecs = " << std::endl << evecs[idx] << std::endl;

    int n = ppGf.nflavour();
    mam::dynamic_matrix_type tmp(mam::dynamic_matrix_type::Zero(n, n));
    tmp.diagonal() = evals[idx].cast<value_type>();
    
    std::cout << "H_sector (reconstr) = " << std::endl
	      << evecs[idx] * tmp * evecs[idx].adjoint() << std::endl;
    */
 
    emin = std::min(emin, evals[idx].minCoeff());
  }
  
  // determine pseudo particle chemical potential pp_mu
  double pp_dens = 0.0;
  for( auto sector : range(0, evals.size()) ) {
    vector_type e = evals[sector];
    
    for( auto idx : range(0, e.rows()) ) {
      pp_dens += std::exp(-beta*(e(idx) - emin));
    }
  }

  pp_mu = emin - std::log(pp_dens)/beta;
  std::cout << "pp_dens, pp_mu = " << pp_dens << ", " << pp_mu << std::endl;

  for( auto sector : range(0, ppGfs.size()) ) {
    cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> ppGf(ppGfs[sector]);

    double dtau = beta / (ppGf.ntau()-1);

    for( auto m : range(0, ppGf.ntau()) ) {

      Eigen::DiagonalMatrix<value_type, Eigen::Dynamic> diag_mat(ppGf.nflavour());      

      double tau = dtau * m;
      auto exp_evals = (-tau * (evals[sector].array() - pp_mu)).exp();
      diag_mat.diagonal() = exp_evals.matrix().cast<value_type>();

      ppGf.mat(m) = - evecs[sector] * diag_mat * evecs[sector].adjoint();

      /*
      if(m == ppGf.ntau()-1) {
	//Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> diag_mat_print(diag_mat);
	std::cout << "sector = " << sector << std::endl;
	//std::cout << "diag_mat = " << std::endl << diag_mat_print << std::endl;
	std::cout << "evecs[sector] = " << std::endl << evecs[sector] << std::endl;
	std::cout << "partial product = " << std::endl
		  << evecs[sector].adjoint() * diag_mat << std::endl;
	std::cout << "full product = " << std::endl
		  << evecs[sector].adjoint() * diag_mat * evecs[sector] << std::endl;
	std::cout << "ppGf.mat(m) = " << std::endl << ppGf.mat(m) << std::endl;
      }
      */
    }
  }
}

// -----------------------------------------------------------------------
std::vector<double> get_kinetic_energy(gf_type & Gloc, gf_type & Delta,
				       double beta, double h, int kt) {

  std::vector<double> Ekin;

  if(Gloc.nt() == 0) return Ekin; // Do not attempt to convolve equilibrium only.

  gf_type tmp(Gloc.nt(), Gloc.ntau(), Gloc.size1(), Gloc.sig());
  
  ::cntr::convolution(tmp, Gloc, Gloc, Delta, Delta,
		      integration::I<double>(kt), beta, h);

  cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> tmp_ref(tmp);
  
  Ekin.resize(Gloc.nt() + 2);
  Ekin[0] = -tmp_ref.mat(tmp_ref.ntau() - 1).trace().real();
  for( auto t : range(0, Gloc.nt()+1) )
    Ekin[t + 1] = -Gloc.sig() * tmp_ref.les(t, t).trace().imag();

  return Ekin;
}

// -----------------------------------------------------------------------
void force_matsubara_hermitian(ppgf_type & G_in) {
  cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> G(G_in);
  for( auto tau : range(0, G.ntau()) ) {
    G.mat(tau) = 0.5 * (G.mat(tau) + G.mat(tau).adjoint());
  }
}
  
// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_CPP
