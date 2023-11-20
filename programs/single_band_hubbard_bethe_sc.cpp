// -----------------------------------------------------------------------

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

// -----------------------------------------------------------------------

#ifndef CNTR_USE_OMP
#define CNTR_USE_OMP
#endif

#ifndef CNTR_USE_MPI
#define CNTR_USE_MPI
#endif

#define NCA_SOLVER_ASSERT_0 0
#define NCA_SOLVER_ASSERT_1 0

#include <cntr/cntr.hpp>
#include <cntr/utils/read_inputfile.hpp>

// -----------------------------------------------------------------------

#include "./ppsc/ppsc.hpp"
#include "./ppsc/solver.hpp"

#include "./ppsc/hilbert_spaces/single_band_fermi_sc.hpp"
#include "./ppsc/hamiltonians/single_band_hubbard_sc.hpp"

// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & Delta,
			       ppsc::gf_type & Delta_cc,
			       HILB & hil_) {

  // spin (u)p/(d)own and (c)reation/(a)nihilation operators

  ppsc::operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];
  
  int fermion=-1, fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;

  pp_ints.push_back(ppsc::pp_int_type(Delta,    0, 0, cuc, cua, fermion, fwd)); 
  pp_ints.push_back(ppsc::pp_int_type(Delta_cc, 0, 0, cua, cuc, fermion, bwd)); 

  pp_ints.push_back(ppsc::pp_int_type(Delta,    0, 1, cuc, cdc, fermion, fwd)); 
  pp_ints.push_back(ppsc::pp_int_type(Delta_cc, 0, 1, cdc, cuc, fermion, bwd)); 

  pp_ints.push_back(ppsc::pp_int_type(Delta,    1, 0, cda, cua, fermion, fwd)); 
  pp_ints.push_back(ppsc::pp_int_type(Delta_cc, 1, 0, cua, cda, fermion, bwd)); 
  
  pp_ints.push_back(ppsc::pp_int_type(Delta,    1, 1, cda, cdc, fermion, fwd)); 
  pp_ints.push_back(ppsc::pp_int_type(Delta_cc, 1, 1, cdc, cda, fermion, bwd)); 

  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::gf_verts_type get_gf_verts(HILB & hil_) {

  ppsc::operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  ppsc::gf_verts_type gf_verts;

  // single-band superconducting Nambu notation \Psi^+ = (cuc, cda)
  // g = -i<\Psi \Psi^+>

  gf_verts.push_back(ppsc::gf_vert_type(0, 0, cua, cuc)); 
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, cua, cda)); 
  gf_verts.push_back(ppsc::gf_vert_type(1, 0, cdc, cuc)); 
  gf_verts.push_back(ppsc::gf_vert_type(1, 1, cdc, cda));
  
  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int itermax, iter_rtime, nt, ntau, kt, iter, tstp, order, nomp;

  int store_pp, read_eq_sym, read_rt_sym, read_state_from_file;
  double linear_mixing;
  
  double beta, h, errmax, dmfterr, dmfterr_equil;
  bool matsubara_converged = false;

  try {
    // ---------------------------------------------------------------------
    // READ GENERAL INPUT (NOT YET MICROSCOPIC PARAMETERS)
    {
      if (argc < 2)
        throw("COMMAND LINE ARGUMENT MISSING");
      // scan the input file, double underscores to avoids mismatch
      find_param(argv[1], "__nt=", nt);
      find_param(argv[1], "__ntau=", ntau);
      find_param(argv[1], "__beta=", beta);
      find_param(argv[1], "__h=", h);
      find_param(argv[1], "__itermax=", itermax);
      find_param(argv[1], "__errmax=", errmax);
      find_param(argv[1], "__iter_rtime=", iter_rtime);
      find_param(argv[1], "__kt=", kt);
      find_param(argv[1], "__order=", order);
      find_param(argv[1], "__store_pp=", store_pp);
      find_param(argv[1], "__linear_mixing=", linear_mixing);
      find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__read_state_from_file=", read_state_from_file);
      find_param(argv[1], "__nomp=", nomp);
    }

    // ---------------------------------------------------------------------
    // -- Setup pp calculator

    typedef ppsc::hilbert_spaces::single_band_fermi_sc hilbert_space_type;
    typedef ppsc::hamiltonians::single_band_hubbard_sc<hilbert_space_type> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    hilbert_space.init();

    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);

    { // -- Setup local Hamiltonian

      find_param_tvector(argv[1], "__mu=", imp.hamiltonian.mu, nt);
      find_param(argv[1], "__mu0=", imp.hamiltonian.mu[0]);

      find_param_tvector(argv[1], "__U=", imp.hamiltonian.U, nt);
      find_param(argv[1], "__U0=", imp.hamiltonian.U[0]);
 
      find_param_tvector(argv[1], "__Bz=", imp.hamiltonian.Bz, nt);
      find_param(argv[1], "__Bz0=", imp.hamiltonian.Bz[0]);
      
      find_param_tvector(argv[1], "__eps=", imp.hamiltonian.epsilon, nt);
      find_param_tvector(argv[1], "__pair_real=", imp.hamiltonian.pair_real, nt);
      find_param_tvector(argv[1], "__pair_imag=", imp.hamiltonian.pair_imag, nt);
      
      imp.update_hamiltonian();
    }
    
    // ---------------------------------------------------------------------
    // -- setup single particle greens functions and hybridizations

    int dim = 2;
    int fermion = -1;
    cntr::herm_matrix<double> Gloc(nt, ntau, dim, fermion);    
    cntr::herm_matrix<double> Delta(nt, ntau, dim, fermion);
    cntr::herm_matrix<double> Delta_cc(nt, ntau, dim, fermion);

    cntr::function<double> sigma_z(nt, 2);
    {
      ppsc::mam::dynamic_matrix_type sigma_z_mat(2, 2);
      sigma_z_mat *= 0.0;
      sigma_z_mat(0, 0) = 1.0;
      sigma_z_mat(1, 1) = -1.0;

      for(int tstp = -1; tstp <= nt; tstp++) {
	sigma_z.set_value(tstp, sigma_z_mat);
      }
    }
    
    // ---------------------------------------------------------------------
    
    if(read_state_from_file) {

      std::string filename = "data_ppsc.h5";
      std::cout << "--> Reading state from file: " << filename << std::endl;

      Gloc.read_from_hdf5(filename.c_str(), "g");
      Delta.read_from_hdf5(filename.c_str(), "d");
      Delta_cc.read_from_hdf5(filename.c_str(), "dcc");

      imp.load(filename);
      
      if(read_eq_sym) {
	std::string filename = "sym_eq.txt";
	imp.read_symmetries(filename);
      }

      imp.update_density_matrix(-1);
      imp.hamiltonian.update_exp_vals(-1, imp.rho);

      std::cout << "pp_mu = " << imp.pp_mu << std::endl;
      
    }

    // ---------------------------------------------------------------------
    // MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
    {
      cntr::herm_matrix<double> gtmp(-1, ntau, dim, fermion);

      if(!read_state_from_file) {
	imp.solve_atomic();
	imp.update_density_matrix(-1);
	imp.hamiltonian.update_exp_vals(-1, imp.rho);	
      } else {
	gtmp.set_timestep(-1, Gloc); // this avoids one iteration
      }

      for (iter = 1; iter <= itermax; iter++) {

	// -- Construct interactions and verticies
	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);

	// -- Solve pseudo particle problem
	imp.update_diagrams(pp_ints, gf_verts);

	if(iter == 1 && read_eq_sym && !imp.has_symmetries()) {
	  std::string filename = "sym_eq.txt";
	  imp.read_symmetries(filename);
	}

	imp.pp_step(-1);

	// -- Update Gloc
	ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);

	for( auto idx : ppsc::range(0, gf_tstps.size()) ) {
	  Gloc.set_matrixelement(
	    -1, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps[idx]);
	}
	
	// -- Check error
	dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc);
        gtmp.set_timestep(-1, Gloc);

	// -- Update Hybridization Delta = \sigma_z Gloc \sigma_z
	cntr::herm_matrix_timestep<double> tmp;
	Gloc.get_timestep(-1, tmp);
	Delta.set_timestep(-1, tmp);
	Delta.left_multiply(-1, sigma_z);
	Delta.right_multiply(-1, sigma_z);
	ppsc::set_bwd_from_fwd(-1, Delta_cc, Delta);
	
	cout << "iter:  " << iter << " err: " << dmfterr_equil << endl;
	
        if (dmfterr_equil < errmax) {

	  /*
	  if(!imp.has_symmetries()) {
	    imp.symmetry_reduction(-1);	    
	    std::string filename = "sym_eq.txt";
	    imp.write_symmetries(filename);
	  }
	  */
	  
	  matsubara_converged = true;
          break;
        }

      }
      if (iter > itermax) {
        cerr << "WARNING: Matsubara not converged  after " << itermax
             << "steps ... abort" << endl;
        cerr << "skip real-time calculation " << endl;
      }
    }

    // ---------------------------------------------------------------------

    /*
    if(nt > 0) imp.clear_symmetries();

    if(read_rt_sym) {
      std::string filename = "sym_rt.txt";
      imp.read_symmetries(filename);
    }
    */

    // ---------------------------------------------------------------------
    // 	START ... same iteration
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      cntr::herm_matrix<double> gtmp(kt, ntau, dim, fermion);
      imp.init_real_time();
      
      for (iter = 1; iter <= itermax; iter++) {

	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	
	imp.update_diagrams(pp_ints, gf_verts);
	imp.pp_step(kt);

	for (int n = 0; n <= kt; n++) {

	  // -- Update Gloc
	  ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);

	  for( auto idx : ppsc::range(0, gf_tstps.size()) ) {
	    Gloc.set_matrixelement(
	      n, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps[idx]);
	  }

	  // -- Update Hybridization Delta = \sigma_z Gloc \sigma_z
	  cntr::herm_matrix_timestep<double> tmp;
	  Gloc.get_timestep(n, tmp);
	  Delta.set_timestep(n, tmp);
	  Delta.left_multiply(n, sigma_z);
	  Delta.right_multiply(n, sigma_z);
	  ppsc::set_bwd_from_fwd(n, Delta_cc, Delta);
	}
	
	// -- Check error
        dmfterr = cntr::distance_norm2(kt, gtmp, Gloc);
        gtmp.set_timestep(kt, Gloc);
	
        cout << "START: iter:  " << iter << " err: " << dmfterr << endl;
        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }
      }
    }

    // ---------------------------------------------------------------------
    // 	REALTIME: ONLY FOR NT>0
    for (tstp = kt + 1; tstp <= nt; tstp++) {
      cout << "tstp:  " << tstp << endl;
      imp.extrapolate_timestep(tstp - 1);
      
      for (iter = 1; iter <= iter_rtime; iter++) {

	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);

	imp.update_diagrams(pp_ints, gf_verts);
	imp.pp_step(tstp);

	// -- Update Gloc
	ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);

	for( auto idx : ppsc::range(0, gf_tstps.size()) ) {
	  Gloc.set_matrixelement(
	    tstp, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps[idx]);
	}

	// -- Update Hybridization Delta = \sigma_z Gloc \sigma_z
	cntr::herm_matrix_timestep<double> tmp;
	Gloc.get_timestep(tstp, tmp);
	Delta.set_timestep(tstp, tmp);
	Delta.left_multiply(tstp, sigma_z);
	Delta.right_multiply(tstp, sigma_z);
	ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);	
      }
    }

    // ---------------------------------------------------------------------
    // -- Reduce number of diagrams using symmetries
    /*
    if(matsubara_converged && nt > 0 && !imp.has_symmetries()) {
      imp.symmetry_reduction(kt);
      std::string filename = "sym_rt.txt";
      imp.write_symmetries(filename);
    }
    */

    // ---------------------------------------------------------------------
    // Kinetic energy

    std::vector<double> Ekin = ppsc::get_kinetic_energy(Gloc, Delta, beta, h, kt);

    // ---------------------------------------------------------------------
    // OBSERVABLES
    {
      std:string filename = "data_ppsc.h5";
      hid_t file_id = open_hdf5_file(filename);
      hid_t group_id;

      imp.store(file_id, store_pp);
      
      group_id = create_group(file_id, "g");
      store_herm_greens_function(group_id, Gloc);
      close_group(group_id);

      group_id = create_group(file_id, "d");
      store_herm_greens_function(group_id, Delta);
      close_group(group_id);

      group_id = create_group(file_id, "dcc");
      store_herm_greens_function(group_id, Delta_cc);
      close_group(group_id);

      group_id = create_group(file_id, "bethe");
      store_real_data_to_hid(group_id, "Ekin", Ekin.data(), Ekin.size());
      store_double_attribute_to_hid(group_id, "dmfterr_equil", dmfterr_equil);
      close_group(group_id);      
      
      close_hdf5_file(file_id);
    }
    
  } // try
  catch (char *message) {
    cerr << "exception\n**** " << message << " ****" << endl;
    cerr << "CDMFT input_file [ --test ]\n" << endl;
  } catch (...) {
    cerr << "unspecified exception " << endl;
    cerr << "\nCDMFT input_file [ --test ]\n" << endl;
  }
  return 0;
}
