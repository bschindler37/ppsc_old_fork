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

//#include "./ppsc/hilbert_spaces/two_band_fermi_diag.hpp"
#include "./ppsc/hilbert_spaces/two_band_fermi_densdens.hpp"
#include "./ppsc/hamiltonians/two_band_hubbard.hpp"

// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------

typedef ppsc::operator_type operator_type;
typedef ppsc::mam::dynamic_matrix_type matrix_type; 

// -----------------------------------------------------------------------

class single_particle_greensfunction_type {

public:

  single_particle_greensfunction_type(int nt, int ntau) : nt(nt), ntau(ntau),
    o1up(nt, ntau, 1, -1), o1do(nt, ntau, 1, -1),
    o2up(nt, ntau, 1, -1), o2do(nt, ntau, 1, -1) {}

  void update(int tstp, ppsc::gf_tstps_type & gf_tstps, double linear_mixing=0.0) {

    if(tstp == -1) {
      ppsc::gf_tstp_type gloc_old(-1, ntau, 1);
      ppsc::gf_tstp_type gloc_mix(-1, ntau, 1);

      gloc_mix.clear();
      this->o1up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[0], 1.0 - linear_mixing);
      this->o1up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o1do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[1], 1.0 - linear_mixing);
      this->o1do.set_timestep(-1, gloc_mix);
      
      gloc_mix.clear();
      this->o2up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[2], 1.0 - linear_mixing);
      this->o2up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o2do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[3], 1.0 - linear_mixing);
      this->o2do.set_timestep(-1, gloc_mix);
      
    } else {

      this->o1up.set_timestep(tstp, gf_tstps[0]);
      this->o1do.set_timestep(tstp, gf_tstps[1]);
      this->o2up.set_timestep(tstp, gf_tstps[2]);
      this->o2do.set_timestep(tstp, gf_tstps[3]);
      
    }
  }

  void store(hid_t file_id) {
    hid_t group_id;
    
    group_id = create_group(file_id, "g1u");
    store_herm_greens_function(group_id, this->o1up);
    close_group(group_id);
    
    group_id = create_group(file_id, "g1d");
    store_herm_greens_function(group_id, this->o1do);
    close_group(group_id);

    group_id = create_group(file_id, "g2u");
    store_herm_greens_function(group_id, this->o2up);
    close_group(group_id);
    
    group_id = create_group(file_id, "g2d");
    store_herm_greens_function(group_id, this->o2do);
    close_group(group_id);
    
  }

  void load(std::string filename) {
    this->o1up.read_from_hdf5(filename.c_str(), "g1u");
    this->o1do.read_from_hdf5(filename.c_str(), "g1d");
    this->o2up.read_from_hdf5(filename.c_str(), "g2u");
    this->o2do.read_from_hdf5(filename.c_str(), "g2d");
  }

  int nt, ntau;
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do;
    
};

// -----------------------------------------------------------------------

class hybridization_function_type {

public:
  
  hybridization_function_type(int nt, int ntau) :
    nt(nt), ntau(ntau),
    o1up(nt, ntau, 1, -1),
    o1do(nt, ntau, 1, -1),
    o2up(nt, ntau, 1, -1),
    o2do(nt, ntau, 1, -1),
    // -- bwd hybridizations
    o1up_cc(nt, ntau, 1, -1),
    o1do_cc(nt, ntau, 1, -1),
    o2up_cc(nt, ntau, 1, -1),
    o2do_cc(nt, ntau, 1, -1)
  {}    

  void update(int tstp, single_particle_greensfunction_type & Gloc,
	      cntr::function<double> & tfunc, int af_bethe_flag=0) {

    cntr::herm_matrix_timestep<double> tmp;
    
    // -- Update Delta.o1up
      
    if(af_bethe_flag == 1)
      Gloc.o1do.get_timestep(tstp, tmp);
    else
      Gloc.o1up.get_timestep(tstp, tmp);
    
    this->o1up.set_timestep(tstp, tmp);
    this->o1up.left_multiply( tstp, tfunc);
    this->o1up.right_multiply(tstp, tfunc);	
    ppsc::set_bwd_from_fwd(tstp, this->o1up_cc, this->o1up);
    
    // -- Update Delta.o1do
    
    if(af_bethe_flag == 1)
      Gloc.o1up.get_timestep(tstp, tmp);
    else
      Gloc.o1do.get_timestep(tstp, tmp);
    
    this->o1do.set_timestep(tstp, tmp);
    this->o1do.left_multiply( tstp, tfunc);
    this->o1do.right_multiply(tstp, tfunc);	
    ppsc::set_bwd_from_fwd(tstp, this->o1do_cc, this->o1do);
    
    // -- Update Delta.o2up
    
    if(af_bethe_flag == 1)
      Gloc.o2do.get_timestep(tstp, tmp);
    else
      Gloc.o2up.get_timestep(tstp, tmp);
    
    this->o2up.set_timestep(tstp, tmp);
    this->o2up.left_multiply( tstp, tfunc);
    this->o2up.right_multiply(tstp, tfunc);	
    ppsc::set_bwd_from_fwd(tstp, this->o2up_cc, this->o2up);
    
    // -- Update Delta.o2do
    
    if(af_bethe_flag == 1)
      Gloc.o2up.get_timestep(tstp, tmp);
    else
      Gloc.o2do.get_timestep(tstp, tmp);
    
    this->o2do.set_timestep(tstp, tmp);
    this->o2do.left_multiply( tstp, tfunc);
    this->o2do.right_multiply(tstp, tfunc);	
    ppsc::set_bwd_from_fwd(tstp, this->o2do_cc, this->o2do);
    
  }

  void store(hid_t file_id) {
    hid_t group_id;
    
    group_id = create_group(file_id, "d1u");
    store_herm_greens_function(group_id, this->o1up);
    close_group(group_id);
    
    group_id = create_group(file_id, "d1d");
    store_herm_greens_function(group_id, this->o1do);
    close_group(group_id);

    group_id = create_group(file_id, "d2u");
    store_herm_greens_function(group_id, this->o2up);
    close_group(group_id);
    
    group_id = create_group(file_id, "d2d");
    store_herm_greens_function(group_id, this->o2do);
    close_group(group_id);

    group_id = create_group(file_id, "d1ucc");
    store_herm_greens_function(group_id, this->o1up_cc);
    close_group(group_id);
    
    group_id = create_group(file_id, "d1dcc");
    store_herm_greens_function(group_id, this->o1do_cc);
    close_group(group_id);

    group_id = create_group(file_id, "d2ucc");
    store_herm_greens_function(group_id, this->o2up_cc);
    close_group(group_id);
    
    group_id = create_group(file_id, "d2dcc");
    store_herm_greens_function(group_id, this->o2do_cc);
    close_group(group_id);
    
  }

  void load(std::string filename) {
    this->o1up.read_from_hdf5(filename.c_str(), "d1u");
    this->o1do.read_from_hdf5(filename.c_str(), "d1d");
    this->o2up.read_from_hdf5(filename.c_str(), "d2u");
    this->o2do.read_from_hdf5(filename.c_str(), "d2d");

    this->o1up_cc.read_from_hdf5(filename.c_str(), "d1ucc");
    this->o1do_cc.read_from_hdf5(filename.c_str(), "d1dcc");
    this->o2up_cc.read_from_hdf5(filename.c_str(), "d2ucc");
    this->o2do_cc.read_from_hdf5(filename.c_str(), "d2dcc");
  }
  
  int nt, ntau;
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do;
  cntr::herm_matrix<double> o1up_cc, o1do_cc, o2up_cc, o2do_cc;
  
};

// -----------------------------------------------------------------------
template<class HAM>
ppsc::pp_ints_type get_pp_ints(hybridization_function_type & Delta, HAM & h) {

  // spin (u)p/(d)own and (c)reation/(a)nihilation operators

  int boson=+1, fermion=-1, fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  
  // orbital no. 1

  pp_ints.push_back(ppsc::pp_int_type(Delta.o1up,    h.c1uc,  h.c1ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o1up_cc, h.c1ua,  h.c1uc,  fermion, bwd)); // spin up bwd
  
  pp_ints.push_back(ppsc::pp_int_type(Delta.o1do,    h.c1dc,  h.c1da,  fermion, fwd)); // spin do fwd 
  pp_ints.push_back(ppsc::pp_int_type(Delta.o1do_cc, h.c1da,  h.c1dc,  fermion, bwd)); // spin do bwd 

  // orbital no. 2
  
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2up,    h.c2uc,  h.c2ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2up_cc, h.c2ua,  h.c2uc,  fermion, bwd)); // spin up bwd
  
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2do,    h.c2dc,  h.c2da,  fermion, fwd)); // spin do fwd 
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2do_cc, h.c2da,  h.c2dc,  fermion, bwd)); // spin do bwd 
  
  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HAM>
ppsc::gf_verts_type get_gf_verts(HAM & h) {

  ppsc::gf_verts_type gf_verts;

  gf_verts.push_back(ppsc::gf_vert_type(0, 0, h.c1ua, h.c1uc)); // spin up orb 1
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, h.c1da, h.c1dc)); // spin do orb 1

  gf_verts.push_back(ppsc::gf_vert_type(1, 1, h.c2ua, h.c2uc)); // spin up orb 2
  gf_verts.push_back(ppsc::gf_vert_type(1, 1, h.c2da, h.c2dc)); // spin do orb 2

  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int itermax, iter_rtime, nt, ntau, kt, iter_rt, iter_equil, iter_warmup, tstp, order, nomp, store_pp;
  double beta, h, errmax, dmfterr, V, mu, omega, dmfterr_equil, linear_mixing, sym_tol, Bz_seed;
  int read_eq_sym, read_rt_sym, read_state_from_file, af_bethe_flag;
  bool matsubara_converged = false;

  std::vector<double> t_vec;

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
      find_param(argv[1], "__sym_tol=", sym_tol);
      find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__read_state_from_file=", read_state_from_file);

      find_param(argv[1], "__af_bethe_flag=", af_bethe_flag);
      
      if (order > 1) {
        cout << "using an OCA impurity solver" << endl;
        find_param(argv[1], "__nomp=", nomp);
      } else {
        nomp = 1;
      }
    }
    // ---------------------------------------------------------------------
    // -- Setup pp calculator

    typedef ppsc::hilbert_spaces::two_band_fermi_densdens hilbert_space_type;
    //typedef ppsc::hilbert_spaces::two_band_fermi_diag hilbert_space_type;
    typedef ppsc::hamiltonians::two_band_hubbard<hilbert_space_type> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    hilbert_space.init();
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    
    { // -- Setup local Hamiltonian
      find_param(argv[1], "__mu=", imp.hamiltonian.mu);

      find_param_tvector(argv[1], "__U=", imp.hamiltonian.U, nt);
      find_param(argv[1], "__U0=", imp.hamiltonian.U[0]);

      find_param_tvector(argv[1], "__J=", imp.hamiltonian.J, nt);
      find_param(argv[1], "__J0=", imp.hamiltonian.J[0]);
      
      find_param_tvector(argv[1], "__delta=", imp.hamiltonian.delta, nt);
      find_param(argv[1], "__delta0=", imp.hamiltonian.delta[0]);

      find_param_tvector(argv[1], "__t=", t_vec, nt);
      find_param(argv[1], "__t0=", t_vec[0]);

      find_param_tvector(argv[1], "__Bz=", imp.hamiltonian.Bz, nt);
      find_param(argv[1], "__Bz0=", imp.hamiltonian.Bz[0]);
      find_param(argv[1], "__Bz_seed=", Bz_seed);
      
      imp.update_hamiltonian();
    }
    // ---------------------------------------------------------------------
    // -- init single particle gf, hyb, and hopping

    single_particle_greensfunction_type Gloc(nt, ntau);
    hybridization_function_type Delta(nt, ntau);
    
    cntr::function<double> tfunc(nt);
    // set tfunc from t_vec
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc[tstp] = t_vec[tstp+1];
    }

    // ---------------------------------------------------------------------
    
    if(read_state_from_file) {

      std::string filename = "data_ppsc.h5";
      std::cout << "--> Reading state from file: " << filename << std::endl;

      Gloc.load(filename);
      Delta.load(filename);
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
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);

      if(!read_state_from_file) {
	imp.solve_atomic();
      } else {
	gtmp.set_timestep(-1, Gloc.o1up); // this avoids one iteration
      }
      
      for (iter_equil = 1; iter_equil <= itermax; iter_equil++) {

        // -- Construct interactions and verticies
        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);
        
        // -- Solve pseudo particle problem
        imp.update_diagrams(pp_ints, gf_verts);

	if(iter_equil == 1 && read_eq_sym && !imp.has_symmetries()) {
	  std::string filename = "sym_eq.txt";
	  imp.read_symmetries(filename);
	}

	if(iter_equil == 1) {
	  std::cout << "--> Applying Bz_seed = " << Bz_seed << std::endl;
	  double tmp = imp.hamiltonian.Bz[0]; // store equil value
	  imp.hamiltonian.Bz[0] = Bz_seed;
	  Bz_seed = tmp;
	} else if(iter_equil == 2) {
	  imp.hamiltonian.Bz[0] = Bz_seed; // Restore equil value	  
	}

	imp.pp_step(-1);

        // -- get spgf
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);
	Gloc.update(-1, gf_tstps, linear_mixing);

        // -- Check error
        dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc.o1up);
        gtmp.set_timestep(-1, Gloc.o1up);
	
        // -- Update Hybridization
	Delta.update(-1, Gloc, tfunc, af_bethe_flag);
	
        cout << "iter:  " << iter_equil << " err: " << dmfterr_equil << endl;

        if(dmfterr_equil < errmax) {
          matsubara_converged = true;
          break;
        }

      }
      if (iter_equil > itermax) {
        cerr << "WARNING: Matsubara not converged  after " << itermax
             << "steps ... abort" << endl;
        cerr << "skip real-time calculation " << endl;
      }
    }

    // ---------------------------------------------------------------------

    if(!imp.has_symmetries()) {
      imp.symmetry_reduction(-1);
      std::string filename = "sym_eq.txt";
      imp.write_symmetries(filename);
    }
	
    if(nt > 0) imp.clear_symmetries();

    if(read_rt_sym) {
      std::string filename = "sym_rt.txt";
      imp.read_symmetries(filename);
    }
    
    // ---------------------------------------------------------------------
    //  START ... same iteration
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      
      for (iter_warmup = 1; iter_warmup <= itermax; iter_warmup++) {

        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);
        
        imp.update_diagrams(pp_ints, gf_verts);
        imp.pp_step(kt);

        for (int n = 0; n <= kt; n++) {
          ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
	  Gloc.update(n, gf_tstps);
	  Delta.update(n, Gloc, tfunc, af_bethe_flag);
        }
        
        dmfterr = cntr::distance_norm2(kt, gtmp, Gloc.o1up);
        gtmp.set_timestep(kt, Gloc.o1up);
        
        cout << "WARMUP: iter:  " << iter_warmup << " err: " << dmfterr << endl;
        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }
      }
    }

    // ---------------------------------------------------------------------
    // -- Reduce number of diagrams using symmetries

    if(matsubara_converged && nt > 0 && !imp.has_symmetries()) {
      imp.symmetry_reduction(kt);
      std::string filename = "sym_rt.txt";
      imp.write_symmetries(filename);
    }
    
    // ---------------------------------------------------------------------
    //  REALTIME: ONLY FOR NT>0
    if(matsubara_converged)
    for (tstp = kt + 1; tstp <= nt; tstp++) {

      imp.extrapolate_timestep(tstp - 1);
      
      for (iter_rt = 1; iter_rt <= iter_rtime; iter_rt++) {

        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);

        imp.update_diagrams(pp_ints, gf_verts);
        imp.pp_step(tstp);

        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
	Gloc.update(tstp, gf_tstps);
	Delta.update(tstp, Gloc, tfunc, af_bethe_flag);
	
      }
    }

    // ---------------------------------------------------------------------
    // Kinetic energy

    //imp.calculate_pp_interaction_energy();
      
    std::vector<double> Ekin_o1up = ppsc::get_kinetic_energy(Gloc.o1up, Delta.o1up, beta, h, kt);
    std::vector<double> Ekin_o1do = ppsc::get_kinetic_energy(Gloc.o1do, Delta.o1do, beta, h, kt);
    std::vector<double> Ekin_o2up = ppsc::get_kinetic_energy(Gloc.o2up, Delta.o2up, beta, h, kt);
    std::vector<double> Ekin_o2do = ppsc::get_kinetic_energy(Gloc.o2do, Delta.o2do, beta, h, kt);

    std::vector<double> Ekin;
    Ekin.resize(Ekin_o1up.size());
    
    for( auto t : ppsc::range(0, Ekin.size()) )
      Ekin[t] = Ekin_o1up[t] + Ekin_o1do[t] + Ekin_o2up[t] + Ekin_o2do[t];

    // ---------------------------------------------------------------------
    // OBSERVABLES
    {
      std:string filename = "data_ppsc.h5";
      hid_t file_id = open_hdf5_file(filename);
      hid_t group_id;

      imp.store(file_id, store_pp);
      Gloc.store(file_id);
      Delta.store(file_id);

      // -- Bethe general properties
      group_id = create_group(file_id, "bethe");

      store_real_data_to_hid(group_id, "t_vec", t_vec.data(), t_vec.size());

      store_real_data_to_hid(group_id, "Ekin", Ekin.data(), Ekin.size());

      store_real_data_to_hid(group_id, "Ekin_o1up", Ekin_o1up.data(), Ekin_o1up.size());
      store_real_data_to_hid(group_id, "Ekin_o1do", Ekin_o1do.data(), Ekin_o1do.size());
      store_real_data_to_hid(group_id, "Ekin_o2up", Ekin_o2up.data(), Ekin_o2up.size());
      store_real_data_to_hid(group_id, "Ekin_o2do", Ekin_o2do.data(), Ekin_o2do.size());

      store_double_attribute_to_hid(group_id, "dmfterr_equil", dmfterr_equil);
      store_int_attribute_to_hid(group_id, "iter_equil", iter_equil);

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

