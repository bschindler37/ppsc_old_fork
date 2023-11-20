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

#include "ppsc/hilbert_spaces/two_band_fermi_tJ.hpp"
#include "ppsc/hamiltonians/two_band_tJ.hpp"

#include "./ppsc/baths/non_int_boson_propagator.hpp"

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
//    o3up(nt, ntau, 1, -1), o3do(nt, ntau, 1, -1) {}

  void update(int tstp, ppsc::gf_tstps_type & gf_tstps, double linear_mixing=0.0) {

    if(tstp == -1) {
      ppsc::gf_tstp_type gloc_old(-1, ntau, 1);
      ppsc::gf_tstp_type gloc_mix(-1, ntau, 1);

      gloc_mix.clear();
      this->o1up.get_timestep(-1, gloc_old);
      //this->o1do.get_timestep(-1, gloc_old); // afm
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[0], 1.0 - linear_mixing);
      this->o1up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o1do.get_timestep(-1, gloc_old);
      //this->o1up.get_timestep(-1, gloc_old); // afm
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[1], 1.0 - linear_mixing);
      this->o1do.set_timestep(-1, gloc_mix);
      
      gloc_mix.clear();
      this->o2up.get_timestep(-1, gloc_old);
      //this->o2do.get_timestep(-1, gloc_old); // afm
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[2], 1.0 - linear_mixing);
      this->o2up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o2do.get_timestep(-1, gloc_old);
      //this->o2up.get_timestep(-1, gloc_old); // afm
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
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do; // o3up, o3do;
    
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
	      cntr::function<double> & tfunc1sin, cntr::function<double> & tfunc2sin, //cntr::function<double> & tfunc3sin,
              cntr::function<double> & tfunc1cos, cntr::function<double> & tfunc2cos, //cntr::function<double> & tfunc3cos,
	      cntr::herm_matrix<double> & D0, // boson bath
	      int af_bethe_flag=0, int bath_flag=0) {

    // -- Update Delta.o1up
    if(af_bethe_flag == 1) _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);
    else                   _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);

    // -- Update Delta.o1do
    if(af_bethe_flag == 1) _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);
    else                   _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);

    // -- Update Delta.o2up
    if(af_bethe_flag == 1) _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o2do, tfunc2sin, tfunc2cos, D0, bath_flag);
    else                   _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o2up, tfunc2sin, tfunc2cos, D0, bath_flag);

    // -- Update Delta.o2do
    if(af_bethe_flag == 1) _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o2up, tfunc2sin, tfunc2cos, D0, bath_flag);
    else                   _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o2do, tfunc2sin, tfunc2cos, D0, bath_flag);

  }

  void _component_update(int tstp,
    cntr::herm_matrix<double> & Delta_comp,
    cntr::herm_matrix<double> & Delta_comp_cc,
    cntr::herm_matrix<double> & Gloc_comp,
    cntr::function<double> & tsin,
    cntr::function<double> & tcos,
    cntr::herm_matrix<double> & D0, // boson bath
    int bath_flag=0) {

    cntr::herm_matrix_timestep<double> tmp;
    
    Gloc_comp.get_timestep(tstp, tmp);
    
    Delta_comp.set_timestep(tstp, tmp);
    Delta_comp.left_multiply(tstp, tsin);
    Delta_comp.right_multiply(tstp, tsin);

    cntr::herm_matrix_timestep<double> tmpsin;
    Delta_comp.get_timestep(tstp, tmpsin);

    Delta_comp.set_timestep(tstp, tmp);
    Delta_comp.left_multiply(tstp, tcos);
    Delta_comp.right_multiply(tstp, tcos);

    Delta_comp.incr_timestep(tstp, tmpsin);

    if(bath_flag) { // add boson bath
      // Bubble2(tstp, tmp, Gloc_comp, D0);
      Delta_comp.incr_timestep(tstp, D0);
    }
    
    ppsc::set_bwd_from_fwd(tstp, Delta_comp_cc, Delta_comp);  
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
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do; // o3up, o3do;
  cntr::herm_matrix<double> o1up_cc, o1do_cc, o2up_cc, o2do_cc; // o3up_cc, o3do_cc;
  
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

  // Add the Sz Sz vertex in the end.
  // ppsc::operator_type Sz = 0.5 * h.m; // Sz is 1/2 * magnetization
  // gf_verts.push_back(ppsc::gf_vert_type(0, 0, Sz, Sz)); 

  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  int itermax, iter_rtime, nt, ntau, kt, iter_rt,
    iter_equil, iter_warmup, tstp, order, nomp, store_pp;

  double beta, h, errmax, dmfterr, V, omega, dmfterr_equil,
    linear_mixing, sym_tol, Bz_seed, Bz_org_val,bath_low,bath_high,mu_bath;


  int read_eq_sym, read_rt_sym, read_state_from_file, af_bethe_flag;

  bool matsubara_converged = false;
  int bath_flag, imp_sym;
  cntr::function<double> ecoup;

  std::vector<double> U,Jh,Jex,Bz,eps,mu,deltaImp;
  cntr::function<double> Jfunc;
  std::vector<double> t_vec,phi_vec;
  try {
    // -------------------------------------------------------------------
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
      find_param(argv[1], "__imp_sym=", imp_sym);
      find_param(argv[1], "__sym_tol=", sym_tol);
      find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__read_state_from_file=", read_state_from_file);
      find_param(argv[1], "__af_bethe_flag=", af_bethe_flag);
      find_param(argv[1], "__bath_flag=", bath_flag);
      find_param(argv[1], "__bath_low=", bath_low);
      find_param(argv[1], "__bath_high=", bath_high);
      find_param(argv[1], "__mu_bath=", mu_bath);
      find_param(argv[1], "__nomp=", nomp);
    }

    // -------------------------------------------------------------------
    // boson bath
      cntr::herm_matrix<double> D0(nt, ntau, 1, -1);
      cntr::herm_matrix<double> D0_cc(nt, ntau, 1, -1);
      double omega;
      std::vector<double> g;
      find_param_tvector(argv[1], "__g=", g, nt);
      find_param(argv[1], "__omega=", omega);

      cntr::smooth_box dos1(bath_low,bath_high,30);
      green_equilibrium_doped(D0,dos1,beta,h,30,30,mu_bath);
      cntr::function<double> gfunc(nt);
      for(int tstp=-1;tstp<=nt;tstp++){
        cdmatrix tmp(1,1);
        tmp(0,0)=g[tstp+1];
        gfunc.set_value(tstp,tmp);
      }

      for(int tstp=-1; tstp <= nt; tstp++) {
        D0.left_multiply(tstp, gfunc);
        D0.right_multiply(tstp, gfunc);
      }
      for(int tstp = -1; tstp <= nt; tstp++) ppsc::set_bwd_from_fwd(tstp, D0_cc, D0);      
      hid_t file_id = open_hdf5_file("bath.h5");
      hid_t group_id1 = create_group(file_id, "eBath");
      store_herm_greens_function(group_id1, D0);
      close_group(group_id1);
      close_hdf5_file(file_id);

    // -------------------------------------------------------------------
    // -- Setup pp calculator
    typedef ppsc::hilbert_spaces::two_band_fermi_tJ hilbert_space_type;
    typedef ppsc::hamiltonians::two_band_tJ<hilbert_space_type> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
        
    hilbert_space_type hilbert_space;
    hilbert_space.init();
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    { // -- Setup local Hamiltonian
      find_param_tvector(argv[1], "__t=", t_vec, nt);
      find_param_tvector(argv[1], "__phi=", phi_vec, nt);
      find_param_tvector(argv[1], "__U=", U, nt);     // Hubbard U
      find_param_tvector(argv[1],"__Jh=",Jh,nt);        // Hund J
      find_param_tvector(argv[1],"__Jex=",Jex,nt);        // exchange coupling strength J
      find_param_tvector(argv[1], "__delta=", deltaImp, nt);
      find_param_tvector(argv[1], "__Bz=", Bz, nt);
      find_param_tvector(argv[1], "__eps=", eps, nt);
      find_param_tvector(argv[1],"__mu=",mu,nt);      // chemical potential
      // assume that file is actually electric field -> integrate
      phi_vec[0]=0.;
      for (int i=1;i<nt+1;i++) {
        phi_vec[i]*=h;
        phi_vec[i]+=phi_vec[i-1]; 
      }
      find_param(argv[1], "__Bz_seed=", Bz_seed);

      imp.hamiltonian.mu = mu[0];
      imp.hamiltonian.U = U;
      imp.hamiltonian.J = Jh;
      imp.hamiltonian.Bz = Bz;
      imp.hamiltonian.delta = deltaImp;
      Bz_org_val = imp.hamiltonian.Bz[0];
      imp.update_hamiltonian();
      imp.print_hamiltonian(-1);
      // -- Hamiltonian Hermicity check
      ppsc::operator_type Htemp;
      imp.hamiltonian.get_hamiltonian(-1, Htemp);
      ppsc::mam::dynamic_matrix_type H = Htemp.to_dense();
      float herm_diff = (H - H.transpose()).cwiseAbs().maxCoeff();
      std::cout << "check hermiticity " << herm_diff << "\n";
      Jfunc=cntr::function<double>(nt);
      for(int tstp=-1;tstp<=nt;tstp++){
        cdmatrix tmp(1,1);
        tmp(0,0)=Jex[tstp+1];
        Jfunc.set_value(tstp,tmp);
      }
    }

    // -------------------------------------------------------------------
    // -- init single particle gf, hyb, and hopping

    single_particle_greensfunction_type Gloc(nt, ntau);
    hybridization_function_type Delta(nt, ntau);
    cntr::function<double> tfunc(nt);
    // set tfunc from t_vec
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc[tstp] = t_vec[tstp+1];
    }

    cntr::function<double> tfunc1sin(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc1sin[tstp] = t_vec[tstp+1]*sin(phi_vec[tstp+1]);
    }
    cntr::function<double> tfunc2sin(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc2sin[tstp] = t_vec[tstp+1]*sin(phi_vec[tstp+1]);
    }

    cntr::function<double> tfunc1cos(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc1cos[tstp] = t_vec[tstp+1]*cos(phi_vec[tstp+1]);
    }
    cntr::function<double> tfunc2cos(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc2cos[tstp] = t_vec[tstp+1]*cos(phi_vec[tstp+1]);
    }
    // -------------------------------------------------------------------
    
    if(read_state_from_file) {

      std::string filename = "data_ppsc.h5";
      std::cout << "--> Reading state from file: " << filename << std::endl;

      Gloc.load(filename);
      Delta.load(filename);
      imp.load(filename);

      if(imp_sym && read_eq_sym) {
	std::string filename = "sym_eq.txt";
	imp.read_symmetries(filename);
      }

      imp.update_density_matrix(-1);
      imp.hamiltonian.update_exp_vals(-1, imp.rho);
      std::cout << "pp_mu = " << imp.pp_mu << std::endl;
    }
    // -------------------------------------------------------------------
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
        imp.update_diagrams(pp_ints, gf_verts);

	      if(iter_equil == 1 && read_eq_sym && !imp.has_symmetries()) {
	         std::string filename = "sym_eq.txt";
	         imp.read_symmetries(filename);
	      }
	      // -- Seeding of spontaneous symmetry breaking
	      int n_seed_iterations = 4;
        // std::cout << "Density matrix " << Gloc.o1up.density_matrix(-1) << " " << Gloc.o2up.density_matrix(-1) << ""  << Gloc.o1do.density_matrix(-1) <<  " " <<  Gloc.o2do.density_matrix(-1) << std::endl;
        double magn=std::real(Gloc.o1up.density_matrix(-1) + Gloc.o2up.density_matrix(-1) - Gloc.o1do.density_matrix(-1) -  Gloc.o2do.density_matrix(-1));
	      if(iter_equil < n_seed_iterations) {
	       std::cout << "--> Applying Bz_seed = " << Bz_seed << std::endl;
	       imp.hamiltonian.Bz[0] = Bz_seed+magn*Jex[0]*0.5;
          imp.update_hamiltonian(-1); // update only tstp -1

	      } else if(iter_equil == n_seed_iterations) {
	        imp.hamiltonian.Bz[0] = Bz[0]+magn*Jex[0]*0.5; // Restore equil value	  
          imp.update_hamiltonian(-1);
	      }

        // -- Solve pseudo particle problem
	      imp.pp_step(-1);

        // -- get spgf
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);
	      Gloc.update(-1, gf_tstps, linear_mixing);


        // -- Check error
        dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc.o1up);
        gtmp.set_timestep(-1, Gloc.o1up);
        Delta.update(-1, Gloc, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos,
                     D0, af_bethe_flag, bath_flag);
	

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

    // -------------------------------------------------------------------

    if(imp_sym && !imp.has_symmetries()) {
      imp.symmetry_reduction(-1);
      std::string filename = "sym_eq.txt";
      imp.write_symmetries(filename);
    }
	
    if(imp_sym && nt > 0) imp.clear_symmetries();

    if(imp_sym && read_rt_sym) {
      std::string filename = "sym_rt.txt";
      imp.read_symmetries(filename);
    }
    // -------------------------------------------------------------------
    //  START ... same iteration
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      
      for (iter_warmup = 1; iter_warmup <= itermax; iter_warmup++) {

        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);
        
        imp.update_diagrams(pp_ints, gf_verts);

        double magn=std::real(Gloc.o1up.density_matrix(-1) + Gloc.o2up.density_matrix(-1) - Gloc.o1do.density_matrix(-1) -  Gloc.o2do.density_matrix(-1));
        for (int n = 0; n <= kt; n++) {
          imp.hamiltonian.Bz[n+1] = Bz[n+1] + magn*Jex[0]*0.5; // Using eq. values should improve the convergence
        }
        imp.update_hamiltonian();
        imp.print_hamiltonian(kt);
        imp.pp_step(kt);

        for (int n = 0; n <= kt; n++) {
          ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
	  Gloc.update(n, gf_tstps);

//	  Delta.update(n, Gloc, tfunc1sin, tfunc2sin, tfunc3sin, tfunc1cos, tfunc2cos, tfunc3cos,
//		       D0, af_bethe_flag, bath_flag);
          Delta.update(n, Gloc, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, 
                       D0, af_bethe_flag, bath_flag);
        }
        
        dmfterr = cntr::distance_norm2(kt, gtmp, Gloc.o1up);
        gtmp.set_timestep(kt, Gloc.o1up);
        
        cout << "WARMUP: iter:  " << iter_warmup
	     << " err: " << dmfterr << endl;
        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }
      }
    }

    // -------------------------------------------------------------------
    // -- Reduce number of diagrams using symmetries

    if(imp_sym && matsubara_converged && nt > 0 && !imp.has_symmetries()) {
      imp.symmetry_reduction(kt);
      std::string filename = "sym_rt.txt";
      imp.write_symmetries(filename);
    }
    
    // -------------------------------------------------------------------
    //  REALTIME: ONLY FOR NT>0
    if(matsubara_converged)
    for (tstp = kt + 1; tstp <= nt; tstp++) {

      imp.extrapolate_timestep(tstp - 1);
      
      for (iter_rt = 1; iter_rt <= iter_rtime; iter_rt++) {

        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);

        imp.update_diagrams(pp_ints, gf_verts);
        double magn;
        if(iter_rt==1){
          magn=std::real(Gloc.o1up.density_matrix(tstp-1) + Gloc.o2up.density_matrix(tstp-1) - Gloc.o1do.density_matrix(tstp-1) -  Gloc.o2do.density_matrix(tstp-1));  
        }else{
          magn=std::real(Gloc.o1up.density_matrix(tstp) + Gloc.o2up.density_matrix(tstp) - Gloc.o1do.density_matrix(tstp) -  Gloc.o2do.density_matrix(tstp));
        }
        
        imp.hamiltonian.Bz[tstp+1] = Bz[tstp+1] +  magn*Jex[tstp+1]*0.5;
        imp.update_hamiltonian(tstp);
        imp.print_hamiltonian(tstp);

        imp.pp_step(tstp);

        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
	      Gloc.update(tstp, gf_tstps);

        Delta.update(tstp, Gloc, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, 
                     D0, af_bethe_flag, bath_flag);

	
      }
    }

    // -------------------------------------------------------------------
    // Kinetic energy


// ? move symmetry analysis here ?


    //imp.calculate_pp_interaction_energy();

//    std::vector<double> Ekin_o1up = ppsc::get_kinetic_energy(Gloc.o1up, Delta.o1up, beta, h, kt);
//    std::vector<double> Ekin_o1do = ppsc::get_kinetic_energy(Gloc.o1do, Delta.o1do, beta, h, kt);
//    std::vector<double> Ekin_o2up = ppsc::get_kinetic_energy(Gloc.o2up, Delta.o2up, beta, h, kt);
//    std::vector<double> Ekin_o2do = ppsc::get_kinetic_energy(Gloc.o2do, Delta.o2do, beta, h, kt);
//    std::vector<double> Ekin_o3up = ppsc::get_kinetic_energy(Gloc.o3up, Delta.o3up, beta, h, kt);
//    std::vector<double> Ekin_o3do = ppsc::get_kinetic_energy(Gloc.o3do, Delta.o3do, beta, h, kt);

// AFM:
    std::vector<double> Ekin_o1up = ppsc::get_kinetic_energy(Gloc.o1up, Delta.o1do, beta, h, kt);
    std::vector<double> Ekin_o1do = ppsc::get_kinetic_energy(Gloc.o1do, Delta.o1up, beta, h, kt);
    std::vector<double> Ekin_o2up = ppsc::get_kinetic_energy(Gloc.o2up, Delta.o2do, beta, h, kt);
    std::vector<double> Ekin_o2do = ppsc::get_kinetic_energy(Gloc.o2do, Delta.o2up, beta, h, kt);


    std::vector<double> Ekin;
    Ekin.resize(Ekin_o1up.size());
    
    for( auto t : ppsc::range(0, Ekin.size()) ) {
      Ekin[t] = Ekin_o1up[t] + Ekin_o1do[t]
	      + Ekin_o2up[t] + Ekin_o2do[t];
//	      + Ekin_o3up[t] + Ekin_o3do[t];

      std::cout << "ekin, epot " << t << " "
		<< Ekin[t] << " " << imp.hamiltonian.Eint_exp[t] << "\n";
    } 

    // -------------------------------------------------------------------
    // OBSERVABLES
    {
      
      // Gloc.o1up.print_to_file("Gloc.o1up.out");
      // Gloc.o2up.print_to_file("Gloc.o2up.out");

      // Gloc.o1do.print_to_file("Gloc.o1do.out");
      // Gloc.o2do.print_to_file("Gloc.o2do.out");

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

