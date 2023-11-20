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

#include "./ppsc/hilbert_spaces/two_band_fermi_diag.hpp"
#include "./ppsc/hilbert_spaces/two_band_fermi_densdens.hpp"
#include "./ppsc/hamiltonians/two_band_hubbard_ising_kanamori_scsusc.hpp"

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

/*
      gloc_mix.clear();
      this->o3up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[4], 1.0 - linear_mixing);
      this->o3up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o3do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[5], 1.0 - linear_mixing);
      this->o3do.set_timestep(-1, gloc_mix);
*/


// HERE WE IMPOSE BY HAND THE SYMMETRIES FOR T8>0, I.E. ORB2=ORB3 AND SPINUP=SPINDOWN
/*
      gloc_mix.clear();
      this->o1up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o1do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o1up.set_timestep(-1, gloc_mix);
      this->o1do.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o2up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.25);
      this->o2do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.25);
      this->o3up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.25);
      this->o3do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.25);
      this->o2up.set_timestep(-1, gloc_mix);
      this->o2do.set_timestep(-1, gloc_mix);
      this->o3up.set_timestep(-1, gloc_mix);
      this->o3do.set_timestep(-1, gloc_mix);
*/
// END T8>0

// IMPOSE SPIN SYMMETRY

/*
      gloc_mix.clear();
      this->o1up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o1do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o1up.set_timestep(-1, gloc_mix);
      this->o1do.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o2up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o2do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o2up.set_timestep(-1, gloc_mix);
      this->o2do.set_timestep(-1, gloc_mix);
*/

/*
      gloc_mix.clear();
      this->o3up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o3do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o3up.set_timestep(-1, gloc_mix);
      this->o3do.set_timestep(-1, gloc_mix);
*/

// END SPIN SYMMETRY

    } else {

      this->o1up.set_timestep(tstp, gf_tstps[0]);
      this->o1do.set_timestep(tstp, gf_tstps[1]);
      this->o2up.set_timestep(tstp, gf_tstps[2]);
      this->o2do.set_timestep(tstp, gf_tstps[3]);
//      this->o3up.set_timestep(tstp, gf_tstps[4]);
//      this->o3do.set_timestep(tstp, gf_tstps[5]);
      
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

/*    
    group_id = create_group(file_id, "g3u");
    store_herm_greens_function(group_id, this->o3up);
    close_group(group_id);

    group_id = create_group(file_id, "g3d");
    store_herm_greens_function(group_id, this->o3do);
    close_group(group_id);
*/

  }

  void load(std::string filename) {
    this->o1up.read_from_hdf5(filename.c_str(), "g1u");
    this->o1do.read_from_hdf5(filename.c_str(), "g1d");
    this->o2up.read_from_hdf5(filename.c_str(), "g2u");
    this->o2do.read_from_hdf5(filename.c_str(), "g2d");
//    this->o3up.read_from_hdf5(filename.c_str(), "g3u");
//    this->o3do.read_from_hdf5(filename.c_str(), "g3d");
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
//    o3up(nt, ntau, 1, -1),
//    o3do(nt, ntau, 1, -1),

    // -- bwd hybridizations
    o1up_cc(nt, ntau, 1, -1),
    o1do_cc(nt, ntau, 1, -1),
    o2up_cc(nt, ntau, 1, -1),
    o2do_cc(nt, ntau, 1, -1)
//    o3up_cc(nt, ntau, 1, -1),
//    o3do_cc(nt, ntau, 1, -1)
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

/*
    // -- Update Delta.o3up
    if(af_bethe_flag == 1) _component_update(tstp, this->o3up, this->o3up_cc, Gloc.o3do, tfunc3sin, tfunc3cos, D0, bath_flag);
    else                   _component_update(tstp, this->o3up, this->o3up_cc, Gloc.o3up, tfunc3sin, tfunc3cos, D0, bath_flag);

    // -- Update Delta.o3do
    if(af_bethe_flag == 1) _component_update(tstp, this->o3do, this->o3do_cc, Gloc.o3up, tfunc3sin, tfunc3cos, D0, bath_flag);
    else                   _component_update(tstp, this->o3do, this->o3do_cc, Gloc.o3do, tfunc3sin, tfunc3cos, D0, bath_flag);
*/    
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
      Bubble2(tstp, tmp, Gloc_comp, D0);
      Delta_comp.incr_timestep(tstp, tmp);
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
/*
    group_id = create_group(file_id, "d3u");
    store_herm_greens_function(group_id, this->o3up);
    close_group(group_id);

    group_id = create_group(file_id, "d3d");
    store_herm_greens_function(group_id, this->o3do);
    close_group(group_id);
*/
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
/*    
    group_id = create_group(file_id, "d3ucc");
    store_herm_greens_function(group_id, this->o3up_cc);
    close_group(group_id);

    group_id = create_group(file_id, "d3dcc");
    store_herm_greens_function(group_id, this->o3do_cc);
    close_group(group_id);
*/
  }

  void load(std::string filename) {
    this->o1up.read_from_hdf5(filename.c_str(), "d1u");
    this->o1do.read_from_hdf5(filename.c_str(), "d1d");
    this->o2up.read_from_hdf5(filename.c_str(), "d2u");
    this->o2do.read_from_hdf5(filename.c_str(), "d2d");
//    this->o3up.read_from_hdf5(filename.c_str(), "d3u");
//    this->o3do.read_from_hdf5(filename.c_str(), "d3d");

    this->o1up_cc.read_from_hdf5(filename.c_str(), "d1ucc");
    this->o1do_cc.read_from_hdf5(filename.c_str(), "d1dcc");
    this->o2up_cc.read_from_hdf5(filename.c_str(), "d2ucc");
    this->o2do_cc.read_from_hdf5(filename.c_str(), "d2dcc");
//    this->o3up_cc.read_from_hdf5(filename.c_str(), "d3ucc");
//    this->o3do_cc.read_from_hdf5(filename.c_str(), "d3dcc");
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

/*  
  // orbital no. 3

  pp_ints.push_back(ppsc::pp_int_type(Delta.o3up,    h.c3uc,  h.c3ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o3up_cc, h.c3ua,  h.c3uc,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.o3do,    h.c3dc,  h.c3da,  fermion, fwd)); // spin do fwd 
  pp_ints.push_back(ppsc::pp_int_type(Delta.o3do_cc, h.c3da,  h.c3dc,  fermion, bwd)); // spin do bwd 
*/

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

/*
  gf_verts.push_back(ppsc::gf_vert_type(2, 2, h.c3ua, h.c3uc)); // spin up orb 3
  gf_verts.push_back(ppsc::gf_vert_type(2, 2, h.c3da, h.c3dc)); // spin do orb 3
*/

// Add the Sz Sz vertex in the end.
  ppsc::operator_type Sz = 0.5 * h.m; // Sz is 1/2 * magnetization
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, Sz, Sz)); 

  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  int itermax, iter_rtime, nt, ntau, kt, iter_rt,
    iter_equil, iter_warmup, tstp, order, nomp, store_pp;

  double beta, h, errmax, dmfterr, V, mu, omega, dmfterr_equil,
    linear_mixing, sym_tol, Bz_seed, Bz_org_val, sosm_seed;

  int read_eq_sym, read_rt_sym, read_state_from_file, af_bethe_flag;

  bool matsubara_converged = false;
  int bath_flag, sym_flag;

//  std::vector<double> t_vec, t_vec1, t_vec2, t_vec3, phi_vec1, phi_vec2, phi_vec3;
  std::vector<double> t_vec, t_vec1, t_vec2, phi_vec1, phi_vec2;

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


      find_param(argv[1], "__sym_flag=", sym_flag);
      find_param(argv[1], "__sym_tol=", sym_tol);
      find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__read_state_from_file=", read_state_from_file);

      find_param(argv[1], "__af_bethe_flag=", af_bethe_flag);
      find_param(argv[1], "__bath_flag=", bath_flag);
     
      find_param(argv[1], "__nomp=", nomp);
    }

    // -------------------------------------------------------------------
    // boson bath

      cntr::herm_matrix<double> D0(nt, ntau, 1, +1);
      cntr::herm_matrix<double> D0_cc(nt, ntau, 1, +1);
      double omega;
      std::vector<double> g;
      find_param_tvector(argv[1], "__g=", g, nt);
      find_param(argv[1], "__omega=", omega);

      ppsc::boson_utils::green_from_eps_phonon(beta, D0, omega, h);
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
    
      for(int tstp = -1; tstp <= nt; tstp++)
        ppsc::set_bwd_from_fwd(tstp, D0_cc, D0);

    // -------------------------------------------------------------------
    // -- Setup pp calculator

    typedef ppsc::hilbert_spaces::two_band_fermi_diag hilbert_space_type;
    //typedef ppsc::hilbert_spaces::two_band_fermi_densdens hilbert_space_type;
    //typedef ppsc::hilbert_spaces::three_band_fermi_diag hilbert_space_type;
    //typedef ppsc::hilbert_spaces::three_band_fermi_densdens hilbert_space_type;

    //typedef ppsc::hamiltonians::three_band_hubbard<hilbert_space_type,
    typedef ppsc::hamiltonians::two_band_hubbard_ising_kanamori<hilbert_space_type,
      ppsc::hamiltonians::interaction_type::ising> hamiltonian_type;
    //  ppsc::hamiltonians::interaction_type::kanamori> hamiltonian_type;   
 
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    hilbert_space.init();
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    
    { // -- Setup local Hamiltonian
      find_param(argv[1], "__mu=", imp.hamiltonian.mu);

      find_param_tvector(argv[1], "__dU1=", imp.hamiltonian.dU1, nt);
      find_param_tvector(argv[1], "__dU2=", imp.hamiltonian.dU2, nt);
//      find_param_tvector(argv[1], "__dU3=", imp.hamiltonian.dU3, nt);

      find_param_tvector(argv[1], "__U=", imp.hamiltonian.U, nt);
      find_param(argv[1], "__U0=", imp.hamiltonian.U[0]);

      find_param_tvector(argv[1], "__J=", imp.hamiltonian.J, nt);
      find_param(argv[1], "__J0=", imp.hamiltonian.J[0]);
      
      find_param_tvector(argv[1], "__delta=", imp.hamiltonian.delta, nt);
      find_param(argv[1], "__delta0=", imp.hamiltonian.delta[0]);

      find_param_tvector(argv[1], "__t=", t_vec, nt);
      find_param(argv[1], "__t0=", t_vec[0]);

      find_param_tvector(argv[1], "__t1=", t_vec1, nt);
      find_param_tvector(argv[1], "__t2=", t_vec2, nt);
//      find_param_tvector(argv[1], "__t3=", t_vec3, nt);

      find_param_tvector(argv[1], "__phi1=", phi_vec1, nt);
      find_param_tvector(argv[1], "__phi2=", phi_vec2, nt);
//      find_param_tvector(argv[1], "__phi3=", phi_vec3, nt);

// assume that file is actually electric field -> integrate
phi_vec1[0]=0.;
for (int i=1;i<nt+1;i++) {
  phi_vec1[i]*=h;
  phi_vec1[i]+=phi_vec1[i-1]; 
}
phi_vec2[0]=0.;
for (int i=1;i<nt+1;i++) {
  phi_vec2[i]*=h;
  phi_vec2[i]+=phi_vec2[i-1];
}
/*
phi_vec3[0]=0.;
for (int i=1;i<nt+1;i++) {
  phi_vec3[i]*=h;
  phi_vec3[i]+=phi_vec3[i-1];
}
*/

      find_param_tvector(argv[1], "__dE1=", imp.hamiltonian.dE1, nt);
      find_param_tvector(argv[1], "__dE2=", imp.hamiltonian.dE2, nt);
//      find_param_tvector(argv[1], "__dE3=", imp.hamiltonian.dE3, nt);

      find_param_tvector(argv[1], "__Bz=", imp.hamiltonian.Bz, nt);
      find_param(argv[1], "__Bz0=", imp.hamiltonian.Bz[0]);
      find_param(argv[1], "__Bz_seed=", Bz_seed);
      Bz_org_val = imp.hamiltonian.Bz[0];

      find_param(argv[1], "__sosm_seed=", sosm_seed);
 
      imp.update_hamiltonian();

      // -- Hamiltonian Hermicity check
      ppsc::operator_type Htemp;
      imp.hamiltonian.get_hamiltonian(-1, Htemp);
      ppsc::mam::dynamic_matrix_type H = Htemp.to_dense();
      float herm_diff = (H - H.transpose()).cwiseAbs().maxCoeff();
      std::cout << "check hermiticity " << herm_diff << "\n";

    }
    // -------------------------------------------------------------------
    // -- init single particle gf, hyb, and hopping

    single_particle_greensfunction_type Gloc(nt, ntau);
    hybridization_function_type Delta(nt, ntau);


int boson=+1, dim=1;
cntr::herm_matrix<double> SzSz(nt, ntau, dim, boson);
    

    cntr::function<double> tfunc(nt);
    // set tfunc from t_vec
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc[tstp] = t_vec[tstp+1];
    }

    cntr::function<double> tfunc1sin(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc1sin[tstp] = t_vec1[tstp+1]*sin(phi_vec1[tstp+1]);
    }
    cntr::function<double> tfunc2sin(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc2sin[tstp] = t_vec2[tstp+1]*sin(phi_vec2[tstp+1]);
    }
/*
    cntr::function<double> tfunc3sin(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc3sin[tstp] = t_vec3[tstp+1]*sin(phi_vec3[tstp+1]);
    }
*/

    cntr::function<double> tfunc1cos(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc1cos[tstp] = t_vec1[tstp+1]*cos(phi_vec1[tstp+1]);
    }
    cntr::function<double> tfunc2cos(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc2cos[tstp] = t_vec2[tstp+1]*cos(phi_vec2[tstp+1]);
    }
/*
    cntr::function<double> tfunc3cos(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc3cos[tstp] = t_vec3[tstp+1]*cos(phi_vec3[tstp+1]);
    }
*/

    // -------------------------------------------------------------------
    
    if(read_state_from_file) {

      std::string filename = "data_ppsc.h5";
      std::cout << "--> Reading state from file: " << filename << std::endl;

      Gloc.load(filename);
      Delta.load(filename);
      imp.load(filename);

      if(sym_flag && read_eq_sym) {
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
	if(iter_equil < n_seed_iterations) {
	  std::cout << "--> Applying Bz_seed = " << Bz_seed << std::endl;
	  imp.hamiltonian.Bz[0] = Bz_seed;

          std::cout << "--> Applying sosm_seed = "
		    << sosm_seed << std::endl;
          imp.hamiltonian.sosm_seed = sosm_seed;
	  
          imp.update_hamiltonian(-1); // update only tstp -1

	} else if(iter_equil == n_seed_iterations) {
	  imp.hamiltonian.Bz[0] = Bz_org_val; // Restore equil value	  
          imp.hamiltonian.sosm_seed = 0.; // set to zero
          imp.update_hamiltonian(-1);
	}

        // -- Solve pseudo particle problem
	imp.pp_step(-1);

        // -- get spgf
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);
	Gloc.update(-1, gf_tstps, linear_mixing);

SzSz.set_timestep(-1, gf_tstps[4]); 

        // -- Check error
        dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc.o1up);
        gtmp.set_timestep(-1, Gloc.o1up);

        // -- Update Hybridization
//	Delta.update(-1, Gloc, tfunc1sin, tfunc2sin, tfunc3sin, tfunc1cos, tfunc2cos, tfunc3cos,
//		     D0, af_bethe_flag, bath_flag);
        Delta.update(-1, Gloc, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos,
                     D0, af_bethe_flag, bath_flag);
	

        cout << "iter:  " << iter_equil
	     << " err: " << dmfterr_equil << endl;

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

    if(sym_flag && !imp.has_symmetries()) {
      imp.symmetry_reduction(-1);
      std::string filename = "sym_eq.txt";
      imp.write_symmetries(filename);
    }
	
    if(sym_flag && nt > 0) imp.clear_symmetries();

    if(sym_flag && read_rt_sym) {
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
        imp.pp_step(kt);

        for (int n = 0; n <= kt; n++) {
          ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
	  Gloc.update(n, gf_tstps);

SzSz.set_timestep(n, gf_tstps[4]); 

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

    if(sym_flag && matsubara_converged && nt > 0 && !imp.has_symmetries()) {
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
        imp.pp_step(tstp);

        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
	Gloc.update(tstp, gf_tstps);

SzSz.set_timestep(tstp, gf_tstps[4]); 

//	Delta.update(tstp, Gloc, tfunc1sin, tfunc2sin, tfunc3sin, tfunc1cos, tfunc2cos, tfunc3cos,
//		     D0, af_bethe_flag, bath_flag);
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
      
      Gloc.o1up.print_to_file("Gloc.o1up.out");
      Gloc.o2up.print_to_file("Gloc.o2up.out");
//      Gloc.o3up.print_to_file("Gloc.o3up.out");

      Gloc.o1do.print_to_file("Gloc.o1do.out");
      Gloc.o2do.print_to_file("Gloc.o2do.out");
//      Gloc.o3do.print_to_file("Gloc.o3do.out");

SzSz.print_to_file("SzSz.out"); 

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
//      store_real_data_to_hid(group_id, "Ekin_o3up", Ekin_o3up.data(), Ekin_o3up.size());
//      store_real_data_to_hid(group_id, "Ekin_o3do", Ekin_o3do.data(), Ekin_o3do.size());

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

