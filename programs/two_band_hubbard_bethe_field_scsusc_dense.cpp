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

#include "./ppsc/hilbert_spaces/n_band_fermi_dense.hpp"
#include "./ppsc/hamiltonians/two_band_hubbard_ising_kanamori_scsusc_dense.hpp"

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
    o2up(nt, ntau, 1, -1), o2do(nt, ntau, 1, -1),
    otaup(nt, ntau, 1, -1), otado(nt, ntau, 1, -1),
    otcup(nt, ntau, 1, -1), otcdo(nt, ntau, 1, -1), 
    otaupdo(nt, ntau, 1, -1), otadoup(nt, ntau, 1, -1),
    otcupdo(nt, ntau, 1, -1), otcdoup(nt, ntau, 1, -1) {}

  void update(int tstp, ppsc::gf_tstps_type & gf_tstps, double linear_mixing=0.0) {

    if(tstp == -1) {
      ppsc::gf_tstp_type gloc_old(-1, ntau, 1);
      ppsc::gf_tstp_type gloc_mix(-1, ntau, 1);

      auto mixer = [&](cntr::herm_matrix<double> & gf_component, int idx_component) {
	gloc_mix.clear();
	gf_component.get_timestep(-1, gloc_old);
	gloc_mix.incr(gloc_old, linear_mixing);
	gloc_mix.incr(gf_tstps[idx_component], 1.0 - linear_mixing);
	gf_component.set_timestep(-1, gloc_mix);
      };

      // anomalous components
      // 4: cc up (otaup)
      // 5: cc do (otado)
      // 6: cdagcdag up (otcup)
      // 7: cdagcdag do (otcdo)
      // 8: cc 2up1do (otaupdo)
      // 9: cc 2do1up (otadoup)
      // 10: cdagcdag 1up2do (otcupdo)
      // 11: cdagcdag 1do2up (otcdoup)
 
      mixer(this->o1up, 0);
      mixer(this->o1do, 1);
      mixer(this->o2up, 2);
      mixer(this->o2do, 3);      

      mixer(this->otaup, 4);
      mixer(this->otado, 5);
      mixer(this->otcup, 6);
      mixer(this->otcdo, 7);

      mixer(this->otaupdo, 8);
      mixer(this->otadoup, 9);
      mixer(this->otcupdo, 10);
      mixer(this->otcdoup, 11);

      // IMPOSE SYMMETRIES

      auto average = [&](cntr::herm_matrix<double> & A, cntr::herm_matrix<double> & B) {
	gloc_mix.clear();
	A.get_timestep(-1, gloc_old); gloc_mix.incr(gloc_old, 0.5);
	B.get_timestep(-1, gloc_old); gloc_mix.incr(gloc_old, 0.5);
	A.set_timestep(-1, gloc_mix);
	B.set_timestep(-1, gloc_mix);
      };

      auto anomalous_average = [&](cntr::herm_matrix<double> & A, cntr::herm_matrix<double> & B) {
        gloc_mix.clear();
        A.get_timestep(-1, gloc_old); gloc_mix.incr(gloc_old, 0.5);
        B.get_timestep(-1, gloc_old); gloc_mix.incr(gloc_old, -0.5);
        A.set_timestep(-1, gloc_mix);
        gloc_mix.clear();
        A.get_timestep(-1, gloc_old); gloc_mix.incr(gloc_old, -0.5);
        B.get_timestep(-1, gloc_old); gloc_mix.incr(gloc_old, 0.5);
        B.set_timestep(-1, gloc_mix);
      };

/*
      // SPIN SYMMETRY 

      // normal components
      average(this->o1up, this->o1do);
      average(this->o2up, this->o2do);

      // anomalous components
      average(this->otaup, this->otado);
      average(this->otcup, this->otcdo);

      // END SPIN SYMMETRY
*/     

/*
      // ORBITAL SYMMETRY
	
      // normal components
      average(this->o1up, this->o2up);
      average(this->o1do, this->o2do);

      // anomalous components
      anomalous_average(this->otaup, this->otcup);
      anomalous_average(this->otado, this->otcdo);
      
      // END ORBITAL SYMMETRY
*/

    } else {

      this->o1up.set_timestep(tstp, gf_tstps[0]);
      this->o1do.set_timestep(tstp, gf_tstps[1]);
      this->o2up.set_timestep(tstp, gf_tstps[2]);
      this->o2do.set_timestep(tstp, gf_tstps[3]);

/*
      // IMPOSE ORBITAL SYMMETRIES

      ppsc::gf_tstp_type gloc_old(tstp, ntau, 1);
      ppsc::gf_tstp_type gloc_mix(tstp, ntau, 1);

      auto average = [&](cntr::herm_matrix<double> & A, cntr::herm_matrix<double> & B) {
        gloc_mix.clear();
        A.get_timestep(tstp, gloc_old); gloc_mix.incr(gloc_old, 0.5);
        B.get_timestep(tstp, gloc_old); gloc_mix.incr(gloc_old, 0.5);
        A.set_timestep(tstp, gloc_mix);
        B.set_timestep(tstp, gloc_mix);
      };

      average(this->o1up, this->o2up);
      average(this->o1do, this->o2do);

      // END ORBITAL SYMMETRIES
*/
      // anomalous components
      this->otaup.set_timestep(tstp, gf_tstps[4]);
      this->otado.set_timestep(tstp, gf_tstps[5]);
      this->otcup.set_timestep(tstp, gf_tstps[6]);
      this->otcdo.set_timestep(tstp, gf_tstps[7]);

      this->otaupdo.set_timestep(tstp, gf_tstps[8]);
      this->otadoup.set_timestep(tstp, gf_tstps[9]);
      this->otcupdo.set_timestep(tstp, gf_tstps[10]);
      this->otcdoup.set_timestep(tstp, gf_tstps[11]);

    }
  }

  void store(hid_t file_id) {
    hid_t group_id;

    auto store = [&](cntr::herm_matrix<double> & component, std::string name) {
      group_id = create_group(file_id, name);
      store_herm_greens_function(group_id, component);
      close_group(group_id);
    };

    store(this->o1up, "g1u");
    store(this->o1do, "g1d");
    store(this->o2up, "g2u");
    store(this->o2do, "g2d");

    store(this->otaup, "g12u");
    store(this->otado, "g12d");
    store(this->otcup, "g21u");
    store(this->otcdo, "g21d");
  }

  void load(std::string filename) {
    this->o1up.read_from_hdf5(filename.c_str(), "g1u");
    this->o1do.read_from_hdf5(filename.c_str(), "g1d");
    this->o2up.read_from_hdf5(filename.c_str(), "g2u");
    this->o2do.read_from_hdf5(filename.c_str(), "g2d");

    this->otaup.read_from_hdf5(filename.c_str(), "g12u");
    this->otado.read_from_hdf5(filename.c_str(), "g12d");
    this->otcup.read_from_hdf5(filename.c_str(), "g21u");
    this->otcdo.read_from_hdf5(filename.c_str(), "g21d");
  }

  int nt, ntau;
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do; 
  cntr::herm_matrix<double> otcup, otcdo, otaup, otado; // otaup=<c1up c2up>, otcup=<cdag2up cdag1up> // [CHECK IF IMPLEMENTATION IS CONSISTENT WITH THIS! PROBABLY NOT (RATHER otcup=<cdag1up cdag2up>)]
  cntr::herm_matrix<double> otcupdo, otcdoup, otaupdo, otadoup; // otaup=<c2do c1up>, otcup=<cdag1up cdag2do>    

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
    o2do_cc(nt, ntau, 1, -1),

    // anomalous components
    otaup(nt, ntau, 1, -1),
    otado(nt, ntau, 1, -1),
    otcup(nt, ntau, 1, -1),
    otcdo(nt, ntau, 1, -1),

    otaup_cc(nt, ntau, 1, -1),
    otado_cc(nt, ntau, 1, -1),
    otcup_cc(nt, ntau, 1, -1),
    otcdo_cc(nt, ntau, 1, -1),

    // anomalous spin-flip components
    otaupdo(nt, ntau, 1, -1),
    otadoup(nt, ntau, 1, -1),
    otcupdo(nt, ntau, 1, -1),
    otcdoup(nt, ntau, 1, -1),

    otaupdo_cc(nt, ntau, 1, -1),
    otadoup_cc(nt, ntau, 1, -1),
    otcupdo_cc(nt, ntau, 1, -1),
    otcdoup_cc(nt, ntau, 1, -1)

  {}    

  void update(int tstp, single_particle_greensfunction_type & Gloc,
	      cntr::function<double> & tfunc1sin, cntr::function<double> & tfunc2sin, 
              cntr::function<double> & tfunc1cos, cntr::function<double> & tfunc2cos, 
	      cntr::herm_matrix<double> & D0, // boson bath
	      int af_bethe_flag=0, int bath_flag=0) {

    // -- Update Delta.o1up
    //if(af_bethe_flag == 1) _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);
    //else                   _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);
    if(af_bethe_flag == 1) _component_update(tstp, this->o1up, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);
    else                   _component_update(tstp, this->o1up, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);

    // -- Update Delta.o1do
    //if(af_bethe_flag == 1) _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);
    //else                   _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);
    if(af_bethe_flag == 1) _component_update(tstp, this->o1do, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);
    else                   _component_update(tstp, this->o1do, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);

    // -- Update Delta.o2up
    //if(af_bethe_flag == 1) _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o2do, tfunc2sin, tfunc2cos, D0, bath_flag);
    //else                   _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o2up, tfunc2sin, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _component_update(tstp, this->o2up, Gloc.o2do, tfunc2sin, tfunc2cos, D0, bath_flag);
    else                   _component_update(tstp, this->o2up, Gloc.o2up, tfunc2sin, tfunc2cos, D0, bath_flag);

    // -- Update Delta.o2do
    //if(af_bethe_flag == 1) _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o2up, tfunc2sin, tfunc2cos, D0, bath_flag);
    //else                   _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o2do, tfunc2sin, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _component_update(tstp, this->o2do, Gloc.o2up, tfunc2sin, tfunc2cos, D0, bath_flag);
    else                   _component_update(tstp, this->o2do, Gloc.o2do, tfunc2sin, tfunc2cos, D0, bath_flag);

    ppsc::set_bwd_from_fwd(tstp, this->o1up_cc, this->o1up);
    ppsc::set_bwd_from_fwd(tstp, this->o1do_cc, this->o1do);
    ppsc::set_bwd_from_fwd(tstp, this->o2up_cc, this->o2up);
    ppsc::set_bwd_from_fwd(tstp, this->o2do_cc, this->o2do);
    
    // anomalous components

    // HS: NB! the anti-hermitian conjugate of the anomalous-components are transposed
    // HS: i.e. [otaup_cc] is given by the anti-hermitian conjugate of [otcup]
    // HS: thus we pass the transposed _cc component below


    // -- Update Delta.otaup
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otaup, this->otcup_cc, Gloc.otcdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otaup, this->otcup_cc, Gloc.otcup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otaup, Gloc.otcdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otaup, Gloc.otcup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);


    // -- Update Delta.otado
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otado, this->otcdo_cc, Gloc.otcup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otado, this->otcdo_cc, Gloc.otcdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otado, Gloc.otcup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otado, Gloc.otcdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);

    // -- Update Delta.otcup
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcup, this->otaup_cc, Gloc.otado, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otcup, this->otaup_cc, Gloc.otaup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcup, Gloc.otado, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otcup, Gloc.otaup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);


    // -- Update Delta.otcdo
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcdo, this->otado_cc, Gloc.otaup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otcdo, this->otado_cc, Gloc.otado, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcdo, Gloc.otaup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otcdo, Gloc.otado, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
  

    ppsc::set_bwd_from_fwd_single_component(tstp, this->otaup_cc, this->otaup, this->otcup);
    ppsc::set_bwd_from_fwd_single_component(tstp, this->otado_cc, this->otado, this->otcdo);
    ppsc::set_bwd_from_fwd_single_component(tstp, this->otcup_cc, this->otcup, this->otaup);
    ppsc::set_bwd_from_fwd_single_component(tstp, this->otcdo_cc, this->otcdo, this->otado);

    
    // anomalous spin-flip components
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otaupdo, this->otcdoup_cc, Gloc.otcdoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otaupdo, this->otcdoup_cc, Gloc.otcupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otaupdo, Gloc.otcdoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otaupdo, Gloc.otcupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    // -- Update Delta.otado
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otadoup, this->otcupdo_cc, Gloc.otcupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otadoup, this->otcupdo_cc, Gloc.otcdoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otadoup, Gloc.otcupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otadoup, Gloc.otcdoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    // -- Update Delta.otcup
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcupdo, this->otadoup_cc, Gloc.otadoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otcupdo, this->otadoup_cc, Gloc.otaupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcupdo, Gloc.otadoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otcupdo, Gloc.otaupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    // -- Update Delta.otcdo
    //if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcdoup, this->otaupdo_cc, Gloc.otaupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    //else                   _anomalous_component_update(tstp, this->otcdoup, this->otaupdo_cc, Gloc.otadoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    if(af_bethe_flag == 1) _anomalous_component_update(tstp, this->otcdoup, Gloc.otaupdo, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);
    else                   _anomalous_component_update(tstp, this->otcdoup, Gloc.otadoup, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, D0, bath_flag);

    ppsc::set_bwd_from_fwd_single_component(tstp, this->otaupdo_cc, this->otaupdo, this->otcdoup);
    ppsc::set_bwd_from_fwd_single_component(tstp, this->otadoup_cc, this->otadoup, this->otcupdo);
    ppsc::set_bwd_from_fwd_single_component(tstp, this->otcupdo_cc, this->otcupdo, this->otadoup);
    ppsc::set_bwd_from_fwd_single_component(tstp, this->otcdoup_cc, this->otcdoup, this->otaupdo);
 
  }

// ATTENTION: HERE WE DO NOT USE THE FACTOR 2 IN
//
//  Delta=2vcos*G*vcos+2vsin*G*vsin,
//
// SO THAT THE FORMULA REDUCES TO THE USUAL
//
// Delta=v^2*G
//
// OF THE BETHE LATTICE WITH W=4v

  void _component_update(int tstp,
    cntr::herm_matrix<double> & Delta_comp,
    cntr::herm_matrix<double> & Gloc_comp,
    cntr::function<double> & tsin,
    cntr::function<double> & tcos,
    cntr::herm_matrix<double> & D0, // boson bath
    int bath_flag=0) {

    // this is implemented below
    // Delta = tsin * Gloc * tsin + tcos * Gloc * tcos

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
  }


  void _anomalous_component_update(int tstp,
    cntr::herm_matrix<double> & Delta_comp,
    cntr::herm_matrix<double> & Gloc_comp,
    cntr::function<double> & t1sin,
    cntr::function<double> & t2sin,
    cntr::function<double> & t1cos,
    cntr::function<double> & t2cos,
    cntr::herm_matrix<double> & D0, // boson bath
    int bath_flag=0) {

    cntr::herm_matrix_timestep<double> tmp;

    // special form for anomalous terms, with additional minus sign!
    // Delta = -( -t1sin * Gloc * t2sin + t1cos * Gloc * t2cos ) =  t1sin * Gloc * t2sin - t1cos * Gloc * t2cos
    
    Gloc_comp.get_timestep(tstp, tmp); // store to tmp

    Delta_comp.set_timestep(tstp, tmp); // set from tmp
    Delta_comp.left_multiply(tstp, t1cos); 
    Delta_comp.right_multiply(tstp, t2cos); 

    cntr::herm_matrix_timestep<double> tmpcos;
    Delta_comp.get_timestep(tstp, tmpcos); // store what we did in tmpcos = t1cos * Gloc * t2cos

    Delta_comp.set_timestep(tstp, tmp); // set (again) to tmp = Gloc_comp
    Delta_comp.left_multiply(tstp, t1sin);
    Delta_comp.right_multiply(tstp, t2sin); // Delta_comp is now = t1sin * Gloc * t2sin

    Delta_comp.incr_timestep(tstp, tmpcos, -1); // subtract t1cos * Gloc * t2cos


// no phonons for off-diagonal hybridizations
/*
    if(bath_flag) { // add boson bath
      Bubble2(tstp, tmp, Gloc_comp, D0);
      Delta_comp.incr_timestep(tstp, tmp);
    }
*/
}



  void store(hid_t file_id) {
    hid_t group_id;
    
    auto store = [&](cntr::herm_matrix<double> & component, std::string name) {
      group_id = create_group(file_id, name);
      store_herm_greens_function(group_id, component);
      close_group(group_id);
    };

    store(this->o1up, "d1u");
    store(this->o1do, "d1d");
    store(this->o2up, "d2u");
    store(this->o2do, "d2d");

    store(this->o1up_cc, "d1ucc");
    store(this->o1do_cc, "d1dcc");
    store(this->o2up_cc, "d2ucc");
    store(this->o2do_cc, "d2dcc");

    store(this->otaup, "d12u");
    store(this->otado, "d12d");
    store(this->otcup, "d21u");
    store(this->otcdo, "d21d");

    store(this->otaup_cc, "d12ucc");
    store(this->otado_cc, "d12dcc");
    store(this->otcup_cc, "d21ucc");
    store(this->otcdo_cc, "d21dcc");
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

    this->otaup.read_from_hdf5(filename.c_str(), "d12u");
    this->otado.read_from_hdf5(filename.c_str(), "d12d");
    this->otcup.read_from_hdf5(filename.c_str(), "d21u");
    this->otcdo.read_from_hdf5(filename.c_str(), "d21d");

    this->otaup_cc.read_from_hdf5(filename.c_str(), "d12ucc");
    this->otado_cc.read_from_hdf5(filename.c_str(), "d12dcc");
    this->otcup_cc.read_from_hdf5(filename.c_str(), "d21ucc");
    this->otcdo_cc.read_from_hdf5(filename.c_str(), "d21dcc");    
  }
  
  int nt, ntau;
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do; 
  cntr::herm_matrix<double> o1up_cc, o1do_cc, o2up_cc, o2do_cc; 
  cntr::herm_matrix<double> otaup, otaup_cc, otado, otado_cc, otcup, otcup_cc, otcdo, otcdo_cc;
  cntr::herm_matrix<double> otaupdo, otaupdo_cc, otadoup, otadoup_cc, otcupdo, otcupdo_cc, otcdoup, otcdoup_cc;

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

  // spin triplet pairing

  // CHECK! NEED TO SWITCH CREATION/ANNIHILATION OPERATORS UNDER _CC?

/*
  pp_ints.push_back(ppsc::pp_int_type(Delta.otaup,    h.c1ua,  h.c2ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otaup_cc, h.c2ua,  h.c1ua,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otado,    h.c1da,  h.c2da,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otado_cc, h.c2da,  h.c1da,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otcup,    h.c1uc,  h.c2uc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcup_cc, h.c2uc,  h.c1uc,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdo,    h.c1dc,  h.c2dc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdo_cc, h.c2dc,  h.c1dc,  fermion, bwd)); // spin up bwd
*/

// new convention with switched indices: Delta11, Delta22, Delta12^cdagcdag, Delta21cc

  pp_ints.push_back(ppsc::pp_int_type(Delta.otaup,    h.c2ua,  h.c1ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otaup_cc, h.c1ua,  h.c2ua,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otado,    h.c2da,  h.c1da,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otado_cc, h.c1da,  h.c2da,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otcup,    h.c1uc,  h.c2uc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcup_cc, h.c2uc,  h.c1uc,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdo,    h.c1dc,  h.c2dc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdo_cc, h.c2dc,  h.c1dc,  fermion, bwd)); // spin up bwd

// anomalous spin-flip components

/*
  pp_ints.push_back(ppsc::pp_int_type(Delta.otaupdo,    h.c1ua,  h.c2da,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otaupdo_cc, h.c2da,  h.c1ua,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otadoup,    h.c1da,  h.c2ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otadoup_cc, h.c2ua,  h.c1da,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otcupdo,    h.c2uc,  h.c1dc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcupdo_cc, h.c1dc,  h.c2uc,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdoup,    h.c2dc,  h.c1uc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdoup_cc, h.c1uc,  h.c2dc,  fermion, bwd)); // spin up bwd
*/


  pp_ints.push_back(ppsc::pp_int_type(Delta.otaupdo,    h.c2ua,  h.c1da,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otaupdo_cc, h.c1da,  h.c2ua,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otadoup,    h.c2da,  h.c1ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otadoup_cc, h.c1ua,  h.c2da,  fermion, bwd)); // spin up bwd


  pp_ints.push_back(ppsc::pp_int_type(Delta.otcupdo,    h.c1uc,  h.c2dc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcupdo_cc, h.c2dc,  h.c1uc,  fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdoup,    h.c1dc,  h.c2uc,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.otcdoup_cc, h.c2uc,  h.c1dc,  fermion, bwd)); // spin up bwd


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

// spin triplet pairs

/*
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1ua, h.c2ua)); // 4: cc_up
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1da, h.c2da)); // 5: cc_do
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1uc, h.c2uc)); // 6: cdagcdag_up
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1dc, h.c2dc)); // 7: cdagcdag_do
*/

// new convention with switched indices: G11, G22, G21^cdagcdag, G12cc
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1ua, h.c2ua)); // 4: cc_up
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1da, h.c2da)); // 5: cc_do
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c2uc, h.c1uc)); // 6: cdagcdag_up
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c2dc, h.c1dc)); // 7: cdagcdag_do

/*
// spin-flip components: G11, G22, G21^cdagcdag, G12cc
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c2ua, h.c1da)); // 8: cc 2up1do (otaupdo)
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c2da, h.c1ua)); // 9: cc 2do1up (otadoup)
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1uc, h.c2dc)); // 10: cdagcdag 1up2do (otcupdo)
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1dc, h.c2ua)); // 11: cdagcdag 1do2up (otcdoup)
*/

// spin-flip components: G11, G22, G21^cdagcdag, G12cc
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1ua, h.c2da)); // 8: cc 1up2do (otaupdo)
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c1da, h.c2ua)); // 9: cc 1do2up (otadoup)
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c2uc, h.c1dc)); // 10: cdagcdag 2up1do (otcupdo)
  gf_verts.push_back(ppsc::gf_vert_type(0, 1, h.c2dc, h.c1uc)); // 11: cdagcdag 2do1up (otcdoup)



// Add the Sz Sz vertex in the end.
  ppsc::operator_type Sz = 0.5 * h.m; // Sz is 1/2 * magnetization
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, Sz, Sz)); // 12: SzSz 

  ppsc::operator_type pair = 0.5 * (h.cdagcdag - h.cc); 
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, pair, pair)); // 13: pairpair // DOES ORBITAL 0,0 WORK FOR TRIPLET PAIRS?
//  gf_verts.push_back(ppsc::gf_vert_type(0, 0, h.cdagcdag, h.cc));

// triplet pairing
  ppsc::operator_type tp = (h.ucuc - h.uaua);
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, tp, tp)); // 14: tptp

  ppsc::operator_type t0 = 1/sqrt(2) * (h.ucdc + h.dcuc - h.uada - h.daua);
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, t0, t0)); // 15: t0t0

  ppsc::operator_type tm = (h.dcdc - h.dada);
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, tm, tm)); // 16: tmtm


// xyz

  ppsc::operator_type tz = h.tz;
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, tz, tz)); // 17: tztz

  ppsc::operator_type tx = h.tx;
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, tx, tx)); // 18: txtx

  ppsc::operator_type ty = h.ty;
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, ty, ty)); // 19: tyty


  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  int itermax, iter_rtime, nt, ntau, kt, iter_rt,
    iter_equil, iter_warmup, tstp, order, nomp, store_pp;

  double beta, h, errmax, dmfterr, V, mu, omega, dmfterr_equil,
    linear_mixing, sym_tol, Bz_seed, Bz_org_val, t0pair_seed, tppair_seed;
  double txpair_seed, typair_seed, tzpair_seed;

  int read_eq_sym, read_rt_sym, read_state_from_file, af_bethe_flag;

  bool matsubara_converged = false;
  int bath_flag, sym_flag;

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

typedef ppsc::hilbert_spaces::n_band_fermi_dense hilbert_space_type;
    //typedef ppsc::hilbert_spaces::two_band_fermi_diag hilbert_space_type;
    //typedef ppsc::hilbert_spaces::two_band_fermi_densdens hilbert_space_type;
    //typedef ppsc::hilbert_spaces::three_band_fermi_diag hilbert_space_type;
    //typedef ppsc::hilbert_spaces::three_band_fermi_densdens hilbert_space_type;

    //typedef ppsc::hamiltonians::three_band_hubbard<hilbert_space_type,
    //typedef ppsc::hamiltonians::two_band_hubbard_ising_kanamori<hilbert_space_type,
typedef ppsc::hamiltonians::two_band_hubbard_ising_kanamori_scsusc_dense<hilbert_space_type,
      ppsc::hamiltonians::interaction_type::ising> hamiltonian_type;
  //  ppsc::hamiltonians::interaction_type::kanamori> hamiltonian_type;   
 
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    // hilbert_space.init();
hilbert_space.init(2); // norb=2
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    
    { // -- Setup local Hamiltonian
      find_param(argv[1], "__mu=", imp.hamiltonian.mu);

      find_param_tvector(argv[1], "__dU1=", imp.hamiltonian.dU1, nt);
      find_param_tvector(argv[1], "__dU2=", imp.hamiltonian.dU2, nt);

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

      find_param_tvector(argv[1], "__phi1=", phi_vec1, nt);
      find_param_tvector(argv[1], "__phi2=", phi_vec2, nt);

      find_param_tvector(argv[1], "__t0pair=", imp.hamiltonian.t0pair, nt);
      find_param_tvector(argv[1], "__tppair=", imp.hamiltonian.tppair, nt);

      find_param_tvector(argv[1], "__txpair=", imp.hamiltonian.txpair, nt);
      find_param_tvector(argv[1], "__typair=", imp.hamiltonian.typair, nt);
      find_param_tvector(argv[1], "__tzpair=", imp.hamiltonian.tzpair, nt);

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

      find_param_tvector(argv[1], "__dE1=", imp.hamiltonian.dE1, nt);
      find_param_tvector(argv[1], "__dE2=", imp.hamiltonian.dE2, nt);

      find_param_tvector(argv[1], "__Bz=", imp.hamiltonian.Bz, nt);
      find_param(argv[1], "__Bz0=", imp.hamiltonian.Bz[0]);
      find_param(argv[1], "__Bz_seed=", Bz_seed);
      Bz_org_val = imp.hamiltonian.Bz[0];

      find_param(argv[1], "__t0pair_seed=", t0pair_seed);
      find_param(argv[1], "__tppair_seed=", tppair_seed); 
 
      find_param(argv[1], "__txpair_seed=", txpair_seed);
      find_param(argv[1], "__typair_seed=", typair_seed);
      find_param(argv[1], "__tzpair_seed=", tzpair_seed);


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
cntr::herm_matrix<double> pairpair(nt, ntau, dim, boson);    
cntr::herm_matrix<double> tptp(nt, ntau, dim, boson);
cntr::herm_matrix<double> t0t0(nt, ntau, dim, boson);
cntr::herm_matrix<double> tmtm(nt, ntau, dim, boson);

cntr::herm_matrix<double> tztz(nt, ntau, dim, boson);
cntr::herm_matrix<double> txtx(nt, ntau, dim, boson);
cntr::herm_matrix<double> tyty(nt, ntau, dim, boson);

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

    cntr::function<double> tfunc1cos(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc1cos[tstp] = t_vec1[tstp+1]*cos(phi_vec1[tstp+1]);
    }
    cntr::function<double> tfunc2cos(nt);
    for(int tstp=-1; tstp <= nt; tstp++) {
      tfunc2cos[tstp] = t_vec2[tstp+1]*cos(phi_vec2[tstp+1]);
    }

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

          std::cout << "--> Applying pair_seed\n";
		    // << sosm_seed << std::endl;
          imp.hamiltonian.t0pair_seed = t0pair_seed;
          imp.hamiltonian.tppair_seed = tppair_seed;

          imp.hamiltonian.txpair_seed = txpair_seed;
          imp.hamiltonian.typair_seed = typair_seed;
          imp.hamiltonian.tzpair_seed = tzpair_seed;	  

          imp.update_hamiltonian(-1); // update only tstp -1

	} else if(iter_equil == n_seed_iterations) {
	  imp.hamiltonian.Bz[0] = Bz_org_val; // Restore equil value	  
          imp.hamiltonian.t0pair_seed = 0.; // set to zero
          imp.hamiltonian.tppair_seed = 0.; // set to zero

          imp.hamiltonian.txpair_seed = 0.; // set to zero
          imp.hamiltonian.typair_seed = 0.; // set to zero
          imp.hamiltonian.tzpair_seed = 0.; // set to zero

          imp.update_hamiltonian(-1);
	}

        // -- Solve pseudo particle problem
	imp.pp_step(-1);

        // -- get spgf
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);
	Gloc.update(-1, gf_tstps, linear_mixing);


SzSz.set_timestep(-1, gf_tstps[12]);
pairpair.set_timestep(-1, gf_tstps[13]);

tptp.set_timestep(-1, gf_tstps[14]);
t0t0.set_timestep(-1, gf_tstps[15]);
tmtm.set_timestep(-1, gf_tstps[16]);

tztz.set_timestep(-1, gf_tstps[17]);
txtx.set_timestep(-1, gf_tstps[18]);
tyty.set_timestep(-1, gf_tstps[19]);

        // -- Check error
        dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc.o1up);
        gtmp.set_timestep(-1, Gloc.o1up);

        // -- Update Hybridization
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

SzSz.set_timestep(n, gf_tstps[12]);
pairpair.set_timestep(n, gf_tstps[13]);

tptp.set_timestep(n, gf_tstps[14]);
t0t0.set_timestep(n, gf_tstps[15]);
tmtm.set_timestep(n, gf_tstps[16]);

tztz.set_timestep(n, gf_tstps[17]);
txtx.set_timestep(n, gf_tstps[18]);
tyty.set_timestep(n, gf_tstps[19]);


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

SzSz.set_timestep(tstp, gf_tstps[12]);
pairpair.set_timestep(tstp, gf_tstps[3]);

tptp.set_timestep(tstp, gf_tstps[14]);
t0t0.set_timestep(tstp, gf_tstps[15]);
tmtm.set_timestep(tstp, gf_tstps[16]);

tztz.set_timestep(tstp, gf_tstps[17]);
txtx.set_timestep(tstp, gf_tstps[18]);
tyty.set_timestep(tstp, gf_tstps[19]);

        Delta.update(tstp, Gloc, tfunc1sin, tfunc2sin, tfunc1cos, tfunc2cos, 
                     D0, af_bethe_flag, bath_flag);
	
      }
    }

    // -------------------------------------------------------------------
    // Kinetic energy


// ? move symmetry analysis here ?
//
    //imp.calculate_pp_interaction_energy();

    std::vector<double> Ekin_o1up = ppsc::get_kinetic_energy(Gloc.o1up, Delta.o1up, beta, h, kt);
    std::vector<double> Ekin_o1do = ppsc::get_kinetic_energy(Gloc.o1do, Delta.o1do, beta, h, kt);
    std::vector<double> Ekin_o2up = ppsc::get_kinetic_energy(Gloc.o2up, Delta.o2up, beta, h, kt);
    std::vector<double> Ekin_o2do = ppsc::get_kinetic_energy(Gloc.o2do, Delta.o2do, beta, h, kt);

    std::vector<double> Ekin;
    Ekin.resize(Ekin_o1up.size());
    
    for( auto t : ppsc::range(0, Ekin.size()) ) {
      Ekin[t] = Ekin_o1up[t] + Ekin_o1do[t]
	      + Ekin_o2up[t] + Ekin_o2do[t];

      std::cout << "ekin, epot: " << t << " "
		<< Ekin[t] << " " << imp.hamiltonian.Eint_exp[t] //<< "\n";
		<< " n1, n2: " << imp.hamiltonian.n1_exp[t] << " " << imp.hamiltonian.n2_exp[t]
                << " n1u, n1d, n2u, n2d: " << imp.hamiltonian.n1u_exp[t] << " " << imp.hamiltonian.n1d_exp[t] << " " << imp.hamiltonian.n2u_exp[t] << " " << imp.hamiltonian.n2d_exp[t]
                << " d1, d2: " << imp.hamiltonian.d1_exp[t] << " " << imp.hamiltonian.d2_exp[t]
                << " m: " << imp.hamiltonian.m_exp[t]
                << " tp, t0, tm: " << imp.hamiltonian.tp_exp[t] << " " << imp.hamiltonian.t0_exp[t] << " " << imp.hamiltonian.tm_exp[t]
                << " tz, tx, ty: " << imp.hamiltonian.tz_exp[t] << " " << imp.hamiltonian.tx_exp[t] << " " << imp.hamiltonian.ty_exp[t]
                << " triplon: " << imp.hamiltonian.triplon_exp[t]
                << " Q: " << imp.hamiltonian.Q_exp[t] 
<< "ekin1u,1do,2u,2do " << Ekin_o1up[t] << " " << Ekin_o1do[t] << " " << Ekin_o2up[t] << " " << Ekin_o2do[t] << " " << "\n";
    } 

    // -------------------------------------------------------------------
    // OBSERVABLES
    {
      
      Gloc.o1up.print_to_file("Gloc.o1up.out");
      Gloc.o2up.print_to_file("Gloc.o2up.out");

      Gloc.o1do.print_to_file("Gloc.o1do.out");
      Gloc.o2do.print_to_file("Gloc.o2do.out");

      Gloc.otaup.print_to_file("Gloc.otaup.out");
      Gloc.otcup.print_to_file("Gloc.otcup.out");

      Gloc.otaupdo.print_to_file("Gloc.otaupdo.out");
      Gloc.otadoup.print_to_file("Gloc.otadoup.out");
      Gloc.otcupdo.print_to_file("Gloc.otcupdo.out");
      Gloc.otcdoup.print_to_file("Gloc.otcdoup.out");

SzSz.print_to_file("SzSz.out"); 
pairpair.print_to_file("pairpair.out");

tptp.print_to_file("tptp.out");
t0t0.print_to_file("t0t0.out");
tmtm.print_to_file("tmtm.out");

tztz.print_to_file("tztz.out");
txtx.print_to_file("txtx.out");
tyty.print_to_file("tyty.out");

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

