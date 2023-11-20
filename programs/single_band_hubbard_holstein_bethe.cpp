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

#include "./ppsc/hilbert_spaces/single_band_fermi_diag.hpp"
#include "./ppsc/hilbert_spaces/single_band_fermi_spin.hpp"
#include "./ppsc/hilbert_spaces/single_band_fermi_spin_sc.hpp"

#include "./ppsc/hamiltonians/single_band_hubbard.hpp"

#include "./ppsc/baths/non_int_boson_propagator.hpp"

// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------

typedef ppsc::operator_type operator_type;
typedef ppsc::mam::dynamic_matrix_type matrix_type;

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & Delta_up,
                               ppsc::gf_type & Delta_up_cc,
                               ppsc::gf_type & D0,
                               ppsc::gf_type & D0_cc,
                               HILB & hil_) {

  // spin (u)p/(d)own and (c)reation/(a)nihilation operators

  operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  operator_type nu = cuc * cua;
  operator_type nd = cdc * cda;
  operator_type n = nu + nd;
  operator_type I = operator_type::Identity(hil_);
  operator_type nbar = n - I;

  //std::cout << "nu =" << std::endl << nu << std::endl;
  //std::cout << "nd =" << std::endl << nd << std::endl;
  //std::cout << "n =" << std::endl << n << std::endl;
  //std::cout << "identity op" << std::endl << I << std::endl;
  //std::cout << "nbar = " << std::endl << nbar << std::endl;

  int boson=+1, fermion=-1, fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  pp_ints.push_back(ppsc::pp_int_type(Delta_up,    cuc,  cua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta_up_cc, cua,  cuc,  fermion, bwd)); // spin up bwd
  
  pp_ints.push_back(ppsc::pp_int_type(Delta_up,    cdc,  cda,  fermion, fwd)); // spin do fwd // assuming spin sym
  pp_ints.push_back(ppsc::pp_int_type(Delta_up_cc, cda,  cdc,  fermion, bwd)); // spin do bwd // assuming spin sym

  pp_ints.push_back(ppsc::pp_int_type(D0,          nbar, nbar, boson,   fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(D0_cc,       nbar, nbar, boson,   bwd)); // spin up fwd

  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::gf_verts_type get_gf_verts(HILB & hil_) {

  ppsc::operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  operator_type nu = cuc * cua;
  operator_type nd = cdc * cda;
  operator_type n = nu + nd;
  operator_type I = operator_type::Identity(hil_);
  operator_type nbar = n - I;

  ppsc::gf_verts_type gf_verts;
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, cua, cuc)); // spin up
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, nbar, nbar)); // n - n
  return gf_verts;
}

 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int itermax, iter_rtime, nt, ntau, kt, iter, tstp, order, nomp, store_pp;
  double beta, h, errmax, dmfterr, V, mu, omega;
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
      find_param(argv[1], "__mu=", mu);
      find_param(argv[1], "__itermax=", itermax);
      find_param(argv[1], "__errmax=", errmax);
      find_param(argv[1], "__iter_rtime=", iter_rtime);
      find_param(argv[1], "__kt=", kt);
      find_param(argv[1], "__order=", order);
      find_param(argv[1], "__store_pp=", store_pp);
      if (order > 1) {
        cout << "using an OCA impurity solver" << endl;
        find_param(argv[1], "__nomp=", nomp);
      } else {
        nomp = 1;
      }
    }
    cntr::herm_matrix<double> gtmp(nt, ntau, 1, -1);
    // ---------------------------------------------------------------------
    // -- Setup pp calculator

    //typedef ppsc::hilbert_spaces::single_band_fermi_spin_sc hilbert_space_type;
    typedef ppsc::hilbert_spaces::single_band_fermi_diag hilbert_space_type;
    //typedef ppsc::hilbert_spaces::single_band_fermi_spin hilbert_space_type;
    typedef ppsc::hamiltonians::single_band_hubbard<hilbert_space_type> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    hilbert_space.init();
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    
    { // -- Setup local Hamiltonian
      imp.hamiltonian.mu = mu;
      find_param_tvector(argv[1], "__U=", imp.hamiltonian.U, nt);
      find_param_tvector(argv[1], "__eps=", imp.hamiltonian.eps_up, nt);
      find_param_tvector(argv[1], "__eps=", imp.hamiltonian.eps_do, nt);

      imp.update_hamiltonian();
    }
    // ---------------------------------------------------------------------
    // -- init the cluster and the "electric field"

    cntr::herm_matrix<double> Gloc_up(nt, ntau, 1, -1);
    cntr::herm_matrix<double> chi(nt, ntau, 1, -1);    
    cntr::herm_matrix<double> Delta_up(nt, ntau, 1, -1);
    cntr::herm_matrix<double> Delta_up_cc(nt, ntau, 1, -1);

    // ---------------------------------------------------------------------
    // -- Boson part
    
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

    // for(int tstp=-1;tstp<=nt;tstp++){
    //     cdmatrix tmp;
    //     gfunc.get_value(tstp,tmp);
    //     std::cout << "Boson param: g = " << tstp << " " << tmp  << std::endl;
    // }
    
    // std::cout << "Boson param: omega0 = " << tstp << " " << omega << std::endl;

    for(int tstp=-1; tstp <= nt; tstp++) {
      D0.left_multiply(tstp, gfunc);
      D0.right_multiply(tstp, gfunc);
    }

    for(int tstp = -1; tstp <= nt; tstp++)
      ppsc::set_bwd_from_fwd(tstp, D0_cc, D0);

    // ---------------------------------------------------------------------
    // MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
    {
      
      imp.solve_atomic();
      for(int s=0;s<imp.ppGfs.size();s++){
        std::cout  << "INIT: n" << s << ": " << imp.ppGfs[s].density_matrix(-1) << " " ;
      }
      std::cout << std::endl;

      for (iter = 1; iter <= itermax; iter++) {

        // -- Construct interactions and verticies
        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta_up, Delta_up_cc, D0, D0_cc, imp.hamiltonian.hil_);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
        
        // -- Solve pseudo particle problem
        imp.update_diagrams(pp_ints, gf_verts);
        imp.pp_step(-1);

        // -- get spgf
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);
        Gloc_up.set_timestep(-1, gf_tstps[0]);  
        chi.set_timestep(-1, gf_tstps[1]); 

        // -- Check error
        dmfterr = cntr::distance_norm2(-1, gtmp, Gloc_up);
        gtmp.set_timestep(-1, Gloc_up);

        // -- Update Hybridization
        cntr::herm_matrix_timestep<double> tmp;
        Gloc_up.get_timestep(-1, tmp);
        Delta_up.set_timestep(-1, tmp);
        ppsc::set_bwd_from_fwd(-1, Delta_up_cc, Delta_up);
        
        cout << "iter:  " << iter << " err: " << dmfterr << endl;
        for(int s=0;s<imp.ppGfs.size();s++){
              std::cout  << "n" << s << ": " << imp.ppGfs[s].density_matrix(-1) << " " ;
        }
        std::cout << std::endl;

        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }

        if(iter == 2) imp.symmetry_reduction(-1);

      }
      if (iter > itermax) {
        cerr << "WARNING: Matsubara not converged  after " << itermax
             << "steps ... abort" << endl;
        cerr << "skip real-time calculation " << endl;
      }
    }
    
    // ---------------------------------------------------------------------

    if(nt > 0) imp.clear_symmetries();

    // ---------------------------------------------------------------------
    //  START ... same iteration
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      // cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      
      for (iter = 1; iter <= itermax; iter++) {

        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta_up, Delta_up_cc, D0, D0_cc, imp.hamiltonian.hil_);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
        
        imp.update_diagrams(pp_ints, gf_verts);
        imp.pp_step(kt);

        for (int n = 0; n <= kt; n++) {
          ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
          Gloc_up.set_timestep(n, gf_tstps[0]);
          chi.set_timestep(n, gf_tstps[1]); 

          cntr::herm_matrix_timestep<double> tmp;
          Gloc_up.get_timestep(n, tmp);
          Delta_up.set_timestep(n, tmp);
          ppsc::set_bwd_from_fwd(n, Delta_up_cc, Delta_up);
        }
        
        for (int n = 0; n <= kt; n++){
          dmfterr = cntr::distance_norm2(n, gtmp, Gloc_up);
          gtmp.set_timestep(n, Gloc_up);
          cout << "START time: " << n  << " iter:  " << iter << " err: " << dmfterr << endl;
        }
        std::cout << std::endl;
        
        
        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }
      }
    }

    // ---------------------------------------------------------------------
    // -- Reduce number of diagrams using symmetries

    if(matsubara_converged)
      imp.symmetry_reduction(kt);
    
    // ---------------------------------------------------------------------
    //  REALTIME: ONLY FOR NT>0
    if(matsubara_converged)
    for (tstp = kt + 1; tstp <= nt; tstp++) {
      // cout << "tstp:  " << tstp << endl;
      imp.extrapolate_timestep(tstp - 1);
      
      for (iter = 1; iter <= iter_rtime; iter++) {

        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta_up, Delta_up_cc, D0, D0_cc, imp.hamiltonian.hil_);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);

        imp.update_diagrams(pp_ints, gf_verts);
        imp.pp_step(tstp);

        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
        Gloc_up.set_timestep(tstp, gf_tstps[0]);        

        cntr::herm_matrix_timestep<double> tmp;
        Gloc_up.get_timestep(tstp, tmp);
        Delta_up.set_timestep(tstp, tmp);
        ppsc::set_bwd_from_fwd(tstp, Delta_up_cc, Delta_up);

        dmfterr = cntr::distance_norm2(tstp, gtmp, Gloc_up);
        gtmp.set_timestep(tstp, Gloc_up);
        cout << "time: " << tstp  << " iter  " << iter << " err: " << dmfterr << endl;
        if (dmfterr < errmax) {
          break;
        }else if(iter==iter_rtime){
          std::cout << "AT TIME " <<tstp <<  " NOT CONVERGED AFTER " << iter_rtime << " ITERATIONS " << std::endl;
        }
      }
      std::cout << std::endl;
    }

    // ---------------------------------------------------------------------
    // Kinetic energy

    // imp.calculate_pp_interaction_energy();
      
    std::vector<double> Ekin = ppsc::get_kinetic_energy(Gloc_up, Delta_up, beta, h, kt);
    for( auto t : ppsc::range(0, Ekin.size()) ) Ekin[t] *= 2; // two spin species

    // ---------------------------------------------------------------------
    // OBSERVABLES
    {
      std:string filename = "data_ppsc.h5";
      hid_t file_id = open_hdf5_file(filename);
      hid_t group_id;

      imp.store(file_id, store_pp);
      
      group_id = create_group(file_id, "g");
      store_herm_greens_function(group_id, Gloc_up);
      close_group(group_id);

      group_id = create_group(file_id, "chi");
      store_herm_greens_function(group_id, chi);
      close_group(group_id);

      group_id = create_group(file_id, "D0");
      store_herm_greens_function(group_id, D0);
      close_group(group_id);
      
      group_id = create_group(file_id, "holstein");
      store_real_data_to_hid(group_id, "Ekin", Ekin.data(), Ekin.size());
      store_real_data_to_hid(group_id, "g_int", g.data(),g.size());
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
