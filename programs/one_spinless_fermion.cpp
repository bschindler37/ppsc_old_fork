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

#include "./ppsc/hamiltonians/single_band_hubbard.hpp"

// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------

#include "./ppsc/hilbert_spaces/hilbert_space_base.hpp"

namespace ppsc {
namespace hilbert_spaces {
  
class one_spinless_fermion : public hilbert_space_base {

public:
  
  typedef hilbert_space_base base_type;

  void init() {
    int norb = 1;
    int spin_degeneracy = 1;
    bool spin_conservation = true;
    base_type::init(norb, spin_degeneracy, spin_conservation);
  }
  
};

} // namespace hilbert_spaces
} // namespace ppsc

// -----------------------------------------------------------------------

namespace ppsc {
namespace hamiltonians {

// -----------------------------------------------------------------------
template<class HILB> class one_spinless_fermion {

public:

  // ---------------------------------------------------------------------
  one_spinless_fermion(int nt, HILB & hilbert_space) :
    nt(nt), Q_exp(nt+2),
    n_exp(nt+2), Eint_exp(nt+2), hil_(hilbert_space) {

    // -- Construct basic operators

    ca = hil_.c_op_[hil_.flavor(0, 0, 0)];
    cc = hil_.c_op_[hil_.flavor(0, 0, 1)];

    n = cc * ca;

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double eps1 = eps[tstp+1] - mu;
    Htemp = eps1 * n;
  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1] = expectation_value(tstp, rho, Q).real();
    n_exp[tstp + 1] = expectation_value(tstp, rho, n).real();

    operator_type Ht;
    get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();
  }

  // ---------------------------------------------------------------------
  void local_obs(int tstp, double & n_out) {
    n_out = n_exp[tstp + 1];
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    store_double_attribute_to_hid(group_id, "mu", mu);
    store_real_data_to_hid(group_id, "eps", eps.data(), eps.size());

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "n_exp", n_exp.data(), n_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());
  }

  int nt;
  HILB hil_;

  double mu;
  std::vector<double> U, eps;
  std::vector<double> Q_exp, n_exp, Eint_exp;

  operator_type ca, cc, n, Q;
};

} // end namespace hamiltonians
} // end namespace ppsc

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & Delta,
			       ppsc::gf_type & Delta_cc,
			       HILB & hil_) {

  // (c)reation/(a)nihilation operators

  ppsc::operator_type ca = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  
  int fermion=-1, fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  
  pp_ints.push_back(ppsc::pp_int_type(Delta,    cc, ca, fermion, fwd));
  pp_ints.push_back(ppsc::pp_int_type(Delta_cc, ca, cc, fermion, bwd));
  
  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::gf_verts_type get_gf_verts(HILB & hil_) {

  ppsc::operator_type ca = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cc = hil_.c_op_[hil_.flavor(0, 0, 1)];

  ppsc::gf_verts_type gf_verts;
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, ca, cc));
  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int itermax, iter_rtime, nt, ntau, kt, iter, tstp, order, nomp;

  int store_pp, read_eq_sym, read_rt_sym, read_state_from_file;
  double linear_mixing;
  
  double beta, h, errmax, dmfterr, dmfterr_equil;
  bool matsubara_converged = false;

  double V, eps2;

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

      if (order > 1) {
        cout << "using an OCA impurity solver" << endl;
        find_param(argv[1], "__nomp=", nomp);
      } else {
	nomp = 1;
      }
    }
    // ---------------------------------------------------------------------
    // -- Setup pp calculator

    typedef ppsc::hilbert_spaces::one_spinless_fermion hilbert_space_type;
    typedef ppsc::hamiltonians::one_spinless_fermion<hilbert_space_type> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    hilbert_space.init();
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    
    { // -- Setup local Hamiltonian

      find_param(argv[1], "__mu=", imp.hamiltonian.mu);
      find_param_tvector(argv[1], "__eps1=", imp.hamiltonian.eps, nt);
      find_param(argv[1], "__eps2=", eps2);
      find_param(argv[1], "__V=", V);

      imp.update_hamiltonian();
    }
    // ---------------------------------------------------------------------
    // -- setup single particle greens functions and hybridizations

    cntr::herm_matrix<double> Gloc(nt, ntau, 1, -1);    
    cntr::herm_matrix<double> Delta(nt, ntau, 1, -1);
    cntr::herm_matrix<double> Delta_cc(nt, ntau, 1, -1);
    
    // ---------------------------------------------------------------------
    // -- Setup Hybridization function

    cntr::function<double> eps2_func(nt, 1);
    eps2_func.set_constant(eps2 * MatrixXcd::Identity(1, 1));

    double mu0 = 0.;
    bool fixHam = false;
    cntr::green_from_H(Delta, mu0, eps2_func, beta, h, fixHam);
    
    cntr::function<double> V_func(nt, 1);
    for(int tstp=-1;tstp<=nt;tstp++){
      cdmatrix tmp(1,1);
      tmp(0,0)=V;
      V_func.set_value(tstp,tmp);
    }

    for(int tstp=-1; tstp <= nt; tstp++) {
      Delta.left_multiply(tstp, V_func);
      Delta.right_multiply(tstp, V_func);
    }
    
    for(int tstp = -1; tstp <= nt; tstp++)
      ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);
    
    // ---------------------------------------------------------------------
    // MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
    {
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);

      if(!read_state_from_file) {
	imp.solve_atomic();
      } else {
	gtmp.set_timestep(-1, Gloc); // this avoids one iteration
      }

      for (iter = 1; iter <= itermax; iter++) {

	// -- Construct interactions and verticies
	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	
	// -- Solve pseudo particle problem
	imp.update_diagrams(pp_ints, gf_verts);
	imp.pp_step(-1);

	// -- get spgf
	ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);

	//Gloc.set_timestep(-1, gf_tstps[0]);	

	// -- linear mixing in Gloc
	ppsc::gf_tstp_type gloc_old(-1, ntau, 1);
	ppsc::gf_tstp_type gloc_mix(-1, ntau, 1);

	gloc_mix.clear();
	Gloc.get_timestep(-1, gloc_old);
	gloc_mix.incr(gloc_old, linear_mixing);
	gloc_mix.incr(gf_tstps[0], 1.0 - linear_mixing);
	Gloc.set_timestep(-1, gloc_mix);
	
	// -- Check error
	dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc);
        gtmp.set_timestep(-1, Gloc);

	// -- Update Hybridization
	//cntr::herm_matrix_timestep<double> tmp;
	//Gloc.get_timestep(-1, tmp);
	//Delta.set_timestep(-1, tmp);
	//ppsc::set_bwd_from_fwd(-1, Delta_cc, Delta);
	
	cout << "iter:  " << iter << " err: " << dmfterr_equil << endl;

        if (dmfterr_equil < errmax) {
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
    // 	START ... same iteration
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      
      for (iter = 1; iter <= itermax; iter++) {

	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	
	imp.update_diagrams(pp_ints, gf_verts);
	imp.pp_step(kt);

	for (int n = 0; n <= kt; n++) {
	  ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
	  Gloc.set_timestep(n, gf_tstps[0]);	

	  //cntr::herm_matrix_timestep<double> tmp;
	  //Gloc.get_timestep(n, tmp);
	  //Delta.set_timestep(n, tmp);
	  //ppsc::set_bwd_from_fwd(n, Delta_cc, Delta);
	}
	
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

	ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
	Gloc.set_timestep(tstp, gf_tstps[0]);	

	//cntr::herm_matrix_timestep<double> tmp;
	//Gloc.get_timestep(tstp, tmp);
	//Delta.set_timestep(tstp, tmp);
	//ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);	
      }
    }

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
