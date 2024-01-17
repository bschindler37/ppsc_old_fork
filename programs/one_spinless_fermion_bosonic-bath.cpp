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

    Htemp = (eps-mu) * n;
    
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
    store_double_attribute_to_hid(group_id, "eps", eps);

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "n_exp", n_exp.data(), n_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());
  }

  int nt;
  HILB hil_;

  double mu, eps; 	// need to be initialized manually (!)
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
  
  ppsc::operator_type n = cc * ca;
  
  int fermion=-1, boson=+1, fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  
  pp_ints.push_back(ppsc::pp_int_type(Delta,    n, n, boson, fwd));
  pp_ints.push_back(ppsc::pp_int_type(Delta_cc, n, n, boson, bwd));
  
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

  int nt, ntau, kt=5, order = 2, nomp=1, itermax = 400, iter;
  double beta, h, mu, eps, gamma, gammatilde, Omega0, err, errmax = 1e-7;
  
  bool store_pp=true;	// enables storing of pp-gf and pp-selfenergy 


    // ---------------------------------------------------------------------
    // READ GENERAL INPUT (NOT YET MICROSCOPIC PARAMETERS)
    
    if (argc < 2)
      throw("COMMAND LINE ARGUMENT (INPUT FILE) MISSING");
      // scan the input file, double underscores to avoids mismatch
    find_param(argv[1], "__nt=", nt);
    find_param(argv[1], "__ntau=", ntau);
    find_param(argv[1], "__beta=", beta);
    find_param(argv[1], "__h=", h);
    find_param(argv[1], "__Omega0=", Omega0);
    find_param(argv[1], "__gamma=", gamma);
      
    gammatilde = gamma*std::sqrt(2*Omega0);
      
    
    // ---------------------------------------------------------------------
    // -- Setup pp calculator

    typedef ppsc::hilbert_spaces::one_spinless_fermion hilbert_space_type;
    typedef ppsc::hamiltonians::one_spinless_fermion<hilbert_space_type> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    hilbert_space.init();
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    
    // -- Setup local Hamiltonian

    find_param(argv[1], "__mu=", imp.hamiltonian.mu);
    find_param(argv[1], "__eps=", imp.hamiltonian.eps);

    imp.update_hamiltonian();
    
    
    // ---------------------------------------------------------------------
    // -- output directory
    
    std::ostringstream output_dir;	// output directory
    output_dir << "Omega" << Omega0 << "_gamma" << gamma << "_eps" << imp.hamiltonian.eps << "_mu" << imp.hamiltonian.mu << "_ntau" << ntau << "_nt" << nt << "_beta" << beta << "_h" << h;
    std::string tmp = output_dir.str();
    const char * output_dir_str = tmp.c_str();
    struct stat sb;
      if ( stat(output_dir_str, &sb) != 0) {
	// create directory 
	const int create_dir = mkdir( output_dir_str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);	// S_ISGID
	if (create_dir == -1) 	std::cerr << "Error :  " << std::strerror(errno) << std::endl;
	else 	std::cout << "Directory Created: " << output_dir.str() << std::endl;
      }
    
    // ---------------------------------------------------------------------
    // -- setup single particle Greens functions

    cntr::herm_matrix<double> Gloc(nt, ntau, 1, -1);    
    cntr::herm_matrix<double> Delta(nt, ntau, 1, +1);
    cntr::herm_matrix<double> Delta_cc(nt, ntau, 1, +1);
    
    // ---------------------------------------------------------------------
    // -- Setup Hybridization function, interaction vertices and green's function vertices

    for(int tstp = -1; tstp <= nt; tstp++) {	
      cntr::green_single_pole_XX_timestep(tstp, Delta, Omega0, beta, h);	// free bosonic propagator [-i < X(t) X(t') >] for one timestep
      Delta.smul(tstp, (gammatilde*gammatilde) / (2.0*Omega0) );
      ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);
    }
      
    ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
    ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
   
    imp.update_diagrams(pp_ints, gf_verts);
       
    // ---------------------------------------------------------------------
    // -- Initialization (tstp = -1)
    
    ppsc::gf_tstps_type gf_tstps;
    
    imp.solve_atomic();
    imp.update_density_matrix(-1);
    imp.hamiltonian.update_exp_vals(-1, imp.rho);
 
    gf_tstps = imp.get_spgf(-1);	// extract local GF without hybridization (because routine 'pp_step' has not been called yet)
    Gloc.set_timestep(-1, gf_tstps[0]); 
	
    // ---------------------------------------------------------------------
    // MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
    {
      
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);

      for (iter = 1; iter <= itermax; iter++) {

	// -- do pseudo-particle step (also updates expectation values)
	imp.pp_step(-1);

	// -- get spgf
	gf_tstps = imp.get_spgf(-1);
	Gloc.set_timestep(-1, gf_tstps[0]);
	
	// -- check error (difference to previous solution 'gtmp')
	err = cntr::distance_norm2(-1, gtmp, Gloc);
	
	// -- update 'gtmp'
        gtmp.set_timestep(-1, Gloc);

	// -- command line output
	std::cout << "iter:  " << iter << " n_exp: " << imp.hamiltonian.n_exp[0] << " err: " << err << std::endl;

        if (err < errmax)
          break;

      }
      
      // -- error handling
      if (iter > itermax) {
        std::cerr << "WARNING: Matsubara not converged after " << itermax
             << "steps ... abort" << std::endl;
        std::cerr << "skip real-time calculation " << std::endl;
        return 1;
      }
      
    }

    // ---------------------------------------------------------------------
    // 	START: ONLY FOR NT>0
    if (nt > 0) {

      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      
      for (iter = 1; iter <= itermax; iter++) {

	imp.pp_step(kt);

	for (int tstp = 0; tstp <= kt; tstp++) {
	  gf_tstps = imp.get_spgf(tstp);
	  Gloc.set_timestep(tstp, gf_tstps[0]);	

	}
	
        err = cntr::distance_norm2(kt, gtmp, Gloc);
        gtmp.set_timestep(kt, Gloc);	// copy only timeslice with tstp = kt to 'gtmp' as it is the latest time computed and the only one compared in the routine above
	
        std::cout << "START: iter:  " << iter << " n_exp(kt): " << imp.hamiltonian.n_exp[kt+1] << " err: " << err << std::endl;
        
        if (err < errmax)
          break;
    
      }
      
      // -- error handling
      if (iter > itermax) {
        std::cerr << "WARNING: Bootstrapping not converged after " << itermax
             << "steps ... abort" << std::endl;
        std::cerr << "skip real-time calculation for t>kt" << std::endl;
        return 1;
      }


    // ---------------------------------------------------------------------
    // 	REALTIME: ONLY FOR NT>0
      for (int tstp = kt + 1; tstp <= nt; tstp++) {
    
        cntr::herm_matrix<double> gtmp(tstp, ntau, 1, -1);

        imp.extrapolate_timestep(tstp - 1);
      
        for (iter = 1; iter <= itermax; iter++) {

	  imp.pp_step(tstp);

	  gf_tstps = imp.get_spgf(tstp);
	  Gloc.set_timestep(tstp, gf_tstps[0]);	
	
          err = cntr::distance_norm2(tstp, gtmp, Gloc);
          gtmp.set_timestep(tstp, Gloc);	
	  
	  std::cout << "tstp:  " << tstp << " n_exp: " << imp.hamiltonian.n_exp[tstp+1] << " err: " << err << std::endl;
	        
	  if (err < errmax)
            break;
          
        }
      
      // -- error handling
        if (iter > itermax) {
          std::cerr << "WARNING: Real-time evolution not converged after " << itermax
               << "steps ... abort" << std::endl;
          std::cerr << "skip remaining real-time calculation" << std::endl;
          return 1;
        }
      
      
      }

    }	// if nt > 0
    
    
    // ---------------------------------------------------------------------
    // GREEN'S FUNCTION
    std::ostringstream gf_outfile;
    gf_outfile << "/GF_oca.txt";
    Gloc.print_to_file((output_dir.str() + gf_outfile.str()).c_str());

    // ---------------------------------------------------------------------
    // OBSERVABLES
    {
      std::string filename = output_dir.str() + "/data_ppsc_oca.h5";
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
      
      close_hdf5_file(file_id);
    }
    

  return 0;
}



