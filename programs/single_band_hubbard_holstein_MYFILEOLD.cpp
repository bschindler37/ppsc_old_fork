// -----------------------------------------------------------------------

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <fstream>

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

#include "./ppsc/hilbert_spaces/single_band_fermi_diag.hpp"	// Hilbert space for spin-conserving dynamics, 1 orbital, spin up & down

#include "./ppsc/hamiltonians/single_band_hubbard.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & Delta_up,
			       ppsc::gf_type & Delta_up_cc,
			       HILB & hil_) {

  // -- spin (u)p/(d)own and (c)reation/(a)nihilation operators
  // -- flavor(orbital, spin, annihilation/creation)
  ppsc::operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];
  
  ppsc::operator_type nd = cdc * cda;	// spin down density
  ppsc::operator_type nu = cuc * cua;	// spin up density
  ppsc::operator_type n = nu + nd;
  
  int fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  pp_ints.push_back(ppsc::pp_int_type(Delta_up, n, n, BOSON, fwd)); // forward
  pp_ints.push_back(ppsc::pp_int_type(Delta_up_cc, n, n, BOSON, bwd)); // backward
  
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
  // -- remark: indices - here (0,0) - are just for bookkeeping -> check 'gf_vert_type' Constructor in ppsc.hpp file
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, cua, cuc)); // <c_up c_up_dag>
  gf_verts.push_back(ppsc::gf_vert_type(1, 0, cda, cuc)); // <c_down c_up_dag>
    
  return gf_verts;
}

// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

	// ---------------------------------------------------------------
	// SET PARAMETERS

	int ntau = 200, nt = 10, order = 1; 	// order = 1 (nca) or 2 (oca) for impurity solver
	int kt = 5; 	// kt = SolveOrder of time-stepping in NESSi
	int nomp = 1; 	// number of threads
	double beta = 1;	// inverse temperature
	double h = 0.02;	// time step width for _real_ time stepping procedure
	
	double Omega_0 = 10.0;	// bosonic frequency
	double gamma = 2.0;	// screening mode interaction with local electron density
	double Hubbard_U = 1.0;
	
	int dim = 1;	// dimension of _local_ Hamiltonian
	
	int itermax = 200;	// maximum number of iterations (for pp & DMFT loop) at every time step
	double errmax = 1e-10;		// maximal error for 2-norm between two successive GFs in iteration
	double itererr0 = 0.0, itererr1 = 0.0, itererr_tot = 0.0;
	int iter;

	// ---------------------------------------------------------------
	// SETUP pp calculator
	
	typedef ppsc::hilbert_spaces::single_band_fermi_diag hilbert_space_type;
	typedef ppsc::hamiltonians::single_band_hubbard<hilbert_space_type> hamiltonian_type;
	typedef ppsc::solver<hamiltonian_type> solver_type;

	hilbert_space_type hilbert_space;
	hilbert_space.init();

	solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
	
	// -- Setup local Hamiltonian    
	imp.hamiltonian.mu = - Hubbard_U / 2.0;	// chemical potential
	
	imp.hamiltonian.eps_up.resize(nt+2);
	imp.hamiltonian.eps_do.resize(nt+2);
	imp.hamiltonian.U.resize(nt+2);
	for (int tstp = -1; tstp <= nt; tstp++) {
		imp.hamiltonian.eps_up[tstp+1] = 0.0;	// energy level of up electrons; zeroth timestep is for Matsubara branch (-> set differently for quench)
		imp.hamiltonian.eps_do[tstp+1] = 0.0;	// energy level of down electrons; zeroth timestep is for Matsubara branch (-> set differently for quench)
		imp.hamiltonian.U[tstp+1] = Hubbard_U;	// Hubbard U interaction; zeroth timestep is for Matsubara branch (-> set differently for quench)
	}
	
	imp.update_hamiltonian(); 
	
	// ---------------------------------------------------------------------
	// SETUP single particle Green's functions and hybridizations	

	GREEN G_upup(nt, ntau, dim, FERMION), G_downup(nt, ntau, dim, FERMION), Delta(nt, ntau, dim, BOSON), Delta_cc(nt, ntau, dim, BOSON);
	ppsc::gf_tstps_type gf_tstps; // vector of single particle GFs G_upup & G_downup
	
	// ---------------------------------------------------------------------
	// INITIALIZATION
	
	for(int tstp = -1; tstp <= nt; tstp++) {
		cntr::green_single_pole_XX_timestep(tstp, Delta, Omega_0, beta, h);	// calculate free bosonic propagator [-i < X(t) X(t') >] for one timestep
		Delta.smul(tstp, (gamma*gamma/2.0)); // see my handwritten calculation on induced interaction by Holstein interaction U(t,t') = + (gamma^2 / 2) * D_free (t,t'), where D_free (t,t') = -i <X(t) X(t')> for real times and D_free(-i tau,-i tau') = - <X(-i tau) X(-i tau')> for imaginary times.
		
		// -- Use one of the two lines below (*) ; Delta is hermitian 
		Delta_cc.set_timestep(tstp,Delta); // (*)
		// ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta); // (*)
	}
	
	// -- Solve atomic problem (here, Delta = Delta_cc = 0)
	imp.solve_atomic();
	imp.update_density_matrix(-1); 		// this routine is calles when performing 'pp_step' function as well
	imp.hamiltonian.update_exp_vals(-1, imp.rho);	// this routine is calles when performing 'pp_step' function as well
	
	// -- Construct interactions and vertices (aka impurity Green's function of interest) for Delta != 0, Delta_cc != 0
	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	
	imp.update_diagrams(pp_ints, gf_verts); 	// -- tell solver which diagrams to evaluate when performing 'get_spgf' or when doing 'pp_step'
	
	// -- Get single-particle GFs ('spgf')
	gf_tstps = imp.get_spgf(-1);
	
	// ---------------------------------------------------------------------
	// MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
	
	GREEN G_temp_upup(nt, ntau, dim, FERMION), G_temp_downup(nt, ntau, dim, FERMION);	// temporary Green's function
	// -- Update temporary GFs according to atomic/initialized solution above
	G_temp_upup.set_timestep(-1, gf_tstps[0]);
	G_temp_downup.set_timestep(-1, gf_tstps[1]);

	// -- TEST OUTPUT
	G_temp_upup.print_to_file("./test_output/test_upup_output0.txt");		

	// -- Iteration loop
	for ( iter = 1 ; iter <= itermax; iter++) {
		
		// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
		//Delta.smul(-1, -0.4*(gamma*gamma/2.0));
		//ppsc::set_bwd_from_fwd(-1, Delta_cc, Delta); 
			
		// -- Solve pseudo particle problem --> this is probably only needed if we change Delta (see above); otherwise just has to be called once before performing 'pp_step' or 'get_spgf' (see 'Initialization' above)
		//imp.update_diagrams(pp_ints, gf_verts);
		
		// -- Do the (Matsubara) time step
		imp.pp_step(-1);
		
		// -- Get single-particle GFs ('spgf')
		gf_tstps = imp.get_spgf(-1);
		
		// -- Update G_upup & G_downup
		G_upup.set_timestep(-1, gf_tstps[0]);
		G_downup.set_timestep(-1, gf_tstps[1]);
		
		// -- Calculate change with last solution
		itererr0 = cntr::distance_norm2(-1, G_upup, G_temp_upup);
		itererr1 = cntr::distance_norm2(-1, G_downup, G_temp_downup);
		itererr_tot = std::sqrt(itererr0*itererr0 + itererr1*itererr1);
		//std::cout << "iter:  " << iter << " err_upup: " << itererr0 << " err_downup: " << itererr1 << " ERROR_TOT: " << itererr_tot << std::endl;		
		
		// -- Update G_temp_upup & G_temp_downup
		G_temp_upup.set_timestep(-1, G_upup);
		G_temp_downup.set_timestep(-1, G_downup);
				
		// -- Check for convergence
		if ( itererr_tot <= errmax ) {
			std::cout << "\nMATSUBARA BRANCH -- Strong coupling self-consistency has converged! Exiting successfully... \n " << std::endl;
			//G_temp_upup.print_to_file("./test_output/test_upup_output1.txt"); 		// -- TEST OUTPUT
			break;
		}
			
	}
	
	// -- Handling non-convergence
	if (iter > itermax) {
        	std::cerr << "WARNING: Matsubara not converged after " << itermax
        	     << " steps ... ABORT" << std::endl;
        	std::cerr << "skip real-time calculation " << std::endl;
        	return 1;
      	}

	// ---------------------------------------------------------------------
	// REAL-TIME PART ('BOOTSTRAPPING PHASE': n <= kt)
	
	if (nt > 0){
		
		// cntr routine to set t=0 components of herm_matrix (ret,lesser, left-mixing) from Matsubara component
		// (defined in line 400ff in 'cntr_utilities_impl.hpp')
		imp.init_real_time();	// should do the same as 'imp.extrapolate_timestep(-1)'
	
		// -- Iteration loop
		for ( iter = 1 ; iter <= itermax; iter++) {
			
			// -- Solve pseudo particle problem
			//	 --> this is probably only needed if we change Delta (see below);
			//	otherwise just has to be called once before performing 'pp_step' or 'get_spgf' (see 'Initialization' above)
			//imp.update_diagrams(pp_ints, gf_verts);
			
			// -- Do the (Bootstrapping) time step; for kt = SolveOrder, the Dyson equation is solved for all time steps 0, 1, ... kt together
			imp.pp_step(kt);	// kt = SolveOrder
		
			for (int tstp = 0; tstp <= kt; tstp++) {
			
				// -- Get single-particle GFs ('spgf')
				gf_tstps = imp.get_spgf(tstp);
				
				// -- Update G_upup & G_downup
				G_upup.set_timestep(tstp, gf_tstps[0]);
				G_downup.set_timestep(tstp, gf_tstps[1]);
				
				// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
				//Delta.smul(tstp, -0.4*(gamma*gamma/2.0));
				//ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta); 
					
			}
			
			// -- Calculate change with last solution for time t = kt only !!!!
			itererr0 = cntr::distance_norm2(kt, G_upup, G_temp_upup);
			itererr1 = cntr::distance_norm2(kt, G_downup, G_temp_downup);
			itererr_tot = std::sqrt(itererr0*itererr0 + itererr1*itererr1);
			//std::cout << "iter:  " << iter << " err_upup: " << itererr0 << " err_downup: " << itererr1 << " ERROR_TOT: " << itererr_tot << std::endl;		
			
			// -- Update G_temp_upup & G_temp_downup
			// (here, in the bootrstrapping phase, we only need kt timeslice of G_temp_upup & G_temp_downup 
			//	--> one could reduce memory here or define only kt time slice for temporary Green's functions)
			for (int tstp = 0; tstp <= kt; tstp++) {
				G_temp_upup.set_timestep(tstp, G_upup);
				G_temp_downup.set_timestep(tstp, G_downup);
			}
			
			// -- Check for convergence
			if ( itererr_tot <= errmax ) {
				std::cout << "\nREAL-TIME BRANCH BOOTSTRAPPING -- Strong coupling self-consistency has converged! Exiting successfully... \n "
					<< std::endl;
				break;
			}
		
		}
	}
	
	// -- Handling non-convergence
	if (iter > itermax) {
        	std::cerr << "WARNING: Bootstrapping not converged after " << itermax
        	     << " steps ... ABORT" << std::endl;
        	std::cerr << "skip real-time calculation for t > kt" << std::endl;
        	return 1;
      	}
	
	// ---------------------------------------------------------------------
	// REAL-TIME PART (TIMESTEPPING: nt > kt)
	
	if (nt > kt) {
		
		// -- Time stepping loop
		for (int tstp = kt+1; tstp <= nt; tstp++) {
			
			// -- First guess for GF via extrapolation from previous, solved GF time-slices; see also NESSi example programs
			imp.extrapolate_timestep(tstp-1);
			
			// -- Iteration loop
			for (iter = 1; iter <= itermax ; iter++) {
			
				// -- Solve pseudo particle problem
				//	 --> this is probably only needed if we change Delta (see below);
				//	otherwise just has to be called once before performing 'pp_step' or 'get_spgf' (see 'Initialization' above)
				//imp.update_diagrams(pp_ints, gf_verts);
				
				// -- Do the real-time step
				imp.pp_step(tstp);
				
				// -- Get single-particle GFs ('spgf')
				gf_tstps = imp.get_spgf(tstp);
				
				// -- Update G_upup & G_downup
				G_upup.set_timestep(tstp, gf_tstps[0]);
				G_downup.set_timestep(tstp, gf_tstps[1]);
				
				// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
				//Delta.smul(tstp, -0.4*(gamma*gamma/2.0));
				//ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);
				
				// -- Calculate change with last solution
				itererr0 = cntr::distance_norm2(tstp, G_upup, G_temp_upup);
				itererr1 = cntr::distance_norm2(tstp, G_downup, G_temp_downup);
				itererr_tot = std::sqrt(itererr0*itererr0 + itererr1*itererr1);
				//std::cout << "iter:  " << iter << " err_upup: " << itererr0 << " err_downup: " << itererr1 << " ERROR_TOT: " << itererr_tot << std::endl;		
				
				// -- Update G_temp_upup & G_temp_downup
				G_temp_upup.set_timestep(tstp, G_upup);
				G_temp_downup.set_timestep(tstp, G_downup);
				
				// -- Check for convergence
				if ( itererr_tot <= errmax ) {
					std::cout << "\nREAL-TIME BRANCH (tstp = " << tstp << " > SolveOrder)"
						<< " -- Strong coupling self-consistency has converged! Exiting successfully... \n " << std::endl;
					break;	// breaks iteration loop, not timetepping loop over 'tstp'
				}
				
			}
			
		}
		
	}
	
	// -- Handling non-convergence
	if (iter > itermax) {
        	std::cerr << "WARNING: real time evolution not converged after " << itermax
        	     << " steps ... ABORT" << std::endl;
        	return 1;
      	}
      	
      	
	// ---------------------------------------------------------------------
	// PRINT OUT / SAVING RESULTS
	
	// -- TEST OUTPUT
	G_temp_upup.print_to_file("./test_output/test_upup_output_final.txt");
	
	// -- Save one-particle expectation values
	std::ofstream one_ptl_expvals ("./test_output/test_one_particle_expvals.txt");
	if (one_ptl_expvals.is_open()) {
		one_ptl_expvals << "# One-particle expectation values\n";
		one_ptl_expvals << "# time \t pseudo particle number Q \t n_up \t n_down \t double occ. n_updown \t local energy \n";
		for (int tstp = -1; tstp <= nt; tstp ++) {
			one_ptl_expvals << tstp << "\t" << imp.hamiltonian.Q_exp[tstp+1] << "\t" << imp.hamiltonian.nu_exp[tstp+1] << "\t" << imp.hamiltonian.nd_exp[tstp+1] << "\t" << imp.hamiltonian.docc_exp[tstp+1] << "\t" << imp.hamiltonian.Eint_exp[tstp+1] << "\n";
		}
		one_ptl_expvals.close();
	}
	else std::cout << "Unable to open file..." << std::endl;
	

	
	
	// ...
	
	return 0;

}


