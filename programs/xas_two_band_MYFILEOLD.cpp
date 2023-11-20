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

#define SCREEN_HYB 1

#include <cntr/cntr.hpp>
#include <cntr/utils/read_inputfile.hpp>

// -----------------------------------------------------------------------

#include "./ppsc/ppsc.hpp"
#include "./ppsc/solver.hpp"

#include "./ppsc/hilbert_spaces/two_band_fermi_diag.hpp"	// Hilbert space for spin-conserving dynamics, 2 orbitals, spin up & down

#include "./ppsc/hamiltonians/two_band_hubbard_xas.hpp"

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------


bool file_exists(const std::string &filename) {
  struct stat buffer;
  
  return (stat(filename.c_str(), &buffer) == 0);
}

// -----------------------------------------------------------------------

void init_core_hybridization(cdvector &result, double width, double mean, double Gamma, double h, int nt ) {
	
	double t;
	result.resize(nt+1);
	for (int time_diff = 0; time_diff <= nt; time_diff++) {
		t = time_diff*h;
		result[time_diff] = -II * Gamma * std::exp(-II*t*mean) * std::exp(-t*t*width*width/2);
	}

}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & Delta_up,
			       ppsc::gf_type & Delta_up_cc,
			       HILB & hil_) {

  // -- spin (u)p/(d)own and (c)reation/(a)nnihilation operators
  // -- flavor(orbital, spin, annihilation/creation)
  // -- orbital 1 = valence (d), orbital 2 = core (c)
  
  ppsc::operator_type c1ua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type c1uc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type c1da = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type c1dc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  ppsc::operator_type c2ua = hil_.c_op_[hil_.flavor(1, 0, 0)];
  ppsc::operator_type c2uc = hil_.c_op_[hil_.flavor(1, 0, 1)];
  ppsc::operator_type c2da = hil_.c_op_[hil_.flavor(1, 1, 0)];
  ppsc::operator_type c2dc = hil_.c_op_[hil_.flavor(1, 1, 1)];
  
  ppsc::operator_type n1 = c1uc*c1ua + c1dc*c1da;  // density on orbital 1
  
  int fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  pp_ints.push_back(ppsc::pp_int_type(Delta_up, n1, n1, BOSON, fwd)); // forward
  pp_ints.push_back(ppsc::pp_int_type(Delta_up_cc, n1, n1, BOSON, bwd)); // backward
  
  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & K,	// K = Screening Hybridization, D = Core-hole bath
			       ppsc::gf_type & K_cc,
			       ppsc::gf_type & D,
			       ppsc::gf_type & D_cc,
			       HILB & hil_) {

  // -- spin (u)p/(d)own and (c)reation/(a)nnihilation operators
  // -- flavor(orbital, spin, annihilation/creation)
  // -- orbital 1 = valence (d), orbital 2 = core (c)
  
  ppsc::operator_type c1ua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type c1uc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type c1da = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type c1dc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  ppsc::operator_type c2ua = hil_.c_op_[hil_.flavor(1, 0, 0)];
  ppsc::operator_type c2uc = hil_.c_op_[hil_.flavor(1, 0, 1)];
  ppsc::operator_type c2da = hil_.c_op_[hil_.flavor(1, 1, 0)];
  ppsc::operator_type c2dc = hil_.c_op_[hil_.flavor(1, 1, 1)];
  
  ppsc::operator_type n1 = c1uc*c1ua + c1dc*c1da;  // density on orbital 1
  
  int fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  
  // -- CORE-HOLE BATH: spin up 
  pp_ints.push_back(ppsc::pp_int_type(D, c2uc, c2ua, FERMION, fwd)); // forward
  pp_ints.push_back(ppsc::pp_int_type(D_cc, c2ua, c2uc, FERMION, bwd)); // backward
  
  // -- CORE-HOLE BATH: spin down
  pp_ints.push_back(ppsc::pp_int_type(D, c2dc, c2da, FERMION, fwd)); // forward
  pp_ints.push_back(ppsc::pp_int_type(D_cc, c2da, c2dc, FERMION, bwd)); // backward
  
  // SCREENING in valence orbital
  pp_ints.push_back(ppsc::pp_int_type(K, n1, n1, BOSON, fwd));  // forward
  pp_ints.push_back(ppsc::pp_int_type(K_cc, n1, n1, BOSON, bwd));  // backward
  
  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::gf_verts_type get_gf_verts(HILB & hil_) {

  ppsc::operator_type c1ua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type c1uc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type c1da = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type c1dc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  ppsc::operator_type c2ua = hil_.c_op_[hil_.flavor(1, 0, 0)];
  ppsc::operator_type c2uc = hil_.c_op_[hil_.flavor(1, 0, 1)];
  ppsc::operator_type c2da = hil_.c_op_[hil_.flavor(1, 1, 0)];
  ppsc::operator_type c2dc = hil_.c_op_[hil_.flavor(1, 1, 1)];
  

  ppsc::gf_verts_type gf_verts;
  // -- remark: indices - here (0,0), (0,1) etc. - are just for bookkeeping -> check 'gf_vert_type' Constructor in ppsc.hpp file

  gf_verts.push_back(ppsc::gf_vert_type(0, 0, c1ua, c1uc)); // <d_up d_up_dag>
  gf_verts.push_back(ppsc::gf_vert_type(1, 1, c1da, c1dc)); // <d_down d_down_dag>
  gf_verts.push_back(ppsc::gf_vert_type(2, 2, c2ua, c2uc)); // <c_up c_up_dag>
  gf_verts.push_back(ppsc::gf_vert_type(3, 3, c2da, c2dc)); // <c_down c_down_dag>    
  
  gf_verts.push_back(ppsc::gf_vert_type(0, 2, c1ua, c2uc)); // <d_up c_up_dag>
  gf_verts.push_back(ppsc::gf_vert_type(1, 3, c1da, c2dc)); // <d_down c_down_dag> 
  
  // #################################### ??? ########################################################
  gf_verts.push_back(ppsc::gf_vert_type(2, 0, c2ua, c1uc)); // <c_up d_up_dag>
  gf_verts.push_back(ppsc::gf_vert_type(3, 1, c2da, c1dc)); // <c_down d_down_dag>  
  // #################################### ??? ########################################################
  
  
  return gf_verts;
}

// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

	// ---------------------------------------------------------------
	// SET PARAMETERS

	int ntau = 100, nt = 400, order = 1; 	// order = 1 (nca) or 2 (oca) for impurity solver
	int kt = 5; 			// kt = SolveOrder of time-stepping in NESSi
	int nomp = 1; 			// number of threads
	double beta = 40.0;		// inverse temperature
	double h = 0.02;		// time step width for _real_ time stepping procedure
	
	double e_c = 0.0; 		// energy core level in ROTATING WAVE APPROXIMATION (!)
	double e_v = 3.0; 		// energy valence level
	double Omega_0 = 1.0;		// bosonic frequency (screening mode)
	double gamma = 2.0;		// screening mode interaction with local valence electron density
	double Hubbard_U = 0.0; 	// Hubbard interaction of valence orbital
	
	double core_bandwidth = 1.0/h;	// core bath bandwidth energy \in [-core_bandwidth + e_c, + corebandwidth + e_c]
	double Gamma = 0.1;		// inverse lifetime of core hole
		
	double omega_in = atof(argv[1]);   	// energy of incoming X-Ray photon in ROTATING WAVE APPROXIMATION (!)
	double mean_pulse = (double(nt)*h / 2.0) , sigma_pulse = ((double(nt)*h - mean_pulse) / 3.0);	// X-ray probe pulse envelope
		
	double g_num = 0.001;		// numerical small factor to do first order expansion in light-matter interaction.
	
	int itermax = 400;		// maximum number of iterations (for pp & DMFT loop) at every time step
	double errmax = 1e-7;		// maximal error for 2-norm between two successive GFs in iteration
	double itererr;
	int iter;

	std::cout << "omega_in = " << omega_in << std::endl;

	// ---------------------------------------------------------------
	// SETUP pp calculator
	
	typedef ppsc::hilbert_spaces::two_band_fermi_diag hilbert_space_type;
	typedef ppsc::hamiltonians::two_band_hubbard_xas<hilbert_space_type> hamiltonian_type;
	typedef ppsc::solver<hamiltonian_type> solver_type;

	hilbert_space_type hilbert_space;
	hilbert_space.init();

	solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
	
	// -- Setup local Hamiltonian    
	imp.hamiltonian.U = Hubbard_U;
	imp.hamiltonian.e_valence = e_v;
	imp.hamiltonian.g = g_num;
			
	imp.hamiltonian.probe_pulse_env.resize(nt+2);
	imp.hamiltonian.omega_in_osc.resize(nt+2);
	imp.hamiltonian.e_core.resize(nt+2);
	double time_diff_aux = 0.0;
	for (int tstp = -1; tstp <= nt; tstp++) {
		if (tstp >= 0) {
			time_diff_aux = double(tstp) * h - mean_pulse; 
			imp.hamiltonian.probe_pulse_env[tstp+1] = std::exp(- (time_diff_aux * time_diff_aux) / (2*sigma_pulse*sigma_pulse) );
			imp.hamiltonian.e_core[tstp+1] = e_c;
		}
		else {	
			imp.hamiltonian.probe_pulse_env[tstp+1] = 0.0;
			imp.hamiltonian.e_core[tstp+1] = e_c - 20.0/beta;	// core level energy; zeroth timestep is for Matsubara branch (-> set differently because we prepare initial state)

		}
		imp.hamiltonian.omega_in_osc[tstp+1] = std::exp(-II*omega_in*h*double(tstp));	
	}
	
	imp.update_hamiltonian();
	
	// -- Setup core bath
	GREEN Hyb_cc(nt, ntau, 1, FERMION);	// initialized as zero
	cdvector core_hybridization;
	init_core_hybridization(core_hybridization, core_bandwidth, e_c, Gamma, h, nt);
	for (int tstp1 = 0; tstp1 <= nt; tstp1++){
		for (int tstp2 = tstp1; tstp2 <= nt; tstp2++) {
			Hyb_cc.set_ret( tstp2 , tstp1 , core_hybridization[tstp2-tstp1] ); 	// set only _retarded_ component non-zero
		}
	}
	
	//***** TEST OUTPUT *****/
	//Hyb_cc.print_to_file("./test_output_xas/test_Delta_cc.txt"); 
	
	// ---------------------------------------------------------------------
	// SETUP single particle Green's functions and hybridizations	

	GREEN Delta(nt, ntau, 1, BOSON), Delta_cc(nt, ntau, 1, BOSON);	// initialized as zero
	GREEN G_loc(nt, ntau, 4, FERMION);
	ppsc::gf_tstps_type gf_tstps; // vector of single particle GFs
	
	// ---------------------------------------------------------------------
	// INITIALIZATION
	
	if (SCREEN_HYB) {
		for(int tstp = -1; tstp <= nt; tstp++) {
			cntr::green_single_pole_XX_timestep(tstp, Delta, Omega_0, beta, h);	// calculate free bosonic propagator [-i < X(t) X(t') >] for one timestep
			Delta.smul(tstp, (gamma*gamma)/(2.0*Omega_0) ); // see my handwritten calculation on induced interaction by Holstein interaction U(t,t') = + (gamma^2 / 2) * D_free (t,t'), where D_free (t,t') = -i <X(t) X(t')> for real times and D_free(-i tau,-i tau') = - <X(-i tau) X(-i tau')> for imaginary times.
				// (1 / Omega_0) factor comes from def. of X = ( a + a^dag ) / sqrt(2) in NESSi [ vs. X = ( a + a^dag ) / sqrt(2*Omega_0) for physical displacement X ]
			// -- Use one of the two lines below (*) ; Delta is hermitian 
			Delta_cc.set_timestep(tstp,Delta); // (*)
			// ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta); // (*)
		}
	}
	
	// -- Solve atomic problem (here, Delta = Delta_cc = 0)
	imp.solve_atomic();
	imp.update_density_matrix(-1); 		// this routine is called when performing 'pp_step' function as well
	imp.hamiltonian.update_exp_vals(-1, imp.rho);	// this routine is called when performing 'pp_step' function as well
		
	// -- Construct interactions and vertices (aka impurity Green's function of interest) for Delta != 0, Delta_cc != 0
	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	
	imp.update_diagrams(pp_ints, gf_verts); 	// -- tell solver which diagrams to evaluate when performing 'get_spgf' or when doing 'pp_step'
	
	// -- Get single-particle GFs ('spgf')
	gf_tstps = imp.get_spgf(-1);
	
	// ---------------------------------------------------------------------
	// MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
		
	GREEN G_loc_temp(nt, ntau, 4, FERMION);	// temporary Green's function
	// -- Update temporary GFs according to atomic/initialized solution above
	for( auto idx : ppsc::range(0, gf_tstps.size()) ) {
	  G_loc_temp.set_matrixelement(-1, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps[idx]);
	}

	// -- TEST OUTPUT
	G_loc_temp.print_to_file("./test_output_xas/test_G_loc0.txt");		
	
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
		
		// -- Update G_loc
		for( auto idx : ppsc::range(0, gf_tstps.size()) ) {
		  G_loc.set_matrixelement(-1, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps[idx]);
		}
		
		// -- Calculate change with last solution
		itererr = cntr::distance_norm2(-1, G_loc, G_loc_temp);
		//std::cout << "iter:  " << iter << " err_upup: " << itererr0 << " err_downup: " << itererr1 << " ERROR_TOT: " << itererr_tot << std::endl;		
		
		// -- Update G_temp's
		G_loc_temp.set_timestep(-1, G_loc);
				
		// -- Check for convergence
		if ( itererr <= errmax ) {
			//std::cout << "\nMATSUBARA BRANCH -- Strong coupling self-consistency has converged! Exiting successfully... \n " << std::endl;
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
				
				// -- Update G_loc
				for( auto idx : ppsc::range(0, gf_tstps.size()) ) {
				  G_loc.set_matrixelement(tstp, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps[idx]);
				}
				
				// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
				//Delta.smul(tstp, -0.4*(gamma*gamma/2.0));
				//ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta); 
					
			}
			
			// -- Calculate change with last solution for time t = kt only !!!!
			itererr = cntr::distance_norm2(kt, G_loc, G_loc_temp);
			//std::cout << "iter:  " << iter << " err_upup: " << itererr0 << " err_downup: " << itererr1 << " ERROR_TOT: " << itererr_tot << std::endl;		
			
			// -- Update G_temp's
			// (here, in the bootrstrapping phase, we only need kt timeslice of G_temp's
			//	--> one could reduce memory here or define only kt time slice for temporary Green's functions)
			for (int tstp = 0; tstp <= kt; tstp++) {
				G_loc_temp.set_timestep(tstp, G_loc);
			}
			
			// -- Check for convergence
			if ( itererr <= errmax ) {
				//std::cout << "\nREAL-TIME BRANCH BOOTSTRAPPING -- Strong coupling self-consistency has converged! Exiting successfully... \n " << std::endl;
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
				
				// -- Update G_loc
				for( auto idx : ppsc::range(0, gf_tstps.size()) ) {
				  G_loc.set_matrixelement(tstp, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps[idx]);
				}
				
				// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
				//Delta.smul(tstp, -0.4*(gamma*gamma/2.0));
				//ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);
				
				// -- Calculate change with last solution
				itererr = cntr::distance_norm2(tstp, G_loc, G_loc_temp);
				//std::cout << "iter:  " << iter << " err_upup: " << itererr0 << " err_downup: " << itererr1 << " ERROR_TOT: " << itererr_tot << std::endl;		
				
				// -- Update G_temp's
				G_loc_temp.set_timestep(tstp, G_loc);

				// -- Check for convergence
				if ( itererr <= errmax ) {
					//std::cout << "\nREAL-TIME BRANCH (tstp = " << tstp << " > SolveOrder)"
					//	<< " -- Strong coupling self-consistency has converged! Exiting successfully... \n " << std::endl;
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
      	
      	// -- TEST OUTPUT
	G_loc.print_to_file("./test_output_xas/test_G_loc.txt");
	
      	// ---------------------------------------------------------------------
      	// ADD CORE BATH TO SOLUTION
      	
      	GREEN G_loc_new(nt, ntau, 4, FERMION);
      	
      	// -- Set up integration kernel F = - G_loc * Hyb_cc_4 (F_cc = - Hyb_cc_4 * G_loc) ; Q is equal to G_loc
      	GREEN F(nt, ntau, 4, BOSON), F_cc(nt, ntau, 4, BOSON);
      	GREEN Hyb_cc_4(nt, ntau, 4, FERMION);
	Hyb_cc_4.set_matrixelement(3, 3, Hyb_cc);	// core bath hybridization is spin symmetric
	Hyb_cc_4.set_matrixelement(2, 2, Hyb_cc);	// core bath hybridization is spin symmetric
	
      	cntr::convolution(F, G_loc, G_loc, Hyb_cc_4, Hyb_cc_4, beta, h, kt);
      	cntr::convolution(F_cc, Hyb_cc_4, Hyb_cc_4, G_loc, G_loc, beta, h, kt);
      	for (int tstp = -1; tstp <= nt; tstp++ ) {
      		F.smul(tstp, -1.0);
      		F_cc.smul(tstp, -1.0);
      	}
      	
      	// -- Matsubara branch
      	cntr::vie2_mat(G_loc_new, F, F_cc, G_loc, beta, kt, CNTR_MAT_FIXPOINT);
      	     	
      	// -- Bootstrapping (time <= kt)
	cntr::vie2_start(G_loc_new, F, F_cc, G_loc, beta, h, kt);

      	// -- Real time propagation
      	for(int tstp = kt; tstp <= nt; tstp++ ) {
      		cntr::vie2_timestep(tstp, G_loc_new, F, F_cc, G_loc, beta, h, kt);
      	}
      	
      	// -- Update G_loc
      	for (int tstp = -1; tstp <= nt; tstp++ ) {
		G_loc.set_timestep(tstp, G_loc_new);
      	}
      	
      	// ---------------------------------------------------------------------
      	// CALCULATE XAS SIGNAL
      	
      	// -- Update density matrix rho and expectation value < H_in^mat (t) > . xas_integrand = II * exp(-II * omega_in * t) * s(t) * < H_in^mat (t) > / g
      	cdmatrix rho(4,4);
      	std::vector<cdouble> xas_integrand(nt+1);	// only for real time arguments
      	std::vector<cdouble> xas_integrand_no_corebath(nt+1);	// only for real time arguments
      	for( int tstp = 0; tstp <= nt; tstp++ ) {
      		G_loc.density_matrix(tstp, rho);
      		xas_integrand[tstp] = h * II * std::exp(-II * omega_in * double(tstp) * h) * ( rho(2,0) + rho(3,1) ) * (imp.hamiltonian.probe_pulse_env[tstp+1]) / g_num ;
 		xas_integrand_no_corebath[tstp] =  h * II * std::exp(-II * omega_in * double(tstp) * h) * ( imp.hamiltonian.P_in_exp[tstp + 1]  ) * (imp.hamiltonian.probe_pulse_env[tstp+1]) / g_num;
      		
      	}
      	
      	// -- Integrate xas_integrand over [0,nt*h] to find I_XAS(omega_in)
      	integration::Integrator<double> gregory_integration(kt);
      	
      	cdouble I_XAS = gregory_integration.integrate(xas_integrand, nt);
      	cdouble I_XAS_no_corebath = gregory_integration.integrate(xas_integrand_no_corebath, nt);
      	

	// ---------------------------------------------------------------------
	// PRINT OUT / SAVING _REAL-TIME_ RESULTS FOR I_XAS
	
	std::ofstream I_XAS_file;

	std::ostringstream output_filename_stream;
	//output_filename_stream << "./test_output_xas/" << "Omega" << Omega_0 << "_gammatilde" << gamma << "_Gamma" << Gamma << "_ev" << e_v << "_ec" << e_c << "_U" << Hubbard_U << "/I_XAS_ntau" << ntau << "_nt" << nt << "_beta" << beta << "_h" << h << ".txt" ;
	output_filename_stream << "./test_output_xas/" << "NOSCREEN" << "_Gamma" << Gamma << "_ev" << e_v << "_ec" << e_c << "_U" << Hubbard_U << "/I_XAS_ntau" << ntau << "_nt" << nt << "_beta" << beta << "_h" << h << "_pulse-mean" << mean_pulse << "_pulse-sigma" << sigma_pulse << "_corebandwidth" << core_bandwidth << ".txt" ;
	std::string output_filename = output_filename_stream.str();
	//std::cout << "output_filename " << output_filename << std::endl; 
	
	// -- Check whether output file already exists ...
	if ( ! file_exists(output_filename) ) {
		I_XAS_file.open(output_filename);
		I_XAS_file << "############ XAS signal ############ \n##\n";
		I_XAS_file << "## Omega = " << Omega_0 << ", gammatilde = " << gamma << ", Gamma =" << Gamma << ", e_v = " << e_v << ", e_c = " << e_c << ", U = " << Hubbard_U << "\n";
		I_XAS_file << "## ntau = " << ntau << ", nt = " << nt << ", beta = " << beta << ", h = " << h << "\n##\n";
		I_XAS_file << "# omega_in" << "\t" << "I_XAS.real" << "\t" << "I_XAS.imag" << "\t" << "I_XAS_no_corebath.real" << "\t" << "I_XAS_no_corebath.imag";
	}
	else 	I_XAS_file.open(output_filename, std::ios_base::app); 
	
	// -- Write data to file
	if (I_XAS_file.is_open()) {
		I_XAS_file << "\n" << omega_in << "\t" << I_XAS.real() << "\t" << I_XAS.imag() << "\t" << I_XAS_no_corebath.real() << "\t" << I_XAS_no_corebath.imag();
		I_XAS_file.close();
	}
	else 	std::cout << "Unable to open file..." << std::endl;
	
	
	/*
	cdouble aux;
	// -- Save one-particle expectation values
	std::ofstream one_ptl_expvals ("./test_output_xas/test_one_particle_expvals.txt");
	if (one_ptl_expvals.is_open()) {
		one_ptl_expvals << "# I_XAS (omega) ..... SPECIFY!!!! ....\n";
		one_ptl_expvals << "# time \t xas_integrand (no core bath) REAL IMAG \t xas_integrand (incl. core bath) REAL IMAG \n";
		for (int tstp = -1; tstp <= nt; tstp ++) {
			aux = ( imp.hamiltonian.P_in_exp[tstp + 1]  ) * (imp.hamiltonian.probe_pulse_env[tstp+1]) / g_num;
			one_ptl_expvals << tstp << "\t" << aux.real() << "\t" << aux.imag() << "\t" << (xas_integrand(tstp+1).real() << "\t" << (xas_integrand(tstp+1)).imag() << "\n";
		}
		one_ptl_expvals.close();
	}
	else std::cout << "Unable to open file..." << std::endl;

	*/
	
	// ...
	
	return 0;

}


