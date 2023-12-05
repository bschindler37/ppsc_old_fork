// -----------------------------------------------------------------------

#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <fstream>
#include <hdf5.h>


// -----------------------------------------------------------------------

#ifndef CNTR_USE_OMP
#define CNTR_USE_OMP
#endif

#ifndef CNTR_USE_MPI
#define CNTR_USE_MPI
#endif

#define NCA_SOLVER_ASSERT_0 0
#define NCA_SOLVER_ASSERT_1 0

#define SCREEN_HYB 0
#define DEBUG 0

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
template<class HAMILTONIAN>
void initialize_hamiltonian(HAMILTONIAN & ham_, double Hubbard_U, double e_v, double e_c, double g_num, int nt, double h, double beta, double omega_in, double sigma_pulse, double mean_pulse ){

	ham_.U = Hubbard_U;
	ham_.e_valence = e_v;
	ham_.g = g_num;
	ham_.h = h;
			
	ham_.probe_pulse_env.resize(nt+2);
	ham_.omega_in_osc.resize(nt+2);
	ham_.e_core.resize(nt+2);
	double time_diff_aux = 0.0;
	for (int tstp = -1; tstp <= nt; tstp++) {
		if (tstp >= 0) {
			time_diff_aux = double(tstp)*h - mean_pulse;
			ham_.probe_pulse_env[tstp+1] = std::exp(- (time_diff_aux * time_diff_aux) / (2*sigma_pulse*sigma_pulse) ) / (sigma_pulse*sqrt(2*PI)) ;
			ham_.e_core[tstp+1] = e_c;
		}
		else {	
			ham_.probe_pulse_env[tstp+1] = 0.0;
			ham_.e_core[tstp+1] = e_c - 20.0/beta;	// core level energy; zeroth timestep is for Matsubara branch (-> set differently because we prepare initial state)
		}
		ham_.omega_in_osc[tstp+1] = std::exp(-II*omega_in*h*double(tstp));	
	}
	


}

// -----------------------------------------------------------------------
void init_core_hybridization(cdvector &result, double width, double mean, double Gamma, double h, int nt ) {
	
	double t;
	result.resize(nt+1);
	for (int time_diff = 0; time_diff <= nt; time_diff++) {
		t = time_diff*h;
		result[time_diff] = -II * Gamma * ( width * sqrt(2/PI) ) * std::exp(-II*t*mean) * std::exp(-t*t*width*width/2);
	} 

}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & D,	// D = Core-hole bath
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
ppsc::pp_ints_type get_pp_ints_CBCON(ppsc::gf_type & K,	// K = Screening Hybridization
			       ppsc::gf_type & K_cc,
			       HILB & hil_) {

  // -- spin (u)p/(d)own and (c)reation/(a)nnihilation operators
  // -- flavor(orbital, spin, annihilation/creation)
  // -- orbital 1 = valence (d), orbital 2 = core (c)
  
  ppsc::operator_type c1ua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type c1uc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type c1da = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type c1dc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  /*
  ppsc::operator_type c2ua = hil_.c_op_[hil_.flavor(1, 0, 0)];
  ppsc::operator_type c2uc = hil_.c_op_[hil_.flavor(1, 0, 1)];
  ppsc::operator_type c2da = hil_.c_op_[hil_.flavor(1, 1, 0)];
  ppsc::operator_type c2dc = hil_.c_op_[hil_.flavor(1, 1, 1)];
  */ 
  
  ppsc::operator_type n1 = c1uc*c1ua + c1dc*c1da;  // density on orbital 1 = valence orbital
  
  int fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  
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

	if (argc < 3) { std::cerr << " \n !!! GIVE 	'omega_in' 'input_file.txt' 	AS COMMAND LINE ARGUMENNTS !!! \n " << std::endl; return 1; }
	// ---------------------------------------------------------------
	// SET PARAMETERS

	int ntau, nt, order = 1; 	// order = 1 (nca) or 2 (oca) for impurity solver
	int kt = 5; 			// kt = SolveOrder of time-stepping in NESSi
	int nomp = 1; 			// number of threads
	double beta;		// inverse temperature
	double h;		// time step width for _real_ time stepping procedure
	
	double e_c; 		// energy core level in ROTATING WAVE APPROXIMATION (!)
	double e_v; 		// energy valence level
	double Omega_0;		// bosonic frequency (screening mode)
	double gammatilde;		// screening mode interaction with local valence electron density
	double Hubbard_U; 	// Hubbard interaction of valence orbital
	
	double core_bandwidth;	// core bath bandwidth energy \in [-core_bandwidth + e_c, + corebandwidth + e_c]
	double Gamma;		// inverse lifetime of core hole
		
	double omega_in = atof(argv[1]);   	// energy of incoming X-Ray photon in ROTATING WAVE APPROXIMATION (!)
	double mean_pulse, sigma_pulse; 	// = (double(nt)*h / 2.0) , sigma_pulse = ((double(nt)*h - mean_pulse) / 3.0);	// X-ray probe pulse envelope

	
	double g_num;		// numerical small factor to do first order expansion in light-matter interaction.
	
	int itermax = 400;		// maximum number of iterations (for pp & DMFT loop) at every time step
	double errmax = 1e-7;		// maximal error for 2-norm between two successive GFs in iteration
	double itererr_CBEXP, itererr_CBCON;
	int iter;

	std::cout << "omega_in = " << omega_in << std::endl; //  << omega_in << " , sigma_pulse = " << sigma_pulse << " , mean_pulse = " << mean_pulse << std::endl;
	
	// OUTPUT DIRECTORY
	std::string output_dir_read_in = "./xas_data/";
	
	// ---------------------------------------------------------------
	// MEMORY REQUIREMENTS
	/*
        std::cout << std::string(50, '-') << "\n" << "     MEMORY requirements" << "\n" << std::string(50, '-') << std::endl;

        const size_t size_MB=1024*1024;
        size_t mem_time=0;

        mem_time += 4*cntr::mem_function<cdouble>(nt,4); 	 // GFs G_loc for CBEXP, CBCON (bith temporary & non-temporary)
        mem_time += 4*cntr::mem_function<cdouble>(nt,1); 	 // Hyb & Hyb_cc and Delta & Delta_cc
        
        mem_time = ceil(mem_time/(double)size_MB);	        // convert to MB

	std::cout << "\t" << mem_time << " MB \n" << std::string(50, '-') << "\n\n" << std::endl;
	*/
	
	// ---------------------------------------------------------------
	// READ PARAMETER VALUES FROM INPUT FILE
	
	// -- scan the input file, double underscores to avoids mismatch
	find_param(argv[2], "__nt=", nt);
	find_param(argv[2], "__ntau=", ntau);
	find_param(argv[2], "__beta=", beta);
	find_param(argv[2], "__h=", h);
	find_param(argv[2], "__e_c=", e_c);
	find_param(argv[2], "__e_v=", e_v);
	find_param(argv[2], "__Omega_0=", Omega_0);
	find_param(argv[2], "__gammatilde=", gammatilde);
	find_param(argv[2], "__mean_pulse=", mean_pulse);
	find_param(argv[2], "__sigma_pulse=", sigma_pulse);
	find_param(argv[2], "__g_num=", g_num);
	find_param(argv[2], "__Hubbard_U=", Hubbard_U);
	find_param(argv[2], "__core_bandwidth=", core_bandwidth);
	find_param(argv[2], "__Gamma=", Gamma);


	// ---------------------------------------------------------------
	// SETUP pp calculator
	
	typedef ppsc::hilbert_spaces::two_band_fermi_diag hilbert_space_type;
	typedef ppsc::hamiltonians::two_band_hubbard_xas<hilbert_space_type> hamiltonian_type;
	typedef ppsc::solver<hamiltonian_type> solver_type;

	hilbert_space_type hilbert_space;
	hilbert_space.init();

	solver_type imp_CBEXP(nt, ntau, beta, h, kt, nomp, hilbert_space, order); 	// core bath hybridization treated with hybridization expansion
	solver_type imp_CBCON(nt, ntau, beta, h, kt, nomp, hilbert_space, order);	// core bath hybridization included by VIE2 / convolution after time-stepping
	
	// -- Setup local Hamiltonian    
	initialize_hamiltonian(imp_CBEXP.hamiltonian, Hubbard_U, e_v, e_c, g_num, nt, h, beta, omega_in, sigma_pulse, mean_pulse);
	initialize_hamiltonian(imp_CBCON.hamiltonian, Hubbard_U, e_v, e_c, g_num, nt, h, beta, omega_in, sigma_pulse, mean_pulse);
	imp_CBEXP.update_hamiltonian();
	imp_CBCON.update_hamiltonian();
	
	// -- Setup core bath
	GREEN Hyb_core(nt, ntau, 1, FERMION), Hyb_core_cc(nt, ntau, 1, FERMION);	// initialized as zero
	cdvector core_hybridization;
	init_core_hybridization(core_hybridization, core_bandwidth, e_c, Gamma, h, nt);
	
	for (int tstp1 = 0; tstp1 <= nt; tstp1++){
		for (int tstp2 = tstp1; tstp2 <= nt; tstp2++) {
			Hyb_core.set_ret( tstp2 , tstp1 , core_hybridization[tstp2-tstp1] ); 	// set _retarded_ component non-zero
			Hyb_core.set_les( tstp1 , tstp2 , core_hybridization[tstp2-tstp1] );	// set _lesser_ component non-zero
		}
	}
	
	for ( int tstp = -1; tstp <= nt; tstp++ ) {	// set hermitian conjugate Hyb_core_cc
		ppsc::set_bwd_from_fwd(tstp, Hyb_core_cc, Hyb_core); 
	}
	
	//***** TEST OUTPUT *****/
	//Hyb_cc.print_to_file("./test_output_xas/test_Hyb_cc.txt"); 
	
	// ---------------------------------------------------------------------
	// SETUP single particle Green's functions and hybridizations	

	GREEN Delta(nt, ntau, 1, BOSON), Delta_cc(nt, ntau, 1, BOSON);	// initialized as zero
	GREEN G_loc_CBEXP(nt, ntau, 4, FERMION), G_loc_CBCON(nt, ntau, 4, FERMION);
	ppsc::gf_tstps_type gf_tstps_CBEXP, gf_tstps_CBCON; // vector of single particle GFs
	
	// ---------------------------------------------------------------------
	// INITIALIZATION
	
	if (SCREEN_HYB) {
		for(int tstp = -1; tstp <= nt; tstp++) {
			cntr::green_single_pole_XX_timestep(tstp, Delta, Omega_0, beta, h);	// calculate free bosonic propagator [-i < X(t) X(t') >] for one timestep
			Delta.smul(tstp, (gammatilde*gammatilde) / (2.0*Omega_0) ); // see my handwritten calculation on induced interaction by Holstein interaction U(t,t') = + (gamma^2 / 2) * D_free (t,t'), where D_free (t,t') = -i <X(t) X(t')> for real times and D_free(-i tau,-i tau') = - <X(-i tau) X(-i tau')> for imaginary times.
				// (1 / Omega_0) factor comes from def. of X = ( a + a^dag ) / sqrt(2) in NESSi [ vs. X = ( a + a^dag ) / sqrt(2*Omega_0) for physical displacement X ]
			
			// -- Use one of the two lines below (*) ; Delta is hermitian 
			Delta_cc.set_timestep(tstp,Delta); // (*)
			// ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta); // (*)
		}
	}

	
	// -- Solve atomic problem (here, Hyb_core = Hyb_core_cc = 0)
	imp_CBEXP.solve_atomic(); imp_CBCON.solve_atomic();
	imp_CBEXP.update_density_matrix(-1); imp_CBCON.update_density_matrix(-1);	// this routine is called when performing 'pp_step' function as well
	imp_CBEXP.hamiltonian.update_exp_vals(-1, imp_CBEXP.rho); imp_CBCON.hamiltonian.update_exp_vals(-1, imp_CBCON.rho);	// this routine is called when performing 'pp_step' function as well
	
	// -- Construct interactions and vertices (aka impurity Green's function of interest)
	ppsc::pp_ints_type pp_ints_CBEXP, pp_ints_CBCON;
	ppsc::gf_verts_type gf_verts;	// is the same for core-bath-hybridization-expansion & core-bath-convolution at the end
	if (SCREEN_HYB) {
		pp_ints_CBEXP = get_pp_ints(Delta, Delta_cc, Hyb_core, Hyb_core_cc, imp_CBEXP.hamiltonian.hil_);
	}
	else {
		//std::cout << "I am here - no screening!" << std::endl;
		pp_ints_CBEXP = get_pp_ints(Hyb_core, Hyb_core_cc, imp_CBEXP.hamiltonian.hil_);
	}
	pp_ints_CBCON = get_pp_ints_CBCON(Delta, Delta_cc, imp_CBCON.hamiltonian.hil_);
	
	gf_verts = get_gf_verts(hilbert_space);
	imp_CBEXP.update_diagrams(pp_ints_CBEXP, gf_verts); imp_CBCON.update_diagrams(pp_ints_CBCON, gf_verts);		// -- tell solver which diagrams to evaluate when performing 'get_spgf' or when doing 'pp_step'

	// -- Get single-particle GFs ('spgf')
	gf_tstps_CBEXP = imp_CBEXP.get_spgf(-1); gf_tstps_CBCON = imp_CBCON.get_spgf(-1);
	
		
	// ---------------------------------------------------------------------
	// MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
	
	
	GREEN G_loc_CBEXP_temp(nt, ntau, 4, FERMION), G_loc_CBCON_temp(nt, ntau, 4, FERMION);	// temporary Green's function
	// -- Update temporary GFs according to atomic/initialized solution aboset_t0_fromve
	for( auto idx : ppsc::range(0, gf_tstps_CBEXP.size()) ) {
	 	G_loc_CBEXP_temp.set_matrixelement(-1, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBEXP[idx]);
		G_loc_CBCON_temp.set_matrixelement(-1, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBCON[idx]);
	}

	// -- TEST OUTPUT
	// G_loc_temp.print_to_file("./test_output_xas/test_G_loc0.txt");		
	
	// -- Iteration loop
	for ( iter = 1 ; iter <= itermax; iter++) {
		
		// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
		//Delta.smul(-1, -0.4*(gamma*gamma/2.0));
		//ppsc::set_bwd_from_fwd(-1, Delta_cc, Delta); 
			
		// -- Solve pseudo particle problem --> this is probably only needed if we change Delta (see above); otherwise just has to be called once before performing 'pp_step' or 'get_spgf' (see 'Initialization' above)
		//imp_CBEXP.update_diagrams(pp_ints_CBEXP, gf_verts); imp_CBCON.update_diagrams(pp_ints_CBCON, gf_verts);
		
		// -- Do the (Matsubara) time step
		imp_CBEXP.pp_step(-1); imp_CBCON.pp_step(-1);
		
		// -- Get single-particle GFs ('spgf')
		gf_tstps_CBEXP = imp_CBEXP.get_spgf(-1); gf_tstps_CBCON = imp_CBCON.get_spgf(-1);
		
		// -- Update G_loc
		for( auto idx : ppsc::range(0, gf_tstps_CBEXP.size()) ) {
		 	G_loc_CBEXP.set_matrixelement(-1, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBEXP[idx]);
			G_loc_CBCON.set_matrixelement(-1, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBCON[idx]);
		}
		
		// -- Calculate change with last solution
		itererr_CBEXP = cntr::distance_norm2(-1, G_loc_CBEXP, G_loc_CBEXP_temp); itererr_CBCON = cntr::distance_norm2(-1, G_loc_CBCON, G_loc_CBCON_temp);
		if (DEBUG) { std::cout << "Matsubara-iter:  " << iter << " err_CBEXP: " << itererr_CBEXP << " err_CBCON: " << itererr_CBCON << std::endl; }		
		
		// -- Update G_temp's
		G_loc_CBEXP_temp.set_timestep(-1, G_loc_CBEXP); G_loc_CBCON_temp.set_timestep(-1, G_loc_CBCON);
				
		// -- Check for convergence
		if ( std::sqrt(itererr_CBEXP*itererr_CBEXP + itererr_CBCON*itererr_CBCON) <= errmax ) {
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
		imp_CBEXP.init_real_time(); imp_CBCON.init_real_time();	// should do the same as 'imp.extrapolate_timestep(-1)'
	
		// -- Iteration loop
		for ( iter = 1 ; iter <= itermax; iter++) {
			
			// -- Solve pseudo particle problem
			//	 --> this is probably only needed if we change Delta (see below);
			//	otherwise just has to be called once before performing 'pp_step' or 'get_spgf' (see 'Initialization' above)
			//imp_CBEXP.update_diagrams(pp_ints_CBEXP, gf_verts); imp_CBCON.update_diagrams(pp_ints_CBCON, gf_verts);
			
			// -- Do the (Bootstrapping) time step; for kt = SolveOrder, the Dyson equation is solved for all time steps 0, 1, ... kt together
			imp_CBEXP.pp_step(kt); imp_CBCON.pp_step(kt); 		// kt = SolveOrder
		
			for (int tstp = 0; tstp <= kt; tstp++) {
			
				// -- Get single-particle GFs ('spgf')
				gf_tstps_CBEXP = imp_CBEXP.get_spgf(tstp); gf_tstps_CBCON = imp_CBCON.get_spgf(tstp);
				
				// -- Update G_loc
				for( auto idx : ppsc::range(0, gf_tstps_CBEXP.size()) ) {
				 	G_loc_CBEXP.set_matrixelement(tstp, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBEXP[idx]);
					G_loc_CBCON.set_matrixelement(tstp, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBCON[idx]);
				}
				
				// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
				//Delta.smul(tstp, -0.4*(gamma*gamma/2.0));
				//ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta); 
					
			}
			
			// -- Calculate change with last solution for time t = kt only !!!!
			itererr_CBEXP = cntr::distance_norm2(kt, G_loc_CBEXP, G_loc_CBEXP_temp); itererr_CBCON = cntr::distance_norm2(kt, G_loc_CBCON, G_loc_CBCON_temp);
			if (DEBUG) { std::cout << "Bootstrapping-iter:  " << iter << " err_CBEXP: " << itererr_CBEXP << " err_CBCON: " << itererr_CBCON << std::endl; }
			
			// -- Update G_temp's
			// (here, in the bootrstrapping phase, we only need kt timeslice of G_temp's
			//	--> one could reduce memory here or define only kt time slice for temporary Green's functions)
			for (int tstp = 0; tstp <= kt; tstp++) {
				G_loc_CBEXP_temp.set_timestep(tstp, G_loc_CBEXP); G_loc_CBCON_temp.set_timestep(tstp, G_loc_CBCON);
			}
			
			// -- Check for convergence
			if ( std::sqrt(itererr_CBEXP*itererr_CBEXP + itererr_CBCON*itererr_CBCON)  <= errmax ) {
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
			imp_CBEXP.extrapolate_timestep(tstp-1); imp_CBCON.extrapolate_timestep(tstp-1);
			
			// -- Iteration loop
			for (iter = 1; iter <= itermax ; iter++) {
			
				// -- Solve pseudo particle problem
				//	 --> this is probably only needed if we change Delta (see below);
				//	otherwise just has to be called once before performing 'pp_step' or 'get_spgf' (see 'Initialization' above)
				//imp_CBEXP.update_diagrams(pp_ints_CBEXP, gf_verts); imp_CBCON.update_diagrams(pp_ints_CBCON, gf_verts);
				
				// -- Do the real-time step
				imp_CBEXP.pp_step(tstp); imp_CBCON.pp_step(tstp);
				
				// -- Get single-particle GFs ('spgf')
				gf_tstps_CBEXP = imp_CBEXP.get_spgf(tstp); gf_tstps_CBCON = imp_CBCON.get_spgf(tstp);
				
				// -- Update G_loc
				for( auto idx : ppsc::range(0, gf_tstps_CBEXP.size()) ) {
				 	G_loc_CBEXP.set_matrixelement(tstp, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBEXP[idx]);
					G_loc_CBCON.set_matrixelement(tstp, gf_verts[idx].idx1, gf_verts[idx].idx2, gf_tstps_CBCON[idx]);
				}
				
				// -- !!! WE COULD CHANGE DELTA HERE WITHOUT INVOKING 'get_pp_ints' AGAIN SINCE REF. GIVEN TO 'pp_ints' ABOVE IS STILL ACTIVE... !!!
				//Delta.smul(tstp, -0.4*(gamma*gamma/2.0));
				//ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);
				
				// -- Calculate change with last solution
				itererr_CBEXP = cntr::distance_norm2(tstp, G_loc_CBEXP, G_loc_CBEXP_temp); itererr_CBCON = cntr::distance_norm2(tstp, G_loc_CBCON, G_loc_CBCON_temp);
				if (DEBUG) { std::cout << "Time-stepping-iter:  " << iter << " err_CBEXP: " << itererr_CBEXP << " err_CBCON: " << itererr_CBCON << std::endl; }
		
				// -- Update G_temp's
				G_loc_CBEXP_temp.set_timestep(tstp, G_loc_CBEXP); G_loc_CBCON_temp.set_timestep(tstp, G_loc_CBCON);

				// -- Check for convergence
				if ( (std::sqrt(itererr_CBEXP*itererr_CBEXP + itererr_CBCON*itererr_CBCON) <= errmax) && iter > 2) {
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
	//G_loc.print_to_file("./test_output_xas/test_G_loc.txt");
	
	
	// ---------------------------------------------------------------------
      	// ADD CORE BATH TO SOLUTION (G_loc_CBCON)
      	
      	GREEN G_loc_CBCON_new(nt, ntau, 4, FERMION);
      	
      	// -- Set up integration kernel F = - G_loc_CBCON * Hyb_cc_4 (F_cc = - Hyb_cc_4 * G_loc_CBCON) ; Q is equal to G_loc_CBCON
      	GREEN F(nt, ntau, 4, BOSON), F_cc(nt, ntau, 4, BOSON);
      	GREEN Hyb_core_4(nt, ntau, 4, FERMION);
	Hyb_core_4.set_matrixelement(3, 3, Hyb_core);	// core bath hybridization is spin symmetric
	Hyb_core_4.set_matrixelement(2, 2, Hyb_core);	// core bath hybridization is spin symmetric
	
      	cntr::convolution(F, G_loc_CBCON, G_loc_CBCON, Hyb_core_4, Hyb_core_4, beta, h, kt);
      	cntr::convolution(F_cc, Hyb_core_4, Hyb_core_4, G_loc_CBCON, G_loc_CBCON, beta, h, kt);
      	for (int tstp = -1; tstp <= nt; tstp++ ) {
      		F.smul(tstp, -1.0);
      		F_cc.smul(tstp, -1.0);
      	}
      	
      	// -- Matsubara branch
      	cntr::vie2_mat(G_loc_CBCON_new, F, F_cc, G_loc_CBCON, beta, CNTR_MAT_FIXPOINT, kt);
      	     	
      	// -- Bootstrapping (time <= kt)
	cntr::vie2_start(G_loc_CBCON_new, F, F_cc, G_loc_CBCON, beta, h, kt);

      	// -- Real time propagation
      	for(int tstp = kt; tstp <= nt; tstp++ ) {
      		cntr::vie2_timestep(tstp, G_loc_CBCON_new, F, F_cc, G_loc_CBCON, beta, h, kt);
      	}
      	
      	// -- Update G_loc_CBCON
      	for (int tstp = -1; tstp <= nt; tstp++ ) {
		G_loc_CBCON.set_timestep(tstp, G_loc_CBCON_new);
      	}
      	
      	
      	// ---------------------------------------------------------------------
      	// CALCULATE XAS SIGNAL
      	
      	// -- Calculate xas_integrand = 2 * II * exp(-II * omega_in * t) * s(t) * < P_in (t) > / g 
      	// -- with P_in(t) , a quantity depending on omega_in (!!!)
      	
      	cdmatrix rho(4,4);
      	std::vector<cdouble> xas_integrand_CBEXP(nt+1);	// only for real time arguments
      	std::vector<cdouble> xas_integrand_CBCON_no_corebath(nt+1);	// only for real time arguments
      	std::vector<cdouble> xas_integrand_CBCON(nt+1);	// only for real time arguments
      	std::vector<cdouble> P_in_CBCON(nt+1);	// only for real time arguments
      	      	      	
      	for( int tstp = 0; tstp <= nt; tstp++ ) {
      	
      	      	G_loc_CBCON.density_matrix(tstp, rho);
      	      	P_in_CBCON[tstp] = rho(2,0) + rho(3,1);
      	      	
      		xas_integrand_CBEXP[tstp] = h * 2 * II * std::exp(-II * omega_in * double(tstp) * h) * ( imp_CBEXP.hamiltonian.P_in_exp[tstp+1] ) * (imp_CBEXP.hamiltonian.probe_pulse_env[tstp+1]) / g_num ;
      		xas_integrand_CBCON[tstp] = h * 2 * II * std::exp(-II * omega_in * double(tstp) * h) * ( P_in_CBCON[tstp] ) * (imp_CBCON.hamiltonian.probe_pulse_env[tstp+1]) / g_num ;
 		xas_integrand_CBCON_no_corebath[tstp] =  h * 2 * II * std::exp(-II * omega_in * double(tstp) * h) * ( imp_CBCON.hamiltonian.P_in_exp[tstp+1]  ) * (imp_CBCON.hamiltonian.probe_pulse_env[tstp+1]) / g_num;
      	}
      	
      	// -- Integrate xas_integrand over [0,nt*h] to find I_XAS(omega_in)
      	integration::Integrator<double> gregory_integration(kt);
      	
      	cdouble I_XAS_CBEXP = gregory_integration.integrate(xas_integrand_CBEXP, nt);
       	cdouble I_XAS_CBCON = gregory_integration.integrate(xas_integrand_CBCON, nt);
      	cdouble I_XAS_CBCON_no_corebath = gregory_integration.integrate(xas_integrand_CBCON_no_corebath, nt);    	

	// ---------------------------------------------------------------------
	// PRINT OUT / SAVING RESULTS FOR I_XAS(omega_in)  &  < P_in >_{g_num} / g_num 
	
	std::ofstream I_XAS_file, P_in_CBEXP_file, P_in_CBCON_file, P_in_CBCON_no_corebath_file;
	
	// -- create directories
	std::ostringstream output_dir, output_file_xas, output_file_time_CBEXP, output_file_time_CBCON, output_file_time_CBCON_no_corebath;
	if (SCREEN_HYB) {
		output_dir << output_dir_read_in << "Omega" << Omega_0 << "_gammatilde" << gammatilde << "_Gamma" << Gamma << "_ev" << e_v << "_ec" << e_c << "_U" << Hubbard_U << "_ntau" << ntau << "_nt" << nt << "_beta" << beta << "_h" << h << "_pulse-mean" << mean_pulse << "_pulse-sigma" << sigma_pulse << "_corebandwidth" << core_bandwidth << "_gnum" << g_num; 
	}
	else {
		output_dir << output_dir_read_in << "NOSCREEN" << "_Gamma" << Gamma << "_ev" << e_v << "_ec" << e_c << "_U" << Hubbard_U << "_ntau" << ntau << "_nt" << nt << "_beta" << beta << "_h" << h << "_pulse-mean" << mean_pulse << "_pulse-sigma" << sigma_pulse << "_corebandwidth" << core_bandwidth << "_gnum" << g_num ;
	}
	
	output_file_xas << "/" << "I_XAS_34.txt" ;
	output_file_time_CBEXP << "/" << "P_in_exp_CBEXP_34.txt" ;
	output_file_time_CBCON << "/" << "P_in_exp_CBCON_34.txt" ;
	output_file_time_CBCON_no_corebath << "/" << "P_in_exp_CBCON_no_corebath_34.txt" ;
		
	std::string tmp = output_dir.str();
	const char * output_dir_str = tmp.c_str();
	std::cout << "output_dir_str: " << output_dir_str << std::endl;
	struct stat sb;
	if ( stat(output_dir_str, &sb) != 0) {
		// create directory 
		const int create_dir = mkdir( output_dir_str, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);	// S_ISGID
		if (create_dir == -1) 	std::cerr << "Error :  " << std::strerror(errno) << std::endl;
    		else 	std::cout << "Directory Created: " << output_dir.str() << std::endl;
	}
	
	std::string output_filename_xas = output_dir.str() + (output_file_xas.str());
	std::string output_filename_CBEXP = output_dir.str() + (output_file_time_CBEXP.str());
	std::string output_filename_CBCON = output_dir.str() + (output_file_time_CBCON.str());
	std::string output_filename_CBCON_no_corebath = output_dir.str() + (output_file_time_CBCON_no_corebath.str());
	
	
	// -- Check whether output file already exists ...
	if ( ! file_exists(output_filename_xas) && ! file_exists(output_filename_CBEXP) && ! file_exists(output_filename_CBCON) && ! file_exists(output_filename_CBCON_no_corebath) ) {
	
		// -- XAS signal
		I_XAS_file.open(output_filename_xas);
		I_XAS_file << "############ XAS signal ############ \n##\n";
		if (!SCREEN_HYB) { I_XAS_file << "## NO SCREENING \n##\n"; }
		I_XAS_file << "## Omega = " << Omega_0 << ", gammatilde = " << gammatilde << ", Gamma =" << Gamma << ", e_v = " << e_v << ", e_c = " << e_c << ", U = " << Hubbard_U << "\n";
		I_XAS_file << "## ntau = " << ntau << ", nt = " << nt << ", beta = " << beta << ", h = " << h << ", g_num = " << g_num << "\n";
		I_XAS_file << "## pulse-mean = " << mean_pulse << ", pulse-sigma = " << sigma_pulse << ", core-bandwidth = " << core_bandwidth  << "\n##\n";
		I_XAS_file << "# omega_in" << "\t" << "I_XAS_CBEXP" << "\t" << "I_XAS_CBCON" << "\t" << "I_XAS_CBCON_no_corebath" ;
		
		// -- P_in (t)
		//
		P_in_CBEXP_file.open(output_filename_CBEXP);
		P_in_CBEXP_file << "############ < P_in_CBEXP >_g / g ############ \n##\n";
		if (!SCREEN_HYB) { P_in_CBEXP_file << "## NO SCREENING \n##\n"; }
		P_in_CBEXP_file << "## Omega = " << Omega_0 << ", gammatilde = " << gammatilde << ", Gamma =" << Gamma << ", e_v = " << e_v << ", e_c = " << e_c << ", U = " << Hubbard_U << "\n";
		P_in_CBEXP_file << "## ntau = " << ntau << ", nt = " << nt << ", beta = " << beta << ", h = " << h << ", g_num = " << g_num << "\n";
		P_in_CBEXP_file << "## pulse-mean = " << mean_pulse << ", pulse-sigma = " << sigma_pulse << ", core-bandwidth = " << core_bandwidth  << "\n##\n";
		P_in_CBEXP_file << "# omega_in" << "\t" << "<P_in>/g (t=0) REAL IMAG" << "\t" << "<P_in>/g (h) REAL IMAG" << "\t" << "<P_in>/g (2*h) REAL IMAG" << "\t" << "..." ;
		//
		P_in_CBCON_file.open(output_filename_CBCON);
		P_in_CBCON_file << "############ < P_in_CBCON >_g / g ############ \n##\n";
		
		if (!SCREEN_HYB) { P_in_CBCON_file << "## NO SCREENING \n##\n"; }
		
		P_in_CBCON_file << "## Omega = " << Omega_0 << ", gammatilde = " << gammatilde << ", Gamma =" << Gamma << ", e_v = " << e_v << ", e_c = " << e_c << ", U = " << Hubbard_U << "\n";
		P_in_CBCON_file << "## ntau = " << ntau << ", nt = " << nt << ", beta = " << beta << ", h = " << h << ", g_num = " << g_num << "\n";
		P_in_CBCON_file << "## pulse-mean = " << mean_pulse << ", pulse-sigma = " << sigma_pulse << ", core-bandwidth = " << core_bandwidth  << "\n##\n";
		P_in_CBCON_file << "# omega_in" << "\t" << "<P_in>/g (t=0) REAL IMAG" << "\t" << "<P_in>/g (h) REAL IMAG" << "\t" << "<P_in>/g (2*h) REAL IMAG" << "\t" << "..." ;
		//
		P_in_CBCON_no_corebath_file.open(output_filename_CBCON_no_corebath);
		P_in_CBCON_no_corebath_file << "############ < P_in_CBCON_no_corebath >_g / g ############ \n##\n";
		if (!SCREEN_HYB) { P_in_CBCON_no_corebath_file << "## NO SCREENING \n##\n"; }
		P_in_CBCON_no_corebath_file << "## Omega = " << Omega_0 << ", gammatilde = " << gammatilde << ", Gamma =" << Gamma << ", e_v = " << e_v << ", e_c = " << e_c << ", U = " << Hubbard_U << "\n";
		P_in_CBCON_no_corebath_file << "## ntau = " << ntau << ", nt = " << nt << ", beta = " << beta << ", h = " << h << ", g_num = " << g_num << "\n";
		P_in_CBCON_no_corebath_file << "## pulse-mean = " << mean_pulse << ", pulse-sigma = " << sigma_pulse << ", core-bandwidth = " << core_bandwidth  << "\n##\n";
		P_in_CBCON_no_corebath_file << "# omega_in" << "\t" << "<P_in>/g (t=0) REAL IMAG" << "\t" << "<P_in>/g (h) REAL IMAG" << "\t" << "<P_in>/g (2*h) REAL IMAG" << "\t" << "..." ;
		
		
	}
	
	else {
		I_XAS_file.open(output_filename_xas, std::ios_base::app); 
		
		P_in_CBEXP_file.open(output_filename_CBEXP, std::ios_base::app); 
		P_in_CBCON_file.open(output_filename_CBCON, std::ios_base::app); 
		P_in_CBCON_no_corebath_file.open(output_filename_CBCON_no_corebath, std::ios_base::app); 
		
		
	}
	
	// -- Write data to file
	if ( I_XAS_file.is_open() && P_in_CBEXP_file.is_open() && P_in_CBCON_file.is_open() && P_in_CBCON_no_corebath_file.is_open() ) {
		
		// -- XAS signal
		I_XAS_file << "\n" << omega_in << "\t" << I_XAS_CBEXP.imag() << "\t" << I_XAS_CBCON.imag() <<  "\t" << I_XAS_CBCON_no_corebath.imag();
		I_XAS_file.close();
		
		// -- P_in (t)
		P_in_CBEXP_file << "\n" << omega_in;
		P_in_CBCON_file << "\n" << omega_in;
		P_in_CBCON_no_corebath_file << "\n" << omega_in;
		for ( int tstp = 0; tstp<= nt; tstp++ ) {
			P_in_CBEXP_file << "\t" << (imp_CBEXP.hamiltonian.P_in_exp[tstp+1]).real() / g_num  << "\t" << (imp_CBEXP.hamiltonian.P_in_exp[tstp+1]).imag() / g_num ;
			P_in_CBCON_file << "\t" << (P_in_CBCON[tstp]).real() / g_num << "\t" << (P_in_CBCON[tstp]).imag() / g_num ;
			P_in_CBCON_no_corebath_file << "\t" << (imp_CBCON.hamiltonian.P_in_exp[tstp+1]).real() / g_num  << "\t" << (imp_CBCON.hamiltonian.P_in_exp[tstp+1]).imag() / g_num ;
		}
		P_in_CBEXP_file.close();
		P_in_CBCON_file.close();
		P_in_CBCON_no_corebath_file.close();
		
		
	}
	else 	std::cout << "Unable to open file(s)..." << std::endl;
	
	
	// ---------------------------------------------------------------------
	/* PRINT OUT / SAVING RESULTS FOR OBSERVABLES DEFINED IN HAMILTONIAN
	

	std::string observable_file_CBEXP = "/observables_CBEXP_5.h5", observable_file_CBCON = "/observables_CBCON_5.h5";
	hid_t file_id_CBEXP, file_id_CBCON;
	hid_t group_id_CBEXP, group_id_CBCON;
	

	// check whether files exist; if yes, only append/create new groups for omega_in
	if ( ! file_exists( output_dir.str() + observable_file_CBEXP ) && ! file_exists( output_dir.str() + observable_file_CBCON) ) {
			std::cout << "Create new hdf5 files." << std::endl;
			file_id_CBEXP = open_hdf5_file(output_dir.str() + observable_file_CBEXP);
			file_id_CBCON = open_hdf5_file(output_dir.str() + observable_file_CBCON);
	
	}
	else {
			file_id_CBEXP = H5Fopen( (output_dir.str() + observable_file_CBEXP).c_str() , H5F_ACC_RDWR, H5P_DEFAULT);
			file_id_CBCON = H5Fopen( (output_dir.str() + observable_file_CBCON).c_str() , H5F_ACC_RDWR, H5P_DEFAULT);
	}

	
	group_id_CBEXP = create_group(file_id_CBEXP, "omega_in=" + std::to_string(omega_in));
	group_id_CBCON = create_group(file_id_CBCON, "omega_in=" + std::to_string(omega_in));
	
	imp_CBEXP.hamiltonian.store(group_id_CBEXP);
	imp_CBCON.hamiltonian.store(group_id_CBCON);
	
	H5Gclose(group_id_CBEXP); H5Gclose(group_id_CBCON);
	close_hdf5_file(file_id_CBEXP); close_hdf5_file(file_id_CBCON);
	
	*/
	
	
	
	return 0;

}


