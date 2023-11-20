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

#define GREEN cntr::herm_matrix<double>
// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & Delta_up,
			       ppsc::gf_type & Delta_up_cc,
			       HILB & hil_) {

  // spin (u)p/(d)own and (c)reation/(a)nihilation operators

  ppsc::operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];
  
  int fermion=-1, fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  pp_ints.push_back(ppsc::pp_int_type(Delta_up,    cuc, cua, fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta_up_cc, cua, cuc, fermion, bwd)); // spin up bwd
  
  pp_ints.push_back(ppsc::pp_int_type(Delta_up,    cdc, cda, fermion, fwd)); // spin do fwd // assuming spin sym
  pp_ints.push_back(ppsc::pp_int_type(Delta_up_cc, cda, cdc, fermion, bwd)); // spin do bwd // assuming spin sym

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
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, cua, cuc)); // spin up
  return gf_verts;
}
 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// solve (1+G1c)*glatt = G  (G1c=G*Delta)
void get_glatt(int tstp,GREEN &g_latt,GREEN &G1,GREEN &G1c,GREEN &Gloc,double beta,double h,int kt){
	cntr::vie2_timestep(tstp,g_latt,G1c,G1,Gloc,beta,h,kt);
		// if(tstp==-1){
		// 	cntr::vie2_mat_sweep(g_latt,G1c,G1,Gloc,beta,integration::I<double>(kt),NSWEEP);
		// }else if(tstp<=kt){
		// 	cntr::vie2_start(g_latt,G1c,G1,Gloc,integration::I<double>(tstp),beta,h);
		// }else{ 
		// 	cntr::vie2_timestep(tstp,g_latt,G1c,G1,Gloc,integration::I<double>(kt),beta,h);		   
		// }
	}

// first update G1c=G*Delta, G1=Delta*G, then solve (1+G1c)*glatt = G 
void get_glatt(int tstp,GREEN &g_latt,GREEN &Delta,GREEN &G1,
					GREEN &G1c,GREEN &Gloc,double beta,double h,int kt){
	{

	int n;
	int n1=(tstp<=kt && tstp>=0 ? 0 : tstp);
	int n2=(tstp<=kt && tstp>=0 ? kt : tstp);
	std::cout << "Get glatt inside " << n1 << " " <<  n2 << std::endl;
	for(n=n1;n<=n2;n++){
	cntr::convolution_timestep(n,G1c,Gloc,Delta,integration::I<double>(kt),beta,h);
	cntr::convolution_timestep(n,G1,Delta,Gloc,integration::I<double>(kt),beta,h);
	}

	for(n=n1;n<=n2;n++){
		get_glatt(n,g_latt,G1,G1c,Gloc,beta,h,kt);
	}
	}
}

// ---------------------------------------------------------------------
// get Sigma (approxiomatic method)
// ---------------------------------------------------------------------
void get_Sigma(GREEN &Sigma,GREEN &F,GREEN &Fcc,GREEN &Z,GREEN Q,double beta,double h, int kt,int nt){
	int tstp;
	
	// get F=id/dt Z - 1
	cntr::deriv1_matsubara(F,Z,integration::I<double>(kt),beta);
	for(tstp=0;tstp<=nt;tstp++) cntr::deriv1_timestep(tstp,F,Z,Z,integration::I<double>(kt),beta,h);
	// get Fcc=-id/dt' Z - 1
	cntr::deriv1_matsubara(Fcc,Z,integration::I<double>(kt),beta); // note: matsubara: id/dt=-id/dt'
	for(tstp=0;tstp<=nt;tstp++) cntr::deriv2_timestep(tstp,Fcc,Z,Z,integration::I<double>(kt),beta,h);
	// get Q=-id/dt'F=id/dt Fcc
	cntr::deriv1_matsubara(Q,F,integration::I<double>(kt),beta);
	for(tstp=0;tstp<=nt;tstp++) cntr::deriv2_timestep(tstp,Q,F,Fcc,integration::I<double>(kt),beta,h);
	// solve (1+F)*Sigma=Q
	cntr::vie2(Sigma,F,Fcc,Q,integration::I<double>(kt),beta,h);
	
}

void get_Sigma_tstp(int iter,int tstp,GREEN &Sigma,GREEN &F,GREEN &Fcc,GREEN &Z,GREEN Q,double beta,double h, int kt,int nt){
	int n;
	
	if (tstp==-1){
		// get F=id/dt Z - 1
		cntr::deriv1_matsubara(F,Z,integration::I<double>(kt),beta);
		// get Fcc=-id/dt' Z - 1
		cntr::deriv1_matsubara(Fcc,Z,integration::I<double>(kt),beta); // note: matsubara: id/dt=-id/dt'
		// get Q=-id/dt'F=id/dt Fcc
		cntr::deriv1_matsubara(Q,F,integration::I<double>(kt),beta);
		//cntr::deriv1_matsubara(Q,F,integration::I<double>(kt),beta);
		if(iter<2){	
		std::ostringstream nameQ;
		   nameQ << "ppsc_data_Q_" << -1 <<"_"<<iter<< ".h5";
		   Q.write_to_hdf5(nameQ.str().c_str(),"Q");
		}
	}else{
		// get F=id/dt Z - 1
		for(n=0;n<=tstp;n++) cntr::deriv1_timestep(n,F,Z,Z,integration::I<double>(kt),beta,h);
		// get Fcc=-id/dt' Z - 1
		for(n=0;n<=tstp;n++) cntr::deriv2_timestep(n,Fcc,Z,Z,integration::I<double>(kt),beta,h);
		// get Q=-id/dt'F=id/dt Fcc
		for(n=0;n<=tstp;n++) cntr::deriv2_timestep(n,Q,F,Fcc,integration::I<double>(kt),beta,h);
	}	
	// solve (1+F)*Sigma=Q
	cntr::vie2(Sigma,F,Fcc,Q,integration::I<double>(kt),beta,h);
}

void get_Sigma_wHartree(GREEN &Sigma,GREEN &F,GREEN &Fcc,GREEN &Z,GREEN Q,double beta,double h, int kt,int nt){
	int tstp;
    cntr::function<double> Sigma_H(nt);
	
	// get F=id/dt Z - 1
	cntr::deriv1_matsubara(F,Z,integration::I<double>(kt),beta);
	for(tstp=0;tstp<=nt;tstp++) cntr::deriv1_timestep(tstp,F,Z,Z,integration::I<double>(kt),beta,h);
	// get Fcc=-id/dt' Z - 1
	cntr::deriv1_matsubara(Fcc,Z,integration::I<double>(kt),beta); // note: matsubara: id/dt=-id/dt'
	for(tstp=0;tstp<=nt;tstp++) cntr::deriv2_timestep(tstp,Fcc,Z,Z,integration::I<double>(kt),beta,h);
	// get Q=-id/dt'F=id/dt Fcc
	cntr::deriv1_matsubara(Q,F,integration::I<double>(kt),beta);
	for(tstp=0;tstp<=nt;tstp++) cntr::deriv2_timestep(tstp,Q,F,Fcc,integration::I<double>(kt),beta,h);
	// solve (1+F)*Sigma=Q
	cntr::vie2(Sigma,F,Fcc,Q,integration::I<double>(kt),beta,h);
	
}


// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int itermax, iter_rtime, nt, ntau, kt, iter, iter_n, tstp, order, nomp;

  int store_pp, read_eq_sym, read_rt_sym, read_state_from_file;
  double linear_mixing;
  
  double beta, h, errmax, dmfterr, hyberr, dmfterr_equil, hyberr_eq;
  bool matsubara_converged = false, hartree=false;

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
      find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__read_state_from_file=", read_state_from_file);
      find_param(argv[1], "__hartree=", hartree);

      if (order > 1) {
        cout << "using an OCA impurity solver" << endl;
        find_param(argv[1], "__nomp=", nomp);
      } else {
	nomp = 1;
      }
    }
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

      find_param(argv[1], "__mu=", imp.hamiltonian.mu);

      find_param_tvector(argv[1], "__U=", imp.hamiltonian.U, nt);
      find_param(argv[1], "__U0=", imp.hamiltonian.U[0]);

      find_param_tvector(argv[1], "__eps=", imp.hamiltonian.eps_up, nt);
      find_param_tvector(argv[1], "__eps=", imp.hamiltonian.eps_do, nt);

      imp.update_hamiltonian();
    }
    // ---------------------------------------------------------------------
    // -- setup single particle greens functions and hybridizations

    GREEN Gloc_up(nt, ntau, 1, -1);    
    GREEN Delta_up(nt, ntau, 1, -1);
    GREEN Delta_up_cc(nt, ntau, 1, -1);

    // -- setup self-energy and auxilary functions
    GREEN Sigma(nt, ntau, 1, -1);    
    GREEN Z(nt, ntau, 1, -1);    
    GREEN G1(nt, ntau, 1, -1);    
    GREEN G1c(nt, ntau, 1, -1);    

	GREEN Q(nt,ntau,1,-1);
	GREEN F(nt,ntau,1,-1);
	GREEN Fcc(nt,ntau,1,-1);

    GREEN Delta_n(nt,ntau,1,-1);
    GREEN A(nt,ntau,1,-1);
    GREEN Acc(nt,ntau,1,-1);
    GREEN Gn(nt,ntau,1,-1);
    GREEN dGn(nt,ntau,1,-1);
    GREEN dGn2(nt,ntau,1,-1);
    GREEN dGncc(nt,ntau,1,-1);
    GREEN tempM(nt,ntau,1,-1);
    GREEN tempMcc(nt,ntau,1,-1);
    GREEN DdGnD(nt,ntau,1,-1);
	GREEN M(nt,ntau,1,-1);
    GREEN Mcc(nt,ntau,1,-1);
    GREEN MD(nt,ntau,1,-1);
    GREEN dDelta(nt,ntau,1,-1);
    GREEN tempdDelta(nt,ntau,1,-1);
    GREEN tmp_zero(nt,ntau,1,-1);
    // ---------------------------------------------------------------------
    
    if(read_state_from_file) {

      std::string filename = "data_ppsc.h5";
      std::cout << "--> Reading state from file: " << filename << std::endl;

      Gloc_up.read_from_hdf5(filename.c_str(), "g");
      Delta_up.read_from_hdf5(filename.c_str(), "d");
      Delta_up_cc.read_from_hdf5(filename.c_str(), "dcc");

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
      GREEN gtmp(-1, ntau, 1, -1);

      if(!read_state_from_file) {
		imp.solve_atomic();
      } else {
		gtmp.set_timestep(-1, Gloc_up); // this avoids one iteration
      }


      for (iter = 1; iter <= itermax; iter++) {

	// -- Construct interactions and verticies
	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta_up, Delta_up_cc, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	
	// -- Solve pseudo particle problem
	imp.update_diagrams(pp_ints, gf_verts);

	if(iter == 1 && read_eq_sym && !imp.has_symmetries()) {
	  std::string filename = "sym_eq.txt";
	  imp.read_symmetries(filename);
	}

	imp.pp_step(-1);

	// -- get spgf
	ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);

	//Gloc_up.set_timestep(-1, gf_tstps[0]);	

	// -- linear mixing in Gloc
	ppsc::gf_tstp_type gloc_old(-1, ntau, 1);
	ppsc::gf_tstp_type gloc_mix(-1, ntau, 1);

	gloc_mix.clear();
	Gloc_up.get_timestep(-1, gloc_old);
	gloc_mix.incr(gloc_old, linear_mixing);
	gloc_mix.incr(gf_tstps[0], 1.0 - linear_mixing);
	Gloc_up.set_timestep(-1, gloc_mix);
	
	// -- Check error
	dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc_up);
        gtmp.set_timestep(-1, Gloc_up);

	// -- Update Hybridization
	cntr::herm_matrix_timestep<double> tmp;
	Gloc_up.get_timestep(-1, tmp);
	Delta_up.set_timestep(-1, tmp);
	ppsc::set_bwd_from_fwd(-1, Delta_up_cc, Delta_up);
	
	cout << "ITER:  " << iter << " err: " << dmfterr_equil << endl;
	
	// calculate Z
	Z.clear();
	get_glatt(-1,Z,Delta_up,G1,G1c,Gloc_up,beta,h,kt);
	// approximate Sigma, self-energy
	Sigma.clear();
	get_Sigma_tstp(iter,-1,Sigma,F,Fcc,Z,Q,beta,h,kt,nt);
	
	// -- make a guess for the hybridization function
	Delta_n.clear();
	Delta_n.set_timestep(-1,Delta_up);

	Delta_n.smul(-1,0.8);
	
	// -- calculate hybridization approximately
	iter_n=1;
	hyberr_eq=1;
	int iter_n_max=7;
	int kt_mat=0;

	while(hyberr_eq>errmax*100000){
		
		// -- calculate Gn=Z+Z*Delta*Glatt
		cntr::convolution_timestep(-1,A,Z,Delta_n,integration::I<double>(kt),beta,h);
		cntr::convolution_timestep(-1,Acc,Delta_n,Z,integration::I<double>(kt),beta,h);
		cntr::convolution_timestep(-1,Gn,A,Acc,Gloc_up,Gloc_up,integration::I<double>(kt),beta,h);

		Gn.incr_timestep(-1,Z,1.0);

		dGn.clear();
		// -- define dG_n = G_loc - G_n
		dGn.set_timestep(-1,Gloc_up);
		dGn.incr_timestep(-1,Gn,-1.0);
		
		// -----------test---------------
		tmp_zero.clear();
		std::cout<<"dGn err total "<< cntr::distance_norm2(-1,dGn,tmp_zero)<<std::endl;
/*	    if((iter_n<=iter_n_max)&&(iter<=2)){	
	    std::ostringstream namedGn;
 	       namedGn << "ppsc_data_dGn_" << -1 <<"_"<<iter<<"_" << iter_n <<  ".h5";
		   dGn.write_to_hdf5(namedGn.str().c_str(),"dGn");
	    }
*/
		// ----------end test-------------

/*		// -- get DdGn=id/dt dGn
		cntr::deriv1_matsubara(tempM,dGn,integration::I<double>(kt_mat),beta);
*/
		// ---------test------------------
		tmp_zero.clear();
		std::cout<<"DdGn err total "<< cntr::distance_norm2(-1,tempM,tmp_zero)<<std::endl;
/*	    if((iter_n<=iter_n_max)&&(iter<=2)){	
	    std::ostringstream nameDdGn;
 	       nameDdGn << "ppsc_data_DdGn_" << -1 <<"_" << iter<<"_" << iter_n <<  ".h5";
	       tempM.write_to_hdf5(nameDdGn.str().c_str(),"DdGn");
	    }
*/
		// ----------end test-------------

//		// -- get dGnD=-id/dt' dGn
//		cntr::deriv1_matsubara(tempMcc,dGn,integration::I<double>(kt),beta);
//		tempMcc.set_mat_herm();
//		// ---------test------------------
//		tmp_zero.clear();
//		std::cout<<"dGnD err total "<< cntr::distance_norm2(-1,tempMcc,tmp_zero)<<std::endl;
//		// ----------end test-------------
//

		// -- get M = Z^-1 dGn
		cntr::convolution_timestep(-1,M,Sigma,dGn,integration::I<double>(kt),beta,h);
		M.smul(-1,-1.0);
		tmp_zero.clear();
		std::cout<<"Sigma*dGn err total "<< cntr::distance_norm2(-1,M,tmp_zero)<<std::endl;

		M.incr_timestep(-1,tempM,1.0);

		// -- get MD = -id/dt' M
		cntr::deriv1_matsubara(MD,M,integration::I<double>(kt_mat),beta);
		
		// -- get dDelta
		tempdDelta.set_timestep(-1,Delta_n);
		tempdDelta.incr_timestep(-1,Sigma,1.0);

		cntr::convolution_timestep(-1,dDelta,M,tempdDelta,integration::I<double>(kt),beta,h);
		dDelta.smul(-1,-1.0);

		dDelta.incr_timestep(-1,MD,1.0);

		tmp_zero.clear();
		hyberr_eq=cntr::distance_norm2(-1,dDelta,tmp_zero);
        cout << "update hybridization: iter_n:  " << iter_n << " err: " << hyberr_eq << endl;

		Delta_n.incr_timestep(-1,dDelta,1.0);


		iter_n++;
		if (iter_n > 10*itermax+1) {
		  cerr << "(HYBRIDIZATION) WARNING: Matsubara not converged  after " << itermax
		  	 << "steps ... abort" << endl;
		  cerr << "skip real-time calculation " << endl;

		  break;
		}
	}

	std::ostringstream nameDelta_n;
 	nameDelta_n << "ppsc_data_Delta_n_" << -1 <<"_" << iter<<".h5";
	Delta_n.write_to_hdf5(nameDelta_n.str().c_str(),"Delta_n");

	if (dmfterr_equil < errmax) {

		if(!imp.has_symmetries()) {
		  imp.symmetry_reduction(-1);	    
		  std::string filename = "sym_eq.txt";
		  imp.write_symmetries(filename);
		}
		
		matsubara_converged = true;
		cout << "converged after:  " << iter << endl;
		
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

    if(nt > 0) imp.clear_symmetries();

    if(read_rt_sym) {
      std::string filename = "sym_rt.txt";
      imp.read_symmetries(filename);
    }

    // ---------------------------------------------------------------------
    // 	START ... same iteration
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      GREEN gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      
      for (iter = 1; iter <= itermax; iter++) {

		ppsc::pp_ints_type pp_ints = get_pp_ints(Delta_up, Delta_up_cc, imp.hamiltonian.hil_);
		ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
		
		imp.update_diagrams(pp_ints, gf_verts);
		imp.pp_step(kt);

		for (int n = 0; n <= kt; n++) {
			ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
			Gloc_up.set_timestep(n, gf_tstps[0]);	

			cntr::herm_matrix_timestep<double> tmp;
			Gloc_up.get_timestep(n, tmp);
			Delta_up.set_timestep(n, tmp);
			ppsc::set_bwd_from_fwd(n, Delta_up_cc, Delta_up);

//			get_glatt(n,Z,Delta_up,G1,G1c,Gloc_up,beta,h,kt);
//			get_Sigma_tstp(iter,n,Sigma,F,Fcc,Z,Q,beta,h,kt,nt);

		}

		dmfterr = cntr::distance_norm2(kt, gtmp, Gloc_up);
		gtmp.set_timestep(kt, Gloc_up);

/*		iter_n=1;
		hyberr=1;
		//while(hyberr>errmax*1000){
			// update hybridization
			for (int n=0; n<=kt;n++){
				Delta_n.set_timestep(n,Delta_up);	
			}
			// -- calculate Gn=Z+Z*Delta*Glatt
			// convolution on warm-up steps:
			for (int n=0; n<=kt;n++){
				cntr::convolution_timestep(n,A,Z,Delta_n,integration::I<double>(kt),beta,h);
				cntr::convolution_timestep(n,Acc,Delta_n,Z,integration::I<double>(kt),beta,h);
			}
			for(int tstp=0;tstp<=kt;tstp++){
				cntr::convolution_timestep(tstp,Gn,A,Acc,Gloc_up,Gloc_up,integration::I<double>(kt),beta,h);
				Gn.incr_timestep(tstp,Z,1.0);
			}


			// -- define dG_n = G_loc - G_n
			for (int n=0; n<=kt;n++){
				dGn.set_timestep(n,Gloc_up);
				dGn.incr_timestep(n,Gn,-1.0);
			}
	
			// -- get DdGn=id/dt dGn
			for (int n=0; n<=kt;n++){
				cntr::deriv1_timestep(n,tempM,dGn,dGn,integration::I<double>(kt),beta,h);
			}
			// -- get dGnD=-id/dt' dGn
			for (int n=0; n<=kt;n++){
				cntr::deriv2_timestep(n,tempMcc,dGn,dGn,integration::I<double>(kt),beta,h);
			}

			// -- get DdGnD
			for (int n=0; n<=kt;n++){
				cntr::deriv1_timestep(n,DdGnD,tempMcc,tempM,integration::I<double>(kt),beta,h);
			}

			// -- get M = Z^-1 dGn
			for (int n=0; n<=kt;n++){
				cntr::convolution_timestep(n,M,Sigma,dGn,integration::I<double>(kt),beta,h);
				M.smul(n,-1.0);

				M.incr_timestep(n,tempM,1.0);
			}

			// -- get Mcc = dGn Z^-1
			for (int n=0; n<=kt;n++){
				cntr::convolution_timestep(n,Mcc,dGn,Sigma,integration::I<double>(kt),beta,h);
				Mcc.smul(n,-1.0);

				Mcc.incr_timestep(n,tempMcc,1.0);
			}

			// -- get MD = -id/dt' M
			for (int n=0; n<=kt;n++){
				cntr::deriv2_timestep(n,MD,M,Mcc,integration::I<double>(kt),beta,h);
			}

			// -- get dDelta
			for (int n=0; n<=kt;n++){
				tempdDelta.set_timestep(n,Delta_n);
				tempdDelta.incr_timestep(n,Sigma,1.0);
				cntr::convolution_timestep(n,dDelta,M,tempdDelta,integration::I<double>(kt),beta,h);
				dDelta.smul(n,-1.0);

				dDelta.incr_timestep(n,MD,1.0);
			}
			hyberr=cntr::distance_norm2(kt,dDelta,tmp_zero);
        	cout << "update hybridization: iter_n:  " << iter_n << " err: " << hyberr << endl;
			iter_n++;*/
//		}

        cout << "START: iter:  " << iter << " err: " << dmfterr << endl;
        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }
      }
	  
	  for (int n=0;n<=kt;n++){
		  get_glatt(n,Z,Delta_up,G1,G1c,Gloc_up,beta,h,kt);
		  //get_Sigma_tstp(n,Sigma,F,Fcc,Z,Q,beta,h,kt,nt);
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
    // 	REALTIME: ONLY FOR NT>0
    for (tstp = kt + 1; tstp <= nt; tstp++) {
      cout << "tstp:  " << tstp << endl;
      imp.extrapolate_timestep(tstp - 1);
      
      for (iter = 1; iter <= iter_rtime; iter++) {

		ppsc::pp_ints_type pp_ints = get_pp_ints(Delta_up, Delta_up_cc, imp.hamiltonian.hil_);
		ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);

		imp.update_diagrams(pp_ints, gf_verts);
		imp.pp_step(tstp);

		ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
		Gloc_up.set_timestep(tstp, gf_tstps[0]);	

		cntr::herm_matrix_timestep<double> tmp;
		Gloc_up.get_timestep(tstp, tmp);
		Delta_up.set_timestep(tstp, tmp);
		ppsc::set_bwd_from_fwd(tstp, Delta_up_cc, Delta_up);	
/*		get_glatt(tstp,Z,Delta_up,G1,G1c,Gloc_up,beta,h,kt);
		get_Sigma_tstp(tstp,Sigma,F,Fcc,Z,Q,beta,h,kt,nt);

		iter_n=1;
		hyberr=1;
//		while(hyberr>errmax*1000){
			// make a guess for the hybridization function
			Delta_n.set_timestep(tstp,Delta_up);
			// -- calculate Gn=Z+Z*Delta_n*Glatt
			// convolution on all other steps:
			cntr::convolution_timestep(tstp,A,Z,Delta_n,integration::I<double>(kt),beta,h);
			cntr::convolution_timestep(tstp,Acc,Delta_n,Z,integration::I<double>(kt),beta,h);
			cntr::convolution_timestep(tstp,Gn,A,Acc,Gloc_up,Gloc_up,integration::I<double>(kt),beta,h);
			Gn.incr_timestep(tstp,Z,1.0);

			// -- define dG_n = G_loc - G_n
			dGn.set_timestep(tstp,Gloc_up);
			dGn.incr_timestep(tstp,Gn,-1.0);

			// -- get DdGn=id/dt dGn
			cntr::deriv1_timestep(tstp,tempM,dGn,dGn,integration::I<double>(kt),beta,h);

			// -- get dGnD=-id/dt' dGn
			cntr::deriv2_timestep(tstp,tempMcc,dGn,dGn,integration::I<double>(kt),beta,h);

			// -- get DdGnD
			cntr::deriv1_timestep(tstp,DdGnD,tempMcc,tempM,integration::I<double>(kt),beta,h);

			// -- get M = Z^-1 dGn
			cntr::convolution_timestep(tstp,M,Sigma,dGn,integration::I<double>(kt),beta,h);
			M.smul(tstp,-1.0);

			M.incr_timestep(tstp,tempM,1.0);

			// -- get Mcc = dGn Z^-1
			cntr::convolution_timestep(tstp,Mcc,dGn,Sigma,integration::I<double>(kt),beta,h);
			Mcc.smul(tstp,-1.0);

			Mcc.incr_timestep(tstp,tempMcc,1.0);

			// -- get MD = -id/dt' M
			cntr::deriv2_timestep(tstp,MD,M,Mcc,integration::I<double>(kt),beta,h);

			//get dDelta
			tempdDelta.set_timestep(tstp,Delta_n);
			tempdDelta.incr_timestep(tstp,Sigma,1.0);
			cntr::convolution_timestep(tstp,dDelta,M,tempdDelta,integration::I<double>(kt),beta,h);
			dDelta.smul(tstp,-1.0);

			dDelta.incr_timestep(tstp,MD,1.0);

			hyberr=cntr::distance_norm2(tstp,dDelta,tmp_zero);
        	cout << "update hybridization: iter_n:  " << iter_n << " err: " << hyberr << endl;
			iter_n++;*/
//	  	}
      }


		get_glatt(tstp,Z,Delta_up,G1,G1c,Gloc_up,beta,h,kt);
		//get_Sigma_tstp(tstp,Sigma,F,Fcc,Z,Q,beta,h,kt,nt);

    }

	
//	get Sigma (approxiomatic method)
	if (hartree){
		get_Sigma_wHartree(Sigma,F,Fcc,Z,Q,beta,h,kt,nt);
	}else get_Sigma(Sigma,F,Fcc,Z,Q,beta,h,kt,nt);

    // ---------------------------------------------------------------------
    //test0: calculate hybridization function

//	tstp=-1;
//	cntr::convolution_timestep(tstp,A,Z,Delta_up,integration::I<double>(kt),beta,h);
//	cntr::convolution_timestep(tstp,Acc,Delta_up,Z,integration::I<double>(kt),beta,h);
//	cntr::convolution_timestep(tstp,Gn,A,Acc,Gloc_up,Gloc_up,integration::I<double>(kt),beta,h);
//	Gn.incr_timestep(tstp,Z,1.0);


	// perform calculation on Matsubara axis
//	std::cout<<"equilibrium "<<"\n"<<std::endl;

	// define dG_n = G_loc - G_n
//	dGn.set_timestep(-1,Gloc_up);
//	dGn.incr_timestep(-1,Gn,-1.0);
  //
//	std::cout<<"dGn err total "<< cntr::distance_norm2(-1,dGn,tmp_zero)<<std::endl;
//	std::cout<<"ret "<< cntr::distance_norm2_ret(-1,dGn,tmp_zero)<<std::endl;
//	std::cout<<"les "<< cntr::distance_norm2_les(-1,dGn,tmp_zero)<<std::endl;
//	std::cout<<"tv "<< cntr::distance_norm2_tv(-1,dGn,tmp_zero)<<"\n"<<std::endl;

	// define dG_ncc = G_loc - G_ncc
/*  dGncc.set_timestep(-1,Gloc_up);
  dGncc.incr_timestep(-1,Gncc,-1.0);*/

  // get DdGn=id/dt dGn
//	cntr::deriv1_matsubara(tempM,dGn,integration::I<double>(kt),beta);
//
//	std::cout<<"DdGn err total "<< cntr::distance_norm2(-1,tempM,tmp_zero)<<std::endl;
//	std::cout<<"ret "<< cntr::distance_norm2_ret(-1,tempM,tmp_zero)<<std::endl;
//	std::cout<<"les "<< cntr::distance_norm2_les(-1,tempM,tmp_zero)<<std::endl;
//	std::cout<<"tv "<< cntr::distance_norm2_tv(-1,tempM,tmp_zero)<<"\n"<<std::endl;

//	cntr::deriv1_matsubara(tempMcc,dGn,integration::I<double>(kt),beta);
//
//	std::cout<<"dGnD err total "<< cntr::distance_norm2(-1,tempMcc,tmp_zero)<<std::endl;
//	std::cout<<"ret "<< cntr::distance_norm2_ret(-1,tempMcc,tmp_zero)<<std::endl;
//	std::cout<<"les "<< cntr::distance_norm2_les(-1,tempMcc,tmp_zero)<<std::endl;
//	std::cout<<"tv "<< cntr::distance_norm2_tv(-1,tempMcc,tmp_zero)<<"\n"<<std::endl;
//
//  // get DdGnD
//	cntr::deriv1_matsubara(DdGnD,tempM,integration::I<double>(kt),beta);
//
//	std::cout<<"DdGnD err total "<< cntr::distance_norm2(-1,DdGnD,tmp_zero)<<std::endl;
//	std::cout<<"ret "<< cntr::distance_norm2_ret(-1,DdGnD,tmp_zero)<<std::endl;
//	std::cout<<"les "<< cntr::distance_norm2_les(-1,DdGnD,tmp_zero)<<std::endl;
//	std::cout<<"tv "<< cntr::distance_norm2_tv(-1,DdGnD,tmp_zero)<<"\n"<<std::endl;
//
//  // get M = Z^-1 dGn
//	cntr::convolution_timestep(-1,M,Sigma,dGn,integration::I<double>(kt),beta,h);
//	M.smul(-1,-1.0);

//	M.incr_timestep(-1,tempM,1.0);

//	std::cout<<"Sigma*dGn err total "<< cntr::distance_norm2(-1,M,tmp_zero)<<std::endl;
//	std::cout<<"ret "<< cntr::distance_norm2_ret(-1,M,tmp_zero)<<std::endl;
//	std::cout<<"les "<< cntr::distance_norm2_les(-1,M,tmp_zero)<<std::endl;
//	std::cout<<"tv "<< cntr::distance_norm2_tv(-1,M,tmp_zero)<<"\n"<<std::endl;

//	// get Mcc = dGn Z^-1
//	cntr::convolution_timestep(-1,Mcc,dGn,Sigma,integration::I<double>(kt),beta,h);
//	Mcc.smul(-1,-1.0);
//
//	Mcc.incr_timestep(-1,tempMcc,1.0);
//  
//	std::cout<<"dGn*Sigma err total "<< cntr::distance_norm2(-1,Mcc,tmp_zero)<<std::endl;
//	std::cout<<"ret "<< cntr::distance_norm2_ret(-1,Mcc,tmp_zero)<<std::endl;
//	std::cout<<"les "<< cntr::distance_norm2_les(-1,Mcc,tmp_zero)<<std::endl;
//	std::cout<<"tv "<< cntr::distance_norm2_tv(-1,Mcc,tmp_zero)<<"\n"<<std::endl;

//	// get MD
//	cntr::deriv1_matsubara(MD,M,integration::I<double>(kt),beta);
//
//	//get dDelta
//	tempdDelta.set_timestep(-1,Delta_up);
//	tempdDelta.incr_timestep(-1,Sigma,1.0);
//	cntr::convolution_timestep(-1,dDelta,M,tempdDelta,integration::I<double>(kt),beta,h);
//	dDelta.smul(-1,-1.0);
//
//	dDelta.incr_timestep(-1,MD,1.0);
//  
//	tmp_zero.clear();
//	hyberr=cntr::distance_norm2(-1,dDelta,tmp_zero);
//    cout << "update hybridization (post process): err: " << hyberr << endl;
//	std::cout<<"ret "<< cntr::distance_norm2_ret(-1,dDelta,tmp_zero)<<std::endl;
//	std::cout<<"les "<< cntr::distance_norm2_les(-1,dDelta,tmp_zero)<<std::endl;
//	std::cout<<"tv "<< cntr::distance_norm2_tv(-1,dDelta,tmp_zero)<<"\n"<<std::endl;

	// approximate hybridization function
	for(int tstp=0;tstp<=nt;tstp++){
		Delta_n.set_timestep(tstp,Delta_up);
	}
	// convolution on warm-up steps:
	for(int tstp=0;tstp<=kt;tstp++){
		cntr::convolution_timestep(tstp,A,Z,Delta_up,integration::I<double>(kt),beta,h);
		cntr::convolution_timestep(tstp,Acc,Delta_up,Z,integration::I<double>(kt),beta,h);
    }

	for(int tstp=0;tstp<=kt;tstp++){
		cntr::convolution_timestep(tstp,Gn,A,Acc,Gloc_up,Gloc_up,integration::I<double>(kt),beta,h);
		Gn.incr_timestep(tstp,Z,1.0);
	}

	// convolution on all other steps:
	for(int tstp=kt+1;tstp<=nt;tstp++){
		cntr::convolution_timestep(tstp,A,Z,Delta_up,integration::I<double>(kt),beta,h);
		cntr::convolution_timestep(tstp,Acc,Delta_up,Z,integration::I<double>(kt),beta,h);
		cntr::convolution_timestep(tstp,Gn,A,Acc,Gloc_up,Gloc_up,integration::I<double>(kt),beta,h);
		Gn.incr_timestep(tstp,Z,1.0);
	}
	for (int n=0;n<=nt;n++){

		std::cout<<"tstp= "<<n<<std::endl;

		// define dG_n = G_loc - G_n
		dGn.set_timestep(n,Gloc_up);
		dGn.incr_timestep(n,Gn,-1.0);
		//
		std::cout<<"dGn err total "<< cntr::distance_norm2(n,dGn,tmp_zero)<<std::endl;
		std::cout<<"dGn err ret "<< cntr::distance_norm2_ret(n,dGn,tmp_zero)<<std::endl;
		std::cout<<"dGn err les "<< cntr::distance_norm2_les(n,dGn,tmp_zero)<<std::endl;
		std::cout<<"dGn err tv "<< cntr::distance_norm2_tv(n,dGn,tmp_zero)<<"\n"<<std::endl;

/*		// define dG_ncc = G_loc - G_ncc
		dGncc.set_timestep(n,Gloc_up);
		dGncc.incr_timestep(n,Gncc,-1.0);
*/
		// get DdGn=id/dt dGn
		//cntr::deriv1_timestep(n,tempM,dGn,dGncc,integration::I<double>(kt),beta,h);
		cntr::deriv1_timestep(n,tempM,dGn,dGn,integration::I<double>(kt),beta,h);

		std::cout<<"DdGn err total "<< cntr::distance_norm2(n,tempM,tmp_zero)<<std::endl;
		std::cout<<"DdGn err ret "<< cntr::distance_norm2_ret(n,tempM,tmp_zero)<<std::endl;
		std::cout<<"DdGn err les "<< cntr::distance_norm2_les(n,tempM,tmp_zero)<<std::endl;
		std::cout<<"DdGn err tv "<< cntr::distance_norm2_tv(n,tempM,tmp_zero)<<"\n"<<std::endl;

		cntr::deriv2_timestep(n,tempMcc,dGn,dGn,integration::I<double>(kt),beta,h);

		std::cout<<"dGnD err total "<< cntr::distance_norm2(n,tempMcc,tmp_zero)<<std::endl;
		std::cout<<"dGnD err ret "<< cntr::distance_norm2_ret(n,tempMcc,tmp_zero)<<std::endl;
		std::cout<<"dGnD err les "<< cntr::distance_norm2_les(n,tempMcc,tmp_zero)<<std::endl;
		std::cout<<"dGnD err tv "<< cntr::distance_norm2_tv(n,tempMcc,tmp_zero)<<"\n"<<std::endl;

		// get DdGnD
		cntr::deriv1_timestep(n,DdGnD,tempMcc,tempM,integration::I<double>(kt),beta,h);
		cntr::deriv2_timestep(n,DdGnD,tempM,tempM,integration::I<double>(kt),beta,h);

		std::cout<<"DdGnD err total "<< cntr::distance_norm2(n,DdGnD,tmp_zero)<<std::endl;
		std::cout<<"DdGnD err ret "<< cntr::distance_norm2_ret(n,DdGnD,tmp_zero)<<std::endl;
		std::cout<<"DdGnD err les "<< cntr::distance_norm2_les(n,DdGnD,tmp_zero)<<std::endl;
		std::cout<<"DdGnD err tv "<< cntr::distance_norm2_tv(n,DdGnD,tmp_zero)<<"\n"<<std::endl;

		// get M = Z^-1 dGn
		cntr::convolution_timestep(n,M,Sigma,dGn,integration::I<double>(kt),beta,h);
		M.smul(n,-1.0);

		M.incr_timestep(n,tempM,1.0);

		std::cout<<"Sigma*dGn err total "<< cntr::distance_norm2(n,M,tmp_zero)<<std::endl;
		std::cout<<"Sigma*dGn err ret "<< cntr::distance_norm2_ret(n,M,tmp_zero)<<std::endl;
		std::cout<<"Sigma*dGn err les "<< cntr::distance_norm2_les(n,M,tmp_zero)<<std::endl;
		std::cout<<"Sigma*dGn err tv "<< cntr::distance_norm2_tv(n,M,tmp_zero)<<"\n"<<std::endl;

		// get Mcc = dGn Z^-1
		cntr::convolution_timestep(n,Mcc,dGn,Sigma,integration::I<double>(kt),beta,h);
		Mcc.smul(n,-1.0);

		Mcc.incr_timestep(n,tempMcc,1.0);
		
		std::cout<<"dGn*Sigma err total "<< cntr::distance_norm2(n,Mcc,tmp_zero)<<std::endl;
		std::cout<<"dGn*Sigma err ret "<< cntr::distance_norm2_ret(n,Mcc,tmp_zero)<<std::endl;
		std::cout<<"dGn*Sigma err les "<< cntr::distance_norm2_les(n,Mcc,tmp_zero)<<std::endl;
		std::cout<<"dGn*Sigma err tv "<< cntr::distance_norm2_tv(n,Mcc,tmp_zero)<<"\n"<<std::endl;

		// get MD
		cntr::deriv2_timestep(n,MD,M,Mcc,integration::I<double>(kt),beta,h);
		cntr::deriv2_timestep(n,MD,M,M,integration::I<double>(kt),beta,h);

		//get dDelta
		tempdDelta.set_timestep(n,Delta_up);
		tempdDelta.incr_timestep(n,Sigma,1.0);
		cntr::convolution_timestep(n,dDelta,M,tempdDelta,integration::I<double>(kt),beta,h);
		dDelta.smul(n,-1.0);

		dDelta.incr_timestep(n,MD,1.0);
		
		tmp_zero.clear();
		hyberr=cntr::distance_norm2(n,dDelta,tmp_zero);
        cout << "update hybridization (post process): err: " << hyberr << endl;
		std::cout<<"test3 ret "<< cntr::distance_norm2_ret(n,dDelta,tmp_zero)<<std::endl;
		std::cout<<"test3 les "<< cntr::distance_norm2_les(n,dDelta,tmp_zero)<<std::endl;
		std::cout<<"test3 tv "<< cntr::distance_norm2_tv(n,dDelta,tmp_zero)<<std::endl;

}

    // ---------------------------------------------------------------------
    //test1: set Sigma into the Dyson equation (id/dt - Sigma)Z=1, in order to compare with the initial Z.
    GREEN newZ(nt,ntau,1,-1);
    double mu=0.0;
    cntr::function<double> tempF(nt);
    tempF.set_zero();
    cntr::dyson(newZ,mu,tempF,Sigma,integration::I<double>(kt), beta,h);
    // ---------------------------------------------------------------------
    //test2: set Sigma into the Dyson equation (id/dt - Sigma-Delta)G=1, in order to compare with the initial G.
    GREEN newG(nt,ntau,1,-1);
    GREEN tempSD(nt,ntau,1,-1);
    mu=0.0;
    tempF.set_zero();
	for (int n=-1;n<=nt;n++){
		tempSD.set_timestep(n,Delta_up);
		tempSD.incr_timestep(n,Sigma,1.0);
		}
    cntr::dyson(newG,mu,tempF,tempSD,integration::I<double>(kt), beta,h);
    // ---------------------------------------------------------------------
    // Kinetic energy

    //imp.calculate_pp_interaction_energy();
      
    std::vector<double> Ekin_up = ppsc::get_kinetic_energy(Gloc_up, Delta_up, beta, h, kt);

    std::vector<double> Ekin;
    Ekin.resize(Ekin_up.size());
    
    for( auto t : ppsc::range(0, Ekin.size()) ) Ekin[t] = 2*Ekin_up[t]; // two spin species

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

      group_id = create_group(file_id, "d");
      store_herm_greens_function(group_id, Delta_up);
      close_group(group_id);

      group_id = create_group(file_id, "dn");
      store_herm_greens_function(group_id, Delta_n);
      close_group(group_id);

      group_id = create_group(file_id, "dcc");
      store_herm_greens_function(group_id, Delta_up_cc);
      close_group(group_id);

      group_id = create_group(file_id, "Z");
      store_herm_greens_function(group_id, Z);
      close_group(group_id);

      group_id = create_group(file_id, "newZ");
      store_herm_greens_function(group_id, newZ);
      close_group(group_id);

      group_id = create_group(file_id, "newG");
      store_herm_greens_function(group_id, newG);
      close_group(group_id);

      group_id = create_group(file_id, "Gn");
      store_herm_greens_function(group_id, Gn);
      close_group(group_id);

      group_id = create_group(file_id, "dGn");
      store_herm_greens_function(group_id, dGn);
      close_group(group_id);

/*   
 *    group_id = create_group(file_id, "dGncc");
      store_herm_greens_function(group_id, dGncc);
      close_group(group_id);
*/
      group_id = create_group(file_id, "tempM");
      store_herm_greens_function(group_id, tempM);
      close_group(group_id);

      group_id = create_group(file_id, "tempMcc");
      store_herm_greens_function(group_id, tempMcc);
      close_group(group_id);

      group_id = create_group(file_id, "M");
      store_herm_greens_function(group_id, M);
      close_group(group_id);

      group_id = create_group(file_id, "Mcc");
      store_herm_greens_function(group_id, Mcc);
      close_group(group_id);

      group_id = create_group(file_id, "MD");
      store_herm_greens_function(group_id, MD);
      close_group(group_id);

      group_id = create_group(file_id, "dDelta");
      store_herm_greens_function(group_id, dDelta);
      close_group(group_id);

      group_id = create_group(file_id, "DdGnD");
      store_herm_greens_function(group_id, DdGnD);
      close_group(group_id);

      group_id = create_group(file_id, "MSDelta");
      store_herm_greens_function(group_id, tempdDelta);
      close_group(group_id);

      group_id = create_group(file_id, "Sigma");
      store_herm_greens_function(group_id, Sigma);
      close_group(group_id);

      group_id = create_group(file_id, "bethe");
      store_real_data_to_hid(group_id, "Ekin", Ekin.data(), Ekin.size());
      store_real_data_to_hid(group_id, "Ekin_up", Ekin_up.data(), Ekin_up.size());
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
