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

// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------

 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  int itermax=100, iter_rtime=10, nt=80, ntau=300, kt=5, iter, tstp, nomp=1;

  double mu = 0.5, beta = 10., h = 0.05, errmax=1e-10, dmfterr=1e-10, dmfterr_equil=1e-10;
  bool matsubara_converged = false;

  cntr::function<double> h_func(nt);
  for(int tstp=-1;tstp<=nt;tstp++){
    cdmatrix tmp(1,1);
    tmp(0,0) = 0.;
    h_func.set_value(tstp, tmp);
  }
  
  // ---------------------------------------------------------------------
    // -- setup single particle greens functions and hybridizations

    cntr::herm_matrix<double> G(nt, ntau, 1, -1);    
    cntr::herm_matrix<double> Delta(nt, ntau, 1, -1);

    // ---------------------------------------------------------------------
    // MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
    {
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);

      gtmp.set_timestep(-1, G); // this avoids one iteration

      for (iter = 1; iter <= itermax; iter++) {
		
	cntr::dyson_mat(G, mu, h_func, Delta, beta, kt);
	
	// -- Check error
	dmfterr_equil = cntr::distance_norm2(-1, gtmp, G);
        gtmp.set_timestep(-1, G);

	// -- Update Hybridization
	cntr::herm_matrix_timestep<double> tmp;
	G.get_timestep(-1, tmp);
	Delta.set_timestep(-1, tmp);
	
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
      
      for (iter = 1; iter <= itermax; iter++) {

	cntr::dyson_start(G, mu, h_func, Delta, beta, h, kt);

	for (int n = 0; n <= kt; n++) {
	  cntr::herm_matrix_timestep<double> tmp;
	  G.get_timestep(n, tmp);
	  Delta.set_timestep(n, tmp);
	}
	
        dmfterr = cntr::distance_norm2(kt, gtmp, G);
        gtmp.set_timestep(kt, G);
	
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
      
      for (iter = 1; iter <= iter_rtime; iter++) {

        cntr::dyson_timestep(tstp, G, mu, h_func, Delta, beta, h, kt);

	cntr::herm_matrix_timestep<double> tmp;
	G.get_timestep(tstp, tmp);
	Delta.set_timestep(tstp, tmp);
      }
    }

    // ---------------------------------------------------------------------
    // Kinetic energy

    //std::vector<double> Ekin = ppsc::get_kinetic_energy(G, Delta, beta, h, kt);

    // ---------------------------------------------------------------------
    // OBSERVABLES
    {
      std:string filename = "data_bethe.h5";
      hid_t file_id = open_hdf5_file(filename);
      hid_t group_id;
      
      group_id = create_group(file_id, "g");
      store_herm_greens_function(group_id, G);
      close_group(group_id);

      /*
      group_id = create_group(file_id, "bethe");
      store_real_data_to_hid(group_id, "Ekin", Ekin.data(), Ekin.size());
      store_double_attribute_to_hid(group_id, "dmfterr_equil", dmfterr_equil);
      close_group(group_id);      
      */
      
      close_hdf5_file(file_id);
    }
    
  return 0;
}
