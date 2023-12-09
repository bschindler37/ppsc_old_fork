
#ifndef _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_XAS2_HPP
#define _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_XAS2_HPP

// -----------------------------------------------------------------------
//
// Two band Hubbard Hamiltonian builder
// with light-matter interaction hopping between the two orbitals for t >= 0
//
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#include "ppsc/ppsc.hpp"
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hamiltonians {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class HILB> class two_band_hubbard_xas2 {

public:

  // ---------------------------------------------------------------------
  two_band_hubbard_xas2(int nt, HILB & hilbert_space) :
    nt(nt),
    Q_exp(nt+2), Eint_exp(nt+2), time(nt+2),
    n1_exp(nt+2), n2_exp(nt+2), P_in_exp(nt+2), P_in_up_exp(nt+2), P_in_down_exp(nt+2), P_out_exp(nt+2),
    hil_(hilbert_space) {

    // -- Construct basic operators 
    // orbital 1 = valence, orbital 2 = core orbital

    c1ua = hil_.c_op_[hil_.flavor(0, 0, 0)];
    c1uc = hil_.c_op_[hil_.flavor(0, 0, 1)];
    c1da = hil_.c_op_[hil_.flavor(0, 1, 0)];
    c1dc = hil_.c_op_[hil_.flavor(0, 1, 1)];

    c2ua = hil_.c_op_[hil_.flavor(1, 0, 0)];
    c2uc = hil_.c_op_[hil_.flavor(1, 0, 1)];
    c2da = hil_.c_op_[hil_.flavor(1, 1, 0)];
    c2dc = hil_.c_op_[hil_.flavor(1, 1, 1)];

    n1u = c1uc * c1ua;
    n1d = c1dc * c1da;
    n1 = n1u + n1d;

    n2u = c2uc * c2ua;
    n2d = c2dc * c2da;
    n2 = n2u + n2d;

    n = n1 + n2;

    // m = n1u - n1d + n2u - n2d; // magnetization

    dens_updown_1 = n1u*n1d;
    
    P_in_up = c1uc * c2ua;	// (spin up) d^dag_up c_up
    P_in_down = c1dc * c2da;	// (spin down) d^dag_down c_down
    P_in = P_in_up + P_in_down;
    P_out = c2uc*c1ua + c2dc*c1da;

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double e_c_temp = e_core[tstp + 1];
    std::complex<double> probe_env_exp = probe_pulse_env[tstp + 1] * omega_in_osc[tstp + 1];
    

    Htemp.InitDiagZero(hil_.ssdim_);

    Htemp = Htemp
    		+ e_c_temp * n2 + e_valence * n1	// core & valence energy level
    	 	+ U * dens_updown_1	// valence Hubbard U interaction
    	 	+ g * ( P_in * probe_env_exp + P_out * std::conj(probe_env_exp) ) ;   // light-induced electron-hole excitation, for tstp >= 0
    
     // - mu * n ???


    /*
    std::cout << "tstp = " << tstp << std::endl;
    std::cout << "mu, d, u, j = " << mu << ", " << d << ", " << u << ", " << j << std::endl;
    std::cout << "H[tstp] = " << Htemp << std::endl;
    */

  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {
	
    time[tstp + 1] = tstp*h;	
	
    Q_exp[tstp + 1]    = expectation_value(tstp, rho, Q ).real();
    operator_type Ht; get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht ).real();

    n1_exp[tstp + 1]   = expectation_value(tstp, rho, n1 ).real();
    n2_exp[tstp + 1]   = expectation_value(tstp, rho, n2 ).real();

    P_in_up_exp[tstp + 1]   = expectation_value(tstp, rho, P_in_up );
    P_in_down_exp[tstp + 1]   = expectation_value(tstp, rho, P_in_down );
    
    P_in_exp[tstp + 1] = P_in_up_exp[tstp + 1] + P_in_down_exp[tstp + 1];
    P_out_exp[tstp + 1] = expectation_value(tstp, rho, P_out );
    
    
    /*
	    std::cout << "tstp, n_d, n_c, Q = "
		      << tstp << ", "
		      << n1_exp[tstp + 1] << ", " << n2_exp[tstp + 1] << ", "
		      << Q_exp[tstp + 1] << std::endl;      
    */
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    // functions like 'store_real_data_to_hid' are defined in hdf scope of cntr (NESSi) library in the file 'nessi/libcntr/cntr/hdf5/hdf_interface.hpp'
    store_real_data_to_hid(group_id, "time", time.data(), time.size());
    
    store_real_data_to_hid(group_id, "e_core", e_core.data(), e_core.size());
    store_cplx_data_to_hid(group_id, "probe pulse envelope s(t)", probe_pulse_env.data(), probe_pulse_env.size());
        
    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());

    store_real_data_to_hid(group_id, "n1_exp", n1_exp.data(), n1_exp.size());
    store_real_data_to_hid(group_id, "n2_exp", n2_exp.data(), n2_exp.size());
    
    store_cplx_data_to_hid(group_id, "P_out_exp", P_out_exp.data(), P_out_exp.size());
    store_cplx_data_to_hid(group_id, "P_in_up_exp", P_in_up_exp.data(), P_in_up_exp.size());
    store_cplx_data_to_hid(group_id, "P_in_down_exp", P_in_down_exp.data(), P_in_down_exp.size());
    store_cplx_data_to_hid(group_id, "P_in_exp", P_in_exp.data(), P_in_exp.size());
        
    /* -- Sub group operators:

    hid_t sub_group_id = create_group(group_id, "op");

    store_operator_type(sub_group_id, c1uc, std::string("c1uc"));
    store_operator_type(sub_group_id, c1ua, std::string("c1ua"));
    store_operator_type(sub_group_id, c1dc, std::string("c1dc"));
    store_operator_type(sub_group_id, c1da, std::string("c1da"));

    store_operator_type(sub_group_id, c2uc, std::string("c2uc"));
    store_operator_type(sub_group_id, c2ua, std::string("c2ua"));
    store_operator_type(sub_group_id, c2dc, std::string("c2dc"));
    store_operator_type(sub_group_id, c2da, std::string("c2da"));

    store_operator_type(sub_group_id, n1u, std::string("n1u"));
    store_operator_type(sub_group_id, n1d, std::string("n1d"));

    store_operator_type(sub_group_id, n2u, std::string("n2u"));
    store_operator_type(sub_group_id, n2d, std::string("n2d"));

    store_operator_type(sub_group_id, n1, std::string("n1"));
    store_operator_type(sub_group_id, n2, std::string("n2"));

    store_operator_type(sub_group_id, n, std::string("n"));
    store_operator_type(sub_group_id, Q, std::string("Q"));
    
    store_operator_type(sub_group_id, P_in_up, std::string("P_in_up"));
    store_operator_type(sub_group_id, P_in_down, std::string("P_in_down"));
    store_operator_type(sub_group_id, P_in, std::string("P_in"));
    
    close_group(sub_group_id);
    */
    
  }

  int nt;	// number of real time steps
  double h;	// time difference/step size
  HILB hil_;
  
  // need to be initialized:
  double e_valence, U;
  double g, omega_in;
  std::vector<double> e_core;
  std::vector< std::complex<double> > probe_pulse_env;	// s(t) for t>= 0; 0 else (no light-matter induced core-valence interaction on Matsubara branch)
  std::vector< std::complex<double> > omega_in_osc;	// exp(-II*omega_in*h*tstp)
  
  // time array
  std::vector<double> time;
  // save expectation values:
  std::vector<double> Q_exp, Eint_exp;
  std::vector<double> n1_exp, n2_exp;
  std::vector<cdouble> P_in_up_exp, P_in_down_exp, P_in_exp, P_out_exp;
  
  operator_type c1ua, c1uc, c1da, c1dc, c2ua, c2uc, c2da, c2dc;
  operator_type dens_updown_1, n1u, n1d, n2u, n2d, n1, n2, n, Q;
  operator_type P_in_up, P_in_down, P_in, P_out;
  
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_XAS2_HPP



