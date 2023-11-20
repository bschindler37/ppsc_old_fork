
#ifndef _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_HPP
#define _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_HPP

// -----------------------------------------------------------------------
//
// Single band Hubbard Hamiltonian builder
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// H = U * n_{up} * n_{do} + (epsup-U/2-mu) * n_{up} + (epsdo-U/2-mu) * n_{do}
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
template<class HILB> class single_band_hubbard {

public:

  // --------------------------------------------------------------------- 
  // -- Constructor (class members are defined below)
  single_band_hubbard(int nt, HILB & hilbert_space) :
    nt(nt), Q_exp(nt+2),
    nu_exp(nt+2), nd_exp(nt+2),
    docc_exp(nt+2), Eint_exp(nt+2), hil_(hilbert_space) {

    // -- Construct basic operators, c_op_ is a list of all c and c_dag, ordered by flavor index
    // -- flavor = ( orb, spin, creation(1)/annihilation(0) )
    cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
    cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
    cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
    cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

    nu = cuc * cua;
    nd = cdc * cda;
    n = nu + nd;

    docc = nu * nd;	// double occupation n_{up down}

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double epsup1 = eps_up[tstp+1] - mu - U[tstp+1] * 0.5;
    double epsdo1 = eps_do[tstp+1] - mu - U[tstp+1] * 0.5;
    Htemp = U[tstp+1] * docc + epsup1 * nu + epsdo1 * nd;
  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1]    = expectation_value(tstp, rho, Q   ).real();
    nu_exp[tstp + 1]   = expectation_value(tstp, rho, nu  ).real();
    nd_exp[tstp + 1]   = expectation_value(tstp, rho, nd  ).real();
    docc_exp[tstp + 1] = expectation_value(tstp, rho, docc).real();

    operator_type Ht;
    get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();

    /*
    std::cout << "tstp, nup, ndo, docc = "
	      << tstp << ", "
	      << nu_exp[tstp + 1] << ", "
	      << nd_exp[tstp + 1] << ", "
	      << docc_exp[tstp + 1] << std::endl;
    */
  }

  // ---------------------------------------------------------------------
  void local_obs(int tstp, double & nup_out, double & ndo_out,
		 double & docc_out) {
    nup_out = nu_exp[tstp + 1];
    ndo_out = nd_exp[tstp + 1];
    docc_out = docc_exp[tstp + 1];
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    store_double_attribute_to_hid(group_id, "mu", mu);
    store_real_data_to_hid(group_id, "U", U.data(), U.size());
    store_real_data_to_hid(group_id, "eps_up", eps_up.data(), eps_up.size());
    store_real_data_to_hid(group_id, "eps_do", eps_do.data(), eps_do.size());

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "nu_exp", nu_exp.data(), nu_exp.size());
    store_real_data_to_hid(group_id, "nd_exp", nd_exp.data(), nd_exp.size());
    store_real_data_to_hid(group_id, "docc_exp", docc_exp.data(), docc_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());
  }

  int nt;
  HILB hil_;

  double mu;
  std::vector<double> U, eps_up, eps_do;
  std::vector<double> Q_exp, nu_exp, nd_exp, docc_exp, Eint_exp;

  operator_type cua, cuc, cda, cdc, nu, nd, n, docc, Q;
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_HPP
