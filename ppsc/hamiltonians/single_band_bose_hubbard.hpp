
#ifndef _PPSC_HAMILTONIANS_SINGLE_BAND_BOSE_HUBBARD_HPP
#define _PPSC_HAMILTONIANS_SINGLE_BAND_BOSE_HUBBARD_HPP

// -----------------------------------------------------------------------
//
// Single band Bose-Hubbard Hamiltonian builder
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// H = U * n * (n + 1) - mu * n
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
template<class HILB> class single_band_bose_hubbard {

public:

  // ---------------------------------------------------------------------
  single_band_bose_hubbard(int nt, HILB & hilbert_space) :
    nt(nt), Q_exp(nt+2),
    n_exp(nt+2), docc_exp(nt+2), Eint_exp(nt+2),
    hil_(hilbert_space) {

    // -- Construct basic operators

    ca = hil_.c_op_[hil_.flavor(0, 0, 0)];
    cc = hil_.c_op_[hil_.flavor(0, 0, 1)];
    Q = operator_type::Identity(hil_);

    n = cc * ca;
    docc = n * (n - Q);

  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {
    Htemp = (0.5 * U[tstp+1]) * docc - mu * n;
  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1]    = expectation_value(tstp, rho, Q   ).real();
    n_exp[tstp + 1]    = expectation_value(tstp, rho, n   ).real();
    docc_exp[tstp + 1] = expectation_value(tstp, rho, docc).real();

    operator_type Ht;
    get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();

    std::cout << "tstp, Q, n, docc = "
	      << tstp << ", "
              << Q_exp[tstp + 1] << ", "
	      << n_exp[tstp + 1] << ", "
	      << docc_exp[tstp + 1] << std::endl;
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    store_double_attribute_to_hid(group_id, "mu", mu);
    store_real_data_to_hid(group_id, "U", U.data(), U.size());

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "n_exp", n_exp.data(), n_exp.size());
    store_real_data_to_hid(group_id, "docc_exp", docc_exp.data(), docc_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());
  }

  int nt;
  HILB hil_;

  double mu;
  std::vector<double> U;
  std::vector<double> Q_exp, n_exp, docc_exp, Eint_exp;

  operator_type ca, cc, n, docc, Q;
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_SINGLE_BAND_BOSE_HUBBARD_HPP
