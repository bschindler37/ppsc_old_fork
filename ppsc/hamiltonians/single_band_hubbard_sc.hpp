
#ifndef _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_SC_HPP
#define _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_SC_HPP

// -----------------------------------------------------------------------
//
// Single band Hubbard Hamiltonian builder with pairing fields
//
// Author: P. Werner, philipp.werner@gmail.com (2016)
//         H. U.R. Strand, hugo.strand@gmail.com (2017)
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
template<class HILB> class single_band_hubbard_sc {

public:

  // ---------------------------------------------------------------------
  single_band_hubbard_sc(int nt, HILB & hilbert_space) :
    nt(nt), pair_seed(0.0),
    // -- parameter vectors
    mu(nt+2, 0.0), U(nt+2, 0.0), epsilon(nt+2, 0.0), Bz(nt+2, 0.0),
    pair_real(nt+2, 0.0), pair_imag(nt+2, 0.0),
    // -- expectation value vectors
    Q_exp(nt+2), n_exp(nt+2), m_exp(nt+2), nu_exp(nt+2), nd_exp(nt+2),
    docc_exp(nt+2), pair_exp(nt+2), Eint_exp(nt+2),
    hil_(hilbert_space)  {

    // -- Construct basic operators

    cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
    cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
    cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
    cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

    nu = cuc * cua;
    nd = cdc * cda;
    n = nu + nd;

    m = nu - nd;

    docc = nu * nd;

    cuc_cdc = cuc * cdc;
    cda_cua = cda * cua;
    cua_cda = cua * cda;

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double mu_val = mu[tstp + 1];
    double u = U[tstp + 1];
    double bz = Bz[tstp + 1];
    double eps = epsilon[tstp + 1];
    std::complex<double> pair_field =
      std::complex<double>(pair_real[tstp + 1], pair_imag[tstp + 1]);
    double half_filling_shift = 0.5*u;

    Htemp = u * docc + bz * m + eps * n - (mu_val + half_filling_shift) * n
      +          (pair_field + pair_seed) * cuc_cdc
      + std::conj(pair_field + pair_seed) * cda_cua;
  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1]    = expectation_value(tstp, rho, Q   ).real();

    m_exp[tstp + 1]    = expectation_value(tstp, rho, m   ).real();
    n_exp[tstp + 1]    = expectation_value(tstp, rho, n   ).real();
    nu_exp[tstp + 1]   = expectation_value(tstp, rho, nu  ).real();
    nd_exp[tstp + 1]   = expectation_value(tstp, rho, nd  ).real();
    docc_exp[tstp + 1] = expectation_value(tstp, rho, docc).real();

    pair_exp[tstp + 1] = expectation_value(tstp, rho, cua_cda);

    operator_type Ht; get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();

    std::cout << "tstp " << tstp << " "
	      << "Q,n,m,nu,nd,d,pair: "
	      << Q_exp[tstp + 1] << " "
	      << n_exp[tstp + 1] << " "
	      << m_exp[tstp + 1] << " "
	      << nu_exp[tstp + 1] << " "
	      << nd_exp[tstp + 1] << " "
	      << docc_exp[tstp + 1] << " "
	      << pair_exp[tstp + 1] << " "
	      << std::endl;
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    store_real_data_to_hid(group_id, "mu", mu.data(), mu.size());
    store_real_data_to_hid(group_id, "U", U.data(), U.size());

    store_real_data_to_hid(group_id, "Bz", Bz.data(), Bz.size());
    store_real_data_to_hid(group_id, "eps", epsilon.data(), epsilon.size());

    store_real_data_to_hid(group_id, "pair_real", pair_real.data(), pair_real.size());
    store_real_data_to_hid(group_id, "pair_imag", pair_imag.data(), pair_imag.size());

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "n_exp", n_exp.data(), n_exp.size());
    store_real_data_to_hid(group_id, "m_exp", m_exp.data(), m_exp.size());
    store_real_data_to_hid(group_id, "nu_exp", nu_exp.data(), nu_exp.size());
    store_real_data_to_hid(group_id, "nd_exp", nd_exp.data(), nd_exp.size());
    store_real_data_to_hid(group_id, "docc_exp", docc_exp.data(), docc_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());

    store_cplx_data_to_hid(group_id, "pair_exp", pair_exp.data(), pair_exp.size());
  }

  int nt;
  HILB hil_;

  std::complex<double> pair_seed;
  std::vector<double> mu, U, epsilon, Bz, pair_real, pair_imag;

  std::vector<double> Q_exp, n_exp, nu_exp, nd_exp, m_exp, docc_exp, Eint_exp;
  std::vector<std::complex<double> > pair_exp;

  operator_type cua, cuc, cda, cdc, nu, nd, n, m, docc, Q;
  operator_type cuc_cdc, cda_cua, cua_cda;
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_SC_HPP
