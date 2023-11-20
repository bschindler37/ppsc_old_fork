
#ifndef _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_HPP
#define _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_HPP

// -----------------------------------------------------------------------
//
// Two band Hubbard Hamiltonian builder
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
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
template<class HILB> class two_band_hubbard {

public:

  // ---------------------------------------------------------------------
  two_band_hubbard(int nt, HILB & hilbert_space) :
    nt(nt),
    Q_exp(nt+2), Eint_exp(nt+2),
    n1_exp(nt+2), n2_exp(nt+2), m_exp(nt+2),
    densdens_intraband_exp(nt+2),
    densdens_interband_exp(nt+2),
    densdens_interband_equalspin_exp(nt+2), hil_(hilbert_space) {

    // -- Construct basic operators

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

    m = n1u - n1d + n2u - n2d; // magnetization

    densdens_intraband = n1u * n1d + n2u * n2d;
    densdens_interband = n1 * n2;
    densdens_interband_equalspin = n1u * n2u + n1d * n2d;

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double d = delta[tstp + 1];
    double u = U[tstp + 1];
    double j = J[tstp + 1];
    double bz = Bz[tstp + 1];

    double half_filling_shift = 0.5*(3*u - 5*j);

    Htemp.InitDiagZero(hil_.ssdim_);

    Htemp = Htemp
      + bz * m
      + 0.5 * d * (n1 - n2)
      - mu * n
      - half_filling_shift * n
      + (u      ) * densdens_intraband
      + (u - 2*j) * densdens_interband
      + (     -j) * densdens_interband_equalspin;

    /*
    std::cout << "tstp = " << tstp << std::endl;
    std::cout << "mu, d, u, j = " << mu << ", " << d << ", " << u << ", " << j << std::endl;
    std::cout << "H[tstp] = " << Htemp << std::endl;
    */

  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1]    = expectation_value(tstp, rho, Q   ).real();
    operator_type Ht; get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();

    m_exp[tstp + 1]    = expectation_value(tstp, rho, m   ).real();
    n1_exp[tstp + 1]   = expectation_value(tstp, rho, n1  ).real();
    n2_exp[tstp + 1]   = expectation_value(tstp, rho, n2  ).real();

    densdens_intraband_exp[tstp + 1] = expectation_value(tstp, rho, densdens_intraband).real();
    densdens_interband_exp[tstp + 1] = expectation_value(tstp, rho, densdens_interband).real();
    densdens_interband_equalspin_exp[tstp + 1] = expectation_value(tstp, rho, densdens_interband_equalspin).real();

    std::cout << "tstp, n, m, Q = "
	      << tstp << ", "
	      << n1_exp[tstp + 1]+n2_exp[tstp + 1] << ", "
	      << m_exp[tstp + 1] << ", "
	      << Q_exp[tstp + 1] << std::endl;
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    store_double_attribute_to_hid(group_id, "mu", mu);
    store_real_data_to_hid(group_id, "U", U.data(), U.size());
    store_real_data_to_hid(group_id, "J", J.data(), J.size());
    store_real_data_to_hid(group_id, "delta", delta.data(), delta.size());

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());

    store_real_data_to_hid(group_id, "m_exp", m_exp.data(), m_exp.size());
    store_real_data_to_hid(group_id, "n1_exp", n1_exp.data(), n1_exp.size());
    store_real_data_to_hid(group_id, "n2_exp", n2_exp.data(), n2_exp.size());

    store_real_data_to_hid(group_id, "densdens_intraband_exp",
			   densdens_intraband_exp.data(), densdens_intraband_exp.size());
    store_real_data_to_hid(group_id, "densdens_interband_exp",
			   densdens_interband_exp.data(), densdens_interband_exp.size());
    store_real_data_to_hid(group_id, "densdens_interband_equalspin_exp",
			   densdens_interband_equalspin_exp.data(), densdens_interband_equalspin_exp.size());

    // -- Sub group operators:

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
    store_operator_type(sub_group_id, m, std::string("m"));

    store_operator_type(sub_group_id, densdens_intraband, std::string("densdens_intraband"));
    store_operator_type(sub_group_id, densdens_interband, std::string("densdens_interband"));
    store_operator_type(sub_group_id, densdens_interband_equalspin, std::string("densdens_interband_equalspin"));

    close_group(sub_group_id);

  }

  int nt;
  HILB hil_;

  double mu;
  std::vector<double> U, J, delta, Bz;
  std::vector<double> Q_exp, Eint_exp;
  std::vector<double> n1_exp, n2_exp, m_exp;
  std::vector<double> densdens_intraband_exp, densdens_interband_exp, densdens_interband_equalspin_exp;

  operator_type c1ua, c1uc, c1da, c1dc, c2ua, c2uc, c2da, c2dc;
  operator_type n1u, n1d, n2u, n2d, n1, n2, n, Q, m;
  operator_type densdens_intraband, densdens_interband, densdens_interband_equalspin;
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_HPP
