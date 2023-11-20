
#ifndef _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_SPIN_HPP
#define _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_SPIN_HPP

// -----------------------------------------------------------------------
//
// Single band Hubbard Hamiltonian builder with magnetic fields x,y,z
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// H = U * n_{up} * n_{do} + (epsup-U/2-mu) * n_{up} + (epsdo-U/2-mu) * n_{do}
//     + Bx * Sx + By * Sy + Bz * Sz
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
template<class HILB> class single_band_hubbard_spin {

public:

  // ---------------------------------------------------------------------
  single_band_hubbard_spin(int nt, HILB & hilbert_space) :
    nt(nt), Q_exp(nt+2),
    nu_exp(nt+2), nd_exp(nt+2),
    docc_exp(nt+2), Eint_exp(nt+2),
    Sx_exp(nt+2), Sy_exp(nt+2), Sz_exp(nt+2),
    hil_(hilbert_space) {

    // -- Construct basic operators

    cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
    cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
    cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
    cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

    nu = cuc * cua;
    nd = cdc * cda;
    n = nu + nd;

    docc = nu * nd;

    std::complex<double> I(0., 1.);

    operator_type Sp = cuc * cda;
    operator_type Sm = cdc * cua;

    Sx =    0.5*(Sp + Sm);
    Sy = -I*0.5*(Sp - Sm);
    Sz =    0.5*(nu - nd);

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    std::complex<double> epsup1 = eps_up[tstp+1] - mu[tstp+1] - U[tstp+1] * 0.5;
    std::complex<double> epsdo1 = eps_do[tstp+1] - mu[tstp+1] - U[tstp+1] * 0.5;

    Htemp = U[tstp+1] * docc + epsup1 * nu + epsdo1 * nd
      + Bx[tstp+1]*Sx + By[tstp+1]*Sy + Bz[tstp+1]*Sz;

    if(tstp == -1) std::cout << "Hamiltonian: at tstp " << tstp << "\n" << Htemp << "\n";
  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1]    = expectation_value(tstp, rho, Q   ).real();
    nu_exp[tstp + 1]   = expectation_value(tstp, rho, nu  ).real();
    nd_exp[tstp + 1]   = expectation_value(tstp, rho, nd  ).real();
    docc_exp[tstp + 1] = expectation_value(tstp, rho, docc).real();

    Sx_exp[tstp + 1] = expectation_value(tstp, rho, Sx).real();
    Sy_exp[tstp + 1] = expectation_value(tstp, rho, Sy).real();
    Sz_exp[tstp + 1] = expectation_value(tstp, rho, Sz).real();

    operator_type Ht;
    get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();
    Eint_exp[tstp + 1] += mu[tstp+1] * (nu_exp[tstp+1] + nd_exp[tstp+1]);

    // -------------------------------------------------------------------
    // Since the ppdyson uses the transposed hamiltonian the Sy component
    // has a sign error.. here we "patch" this by signs...

    //Sy_exp[tstp + 1] *= -1.0; // Hermicity bug in PPDyson
    //Eint_exp[tstp + 1] += 2.0 * Sy_exp[tstp + 1] * By[tstp + 1]; // Hermicity bug in PPDyson
    // -------------------------------------------------------------------

    double Stot = std::sqrt(
      Sx_exp[tstp+1] * Sx_exp[tstp+1] +
      Sy_exp[tstp+1] * Sy_exp[tstp+1] +
      Sz_exp[tstp+1] * Sz_exp[tstp+1] );

    std::cout << "Q, n, d, St,x,y,z = "
	      << Q_exp[tstp+1] << ", "
	      << nu_exp[tstp+1] + nd_exp[tstp+1] << ", "
	      << docc_exp[tstp+1] << ", "
	      << Stot << ", "
	      << Sx_exp[tstp+1] << ", "
	      << Sy_exp[tstp+1] << ", "
	      << Sz_exp[tstp+1] << std::endl;
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

    //store_double_attribute_to_hid(group_id, "mu", mu);
    store_real_data_to_hid(group_id, "mu", mu.data(), mu.size());
    store_real_data_to_hid(group_id, "U", U.data(), U.size());
    store_real_data_to_hid(group_id, "Bx", Bx.data(), Bx.size());
    store_real_data_to_hid(group_id, "By", By.data(), By.size());
    store_real_data_to_hid(group_id, "Bz", Bz.data(), Bz.size());
    store_real_data_to_hid(group_id, "eps_up", eps_up.data(), eps_up.size());
    store_real_data_to_hid(group_id, "eps_do", eps_do.data(), eps_do.size());

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "nu_exp", nu_exp.data(), nu_exp.size());
    store_real_data_to_hid(group_id, "nd_exp", nd_exp.data(), nd_exp.size());
    store_real_data_to_hid(group_id, "docc_exp", docc_exp.data(), docc_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());

    store_real_data_to_hid(group_id, "Sx_exp", Sx_exp.data(), Sx_exp.size());
    store_real_data_to_hid(group_id, "Sy_exp", Sy_exp.data(), Sy_exp.size());
    store_real_data_to_hid(group_id, "Sz_exp", Sz_exp.data(), Sz_exp.size());
  }

  int nt;
  HILB hil_;

  std::vector<double> mu;
  std::vector<double> U, eps_up, eps_do, Bx, By, Bz;
  std::vector<double> Q_exp, nu_exp, nd_exp, docc_exp, Eint_exp;
  std::vector<double> Sx_exp, Sy_exp, Sz_exp;

  operator_type cua, cuc, cda, cdc, nu, nd, n, docc, Q, Sx, Sy, Sz;
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_SINGLE_BAND_HUBBARD_SPIN_HPP
