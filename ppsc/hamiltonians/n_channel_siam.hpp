
#ifndef _PPSC_HAMILTONIANS_TWO_CHANNEL_KONDO_HPP
#define _PPSC_HAMILTONIANS_TWO_CHANNEL_KONDO_HPP

/* -----------------------------------------------------------------------
//
// N-Channel Single-impurity Anderson Hamiltonian builder, should be based on a
// n_channel_one_spin hilbert space. This is the local hamiltonian for the following impurity model with n bath channels

H = sum_{s=up,do} sum_m={0}^{nchannel-1} sum_k  psi_{m,s}^dag c_{k,m,s} + h.c.
- epsf*(f_{up}^dagger f_{up} + f_{do}^dagger f_{do})
+ bz * (f_{up}^dagger f_{up} - f_{do}^dagger f_{do})
+ sum_{k,s,m} c_{k,m,s}^dag c_{k,m,s} eps_k

with psi_{m,s} = f_{s}b_m^dag build from a fermion and a pseudo-boson

see Cox and Ruckenstein, PRL 71, 1613 (1993)

one may use that

f_{s}^dagger f_{s} =
= psi_{m,s}^dag psi_{m,s}   (for each m)
= 1/nchannel sum_m psi_{m,s}^dag psi_{m,s},

// Author: Martin Eckstein, martin.eckstein@mpsd.cfel.de (2016)

//---------------------------------------------------------------------*/



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
template<class HILB,int NC> class n_channel_siam {

public:

  // ---------------------------------------------------------------------
  n_channel_siam(int nt) :
    nt(nt),
    Q_exp(nt+2), Eint_exp(nt+2),
    n_exp(nt+2), m_exp(nt+2)
    {
    nchannel=NC;
    hil_.init(nchannel);

    // -- Construct basic operators
    // f_{s}^dagger f_s
    int spin=0;
    nu = hil_.c_op_[hil_.flavor(0,spin, 1)] * hil_.c_op_[hil_.flavor(0,spin, 0)];
    for(int m=1;m<nchannel;m++){
       nu = nu + hil_.c_op_[hil_.flavor(m,spin, 1)] * hil_.c_op_[hil_.flavor(m,spin, 0)];
    }
    nu = (1.0/nchannel)*nu;
    spin=1;
    nd = hil_.c_op_[hil_.flavor(0,spin, 1)] * hil_.c_op_[hil_.flavor(0,spin, 0)];
    for(int m=1;m<nchannel;m++){
       nd = nd + hil_.c_op_[hil_.flavor(m,spin, 1)] * hil_.c_op_[hil_.flavor(m,spin, 0)];
    }
    nd = (1.0/nchannel)*nd;
    n=nu+nd;
    m = nu - nd; // magnetization
    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double epsf = Epsf[tstp + 1];
    double bz = Bz[tstp + 1];

    Htemp.InitDiagZero(hil_.ssdim_);

    Htemp = Htemp + bz * m - epsf * n;

  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1]    = expectation_value(tstp, rho, Q   ).real();
    operator_type Ht; get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();

    m_exp[tstp + 1]    = expectation_value(tstp, rho, m   ).real();
    n_exp[tstp + 1]   = expectation_value(tstp, rho, n  ).real();

    std::cout << "tstp, n, m, Q = "
	      << tstp << ", "
	      << n_exp[tstp + 1] << ", "
	      << m_exp[tstp + 1] << ", "
	      << Q_exp[tstp + 1] << std::endl;
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    store_real_data_to_hid(group_id, "Epsf", Epsf.data(), Epsf.size());
    store_real_data_to_hid(group_id, "Bz", Bz.data(), Bz.size());
    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());
    store_real_data_to_hid(group_id, "m_exp", m_exp.data(), m_exp.size());
    store_real_data_to_hid(group_id, "n_exp", n_exp.data(), n_exp.size());

    // -- Sub group operators:

    hid_t sub_group_id = create_group(group_id, "op");

    store_operator_type(sub_group_id, nu, std::string("nu"));
    store_operator_type(sub_group_id, nd, std::string("nd"));
    store_operator_type(sub_group_id, n, std::string("n"));
    store_operator_type(sub_group_id, Q, std::string("Q"));
    store_operator_type(sub_group_id, m, std::string("m"));

    close_group(sub_group_id);

  }

  int nt,nchannel;
  HILB hil_;
  std::vector<double> Epsf, Bz;
  std::vector<double> Q_exp, Eint_exp;
  std::vector<double> n_exp, m_exp;
  operator_type nu, nd,n, Q, m;
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_HPP
