
#ifndef _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_ISING_KANAMORI_HPP
#define _PPSC_HAMILTONIANS_TWO_BAND_HUBBARD_ISING_KANAMORI_HPP

// -----------------------------------------------------------------------
//
// Three band Hubbard Hamiltonian builder
//
// Author: P. Werner, philipp.werner@gmail.com (2016)
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
enum interaction_type { ising, kanamori };

// -----------------------------------------------------------------------
template<class HILB, interaction_type INTERACTION>
class two_band_hubbard_ising_kanamori_scsusc_dense {

public:

  const interaction_type interaction = INTERACTION;

  // ---------------------------------------------------------------------
  two_band_hubbard_ising_kanamori_scsusc_dense(int nt, HILB & hilbert_space) :
    nt(nt),
    mu(0.0), t0pair_seed(0.0), tppair_seed(0.0),
    txpair_seed(0.0), typair_seed(0.0), tzpair_seed(0.0),
    // -- parameter vectors
    U(nt+2, 0.0), J(nt+2, 0.0), delta(nt+2, 0.0), Bz(nt+2, 0.0),
    dU1(nt+2, 0.0), dU2(nt+2, 0.0), //dU3(nt+2, 0.0),
    dE1(nt+2, 0.0), dE2(nt+2, 0.0), //dE3(nt+2, 0.0),
    t0pair(nt+2, 0.0), tppair(nt+2, 0.0),
    txpair(nt+2, 0.0), typair(nt+2, 0.0), tzpair(nt+2, 0.0),
    // -- expectation value vectors
    Q_exp(nt+2), Eint_exp(nt+2),
    n1_exp(nt+2), n2_exp(nt+2), //n3_exp(nt+2),
    m_exp(nt+2),
    n1u_exp(nt+2), n1d_exp(nt+2), n2u_exp(nt+2), n2d_exp(nt+2),
    d1_exp(nt+2), d2_exp(nt+2), //d3_exp(nt+2),
    doublonS0sameorb_exp(nt+2), doublonS0oppositeorb_exp(nt+2), doublonS1_exp(nt+2), triplon_exp(nt+2), quadruplon_exp(nt+2),
    triplet_exp(nt+2),
    tp_exp(nt+2), t0_exp(nt+2), tm_exp(nt+2),
    tz_exp(nt+2), tx_exp(nt+2), ty_exp(nt+2),
    densdens_intraband_exp(nt+2),
    densdens_interband_exp(nt+2),
    densdens_interband_equalspin_exp(nt+2), hil_(hilbert_space)  {

    // -- Construct basic operators

    c1ua = hil_.c_op_[hil_.flavor(0, 0, 0)];
    c1uc = hil_.c_op_[hil_.flavor(0, 0, 1)];
    c1da = hil_.c_op_[hil_.flavor(0, 1, 0)];
    c1dc = hil_.c_op_[hil_.flavor(0, 1, 1)];

    c2ua = hil_.c_op_[hil_.flavor(1, 0, 0)];
    c2uc = hil_.c_op_[hil_.flavor(1, 0, 1)];
    c2da = hil_.c_op_[hil_.flavor(1, 1, 0)];
    c2dc = hil_.c_op_[hil_.flavor(1, 1, 1)];

/*
    c3ua = hil_.c_op_[hil_.flavor(2, 0, 0)];
    c3uc = hil_.c_op_[hil_.flavor(2, 0, 1)];
    c3da = hil_.c_op_[hil_.flavor(2, 1, 0)];
    c3dc = hil_.c_op_[hil_.flavor(2, 1, 1)];
*/
    n1u = c1uc * c1ua;
    n1d = c1dc * c1da;
    n1 = n1u + n1d;

    n2u = c2uc * c2ua;
    n2d = c2dc * c2da;
    n2 = n2u + n2d;
/*
    n3u = c3uc * c3ua;
    n3d = c3dc * c3da;
    n3 = n3u + n3d;
*/
    n = n1 + n2; // + n3;

    m = n1u - n1d + n2u - n2d; // + n3u - n3d; // magnetization

    cdagcdag = c1uc*c2uc;
    cc = c1ua*c2ua;

ucuc = c1uc*c2uc;
uaua = c1ua*c2ua;
dcdc = c1dc*c2dc;
dada = c1da*c2da;

ucdc = c1uc*c2dc;
dcuc = c1dc*c2uc;
uada = c1ua*c2da;
daua = c1da*c2ua;

//triplet = cdagcdag - cc;
triplet = ucdc + dcuc - uada - daua;
tp = ucuc - uaua;
t0 =  1/sqrt(2) * (ucdc + dcuc - uada - daua);
tm = dcdc - dada;

// bracket contains c.c.
tz = (ucuc - uaua) - (dcdc - dada);
tx = -1.*(ucdc - uada) -1.*(dcuc - daua);
ty = (ucuc - uaua) + (dcdc - dada);

    densdens_intraband = n1u * n1d + n2u * n2d; // + n3u * n3d;
    densdens_interband = n1 * n2; // + n1 * n3 + n2 * n3;
    densdens_interband_equalspin =
        n1u * n2u + n1d * n2d;
    //  + n1u * n3u + n1d * n3d
    //  + n2u * n3u + n2d * n3d;

    d1 = n1u * n1d;
    d2 = n2u * n2d;
    //d3 = n3u * n3d;

    triplon = d1*n2 + d2*n1 - 4*d1*d2;
    doublonS0sameorb = d1 + d2 - d1*n2 - d2*n1 + 2*d1*d2; // d1*(1-n2u)*(1-n2d) + d2*(1-n1u)*(1-n1d);
    doublonS0oppositeorb = (n1u-d1)*(n2d-d2) + (n1d-d1)*(n2u-d2); // n1u*(1-n1d)*n2d*(1-n2u) + n1d*(1-n1u)*n2u*(1-n2d);
    doublonS1 = (n1u-d1)*(n2u-d2) + (n1d-d1)*(n2d-d2); // n1u*(1-n1d)*n2u*(1-n2d) + n1d*(1-n1u)*n2d*(1-n2u);
    quadruplon = d1*d2;

    if(interaction == interaction_type::kanamori) {

      pair_hop =
	  c2uc * c2dc * c1ua * c1da
	+ c1dc * c1uc * c2da * c2ua;
//	+ c3uc * c3dc * c2ua * c2da
//	+ c2dc * c2uc * c3da * c3ua
//	+ c1uc * c1dc * c3ua * c3da
//	+ c3dc * c3uc * c1da * c1ua;

      spin_flip =
	  c1dc * c2uc * c2da * c1ua
	+ c1uc * c2dc * c2ua * c1da;
//	+ c2dc * c3uc * c3da * c2ua
//	+ c2uc * c3dc * c3ua * c2da
//	+ c3dc * c1uc * c1da * c3ua
//	+ c3uc * c1dc * c1ua * c3da;
    }

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double d = delta[tstp + 1];
    double u = U[tstp + 1];
    double j = J[tstp + 1];
    double bz = Bz[tstp + 1];

    double du1 = dU1[tstp + 1];
    double du2 = dU2[tstp + 1];
//    double du3 = dU3[tstp + 1];

    double de1 = dE1[tstp + 1];
    double de2 = dE2[tstp + 1];
//    double de3 = dE3[tstp + 1];

    double t0pair_field = t0pair[tstp + 1];
    double tppair_field = tppair[tstp + 1];

    double txpair_field = txpair[tstp + 1];
    double typair_field = typair[tstp + 1];
    double tzpair_field = tzpair[tstp + 1];

//    double half_filling_shift = 0.5*(5*u - 10*j);
    double half_filling_shift = 0.5*(3*u - 5*j);

    Htemp.InitDiagZero(hil_.ssdim_);

    Htemp = Htemp
      + bz * m
//      + 0.5 * d * (n1 - n2 - n3) // "one up, two down" type crystal field splitting
      + 0.5 * d * (n1 - n2) // "crystal field splittin
//      + sosm_seed * d3 // seed for sosm symmetry breaking
//+ pair_seed*( cdagcdag - cc ) // triplet pairing, spin up
//+ pair_field*( cdagcdag - cc )
//
+ tppair_seed*( ucuc - uaua ) // tp triplet pairing
+ t0pair_seed/sqrt(2)*(  ucdc + dcuc - uada - daua ) // t0 triple pairing
+ tppair_field*(  ucuc - uaua )
+ t0pair_field/sqrt(2)*(  ucdc + dcuc - uada - daua )
+ tzpair_seed*tz
+ txpair_seed*tx
+ typair_seed*ty
+ tzpair_field*tz
+ txpair_field*tx
+ typair_field*ty
      - mu * n
      - half_filling_shift * n
      + de1 * n1 + de2 * n2 //+ de3 * n3
      + du1 * n1u * n1d + du2 * n2u * n2d //+ du3 * n3u * n3d
      + (u      ) * densdens_intraband
      + (u - 2*j) * densdens_interband
      + (     -j) * densdens_interband_equalspin;

    if(interaction == interaction_type::kanamori) {
      Htemp = Htemp
	- j * (spin_flip + pair_hop);
    }

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
//    n3_exp[tstp + 1]   = expectation_value(tstp, rho, n3  ).real();

    n1u_exp[tstp + 1]   = expectation_value(tstp, rho, n1u  ).real();
    n1d_exp[tstp + 1]   = expectation_value(tstp, rho, n1d  ).real();
    n2u_exp[tstp + 1]   = expectation_value(tstp, rho, n2u  ).real();
    n2d_exp[tstp + 1]   = expectation_value(tstp, rho, n2d  ).real();


    densdens_intraband_exp[tstp + 1] =
      expectation_value(tstp, rho, densdens_intraband).real();

    densdens_interband_exp[tstp + 1] =
      expectation_value(tstp, rho, densdens_interband).real();

    densdens_interband_equalspin_exp[tstp + 1] =
      expectation_value(tstp, rho, densdens_interband_equalspin).real();

    d1_exp[tstp + 1]   = expectation_value(tstp, rho, d1  ).real();
    d2_exp[tstp + 1]   = expectation_value(tstp, rho, d2  ).real();
//    d3_exp[tstp + 1]   = expectation_value(tstp, rho, d3  ).real();

    triplon_exp[tstp + 1] = expectation_value(tstp, rho, triplon  ).real();
    doublonS0sameorb_exp[tstp + 1] = expectation_value(tstp, rho, doublonS0sameorb  ).real();
    doublonS0oppositeorb_exp[tstp + 1] = expectation_value(tstp, rho, doublonS0oppositeorb  ).real();
    doublonS1_exp[tstp + 1] = expectation_value(tstp, rho, doublonS1  ).real();
    quadruplon_exp[tstp + 1] = expectation_value(tstp, rho, quadruplon  ).real();
    triplet_exp[tstp + 1] = expectation_value(tstp, rho, triplet).real();
    tp_exp[tstp + 1] = expectation_value(tstp, rho, tp).real();
    t0_exp[tstp + 1] = expectation_value(tstp, rho, t0).real();
    tm_exp[tstp + 1] = expectation_value(tstp, rho, tm).real();
    tx_exp[tstp + 1] = expectation_value(tstp, rho, tx).real();
    ty_exp[tstp + 1] = expectation_value(tstp, rho, ty).real();
    tz_exp[tstp + 1] = expectation_value(tstp, rho, tz).real();


    //triplon_exp[tstp + 1] = expectation_value(tstp, rho, d1*n2 + d2*n1 - 4*d1*d2  ).real();
    //doublonS0sameorb_exp[tstp + 1] = expectation_value(tstp, rho, d1*(1-n2u)*(1-n2d) + d2*(1-n1u)*(1-n1d)  ).real();
    //doublonS0oppositeorb_exp[tstp + 1] = expectation_value(tstp, rho, n1u*(1-n1d)*n2d*(1-n2u) + n1d*(1-n1u)*n2u*(1-n2d)  ).real();
    //doublonS1_exp[tstp + 1] = expectation_value(tstp, rho, n1u*(1-n1d)*n2u*(1-n2d) + n1d*(1-n1u)*n2d*(1-n2u)  ).real();
    //quadruplon_exp[tstp + 1] = expectation_value(tstp, rho, d1*d2  ).real();

//    std::cout << "tstp, n, d1, d2, d3, m, Q = "
    std::cout << "tstp " << tstp << " "
	      << "n: " << n1_exp[tstp + 1]+n2_exp[tstp + 1] << " "
              << "d1: " << d1_exp[tstp + 1] << " "
              << "d2: " << d2_exp[tstp + 1] << " "
              << "2S0same: " << doublonS0sameorb_exp[tstp + 1] << " "
              << "2S0opp: " << doublonS0oppositeorb_exp[tstp + 1] << " "
              << "2S1: " << doublonS1_exp[tstp + 1] << " "
              << "3: " << triplon_exp[tstp + 1] << " "
    	      << "4: " << quadruplon_exp[tstp + 1] << " "
	      << "mag: " << m_exp[tstp + 1] << " "
	      << "Q: " << Q_exp[tstp + 1] << " "
	      << "n1: " << n1_exp[tstp + 1] << " "
	      << "n2: " << n2_exp[tstp + 1] << " "
              << "triplet: " << triplet_exp[tstp + 1] << " "
              << "tp: " << tp_exp[tstp + 1] << " "
              << "t0: " << t0_exp[tstp + 1] << " "
              << "tm: " << tm_exp[tstp + 1] << " "
              << "tz: " << tz_exp[tstp + 1] << " "
              << "tx: " << tx_exp[tstp + 1] << " "
              << "ty: " << ty_exp[tstp + 1] << " "
            << std::endl;
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
//    store_real_data_to_hid(group_id, "n3_exp", n3_exp.data(), n3_exp.size());

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
/*
    store_operator_type(sub_group_id, c3uc, std::string("c3uc"));
    store_operator_type(sub_group_id, c3ua, std::string("c3ua"));
    store_operator_type(sub_group_id, c3dc, std::string("c3dc"));
    store_operator_type(sub_group_id, c3da, std::string("c3da"));
*/
    store_operator_type(sub_group_id, n1u, std::string("n1u"));
    store_operator_type(sub_group_id, n1d, std::string("n1d"));

    store_operator_type(sub_group_id, n2u, std::string("n2u"));
    store_operator_type(sub_group_id, n2d, std::string("n2d"));
/*
    store_operator_type(sub_group_id, n3u, std::string("n3u"));
    store_operator_type(sub_group_id, n3d, std::string("n3d"));
*/
    store_operator_type(sub_group_id, n1, std::string("n1"));
    store_operator_type(sub_group_id, n2, std::string("n2"));
//    store_operator_type(sub_group_id, n3, std::string("n3"));

    store_operator_type(sub_group_id, d1, std::string("d1"));
    store_operator_type(sub_group_id, d2, std::string("d2"));
//    store_operator_type(sub_group_id, d3, std::string("d3"));

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

  double mu, tppair_seed, t0pair_seed, tppair_field, t0pair_field;
  double txpair_seed, typair_seed, tzpair_seed, txpair_field, typair_field, tzpair_field;
  std::vector<double> U, J, delta, Bz;
//  std::vector<double> dU1, dU2, dU3, dE1, dE2, dE3;
  std::vector<double> dU1, dU2, dE1, dE2, tppair, t0pair;
  std::vector<double> txpair, typair, tzpair;
  std::vector<double> Q_exp, Eint_exp;
//  std::vector<double> n1_exp, n2_exp, n3_exp, m_exp;
  std::vector<double> n1_exp, n2_exp, m_exp;
  std::vector<double> n1u_exp, n1d_exp, n2u_exp, n2d_exp;
  std::vector<double> densdens_intraband_exp,
		      densdens_interband_exp,
		      densdens_interband_equalspin_exp;
  std::vector<double> d1_exp, d2_exp; //, d3_exp;
  std::vector<double> doublonS0sameorb_exp, doublonS0oppositeorb_exp, doublonS1_exp, triplon_exp, quadruplon_exp;
  std::vector<double> triplet_exp;
  std::vector<double> tp_exp, t0_exp, tm_exp;
  std::vector<double> tz_exp, tx_exp, ty_exp;

  operator_type c1ua, c1uc, c1da, c1dc,
		c2ua, c2uc, c2da, c2dc;
//		c3ua, c3uc, c3da, c3dc;

//  operator_type n1u, n1d, n2u, n2d, n3u, n3d, n1, n2, n3, n, Q, m;
  operator_type n1u, n1d, n2u, n2d, n1, n2, n, Q, m, cdagcdag, cc, ucuc, dcdc, uaua, dada, ucdc, dcuc, uada, daua;

  operator_type d1, d2; // d3;
  operator_type doublonS0sameorb, doublonS0oppositeorb, doublonS1, triplon, quadruplon;
  operator_type triplet;
  operator_type tp, t0, tm;
  operator_type tz, tx, ty;

  operator_type densdens_intraband,
		densdens_interband,
		densdens_interband_equalspin;

  operator_type spin_flip, pair_hop;
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_THREE_BAND_HUBBARD_HPP
