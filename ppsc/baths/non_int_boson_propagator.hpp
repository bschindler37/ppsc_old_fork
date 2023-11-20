
#ifndef _NON_INT_BOSON_PROPAGATOR_HPP
#define _NON_INT_BOSON_PROPAGATOR_HPP

// -----------------------------------------------------------------------
//
// Real-time non-interacting boson <XX> propagator
//
// Author: Denis, Golez (2015)
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace boson_utils {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
double bose(double beta, double omega) {
  const double EXPMAX = 100.0;
  double arg = omega * beta;
  if (fabs(arg) > EXPMAX) {
    return (arg > 0.0 ? 0.0 : -1.0);
  } else {
    return 1.0 / (exp(arg) - 1.0);
  }
}
  //#undef EXPMAX

// -----------------------------------------------------------------------
/*! Auxiliary function for phonon propagator

        bose_exp_1 = Exp(-omega tau)*(1-exp(-omega beta))=-Exp(-omega
   tau)*f_B(-omega)

        @param double beta - inverse temperature
        @param double tau -  imaginary time
        @param double omega - frequency
        @return double - bose_exp_1
*/
double bose_exp_1(double beta, double tau, double omega) {
  if (omega > 0) {
    return -exp(-omega * tau) *
           bose(beta, -omega); // -exp(-w*t)/(1-exp(-b*w)) always OK for w>0
  } else {
    return exp(-(tau - beta) * omega) *
           bose(beta, omega); // exp(-(t-b)*w)/(exp(w*b)-1)
  }
}

// -----------------------------------------------------------------------
/*! Auxiliary function for phonon propagator

bose_exp_2 = Exp(omega tau)*(1-exp(-omega beta))=-Exp(-omega tau)*f_B(-omega)

@param T beta - inverse temperature
@param T tau -  imaginary time
@param T omega - frequency
@return double - bose_exp_2
*/

double bose_exp_2(double beta, double tau, double omega) {
  return -exp(omega * tau) * bose(beta, -omega);
}

// -----------------------------------------------------------------------
/*! Noninteracting phononic propagator  <<X,X>>

@param T beta -  inverse temperature
@param herm_matrix<T> &G - reference to set the noninteracting phononic
propagator
*/

template <int SIZE>
void green_from_eps_phonon_dispatch(double beta, gf_type & G,
                                    double omega, double h) {

  typedef std::complex<double> CPLX;

  int nt = G.nt(), ntau = G.ntau(), size1 = G.size1(), i, i0, n, m;

  double eps1, fB, fBm, dtau, tau;
  std::complex<double> x;
  std::vector<std::complex<double>> expp;
  G.clear();
  expp.resize(nt + 1);

  for (i = 0; i < size1; i++) {
    i0 = i * size1 + i;
    for (m = 0; m <= nt; m++)
      expp[m] = std::complex<double>(cos(omega * m * h), -sin(omega * m * h));
    fB = bose(beta, omega);
    fBm = -bose(beta, -omega);
    dtau = beta / ntau;
    for (m = 0; m <= ntau; m++) {
      tau = m * dtau;

      G.matptr(m)[i0] =
          (-1.0) *
          (bose_exp_1(beta, tau, omega) +
           (-1.0) *
               bose_exp_1(beta, tau,
                          -omega)); // Minus due to conversion from bose_exp_1;
      for (n = 0; n <= nt; n++) {
        G.tvptr(n, m)[i0] = std::complex<double>(0, -1.0) *
                            (bose_exp_1(beta, tau, omega) * std::conj(expp[n]) +
                             expp[n] * (-1.0) * bose_exp_1(beta, tau, -omega));
      }
    }

    for (m = 0; m <= nt; m++) {
      for (n = 0; n <= m; n++) {
        G.retptr(m, n)[i0] = std::complex<double>(0, 1.0) *
                             (-expp[m - n] + std::conj(expp[m - n]));
        G.lesptr(n, m)[i0] = std::complex<double>(0, -1.0) *
                             (CPLX(fB, 0.0) * std::conj(expp[m - n]) +
                              CPLX(fBm, 0.0) * expp[m - n]);
      }
    }
  }
}

// -----------------------------------------------------------------------
void green_from_eps_phonon(double beta, gf_type & D0,
                           double omega, double h) {
  if (D0.size1() == 1)
    green_from_eps_phonon_dispatch<1>(beta, D0, omega, h);
  else
    green_from_eps_phonon_dispatch<LARGESIZE>(beta, D0, omega, h);
}

// -----------------------------------------------------------------------
} // namespace boson_utils
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------
  
#endif // NON_INT_BOSON_PROPAGATOR_HPP

