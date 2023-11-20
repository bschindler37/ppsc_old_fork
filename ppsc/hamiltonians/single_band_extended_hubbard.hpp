
#ifndef _PPSC_HAMILTONIANS_SINGLE_BAND_EXTENDED_HUBBARD_HPP
#define _PPSC_HAMILTONIANS_SINGLE_BAND_EXTENDED_HUBBARD_HPP

// -----------------------------------------------------------------------
//
// Single band Hubbard Hamiltonian builder for the use of EDMFT+GW 
// (just additional shifts with respect to the SINGLE_BAND_HUBBARD.HPP )
// Author: D. Golez, denis.golez@gmail.com (2016)
//
// H = U * n_{up} * n_{do} + (epsup+U/2+Delta_sin-mu) * n_{up} + (epsdo+U/2+Delta_sin-mu) * n_{do}
//
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#include "ppsc/ppsc.hpp"
// -----------------------------------------------------------------------

#include "ppsc/hamiltonians/single_band_hubbard.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hamiltonians {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class HILB> class single_band_extended_hubbard :
    public single_band_hubbard<HILB> {

public:

  typedef single_band_hubbard<HILB> base_type;

  // ---------------------------------------------------------------------
  single_band_extended_hubbard(int nt, HILB & hilbert_space) :
    base_type::single_band_hubbard(nt,hilbert_space), Delta_sin(nt+2) {}

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp,double nexp=1.0) {

    double epsup1 = base_type::eps_up[tstp+1] - base_type::mu + base_type::U[tstp+1]*nexp/2  + Delta_sin[tstp+1];
    double epsdo1 = base_type::eps_do[tstp+1] - base_type::mu + base_type::U[tstp+1]*nexp/2  + Delta_sin[tstp+1];
    Htemp = base_type::U[tstp+1] * base_type::docc + epsup1 * base_type::nu + epsdo1 * base_type::nd;
    std::cout << "Set Hamiltonian " << tstp << " " << epsup1 << " " <<epsdo1 << " " << base_type::mu << " " << Delta_sin[tstp+1] << " " << base_type::U[tstp+1]  << " " << nexp << " " <<  base_type::U[tstp+1]*nexp <<  " " << base_type::eps_up[tstp+1]  << std::endl; 

  }

  std::vector<double> Delta_sin;

};

// -----------------------------------------------------------------------
template<class HILB> class single_band_extended_hubbard_spin :
    public single_band_hubbard<HILB> {

public:

  typedef single_band_hubbard<HILB> base_type;

  // ---------------------------------------------------------------------
  single_band_extended_hubbard_spin(int nt, HILB & hilbert_space) :
    base_type::single_band_hubbard(nt,hilbert_space), Delta_sin_up(nt+2),Delta_sin_do(nt+2) {}

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp,double nd_exp=0.5,double nu_exp=0.5) {

    // double epsup1 = base_type::eps_up[tstp+1] - base_type::mu + base_type::U[tstp+1]*nd_exp  + Delta_sin_up[tstp+1];
    // double epsdo1 = base_type::eps_do[tstp+1] - base_type::mu + base_type::U[tstp+1]*nu_exp  + Delta_sin_do[tstp+1];
    double epsup1 = base_type::eps_up[tstp+1] - base_type::mu + Delta_sin_up[tstp+1];
    double epsdo1 = base_type::eps_do[tstp+1] - base_type::mu + Delta_sin_do[tstp+1];
    Htemp = base_type::U[tstp+1] * base_type::docc + epsup1 * base_type::nu + epsdo1 * base_type::nd;
    std::cout << "Set Hamiltonian " << tstp << " " << epsup1 << " " <<epsdo1 << " " << base_type::mu << " " << Delta_sin_up[tstp+1] << " " << Delta_sin_up[tstp+1] << " " << base_type::U[tstp+1]  << " " << nd_exp << " " <<  base_type::U[tstp+1]*nd_exp << " " << nu_exp << " " <<  base_type::U[tstp+1]*nu_exp <<  " " << base_type::eps_up[tstp+1]  << std::endl;
  }

  std::vector<double> Delta_sin_up,Delta_sin_do;

};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_HAMILTONIANS_SINGLE_BAND_EXTENDED_HUBBARD_HPP
