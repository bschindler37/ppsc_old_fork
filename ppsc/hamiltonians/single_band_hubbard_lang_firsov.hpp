// -----------------------------------------------------------------------
//
// Single band Hubbard Hamiltonian builder, with Lang-Firsov U_eff
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2017)
//
// H = (U - 2g^2/omega) * n_{up} * n_{do}
//      + (epsup-U/2-mu) * n_{up} + (epsdo-U/2-mu) * n_{do}
//
// -----------------------------------------------------------------------
#pragma once

#include "single_band_hubbard.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hamiltonians {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class HILB> class single_band_hubbard_lang_firsov :
  public single_band_hubbard<HILB> {

public:

  typedef single_band_hubbard<HILB> B;
  
  // ---------------------------------------------------------------------
  single_band_hubbard_lang_firsov(int nt, HILB & hilbert_space) :
    B::single_band_hubbard(nt, hilbert_space), g(0.), omega(1.) {}

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {

    double U_eff = B::U[tstp+1] - 2*g*g / omega;
    double epsup1 = B::eps_up[tstp+1] - B::mu - U_eff * 0.5;
    double epsdo1 = B::eps_do[tstp+1] - B::mu - U_eff * 0.5;
    Htemp = U_eff * B::docc + epsup1 * B::nu + epsdo1 * B::nd;
  }
  
  // ---------------------------------------------------------------------
  void store(hid_t group_id) {
    store_double_attribute_to_hid(group_id, "g", g);
    store_double_attribute_to_hid(group_id, "omega", omega);
    B::store(group_id);
  }
  
  double g; // Fermion-Boson coupling
  double omega; // Boson frequency
  
};

// -----------------------------------------------------------------------
} // end namespace hamiltonians
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------
