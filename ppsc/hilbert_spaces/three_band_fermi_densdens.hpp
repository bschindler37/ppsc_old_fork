
#ifndef _PPSC_HILBERT_SPACE_THREE_BAND_FERMI_DENSDENS_HPP
#define _PPSC_HILBERT_SPACE_THREE_BAND_FERMI_DENSDENS_HPP

// -----------------------------------------------------------------------
//
// Three band fermionic hilbert space enabling density-density interactions
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "hilbert_space_base.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hilbert_spaces {
// -----------------------------------------------------------------------

  /*

    only the Jz term of the Hunds coupling is allowed, 
    such that there are additional symmeties, 
    and each block is one-dimensional

    The following type of terms in the Hamiltonian are allowed:
    *  density-density interaction terms
    *  hybridization diagonal in spin and orbital

  */

// -----------------------------------------------------------------------
class three_band_fermi_densdens : public hilbert_space_base {

public:

  // ------------------------------------------------------------------
  void init(void) {

    norb_ = 3;
    spin_degeneracy_ = 2;
    nh_ = 64;
    ns_ = 64;

    // ----------------------------------------------------------------

    init_storage();    
    
    ssdim_ = std::vector<int>(ns_, 1);
    
    for (int sector = 0; sector < ns_; sector++)
      nstock_[sector].resize(ssdim_[sector]);

    // ----------------------------------------------------------------
    // determine all states by a loop

    for(int state = 0; state < ns_; state++) {
      int sector = state;
      nstock_[sector][0] = state;
      ninv_[state] = 0;
      sector_[state] = sector;
    }
    
    // ----------------------------------------------------------------

    init_npart_sig();
    init_c_cdag();
    init_vacuum();
  }

  // ------------------------------------------------------------------
  void init_npart_sig() {

    // There is an assumption that the number of particles
    // in a sector is constant! npart_.size() == ns_ !!
    // This is not the case for symmetry breaking,
    // i.e. |0> and |up,down> are in the same sector.

    for (int sector = 0; sector < ns_; sector++) {
      // count particles:
      int psi = nstock_[sector][0]; // Nb! using only one state from the sector
      int np;
      for (np = 0; psi; psi >>= 1)
        np += psi & 1; // Count the number of bits set (destructive on psi)

      npart_[sector] = np;
      sig_[sector] = 1 - 2 * (np % 2);
    }
  }  
};

// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_THREE_BAND_FERMI_DENSDENS_HPP
// -----------------------------------------------------------------------
