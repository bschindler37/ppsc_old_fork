
#ifndef _PPSC_HILBERT_SPACE_SINGLE_BAND_FERMI_SPIN_SC_HPP
#define _PPSC_HILBERT_SPACE_SINGLE_BAND_FERMI_SPIN_SC_HPP

// -----------------------------------------------------------------------
//
// Single band fermionic hilbert space with non-colinear spin and pairing
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

// -----------------------------------------------------------------------
class single_band_fermi_spin_sc : public hilbert_space_base {

public:
  void init() {

    norb_ = 1;
    spin_degeneracy_ = 2;
    nh_ = 4;
    ns_ = 2;

    // ----------------------------------------------------------------

    init_storage();

    ssdim_[0] = 2;
    ssdim_[1] = 2;

    for (int sector = 0; sector < ns_; sector++)
      nstock_[sector].resize(ssdim_[sector]);

    // ----------------------------------------------------------------

    int sector;

    sector = 0;
    nstock_[sector][0] = 0;
    ninv_[0] = 0;
    sector_[0] = sector; // 00
    nstock_[sector][1] = 3;
    ninv_[3] = 1;
    sector_[3] = sector; // 11

    sector = 1;
    nstock_[sector][0] = 1;
    ninv_[1] = 0;
    sector_[1] = sector; // 01
    nstock_[sector][1] = 2;
    ninv_[2] = 1;
    sector_[2] = sector; // 10

    // ----------------------------------------------------------------

    init_npart_sig();
    init_c_cdag();
    init_vacuum();
  }

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

  void init_vacuum() {
    // sector 0 is the vacuum
    cdvector v1(1);
    vacuum_.InitZero(ssdim_);
    v1(0) = 1.0;
    vacuum_.M_[0] = v1;
  }
};
  
// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_SINGLE_BAND_FERMI_SPIN_SC_HPP
// -----------------------------------------------------------------------
