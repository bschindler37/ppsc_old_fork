
#ifndef _PPSC_HILBERT_SPACE_TWO_BAND_FERMI_DENSDENS_HPP
#define _PPSC_HILBERT_SPACE_TWO_BAND_FERMI_DENSDENS_HPP

// -----------------------------------------------------------------------
//
// Two band fermionic hilbert space enabling density-density interactions
//
// Author: Martin Eckstein (201?)
//         Hugo U. R. Strand, hugo.strand@gmail.com (2016)
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

    Almost identical to HilbertSpaceFermi_twoorbital,
    but there are additional symmeries because only the Jz term
    of the Hunds coupling is allowed, such that there are additional
    symmeties, and each block is one-dimensional

    The following type of terms in the Hamiltonian are allowed:
    *  density-density interaction terms
    *  hybridization diagonal in spin and orbital


    do  up
    B A B A    Ndo,Nup    Sector
    ---------------------------------------------
    0 0 0 0     0,0         0
    ---------------------------------------------
    0 0 0 1                 1
                0,1  ----------------------------
    0 0 1 0                 2
    ---------------------------------------------
    0 1 0 0                 3
                1,0  ----------------------------
    1 0 0 0                 4
    ---------------------------------------------
    0 0 1 1     0,2         5
    ---------------------------------------------
    1 1 0 0     2,0         6
    ---------------------------------------------
    0 1 0 1                 7
                     ----------------------------
    1 0 1 0                 8
                1,1  ----------------------------
    0 1 1 0                 9
                     ----------------------------
    1 0 0 1                10
    ---------------------------------------------
    1 1 1 0                11
                2,1  ----------------------------
    1 1 0 1                12
    ---------------------------------------------
    1 0 1 1                13
                1,2  ----------------------------
    0 1 1 1                14
    ---------------------------------------------
    1 1 1 1     2,2        15
    ---------------------------------------------

  */

// -----------------------------------------------------------------------
class two_band_fermi_densdens : public hilbert_space_base {

public:

  // ------------------------------------------------------------------
  void init(void) {

    norb_ = 2;
    spin_degeneracy_ = 2;
    nh_ = 16;
    ns_ = 16;

    // ----------------------------------------------------------------

    init_storage();    
    
    ssdim_ = std::vector<int>(ns_, 1);
    
    for (int sector = 0; sector < ns_; sector++)
      nstock_[sector].resize(ssdim_[sector]);

    // ----------------------------------------------------------------
    // determine all states by hand

    int sector;
    
    sector = 0;
    nstock_[sector][0] = 0;
    ninv_[0] = 0;
    sector_[0] = sector; // 0000

    sector = 1;
    nstock_[sector][0] = 1;
    ninv_[1] = 0;
    sector_[1] = sector; // 0001

    sector = 2;
    nstock_[sector][0] = 2;
    ninv_[2] = 0;
    sector_[2] = sector; // 0010

    sector = 3;
    nstock_[sector][0] = 4;
    ninv_[4] = 0;
    sector_[4] = sector; // 0100

    sector = 4;
    nstock_[sector][0] = 8;
    ninv_[8] = 0;
    sector_[8] = sector; // 1000

    sector = 5;
    nstock_[sector][0] = 3;
    ninv_[3] = 0;
    sector_[3] = sector; // 0011

    sector = 6;
    nstock_[sector][0] = 12;
    ninv_[12] = 0;
    sector_[12] = sector; // 1100

    sector = 7;
    nstock_[sector][0] = 5;
    ninv_[5] = 0;
    sector_[5] = sector; // 0101

    sector = 8;
    nstock_[sector][0] = 10;
    ninv_[10] = 0;
    sector_[10] = sector; // 1010

    sector = 9;
    nstock_[sector][0] = 6;
    ninv_[6] = 0;
    sector_[6] = sector; // 0110

    sector = 10;
    nstock_[sector][0] = 9;
    ninv_[9] = 0;
    sector_[9] = sector; // 1001

    sector = 11;
    nstock_[sector][0] = 14;
    ninv_[14] = 0;
    sector_[14] = sector; // 1110

    sector = 12;
    nstock_[sector][0] = 13;
    ninv_[13] = 0;
    sector_[13] = sector; // 1101

    sector = 13;
    nstock_[sector][0] = 11;
    ninv_[11] = 0;
    sector_[11] = sector; // 1011

    sector = 14;
    nstock_[sector][0] = 7;
    ninv_[7] = 0;
    sector_[7] = sector; // 0111

    sector = 15;
    nstock_[sector][0] = 15;
    ninv_[15] = 0;
    sector_[15] = sector; // 1111

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
#endif // _PPSC_HILBERT_SPACE_TWO_BAND_FERMI_DENSDENS_HPP
// -----------------------------------------------------------------------
