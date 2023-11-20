
#ifndef _PPSC_HILBERT_SPACE_N_BAND_FERMI_DENSE_HPP
#define _PPSC_HILBERT_SPACE_N_BAND_FERMI_DENSE_HPP

// -----------------------------------------------------------------------
//
// N-band fermionic hilbert space with only even/odd fermion number sectors
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2017)
//
// -----------------------------------------------------------------------

#include <bitset>

#include "hilbert_space_base.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hilbert_spaces {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
class n_band_fermi_dense : public hilbert_space_base {

public:
  void init(int norb) {

    norb_ = norb;
    spin_degeneracy_ = 2;
    nh_ = std::pow(2, norb_*spin_degeneracy_);
    ns_ = 2;

    // ----------------------------------------------------------------

    init_storage();

    ssdim_[0] = nh_/2; // even number of fermions sector
    ssdim_[1] = nh_/2; // odd number of fermions sector

    for (int sector = 0; sector < ns_; sector++)
      nstock_[sector].resize(ssdim_[sector]);

    // ----------------------------------------------------------------

    int even_state_idx = 0;
    int odd_state_idx = 0;

    // iterate over all states, and sort them on
    // even/odd numbers of fermions
    
    for(unsigned long state = 0; state < nh_; state++) {
      int nf = std::bitset<16>(state).count();
      int sector = nf % 2; // even or odd no of fermions
      if( sector == 0 ) {

	nstock_[sector][even_state_idx] = state;
	ninv_[state] = even_state_idx;
	sector_[state] = sector;
	even_state_idx++;
	
      } else { // sector == 1

	nstock_[sector][odd_state_idx] = state;
	ninv_[state] = odd_state_idx;
	sector_[state] = sector;
	odd_state_idx++;
      }
	
      
    }

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
#endif // _PPSC_HILBERT_SPACE_N_BAND_FERMI_DENSE_HPP
// -----------------------------------------------------------------------
