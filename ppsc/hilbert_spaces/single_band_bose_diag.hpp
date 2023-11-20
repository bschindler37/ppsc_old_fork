#ifndef _PPSC_HILBERT_SPACE_SINGLE_BAND_BOSE_DIAG_HPP
#define _PPSC_HILBERT_SPACE_SINGLE_BAND_BOSE_DIAG_HPP

// -----------------------------------------------------------------------
//
// Single band bosonic hilbert space (diagonal basis, no symmetry break)
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "single_band_bose.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hilbert_spaces {
// -----------------------------------------------------------------------
  
// -----------------------------------------------------------------------
class single_band_bose_diag : public single_band_bose {

public:
  
  typedef single_band_bose base_type;

  // ---------------------------------------------------------------------
  void init(int Nmax) {

    Nmax_ = Nmax;
    norb_ = 1;
    spin_degeneracy_ = 1;
    // spin_conservation_ = false;
    xi_ = 1;
    
    //base_type::init(norb, spin_degeneracy, spin_conservation);

    // hilbert space dimension Nmax^{nflavour}
    nflavor = norb_ * spin_degeneracy_;
    nh_ = std::pow(Nmax, nflavor);

    // diagonal version all states have a sector each
    ns_ = nh_;

    // Initialize local class vectors
    init_storage();
    init_states();
    init_c_cdag();
    init_vacuum(); 
  }
  
  // ---------------------------------------------------------------------
  void init_states() {
    int np, i, sector;

    // Loop over all states
    for(state psi = 0; psi < nh_; psi++){

      // count number of each spin flavor in psi, and find the sector
      // sector = (nspin(0)*norb + spin(1))*norb + nspin(2) ....  
      // (spin-conservation)
      // sector = nparticles (no spin-conservation)

      // Compute occupation numbers for each flavour
      size_t Ntot = 0;
      state psi1 = psi;
      std::vector<size_t> state_vec(nflavor);

      for(size_t flavor = 0; flavor < (size_t) nflavor; flavor++) {

	// <n> = numeric_cast<int>(psi / Nmax^flavor) % Nmax
	size_t n = (psi1 / TemplatePow(Nmax_, flavor)) % Nmax_;
	state_vec[flavor] = n;
	Ntot += n;

      }

      np = Ntot; // Number of particles
      sector = psi; // diagonal version has one sector per state

      i = ssdim_[sector];

      npart_[sector] = np;
      nstock_[sector][i] = psi;
      ninv_[psi] = i;
      sector_[psi] = sector;
      ssdim_[sector]++;
      sig_[sector] = 1;
	    
    }

    // each element in nstock_ was original a vector of length nh_
    // resize to final (actually needed) size
    for(sector = 0; sector < ns_; sector++) {
      nstock_[sector].resize(ssdim_[sector]);
    }
  }
  // ---------------------------------------------------------------------
  
private:
  int nflavor;
  
};

// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_SINGLE_BAND_BOSON_DIAG_HPP
// -----------------------------------------------------------------------
